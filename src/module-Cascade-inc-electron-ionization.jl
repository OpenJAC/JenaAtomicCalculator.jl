
# Functions and methods for scheme::Cascade.ElectronIonizationScheme computations
#
# This scheme implements the INDIRECT contribution to electron-impact ionization, i.e. excitation-autoionization
# (EA): an inner-shell electron is excited above the ionization threshold by electron impact, and the resulting
# level then autoionizes.  It is the mirror image of Cascade.DielectronicRecombinationScheme and is built on the
# same two-step pattern; this file follows that one closely and deliberately.
#
#      DR :  Auger (capture, by detailed balance)  +  Radiative (stabilization)      ... net recombination
#      EA :  ImpactExc (excitation)                +  Auger (autoionization)         ... net ionization
#
# NOT covered here: the DIRECT channel (Cascade.ImpactIonizationScheme, not implemented -- see the note at that
# struct) and the two resonant-electron-capture channels, in which the incident electron is CAPTURED into a
# doubly-excited resonance that then sheds two electrons, sequentially or simultaneously.  The electron-capture
# scheme they need is NOT outstanding -- it has run since 16-Aug-2026; what is missing is the second Auger
# generation and the strengths.  See the note at Cascade.ElectronIonizationScheme.


"""
`Cascade.computeSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})`
    ... computes in turn all the requested transition amplitudes for all pre-specified steps of the cascade, i.e. the
        ImpactExcitation.Line's of an impact excitation, the AutoIonization.Line's of an autoionization or a capture, and
        the PhotoEmission.Line's of a radiative stabilization.  As for the other cascade schemes, the amount of printout
        is reduced and sent to the summary file.  A set of data::Array{Cascade.Data,1} is returned, always of length three
        and always in the order (impact-excitation, autoionization, radiative), so that a caller may index it without
        first asking which channels were requested; an unrequested channel contributes an empty list.
"""
function computeSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesE = ImpactExcitation.Line[];    linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc lines: "
        println(sa);    if  printSummary   println(iostream, sa)   end
        #
        if      step.process == Basics.ImpactExc()
            newLines = ImpactExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel,
                                                            comp.grid, step.settings, output=true, printout=false)
            append!(linesE, newLines);    nt = length(linesE)
        elseif  step.process == Basics.Auger()
            newLines = AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel,
                                                          comp.grid, step.settings, output=true, printout=false)
            append!(linesA, newLines);    nt = length(linesA)
        elseif  step.process == Basics.Radiative()
            newLines = PhotoEmission.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid,
                                                         step.settings, output=true, printout=false)
            append!(linesR, newLines);    nt = length(linesR)
        else    error("Unsupported atomic process for electron-ionization computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
             "to a total of $nt $(string(step.process)) lines."
        println(sa);    if  printSummary   println(iostream, sa)   end
    end
    #
    data = [ Cascade.Data{ImpactExcitation.Line}(linesE), Cascade.Data{AutoIonization.Line}(linesA),
             Cascade.Data{PhotoEmission.Line}(linesR) ]
    return( data )
end


"""
`Cascade.determineSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation,
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1},
                        ionizedList::Array{Cascade.Block,1})`
    ... determines all step::Cascade.Step's of this excitation-autoionization cascade.  It pairs the initial blocks
        with the excited ones for the electron-impact excitation, and the excited blocks with the ionized ones for
        the subsequent autoionization.  An Auger step is only set up if it is energetically possible at all, i.e.
        if the excited block lies above the ionized one; that test is what restricts the cascade to the
        AUTOIONIZING excited levels, which are the only ones that contribute to ionization.
        A stepList::Array{Cascade.Step,1} is returned.
"""
function determineSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation,
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1},
                        ionizedList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        maxKappa = length(scheme.lValues) > 0 ? maximum(scheme.lValues) + 1 : 300
        # (1) the electron-impact excitation of an inner-shell electron
        for  initialBlock in initialList
            for  excitedBlock in excitedList
                if  initialBlock.NoElectrons == excitedBlock.NoElectrons
                    settings = ImpactExcitation.Settings(ImpactExcitation.Settings(),
                                                         electronEnergies = scheme.electronEnergies,
                                                         energyShift      = scheme.electronEnergyShift,
                                                         maxKappa         = maxKappa)
                    push!( stepList, Cascade.Step(Basics.ImpactExc(), settings, initialBlock.confs, excitedBlock.confs,
                                                                                initialBlock.multiplet, excitedBlock.multiplet) )
                end
            end
        end
        # (2) the autoionization of the excited levels; only those that lie above the ionized block can decay
        for  excitedBlock in excitedList
            for  ionizedBlock in ionizedList
                if  excitedBlock.NoElectrons == ionizedBlock.NoElectrons + 1   &&
                    Basics.determineMeanEnergy(excitedBlock.multiplet) - Basics.determineMeanEnergy(ionizedBlock.multiplet) > 0.
                    settings = AutoIonization.Settings(AutoIonization.Settings(), maxKappa=maxKappa)
                    push!( stepList, Cascade.Step(Basics.Auger(), settings, excitedBlock.confs, ionizedBlock.confs,
                                                                            excitedBlock.multiplet, ionizedBlock.multiplet) )
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`
    ... generates all block::Cascade.Block's of this excitation-autoionization cascade and computes the corresponding
        multiplets.  A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        if  printout
            sa = "\n* Generate blocks for electron-ionization (excitation-autoionization) computations: " *
                 "\n\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                 "\n    + orbitals are generated independently for each block in a Dirac-Fock-Slater potential; " *
                 "\n    + all blocks (multiplets) are generated from single-CSF levels and without configuration mixing; " *
                 "\n    + only the Coulomb interaction is considered for the autoionization. \n"
            println(sa)
            if  printSummary   println(iostream, sa)    end
        end
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   " *
                                                 "with $(confa.NoElectrons) electrons ... ")   end
            scfMultiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet    = Hamiltonian.performCIwithFrozenOrbitals([confa], scfMultiplet.levels[1].basis.orbitals, comp.nuclearModel,
                                                                   comp.grid, Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings);
                                                                   printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif  comp.approach == SCA()   error("Not yet implemented.")
    else    error("Unsupported cascade approach.")
    end
    #
    return( blockList )
end


"""
`Cascade.generateConfigurationsForElectronIonization(multiplets::Array{Multiplet,1}, scheme::ElectronIonizationScheme,
                                                     nm::Nuclear.Model)`
    ... generates the three sets of configurations of an excitation-autoionization cascade: the initial ones, the
        (inner-shell) excited ones with the same number of electrons, and the ionized ones with one electron less.
        A Tuple(initialConfList, excitedConfList, ionizedConfList) is returned.
"""
function generateConfigurationsForElectronIonization(multiplets::Array{Multiplet,1}, scheme::ElectronIonizationScheme,
                                                     nm::Nuclear.Model)
    initialConfList = Configuration[]
    for  mp  in  multiplets
        confList = Basics.extractConfigurations(Basics.FromBasis(), mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)   end   end
    end
    # The inner-shell excited configurations; these keep the number of electrons
    excitedConfList = Basics.generateConfigurations(initialConfList, scheme.excitationFromShells, scheme.excitationToShells,
                                                    scheme.NoExcitations)
    excitedConfList = unique(excitedConfList)
    # The ionized configurations, i.e. one electron removed from any shell that is occupied in the initial ones
    shells          = Basics.extractShellList(initialConfList)
    ionizedConfList = Basics.generateConfigurations(Basics.RemoveElectrons(1, shells), initialConfList)
    ionizedConfList = unique(ionizedConfList)
    #
    return( initialConfList, excitedConfList, ionizedConfList )
end


"""
`Cascade.perform(scheme::ElectronIonizationScheme, comp::Cascade.Computation)`
    ... to set up and perform an excitation-autoionization cascade computation for the given initial configurations
        xor initial multiplets.  Nothing is returned unless output::Bool=false.

`Cascade.perform(scheme::ElectronIonizationScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`
    ... to perform the same but to return the complete output in a dictionary; this dictionary can be saved to disk
        and is the input of a subsequent Cascade.Simulation.
"""
function perform(scheme::ElectronIonizationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true,
                 outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Which channels were requested?  Tested FIRST, before any SCF, so that an empty request costs nothing.
    calcImpact   = Basics.ImpactExcAuto()  in  scheme.processes
    calcResonant = Cascade.hasResonantChannel(scheme)
    if  !calcImpact  &&  !calcResonant
        error("Cascade.perform(::ElectronIonizationScheme): the scheme requests no channel at all.  Put " *
              "Basics.ImpactExcAuto() into its processes for impact-excitation with subsequent autoionization, " *
              "and/or ResonantImpactIonization.SequentialAuger() / .SimultaneousAuger() for the resonant capture " *
              "channels.")
    end
    #
    # Perform the SCF and CI computation for the initial-state multiplets if initial configurations are given
    if  comp.initialConfigs != Configuration[]
        multiplet  = SelfConsistent.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        multiplets = [Multiplet("initial states", multiplet.levels)]
    else
        multiplets = comp.initialMultiplets
    end
    Cascade.displayLevels(stdout, multiplets, sa="initial ")
    if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")    end
    #
    # Generate the three sets of configurations and display them
    wa  = Cascade.generateConfigurationsForElectronIonization(multiplets, scheme, comp.nuclearModel)
    wb1 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[1], sa="(initial part of the) electron-ionization ")
    wb2 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[2], sa="(inner-shell excited part of the) electron-ionization ")
    wb3 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[3], sa="(ionized part of the) electron-ionization ")
    #
    # Determine the blocks and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp, wb1)
    wc2 = Cascade.generateBlocks(scheme, comp, wb2, printout=false)
    wc3 = Cascade.generateBlocks(scheme, comp, wb3, printout=false)
    Cascade.displayBlocks(stdout, wc1, sa="from the initial configurations of the electron-ionization cascade ")
    Cascade.displayBlocks(stdout, wc2, sa="from the inner-shell excited configurations of the electron-ionization cascade ")
    Cascade.displayBlocks(stdout, wc3, sa="from the ionized configurations of the electron-ionization cascade ")
    if  printSummary
        Cascade.displayBlocks(iostream, wc1, sa="from the initial configurations of the electron-ionization cascade ")
        Cascade.displayBlocks(iostream, wc2, sa="from the inner-shell excited configurations of the electron-ionization cascade ")
        Cascade.displayBlocks(iostream, wc3, sa="from the ionized configurations of the electron-ionization cascade ")
    end
    #
    gMultiplets = Multiplet[];     for block in wc2   push!(gMultiplets, block.multiplet)    end
    #
    we = Cascade.Step[]
    if  calcImpact      append!(we, Cascade.determineSteps(scheme, comp, wc1, wc2, wc3))     end
    #
    # The resonant channels need one further configuration set -- the doubly-excited resonances, with one electron MORE
    # than the initial ones -- and the (N+1)-electron configurations they radiate into.  The intermediate and final sets
    # are the very same excited and ionized blocks the impact-excitation channel uses, which is what lets both channels
    # live in one scheme.
    wc4 = Cascade.Block[];   wc5 = Cascade.Block[]
    if  calcResonant
        wr  = Cascade.generateConfigurationsForResonantIonization(multiplets, scheme, comp.nuclearModel)
        wb4 = Basics.displayConfigurations(comp.nuclearModel.Z, wr[1], sa="(doubly-excited resonances of the) electron-ionization ")
        wb5 = Basics.displayConfigurations(comp.nuclearModel.Z, wr[2], sa="(radiative-stabilization part of the) electron-ionization ")
        wc4 = Cascade.generateBlocks(scheme, comp, wb4, printout=false)
        wc5 = isempty(wb5) ? Cascade.Block[] : Cascade.generateBlocks(scheme, comp, wb5, printout=false)
        Cascade.displayBlocks(stdout, wc4, sa="from the doubly-excited resonances of the electron-ionization cascade ")
        if  printSummary
            Cascade.displayBlocks(iostream, wc4, sa="from the doubly-excited resonances of the electron-ionization cascade ")
        end
        for block in wc4   push!(gMultiplets, block.multiplet)    end
        append!(we, Cascade.determineResonantSteps(scheme, comp, wc1, wc4, wc2, wc3, wc5;
                                                   withSecondAuger = !calcImpact))
    end
    #
    Cascade.displaySteps(stdout, we, sa="electron-ionization ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="electron-ionization ")    end
    wf   = Cascade.modifySteps(we)
    data = Cascade.computeSteps(scheme, comp, wf)
    #
    if  output
        results = Base.merge( results, Dict("name"                          => comp.name) )
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) )
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )
        results = Base.merge( results, Dict("impact-excitation lines:"      => data[1].lines) )
        results = Base.merge( results, Dict("autoionization lines:"         => data[2].lines) )
        results = Base.merge( results, Dict("radiative lines:"              => data[3].lines) )
        results = Base.merge( results, Dict("cascade data:"                 => data) )
        #
        if  outputToFile
            filename = "zzz-cascade-electron-ionization-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n" *
                                                 "   using JLD2 \n   results = JLD2.load(''$filename'')    ... to load them back." )   end
            Cascade.writeDataFile(filename, results)
        end
    end
    #
    return( results )
end
