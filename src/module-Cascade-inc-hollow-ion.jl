
# Functions and methods for scheme::Cascade.HollowIonScheme computations


"""
`Cascade.computeSteps(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested (decay) transition amplitudes as well as AutoIonization.Line's, etc. for all 
        pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
        processes, however, the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Array{Cascade.Data,1}  is returned.
"""
function computeSteps(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    linesC = ElectronCapture.Line[];    cOrbitals = Dict{Subshell, Orbital}()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): "
        println(sa);    if  printSummary   println(iostream, sa)   end 
                                                
        if      step.process == Basics.ElecCapture() 
            # Compute continuum orbitals independently for all transitions in the given block.
            newLines = DielectronicCapture.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                            step.settings, output=true, printout=false) 
            append!(linesC, newLines);    nt = length(linesC)
        elseif  step.process == Basics.Auger()  &&   comp.approach == Cascade.AverageSCA()
            # First determine the `mean' free-electron energy for this Auger block and calculate a common set of continuum orbital
            meanEn = 0.;    NoEn = 0
            for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                if  en > 0.1    meanEn = meanEn + en;    NoEn = NoEn + 1    end
            end
            if  NoEn > 0     meanEn = meanEn/ NoEn     else     meanEn = 0.1    end
            if  abs(meanEn - previousMeanEn) / meanEn  < 0.15   # no new continuum orbitals
                println(">> No new continum orbitals are generated for $(keys(cOrbitals)) and for the energy $meanEn ")
            else
                previousMeanEn = meanEn;    cOrbitals = Dict{Subshell, Orbital}()
                for kappa = -step.settings.maxKappa-1:step.settings.maxKappa        if   kappa == 0     continue    end
                    sh           = Subshell(101, kappa);     nrContinuum = Continuum.gridConsistency(meanEn, comp.grid)
                    contSettings = Continuum.Settings(false, nrContinuum);   
                    npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
                    ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
                    ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
                    ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
                    ## wp           = Basics.compute("radial potential: Dirac-Fock-Slater", comp.grid, step.finalMultiplet.levels[1].basis)
                    ## Basics.computePotentialDFS does not exist.  The working pattern is the one of
                    ## module-Cascade-inc-dielectronic-recombination.jl: dispatch on the self-consistent field of
                    ## the settings, which also honours a non-DFS choice instead of hard-wiring one.
                    ## The same dead call still sits in module-Cascade-inc-stepwise-decay.jl:174 and
                    ## module-PhotoDoubleIonization.jl:356; both are left for their own tasks.
                    wp           = Basics.computePotential(comp.asfSettings.scField, comp.grid, step.finalMultiplet.levels[1])
                    pot          = Basics.add(npot, wp)
                    cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(meanEn, sh, pot, contSettings)
                    cOrbitals[sh] = cOrbital
                end
                println(">> Generate continum orbitals in DFS potential for $(keys(cOrbitals)) and for the energy $meanEn ")
            end
        
            newLines = AutoIonization.computeLinesFromOrbitals(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                                step.settings, cOrbitals, output=true, printout=false) 
            append!(linesA, newLines);    nt = length(linesA)
        elseif  step.process == Basics.Auger() 
            # Compute continuum orbitals independently for all transitions in the given block.
            newLines = AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                            step.settings, output=true, printout=false) 
                append!(linesA, newLines);    nt = length(linesA)
        elseif  step.process == Basics.Radiative()
            newLines = PhotoEmission.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                            step.settings, output=true, printout=false) 
            append!(linesR, newLines);    nt = length(linesR)
        else   error("Unsupported atomic process for cascade computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) decay lines."
        println(sa);    if  printSummary   println(iostream, sa)   end 
    end
    #
    ## Cascade.DecayData no longer exists; the generalised Cascade.Data{T} replaced the per-scheme data types,
    ## exactly as in module-Cascade-inc-stepwise-decay.jl.
    data = [ Cascade.Data{PhotoEmission.Line}(linesR), Cascade.Data{AutoIonization.Line}(linesA) ]
end

"""
`Cascade.determineSteps(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, capturedList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all decay processes 
        of the given scheme and selects all pairs of blocks due to the selected cascade approach. It checks that at least 
        on pair of levels supports a decay within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently 
        all required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, capturedList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        @warn("Transition within a cascade block are not yet considered.")
        for  (ia, blocka) in enumerate(capturedList)
            for  (ib, blockb) in enumerate(capturedList)
                if   ia == ib   continue     end
                for  process  in  comp.scheme.processes
                    #
                    if      process == Basics.Radiative()  
                        if  blocka.NoElectrons == blockb.NoElectrons   &&
                            Basics.determineMeanEnergy(blocka.multiplet) - Basics.determineMeanEnergy(blockb.multiplet) > 0.
                            settings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=scheme.multipoles)
                            push!( stepList, Cascade.Step(process, settings, blocka.confs, blockb.confs, blocka.multiplet, blockb.multiplet) )
                        end
                    #
                    elseif  process == Basics.Auger()  
                        if  blocka.NoElectrons == blockb.NoElectrons + 1   &&
                            Basics.determineMeanEnergy(blocka.multiplet) - Basics.determineMeanEnergy(blockb.multiplet) > 0.
                            settings = AutoIonization.Settings(AutoIonization.Settings(); maxKappa=4)
                            push!( stepList, Cascade.Step(process, settings, blocka.confs, blockb.confs, blocka.multiplet, blockb.multiplet) )
                        end
                    end
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this electron-capture and subsequent stabilization (DR) cascade, 
        and compute also the corresponding multiplets. The different cascade approches realized different strategies how these block are 
        selected and computed. A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.HollowIonScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        sa = "\n* Generate blocks for hollow-ion cascade computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; " *
                "\n    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; " *
                "\n    + only the Coulomb interaction is considered for the electron capture. " *
                "\n    + only E1 excitations are considered for the stabilization. \n"
        if  printout       println(sa)              end
        if  printSummary   println(iostream, sa)    end
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
                ## SelfConsistent.performSCF returns a Multiplet, not a Basis; the orbitals sit one level deeper.
                ## The same stale assumption was found in module-Cascade-inc-impact-excitation.jl.
                scfMultiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                multiplet    = Hamiltonian.performCIwithFrozenOrbitals([confa], scfMultiplet.levels[1].basis.orbitals, comp.nuclearModel,
                                                                       comp.grid, Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings);
                                                                       printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()  
        sa = "\n* Generate blocks for hollow-ion cascade computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + each single configuration forms an individual cascade block; " *
                "\n    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; " *
                "\n    + configuration mixing is included for each block, based on H^(DC); " *
                "\n    + all requested multipoles are considered for the stabilization. \n"
        if  printout       println(sa)              end
        if  printSummary   println(iostream, sa)    end
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
            multiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    else  error("Unsupported cascade approach.")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationsForHollowIons(initialConfigs::Array{Configurations,1}, intoShells::Array{Shell,1}, 
                                                decayShells::Array{Shell,1}, noElectrons::Int64)`  
    ... generates all possible configurations as obtained by the capture of noElectrons into the intoShells. 
        A confList::Array{Configuration,1} is returned.
"""
function generateConfigurationsForHollowIons(initialConfigs::Array{Configuration,1}, intoShells::Array{Shell,1}, 
                                                decayShells::Array{Shell,1}, noElectrons::Int64)
    # Generate all configurations with additional noElectrons in the intoShells
    newConfigs = copy(initialConfigs)
    for  ne = 1:noElectrons
        newConfigs = Basics.generateConfigurationsForExcitationScheme(newConfigs, Basics.ExciteByCapture(), Shell[], Shell[], intoShells, 0)
    end
    #
    # Build configurations with all decayShells 'in between', even if zero occupation
    configs = newConfigs;    newConfigs = Configuration[]
    for  conf  in  configs  
        nshells = copy(conf.shells)
        for shell in decayShells    if  haskey(nshells, shell)   else    nshells[shell] = 0    end   end
        push!(newConfigs, Configuration(nshells, conf.NoElectrons))
    end
    #
    # Generate the configurations that are REACHED BY DECAY, i.e. with one or several electrons moved from the
    # captured (into) shells DOWN into the decayShells.  Without this step the decayShells were only ever
    # emptied by RemoveElectrons below, never filled, so `decayShells` did not mean what its name says: for a
    # hollow ion started with an empty K shell no configuration with an occupied 1s was generated at all, and
    # neither the dominant 2p -> 1s radiative channel nor the Auger decay into the ground state of the next
    # ion could appear.  The transfer is iterated, since a cascade moves the electrons down one at a time, and
    # it terminates because each transfer strictly reduces the occupation of the into-shells.
    ## NB the electron number has to be enforced explicitly.  intoShells and decayShells generally OVERLAP --
    ## for a K-empty hollow ion one naturally writes intoShells = [2s,2p], decayShells = [1s,2s,2p] -- and
    ## Basics.generateConfigurations then also returns configurations with an electron ADDED rather than moved
    ## (e.g. 1s^1 2s^2 or 2p^3 from a two-electron parent).  A radiative transition conserves the electron
    ## number, so only the configurations that keep it are decay targets.
    ## Count the electrons from the OCCUPATIONS, not from Configuration.NoElectrons: that field is not kept in
    ## step with the shells by the generators used here, which is why the final loop of this function recomputes
    ## it as well.  Filtering on the field silently lets three-electron configurations through.
    neOf(c)     = sum(values(c.shells))
    neCaptured  = length(newConfigs) > 0 ? neOf(newConfigs[1]) : 0
    transferred = Configuration[];    tConfigs = copy(newConfigs)
    for  nt = 1:noElectrons + length(decayShells)
        tConfigs = Basics.generateConfigurations(tConfigs, intoShells, decayShells, 1)
        tConfigs = [c for c in unique(tConfigs)  if neOf(c) == neCaptured]
        tConfigs = setdiff(tConfigs, unique(vcat(newConfigs, transferred)))
        if  length(tConfigs) == 0    break    end
        append!(transferred, tConfigs)
    end
    append!(newConfigs, transferred);    newConfigs = unique(newConfigs)
    #
    # Now add all decay configurations
    decayConfigs = Configuration[];    dConfigs = copy(newConfigs)
    further = true
    while  further
        dConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(2, decayShells), dConfigs)
        if length(dConfigs) > 0     append!(decayConfigs, dConfigs)     else   further = false      end
    end
    decayConfigs = unique(decayConfigs)
    #
    dConfigs = copy(newConfigs);   append!(dConfigs, decayConfigs)
    further = true
    while  further
        dConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), dConfigs)
        if length(dConfigs) > 0     append!(decayConfigs, dConfigs)     else   further = false      end
    end
    decayConfigs = unique(decayConfigs)
    
    dConfigs = copy(newConfigs);   append!(dConfigs, decayConfigs)
    dConfigs = unique(dConfigs)
    
    #
    # Remove obsolete shells and double configurations
    nconfList = Configuration[];   nshells = Dict{Shell,Int64}();    ne = 0
    for  conf  in  dConfigs  
        nshells = Dict{Shell,Int64}();      ne = 0
        for  (sh,occ) in conf.shells   if  occ > 0     nshells[sh] = occ;  ne = ne + occ    end     end
        ## Skip the fully-stripped ion: repeated RemoveElectrons eventually empties the configuration, and a
        ## bare nucleus is not a cascade block -- it has no levels.  Left in, it reaches Hamiltonian.performCI
        ## and aborts there in Basics.merge with a bare error("stop a").
        if  ne == 0    continue    end
        push!(nconfList, Configuration( nshells, ne))
    end
    #
    # Discard configurations with energies higher than initial configurations after electron capture;
    # this may be required for a very large number of configurations (not yet)

    return( nconfList )
end


#==
"""
`Cascade.generateConfigurationsForDielectronicCapture(multiplets::Array{Multiplet,1},  scheme::HollowIonScheme, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)`  
    ... generates all possible doubly-excited configurations due to (dielectronic) electron capture into the given multiplets.
        The number and type of such doubly-generated configurations depend on (1) the maximum (electron) energy for capturing an electron
        that is closely related to the (maximum) temperature of the plasma; (2) the fromShells from which (and how many displacements)
        are accepted as well as (3) the maximum principle and orbital angular quantum number of the additional (to-) shellsthe fromShells
        into which electrons excited and/or captured. A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) 
        is returned.
"""
function generateConfigurationsForDielectronicCapture(multiplets::Array{Multiplet,1},  scheme::HollowIonScheme, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)
    # Determine all (reference) configurations from multiplets and generate the 'excited' configurations due to the specificed excitations
    initialConfList = Configuration[]
    for mp  in  multiplets   
        confList = Basics.extractNonrelativisticConfigurations(mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
    end
    captureConfList = copy(initialConfList)
    for nc = 1:scheme.NoCapturedElectrons
        captureConfList = Basics.generateConfigurationsWithElectronCapture(captureConfList, Shell[], scheme.intoShells, 0)
    end
    shellList       = Basics.extractNonrelativisticShellList(multiplets)
    shellList       = Basics.merge(scheme.intoShells, scheme.decayShells)
    blockConfList   = copy(initialConfList)
    for nc = 1:scheme.NoCapturedElectrons
        blockConfList   = Basics.generateConfigurationsWithElectronCapture(blockConfList, Shell[], shellList, 0)
    end
    @warn("blockConfList does not yet support the autoionization of the high-n electrons.")
    #
    # Determine first a hydrogenic spectrum for all subshells of the initial and doubly-excited states
    allConfList   = Configuration[];      append!(allConfList, initialConfList);      append!(allConfList, blockConfList)
    allSubshells  = Basics.extractRelativisticSubshellList(allConfList)
    primitives    = Bsplines.generatePrimitives(grid)
    ## The argument order had drifted here as well: Bsplines.generateOrbitalsHydrogenic takes
    ## (subshells, nuclearModel, primitives).  The same defect was found and fixed in
    ## module-Cascade-inc-expansion-opacity.jl; this call sits in a path that example-Fj.jl does not reach,
    ## so it would have raised a MethodError only for a hollow-ion cascade with high-n spectator shells.
    orbitals      = Bsplines.generateOrbitalsHydrogenic(allSubshells, nm, primitives, printout=true)
    # Exclude configurations with too high mean energies
    en            = Float64[];   
    for conf in initialConfList
        wen = Basics.determineMeanEnergy(conf, orbitals, nm, grid)
        push!(en, wen)
    end
    initialMean = sum(en) / length(en)
    println(">>> initial configuration(s) have mena energies  $initialMean  [a.u.].")

    return( (initialConfList, blockConfList)  )
end  ==#


"""
`Cascade.perform(scheme::HollowIonScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a dielectronic-recombination (DR) plasma rate coefficient computation that combines the electron capture
        and the subsequent radiative stabilization steps. Such a computation starts from a given set of initial configurations xor 
        initial multiplets and (1) generates all doubly-excited configurations due to the capture of an electron with a given maximum
        electron energy; (2) selects all electron capture (inverse Auger) and re-autoionization (Auger) steps and (3) selects
        all steps for radiative stabilization due to given parameters of the scheme::HollowIonScheme. The results of 
        these DR plasma rate computation are comprised into (output) data::ExcitationData, while these data are only printed during 
        the generation and nothing is returned.

`Cascade.perform(scheme::HollowIonScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::HollowIonScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the SCF and CI computation for the intial-state multiplets if initial configurations are given
    if  comp.initialConfigs != Configuration[]
        multiplet  = SelfConsistent.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        multiplets = [Multiplet("initial states", multiplet.levels)]
    else
        multiplets = comp.initialMultiplets
    end
    # Print out initial configurations and levels 
    Cascade.displayLevels(stdout, multiplets, sa="initial ")
    if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")                            end
    #
    # Generate all configurations with NoCapturedElectrons in the intoShells 
    wa = Cascade.generateConfigurationsForHollowIons(comp.initialConfigs, comp.scheme.intoShells, comp.scheme.decayShells, 
                                                        comp.scheme.NoCapturedElectrons)
    # Display and group all configuration together
    wb = Basics.displayConfigurations(comp.nuclearModel.Z, wa, sa="hollow ion configurations ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb, printout=false)
    Cascade.displayBlocks(stdout, wc, sa="for the hollow ion cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc, sa="for the hollow ion cascade ")    end      
    #
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     for block in wc  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc)
    Cascade.displaySteps(stdout, we, sa="hollow ion decay ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="hollow ion decay ")    end      
    wf   = Cascade.modifySteps(we)
    #
    data = Cascade.computeSteps(scheme, comp, wf)
    if output    
        results = Base.merge( results, Dict("name"                                  => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                        => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"                   => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"                 => gMultiplets) )    
        ## The simulations read their line data from "cascade data:" (Cascade.reviewData); storing it only under
        ## "hollow-ion line data:" meant that NO Cascade.Simulation could consume a hollow-ion cascade at all.
        ## Both keys are written: the descriptive one for direct access, the standard one for the simulations.
        results = Base.merge( results, Dict("hollow-ion line data:" => data) )
        results = Base.merge( results, Dict("cascade data:"         => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-hollow-ion-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        JLD2.@save filename results
    end
    ## return( results )
    return( results )
end
