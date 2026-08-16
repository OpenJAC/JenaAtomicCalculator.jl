
# Functions and methods for scheme::Cascade.DielectronicRecombinationScheme computations


"""
`Cascade.computeSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested capture & transition amplitudes as well as DielectronicCapture.Line's, AutoIonization.Line's, 
        etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
        processes, however, the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.CaptureData  is returned.
"""
function computeSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    cOrbitals = Dict{Subshell, Orbital}();    cPhases = Dict{Subshell, Float64}()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): "
        println(sa);    if  printSummary   println(iostream, sa)   end 
                                                
        if      step.process == Basics.Auger()  &&   comp.approach == Cascade.AverageSCA()
            meanEn = 0.;    NoEn = 0
            for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                if  en > 0.03    meanEn = meanEn + en;    NoEn = NoEn + 1    end
            end
            if  NoEn > 0     meanEn = meanEn/ NoEn     
            else     
                println(">> No transition with positive energy.")
                continue    
            end
        
            # Generate potential for continuum orbitals for this step
            nrContinuum  = Continuum.gridConsistency(meanEn + 0.1, comp.grid)
            contSettings = Continuum.Settings(false, nrContinuum);   
            npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
            ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
            ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
            wp           = Basics.computePotential(comp.asfSettings.scField, comp.grid, step.finalMultiplet.levels[1])         
            pot          = Basics.add(npot, wp)
            #  Generate continuum if not yet available
            generatedKappas = Int64[]
            for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
                en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
                if  en > 0.03    
                    symi      = LevelSymmetry(step.initialMultiplet.levels[p].J, step.initialMultiplet.levels[p].parity)
                    symf      = LevelSymmetry(step.finalMultiplet.levels[q].J, step.finalMultiplet.levels[q].parity)
                    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symf)
                    for  kappa in kappaList
                        sh    = Subshell(101, kappa)
                        if      abs(kappa) > step.settings.maxKappa                               ## continuum orbital not needed
                        elseif  kappa in generatedKappas                                          ## already generated for this step
                        elseif  haskey(cOrbitals,sh)   &&   cOrbitals[sh].energy - en / en < 0.15 ## use previous one
                            println(">> No new continum orbital generated for $sh and energy $en ")
                        else                                                                      ## generate new continuum orbital
                            push!(generatedKappas, kappa)
                            cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(en, sh, pot, contSettings)
                            cOrbitals[sh]           = cOrbital;    cPhases[sh] = phase
                            println(">> New continum orbital generated for $sh and energy $en ")
                        end
                    end
                end
            end
        
            newLines = AutoIonization.computeLinesFromOrbitals(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                                step.settings, cOrbitals, cPhases, output=true, printout=false) 
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
    data = [ Cascade.Data{PhotoEmission.Line}(linesR), Cascade.Data{AutoIonization.Line}(linesA) ]
end

"""
`Cascade.determineSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1}, decayList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this dielectronic cascade. It considers the autoionization
        between the blocks from the captureList and initialList as well as the radiative stabilization between the blocks from the 
        captureList and decayList. It checks that at least on pair of levels supports either an `electron-capture' or 
        `radiative stabilization' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently all 
        required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1}, decayList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  capturedBlock in capturedList
            for  initialBlock in initialList
                if  initialBlock.NoElectrons + 1 != capturedBlock.NoElectrons   error("stop a")     end
                # Check that at least one energy supports autoionization
                settings = AutoIonization.Settings(AutoIonization.Settings(), maxKappa=7)
                push!( stepList, Cascade.Step(Basics.Auger(), settings, capturedBlock.confs,     initialBlock.confs,
                                                                        capturedBlock.multiplet, initialBlock.multiplet) )
            end
        end
        #
        for  capturedBlock in capturedList
            for  decayBlock in decayList
                if  decayBlock.NoElectrons != capturedBlock.NoElectrons   error("stop b")     end
                # Check that at least one energy supports radiative stabilization
                if   Basics.determineMeanEnergy(capturedBlock.multiplet) - Basics.determineMeanEnergy(decayBlock.multiplet) > scheme.minPhotonEnergy
                    ## The multipoles must come from the scheme.  PhotoEmission.Settings() defaults to [E1], so the
                    ## stabilization was always computed in E1 alone -- while generateBlocks() prints "all requested
                    ## multipoles are considered for the stabilization".  scheme.multipoles was read nowhere in this
                    ## file; this is the same defect that example-Fc.jl uncovered for the photo-excitation scheme.
                    settings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=scheme.multipoles,
                                                      gauges=[UseCoulomb, UseBabushkin])
                    push!( stepList, Cascade.Step(Basics.Radiative(), settings, capturedBlock.confs,     decayBlock.confs,
                                                                                capturedBlock.multiplet, decayBlock.multiplet) )
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this electron-capture and subsequent stabilization (DR) cascade, 
        and compute also the corresponding multiplets. The different cascade approches realized different strategies how these blocks are 
        selected and computed. A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.DielectronicRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[];    basis = Basis()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        sa = "\n* Generate blocks for DR plasma rate coefficient computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + orbitals are generated independently for the first block in a Dirac-Fock-Slater potential; " *
                "\n    + these orbitals are re-used for all other block, together with hydrogenic orbitals for the outer part; " *
                "\n    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; " *
                "\n    + only the Coulomb interaction is considered for the electron capture. " *
                "\n    + only E1 excitations are considered for the stabilization. \n"
        if  printout       println(sa)              end
        if  printSummary   println(iostream, sa)    end
        #
        if  length(confs) > 1
            # Determine a list of hydrogenic orbitals for later use 
            relconfList = ConfigurationR[]
            for  confa in confs
                wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confa)
                append!( relconfList, wa)
            end
            subshellList = Basics.generateSubshellList(relconfList)
            Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)
            wa                 = Bsplines.generatePrimitives(comp.grid)
            hydrogenicOrbitals = Bsplines.generateOrbitalsHydrogenic(subshellList, comp.nuclearModel, wa; printout=printout)
        end
        
        for  (ia, confa)  in  enumerate(confs)
            sa = "  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... "
            print(sa);      if  printSummary   println(iostream, sa)   end
            # Now distinguish between the first and all other blocks; for the first block, a SCF is generated and the occupied orbital
            # used also for all other blocks. In addition, a set of hydrogenic orbitals generated for later use
            if  ia == 1
                multiplet     = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                basis         = multiplet.levels[1].basis
            else
                # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
                relconfList  = ConfigurationR[]
                wa           = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confa)
                append!( relconfList, wa)
                subshellList = Basics.generateSubshellList(relconfList)
                Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
                # Generate the relativistic CSF's for the given subshell list
                csfList = CsfR[]
                for  relconf in relconfList
                    newCsfs = Basics.generateCsfRs(relconf, subshellList)
                    append!( csfList, newCsfs)
                end
                # Determine the list of coreSubshells
                coreSubshellList = Subshell[]
                for  k in 1:length(subshellList)
                    mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
                    for  csf in csfList
                        if  csf.occupation[k] != mocc    is_filled = false;    break   end
                    end
                    if   is_filled    push!( coreSubshellList, subshellList[k])    end
                end
                # Add all missing orbitals as hydrogenic
                orbitals      = copy(basis.orbitals)
                for  subsh in subshellList
                    if haskey(orbitals, subsh)   ## do nothing
                    else      orbitals[subsh] = hydrogenicOrbitals[subsh];   print("hydrogenic $subsh ...")
                    end
                end
                
                basis         = Basis(true, confa.NoElectrons, subshellList, csfList, coreSubshellList, orbitals)
            end
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa],  basis.orbitals, comp.nuclearModel, comp.grid, Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings); printout=false)
            
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()
        sa = "\n* Generate blocks for electron-capture & decay cascade computations: \n" *
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
`Cascade.generateCaptureConfigurations(multiplets::Array{Multiplet,1},  coreConfList::Array{Configuration,1}, 
                                       scheme::DielectronicRecombinationScheme, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... checks and compiles are doubly-excited configurations due to an (dielectronic) electron capture into the given multiplets,
        for which the total energy is expected to be in the interval E_initial < E_Captured < E_initial + maxExcitationEnergy.
        A ConfList::Array{Configuration,1} is returned that contains all the doubly-excited configurations. Moreover, a neat table 
        of all configurations is printed to check and control for the correct input.
            
        In practice, the generation of useful doubly-excited configurations is not simple as long as no SCF computations
        are performed explicitly. We here apply the following strategy to keep the computations feasible:
        (1) Compute a mean-field basis for the coreConfList to obtain a proper set of orbitals.
        (2) Generate a set of hydrogenic orbitals for all (n,l) shells with n <= nMax and l <= lMax; the orbitals
            are generated for the nuclear charge Zeff = Z - NoElectrons(core) + 1
        (3) Set-up and diagonalize a multiplet for each doubly-excited configuration; determine eLowest, eHighest
            as the lowest and highest level energies.
        (4) Take the energy conditions in the sloppy form:
                eMin - 1. <= eLowest <= eMax + 1. || eMin - 1. <= eHighest <= eMax + 1.
            to make sure that no relevant configurations are missing
        (5) Report about energies   eMin, eLowest, eHighest eMax  
            to modify the energy condition if appropriate        
"""
function generateCaptureConfigurations(multiplets::Array{Multiplet,1},  coreConfList::Array{Configuration,1}, 
                                       scheme::DielectronicRecombinationScheme, nm::Nuclear.Model, grid::Radial.Grid)
    # Run through all configurations with an additional electron and check whether the total energy is within the given range
    captureConfList = Configuration[];    allSubshells = Subshell[] 
    nMax = scheme.maxIntoShell.n;   lMax = scheme.maxIntoShell.l;    eMin = 0.
    for  mp  in  multiplets    if  mp.levels[1].energy < eMin    eMin = mp.levels[1].energy   end   end 
    for  n = 1:nMax,   l = 0:lMax     # Determine all required subshells
        if  lMax > nMax - 1    continue   end
        subshells = Basics.shellSplitIntoSubshells(Shell(n,l))
        for subsh in subshells   push!(allSubshells, subsh)   end
    end

    eMax = eMin + scheme.maxExcitationEnergy
    
    # Compute a mean-field basis for coreConfList
    println(">> Generate mean-field multiplet for $(length(coreConfList)) excited core configurations \n $coreConfList \n ...  ")
    asfSettings = AsfSettings(AsfSettings(), scField=Basics.DFSField(1.0))
    mp          = SelfConsistent.performSCF(coreConfList, nm, grid, asfSettings, printout=true)
    println(" ... mean-field multiplet done")
    
    # Generate hydrogenic spectrum for nuclear charge Zx 
    nmx = Nuclear.Model(nm; Z = nm.Z - coreConfList[1].NoElectrons + 1.)
    print(">> Generate hydrogenic spectrum for nuclear Z_eff = $(nmx.Z), nMax = $nMax and lMax = $lMax ...  ")
    Defaults.setDefaults("standard grid", grid)
    primitives          = Bsplines.generatePrimitives(grid)  
    hydrogenicOrbitals  = Bsplines.generateOrbitalsHydrogenic(allSubshells, nm, primitives; printout=false)
    
    println(" ... hydrogenic spectrum done")
    
    # Generate doubly-excited configurations, calculate the multiplet with given basis and compare energies
    # for the applied energy condition(s); report about energies; accept or refuse additional configuration
    println(">> Loop through all -- maximally $(nMax*(lMax+1)*length(coreConfList)) -- doubly-excited configurations " *
            "and compare energies ...")
    nCount = 0
    for  n = 1:nMax,   l = 0:lMax 
        if  l > n - 1    continue   end
        newConfList = Basics.generateConfigurations(Basics.AddElectrons(1, [Shell(n,l)]), coreConfList)
        for  conf in newConfList
            orbitals  = copy(mp.levels[1].basis.orbitals)
            subshells = Basics.extractSubshellList(conf, orbitals)
            for  subsh in subshells   orbitals[subsh] = hydrogenicOrbitals[subsh]   end
            asfSettings = AsfSettings(AsfSettings(), scField = Basics.DFSField(1.0), 
                                      startScfFrom = ManyElectron.StartFromPrevious(orbitals))
            if true  
                # Accept all configurations
                eLowest   = eMin;     eHighest = eMin;    accepted = false
            elseif false
                basis   = SelfConsistent.initializeBasis([conf], nm, primitives, asfSettings)
                mp      = Hamiltonian.performCI(basis, nm, grid, AsfSettings(), printout=false)
                eLowest   = mp.levels[1].energy;     eHighest = mp.levels[end].energy;    accepted = false
            elseif true
                # Screening step only: the multiplet is used solely for its lowest/highest level energy, to decide
                # whether this capture configuration falls into the requested window. Single-CSF levels are
                # sufficient here (and were the behaviour before commit 7cc164b), so the cheap representation is
                # requested explicitly rather than inherited -- this is deliberately independent of comp.approach,
                # which is not available in this function anyway.
                mp      = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, nm, grid,
                                                                  AsfSettings(asfSettings; eeInteractionCI=DiagonalCoulomb()), printout=false)
                eLowest   = mp.levels[1].energy;     eHighest = mp.levels[end].energy;    accepted = false
            end
            if  eMin - 4.0 <= eLowest <=  eMax + 4.0   ||   eMin - 4.0 <= eHighest <=  eMax + 4.0
                push!(captureConfList, conf);   accepted = true
            end
            nCount = nCount + 1
            println("  $nCount)  $conf .. $accepted   with  eMin = $eMin  [eL, eH] = [$eLowest, $eHighest]  eMax = $eMax ")
        end
    end
    println("   ... loop through all doubly-excited configurations done")

    return( captureConfList )
end





"""
`Cascade.perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a dielectronic-recombination (DR) plasma rate coefficient computation that combines the electron capture
        and the subsequent radiative stabilization steps. Such a computation starts from a given set of initial configurations xor 
        initial multiplets and (1) generates all doubly-excited configurations due to the capture of an electron with a given maximum
        electron energy; (2) selects all electron capture (inverse Auger) and re-autoionization (Auger) steps and (3) selects
        all steps for radiative stabilization due to given parameters of the scheme::DielectronicRecombinationScheme. The results of 
        these DR plasma rate computation are comprised into (output) data::ExcitationData, while these data are only printed during 
        the generation and nothing is returned.

`Cascade.perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::DielectronicRecombinationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
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
    # Generate subsequent cascade configurations as well as display and group them together
    
    
    initialConfigs  = Basics.extractConfigurations(Basics.FromMultiplet(), multiplets)
    Basics.displayConfigurations(stdout, initialConfigs, details="initial configurations of the DR cascade ")
    theme           = Basics.ForDielectronicRecombination(comp.scheme.excitationFromShells, comp.scheme.excitationToShells,
                                                          comp.scheme.intoShells, comp.scheme.decayShells)
    capturedConfigs, decayConfigs = Basics.generateConfigurations(theme, initialConfigs)
    Basics.displayConfigurations(stdout, capturedConfigs, details="doubly-excited capture configurations of the DR cascade ")
    Basics.displayConfigurations(stdout, decayConfigs, details="decay configurations of the DR cascade ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, initialConfigs)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, capturedConfigs, printout=false)
    wc3 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, decayConfigs,    printout=false)
    # Shift the initial level energy by -electronEnergyShift
    if  scheme.electronEnergyShift != 0. 
        wc1old = wc1;   wc1 = Cascade.Block[]
        for  block in wc1old
            newMultiplet = Basics.shiftTotalEnergies(multiplet, Defaults.convertUnits("energy: to atomic", -scheme.electronEnergyShift))
            push!(wc1, Cascade.Block(block.NoElectrons, block.confs, block.hasMultiplet, newMultiplet))
        end
        println(">> Shift all initial level energies by $(-scheme.electronEnergyShift) $(Defaults.getDefaults("unit: energy"))")
    end
    #
    Cascade.displayBlocks(stdout, wc1, sa="from the initial configurations of the DR cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="from the doubly-excited capture configurations of the DR cascade ")
    Cascade.displayBlocks(stdout, wc3, sa="from the decay configurations of the DR cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="from the initial configurations of the DR cascade ");      
                        Cascade.displayBlocks(iostream, wc2, sa="from the doubly-excited capture configurations of the DR cascade ")
                        Cascade.displayBlocks(iostream, wc3, sa="from the decay configurations of the DR cascade ")            end      
    #
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     
    for block in wc1  push!(gMultiplets, block.multiplet)    end
    for block in wc2  push!(gMultiplets, block.multiplet)    end
    for block in wc3  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc1, wc2, wc3)
    Cascade.displaySteps(stdout, we, sa="electron capture and stabilization ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="electron capture and stabilization ")    end      
    wf   = Cascade.modifySteps(we)
    #
    data = Cascade.computeSteps(scheme, comp, wf)
    if output    
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )    
        results = Base.merge( results, Dict("cascade data:"                 => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-dr-rate-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        Cascade.writeDataFile(filename, results)
    end
    ## return( results )
    return( results )
end


"""
`Cascade.perform(scheme::Cascade.DielectronicCaptureScheme, comp::Cascade.Computation; output::Bool=false,
                 outputToFile::Bool=true, outputDirectory::String="")`
    ... sets up and performs a dielectronic-CAPTURE cascade: the capture of a free electron into doubly-excited levels, together with the
        autoionization of those levels back into the target -- both into its GROUND state, which is the capture channel itself, and into
        its EXCITED states, which is what turns a capture into a contribution to electron-impact EXCITATION. No radiative stabilization is
        computed; that path is dielectronic recombination and belongs to Cascade.DielectronicRecombinationScheme. A
        results::Dict{String,Any} is returned if output=true, and nothing otherwise.

        THREE GROUPS OF CONFIGURATIONS are built, and the third is what distinguishes this from a plain capture calculation:
          * the initial target, N electrons, from comp.initialConfigs;
          * the EXCITED target, also N electrons, by exciting one electron from scheme.excitationFromShells into
            scheme.excitationToShells.  These are the exit channels of a resonant excitation;
          * the doubly-excited resonances, N+1 electrons, from Basics.ForDielectronicRecombination.
        Auger steps are then generated from every resonance block to every N-electron block, ground and excited alike, so that each
        resonance carries both its capture width (the channel back to the ground state) and its excitation widths. A TOTAL
        impact-excitation cross section still has to be assembled from these: the resonant contribution is the capture strength times the
        branching ratio into the excited exit channel, and it adds to the direct cross section of Cascade.ImpactExcitationScheme.

        WHICH RESONANCES CONTRIBUTE TO EXCITATION, and this is easy to get wrong. A resonance can only autoionize into a target state that
        lies BELOW it, and the 1s 2l nl' series converges to the 1s2l threshold FROM BELOW however large n is taken. Capturing into ever
        higher n therefore never opens the 1s->2l excitation channel. What opens it is a resonance built on the SAME core excitation as the
        capture: 1s 3l 3l' lies above 1s2l and autoionizes into it. In practice scheme.intoShells must reach the same principal quantum
        number as scheme.excitationToShells.

        Measured for He-like carbon, which is the case to reason from. With excitationToShells = 2l and intoShells = 2l, all 16 K-LL
        resonances (229-255 eV) autoionize ONLY to 1s^2: pure dielectronic capture, no excitation. With excitationToShells = {2l,3l} and
        intoShells = 3l, the 507 lines split into 129 back to 1s^2 (268-335 eV) and 378 into 1s2s and 1s2p, the latter ejecting electrons
        of 12-37 eV. The 1s2s and 1s2p ejection energies differ by the 2s-2p separation of the target, which is a free consistency check.

        The capture and excitation lines are AutoIonization.Line's -- a capture is the time reverse of an Auger transition -- and are
        returned under "dielectronic-capture lines:" as well as in a Cascade.Data{AutoIonization.Line}.
"""
function perform(scheme::Cascade.DielectronicCaptureScheme, comp::Cascade.Computation; output::Bool=false,
                 outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if  comp.initialConfigs != Configuration[]
        multiplet  = SelfConsistent.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        multiplets = [Multiplet("initial states", multiplet.levels)]
    else
        multiplets = comp.initialMultiplets
    end
    Cascade.displayLevels(stdout, multiplets, sa="initial ")
    #
    initialConfigs = Basics.extractConfigurations(Basics.FromMultiplet(), multiplets)
    Basics.displayConfigurations(stdout, initialConfigs, details="initial configurations of the capture cascade ")
    # The EXCITED target: the exit channels of a resonant excitation, same electron number as the initial target.
    excitedConfigs = Basics.generateConfigurations(Basics.ExciteElectrons(1, scheme.excitationFromShells,
                                                                         scheme.excitationToShells), initialConfigs)
    excitedConfigs = setdiff(excitedConfigs, initialConfigs)
    Basics.displayConfigurations(stdout, excitedConfigs, details="excited target configurations (excitation exit channels) ")
    # The doubly-excited resonances, one electron more.
    theme = Basics.ForDielectronicRecombination(scheme.excitationFromShells, scheme.excitationToShells,
                                                scheme.intoShells, Shell[])
    capturedConfigs, _ = Basics.generateConfigurations(theme, initialConfigs)
    Basics.displayConfigurations(stdout, capturedConfigs, details="doubly-excited capture configurations ")
    #
    drScheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell(0,0), scheme.maxExcitationEnergy,
                                                        scheme.electronEnergyShift, 0., scheme.NoExcitations,
                                                        scheme.excitationFromShells, scheme.excitationToShells,
                                                        scheme.intoShells, Shell[])
    wc1 = Cascade.generateBlocks(drScheme, comp, initialConfigs)
    wc1x = isempty(excitedConfigs) ? Cascade.Block[] :
           Cascade.generateBlocks(drScheme, comp, excitedConfigs, printout=false)
    wc2 = Cascade.generateBlocks(drScheme, comp, capturedConfigs, printout=false)
    if  scheme.electronEnergyShift != 0.
        wc1old = wc1;   wc1 = Cascade.Block[]
        for  block in wc1old
            newMultiplet = Basics.shiftTotalEnergies(block.multiplet,
                                Defaults.convertUnits("energy: to atomic", -scheme.electronEnergyShift))
            push!(wc1, Cascade.Block(block.NoElectrons, block.confs, block.hasMultiplet, newMultiplet))
        end
    end
    Cascade.displayBlocks(stdout, wc1,  sa="from the initial configurations ")
    Cascade.displayBlocks(stdout, wc1x, sa="from the excited target configurations ")
    Cascade.displayBlocks(stdout, wc2,  sa="from the doubly-excited capture configurations ")
    #
    gMultiplets = Multiplet[]
    for block in wc1   push!(gMultiplets, block.multiplet)    end
    for block in wc1x  push!(gMultiplets, block.multiplet)    end
    for block in wc2   push!(gMultiplets, block.multiplet)    end
    # Every resonance autoionizes into every N-electron block: the ground one gives the capture width, the
    # excited ones the resonant-excitation widths.
    we = Cascade.determineSteps(drScheme, comp, vcat(wc1, wc1x), wc2, Cascade.Block[])
    Cascade.displaySteps(stdout, we, sa="electron capture and re-autoionization ")
    wf   = Cascade.modifySteps(we)
    data = Cascade.computeSteps(drScheme, comp, wf)
    #
    if  output
        linesC = AutoIonization.Line[]
        for  cData in data
            if  typeof(cData) == Cascade.Data{AutoIonization.Line}   append!(linesC, cData.lines)    end
        end
        results = Base.merge( results, Dict("name"                          => comp.name) )
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) )
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )
        results = Base.merge( results, Dict("dielectronic-capture lines:"   => linesC) )
        results = Base.merge( results, Dict("cascade data:"                 => data) )
    end
    #
    if  outputToFile
        filename = "zzz-cascade-dielectronic-capture-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:  JLD2.load(''$filename'')")
        Cascade.writeDataFile(filename, results)
    end

    return( results )
end


"""
`Cascade.resonantExcitationStrengths(lines::Array{AutoIonization.Line,1}; totalPhotonRates::Dict{Int64,Float64}=Dict{Int64,Float64}())`
    ... assembles the RESONANT contribution to electron-impact excitation from the autoionization lines of a dielectronic-capture cascade,
        i.e. from the output of Cascade.perform(::Cascade.DielectronicCaptureScheme). Each doubly-excited resonance m is formed by capture
        from the ground level i and then decays; the part that ends on an EXCITED target level f is an excitation that went through a
        resonance. Its integrated cross section (resonance strength) follows JAC's dielectronic-recombination convention,

            C(i,m) = pi^2 / k^2 * A_a(m->i) * (2J_m+1)/(2J_i+1),      k^2 = 2 E_res
            S(i->m->f) = C(i,m) * A_a(m->f) / (Gamma_a + Gamma_r),    Gamma_a = sum_k A_a(m->k)

        which is the same expression the DR strength uses with the radiative width Gamma_r in the numerator replaced by the Auger width
        into the excited channel. A vector of named tuples (resonance, finalLevel, energy, strength) is returned, with the energy in atomic
        units and the strength in a_0^2 * Hartree.

        Gamma_r DEFAULTS TO ZERO, because a capture cascade computes no radiative rates. Neglecting it OVERESTIMATES every branching ratio,
        and the error grows with the nuclear charge, where radiative stabilization competes with autoionization. Pass totalPhotonRates,
        keyed on the resonance level index, whenever those rates are available.

        The ground level is identified as the target level of LOWEST energy among the lines, so the capture channel is the line that ends
        on it; every other final level is treated as an excitation channel.
"""
function resonantExcitationStrengths(lines::Array{AutoIonization.Line,1};
                                     totalPhotonRates::Dict{Int64,Float64}=Dict{Int64,Float64}())
    isempty(lines)  &&  return( NamedTuple[] )
    groundEnergy = minimum(l.finalLevel.energy for l in lines)

    # Group the lines by the resonance they decay from; a resonance is identified by its level energy.
    byRes = Dict{Float64, Vector{AutoIonization.Line}}()
    for  l in lines    push!( get!(byRes, l.initialLevel.energy, AutoIonization.Line[]), l )    end

    strengths = NamedTuple[]
    for  (enRes, ls) in byRes
        gammaA = sum(l.totalRate for l in ls)
        gammaR = get(totalPhotonRates, ls[1].initialLevel.index, 0.)
        gammaA + gammaR <= 0.  &&  continue
        capture = ls[findmin([abs(l.finalLevel.energy - groundEnergy) for l in ls])[2]]
        eRes    = capture.electronEnergy
        eRes <= 0.  &&  continue
        cFactor = pi*pi / (2 * eRes) * capture.totalRate *
                  ((Basics.twice(capture.initialLevel.J) + 1) / (Basics.twice(capture.finalLevel.J) + 1))
        for  l in ls
            l === capture  &&  continue                      # the capture channel itself is not an excitation
            push!(strengths, (resonance = l.initialLevel, finalLevel = l.finalLevel, energy = eRes,
                              strength = cFactor * l.totalRate / (gammaA + gammaR)))
        end
    end

    return( sort(strengths, by = x -> x.energy) )
end


"""
`Cascade.resonantExcitationRateCoefficient(strengths::Array{<:NamedTuple,1}, temp::Float64)`
    ... converts the resonance strengths of Cascade.resonantExcitationStrengths into a plasma rate coefficient at the electron temperature
        temp [K], in the isolated-resonance (delta-like) approximation and with the same prefactor that
        DielectronicRecombination.computeRateCoefficient uses,

            alpha(T) = 4/sqrt(2 pi) * T^(-3/2) * sum_d E_d * exp(-E_d/T) * S_d

        An alpha::Float64 is returned in cm^3/s. Summing over all entries gives the total resonant contribution; filter the vector by
        finalLevel first to obtain the contribution to one excitation channel.
"""
function resonantExcitationRateCoefficient(strengths::Array{<:NamedTuple,1}, temp::Float64)
    temp_au = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
    alpha   = 0.
    for  s in strengths
        alpha = alpha + 4 / sqrt(2pi) * temp_au^(-3/2) * s.energy * exp(-s.energy/temp_au) * s.strength
    end
    # atomic units -> cm^3/s, exactly as in DielectronicRecombination.computeRateCoefficient
    factor = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 *
             Defaults.convertUnits("rate: from atomic to 1/s", 1.0)

    return( factor * alpha )
end
