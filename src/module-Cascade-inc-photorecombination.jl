
# Functions and methods for scheme::Cascade.RadiativeRecombinationScheme computations


"""
`Cascade.computeContinuumOrbitals(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, level::ManyElectron.Level)` 
    ... computes in turn all the necessary continuum orbitals for the given energy grid and maxKappa value; a set of orbitals
        cOrbitals::Dict{Subshell, Orbital} is returned.
"""
function computeContinuumOrbitals(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, level::ManyElectron.Level)
    cOrbitals = Dict{Subshell, Orbital}()
    # Generate potential for continuum orbitals for this step
    maxFreeElectronEnergy_au  = Defaults.convertUnits("energy: to atomic", scheme.maxFreeElectronEnergy)
    enGrid       = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.0, maxFreeElectronEnergy_au, scheme.NoFreeElectronEnergies, printout=false)
    maxKappa     = maximum(scheme.lValues) + 1
    nrContinuum  = Continuum.gridConsistency(maximum(enGrid.t) + 0.1, comp.grid)
    contSettings = Continuum.Settings(false, nrContinuum);   
    npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
    ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
    ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
    wp           = Basics.computePotential(Basics.DFSField(1.0), comp.grid, level)
    pot          = Basics.add(npot, wp)
    #  Generate continuum orbitals
    for  (ie, en)  in  enumerate(enGrid.t)
        for  kappa = -maxKappa:maxKappa     
            if  kappa == 0      continue    end
            sh    = Subshell(100+ie, kappa)
            cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(en, sh, pot, contSettings)
            cOrbitals[sh]           = cOrbital
            println(">> New continum orbital generated for $sh and energy $en ")
        end
    end
    
    return( cOrbitals )
end


"""
`Cascade.computeSteps(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested capture amplitudes as well as PhotoRecombination.Line's, etc. for all pre-specified 
        radiative recombination steps of the cascade. When compared with standard computations of these atomic processes, however, 
        the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.CaptureData  is returned.
"""
function computeSteps(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesRR = PhotoRecombination.Line[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    maxFreeElectronEnergy_au  = Defaults.convertUnits("energy: to atomic", scheme.maxFreeElectronEnergy)
    enGrid       = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.0, maxFreeElectronEnergy_au, scheme.NoFreeElectronEnergies, printout=false)
    nt = 0;   st = 0;   previousMeanEn = 0.
    # First compute all necessary continuum orbitals
    cOrbitals = Cascade.computeContinuumOrbitals(scheme, comp, stepList[1].initialMultiplet.levels[1])
    
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): "
        println(sa);    if  printSummary   println(iostream, sa)   end 
                                                
        if      step.process == Basics.Rec()
            newLines = PhotoRecombination.computeLinesWithContinuumOrbital(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, 
                                                                            comp.grid, cOrbitals, enGrid, step.settings, output=true) 
            append!(linesRR, newLines);    nt = length(linesRR)
        else   error("Unsupported atomic process for cascade computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) decay lines."
        println(sa);    if  printSummary   println(iostream, sa)   end 
    end
    #
    data = Cascade.Data{PhotoRecombination.Line}(linesRR)
end



"""
`Cascade.determineSteps(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this radiative recombination cascade. It considers the 
        radiative recombination between the blocks from the initialList and captureList. It checks that at least on pair of
        levels supports such a `radiative stabilization' within the step. A stepList::Array{Cascade.Step,1} is returned, and 
        for which subsequently all required capture amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1})
    # Determine a (free-electron) energy grid 
    println(" ")
    ## NB: settings.electronEnergies are in USER units -- PhotoRecombination.determineLines converts them with
    ## convertUnits("energy: to atomic", ...).  This grid must therefore be spanned with the RAW
    ## scheme.maxFreeElectronEnergy, not the converted one; the *_au value below is only used for the
    ## continuum orbitals in computeSteps/computeContinuumOrbitals, which work in atomic units.
    maxFreeElectronEnergy_au  = Defaults.convertUnits("energy: to atomic", scheme.maxFreeElectronEnergy)
    enGrid       = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.0, scheme.maxFreeElectronEnergy,
                                 scheme.NoFreeElectronEnergies, printout=true)
    println(">> Energy grid points:  $(enGrid.t[1:3])  ...  $(enGrid.t[end-2:end])")
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  capturedBlock in capturedList
            for  initialBlock in initialList
                if  initialBlock.NoElectrons + 1 != capturedBlock.NoElectrons   error("stop a")     end
                # Check that at least one energy supports radiative stabilization
                maxKappa = maximum(scheme.lValues) + 1
                settings = PhotoRecombination.Settings(scheme.multipoles, [UseCoulomb, UseBabushkin], enGrid.t, 
                                                        Float64[], false, false, false, false, true, maxKappa, LineSelection())
                push!( stepList, Cascade.Step(Basics.Rec(), settings, initialBlock.confs,     capturedBlock.confs,     
                                                                        initialBlock.multiplet, capturedBlock.multiplet ) )
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this radiative recombination cascade, and compute also the 
        corresponding multiplets. The different cascade approaches help realize different strategies how these blocks are 
        selected and computed. A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.RadiativeRecombinationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[];    basis = Basis()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        sa = "\n* Generate blocks for RR plasma rate coefficient computations: \n" *
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
                multiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                basis     = multiplet.levels[1].basis
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
            multiplet = Hamiltonian.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings, printout=false)
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
`Cascade.perform(scheme::RadiativeRecombinationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a radiative-recombination (RR) plasma rate coefficient computation basedm on all selected electron 
        capture steps. Such a computation starts from a given set of initial configurations xor initial multiplets and 
        (1) generates all singly-excited configurations due to the capture of an electron with a given maximum electron energy; 
        (2) selects all steps for radiative recombination due to given parameters of the scheme::RadiativeRecombinationScheme. 
        The results of these RR plasma rate computation are comprised into (output) data::DecayData, while these data are only 
        printed during the generation and nothing is returned.

`Cascade.perform(scheme::RadiativeRecombinationScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::RadiativeRecombinationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
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
    Basics.displayConfigurations(stdout, initialConfigs, details="initial configurations of the RR cascade ")
    capturedConfigs = Basics.generateConfigurations(Basics.ForPhotoRecombination(comp.scheme.intoShells), initialConfigs)
    Basics.displayConfigurations(stdout, capturedConfigs, details="final configurations of the RR cascade ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, initialConfigs)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, capturedConfigs, printout=false)
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
    Cascade.displayBlocks(stdout, wc1, sa="from the initial configurations of the RR cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="from the final configurations of the RR cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="from the initial configurations of the DR cascade ");      
                        Cascade.displayBlocks(iostream, wc2, sa="from the final configurations of the RR cascade ")           end      
    #
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     
    for block in wc1  push!(gMultiplets, block.multiplet)    end
    for block in wc2  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc1, wc2)
    Cascade.displaySteps(stdout, we, sa="radiative recombination ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="radiative recombination ")    end      
    wf   = Cascade.modifySteps(we)
    #
    data = Cascade.computeSteps(scheme, comp, wf)
    if output    
        ## results = Base.merge( results, Dict("name"                              => comp.name) ) 
        ## results = Base.merge( results, Dict("cascade scheme"                    => comp.scheme) ) 
        ## results = Base.merge( results, Dict("initial multiplets:"               => multiplets) )    
        ## results = Base.merge( results, Dict("generated multiplets:"             => gMultiplets) )    
        ## results = Base.merge( results, Dict("photo-recombination line data:"    => data) )
        #
        results["name"]                              = comp.name
        results["cascade scheme"]                    = comp.scheme
        results["initial multiplets:"]               = multiplets  
        results["generated multiplets:"]             = gMultiplets 
        results["photo-recombination line data:"]    = data
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-rr-rate-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        JLD2.@save filename results
    end
    ## return( results )
    return( results )
end
