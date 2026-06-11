
#== August 2025, the following replacements need to be made and tested properly:
++ module-Cascade-inc-photoionization.jl:    ionConfList = Basics.generateConfigurationsWithElectronLoss(initialConfList, scheme.excitationFromShells)
++ Basics.generateConfigurationsWithElectronLoss(confs::Array{Configuration,1}, fromShells::Array{Shell,1}) 
   ... generates a list of non-relativistic configurations for the given (reference) confs and with one 
       removed (ionized) electron fromShells.
++ See Basics.generateConfigurations(RemoveElectrons(), confs)
==#


# Functions and methods for photoionization (cascade) computations

"""
`Cascade.computeSteps(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested transition amplitudes as well as PhotoExcitation.Line's, etc. for all pre-specified 
        excitation steps of the cascade. When compared with standard excitation computations, however, the amount of output is 
        largely reduced and often just printed into the summary file. A set of  data::Cascade.ExcitationData  is returned.
"""
function computeSteps(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesP = PhotoIonization.Line[]    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0;   previousMeanEn = 0.
    @show scheme.initialLevelSelection
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels) * length(scheme.photonEnergies)
        println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc photoionization lines (without selection rules): ")
        if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                "up to $nc photoionization lines (without selection rules): ")   end 
                                                
        if  step.process == Basics.Photo() 
            newLines = PhotoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                           ## step.settings, scheme.initialLevelSelection, output=true, printout=false) 
                                                           step.settings, LevelSelection(true, indices=[1]), output=true, printout=false) 
            append!(linesP, newLines);    nt = length(linesP)
        else   error("Unsupported atomic process for photoionization computations.")
        end
        println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) photoionization lines." )
        if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                "giving now rise to a total of $nt $(string(step.process)) photoionization lines." )   end      
    end
    #
    return( linesP )
end


"""
`Cascade.determineSteps(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
        decay scheme and selects all pairs of blocks due to the selected processes and cascade approach. It checks that at least on 
        pair of levels supports either a `photo-ionization' or `photo-excitation' within the step. A stepList::Array{Cascade.Step,1} 
        is returned, and for which subsequently all required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1})
    #
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  initialBlock in initialList
            for  ionizedBlock in ionizedList
                if  initialBlock.NoElectrons == ionizedBlock.NoElectrons + 1
                    photonEnergies = scheme.photonEnergies
                    println(">>> Photon energies must still be given in user-selected units: $(photonEnergies)")
                    settings = PhotoIonization.Settings(PhotoIonization.Settings(), multipoles = scheme.multipoles, 
                                                        gauges = [UseCoulomb, UseBabushkin], photonEnergies = photonEnergies, 
                                                        lValues = scheme.lValues)
                    push!( stepList, Cascade.Step(Basics.Photo(), settings, initialBlock.confs, ionizedBlock.confs, 
                                                                            initialBlock.multiplet, ionizedBlock.multiplet) )
                end
            end
        end
            #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this excitation cascade, and compute also the corresponding multiplets.
        The different cascade approaches realizes different strategies how these block are selected and computed. 
        A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.PhotoIonizationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        if  printout
        println("\n* Generate blocks for photoabsorption computations:")
        println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
        println("    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; ")
        println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
        println("    + only E1 excitations are considered. \n")
        if  printSummary   
        println(iostream, "\n* Generate blocks for photoabsorption computations:")
        println(iostream, "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
        println(iostream, "    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; ")
        println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
        println(iostream, "    + only E1 excitations are considered. \n")
        end
        end     # printout
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
            multiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa], multiplet.levels[1].basis.orbitals, comp.nuclearModel, 
                                                                comp.grid, comp.asfSettings; printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()  error("Not yet implemented.")
    else  error("Unsupported cascade approach.")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationsForPhotoionization(multiplets::Array{Multiplet,1},  scheme::PhotoIonizationScheme, nm::Nuclear.Model)`  
    ... generates all possible (photo-absorption) configurations with a single displacements of an electron from
        scheme.excitationFromShells  to either scheme.excitationToShells or into partial waves with scheme.lValues.
        A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) is returned.
"""
function generateConfigurationsForPhotoionization(multiplets::Array{Multiplet,1},  scheme::PhotoIonizationScheme, nm::Nuclear.Model)
    # Determine all (reference) configurations from multiplets and generate the 'excited/ionized' configurations due to the 
    # specificed excitation/ionization
    initialConfList = Configuration[];   ionConfList = Configuration[];   excConfList = Configuration[];
    for mp  in  multiplets   
        confList = Basics.extractConfigurations(Basics.FromBasis(), mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
    end
    # 
    # Generate all photoionized configurations if photoionization is to be considered; no configuration need to be excluded
    # since parity is give by the partial waves.
    ionConfList = Basics.generateConfigurations(Basics.RemoveElectrons(1, scheme.excitationFromShells), initialConfList)
    ionConfList = unique(ionConfList)
    #
    return( initialConfList, ionConfList )
end


"""
`Cascade.perform(scheme::PhotoIonizationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a photo-excitation computation that starts from a given set of initial configurations xor initial multiplets
        and comprises (various) photoexcitation processes into configurations with up-to NoExcitations single-electron excitations with 
        regard to the initial multiplets. The results of these excitation are comprised into (output) data::PhotoExcData, while these 
        data are only printed during the generation and nothing is returned.

`Cascade.perform(scheme::PhotoIonizationScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::PhotoIonizationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
    ##
    ##
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
    wa  = Cascade.generateConfigurationsForPhotoionization(multiplets, comp.scheme, comp.nuclearModel)
    wb1 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[1], sa="(initial part of the) photoionization ")
    wb2 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[2], sa="(generated part of the) photoionization ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb1)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb2, printout=false)
    Cascade.displayBlocks(stdout, wc1, sa="for the (initial part of the) photoabsorption cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="for the (photo-ionized part of the) photoabsorption cascade ")
    if  printSummary    Cascade.displayBlocks(iostream, wc1, sa="for the (initial part of the) photoionization cascade ");      
                        Cascade.displayBlocks(iostream, wc2, sa="for the (photo-ionized part of the) photoionization cascade ")    end      
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     for block in wc2  push!(gMultiplets, block.multiplet)    end
    #
    we = Cascade.determineSteps(scheme, comp, wc1, wc2)
    Cascade.displaySteps(stdout, we, sa="ionized ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="ionized ")    end      
    wf     = Cascade.modifySteps(we)
    linesP = Cascade.computeSteps(scheme, comp, wf)
    data   = Cascade.Data[];    push!(data, Cascade.Data{PhotoIonization.Line}(linesP) )
    if output    
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )    
        results = Base.merge( results, Dict("photoionization lines:"        => linesP) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        if  outputToFile
            filename = "zzz-cascade-photoionization-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
            JLD2.@save filename results
        end
    end

    return( results )
end
