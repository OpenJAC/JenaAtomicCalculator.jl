
# Functions and methods for cascade computation


"""
`Cascade.computeSteps(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested transition amplitudes as well as PhotoExcitation.Line's, etc. for all pre-specified 
        excitation steps of the cascade. When compared with standard excitation computations, however, the amount of output is 
        largely reduced and often just printed into the summary file. A set of  data::Cascade.ExcitationData  is returned.
"""
function computeSteps(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesE = PhotoExcitation.Line[]    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc excitation lines (without selection rules): ")
        if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                "up to $nc excitation lines (without selection rules): ")   end 
                                                
        if  step.process == Basics.PhotoExc() 
            newLines = PhotoExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                           ## step.settings, scheme.initialLevelSelection, output=true, printout=false) 
                                                           step.settings, LevelSelection(true, indices=[1,2,3,4]), output=true, printout=false) 
            append!(linesE, newLines);    nt = length(linesE)
        else   error("Unsupported atomic process for excitation computations.")
        end
        println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) decay lines." )
        if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                "giving now rise to a total of $nt $(string(step.process)) decay lines." )   end      
    end
    #
    return( linesE )
end


"""
`Cascade.determineSteps(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
        decay scheme and selects all pairs of blocks due to the selected cascade approach. It checks that at least on pair of levels
        supports a `photo-excitation' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently all 
        required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  initialBlock in initialList
            for  excitedBlock in excitedList
                if  initialBlock.NoElectrons == excitedBlock.NoElectrons
                    ## The multipoles must be taken from the scheme; a hard-wired [E1] made the multipole list of
                    ## PhotoExcitationScheme act on the configuration filter (generateConfigurationsForPhotoexcitation)
                    ## only, so that M1/E2/M2 enlarged the cascade by same-parity blocks that could never yield a line.
                    settings = PhotoExcitation.Settings(scheme.multipoles, [UseCoulomb, UseBabushkin], false, false, false, false, false,
                                                        LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes() )
                    push!( stepList, Cascade.Step(Basics.PhotoExc(), settings, initialBlock.confs, excitedBlock.confs, 
                                                                                initialBlock.multiplet, excitedBlock.multiplet) )
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this excitation cascade, and compute also the corresponding multiplets.
        The different cascade approches realized different strategies how these block are selected and computed. A blockList::Array{Cascade.Block,1} 
        is returned.
"""
function generateBlocks(scheme::Cascade.PhotoExcitationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        if  printout
            println("\n* Generate blocks for excitation computations:")
            Cascade.displayApproach(stdout, comp.approach, comp.asfSettings)
        end
        ## Cascade.generateBoundOrbitals prints the bound-orbital assumption itself, to stdout and to the
        ## summary file, so that what is announced can never drift from what is done.
        ## It takes the APPROACH, not its conditions: that signature changed in 46dd775 and this call site was
        ## missed, so every photoexcitation cascade has raised a MethodError since. Nothing caught it because
        ## testModule_Cascade_PhotonExcitation is hollow -- it computes nothing and hardcodes success=true.
        orbitals = Cascade.generateBoundOrbitals(comp.approach, comp, confs; printout=printout)
        if  printout
        println("    + only E1 excitations are considered. \n")
        if  printSummary
        println(iostream, "\n* Generate blocks for excitation computations:")
        println(iostream, "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
        println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
        println(iostream, "    + only E1 excitations are considered. \n")
        end
        end     # printout
        #
        for  confa  in confs
            print("  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")
            if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... ")   end
            ## generateBoundOrbitals returns one orbital set PER CHARGE STATE, keyed by the electron count,
            ## so it has to be indexed here; passing the whole dictionary was the second half of the same
            ## signature change missed in 46dd775.
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa], orbitals[confa.NoElectrons], comp.nuclearModel, comp.grid,
                                                                Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings); printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()  error("Not yet implemented.")
    else  error("Unsupported cascade approach.")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationsForPhotoexcitation(multiplets::Array{Multiplet,1},  scheme::PhotoExcitationScheme, nm::Nuclear.Model)`  
    ... generates all possible (photo-excited) configurations with upto  scheme.NoExcitations of displacements of electron from
        scheme.excitationFromShells  to  scheme.excitationToShells  with regard to the given mutiplets.
        A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) is returned.
"""
function generateConfigurationsForPhotoexcitation(multiplets::Array{Multiplet,1},  scheme::PhotoExcitationScheme, nm::Nuclear.Model)
    # Determine all (reference) configurations from multiplets and generate the 'excited' configurations due to the specificed excitations
    initialConfList = Configuration[]
    for mp  in  multiplets   
        confList = Basics.extractConfigurations(Basics.FromBasis(), mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
    end
    blockConfList = Basics.generateConfigurations(initialConfList, scheme.excitationFromShells, scheme.excitationToShells, scheme.NoExcitations)
    # Exclude configurations with too low or too high mean energies as well as those that are parity forbidden for the given multipoles
    hasPlus = false;   hasMinus = false
    for  conf  in  initialConfList   
        if  Basics.extractFromConfiguration(Basics.GetParity(), conf) == Basics.plus   hasPlus = true    else   hasMinus = true      end
    end
    en     = Float64[];   
    for conf in initialConfList    push!(en, Empirical.totalEnergy(round(Int64, nm.Z), conf, data = PeriodicTable.XrayDataBooklet() ))    end
    maxen  = maximum(en);    minen  = minimum(en);  
    println(">>> initial configuration(s) have energies from $minen  to  $maxen  [a.u.].")
    #
    newBlockConfList = Configuration[]
    for  conf  in  blockConfList    meanEnergy = Empirical.totalEnergy(round(Int64, nm.Z), conf, data = PeriodicTable.XrayDataBooklet() )
        minEn = 0.2* Defaults.convertUnits("energy: from predefined to atomic unit", scheme.minPhotonEnergy) 
        maxEn = 5*   Defaults.convertUnits("energy: from predefined to atomic unit", scheme.maxPhotonEnergy) 
        if  minen + minEn  <= meanEnergy <= maxen + maxEn
            if      hasPlus  &&  hasMinus                                                                        push!(newBlockConfList, conf)   
            elseif  M1 in scheme.multipoles  ||  E2 in scheme.multipoles  ||  M2 in scheme.multipoles            push!(newBlockConfList, conf)  
            elseif  hasPlus  &&  E1 in scheme.multipoles  &&   Basics.extractFromConfiguration(Basics.GetParity(), conf) == Basics.minus 
                push!(newBlockConfList, conf)  
            elseif  hasMinus &&  E1 in scheme.multipoles  &&   Basics.extractFromConfiguration(Basics.GetParity(), conf) == Basics.plus 
                push!(newBlockConfList, conf) 
            end
        else    println(">>> exclude $conf with energy $meanEnergy [a.u.] because of energy reasons." *
                        "\n min=$(minen + minEn) ... max=$(maxen + maxEn) ")
        end
    end

    return( (initialConfList, newBlockConfList)  )
end


"""
`Cascade.perform(scheme::PhotoExcitationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a photo-excitation computation that starts from a given set of initial configurations xor initial multiplets
        and comprises (various) photoexcitation processes into configurations with up-to NoExcitations single-electron excitations with 
        regard to the initial multiplets. The results of these excitation are comprised into (output) data::PhotoExcData, while these 
        data are only printed during the generation and nothing is returned.

`Cascade.perform(scheme::PhotoExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::PhotoExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
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
    wa  = Cascade.generateConfigurationsForPhotoexcitation(multiplets, comp.scheme, comp.nuclearModel)
    wb1 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[1]; sa="(initial part of the) photo-excited ")
    wb2 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[2]; sa="(generated part of the) photo-excited ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb1)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb2, printout=false)
    Cascade.displayBlocks(stdout, wc1, sa="for the (initial part of the) excited cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="for the (generated part of the) excited cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="for the (initial part of the) excited cascade ")
                        Cascade.displayBlocks(iostream, wc2, sa="for the (generated part of the) excited cascade ")    end      
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     for block in wc2  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc1, wc2)
    Cascade.displaySteps(stdout, we, sa="excited ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="ionizing ")    end      
    wf     = Cascade.modifySteps(we)
    linesE = Cascade.computeSteps(scheme, comp, wf)
    data   = Cascade.Data[];    push!(data, Cascade.Data{PhotoExcitation.Line}(linesE) )
    if output    
        results = Base.merge( results, Dict("name"                       => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"             => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"        => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"      => gMultiplets) )    
        results = Base.merge( results, Dict("photoexcitation lines:"     => linesE) )
        results = Base.merge( results, Dict("cascade data:"              => data ) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        if  outputToFile
            filename = "zzz-cascade-photoexcitation-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
            JLD2.@save filename results
        end
    end
    ## return( results )
    return( results )
end
