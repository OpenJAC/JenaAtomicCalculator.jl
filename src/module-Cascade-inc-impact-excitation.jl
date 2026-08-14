
# Functions and methods for impact-excitation (cascade) computation

"""
`Cascade.computeSteps(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested transition amplitudes as well as ImpactExcitation.Line's, etc. for all pre-specified 
        excitation steps of the cascade. When compared with standard excitation computations, however, the amount of output is 
        largely reduced and often just printed into the summary file. A set of  data::Cascade.Data  is returned.
"""
function computeSteps(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesE = ImpactExcitation.Line[]    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc impact-excitation lines: ")
        if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                "up to $nc impact-excitation lines: ")   end 
                                                
        if  step.process == Basics.ImpactExc() 
            ## The nuclear model has to be passed on: the continuum orbital of the scattered electron is generated
            ## in its field.  The former call omitted it, but that could not be noticed because
            ## ImpactExcitation.computeLinesCascade did not exist at all.
            newLines = ImpactExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel,
                                                            comp.grid, step.settings, output=true, printout=false)
            append!(linesE, newLines);    nt = length(linesE)
        else   error("Unsupported atomic process for excitation computations.")
        end
        println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) impact-excitation lines." )
        if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                "giving now rise to a total of $nt $(string(step.process)) impact-excitation lines." )   end      
    end
    #
    return( linesE )
end


"""
`Cascade.determineSteps(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
        decay scheme and selects all pairs of blocks due to the selected cascade approach. It checks that at least on pair of levels
        supports a `electron-impact excitation' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which 
        subsequently all required transition amplitudes and rates/cross sections are computed.
"""
function determineSteps(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, 
                        initialList::Array{Cascade.Block,1}, excitedList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  initialBlock in initialList
            for  excitedBlock in excitedList
                if  initialBlock.NoElectrons == excitedBlock.NoElectrons
                    ## The former call was ImpactExcitation.Settings(false, false, false, false, LineSelection()),
                    ## i.e. five positional arguments of the wrong types against a ten-field struct -- a MethodError
                    ## before any physics happened, and it also passed no electron energies at all.  Take the
                    ## energies, the partial-wave cap and the energy shift from the scheme, as the other schemes do.
                    maxKappa = length(scheme.lValues) > 0 ? maximum(scheme.lValues) + 1 : 300
                    ## calcRateCoefficient is deliberately NOT set here.  A cascade COMPUTATION produces lines and
                    ## nothing else; turning lines into collision strengths and rate coefficients is the job of a
                    ## cascade SIMULATION, which is the distinction drawn in Symmetry 13, 520 (2021) and
                    ## Eur. Phys. J. D 78, 75 (2024).  See Cascade.EieRateCoefficients for the simulation side.
                    settings = ImpactExcitation.Settings(ImpactExcitation.Settings(),
                                                         electronEnergies = scheme.electronEnergies,
                                                         energyShift      = scheme.electronEnergyShift,
                                                         maxKappa         = maxKappa)
                    push!( stepList, Cascade.Step(Basics.ImpactExc(), settings, initialBlock.confs, excitedBlock.confs, 
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
`Cascade.generateBlocks(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this impact-excitation cascade, and compute also 
        the corresponding multiplets. The different cascade approches realized different strategies how these block are selected 
        and computed. A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.ImpactExcitationScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        if  printout
        println("\n* Generate blocks for excitation computations:")
        println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
        println("    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; ")
        println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
        println("    + only E1 excitations are considered. \n")
        if  printSummary   
        println(iostream, "\n* Generate blocks for excitation computations:")
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
            ## SelfConsistent.performSCF returns a Multiplet, not a Basis; the orbitals sit one level deeper.
            ## The variable was still named and used as if it were a Basis, which raised a FieldError here.
            scfMultiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            multiplet    = Hamiltonian.performCIwithFrozenOrbitals([confa], scfMultiplet.levels[1].basis.orbitals, comp.nuclearModel,
                                                                   comp.grid, Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings); printout=false)
            push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
            println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
        end
    elseif    comp.approach == SCA()  error("Not yet implemented.")
    else  error("Unsupported cascade approach.")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationsForImpactExcitation(multiplets::Array{Multiplet,1},  scheme::ImpactExcitationScheme, nm::Nuclear.Model)`  
    ... generates all possible (electron-impact excited) configurations due to the single excitation scheme.fromShells to
        scheme.toShells with regard to the given mutiplets. 
        A Tuple(initialConfList::Array{Configuration,1}, confList::Array{Configuration,1}) is returned.
"""
function generateConfigurationsForImpactExcitation(multiplets::Array{Multiplet,1},  scheme::ImpactExcitationScheme, nm::Nuclear.Model)
    # Determine all (reference) configurations from multiplets and generate the 'excited' configurations due to the specificed excitations
    initialConfList = Configuration[]
    for mp  in  multiplets   
        ## Basics.extractNonrelativisticConfigurations no longer exists; every other cascade scheme uses the
        ## FromBasis theme of Basics.extractConfigurations for exactly this purpose.
        confList = Basics.extractConfigurations(Basics.FromBasis(), mp.levels[1].basis)
        for  conf in confList   if  conf in initialConfList   nothing   else   push!(initialConfList, conf)      end      end
    end
    blockConfList = Basics.generateConfigurations(initialConfList, scheme.fromShells, scheme.toShells, 1)
    # Exclude configurations with too low or too high mean energies 
    en     = Float64[];   for conf in initialConfList    push!(en, -Semiempirical.estimate(EstimateBindingEnergyWilliams2000(), round(Int64, nm.Z), conf))   end
    maxen  = maximum(en);    minen  = minimum(en);  
    println(">>> initial configuration(s) have energies from $minen  to  $maxen  [a.u.].")
    #
    if      scheme.NoFreeElectronEnergies > 0  mx = scheme.maxFreeElectronEnergy 
    else                                       mx = maximum( scheme.electronEnergies )
    end
    maxen = maxen + Defaults.convertUnits("energy: from predefined to atomic unit", mx) 
    #
    newBlockConfList = Configuration[]
    for  conf  in  blockConfList    meanEnergy = -Semiempirical.estimate(EstimateBindingEnergyWilliams2000(), round(Int64, nm.Z), conf)
        if  meanEnergy <= maxen     push!(newBlockConfList, conf) 
        else    println(">>> exclude $conf with energy $meanEnergy [a.u.] because of energy reasons; maximum total energy = $maxen ")
        end
    end

    return( (initialConfList, newBlockConfList)  )
end


"""
`Cascade.perform(scheme::ImpactExcitationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a electron-impact excitation computation that starts from a given set of initial configurations 
        xor initial multiplets and comprises the impact-excitation of different fromShells to the toShells. All configurations 
        are constructed with regard to the initial configurations. The results of these impact-excitation computations 
        are comprised into (output) data::??? but nothing is returned as default.

`Cascade.perform(scheme::ImpactExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::ImpactExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
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
    wa  = Cascade.generateConfigurationsForImpactExcitation(multiplets, comp.scheme, comp.nuclearModel)
    wb1 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[1], sa="(initial part of the) impact-excited ")
    wb2 = Basics.displayConfigurations(comp.nuclearModel.Z, wa[2], sa="(generated part of the) impact-excited ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc1 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb1)
    wc2 = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb2, printout=false)
    Cascade.displayBlocks(stdout, wc1, sa="for the (initial part of the) impact-excited cascade ");      
    Cascade.displayBlocks(stdout, wc2, sa="for the (generated part of the) impact-excited cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc1, sa="for the (initial part of the) impact-excited cascade ")
                        Cascade.displayBlocks(iostream, wc2, sa="for the (generated part of the) impact-excited cascade ")    end   
    # Determine, modify and compute the transition data for all steps, ie. the PhotoIonization.Line's, etc.
    gMultiplets = Multiplet[];     for block in wc2  push!(gMultiplets, block.multiplet)    end
    we = Cascade.determineSteps(scheme, comp, wc1, wc2)
    Cascade.displaySteps(stdout, we, sa="excited ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="impact-exciting ")    end      
    wf      = Cascade.modifySteps(we)
    linesE  = Cascade.computeSteps(scheme, comp, wf)
    data    = Cascade.Data[];    push!(data, Cascade.Data{ImpactExcitation.Line}(linesE) )
    if output    
        results = Base.merge( results, Dict("name"                       => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"             => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"        => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"      => gMultiplets) )    
        results = Base.merge( results, Dict("impact-excitation lines:"   => linesE) )
        results = Base.merge( results, Dict("cascade data:"              => data ) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        if  outputToFile
            ## Was "zzz-cascade-photoexcitation-computations-...", a copy-paste artefact: this scheme writes
            ## impact-excitation lines, not photo-excitation ones.
            filename = "zzz-cascade-impact-excitation-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
            Cascade.writeDataFile(filename, results)
        end
    end
    ## return( results )
    return( results )
end
