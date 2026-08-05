
# Functions and methods for scheme::Cascade.ExpansionOpacityScheme computations


"""
`Cascade.computeSteps(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested photon excitation/absorption amplitudes as well as PhotoExcitation.Line's for all 
        pre-specified decay steps of the cascade. When compared with standard computations of photoexcitation, however, the amount 
        of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.ExcitationData  is returned.
"""
function computeSteps(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesE = PhotoExcitation.Line[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    nt = 0;   st = 0
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        sa = "\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): "
        println(sa);    if  printSummary   println(iostream, sa)   end 
                                                
        if      step.process == Basics.PhotoExc()
            newLines = PhotoExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                            step.settings, output=true, printout=false)
            if  scheme.printTransitions 
                println("\n ****** Prepare a printout of the emission rates for photoexcitation lines (not yet !!). ****** \n")
                ## PhotoExcitation.displayEmissionRates(stdout, newLines, PhotoExcitation.Settings())
            end 
            append!(linesE, newLines);    nt = length(linesE)
        else   error("Unsupported atomic process for cascade computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) excitation lines."
        println(sa);    if  printSummary   println(iostream, sa)   end 
    end
    #
    data = Cascade.ExcitationData(linesE)
end


"""
`Cascade.determineSteps(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this expansion opacity cascade. 
        It considers the pairwise photoexcitation between all blocks and checks that at least one transition is allowed for these blocks.
        A stepList::Array{Cascade.Step,1} is returned, and for which subsequently all required transition amplitudes and oscillator strengths
        are computed.
"""
function determineSteps(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  blocka  in  blockList
            for  blockb  in  blockList
                if  blocka.NoElectrons  !=  blockb.NoElectrons   error("stop a")     end
                # Check that at least one energy supports photoexcitation/photoabsortion
                settings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=scheme.multipoles, gauges=[UseCoulomb, UseBabushkin],
                                                    mimimumPhotonEnergy=scheme.minPhotonEnergy, maximumPhotonEnergy=scheme.maxPhotonEnergy)
                push!( stepList, Cascade.Step(Basics.PhotoExc(), settings, blocka.confs, blockb.confs, blocka.multiplet, blockb.multiplet) )
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


"""
`Cascade.generateBlocks(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this expansion opacity cascade, and compute also the corresponding 
        multiplets. The different cascade approches realizes different strategies how these blocks are selected and computed. 
        A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(scheme::Cascade.ExpansionOpacityScheme, comp::Cascade.Computation, confs::Array{Configuration,1}; printout::Bool=true)
    blockList = Cascade.Block[];    basis = Basis()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    if    comp.approach == AverageSCA()
        sa = "\n* Generate blocks for expansion opacity computations: \n" *
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
            hydrogenicOrbitals = Bsplines.generateOrbitalsHydrogenic(wa, comp.nuclearModel, subshellList; printout=printout)
        end
        
        for  (ia, confa)  in  enumerate(confs)
            sa = "  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... "
            print(sa);      if  printSummary   println(iostream, sa)   end
            # Now distinguish between the first and all other blocks; for the first block, a SCF is generated and the occupied orbital
            # used also for all other blocks. In addition, a set of hydrogenic orbitals generated for later use
            if  ia == 1
                basis     = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
            else
                # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
                relconfList  = ConfigurationR[]
                #
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
`Cascade.generateConfigurationsForExpansionOpacity(initialConfigs::Array{Configuration,1}, scheme::ExpansionOpacityScheme, 
                                                    nm::Nuclear.Model, grid::Radial.Grid)`  
    ... generates all excited configurations for the expansion opacity computations due the given fromShells, toShells and number
        of excited electrons. A confList::Array{Configurations,1} is returned.
"""
function generateConfigurationsForExpansionOpacity(initialConfigs::Array{Configuration,1}, scheme::ExpansionOpacityScheme,
                                                    nm::Nuclear.Model, grid::Radial.Grid)
    newConfigs = Basics.generateConfigurations(initialConfigs, scheme.excitationFromShells, scheme.excitationToShells, scheme.NoExcitations)  

    return( newConfigs )
end


"""
`Cascade.perform(scheme::ExpansionOpacityScheme, comp::Cascade.Computation)`  
    ... to set-up and perform an expansion opacity computation; it starts from a given set of initial configurations xor initial 
        multiplets and (1) generates all excited configurations with regard to the initial configuration, (2) selects all 
        radiative photoabsorption steps and (3) computes the corresponding transition amplitudes. The results of these expansion opacity 
        computations are comprised into (output) data::ExcitationData, while these data are only printed during the generation and 
        nothing is returned.

`Cascade.perform(scheme::ExpansionOpacityScheme, comp::Cascade.Computation; output=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::ExpansionOpacityScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
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
    wa  = Cascade.generateConfigurationsForExpansionOpacity(comp.initialConfigs, comp.scheme, comp.nuclearModel, comp.grid)
    wb  = Basics.displayConfigurations(comp.nuclearModel.Z, wa, sa="excited configurations of the expansion opacity ")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc  = Cascade.generateBlocks(scheme, comp::Cascade.Computation, wb)
    #
    Cascade.displayBlocks(stdout, wc, sa="from the expansion opacity cascade ")      
    if  printSummary   Cascade.displayBlocks(iostream, wc, sa="from the expansion opacity cascade ")           end      
    #
    # Determine, modify and compute the transition data for all steps, ie. the PhotoExcitation.Line's, etc.
    gMultiplets = Multiplet[];     
    for block in wc  push!(gMultiplets, block.multiplet)    end
    #
    we = Cascade.determineSteps(scheme, comp, wc)
    Cascade.displaySteps(stdout, we, sa="expansion opacity ")
    if  printSummary   Cascade.displaySteps(iostream, we, sa="expansion opacity ")    end      
    wf = Cascade.modifySteps(we)
    #
    data = Cascade.computeSteps(scheme, comp, wf)
    if output    
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:"         => gMultiplets) )    
        results = Base.merge( results, Dict("photoexcitation line data:"    => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-expansion-opacity-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        JLD2.@save filename results
    end
    ## return( results )
    return( results )
end
