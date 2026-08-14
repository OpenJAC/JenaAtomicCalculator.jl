
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
            ## computeLinesCascade gained a LevelSelection argument. Here EVERY level of the lower block must be
            ## allowed to absorb -- a Boltzmann distribution puts weight on all of them, and the metastable
            ## levels are precisely the ones that carry the interesting lines -- so no selection is applied.
            newLines = PhotoExcitation.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                            step.settings, LevelSelection(false), output=true, printout=false)
            append!(linesE, newLines);    nt = length(linesE)
        else   error("Unsupported atomic process for cascade computations.")
        end
        sa = "     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) excitation lines."
        println(sa);    if  printSummary   println(iostream, sa)   end 
    end
    #
    ## The line list is printed here, once and sorted by wavelength, rather than step by step: it is the
    ## bound-bound list the expansion opacity is built from, and one wants to see which lines dominate a
    ## wavelength bin.  Previously this branch only announced itself with "(not yet !!)".
    if  scheme.printTransitions
        Cascade.displayExpansionOpacityTransitions(stdout, linesE)
        if  printSummary    Cascade.displayExpansionOpacityTransitions(iostream, linesE)    end
    end
    #
    ## Return the generalised Cascade.Data{T}, not the older per-scheme Cascade.ExcitationData: the rest of the
    ## module -- Cascade.extractPhotoExcitationData in particular -- has moved to Data{T} (field `lines`), and
    ## the half-finished migration is what broke the computation -> simulation chain of this scheme.
    data = [ Cascade.Data{PhotoExcitation.Line}(linesE) ]
end


"""
`Cascade.displayExpansionOpacityTransitions(stream::IO, lines::Array{PhotoExcitation.Line,1})`  
    ... displays the bound-bound transitions from which the expansion opacity is built, sorted by wavelength.
        For each line, the wavelength, the transition energy, the absorption oscillator strength in both gauges
        and the excitation energy of the LOWER level are shown; the latter decides, together with the temperature,
        how strongly the line is populated. Nothing is returned.
"""
function displayExpansionOpacityTransitions(stream::IO, lines::Array{PhotoExcitation.Line,1})
    if  length(lines) == 0    println(stream, "\n  No bound-bound transitions have been calculated. \n");    return( nothing )    end
    ## The lower level of the lowest-lying line is taken as the reference for the excitation energies.
    eGround = minimum( [ line.initialLevel.energy  for line in lines ] )
    sortedLines = sort( lines, by = line -> -line.omega )
    #
    nx = 116
    println(stream, " ")
    println(stream, "  Bound-bound transitions that enter the expansion opacity:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "lambda"; na=3)
    sb = sb * TableStrings.center(14, "[Angstrom]"; na=3)
    sa = sa * TableStrings.center(14, "omega"; na=3)
    sb = sb * TableStrings.center(14, "[eV]"; na=3)
    sa = sa * TableStrings.center(26, "f_ik (absorption)"; na=2)
    sb = sb * TableStrings.center(26, "Coulomb -- Babushkin"; na=2)
    sa = sa * TableStrings.center(14, "E(lower)"; na=2)
    sb = sb * TableStrings.center(14, "[eV]"; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  line in sortedLines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym);  na=4)
        sa = sa * @sprintf("%.5e", convertUnits("energy: from atomic to Angstrom", line.omega))         * "    "
        sa = sa * @sprintf("%.5e", convertUnits("energy: from atomic to eV", line.omega))      * "    "
        sa = sa * @sprintf("%.5e", line.oscStrength.Coulomb)                                            * "    "
        sa = sa * @sprintf("%.5e", line.oscStrength.Babushkin)                                          * "    "
        sa = sa * @sprintf("%.5e", convertUnits("energy: from atomic to eV", line.initialLevel.energy - eGround))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
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
            ## The argument order had drifted:  Bsplines.generateOrbitalsHydrogenic takes
            ## (subshells, nuclearModel, primitives), as the dielectronic-recombination and photorecombination
            ## schemes call it.  This line raised a MethodError on the first block of any cascade.
            hydrogenicOrbitals = Bsplines.generateOrbitalsHydrogenic(subshellList, comp.nuclearModel, wa; printout=printout)
        end
        
        for  (ia, confa)  in  enumerate(confs)
            sa = "  Multiplet computations for $(string(confa)[1:end])   with $(confa.NoElectrons) electrons ... "
            print(sa);      if  printSummary   println(iostream, sa)   end
            # Now distinguish between the first and all other blocks; for the first block, a SCF is generated and the occupied orbital
            # used also for all other blocks. In addition, a set of hydrogenic orbitals generated for later use
            if  ia == 1
                ## SelfConsistent.performSCF returns a Multiplet, not a Basis; basis.orbitals below would fail.
                scfMultiplet  = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                basis         = scfMultiplet.levels[1].basis
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
        sa = "\n* Generate blocks for expansion opacity computations: \n" *
                "\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: " *
                "\n    + each single configuration forms an individual cascade block; " *
                "\n    + orbitals are generated independently for each block for a Dirac-Fock-Slater potential; " *
                "\n    + configuration mixing is included for each block, based on H^(DC); " *
                "\n    + all requested multipoles are considered for the bound-bound transitions. \n"
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
    ## Basics.generateConfigurations returns the EXCITED configurations only. Without the initial ones the
    ## ground configuration never becomes a cascade block, and no line can start from the level that carries
    ## almost the whole population -- the resonance lines, i.e. the strongest contributors to any opacity,
    ## were therefore missing altogether.
    for  conf in initialConfigs
        if  !(conf in newConfigs)    push!(newConfigs, conf)    end
    end

    return( newConfigs )
end


"""
`Cascade.perform(scheme::ExpansionOpacityScheme, comp::Cascade.Computation)`  
    ... to set-up and perform an expansion opacity computation; it starts from a given set of initial configurations xor initial 
        multiplets and (1) generates all excited configurations with regard to the initial configuration, (2) selects all 
        radiative photoabsorption steps and (3) computes the corresponding transition amplitudes. The results of these expansion opacity 
        computations are comprised into (output) data::ExcitationData, while these data are only printed during the generation and 
        nothing is returned.

`Cascade.perform(scheme::ExpansionOpacityScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
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
        ## "photoexcitation lines:" is the key Cascade.extractPhotoExcitationData reads; storing only the
        ## descriptive "photoexcitation line data:" meant the extractor returned an empty array and every
        ## opacity simulation aborted with "No photoexcitationData provided."  Both keys are written: the raw
        ## line list under the standard key, the Data{T} container under the descriptive one.
        results = Base.merge( results, Dict("photoexcitation line data:"    => data) )
        results = Base.merge( results, Dict("photoexcitation lines:"        => data[1].lines) )
        results = Base.merge( results, Dict("cascade data:"                 => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        filename = "zzz-cascade-expansion-opacity-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        Cascade.writeDataFile(filename, results)
    end
    ## return( results )
    return( results )
end
