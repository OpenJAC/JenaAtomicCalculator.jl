
# Functions and methods for scheme::Cascade.StepwiseDecayScheme computations

"""
`Cascade.computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested transition amplitudes as well as PhotoEmission.Line's, AutoIonization.Line's, 
        etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
        processes, however, the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.DecayData  is returned.
"""
function computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    cOrbitals = Dict{Subshell, Orbital}()    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
                                                        
    nt = 0;   st = 0;   previousMeanEn = 0.
    for  step  in  stepList
        st = st + 1
        nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels)
        println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc decay lines (without selection rules): ")
        if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                "up to $nc decay lines (without selection rules): ")   end 
                                                
        if      step.process == Basics.Auger()  &&   comp.approach == Cascade.AverageSCA()
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
                    wp           = Basics.computePotential(Basics.DFSField(), comp.grid, step.finalMultiplet.levels[1])
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
        println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                "to a total of $nt $(string(step.process)) decay lines." )
        if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                "giving now rise to a total of $nt $(string(step.process)) decay lines." )   end      
    end
    #
    data = [ Cascade.Data{PhotoEmission.Line}(linesR), Cascade.Data{AutoIonization.Line}(linesA) ]
end


#=====================================================  
    ... Beginning of Aloka's part ... for parallel computing via multiprocessing and multithreading
    
"""
`Cascade.computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
    ... computes in turn all the requested transition amplitudes as well as PhotoEmission.Line's, AutoIonization.Line's, 
        etc. for all pre-specified decay steps of the cascade. When compared with standard computations of these atomic 
        processes, however, the amount of output is largely reduced and often just printed into the summary file. 
        A set of  data::Cascade.DecayData  is returned.
"""
function computeSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
    #linesA = AutoIonization.Line[];    #linesR = PhotoEmission.Line[];    cOrbitals = Dict{Subshell, Orbital}()    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")

    stepsA = [ step for step in stepList if step.process ==  Basics.Auger() ]
    stepsR = [ step for step in stepList if step.process ==  Basics.Radiative() ]

    println("No of AutoIonization steps", length(stepsA))
    println("No of Radiative steps", length(stepsR))

    @assert length(stepList) == ( length(stepsA) + length( stepsR) )
                                                
    if  comp.approach == Cascade.AverageSCA()
        
        if nworkers() > 1
            linesA = @showprogress desc="Comp. AI Lines..." pmap( step -> computeStepAugerAverageSCA(step, comp), stepsA )
        else
            pA = Progress(length(stepsA); desc="Comp. AI Lines...")
            linesA = Vector{Vector{AutoIonization.Line}}(undef, length(stepsA))

            Threads.@threads for i = 1:length(linesA)
                linesA[i] = computeStepAugerAverageSCA(stepsA[i], comp)
                println("the AutoIonization lines from step $i is ", length(linesA[i]) )
                next!(pA)
            end
            finish!(pA)
        end

    else
        # Compute continuum orbitals independently for all transitions in the given block.
        linesA = @showprogress pmap( step -> AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet,
                                        comp.nuclearModel, comp.grid, step.settings, output=true, printout=false), stepsA )
    end
    
    linesA = collect(Iterators.flatten(linesA))

    
    if nworkers() > 1
        linesR = @showprogress desc="Comp. PE Lines..." pmap( step -> PhotoEmission.computeLinesCascade(step.finalMultiplet, 
                                    step.initialMultiplet, comp.grid, step.settings, output=true, printout=false), stepsR )
    else
        pR = Progress(length(stepsR); desc="Comp. PE Lines...")
        linesR = Vector{Vector{PhotoEmission.Line}}(undef, length(stepsR))

        Threads.@threads for i = 1:length(linesR)
            linesR[i] = PhotoEmission.computeLinesCascade(stepsR[i].finalMultiplet, stepsR[i].initialMultiplet, comp.grid,
                                                                            stepsR[i].settings, output=true, printout=false)
            println("the Radiative lines from step $i is ", length(linesR[i]) )
            next!(pR)
        end
        finish!(pR)
    end
    
    linesR = collect(Iterators.flatten(linesR))
    
    #
    data = [ Cascade.Data{PhotoEmission.Line}(linesR), Cascade.Data{AutoIonization.Line}(linesA) ]
end


"""
`Cascade.computeStepAugerAverageSCA(step::Cascade.Step, comp::Cascade.Computation)` 
    ... computes the AutoIonization.Line for AverageSCA approach for one pre-specified decay step of the cascade. 
        An Array{AutoIonization.Line,1} is returned.
"""
function computeStepAugerAverageSCA(step::Cascade.Step, comp::Cascade.Computation)
    # First determine the `mean' free-electron energy for this Auger block and calculate a common set of continuum orbital
    meanEn = 0.;    NoEn = 0
    for  p = 1:length(step.initialMultiplet.levels),  q = 1:length(step.finalMultiplet.levels)
        en = step.initialMultiplet.levels[p].energy - step.finalMultiplet.levels[q].energy
        if  en > 0.1    meanEn = meanEn + en;    NoEn = NoEn + 1    end
    end
    if  NoEn > 0     meanEn = meanEn/ NoEn     else     meanEn = 0.1    end
    #if  abs(meanEn - previousMeanEn) / meanEn  < 0.15   # no new continuum orbitals
    #    println(">> No new continum orbitals are generated for $(keys(cOrbitals)) and for the energy $meanEn ")
    #else
        previousMeanEn = meanEn;    cOrbitals = Dict{Subshell, Orbital}()
        for kappa = -step.settings.maxKappa-1:step.settings.maxKappa        if   kappa == 0     continue    end
            sh           = Subshell(101, kappa);     nrContinuum = Continuum.gridConsistency(meanEn, comp.grid)
            contSettings = Continuum.Settings(false, nrContinuum);   
            npot         = Nuclear.nuclearPotential(comp.nuclearModel, comp.grid)
            ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
            ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
            ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
            ## wp           = Basics.compute("radial potential: Dirac-Fock-Slater", comp.grid, step.finalMultiplet.levels[1].basis)
            wp           = Basics.computePotentialDFS(comp.grid, step.finalMultiplet.levels[1])
            pot          = Basics.add(npot, wp)
            cOrbital, phase, normF  = Continuum.generateOrbitalLocalPotential(meanEn, sh, pot, contSettings)
            cOrbitals[sh] = cOrbital
        end
        #println(">> Generate continum orbitals in DFS potential for $(keys(cOrbitals)) and for the energy $meanEn ")
    #end

    newLines = AutoIonization.computeLinesFromOrbitals(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                        step.settings, cOrbitals, output=true, printout=false)
    # Deallocating
    cOrbitals = Dict{Subshell, Orbital}()
    return newLines
end
... End of Aloka's part   ===============================================================#  


"""
`Cascade.determineSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})`  
    ... determines all step::Cascade.Step's that need to be computed for this decay cascade. It cycles through all processes of the given
        decay scheme and selects all pairs of blocks due to the selected cascade approach. It is checked that the (averaged) energies 
        each block or level supports a `decay' within the step. A stepList::Array{Cascade.Step,1} is returned, and for which subsequently
        all required transition amplitudes and rates are computed.
"""
function determineSteps(scheme::Cascade.StepwiseDecayScheme, comp::Cascade.Computation, blockList::Array{Cascade.Block,1})
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        for  a = 1:length(blockList)
            for  b = 1:length(blockList)
                minEn = 100000.;   maxEn = -100000.
                for  p = 1:length(blockList[a].multiplet.levels),  q = 1:length(blockList[b].multiplet.levels)
                    minEn = min(minEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                    maxEn = max(maxEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                end
                for  process  in  comp.scheme.processes
                    if      process == Basics.Radiative()   
                        if  a == b   ||   maxEn < 0.    continue   end
                        if  blockList[a].NoElectrons == blockList[b].NoElectrons
                            settings = PhotoEmission.Settings([E1], [UseBabushkin,UseCoulomb], false, false, CorePolarization(), LineSelection(), 0., 0., 1.0e6, false)
                            push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                            blockList[a].multiplet, blockList[b].multiplet) )
                        end
                    elseif  process == Basics.Auger()       
                        if  a == b   ||   maxEn < 0.    continue   end
                        if  blockList[a].NoElectrons == blockList[b].NoElectrons + 1
                            settings = AutoIonization.Settings(AutoIonization.Settings(), maxKappa=3)
                            push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                            blockList[a].multiplet, blockList[b].multiplet) )
                        end
                    else    error("stop a")
                    end
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end
    return( stepList )
end


#==  October 2025, replaced by Basics.ForStepwiseDecay(), ...        
"""
`Cascade.generateConfigurationsForStepwiseDecay(scheme::Cascade.StepwiseDecayScheme, initialConfigs::Array{Configurations,1})`  
    ... generates all possible configurations as obtained by a stepwise loss of maxElectronLoss electrons and by considering
        NoShakeDisplacements displacements. A confList::Array{Configuration,1} is returned.
"""
function generateConfigurationsForStepwiseDecay(scheme::Cascade.StepwiseDecayScheme, initialConfigs::Array{Configuration,1})
    initialNo = initialConfigs[1].NoElectrons
    
    #
    decayShells = scheme.decayShells;   maxElectronLoss = scheme.maxElectronLoss
    # Build configurations with all decayShells 'in between', even if zero occupation
    newConfigs = Configuration[]
    for  conf  in  initialConfigs  
        nshells = copy(conf.shells)
        for shell in decayShells    if  haskey(nshells, shell)   else    nshells[shell] = 0    end   end
        push!(newConfigs, Configuration(nshells, conf.NoElectrons))
    end
    newConfigs = Basics.extractConfigurations(Basics.ByNumber([initialNo-maxElectronLoss]), newConfigs)
    @show "a", newConfigs
    # Now add all decay configurations
    allConfigs   = copy(newConfigs)
    decayConfigs = Configuration[];    
    #
    if      scheme.maxElectronLoss == 0
        # Just move the hole outwards into the decayShells
        decayConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), allConfigs)
        decayConfigs = unique(decayConfigs)
        allConfigs   = append!(allConfigs, decayConfigs)
        return( allConfigs )
        #
    elseif  scheme.maxElectronLoss == 1
        decayConfigs1 = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), allConfigs)
        decayConfigs1 = unique(decayConfigs1)
        decayConfigs2 = Basics.generateConfigurations(Basics.RemoveElectrons(2, decayShells), allConfigs)
        decayConfigs2 = unique(decayConfigs2)
        decayConfigs3 = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), decayConfigs2)
        decayConfigs3 = unique(decayConfigs3)
        allConfigs    = append!(allConfigs, decayConfigs1, decayConfigs2, decayConfigs3)
        allConfigs    = unique(allConfigs)
        return( allConfigs )
        #
    else
        # The general case is quite complex and follows below
    end        
    # 
    # This earlier code for scheme.maxElectronLoss > 1 is highly inefficient if medium and heavy ions are considered.
    further = true;   dConfigs = copy(newConfigs)
    while  further
        dConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), dConfigs)
        if length(dConfigs) > 0     append!(decayConfigs, dConfigs)     else   further = false      end
    end
    decayConfigs = unique(decayConfigs)
    decayConfigs = Basics.extractConfigurations(Basics.ByNumber([initialNo-maxElectronLoss]), decayConfigs)
    allConfigs   = append!(allConfigs, decayConfigs)
    dConfigs     = copy(newConfigs)
    #
    further = true
    while  further
        dConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(2, decayShells), dConfigs)
        if length(dConfigs) > 0     append!(decayConfigs, dConfigs)     else   further = false      end
    end
    decayConfigs = unique(decayConfigs)
    decayConfigs = Basics.extractConfigurations(Basics.ByNumber([initialNo-maxElectronLoss]), decayConfigs)
    allConfigs   = append!(allConfigs, decayConfigs)
    dConfigs     = copy(decayConfigs)
    #
    further = true
    while  further
        dConfigs = Basics.generateConfigurations(Basics.RemoveElectrons(1, decayShells), dConfigs)
        if length(dConfigs) > 0     append!(decayConfigs, dConfigs)     else   further = false      end
    end
    allConfigs   = append!(allConfigs, decayConfigs)
    allConfigs   = unique(allConfigs)
    allConfigs   = Basics.extractConfigurations(Basics.ByNumber([initialNo-maxElectronLoss]), allConfigs)
    #
    # Discard configurations with energies higher than initial configurations after electron capture;
    # this may be required for a very large number of configurations (not yet)

    return( allConfigs )
end
==#


"""
`Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
        various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
        to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
        otherwise.

`Cascade.perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function perform(scheme::StepwiseDecayScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the SCF and CI computation for the intial-state multiplets if initial configurations are given
    if  comp.initialConfigs != Configuration[]
        multiplet  = SelfConsistent.performSCF(comp.initialConfigs, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        basis      = multiplet.levels[1].basis
        multiplets = [Multiplet("initial states", multiplet.levels)]
    else
        basis      = comp.initialMultiplets[1].levels[1].basis
        multiplets = comp.initialMultiplets
    end
    # Print out initial configurations and levels 
    Cascade.displayLevels(stdout, multiplets, sa="initial ")
    if  printSummary   Cascade.displayLevels(iostream, multiplets, sa="initial ")          end
    #
    # Generate subsequent cascade configurations as well as display and group them together
    initialConfigs  = Basics.extractConfigurations(Basics.FromMultiplet(), multiplets)
    @show initialConfigs
    ## wax = Cascade.generateConfigurationList(multiplets, comp.scheme.maxElectronLoss, comp.scheme.NoShakeDisplacements)
    ## @show wax
    ## wbx = Basics.displayConfigurations(comp.nuclearModel.Z, wax, sa="OLD decay ")
    ## @show wax
    decayConfigs    = Basics.generateConfigurations(Basics.ForStepwiseDecay(comp.scheme.maxElectronLoss), initialConfigs)
    append!(decayConfigs, initialConfigs)
    NoElectrons     = unique( Basics.extractFromConfigurations(Basics.NumberOfElectrons(), decayConfigs) )
    Basics.displayConfigurations(stdout, Basics.ByNumber(NoElectrons), decayConfigs, details="decay configurations")
    #
    # Determine first all configuration 'blocks' and from them the individual steps of the cascade
    wc = Cascade.generateBlocks(comp, decayConfigs, basis.orbitals, sa="for the decay cascade:")
    Cascade.displayBlocks(stdout, wc, sa="for the decay cascade ")
    if  printSummary   Cascade.displayBlocks(iostream, wc, sa="for the decay cascade ")    end 
    gMultiplets = Multiplet[];     for block in wc  push!(gMultiplets, block.multiplet)    end
    # Determine, modify and compute the transition data for all steps, ie. the PhotoEmission.Line's, the AutoIonization.Line's, etc.
    wd = Cascade.determineSteps(scheme, comp, wc)
    Cascade.displaySteps(stdout, wd, sa="decay ")
    if  printSummary   Cascade.displaySteps(iostream, wd, sa="decay ")    end      
    we   = Cascade.modifySteps(wd)
    @time data = Cascade.computeSteps(scheme, comp, we)
    @show typeof(data)
    if output    
        results = Base.merge( results, Dict("name"                  => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"        => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"   => multiplets) )    
        results = Base.merge( results, Dict("generated multiplets:" => gMultiplets) ) 
        results = Base.merge( results, Dict("photoemission lines:"  => data[1].lines) ) 
        results = Base.merge( results, Dict("autoionization lines:" => data[2].lines) ) 
        results = Base.merge( results, Dict("cascade data:"         => data) )
        #
        #  Write out the result to file to later continue with simulations on the cascade data
        if outputToFile
            filename = "zzz-cascade-decay-computations-" * string(Dates.now())[1:13] * ".jld"
            println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
            if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                    "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
            JLD2.@save filename results
        end
    else    results = nothing
    end
    
    return( results )
end
