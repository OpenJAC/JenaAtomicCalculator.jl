
"""
`DielectronicRecombination.addEmpiricalPassages!(passages::Array{DielectronicRecombination.Passage,1},
                                                 empTreatment::EmpiricalTreatment)`  
    ... to add further (empirical) passages due to the capture of electrons into higher n-shells; information
        about the last "captured" shell as well as the interval [nLowerEmpirical, nUpperEmpirical] are taken from
        empTreatment. A simple scaling rules of energies and rates are presently applied bu could be improved
        if needed. The array passages::Array{DielectronicRecombination.Passage,1} is modified and nothing is
        returned otherwise.
"""
function  addEmpiricalPassages!(passages::Array{DielectronicRecombination.Passage,1}, empTreatment::EmpiricalTreatment)
    nLastCaptured = empTreatment.nLowerEmpirical-1;   npIn = length(passages);   Zeff = empTreatment.empiricalEffectiveZ
    
    @show empTreatment.doEmpiricalCorrections, empTreatment.nLowerEmpirical, empTreatment.nUpperEmpirical,
          empTreatment.empiricalEffectiveZ
          
    if  empTreatment.nUpperEmpirical == 0   return(nothing)  end
    #
    for   n = empTreatment.nLowerEmpirical:empTreatment.nUpperEmpirical
        for  np = 1:npIn
            passage = passages[np]
            # Figure out that passage belongs to nLastCaptured ... otherwise skip
            electronEnergy  = passage.electronEnergy - (Zeff/n)^2 / 2.  +  (Zeff/nLastCaptured)^2 / 2.
            captureRate     = passage.captureRate    * (nLastCaptured/n)^(0.8)
            photonRate      = passage.photonRate     * (nLastCaptured/n)^3
            #  Factor due to UserGuide
            wb              = Defaults.convertUnits("kinetic energy to wave number: atomic units", electronEnergy)
            wb              = pi*pi / (wb*wb) * captureRate
                                ((Basics.twice(passage.intermediateLevel.J) + 1) / (Basics.twice(passage.initialLevel.J) + 1))
            reducedStrength = EmProperty(wb * photonRate.Coulomb, wb * photonRate.Babushkin)

            # Check that the resonance is in the desired "window of resonances"
            if  empTreatment.doResonanceWindowCorrection   &&  (electronEnergy <  empTreatment.resonanceEnergyMin  ||   
                                                                electronEnergy >  empTreatment.resonanceEnergyMax)  continue    end
            
            newPassage = DielectronicRecombination.Passage(passage.initialLevel, passage.intermediateLevel, electronEnergy, 
                                                           captureRate, photonRate, reducedStrength, AutoIonization.Channel[])
            push!(passages, newPassage)
        end
    end
    
    return(nothing)
end




"""
`DielectronicRecombination.checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                                     initialMultiplet::Multiplet)`  
    ... to check that the given initial-, intermediate- and final-state levels and multiplets are consistent to each other and
        to avoid later problems with the computations. An error message is issued if an inconsistency occurs,
        and nothing is returned otherwise.
"""
function  checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    initialSubshells      = initialMultiplet.levels[1].basis.subshells;             ni = length(initialSubshells)
    intermediateSubshells = intermediateMultiplet.levels[1].basis.subshells
    finalSubshells        = finalMultiplet.levels[1].basis.subshells;               nf = length(finalSubshells)
    
    if initialSubshells[1:end] == intermediateSubshells[1:ni]   &&
        intermediateSubshells[1:nf]  == finalSubshells
    else
        @show initialSubshells
        @show intermediateSubshells
        @show finalSubshells
        error("\nThe order of subshells must be equal for the initial-, intermediate and final states, \n"     *
                "and the same subshells must occur in the definition of the intermediate and final states. \n" *
                "Only the initial states can have less subshells; this limitation arises from the angular coefficients.")
    end
        
    return( nothing )
end


"""
`DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                                      initialMultiplet::Multiplet)`  
    ... to check (and analyze) that all high nl orbitals in these multiplets are properly represented on the given grid.
        The function prints for each symmetry block kappa the high-nl orbitals and checks that they are all bound.
        An error message is issued if this is not the case, and nothing is returned otherwise.
"""
function  checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    
    println("\n  Check the energies and orbital representation:")
    println("\n  ----------------------------------------------")
    subshells = initialMultiplet.levels[1].basis.subshells;     basis = initialMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Initial occupied subshells:      $(subshells)")
    
    subshells = intermediateMultiplet.levels[1].basis.subshells;     basis = intermediateMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    @show sub, basis.orbitals[sub].energy
            error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Intermediate occupied subshells: $(subshells)")
    
    subshells = finalMultiplet.levels[1].basis.subshells;     basis = finalMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Final occupied subshells:        $(subshells)")
        
    return( nothing )
end


"""
`DielectronicRecombination.computeAmplitudesProperties(passage::DielectronicRecombination.Passage, 
                           finalMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                           empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)` 
    ... to compute all amplitudes and properties of the given line; a line::DielectronicRecombination.Pathway is returned 
        for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(passage::DielectronicRecombination.Passage, 
                                      finalMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                      empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)
    rateA = 0.
    println(">> passage: $(passage.initialLevel.index)--$(passage.intermediateLevel.index) ...")
    newcChannels      = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)  
    initialLevel      = deepcopy(passage.initialLevel)
    intermediateLevel = deepcopy(passage.intermediateLevel)
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    for cChannel in passage.captureChannels
        newnLevel = Basics.generateLevelWithSymmetryReducedBasis(intermediateLevel, intermediateLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, newnLevel.basis.subshells)
        newnLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, cChannel.kappa), newnLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(passage.electronEnergy, Subshell(101, cChannel.kappa), 
                                                             newiLevel, nm, grid, contSettings)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, cChannel.symmetry, newiLevel)
        newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, Complex(0.))
        amplitude   = AutoIonization.amplitude(settings.augerOperator, cChannel, newcLevel, newnLevel, grid)
        rateA       = rateA + conj(amplitude) * amplitude
        newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, amplitude)
        push!( newcChannels, newcChannel)
    end
    #
    # Compute and sum-up the photon rates of all f levels
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    photonEnergyShift = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    photonRate = Basics.EmProperty(0.)
    #
    for  fLevel  in  finalMultiplet.levels
        pEnergy = intermediateLevel.energy - fLevel.energy + photonEnergyShift
        if  pEnergy < 0.    continue    end
        pChannels = DielectronicRecombination.determinePhotonChannels( fLevel, intermediateLevel, settings)
        #
        newpChannels = PhotoEmission.Channel[];    rateC = 0.;    rateB = 0.
        for pChannel in pChannels
            amplitude   = PhotoEmission.amplitude("emission", pChannel.multipole, pChannel.gauge, pEnergy,
                                                  fLevel, intermediateLevel, grid, display=false, printout=false)
            newpChannel = PhotoEmission.Channel( pChannel.multipole, pChannel.gauge, amplitude)
            push!( newpChannels, newpChannel)
            if       newpChannel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
            elseif   newpChannel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
            elseif   newpChannel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
            end
        end
        #
        wa          = 8.0pi * Defaults.getDefaults("alpha") * pEnergy / (Basics.twice(intermediateLevel.J) + 1) *
                                                                        (Basics.twice(fLevel.J) + 1)
        wa          = wa / pi  ## modified for test with Xe^53+ (March/2024)
        photonRate  = photonRate + EmProperty(wa * rateC, wa * rateB)  
        #
    end
    #
    #
    # Now add to the photonRate, if HydrogenicCorrections() are needed
    if  empTreatment.doHydrogenicCorrections
        lmax = empTreatment.maximum_l
        Zeff          = empTreatment.hydrogenicEffectiveZ
        #            
        # Determine nRydberg, lRydberg for the given resonance
        rydbergSubshs = Basics.extractRydbergSubshellList(intermediateLevel, empTreatment.nLowestCaptured-1, 1.0e-1)
        rydbergShells = Basics.extractNonrelativisticShellList(rydbergSubshs)
        #
        if  length(rydbergShells) == 1  rShell = rydbergShells[1];   ni = rShell.n;   li = rShell.l
        else   error("Inappropriate number of Rydberg shells = $rydbergShells ")
        end
        #
        # Compute and add hydrogenic rates A(ni,li --> nDetailed < n <= ni-1, li +- 1)
        hydrogenicRate = 0.
        for  nf = empTreatment.nFinal+1:empTreatment.nHydrogenic
            if  li > lmax
            else     hydrogenicRate = hydrogenicRate + computeHydrogenicRate(ni, li, nf, li-1, Zeff) + 
                                                       computeHydrogenicRate(ni, li, nf, li+1, Zeff)
            end
        end
        photonRate = photonRate + hydrogenicRate * empTreatment.hydrogenicRateScaling
        println(">>> Add hydrogenic corrections from nLower:nUpper = $(empTreatment.nFinal+1):$(empTreatment.nHydrogenic)")
        println(">>> Total photon rate & total hydrogenic rate [a.u.] for A($ni, $li;  -> ...) = $photonRate,  " * 
                "$(hydrogenicRate * empTreatment.hydrogenicRateScaling)")
    end
    #
    #
    captureRate     = 2pi * rateA
    #  Factor due to UserGuide
    wb              = Defaults.convertUnits("kinetic energy to wave number: atomic units", passage.electronEnergy)
    wb              = pi*pi / (wb*wb) * captureRate  * ##  2 * # factor 2 is not really clear.
                        ((Basics.twice(intermediateLevel.J) + 1) / (Basics.twice(initialLevel.J) + 1))
    reducedStrength = EmProperty(wb * photonRate.Coulomb, wb * photonRate.Babushkin)
    #
    #  Factor due to Tu et al. (Plasma Phys., 2016)
    ## wb              = pi*pi / 2. / pathway.electronEnergy * captureRate * 
    ##                   (Basics.twice(pathway.intermediateLevel.J) + 1) / (Basics.twice(pathway.initialLevel.J) + 1)
    #
    passage = DielectronicRecombination.Passage(passage.initialLevel, passage.intermediateLevel, passage.electronEnergy,
                                                captureRate.re, photonRate, reducedStrength, newcChannels)
    # 
    return( passage )
end


"""
` (pathway::DielectronicRecombination.Pathway, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
   settings::DielectronicRecombination.Settings, hasCaptureChannels::Bool, 
   lastCaptureChannels::Array{AutoIonization.Channel,1})` 
    ... to compute all amplitudes and properties of the given line; a line::DielectronicRecombination.Pathway is returned 
        for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(pathway::DielectronicRecombination.Pathway, nm::Nuclear.Model, grid::Radial.Grid, 
                                      nrContinuum::Int64, settings::DielectronicRecombination.Settings, 
                                      hasCaptureChannels::Bool, lastCaptureChannels::Array{AutoIonization.Channel,1})
    rateA = 0.
    # println(">> pathway: $(pathway.initialLevel.index)--$(pathway.intermediateLevel.index)--$(pathway.finalLevel.index) ...")
    
    if hasCaptureChannels
        # Simply copy the results from previous computation of the same channels
        newcChannels = deepcopy(lastCaptureChannels)
        for cChannel in newcChannels
            rateA     = rateA + conj(cChannel.amplitude) * cChannel.amplitude
        end
    else
        newcChannels      = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)  
        initialLevel      = deepcopy(pathway.initialLevel)
        intermediateLevel = deepcopy(pathway.intermediateLevel)
        Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
        for cChannel in pathway.captureChannels
            newnLevel = Basics.generateLevelWithSymmetryReducedBasis(intermediateLevel, intermediateLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, newnLevel.basis.subshells)
            newnLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, cChannel.kappa), newnLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(pathway.electronEnergy, Subshell(101, cChannel.kappa), 
                                                                    newiLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, cChannel.symmetry, newiLevel)
            newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, Complex(0.))
            amplitude   = AutoIonization.amplitude(settings.augerOperator, cChannel, newcLevel, newnLevel, grid)
            rateA       = rateA + conj(amplitude) * amplitude
            newcChannel = AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, amplitude)
            push!( newcChannels, newcChannel)
        end
        println(">> pathway: $(pathway.initialLevel.index)--$(pathway.intermediateLevel.index)--$(pathway.finalLevel.index) ...")
    end
    #
    finalLevel        = deepcopy(pathway.finalLevel)
    intermediateLevel = deepcopy(pathway.intermediateLevel)
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    
    newpChannels = PhotoEmission.Channel[];    rateC = 0.;    rateB = 0.
    for pChannel in pathway.photonChannels
        amplitude   = PhotoEmission.amplitude("emission", pChannel.multipole, pChannel.gauge, pathway.photonEnergy, 
                                                finalLevel, intermediateLevel, grid, display=false, printout=false)
        newpChannel = PhotoEmission.Channel( pChannel.multipole, pChannel.gauge, amplitude)
        push!( newpChannels, newpChannel)
        if       newpChannel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   newpChannel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   newpChannel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    captureRate     = 2pi * rateA
    wa              = 8.0pi * Defaults.getDefaults("alpha") * pathway.photonEnergy / (Basics.twice(pathway.intermediateLevel.J) + 1) * 
                                                                                        (Basics.twice(pathway.finalLevel.J) + 1)
    wa              = wa / pi  ## modified for test with Xe^53+ (March/2024)
    photonRate      = EmProperty(wa * rateC, wa * rateB)  
    angularBeta     = EmProperty(-9., -9.)
    #  Factor due to UserGuide
    wa              = Defaults.convertUnits("kinetic energy to wave number: atomic units", pathway.electronEnergy)
    wa              = pi*pi / (wa*wa) * captureRate  * ##  2 * # factor 2 is not really clear.
                        ((Basics.twice(pathway.intermediateLevel.J) + 1) / (Basics.twice(pathway.initialLevel.J) + 1))
    #  Factor due to Tu et al. (Plasma Phys., 2016)
    ## wa              = pi*pi / 2. / pathway.electronEnergy * captureRate * 
    ##                   (Basics.twice(pathway.intermediateLevel.J) + 1) / (Basics.twice(pathway.initialLevel.J) + 1)
    reducedStrength = EmProperty(wa * photonRate.Coulomb, wa * photonRate.Babushkin)
    pathway = DielectronicRecombination.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, 
                                                 pathway.electronEnergy, pathway.photonEnergy, captureRate, photonRate, 
                                                 angularBeta, reducedStrength, newcChannels, newpChannels)
    return( pathway )
end


"""
`DielectronicRecombination.computePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                           nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment, 
                                           settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all resonances (resonance lines) directly from the given multiplets of the initial-, intermediate- 
        and final states. It also enables one to (successively) include a set of corrections to the resonance strength to incorporate
        the contributions of shells that were not considered explicitly. 
        A list of resonances::Array{DielectronicRecombination.Resonance,1} is returned.
"""
function  computePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                          nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment,
                          settings::DielectronicRecombination.Settings)
    
    passages = DielectronicRecombination.determinePassages(intermediateMultiplet, initialMultiplet, empTreatment, settings)
    # Display all selected resonances before the computations start
    if  settings.printBefore    DielectronicRecombination.displayPassages(stdout, passages)           end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  passage in passages   maxEnergy = max(maxEnergy, passage.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    # Calculate all amplitudes and requested properties; simply copy if the captureChannels have been computed before
    # Here, the selected set of "corrections" can also be considered for each passage.
    #
    if Distributed.nworkers() > 1
        # Distributed loop
        newPassages_ = @showprogress desc="Computing Passages ..." pmap(passage -> 
                        DielectronicRecombination.computeAmplitudesProperties(passage, finalMultiplet, nm, grid, nrContinuum, empTreatment, settings), passages)
        newPassages = convert(Vector{DielectronicRecombination.Passage}, newPassages_)
    else
    # Multithreading loop
        localPassages = [DielectronicRecombination.Passage[] for _ in 1:nthreads()]
        @threads for  p in eachindex(passages)
            newPassage = DielectronicRecombination.computeAmplitudesProperties(passages[p], finalMultiplet, nm, grid, nrContinuum, 
                                                                            empTreatment, settings) 
            push!(localPassages[threadid()], newPassage)
        end 
        newPassages = vcat(localPassages...)
    end
    #
    # Original implementation
    # newPassages = DielectronicRecombination.Passage[]; 
    # for  passage in passages
    #     newPassage = DielectronicRecombination.computeAmplitudesProperties(passage, finalMultiplet, nm, grid, nrContinuum, 
    #                                                                        empTreatment, settings) 
    #     push!( newPassages, newPassage)
    # end 

    # Add empirical passages to newPassages, if requested as correction
    if  empTreatment.doEmpiricalCorrections   &&   empTreatment.nUpperEmpirical > 0
        DielectronicRecombination.addEmpiricalPassages!(newPassages, empTreatment)      end
    # 
    # Calculate all corresponding resonance
    resonances = DielectronicRecombination.computeResonances(newPassages, settings)
    # Print all results to screen
    DielectronicRecombination.displayResults(stdout, resonances,  settings)
    DielectronicRecombination.displayRateCoefficients(stdout, resonances,  settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   DielectronicRecombination.displayResults(iostream, resonances,  settings)
                       DielectronicRecombination.displayRateCoefficients(iostream, resonances,  settings)    end
                
    return( newPassages )
end


"""
`DielectronicRecombination.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                           nm::Nuclear.Model, grid::Radial.Grid, settings::DielectronicRecombination.Settings; output=true)`  
    ... to compute the dielectronic recombination amplitudes and all properties as requested by the given settings. 
        A list of pathways::Array{DielectronicRecombination.Pathway,1} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                          nm::Nuclear.Model, grid::Radial.Grid, settings::DielectronicRecombination.Settings; output=true)
    println("")
    printstyled("DielectronicRecombination.computePathways(): The computation of dielectronic resonance strength, etc. starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", 
                color=:light_green)
    println("")
    # First check that the initial, intermediate and final-state multiplets are consistent with each other to allow all requested computations
    DielectronicRecombination.checkConsistentMultiplets(finalMultiplet, intermediateMultiplet, initialMultiplet)
    # Second, analyze that the high nl orbitals are properly represented with the given grid
    DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet, intermediateMultiplet, initialMultiplet)
    # Define all parameters for the empirical treatment of the ..
    empTreatment = DielectronicRecombination.determineEmpiricalTreatment(finalMultiplet, intermediateMultiplet, 
                                                                         nm, initialMultiplet, settings)
    #
    #
    if     settings.calcHyperfineResolved  &&  settings.calcOnlyPassages
        # If settings.calcHyperfineResolved, the calculations are performed throughout with hyperfine-resolved levels and by collecting
        # the contributions of individual resonances directly in terms of hyperfine-resolved passages. This follows similar concepts
        # as if high-n Rydberg shells are involved in the intermediate and final levels, and if additional (semi-empirical) corrections 
        # need to be considered.
        hfPassages = DielectronicRecombination.computeHyperfinePassages(finalMultiplet, intermediateMultiplet, initialMultiplet, nm, 
                                                                        grid, empTreatment, settings)
        if    output    return( hfPassages )
        else            return( nothing )
        end
    elseif settings.calcOnlyPassages
        # If settings.calcOnlyPassages, only the contributions of individual resonances are collected and the results
        # are returned directly. This branch if usually considered if high-n Rydberg shells are involved in the intermediate and 
        # final levels, and if additional (semi-empirical) corrections need to be considered.
        passages = DielectronicRecombination.computePassages(finalMultiplet, intermediateMultiplet, initialMultiplet, nm, 
                                                             grid, empTreatment, settings)
        if    output    return( passages )
        else            return( nothing )
        end
    end
    #
    #
    if  empTreatment.doEmpiricalCorrections
        @warn("EmpiricalCorrections() are not yet properly implemented for DR pathway computations.")       end
    if  empTreatment.doHydrogenicCorrections
        @warn("HydrogenicCorrections() are not yet properly implemented for DR pathway computations.")      end
    if  empTreatment.doMaximumlCorrection
        @warn("MaximumlCorrection() are not yet properly implemented for DR pathway computations.")         end
    if  empTreatment.doResonanceWindowCorrection
        @warn("ResonanceWindowCorrection() are not yet properly implemented for DR pathway computations.")  end
    #
    pathways = DielectronicRecombination.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    # Display all selected pathways before the computations start
    if  settings.printBefore    DielectronicRecombination.displayPathways(stdout, pathways)    end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    # Calculate all amplitudes and requested properties; simply copy if the captureChannels have been computed before
    newPathways = DielectronicRecombination.Pathway[];    
    hasCaptureChannels = false;    lastCaptureIndices = (0,0);   lastCaptureChannels = AutoIonization.Channel[]
    for  pathway in pathways
        if  (pathway.initialLevel.index, pathway.intermediateLevel.index) == lastCaptureIndices   hasCaptureChannels = true
        else                                                                                      hasCaptureChannels = false
        end
        newPathway = DielectronicRecombination.computeAmplitudesProperties(pathway, nm, grid, nrContinuum, settings, 
                                                                           hasCaptureChannels, lastCaptureChannels) 
        push!( newPathways, newPathway)
        #
        lastCaptureIndices  = (newPathway.initialLevel.index, newPathway.intermediateLevel.index)
        lastCaptureChannels = newPathway.captureChannels
    end
    # 
    # Calculate all corresponding resonance
    resonances = DielectronicRecombination.computeResonances(newPathways, settings)
    # Print all results to screen
    DielectronicRecombination.displayResults(stdout, newPathways, settings)
    DielectronicRecombination.displayResults(stdout, resonances,  settings)
    DielectronicRecombination.displayRateCoefficients(stdout, resonances,  settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   DielectronicRecombination.displayResults(iostream, newPathways, settings)
                       DielectronicRecombination.displayResults(iostream, resonances,  settings)
                       DielectronicRecombination.displayRateCoefficients(iostream, resonances,  settings)    end
    #
    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`DielectronicRecombination.computeRateCoefficient(resonance::DielectronicRecombination.Resonance, temp::Float64)`  
    ... computes for a delta-like resonance the DR rate coefficient alpha_d (i, Te) from the given resonance strength
        and temperature [K], and for both, Coulomb and Babushkin gauge. All values are directly returned in [cm^3/s].
        An alphaDR::EmProperty is returned.
"""
function computeRateCoefficient(resonance::DielectronicRecombination.Resonance, temp::Float64)
    temp_au = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
    factor  = 4 / sqrt(2pi) * temp_au^(-3/2) * resonance.resonanceEnergy * exp(- resonance.resonanceEnergy/temp_au)
    alphaDR = factor * resonance.resonanceStrength 
    # Convert units into cm^3 / s
    factor  = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 * 
                Defaults.convertUnits("rate: from atomic to 1/s", 1.0) 
    alphaDR = factor * alphaDR
                
    return( alphaDR )
end


"""
`DielectronicRecombination.computeResonances(passages::Array{DielectronicRecombination.Passage,1}, 
                                             settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all resonances (resonance lines) as defined by the given passages and and settings. 
        In fact, passages and resonances are treated rather similar to each other in JAC and can be readily
        extracted from each other. A list of resonances::Array{DielectronicRecombination.Resonance,1} is returned.
"""
function  computeResonances(passages::Array{DielectronicRecombination.Passage,1}, settings::DielectronicRecombination.Settings)
    resonances = DielectronicRecombination.Resonance[]
    for  passage in passages
        augerRate = 0.;   
        for  ps in passages   
            if  ps.intermediateLevel.index == passage.intermediateLevel.index   augerRate = augerRate + ps.captureRate   end   
        end
        resonanceStrength = passage.reducedStrength / (augerRate + passage.photonRate)
        #
        if (isequal(resonanceStrength.Babushkin, NaN) || isequal(resonanceStrength.Coulomb, NaN)) continue end
        #
        push!( resonances, DielectronicRecombination.Resonance(passage.initialLevel, passage.intermediateLevel, 
                                                               passage.electronEnergy, resonanceStrength, 
                                                               passage.captureRate, augerRate, passage.photonRate) )
    end
    
    return( resonances )
end


"""
` (pathways::Array{DielectronicRecombination.Pathway,1}, settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all resonances (resonance lines) as defined by the given pathways and and settings. 
        A list of resonances::Array{DielectronicRecombination.Resonance,1} is returned.
"""
function  computeResonances(pathways::Array{DielectronicRecombination.Pathway,1}, settings::DielectronicRecombination.Settings)
    # Determine all pre-defined resonances in pathways
    resonances = DielectronicRecombination.Resonance[];    idxTuples = Tuple{Int64,Int64}[]
    for  pathway in pathways    idxTuples = union(idxTuples, [(pathway.initialLevel.index, pathway.intermediateLevel.index)] )    end
    # Determine the capture rate for each idxTuple
    captureRates = zeros( length(idxTuples) )
    for  (idx, idxTuple)  in  enumerate(idxTuples)
        for  pathway in pathways    
            if  idxTuple == (pathway.initialLevel.index, pathway.intermediateLevel.index)
                captureRates[idx] = pathway.captureRate
            end
        end
    end
    #
    # Determine the resonances
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    for  idxTuple in idxTuples
        iLevel = Level();    nLevel = Level();    resonanceEnergy = 0.;    reducedStrength = EmProperty(0., 0.)
        captureRate = 0.;    augerRate = 0.;      photonRate =  EmProperty(0., 0.)
        for  pathway in pathways    
            if  idxTuple == (pathway.initialLevel.index, pathway.intermediateLevel.index)
                iLevel            = pathway.initialLevel
                nLevel            = pathway.intermediateLevel
                resonanceEnergy   = pathway.intermediateLevel.energy - pathway.initialLevel.energy + electronEnergyShift
                captureRate       = pathway.captureRate 
                photonRate        = photonRate + pathway.photonRate
                reducedStrength   = reducedStrength + pathway.reducedStrength
            end
        end
        # Determine total Auger rates by summation over all idxTupels with the same m
        for  (idx, idxxTuple)  in  enumerate(idxTuples)
            if idxTuple[2] == idxxTuple[2]   augerRate = captureRates[idx]   end
        end
        #
        #== This code need to be adapted: Correct the photon rate if requested
        if  typeof(settings.corrections) == DielectronicRecombination.HydrogenicCorrections
            println(">>> Add hydrogenic corrections from n_low = $(settings.corrections.nDetailed+1) upwards")
            nDetailed     = settings.corrections.nDetailed
            Zeff          = settings.corrections.effectiveZ
            # Determine ni, li for the given resonance
            rydbergSubshs = Basics.extractRydbergSubshellList(nLevel, nDetailed+1, 1.0e-1)
            rydbergShells = Basics.extractNonrelativisticShellList(rydbergSubshs)
            if  length(rydbergShells) == 1  rShell = rydbergShells[1];   ni = rShell.n;   li = rShell.l
            else   error("Inappropriate number of Rydberg shells = $rydbergShells ")
            end
            # Compute and add hydrogenic rates A(ni,li --> nDetailed < n <= ni-1, li +- 1)
            hydrogenicRate = 0.
            for  nf = nDetailed+1:ni-1
                hydrogenicRate = hydrogenicRate + computeHydrogenicRate(ni, li, nf, li-1, Zeff) + 
                                                    computeHydrogenicRate(ni, li, nf, li+1, Zeff)
            end
            photonRate = photonRate + hydrogenicRate * settings.corrections.rateScaling
            println(">>> Total photon rate & total hydrogenic rate for A($ni, $li;  -> ...) = $photonRate, " * 
                    "$(hydrogenicRate * settings.corrections.rateScaling)")

        end  ==#
        
        resonanceStrength = reducedStrength / (augerRate + photonRate)
        push!( resonances, DielectronicRecombination.Resonance( iLevel, nLevel, resonanceEnergy, resonanceStrength, captureRate, augerRate, photonRate) )
    end
    
    return( resonances )
end


"""
`DielectronicRecombination.computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64,  Zeff::Float64)`
    ... to compute the non-relativistic electric-dipole rate for the transition from shell ni,li --> nf,lf of a
        hydrogenic ion with effective charge Zeff. The recursion formulas by Infeld and Hull (1951) are used
        together with the absorption oscillator strength. This makes the overall formulation/computation rather
        obscure, unfortunately. Uses SpecialFunctions.logfactorial. A rate::Float64 [a.u.] is returned.
        This procedure has been worked out by Stefan Schippers (2023).
"""
function  computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64,  Zeff::Float64)
    # Compute A(n,l) coeffient in the recursion formulas
    function computeA(n::Int64, l::Int64)
        A = 0.0
        if n>l  &&  n*l > 0     A = sqrt( (n+l) * (n-l)) / (n*l)     end
        return ( A )
    end
    # Compute I(n,l; n',l') integral in the recursion formulas
    function computeI(n::Int64, l::Int64, np::Int64, lp::Int64)
        wi = 0.0
        if       l>=n  ||  l<0  || lp>n  ||  lp>=np  ||  lp<0  || abs(l-lp)!=1    wi = 0.0
        elseif   l == lp-1
            if  lp == n
                wi = (n+2)*log(4*n*np)+(np-n-2)*log(np-n)-(np+n+2)*log(np+n) + 
                        0.5 * ( SpecialFunctions.logfactorial(np+n) - SpecialFunctions.logfactorial(np-n-1) -
                                SpecialFunctions.logfactorial(2*n-1) )
                wi = 0.25* exp(wi)
            else
                wi = (2*lp+1)   * computeA(np, lp+1) * computeI(n,lp, np,lp+1) + computeA(n, lp+1) * computeI(n,lp+1, np,lp)
                wi = wi / (2*lp * computeA(n, lp))
            end
        elseif   l == lp+1
            wi = computeA(np, l+1) * computeI(n,l, np,l+1) + (2*l+1) * computeA(n, l+1) * computeI(n,l+1, np,l)
            wi = wi / (2*l * computeA(np, l))
        else
            error("Unexpected set of quantum number n=$n l=$l  np=$np  lp=$lp ")
        end
        return ( wi )
    end
    # Compute absorption oscillator strength
    function computeOsc(n::Int64, l::Int64, np::Int64, lp::Int64)
        # This oscillator strength is used in absorption
        wx = 0.0
        if       l == lp+1   wx = (1/n^2 - 1/np^2) * (lp+1) / (2*(lp+1) +1) * computeI(n,lp+1, np,lp)^2
        elseif   l == lp-1   wx = (1/n^2 - 1/np^2) *  lp    / (2*(lp-1) +1) * computeI(n,lp-1, np,lp)^2
        else     error("stop a")
        end
        return( wx / 3.0 )
    end
    #
    rate = 0.0;
    if  abs(li-lf)!=1  ||  ni <= nf  ||  li >= ni  ||  li<0  ||  lf >= nf  ||  lf<0   return( rate )   
    elseif  nf > 40    error("Don't use a recursive scheme ... but make a new implementation for ni = $ni ")
    end
    #
    if      lf == (li + 1)     rate = (2*(li+1) + 1) / (2*li+1) * computeOsc(nf, li+1, ni, li)
    elseif  lf == (li - 1)     rate = (2*(li-1) + 1) / (2*li+1) * computeOsc(nf, li-1, ni, li)
    else    error("stop b")
    end
    rate = rate * Zeff^4 / 2.0 * Defaults.getDefaults("alpha")^3 * (1/nf^2 - 1/ni^2)^2
    ## println(">>>> Hydrogenic rate for A($ni, $li;  -> $nf, $lf) = $rate")

    return( rate )
end


"""
`DielectronicRecombination.determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, 
                                                    settings::DielectronicRecombination.Settings)` 
    ... to determine a list of AutoIonization.Channel for a (Auger) capture transitions from the initial to an 
        intermediate level, and by taking into account the particular settings of for this computation;  
        an Array{AutoIonization.Channel,1} is returned.
"""
function determineCaptureChannels(intermediateLevel::Level, initialLevel::Level, settings::DielectronicRecombination.Settings)
    channels = AutoIonization.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symn)
    for  kappa in kappaList
        push!( channels, AutoIonization.Channel(kappa, symn, 0., Complex(0.)) )
    end

    return( channels )  
end


"""
`DielectronicRecombination.determineEmpiricalTreatment(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                    nm::Nuclear.Model, initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)` 
    ... to determine an instance of empiricalTreatment::EmpiricalTreatment that is (internally) applied to simplify
        the use of "corrections" to the DR strenghts. This data structure summarizes all parameters that help introduce
        several empirical corrections. The procedure is simple but slightly sophisticated as we wish to support "missing"
        parameters in the individual corrections as well as the knowledge that can be derived internally. 
        The definition of EmpiricalTreatment() can readily be extended as the need arises from the user side.
"""
function determineEmpiricalTreatment(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                  nm::Nuclear.Model, initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)
    # First specify all parameters in turn 
    doEmpiricalCorrections = doHydrogenicCorrections = doMaximumlCorrection = doResonanceWindowCorrection = false
    nHydrogenic            = nUpperEmpirical         = 0
    maximum_l              = 1000
    empiricalRateScaling   = hydrogenicEffectiveZ = hydrogenicRateScaling = 1.0
    resonanceEnergyMin     = 0.;    resonanceEnergyMax = 1.0e8;    empiricalEffectiveZ = 0.
    
    for correction in settings.corrections
        if  typeof(correction) == EmpiricalCorrections    
            doEmpiricalCorrections  = true;    empiricalEffectiveZ  = correction.effectiveZ
                                               empiricalRateScaling = correction.rateScaling
                                               nUpperEmpirical      = correction.nUpperEmpirical
        end
        if  typeof(correction) == HydrogenicCorrections   
            doHydrogenicCorrections = true;    nHydrogenic            = correction.nHydrogenic
                                               hydrogenicEffectiveZ   = correction.effectiveZ
                                               hydrogenicRateScaling  = correction.rateScaling  
        end
        if  typeof(correction) == MaximumlCorrection      
            doMaximumlCorrection = true;       maximum_l = correction.maximum_l     
        end
        if  typeof(correction) == ResonanceWindowCorrection   
            doResonanceWindowCorrection = true;     resonanceEnergyMin = correction.energyMin    
                                                    resonanceEnergyMax = correction.energyMax  
        end
    end
    
    # Extract the n-shell quantum numbers from the various multiplets
    NoCoreElectrons          = initialMultiplet.levels[1].basis.NoElectrons
    basis                    = intermediateMultiplet.levels[1].basis
    intermediateConfs        = Basics.extractConfigurations(Basics.FromBasis(), basis)
    intermediateShellList    = Basics.extractNonrelativisticShellList(intermediateMultiplet.levels[1].basis.subshells)
    finalShellList           = Basics.extractNonrelativisticShellList(finalMultiplet.levels[1].basis.subshells)
    ## Basics.extractNonrelativisticShellList(finalMultiplet)
    @show NoCoreElectrons, intermediateShellList, finalShellList
    
    intermediateNs = Int64[];   for  shell  in  intermediateShellList    push!(intermediateNs, shell.n)   end
    finalNs        = Int64[];   for  shell  in  finalShellList           push!(finalNs,        shell.n)   end
    coreNs         = Int64[]
    for  conf  in  intermediateConfs     ne = 0;   occ = 0
        for  shell in intermediateShellList
            if haskey(conf.shells, shell)   
                ne = ne + conf.shells[shell]   
                ne = ne + occ;    push!(coreNs, shell.n)
            end
            if  ne >= NoCoreElectrons    break   end
        end
    end
    #
    nCore          = maximum(coreNs)
    nFinal         = maximum(finalNs)
    if  typeof(nHydrogenic)     == Missing    nHydrogenic = 0   end
    afterGap = false;   nLowestCaptured = 0
    for  i = 1:100    
        if     i in intermediateNs && afterGap   nLowestCaptured = i;   break
        elseif !(i in intermediateNs)            afterGap = true
        end 
    end
    
    if  typeof(nUpperEmpirical)         == Missing    nLowerEmpirical = 0;    nUpperEmpirical = 0  
    else                                              nLowerEmpirical = maximum(intermediateNs) + 1      end
    if  typeof(empiricalEffectiveZ)     == Missing    empiricalEffectiveZ   = nm.Z - nCore               end
    if  typeof(empiricalRateScaling)    == Missing    empiricalRateScaling  = 1.0                        end
    if  typeof(hydrogenicEffectiveZ)    == Missing    hydrogenicEffectiveZ  = nm.Z - nCore               end
    if  typeof(hydrogenicRateScaling)   == Missing    hydrogenicRateScaling = 1.0                        end
    if  typeof(maximum_l)               == Missing    maximum_l             = 1                          end
    if  typeof(resonanceEnergyMin)      == Missing    resonanceEnergyMin    = 0.                         end
    if  typeof(resonanceEnergyMax)      == Missing    resonanceEnergyMax    = 1.0e7                      end
    
   
    empTreatment = DielectronicRecombination.EmpiricalTreatment(doEmpiricalCorrections, doHydrogenicCorrections,
                                             doMaximumlCorrection, doResonanceWindowCorrection, nCore, nFinal, nHydrogenic, 
                                             nLowestCaptured, nLowerEmpirical, nUpperEmpirical, maximum_l, hydrogenicEffectiveZ, 
                                             hydrogenicRateScaling, empiricalEffectiveZ, empiricalRateScaling,
                                             resonanceEnergyMin, resonanceEnergyMax)
    
    println("  > An EmpiricalTreatment is defined with the following parameters:\n$empTreatment")
    
    @show length(intermediateMultiplet.levels)
    
    return( empTreatment )
end


"""
`DielectronicRecombination.determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, 
                                                   settings::DielectronicRecombination.Settings)` 
    ... to determine a list of PhotoEmission.Channel for the photon transitions from the intermediate and to a final level, and by 
        taking into account the particular settings of for this computation;  an Array{PhotoEmission.Channel,1} is returned.
"""
function determinePhotonChannels(finalLevel::Level, intermediateLevel::Level, settings::DielectronicRecombination.Settings)
    channels = PhotoEmission.Channel[];   
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symn, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)           push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end

    return( channels )  
end


"""
`DielectronicRecombination.determinePassages(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                             empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)`  
    ... to determine a list of dielectronic-recombination resonances between the levels from the given initial- and intermediate- 
        states, whereas the final states are considered "on-fly"; the particular selections and settings for this computation
        are taken into account; an Array{DielectronicRecombination.Passsage,1} is returned. Apart from the level specification, 
        all physical properties are set to zero during the initialization process.  
"""
function  determinePassages(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                            empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)
    passages = DielectronicRecombination.Passage[]
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    @warn("No pathway selection is considered, if settings.calcOnlyPassages=true.")
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            eEnergy = nLevel.energy - iLevel.energy + electronEnergyShift
            if  eEnergy < 0.  ||   eEnergy <  empTreatment.resonanceEnergyMin  ||   
                                   eEnergy >  empTreatment.resonanceEnergyMax  continue    end
            cChannels = DielectronicRecombination.determineCaptureChannels(nLevel, iLevel, settings) 
            push!( passages, DielectronicRecombination.Passage(iLevel, nLevel, eEnergy, 0., EmProperty(0., 0.), 
                                                                EmProperty(0., 0.), cChannels) )
        end
    end
    return( passages )
end


"""
`DielectronicRecombination.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                             settings::DielectronicRecombination.Settings)`  
    ... to determine a list of dielectronic-recombination pathways between the levels from the given initial-, intermediate- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{DielectronicRecombination.Pathway,1} is returned. Apart from the level specification, all physical properties 
        are set to zero during the initialization process.  
"""
function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                            settings::DielectronicRecombination.Settings)
    pathways = DielectronicRecombination.Pathway[]
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    photonEnergyShift   = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                    eEnergy = nLevel.energy - iLevel.energy + electronEnergyShift
                    pEnergy = nLevel.energy - fLevel.energy + photonEnergyShift
                    if  pEnergy < 0.   ||   eEnergy < 0.    continue    end
                    cChannels = DielectronicRecombination.determineCaptureChannels(nLevel, iLevel, settings) 
                    pChannels = DielectronicRecombination.determinePhotonChannels( fLevel, nLevel, settings) 
                    push!( pathways, DielectronicRecombination.Pathway(iLevel, nLevel, fLevel, eEnergy, pEnergy, 0., 
                                                               EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), 
                                                               cChannels, pChannels) )
                end
            end
        end
    end
    return( pathways )
end


"""
`DielectronicRecombination.displayPassages(stream::IO, passages::Array{DielectronicRecombination.Passage,1})`  
    ... to display a list of passages and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPassages(stream::IO, passages::Array{DielectronicRecombination.Passage,1})
    nx = 120
    println(stream, " ")
    println(stream, "  Selected dielectronic-recombination passages:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(16, "Levels"; na=4);            sb = sb * TableStrings.center(16, "i  --  m"; na=4);          
    sa = sa * TableStrings.center(16, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(16, "i  --  m"; na=3);
    sa = sa * TableStrings.center(18, "Energies  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(18, "electron     "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  passage in passages
        sa  = "  ";     isym = LevelSymmetry( passage.initialLevel.J,      passage.initialLevel.parity)
                        msym = LevelSymmetry( passage.intermediateLevel.J, passage.intermediateLevel.parity)
        sa = sa * TableStrings.center(17, TableStrings.levels_if(passage.initialLevel.index, passage.intermediateLevel.index); na=7)
        sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", passage.electronEnergy)) * "        "
        kappaSymmetryList = Tuple{Int64,LevelSymmetry}[]
        for  cChannel in passage.captureChannels
            push!( kappaSymmetryList, (cChannel.kappa, cChannel.symmetry) )
        end
        wa = TableStrings.kappaSymmetryTupelList(85, kappaSymmetryList)
        if  length(wa) > 0    sb = sa * wa[1];    println(stream,  sb )    end  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "\n>> A total of $(length(passages)) dielectronic-recombination passages will be calculated. \n")
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayPathways(stream::IO, pathways::Array{DielectronicRecombination.Pathway,1})`  
    ... to display a list of pathways and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPathways(stream::IO, pathways::Array{DielectronicRecombination.Pathway,1})
    nx = 180
    println(stream, " ")
    println(stream, "  Selected dielectronic-recombination pathways:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(23, "Levels"; na=4);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=4);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=3);
    sa = sa * TableStrings.center(26, "Energies  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(26, "electron        photon "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=5)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  cChannel in pathway.captureChannels
            for  pChannel in pathway.photonChannels
                push!( kappaMultipoleSymmetryList, (cChannel.kappa, pChannel.multipole, pChannel.gauge, cChannel.symmetry) )
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        if  length(wa) > 0    sb = sa * wa[1];    println(stream,  sb )    end  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "\n>> A total of $(length(pathways)) dielectronic-recombination pathways will be calculated. \n")
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayResults(stream::IO, pathways::Array{DielectronicRecombination.Pathway,1},
                                          settings::DielectronicRecombination.Settings)`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, pathways::Array{DielectronicRecombination.Pathway,1}, settings::DielectronicRecombination.Settings)
    function  computeTotalRate(mLevel::Level)
        tRate = EmProperty(0.);   captureRateAdded = false 
        # Add all capture and photon rate for the given level
        for  pathway in pathways
            if  mLevel.index == pathway.intermediateLevel.index   
                if  captureRateAdded              tRate = tRate + pathway.photonRate   
                else   captureRateAdded = true;   tRate = tRate + pathway.captureRate + pathway.photonRate   
                end
            end
        end
        
        return( tRate )
    end
    
    nx = 150
    println(stream, " ")
    println(stream, "  Partial (Auger) capture and radiative decay rates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);
    sa = sa * TableStrings.center(38, "Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(38, "electron        m--i       photon  "; na=1)
    sa = sa * TableStrings.center(10, "Multipoles"; na=6);        sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(36, "Rates  " * TableStrings.inUnits("rate"); na=2);   
    sb = sb * TableStrings.center(36, "(Auger) capture    Cou--photon--Bab";        na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        if  pathway.photonRate.Coulomb == 0.  &&   pathway.photonRate.Babushkin == 0.   continue    end
        sa  = " ";      isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", en_mi))                  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        multipoles = EmMultipole[]
        for  pch in pathway.photonChannels
            multipoles = push!( multipoles, pch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)     * "                   "
        sa = sa * TableStrings.flushleft(16, mpString[1:15];  na=0)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.captureRate))            * "     "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.photonRate.Coulomb))     * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pathway.photonRate.Babushkin))   * "  "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    #
    nx = 170
    println(stream, " ")
    println(stream, "  Partial resonance strength:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);
    sa = sa * TableStrings.center(38, "Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(38, "electron        m--i       photon  "; na=1)
    sa = sa * TableStrings.center(10, "Multipoles"; na=5);        sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(26, "S * Gamma_m  " * TableStrings.inUnits("reduced strength"); na=2);   
    sb = sb * TableStrings.center(26, " Cou -- photon -- Bab";        na=5)
    sa = sa * TableStrings.center(26, "   Width_m   " * TableStrings.inUnits("energy"); na=2);   
    sb = sb * TableStrings.center(26, " Cou --  width  -- Bab";       na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        if  pathway.reducedStrength.Coulomb == 0.  &&   pathway.reducedStrength.Babushkin == 0.     continue   end
        sa  = " ";      isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", en_mi))                  * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "    "
        multipoles = EmMultipole[]
        for  pch in pathway.photonChannels
            multipoles = push!( multipoles, pch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)                          * "                   "
        sa = sa * TableStrings.flushleft(15, mpString[1:15];  na=0)
        wa = Defaults.convertUnits("strength: from atomic", 1.0) * Defaults.convertUnits("energy: from atomic", 1.0)
        sa = sa * @sprintf("%.4e", wa * pathway.reducedStrength.Coulomb)     * "     "
        sa = sa * @sprintf("%.4e", wa * pathway.reducedStrength.Babushkin)   * "     "
        wb = computeTotalRate(pathway.intermediateLevel)
        wa = Defaults.convertUnits("energy: from atomic", 1.0)
        sa = sa * @sprintf("%.4e", wa * wb.Coulomb)     * "     "
        sa = sa * @sprintf("%.4e", wa * wb.Babushkin)   * "     "
        
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
    + (stream::IO, resonances::Array{DielectronicRecombination.Resonance,1}, settings::DielectronicRecombination.Settings)`  
    ... to list all results for the resonances. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, resonances::Array{DielectronicRecombination.Resonance,1}, settings::DielectronicRecombination.Settings)
    nx = 160
    println(stream, " ")
    println(stream, "  Total Auger rates, radiative rates and resonance strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(42, "Auger rate     Cou -- rad. rates -- Bab"; na=1);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=1)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=6)
    sa = sa * TableStrings.center(30, "Cou -- res. strength -- Bab"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=2)
    sa = sa * TableStrings.center(18, "Widths Gamma_m"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=6)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                       msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))          * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.augerRate))                  * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Coulomb))         * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Babushkin))       * "        "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Coulomb))    * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Babushkin))  * "     "
        wa = resonance.augerRate + resonance.photonRate.Coulomb 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa))                                 * "   "
        wa = resonance.augerRate + resonance.photonRate.Babushkin 
        sa = sa * "(" * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa)) * ")"                     * "   "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayRateCoefficients(stream::IO, resonances::Array{DielectronicRecombination.Resonance,1},
                                                   settings::DielectronicRecombination.Settings)`  
    ... to list, if settings.calcRateAlpha, all rate coefficients for the selected temperatures. Both, the individual as well as
        the total DR plasma rate coefficients are printed in neat tables, though nothing is returned otherwise.
"""
function  displayRateCoefficients(stream::IO, resonances::Array{DielectronicRecombination.Resonance,1},
                                  settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures)
    if  !settings.calcRateAlpha  ||  ntemps == 0     return(nothing)     end
    #
    nx = 54 + 17 * min(ntemps, 7)
    println(stream, " ")
    println(stream, "  Rate coefficients for delta-like resonances [cm^3/s]:        ... all results in Babushkin gauge")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    for  nt = 1:min(ntemps, 7)
        sa = sa * TableStrings.center(14, "T = " * @sprintf("%.2e", settings.temperatures[nt]); na=3);       
        sb = sb * TableStrings.center(14, "[K]"; na=3)
    end
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                        msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))    * "      "
        for  nt = 1:min(ntemps, 7)
            alphaDR      = DielectronicRecombination.computeRateCoefficient(resonance, settings.temperatures[nt])
            sa = sa * @sprintf("%.4e", alphaDR.Babushkin)  * "       "
        end
        println(stream, sa)
    end
    #
    println(stream, "  ")
    sa = "       alpha^DR (T, i; Coulomb gauge):                      " 
    sb = "       alpha^DR (T, i; Babushkin gauge):                    " 
    for  nt = 1:min(ntemps, 7) 
        alphaDRtotal = EmProperty(0.)
        for  resonance in resonances
            alphaDR      = DielectronicRecombination.computeRateCoefficient(resonance, settings.temperatures[nt])
            alphaDRtotal = alphaDRtotal + alphaDR
        end
        sa = sa * @sprintf("%.4e", alphaDRtotal.Coulomb)    * "       "
        sb = sb * @sprintf("%.4e", alphaDRtotal.Babushkin)  * "       "
    end
    println(stream, sa);    println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`DielectronicRecombination.extractRateCoefficients(resonances::Array{DielectronicRecombination.Resonance,1}, settings::DielectronicRecombination.Settings)`  
    ... to extract, if settings.calcRateAlpha, the total DR rate coefficients for all temperatures. 
        A list of total rate coefficients [cm^3/s] alphaDR::Array{EmProperty,1} is returned.
"""
function  extractRateCoefficients(resonances::Array{DielectronicRecombination.Resonance,1}, settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures);     alphaDRs = EmProperty[]
    if  !settings.calcRateAlpha  ||  ntemps == 0     return( EmProperty[] )     end
    #
    for  nt = 1:ntemps
        alphaDRtotal = EmProperty(0.)
        for  resonance in resonances
            alphaDR      = DielectronicRecombination.computeRateCoefficient(resonance, settings.temperatures[nt])
            alphaDRtotal = alphaDRtotal + alphaDR
        end
        push!(alphaDRs, alphaDRtotal)
    end
    #
    return( alphaDRs )
end


"""
`DielectronicRecombination.isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)`  
    returns true, if level is to be excluded from the valid resonances, and false otherwise.
    It returns false if the ResonanceSelection() is inactive or if level belongs to the selected resoances.
    It is true only of ResonanceSelection() is active but the level does not belong to the selected resonances.
"""
function  isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)
    wa = false
    if !rSelection.active    return( wa )
        # This is the standard case if no additional limitations are specified by the user
    else
        # Analyze of whether level belongs to the selected resonances; it determines of whether level has 
        # electrons in either toShells or intoShells, and if there is one electron less in the fromShells
        fromwb    = false;    towb    = false;    intowb    = true
        confList  = Basics.extractConfigurations(Basics.FromBasis(), level.basis)
        occShells = Basics.extractShellList(confList) 
        for  shell in rSelection.toShells   
            if  shell in occShells   towb   = true;   break    end
        end
        for  shell in occShells   
            if  !(shell in rSelection.intoShells)   intowb = false;   break    end
        end
        # Compare the occupation of the given resonance level with those of the refLevel; there should be (at least)
        # one electron more in the shells of leadingConfig than in the same shells of level
        refConfig = Basics.extractConfiguration(Basics.LeadingConfiguration(), refLevel)
        levConfig = confList[1];   NoElectrons = 0
        for  (k,v) in refConfig.shells
            NoElectrons = NoElectrons + levConfig.shells[k]
        end
        @show NoElectrons, refConfig.NoElectrons
        if  NoElectrons < refConfig.NoElectrons     fromwb = true   end
        
        wa = !(fromwb &&  towb  && intowb)
    end
    
    return( wa )
end
