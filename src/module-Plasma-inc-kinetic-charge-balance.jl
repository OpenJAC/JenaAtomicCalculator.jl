
# Functions and methods for scheme::Plasma.KineticChargeBalanceScheme computations


"""
`Plasma.generateChargeStageMultiplet(stage::Plasma.ChargeStage, nm::Nuclear.Model, grid::Radial.Grid, asfSettings::AsfSettings)`
    ... generates the small 2-level reference Multiplet of one charge stage: one shared mean-field orbital basis for
        [stage.groundConfig, stage.excitedConfig], followed by Hamiltonian.performCIwithFrozenOrbitals for each
        configuration separately and Basics.sortByEnergy over the combined pair -- the same generation policy as
        Plasma.generateCollisionalRadiativeLevels, specialized to exactly two configurations. After sorting, level 1 is
        always the lower-energy (ground) level and level 2 the higher-energy (excited) level, regardless of which
        configuration happened to be passed as groundConfig/excitedConfig; a Multiplet is returned.
"""
function generateChargeStageMultiplet(stage::Plasma.ChargeStage, nm::Nuclear.Model, grid::Radial.Grid, asfSettings::AsfSettings)
    allConfigs = [stage.groundConfig, stage.excitedConfig]
    name       = "KineticChargeBalanceScheme stage multiplet"
    mfSettings = AtomicState.MeanFieldSettings()
    repBasis   = Representation(name, nm, grid, allConfigs, MeanFieldBasis(mfSettings))
    repOutput  = generate(repBasis, output=true)
    orbitals   = repOutput[ResultKeys.MeanFieldBasis].orbitals

    levels = Level[]
    for  conf in allConfigs
        mp = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, nm, grid, asfSettings; printout=false)
        mp = Basics.sortByEnergy(mp)
        push!(levels, mp.levels[1])
    end

    return( Basics.sortByEnergy( Multiplet(name, levels) ) )
end


"""
`Plasma.structureBasedRrRateCoefficient(groundLevel::Level, higherLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                        Te::Float64)`
    ... computes the structure-based photorecombination plasma rate coefficient alpha^(RR)(T_e; higherLevel -> groundLevel)
        by a Gauss-Legendre Maxwellian average of JAC's own PhotoRecombination.computeLines cross section, mirroring
        Cascade.simulateRrRateCoefficients but as a plain-argument local helper (that function is coupled to
        Cascade.Simulation's own types and is not reused here, to keep this task to one module). Te is given in Kelvin.
        A named tuple (Coulomb::Float64, Babushkin::Float64) [cm^3/s] is returned.
"""
function structureBasedRrRateCoefficient(groundLevel::Level, higherLevel::Level, nm::Nuclear.Model, grid::Radial.Grid, Te::Float64)
    Te_au     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", Te)
    finalMp   = Multiplet("RR final",   [groundLevel])
    initialMp = Multiplet("RR initial", [higherLevel])

    gridGL     = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., 30*Te_au, 16; printout=false)
    prSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(); multipoles=[Basics.E1], gauges=[UseCoulomb, UseBabushkin],
                     electronEnergies=gridGL.t, printBefore=false)
    lines = PhotoRecombination.computeLines(finalMp, initialMp, nm, grid, prSettings; output=true)

    # Energies are matched by VALUE, not by list position: PhotoRecombination.computeLines is not documented to
    # preserve the order of settings.electronEnergies, and a single level pair makes an exact-value lookup safe.
    weightOf = Dict{Float64,Float64}( gridGL.t[n] => gridGL.wt[n]  for n = 1:gridGL.nt )

    factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    alphaC = 0.;   alphaB = 0.
    for  line  in  lines
        wt     = weightOf[line.electronEnergy]
        fe     = Distribution.electronEnergyDistribution(Distribution.ElectronMaxwell(Te_au), line.electronEnergy)
        alphaC = alphaC + line.crossSection.Coulomb   * sqrt(2*line.electronEnergy) * fe * wt
        alphaB = alphaB + line.crossSection.Babushkin * sqrt(2*line.electronEnergy) * fe * wt
    end

    return( (Coulomb = factor*alphaC, Babushkin = factor*alphaB) )
end


"""
`Plasma.displayKineticChargeBalanceTreatments(stream::IO, scheme::Plasma.KineticChargeBalanceScheme)`
    ... prints, once per leg, which Plasma.AbstractRateTreatment was used -- the user-facing provenance table that
        lets a reader recognize, at a glance, at which (empirical or structure-based) level each process was treated,
        without having to read the scheme's fields back out of the code that built it.
"""
function displayKineticChargeBalanceTreatments(stream::IO, scheme::Plasma.KineticChargeBalanceScheme)
    println(stream, "\n  Plasma.KineticChargeBalanceScheme -- rate treatment per leg:")
    println(stream, "  ", TableStrings.hLine(80))
    sa = "    " * TableStrings.flushleft(28, "Leg") * TableStrings.flushleft(13, "Within/") * TableStrings.flushleft(40, "Treatment")
    sb = "    " * TableStrings.flushleft(28, "")    * TableStrings.flushleft(13, "cross-stage")
    println(stream, sa);   println(stream, sb)
    println(stream, "  ", TableStrings.hLine(80))
    rows = [ ("PhotoEmission",          "within",      scheme.peTreatment),
             ("ImpactExcitation",       "within",      scheme.ieTreatment),
             ("ImpactIonization (up)",  "cross-stage", scheme.eiiTreatment),
             ("PhotoRecombination (dn)","cross-stage", scheme.rrTreatment),
             ("ThreeBodyRecomb. (dn)",  "cross-stage", scheme.tbrTreatment) ]
    for  (name, kind, treatment)  in  rows
        println(stream, "    " * TableStrings.flushleft(28, name) * TableStrings.flushleft(13, kind) * TableStrings.flushleft(40, string(treatment)))
    end
    println(stream, "  ", TableStrings.hLine(80))

    return( nothing )
end


"""
`Plasma.displayKineticChargeBalanceResults(stream::IO, scheme::Plasma.KineticChargeBalanceScheme, Ne::Float64, popVectors::Dict)`
    ... prints a table of the relative populations (Coulomb/Babushkin) of every (stage, ground/excited) slot, one row
        per (slot, temperature) pair, following JAC's usual table style.
"""
function displayKineticChargeBalanceResults(stream::IO, scheme::Plasma.KineticChargeBalanceScheme, Ne::Float64, popVectors::Dict)
    nx = 100
    println(stream, "\n  Kinetic charge-balance relative populations (Ne = $Ne cm^-3):")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    " * TableStrings.flushleft(20, "Stage") * TableStrings.flushleft(10, "Slot") *
                  TableStrings.flushright(14, "Temperature") * TableStrings.flushright(22, "Population") * TableStrings.flushright(22, "Population")
    sb = "    " * TableStrings.flushleft(20, "") * TableStrings.flushleft(10, "") *
                  TableStrings.flushright(14, "[K]") * TableStrings.flushright(22, "(Coulomb)") * TableStrings.flushright(22, "(Babushkin)")
    println(stream, sa);   println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    for  (k, stage)  in  enumerate(scheme.chargeStages)
        stageLabel = string(stage.groundConfig)
        for  slot  in  ("ground", "excited")
            idx = slot == "ground" ? 2k-1 : 2k
            for  temp  in  sort(collect(keys(popVectors)))
                pv = popVectors[temp]
                sc = "    " * TableStrings.flushleft(20, stageLabel) * TableStrings.flushleft(10, slot) *
                              TableStrings.flushright(14, @sprintf("%.4e", temp)) *
                              TableStrings.flushright(22, @sprintf("%.4e", pv.Coulomb[idx])) *
                              TableStrings.flushright(22, @sprintf("%.4e", pv.Babushkin[idx]))
                println(stream, sc)
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Plasma.perform(scheme::Plasma.KineticChargeBalanceScheme, computation::Plasma.Computation; output::Bool=true)`
    ... to compute the kinetic (rate-based, non-LTE) population balance across scheme.chargeStages. For output=true,
        a dictionary is returned from which the relevant results can be accessed by proper keys.
"""
function  perform(scheme::Plasma.KineticChargeBalanceScheme, computation::Plasma.Computation; output::Bool=true)
    nm  = computation.nuclearModel;   grid = computation.grid;   asf = computation.asfSettings
    nStages = length(scheme.chargeStages)
    if  nStages < 2
        error("Plasma.KineticChargeBalanceScheme needs at least 2 chargeStages to have anything to balance.")
    end
    if  !computation.settings.useNumberDensity
        error("Plasma.KineticChargeBalanceScheme: computation.settings.useNumberDensity must be true, with " *
              "settings.density interpreted as the ELECTRON density Ne [cm^-3], exactly as for Plasma.CollisionalRadiativeScheme.")
    end
    Ne = computation.settings.density

    # Empirical's plasma-rate functions read the nuclear charge from the global defaults rather than taking it as
    # an argument (cf. Cascade.simulateEiiRateCoefficients); it is set here and reliably restored afterwards, even on
    # an error, so that a KineticChargeBalanceScheme run never leaves a wrong "nuclear: charge" default behind for
    # whatever the caller does next in the same session.
    oldZ = Defaults.getDefaults("nuclear: charge")
    Defaults.setDefaults("nuclear: charge", nm.Z)
    try
        return  Plasma.computeKineticChargeBalance(scheme, nm, grid, asf, Ne, output)
    finally
        Defaults.setDefaults("nuclear: charge", oldZ)
    end
end


"""
`Plasma.computeKineticChargeBalance(scheme::Plasma.KineticChargeBalanceScheme, nm::Nuclear.Model, grid::Radial.Grid,
                                    asf::AsfSettings, Ne::Float64, output::Bool)`
    ... the body of Plasma.perform(scheme::Plasma.KineticChargeBalanceScheme, ...), factored out only so that the
        global "nuclear: charge" default can be set and reliably restored around it via a try/finally in the caller.
"""
function computeKineticChargeBalance(scheme::Plasma.KineticChargeBalanceScheme, nm::Nuclear.Model, grid::Radial.Grid,
                                     asf::AsfSettings, Ne::Float64, output::Bool)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nStages = length(scheme.chargeStages)

    @warn("Plasma.KineticChargeBalanceScheme: DielectronicRecombination is deliberately not included (see the scheme's " *
          "own docstring), and the cross-stage legs (ImpactIonization, PhotoRecombination, ThreeBodyRecombination) " *
          "connect GROUND levels of adjacent stages only -- the excited slot of a stage is populated/depopulated by " *
          "PhotoEmission/ImpactExcitation alone, never by ionization or recombination directly.")

    treatments = (scheme.peTreatment, scheme.ieTreatment, scheme.eiiTreatment, scheme.rrTreatment, scheme.tbrTreatment)
    anyStructureBased = any( t -> t isa Plasma.StructureBasedTreatment, treatments )

    stageMultiplets = Vector{Union{Nothing,Multiplet}}(nothing, nStages)
    if  anyStructureBased
        println(">> Plasma.KineticChargeBalanceScheme: at least one leg is structure-based; generating a small " *
                "reference Multiplet for every one of the $nStages charge stages.")
        for  k = 1:nStages    stageMultiplets[k] = Plasma.generateChargeStageMultiplet(scheme.chargeStages[k], nm, grid, asf)    end
    else
        println(">> Plasma.KineticChargeBalanceScheme: every leg is Plasma.EmpiricalTreatment; no SCF/CI is performed.")
    end

    Plasma.displayKineticChargeBalanceTreatments(stdout, scheme)

    N             = 2*nStages
    groundIndex(k)  = 2k-1
    excitedIndex(k) = 2k
    temperatures  = scheme.electronTemperatures
    rateFactor    = Defaults.convertUnits("rate: from atomic to 1/s", 1.0)
    coeffFactor   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)

    # -------------------------------------------------------------------------------------------------------
    # (1) Within-stage PhotoEmission, excited -> ground: temperature-independent, computed once.
    # -------------------------------------------------------------------------------------------------------
    peRatesC = zeros(nStages);   peRatesB = zeros(nStages)
    einsteinAs = Vector{Any}(undef, nStages)                                 # keep (multipole, energy, rate) for the IE detailed balance below
    for  k = 1:nStages
        stage = scheme.chargeStages[k]
        if  scheme.peTreatment isa Plasma.StructureBasedTreatment
            peSettings = PhotoEmission.Settings(PhotoEmission.Settings(); gauges=[UseCoulomb, UseBabushkin], printBefore=false)
            peLines    = PhotoEmission.computeLines(stageMultiplets[k], stageMultiplets[k], grid, peSettings; output=true)
            for  line  in  peLines
                if  line.initialLevel.index == 2  &&  line.finalLevel.index == 1
                    peRatesC[k] = rateFactor * line.photonRate.Coulomb;   peRatesB[k] = rateFactor * line.photonRate.Babushkin
                end
            end
            einsteinAs[k] = Empirical.photoemissionEinsteinA(stage.excitedConfig, stage.groundConfig, Empirical.ScaledHydrogenic(); printout=false)
        else
            ea = Empirical.photoemissionEinsteinA(stage.excitedConfig, stage.groundConfig, scheme.peTreatment.approx; printout=false)
            peRatesC[k] = peRatesB[k] = rateFactor * ea.rate
            einsteinAs[k] = ea
        end
    end

    popVectors = Dict{Float64, NamedTuple{(:Coulomb,:Babushkin), Tuple{Array{Float64,1},Array{Float64,1}}}}()

    for  TeK  in  temperatures
        Te_au   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", TeK)
        cMatrix = zeros(N, N);   bMatrix = zeros(N, N)

        # ---------------------------------------------------------------------------------------------------
        # (2) Within-stage legs: PhotoEmission (precomputed above) and ImpactExcitation (both directions).
        # ---------------------------------------------------------------------------------------------------
        for  k = 1:nStages
            stage = scheme.chargeStages[k]
            gIdx  = groundIndex(k);   eIdx = excitedIndex(k)

            if  isfinite(peRatesC[k])   cMatrix[eIdx,gIdx] += peRatesC[k];   cMatrix[eIdx,eIdx] -= peRatesC[k]   end
            if  isfinite(peRatesB[k])   bMatrix[eIdx,gIdx] += peRatesB[k];   bMatrix[eIdx,eIdx] -= peRatesB[k]   end

            if  scheme.ieTreatment isa Plasma.StructureBasedTreatment
                ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings(); calcRateCoefficient=true, printBefore=false,
                                 lineSelection=LineSelection(true; indexPairs=[(1,2)]), temperatures=[TeK])
                ieLines, _ = ImpactExcitation.computeLines(stageMultiplets[k], stageMultiplets[k], nm, grid, ieSettings; output=true)
                ieRates    = ImpactExcitation.computeRateCoefficients(ieLines, ieSettings)
                if  length(ieRates) > 0
                    rc        = ieRates[1]
                    excRate   = rc.alphas[1]
                    gUpper    = Float64(Basics.twice(rc.finalLevel.J) + 1)
                    deExcRate = 8.6291269e-6 * rc.effOmegas[1] / (gUpper * sqrt(TeK))
                else
                    excRate = deExcRate = 0.
                end
            else
                eDist   = Distribution.ElectronMaxwell(Te_au)
                excRate = coeffFactor * Empirical.impactExcitationPlasmaAlpha(eDist, stage.groundConfig, stage.excitedConfig;
                                            approx=scheme.ieTreatment.approx, aSource=Empirical.ScaledHydrogenic(), printout=false)
                # Detailed balance for the (empirical) de-excitation rate coefficient:  g_g C_exc(Te) = g_e C_deexc(Te) exp(-dE/Te).
                gi        = Basics.extractFromConfiguration(Basics.Multiplicity(), stage.groundConfig)
                ge        = Basics.extractFromConfiguration(Basics.Multiplicity(), stage.excitedConfig)
                deExcRate = (gi/ge) * excRate * exp( einsteinAs[k].energy / Te_au )
            end

            if  isfinite(excRate)
                cMatrix[gIdx,eIdx] += excRate*Ne;   cMatrix[gIdx,gIdx] -= excRate*Ne
                bMatrix[gIdx,eIdx] += excRate*Ne;   bMatrix[gIdx,gIdx] -= excRate*Ne
            end
            if  isfinite(deExcRate)
                cMatrix[eIdx,gIdx] += deExcRate*Ne;   cMatrix[eIdx,eIdx] -= deExcRate*Ne
                bMatrix[eIdx,gIdx] += deExcRate*Ne;   bMatrix[eIdx,eIdx] -= deExcRate*Ne
            end
        end

        # ---------------------------------------------------------------------------------------------------
        # (3) Cross-stage legs, ground(k) <-> ground(k+1): ImpactIonization (up), PhotoRecombination +
        #     ThreeBodyRecombination (down, summed into the same matrix entry).
        # ---------------------------------------------------------------------------------------------------
        for  k = 1:nStages-1
            stageLo = scheme.chargeStages[k];   stageHi = scheme.chargeStages[k+1]
            iLo = groundIndex(k);   iHi = groundIndex(k+1)

            if  scheme.eiiTreatment isa Plasma.StructureBasedTreatment
                groundLo = stageMultiplets[k].levels[1];   groundHi = stageMultiplets[k+1].levels[1]
                sh          = ThreeBodyRecombination.capturedSubshell(groundHi, groundLo)
                eiiSettings = ThreeBodyRecombination.Settings(ThreeBodyRecombination.Settings(); eiiModel=scheme.eiiTreatment.eiiModel,
                                                              electronTemperature=TeK, printBefore=false)
                eiiRate     = ThreeBodyRecombination.eiiRateCoefficient(sh, groundLo, nm, grid, eiiSettings)
            else
                eDist   = Distribution.ElectronMaxwell(Te_au)
                eiiRate = coeffFactor * Empirical.impactIonizationPlasmaAlpha(eDist, stageLo.groundConfig, stageHi.groundConfig;
                                            approx=scheme.eiiTreatment.approx, printout=false)
            end

            if  scheme.rrTreatment isa Plasma.StructureBasedTreatment
                groundLo = stageMultiplets[k].levels[1];   groundHi = stageMultiplets[k+1].levels[1]
                rr      = Plasma.structureBasedRrRateCoefficient(groundLo, groundHi, nm, grid, TeK)
                rrRateC = rr.Coulomb;   rrRateB = rr.Babushkin
            else
                eDist    = Distribution.ElectronMaxwell(Te_au)
                rrRateAU = coeffFactor * Empirical.photorecombinationPlasmaAlpha(eDist, stageHi.groundConfig, stageLo.groundConfig;
                                            approx=scheme.rrTreatment.approx, printout=false)
                rrRateC = rrRateB = rrRateAU
            end

            if  scheme.tbrTreatment isa Plasma.StructureBasedTreatment
                groundLo = stageMultiplets[k].levels[1];   groundHi = stageMultiplets[k+1].levels[1]
                tbrSettings = ThreeBodyRecombination.Settings(ThreeBodyRecombination.Settings(); eiiModel=scheme.tbrTreatment.eiiModel,
                                                              electronTemperature=TeK, printBefore=false)
                finalMp     = Multiplet("TBR final",   [groundLo]);   initialMp = Multiplet("TBR initial", [groundHi])
                tbrLines    = ThreeBodyRecombination.computeLines(finalMp, initialMp, nm, grid, tbrSettings; output=true)
                tbrRate     = length(tbrLines) > 0 ? tbrLines[1].rateCoefficient : 0.
            else
                eDist   = Distribution.ElectronMaxwell(Te_au)
                tbrRate = coeffFactor * Empirical.threeBodyRecombinationPlasmaAlpha(eDist, stageHi.groundConfig, stageLo.groundConfig;
                                            approx=scheme.tbrTreatment.approx, printout=false)
            end

            eiiTerm   = eiiRate * Ne
            downTermC = (rrRateC + tbrRate) * Ne;   downTermB = (rrRateB + tbrRate) * Ne

            if  isfinite(eiiTerm)
                cMatrix[iLo,iHi] += eiiTerm;   cMatrix[iLo,iLo] -= eiiTerm
                bMatrix[iLo,iHi] += eiiTerm;   bMatrix[iLo,iLo] -= eiiTerm
            end
            if  isfinite(downTermC)   cMatrix[iHi,iLo] += downTermC;   cMatrix[iHi,iHi] -= downTermC   end
            if  isfinite(downTermB)   bMatrix[iHi,iLo] += downTermB;   bMatrix[iHi,iHi] -= downTermB   end
        end

        popVectors[TeK] = (Coulomb = Plasma.solveCollisionalRadiativeBalance(cMatrix), Babushkin = Plasma.solveCollisionalRadiativeBalance(bMatrix))
    end

    Plasma.displayKineticChargeBalanceResults(stdout, scheme, Ne, popVectors)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayKineticChargeBalanceResults(iostream, scheme, Ne, popVectors)   end

    if output    results = Base.merge( results, Dict("populations:" => popVectors) )   end

    println("KineticChargeBalanceScheme computation complete ...")
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
