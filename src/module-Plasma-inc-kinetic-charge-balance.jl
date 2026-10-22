
# Functions and methods for scheme::Plasma.KineticChargeBalanceScheme computations


"""
`Plasma.generateChargeStageMultiplet(stage::Plasma.ChargeStage, nm::Nuclear.Model, grid::Radial.Grid, asfSettings::AsfSettings)`
    ... generates the reference Multiplet of one charge stage: one shared mean-field orbital basis for stage.configs,
        followed by Hamiltonian.performCIwithFrozenOrbitals for each configuration separately, retaining its own
        stage.nLevels[i] lowest CI levels (clipped if that configuration has fewer), and a single energy sort across
        the combined set -- the same generation policy as Plasma.generateCollisionalRadiativeLevels, applied once per
        configuration. A named tuple (multiplet::Multiplet=, configIndex::Array{Int64,1}=) is returned, where
        configIndex[j] gives the index into stage.configs that multiplet.levels[j] was generated from (levels sorted
        by energy, lowest first, so for a CLASSIC stage multiplet.levels[1] is always the ground level).
"""
function generateChargeStageMultiplet(stage::Plasma.ChargeStage, nm::Nuclear.Model, grid::Radial.Grid, asfSettings::AsfSettings)
    allConfigs = stage.configs
    name       = "KineticChargeBalanceScheme stage multiplet"
    mfSettings = AtomicState.MeanFieldSettings()
    repBasis   = Representation(name, nm, grid, allConfigs, MeanFieldBasis(mfSettings))
    repOutput  = generate(repBasis, output=true)
    orbitals   = repOutput[ResultKeys.MeanFieldBasis].orbitals

    levels      = Level[]
    configIndex = Int64[]
    for  (ci, conf)  in  enumerate(allConfigs)
        mp = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, nm, grid, asfSettings; printout=false)
        mp = Basics.sortByEnergy(mp)
        for  lv = 1:min(stage.nLevels[ci], length(mp.levels))
            push!(levels, mp.levels[lv]);   push!(configIndex, ci)
        end
    end

    order = sortperm([lv.energy for lv in levels])

    return( (multiplet = Multiplet(name, levels[order]), configIndex = configIndex[order]) )
end


"""
`Plasma.validateKineticChargeBalanceSupport(scheme::Plasma.KineticChargeBalanceScheme)`
    ... checks every leg of scheme against the chargeStages it touches and raises -- listing every violation found,
        each naming the leg, the stage(s) involved, and the nearest supported alternative -- if an unsupported
        (treatment tier, stage shape) combination is requested. Called BEFORE any Multiplet is generated (cf.
        Plasma.KineticChargeBalanceScheme), so an unsupported request fails fast and cheaply rather than part-way
        through real computation. Nothing is returned if every leg/stage combination is supported; an error is
        raised otherwise.
"""
function validateKineticChargeBalanceSupport(scheme::Plasma.KineticChargeBalanceScheme)
    problems = String[]
    nStages  = length(scheme.chargeStages)

    # --- within-stage legs: PhotoEmission, ImpactExcitation. Empirical/SemiEmpirical need a CLASSIC stage. --------
    for  (k, stage)  in  enumerate(scheme.chargeStages)
        Plasma.isClassicChargeStage(stage)   &&   continue
        for  (legName, treatment)  in  (("peTreatment", scheme.peTreatment), ("ieTreatment", scheme.ieTreatment))
            if  !(treatment isa Plasma.StructureBasedTreatment)
                push!(problems, "chargeStages[$k] is a GENERAL (non-classic) stage, but $legName = $treatment resolves " *
                                 "rates at configuration level; either make chargeStages[$k] classic (2 configs, 1 " *
                                 "level each), or set $legName = Plasma.StructureBasedTreatment().")
            end
        end
    end

    # --- cross-stage legs: ImpactIonization, PhotoRecombination, ThreeBodyRecombination. ---------------------------
    for  k = 1:nStages-1
        bothClassic = Plasma.isClassicChargeStage(scheme.chargeStages[k])  &&  Plasma.isClassicChargeStage(scheme.chargeStages[k+1])
        bothClassic   &&   continue
        for  (legName, treatment)  in  (("eiiTreatment", scheme.eiiTreatment), ("rrTreatment", scheme.rrTreatment), ("tbrTreatment", scheme.tbrTreatment))
            if  !(treatment isa Plasma.StructureBasedTreatment)
                push!(problems, "chargeStages[$k]/[$(k+1)] are not both classic stages, but $legName = $treatment resolves " *
                                 "rates at configuration level; either keep both stages classic, or set $legName = " *
                                 "Plasma.StructureBasedTreatment().")
            end
        end
    end

    # --- SemiEmpiricalTreatment's own approx requirement, one check per leg. ---------------------------------------
    if  scheme.peTreatment isa Plasma.SemiEmpiricalTreatment  &&  !(scheme.peTreatment.approx isa Empirical.UsingJAC)
        push!(problems, "peTreatment = Plasma.SemiEmpiricalTreatment requires approx = Empirical.UsingJAC(), got " *
                         "$(typeof(scheme.peTreatment.approx)).")
    end
    if  scheme.ieTreatment isa Plasma.SemiEmpiricalTreatment  &&  !(scheme.ieTreatment.approx isa Empirical.VanRegemorter1962)
        push!(problems, "ieTreatment = Plasma.SemiEmpiricalTreatment requires approx = Empirical.VanRegemorter1962(), got " *
                         "$(typeof(scheme.ieTreatment.approx)).")
    end
    if  scheme.rrTreatment isa Plasma.SemiEmpiricalTreatment  &&  !(scheme.rrTreatment.approx isa Empirical.UsingJAC)
        push!(problems, "rrTreatment = Plasma.SemiEmpiricalTreatment requires approx = Empirical.UsingJAC(), got " *
                         "$(typeof(scheme.rrTreatment.approx)).")
    end
    if  scheme.eiiTreatment isa Plasma.SemiEmpiricalTreatment  &&  !(scheme.eiiTreatment.approx isa Empirical.Lotz1967)
        push!(problems, "eiiTreatment = Plasma.SemiEmpiricalTreatment requires approx = Empirical.Lotz1967(), got " *
                         "$(typeof(scheme.eiiTreatment.approx)) -- SemiEmpirical feeds a mean-field binding energy into " *
                         "the SAME Lotz formula, cf. Empirical.impactIonizationBindingEnergy.")
    end
    if  scheme.tbrTreatment isa Plasma.SemiEmpiricalTreatment  &&  !(scheme.tbrTreatment.approx isa Empirical.Lotz1967)
        push!(problems, "tbrTreatment = Plasma.SemiEmpiricalTreatment requires approx = Empirical.Lotz1967(), got " *
                         "$(typeof(scheme.tbrTreatment.approx)) -- SemiEmpirical feeds a mean-field binding energy into " *
                         "the SAME Saha/Lotz formula, cf. Empirical.impactIonizationBindingEnergy.")
    end

    if  length(problems) > 0
        error("Plasma.KineticChargeBalanceScheme: $(length(problems)) unsupported treatment/stage combination(s) requested:\n  - " *
              join(problems, "\n  - "))
    end

    return( nothing )
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
    # PhotoRecombination.Settings() leaves maxKappa = 0, which admits NO partial wave and would raise outright
    # (cf. PhotoRecombination.checkPartialWavesAreAdmitted); 4 matches what the PhotoRecombination example branches use.
    prSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(); multipoles=[Basics.E1], gauges=[UseCoulomb, UseBabushkin],
                     electronEnergies=gridGL.t, maxKappa=4, printBefore=false)
    lines = PhotoRecombination.computeLines(finalMp, initialMp, nm, grid, prSettings; output=true)

    # Energies are matched by NEAREST VALUE, not exact equality: PhotoRecombination.computeLines does not return
    # settings.electronEnergies unchanged (it passes through its own internal unit handling first), so an exact
    # Dict lookup on the returned line.electronEnergy misses -- found by actually running this function for the
    # first time (it was committed but never exercised until a general multi-level stage called it).
    factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    alphaC = 0.;   alphaB = 0.
    for  line  in  lines
        n      = argmin(abs.(gridGL.t .- line.electronEnergy))
        wt     = gridGL.wt[n]
        fe     = Distribution.electronEnergyDistribution(Distribution.ElectronMaxwell(Te_au), line.electronEnergy)
        alphaC = alphaC + line.crossSection.Coulomb   * sqrt(2*line.electronEnergy) * fe * wt
        alphaB = alphaB + line.crossSection.Babushkin * sqrt(2*line.electronEnergy) * fe * wt
    end

    return( (Coulomb = factor*alphaC, Babushkin = factor*alphaB) )
end


"""
`Plasma.displayKineticChargeBalanceTreatments(stream::IO, scheme::Plasma.KineticChargeBalanceScheme)`
    ... prints, once per leg, which Plasma.AbstractRateTreatment was used -- the user-facing provenance table that
        lets a reader recognize, at a glance, at which (empirical, semi-empirical or structure-based) tier each
        process was treated, without having to read the scheme's fields back out of the code that built it.
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
    println(stream, "    chargeStages: " * join( ["$(length(s.configs)) config(s)/" * string(sum(s.nLevels)) * " level(s)"
                                                    for s in scheme.chargeStages], ", " ))

    return( nothing )
end


"""
`Plasma.displayKineticChargeBalanceResults(stream::IO, scheme::Plasma.KineticChargeBalanceScheme, Ne::Float64,
                                           popVectors::Dict, levelCounts::Array{Int64,1})`
    ... prints a table of the relative populations (Coulomb/Babushkin) of every (stage, local level) slot, one row
        per (slot, temperature) pair, following JAC's usual table style.
"""
function displayKineticChargeBalanceResults(stream::IO, scheme::Plasma.KineticChargeBalanceScheme, Ne::Float64, popVectors::Dict,
                                            levelCounts::Array{Int64,1})
    nx = 100
    println(stream, "\n  Kinetic charge-balance relative populations (Ne = $Ne cm^-3):")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    " * TableStrings.flushleft(20, "Stage") * TableStrings.flushleft(10, "Level") *
                  TableStrings.flushright(14, "Temperature") * TableStrings.flushright(22, "Population") * TableStrings.flushright(22, "Population")
    sb = "    " * TableStrings.flushleft(20, "") * TableStrings.flushleft(10, "") *
                  TableStrings.flushright(14, "[K]") * TableStrings.flushright(22, "(Coulomb)") * TableStrings.flushright(22, "(Babushkin)")
    println(stream, sa);   println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    offset = cumsum([0; levelCounts[1:end-1]])
    for  (k, stage)  in  enumerate(scheme.chargeStages)
        stageLabel = string(stage.configs[1])
        for  loc = 1:levelCounts[k]
            idx = offset[k] + loc
            for  temp  in  sort(collect(keys(popVectors)))
                pv = popVectors[temp]
                sc = "    " * TableStrings.flushleft(20, stageLabel) * TableStrings.flushleft(10, "#$loc") *
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

    Plasma.validateKineticChargeBalanceSupport(scheme)

    @warn("Plasma.KineticChargeBalanceScheme: DielectronicRecombination is deliberately not included (see the scheme's " *
          "own docstring), and the cross-stage legs (ImpactIonization, PhotoRecombination, ThreeBodyRecombination) " *
          "connect only the level pairs they are configured for -- on a classic stage pair that is the ground levels " *
          "only; a non-ground level is populated/depopulated by PhotoEmission/ImpactExcitation alone, never directly " *
          "by ionization or recombination, unless the corresponding leg is Plasma.StructureBasedTreatment.")

    treatments = (scheme.peTreatment, scheme.ieTreatment, scheme.eiiTreatment, scheme.rrTreatment, scheme.tbrTreatment)
    anyNonEmpirical = any( t -> !(t isa Plasma.EmpiricalTreatment), treatments )

    multipletData = Vector{Any}(nothing, nStages)
    levelCounts   = fill(2, nStages)
    if  anyNonEmpirical
        println(">> Plasma.KineticChargeBalanceScheme: at least one leg is not Plasma.EmpiricalTreatment; generating " *
                "a reference Multiplet for every one of the $nStages charge stages.")
        for  k = 1:nStages
            multipletData[k] = Plasma.generateChargeStageMultiplet(scheme.chargeStages[k], nm, grid, asf)
            levelCounts[k]   = length(multipletData[k].multiplet.levels)
        end
    else
        println(">> Plasma.KineticChargeBalanceScheme: every leg is Plasma.EmpiricalTreatment; no SCF/CI is performed.")
    end

    Plasma.displayKineticChargeBalanceTreatments(stdout, scheme)

    offset        = cumsum([0; levelCounts[1:end-1]])
    globalIndex(k, loc) = offset[k] + loc
    N             = sum(levelCounts)
    temperatures  = scheme.electronTemperatures
    rateFactor    = Defaults.convertUnits("rate: from atomic to 1/s", 1.0)
    coeffFactor   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)

    # -------------------------------------------------------------------------------------------------------
    # (1) Within-stage PhotoEmission: temperature-independent, computed once. peTransitions[k] lists every
    #     (iLoc, jLoc, rateC, rateB) radiative-decay pair found for stage k, in its OWN local level numbering.
    # -------------------------------------------------------------------------------------------------------
    TransitionT = NamedTuple{(:iLoc,:jLoc,:rateC,:rateB), Tuple{Int64,Int64,Float64,Float64}}
    peTransitions = Vector{Vector{TransitionT}}(undef, nStages)
    einsteinAs    = Vector{Any}(undef, nStages)                                  # classic-stage IE detailed balance only
    for  k = 1:nStages
        stage       = scheme.chargeStages[k]
        transitions = TransitionT[]
        if  scheme.peTreatment isa Plasma.StructureBasedTreatment
            mp         = multipletData[k].multiplet
            peSettings = PhotoEmission.Settings(PhotoEmission.Settings(); gauges=[UseCoulomb, UseBabushkin], printBefore=false)
            peLines    = PhotoEmission.computeLines(mp, mp, grid, peSettings; output=true)
            for  line  in  peLines
                push!(transitions, (iLoc=line.initialLevel.index, jLoc=line.finalLevel.index,
                                     rateC=rateFactor*line.photonRate.Coulomb, rateB=rateFactor*line.photonRate.Babushkin))
            end
            if  Plasma.isClassicChargeStage(stage)
                einsteinAs[k] = Empirical.photoemissionEinsteinA(stage.configs[2], stage.configs[1], Empirical.ScaledHydrogenic(); printout=false)
            end
        elseif  scheme.peTreatment isa Plasma.SemiEmpiricalTreatment
            ea = Empirical.photoemissionEinsteinA(stage.configs[2], stage.configs[1], Empirical.UsingJAC(); printout=false)
            push!(transitions, (iLoc=2, jLoc=1, rateC=rateFactor*ea.rate, rateB=rateFactor*ea.rate))
            einsteinAs[k] = ea
        else                                                                     # Plasma.EmpiricalTreatment
            ea = Empirical.photoemissionEinsteinA(stage.configs[2], stage.configs[1], scheme.peTreatment.approx; printout=false)
            push!(transitions, (iLoc=2, jLoc=1, rateC=rateFactor*ea.rate, rateB=rateFactor*ea.rate))
            einsteinAs[k] = ea
        end
        peTransitions[k] = transitions
    end

    popVectors = Dict{Float64, NamedTuple{(:Coulomb,:Babushkin), Tuple{Array{Float64,1},Array{Float64,1}}}}()

    for  TeK  in  temperatures
        Te_au   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", TeK)
        cMatrix = zeros(N, N);   bMatrix = zeros(N, N)

        # -----------------------------------------------------------------------------------------------------
        # (2) Within-stage legs: PhotoEmission (precomputed above) and ImpactExcitation (both directions).
        # -----------------------------------------------------------------------------------------------------
        for  k = 1:nStages
            stage = scheme.chargeStages[k]

            for  t  in  peTransitions[k]
                gi = globalIndex(k, t.iLoc);   gj = globalIndex(k, t.jLoc)       # gi = source (initial/excited), gj = destination (final/ground)
                cMatrix[gi,gj] += t.rateC;   cMatrix[gi,gi] -= t.rateC
                bMatrix[gi,gj] += t.rateB;   bMatrix[gi,gi] -= t.rateB
            end

            if  scheme.ieTreatment isa Plasma.StructureBasedTreatment
                mp    = multipletData[k].multiplet
                L     = levelCounts[k]
                pairs = [ (i,j)  for i = 1:L  for j = i+1:L ]
                if  length(pairs) > 0
                    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings(); calcRateCoefficient=true, printBefore=false,
                                     lineSelection=LineSelection(true; indexPairs=pairs), temperatures=[TeK])
                    ieLines, _ = ImpactExcitation.computeLines(mp, mp, nm, grid, ieSettings; output=true)
                    ieRates    = ImpactExcitation.computeRateCoefficients(ieLines, ieSettings)
                    for  rc  in  ieRates
                        gLower    = globalIndex(k, rc.initialLevel.index);   gUpperIdx = globalIndex(k, rc.finalLevel.index)
                        excRate   = rc.alphas[1]
                        gUpper    = Float64(Basics.twice(rc.finalLevel.J) + 1)
                        deExcRate = 8.6291269e-6 * rc.effOmegas[1] / (gUpper * sqrt(TeK))
                        if  isfinite(excRate)                                    # source = gLower (ground), destination = gUpperIdx (excited)
                            cMatrix[gLower,gUpperIdx] += excRate*Ne;   cMatrix[gLower,gLower] -= excRate*Ne
                            bMatrix[gLower,gUpperIdx] += excRate*Ne;   bMatrix[gLower,gLower] -= excRate*Ne
                        end
                        if  isfinite(deExcRate)                                  # source = gUpperIdx (excited), destination = gLower (ground)
                            cMatrix[gUpperIdx,gLower] += deExcRate*Ne;   cMatrix[gUpperIdx,gUpperIdx] -= deExcRate*Ne
                            bMatrix[gUpperIdx,gLower] += deExcRate*Ne;   bMatrix[gUpperIdx,gUpperIdx] -= deExcRate*Ne
                        end
                    end
                end
            else                                                                 # Empirical or SemiEmpirical, classic stage
                gIdx = globalIndex(k, 1);   eIdx = globalIndex(k, 2)
                if  scheme.ieTreatment isa Plasma.SemiEmpiricalTreatment
                    bridge  = Empirical.photoemissionEinsteinA(stage.configs[2], stage.configs[1], Empirical.UsingJAC(); printout=false)
                    aSource = Empirical.GivenEinsteinA(bridge.multipole, bridge.energy, bridge.rate)
                    eDist   = Distribution.ElectronMaxwell(Te_au)
                    excRate = coeffFactor * Empirical.impactExcitationPlasmaAlpha(eDist, stage.configs[1], stage.configs[2];
                                                approx=scheme.ieTreatment.approx, aSource=aSource, printout=false)
                    deltaE  = bridge.energy
                else
                    eDist   = Distribution.ElectronMaxwell(Te_au)
                    excRate = coeffFactor * Empirical.impactExcitationPlasmaAlpha(eDist, stage.configs[1], stage.configs[2];
                                                approx=scheme.ieTreatment.approx, aSource=Empirical.ScaledHydrogenic(), printout=false)
                    deltaE  = einsteinAs[k].energy
                end
                # Detailed balance for the (empirical/semi-empirical) de-excitation rate coefficient:
                #   g_ground C_exc(Te) = g_excited C_deexc(Te) exp(-dE/Te).
                gi        = Basics.extractFromConfiguration(Basics.Multiplicity(), stage.configs[1])
                ge        = Basics.extractFromConfiguration(Basics.Multiplicity(), stage.configs[2])
                deExcRate = (gi/ge) * excRate * exp( deltaE / Te_au )

                if  isfinite(excRate)                                           # source = gIdx (ground), destination = eIdx (excited)
                    cMatrix[gIdx,eIdx] += excRate*Ne;   cMatrix[gIdx,gIdx] -= excRate*Ne
                    bMatrix[gIdx,eIdx] += excRate*Ne;   bMatrix[gIdx,gIdx] -= excRate*Ne
                end
                if  isfinite(deExcRate)                                         # source = eIdx (excited), destination = gIdx (ground)
                    cMatrix[eIdx,gIdx] += deExcRate*Ne;   cMatrix[eIdx,eIdx] -= deExcRate*Ne
                    bMatrix[eIdx,gIdx] += deExcRate*Ne;   bMatrix[eIdx,eIdx] -= deExcRate*Ne
                end
            end
        end

        # ---------------------------------------------------------------------------------------------------
        # (3) Cross-stage legs, stage(k) <-> stage(k+1): ImpactIonization (up), PhotoRecombination +
        #     ThreeBodyRecombination (down, summed into the same matrix entry). On a classic stage pair with
        #     an Empirical leg this connects ground levels only, exactly as before; a StructureBasedTreatment
        #     leg connects EVERY level pair between the two stages' multiplets.
        # ---------------------------------------------------------------------------------------------------
        for  k = 1:nStages-1
            stageLo = scheme.chargeStages[k];   stageHi = scheme.chargeStages[k+1]

            # A CLASSIC stage pair connects ground levels only, regardless of treatment -- this is the scheme's own
            # documented scope decision (cf. Plasma.KineticChargeBalanceScheme), not a side effect of which
            # treatment happens to be selected; branch (b) of example-Ji.jl relies on exactly this to keep its
            # already-published ~135x finding reproducible. Only a GENERAL (non-classic) stage, which the guard
            # already restricts to Plasma.StructureBasedTreatment on every leg that touches it, uses all level pairs.
            bothClassicPair = Plasma.isClassicChargeStage(stageLo)  &&  Plasma.isClassicChargeStage(stageHi)
            eiiPairs = bothClassicPair ? [ (1,1) ] : [ (i,j) for i=1:levelCounts[k] for j=1:levelCounts[k+1] ]
            for  (iLo, iHi)  in  eiiPairs
                gLo = globalIndex(k, iLo);   gHi = globalIndex(k+1, iHi)
                if  scheme.eiiTreatment isa Plasma.StructureBasedTreatment
                    levelLo     = multipletData[k].multiplet.levels[iLo];   levelHi = multipletData[k+1].multiplet.levels[iHi]
                    sh          = ThreeBodyRecombination.capturedSubshell(levelHi, levelLo)
                    eiiSettings = ThreeBodyRecombination.Settings(ThreeBodyRecombination.Settings(); eiiModel=scheme.eiiTreatment.eiiModel,
                                                                  electronTemperature=TeK, printBefore=false)
                    eiiRate     = ThreeBodyRecombination.eiiRateCoefficient(sh, levelLo, nm, grid, eiiSettings)
                else
                    eDist    = Distribution.ElectronMaxwell(Te_au)
                    pSource  = scheme.eiiTreatment isa Plasma.SemiEmpiricalTreatment ? Empirical.UsingJAC() : Empirical.ScaledHydrogenic()
                    eiiRate = coeffFactor * Empirical.impactIonizationPlasmaAlpha(eDist, stageLo.configs[1], stageHi.configs[1];
                                                approx=scheme.eiiTreatment.approx, pSource=pSource, printout=false)
                end
                eiiTerm = eiiRate * Ne
                if  isfinite(eiiTerm)
                    cMatrix[gLo,gHi] += eiiTerm;   cMatrix[gLo,gLo] -= eiiTerm
                    bMatrix[gLo,gHi] += eiiTerm;   bMatrix[gLo,gLo] -= eiiTerm
                end
            end

            downPairs = bothClassicPair ? [ (1,1) ] : [ (i,j) for i=1:levelCounts[k] for j=1:levelCounts[k+1] ]
            for  (iLo, iHi)  in  downPairs
                gLo = globalIndex(k, iLo);   gHi = globalIndex(k+1, iHi)

                if  scheme.rrTreatment isa Plasma.StructureBasedTreatment
                    levelLo = multipletData[k].multiplet.levels[iLo];   levelHi = multipletData[k+1].multiplet.levels[iHi]
                    rr      = Plasma.structureBasedRrRateCoefficient(levelLo, levelHi, nm, grid, TeK)
                    rrRateC = rr.Coulomb;   rrRateB = rr.Babushkin
                else
                    eDist    = Distribution.ElectronMaxwell(Te_au)
                    rrRateAU = coeffFactor * Empirical.photorecombinationPlasmaAlpha(eDist, stageHi.configs[1], stageLo.configs[1];
                                                approx=scheme.rrTreatment.approx, printout=false)
                    rrRateC = rrRateB = rrRateAU
                end

                if  scheme.tbrTreatment isa Plasma.StructureBasedTreatment
                    levelLo     = multipletData[k].multiplet.levels[iLo];   levelHi = multipletData[k+1].multiplet.levels[iHi]
                    tbrSettings = ThreeBodyRecombination.Settings(ThreeBodyRecombination.Settings(); eiiModel=scheme.tbrTreatment.eiiModel,
                                                                  electronTemperature=TeK, printBefore=false)
                    finalMp     = Multiplet("TBR final",   [levelLo]);   initialMp = Multiplet("TBR initial", [levelHi])
                    tbrLines    = ThreeBodyRecombination.computeLines(finalMp, initialMp, nm, grid, tbrSettings; output=true)
                    tbrRate     = length(tbrLines) > 0 ? tbrLines[1].rateCoefficient : 0.
                else
                    eDist    = Distribution.ElectronMaxwell(Te_au)
                    pSource  = scheme.tbrTreatment isa Plasma.SemiEmpiricalTreatment ? Empirical.UsingJAC() : Empirical.ScaledHydrogenic()
                    tbrRate = coeffFactor * Empirical.threeBodyRecombinationPlasmaAlpha(eDist, stageHi.configs[1], stageLo.configs[1];
                                                approx=scheme.tbrTreatment.approx, pSource=pSource, printout=false)
                end

                downTermC = (rrRateC + tbrRate) * Ne;   downTermB = (rrRateB + tbrRate) * Ne
                if  isfinite(downTermC)   cMatrix[gHi,gLo] += downTermC;   cMatrix[gHi,gHi] -= downTermC   end
                if  isfinite(downTermB)   bMatrix[gHi,gLo] += downTermB;   bMatrix[gHi,gHi] -= downTermB   end
            end
        end

        popVectors[TeK] = (Coulomb = Plasma.solveCollisionalRadiativeBalance(cMatrix), Babushkin = Plasma.solveCollisionalRadiativeBalance(bMatrix))
    end

    Plasma.displayKineticChargeBalanceResults(stdout, scheme, Ne, popVectors, levelCounts)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayKineticChargeBalanceResults(iostream, scheme, Ne, popVectors, levelCounts)   end

    if output    results = Base.merge( results, Dict("populations:" => popVectors) )   end

    println("KineticChargeBalanceScheme computation complete ...")
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
