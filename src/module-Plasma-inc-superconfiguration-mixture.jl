# Functions and methods for scheme::Plasma.SuperConfigurationMixtureOpacityScheme computations.
#
# EXTENDED 18-Oct-2026 to cover a species' own STAGE axis, not just the species axis: every (species,stage)
# pair gets its own Group (via SuperConfigurations.buildGroup, for the coupled Saha-Boltzmann balance) AND
# is run through Plasma.perform(::SuperConfigurationOpacityScheme, ...) UNCHANGED (for its own Stage-2/3
# kappa). The NEW physics is the coupled multi-species charge-neutrality solve; the mixing and mean-opacity
# machinery below is otherwise the SAME as the original (17-Oct-2026) version. See
# Plasma.SuperConfigurationMixtureOpacityScheme's own docstring for the coupled-balance equation, the
# two-stage mixing rule, and the weight functions.


"""
`Plasma.planckWeight(x::Float64)`
    ... the normalized Planck mean-opacity weight function at x = omega/kT: `(15/pi^4) * x^3/(e^x-1)`,
        INTEGRAL over x in (0,infinity) of this function is exactly 1 (verified numerically to 1e-10 over
        x in [0,60] before this formula was written into the code, not assumed from memory -- see this
        module's own docstring for why that check matters here). A value::Float64 is returned.
"""
function  planckWeight(x::Float64)
    x <= 0.0   &&   return( 0.0 )
    ex = exp(-x)

    return( (15/pi^4) * x^3 * ex / (1-ex) )
end


"""
`Plasma.rosselandWeight(x::Float64)`
    ... the normalized Rosseland mean-opacity weight function at x = omega/kT: `(15/(4*pi^4)) * x^4*e^x/(e^x-1)^2`,
        INTEGRAL over x in (0,infinity) of this function is exactly 1 (verified the same way as
        `Plasma.planckWeight`). A value::Float64 is returned.
"""
function  rosselandWeight(x::Float64)
    x <= 0.0   &&   return( 0.0 )
    ex = exp(x)

    return( (15/(4*pi^4)) * x^4 * ex / (ex-1)^2 )
end


"""
`Plasma.planckMeanOpacity(xGrid::Array{Float64,1}, kappa::Array{Float64,1})`
    ... the Planck mean of `kappa` over `xGrid` (x = omega/kT), `INTEGRAL kappa(x)*planckWeight(x) dx`,
        evaluated by the trapezoidal rule over the GIVEN grid -- the caller's own `opacityOmegas` choice,
        not a refined quadrature, so the grid must already resolve both `planckWeight`'s own peak (x ~ 3-4)
        and every line in `kappa`. A value::Float64 is returned.
"""
function  planckMeanOpacity(xGrid::Array{Float64,1}, kappa::Array{Float64,1})
    n = length(xGrid)
    n < 2   &&   error("Plasma.planckMeanOpacity(): xGrid needs at least two points.")

    integral = 0.0
    for  k = 1:n-1
        dx       = xGrid[k+1] - xGrid[k]
        fk       = kappa[k]   * planckWeight(xGrid[k])
        fk1      = kappa[k+1] * planckWeight(xGrid[k+1])
        integral = integral + 0.5*(fk+fk1)*dx
    end

    return( integral )
end


"""
`Plasma.rosselandMeanOpacity(xGrid::Array{Float64,1}, kappa::Array{Float64,1})`
    ... the Rosseland mean of `kappa` over `xGrid` (x = omega/kT), the HARMONIC mean
        `1 / INTEGRAL( rosselandWeight(x)/kappa(x) dx )`, evaluated by the trapezoidal rule over the GIVEN
        grid.

        **NOT CALLED BY `Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...)`, AND NAMED HERE
        RATHER THAN SILENTLY SHIPPED, because it was MEASURED, not assumed, to have no stable value for a
        LINE-ONLY spectrum.** A Rosseland mean is a harmonic average, dominated by whichever frequencies
        are LEAST opaque -- and a `kappa(omega)` built only from Gaussian E1 lines (no bound-free, no
        free-free continuum) never reaches a true zero between lines, only an exponentially small one, so
        the function needs a FLOOR below which `kappa` is treated as the model's own zero rather than
        divided into. Tried on a real two-configuration case, 17-Oct-2026: the returned value moved by
        **more than eleven orders of magnitude** as the floor was swept from `1e-3` to `1e-15` of the peak
        (0.049 down to 7.7e-13) -- proof that the number this function returns is set by an ARBITRARY
        internal threshold, not by the physics of the lines. A quantity that moves eleven orders of
        magnitude with an implementation detail has no content worth reporting. Kept here, floored at
        `max(kappa)*1e-12` purely so the function itself does not divide by zero or overflow, as a
        BUILDING BLOCK for the day a continuum opacity exists to fill the gaps a line-only model cannot --
        at that point the harmonic mean becomes well-posed again, because the continuum (not an arbitrary
        floor) sets what "least opaque" means. Until then, `Plasma.planckMeanOpacity` -- an ARITHMETIC,
        kappa-weighted average with no reciprocal and therefore no gap-sensitivity -- is the only mean
        opacity this scheme reports. A value::Float64 is returned.
"""
function  rosselandMeanOpacity(xGrid::Array{Float64,1}, kappa::Array{Float64,1})
    n = length(xGrid)
    n < 2   &&   error("Plasma.rosselandMeanOpacity(): xGrid needs at least two points.")
    kappaMax = maximum(kappa)
    kappaMax <= 0.0   &&
        error("Plasma.rosselandMeanOpacity(): kappa is zero everywhere -- no opacity to average.")
    floorVal = kappaMax * 1.0e-12

    integral = 0.0
    for  k = 1:n-1
        dx  = xGrid[k+1] - xGrid[k]
        fk  = kappa[k]   < floorVal   ?   0.0   :   rosselandWeight(xGrid[k])   / kappa[k]
        fk1 = kappa[k+1] < floorVal   ?   0.0   :   rosselandWeight(xGrid[k+1]) / kappa[k+1]
        integral = integral + 0.5*(fk+fk1)*dx
    end

    integral <= 0.0   &&
        error("Plasma.rosselandMeanOpacity(): the Rosseland integral is zero -- kappa never reaches the " *
              "floor anywhere the weight function has support; widen opacityOmegas or check the species' arrays.")

    return( 1.0 / integral )
end


"""
`Plasma.displayMixtureSpectrum(stream::IO, omegas::Array{Float64,1}, kappaMix::Array{Float64,1}, kappaP::Float64)`
    ... prints the mixture's own summed per-absorber cross section kappaMix(omega) at every `omegas` entry,
        followed by its Planck mean opacity `kappaP`. The Rosseland mean is DELIBERATELY not printed here --
        see `Plasma.rosselandMeanOpacity`'s own docstring for why it is not yet a trustworthy number for a
        line-only spectrum. Nothing is returned.
"""
function  displayMixtureSpectrum(stream::IO, omegas::Array{Float64,1}, kappaMix::Array{Float64,1}, kappaP::Float64)
    println(stream, "\n  Superconfiguration MIXTURE opacity spectrum (per-absorber cross section):\n")
    println(stream, "  ", TableStrings.hLine(46))
    sa = "    " * TableStrings.center(18, "omega [a.u.]"; na=2)
    sa = sa     * TableStrings.center(20, "kappa_mix(omega) [a0^2]"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(46))
    for  (om, kp) in zip(omegas, kappaMix)
        sa = "    " * @sprintf("%18.6f", om)
        sa = sa     * @sprintf("%20.6e", kp)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(46))
    println(stream, "\n  Planck mean opacity:     $(@sprintf("%.6e", kappaP)) [a0^2]")
    println(stream, "  (Rosseland mean not reported -- not yet well-defined for a line-only spectrum; " *
                     "see Plasma.rosselandMeanOpacity's own docstring)")

    return( nothing )
end


"""
`Plasma.displayMixtureBalance(stream::IO, charges::Array{Array{Int64,1},1},
                                 stagePopulations::Array{Array{Float64,1},1}, numberFractions::Array{Float64,1},
                                 electronDensity::Float64)`
    ... prints, for every species, its own stages' charges and population fractions (summing to 1 WITHIN
        that species), then the ONE electron density shared by the whole mixture. Nothing is returned.
"""
function  displayMixtureBalance(stream::IO, charges::Array{Array{Int64,1},1},
                                   stagePopulations::Array{Array{Float64,1},1}, numberFractions::Array{Float64,1},
                                   electronDensity::Float64)
    println(stream, "\n  Mixture ionization balance, species by species:\n")
    for  s = 1:length(charges)
        println(stream, "    species $s, number fraction $(numberFractions[s]):")
        println(stream, "    ", TableStrings.hLine(30))
        sa = "      " * TableStrings.center(10, "Charge q"; na=2)
        sa = sa       * TableStrings.center(14, "Population"; na=2)
        println(stream, sa)
        println(stream, "    ", TableStrings.hLine(30))
        for  (q, pop) in zip(charges[s], stagePopulations[s])
            sa = "      " * TableStrings.center(10, string(q); na=2)
            sa = sa       * @sprintf("%14.6f", pop)
            println(stream, sa)
        end
        println(stream, "    ", TableStrings.hLine(30))
    end
    println(stream, "\n    shared converged electron number density = ", @sprintf("%.6e", electronDensity), " [a.u.]")

    return( nothing )
end


"""
`Plasma.perform(scheme::Plasma.SuperConfigurationMixtureOpacityScheme, computation::Plasma.Computation; output::Bool=true)`
    ... STAGE 0: for every species and every one of its OWN stages in `scheme.refConfigsPerSpeciesStages`,
        builds one `SuperConfigurations.Group` per configuration, each on that (species,stage)'s own
        automatically-sized grid. STAGE 1: solves the COUPLED multi-species grand-canonical Saha-Boltzmann
        balance (see this scheme's own docstring for the equation) for the ONE electron density shared by
        the whole mixture at `computation.settings.temperature` [Hartree] and `computation.settings.density`
        [a.u., the TOTAL ION number density of the WHOLE mixture; needs `computation.settings.useNumberDensity`
        true] -- giving each species its OWN stage populations, summing to 1 within that species. **SKIPPED
        ENTIRELY, with every stage population set to 1 directly, when NO species has more than one stage**
        (the balance is then trivial by construction) -- which is exactly why the backward-compatible,
        single-stage-per-species constructor needs no density settings at all, same as the ORIGINAL
        (17-Oct-2026) scheme never did. A MIXTURE WITH ANY MULTI-STAGE SPECIES, however, DOES need a valid
        `computation.settings.density`/`useNumberDensity`, and raises if given none. STAGE 2/3: runs
        `Plasma.perform(::SuperConfigurationOpacityScheme, ...)` UNCHANGED for every
        (species,stage) pair, population-weights each stage's own `kappa` WITHIN its species by that stage's
        balance population, then weights ACROSS species by `numberFractions` into `kappaMix`, and reduces
        that to its Planck mean opacity. The Rosseland mean is DELIBERATELY NOT computed -- see
        `Plasma.rosselandMeanOpacity`'s own docstring. `computation.nuclearModel`/`computation.grid`/
        `computation.refConfigs` are UNUSED; only `computation.asfSettings` and `computation.settings` are
        shared across every (species,stage). A MEASURED, ACCEPTED COST of composing two already-separate
        schemes this way: every (species,stage)'s SCF and grid are built TWICE -- once here for its `Group`
        (Stage 0) and again inside `SuperConfigurationOpacityScheme`'s own `perform` call (Stage 2/3) --
        rather than changing THAT scheme's interface to accept externally-built orbitals, which would be a
        separate module's change. For output=true, a dictionary is returned under the keys
        "superconfiguration mixture species opacities" (`Array{Array{Float64,1},1}`, one spectrum per
        species), "superconfiguration mixture opacity spectrum" (`kappaMix`, `Array{Float64,1}`),
        "superconfiguration mixture planck opacity" (`Float64`), "superconfiguration mixture charges"
        (`Array{Array{Int64,1},1}`, one array per species), "superconfiguration mixture stage populations"
        (`Array{Array{Float64,1},1}`, matching, summing to 1 WITHIN each species), and "superconfiguration
        mixture electron density" (`Float64`, the ONE value shared by the whole mixture).
"""
function  perform(scheme::Plasma.SuperConfigurationMixtureOpacityScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nSpecies = length(scheme.refConfigsPerSpeciesStages)
    nSpecies == 0   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): refConfigsPerSpeciesStages is empty.")
    ( length(scheme.nuclearModels) == nSpecies  &&  length(scheme.numberFractions) == nSpecies )   ||
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): nuclearModels, " *
              "refConfigsPerSpeciesStages and numberFractions must all have the same length.")
    isempty(scheme.opacityOmegas)   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): opacityOmegas is empty -- " *
              "a mixture needs a shared frequency grid to combine species on.")
    fracSum = sum(scheme.numberFractions)
    abs(fracSum-1.0) > 1.0e-6   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): numberFractions sum to " *
              "$fracSum, not 1.")
    kT = computation.settings.temperature
    kT <= 0.0   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): temperature must be positive.")

    # A species with exactly ONE stage has a TRIVIAL balance (population = 1, whatever mu_e is), so a
    # number density is only genuinely NEEDED when some species has more than one -- exactly the condition
    # under which the ORIGINAL (17-Oct-2026) single-stage-per-species call sites, which never supplied one,
    # must keep working unchanged.
    nStages      = [ length(stages)   for stages in scheme.refConfigsPerSpeciesStages ]
    anyMultiStage = any( n -> n > 1, nStages )
    haveDensity   = computation.settings.useNumberDensity  &&  computation.settings.density > 0.0
    anyMultiStage  &&  !haveDensity   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): some species has more than " *
              "one stage, so settings.useNumberDensity must be true with settings.density the TOTAL ION " *
              "number density of the WHOLE mixture, in a.u.")
    nIonTotal = haveDensity ? computation.settings.density : NaN

    SC = superConfigurationsModule()

    # STAGE 0: one Group per configuration, for every stage of every species, each on its OWN grid/orbitals.
    groupsPerSpeciesStage = Vector{Vector{Any}}(undef, nSpecies)
    chargesPerSpecies     = Vector{Vector{Int64}}(undef, nSpecies)
    ErefPerSpecies        = zeros(Float64, nSpecies)
    for  s = 1:nSpecies
        nm     = scheme.nuclearModels[s]
        stages = scheme.refConfigsPerSpeciesStages[s]
        isempty(stages)   &&
            error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): species $s has no stages.")

        groupsPerSpeciesStage[s] = Vector{Any}(undef, nStages[s])
        chargesPerSpecies[s]     = zeros(Int64, nStages[s])
        for  i = 1:nStages[s]
            refConfigs = stages[i]
            isempty(refConfigs)   &&
                error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): species $s, stage $i " *
                      "has no reference configurations.")
            any( c -> c.NoElectrons != refConfigs[1].NoElectrons, refConfigs )   &&
                error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): species $s, stage $i " *
                      "mixes configurations of different electron count; one stage must be ONE electron count.")

            grid     = Basics.recommendedGrid(refConfigs, nm; printout=false)
            wb       = SelfConsistent.performSCF(refConfigs, nm, grid, computation.asfSettings; printout=false)
            orbitals = wb.levels[1].basis.orbitals
            groupsPerSpeciesStage[s][i] = [ SC.buildGroup(conf, orbitals, nm, grid, computation.asfSettings)
                                             for conf in refConfigs ]
            chargesPerSpecies[s][i]     = round(Int64, nm.Z) - refConfigs[1].NoElectrons
        end
        ErefPerSpecies[s] = minimum( minimum(g.meanEnergy for g in groups)  for groups in groupsPerSpeciesStage[s] )
    end

    # STAGE 1: the COUPLED fixed-point for mu_e/n_e -- SuperConfigurationIonizationScheme's own equations,
    # with the charge-neutrality sum lifted to range over every (species,stage) pair, each species weighted
    # by its OWN numberFractions entry. See this scheme's own docstring for the equation. SKIPPED entirely
    # when no density was given AND none is needed (every species single-stage, guaranteed by the earlier
    # check): the balance is then trivial, population = 1, by construction rather than by iterating to it.
    local ne, stagePopulations, converged
    if  haveDensity
        ne               = nIonTotal
        stagePopulations = [ zeros(Float64, nStages[s])   for s = 1:nSpecies ]
        converged        = false
        for  iteration = 1:100
            eLambda3 = (2*pi/kT)^1.5
            muE      = -kT * log( 2.0 / (ne*eLambda3) )

            for  s = 1:nSpecies
                weights = [ [ g.degeneracy * exp( -(muE*chargesPerSpecies[s][i] + g.meanEnergy - ErefPerSpecies[s])/kT )
                              for g in groupsPerSpeciesStage[s][i] ]   for i = 1:nStages[s] ]
                Zspecies = sum( sum(w)  for w in weights )
                for  i = 1:nStages[s]
                    stagePopulations[s][i] = sum(weights[i]) / Zspecies
                end
            end

            neNew = nIonTotal * sum( scheme.numberFractions[s] *
                                       sum( chargesPerSpecies[s][i]*stagePopulations[s][i]   for i = 1:nStages[s] )
                                       for s = 1:nSpecies )
            if  abs(neNew-ne) < 1.0e-10*ne
                converged = true;   ne = neNew;   break
            end
            ne = neNew
        end
        converged   ||
            @warn("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): the coupled electron-density " *
                  "self-consistency did NOT converge in 100 iterations; ne = $ne a.u. is the last iterate, not a " *
                  "solution.")
    else
        ne               = NaN    # no density was given, and none is needed -- see this scheme's own docstring
        stagePopulations = [ [1.0]   for s = 1:nSpecies ]
        converged        = true
    end

    Plasma.displayMixtureBalance(stdout, chargesPerSpecies, stagePopulations, scheme.numberFractions, ne)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayMixtureBalance(iostream, chargesPerSpecies, stagePopulations, scheme.numberFractions, ne)   end

    # STAGE 2/3, per (species,stage): SuperConfigurationOpacityScheme's own machinery, UNCHANGED; weight
    # WITHIN a species by that stage's own balance population, then ACROSS species by numberFractions.
    kappaMix      = zeros(Float64, length(scheme.opacityOmegas))
    speciesKappas = Array{Float64,1}[]
    for  s = 1:nSpecies
        nm           = scheme.nuclearModels[s]
        stages       = scheme.refConfigsPerSpeciesStages[s]
        kappaSpecies = zeros(Float64, length(scheme.opacityOmegas))
        for  i = 1:nStages[s]
            refConfigs = stages[i]
            println("\n  -- species $s of $nSpecies, stage $i of $(nStages[s]) " *
                     "(balance population $(stagePopulations[s][i])) --")

            grid               = Basics.recommendedGrid(refConfigs, nm; printout=false)
            stageScheme        = SuperConfigurationOpacityScheme(refConfigs, scheme.opacityOmegas)
            stageComputation   = Plasma.Computation(computation; scheme=stageScheme, nuclearModel=nm, grid=grid,
                                                       refConfigs=Configuration[])
            stageResults       = Plasma.perform(stageScheme, stageComputation; output=true)
            kappaStage         = stageResults[ResultKeys.SuperConfigurationOpacitySpectrum]
            kappaSpecies .+= stagePopulations[s][i] .* kappaStage
        end
        push!(speciesKappas, kappaSpecies)
        kappaMix .+= scheme.numberFractions[s] .* kappaSpecies
    end

    xGrid  = scheme.opacityOmegas ./ kT
    kappaP = Plasma.planckMeanOpacity(xGrid, kappaMix)

    Plasma.displayMixtureSpectrum(stdout, scheme.opacityOmegas, kappaMix, kappaP)
    if  printSummary   Plasma.displayMixtureSpectrum(iostream, scheme.opacityOmegas, kappaMix, kappaP)   end

    println("SuperConfigurationMixtureOpacityScheme computation complete ($nSpecies species, electron " *
            "density " * (converged ? "converged" : "NOT CONVERGED") * ", Planck $(@sprintf("%.4e", kappaP)) " *
            "[a0^2]; Rosseland not reported, see this scheme's own docstring) ...")

    if  output
        results = Base.merge( results, Dict(
            "superconfiguration mixture species opacities" => speciesKappas,
            "superconfiguration mixture opacity spectrum"  => kappaMix,
            "superconfiguration mixture planck opacity"    => kappaP,
            "superconfiguration mixture charges"           => chargesPerSpecies,
            "superconfiguration mixture stage populations" => stagePopulations,
            "superconfiguration mixture electron density"  => ne) )
    end

    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
