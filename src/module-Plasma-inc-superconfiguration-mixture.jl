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
    ... prints the mixture's own summed BOUND-BOUND per-absorber cross section kappaMix(omega) at every
        `omegas` entry, followed by its Planck mean opacity `kappaP`. Used UNCHANGED whenever
        `scheme.boundFreeApprox === nothing` (the default); see `Plasma.displayMixtureSpectrumBreakdown` for
        the richer, three-column table used once bound-free is turned on. The Rosseland mean is DELIBERATELY
        not printed here -- see `Plasma.rosselandMeanOpacity`'s own docstring for why it is not yet a
        trustworthy number for a line-only spectrum. Nothing is returned.
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
`Plasma.boundFreeOpacityContribution(stages::Array{Array{Configuration,1},1}, popsPerStage::Array{Array{Float64,1},1},
                                        stagePopulations::Array{Float64,1}, omegas::Array{Float64,1},
                                        approx::Empirical.AbstractEmpiricalApproximation)`
    ... sums the bound-free (photoionization) contribution to kappa(omega), for ONE species, over every pair
        of groups in every pair of ADJACENT stages whose configurations differ by exactly one electron in
        exactly one shell -- checked by `Basics.extractFromConfigurations(Basics.OccupationDifference(),...)`,
        the SAME compatibility test `Empirical.photoionizationCrossSection` applies internally. A pair that is
        not a genuine single-electron ionization channel is simply SKIPPED, not an error, the same way a
        same-parity pair is skipped (not an error) for the bound-bound Stage 2 in
        `Plasma.perform(::SuperConfigurationOpacityScheme,...)`. Each valid channel's cross section -- from
        `Empirical.photoionizationCrossSection`, Kramers' (1923) formula scaled by the ionized shell's own
        binding energy and a bound-free Gaunt factor -- is weighted by the INITIAL group's own ABSOLUTE
        population within the species, `stagePopulations[i] * popsPerStage[i][g]`; the FINAL state's own
        population plays no role, since a photoionization event is driven by the absorber, not by how occupied
        the state it is ejected into happens to be. A `(kappaBF::Array{Float64,1},
        channels::Array{Tuple{Configuration,Configuration,Float64,Array{Float64,1}},1})` tuple is returned --
        `channels` holding `(initialConf, finalConf, population, crossSection(omegas))` for every channel
        actually used, for display.
"""
function  boundFreeOpacityContribution(stages::Array{Array{Configuration,1},1}, popsPerStage::Array{Array{Float64,1},1},
                                          stagePopulations::Array{Float64,1}, omegas::Array{Float64,1},
                                          approx::Empirical.AbstractEmpiricalApproximation)
    kappaBF  = zeros(Float64, length(omegas))
    channels = Tuple{Configuration,Configuration,Float64,Array{Float64,1}}[]
    for  i = 1:length(stages)-1
        for  (g, iConf)  in  enumerate(stages[i])
            popAbs = stagePopulations[i] * popsPerStage[i][g]
            popAbs <= 0.0   &&   continue
            for  fConf  in  stages[i+1]
                wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
                ( length(wa) == 1  &&  wa[1][2] == 1 )   ||   continue
                css = Empirical.photoionizationCrossSection(omegas, iConf, fConf, approx)
                kappaBF .+= popAbs .* css
                push!(channels, (iConf, fConf, popAbs, css))
            end
        end
    end

    return( kappaBF, channels )
end


"""
`Plasma.displayMixtureBoundFreeChannels(stream::IO, channels::Array{Tuple{Configuration,Configuration,Float64,Array{Float64,1}},1})`
    ... prints, for every bound-free channel `Plasma.boundFreeOpacityContribution` found, its initial and
        final configuration, the initial group's own absolute population, and the PEAK cross section (the
        maximum of `crossSection(omegas)` over the grid it was evaluated on -- not the full spectrum, which is
        in the opacity table that follows). Nothing is returned.
"""
function  displayMixtureBoundFreeChannels(stream::IO,
                                             channels::Array{Tuple{Configuration,Configuration,Float64,Array{Float64,1}},1})
    println(stream, "\n  Bound-free (photoionization) channels contributing to the mixture opacity:\n")
    println(stream, "  ", TableStrings.hLine(92))
    sa = "    " * TableStrings.center(22, "Initial"; na=2)
    sa = sa     * TableStrings.center(22, "Final"; na=2)
    sa = sa     * TableStrings.center(16, "Population"; na=2)
    sa = sa     * TableStrings.center(20, "Peak cross sec. [a0^2]"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(92))
    if  isempty(channels)
        println(stream, "    (none -- no pair of configurations in adjacent stages differs by exactly one " *
                         "electron in exactly one shell)")
    end
    for  (iConf, fConf, pop, css)  in  channels
        sa = "    " * TableStrings.center(22, string(iConf); na=2)
        sa = sa     * TableStrings.center(22, string(fConf); na=2)
        sa = sa     * @sprintf("%16.6f", pop)
        sa = sa     * @sprintf("%20.6e", maximum(css))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(92))

    return( nothing )
end


"""
`Plasma.classifyResonantGroups(groupsPerStage::Array)`
    ... classifies every `SuperConfigurations.Group` (typed loosely here for the SAME include-order reason
        `superConfigurationsModule()` exists for, see `module-Plasma-inc-superconfiguration-opacity.jl`) in
        every stage but the LAST as BOUND or RESONANT (autoionizing), by comparing its own `meanEnergy`
        against the GROUND energy of the adjacent, one-electron-more-ionized stage --
        `minimum(g.meanEnergy for g in groupsPerStage[i+1])`, the SAME stage-to-stage pairing
        `Plasma.boundFreeOpacityContribution` already uses for its own configuration-level connectivity
        test. A group whose mean energy EXCEEDS that threshold has, on average, enough energy to
        autoionize into the adjacent stage's ground level and is flagged `true`; one below it is an
        ordinary bound excited configuration and is flagged `false`. The LAST stage has no stage above it
        to compare against and is returned all-`false`, by convention rather than by any physics claim.
        Classification ALONE: no width, no line profile, and nothing here is wired into any `kappa` -- an
        `Array{Array{Bool,1},1}` is returned, matching `groupsPerStage`'s own shape.

        **ONE CAVEAT STATED RATHER THAN HIDDEN**, since it bears directly on whether the comparison means
        anything: `Plasma.perform(::SuperConfigurationMixtureOpacityScheme,...)`'s own Stage 0 builds each
        stage's groups from an INDEPENDENT SCF, on that stage's OWN grid and orbitals. The comparison here
        assumes the two stages' absolute energies already sit on a mutually consistent scale -- true for a
        shared nuclear model and adequately converged orbitals (Rule 12), but NOT independently checked by
        this function.
"""
function  classifyResonantGroups(groupsPerStage::Array)
    nStages = length(groupsPerStage)
    flags   = [ zeros(Bool, length(groupsPerStage[i]))   for i = 1:nStages ]
    for  i = 1:nStages-1
        isempty(groupsPerStage[i+1])   &&
            error("Plasma.classifyResonantGroups(): stage $(i+1) has no groups to set a threshold from.")
        threshold = minimum( g.meanEnergy   for g in groupsPerStage[i+1] )
        for  (k, g)  in  enumerate(groupsPerStage[i])
            flags[i][k] = g.meanEnergy > threshold
        end
    end

    return( flags )
end


"""
`Plasma.displayResonantClassification(stream::IO, groupsPerStage::Array, flags::Array{Array{Bool,1},1})`
    ... prints, for every stage but the last, each of its groups' configuration and mean energy, the
        adjacent higher stage's own threshold energy, and whether `Plasma.classifyResonantGroups` flagged
        it bound or resonant. Nothing is returned.
"""
function  displayResonantClassification(stream::IO, groupsPerStage::Array, flags::Array{Array{Bool,1},1})
    println(stream, "\n  Bound/resonant classification of each stage's excited configurations:\n")
    println(stream, "  ", TableStrings.hLine(100))
    sa = "    " * TableStrings.center(10, "Stage"; na=2)
    sa = sa     * TableStrings.center(22, "Configuration"; na=2)
    sa = sa     * TableStrings.center(18, "Mean energy [a.u.]"; na=2)
    sa = sa     * TableStrings.center(18, "Threshold [a.u.]"; na=2)
    sa = sa     * TableStrings.center(14, "Classification"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(100))
    nStages = length(groupsPerStage)
    for  i = 1:nStages-1
        threshold = minimum( g.meanEnergy   for g in groupsPerStage[i+1] )
        for  (k, g)  in  enumerate(groupsPerStage[i])
            sa = "    " * @sprintf("%10i", i)
            sa = sa     * TableStrings.center(22, string(g.conf); na=2)
            sa = sa     * @sprintf("%18.6f", g.meanEnergy)
            sa = sa     * @sprintf("%18.6f", threshold)
            sa = sa     * TableStrings.center(14, flags[i][k] ? "RESONANT" : "bound"; na=2)
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(100))

    return( nothing )
end


"""
`Plasma.freeFreeOpacityContribution(chargesPerStage::Array{Int64,1}, stagePopulations::Array{Float64,1},
                                       omegas::Array{Float64,1}, ne::Float64, kT::Float64)`
    ... sums the free-free (inverse-bremsstrahlung) contribution to kappa(omega), for ONE species, over
        every one of its OWN stages. Unlike `Plasma.boundFreeOpacityContribution` this needs no PAIR of
        adjacent stages: a stage's own ion charge alone sets the Coulomb field a passing free electron
        feels, so each stage contributes independently through
        `Empirical.freeFreeAbsorptionCoefficient(omegas, chargesPerStage[i], ne, kT)`, weighted by that
        stage's own ABSOLUTE population within the species, `stagePopulations[i]`. A neutral stage
        (charge 0) contributes exactly zero, through that function's own Z^2 factor, not as a special case
        here. A kappaFF::Array{Float64,1} is returned, matching `omegas` in length.
"""
function  freeFreeOpacityContribution(chargesPerStage::Array{Int64,1}, stagePopulations::Array{Float64,1},
                                         omegas::Array{Float64,1}, ne::Float64, kT::Float64)
    kappaFF = zeros(Float64, length(omegas))
    for  i = 1:length(chargesPerStage)
        stagePopulations[i] <= 0.0   &&   continue
        kappaFF .+= stagePopulations[i] .* Empirical.freeFreeAbsorptionCoefficient(omegas, chargesPerStage[i], ne, kT)
    end

    return( kappaFF )
end


"""
`Plasma.displayMixtureFreeFreeStages(stream::IO, stages::Array{Tuple{Int64,Int64,Float64,Array{Float64,1}},1})`
    ... prints, for every (species,stage) `Plasma.freeFreeOpacityContribution` actually used (charge > 0),
        its species index, ion charge, absolute population, and PEAK free-free absorption coefficient (the
        maximum of `kappaFF(omegas)` over the grid it was evaluated on -- not the full spectrum, which is in
        the opacity table that follows). Nothing is returned.
"""
function  displayMixtureFreeFreeStages(stream::IO, stages::Array{Tuple{Int64,Int64,Float64,Array{Float64,1}},1})
    println(stream, "\n  Free-free (inverse-bremsstrahlung) stages contributing to the mixture opacity:\n")
    println(stream, "  ", TableStrings.hLine(70))
    sa = "    " * TableStrings.center(10, "Species"; na=2)
    sa = sa     * TableStrings.center(10, "Charge q"; na=2)
    sa = sa     * TableStrings.center(16, "Population"; na=2)
    sa = sa     * TableStrings.center(20, "Peak kappa_ff [a0^2]"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(70))
    if  isempty(stages)
        println(stream, "    (none -- every stage of every species has charge 0)")
    end
    for  (s, q, pop, kff)  in  stages
        sa = "    " * TableStrings.center(10, string(s); na=2)
        sa = sa     * TableStrings.center(10, string(q); na=2)
        sa = sa     * @sprintf("%16.6f", pop)
        sa = sa     * @sprintf("%20.6e", maximum(kff))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(70))

    return( nothing )
end


"""
`Plasma.displayMixtureSpectrumBreakdown(stream::IO, omegas::Array{Float64,1}, kappaBB::Array{Float64,1},
                                           kappaBF::Union{Nothing,Array{Float64,1}},
                                           kappaFF::Union{Nothing,Array{Float64,1}}, kappaTotal::Array{Float64,1},
                                           kappaPbb::Float64, kappaPtotal::Float64, kappaR::Float64)`
    ... prints the mixture's BOUND-BOUND, and whichever of BOUND-FREE/FREE-FREE are active (`kappaBF`/
        `kappaFF` not `nothing`), and TOTAL per-absorber cross sections SIDE BY SIDE at every `omegas`
        entry, each column clearly labelled so no contribution is ever read as one undifferentiated number,
        followed by the Planck mean of the bound-bound-only spectrum (`kappaPbb`, UNCHANGED in meaning from
        `Plasma.displayMixtureSpectrum`'s own number), the Planck mean of the TOTAL spectrum (`kappaPtotal`),
        and the Rosseland mean of the TOTAL spectrum (`kappaR` -- MUCH LESS floor-sensitive now that a
        continuum contribution supplies a smooth, physically decaying tail in place of a Gaussian line's
        literal float64 underflow, but NOT UNCONDITIONALLY floor-free: measured, not assumed, 18-Oct-2026 for
        bound-free alone, a floor sweep on a real case stayed CONSTANT across 9 orders of magnitude of floor
        choice -- `kappaMax*1e-6` down to `kappaMax*1e-15` all agreed to 6 figures on a moderate
        `opacityOmegas` grid -- but widening that SAME grid further past threshold moved the `1e-6` point away
        from the (still-stable) `1e-9`-to`-1e-15` plateau, since a power-law tail can still underflow a FIXED
        multiplicative floor given enough span. `Plasma.rosselandMeanOpacity`'s hardcoded `1e-12` sat inside
        the stable plateau in both tests run, but this is a measured observation on those cases, not a proof
        for every grid). Ends with an EXPLICIT reminder naming whichever of bound-free/free-free is still
        `nothing`, so neither mean is mistaken for a complete continuum opacity unless BOTH are active.
        Nothing is returned.
"""
function  displayMixtureSpectrumBreakdown(stream::IO, omegas::Array{Float64,1}, kappaBB::Array{Float64,1},
                                             kappaBF::Union{Nothing,Array{Float64,1}},
                                             kappaFF::Union{Nothing,Array{Float64,1}}, kappaTotal::Array{Float64,1},
                                             kappaPbb::Float64, kappaPtotal::Float64, kappaR::Float64)
    println(stream, "\n  Superconfiguration MIXTURE opacity spectrum, by CONTRIBUTION (per-absorber cross section):\n")
    width = 14 + 18 + (isnothing(kappaBF) ? 0 : 18) + (isnothing(kappaFF) ? 0 : 18) + 16
    println(stream, "  ", TableStrings.hLine(width))
    sa = "    " * TableStrings.center(14, "omega [a.u.]"; na=2)
    sa = sa     * TableStrings.center(18, "bound-bound [a0^2]"; na=2)
    isnothing(kappaBF)   ||   (sa = sa * TableStrings.center(18, "bound-free [a0^2]"; na=2))
    isnothing(kappaFF)   ||   (sa = sa * TableStrings.center(18, "free-free [a0^2]"; na=2))
    sa = sa     * TableStrings.center(16, "TOTAL [a0^2]"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(width))
    for  k = 1:length(omegas)
        sa = "    " * @sprintf("%14.6f", omegas[k])
        sa = sa     * @sprintf("%18.6e", kappaBB[k])
        isnothing(kappaBF)   ||   (sa = sa * @sprintf("%18.6e", kappaBF[k]))
        isnothing(kappaFF)   ||   (sa = sa * @sprintf("%18.6e", kappaFF[k]))
        sa = sa     * @sprintf("%16.6e", kappaTotal[k])
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(width))
    label = "bound-bound" * (isnothing(kappaBF) ? "" : " + bound-free") * (isnothing(kappaFF) ? "" : " + free-free")
    println(stream, "\n  Planck mean opacity    (bound-bound only)  = $(@sprintf("%.6e", kappaPbb)) [a0^2]")
    println(stream, "  Planck mean opacity    ($label)  = $(@sprintf("%.6e", kappaPtotal)) [a0^2]")
    println(stream, "  Rosseland mean opacity ($label)  = $(@sprintf("%.6e", kappaR)) [a0^2]")
    missingParts = String[]
    isnothing(kappaBF)   &&   push!(missingParts, "bound-free")
    isnothing(kappaFF)   &&   push!(missingParts, "free-free")
    if  !isempty(missingParts)
        println(stream, "\n  " * join(missingParts, " AND ") * " NOT INCLUDED IN ANY COLUMN OR MEAN ABOVE. The")
        println(stream, "  means above are therefore LOWER BOUNDS on the true continuum-inclusive opacity.")
    end

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
        separate module's change.

        STAGE 2.5, ONLY WHEN `scheme.boundFreeApprox !== nothing` (18-Oct-2026): for every species with more
        than one stage, `Plasma.boundFreeOpacityContribution` sums a bound-free (photoionization) contribution
        between every pair of ADJACENT stages, reusing the within-stage populations `SuperConfigurationOpacityScheme`
        already returned in Stage 2/3 -- no extra SCF, no extra grid.

        STAGE 2.6, ONLY WHEN `scheme.freeFreeIncluded` (19-Oct-2026): for every species, EVERY one of its
        stages (no adjacency needed -- a stage's own ion charge alone sets the free-free cross section),
        `Plasma.freeFreeOpacityContribution` sums the free-free (inverse-bremsstrahlung) contribution,
        weighted the same way. This needs a real electron density: `scheme.freeFreeIncluded` with no
        density given (STAGE 1 skipped) raises, naming what is missing.

        The bound-bound-only `ResultKeys` (`SuperConfigurationMixtureOpacitySpectrum`, `...PlanckOpacity`)
        keep their EXACT, UNCHANGED pre-18-Oct meaning regardless of either switch; the bound-free spectrum
        and the free-free spectrum are each reported under their OWN key (present only when their own switch
        is on); and the TOTAL spectrum, its Planck mean, and its Rosseland mean (MUCH LESS floor-sensitive
        than bound-bound alone once EITHER continuum contribution supplies a floor -- see
        `Plasma.displayMixtureSpectrumBreakdown`'s own docstring for the measured, not-unconditional extent
        of that improvement) sum bound-bound plus WHICHEVER of bound-free/free-free are active, and are
        present whenever EITHER is. No contribution is ever merged into another's own key; see
        `Plasma.displayMixtureSpectrumBreakdown`. When NEITHER switch is on (the default), both new stages
        are skipped entirely and every printout/result is BYTE-IDENTICAL to the pre-18-Oct scheme.

        For output=true, a dictionary is returned under the keys "superconfiguration mixture species
        opacities" (`Array{Array{Float64,1},1}`, one BOUND-BOUND spectrum per species), "superconfiguration
        mixture opacity spectrum" (`kappaMix`, BOUND-BOUND ONLY, `Array{Float64,1}`), "superconfiguration
        mixture planck opacity" (`Float64`, of that bound-bound-only spectrum), "superconfiguration mixture
        charges" (`Array{Array{Int64,1},1}`, one array per species), "superconfiguration mixture stage
        populations" (`Array{Array{Float64,1},1}`, matching, summing to 1 WITHIN each species), and
        "superconfiguration mixture electron density" (`Float64`, the ONE value shared by the whole mixture)
        -- ALWAYS present -- plus, ONLY when `scheme.boundFreeApprox !== nothing`: "superconfiguration mixture
        bound-free opacity" (`Array{Float64,1}`) -- plus, ONLY when `scheme.freeFreeIncluded`:
        "superconfiguration mixture free-free opacity" (`Array{Float64,1}`) -- plus, WHENEVER EITHER is
        active: "superconfiguration mixture total opacity spectrum" (bound-bound plus whichever of
        bound-free/free-free are on, `Array{Float64,1}`), "superconfiguration mixture total planck opacity"
        (`Float64`), and "superconfiguration mixture rosseland opacity" (`Float64`, of that same TOTAL
        spectrum).
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
    scheme.freeFreeIncluded  &&  !haveDensity   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): freeFreeIncluded is true, so " *
              "settings.useNumberDensity must be true with settings.density the TOTAL ION number density of " *
              "the WHOLE mixture, in a.u. -- free-free needs a real electron density.")
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
    haveBoundFree      = !isnothing(scheme.boundFreeApprox)
    haveFreeFree       = scheme.freeFreeIncluded
    kappaMix           = zeros(Float64, length(scheme.opacityOmegas))
    kappaMixBF         = zeros(Float64, length(scheme.opacityOmegas))
    kappaMixFF         = zeros(Float64, length(scheme.opacityOmegas))
    speciesKappas      = Array{Float64,1}[]
    allChannels        = Tuple{Configuration,Configuration,Float64,Array{Float64,1}}[]
    allFreeFreeStages  = Tuple{Int64,Int64,Float64,Array{Float64,1}}[]
    for  s = 1:nSpecies
        nm           = scheme.nuclearModels[s]
        stages       = scheme.refConfigsPerSpeciesStages[s]
        kappaSpecies = zeros(Float64, length(scheme.opacityOmegas))
        popsPerStage = Vector{Array{Float64,1}}(undef, nStages[s])
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
            popsPerStage[i]    = stageResults[ResultKeys.SuperConfigurationPopulations]
            kappaSpecies .+= stagePopulations[s][i] .* kappaStage
        end
        push!(speciesKappas, kappaSpecies)
        kappaMix .+= scheme.numberFractions[s] .* kappaSpecies

        # STAGE 2.5, bound-free -- see Plasma.boundFreeOpacityContribution's own docstring. Skipped for a
        # single-stage species (no adjacent stage to ionize into) even when boundFreeApprox is set.
        # Empirical.photoionizationCrossSection reads the nuclear charge from the GLOBAL Defaults (1.0,
        # hydrogen, unless set) -- see Cascade.spectralOpacityContribution(::Cascade.BoundFreeOpacity,...)
        # for the SAME trap met there first. Set from THIS species' own nm.Z, restored afterwards, so a
        # multi-species mixture never silently scores one species' channels against another's Z.
        if  haveBoundFree  &&  nStages[s] > 1
            zOld = Defaults.getDefaults("nuclear: charge")
            Defaults.setDefaults("nuclear: charge", nm.Z)
            try
                kappaSpeciesBF, channelsS = Plasma.boundFreeOpacityContribution(stages, popsPerStage, stagePopulations[s],
                                                                                   scheme.opacityOmegas, scheme.boundFreeApprox)
                kappaMixBF .+= scheme.numberFractions[s] .* kappaSpeciesBF
                append!(allChannels, channelsS)
            finally
                Defaults.setDefaults("nuclear: charge", zOld)
            end
        end

        # STAGE 2.6, free-free -- see Plasma.freeFreeOpacityContribution's own docstring. Unlike bound-free,
        # runs for EVERY species including a single-stage one, since a stage's own charge alone suffices.
        if  haveFreeFree
            kappaSpeciesFF = Plasma.freeFreeOpacityContribution(chargesPerSpecies[s], stagePopulations[s],
                                                                   scheme.opacityOmegas, ne, kT)
            kappaMixFF .+= scheme.numberFractions[s] .* kappaSpeciesFF
            for  i = 1:nStages[s]
                chargesPerSpecies[s][i] <= 0   &&   continue
                push!(allFreeFreeStages, (s, chargesPerSpecies[s][i], stagePopulations[s][i],
                        Empirical.freeFreeAbsorptionCoefficient(scheme.opacityOmegas, chargesPerSpecies[s][i], ne, kT)))
            end
        end
    end

    xGrid  = scheme.opacityOmegas ./ kT
    kappaP = Plasma.planckMeanOpacity(xGrid, kappaMix)

    haveContinuum = haveBoundFree  ||  haveFreeFree
    if  haveContinuum
        kappaMixTotal = kappaMix .+ (haveBoundFree ? kappaMixBF : 0.0) .+ (haveFreeFree ? kappaMixFF : 0.0)
        kappaPtotal   = Plasma.planckMeanOpacity(xGrid, kappaMixTotal)
        kappaR        = Plasma.rosselandMeanOpacity(xGrid, kappaMixTotal)

        if  haveBoundFree
            Plasma.displayMixtureBoundFreeChannels(stdout, allChannels)
            if  printSummary   Plasma.displayMixtureBoundFreeChannels(iostream, allChannels)   end
        end
        if  haveFreeFree
            Plasma.displayMixtureFreeFreeStages(stdout, allFreeFreeStages)
            if  printSummary   Plasma.displayMixtureFreeFreeStages(iostream, allFreeFreeStages)   end
        end
        Plasma.displayMixtureSpectrumBreakdown(stdout, scheme.opacityOmegas, kappaMix, haveBoundFree ? kappaMixBF : nothing,
                                                  haveFreeFree ? kappaMixFF : nothing, kappaMixTotal,
                                                  kappaP, kappaPtotal, kappaR)
        if  printSummary
            Plasma.displayMixtureSpectrumBreakdown(iostream, scheme.opacityOmegas, kappaMix, haveBoundFree ? kappaMixBF : nothing,
                                                      haveFreeFree ? kappaMixFF : nothing, kappaMixTotal,
                                                      kappaP, kappaPtotal, kappaR)
        end

        println("SuperConfigurationMixtureOpacityScheme computation complete ($nSpecies species, electron " *
                "density " * (converged ? "converged" : "NOT CONVERGED") * "; Planck (bound-bound) " *
                "$(@sprintf("%.4e", kappaP)), Planck (total) $(@sprintf("%.4e", kappaPtotal)), Rosseland " *
                "(total) $(@sprintf("%.4e", kappaR)) [a0^2]; bound-free $(haveBoundFree ? "included" : "NOT included"), " *
                "free-free $(haveFreeFree ? "included" : "NOT included")) ...")
    else
        Plasma.displayMixtureSpectrum(stdout, scheme.opacityOmegas, kappaMix, kappaP)
        if  printSummary   Plasma.displayMixtureSpectrum(iostream, scheme.opacityOmegas, kappaMix, kappaP)   end

        println("SuperConfigurationMixtureOpacityScheme computation complete ($nSpecies species, electron " *
                "density " * (converged ? "converged" : "NOT CONVERGED") * ", Planck $(@sprintf("%.4e", kappaP)) " *
                "[a0^2]; Rosseland not reported, see this scheme's own docstring) ...")
    end

    if  output
        results = Base.merge( results, Dict(
            "superconfiguration mixture species opacities" => speciesKappas,
            "superconfiguration mixture opacity spectrum"  => kappaMix,
            "superconfiguration mixture planck opacity"    => kappaP,
            "superconfiguration mixture charges"           => chargesPerSpecies,
            "superconfiguration mixture stage populations" => stagePopulations,
            "superconfiguration mixture electron density"  => ne) )
        if  haveBoundFree
            results["superconfiguration mixture bound-free opacity"]       = kappaMixBF
        end
        if  haveFreeFree
            results["superconfiguration mixture free-free opacity"]       = kappaMixFF
        end
        if  haveContinuum
            results["superconfiguration mixture total opacity spectrum"]   = kappaMixTotal
            results["superconfiguration mixture total planck opacity"]    = kappaPtotal
            results["superconfiguration mixture rosseland opacity"]       = kappaR
        end
    end

    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
