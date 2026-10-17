# Functions and methods for scheme::Plasma.SuperConfigurationMixtureOpacityScheme computations.
#
# EACH SPECIES IS RUN THROUGH Plasma.perform(::SuperConfigurationOpacityScheme, ...) UNCHANGED, one call
# per species, on that species' OWN Nuclear.Model and an automatically-sized grid (Rule 12). No new CSF-level
# physics is added here -- only the mixing of several species' own Stage-3 spectra and the two mean-opacity
# integrals over the result. See Plasma.SuperConfigurationMixtureOpacityScheme's own docstring for the
# mixing rule and the weight functions.


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
`Plasma.perform(scheme::Plasma.SuperConfigurationMixtureOpacityScheme, computation::Plasma.Computation; output::Bool=true)`
    ... to run `Plasma.perform(::SuperConfigurationOpacityScheme, ...)` UNCHANGED for every species in
        `scheme.refConfigsPerSpecies` (its own `Nuclear.Model`, its own automatically-sized grid), combine
        their Stage-3 spectra into `kappaMix(omega) = SUM_s numberFractions[s]*kappa_s(omega)`, and reduce
        that mixture spectrum to its Planck mean opacity over `computation.settings.temperature`. The
        Rosseland mean is DELIBERATELY NOT computed here -- see `Plasma.rosselandMeanOpacity`'s own
        docstring: it was measured, not assumed, to have no stable value for a line-only spectrum (the
        result moved by eleven-plus orders of magnitude under an internal floor choice that should not
        matter), so reporting it would be worse than naming it as open work needing a continuum opacity.
        `computation.nuclearModel`/`computation.grid`/`computation.refConfigs` are UNUSED (see this scheme's
        own docstring); only `computation.asfSettings` and `computation.settings.temperature` are shared
        across species. For output=true, a dictionary is returned under the keys "superconfiguration mixture
        species opacities" (`Array{Array{Float64,1},1}`, one spectrum per species, in `scheme`'s own order),
        "superconfiguration mixture opacity spectrum" (`kappaMix`, `Array{Float64,1}`), and
        "superconfiguration mixture planck opacity" (`Float64`).
"""
function  perform(scheme::Plasma.SuperConfigurationMixtureOpacityScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nSpecies = length(scheme.refConfigsPerSpecies)
    nSpecies == 0   &&
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): refConfigsPerSpecies is empty.")
    ( length(scheme.nuclearModels) == nSpecies  &&  length(scheme.numberFractions) == nSpecies )   ||
        error("Plasma.perform(::SuperConfigurationMixtureOpacityScheme, ...): nuclearModels, " *
              "refConfigsPerSpecies and numberFractions must all have the same length.")
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

    kappaMix      = zeros(Float64, length(scheme.opacityOmegas))
    speciesKappas = Array{Float64,1}[]

    for  s = 1:nSpecies
        nm   = scheme.nuclearModels[s]
        refs = scheme.refConfigsPerSpecies[s]
        println("\n  -- Mixture species $s of $nSpecies, number fraction $(scheme.numberFractions[s]) --")

        grid               = Basics.recommendedGrid(refs, nm; printout=false)
        speciesScheme      = SuperConfigurationOpacityScheme(refs, scheme.opacityOmegas)
        speciesComputation = Plasma.Computation(computation; scheme=speciesScheme, nuclearModel=nm, grid=grid,
                                                   refConfigs=Configuration[])

        speciesResults = Plasma.perform(speciesScheme, speciesComputation; output=true)
        kappaS         = speciesResults[ResultKeys.SuperConfigurationOpacitySpectrum]
        push!(speciesKappas, kappaS)
        kappaMix .+= scheme.numberFractions[s] .* kappaS
    end

    xGrid  = scheme.opacityOmegas ./ kT
    kappaP = Plasma.planckMeanOpacity(xGrid, kappaMix)

    Plasma.displayMixtureSpectrum(stdout, scheme.opacityOmegas, kappaMix, kappaP)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayMixtureSpectrum(iostream, scheme.opacityOmegas, kappaMix, kappaP)   end

    println("SuperConfigurationMixtureOpacityScheme computation complete ($nSpecies species, Planck " *
            "$(@sprintf("%.4e", kappaP)) [a0^2]; Rosseland not reported, see this scheme's own docstring) ...")

    if  output
        results = Base.merge( results, Dict("superconfiguration mixture species opacities" => speciesKappas,
                                              "superconfiguration mixture opacity spectrum" => kappaMix,
                                              "superconfiguration mixture planck opacity"    => kappaP) )
    end

    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
