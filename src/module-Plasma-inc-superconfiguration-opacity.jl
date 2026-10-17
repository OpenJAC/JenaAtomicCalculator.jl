# Functions and methods for scheme::Plasma.SuperConfigurationOpacityScheme computations.
#
# THIS FILE DOES NOT `using ..SuperConfigurations` AT MODULE TOP LEVEL, unlike most sibling -inc- files'
# dependencies. `SuperConfigurations` is included in JenaAtomicCalculator.jl AFTER module-Plasma.jl -- it
# depends on BasicsAZ.meanEnergyInClosedForm, and BasicsAZ is deliberately included near the very end of
# JenaAtomicCalculator.jl, after Cascade/ForPedestrians/RacahAlgebra and more, since its own docstring says
# it "contains methods from the Basics module but with reference to many other modules" -- so moving it (or
# SuperConfigurations) earlier to let Plasma `using` it directly risks breaking dependencies this session
# has not audited. `superConfigurationsModule(...)` below resolves the sibling module by RUNTIME reflection
# instead, which only runs when `perform` is actually CALLED -- by which point every module has long since
# been included -- so no existing include order needs to change.
function  superConfigurationsModule()
    return( getfield(parentmodule(@__MODULE__), :SuperConfigurations) )
end


"""
`Plasma.displayResults(stream::IO, groups::Array{Any,1}, weights::Array{Float64,1}, Z::Float64)`
    ... prints, for each group (each a `SuperConfigurations.Group`, typed loosely here for the SAME
        include-order reason `superConfigurationsModule()` exists above), its configuration, mean energy,
        width (sqrt variance), degeneracy, and Boltzmann population fraction weights[i]/Z. Nothing is
        returned.
"""
function  displayResults(stream::IO, groups::Array, weights::Array{Float64,1}, Z::Float64)
    println(stream, "\n  Superconfiguration groups and their Boltzmann population fractions:\n")
    println(stream, "  ", TableStrings.hLine(100))
    sa = "    " * TableStrings.center(22, "Configuration"; na=2)
    sa = sa     * TableStrings.center(18, "Mean energy [a.u.]"; na=2)
    sa = sa     * TableStrings.center(16, "Width [a.u.]"; na=2)
    sa = sa     * TableStrings.center(12, "Degeneracy"; na=2)
    sa = sa     * TableStrings.center(14, "Population"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(100))
    for  (group, w) in zip(groups, weights)
        sa = "    " * TableStrings.center(22, string(group.conf); na=2)
        sa = sa     * @sprintf("%18.6f", group.meanEnergy)
        sa = sa     * @sprintf("%16.6f", sqrt(max(group.variance, 0.0)))
        sa = sa     * @sprintf("%12.3f", group.degeneracy)
        sa = sa     * @sprintf("%14.6f", w/Z)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(100))

    return( nothing )
end


"""
`Plasma.displayTransitionArrays(stream::IO, arrays::Array, gfValues::Array{Float64,1},
                                   gfEffValues::Array{Float64,1})`
    ... prints, for each `SuperConfigurations.TransitionArray` found among `scheme.refConfigs` (typed
        loosely here for the SAME include-order reason `superConfigurationsModule()` exists above), the
        initial and final configuration, the mean transition energy, its width (sqrt variance), the array's
        OWN total oscillator-strength sum `gf = (2/3) * meanEnergy * totalStrength`, and that same `gf`
        weighted by the initial configuration's Boltzmann population fraction -- the quantity that actually
        enters a plasma opacity sum. Nothing is returned.
"""
function  displayTransitionArrays(stream::IO, arrays::Array, gfValues::Array{Float64,1}, gfEffValues::Array{Float64,1})
    println(stream, "\n  Superconfiguration E1 transition arrays and their oscillator-strength sums:\n")
    println(stream, "  ", TableStrings.hLine(118))
    sa = "    " * TableStrings.center(20, "Initial"; na=2)
    sa = sa     * TableStrings.center(20, "Final"; na=2)
    sa = sa     * TableStrings.center(18, "Mean energy [a.u.]"; na=2)
    sa = sa     * TableStrings.center(16, "Width [a.u.]"; na=2)
    sa = sa     * TableStrings.center(14, "gf (array)"; na=2)
    sa = sa     * TableStrings.center(16, "Population x gf"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(118))
    for  (ta, gf, gfEff) in zip(arrays, gfValues, gfEffValues)
        sa = "    " * TableStrings.center(20, string(ta.confInitial); na=2)
        sa = sa     * TableStrings.center(20, string(ta.confFinal); na=2)
        sa = sa     * @sprintf("%18.6f", ta.meanEnergy)
        sa = sa     * @sprintf("%16.6f", sqrt(max(ta.variance, 0.0)))
        sa = sa     * @sprintf("%14.6f", gf)
        sa = sa     * @sprintf("%16.6f", gfEff)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(118))

    return( nothing )
end


"""
`Plasma.perform(scheme::Plasma.SuperConfigurationOpacityScheme, computation::Plasma.Computation; output::Bool=true)`
    ... to build one `SuperConfigurations.Group` per configuration in `scheme.refConfigs`, on a SHARED
        orbital set obtained from `SelfConsistent.performSCF` over all of `refConfigs` together, to report
        each group's Boltzmann population fraction at `computation.settings.temperature`, and -- STAGE 2 --
        to build a `SuperConfigurations.TransitionArray` for every PAIR of `refConfigs` connected by a
        single-electron E1 transition (opposite parity), reporting each array's own oscillator-strength sum
        `gf` and that `gf` weighted by the lower configuration's population. See
        `Plasma.SuperConfigurationOpacityScheme`'s own docstring for what this does and does not yet
        compute -- in particular, these are per-array Gaussian line-shape PARAMETERS, not yet a single
        summed opacity as a function of photon energy. For output=true, a dictionary is returned under the
        keys "superconfiguration groups" (`Array{SuperConfigurations.Group,1}`), "superconfiguration
        populations" (the matching Boltzmann population fractions, `Array{Float64,1}`, summing to 1),
        "superconfiguration transition arrays" (`Array{SuperConfigurations.TransitionArray,1}`, one per
        E1-connected pair), "superconfiguration transition gf" (each array's own `gf`, `Array{Float64,1}`),
        and "superconfiguration opacity contributions" (that `gf` weighted by the lower configuration's
        population, `Array{Float64,1}` -- the quantity a later summed spectrum would be built from).
"""
function  perform(scheme::Plasma.SuperConfigurationOpacityScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    isempty(scheme.refConfigs)   &&
        error("Plasma.perform(::SuperConfigurationOpacityScheme, ...): scheme.refConfigs is empty.")

    nm   = computation.nuclearModel
    grid = computation.grid

    wb       = SelfConsistent.performSCF(scheme.refConfigs, nm, grid, computation.asfSettings; printout=false)
    orbitals = wb.levels[1].basis.orbitals

    SC     = superConfigurationsModule()
    groups = [ SC.buildGroup(conf, orbitals, nm, grid, computation.asfSettings)   for conf in scheme.refConfigs ]

    kT = computation.settings.temperature
    kT <= 0.0   &&   error("Plasma.perform(::SuperConfigurationOpacityScheme, ...): temperature must be positive.")

    minEnergy = minimum( group.meanEnergy  for group in groups )
    weights   = [ group.degeneracy * exp( -(group.meanEnergy - minEnergy)/kT )   for group in groups ]
    Z           = sum(weights)
    populations = weights ./ Z

    Plasma.displayResults(stdout, groups, weights, Z)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayResults(iostream, groups, weights, Z)   end

    # STAGE 2: one TransitionArray for every pair connected by a single-electron E1 transition (opposite
    # parity); a same-parity pair is simply not an E1 array and is skipped, not an error here.
    n          = length(scheme.refConfigs)
    tArrays    = Any[];   gfValues = Float64[];   gfEffValues = Float64[]
    for  i = 1:n
        for  j = i+1:n
            parityI = Basics.extractFromConfiguration(Basics.GetParity(), scheme.refConfigs[i])
            parityJ = Basics.extractFromConfiguration(Basics.GetParity(), scheme.refConfigs[j])
            parityI == parityJ   &&   continue

            if  groups[i].meanEnergy <= groups[j].meanEnergy
                  idxInitial, idxFinal = i, j
            else  idxInitial, idxFinal = j, i
            end
            ta = SC.buildTransitionArray(scheme.refConfigs[idxFinal], scheme.refConfigs[idxInitial], orbitals,
                                           nm, grid, computation.asfSettings)
            gf    = 2/3 * ta.meanEnergy * ta.totalStrength
            gfEff = populations[idxInitial] * gf
            push!(tArrays, ta);   push!(gfValues, gf);   push!(gfEffValues, gfEff)
        end
    end

    Plasma.displayTransitionArrays(stdout, tArrays, gfValues, gfEffValues)
    if  printSummary   Plasma.displayTransitionArrays(iostream, tArrays, gfValues, gfEffValues)   end

    println("SuperConfigurationOpacityScheme computation complete (Stage 2: population weighting plus " *
            "per-array oscillator-strength sums for every E1-connected pair -- see this scheme's own " *
            "docstring for what is not yet built from them) ...")

    if  output
        results = Base.merge( results, Dict("superconfiguration groups"                => groups,
                                              "superconfiguration populations"           => populations,
                                              "superconfiguration transition arrays"     => tArrays,
                                              "superconfiguration transition gf"         => gfValues,
                                              "superconfiguration opacity contributions" => gfEffValues) )
    end

    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
