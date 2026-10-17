# Functions and methods for scheme::Plasma.SuperConfigurationIonizationScheme computations.
#
# SHARES `superConfigurationsModule()`, defined once in module-Plasma-inc-superconfiguration-opacity.jl for
# the SAME include-order reason explained there -- visible here without redefinition, since both files are
# `include`d into the one `module Plasma ... end` block and share its namespace.


"""
`Plasma.displayIonizationBalance(stream::IO, charges::Array{Int64,1}, groupsPerStage::Array,
                                    stagePopulations::Array{Float64,1}, electronDensity::Float64)`
    ... prints, for every ionization stage, its charge, the configuration of its lowest-meanEnergy group,
        how many groups it carries, and its own population fraction (summed over every group in that
        stage); then the converged electron number density. Nothing is returned.
"""
function  displayIonizationBalance(stream::IO, charges::Array{Int64,1}, groupsPerStage::Array,
                                     stagePopulations::Array{Float64,1}, electronDensity::Float64)
    println(stream, "\n  Superconfiguration ionization balance:\n")
    println(stream, "  ", TableStrings.hLine(80))
    sa = "    " * TableStrings.center(10, "Charge q"; na=2)
    sa = sa     * TableStrings.center(26, "Lowest configuration"; na=2)
    sa = sa     * TableStrings.center(12, "# groups"; na=2)
    sa = sa     * TableStrings.center(16, "Population"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(80))
    for  (q, groups, pop) in zip(charges, groupsPerStage, stagePopulations)
        minIdx = argmin( [ g.meanEnergy  for g in groups ] )
        sa = "    " * TableStrings.center(10, string(q); na=2)
        sa = sa     * TableStrings.center(26, string(groups[minIdx].conf); na=2)
        sa = sa     * TableStrings.center(12, string(length(groups)); na=2)
        sa = sa     * @sprintf("%16.6f", pop)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(80))
    println(stream, "    converged electron number density = ", @sprintf("%.6e", electronDensity), " [a.u.]")

    return( nothing )
end


"""
`Plasma.perform(scheme::Plasma.SuperConfigurationIonizationScheme, computation::Plasma.Computation;
                  output::Bool=true)`
    ... to build one `SuperConfigurations.Group` per configuration in EVERY stage of
        `scheme.refConfigsPerStage` (one SHARED orbital set PER STAGE, since different stages have
        different electron counts and cannot share orbitals the way configurations WITHIN one stage do),
        and to solve the grand-canonical Saha-Boltzmann equilibrium self-consistently for the electron
        density at `computation.settings.temperature` [Hartree] and `computation.settings.density` [a.u.,
        the TOTAL ION number density; `computation.settings.useNumberDensity` must be true]. See this
        scheme's own docstring for the equations and their scope. For output=true, a dictionary is returned
        under the keys "superconfiguration ionization groups" (`Array{Array{SuperConfigurations.Group,1},1}`,
        one inner array per stage), "superconfiguration ionization charges" (`Array{Int64,1}`, matching),
        "superconfiguration ionization stage populations" (`Array{Float64,1}`, summing to 1),
        "superconfiguration ionization group populations" (`Array{Array{Float64,1},1}`, matching the groups,
        summing to 1 OVER ALL stages together), and "superconfiguration ionization electron density"
        (`Float64`, the converged value).
"""
function  perform(scheme::Plasma.SuperConfigurationIonizationScheme, computation::Plasma.Computation;
                    output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    isempty(scheme.refConfigsPerStage)   &&
        error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): scheme.refConfigsPerStage is empty.")
    computation.settings.useNumberDensity   ||
        error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): settings.useNumberDensity must be " *
              "true -- settings.density is read as the TOTAL ION number density, in a.u.")

    nm   = computation.nuclearModel
    grid = computation.grid
    kT   = computation.settings.temperature
    kT <= 0.0   &&
        error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): temperature must be positive.")
    nIon = computation.settings.density
    nIon <= 0.0   &&
        error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): density must be positive.")

    SC = superConfigurationsModule()

    nStages        = length(scheme.refConfigsPerStage)
    groupsPerStage = Array{Any,1}(undef, nStages)
    charges        = zeros(Int64, nStages)
    for  i = 1:nStages
        refConfigs = scheme.refConfigsPerStage[i]
        isempty(refConfigs)   &&
            error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): stage $i has no reference " *
                  "configurations.")
        any( c -> c.NoElectrons != refConfigs[1].NoElectrons, refConfigs )   &&
            error("Plasma.perform(::SuperConfigurationIonizationScheme, ...): stage $i mixes configurations " *
                  "of different electron count; one stage must be ONE electron count.")

        wb       = SelfConsistent.performSCF(refConfigs, nm, grid, computation.asfSettings; printout=false)
        orbitals = wb.levels[1].basis.orbitals
        groupsPerStage[i] = [ SC.buildGroup(conf, orbitals, nm, grid, computation.asfSettings)
                              for conf in refConfigs ]
        charges[i]        = round(Int64, nm.Z) - refConfigs[1].NoElectrons
    end

    # ONE shared energy reference for every stage and group, so that charge-state differences come out
    # directly from the raw SCF energies -- no separate "ionization energy" is ever extracted, exactly as
    # Plasma.SahaBoltzmannScheme's own grand-canonical formula needs none either (see this scheme's own
    # docstring).
    Eref = minimum( minimum(g.meanEnergy for g in groups)  for groups in groupsPerStage )

    # SELF-CONSISTENT fixed-point iteration on the electron chemical potential / density, via charge
    # neutrality -- the SAME grand-canonical Saha-Boltzmann equations Plasma.SahaBoltzmannScheme already
    # uses, reimplemented here (not called there) since that scheme carries isotopic-mixture and
    # ionization-potential-depression machinery this one deliberately does not need, and is itself not in
    # the test suite.
    ne               = nIon          # a defensible first guess: at most one electron per ion on average
    stagePopulations = zeros(Float64, nStages)
    groupPopulations = [ zeros(Float64, length(groups))  for groups in groupsPerStage ]
    converged        = false
    for  iteration = 1:100
        eLambda3 = (2*pi/kT)^1.5
        muE      = -kT * log( 2.0 / (ne * eLambda3) )

        weights = [ [ g.degeneracy * exp( -(muE*charges[i] + g.meanEnergy - Eref)/kT )   for g in groupsPerStage[i] ]
                    for i = 1:nStages ]
        Zgrand  = sum( sum(w)  for w in weights )

        for  i = 1:nStages
            groupPopulations[i] = weights[i] ./ Zgrand
            stagePopulations[i] = sum(weights[i]) / Zgrand
        end

        neNew = nIon * sum( charges[i]*stagePopulations[i]   for i = 1:nStages )
        if  abs(neNew - ne) < 1.0e-10 * ne
            converged = true;   ne = neNew;   break
        end
        ne = neNew
    end
    converged   ||
        @warn("Plasma.perform(::SuperConfigurationIonizationScheme, ...): the electron-density self-" *
              "consistency did NOT converge in 100 iterations; ne = $ne a.u. is the last iterate, not a " *
              "solution.")

    Plasma.displayIonizationBalance(stdout, charges, groupsPerStage, stagePopulations, ne)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayIonizationBalance(iostream, charges, groupsPerStage, stagePopulations, ne)   end

    println("SuperConfigurationIonizationScheme computation complete (equilibrium ionization balance " *
            "across $nStages stages, electron density " * (converged ? "converged" : "NOT CONVERGED") * ") ...")

    if  output
        results = Base.merge( results, Dict(
            "superconfiguration ionization groups"            => groupsPerStage,
            "superconfiguration ionization charges"           => charges,
            "superconfiguration ionization stage populations" => stagePopulations,
            "superconfiguration ionization group populations" => groupPopulations,
            "superconfiguration ionization electron density"  => ne) )
    end

    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
