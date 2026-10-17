
"""
`module  JAC.SuperConfigurations`
    ... a submodel of JAC that gives the Gaussian (two-moment) statistical description of ONE configuration's
        UNRESOLVED level spectrum -- a mean energy, a variance, and a total statistical weight -- reported
        WITHOUT diagonalizing the configuration's CI matrix. This is the building block a superconfiguration
        opacity (and, later, ionization-balance or NLTE-kinetics) scheme needs for an ion whose open-shell
        CSF count is too large for ordinary, level-by-level JAC machinery to resolve one level at a time.

        ONE-WAY DEPENDENCY, the same discipline every statistical/empirical module this project has built
        keeps: no process or `Plasma` module is called FROM here, and this module is never called FROM any
        process module either. A caller -- a future `Plasma.SuperConfigurationOpacityScheme` -- builds a
        `Group` and reads its `meanEnergy`/`variance`/`degeneracy` fields directly.

        THE MEAN IS REUSED, NOT REIMPLEMENTED. `BasicsAZ.meanEnergyInClosedForm` already gives a
        configuration's (2J+1)-weighted average energy in closed form, touching no CSF at all, independently
        validated to 1.6e-9 relative against the exact trace on nine configurations through `[Xe] 4f^7`
        (`module-BasicsAZ-inc-AG.jl`). This module calls it directly for `Group.meanEnergy` -- which is also
        why this module is included AFTER `module-BasicsAZ.jl` in `JenaAtomicCalculator.jl`, not alongside
        the other statistical modules built earlier.

        THE VARIANCE IS DELIBERATELY NOT YET THE SAME KIND OF CLOSED FORM. A true shell-pair-level closed
        form for the SECOND moment needs genuine Racah-algebra sum-rule reduction -- 6-j-symbol-level
        identities for sums of PRODUCTS of two-body matrix elements -- which this module's own first version
        chose not to derive from memory without a reference to check it against: the risk of a subtly wrong
        combinatorial factor, one that passes every test that does not happen to probe it, was judged higher
        than the value of rushing it. WHAT IS IMPLEMENTED INSTEAD IS EXACT AND STILL USEFUL:
        `varianceFromCsfMatrix` builds the full CSF-pair Hamiltonian matrix (the SAME basis
        `Basics.determineMeanEnergy` builds) and reads the variance from Tr(H^2), avoiding the O(n^3)
        DIAGONALIZATION step but not the O(n^2) matrix CONSTRUCTION -- exactly the role
        `Basics.determineMeanEnergy` played before `meanEnergyInClosedForm` replaced it, and the natural
        internal validation target for whatever shell-pair closed form is derived later. For a configuration
        too large even for this, the module is not yet usable; that is Stage 2, named here rather than hidden.

        STAGE 2 STARTS HERE. `SuperConfigurations.TransitionArray` gives the SAME two-moment description for
        the UNRESOLVED electric-dipole (E1) transition array BETWEEN two configurations -- the building
        block `Plasma.SuperConfigurationOpacityScheme` needs to turn a population weighting into an actual
        opacity. It is built with the SAME discipline as `Group.variance`: EXACT, at CSF-pair level (no
        diagonalization of either configuration, and no closed-form Bauche-Arnoult moment formula reached
        for without a reference to check it against), via `SuperConfigurations.buildTransitionArray`. See
        that function's own docstring for the trace identities it uses and why no extra statistical-weight
        factor is needed beyond the dipole reduced matrix element's own Wigner-Eckart normalization.

        PRIORITY ITEM 79, CLOSED 17-Oct-2026. The blind double loop over CSF pairs in
        `buildHamiltonianMatrix`/`buildDipoleMatrix` was measured to need ~6700 s at 2725 CSFs (Dy3+
        `4f^9`/`4f^8 5d`) -- the regime real superconfiguration work is FOR. The fix is NOT skipping
        zero-coupling CSF pairs, which a related register entry for the EOL solver already found plateaus
        at a constant factor: it is memoizing the one-particle RADIAL integral (`RadialIntegrals.GrantIab`,
        `InteractionStrength.eMultipole`) BY ORBITAL PAIR -- `oneParticleIntegralCached`,
        `matrixElementCached`, and the matching cache inside `buildDipoleMatrix`. Measured on Co2+
        `3d^7`/`3d^6 4p`: the same orbital pair recurs about 12 times per 400 same-block CSF pairs, and that
        redundancy GROWS with the configuration (CSF pairs scale as n^2, distinct orbital pairs do not) --
        the opposite of the EOL's constant-factor finding, and why this lever actually changes the order of
        the cost rather than a prefactor.

        SCOPE, DELIBERATELY NARROW, the same choice `AdjustedLevels` made for the same reason: ONE ordinary
        (non-relativistic) `Configuration` at a time, pure Coulomb (no Breit/QED).
"""
module SuperConfigurations

using ..Basics, ..BasicsAZ, ..ManyElectron, ..Nuclear, ..Radial, ..InteractionStrength, ..InteractionStrengthQED,
      ..RadialIntegrals, ..Hamiltonian, ..SpinAngular, LinearAlgebra


"""
`struct  SuperConfigurations.Group`
    ... the Gaussian (two-moment) statistical description of one configuration's UNRESOLVED level spectrum:
        its mean energy, variance, and total (2J+1)-weighted degeneracy (the statistical weight the group
        represents), computed without diagonalizing its CI matrix. See this module's own docstring for
        exactly how `meanEnergy` and `variance` are obtained, and what each still costs.

    + conf          ::Configuration   ... the configuration this group describes.
    + meanEnergy    ::Float64          ... (2J+1)-weighted mean energy, in Hartree.
    + variance      ::Float64          ... (2J+1)-weighted variance of the level energies, in Hartree^2.
    + degeneracy    ::Float64          ... total (2J+1)-weighted statistical weight, SUM_csf (2*J_csf+1).
"""
struct  Group
    conf            ::Configuration
    meanEnergy      ::Float64
    variance        ::Float64
    degeneracy      ::Float64
end


"""
`Base.show(io::IO, group::SuperConfigurations.Group)`
    ... prints the configuration, its mean energy, its width (sqrt(variance), the natural quantity to quote
        for a spread) and its degeneracy. Nothing is returned.
"""
function  Base.show(io::IO, group::Group)
    println(io, "SuperConfigurations.Group: ", string(group.conf))
    println(io, "  mean energy           = $(group.meanEnergy) Hartree")
    println(io, "  width (sqrt variance) = $(sqrt(max(group.variance, 0.0))) Hartree")
    println(io, "  degeneracy (SUM 2J+1) = $(group.degeneracy)")
end


"""
`struct  SuperConfigurations.TransitionArray`
    ... the Gaussian (two-moment) statistical description of the UNRESOLVED electric-dipole (E1) transition
        array between two configurations' level spectra: the line-strength-weighted mean transition energy,
        its variance, and the array's total line strength, computed without diagonalizing either
        configuration's CI matrix. See `SuperConfigurations.buildTransitionArray` for exactly how these are
        obtained, and this module's own docstring for why they are EXACT rather than a closed form.

        WHETHER A GAUSSIAN IS A GOOD MODEL FOR THE TRUE DISTRIBUTION -- not just whether these two moments
        are computed correctly -- was CALIBRATED, not assumed, 17-Oct-2026, on Co2+ [Ar] 3d^7 -> 3d^6 4p
        (3420 individual CSF-pair lines, exact eigenbasis via direct diagonalization, used only for this
        one-time check and nowhere in the production path): skewness +0.28, excess kurtosis +1.11 against
        the Gaussian's 0 and 0, with no single line carrying more than 2.2 % of the total strength (the top
        five together 8.4 %). So at this CSF count a Gaussian is a reasonable, mildly heavy-tailed and
        mildly asymmetric approximation, not an exact one, and NOT an artefact of a few outlier lines
        dominating the sum. This is ONE measurement on ONE array, not a theorem -- re-check before trusting
        it at a very different size or shell structure, in particular before building Stage 3 (a single
        summed opacity spectrum) on top of it.

    + confFinal     ::Configuration   ... the configuration the transition ends in.
    + confInitial   ::Configuration   ... the configuration the transition starts from.
    + meanEnergy    ::Float64          ... line-strength-weighted mean of (E_final - E_initial), in Hartree;
                                            POSITIVE for an allowed absorption transition to a higher complex.
    + variance      ::Float64          ... line-strength-weighted variance of the transition energy, in
                                            Hartree^2.
    + totalStrength ::Float64          ... the array's total E1 line strength, SUM over every CSF pair of
                                            |<confFinal CSF || D^(1) || confInitial CSF>|^2, in atomic units.
"""
struct  TransitionArray
    confFinal       ::Configuration
    confInitial     ::Configuration
    meanEnergy      ::Float64
    variance        ::Float64
    totalStrength   ::Float64
end


"""
`Base.show(io::IO, ta::SuperConfigurations.TransitionArray)`
    ... prints the two configurations, the mean transition energy, its width (sqrt(variance)) and the
        array's total line strength. Nothing is returned.
"""
function  Base.show(io::IO, ta::TransitionArray)
    println(io, "SuperConfigurations.TransitionArray:  ", string(ta.confInitial), "  ->  ", string(ta.confFinal))
    println(io, "  mean transition energy = $(ta.meanEnergy) Hartree")
    println(io, "  width (sqrt variance)  = $(sqrt(max(ta.variance, 0.0))) Hartree")
    println(io, "  total line strength    = $(ta.totalStrength)")
end


"""
`SuperConfigurations.buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital};
                                         subshellList::Union{Nothing,Array{Subshell,1}}=nothing)`
    ... a private helper that builds the SAME relativistic-CSF basis `Basics.determineMeanEnergy` builds, so
        that `varianceFromCsfMatrix` and `BasicsAZ.meanEnergyInClosedForm` describe EXACTLY the same set of
        states. If `subshellList` is given, `conf`'s CSFs are built on THAT list instead of generating one
        from `conf` alone -- the one case this is needed is `buildTransitionArray`, where two DIFFERENT
        configurations' CSFs must be indexed against their COMBINED subshell list before a transition
        operator can connect them. A (basis::Basis, csfList::Array{CsfR,1}) tuple is returned.
"""
function  buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital};
                              subshellList::Union{Nothing,Array{Subshell,1}}=nothing)
    relconfList = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
    subshells   = isnothing(subshellList) ? Basics.generateSubshellList(relconfList) : subshellList
    csfList     = CsfR[]
    for  relconf in relconfList    append!(csfList, Basics.generateCsfRs(relconf, subshells))    end
    NoElectrons = sum( csfList[1].occupation )

    coreSubshellList = Subshell[]
    for  k = 1:length(subshells)
        mocc = Basics.subshell_2j(subshells[k]) + 1;    isFilled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    isFilled = false;    break    end
        end
        isFilled   &&   push!(coreSubshellList, subshells[k])
    end
    basis = Basis(true, NoElectrons, subshells, csfList, coreSubshellList, orbitals)

    return( (basis, csfList) )
end


"""
`SuperConfigurations.oneParticleIntegralCached(a::Subshell, b::Subshell, orbitals::Dict{Subshell,Orbital},
                                                  nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                                                  potential::Radial.Potential,
                                                  cache::Dict{Tuple{Subshell,Subshell},Float64})`
    ... a private helper returning the SAME one-particle quantity `Hamiltonian.matrixElement`'s own inner
        loop computes per coefficient -- `RadialIntegrals.GrantIab` plus the QED local term when
        `settings.qedModel` asks for one -- but MEMOIZED by the orbital-pair key `(a,b)`. THIS IS THE FIX
        FOR PRIORITY ITEM 79. Measured 17-Oct-2026 on Co2+ `3d^7`/`3d^6 4p`: across 400 sampled same-block
        CSF pairs, the SAME orbital pair recurred about 12 TIMES, because every CSF pair's one-particle
        coefficients can only ever name orbitals from the small, SHARED active set, while the number of CSF
        PAIRS grows as the square of the CSF count -- so the redundancy grows with the configuration, not a
        constant. `GrantIab` alone was 82 % of the sampled cost (measured ~199 microseconds per call
        against ~12.5 for the angular coefficients feeding it), so caching IT specifically, by orbital pair
        rather than by CSF pair, is what changes the ORDER of the cost -- unlike skipping zero-coupling CSF
        pairs, which the EOL solver's own register entry already found plateaus at a constant factor rather
        than vanishing. A value::Float64 is returned.
"""
function  oneParticleIntegralCached(a::Subshell, b::Subshell, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                      grid::Radial.Grid, settings::AsfSettings, potential::Radial.Potential,
                                      cache::Dict{Tuple{Subshell,Subshell},Float64})
    key = (a, b)
    haskey(cache, key)   &&   return( cache[key] )

    val = RadialIntegrals.GrantIab(orbitals[a], orbitals[b], grid, potential)
    if  settings.qedModel != NoneQed()
        val = val + InteractionStrengthQED.qedLocal(orbitals[a], orbitals[b], nm, settings.qedModel, potential, grid)
    end
    cache[key] = val

    return( val )
end


"""
`SuperConfigurations.matrixElementCached(basis::Basis, r::Int64, s::Int64, nm::Nuclear.Model, grid::Radial.Grid,
                                            settings::AsfSettings, potential::Radial.Potential,
                                            cache::InteractionStrength.XLCache,
                                            grantCache::Dict{Tuple{Subshell,Subshell},Float64})`
    ... a private helper, EXACTLY `Hamiltonian.matrixElement`'s own logic -- the two-particle
        `XL_Coulomb`/`XL_Breit` terms already reuse `cache` across every CSF pair and are unchanged here --
        with its one-particle term routed through `oneParticleIntegralCached` instead of calling
        `RadialIntegrals.GrantIab` fresh every time. Kept as a LOCAL copy rather than a change to
        `Hamiltonian.matrixElement` itself, since that function is shared by the whole RAS/CI machinery and
        extending it is a larger, separate decision this module does not make for it -- the fix here stays
        inside this module's own file. A value::Float64 is returned.
"""
function  matrixElementCached(basis::Basis, r::Int64, s::Int64, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::AsfSettings, potential::Radial.Potential, cache::InteractionStrength.XLCache,
                                grantCache::Dict{Tuple{Subshell,Subshell},Float64})
    if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    return( 0. )    end
    subshellList = basis.subshells
    opa  = SpinAngular.OneParticleOperator(0, plus)
    waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)
    opa  = SpinAngular.TwoParticleOperator(0, plus)
    waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)

    me = 0.
    for  coeff in waG1
        me = me + coeff.T * oneParticleIntegralCached(coeff.a, coeff.b, basis.orbitals, nm, grid, settings,
                                                         potential, grantCache)
    end
    for  coeff in waG2
        if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                        basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, cache)
        end
        if  typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                      basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                      settings.eeInteractionCI, cache)
        end
    end

    return( me )
end


"""
`SuperConfigurations.varianceFromCsfMatrix(conf::Configuration, orbitals::Dict{Subshell,Orbital},
                                             nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                                             meanEnergy::Float64)`
    ... to determine the (2J+1)-weighted variance of a configuration's level energies from the full CSF-pair
        Hamiltonian matrix, via Tr(H^2), WITHOUT diagonalizing it. THIS IS NOT YET THE CLOSED FORM this
        module's own docstring names as future work: it still touches every CSF PAIR (O(n^2) matrix
        elements), only the O(n^3) diagonalization is avoided. Pairs of CSFs in different J^P blocks are
        skipped outright rather than relied on to return zero, since the Hamiltonian's own J/parity
        selection rule makes them exactly zero by construction. A value::Float64 (in Hartree^2) is returned.
"""
function  varianceFromCsfMatrix(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                  grid::Radial.Grid, settings::AsfSettings, meanEnergy::Float64)
    basis, csfList = buildBasisAndCsfs(conf, orbitals)
    potential  = Nuclear.nuclearPotential(nm, grid)
    cache      = InteractionStrength.XLCache()
    grantCache = Dict{Tuple{Subshell,Subshell},Float64}()

    num2 = 0.;   den = 0.
    for  r = 1:length(csfList)
        wr  = Basics.twice(csfList[r].J) + 1
        den = den + wr
        for  s = 1:length(csfList)
            ( csfList[s].J != csfList[r].J   ||   csfList[s].parity != csfList[r].parity )   &&   continue
            Hrs  = matrixElementCached(basis, r, s, nm, grid, settings, potential, cache, grantCache)
            num2 = num2 + wr * Hrs^2
        end
    end

    return( num2/den - meanEnergy^2 )
end


"""
`SuperConfigurations.buildGroup(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                  grid::Radial.Grid, settings::AsfSettings)`
    ... to build the Gaussian statistical description of `conf`: the mean energy from
        `BasicsAZ.meanEnergyInClosedForm` (closed form, no CSF touched) and the variance from
        `SuperConfigurations.varianceFromCsfMatrix` (exact, CSF-pair-level, not yet closed form -- see this
        module's own docstring for why). A group::SuperConfigurations.Group is returned.
"""
function  buildGroup(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                       grid::Radial.Grid, settings::AsfSettings)
    meanE          = BasicsAZ.meanEnergyInClosedForm(conf, orbitals, nm, grid)
    varE           = varianceFromCsfMatrix(conf, orbitals, nm, grid, settings, meanE)
    basis, csfList = buildBasisAndCsfs(conf, orbitals)
    degeneracy     = sum( Basics.twice(csf.J) + 1   for csf in csfList )

    return( Group(conf, meanE, varE, degeneracy) )
end


"""
`SuperConfigurations.buildHamiltonianMatrix(csfList::Array{CsfR,1}, basis::Basis, nm::Nuclear.Model,
                                              grid::Radial.Grid, settings::AsfSettings)`
    ... a private helper that fills the FULL CSF-pair Hamiltonian matrix of `basis` -- zero between different
        J^P blocks, exactly as `varianceFromCsfMatrix`'s own loop skips them -- needed here because
        `buildTransitionArray` contracts this matrix with the transition operator rather than only reading
        its trace. Uses `matrixElementCached`, NOT `Hamiltonian.matrixElement` directly, so the one-particle
        radial integral is memoized by orbital pair across every CSF pair in this build (priority item 79;
        see `oneParticleIntegralCached`'s own docstring for the measurement). A matrix::Array{Float64,2} is
        returned.
"""
function  buildHamiltonianMatrix(csfList::Array{CsfR,1}, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                   settings::AsfSettings)
    potential  = Nuclear.nuclearPotential(nm, grid)
    cache      = InteractionStrength.XLCache()
    grantCache = Dict{Tuple{Subshell,Subshell},Float64}()
    n          = length(csfList)
    H          = zeros(Float64, n, n)
    for  r = 1:n
        for  s = 1:n
            ( csfList[s].J != csfList[r].J   ||   csfList[s].parity != csfList[r].parity )   &&   continue
            H[r,s] = matrixElementCached(basis, r, s, nm, grid, settings, potential, cache, grantCache)
        end
    end

    return( H )
end


"""
`SuperConfigurations.buildDipoleMatrix(csfListFinal::Array{CsfR,1}, basisFinal::Basis,
                                         csfListInitial::Array{CsfR,1}, basisInitial::Basis,
                                         subshellList::Array{Subshell,1}, grid::Radial.Grid)`
    ... a private helper that fills the CSF-pair electric-dipole (E1) REDUCED matrix element
        <csfListFinal[s] || D^(1) || csfListInitial[r]> for every pair, the two-configuration analogue of
        `buildHamiltonianMatrix` for a rank-1 operator instead of the rank-0 Hamiltonian. Reuses the SAME
        `InteractionStrength.eMultipole(1, ...)` radial integral `MultipoleMoment.emmStaticAmplitude` pairs
        with `coeff.T` (GRASP-dressed, Rule 18) -- NOT `InteractionStrength.MabEmission`, which is the
        retarded photon-emission operator and needs a transition energy this module does not yet have. The
        residual sqrt(2J_f+1) factor Rule 18 places AFTER the CSF loop is applied HERE PER ROW, using that
        row's OWN csf.J, since `csfListFinal` spans several J blocks at once rather than one level's fixed J
        -- the direct generalization of `varianceFromCsfMatrix`'s per-CSF `wr` weight to a rank-1 operator.
        `eMultipole` is MEMOIZED by orbital pair here too, the SAME fix and the SAME reason as
        `oneParticleIntegralCached` (priority item 79) -- the number of distinct orbital pairs is bounded
        by the two configurations' active shells, not by `nF x nI`. A matrix::Array{Float64,2} of size
        (length(csfListFinal), length(csfListInitial)) is returned.
"""
function  buildDipoleMatrix(csfListFinal::Array{CsfR,1}, basisFinal::Basis, csfListInitial::Array{CsfR,1},
                              basisInitial::Basis, subshellList::Array{Subshell,1}, grid::Radial.Grid)
    nF  = length(csfListFinal);   nI = length(csfListInitial)
    D   = zeros(Float64, nF, nI)
    opa = SpinAngular.OneParticleOperator(1, Basics.multipoleParity(EmMultipole(1, true)))
    dipoleCache = Dict{Tuple{Subshell,Subshell},Float64}()
    for  s = 1:nF
        wf = sqrt( Basics.twice(csfListFinal[s].J) + 1.0 )
        for  r = 1:nI
            wa = SpinAngular.computeCoefficients(opa, csfListFinal[s], csfListInitial[r], subshellList)
            me = 0.
            for  coeff in wa
                key = (coeff.a, coeff.b)
                if  haskey(dipoleCache, key)
                      tamp = dipoleCache[key]
                else  tamp = InteractionStrength.eMultipole(1, basisFinal.orbitals[coeff.a], basisInitial.orbitals[coeff.b], grid)
                      dipoleCache[key] = tamp
                end
                me = me + coeff.T * tamp
            end
            D[s,r] = me * wf
        end
    end

    return( D )
end


"""
`SuperConfigurations.buildTransitionArray(confFinal::Configuration, confInitial::Configuration,
                                            orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                            grid::Radial.Grid, settings::AsfSettings)`
    ... to build the Gaussian statistical description of the UNRESOLVED E1 transition array from
        `confInitial` to `confFinal`: the line-strength-weighted mean energy, variance and total line
        strength, EXACT at CSF-pair level, with no diagonalization of either configuration's Hamiltonian and
        no closed-form Bauche-Arnoult moment formula. The two configurations' CSFs are built on their
        COMBINED subshell list (`Basics.merge`), exactly as every two-configuration transition amplitude
        elsewhere in JAC already does, so the dipole operator can connect them.

        THE ALGEBRA, so it can be checked rather than trusted. Let H1, H2 be the ordinary CSF-pair
        Hamiltonian matrices of confInitial, confFinal (each block-diagonal in J, parity, with EXACT
        eigenvalues equal to that configuration's true level energies), and D the CSF-pair E1 REDUCED matrix
        element (Rule 18: already Wigner-Eckart-normalized, so |D_sr|^2 IS the standard line strength between
        the levels D's CSF indices expand into -- no extra statistical weight is needed anywhere below).
        Inserting the completeness relation of each configuration's OWN exact eigenvectors on either side of
        D, with NO approximation:

            S_total = Tr(D^T D)                                           -- the array's total line strength
            mean    = [Tr(H2 D D^T) - Tr(H1 D^T D)] / S_total              -- strength-weighted <E_final-E_initial>
            var     = [Tr(H2^2 D D^T) - 2 Tr(H1 D^T H2 D) + Tr(H1^2 D^T D)] / S_total - mean^2

        See `TestFrames.testModule_SuperConfigurations` for the brute-force check -- diagonalizing H1, H2
        directly and summing over the resulting level pairs -- that these traces reproduce the double sum
        rather than being a tidy-looking guess. A transitionArray::SuperConfigurations.TransitionArray is
        returned.
"""
function  buildTransitionArray(confFinal::Configuration, confInitial::Configuration,
                                 orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model, grid::Radial.Grid,
                                 settings::AsfSettings)
    Basics.extractFromConfiguration(Basics.GetParity(), confFinal) ==
        Basics.extractFromConfiguration(Basics.GetParity(), confInitial)   &&
        error("SuperConfigurations.buildTransitionArray(): confFinal and confInitial have the SAME parity; " *
              "an E1 transition array between them is identically zero.")

    relconfFinal   = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confFinal)
    relconfInitial = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confInitial)
    subshellList   = Basics.merge( Basics.generateSubshellList(relconfInitial),
                                     Basics.generateSubshellList(relconfFinal) )

    basisFinal,   csfListFinal   = buildBasisAndCsfs(confFinal,   orbitals; subshellList=subshellList)
    basisInitial, csfListInitial = buildBasisAndCsfs(confInitial, orbitals; subshellList=subshellList)

    H2 = buildHamiltonianMatrix(csfListFinal,   basisFinal,   nm, grid, settings)
    H1 = buildHamiltonianMatrix(csfListInitial, basisInitial, nm, grid, settings)
    D  = buildDipoleMatrix(csfListFinal, basisFinal, csfListInitial, basisInitial, subshellList, grid)

    DDt    = D * transpose(D)                 # (nF x nF)
    DtD    = transpose(D) * D                 # (nI x nI)
    Stotal = LinearAlgebra.tr(DtD)

    Stotal < 1.0e-14   &&
        error("SuperConfigurations.buildTransitionArray(): the total E1 line strength between $confInitial " *
              "and $confFinal is zero to machine precision -- check that they are connected by a genuine " *
              "single-electron E1 transition.")

    meanE = ( LinearAlgebra.tr(H2 * DDt) - LinearAlgebra.tr(H1 * DtD) ) / Stotal

    secondMoment = ( LinearAlgebra.tr(H2^2 * DDt)   -   2*LinearAlgebra.tr(H1 * transpose(D) * H2 * D)   +
                      LinearAlgebra.tr(H1^2 * DtD) ) / Stotal
    varE = secondMoment - meanE^2

    return( TransitionArray(confFinal, confInitial, meanE, varE, Stotal) )
end

end # module
