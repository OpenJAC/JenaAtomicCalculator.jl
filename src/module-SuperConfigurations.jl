
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

        COMPETITION ITEM C9, NARROWED 17-Oct-2026: `SuperConfigurations.SuperGroup` pools several `Group`'s
        into one genuine superconfiguration bin by the STANDARD combined-sample (pooled) moment formulas --
        no new CSF-level computation, and no new approximation beyond what `Group` itself already carries
        (no configuration-interaction mixing between the bin's own configurations). The pooled MEAN is
        exact regardless of mixing (a trace argument); the pooled VARIANCE is in principle only a LOWER
        BOUND, missing a CI cross term the pooled formula cannot see -- but RE-MEASURED after the
        catastrophic-cancellation fix below, that cross term came out NEGLIGIBLE (agreement to 5
        significant figures, not the 19.5 % shortfall an earlier, bug-corrupted measurement reported) on
        the one realistic strongly-mixing pair checked against full CI -- see
        `SuperConfigurations.buildSuperGroup`'s own docstring.

        CATASTROPHIC CANCELLATION, FOUND AND FIXED 17-Oct-2026, while running the f-shell-scale case C9's
        narrowing itself called for as closing evidence (Dy3+ `[Xe] 4f^9`, Z=66). Both `Group.variance`
        (via `varianceFromCsfMatrix`) and `TransitionArray.variance` (via `buildTransitionArray`) computed
        their second moment as `Tr(H^2)/den - mean^2` -- a FINAL subtraction of two numbers that, for a
        heavy ion, agree to as many as 9 of Float64's ~16 significant digits (`mean^2` alone reaches
        ~1.5e8 for Dy3+; the physical variance is order 1). The symptom was not a slightly-off number but
        a NEGATIVE variance, -0.71 Hartree^2, impossible on its face -- caught only because this module's
        own `Base.show` prints a width, and `sqrt` of a negative number is what made it visible rather than
        merely imprecise. Light-ion cases (C I at Z=6, Co2+ at Z=27) never tripped the sign, because the
        cancellation there cost only 1-5 digits, not 9 -- which means this was a LATENT defect since Stage 1
        was first built, silently eating precision in every case run before this one, just not enough to
        be caught by eye.

        THE FIX is "shift before squaring", in both functions, and it is EXACT, not approximate: subtract
        a representative constant (the diagonal mean) from the matrix BEFORE forming the square, so every
        term that is actually summed stays at the PHYSICAL energy scale and no large intermediate number is
        ever created for later cancellation to erase. For `Group`, the constant is `meanEnergy` itself and
        the shift applies to the Hamiltonian's diagonal CSF-pair elements only. For `TransitionArray`, H1
        and H2 are shifted independently by their own trace-averages before `H1^2`/`H2^2`/`H1*D'*H2*D` are
        formed; the algebra (`Tr((H+c*I)*A) = Tr(H*A) + c*Tr(A)`) shows the shift cancels exactly out of the
        variance and leaves only an ordinary additive correction in the mean. Verified three ways: (1) the
        light-ion cases this module's own tests already carry reproduce their old values exactly (diff
        0.000e+00 for `Group`/`SuperGroup`, matching to the last printed digit for `TransitionArray`); (2)
        the Dy3+ `Group.variance` that was -0.71 Hartree^2 is now +0.0266 Hartree^2 (width 0.163 a.u.),
        the physically sensible sign; (3) the pre-fix and post-fix `TransitionArray` numbers for the SAME
        Dy3+ pair agree to ~6 digits (9.5859436521e-03 against 9.5859803569e-03), consistent with the
        9-digit loss measured directly from the un-shifted traces (`Tr(H2^2*DDt)` ~ 3.7e11 collapsing to
        ~2.8e2) -- the old number was not WRONG at this particular case's scale, but one accident away from
        being so, and the margin visibly shrinks with Z and CSF count.

        SCOPE, DELIBERATELY NARROW, the same choice `AdjustedLevels` made for the same reason: ordinary
        (non-relativistic) `Configuration`s, pure Coulomb (no Breit/QED).
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
`struct  SuperConfigurations.SuperGroup`
    ... the Gaussian (two-moment) statistical description of a GENUINE superconfiguration: several
        `Group`'s, each already an unresolved configuration, POOLED into one wider statistical bin --
        competition item C9's own name for this module, completed 17-Oct-2026. The mean, variance and
        degeneracy are the STANDARD combined-sample (pooled) moments: the between-configuration spread
        is folded into the variance via the law of total variance, exactly as `SuperConfigurations.Group`
        folds the within-configuration spread into ITS OWN variance one level down. See
        `SuperConfigurations.buildSuperGroup` for the formula and for the ONE assumption it carries, named
        rather than hidden: NO configuration-interaction MIXING between the bin's own configurations is
        included, only within each one individually (the same scope `Group` itself already has). The MEAN
        is exact regardless -- a trace argument, not an approximation -- while the VARIANCE is, in
        principle, only a LOWER BOUND where configurations genuinely mix; measured on a realistic
        strongly-mixing pair against full CI, that bound came out tight (5 significant figures), not
        substantially low (see `buildSuperGroup`'s own docstring for what an earlier, bug-corrupted
        measurement claimed instead, and why it was wrong).

    + confs         ::Array{Configuration,1}   ... the configurations pooled into this bin.
    + meanEnergy    ::Float64                   ... the pooled (degeneracy-weighted) mean energy, in Hartree.
    + variance      ::Float64                   ... the pooled variance (within-configuration spread,
                                                       degeneracy-weighted, PLUS the between-configuration
                                                       spread of the individual means around the pooled
                                                       mean), in Hartree^2.
    + degeneracy    ::Float64                   ... the total (2J+1)-weighted statistical weight, SUM over
                                                       every configuration in the bin.
"""
struct  SuperGroup
    confs           ::Array{Configuration,1}
    meanEnergy      ::Float64
    variance        ::Float64
    degeneracy      ::Float64
end


"""
`Base.show(io::IO, superGroup::SuperConfigurations.SuperGroup)`
    ... prints every pooled configuration, the bin's mean energy, its width (sqrt(variance)) and its total
        degeneracy. Nothing is returned.
"""
function  Base.show(io::IO, superGroup::SuperGroup)
    println(io, "SuperConfigurations.SuperGroup:  ", join(string.(superGroup.confs), "  +  "))
    println(io, "  mean energy           = $(superGroup.meanEnergy) Hartree")
    println(io, "  width (sqrt variance) = $(sqrt(max(superGroup.variance, 0.0))) Hartree")
    println(io, "  degeneracy (SUM 2J+1) = $(superGroup.degeneracy)")
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
        Hamiltonian matrix, via Tr((H - meanEnergy*I)^2), WITHOUT diagonalizing it. THIS IS NOT YET THE
        CLOSED FORM this module's own docstring names as future work: it still touches every CSF PAIR
        (O(n^2) matrix elements), only the O(n^3) diagonalization is avoided. Pairs of CSFs in different J^P
        blocks are skipped outright rather than relied on to return zero, since the Hamiltonian's own
        J/parity selection rule makes them exactly zero by construction.

        THE MEAN IS SUBTRACTED PER ELEMENT, BEFORE SQUARING -- NOT, as an earlier version did,
        `Tr(H^2) - meanEnergy^2` AFTER the fact. Found 17-Oct-2026 on Dy3+ `[Xe] 4f^9` (meanEnergy ~
        -12160 Hartree): `Tr(H^2)/den` and `meanEnergy^2` both carry 9 integer digits, Float64 keeps about
        16 significant digits total, and their difference -- the true variance, of order 1 Hartree^2 --
        needs 7-8 of those digits just to survive the subtraction. It did not: the result came out
        NEGATIVE, -0.71 Hartree^2, a value no variance can honestly be. Subtracting `meanEnergy` from
        the DIAGONAL of H before squaring (`Tr((H-meanEnergy*I)^2)`, algebraically IDENTICAL to the old
        formula -- only `delta_rs`-weighted terms change) keeps every summand at the PHYSICAL scale of an
        energy deviation from the mean, never forming the ~1.5e8 intermediate the cancellation fed on.
        Verified byte-for-byte unchanged on every light-ion case already validated (where the old formula
        was never wrong, only ever one accident away from being so). A value::Float64 (in Hartree^2) is
        returned.
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
            Hrs  = ( r == s ? Hrs - meanEnergy : Hrs )    # shift the DIAGONAL only, before squaring
            num2 = num2 + wr * Hrs^2
        end
    end

    return( num2/den )
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
`SuperConfigurations.buildSuperGroup(confs::Array{Configuration,1}, orbitals::Dict{Subshell,Orbital},
                                       nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings)`
    ... to build the pooled Gaussian description of a GENUINE superconfiguration: build one
        `SuperConfigurations.Group` per configuration in `confs` (reusing `buildGroup` exactly, no new
        CSF-level computation), then combine them with the STANDARD combined-sample (pooled) moment
        formulas -- the same ones the Bauche-Arnoult/STA literature itself uses at this level of theory:

            G  =  SUM_i g_i                                            -- total degeneracy
            M  =  SUM_i g_i * mean_i  /  G                             -- pooled mean
            V  =  SUM_i g_i * ( var_i + (mean_i - M)^2 )  /  G         -- pooled variance

        THE MEAN `M` IS EXACT, PROVABLY, NOT JUST USUALLY. Configuration-interaction mixing between
        `confs[i]` and `confs[j]` would add OFF-DIAGONAL blocks to the full combined Hamiltonian, and an
        off-diagonal block never contributes to a trace -- so Tr(H_full) = Tr(H_AA) + Tr(H_BB) + ... holds
        with NO approximation, whatever the mixing, and `M` built from the separate `Group`s' own means is
        exactly the mean a full multi-configuration diagonalization would give.

        THE VARIANCE `V` IS NOT EXACT WHERE CONFIGURATIONS GENUINELY MIX, IN PRINCIPLE: Tr(H_full^2) picks
        up an off-diagonal cross term the pooled formula has no way to see, and that term is exactly the
        mixing `SuperGroup` omits -- `V` is therefore a LOWER BOUND on the true variance, never an
        overestimate. HOW LARGE THAT CROSS TERM ACTUALLY IS was measured, not guessed, on a realistic
        strongly-mixing pair, Co2+ `[Ar] 3d^7` + `[Ar] 3d^6 4s` (near-degenerate nd^n/nd^(n-1)(n+1)s mixing,
        a textbook case), against `Hamiltonian.performCI` on the two configurations TOGETHER (the SAME
        validation route `TestFrames.testModule_SuperConfigurations` already uses for one configuration).

        RE-MEASURED 17-Oct-2026 after the catastrophic-cancellation fix to `varianceFromCsfMatrix` (see this
        module's own top docstring): the FIRST measurement, taken through the unfixed formula, reported the
        pooled width 19.5 % too narrow (0.115 Hartree against 0.143 Hartree) -- that number was itself an
        artifact of the bug, not a property of CI mixing, because each `Group`'s own variance feeding the
        pool was already corrupted before pooling began. Through the FIXED formula, the pooled width is
        0.142711 Hartree against a true (full-CI) 0.142714 Hartree -- agreement to 5 significant figures,
        not a 19.5 % shortfall. For THIS pair, in other words, the omitted cross term turns out to be
        negligible despite the strong level mixing, because Tr(H_AB*H_BA) can be small even when individual
        off-diagonal matrix elements are not. The lower-bound property above still holds in general and the
        cross term is not zero by any argument -- only ever measured small here -- so a new pair should
        still be checked rather than assumed, but do not re-use the retracted 19.5 % figure for anything.
        A superGroup::SuperConfigurations.SuperGroup is returned.
"""
function  buildSuperGroup(confs::Array{Configuration,1}, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                            grid::Radial.Grid, settings::AsfSettings)
    isempty(confs)   &&   error("SuperConfigurations.buildSuperGroup(): confs is empty.")
    groups = [ buildGroup(conf, orbitals, nm, grid, settings)   for conf in confs ]

    G = sum( g.degeneracy   for g in groups )
    M = sum( g.degeneracy * g.meanEnergy   for g in groups ) / G
    V = sum( g.degeneracy * ( g.variance + (g.meanEnergy - M)^2 )   for g in groups ) / G

    return( SuperGroup(confs, M, V, G) )
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
        rather than being a tidy-looking guess.

        H1 AND H2 ARE SHIFTED BY A REPRESENTATIVE CONSTANT BEFORE EITHER IS SQUARED, for the same reason as
        `varianceFromCsfMatrix`: found 17-Oct-2026 on Dy3+ `[Xe] 4f^9` -> `[Xe] 4f^8 5d` (absolute energies
        ~ -12160 Hartree) that `Tr(H2^2 D D^T)` and its two partners each reach ~3.7e11 while the combination
        that survives is of order 1 -- 9 of Float64's ~16 significant digits lost in that one subtraction,
        leaving only 5-6 reliable digits in the printed variance, measured by comparing against the shifted
        form below (agreement to 3.7e-8 absolute, i.e. about 6 digits, exactly as the lost-digit count
        predicts). Writing `a = Tr(H1)/n1`, `b = Tr(H2)/n2` and `H1' = H1-a*I`, `H2' = H2-b*I`, the identity
        `Tr((H+c*I)*A) = Tr(H*A) + c*Tr(A)` propagates through every trace above and the shift cancels
        EXACTLY out of `varE` (algebra kept in the module's own validation script, not reproduced here);
        `meanE` picks up the leftover constant `(b-a)`, which is itself an ordinary, unproblematic number.
        The shifted matrices carry only the physically small energies, so no large intermediate ever forms.
        A transitionArray::SuperConfigurations.TransitionArray is returned.
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

    a  = LinearAlgebra.tr(H1) / size(H1,1)       # representative energy of confInitial, confFinal -- any
    b  = LinearAlgebra.tr(H2) / size(H2,1)       # value works exactly (see docstring); the mean is the natural one
    H1s = H1 - a*LinearAlgebra.I
    H2s = H2 - b*LinearAlgebra.I

    meanEs = ( LinearAlgebra.tr(H2s * DDt) - LinearAlgebra.tr(H1s * DtD) ) / Stotal
    meanE  = (b - a) + meanEs

    secondMomentShifted = ( LinearAlgebra.tr(H2s^2 * DDt)   -   2*LinearAlgebra.tr(H1s * transpose(D) * H2s * D)   +
                             LinearAlgebra.tr(H1s^2 * DtD) ) / Stotal
    varE = secondMomentShifted - meanEs^2

    return( TransitionArray(confFinal, confInitial, meanE, varE, Stotal) )
end

end # module
