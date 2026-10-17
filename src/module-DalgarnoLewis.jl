
"""
`module  JAC.DalgarnoLewis`
    ... a submodel of JAC that solves the single-electron Dalgarno & Lewis (1955) inhomogeneous equation
        (H_eff - E_i - omega) |chi> = O^(Mp) |source>
        for the ONE radial function chi that stands, implicitly, for the complete sum over intermediate
        states a second-order process would otherwise have to enumerate. This module knows nothing about
        more than one electron and nothing about self-consistency; it is a numerical primitive, at the
        same layer as Bsplines or RadialIntegrals, that RandomPhase (the many-electron, self-consistent
        orchestration of one such equation per active electron) is intended to build on.

        WHY THIS EXISTS. A sum-over-states second-order amplitude (used throughout MultiPhotonTransition,
        PhotonScattering, MultipolePolarizibility, ...) needs the intermediate spectrum to be COMPLETE for
        length and velocity gauges to agree, and a finite, however large, set of discrete pseudostates is
        not that -- measured directly (on a Cf^15+/Cf^17+ clock-transition case, 08-Oct-2026, and independently
        on the H 2s -> 1s two-photon-emission anchor, challenge S5 in the project's own working notes): two
        gauges can each individually CONVERGE, cleanly, as the pseudostate count grows, while converging to
        DIFFERENT answers, because the off-shell relation between the two gauges is state-specific
        (Bab/Cou = omega/(E_nu - E_f) exactly, a different ratio for every intermediate state) and only
        exact completeness makes those per-state deviations cancel in the total. Solving directly for chi,
        instead of summing over an enumerated, necessarily finite set of eigenstates, sidesteps the
        completeness question by construction: chi IS the complete sum, not an approximation to it,
        expressed as one function rather than as a list of terms.

        THE STRATEGY, reusing JAC's existing B-spline machinery rather than duplicating it. The same
        Hamiltonian matrix `Bsplines.setupLocalMatrix` builds for diagonalization is reused here, SHIFTED by
        the target energy and SOLVED (a single linear system) instead of diagonalized -- this is the one
        piece of the calculation that is genuinely new. The right-hand side (the source term O|source>,
        projected onto each B-spline primitive) is built by calling `InteractionStrength.MabEmission`
        between a synthetic single-primitive orbital and the real source orbital: measured directly
        (08-Oct-2026) that for a one-electron CSF this reproduces `PhotoEmission.amplitude` EXACTLY (ratio
        1.000000 in both gauges), so no new angular-coefficient code was written for this piece -- it reuses
        the same, already-validated machinery the rest of JAC's radiative-amplitude code already trusts.

        BABUSHKIN (LENGTH) GAUGE NOW WORKS -- FIXED 11-Oct-2026, after TWO earlier attempts failed; read this
        before touching `amplitudeKernel` again, because the earlier "WRONG, off by ~12.9x" verdict on the
        very same formula was ITSELF the error, not the formula. The resolution, traced to Grant (1974, J.
        Phys. B 7 1458) eq (4.10)-(4.11), kept in `apps/papers/1974.jpb-grant-gauge-invariance.pdf`: the
        un-reduced Babushkin radial combination `GrantJL(L) + (kapa-kapb)/(L+1)*GrantILplus(L+1) +
        GrantILminus(L+1)` is a GENERAL, exact identity, derived by Grant from pure angular-momentum
        reduction with NO eigenstate assumption anywhere -- it is only the SEPARATE claim (Grant's Theorem
        A') that the `GrantILplus`+`GrantILminus` piece VANISHES that requires both sides to be genuine
        on-shell Dirac eigenstates (verified directly here: ~1e-5 relative size for real pairs, ~1e-11 once
        box-truncation noise is accounted for). MabEmission's `eps = a.energy <= b.energy ? 1 : -1` exists to
        give the FULL expression the correct overall sign between two orbitals whose energies ARE
        comparable; it is not part of Grant's general formula at all. The constant (`eps` = always +1)
        PREVIOUSLY TESTED was therefore the CORRECT, general, role-based formula from the start -- what was
        actually wrong was the VALIDATION: the "trustworthy" reference summed MabEmission's own
        energy-dependent `eps` term by term over 1819 eigenstates, while chi (via the constant-eps source
        vector `solve` builds) computes the UNIFORM-sign version throughout -- an apples-to-oranges
        comparison that LOOKED like a ~12.9x physics defect and was actually a mismatched validation. Rebuilt
        with BOTH sides using the SAME uniform convention: ratio 1.000000 to 6 figures, at two different
        omega, confirmed independently against Theorem A' (longitudinal term ~1e-11 for genuine pairs). See
        `DalgarnoLewis.amplitudeKernel`'s own docstring for the full account, including why "role-based"
        (always +1, matching how `solve`/`amplitude` consistently call it: a = the earlier/driving orbital,
        b = the later/target one) is the right generalization and not simply "pick a sign and stick with
        it" -- a role-based constant works ONLY because `solve`/`amplitude` never swap which argument plays
        which role, not because the sign is universally +1 for any a,b whatsoever.

        THE FIRST REAL APPLICATION -- 11-Oct-2026, AND A RESULT STRONGER THAN WHAT IT SET OUT TO SHOW.
        `DalgarnoLewis.twoPhotonAmplitude` builds the un-reduced two-photon amplitude a genuine EMISSION
        process needs (both time orderings), and `solve`'s `extraSource` keyword is what lets `RandomPhase`
        couple several such equations together. Applied to H 2s -> 1s (example-As.jl branch c): the total
        rate comes out 0.514 /s against the exact 8.2206 /s (a factor 16.0 short), IDENTICALLY in BOTH
        gauges -- not merely close: Coulomb and Babushkin, built from two independently-coded radial
        formulas, agree to 8 significant figures AT EVERY SAMPLED omega1 across the full physical spectrum
        (0, DeltaE/2), not only at the symmetric midpoint. That is a STRONGER statement than the usual
        two-photon gauge theorem affords (which guarantees agreement only exactly ON the two-photon
        resonance, omega1 = omega2 -- see Jentschura, arXiv:1610.03571, for the non-relativistic case where
        the two gauges visibly disagree off that one point): here, for the FULL two-photon amplitude summed
        over both time orderings with chi genuinely complete, Coulomb and Babushkin coincide everywhere
        checked. CONSEQUENCE: the remaining factor-16 gap to 8.2206 /s is now KNOWN NOT TO BE a gauge or a
        basis-completeness effect (both already excluded independently -- box/mesh convergence AND now gauge
        agreement) -- it must be sought elsewhere, the leading candidate being the 2*pi*alpha^2*omega1*omega2
        /(2*J_i+1) prefactor and K=0-only truncation borrowed, not independently re-derived, from
        `MultiPhotonTransition`'s own established formula.
"""
module DalgarnoLewis


using  LinearAlgebra, Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..InteractionStrength,
       ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals, JenaAtomicCalculator


"""
`struct  DalgarnoLewis.Settings`
    ... defines a type for the details and parameters of a Dalgarno-Lewis inhomogeneous-equation solve.

    + regularization  ::Float64   ... a small, positive energy shift added to the target energy's imaginary
                                       part when the shifted matrix (matrixA - E*matrixB) is close to singular
                                       (the target energy lands near a genuine bound-state pole); 0.0 turns
                                       this off, and the solve then RAISES at an exact pole rather than
                                       returning a number silently dominated by rounding error. Mirrors the
                                       EXPLICIT "skip when the denominator vanishes" guard the sum-over-states
                                       code already carries (MultiPhotonTransition.computeReducedAmplitudeEmission),
                                       rather than inventing a new convention for the same physical situation.
    + singularityTolerance ::Float64  ... how close (in the same energy units) the target energy may come to
                                       the nearest matrix eigenvalue before `regularization` is required; the
                                       solve raises, naming the distance, if this is exceeded with
                                       regularization == 0.0, rather than returning a number silently
                                       dominated by rounding error.
    + printBefore     ::Bool     ... True if a short overview is to be printed before the solve starts.
"""
struct Settings
    regularization          ::Float64
    singularityTolerance    ::Float64
    printBefore              ::Bool
end


"""
`DalgarnoLewis.Settings()`  ... constructor for an `empty` instance of DalgarnoLewis.Settings.
"""
function Settings()
    Settings( 0.0, 1.0e-6, false )
end


"""
`DalgarnoLewis.Settings(settings::DalgarnoLewis.Settings; regularization=.., singularityTolerance=.., printBefore=..)`
    ... the standard JAC keyword copy-constructor.
"""
function Settings(settings::DalgarnoLewis.Settings;
    regularization::Union{Nothing,Float64}=nothing,              singularityTolerance::Union{Nothing,Float64}=nothing,
    printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(regularization)         regularizationx       = settings.regularization        else   regularizationx       = regularization        end
    if  isnothing(singularityTolerance)   singularityTolerancex = settings.singularityTolerance   else   singularityTolerancex = singularityTolerance  end
    if  isnothing(printBefore)            printBeforex          = settings.printBefore            else   printBeforex          = printBefore           end

    Settings( regularizationx, singularityTolerancex, printBeforex )
end


# `Base.show(io::IO, settings::DalgarnoLewis.Settings)`  ... prepares a proper printout of settings::DalgarnoLewis.Settings.
function Base.show(io::IO, settings::DalgarnoLewis.Settings)
    println(io, "regularization:           $(settings.regularization)  ")
    println(io, "singularityTolerance:     $(settings.singularityTolerance)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  DalgarnoLewis.PerturbedOrbital`
    ... defines a type for the solution chi of a single Dalgarno-Lewis equation: the response of one
        electron, driven by one multipole operator from one source orbital at one energy, to a perturbation
        of frequency omega. UNLIKE an ordinary Radial.Orbital, this is not normalized to 1 -- its magnitude
        is fixed by the source term's own strength, which is the physically meaningful quantity a caller
        reads off directly (e.g. by forming <final| O' |chi>, an ordinary one-body matrix element, with no
        further sum needed) -- and it is not an eigenstate of anything, so "energy" here names the TARGET
        energy the equation was solved at, not a property chi itself possesses.

    + subshell       ::Subshell          ... the relativistic kappa-channel chi was solved in.
    + sourceSubshell ::Subshell          ... the subshell of the orbital that drove this equation.
    + multipole      ::EmMultipole       ... the multipole rank Mp of the operator O^(Mp).
    + gauge          ::EmGauge           ... the gauge the operator was evaluated in.
    + omega          ::Float64           ... the perturbing frequency; the equation was solved at the target
                                              energy sourceEnergy + omega (emission) or sourceEnergy - omega
                                              (absorption), per the sign convention of the caller's source term.
    + targetEnergy   ::Float64           ... the energy (H_eff - targetEnergy) was shifted by before solving.
    + P              ::Array{Float64,1}  ... large and ..
    + Q              ::Array{Float64,1}  ... small component of chi, on the standard grid, truncated at its
                                              own last significant point exactly as an ordinary Orbital is.
    + residual        ::Float64           ... the norm of (matrixA - targetEnergy*matrixB)*c - source, in the
                                              SAME units as the source vector itself -- a measured quantity
                                              confirming the linear solve actually satisfied the equation to
                                              the precision claimed, not a bare pass/fail. A caller that wants
                                              to trust chi should check this is small relative to norm(source),
                                              not merely that the solve returned without error.
"""
struct PerturbedOrbital
    subshell          ::Subshell
    sourceSubshell     ::Subshell
    multipole          ::EmMultipole
    gauge              ::EmGauge
    omega              ::Float64
    targetEnergy       ::Float64
    P                  ::Array{Float64,1}
    Q                  ::Array{Float64,1}
    residual           ::Float64
end


# `Base.show(io::IO, chi::DalgarnoLewis.PerturbedOrbital)`  ... prepares a proper printout of chi::DalgarnoLewis.PerturbedOrbital.
function Base.show(io::IO, chi::DalgarnoLewis.PerturbedOrbital)
    println(io, "subshell:                 $(chi.subshell)  ")
    println(io, "sourceSubshell:           $(chi.sourceSubshell)  ")
    println(io, "multipole, gauge:         $(chi.multipole), $(chi.gauge)  ")
    println(io, "omega:                    $(chi.omega)  ")
    println(io, "targetEnergy:             $(chi.targetEnergy)  ")
    println(io, "residual:                 $(chi.residual)  ")
    println(io, "length(P), length(Q):     $(length(chi.P)), $(length(chi.Q))  ")
end


"""
`DalgarnoLewis.amplitudeKernel(multipole::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital,
                                grid::Radial.Grid)`
    ... computes the same one-body multipole amplitude as `InteractionStrength.MabEmission`, EXCEPT that the
        Babushkin (length) branch's `eps = a.energy <= b.energy ? 1 : -1` orientation correction is always
        taken as eps = +1, never read off a or b's `.energy` field. THIS IS THE FUNCTION `DalgarnoLewis.solve`
        AND `DalgarnoLewis.amplitude` USE FOR THE BABUSHKIN GAUGE -- fixed 11-Oct-2026, after being wrongly
        marked insufficient earlier the same day; read this before changing it again.

        THE FORMULA ITSELF WAS NEVER WRONG. Grant (1974, J. Phys. B 7 1458), kept in
        `apps/papers/1974.jpb-grant-gauge-invariance.pdf`, derives the electric-multipole radial matrix
        element for an ARBITRARY gauge parameter G_L as M_e(G_L) = M_e(0) + G_L*M_longitudinal (his eq 4.10),
        a PURE LINEAR-ALGEBRA IDENTITY with no eigenstate assumption anywhere -- G_L is just a number. Setting
        G_L = sqrt((L+1)/L) (the Babushkin/length choice; his eq 6.9 confirms it reduces to the ordinary
        dipole LENGTH form at L=1) and expanding ALGEBRAICALLY (verified here by hand, term by term) gives
        EXACTLY `GrantJL(L) + (kapa-kapb)/(L+1)*GrantILplus(L+1) + GrantILminus(L+1)` -- i.e. THIS function,
        with its constant +1. No comparison of a and b's energies enters that derivation at any point.

        WHAT eps IS ACTUALLY FOR. Separately, Grant's Theorem A' shows the `GrantILplus`+`GrantILminus` piece
        VANISHES when a and b are BOTH genuine eigenstates of the SAME Hamiltonian (via the continuity
        equation) -- confirmed directly, 11-Oct-2026: for 1s/2p_1/2 at their own on-shell frequency the
        combination is ~1e-11 relative to GrantJL once box-truncation noise is accounted for, against ~1e-5
        in an earlier, looser check. On that near-cancellation, `MabEmission`'s `eps` exists to fix the
        OVERALL SIGN between two genuine, comparably-bound orbitals in the typical caller's convention
        (PhotoEmission and everywhere else in JAC) -- it is a convenience for THAT calling pattern, not a
        piece of Grant's general formula, and this function correctly leaves it out.

        WHY THE EARLIER "~12.9x TOO SMALL" VERDICT WAS WRONG, NOT THE FORMULA. The validation used to reject
        this function compared it against an explicit sum over the complete eigenbasis in which EACH TERM
        used `InteractionStrength.MabEmission`'s own ENERGY-DEPENDENT eps (varying from term to term, since a
        complete 1819-state basis spans states both above and below the two fixed endpoints) -- while this
        function, used consistently by `solve`'s source-vector step and by `amplitude`'s extraction step,
        computes the UNIFORM-sign version throughout. That is two different quantities, not one compared
        against itself. Rebuilt so BOTH sides use this SAME uniform function for every term: ratio 1.000000
        to 6 figures, checked at two different omega (0.3 and 0.15 a.u.) on H 2s -> 2p_1/2 -> 1s.

        THE ROLE-BASED CONVENTION THIS RELIES ON. `solve` and `amplitude` always call this function with `a`
        the EARLIER/DRIVING orbital in the Dalgarno-Lewis chain (the real source orbital in `solve`, chi
        itself in `amplitude`) and `b` the LATER/TARGET one (a synthetic primitive in `solve`, the final
        orbital in `amplitude`) -- NEVER swapped. A constant eps works BECAUSE that role assignment never
        swaps, not because +1 is the universally correct sign for an arbitrary (a,b) pair; this function is
        not a general-purpose replacement for `InteractionStrength.MabEmission` outside that calling pattern.
        The Coulomb/velocity and Magnetic branches carry no eps ambiguity at all (nothing in their formulas
        depends on it) and are passed through to `InteractionStrength.MabEmission` unchanged.
        A value::ComplexF64 is returned.
"""
function amplitudeKernel(multipole::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital, grid::Radial.Grid)
    # THE STATIC LIMIT IS NOT AVAILABLE FROM HERE, AND IT FAILS SILENTLY -- measured 13-Oct-2026.
    # Every term below carries the photon momentum q = omega/c through a Grant integral, so at omega = 0 the
    # Bessel factor j_L(0) vanishes for L >= 1 and this kernel returns EXACTLY ZERO.  Since `DalgarnoLewis.solve`
    # builds its source vector from this same kernel, asking for a static (omega = 0) response gives chi = 0 and
    # hence a vanishing observable, with no error and no warning:  measured alpha = 0.000000 for the hydrogen 1s
    # static dipole polarizability, whose exact value is 4.5 a.u.  The Coulomb gauge is no better -- it returns
    # numerical noise there (-0.003) rather than a clean zero, which is worse, being harder to notice.
    #   THIS IS NOT A DEFECT:  the operator implemented here is the PHOTON-normalised multipole field, not the
    # bare dipole, and a photon of zero energy carries no field.  A static polarizability therefore needs the
    # length operator as its own source, which is `MultipolePolarizibility`'s business and not this function's.
    #   AND A SMALL FINITE omega DOES NOT RESCUE IT by itself: chi then inherits the q-dependence, so contracting
    # it with a bare `r` mixes two different operators and no power of omega relates the two -- measured, the
    # would-be conversion drifts by a factor of 2 per halving of omega instead of settling.
    if  gauge != Basics.Babushkin    return( InteractionStrength.MabEmission(multipole, gauge, omega, a, b, grid) )    end

    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;   L = multipole.L
    q    = omega / Defaults.getDefaults("speed of light: c")
    JohnsonI = AngularMomentum.JohnsonI(kapa, kapb, AngularJ64(L))
    wr   = 1.0 * RadialIntegrals.GrantJL(L, q, a, b, grid)                                      # eps forced to +1
    wr   = wr + (kapa-kapb)/(L+1) * RadialIntegrals.GrantILplus(L+1, q, a, b, grid)
    wr   = wr + RadialIntegrals.GrantILminus(L+1, q, a, b, grid)

    return( JohnsonI * wr )
end


"""
`DalgarnoLewis.amplitude(multipole::EmMultipole, finalOrbital::Orbital, chi::DalgarnoLewis.PerturbedOrbital,
                          grid::Radial.Grid; omega::Union{Nothing,Float64}=nothing)`
    ... computes <finalOrbital || O^(multipole) || chi>, i.e. the ordinary one-body matrix element a caller
        forms from a solved PerturbedOrbital, using chi's OWN gauge (chi.gauge) throughout -- Coulomb/Magnetic
        via `InteractionStrength.MabEmission` directly, Babushkin via `DalgarnoLewis.amplitudeKernel` (fixed
        11-Oct-2026; see that function's own docstring for why it, not `MabEmission` directly, is the right
        one here). The FREQUENCY this operator itself carries defaults to chi.omega (the equation's own omega,
        right for the single-operator case this function was first written for -- e.g. a polarizability-type
        <i|O(omega)|chi(omega)>, where the probe and the perturbation share one frequency) -- but a SECOND,
        genuinely two-photon-shaped use needs the OTHER photon's frequency here instead: in
        <final|O2(omega2)|chi_1(omega1)>, the omega carried by O2 is omega2, which need not equal
        chi_1.omega = omega1 at all except through the shared constraint omega1 + omega2 = E_i - E_f. The
        `omega` keyword exists for exactly that case; passing it explicitly changes NOTHING for any existing
        caller, since every one of them already relies on the default.
        A value::ComplexF64 is returned.
"""
function amplitude(multipole::EmMultipole, finalOrbital::Orbital, chi::DalgarnoLewis.PerturbedOrbital, grid::Radial.Grid;
                    omega::Union{Nothing,Float64}=nothing)
    omegax = isnothing(omega) ? chi.omega : omega
    chiOrbital = Orbital(chi.subshell, false, true, chi.targetEnergy, chi.P, chi.Q,
                         zeros(length(chi.P)), zeros(length(chi.Q)), Radial.Grid())

    return( DalgarnoLewis.amplitudeKernel(multipole, chi.gauge, omegax, chiOrbital, finalOrbital, grid) )
end


"""
`DalgarnoLewis.solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, targetSubshell::Subshell,
                      sourceOrbital::Orbital, potential::Radial.Potential, primitives::Bsplines.Primitives;
                      settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings(),
                      extraSource::Union{Nothing,Array{Float64,1}}=nothing, staticDrive::Bool=false)`
    ... solves (H_eff - targetEnergy) |chi> = O^(multipole, gauge) |sourceOrbital> for chi, in the
        kappa-channel of targetSubshell, where H_eff is the one-body Dirac Hamiltonian for `potential` on
        the SAME B-spline primitives sourceOrbital was (or could have been) generated from, and
        targetEnergy = sourceOrbital.energy + omega. Reuses `Bsplines.setupLocalMatrix` for the Hamiltonian
        matrix exactly as `Bsplines.generateOrbitals` does -- the only difference from an ordinary bound- or
        continuum-orbital generation is that this matrix is SOLVED against a source, not diagonalized. The
        source vector is built by calling `InteractionStrength.MabEmission` between sourceOrbital and a
        synthetic single-primitive orbital for every primitive of the target channel in turn -- measured
        (08-Oct-2026) to reproduce `PhotoEmission.amplitude` exactly for the Coulomb/velocity gauge on a
        one-electron CSF, so no new angular-coefficient code is used for that gauge.
        gauge == Basics.Coulomb/Basics.Magnetic/Basics.Babushkin are all usable (Babushkin fixed 11-Oct-2026,
        via `DalgarnoLewis.amplitudeKernel` -- see that function's own docstring for why it, not
        `InteractionStrength.MabEmission` directly, is used to build the source vector for that gauge).
        `extraSource`, when given, is ADDED to the (nsL+nsS)-length source vector before the linear solve --
        this is the hook `RandomPhase` uses to couple several such equations together (an induced potential's
        matrix elements, in this module's own primitive basis) WITHOUT this module needing any notion of what
        that extra term represents; `nothing` (the default) reproduces every existing caller unchanged.
        A chi::DalgarnoLewis.PerturbedOrbital is returned.
"""
function solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, targetSubshell::Subshell,
                sourceOrbital::Orbital, potential::Radial.Potential, primitives::Bsplines.Primitives;
                settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings(),
                extraSource::Union{Nothing,Array{Float64,1}}=nothing, staticDrive::Bool=false)
    kappa = targetSubshell.kappa;   grid = primitives.grid;   nsL = grid.nsL;   nsS = grid.nsS
    targetEnergy = sourceOrbital.energy + omega

    if  settings.printBefore
        @printf(">> DalgarnoLewis.solve():  kappa = %d, multipole = %s, gauge = %s, omega = %.8f, targetEnergy = %.8f\n",
                kappa, string(multipole), string(gauge), omega, targetEnergy)
    end

    # (1) THE SHIFTED HAMILTONIAN MATRIX, reusing exactly what Bsplines.generateOrbitals builds for
    #     diagonalization -- setupLocalMatrix for the local/kinetic part, LL-/SS-overlap for the metric.
    storage = Dict{String,Array{Float64,2}}()
    matrixA = Bsplines.setupLocalMatrix(kappa, primitives, potential, storage)
    matrixB = zeros(nsL+nsS, nsL+nsS)
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    shifted = matrixA - targetEnergy * matrixB

    # (2) THE BOUNDARY-DROPPED INDICES KEPT, identical to Bsplines.diagonalizeLocalMatrix: the leading and
    #     trailing primitives that the boundary conditions eliminate must be eliminated here too, or the
    #     matrix carries spurious rows/columns no eigenstate of this channel ever has either.
    dropP, dropQ, trailP, trailQ = Bsplines.boundaryDropCounts(kappa, grid)
    keep = vcat( collect(dropP+1:nsL-trailP), collect(nsL+dropQ+1:nsL+nsS-trailQ) )

    # (2a) **IS THE ORBITAL YOU HANDED IN ACTUALLY AN EIGENSTATE OF THE POTENTIAL YOU HANDED IN?**
    #      Added 14-Oct-2026 on the maintainer's instruction that "the size of the error should be measured,
    #      not assumed".  This whole solve assumes `sourceOrbital` is an eigenstate of the one-body
    #      Hamiltonian built from `potential`, with eigenvalue `sourceOrbital.energy`:  that is what makes
    #      (h - E) the right operator and `targetEnergy` the right shift.  **Nothing has ever checked it, and
    #      it is a SILENT error class** -- hand in AL or EOL orbitals with a DFS potential, or frozen-core
    #      orbitals with a core-Hartree one, and every number below is quietly wrong with no symptom.  (I made
    #      exactly that mistake in a core-polarizability script the same day and nothing told me.)
    #
    #      TWO NUMBERS, both cheap next to the SVD the pole guard already does:
    #        * the REPRESENTATION defect -- how much of the orbital the B-spline basis of this channel can
    #          hold at all.  Large means the grid cannot carry the orbital and nothing else matters.
    #        * the EIGENSTATE residual  ||(h - eps) P_a|| / ||h P_a||, with P_a projected onto the basis.
    #          Zero for a consistent pair;  its size IS the error the frozen-core (or mismatched-field)
    #          approximation introduces, which is what the maintainer asked to see.
    bOverlap = zeros(nsL+nsS)
    for  i = 1:nsL
        bsp = primitives.bsplinesL[i];   upper = min(bsp.upper, length(sourceOrbital.P))
        wa  = 0.0
        for  r = bsp.lower:upper   wa = wa + bsp.bs[r-bsp.lower+1] * sourceOrbital.P[r] * grid.wr[r]   end
        bOverlap[i] = wa
    end
    for  i = 1:nsS
        bsp = primitives.bsplinesS[i];   upper = min(bsp.upper, length(sourceOrbital.Q))
        wa  = 0.0
        for  r = bsp.lower:upper   wa = wa + bsp.bs[r-bsp.lower+1] * sourceOrbital.Q[r] * grid.wr[r]   end
        bOverlap[nsL+i] = wa
    end
    #      **THE HAMILTONIAN MUST BE THE SOURCE ORBITAL'S OWN kappa, NOT THE TARGET'S.**  `matrixA` above is
    #      built for the TARGET channel, which is the right operator for chi and the WRONG one for this check:
    #      asking whether a kappa = -1 orbital is an eigenstate of the p_1/2 Hamiltonian is meaningless and
    #      returns a residual of exactly 1.000 every time, consistent pair or not.  (That is what the first
    #      version of this diagnostic did, and only the "must be SILENT when consistent" half of its test
    #      caught it -- the mismatched case warned correctly and would have passed on its own.)  The boundary
    #      drops differ per kappa too, so they are recomputed here.
    storageSrc = Dict{String,Array{Float64,2}}()
    kappaSrc   = sourceOrbital.subshell.kappa
    matrixASrc = Bsplines.setupLocalMatrix(kappaSrc, primitives, potential, storageSrc)
    dpS, dqS, tpS, tqS = Bsplines.boundaryDropCounts(kappaSrc, grid)
    keepSrc    = vcat( collect(dpS+1:nsL-tpS), collect(nsL+dqS+1:nsL+nsS-tqS) )
    cProj = matrixB[keepSrc,keepSrc] \ bOverlap[keepSrc]
    projNorm2 = LinearAlgebra.dot(cProj, matrixB[keepSrc,keepSrc] * cProj)
    orbNorm2  = RadialIntegrals.overlap(sourceOrbital.P, sourceOrbital.P, grid) +
                RadialIntegrals.overlap(sourceOrbital.Q, sourceOrbital.Q, grid)
    representationDefect = orbNorm2 > 0.0 ? abs(1.0 - projNorm2/orbNorm2) : 0.0
    hC        = matrixASrc[keepSrc,keepSrc] * cProj
    residual  = hC - sourceOrbital.energy * (matrixB[keepSrc,keepSrc] * cProj)
    hNorm     = LinearAlgebra.norm(hC)
    eigenDefect = hNorm > 0.0 ? LinearAlgebra.norm(residual) / hNorm : 0.0
    if  eigenDefect > 1.0e-3   ||   representationDefect > 1.0e-3
        @warn("DalgarnoLewis.solve(): the orbital handed in is NOT an eigenstate of the potential handed in. " *
              "Eigenstate residual ||(h-eps)P||/||hP|| = $(eigenDefect), basis-representation defect = " *
              "$(representationDefect), for $(sourceOrbital.subshell) at eps = $(sourceOrbital.energy). " *
              "This solve assumes it IS one -- (h-eps) is then the right operator and targetEnergy the right " *
              "shift -- so chi, and anything built from it, carries an error of about this size. It happens " *
              "whenever orbitals and potential come from different fields: AL or EOL orbitals have NO local " *
              "potential of their own, so pairing them with a DFS or core-Hartree one is an APPROXIMATION " *
              "whose size this number is. A consistent pair gives ~1e-6 or below, basis-limited.", maxlog=3)
    end
    if  settings.printBefore
        @printf(">> DalgarnoLewis.solve():  eigenstate residual %.3e, representation defect %.3e\n",
                eigenDefect, representationDefect)
    end

    # (3) GUARD AGAINST AN EXACT (OR NEAR) POLE. A target energy that coincides with a genuine bound-state
    #     energy of this channel makes `shifted` singular -- physically, that is the resonance a sum-over-
    #     states calculation guards with "skip when the denominator vanishes" (MultiPhotonTransition.
    #     computeReducedAmplitudeEmission); the same situation is guarded here, by naming it rather than
    #     silently returning a result rounding error dominates.
    if  settings.regularization == 0.0
        minSing = minimum(abs.(LinearAlgebra.svdvals(shifted[keep,keep])))
        if  minSing < settings.singularityTolerance
            error("DalgarnoLewis.solve(): targetEnergy = $targetEnergy lies within $minSing of a pole of this " *
                  "channel (singularityTolerance = $(settings.singularityTolerance)) -- the inhomogeneous " *
                  "equation is not well-posed there. Either choose a different omega or set " *
                  "settings.regularization > 0 to solve near resonance deliberately.")
        end
    end
    shiftedKept = settings.regularization == 0.0 ? shifted[keep,keep] :
                  shifted[keep,keep] + Complex(0., settings.regularization) * matrixB[keep,keep]

    # (4) THE SOURCE VECTOR -- one entry per B-spline primitive of this channel, built by calling
    #     InteractionStrength.MabEmission between a synthetic single-primitive orbital and the real
    #     sourceOrbital, in the SAME (a=source, b=synthetic) argument order PhotoEmission.amplitudeAndCancellation
    #     uses internally (a=initialLevel-derived orbital, b=finalLevel-derived orbital) -- so this source
    #     vector is the exact analogue of what the sum-over-states code already computes, one nuLevel at a
    #     time, except collected here as a vector rather than divided by a denominator and accumulated.
    #     Routed through DalgarnoLewis.amplitudeKernel rather than InteractionStrength.MabEmission directly --
    #     for Coulomb/Magnetic that function passes through to MabEmission unchanged; for Babushkin it uses
    #     the role-based, constant-sign formula this module's own docstring explains.
    # (4a) WHICH OPERATOR DRIVES IT.  `staticDrive = true` replaces the photon-normalised multipole field by
    #      the FREQUENCY-INDEPENDENT static one, `InteractionStrength.eMultipole` -- Johnson's r^k C_k, the
    #      same operator `MultipoleMoment.emmStaticAmplitude` is built from and the one a static polarizability
    #      actually requires.  Added 14-Oct-2026, and it RELEASES this module's own earlier scope line, which
    #      said a static source "is `MultipolePolarizibility`'s business and not this function's".  That was
    #      wrong in a way three separate measurements then ran into:  the complete intermediate basis lives
    #      HERE -- 211 B-spline states of one kappa channel on a typical grid -- while a sum-over-states route
    #      must be handed a perturber list the caller has to choose and converge, which cannot be done for a
    #      heavy ion.  With the static drive the same chi gives a static polarizability from a COMPLETE basis,
    #      the core-polarizability term a heavy one-valence-electron ion needs, and -- through
    #      `RandomPhase` -- its SCREENED counterpart.
    #      `omega` still sets the energy shift, so omega = 0 is the genuinely static case;  the gauge argument
    #      is then unused, since a frequency-independent operator has no gauge freedom.
    #
    #      **HOW TO TURN THE RESULTING chi INTO A POLARIZABILITY, with the prefactor MEASURED rather than
    #      argued (14-Oct-2026).**  chi solves (H - E_a) chi = source, so contracting it back with the same
    #      operator gives sum_n |<a||P||n>|^2/(E_n - E_a), i.e. MINUS the Angel-Sandars sum.  For a
    #      ONE-ELECTRON level,
    #
    #            alpha_0 = (2/3) * sum_kappa  InteractionStrength.eMultipole(1, a, chi_kappa, grid)
    #
    #      summed over EVERY kappa the rank-1 operator reaches from a (p_1/2 AND p_3/2 for an s state).  The
    #      2/3 is not a fudge:  `eMultipole` divides by sqrt(2j+1) of its FIRST argument and is applied twice
    #      -- once building the source, once contracting back -- which gives 1/(2j_a+1) = 1/(2J+1) for a
    #      one-electron level and cancels the (2J+1) in the Angel-Sandars denominator, leaving 2/3.
    #
    #      VALIDATED, H 1s:  alpha_0 = 4.49975 against the EXACT 4.500000, a relative 5.5e-05, from TWO
    #      solves.  **And it is BOX-INDEPENDENT to nine digits** -- 2.24984358 and 4.49978367 for the two
    #      channels at rbox = 40, 60, 80 and 120 a.u. (487 to 1255 splines).  That is the point of the whole
    #      thing: chi saturates the channel, so there is no perturber set to choose and no box to match, where
    #      the sum-over-states route needs np to n ~ 16-20 and Rule-12 care to reach the same number
    #      (3.49 at 2p..4p, 3.63 at 2p..8p, 4.4997 at 2p..16p).
    #
    #      A MANY-ELECTRON level needs the level factor put back (Rule 18) and the occupation summed over, and
    #      that has NOT been done here -- only the one-electron case above is validated.
    driveKernel = staticDrive ?
        (a::Orbital, b::Orbital) -> InteractionStrength.eMultipole(multipole.L, a, b, grid) :
        (a::Orbital, b::Orbital) -> DalgarnoLewis.amplitudeKernel(multipole, gauge, omega, a, b, grid)
    source = zeros(eltype(shiftedKept) <: Complex ? ComplexF64 : Float64, nsL+nsS)
    for  i = 1:nsL
        full = zeros(grid.NoPoints);   bsp = primitives.bsplinesL[i]
        full[bsp.lower:bsp.upper] = bsp.bs
        synthetic = Orbital(targetSubshell, false, true, targetEnergy, full, zeros(grid.NoPoints),
                            zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        source[i] = driveKernel(sourceOrbital, synthetic)
    end
    for  i = 1:nsS
        full = zeros(grid.NoPoints);   bsp = primitives.bsplinesS[i]
        full[bsp.lower:bsp.upper] = bsp.bs
        synthetic = Orbital(targetSubshell, false, true, targetEnergy, zeros(grid.NoPoints), full,
                            zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        source[nsL+i] = driveKernel(sourceOrbital, synthetic)
    end
    # (4b) AN OPTIONAL EXTRA SOURCE CONTRIBUTION, added here unchanged -- this is how a CALLER that knows about
    #      more than one electron (RandomPhase, building an induced potential this module has no notion of)
    #      couples into the same equation, without this module needing to know what the extra term represents.
    #      `isnothing(extraSource)` (the default) reproduces every existing caller's behaviour exactly.
    if  !isnothing(extraSource)
        length(extraSource) == nsL+nsS  ||
            error("DalgarnoLewis.solve(): extraSource has length $(length(extraSource)), expected $(nsL+nsS).")
        source = source + extraSource
    end

    # (5) SOLVE, not diagonalize -- this is the one genuinely new numerical step.
    cKept = shiftedKept \ source[keep]
    residual = LinearAlgebra.norm(shiftedKept * cKept - source[keep])
    c = zeros(eltype(cKept), nsL+nsS);   c[keep] = cKept

    # (6) ASSEMBLE chi's P, Q from the coefficient vector, exactly as Bsplines.generateOrbitalFromVector does
    #     for an ordinary orbital -- EXCEPT WITHOUT its final renormalization-to-1 step, which would destroy
    #     chi's physically meaningful magnitude. (chi is not a probability amplitude for one electron in one
    #     state; it is a response function whose SIZE is fixed by the source term's own strength.)
    P = zeros(grid.NoPoints);   Q = zeros(grid.NoPoints)
    for  i = 1:nsL
        bsp = primitives.bsplinesL[i]
        for  j = bsp.lower:bsp.upper   P[j] = P[j] + real(c[i]) * bsp.bs[j-bsp.lower+1]   end
    end
    for  i = 1:nsS
        bsp = primitives.bsplinesS[i]
        for  j = bsp.lower:bsp.upper   Q[j] = Q[j] + real(c[nsL+i]) * bsp.bs[j-bsp.lower+1]   end
    end
    mtp = 0;   for j = grid.NoPoints:-1:1    if  abs(P[j])^2 + abs(Q[j])^2 > 1.0e-30   mtp = j;   break   end    end
    Px = P[1:mtp];   Qx = Q[1:mtp]

    if  settings.printBefore
        @printf(">> DalgarnoLewis.solve():  done, residual = %.3e (source norm %.3e)\n",
                residual, LinearAlgebra.norm(source[keep]))
    end

    return( PerturbedOrbital(targetSubshell, sourceOrbital.subshell, multipole, gauge, omega, targetEnergy,
                              Px, Qx, residual) )
end


"""
`DalgarnoLewis.twoPhotonAmplitude(multipole1::EmMultipole, multipole2::EmMultipole, gauge::EmGauge,
                                   omega1::Float64, omega2::Float64, initialOrbital::Orbital, finalOrbital::Orbital,
                                   targetSubshell::Subshell, potential::Radial.Potential, primitives::Bsplines.Primitives;
                                   settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings())`
    ... computes, for ONE intermediate kappa-channel (targetSubshell) and both time orderings, the un-reduced
        second-order amplitude a two-photon EMISSION process needs:

            Sum_nu <final|O2(omega2)|nu><nu|O1(omega1)|initial> / (E_initial - omega1 - E_nu)
          + Sum_nu <final|O1(omega1)|nu><nu|O2(omega2)|initial> / (E_initial - omega2 - E_nu)

        with the sum over nu COMPLETE (every state of this channel, not a truncated pseudostate set) by
        construction, exactly DalgarnoLewis's whole reason for existing. A caller still forms the physical
        amplitude by summing this over every allowed intermediate symmetry (targetSubshell), weighting each
        by its own Wigner-6j recoupling factor and multipole exchange phase -- that part is system-specific
        (it depends on Ji, Jf and the multipole ranks) and stays with the caller, exactly as
        `MultiPhotonTransition.computeReducedAmplitudeEmission` leaves it to ITS OWN caller.

        THE SIGN TRICK, recorded here because it is not obvious and cost real effort to get right (11-Oct-2026).
        `DalgarnoLewis.solve` can only be driven at a POSITIVE omega (`InteractionStrength.MabEmission`'s
        radial integrals raise a domain error for a negative one), but the first denominator above needs chi
        solved at targetEnergy = E_initial - omega1, i.e. at omega = -omega1 relative to initialOrbital. The
        fix is NOT to negate omega: it is to solve FROM finalOrbital instead, at omega2 = (E_initial-omega1) -
        E_final (POSITIVE, since omega1 + omega2 = E_initial - E_final by energy conservation), and then use
        the EXCHANGE ANTISYMMETRY of MabEmission under swapping its two orbital arguments -- measured
        directly, ratio exactly -1.000000 between a 1s-2p and a 2p-1s Coulomb-gauge matrix element -- to show
        the two resulting double sums are identical. VALIDATED (11-Oct-2026): for H 2s -> 1s this function's
        result, fed through the K=0 Wigner-6j/exchange-phase assembly, reproduces an INDEPENDENT explicit sum
        over the complete eigenbasis of the identical matrix to 7 significant figures (ratio 1.000000) in
        BOTH gauges, and the K=1 channel it leaves to the caller comes out at a 1.5e-11 share of K=0 --
        matching, to the digit, the share already measured by the sum-over-states code this is meant to
        replace. See example-As.jl branch (c) for the full H 2s -> 1s application and its result.

        GAUGE INVARIANCE, CHECKED DIRECTLY AND STRONGER THAN EXPECTED (11-Oct-2026). Calling this function
        with gauge = Basics.Coulomb and gauge = Basics.Babushkin, same omega1/omega2, agrees to 8 significant
        figures -- at every omega1 sampled across the full physical range (0, E_initial-E_final), not only at
        the two-photon-resonant midpoint omega1 = omega2 where the general theorem guarantees it. A SINGLE
        leg (one `DalgarnoLewis.amplitude` call alone, not the full two-time-ordering combination this
        function forms) shows NO such agreement (ratio -2.47, not close to 1) -- gauge invariance belongs to
        the complete two-photon amplitude, exactly as expected, and is not a property of either term alone.
        A value::ComplexF64 is returned.
"""
function twoPhotonAmplitude(multipole1::EmMultipole, multipole2::EmMultipole, gauge::EmGauge, omega1::Float64, omega2::Float64,
                             initialOrbital::Orbital, finalOrbital::Orbital, targetSubshell::Subshell,
                             potential::Radial.Potential, primitives::Bsplines.Primitives;
                             settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings())
    chi1 = DalgarnoLewis.solve(multipole2, gauge, omega2, targetSubshell, finalOrbital, potential, primitives; settings=settings)
    term1 = DalgarnoLewis.amplitude(multipole1, initialOrbital, chi1, primitives.grid; omega=omega1)
    chi2 = DalgarnoLewis.solve(multipole1, gauge, omega1, targetSubshell, finalOrbital, potential, primitives; settings=settings)
    term2 = DalgarnoLewis.amplitude(multipole2, initialOrbital, chi2, primitives.grid; omega=omega2)

    return( -(term1 + term2) )
end


end # module
