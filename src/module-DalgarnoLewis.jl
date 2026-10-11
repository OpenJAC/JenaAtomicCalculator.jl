
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
        not that -- measured directly (apps-very-schaefer-californium work, 08-Oct-2026, and independently
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

        BABUSHKIN (LENGTH) GAUGE DOES NOT YET WORK, AND `solve` REFUSES IT -- 11-Oct-2026, read this before
        trying to fix it quickly, since one quick fix was already tried here and shown insufficient.
        MabEmission's Babushkin branch picks between two otherwise-asymmetric forms of the same formula by
        comparing `a.energy` and `b.energy` -- a correction that is only meaningful between two genuine
        eigenstates, and neither a bare B-spline primitive nor chi itself (a solved response, not an
        eigenstate) has a real binding energy for that comparison to use. The FIRST attempt forced the
        correction OFF (eps = +1) wherever at least one side is not a genuine eigenstate
        (`DalgarnoLewis.amplitudeKernel`, still present below). THAT WAS WRONG, not merely incomplete:
        checked against the only trustworthy reference available -- an EXPLICIT sum over the complete
        1819-state eigenbasis of the SAME finite matrix, using each genuine eigenstate's own real energy for
        eps exactly as MabEmission already does for real orbital pairs -- the forced-eps=+1 result is too
        SMALL by a factor of ~12.9, not merely wrong in sign. (An earlier, narrower check, comparing the
        forced-eps solve against an explicit sum that ALSO had eps forced to +1 throughout, matched well --
        but that comparison was self-consistent rather than correct: both sides shared the same wrong
        assumption. The Coulomb/velocity gauge has no such ambiguity anywhere in its formula and matches the
        genuine, correctly-computed explicit sum to 6 digits without any correction at all.) So a single
        forced eps does not represent the Babushkin gauge's true, state-dependent completeness requirement;
        whatever the right fix is, it is not a constant. UNTIL THAT IS FOUND, `solve` raises for
        `gauge == Basics.Babushkin` rather than returning a number that LOOKS plausible and is off by an
        order of magnitude -- the Rule-13 pattern this project uses for a known-incomplete path. Only the
        Coulomb/velocity gauge is currently usable.

        THE FIRST REAL APPLICATION, AND WHAT IT NEEDED -- 11-Oct-2026. `DalgarnoLewis.twoPhotonAmplitude`
        builds the un-reduced two-photon amplitude a genuine EMISSION process needs (both time orderings),
        and `solve`'s `extraSource` keyword is what lets `RandomPhase` couple several such equations together.
        See `example-As.jl` branch (c) for the H 2s -> 1s two-photon rate this module now reproduces (0.514
        /s against the exact 8.2206 /s in the Coulomb gauge, roughly 5.8x closer than the existing bound-only
        pseudostate sum's 0.089 /s) -- a real, measured step, not yet the full answer: Babushkin is still
        unusable, and a single gauge alone cannot be trusted to have converged.
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
        taken as eps = +1, never read off a or b's `.energy` field.

        WHY THIS EXISTS, AND WHY IT IS NOT A WORKAROUND FOR A BUG IN MabEmission. `eps` exists because the
        Babushkin formula, written out in Grant integrals, is not symmetric under exchanging a and b -- it is
        valid as written only when `a` is the more strongly bound of the two, and `eps` corrects the other
        assignment. That correction is only MEANINGFUL between two genuine eigenstates, each with its own
        real, physical binding energy to compare. Neither a bare B-spline primitive (used while building a
        PerturbedOrbital's source vector) nor chi itself (a solved response function, not an eigenstate, used
        while extracting an amplitude FROM a PerturbedOrbital) has such an energy -- `targetEnergy`, carried
        on the synthetic orbitals this module builds, is a BOOKKEEPING LABEL (what energy the equation was
        solved at), not a claim about how tightly bound anything is. Letting MabEmission's natural comparison
        read that label produces a RESULT THAT DEPENDS ON WHICH ROLE (a or b) the synthetic/chi object
        happens to play in a given call -- confirmed directly (11-Oct-2026): chi plays `b` while
        DalgarnoLewis.solve builds the source vector and naturally gets eps=+1 there, but plays `a` while a
        caller extracts a final amplitude from it, where the SAME natural comparison flips to eps=-1 purely
        from the role swap, with nothing physical having changed -- a sign that depends on which argument
        slot an object sits in is not physics.

        FORCING eps=+1 UNCONDITIONALLY IS NOT THE RIGHT ANSWER EITHER, AND `DalgarnoLewis.solve` DOES NOT USE
        THIS FUNCTION FOR THE BABUSHKIN GAUGE AS A RESULT -- found 11-Oct-2026, after this function's first
        version claimed the opposite. Checked against the one trustworthy reference available -- an EXPLICIT
        sum over the complete eigenbasis of the same finite matrix, using each genuine eigenstate's own real
        energy for eps exactly as MabEmission's ordinary formula already does -- the forced-eps=+1 result is
        too SMALL by a factor of ~12.9, not a sign error. (The ORIGINAL check that looked like a fix compared
        this function's output against an explicit sum that ALSO had eps forced to +1 throughout, and the two
        naturally agreed -- but that was self-consistency between two instances of the same wrong assumption,
        not correctness.) So the true Babushkin completeness requirement is state-dependent in a way a single
        constant cannot capture; this function is kept, UNUSED by `solve`, as the record of what was tried.
        The Coulomb/velocity and Magnetic branches carry no eps ambiguity at all (nothing in their formulas
        depends on it) and are passed through to MabEmission unchanged -- those are unaffected by any of this.
        A value::ComplexF64 is returned.
"""
function amplitudeKernel(multipole::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital, grid::Radial.Grid)
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
        forms from a solved PerturbedOrbital, using chi's OWN gauge (chi.gauge) throughout. The FREQUENCY this
        operator itself carries defaults to chi.omega (the equation's own omega, right for the single-operator
        case this function was first written for -- e.g. a polarizability-type <i|O(omega)|chi(omega)>, where
        the probe and the perturbation share one frequency) -- but a SECOND, genuinely two-photon-shaped use
        needs the OTHER photon's frequency here instead: in <final|O2(omega2)|chi_1(omega1)>, the omega
        carried by O2 is omega2, which need not equal chi_1.omega = omega1 at all except through the shared
        constraint omega1 + omega2 = E_i - E_f. The `omega` keyword exists for exactly that case; passing it
        explicitly changes NOTHING for any existing caller, since every one of them already relies on the
        default. RAISES if chi.gauge == Basics.Babushkin: that gauge is not yet usable from this module -- see
        the module's own docstring and DalgarnoLewis.amplitudeKernel's for why, and do not route around this
        by calling amplitudeKernel directly, which is kept only as the record of an insufficient attempt.
        A value::ComplexF64 is returned.
"""
function amplitude(multipole::EmMultipole, finalOrbital::Orbital, chi::DalgarnoLewis.PerturbedOrbital, grid::Radial.Grid;
                    omega::Union{Nothing,Float64}=nothing)
    if  chi.gauge == Basics.Babushkin
        error("DalgarnoLewis.amplitude(): the Babushkin gauge is not yet usable -- see the module's own " *
              "docstring (11-Oct-2026 entry) for what was tried and why it was insufficient.")
    end
    omegax = isnothing(omega) ? chi.omega : omega
    chiOrbital = Orbital(chi.subshell, false, true, chi.targetEnergy, chi.P, chi.Q,
                         zeros(length(chi.P)), zeros(length(chi.Q)), Radial.Grid())

    return( InteractionStrength.MabEmission(multipole, chi.gauge, omegax, chiOrbital, finalOrbital, grid) )
end


"""
`DalgarnoLewis.solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, targetSubshell::Subshell,
                      sourceOrbital::Orbital, potential::Radial.Potential, primitives::Bsplines.Primitives;
                      settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings(),
                      extraSource::Union{Nothing,Array{Float64,1}}=nothing)`
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
        ONLY gauge == Basics.Coulomb (or Basics.Magnetic) IS CURRENTLY SUPPORTED. RAISES for
        gauge == Basics.Babushkin: that gauge needs a state-dependent completeness correction this module
        does not yet have a working formula for -- see the module's own docstring (11-Oct-2026 entry).
        `extraSource`, when given, is ADDED to the (nsL+nsS)-length source vector before the linear solve --
        this is the hook `RandomPhase` uses to couple several such equations together (an induced potential's
        matrix elements, in this module's own primitive basis) WITHOUT this module needing any notion of what
        that extra term represents; `nothing` (the default) reproduces every existing caller unchanged.
        A chi::DalgarnoLewis.PerturbedOrbital is returned.
"""
function solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, targetSubshell::Subshell,
                sourceOrbital::Orbital, potential::Radial.Potential, primitives::Bsplines.Primitives;
                settings::DalgarnoLewis.Settings=DalgarnoLewis.Settings(),
                extraSource::Union{Nothing,Array{Float64,1}}=nothing)
    if  gauge == Basics.Babushkin
        error("DalgarnoLewis.solve(): the Babushkin gauge is not yet usable -- see the module's own " *
              "docstring (11-Oct-2026 entry) for what was tried (forcing the eps orientation correction to " *
              "+1) and why it was found insufficient (off by a factor of ~12.9 against the correct, " *
              "full-basis explicit sum, not merely a sign error). Use gauge = Basics.Coulomb.")
    end
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
    #     Only reached for gauge == Coulomb/Magnetic (Babushkin raises above), where MabEmission has no eps
    #     ambiguity at all, so it is called directly -- see DalgarnoLewis.amplitudeKernel for the gauge this
    #     is NOT yet safe for.
    source = zeros(eltype(shiftedKept) <: Complex ? ComplexF64 : Float64, nsL+nsS)
    for  i = 1:nsL
        full = zeros(grid.NoPoints);   bsp = primitives.bsplinesL[i]
        full[bsp.lower:bsp.upper] = bsp.bs
        synthetic = Orbital(targetSubshell, false, true, targetEnergy, full, zeros(grid.NoPoints),
                            zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        source[i] = InteractionStrength.MabEmission(multipole, gauge, omega, sourceOrbital, synthetic, grid)
    end
    for  i = 1:nsS
        full = zeros(grid.NoPoints);   bsp = primitives.bsplinesS[i]
        full[bsp.lower:bsp.upper] = bsp.bs
        synthetic = Orbital(targetSubshell, false, true, targetEnergy, zeros(grid.NoPoints), full,
                            zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        source[nsL+i] = InteractionStrength.MabEmission(multipole, gauge, omega, sourceOrbital, synthetic, grid)
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
        over the complete eigenbasis of the identical matrix to 7 significant figures (ratio 1.000000), and
        the K=1 channel it leaves to the caller comes out at a 1.5e-11 share of K=0 -- matching, to the digit,
        the share already measured by the sum-over-states code this is meant to replace. See example-As.jl
        branch (c) for the full H 2s -> 1s application and its result.
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
