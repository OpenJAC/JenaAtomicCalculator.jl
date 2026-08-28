
"""
`module  JAC.WeakInteractionMoment`
... a submodel of JAC that computes the P-odd and P,T-odd one-electron interaction amplitudes between two bound-state levels: the nuclear
    WEAK CHARGE, the nuclear ANAPOLE MOMENT and the nuclear SCHIFF MOMENT.  These are the bare matrix elements only; every observable built
    from them -- the parity-non-conserving E1 amplitude, an EDM enhancement factor -- needs a sum over intermediate states in addition and
    belongs to `WeakInteractionEnhancement`.

    This module is the P-odd counterpart of `MultipoleMoment` and has the same shape: no `Settings`, called directly, returning reduced
    matrix elements <alpha_f J_f || H || alpha_i J_i>.

    THE THREE OPERATORS, and why each takes the radial form it does.  The Dirac matrices gamma_5 and alpha both MIX the large and small
    components, so the weak-charge and anapole integrands are of the P*Q type; the Schiff operator is a function of position only and so
    keeps the P*P + Q*Q type of an ordinary electric multipole:

        weak charge   H_W = (G_F/2sqrt2) Q_W rho(r) gamma_5                  rank 0, P-odd
                      one-electron ME  =  i (G_F/2sqrt2) Q_W  INT rho(r) [P_a Q_b - Q_a P_b] dr,   kappa_b = -kappa_a

        anapole       H_A = (G_F/sqrt2) kappa_anapole alpha.I rho(r)         rank 1, P-odd
                      one-electron ME  =  i [ <Om(k_a)||sigma||Om(-k_b)> INT rho P_a Q_b dr
                                             - <Om(-k_a)||sigma||Om(k_b)> INT rho Q_a P_b dr ]

        Schiff        H_S = 4 pi S . grad rho(r)                             rank 1, P-odd AND T-odd
                      electric-multipole structure: <kappa_a||C^1||kappa_b> INT (d rho/dr) [P_a P_b + Q_a Q_b] dr

    The weak-charge amplitude is PURELY IMAGINARY between real orbitals -- that is not a convention but a consequence of gamma_5 being
    antisymmetric in the large/small components, and it is one of the things the example branch checks.

    Because `rho(r)` is normalized to INT rho 4 pi r^2 dr = 1, the nuclear charge (or weak charge) sits entirely in the prefactor and never
    in the density, which keeps the two separable.
"""
module WeakInteractionMoment


using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..SpinAngular


"""
`WeakInteractionMoment.GF`  ... the Fermi coupling constant in atomic units, G_F = 2.22254e-14 a.u.
"""
const GF          = 2.22254e-14


"""
`WeakInteractionMoment.sinThetaW2`  ... the weak mixing angle, sin^2(theta_W) = 0.23122 (PDG).
"""
const sinThetaW2  = 0.23122


"""
`WeakInteractionMoment.anapoleAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid;`
                                        `kappaAnapole::Float64=1.0, display::Bool=false)`
    ... to compute the anapole-moment amplitude <alpha_f J_f || H^(anapole) || alpha_i J_i> for the given final and initial level and for
        the nuclear density of the given nuclear model.  The operator is `(G_F/sqrt2) kappa_anapole alpha.I rho(r)`, of rank 1 and odd
        parity.  A value::ComplexF64 is returned, purely imaginary for real radial orbitals.

        THE ANGULAR STRUCTURE IS NOT THE MAGNETIC-MULTIPOLE ONE, and a first implementation that assumed it was had to be withdrawn: the
        template <-kappa_a||C^1||kappa_b> (kappa_a + kappa_b) is parity-EVEN, so it survived for p_1/2 <- p_3/2 and vanished for BOTH
        s_1/2 <- p_1/2 and s_1/2 <- p_3/2 -- exactly backwards for a P-odd operator, and it forbade precisely the s_1/2 <-> p_1/2 mixing
        through which an anapole moment reaches a measured PNC transition.  Branch a of `examples/example-Bb.jl` caught it as a column of
        sixteen identical zeros.

        What is correct follows from alpha = [[0, sigma], [sigma, 0]], which connects the LARGE component of one orbital with the SMALL
        component of the other.  With `Omega_(-kappa) = -(sigma.rhat) Omega_(kappa)` the two cross terms combine, and

            <kappa_a || alpha rho || kappa_b>  =  i [ <Omega_(kappa_a)||sigma||Omega_(-kappa_b)>  INT rho P_a Q_b dr
                                                    - <Omega_(-kappa_a)||sigma||Omega_(kappa_b)>  INT rho Q_a P_b dr ] .

        Since sigma acts only on spin it cannot change the orbital angular momentum, so a term vanishes unless the two kappa it couples
        share the same l -- which is what makes the operator parity-odd, Omega_(-kappa) carrying orbital parity l+-1.

        THE TWO CROSS TERMS MUST BE KEPT APART.  A first version combined them into the symmetric P_a Q_b + Q_a P_b, assuming their angular
        factors differ only by a sign.  That is true for gamma_5, whose angular factors are both simply delta(kappa_a, -kappa_b), and FALSE
        here: for s_1/2 <-> p_1/2 the two are <s||sigma||s> = sqrt(6) and <p_1/2||sigma||p_1/2> = -sqrt(6)/3, a ratio of -3 and not -1.  The
        error was invisible in the selection rules and in the reality of the amplitude, and showed up only in the HERMITICITY check, where
        exchanging the two levels multiplied the amplitude by -3 instead of -1.

        `kappaAnapole` is the dimensionless nuclear anapole constant, left at 1.0 so that the amplitude may be scaled by whatever value the
        user adopts.
"""
function anapoleAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid;
                          kappaAnapole::Float64=1.0, display::Bool=false)
    # rank 1 and odd parity: the triangular rule and a parity change are both required, and both zeros are still displayed
    if      finalLevel.parity == initialLevel.parity                                 amplitude = ComplexF64(0.)
    elseif  !AngularMomentum.isTriangle(finalLevel.J, AngularJ64(1), initialLevel.J)  amplitude = ComplexF64(0.)
    else
        rho     = WeakInteractionMoment.nuclearDensity(nm, grid)
        prefac  = WeakInteractionMoment.GF / sqrt(2.0) * kappaAnapole
        kernel  = (orba, orbb) -> begin
            kapa = orba.subshell.kappa;    kapb = orbb.subshell.kappa
            sgA  = WeakInteractionMoment.sigmaReducedMe( kapa, -kapb)    # weights P_a Q_b
            sgB  = WeakInteractionMoment.sigmaReducedMe(-kapa,  kapb)    # weights Q_a P_b
            if  sgA == 0.  &&  sgB == 0.    return( ComplexF64(0.) )    end
            im / sqrt(Basics.subshell_2j(orba.subshell) + 1.0) *
                ( sgA * WeakInteractionMoment.radialIntegralPQ(rho, orba, orbb, grid) -
                  sgB * WeakInteractionMoment.radialIntegralPQ(rho, orbb, orba, grid) )
        end
        amplitude = prefac * WeakInteractionMoment.oneParticleAmplitude(1, kernel, finalLevel, initialLevel)
    end
    #
    if  display   WeakInteractionMoment.displayAmplitude("Anapole moment amplitude:  ", "H^(anapole)", nm,
                                                          finalLevel, initialLevel, amplitude)    end

    return( amplitude )
end


"""
`WeakInteractionMoment.displayAmplitude(label::String, opName::String, nm::Nuclear.Model, finalLevel::Level, initialLevel::Level,`
                                        `amplitude::ComplexF64)`
    ... to print one amplitude to screen and, if it is open, to the summary stream; nothing is returned.

        The two parts are printed with a separator between them: without it the real and imaginary parts run together into one unparsable
        token, and `TestFrames.testCompareLines` can then only compare that token as a STRING, so that a change in the last digit fails the
        test outright instead of being weighed against its tolerance.
"""
function displayAmplitude(label::String, opName::String, nm::Nuclear.Model, finalLevel::Level, initialLevel::Level,
                          amplitude::ComplexF64)
    sa = @sprintf("%.8e", amplitude.re) * "  " * @sprintf("%.8e", amplitude.im)
    sb = "   " * label *
         "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || $opName ($(Nuclear.name(nm.model))) ||" *
         " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa
    println(sb)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    println(iostream, sb)    end

    return( nothing )
end


"""
`WeakInteractionMoment.nuclearDensity(nm::Nuclear.Model, grid::Radial.Grid)`
    ... to compute the NORMALIZED nuclear density rho(r) on the given grid, i.e. with INT rho(r) 4 pi r^2 dr = 1, so that the nuclear charge
        or weak charge enters only through the prefactor of an operator and never through the density.  An Array{Float64,1} over the grid
        points is returned.

        `Nuclear.fermiDistributedNucleus` cannot serve here: it returns the effective charge Z(r) of the POTENTIAL, not a density, and its
        Fermi shape is a local closure that cannot be reused.  A point nucleus is rejected rather than approximated, since a delta density
        would have every one of these amplitudes depend on the orbital exactly at the origin, where a finite grid says least.
"""
function nuclearDensity(nm::Nuclear.Model, grid::Radial.Grid)
    rho = zeros( grid.NoPoints )
    if      typeof(nm.model) == Nuclear.PointNucleus
        error("A point nucleus has rho(r) = delta(r), for which the weak-interaction amplitudes of this module are not defined on a " *
              "radial grid: they would be fixed by the orbitals exactly at the origin. Use UniformNucleus() or FermiNucleus() instead.")
    elseif  typeof(nm.model) == Nuclear.UniformNucleus
        # a homogeneously charged sphere of radius R, related to the rms radius by R = sqrt(5/3) R_rms
        R = sqrt(5.0/3.0) * Defaults.convertUnits("length: from fm to atomic", nm.radius)
        for  i = 1:grid.NoPoints
            if  grid.r[i] <= R      rho[i] = 1.0    end
        end
    else
        # every Fermi-type model uses the two-parameter shape 1/(1 + exp((r-c)/a))
        c = Defaults.convertUnits("length: from fm to atomic", Nuclear.computeFermiBParameter(nm.radius))
        a = Defaults.convertUnits("length: from fm to atomic", Nuclear.fermiA)
        for  i = 1:grid.NoPoints
            rho[i] = 1.0 / (1.0 + exp( (grid.r[i] - c)/a ))
        end
    end
    #
    wn = 0.
    for  i = 1:grid.NoPoints    wn = wn + rho[i] * 4pi * grid.r[i]^2 * grid.wr[i]    end
    wn == 0.  &&  error("The nuclear density integrates to zero on this grid; the nucleus is smaller than the innermost grid point.")
    for  i = 1:grid.NoPoints    rho[i] = rho[i] / wn    end

    return( rho )
end


"""
`WeakInteractionMoment.oneParticleAmplitude(rank::Int64, kernel::Function, finalLevel::Level, initialLevel::Level)`
    ... to contract a one-electron kernel with the spin-angular coefficients of the two levels and with their mixing coefficients, giving
        the many-electron REDUCED matrix element <alpha_f J_f || H || alpha_i J_i>.  A value::ComplexF64 is returned.

        `SpinAngularGaigalas.computeCoefficients` does NOT return its coefficients in one normalization, and the difference is exactly what this
        function has to absorb.  For rank >= 1 (`computeCoefficientsNonScalar`) the contraction yields the reduced matrix element up to a
        factor sqrt(2J_f+1); for rank 0 (`computeCoefficientsScalar`) it yields the ORDINARY matrix element, as it must, since its main
        client is the one-body Hamiltonian -- and the reduced one then follows from the Wigner-Eckart theorem as sqrt(2J_f+1) times it.
        Both factors were measured against exact one-electron reduced matrix elements rather than assumed; see the module's example branch.
"""
function oneParticleAmplitude(rank::Int64, kernel::Function, finalLevel::Level, initialLevel::Level)
    if  initialLevel.basis.subshells == finalLevel.basis.subshells
        iLevel = initialLevel;   fLevel = finalLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
        iLevel    = Level(initialLevel, subshells);     fLevel = Level(finalLevel, subshells)
    end
    nf = length(fLevel.basis.csfs);     ni = length(iLevel.basis.csfs)
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        for  s = 1:ni
            if  fLevel.mc[r] == 0.  ||  iLevel.mc[s] == 0.    continue    end
            # the parity field of OneParticleOperator is not read by SpinAngular; plus is used, as by every other caller
            opa = SpinAngular.OneParticleOperator(rank, Basics.minus)
            wa  = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], iLevel.basis.subshells)
            for  coeff in wa
                matrix[r,s] = matrix[r,s] + coeff.T * kernel(fLevel.basis.orbitals[coeff.a], iLevel.basis.orbitals[coeff.b])
            end
        end
    end
    # rank 0 needs (2J_f+1) and rank >= 1 only sqrt(2J_f+1); see the docstring above.  Measured, not assumed: a one-electron
    # weak-charge amplitude built here came out 1/sqrt(2) times the exact reduced matrix element while this factor read
    # sqrt(2J_f+1) for both cases (21-Aug-2026).
    if  rank == 0   wNorm = Basics.twice(fLevel.J) + 1.0
    else            wNorm = sqrt( Basics.twice(fLevel.J) + 1.0 )
    end

    return( ComplexF64( wNorm * transpose(fLevel.mc) * matrix * iLevel.mc ) )
end


"""
`WeakInteractionMoment.radialDerivative(f::Array{Float64,1}, grid::Radial.Grid)`
    ... to compute df/dr on the given grid by central differences, with one-sided differences at the two ends; an Array{Float64,1} is
        returned.  The grid is not uniform in r, so the spacing is taken from `grid.r` itself rather than assumed.
"""
function radialDerivative(f::Array{Float64,1}, grid::Radial.Grid)
    n  = min(length(f), grid.NoPoints);    df = zeros(n)
    for  i = 2:n-1    df[i] = (f[i+1] - f[i-1]) / (grid.r[i+1] - grid.r[i-1])    end
    if  n >= 2    df[1] = (f[2] - f[1]) / (grid.r[2] - grid.r[1]);    df[n] = (f[n] - f[n-1]) / (grid.r[n] - grid.r[n-1])    end

    return( df )
end


"""
`WeakInteractionMoment.radialIntegralPPplus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute INT weight(r) [P_a P_b + Q_a Q_b] dr, the combination of an operator that does not mix the large and small components;
        a value::Float64 is returned.
"""
function radialIntegralPPplus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)
    mtp = min( size(a.P, 1), size(b.P, 1), length(weight) )
    wa  = 0.
    for  i = 2:mtp   wa = wa + weight[i] * (a.P[i]*b.P[i] + a.Q[i]*b.Q[i]) * grid.wr[i]    end

    return( wa )
end


"""
`WeakInteractionMoment.radialIntegralPQminus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute INT weight(r) [P_a Q_b - Q_a P_b] dr, the antisymmetric large/small combination that gamma_5 produces; a value::Float64
        is returned.  It is the same shape as `RadialIntegrals.GrantILminus`, with a nuclear density in place of a spherical Bessel
        function.
"""
function radialIntegralPQminus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)
    mtp = min( size(a.P, 1), size(b.P, 1), length(weight) )
    wa  = 0.
    for  i = 2:mtp   wa = wa + weight[i] * (a.P[i]*b.Q[i] - a.Q[i]*b.P[i]) * grid.wr[i]    end

    return( wa )
end


"""
`WeakInteractionMoment.radialIntegralPQ(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute INT weight(r) P_a(r) Q_b(r) dr, i.e. the LARGE component of the first orbital against the SMALL component of the second;
        a value::Float64 is returned.  It is the same shape as `RadialIntegrals.GrantIL0`, with a nuclear density in place of a spherical
        Bessel function.

        This deliberately does NOT symmetrize.  An operator built from the Dirac alpha matrix produces two cross terms, P_a Q_b and Q_a P_b,
        and they carry DIFFERENT angular factors; combining them into P_a Q_b + Q_a P_b is legitimate only when those factors happen to be
        equal up to a sign, which is true for gamma_5 but false for alpha.  Assuming it for the anapole cost a factor of -3 in the
        Hermiticity check.  The caller therefore weights the two terms itself; `radialIntegralPQ(w, b, a, grid)` supplies the second.
"""
function radialIntegralPQ(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)
    mtp = min( size(a.P, 1), size(b.P, 1), length(weight) )
    wa  = 0.
    for  i = 2:mtp   wa = wa + weight[i] * a.P[i] * b.Q[i] * grid.wr[i]    end

    return( wa )
end


"""
`WeakInteractionMoment.schiffMomentAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid;`
                                             `display::Bool=false)`
    ... to compute the Schiff-moment amplitude <alpha_f J_f || H^(Schiff) || alpha_i J_i> for the given final and initial level and for the
        nuclear density of the given nuclear model.  The operator is `4 pi S . grad rho(r)`, of rank 1 and odd parity, and -- being a
        function of position rather than a Dirac matrix -- it carries the ELECTRIC-multipole structure, i.e. <kappa_a||C^1||kappa_b> and the
        P*P + Q*Q radial combination, with the density GRADIENT as the radial weight.  A value::ComplexF64 is returned.

        The amplitude is returned for a unit Schiff moment, so that it may be scaled by whatever nuclear value is adopted.
"""
function schiffMomentAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid; display::Bool=false)
    # rank 1 and odd parity, as for the anapole; the zeros are displayed rather than suppressed
    if      finalLevel.parity == initialLevel.parity                               amplitude = ComplexF64(0.)
    elseif  !AngularMomentum.isTriangle(finalLevel.J, AngularJ64(1), initialLevel.J)  amplitude = ComplexF64(0.)
    else
        rho     = WeakInteractionMoment.nuclearDensity(nm, grid)
        drho    = WeakInteractionMoment.radialDerivative(rho, grid)
        kernel  = (orba, orbb) -> begin
            cl = AngularMomentum.CL_reduced_me(orba.subshell, 1, orbb.subshell)
            if  cl == 0.    return( 0. )    end
            cl / sqrt(Basics.subshell_2j(orba.subshell) + 1.0) *
                WeakInteractionMoment.radialIntegralPPplus(drho, orba, orbb, grid)
        end
        amplitude = 4pi * WeakInteractionMoment.oneParticleAmplitude(1, kernel, finalLevel, initialLevel)
    end
    #
    if  display   WeakInteractionMoment.displayAmplitude("Schiff moment amplitude:   ", "H^(Schiff)", nm,
                                                          finalLevel, initialLevel, amplitude)    end

    return( amplitude )
end


"""
`WeakInteractionMoment.sigmaReducedMe(kappaA::Int64, kappaB::Int64)`
    ... to compute the reduced matrix element <Omega_(kappaA) || sigma || Omega_(kappaB)> of the Pauli spin operator between two
        spin-angular functions; a value::Float64 is returned.

        sigma acts on the spin alone, so it cannot change the orbital angular momentum and the result vanishes unless l_A = l_B.  With that
        understood, the two functions are just |(l 1/2) j> couplings and the Wigner-Eckart theorem for an operator acting on the second
        factor gives

            <(l 1/2)j_A || sigma || (l 1/2)j_B> = (-1)^(l + 1/2 + j_B + 1) sqrt((2j_A+1)(2j_B+1)) {1/2  j_A  l ; j_B  1/2  1} sqrt(6),

        the sqrt(6) being <1/2||sigma||1/2>.  Note that j_A and j_B need NOT be equal: for a given l the two values l +- 1/2 are connected,
        which is why the anapole amplitude reaches p_1/2 <-> d_3/2 as well as s_1/2 <-> p_1/2 once the -kappa_B of the small component is
        taken into account.
"""
function sigmaReducedMe(kappaA::Int64, kappaB::Int64)
    shA = Subshell(9, kappaA);    shB = Subshell(9, kappaB)
    lA  = Basics.subshell_l(shA);    lB = Basics.subshell_l(shB)
    if  lA != lB    return( 0. )    end
    jA  = Basics.subshell_j(shA);    jB = Basics.subshell_j(shB)
    wa  = AngularMomentum.phaseFactor([AngularJ64(lA), +1, AngularJ64(1//2), +1, jB, +1, AngularJ64(1)]) *
          sqrt( (Basics.twice(jA) + 1.0) * (Basics.twice(jB) + 1.0) ) * sqrt(6.0) *
          AngularMomentum.Wigner_6j(AngularJ64(1//2), jA, AngularJ64(lA), jB, AngularJ64(1//2), AngularJ64(1))

    return( wa )
end


"""
`WeakInteractionMoment.weakCharge(nm::Nuclear.Model)`
    ... to compute the nuclear weak charge Q_W = -N + Z(1 - 4 sin^2 theta_W) from the given nuclear model, with the neutron number taken as
        N = round(mass) - Z; a value::Float64 is returned.

        `Nuclear.Model` carries no neutron number and no weak charge, so both are derived here rather than by extending that struct, whose
        eight-field constructor is called positionally throughout src/ and examples/.
"""
function weakCharge(nm::Nuclear.Model)
    N = round(nm.mass) - nm.Z
    N < 0.  &&  error("A negative neutron number, N = $N, follows from mass = $(nm.mass) and Z = $(nm.Z); give the nuclear model a " *
                      "physical mass number before asking for a weak charge.")

    return( -N + nm.Z * (1.0 - 4.0*WeakInteractionMoment.sinThetaW2) )
end


"""
`WeakInteractionMoment.weakChargeAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid;`
                                           `display::Bool=false)`
    ... to compute the weak-charge amplitude <alpha_f J_f || H^(weak-charge) || alpha_i J_i> for the given final and initial level and for
        the nuclear density of the given nuclear model.  A value::ComplexF64 is returned.

        The operator `(G_F/2sqrt2) Q_W rho(r) gamma_5` is a PSEUDOSCALAR: rank 0 and odd parity.  It therefore connects only levels of the
        SAME total angular momentum and OPPOSITE parity, and its one-electron matrix element vanishes unless kappa_b = -kappa_a, i.e.
        between the two orbitals of equal j and opposite parity such as s_1/2 and p_1/2.  All three conditions are enforced here rather
        than left to cancel numerically.

        The result is PURELY IMAGINARY for real radial orbitals, since gamma_5 produces the antisymmetric P_a Q_b - Q_a P_b combination.
"""
function weakChargeAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid; display::Bool=false)
    # A pseudoscalar connects only equal J and opposite parity.  The two guards are kept EXPLICIT rather than left to cancel
    # numerically, and the resulting zero is still displayed, because a selection rule that is visibly satisfied is worth more
    # than one that is silently assumed -- it is what the test suite now asserts in place of the former fabricated numbers.
    if      finalLevel.parity == initialLevel.parity   amplitude = ComplexF64(0.)
    elseif  finalLevel.J      != initialLevel.J        amplitude = ComplexF64(0.)
    else
        rho     = WeakInteractionMoment.nuclearDensity(nm, grid)
        prefac  = WeakInteractionMoment.GF / (2.0*sqrt(2.0)) * WeakInteractionMoment.weakCharge(nm)
        kernel  = (orba, orbb) -> begin
            if  orba.subshell.kappa != -orbb.subshell.kappa    return( ComplexF64(0.) )    end
            im * WeakInteractionMoment.radialIntegralPQminus(rho, orba, orbb, grid)
        end
        amplitude = prefac * WeakInteractionMoment.oneParticleAmplitude(0, kernel, finalLevel, initialLevel)
    end
    #
    if  display   WeakInteractionMoment.displayAmplitude("weak-charge amplitude:     ", "H^(weak-charge)", nm,
                                                          finalLevel, initialLevel, amplitude)    end

    return( amplitude )
end

end # module
