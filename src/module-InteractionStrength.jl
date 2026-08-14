
"""
`module JAC.InteractionStrength`  
... a submodel of JAC that contains all methods for evaluating the interaction strength (reduced matrix elements) 
    for various atomic interactions.
"""
module InteractionStrength


using  GSL, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals


"""
`struct  InteractionStrength.XLCoefficient`  ... defines a type for coefficients of the two-electron (Breit) interaction

    + kind      ::Char       ... Kind of integral, either 'S' or 'T'
    + nu        ::Int64      ... Rank of the integral.
    + a         ::Orbital    ... Orbitals a, b, c, d.
    + b         ::Orbital
    + c         ::Orbital
    + d         ::Orbital
    + coeff     ::Float64    ... corresponding coefficient.
"""
struct XLCoefficient 
    kind        ::Char
    nu          ::Int64
    a           ::Orbital
    b           ::Orbital
    c           ::Orbital
    d           ::Orbital
    coeff       ::Float64 
end 


"""
`InteractionStrength.bosonShift(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the  <a|| h^(boson-field) ||b>  reduced matrix element of the boson-field shift Hamiltonian for orbital 
        functions a, b. This boson-field shift Hamiltonian just refers to the effective potential of the given 
        isotope due to the (assumed) boson mass. A value::Float64 is returned.  
"""
function bosonShift(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)
    wa = RadialIntegrals.isotope_boson(a, b, potential, grid) 
    ## wa = RadialIntegrals.isotope_boson(a, b, potential, grid) 
    ## println("**  <$(a.subshell) || h^(boson-field shift) || $(b.subshell)>  = $wa" )
    return( wa )
end


"""
`InteractionStrength.dipole(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the  <a|| d ||b>  reduced matrix element of the dipole operator for orbital functions a, b. 
        A value::Float64 is returned. 
"""
function dipole(a::Orbital, b::Orbital, grid::Radial.Grid)
    wa = AngularMomentum.CL_reduced_me(a.subshell, 1, b.subshell) * RadialIntegrals.rkDiagonal(1, a, b, grid)
    return( wa )
end


"""
`InteractionStrength.eMultipole(k::Int64, a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the  <a|| t^(Ek) ||b>  reduced matrix element of the dipole operator for orbital functions a, b. 
        A value::Float64 is returned. 
"""
function eMultipole(k::Int64, a::Orbital, b::Orbital, grid::Radial.Grid)
    ## the former CL_reduced_me_rb convention: Grant's reduced matrix element divided by sqrt(2 j_a + 1),
    ## with j_a from the FIRST argument as passed (10-Aug-2026).
    wa = AngularMomentum.CL_reduced_me(a.subshell, k, b.subshell) / sqrt( Basics.subshell_2j(a.subshell) + 1 ) *
            RadialIntegrals.rkDiagonal(k, a, b, grid)
    return( wa )
end


"""
`InteractionStrength.fieldShift(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the  <a|| h^(field-shift) ||b>  reduced matrix element of the field-shift Hamiltonian for orbital 
        functions a, b. This field-shift Hamiltonian just refers to the difference of the nuclear potential 
        deltaPotential of two isotopes, and which is already divided by the difference of the mean-square radii. 
        A value::Float64 is returned.  
"""
function fieldShift(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)
    wa = RadialIntegrals.isotope_field(a, b, deltaPotential, grid) 
    return( wa )
end


"""
`InteractionStrength.hamiltonian_nms(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the  <a|| h_nms ||b>  reduced matrix element of the normal-mass-shift Hamiltonian for orbital 
        functions a, b. A value::Float64 is returned.  For details, see Naze et al., CPC 184 (2013) 2187, Eq. (37).
"""
function hamiltonian_nms(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    if  a.subshell.kappa != b.subshell.kappa   return( 0. )   end
    wa = RadialIntegrals.isotope_nms(a, b, nm.Z, grid)
    ## println("**  <$(a.subshell) || h^nms || $(b.subshell)>  = $wa" )
    return( wa )
end



"""
`InteractionStrength.hfs_tM1(a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... computes the <a|| t^(1) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment
        of the nucleus for orbital functions a, b. A value::Float64 is returned.

        Note (26-Jul-2026): was missing the alpha (fine-structure constant) prefactor required by Andersson &
        Jonsson (2008), CPC 178, Eq. (49): <n_a kappa_a || t^(1) || n_b kappa_b> = -alpha(kappa_a+kappa_b)
        <-kappa_a||C^(1)||kappa_b>[r^-2] -- confirmed by direct re-reading of the paper (p. 161). This alone
        does not fully resolve the long-standing H(1s) HFS A-constant discrepancy; see
        examples/example-Cb.jl and memory project_zeeman_hfs_bugs.md for the residual after this fix.
"""
function hfs_tM1(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - Defaults.getDefaults("alpha") * (a.subshell.kappa + b.subshell.kappa) *
                AngularMomentum.CL_reduced_me(minusa, 1, b.subshell)
    wc =   RadialIntegrals.rkNonDiagonal(-2, a, b, grid)
    wa =   wb * wc
    #
    return( wa )
end

"""
`InteractionStrength.hfs_tM2(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(M2) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment 
        of the nucleus for orbital functions a, b. A value::Float64 is returned.  
"""
function hfs_tM2(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - (a.subshell.kappa + b.subshell.kappa) * AngularMomentum.CL_reduced_me(minusa, 2, b.subshell)
    wc =   RadialIntegrals.rkNonDiagonal(-3, a, b, grid)/2
    wa =   wb * wc
    #
    return( wa )
end

"""
`InteractionStrength.hfs_tM3(a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... computes the <a|| t^(M3) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment
        of the nucleus for orbital functions a, b. A value::Float64 is returned.

        Note (30-Jul-2026): was missing the same alpha (fine-structure constant) prefactor found missing in
        hfs_tM1 on 26-Jul-2026 -- Andersson & Jonsson (2008), CPC 178, Eq. (49) is generic in the multipole
        rank L, so the same alpha(kappa_a+kappa_b) prefactor applies for L=3 (M3) as for L=1 (M1). See
        examples/example-Cb.jl and memory project_zeeman_hfs_bugs.md.
"""
function hfs_tM3(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - Defaults.getDefaults("alpha") * (a.subshell.kappa + b.subshell.kappa) *
                AngularMomentum.CL_reduced_me(minusa, 3, b.subshell)
    wc =   RadialIntegrals.rkNonDiagonal(-4, a, b, grid)/3
    wa =   wb * wc
    #
    return( wa )
end

"""
`InteractionStrength.hfs_tE2(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(2) ||b> reduced matrix element for the HFS coupling to the electric-quadrupole moment of 
        the nucleus for orbital functions a, b. A value::Float64 is returned.  
    """
function hfs_tE2(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    wb = - AngularMomentum.CL_reduced_me(a.subshell, 2, b.subshell)
    ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wa =   wb * wc
    #
    ## println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
    return( wa )
end

"""
`InteractionStrength.hfs_tE1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(E1) ||b> reduced matrix element for the HFS coupling to the electric-quadrupole moment of 
        the nucleus for orbital functions a, b. A value::Float64 is returned.  
    """
function hfs_tE1(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    wb = - AngularMomentum.CL_reduced_me(a.subshell, 1, b.subshell)
    ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wc =   RadialIntegrals.rkDiagonal(-2, a, b, grid)
    wa =   wb * wc
    #
    ## println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
    return( wa )
end

"""
`InteractionStrength.hfs_tE3(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(E3) ||b> reduced matrix element for the HFS coupling to the electric-quadrupole moment of 
        the nucleus for orbital functions a, b. A value::Float64 is returned.  
    """
function hfs_tE3(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    wb = - AngularMomentum.CL_reduced_me(a.subshell, 3, b.subshell)
    ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wc =   RadialIntegrals.rkDiagonal(-4, a, b, grid)
    wa =   wb * wc
    #
    ## println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
    return( wa )
end












"""
`InteractionStrength.MabEmission(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
    ... computes the single-electron reduced matrix element <a || O^(Mp, emission) || b> of the electron-photon
        multipole interaction, for the multipole mp of the radiation field at photon energy omega and in the
        given gauge. A value::Float64 is returned.

        THE ONE VERSION. Introduced 09-Aug-2026 to replace a family of seven variants that had accumulated along
        historical routes -- MabEmissionJohnsony, MabEmissionJohnsony_Wu, MbaEmissionJohnsonx, MbaEmissionCheng,
        MbaAbsorptionCheng, MbaEmissionAndreyOld and MbaEmissionMigdalek -- which differed in normalisation, in
        return type, in sign and in argument order, so that amplitudes computed through different ones could not
        be added.

        ARGUMENT ORDER, stated because the family disagreed about it: the orbitals are passed as (b, a) and the
        matrix element returned is <a || O || b>. That is, the SECOND orbital argument is the one that appears
        on the LEFT of the matrix element. This is the order every live call site already uses.

        THE GAUGES, following Grant, J. Phys. B 7, 12 (1974): `Coulomb` is the velocity form and `Babushkin` the
        length form of the electric multipole operator; `Magnetic` is used for magnetic multipoles. For exact
        one-body eigenfunctions the two electric forms must agree at the physical photon energy, and here they
        do -- H 1s-2p gives Coulomb/Babushkin = 1.000000 on shell. Away from that energy they legitimately
        differ: the length form scales with the omega it is given, while the velocity form carries the level
        difference through <f|p|i> = i m (E_f - E_i) <f|r|i>, so an OFF-SHELL call is not gauge invariant and
        must not be expected to be.

        THE VALUE IS REAL, and this deserves a word because the literature usually writes it as complex. The
        Racah phase of the angular factor is (-1)^(j_b + 1/2); j_b is half-integer, so that exponent is an
        integer and the phase is +-1, alternating with j_b. Until 09-Aug-2026 `AngularMomentum.JohnsonI` raised
        (-1+0im) to `jb.num + 1/2`, i.e. to a HALF-integer, which returned a constant -i for every j_b: it turned
        every matrix element imaginary and, worse, discarded the alternating sign. With that corrected the
        quantity is real, as it should be. Any overall factor i that a given convention prefers is a global
        phase and cancels from every observable; the j-dependent sign does NOT cancel and is carried here.

        WHAT IS INSIDE, AND WHAT THE CALLER MUST SUPPLY. This function returns the SINGLE-ELECTRON reduced
        matrix element and nothing else. A many-electron amplitude is assembled by the caller as

            amplitude = sum_(r,s) c_r c_s sum_coeff  coeff.T * MabEmission(...) / sqrt(2j_a+1) * sqrt(2J_f+1)

        with the angular coefficients from SpinAngular and the CI mixing coefficients. NO FURTHER MULTIPOLARITY
        FACTOR IS TO BE APPLIED: a factor sqrt((2L+1)(L+1)/L) appears commented out at the PhotoEmission call
        site, and it is commented out CORRECTLY -- restoring it would multiply E1 by 2.449 and E2 by 2.739 and
        destroy the agreement documented below.

        VALIDATED, and additive across multipoles. Against Jitrik & Bunge, J. Phys. Chem. Ref. Data 33, 1059
        (2004), hydrogen Z = 1, point-nucleus Dirac -- the same model as the test:

            transition        JAC A(Cou)      JAC A(Bab)      reference      ratio
            E1 1s-2p_1/2      6.268354e+08    6.268020e+08    6.2649e+08     1.000498
            M1 1s-2s_1/2      2.481059e-06    2.481059e-06    2.4946e-06     0.994572
            E2 1s-3d_3/2      5.940766e+02    5.940251e+02    5.937500e+02   1.000463

        The three ratios agree with one another, which is the statement that E1, M1 and E2 amplitudes computed
        through this function stand on a common footing and may be added as they occur in the full
        electron-photon interaction.

        THE COULOMB GAUGE IS STRUCTURALLY ZERO FOR A DIAGONAL ELECTRIC MULTIPOLE, i.e. whenever the two
        orbitals belong to the SAME subshell and hence kapa == kapb. The Coulomb branch below is built only
        from I^+ and I^-, and there both die at once: the prefactor (kapa-kapb) vanishes, and
        I^-(a,a) = int j_L (P_a Q_a - Q_a P_a) dr is identically zero. The Babushkin branch keeps
        J_L = int j_L (P_a^2 + Q_a^2) dr and stays finite. This is not a defect of either branch: the velocity
        form of an electric multipole between two states built on the same orbital vanishes at the one-body
        level, and the amplitude of a transition INSIDE one configuration (2p^4 1S_0 - 1D_2, say) is carried
        entirely by the length form. Use Babushkin for such transitions, and do not read the resulting
        Coulomb/Babushkin ratio as a gauge-consistency check -- there is nothing there to compare.
        The magnetic branch is unaffected, since it carries (kapa+kapb) = 2*kappa and I^+(a,a) != 0; that is
        why M1 fine-structure rates inside a term come out right to a few parts in 1000.
"""
function MabEmission(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c")
    #
    if       gauge == Basics.Magnetic
        JohnsonI = AngularMomentum.JohnsonI(-kapa, kapb, AngularJ64(mp.L))
        wa     = JohnsonI * (kapa + kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
        #
    elseif   gauge == Basics.Babushkin
        JohnsonI = AngularMomentum.JohnsonI(kapa, kapb, AngularJ64(mp.L))
        ## ORIENTATION (corrected 09-Aug-2026). Unlike the velocity form below, the length form is not
        ## antisymmetric term by term under exchange of the two orbitals: GrantJL is symmetric, GrantILminus
        ## is antisymmetric and (kapa-kapb) changes sign. The written expression is therefore valid for ONE
        ## assignment only -- the one with `a` the more strongly bound orbital -- and the mirrored order needs
        ## the GrantJL term negated, the other two flips cancelling against the required overall sign. Call
        ## sites disagreed about the order, so the orientation is settled here from the orbital energies
        ## instead of being left to the caller. For a == b both orders coincide and eps = +1 is the right one.
        eps    = a.energy <= b.energy   ?   1.0   :   -1.0
        wr     = eps * RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapa-kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr +  RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid)
        wa     = JohnsonI * wr
        #
    elseif   gauge == Basics.Coulomb
        JohnsonI = AngularMomentum.JohnsonI(kapa, kapb, AngularJ64(mp.L))
        wr       = (1 - mp.L/(2mp.L+1)) * RadialIntegrals.GrantILplus(mp.L-1, q, a, b, grid::Radial.Grid)  -
                    mp.L/(2mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid)
        wr       = -(kapa-kapb) / (mp.L+1) * wr
        wr       = wr  +  mp.L/(2mp.L+1) * RadialIntegrals.GrantILminus(mp.L-1, q, a, b, grid::Radial.Grid)
        wr       = wr  +  mp.L/(2mp.L+1) * RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid)
        wa       = JohnsonI * wr
    else     error("stop a")
    end

    return( wa )
end







"""
`InteractionStrength.multipoleTransition(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) multipole-transition interaction strength 
        <b || T^(Mp, absorption) || a> due to Johnson (2007) for the interaction with the Mp multipole component of the
        radiation field and the transition frequency omega, and within the given gauge. A value::Float64 is returned.  
"""
function multipoleTransition(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    function besselPrime_jl(L::Int64, x::Float64)    return( GSL.sf_bessel_jl(L-1, x) - (L+1)/x * GSL.sf_bessel_jl(L, x) )       end

    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    mtp  = min(size(a.P, 1), size(b.P, 1))
    #
    if       gauge == Basics.Magnetic
        ChengI = AngularMomentum.ChengI(-kapa, kapb, AngularJ64(mp.L));   if  abs(ChengI) < 1.0e-10  return( 0. )   end
        wa = Complex(0.)
        for  i = 2:mtp
            wa = wa + (kapa+kapb) / (mp.L+1) * GSL.sf_bessel_jl(mp.L, q * grid.r[i]) * (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * grid.wr[i]  
        end
        wa = ChengI * wa
        #
    elseif   gauge == Basics.Velocity
        ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L));    if  abs(ChengI) < 1.0e-10  return( 0. )   end
        wa = Complex(0.)
        for  i = 2:mtp
            wa = wa - (kapa-kapb) / (mp.L+1) * 
                        ( besselPrime_jl(mp.L, q * grid.r[i]) + GSL.sf_bessel_jl(mp.L, q * grid.r[i]) / (q * grid.r[i]) ) * 
                        (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * grid.wr[i]
            wa = wa + mp.L * GSL.sf_bessel_jl(mp.L, q * grid.r[i]) / (q * grid.r[i]) * (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) * grid.wr[i]  
        end
        wa = ChengI * wa
        #
    elseif   gauge == Basics.Length
        ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L));    if  abs(ChengI) < 1.0e-10  return( 0. )   end
        wa = Complex(0.)
        for  i = 2:mtp
            wa = wa + GSL.sf_bessel_jl(mp.L, q * grid.r[i]) * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * grid.wr[i] 
                    + GSL.sf_bessel_jl(mp.L+1, q * grid.r[i]) * 
                        ( (kapa-kapb) / (mp.L+1) * (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) +
                        (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) ) * grid.wr[i]
        end
        wa = ChengI * wa
        #
    else     error("stop a")
    end

    return( wa )
end


"""
`InteractionStrength.schiffMoment(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the  <a|| h^(Schiff-moment) ||b>  reduced matrix element of the Schiff-moment Hamiltonian for orbital 
        functions a, b and for the nuclear density as given by the nuclear model. A value::Float64 is returned.  
"""
function schiffMoment(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    printstyled("\nWarning -- InteractionStrength.schiffMoment():: Not yet implemented.", color=:cyan)
    wb = 1.0 + 2.0im
    return( wb )
end


"""
`InteractionStrength.weakCharge(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the  <a|| h^(weak-charge) ||b>  reduced matrix element of the weak-charge Hamiltonian for orbital functions 
        a, b and for the nuclear density as given by the nuclear model. A value::Float64 is returned.  
"""
function weakCharge(a::Orbital, b::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    printstyled("\nWarning -- InteractionStrength.weakCharge():: Not yet implemented.", color=:cyan)
    wb = 1.0 + 2.0im
    return( wb )
end


##
## ============================================================================================================
##  THE BREIT INTERACTION IN JAC -- THE REFERENCE FORMULATION
##  Written 13-Aug-2026 as Stage 1 of the frequency-dependent Breit work, BEFORE any code was changed, so
##  that the implementation below can be checked against a statement of the physics rather than against
##  itself.
##
##  PROVENANCE, STATED HONESTLY: the standard references are Grant & Pyper, J. Phys. B 9, 761 (1976) and
##  Grant, "Relativistic Quantum Theory of Atoms and Molecules" (2007), ch. 8.  THIS BLOCK WAS NOT
##  TRANSCRIBED FROM THEM -- it is written from recollection of the standard formulation, cross-checked
##  against secondary sources only.  It must therefore be checked against the primary texts by someone who
##  has them before any number produced from it is trusted.  Where it has been confirmed NUMERICALLY, that
##  is said at the point in question.
## ============================================================================================================
##
## 1. THE OPERATOR.  Exchange of one transverse photon of frequency omega between electrons 1 and 2, in the
##    COULOMB GAUGE:
##
##        B(omega) = - (alpha_1 . alpha_2) cos(omega r_12) / r_12
##                   + (alpha_1 . grad_1)(alpha_2 . grad_2) [cos(omega r_12) - 1] / (omega^2 r_12)
##
##    The second term is finite as omega -> 0 only through a cancellation: [cos(x) - 1]/x^2 -> -1/2.
##
## 2. UNITS OF omega -- AND A DEFECT.  omega is a WAVE NUMBER, not an energy.  In atomic units
##
##        omega = |E_a - E_c| / c,        c = 137.036 ,
##
##    so that omega*r is the dimensionless phase the Bessel functions below require.  Retardation is small
##    precisely because omega*r ~ (Delta E) r / c << 1 for atomic transitions.
##
##    XL_Breit_densities below formed  omg_ac = factor * abs(E_a - E_c)  and NEVER DIVIDED BY c until
##    13-Aug-2026, feeding an energy in Hartree into sf_bessel_jl(nu, omega*r) as if it were a wave number
##    and overstating the phase by c = 137.  FIXED.  It never affected a published number, since no branch
##    consumed those values (see point 6).
##
## 3. THE TWO omega -> 0 LIMITS ARE DIFFERENT OPERATORS, AND THE DIFFERENCE IS A GAUGE CHOICE.
##
##        Coulomb gauge, omega -> 0:   B = - 1/(2 r_12) [ alpha_1.alpha_2 + (alpha_1.rhat)(alpha_2.rhat) ]
##                                         ... the BREIT operator            -> JAC: CoulombBreit(0.)
##        Feynman gauge, omega -> 0:   G = - alpha_1.alpha_2 / r_12
##                                         ... the GAUNT operator            -> JAC: CoulombGaunt
##
##    These are NOT the same approximation.  The total one-photon exchange is gauge independent; truncating
##    at omega -> 0 is what breaks that, which is why two gauges leave two different instantaneous operators.
##    The retardation correction to instantaneous Gaunt is O(alpha^2) -- LARGER than the omega-independent
##    Coulomb-gauge term itself.  A user choosing CoulombGaunt over CoulombBreit(0.) is therefore making an
##    undocumented gauge choice, and the frequency dependence does NOT vanish in either gauge.
##
## 4. MULTIPOLE DECOMPOSITION.  XL_Breit_coefficients splits the operator into two families, and this split
##    IS the Gaunt/retardation split (onlyGaunt returns after the 'T' blocks):
##
##        'T'  (magnetic / GAUNT)      nu = L-1, L, L+1,   four mu permutations each
##        'S'  (RETARDATION)           nu = L +- 1
##
##    The many-electron angular factors are NOT computed here: they come from SpinAngular as
##    Coefficient2p(nu, a, b, c, d, V), and the CI matrix forms V * XL_Breit(nu, ...).  The coefficients
##    below belong to the OPERATOR, not to the CSFs.
##
## 5. THE RADIAL KERNELS.  The two kinds enter DIFFERENTLY, which is easy to get wrong:
##
##        Gaunt / 'T':   the static kernel Ubar_nu = r_<^nu / r_>^(nu+1) times a FREQUENCY FACTOR
##                           V_nu = -omega (2nu+1) j_nu(omega r_<) y_nu(omega r_>)  ->  1   as omega -> 0.
##
##        Retardation / 'S':  the kernel W_(L-1,L+1,L) of Grant & Pyper equation (6) IN FULL -- it is not a
##                            factor on Ubar_nu, and for r_1 < r_2 it is a difference of two 1/omega^2-
##                            divergent pieces, which is the whole numerical difficulty.
##
##    Both are evaluated through the normalised phi, psi of besselPhiPsi, so the cancellations are taken
##    ALGEBRAICALLY and each kernel reaches its omega -> 0 limit by construction rather than by a special case.
##    See the closures V() and W() in XL_Breit_densities for the derivations and their verification.
##
##    MEASURED 13-Aug-2026, on the forms that stood here before: V() omitted the leading factor omega
##    (V*omega/static -> 1.000001 at omega = 1e-3 for nu = 1, while V/static -> 1000), and W() was equation (6)
##    with the omega missing from its first term AND that term's sign reversed, so its cancellation failed
##    outright -- neither W, W*omega nor W*omega^2 converged.  Neither could ever have been switched on.
##
## 6. WHAT IS COMPUTED.  Since 14-Aug-2026 both parts are frequency dependent, and the argument `factor` of
##    CoulombBreit(factor) does exactly what its name says: it scales omega, so that
##
##        factor = 0  ->  omega = 0, the standard frequency-independent Breit interaction, recovered as the
##                        EXACT limit of the same expressions rather than by a separate code path;
##        factor = 1  ->  the full frequency-dependent interaction at omega = |E_a - E_c| / c.
##
##    Before that date nothing frequency dependent was computed at all: omg_ac and omg_bd were formed and then
##    discarded, factor == 0 taking wy = 1 and factor == 1 taking wy = 1.05, a placeholder with no derivation
##    that is now gone.  No JAC number published before 14-Aug-2026 contains a retardation correction.
##
##    NOT included: for a given 'S' permutation only the radial ordering carried by that permutation is
##    summed, its transpose supplying the other.  This is exact at omega = 0 and drops a term of relative
##    order omega^2 WITHIN the retardation -- O(omega^4) in the interaction -- since W's r_1 > r_2 branch is
##    itself O(omega^2).  With omega ~ 1e-2 for the transitions of interest that is a ~1e-8 relative effect.
##
##
## 7. AUDIT OF THE ANGULAR COEFFICIENTS, 13-Aug-2026.  XL_Breit_coefficients was compared TERM BY TERM against
##    Grant & Pyper's table 2 and against GRASP2018's rci90/cxk.f90, which implements the same table.
##    RESULT: THE ANGULAR DECOMPOSITION IS CORRECT.  Every weight, every kappa factor and every mu ordering
##    agrees.  In GRASP's notation DK1 = kappa_c - kappa_a, DK2 = kappa_d - kappa_b, F1..F4 = DK -+ K,
##    G1..G4 = DK -+ (K+1), H = the two CL reduced matrix elements with the odd-L sign:
##
##      'T', nu = L    : JAC -(ka+kc)(kb+kd)/(L(L+1)), equal for all mu   ==  GRASP S(1..4) = -(KA+KC)(KD+KB)H/(K(K+1))
##      'T', nu = L+1  : JAC L/((L+1)(2L+1)(2L+3)), four kappa products   ==  GRASP A = H*K/((K+1)(2K+1)(2K+3)),
##                       G1G3, G2G4, G1G4, G2G3
##      'T', nu = L-1  : JAC (L+1)/(L(2L-1)(2L+1))                        ==  GRASP A = H*(K+1)/(K(2K+1)(2K-1)),
##                       F2F4, F1F3, F2F3, F1F4
##      'S', all 8 mu  : JAC 1/(2L+1)^2 times, in order,
##                       F2G3, F4G1, F1G4, F3G2, F2G4, F3G1, F1G3, F4G2   ==  GRASP S(5), S(6), ..., S(12)
##
##    THE ONE STRUCTURAL DIFFERENCE, and it is deliberate rather than a defect: JAC emits each 'S' coefficient
##    TWICE, at nu = L+1 with weight +[L]/2 and at nu = L-1 with -[L]/2, [L] = 2L+1.  That is precisely Grant's
##    STATIC LIMIT of the retardation kernel, their equations (9)-(10),
##
##        W_(nu-1,nu+1,nu) -> -(1/2)[nu] ( Ubar_(nu-1) - Ubar_(nu+1) ),
##
##    with Ubar the one-sided kernel of equation (10) -- which the r/s loop below realises by summing only
##    s <= r and halving the diagonal, since U = Ubar(1,2) + Ubar(2,1).  Signs included, JAC matches.
##
##    CONSEQUENCE FOR THE FREQUENCY-DEPENDENT RETARDATION: it CANNOT be switched on by changing a multiplier.
##    The 'S' entries currently carry the static limit inside their coefficients, so they would have to be
##    re-emitted as a single term in W_(L-1,L+1,L) (Grant equation 6) instead of two terms in Ubar.  That is
##    the remaining work, and it is a change to XL_Breit_coefficients, not to the kernels.
##
## ============================================================================================================
##


"""
`InteractionStrength.XL_Breit_reset_storage(keep::Bool; printout::Bool=false)`  
    ... resets the global storage of XL_Breit interaction strength; nothing is returned.
"""
function XL_Breit_reset_storage(keep::Bool; printout::Bool=false)
    if  keep
        if  printout   println("  reset GBL_Storage_XL_Breit storage ...")    end
        global GBL_Storage_XL_Breit = Dict{String, Float64}()
    else
    end
    return( nothing )      
end


"""
`InteractionStrength.XL_Breit(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                                eeint::Union{BreitInteraction, CoulombBreit, CoulombGaunt}; keep::Bool=false)`  
    ... computes the the effective Breit interaction strengths X^L_Breit (abcd) or Gaunt interaction strengths 
        X^L_Gaunt (abcd) for given rank L and orbital functions a, b, c and d  at the given grid. 
        For keep=true, the procedure looks up the (global) directory GBL_Storage_XL_Coulomb
        and returns the corresponding value without re-calculation of the interaction strength; it also 'stores' the calculated
        value if not yet included. For keep=false, the interaction strength is always computed on-fly. A value::Float64 is returned. 
        At present ONLY THE ZERO-FREQUENCY LIMIT IS RETURNED, whatever factor is given: the
        frequency-dependent kernels below are not consumed by any branch.  See the reference formulation
        heading this section.
"""
function XL_Breit(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                    eeint::Union{BreitInteraction, CoulombBreit, CoulombGaunt}; keep::Bool=false)
    global GBL_Storage_XL_Breit
    ja2 = Basics.subshell_2j(a.subshell)
    jb2 = Basics.subshell_2j(b.subshell)
    jc2 = Basics.subshell_2j(c.subshell)
    jd2 = Basics.subshell_2j(d.subshell)
    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||  L == 0  
        return( 0. )
    end
    
    # Calculate a reduced number of cofficients for the CoulombGaunt() interaction
    if    typeof(eeint) == CoulombGaunt   onlyGaunt = true;    factor = 0.
    else                                  onlyGaunt = false;   factor = eeint.factor   
    end
    
    # Now distiguish due to the optional argument keep
    if  keep
        sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
        if haskey(GBL_Storage_XL_Breit, sa )
            XL_Breit = GBL_Storage_XL_Breit[sa]
        else
            xcList   = XL_Breit_coefficients(L,a,b,c,d, onlyGaunt=onlyGaunt)
            XL_Breit = XL_Breit_densities(xcList, factor, grid)
            global GBL_Storage_XL_Breit = Base.merge(GBL_Storage_XL_Breit, Dict( sa => XL_Breit))
        end
    else
        xcList   = XL_Breit_coefficients(L,a,b,c,d, onlyGaunt=onlyGaunt)
        XL_Breit = XL_Breit_densities(xcList, factor, grid)
    end
    #

    return( XL_Breit )
end


"""
`InteractionStrength.XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                                    eeint::Union{BreitInteraction, CoulombBreit, CoulombGaunt})`
    ... computes the effective Breit interaction strength X^L_Breit (abcd) for given rank L and orbital functions
        a, b, c and d at the given grid, with the radial integrand damped by exp(-tau r_1 - tau r_2) as the
        DampedSpaceCI Green function approach requires. A value::Float64 is returned. This is XL_Breit with
        that one factor added, so it follows the same gauge, the same frequency dependence and the same
        coefficients; see XL_Breit and the reference formulation heading this section.

        UNTIL 14-Aug-2026 the body of this function was `error("stop a")`, and its signature took SEVEN
        arguments while its only caller -- the Breit branch of Basics.computeMultipletForGreenApproach for
        AtomicState.DampedSpaceCI -- passed EIGHT. The path therefore raised a MethodError before it could
        even reach the stub. It had never been noticed because every use of DampedSpaceCI in the examples and
        in the suite selects a Coulomb-only e-e interaction, so that branch was never taken.
"""
function XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                        eeint::Union{BreitInteraction, CoulombBreit, CoulombGaunt})
    ja2 = Basics.subshell_2j(a.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    jc2 = Basics.subshell_2j(c.subshell);    jd2 = Basics.subshell_2j(d.subshell)
    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||  L == 0
        return( 0. )
    end
    #
    if    typeof(eeint) == CoulombGaunt   onlyGaunt = true;    factor = 0.
    else                                  onlyGaunt = false;   factor = eeint.factor
    end
    xcList = XL_Breit_coefficients(L, a, b, c, d, onlyGaunt=onlyGaunt)
    return( XL_Breit_densities(xcList, factor, grid, tau=tau) )
end


"""
`InteractionStrength.XL_Breit_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital; onlyGaunt::Bool=false)`  
    ... evaluates the combinations and pre-coefficients of the Breit interaction X^L_Breit(abcd) for given
        rank L and orbital functions a, b, c and d; a list xcList::Array{XLCoefficient,1} is returned.

        THESE COEFFICIENTS ARE FREQUENCY INDEPENDENT and are correct for both cases: the omega dependence
        enters only through the radial kernels in XL_Breit_densities, not here.  The list splits into
        kind = 'T' (magnetic / GAUNT, nu = L-1, L, L+1) and kind = 'S' (RETARDATION, nu = L+-1); onlyGaunt
        returns after the 'T' blocks, so that split IS the Gaunt/retardation split.
"""
function XL_Breit_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital; onlyGaunt::Bool=false)
    xcList = XLCoefficient[]
    
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    ## TWO angular prefactors are needed, because the blocks below do not all want the same parity.
    ##
    ##   xc      = <kappa_a || C^L || kappa_c> <kappa_b || C^L || kappa_d>,  nonzero for l_a+l_c+L EVEN.
    ##             Used by the 'T' blocks at nu = L-1 and nu = L+1 and by the whole 'S' block.
    ##   xcFlip  = the same with the kappa sign flipped on the first argument, nonzero for l_a+l_c+L ODD.
    ##             Used by the 'T' block at nu = L.
    ##
    ## The nu = L integrand is a.P * c.Q -- LARGE component against SMALL.  The small component of kappa_c
    ## carries lbar_c = l_c +- 1, so the angular factor there is <kappa_a || C^L || -kappa_c> and its parity
    ## condition is the OPPOSITE of the one that governs the other blocks.  CL_reduced_me depends on kappa
    ## only through 2j = 2|kappa|-1, so flipping the sign leaves the VALUE untouched and moves only the
    ## selection rule; xcFlip is therefore exactly the number this code used before 8f0930b.
    ##
    ## THAT COMMIT (10-Aug-2026) gave CL_reduced_me the parity rule it had genuinely been missing -- correct
    ## in itself, and needed for the Coulomb and multipole call sites -- but here it silently zeroed the
    ## nu = L block: its guard demands l_a+l_c+L odd, which the single shared xc could no longer satisfy.
    ## That block carries the DOMINANT magnetic term, so every Breit and Gaunt number computed between
    ## 10-Aug-2026 and 14-Aug-2026 came out far too small (a factor 3.3 in Gaunt for Cl-like Xe), while the
    ## retardation part, which wants even parity, was untouched.  The test suite did not notice.
    ##
    ## GRASP2018 keeps the same division of labour differently: its CLRX depends only on |kappa_a|, |kappa_b|
    ## and K, vetoes on the triangle condition alone, and evaluates BOTH parity branches, leaving the parity
    ## vetoes to callers such as rci90/cxk.f90.  Spelling the -kappa out here says the same thing explicitly.
    xc     = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) *
             AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    xcFlip = AngularMomentum.CL_reduced_me(Subshell(a.subshell.n, -a.subshell.kappa), L, c.subshell) *
             AngularMomentum.CL_reduced_me(Subshell(b.subshell.n, -b.subshell.kappa), L, d.subshell)
    if   rem(L,2) == 1    xc = - xc;   xcFlip = - xcFlip                       end
    if   abs(xc) < 1.0e-10   &&   abs(xcFlip) < 1.0e-10      return( xcList )  end

    # Consider the individual contributions from sum_nu and sum_mu. First, take T^(nu,L)_mu = R^(nu,L)_mu
    nu = L - 1

    if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
        wa = (L+1) / ( L*(L+L-1)*(L+L+1) ) 
        # mu = 1
        xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa+L) * (d.subshell.kappa-b.subshell.kappa+L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, b, c, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, a, d, c, xcc) )   end
        # mu = 2
        xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa-L) * (d.subshell.kappa-b.subshell.kappa-L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, d, a, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, c, b, a, xcc) )   end
        # mu = 3
        xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa+L) * (d.subshell.kappa-b.subshell.kappa-L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, d, c, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, a, b, c, xcc) )   end
        # mu = 4
        xcc = xc * wa * (c.subshell.kappa-a.subshell.kappa-L) * (d.subshell.kappa-b.subshell.kappa+L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, b, a, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, c, d, a, xcc) )   end
    end
    #
    nu = L

    if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
        ## The ONLY block that wants odd parity, and hence the only consumer of xcFlip; see the note above.
        wa = - (a.subshell.kappa + c.subshell.kappa) * (b.subshell.kappa + d.subshell.kappa) / (L*(L+1))
        # mu = 1
        xcc = xcFlip * wa
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, b, c, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, a, d, c, xcc) )   end
        # mu = 2
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, d, a, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, c, b, a, xcc) )   end
        # mu = 3
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, d, c, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, a, b, c, xcc) )   end
        # mu = 4
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, b, a, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, c, d, a, xcc) )   end
    end
    #
    nu = L + 1

    if  rem(la+lc+nu,2) == 1   &&   rem(lb+ld+nu,2) == 1   &&   L != 0
        wa =  L / ( (L+1)*(L+L+1)*(L+L+3) )
        # mu = 1
        xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa - L - 1) * (d.subshell.kappa - b.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, b, c, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, a, d, c, xcc) )   end
        # mu = 2
        xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa + L + 1) * (d.subshell.kappa - b.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, d, a, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, c, b, a, xcc) )   end
        # mu = 3
        xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa - L - 1) * (d.subshell.kappa - b.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, a, d, c, b, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, d, a, b, c, xcc) )   end
        # mu = 4
        xcc = xc * wa * (c.subshell.kappa - a.subshell.kappa + L + 1) * (d.subshell.kappa - b.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, c, b, a, d, xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('T', nu, b, c, d, a, xcc) )   end
    end
    
    # Return here if onlyGaunt = true
    if   onlyGaunt     return( xcList )     end

    ## Add contributions of the S^k_mu integrals -- the RETARDATION part, Grant & Pyper's B^(2), their eq (8).
    ##
    ## Each mu carries the multipole L itself and the bare coefficient xcc; the whole radial kernel is then
    ## W_(L-1,L+1,L) of their equation (6), supplied by XL_Breit_densities.  The eight mu terms form four
    ## TRANSPOSED PAIRS (1<->2, 3<->4, 5<->6, 7<->8), which is exactly the two-term W(1,2) / W(2,1) structure
    ## of equation (8): each partner is supported on one of the two radial orderings.
    ##
    ## UNTIL 14-Aug-2026 each mu emitted a PAIR of entries instead,
    ##     ('S', L+1, perm, +[L]/2 xcc)  and  ('S', L-1, perm, -[L]/2 xcc),   [L] = 2L+1,
    ## whose sum against the static kernels Ubar_(L+1), Ubar_(L-1) is  -[L]/2 xcc (Ubar_(L-1) - Ubar_(L+1)),
    ## i.e. precisely xcc * W_(L-1,L+1,L) in the limit omega -> 0 of their equations (9)-(10).  That is a
    ## correct static Breit interaction, but it FORECLOSES omega /= 0 in the coefficient list, where no factor
    ## applied later in the radial loop can restore it.  The collapse to a single entry is therefore the whole
    ## of what makes the retardation frequency-dependent; the angular algebra is untouched and needs no extra
    ## prefactor, since the omega -> 0 value above is reproduced by coefficient xcc and nothing else.
    if  rem(la+lc+L-1,2) == 1   &&   rem(lb+ld+L+1,2) == 1
        # mu = 1
        wb =  1 / ( (L+L+1)*(L+L+1) )
        xcc = xc * wb * (c.subshell.kappa - a.subshell.kappa + L) * (d.subshell.kappa - b.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   b, a, d, c,   xcc) )   end
        # mu = 2
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   a, b, c, d,   xcc) )   end
        # mu = 3
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa - L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   d, c, b, a,   xcc) )   end
        # mu = 4
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   c, d, a, b,   xcc) )   end
        # mu = 5
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa + L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   d, a, b, c,   xcc) )   end
        # mu = 6
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   a, d, c, b,   xcc) )   end
        # mu = 7
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L - 1) * (c.subshell.kappa - a.subshell.kappa - L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   b, c, d, a,   xcc) )   end
        # mu = 8
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L,   c, b, a, d,   xcc) )   end
    end

    return( xcList )
end


"""
`InteractionStrength.besselPhiPsi(K::Int64, x::Float64; nmax::Int64=12, tol::Float64=1.0e-16)`
    ... returns the pair (phi_K(x), psi_K(x)) of NORMALISED spherical Bessel functions

            phi_K(x) = (2K+1)!! j_K(x) / x^K            psi_K(x) = - x^(K+1) y_K(x) / (2K-1)!!

        both of which tend to 1 as x -> 0.  They are the natural building blocks of the frequency-dependent
        Breit kernels, because those kernels are the static kernel TIMES a product of one phi and one psi:
        the omega -> 0 limit is then exact by construction rather than a special case, and nothing has to be
        recovered from a cancellation between large numbers.

        Evaluated from the ascending power series in -x^2/2 rather than from j_K and y_K separately.  For
        small x -- which is the physical case, since x = omega r ~ (Delta E) r / c -- forming j_K y_K
        directly loses the answer to cancellation, and each factor alone diverges as omega -> 0.

        Verified against (2K+1)!! j_K(x)/x^K and -x^(K+1) y_K(x)/(2K-1)!! from GSL to <= 4.6e-16 for
        K = 0..3 and x = 1e-4..0.5, and to 1.000000000000000 in the limit x -> 0.

        The normalisation follows GRASP2018 (rci90/bessel.f90, which stores phi-1 and psi-1); the
        formulation was read there and re-derived here rather than transcribed.
"""
function besselPhiPsi(K::Int64, x::Float64; nmax::Int64=12, tol::Float64=1.0e-16)
    x == 0.  &&  return( (1.0, 1.0) )
    wa = -0.5 * x^2;    t1 = 1.0;   t2 = 1.0;   s1 = 0.0;   s2 = 0.0
    for  i = 1:nmax
        t1 = t1 * wa / ( i * (2*(K+i) + 1) );    t2 = t2 * wa / ( i * (2*(i-K) - 1) )
        s1 = s1 + t1;                            s2 = s2 + t2
        if  abs(t1) < abs(s1)*tol  &&  abs(t2) < abs(s2)*tol    break    end
    end
    return( (1.0 + s1, 1.0 + s2) )
end


"""
`InteractionStrength.XL_Breit_densities(xcList::Array{XLCoefficient,1}, factor::Float64, grid::Radial.Grid; tau::Float64=0.)`
    ... computes the the effective Breit interaction strengths X^L_Breit (abcd) for given rank L and a list of
        orbital functions a, b, c, d and angular coefficients at the given grid. A value::Float64 is returned.
        The argument factor scales the photon wave number omega = factor * |E_a - E_c| / c, so that factor = 0
        gives the frequency-independent interaction as the exact limit of the same expressions; see the
        reference formulation heading this section.
        A positive tau damps the radial integrand by exp(-tau r_1 - tau r_2), as the DampedSpaceCI Green
        function approach requires; tau = 0 is the undamped interaction and costs nothing extra.
"""
function XL_Breit_densities(xcList::Array{XLCoefficient,1}, factor::Float64, grid::Radial.Grid; tau::Float64=0.)
    ## V is the MULTIPLIER on the static kernel r_<^nu / r_>^(nu+1), so it must equal 1 at omega = 0.
    ##
    ## REWRITTEN 13-Aug-2026.  It read  -(2nu+1) j_nu(omega r_<) y_nu(omega r_>), which is the frequency-
    ## dependent kernel DIVIDED BY THE STATIC ONE AND BY omega -- it diverged as 1/omega and could never be
    ## switched on (measured: V/static = 1000 at omega = 1e-3, V*omega/static = 1.000001).  Two bare catch
    ## blocks then substituted wx = 1.0 whenever the Bessel evaluation failed, which would have silently
    ## replaced the kernel by its static limit.
    ##
    ## The identity that fixes it, verified numerically to 3e-15 over nu = 1..3 and omega = 0.005..0.5:
    ##
    ##     -omega (2nu+1) j_nu(omega r_<) y_nu(omega r_>)  =  phi_nu(omega r_<) psi_nu(omega r_>) *
    ##                                                        [ r_<^nu / r_>^(nu+1) ]
    ##
    ## with the NORMALISED functions phi and psi of Nuclear-style small-argument form
    ##
    ##     phi_K(x) = (2K+1)!! j_K(x) / x^K        psi_K(x) = -x^(K+1) y_K(x) / (2K-1)!!
    ##
    ## both of which tend to 1 as x -> 0.  So the multiplier is simply phi*psi, the static limit is exact by
    ## construction rather than by a special case, and no cancellation is taken between large numbers.
    ##
    ## This normalisation is the one used by GRASP2018 (rci90/bessel.f90, which stores phi-1 and psi-1); the
    ## formulation was read there and re-derived here, not transcribed.
    function V(nu::Int64, r::Float64, s::Float64, omega::Float64)
        if  omega <= 0.  ||  nu < 0     return( 1.0 )    end
        rSmall = min(r,s);    rLarge = max(r,s)
        phi, _ = InteractionStrength.besselPhiPsi(nu, omega*rSmall)
        _, psi = InteractionStrength.besselPhiPsi(nu, omega*rLarge)
        return( phi * psi )
    end
    #
    ## W is Grant & Pyper's RETARDATION kernel.  Their paper carries two objects both written W, and it matters
    ## which one is meant here:
    ##
    ##   * equation (4)  W_nu = ( V_nu - U_nu ) / omega^2  -- an auxiliary of the general equation (5); and
    ##   * equation (6)  W_(nu-1,nu+1,nu)(1,2)             -- the one that actually appears in the B^(2) term
    ##                                                        of their equation (8), and hence the one needed.
    ##
    ## It is equation (6) that this closure implements:
    ##
    ##     W_(nu-1,nu+1,nu)(1,2) = [nu] omega j_(nu-1)(omega r_1) n_(nu+1)(omega r_2)
    ##                             + ([nu]^2/omega^2) r_1^(nu-1) / r_2^(nu+2)              for r_1 < r_2
    ##                           = [nu] omega n_(nu-1)(omega r_1) j_(nu+1)(omega r_2)      for r_1 > r_2
    ##
    ## with [nu] = 2nu+1, and W_(nu+1,nu-1,nu)(1,2) = W_(nu-1,nu+1,nu)(2,1).  The two orderings are DIFFERENT
    ## functions, not one function of r_< and r_>, which is why equation (8) carries both of them.
    ##
    ## For r_1 < r_2 this is a DIFFERENCE OF TWO 1/omega^2-DIVERGENT QUANTITIES, which is the whole numerical
    ## difficulty.  Substituting the normalised phi, psi of besselPhiPsi,
    ##     j_(nu-1)(x) = phi_(nu-1)(x) x^(nu-1)/(2nu-1)!!      n_(nu+1)(x) = -psi_(nu+1)(x) (2nu+1)!!/x^(nu+2)
    ## the first term becomes exactly -([nu]^2/omega^2) phi psi r_1^(nu-1)/r_2^(nu+2), so the divergence
    ## cancels ALGEBRAICALLY against the second and what remains is
    ##
    ##     W = -([nu]^2/omega^2) (r_1^(nu-1)/r_2^(nu+2)) (a + b + a*b),   a = phi_(nu-1)-1,  b = psi_(nu+1)-1
    ##
    ## in which nothing is subtracted numerically: a and b come straight from their own power series.
    ## For r_1 > r_2 no cancellation arises and the same substitution gives an O(omega^2) expression.
    ##
    ## VERIFIED against Grant's analytic omega -> 0 limit, their equations (9)-(10),
    ##     W_(nu-1,nu+1,nu) -> -(1/2) [nu] ( Ubar_(nu-1) - Ubar_(nu+1) ),   Ubar_k = r_1^k/r_2^(k+1), r_1 < r_2
    ##                      -> 0                                                                     r_1 > r_2
    ## the deviation scaling as O(omega^2): 4.9e-3, 4.9e-5, 4.9e-7 at omega = 0.1, 0.01, 0.001 for nu = 1, and
    ## likewise for nu = 2, 3.  That limit is precisely what the 'S' coefficients of XL_Breit_coefficients
    ## already encode in their +-[nu]/2 weights, so at omega = 0 kernel and coefficients agree by construction.
    ##
    ## Below omegaFloor the residual roundoff of a and b, amplified by 1/omega^2, exceeds that deviation, so
    ## the analytic limit is returned instead -- the device of GRASP2018's IF (W < EPSI**2) in rci90/bessel.f90.
    ##
    ## HISTORY.  What stood here before 13-Aug-2026 was equation (6) with the factor omega missing from its
    ## first term and that term's sign reversed, so the cancellation failed outright; measured, neither W,
    ## W*omega nor W*omega^2 converged, and two bare catch blocks then substituted wx = 1.0.  The rewrite of
    ## bd0e9fb replaced it with equation (4) instead, which is a different function -- corrected here.
    function W(nu::Int64, r1::Float64, r2::Float64, omega::Float64)
        nu < 1  &&  return( 0.0 )
        nn = 2nu + 1
        omegaFloor = 1.0e-4
        if  r1 < r2
            if  omega <= omegaFloor
                return( -0.5 * nn * ( r1^(nu-1) / r2^nu - r1^(nu+1) / r2^(nu+2) ) )
            end
            aa = InteractionStrength.besselPhiPsi(nu-1, omega*r1)[1] - 1.0
            bb = InteractionStrength.besselPhiPsi(nu+1, omega*r2)[2] - 1.0
            return( -(nn^2 / omega^2) * (r1^(nu-1) / r2^(nu+2)) * (aa + bb + aa*bb) )
        else
            omega <= omegaFloor  &&  return( 0.0 )
            psi = InteractionStrength.besselPhiPsi(nu-1, omega*r1)[2]
            phi = InteractionStrength.besselPhiPsi(nu+1, omega*r2)[1]
            return( -omega^2 * psi * phi * r2^(nu+1) / ((2nu-1) * (2nu+3) * r1^nu) )
        end
    end
    
    wa = 0.
    for  xc  in  xcList  ## [end:end]
        # Use the minimal extent of any involved orbitals; this need to be improved
        mtp_ac = min(size(xc.a.P, 1), size(xc.c.P, 1));     mtp_bd = min(size(xc.b.P, 1), size(xc.d.P, 1))
        ## omega is a WAVE NUMBER, not an energy: omega = |E_a - E_c| / c in atomic units, so that omega*r is
        ## the dimensionless phase the spherical Bessel functions require.  The division by c was MISSING here
        ## until 13-Aug-2026, which overstated the phase by a factor c = 137 and would have sampled j_nu*y_nu
        ## in their oscillatory regime instead of the small-argument one where retardation actually lives.
        ## It never affected a published number, because no branch below ever consumed these values.
        cLight = Defaults.getDefaults("speed of light: c")
        omg_ac = factor * abs(xc.a.energy - xc.c.energy) / cLight
        omg_bd = factor * abs(xc.b.energy - xc.d.energy) / cLight
        ## Only s <= r contributes: every permutation is emitted together with its TRANSPOSE, and the transposed
        ## partner covers the other radial ordering.  Grant's coordinate 1 is then always the SMALLER radius --
        ## his Ubar_nu(1,2) vanishes for r_1 > r_2, and the static factor here is r_s^nu / r_r^(nu+1) -- so both
        ## kernels are called with grid.r[s] first.  The two frequencies omg_ac, omg_bd coincide when energy is
        ## conserved; their mean is used where it does not, as GRASP2018 likewise works with a single omega.
        for  r = 2:mtp_ac
            for  s = 2:min(r, mtp_bd)
                if      xc.kind == 'T'
                    ## Gaunt (magnetic) part: the static kernel Ubar_nu times the frequency factor V_nu -> 1.
                    wk = ( V(xc.nu, grid.r[s], grid.r[r], omg_ac) + V(xc.nu, grid.r[s], grid.r[r], omg_bd) ) / 2.0 *
                         (grid.r[s]^xc.nu) / (grid.r[r]^(xc.nu+1))
                elseif  xc.kind == 'S'
                    ## Retardation: W_(L-1,L+1,L) is the COMPLETE kernel, not a factor multiplying Ubar_nu.
                    wk = ( W(xc.nu, grid.r[s], grid.r[r], omg_ac) + W(xc.nu, grid.r[s], grid.r[r], omg_bd) ) / 2.0
                else    error("stop a")
                end
                #
                wc = xc.coeff * grid.wr[r] * grid.wr[s]
                if  s == r      wc = wc / 2.0                                          end
                if  tau > 0.    wc = wc * exp( -tau*grid.r[r] - tau*grid.r[s] )         end
                wa = wa + wc * wk * (xc.a.P[r] * xc.c.Q[r]) * (xc.b.P[s] * xc.d.Q[s])
            end
        end
    end
    return( wa )
end


"""
`InteractionStrength.XL_CoulombReference(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
        a, b, c and d at the given grid but without optimization. A value::Float64 is returned.
"""
function XL_CoulombReference(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        return( 0. )
    end
    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end 
    
    XL_Coulomb = xc * RadialIntegrals.SlaterRkReference(L, a, b, c, d, grid)
    return( XL_Coulomb )
end


"""
`InteractionStrength.XL_Coulomb_reset_storage(keep::Bool; printout::Bool=false)`  
    ... resets the global storage of XL_Coulomb interaction strength; nothing is returned.
"""
function XL_Coulomb_reset_storage(keep::Bool; printout::Bool=false)
    if  keep
        if printout     println(">> Reset GBL_Storage_XL_Coulomb storage.")     end
        global GBL_Storage_XL_Coulomb = Dict{String, Float64}()
    else
    end
    return( nothing )      
end


"""
`InteractionStrength.XL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)`  
    ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
        a, b, c and d at the given grid. For keep=true, the procedure looks up the (global) directory GBL_Storage_XL_Coulomb
        and returns the corresponding value without re-calculation of the interaction strength; it also 'stores' the calculated
        value if not yet included. For keep=false, the interaction strength is always computed on-fly. A value::Float64 is 
        returned.
"""
function XL_Coulomb(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)
    global GBL_Storage_XL_Coulomb
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        return( 0. )
    end
    
    # Now distiguish due to the optional argument keep
    if  keep
        sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
        if haskey(GBL_Storage_XL_Coulomb, sa )
            XL_Coulomb = GBL_Storage_XL_Coulomb[sa]
        else
            XL_Coulomb = InteractionStrength.XL_Coulomb(L::Int64, a, b, c, d, grid)
            ## global GBL_Storage_XL_Coulomb = Base.merge(GBL_Storage_XL_Coulomb, Dict( sa => XL_Coulomb))
            global GBL_Storage_XL_Coulomb[sa] = XL_Coulomb
        end
    else
        xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
        if   rem(L,2) == 1    xc = - xc    end
        ## NOT SWITCHED to the kink-aware RadialIntegrals.SlaterRkKinkAware -- attempted 13-Aug-2026 and
        ## REVERTED, with what was measured recorded here so the attempt is not simply repeated.
        ##
        ## The kink-aware integral IS the better quadrature, and that part is settled: against the analytic
        ## F^0(1s,1s) = 5Z/8 it is converged already on the coarsest grid tried (identical to nine digits over
        ## a 10x refinement), whereas this line's rule still drifts; on JAC's DEFAULT exponential grid the
        ## errors are 2.59e-4 against 7.85e-5.  The two agree to 1e-5..2e-4 on direct AND cross terms, both
        ## satisfy R^k(abcd) = R^k(badc), and the discrepancy grows with rank as the r_< / r_> cusp sharpens.
        ##
        ## WHAT STOPPED THE SWITCH is that its effect on the approved references could not be interpreted.
        ## Median changes are 1e-5..1e-3, as a 1e-4 shift in the integrals should give -- but individual
        ## entries move by factors of 10 to 250, and it could NOT be established whether those are real
        ## near-cancellations or merely rows changing places, since test-Cascade-StepwiseDecay is known to
        ## reorder degenerate levels and a line-by-line file comparison cannot tell the two apart.
        ##
        ## Re-approving twelve references on evidence that cannot be read would be the opposite of a
        ## deliberate editorial act.  The switch needs a comparison that matches transitions by their QUANTUM
        ## NUMBERS rather than by line position; until that exists, this stays as it is.
        XL_Coulomb = xc * RadialIntegrals.SlaterRk(L, a, b, c, d, grid)

    end

    return( XL_Coulomb )
end


"""
`InteractionStrength.XL_CoulombKinkAware_reset_storage(keep::Bool; printout::Bool=false)`
    ... resets the global storage of XL_CoulombKinkAware interaction strengths (a SEPARATE cache from
        XL_Coulomb's own GBL_Storage_XL_Coulomb, since the kink-aware and standard quadratures give
        different numeric results for the same subshell labels and must not share a cache namespace);
        nothing is returned.
"""
function XL_CoulombKinkAware_reset_storage(keep::Bool; printout::Bool=false)
    if  keep
        if printout     println(">> Reset GBL_Storage_XL_CoulombKinkAware storage.")     end
        global GBL_Storage_XL_CoulombKinkAware = Dict{String, Float64}()
    else
    end
    return( nothing )
end


"""
`InteractionStrength.XL_CoulombKinkAware(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)`
    ... computes the same effective Coulomb interaction strength as XL_Coulomb(L, a, b, c, d, grid), including the
        same triangular-delta veto and angular reduced-matrix-element prefactor xc, but using the kink-aware
        RadialIntegrals.SlaterRkKinkAware for the underlying radial integral instead of RadialIntegrals.
        SlaterRk. For keep=true, looks up (and stores into) the global GBL_Storage_XL_CoulombKinkAware
        Dict, mirroring XL_Coulomb's own keep/GBL_Storage_XL_Coulomb pattern exactly -- used by
        Hamiltonian.setupMatrixKinkAware (the CI-matrix Coulomb term for ALField/EOLField), where orbitals
        are FIXED for the whole performCIKinkAware call, so caching is unconditionally safe there. NOT enabled
        (keep=false, the default) at SelfConsistent.computeTwoElectronV's own call site: that call sits
        inside the outer SCF iteration, where orbitals change every iteration, so a cache surviving across
        iterations would silently return stale integrals from an earlier orbital shape -- extending caching
        safely into that loop needs its own explicit per-iteration reset wiring, deferred as a separate item.
        Isolated from XL_Coulomb; shared by Hamiltonian.setupMatrixKinkAware and
        SelfConsistent.computeTwoElectronV (their Fock matrix, uncached). A value::Float64 is
        returned.
"""
function XL_CoulombKinkAware(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)
    global GBL_Storage_XL_CoulombKinkAware
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        return( 0. )
    end

    if  keep
        sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
        if haskey(GBL_Storage_XL_CoulombKinkAware, sa)
            return( GBL_Storage_XL_CoulombKinkAware[sa] )
        end
    end

    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end

    XL_CoulombKinkAwareValue = xc * RadialIntegrals.SlaterRkKinkAware(L, a, b, c, d, grid)

    if  keep
        sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
        global GBL_Storage_XL_CoulombKinkAware[sa] = XL_CoulombKinkAwareValue
    end

    return( XL_CoulombKinkAwareValue )
end


"""
`InteractionStrength.XL_Coulomb(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)`  
    ... computes the (direct) Coulomb interaction strengths X^L_Coulomb (.b.d) for given rank L and orbital functions
        as well as the given primitives. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_Coulomb(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wm  = zeros(nsL+nsS, nsL+nsS)
    
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a);             ja2 = Basics.subshell_2j(a)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c);             jc2 = Basics.subshell_2j(c)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        @warn("stop aa")  ## This should not occur.
        return( wm )
    end
    
    xc = AngularMomentum.CL_reduced_me(a, L, c) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end 
    
    # Direct interaction; contract the full interaction array over the orbitals b and d
    wm = zeros(nsL+nsS, nsL+nsS)
    for  i = 1:nsL
        for  k = 1:nsL 
            ## Ba = primitives.bsplinesL[i].bs;    Bc = primitives.bsplinesL[k].bs
            Ba = primitives.bsplinesL[i];    Bc = primitives.bsplinesL[k]
            Pa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Pa[j] = Pa[j] + Ba.bs[j+add]   end
            Pc = zeros(Bc.upper);   add = 1 - Bc.lower;   
            for  j = Bc.lower:Bc.upper  Pc[j] = Pc[j] + Bc.bs[j+add]   end
            wm[i,k] = RadialIntegrals.SlaterRkComponent(L, Pa, b.P, Pc, d.P, grid) + 
                      RadialIntegrals.SlaterRkComponent(L, Pa, b.Q, Pc, d.Q, grid)
        end
    end
    for  i = 1:nsS
        for  k = 1:nsS
            ## Ba = primitives.bsplinesS[i].bs;    Bc = primitives.bsplinesS[k].bs
            Ba = primitives.bsplinesS[i];    Bc = primitives.bsplinesS[k]
            Qa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Qa[j] = Qa[j] + Ba.bs[j+add]   end
            Qc = zeros(Bc.upper);   add = 1 - Bc.lower;   
            for  j = Bc.lower:Bc.upper  Qc[j] = Qc[j] + Bc.bs[j+add]   end            
            wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent(L, Qa, b.P, Qc, d.P, grid) + 
                              RadialIntegrals.SlaterRkComponent(L, Qa, b.Q, Qc, d.Q, grid)
        end
    end
    
    return( xc * wm )
end


"""
`InteractionStrength.XL_CoulombKinkAware(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)`
    ... computes the same (direct) Coulomb interaction strengths X^L_Coulomb (.b.d) as
        XL_Coulomb(L,a::Subshell,b::Orbital,c::Subshell,d::Orbital,primitives), for given rank L and orbital
        functions as well as the given primitives, but using the kink-aware screened-potential construction
        (RadialIntegrals.buildScreenedPotential) instead of RadialIntegrals.SlaterRkComponent's naive
        tensor-product double sum. The screened potential V_L(r), which depends only on the fixed orbital pair
        (b,d), is built ONCE (adaptive quadrature) and then reused cheaply -- via the existing grid quadrature
        weights, since V_L(r) is smooth once built -- for every B-spline pair (i,k) of the L- and S-block.
        Isolated from XL_Coulomb; shared by the ALField/EOLField code lines, cf.
        SelfConsistent.computeTwoElectronV. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_CoulombKinkAware(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wm  = zeros(nsL+nsS, nsL+nsS)

    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a);             ja2 = Basics.subshell_2j(a)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c);             jc2 = Basics.subshell_2j(c)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        @warn("stop aa")  ## This should not occur.
        return( wm )
    end

    xc = AngularMomentum.CL_reduced_me(a, L, c) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end

    # Build the screened potential once for the fixed orbital pair (b,d); reused below for both L- and S-block.
    # mtpOut is forced to the full grid extent: unlike SlaterRkKinkAware's orbital-orbital use (where the
    # OTHER factor in the contraction is also naturally truncated to some orbital's own extent), here Vk gets
    # contracted against B-spline ROW/COLUMN indices that span the FULL basis and can extend well past (b,d)'s
    # own reach -- leaving mtpOut at its default silently drops that tail and was traced to a real bug (Ne's 1s
    # orbital coming out measurably too deeply bound from a missing part of its screening by the 2p shell).
    Vk = RadialIntegrals.buildScreenedPotential(L, b, d, grid; mtpOut=grid.NoPoints)

    # Direct interaction; contract Vk against every B-spline pair of the L- and S-block
    wm = zeros(nsL+nsS, nsL+nsS)
    for  i = 1:nsL
        for  k = 1:nsL
            Ba = primitives.bsplinesL[i];    Bc = primitives.bsplinesL[k]
            Pa = zeros(Ba.upper);   add = 1 - Ba.lower;
            for  j = Ba.lower:Ba.upper  Pa[j] = Pa[j] + Ba.bs[j+add]   end
            Pc = zeros(Bc.upper);   add = 1 - Bc.lower;
            for  j = Bc.lower:Bc.upper  Pc[j] = Pc[j] + Bc.bs[j+add]   end
            mtp = min(Ba.upper, Bc.upper, length(Vk))
            wa  = 0.
            for  r = 2:mtp   wa = wa + Pa[r]*Pc[r] * grid.wr[r] * Vk[r]   end
            wm[i,k] = wa
        end
    end
    for  i = 1:nsS
        for  k = 1:nsS
            Ba = primitives.bsplinesS[i];    Bc = primitives.bsplinesS[k]
            Qa = zeros(Ba.upper);   add = 1 - Ba.lower;
            for  j = Ba.lower:Ba.upper  Qa[j] = Qa[j] + Ba.bs[j+add]   end
            Qc = zeros(Bc.upper);   add = 1 - Bc.lower;
            for  j = Bc.lower:Bc.upper  Qc[j] = Qc[j] + Bc.bs[j+add]   end
            mtp = min(Ba.upper, Bc.upper, length(Vk))
            wa  = 0.
            for  r = 2:mtp   wa = wa + Qa[r]*Qc[r] * grid.wr[r] * Vk[r]   end
            wm[nsL+i,nsL+k] = wa
        end
    end

    return( xc * wm )
end


"""
`InteractionStrength.XL_Coulomb(L::Int64, a::Subshell, b::Orbital, c::Orbital, d::Subshell, primitives::Bsplines.Primitives)`
    ... computes the (exchange) Coulomb interaction strengths X^L_Coulomb (.bc.) for given rank L and orbital functions
        as well as the given primitives. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_Coulomb(L::Int64, a::Subshell, b::Orbital, c::Orbital, d::Subshell, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wm  = zeros(nsL+nsS, nsL+nsS)
    
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a);             ja2 = Basics.subshell_2j(a)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d);             jd2 = Basics.subshell_2j(d)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        @warn("stop ab")  ## This should not occur.
        return( wm )
    end
    
    xc = AngularMomentum.CL_reduced_me(a, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d)
    if   rem(L,2) == 1    xc = - xc    end 
    
    # Exchange interaction; contract the full interaction array over the orbitals b and c
    for  i = 1:nsL
        for  k = 1:nsL 
            ## Ba = primitives.bsplinesL[i].bs;    Bd = primitives.bsplinesL[k].bs
            Ba = primitives.bsplinesL[i];    Bd = primitives.bsplinesL[k]
            Pa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Pa[j] = Pa[j] + Ba.bs[j+add]   end
            Pd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Pd[j] = Pd[j] + Bd.bs[j+add]   end            
            wm[i,k] = RadialIntegrals.SlaterRkComponent(L, Pa, b.P, c.P, Pd, grid)
        end
        for  k = 1:nsS 
            ## Ba = primitives.bsplinesL[i].bs;    Bd = primitives.bsplinesS[k].bs
            Ba = primitives.bsplinesL[i];    Bd = primitives.bsplinesS[k]
            Pa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Pa[j] = Pa[j] + Ba.bs[j+add]   end
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Qd[j] = Qd[j] + Bd.bs[j+add]   end            
            wm[i,nsL+k] = RadialIntegrals.SlaterRkComponent(L, Pa, b.P, c.Q, Qd, grid)
        end
    end
    for  i = 1:nsS
        for  k = 1:nsL 
            ## Ba = primitives.bsplinesS[i].bs;    Bd = primitives.bsplinesL[k].bs
            Ba = primitives.bsplinesS[i];    Bd = primitives.bsplinesL[k]
            Qa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Qa[j] = Qa[j] + Ba.bs[j+add]   end
            Pd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Pd[j] = Pd[j] + Bd.bs[j+add]   end            
            wm[nsL+i,k] = RadialIntegrals.SlaterRkComponent(L, Qa, b.Q, c.P, Pd, grid)
        end
        for  k = 1:nsS 
            ## Ba = primitives.bsplinesS[i].bs;    Bd = primitives.bsplinesS[k].bs
            Ba = primitives.bsplinesS[i];    Bd = primitives.bsplinesS[k]
            Qa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Qa[j] = Qa[j] + Ba.bs[j+add]   end
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Qd[j] = Qd[j] + Bd.bs[j+add]   end            
            wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent(L, Qa, b.Q, c.Q, Qd, grid)
        end
    end

    return( xc * wm )
end


"""
`InteractionStrength.XL_CoulombTensor(L::Int64, a::Subshell, b::Orbital, c::Orbital, cVector::Vector{Float64},
                                            d::Subshell, cacheLL::RadialIntegrals.ScreenedPotentialCache,
                                            cacheLS::RadialIntegrals.ScreenedPotentialCache,
                                            cacheSS::RadialIntegrals.ScreenedPotentialCache,
                                            primitives::Bsplines.Primitives)`
    ... computes the same (exchange) Coulomb interaction strengths X^L_Coulomb (.bc.) as
        XL_Coulomb(L,a::Subshell,b::Orbital,c::Orbital,d::Subshell,primitives) / XL_CoulombKinkAware of the same
        signature, but using PRECOMPUTED RadialIntegrals.ScreenedPotentialCache tensors -- built ONCE, outside
        the SCF iteration, via RadialIntegrals.buildScreenedPotentialCache for rank L -- instead of any per-call
        adaptive quadrature. Re-deriving the original (validated) block structure carefully shows that, in each
        block, the B-spline row index pairs with orbital c's component matching the COLUMN's type (P for an
        L-column, Q for an S-column), while orbital b's component matching the ROW's type pairs with the column
        B-spline index. Consequently it is orbital c -- not b -- whose OWN B-spline expansion coefficients
        (cVector, the B-spline coefficient vector carried alongside the orbital during SCF) are needed, to build, for
        each row i, a "partial potential" Psi_i(s) = sum_m cVector[m] * Phi_(i,m)(s) as a CHEAP weighted sum of
        already-cached, already-kink-aware Phi_(i,m) entries (only ~band terms, no quadrature at all); orbital b
        itself stays a plain tabulated orbital, exactly like the column B-spline, entering only the final cheap
        grid-quadrature dot product. Three caches are needed because the row/column-type combination can be
        LL, LS, SL, or SS -- and, unlike the "direct" signature, a single screened potential cannot be shared
        across the whole (nsL+nsS) x (nsL+nsS) matrix here, since both B-spline indices (row and column) sit on
        opposite sides of the two-electron kernel. cacheLL and cacheSS are ordinary same-basis caches
        (bsplinesL, bsplinesL) and (bsplinesS, bsplinesS); cacheLS is the cross-basis cache
        (bsplinesL, bsplinesS) and is also used, with indices swapped at lookup, for the SL combination -- see
        RadialIntegrals.buildScreenedPotentialCache. This turns what used to be an expensive per-call
        adaptive-quadrature computation into one with NO further quadrature at all, for every SCF iteration and
        every exchange coefficient that shares this rank L, once the caches themselves have been built (see
        SelfConsistent.solveAverageLevelField for how they get built and cached once per SCF run).
        Isolated from XL_Coulomb; only used by the average-level (ALField) code line. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_CoulombTensor(L::Int64, a::Subshell, b::Orbital, c::Orbital, cVector::Vector{Float64}, d::Subshell,
                                cacheLL::RadialIntegrals.ScreenedPotentialCache,
                                cacheLS::RadialIntegrals.ScreenedPotentialCache,
                                cacheSS::RadialIntegrals.ScreenedPotentialCache,
                                primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;  nsS = primitives.grid.nsS;  grid = primitives.grid
    wm  = zeros(nsL+nsS, nsL+nsS)

    la = Basics.subshell_l(a);             ja2 = Basics.subshell_2j(a)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d);             jd2 = Basics.subshell_2j(d)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        @warn("stop ab")  ## This should not occur.
        return( wm )
    end

    xc = AngularMomentum.CL_reduced_me(a, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d)
    if   rem(L,2) == 1    xc = - xc    end

    if  length(cVector) != nsL+nsS    error("stop a; cVector must hold both L- and S-basis expansion coefficients")   end
    cP = cVector[1:nsL];   cQ = cVector[nsL+1:nsL+nsS]

    # L-block rows: row B-spline "i" from bsplinesL
    for  i = 1:nsL
        PsiLL = zeros(grid.NoPoints)   # for L-columns (contracted below against c.P's expansion, cacheLL)
        for  m = max(1,i-cacheLL.band):min(nsL,i+cacheLL.band)
            Phi_im = get(cacheLL.Phi, (i,m), nothing)
            if  Phi_im === nothing || cP[m] == 0.   continue   end
            mtp = min(length(Phi_im), length(PsiLL))
            for  r = 1:mtp   PsiLL[r] += cP[m] * Phi_im[r]   end
        end
        PsiLS = zeros(grid.NoPoints)   # for S-columns (contracted below against c.Q's expansion, cacheLS)
        for  n = 1:nsS
            Phi_in = get(cacheLS.Phi, (i,n), nothing)
            if  Phi_in === nothing || cQ[n] == 0.   continue   end
            mtp = min(length(Phi_in), length(PsiLS))
            for  r = 1:mtp   PsiLS[r] += cQ[n] * Phi_in[r]   end
        end

        for  k = 1:nsL
            Bd = primitives.bsplinesL[k]
            Pd = zeros(Bd.upper);   add = 1 - Bd.lower
            for  j = Bd.lower:Bd.upper   Pd[j] = Pd[j] + Bd.bs[j+add]   end
            mtp = min(length(PsiLL), length(b.P), length(Pd))
            wa  = 0.
            for  s = 2:mtp   wa += PsiLL[s] * b.P[s] * Pd[s] * grid.wr[s]   end
            wm[i,k] = wa
        end
        for  k = 1:nsS
            Bd = primitives.bsplinesS[k]
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower
            for  j = Bd.lower:Bd.upper   Qd[j] = Qd[j] + Bd.bs[j+add]   end
            mtp = min(length(PsiLS), length(b.P), length(Qd))
            wa  = 0.
            for  s = 2:mtp   wa += PsiLS[s] * b.P[s] * Qd[s] * grid.wr[s]   end
            wm[i,nsL+k] = wa
        end
    end

    # S-block rows: row B-spline "i" from bsplinesS
    for  i = 1:nsS
        PsiSL = zeros(grid.NoPoints)   # for L-columns; cacheLS was built as (bsplinesL, bsplinesS), so the
        for  m = 1:nsL                 # row (S-basis) index is now the SECOND slot: look up (m,i), not (i,m)
            Phi_mi = get(cacheLS.Phi, (m,i), nothing)
            if  Phi_mi === nothing || cP[m] == 0.   continue   end
            mtp = min(length(Phi_mi), length(PsiSL))
            for  r = 1:mtp   PsiSL[r] += cP[m] * Phi_mi[r]   end
        end
        PsiSS = zeros(grid.NoPoints)   # for S-columns (cacheSS)
        for  n = max(1,i-cacheSS.band):min(nsS,i+cacheSS.band)
            Phi_in = get(cacheSS.Phi, (i,n), nothing)
            if  Phi_in === nothing || cQ[n] == 0.   continue   end
            mtp = min(length(Phi_in), length(PsiSS))
            for  r = 1:mtp   PsiSS[r] += cQ[n] * Phi_in[r]   end
        end

        for  k = 1:nsL
            Bd = primitives.bsplinesL[k]
            Pd = zeros(Bd.upper);   add = 1 - Bd.lower
            for  j = Bd.lower:Bd.upper   Pd[j] = Pd[j] + Bd.bs[j+add]   end
            mtp = min(length(PsiSL), length(b.Q), length(Pd))
            wa  = 0.
            for  s = 2:mtp   wa += PsiSL[s] * b.Q[s] * Pd[s] * grid.wr[s]   end
            wm[nsL+i,k] = wa
        end
        for  k = 1:nsS
            Bd = primitives.bsplinesS[k]
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower
            for  j = Bd.lower:Bd.upper   Qd[j] = Qd[j] + Bd.bs[j+add]   end
            mtp = min(length(PsiSS), length(b.Q), length(Qd))
            wa  = 0.
            for  s = 2:mtp   wa += PsiSS[s] * b.Q[s] * Qd[s] * grid.wr[s]   end
            wm[nsL+i,nsL+k] = wa
        end
    end

    return( xc * wm )
end


"""
`InteractionStrength.XL_CoulombDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`
    ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
        a, b, c and d at the given grid. A value::Float64 is returned.
"""
function XL_CoulombDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        return( 0. )
    end
    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end 
    
    XL_Coulomb = xc * RadialIntegrals.SlaterRkDamped(tau::Float64, L, a, b, c, d, grid)
    return( XL_Coulomb )
end


"""
`InteractionStrength.XL_Coulomb_DH(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid, lambda::Float64)`  
    ... computes the the effective Coulomb-Debye-Hückel interaction strengths X^L_Coulomb_DH (abcd) for given rank L and 
        orbital functions a, b, c and d at the given grid and for the given screening parameter lambda. A value::Float64 is 
        returned.
"""
function XL_Coulomb_DH(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid, lambda::Float64)
    # Test for the triangular-delta conditions and calculate the reduced matrix elements of the C^L tensors
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
        rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
        return( 0. )
    end
    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end 

    XL_Coulomb_DH = xc * RadialIntegrals.SlaterRkDebyeHueckel(L, a, b, c, d, grid, lambda)
    return( XL_Coulomb_DH )
end


"""
`InteractionStrength.XL_plasma_ionSphere(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, lambda::Float64)`  
    ... computes the effective interaction strengths X^L_ion-sphere (abcd) for given rank L and orbital functions 
        a, b, c and d and for the plasma parameter lambda. A value::Float64 is returned.  **Not yet implemented !**
"""
function XL_plasma_ionSphere(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, lambda::Float64)
    error("Not yet implemented")
end


"""
`InteractionStrength.X_smsA(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the the effective interaction strengths X^1_sms,A (abcd) for fixed rank 1 and orbital functions 
        a, b, c and d at the given grid. A value::Float64 is returned.
"""
function X_smsA(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    wa = AngularMomentum.CL_reduced_me(a.subshell, 1, c.subshell) * 
            AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) *
            RadialIntegrals.Vinti(a, c, grid) * RadialIntegrals.Vinti(b, d, grid) / 2
    ## println("**  <$(a.subshell) || Vinti || $(c.subshell)>  = $(RadialIntegrals.Vinti(a, c, grid)) " )
    ## println("**  <$(b.subshell) || Vinti || $(d.subshell)>  = $(RadialIntegrals.Vinti(b, d, grid)) " )
    return( wa )
end


"""
`InteractionStrength.X_smsB(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... computes the the effective interaction strengths X^1_sms,B (abcd) for fixed rank 1 and orbital functions 
        a, b, c and d at the given grid. A value::Float64 is returned.
"""
function X_smsB(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    wa = - AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) * RadialIntegrals.Vinti(b, d, grid) *
            RadialIntegrals.isotope_smsB(a, c, nm.Z, grid) / 2
    return( wa )
end


"""
`InteractionStrength.X_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)` 
    ... computes the the effective interaction strengths X^1_sms,C (abcd) for fixed rank 1 and orbital functions 
        a, b, c and d at the given grid. A value::Float64 is returned.
"""
function X_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    wa = - AngularMomentum.CL_reduced_me(b.subshell, 1, d.subshell) * 
            AngularMomentum.CL_reduced_me(a.subshell, 1, c.subshell) * 
            RadialIntegrals.Vinti(b, d, grid) * RadialIntegrals.isotope_smsC(a, c, nm.Z, grid) / 2
    return( wa )
end


"""
`InteractionStrength.zeeman_Delta_n1(a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... computes the <a|| Delta n^(1) ||b> reduced matrix element for the Zeeman-Schwinger contribution to the coupling
        to an external magnetic field for orbital functions a, b. A value::Float64 is returned.

        Note (25-Jul-2026): the prefactor is (g_s-2)/4, not the (g_s-2)/2 of Andersson & Jonsson (2008), CPC, Eq. (26)/(52)
        as literally printed -- found and fixed after the printed formula gave a Delta N1 contribution to g_J exactly 2x
        too large, confirmed against two independent references: (i) the standard non-relativistic Lande g_J decomposition
        for H(2p_1/2) and H(2p_3/2) (different kappa, opposite-sign corrections, both matched after the /4 fix to 5-6
        significant figures), and (ii) the paper's own published He 1s2p g_J table (Fig. 10: 1.5011166/1.5011183/0.9999936
        for 3P1/3P2/1P1), matched to 5-6 significant figures only after this fix, not before. Likely root cause (plausible,
        not fully proven): the Sigma operator in Eq. (26) is Pauli-normalized (eigenvalue +/-1), while the physical
        magnetic-moment operator needs S (eigenvalue +/-1/2) -- a factor-of-2 normalization slip, not a JAC transcription
        error, since the code faithfully reproduced what Eq. (52) says.
"""
function zeeman_Delta_n1(a::Orbital, b::Orbital, grid::Radial.Grid)
    ka = a.subshell.kappa
    kb = b.subshell.kappa

    rad = RadialIntegrals.rkDiagonal(0, a.P, b.P, grid) * (ka + kb - 1) + RadialIntegrals.rkDiagonal(0, a.Q, b.Q, grid) * (ka + kb + 1)
    ## the former CL_reduced_me_rb convention: divided by sqrt(2 j + 1) of the FIRST argument, which here is
    ## the SIGN-FLIPPED kappa, i.e. j(-ka) and not j(ka) (10-Aug-2026).
    ang = AngularMomentum.CL_reduced_me(Subshell(1, -ka), 1, b.subshell) /
              sqrt( Basics.subshell_2j(Subshell(1, -ka)) + 1 )

    return ( (Defaults.getDefaults("electron g-factor") - 2)/4 * rad * ang )
end


"""
`InteractionStrength.zeeman_n1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| n^(1) ||b> reduced matrix element for the Zeeman coupling to an external magnetic field for 
        orbital functions a, b. A value::Float64 is returned. 
"""
function zeeman_n1(a::Orbital, b::Orbital, grid::Radial.Grid)
    ka = a.subshell.kappa
    kb = b.subshell.kappa

    rad = RadialIntegrals.rkNonDiagonal(1, a, b, grid)
    ## the former CL_reduced_me_rb convention: divided by sqrt(2 j + 1) of the FIRST argument, which here is
    ## the SIGN-FLIPPED kappa, i.e. j(-ka) and not j(ka) (10-Aug-2026).
    ang = AngularMomentum.CL_reduced_me(Subshell(1, -ka), 1, b.subshell) /
              sqrt( Basics.subshell_2j(Subshell(1, -ka)) + 1 )

    return ( -rad * ang/(2 * Defaults.getDefaults("alpha")) * (ka + kb) )
end


end # module

