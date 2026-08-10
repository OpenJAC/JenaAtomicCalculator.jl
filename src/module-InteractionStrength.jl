
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
    !(grid.meshType == Radial.MeshGL())  &&  error("Only for Radial.MeshGL() implemented so far.")
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
        At present, only the zero-frequency Breit or Gaunt interaction is taken into account.
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
`InteractionStrength.XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the the effective Breit interaction strengths X^L_Breit (abcd) for given rank L and orbital functions 
        a, b, c and d  at the given grid. A value::Float64 is returned. At present, only the zero-frequency Breit 
        interaction is taken into account.
"""
function XL_BreitDamped(tau::Float64, L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
    error("stop a")
end


"""
`InteractionStrength.XL_Breit_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital; onlyGaunt::Bool=false)`  
    ... evaluates the combinations and pre-coefficients for the zero-frequency Breit interaction  
        X^L_Breit (omega=0.; abcd) for given rank L and orbital functions a, b, c and d. A list of coefficients 
        xcList::Array{XLCoefficient,1} is returned.
"""
function XL_Breit_coefficients(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital; onlyGaunt::Bool=false)
    xcList = XLCoefficient[]
    
    la = Basics.subshell_l(a.subshell);    ja2 = Basics.subshell_2j(a.subshell)
    lb = Basics.subshell_l(b.subshell);    jb2 = Basics.subshell_2j(b.subshell)
    lc = Basics.subshell_l(c.subshell);    jc2 = Basics.subshell_2j(c.subshell)
    ld = Basics.subshell_l(d.subshell);    jd2 = Basics.subshell_2j(d.subshell)

    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc                end 
    if   abs(xc)  <  1.0e-10    return( xcList )   end

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
        wa = - (a.subshell.kappa + c.subshell.kappa) * (b.subshell.kappa + d.subshell.kappa) / (L*(L+1))
        # mu = 1
        xcc = xc * wa
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

    # Add contributions of the S^k_mu integrals
    if  rem(la+lc+L-1,2) == 1   &&   rem(lb+ld+L+1,2) == 1
        # mu = 1
        wb =  1 / ( (L+L+1)*(L+L+1) )
        xcc = xc * wb * (c.subshell.kappa - a.subshell.kappa + L) * (d.subshell.kappa - b.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, b, a, d, c,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, b, a, d, c,  -(L+L+1) / 2 * xcc) )   end
        # mu = 2
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, a, b, c, d,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, a, b, c, d,  -(L+L+1) / 2 * xcc) )   end
        # mu = 3
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa - L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, d, c, b, a,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, d, c, b, a,  -(L+L+1) / 2 * xcc) )   end
        # mu = 4
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, c, d, a, b,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, c, d, a, b,  -(L+L+1) / 2 * xcc) )   end
        # mu = 5
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L + 1) * (c.subshell.kappa - a.subshell.kappa + L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, d, a, b, c,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, d, a, b, c,  -(L+L+1) / 2 * xcc) )   end
        # mu = 6
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L) * (c.subshell.kappa - a.subshell.kappa - L - 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, a, d, c, b,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, a, d, c, b,  -(L+L+1) / 2 * xcc) )   end
        # mu = 7
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa - L - 1) * (c.subshell.kappa - a.subshell.kappa - L)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, b, c, d, a,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, b, c, d, a,  -(L+L+1) / 2 * xcc) )   end
        # mu = 8
        xcc = xc * wb * (d.subshell.kappa - b.subshell.kappa + L) * (c.subshell.kappa - a.subshell.kappa + L + 1)
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L+1, c, b, a, d,   (L+L+1) / 2 * xcc) )   end
        if  abs(xcc) > 1.0e-10   push!( xcList, XLCoefficient('S', L-1, c, b, a, d,  -(L+L+1) / 2 * xcc) )   end
    end

    return( xcList )
end


"""
`InteractionStrength.XL_Breit_densities(xcList::Array{XLCoefficient,1}, factor::Float64, grid::Radial.Grid)`  
    ... computes the the effective Breit interaction strengths X^L,0_Breit (abcd) for given rank L and a list of 
        orbital functions a, b, c, d and angular coefficients at the given grid. A value::Float64 is returned. 
        At present, only the zero-frequency Breit interaction is taken into account.
"""
function XL_Breit_densities(xcList::Array{XLCoefficient,1}, factor::Float64, grid::Radial.Grid)
    function V(nu::Int64, r::Float64, s::Float64, omega::Float64)
        wx = 1.0;            
        if      omega <= 0.  
        elseif  nu < 0       @show "V", nu
        elseif   r < s   
            try     wx = -(2*nu+1) * GSL.sf_bessel_jl(nu, omega*r) * GSL.sf_bessel_yl(nu, omega*s)
            catch
                    wx = 1.0;  @show  "Va", nu, omega*r, omega*s
            end
        else             
            try     wx = -(2*nu+1) * GSL.sf_bessel_jl(nu, omega*s) * GSL.sf_bessel_yl(nu, omega*r)
            catch
                    wx = 1.0;  @show  "Vb", nu, omega*r, omega*s
            end
        end
        return(wx)
    end
    #
    function W(nu::Int64, r::Float64, s::Float64, omega::Float64)
        wx = 1.0
        if      omega <= 0.  
        elseif      nu < 1   ## @show "W", nu
        elseif   r < s   
            try     wx = -(2*nu+1) * GSL.sf_bessel_jl(nu-1, omega*r) * GSL.sf_bessel_yl(nu+1, omega*s) +
                            ( (2*nu+1)/omega )^2 * r^(nu-1) / s^(nu+2)
            catch
                    wx = 1.0;  @show  "Wa", nu, omega*r, omega*s
            end
        elseif   r > s
        else             
            try     wx = -(2*nu+1) * GSL.sf_bessel_jl(nu-1, omega*r) * GSL.sf_bessel_yl(nu+1, omega*s) +
                            ( (2*nu+1)/omega )^2 * r^(nu-1) / s^(nu+2)
            catch
                    wx = 1.0;  @show  "Wb", nu, omega*r, omega*s
            end
        end
        return(wx)
    end
    
    if  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  xc  in  xcList  ## [end:end]
            # Use the minimal extent of any involved orbitals; this need to be improved
            mtp_ac = min(size(xc.a.P, 1), size(xc.c.P, 1));     mtp_bd = min(size(xc.b.P, 1), size(xc.d.P, 1))
            omg_ac = factor * abs(xc.a.energy - xc.c.energy);   omg_bd = factor * abs(xc.b.energy - xc.d.energy)
            for  r = 2:mtp_ac
                for  s = 2:mtp_bd
                    if      factor  == 0.    wy = 1.0
                    elseif  factor  == 1.    wy = 1.05
                    elseif  xc.kind == 'S'   wy = 1.0
                                                ## wy = (W(xc.nu, grid.r[r], grid.r[s], omg_ac) + W(xc.nu, grid.r[r], grid.r[s], omg_bd)) / 2.0
                    elseif  xc.kind == 'T'   wy = (V(xc.nu, grid.r[r], grid.r[s], omg_ac) + V(xc.nu, grid.r[r], grid.r[s], omg_bd)) / 2.0
                    else    error("stop a")
                    end
                    #
                    wc = xc.coeff * grid.wr[r] * grid.wr[s] 
                    #
                    if      s > r   continue
                    elseif  s == r
                        wa = wa + wy * wc * (xc.a.P[r] * xc.c.Q[r]) * (grid.r[s]^xc.nu) / (grid.r[r]^(xc.nu+1)) * (xc.b.P[s] * xc.d.Q[s]) / 2.0
                    else 
                        wa = wa + wy * wc * (xc.a.P[r] * xc.c.Q[r]) * (grid.r[s]^xc.nu) / (grid.r[r]^(xc.nu+1)) * (xc.b.P[s] * xc.d.Q[s])
                    end
                end
            end
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`InteractionStrength.XL_Coulomb_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the the effective Coulomb interaction strengths X^L_Coulomb (abcd) for given rank L and orbital functions 
        a, b, c and d at the given grid but without optimization. A value::Float64 is returned.
"""
function XL_Coulomb_WO(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid)
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
    
    XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim_WO(L, a, b, c, d, grid)
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
        XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim(L, a, b, c, d, grid)

    end

    return( XL_Coulomb )
end


"""
`InteractionStrength.XL_CoulombClaude_reset_storage(keep::Bool; printout::Bool=false)`
    ... resets the global storage of XL_CoulombClaude interaction strengths (a SEPARATE cache from
        XL_Coulomb's own GBL_Storage_XL_Coulomb, since the kink-aware and standard quadratures give
        different numeric results for the same subshell labels and must not share a cache namespace);
        nothing is returned.
"""
function XL_CoulombClaude_reset_storage(keep::Bool; printout::Bool=false)
    if  keep
        if printout     println(">> Reset GBL_Storage_XL_CoulombClaude storage.")     end
        global GBL_Storage_XL_CoulombClaude = Dict{String, Float64}()
    else
    end
    return( nothing )
end


"""
`InteractionStrength.XL_CoulombClaude(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)`
    ... computes the same effective Coulomb interaction strength as XL_Coulomb(L, a, b, c, d, grid), including the
        same triangular-delta veto and angular reduced-matrix-element prefactor xc, but using the kink-aware
        RadialIntegrals.SlaterRk_2dimClaude for the underlying radial integral instead of RadialIntegrals.
        SlaterRk_2dim. For keep=true, looks up (and stores into) the global GBL_Storage_XL_CoulombClaude
        Dict, mirroring XL_Coulomb's own keep/GBL_Storage_XL_Coulomb pattern exactly -- used by
        Hamiltonian.setupMatrixClaude (the CI-matrix Coulomb term for ALField/EOLField), where orbitals
        are FIXED for the whole performCIClaude call, so caching is unconditionally safe there. NOT enabled
        (keep=false, the default) at SelfConsistent.computeTwoElectronVClaude2's own call site: that call sits
        inside the outer SCF iteration, where orbitals change every iteration, so a cache surviving across
        iterations would silently return stale integrals from an earlier orbital shape -- extending caching
        safely into that loop needs its own explicit per-iteration reset wiring, deferred as a separate item.
        Isolated from XL_Coulomb; shared by Hamiltonian.setupMatrixClaude and
        SelfConsistent.computeTwoElectronVClaude2 (their Fock matrix, uncached). A value::Float64 is
        returned.
"""
function XL_CoulombClaude(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid; keep::Bool=false)
    global GBL_Storage_XL_CoulombClaude
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
        if haskey(GBL_Storage_XL_CoulombClaude, sa)
            return( GBL_Storage_XL_CoulombClaude[sa] )
        end
    end

    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if   rem(L,2) == 1    xc = - xc    end

    XL_CoulombClaudeValue = xc * RadialIntegrals.SlaterRk_2dimClaude(L, a, b, c, d, grid)

    if  keep
        sa = "XL" * string(L) * " " * string(a.subshell) * string(b.subshell) * string(c.subshell) * string(d.subshell)
        global GBL_Storage_XL_CoulombClaude[sa] = XL_CoulombClaudeValue
    end

    return( XL_CoulombClaudeValue )
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
            wm[i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Pa, b.P, Pc, d.P, grid) + 
                      RadialIntegrals.SlaterRkComponent_2dim(L, Pa, b.Q, Pc, d.Q, grid)
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
            wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Qa, b.P, Qc, d.P, grid) + 
                              RadialIntegrals.SlaterRkComponent_2dim(L, Qa, b.Q, Qc, d.Q, grid)
        end
    end
    
    return( xc * wm )
end


"""
`InteractionStrength.XL_CoulombClaude(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)`
    ... computes the same (direct) Coulomb interaction strengths X^L_Coulomb (.b.d) as
        XL_Coulomb(L,a::Subshell,b::Orbital,c::Subshell,d::Orbital,primitives), for given rank L and orbital
        functions as well as the given primitives, but using the kink-aware screened-potential construction
        (RadialIntegrals.buildScreenedPotentialClaude) instead of RadialIntegrals.SlaterRkComponent_2dim's naive
        tensor-product double sum. The screened potential V_L(r), which depends only on the fixed orbital pair
        (b,d), is built ONCE (adaptive quadrature) and then reused cheaply -- via the existing grid quadrature
        weights, since V_L(r) is smooth once built -- for every B-spline pair (i,k) of the L- and S-block.
        Isolated from XL_Coulomb; shared by the ALField/EOLField code lines, cf.
        SelfConsistent.computeTwoElectronVClaude2. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_CoulombClaude(L::Int64, a::Subshell, b::Orbital, c::Subshell, d::Orbital, primitives::Bsplines.Primitives)
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
    # mtpOut is forced to the full grid extent: unlike SlaterRk_2dimClaude's orbital-orbital use (where the
    # OTHER factor in the contraction is also naturally truncated to some orbital's own extent), here Vk gets
    # contracted against B-spline ROW/COLUMN indices that span the FULL basis and can extend well past (b,d)'s
    # own reach -- leaving mtpOut at its default silently drops that tail and was traced to a real bug (Ne's 1s
    # orbital coming out measurably too deeply bound from a missing part of its screening by the 2p shell).
    Vk = RadialIntegrals.buildScreenedPotentialClaude(L, b, d, grid; mtpOut=grid.NoPoints)

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
            wm[i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Pa, b.P, c.P, Pd, grid)
        end
        for  k = 1:nsS 
            ## Ba = primitives.bsplinesL[i].bs;    Bd = primitives.bsplinesS[k].bs
            Ba = primitives.bsplinesL[i];    Bd = primitives.bsplinesS[k]
            Pa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Pa[j] = Pa[j] + Ba.bs[j+add]   end
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Qd[j] = Qd[j] + Bd.bs[j+add]   end            
            wm[i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Pa, b.P, c.Q, Qd, grid)
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
            wm[nsL+i,k] = RadialIntegrals.SlaterRkComponent_2dim(L, Qa, b.Q, c.P, Pd, grid)
        end
        for  k = 1:nsS 
            ## Ba = primitives.bsplinesS[i].bs;    Bd = primitives.bsplinesS[k].bs
            Ba = primitives.bsplinesS[i];    Bd = primitives.bsplinesS[k]
            Qa = zeros(Ba.upper);   add = 1 - Ba.lower;   
            for  j = Ba.lower:Ba.upper  Qa[j] = Qa[j] + Ba.bs[j+add]   end
            Qd = zeros(Bd.upper);   add = 1 - Bd.lower;   
            for  j = Bd.lower:Bd.upper  Qd[j] = Qd[j] + Bd.bs[j+add]   end            
            wm[nsL+i,nsL+k] = RadialIntegrals.SlaterRkComponent_2dim(L, Qa, b.Q, c.Q, Qd, grid)
        end
    end

    return( xc * wm )
end


"""
`InteractionStrength.XL_CoulombTensorClaude(L::Int64, a::Subshell, b::Orbital, c::Orbital, cVector::Vector{Float64},
                                            d::Subshell, cacheLL::RadialIntegrals.SlaterMomentCacheClaude,
                                            cacheLS::RadialIntegrals.SlaterMomentCacheClaude,
                                            cacheSS::RadialIntegrals.SlaterMomentCacheClaude,
                                            primitives::Bsplines.Primitives)`
    ... computes the same (exchange) Coulomb interaction strengths X^L_Coulomb (.bc.) as
        XL_Coulomb(L,a::Subshell,b::Orbital,c::Orbital,d::Subshell,primitives) / XL_CoulombClaude of the same
        signature, but using PRECOMPUTED RadialIntegrals.SlaterMomentCacheClaude tensors -- built ONCE, outside
        the SCF iteration, via RadialIntegrals.buildSlaterMomentCacheClaude for rank L -- instead of any per-call
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
        RadialIntegrals.buildSlaterMomentCacheClaude. This turns what used to be an expensive per-call
        adaptive-quadrature computation into one with NO further quadrature at all, for every SCF iteration and
        every exchange coefficient that shares this rank L, once the caches themselves have been built (see
        SelfConsistent.solveAverageLevelFieldClaude for how they get built and cached once per SCF run).
        Isolated from XL_Coulomb; only used by the ALFieldClaude code line, cf.
        SelfConsistent.computeDirectExchangeVClaude. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function XL_CoulombTensorClaude(L::Int64, a::Subshell, b::Orbital, c::Orbital, cVector::Vector{Float64}, d::Subshell,
                                cacheLL::RadialIntegrals.SlaterMomentCacheClaude,
                                cacheLS::RadialIntegrals.SlaterMomentCacheClaude,
                                cacheSS::RadialIntegrals.SlaterMomentCacheClaude,
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
    
    XL_Coulomb = xc * RadialIntegrals.SlaterRk_2dim_Damped(tau::Float64, L, a, b, c, d, grid)
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

    XL_Coulomb_DH = xc * RadialIntegrals.SlaterRk_DebyeHueckel_2dim(L, a, b, c, d, grid, lambda)
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

