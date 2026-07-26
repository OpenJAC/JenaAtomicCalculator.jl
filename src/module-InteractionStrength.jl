
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
    wa = AngularMomentum.CL_reduced_me_rb(a.subshell, k, b.subshell) * RadialIntegrals.rkDiagonal(k, a, b, grid)
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

#==
"""
`InteractionStrength.hfs_t1(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(1) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment 
        of the nucleus for orbital functions a, b. A value::Float64 is returned.  
"""
function hfs_t1(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - (a.subshell.kappa + b.subshell.kappa) * AngularMomentum.CL_reduced_me_rb(minusa, 1, b.subshell)
    wc =   RadialIntegrals.rkNonDiagonal(-2, a, b, grid)
    wa =   wb * wc
    #
    return( wa )
end


"""
`InteractionStrength.hfs_t2(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(2) ||b> reduced matrix element for the HFS coupling to the electric-quadrupole moment of 
        the nucleus for orbital functions a, b. A value::Float64 is returned.  
    """
function hfs_t2(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    wb = - AngularMomentum.CL_reduced_me_rb(a.subshell, 2, b.subshell)
    ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wa =   wb * wc
    #
    ## println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
    return( wa )
end
==#


"""
`InteractionStrength.hfs_tM1(a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... computes the <a|| t^(1) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment
        of the nucleus for orbital functions a, b. A value::Float64 is returned.

        Note (26-Jul-2026): was missing the alpha (fine-structure constant) prefactor required by Andersson &
        Jonsson (2008), CPC 178, Eq. (49): <n_a kappa_a || t^(1) || n_b kappa_b> = -alpha(kappa_a+kappa_b)
        <-kappa_a||C^(1)||kappa_b>[r^-2] -- confirmed by direct re-reading of the paper (p. 161). Combined
        with a separate, second fix in module-Hfs.jl (a missing sqrt(2j_a+1) normalization undo), this
        resolves the long-standing H(1s) HFS A-constant discrepancy: verified to 0.19% against the 21 cm
        line and 1.65% against Na(3s)'s well-known hyperfine constant -- see examples/example-Cb.jl and
        memory project_zeeman_hfs_bugs.md.
"""
function hfs_tM1(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - Defaults.getDefaults("alpha") * (a.subshell.kappa + b.subshell.kappa) *
                AngularMomentum.CL_reduced_me_sms(minusa, 1, b.subshell)
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
    wb =   - (a.subshell.kappa + b.subshell.kappa) * AngularMomentum.CL_reduced_me_sms(minusa, 2, b.subshell)
    wc =   RadialIntegrals.rkNonDiagonal(-3, a, b, grid)/2
    wa =   wb * wc
    #
    return( wa )
end

"""
`InteractionStrength.hfs_tM3(a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... computes the <a|| t^(M3) ||b> reduced matrix element for the HFS coupling to the magnetic-dipole moment 
        of the nucleus for orbital functions a, b. A value::Float64 is returned.  
"""
function hfs_tM3(a::Orbital, b::Orbital, grid::Radial.Grid)
    # Use Andersson, Jönson (2008), CPC, Eq. (49) ... test for the proper definition of the C^L tensors.
    minusa = Subshell(1, -a.subshell.kappa)
    wb =   - (a.subshell.kappa + b.subshell.kappa) * AngularMomentum.CL_reduced_me_sms(minusa, 3, b.subshell)
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
    wb = - AngularMomentum.CL_reduced_me_sms(a.subshell, 2, b.subshell)
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
    wb = - AngularMomentum.CL_reduced_me_sms(a.subshell, 1, b.subshell)
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
    wb = - AngularMomentum.CL_reduced_me_sms(a.subshell, 3, b.subshell)
    ## wc =   RadialIntegrals.rkDiagonal(-3, a, b, grid)
    wc =   RadialIntegrals.rkDiagonal(-4, a, b, grid)
    wa =   wb * wc
    #
    ## println("**  <$(a.subshell) || t2 || $(b.subshell)>  = $wa   = $wb * $wc" )
    return( wa )
end

"""
`InteractionStrength.MbaAbsorptionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp, absorption) || a> 
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
        within the given gauge. A value::Float64 is returned.  
"""
function MbaAbsorptionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    wa = MbaEmissionCheng(mp, gauge, omega, b, a, grid)
    wa = conj(wa)
    return( wa )
end



"""
`InteractionStrength.MbaEmissionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp,emission) || a> 
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
        within the given gauge. A value::Float64 is returned.  
        This procedure has been first worked out with Andrey; in this case, however, the phases are not under good control,
        and this gives rise to wrong amplitudes and rates. The procedure is currently not in use.
"""
function MbaEmissionCheng(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    #
    if       gauge == Basics.Magnetic
        ChengI = AngularMomentum.ChengI(-kapa, kapb, AngularJ64(mp.L))
        wa     = -1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * (kapb + kapa) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
    #
    elseif   gauge == Basics.Babushkin
        ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L))
        wr     = (mp.L+1) * RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapb-kapa-mp.L-1) * RadialIntegrals.GrantIL0(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapb-kapa+mp.L+1) * RadialIntegrals.GrantIL0(mp.L+1, q, b, a, grid::Radial.Grid)
        wa     = 1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * wr
    #
    elseif   gauge == Basics.Coulomb
        ChengI = AngularMomentum.ChengI(kapa, kapb, AngularJ64(mp.L))
        wr     = (mp.L+1) * (kapb-kapa+mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, a, b, grid::Radial.Grid)
        wr     = wr +  ((mp.L+1) * (kapb-kapa+mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                            RadialIntegrals.GrantIL0(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr + (mp.L+1) * (kapb-kapa-mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, b, a, grid::Radial.Grid)
        wr     = wr +  ((mp.L+1) * (kapb-kapa-mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                            RadialIntegrals.GrantIL0(mp.L+1, q, b, a, grid::Radial.Grid) 
        wa     = 1.0im / sqrt(mp.L*(mp.L+1)) * ChengI * wr
    else     error("stop a")
    end
    wa = conj(wa)

    return( wa )
end



"""
`InteractionStrength.MbaEmissionJohnsonx(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp,emission) || a> 
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
        within the given gauge. A value::Float64 is returned.  
        This procedure has been adapted from Jiri's work but modified for the Coulomb gauge which was apparently wrong following
        some former implementation with RATIP.
"""
function MbaEmissionJohnsonx(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    #
    if       gauge == Basics.Magnetic
        JohnsonI = AngularMomentum.JohnsonI(-kapb, kapa, AngularJ64(mp.L))
        wa     = JohnsonI * (kapa + kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
    #
    elseif   gauge == Basics.Babushkin
        JohnsonI = AngularMomentum.JohnsonI(kapb, kapa, AngularJ64(mp.L))
        wr     = - RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapb-kapa)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr +  RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid)
        wa     = JohnsonI * wr
    #
    elseif   gauge == Basics.Coulomb
        # This reduced matrix element still follows RATIP apart from JohnsonI
        ChengI = AngularMomentum.JohnsonI(kapb, kapa, AngularJ64(mp.L))
        wr     = (mp.L+1) * (kapb-kapa+mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, a, b, grid::Radial.Grid)
        wr     = wr +  ((mp.L+1) * (kapb-kapa+mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                            RadialIntegrals.GrantIL0(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr + (mp.L+1) * (kapb-kapa-mp.L) / (2mp.L+1) * RadialIntegrals.GrantIL0(mp.L-1, q, b, a, grid::Radial.Grid)
        wr     = wr +  ((mp.L+1) * (kapb-kapa-mp.L) - (2mp.L+1) * (kapb-kapa)) / (2mp.L+1) * 
                                                            RadialIntegrals.GrantIL0(mp.L+1, q, b, a, grid::Radial.Grid) 
        wa     = 1.0 / (mp.L+1) * ChengI * wr
    else     error("stop a")
    end

    return( wa )
end



"""
`InteractionStrength.MabEmissionJohnsony(mp::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) interaction strength <a || O^(Mp,emission) || b> 
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
        within the given gauge. A value::Float64 is returned. This procedure has been re-worked due to the book by Johnson (2007). 
"""
function MabEmissionJohnsony(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    #
    if       gauge == Basics.Magnetic
        JohnsonI = AngularMomentum.JohnsonI(-kapa, kapb, AngularJ64(mp.L))
        wa     = JohnsonI * (kapa + kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
        #
    elseif   gauge == Basics.Babushkin
        JohnsonI = AngularMomentum.JohnsonI(kapa, kapb, AngularJ64(mp.L))
        wr     = RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapa-kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr +  RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid)
        wa     = JohnsonI * wr
        #
    elseif   gauge == Basics.Coulomb
        # This reduced matrix element has been re-worked due to Johnson (2007)
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
`InteractionStrength.MabEmissionJohnsony_Wu(mp::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute the (single-electron reduced matrix element) interaction strength <a || O^(Mp,emission) || b> 
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega, and 
        within the given gauge. The caluclation is performed by using the function AngularMomentum.CL_reduced_me_sms. 
        A value::Float64 is returned. This procedure has been re-worked due to the book by Johnson (2007). 
"""
function MabEmissionJohnsony_Wu(mp::EmMultipole, gauge::EmGauge, omega::Float64, a::Orbital, b::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    #
    if       gauge == Basics.Magnetic
        minusa = Subshell(1, -kapa)
        JohnsonI =AngularMomentum.CL_reduced_me_sms(minusa, mp.L, b.subshell)*sqrt((2*mp.L+1)*(mp.L+1) /(4*pi*mp.L ));
        wa     = JohnsonI * (kapa + kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
        #
    elseif   gauge == Basics.Babushkin
        JohnsonI = AngularMomentum.CL_reduced_me_sms(a.subshell,mp.L, b.subshell)*sqrt((2*mp.L+1)*(mp.L+1) /(4*pi*mp.L ));
        wr     = RadialIntegrals.GrantJL(mp.L, q, a, b, grid::Radial.Grid)
        wr     = wr +  (kapa-kapb)/(mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid)
        wr     = wr +  RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid)
        wa     = -JohnsonI * wr
        #
    elseif   gauge == Basics.Coulomb
        # This reduced matrix element has been re-worked due to Johnson (2007)
        JohnsonI = AngularMomentum.CL_reduced_me_sms(a.subshell, mp.L, b.subshell)*sqrt((2*mp.L+1)*(mp.L+1) /(4*pi*mp.L ));
        wr       = (1 - mp.L/(2mp.L+1)) * RadialIntegrals.GrantILplus(mp.L-1, q, a, b, grid::Radial.Grid)  -
                    mp.L/(2mp.L+1) * RadialIntegrals.GrantILplus(mp.L+1, q, a, b, grid::Radial.Grid) 
        wr       = -(kapa-kapb) / (mp.L+1) * wr 
        wr       = wr  +  mp.L/(2mp.L+1) * RadialIntegrals.GrantILminus(mp.L-1, q, a, b, grid::Radial.Grid) 
        wr       = wr  +  mp.L/(2mp.L+1) * RadialIntegrals.GrantILminus(mp.L+1, q, a, b, grid::Radial.Grid) 
        wa       = -JohnsonI * wr
    else     error("stop a")
    end

    return( wa )
end


"""
`InteractionStrength.MbaEmissionAndrey(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)`  
    ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(Mp, emission) || a>  
        for the interaction with the Mp multipole component of the radiation field and the transition frequency omega,
        and within the given gauge. A value::Float64 is returned. At present, only the magnetic matrix elements are 
        implemented. 
        This procedure has been worked out with Andrey but is currently not in use.
"""
function MbaEmissionAndreyOld(mp::EmMultipole, gauge::EmGauge, omega::Float64, b::Orbital, a::Orbital, grid::Radial.Grid)
    kapa = a.subshell.kappa;   kapb = b.subshell.kappa;    q = omega / Defaults.getDefaults("speed of light: c") 
    ja   = Basics.subshell_j(a.subshell);   jb   = Basics.subshell_j(b.subshell);   
    #
    if       gauge == Basics.Magnetic
        wa = -1.0im * sqrt( (2mp.L+1) * (Basics.subshell_2j(a.subshell)+1) ) / sqrt(4pi) / sqrt(mp.L * (mp.L+1)) * (-1)^mp.L * (kapb + kapa)
        wa = wa * AngularMomentum.ClebschGordan(ja, AngularM64(1//2), AngularJ64(mp.L), AngularM64(0), jb, AngularM64(1//2)) 
        wa = wa * RadialIntegrals.GrantILplus(mp.L, q, a, b, grid::Radial.Grid)
        wa = conj(wa)
    #
    elseif   gauge == Basics.Babushkin
        wa = 0.
    #
    elseif   gauge == Basics.Coulomb
        wa = 0.
    else     error("stop a")
    end
        
    return( wa )
end

"""
`InteractionStrength.MbaEmissionMigdalek(cp::CorePolarization, a::Orbital, b::Orbital, grid::Radial.Grid)`  
    ... to compute the (single-electron reduced matrix element) interaction strength <b || O^(E1, emission with core-polarization) || a>  
        in length gauge. A value::Float64 is returned. 
"""
function MbaEmissionMigdalek(cp::CorePolarization, a::Orbital, b::Orbital, grid::Radial.Grid)
    mp   = E1;    kapa = a.subshell.kappa;   kapb = b.subshell.kappa
    ja   = Basics.subshell_j(a.subshell);    jb   = Basics.subshell_j(b.subshell);   
    #
    JohnsonI = AngularMomentum.JohnsonI(kapa, kapb, AngularJ64(mp.L))
    wr       = RadialIntegrals.GrantJL_cp(mp.L, 0., a, b, grid::Radial.Grid, cp::CorePolarization)
    @show "MbaEmissionMigdalek", wr
    wa       = JohnsonI * wr
        
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
    wa = AngularMomentum.CL_reduced_me_sms(a.subshell, 1, c.subshell) * 
            AngularMomentum.CL_reduced_me_sms(b.subshell, 1, d.subshell) *
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
    wa = - AngularMomentum.CL_reduced_me_sms(b.subshell, 1, d.subshell) * RadialIntegrals.Vinti(b, d, grid) *
            RadialIntegrals.isotope_smsB(a, c, nm.Z, grid) / 2
    return( wa )
end


"""
`InteractionStrength.X_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)` 
    ... computes the the effective interaction strengths X^1_sms,C (abcd) for fixed rank 1 and orbital functions 
        a, b, c and d at the given grid. A value::Float64 is returned.
"""
function X_smsC(a::Orbital, b::Orbital, c::Orbital, d::Orbital, nm::Nuclear.Model, grid::Radial.Grid)
    wa = - AngularMomentum.CL_reduced_me_sms(b.subshell, 1, d.subshell) * 
            AngularMomentum.CL_reduced_me_sms(a.subshell, 1, c.subshell) * 
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
    ang = AngularMomentum.CL_reduced_me_rb(Subshell(1, -ka), 1, b.subshell)

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
    ang = AngularMomentum.CL_reduced_me_rb(Subshell(1, -ka), 1, b.subshell)

    return ( -rad * ang/(2 * Defaults.getDefaults("alpha")) * (ka + kb) )
end


end # module

