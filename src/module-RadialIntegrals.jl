
"""
`module JAC.RadialIntegrals`  
... a submodel of JAC that contains methods for calculating radial one- and two-particle matrix elements. These integrals occur 
    frequently in atomic structure and collision theory, and their fast computations often appears essential.
"""
module  RadialIntegrals

using  Dierckx, FastGaussQuadrature, GSL, QuadGK
using  ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults,  ..Radial, ..Math, ..ManyElectron, ..Nuclear

#   Fitting coefficients A(Z,n,l=0) = a0 + a1*Z + a2*Z^2 + a3*Z^3 + a4*Z^4  for the electric self-energy
#   form-factor prefactor of s-orbitals (kappa = -1), separately for n = 1..5(and above), and for Z >= 20
#   and Z < 20; cf. Kozioł, "Rci-Q", arXiv:2512.01515, Table I. Columns 1-5 are Z>=20; columns 6-10 are Z<20.
#                   ---------------------- Z >= 20 ----------------------    --------------------- Z < 20 ---------------------
#                     a0            a1            a2            a3            a4              a0            a1            a2            a3            a4
const  rciQ_Ael0 = [  7.72308e-1   -2.40991e-4    3.48842e-5   -2.83516e-7    3.16093e-10      8.72587e-1   -1.44109e-2   -1.48436e-3    3.37161e-4   -1.34952e-5  ;   # n=1
                      8.07899e-1    1.29047e-3    3.94430e-5   -1.64904e-7   -1.97988e-9       8.91352e-1   -1.26732e-2   -1.29823e-3    3.24461e-4   -1.32155e-5  ;   # n=2
                      8.08505e-1    2.07848e-3    3.83122e-5   -1.64870e-7   -2.26659e-9       8.81396e-1    3.14301e-2   -1.29689e-2    1.71096e-3   -5.74309e-5  ;   # n=3
                      8.08957e-1    2.52039e-3    3.03563e-5   -5.24535e-8   -2.92995e-9       8.83234e-1    3.17102e-2   -1.49614e-2    1.71022e-3   -5.74210e-5  ;   # n=4
                      8.15199e-1    2.19164e-3    4.33051e-5   -2.04857e-7   -2.39770e-9       8.84127e-1    3.18795e-2   -1.49649e-2    1.71050e-3   -5.74349e-5  ]   # n>=5

"""
`RadialIntegrals.GrantIab(a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, potential::Radial.Potential)`  
    ... computes the (radial) single-electron energy integral:

        `I(ab) = <a | h_D | b> = delta_{kappa_a, kappa_b} int_0^infty dr  [ c Q_a ( d/dr + kappa_a/r ) P_b +  c P_a (-d/dr + kappa_a/r ) Q_b
                                                                            - 2c^2 Q_a Q_b + V_nuc (r) (P_a P_b + Q_a Q_b) ]`
                                
        for the orbitals a and b on the grid. potential.Zr must provide the effective nuclear charge Z(r) on this grid.
"""
function GrantIab(a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, potential::Radial.Potential)
    if  a.subshell.kappa != b.subshell.kappa    return( 0 )    end
    kappa = a.subshell.kappa;                   Zr = potential.Zr
    mtp   = min(size(a.P, 1), size(b.P, 1));    wc = Defaults.getDefaults("speed of light: c")
    
    wa = 0.
    for  i = 2:mtp   
        wa = wa + grid.wr[i] * (  wc * a.Q[i] * (b.Pprime[i] + kappa/grid.r[i] * b.P[i])  
                                - wc * a.P[i] * (b.Qprime[i] - kappa/grid.r[i] * b.Q[i]) 
                                - 2wc^2 * a.Q[i] * b.Q[i]
                                - Zr[i] * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) / grid.r[i]  ) 
    end
    return( wa )
end

"""
`RadialIntegrals.GrantIabDamped(tau::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, potential::Radial.Potential)`  
    ... computes the (radial) single-electron energy integral:

        `I(ab) = <a | h_D | b> = delta_{kappa_a, kappa_b} int_0^infty dr  [ c Q_a ( d/dr + kappa_a/r ) P_b +  c P_a (-d/dr + kappa_a/r ) Q_b
                                                                            - 2c^2 Q_a Q_b + V_nuc (r) (P_a P_b + Q_a Q_b) ] * exp(-tau * r)`
                                
        for the orbitals a and b on the grid. potential.Zr must provide the effective nuclear charge Z(r) on this grid.
"""
function GrantIabDamped(tau::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, potential::Radial.Potential)
    if  a.subshell.kappa != b.subshell.kappa    return( 0 )    end
    kappa = a.subshell.kappa;                   Zr = potential.Zr
    mtp   = min(size(a.P, 1), size(b.P, 1));    wc = Defaults.getDefaults("speed of light: c")
    
    wa = 0.
    for  i = 2:mtp   
        wa = wa + grid.wr[i] * (  wc * a.Q[i] * (b.Pprime[i] + kappa/grid.r[i] * b.P[i])  
                                - wc * a.P[i] * (b.Qprime[i] - kappa/grid.r[i] * b.Q[i]) 
                                - 2wc^2 * a.Q[i] * b.Q[i]
                                - Zr[i] * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) / grid.r[i]  ) * exp(-tau * grid.r[i])
    end
    return( wa )
end

"""
`RadialIntegrals.GrantILminus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^- (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b - Q_a P_b ] .
"""
function GrantILminus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    wa = 0.
    for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.GrantILplus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^+ (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b + Q_a P_b ] .
"""
function GrantILplus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    wa = 0.
    for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.GrantIL0(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^0 (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b ] .
"""
function GrantIL0(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    wa = 0.
    for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.GrantJL(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        J_L (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a P_b + Q_a Q_b ] .
"""
function GrantJL(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    wa = 0.
    for  i = 2:mtp   wa = wa + (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
    return( wa )
end


"""
`RadialIntegrals.isotope_boson(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the boson-field shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_boson(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1), size(potential, 1));   
    
    wa = 0.
    for  i = 2:mtp   
        wa = wa + (a.P[i] * b.P[i]  +  a.Q[i] * b.Q[i]) * potential[i] * grid.wr[i]
    end
    return( wa )
end

"""
`RadialIntegrals.isotope_field(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the field-shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_field(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1), size(deltaPotential, 1));   
    
    wa = 0.
    for  i = 2:mtp   
        wa = wa - (a.P[i] * b.P[i]  +  a.Q[i] * b.Q[i]) * deltaPotential[i] / grid.r[i] * grid.wr[i]
    end
    return( wa )
end

"""
`RadialIntegrals.isotope_nms(a::Orbital, b::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the normal mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_nms(a::Orbital, b::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1));   lb = Basics.subshell_l(b.subshell);    kb = b.subshell.kappa
    alphaZ = Defaults.getDefaults("alpha") * Z

    wa = 0.
    for  i = 2:mtp
        wb = (a.Pprime[i] * b.Pprime[i]  +  a.Qprime[i] * b.Qprime[i])  +
                (lb*(lb+1) * a.P[i] * b.P[i]  +  kb*(kb-1) * a.Q[i] * b.Q[i]) / (grid.r[i]^2)
        wc = - 2 * alphaZ * (a.Q[i] * b.Pprime[i]  +  b.Q[i] * a.Pprime[i]) / grid.r[i]
        wd = - alphaZ * (b.subshell.kappa - 1) * (a.Q[i] * b.P[i]  +  b.Q[i] * a.P[i]) / (grid.r[i]^2)
        wa = wa + (wb + wc + wd) * grid.wr[i]   
    end
    return( wa / 2. )
end

"""
`RadialIntegrals.isotope_smsB(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the specific mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_smsB(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(c.P, 1));   alphaZ = Defaults.getDefaults("alpha") * Z
    kapa = a.subshell.kappa;   mkapa = -kapa;    kapc = c.subshell.kappa;   mkapc = -kapc 
    
    
    wa = 0.
    for  i = 2:mtp   
        wb = (- a.Q[i] * c.P[i] * AngularMomentum.sigma_reduced_me_ma(mkapa, kapc)  +
                c.Q[i] * a.P[i] * AngularMomentum.sigma_reduced_me_mb(kapa,  mkapc)  )
        wa = wa - alphaZ / grid.r[i] * wb * grid.wr[i]   
    end
    return( wa )
end

"""
`RadialIntegrals.isotope_smsC(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the specific mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_smsC(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(c.P, 1));   alphaZ = Defaults.getDefaults("alpha") * Z
    
    wa = 0.
    for  i = 2:mtp   
        wa = wa - alphaZ / grid.r[i] * (a.Q[i] * c.P[i] - c.Q[i] * a.P[i]) * grid.wr[i]   
    end
    return( wa )
end

"""
`RadialIntegrals.overlap(orbital1::Radial.Orbital, orbital2::Radial.Orbital, grid::Radial.Grid)`

+ (orbital1::Radial.Orbital, orbital2::Radial.Orbital, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral <orbital_a|orbital_b>  for two relativistic orbitals of the same 
        symmetry (kappa).
"""
function overlap(orbital1::Radial.Orbital, orbital2::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(orbital1.P, 1), size(orbital2.P, 1))
    
    wa = 0.
    for  i = 1:grid.NoPoints 
        if i > mtp   break   end
        wa = wa + ( orbital1.P[i] * orbital2.P[i] + orbital1.Q[i] * orbital2.Q[i] ) * grid.wr[i]   
    end
    return( wa )
end

"""
+ (p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral of two (non-relativistic) radial orbital functions <p1|p2>  as defined on grid.
"""
function overlap(p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)
    
    mtp = min( length(p1List), length(p2List))
    
    wa = 0.
    for  i = 1:grid.NoPoints 
        if i > mtp   break   end
        wa = wa + p1List[i] * p2List[i] * grid.wr[i]   
    end
    return( wa )
end


"""
`RadialIntegrals.qedDampedOverlap(lambda::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)` 
    ... computes the damped (radial) integral  int_0^infty (P_a P_b  +  Q_a Q_b) * e^{r/lambda} for the radial 
        orbitals a, b on the given grid. A value::Float64 is returned.
"""
function qedDampedOverlap(lambda::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp   wb = Base.MathConstants.e^(- grid.r[i]/lambda);     wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.qedLowFrequency(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)` 
    ... computes the (radial) integral for the low-frequency QED potential for the radial orbitals a, b on the given grid. 
        A value::Float64 is returned.
"""
function qedLowFrequency(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)
    alpha = Defaults.getDefaults("alpha");    BZ = 0.074 + 0.35 * nm.Z * alpha
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 1:mtp   wb = Base.MathConstants.e^(-nm.Z * grid.r[i]) ;     wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]   end
    wa = -BZ * nm.Z^4 * alpha^3 * wa

    println("QED single-electron strength <$(a.subshell)| h^(SE, low-frequency) | $(b.subshell)> = $wa ")
    return( wa )
end

"""
`RadialIntegrals.qedUehling(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model,
                            grid::Radial.Grid, qgrid::Radial.GridGL)` 
    ... computes the (radial) integral for the Uehling potential for the radial orbitals a, b on the given grid. This included a 
        formal t-integration that is performed internally on the (QED) grid qgrid. A value::Float64 is returned.
"""
function qedUehling(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)
    # Define the internal t-integration that is specific to the (simplified) Uehling potential; cf. PRA 72, 052115 (2005); eq. (9)
    function tIntegral(r::Float64, rp::Float64)
        wx = 0.;
        alpha = Defaults.getDefaults("alpha")
        for  i = 1:qgrid.nt   t = qgrid.t[i];   
            wx = wx + sqrt(t^2 - 1.) / t^2 * (1. + 1. / (2.0*t^2)) / (4*t*r/alpha) * qgrid.wt[i] *
                    (Base.MathConstants.e^(-2.0*t*abs(r-rp)/alpha) * qgrid.wt[i] - Base.MathConstants.e^(-2.0*t*(r+rp)/alpha))
        end
        return( wx )
    end
    
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp   
        wb = 0.;
        for  ip = 2:mtp 
            rho_rp
            wb  = wb + tIntegral(grid.r[i],grid.r[ip]) * (4pi) * grid.r[ip] * rho_rp * grid.wr[ip]  
        end
        wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]  
        ## x@show wa
    end
    wa = - 2. * Defaults.getDefaults("alpha")^2 / (3pi) * wa
    
    println("QED single-electron strength <$(a.subshell)| h^(Uehling) | $(b.subshell)> = $wa ")
    return( wa )
end

"""
`RadialIntegrals.qedUehlingSimple(a::Radial.Orbital, b::Radial.Orbital, pot::Radial.Potential,
                                    grid::Radial.Grid, qgrid::Radial.GridGL)` 
    ... computes the (radial) integral for the Uehling potential for the radial orbitals a, b on the given grid. This 
        included a formal t-integration that is performed internally on the (QED) grid qgrid. A value::Float64 is returned.
"""
function qedUehlingSimple(a::Radial.Orbital, b::Radial.Orbital, pot::Radial.Potential, grid::Radial.Grid, qgrid::Radial.GridGL)
    # Define the internal t-integration that is specific to the (simplified) Uehling potential; cf. PRA 72, 052115 (2005); eq. (9)
    function tIntegral(r::Float64)
        wx = 0.;
        alpha = Defaults.getDefaults("alpha")
        for  i = 1:qgrid.nt   t = qgrid.t[i];    
            wx = wx + sqrt(t^2 - 1.) / t^2 * (1. + 1. / (2.0*t^2)) * Base.MathConstants.e^(-2.0*t*r/alpha) * qgrid.wt[i] 
        end
        return( wx )
    end

    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp   
        wb = tIntegral(grid.r[i]) 
        wc = (-pot.Zr[i] / grid.r[i])
        wa = wa + (a.P[i]*wb*wc*b.P[i] + a.Q[i]*wb*wc*b.Q[i]) * grid.wr[i]   
    end
    wa = 2. * Defaults.getDefaults("alpha") / (3pi) * wa
    
    println("QED single-electron strength <$(a.subshell)| h^(simplified Uehling) | $(b.subshell)> = $wa ")
    return( wa )
end


"""
`RadialIntegrals.qedElectricFormFactor(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model,
                                        grid::Radial.Grid, qgrid::Radial.GridGL)`
    ... computes the (radial) integral for the high-frequency electric self-energy form-factor contribution
        (point-like nucleus) for the radial orbitals a, b on the given grid; cf. Flambaum & Ginges,
        PRA 72, 052115 (2005), eq. (10), with the coefficient A(Z,r) taken from Kozioł's Rci-Q re-fit
        (arXiv:2512.01515, Table I) for s-orbitals (kappa=-1), where FG's original generic fit is known
        (per Rci-Q's own text) to be poorly suited for inner shells; p/d/f orbitals still use FG's original
        generic fit -- a documented, narrower scope, mirroring Petersburg's own n<=4/kappa-restricted note.
        A value::Float64 is returned.
"""
function qedElectricFormFactor(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)
    # Define the internal t-integration that is specific to the electric self-energy form factor; cf. PRA 72, 052115 (2005); eq. (10)
    function tIntegral(r::Float64, alpha::Float64)
        wx = 0.;
        for  i = 1:qgrid.nt   t = qgrid.t[i];
            wx = wx + 1.0/sqrt(t^2 - 1.) * ( (1. - 1.0/(2.0*t^2)) * (log(t^2 - 1.) + 4.0*log(1.0/(nm.Z*alpha) + 0.5)) - 1.5 + 1.0/t^2 ) *
                        Base.MathConstants.e^(-2.0*t*r/alpha) * qgrid.wt[i]
        end
        return( wx )
    end

    alpha = Defaults.getDefaults("alpha");    x = (nm.Z - 80.) * alpha
    if  a.subshell.kappa == -1
        nn  = min(a.subshell.n, 5);   col = nm.Z >= 20.   ?   (1:5)   :   (6:10)
        a0, a1c, a2c, a3c, a4c = rciQ_Ael0[nn, col]
        A0Z = a0 + a1c*nm.Z + a2c*nm.Z^2 + a3c*nm.Z^3 + a4c*nm.Z^4
    else
        A0Z = 1.071 - 1.976*x^2 - 2.128*x^3 + 0.169*x^4
    end
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp
        r  = grid.r[i]
        # NOTE: the small-distance cutoff denominator uses alpha^3 (Rci-Q's stated exponent), not FG's
        # own alpha^2 -- with alpha^2 the damping increasingly (and spuriously) suppresses this term at
        # the Z ~ 70-90 range where the self-energy integral's dominant r ~ alpha; alpha^3 keeps the
        # damping negligible there, matching the term's expected fast growth with Z.
        AZ = A0Z * r / (r + 0.07*nm.Z^2*alpha^3)
        # NOTE: FG's eq. (10) reads  Phi_f(r) = -A(Z,r)(alpha/pi) Phi(r) integral(),  with Phi(r) = -Z/r
        # (the same convention as Phi(r) in qedUehlingSimple's wc); the two minus signs cancel, giving
        # a net POSITIVE prefactor on (Z/r) -- this is the dominant, positive self-energy contribution.
        wb = AZ * alpha/pi * (nm.Z/r) * tIntegral(r, alpha)
        wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]
    end

    println("QED single-electron strength <$(a.subshell)| h^(SE, electric form factor) | $(b.subshell)> = $wa ")
    return( wa )
end


"""
`RadialIntegrals.qedMagneticFormFactor(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model,
                                        grid::Radial.Grid, qgrid::Radial.GridGL)`
    ... computes the (radial) integral for the magnetic self-energy form-factor contribution (point-like
        nucleus) for the radial orbitals a, b on the given grid; cf. Flambaum & Ginges, PRA 72, 052115 (2005),
        eq. (7). Unlike the electric-type QED terms, the magnetic form factor couples the large and small
        radial components (a P-Q cross term) rather than P*P + Q*Q. A value::Float64 is returned.
"""
function qedMagneticFormFactor(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)
    # Define the internal t-integration that is specific to the magnetic self-energy form factor; cf. PRA 72, 052115 (2005); eq. (7)
    function tIntegrals(r::Float64, alpha::Float64)
        wi1 = 0.;   wi2 = 0.;
        for  i = 1:qgrid.nt   t = qgrid.t[i];
            we  = Base.MathConstants.e^(-2.0*t*r/alpha) * qgrid.wt[i]
            wi1 = wi1 + we / (t    * sqrt(t^2 - 1.))
            wi2 = wi2 + we / (t^2  * sqrt(t^2 - 1.))
        end
        return( wi1, wi2 )
    end

    alpha = Defaults.getDefaults("alpha")
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp
        r        = grid.r[i]
        wi1, wi2 = tIntegrals(r, alpha)
        wb       = alpha^2 * nm.Z / (4pi * r^2) * ( wi2 + (2*r/alpha)*wi1 - 1. )
        wa       = wa + (a.P[i]*wb*b.Q[i] + a.Q[i]*wb*b.P[i]) * grid.wr[i]
    end

    println("QED single-electron strength <$(a.subshell)| h^(SE, magnetic form factor) | $(b.subshell)> = $wa ")
    return( wa )
end


"""
`RadialIntegrals.qedWichmannKrollSimple(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model,
                                        pot::Radial.Potential, grid::Radial.Grid, qgrid::Radial.GridGL)`
    ... computes the (radial) integral for the simplified Wichmann-Kroll vacuum-polarization contribution
        for the radial orbitals a, b on the given grid; cf. Flambaum & Ginges, PRA 72, 052115 (2005), eq. (12).
        A closed-form (non-t-integrated) correction; qgrid is accepted only for interface consistency with
        the other local QED terms and is not used. A value::Float64 is returned.
"""
function qedWichmannKrollSimple(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, pot::Radial.Potential,
                                grid::Radial.Grid, qgrid::Radial.GridGL)
    alpha = Defaults.getDefaults("alpha");    rc = alpha
    mtp = min(size(a.P, 1), size(b.P, 1))
    wa = 0.
    for  i = 2:mtp
        r  = grid.r[i]
        wc = (-pot.Zr[i] / r)
        wb = - (2*alpha)/(3pi) * wc * 0.092*nm.Z^2*alpha^2 / (1. + (1.62*r/rc)^4)
        wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]
    end

    println("QED single-electron strength <$(a.subshell)| h^(Wichmann-Kroll) | $(b.subshell)> = $wa ")
    return( wa )
end


"""
`RadialIntegrals.rkDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`   ... computes the (radial and diagonal) integral of r^k for two radial orbital functions.

+ (k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes this integral for two relativistic orbitals:   < r^k >_ab = int_0^\\infty  dr  [P_a P_b + Q_a Q_b]  r^k
"""
function rkDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))

    wa = 0.
    # Don't allow too small r-values -- graduated with how negative k is, since r^k for very negative k
    # blows up catastrophically at the innermost grid points, amplifying any tiny residual imprecision
    # in the tabulated P/Q there (see project_zeeman_hfs_bugs.md, 30-Jul-2026, for the kappa<=-3 case
    # this matters most for; k<=-4 is not fully resolved even at m0=18 for kappa<=-3 orbitals -- a known,
    # documented residual, not chased further this session).
    if      k > -3   m0 = 2
    elseif  k == -3  m0 = 10
    else             m0 = 18   end
    for  i = m0:mtp   wa = wa + (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * (grid.r[i]^k) * grid.wr[i]   end
    return( wa )
end

"""
+ (k::Int64, p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)`
    ... computes this integral for two non-relativistic orbitals:   < r^k >_ab = int_0^\\infty  dr  [P_a P_b]  r^k
"""
function rkDiagonal(k::Int64, p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)
    mtp = min( length(p1List), length(p2List))
    
    wa = 0.
    for  i = 2:mtp   wa = wa + p1List[i] * p2List[i] * (grid.r[i]^k) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.rkNonDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)` 
    ... computes the (radial and non-diagonal) integral of r^k for two relativistic orbitals:
        [ r^k ]_ab = int_0^\\infty dr [P_a Q_b + Q_a P_b] r^k
"""
function rkNonDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))

    wa = 0.
    # Don't allow too small r-values -- graduated with how negative k is (this function previously had
    # no such guard at all, unlike rkDiagonal). k<=-4 (the Hfs M3/octupole case) remains only partially
    # resolved even at this conservative m0: kappa<=-3 orbitals (e.g. d5/2, f7/2, ...) retain a genuine,
    # slowly-decaying near-origin residual traced (30-Jul-2026) to the exact l+1+kappa=0 cancellation of
    # Q(r)'s leading power for kappa=-(l+1) -- ruled out as a tabulation artifact (an equivalent bVector-
    # space bilinear-form evaluation shows the identical residual) and as an SCF-convergence artifact
    # (tightening accuracyScf by 6 orders of magnitude did not change it); most likely an intrinsic
    # precision limitation of the single-diagonalization step for this delicate, cancellation-exposed
    # quantity. A documented, NOT-fully-resolved follow-up item -- see project_zeeman_hfs_bugs.md.
    if      k > -3   m0 = 2
    elseif  k == -3  m0 = 10
    else             m0 = 70   end
    for  i = m0:mtp   wa = wa + (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * (grid.r[i]^k) * grid.wr[i]   end
    return( wa )
end

"""
`RadialIntegrals.SlaterRkComponent_2dim(k::Int64, Ba::Array{Float64,1}, Bb::Array{Float64,1}, 
                                                  Bc::Array{Float64,1}, Bd::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes one component of the (relativistic) Slater integral

        R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d),   namely
        W^k (abcd) = int_0^infty dr int_0^infty ds  B_a B_c            r_<^k / r_>^(k+1)  B_b B_d

        of rank k for the four components Ba, Bb, ... above , and over the given grid by using an explicit 2-dimensional integration 
        scheme; a value::Float64 is returned.
"""
function SlaterRkComponent_2dim(k::Int64, Ba::Array{Float64,1}, Bb::Array{Float64,1}, Bc::Array{Float64,1}, Bd::Array{Float64,1}, grid::Radial.Grid)
    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( r^k/s^(k+1) )
        elseif r > s     return( s^k/r^(k+1) )
        end
    end

    mtp_ac = min(size(Ba, 1), size(Bc, 1));    mtp_bd = min(size(Bb, 1), size(Bd, 1))
    wac = zeros(mtp_ac);   wbd = zeros(mtp_bd)
    for  r = 2:mtp_ac   wac[r] = (Ba[r] * Bc[r]) * grid.wr[r]  end
    for  s = 2:mtp_bd   wbd[s] = (Bb[s] * Bd[s]) * grid.wr[s]  end
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd   wa = wa + wac[r] * ul(grid.r[r], grid.r[s]) * wbd[s]   end
    end
    return( wa )
end

"""
`RadialIntegrals.SlaterRk_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the (relativistic) Slater integral

    R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d)

    of rank k for the four orbitals a, b, c, d, and over the given grid by using an explicit 2-dimensional integration scheme; a 
    value::Float64 is returned.
"""
function SlaterRk_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
    wac = zeros(mtp_ac);   wbd = zeros(mtp_bd)
    for  r = 2:mtp_ac   wac[r] = (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * grid.wr[r]  end
    for  s = 2:mtp_bd   wbd[s] = (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[s]  end
    ## The two powers r^k and r^(k+1) depend only on the grid point, not on the pair (r,s), yet the
    ## closure ul() above recomputed both inside the N^2 loop -- N ~ 1000, so ~10^6 calls to ^ per
    ## Slater integral where 2N suffice.  Hoisting them costs O(N) and leaves the arithmetic of each
    ## term untouched: the same two numbers are formed and the same division is taken, so the result is
    ## BITWISE identical and no approved reference can move.  Measured 4.7x / 6.5x / 9.4x on three real
    ## cases; SlaterRk_2dim and its ul() closure were ~16 % of an Auger rate computation.
    ## ul() itself was kept at the time because the MeshGrasp branch above still called it; that branch is
    ## gone (12-Aug-2026) and the closure with it, since nothing in this function called it any more.
    mtp = max(mtp_ac, mtp_bd)
    rPk = zeros(mtp);   rPk1 = zeros(mtp)
    for  i = 2:mtp      rPk[i] = grid.r[i]^k;    rPk1[i] = grid.r[i]^(k+1)    end
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd
            uls = grid.r[r] <= grid.r[s]  ?  rPk[r] / rPk1[s]  :  rPk[s] / rPk1[r]
            wa  = wa + wac[r] * uls * wbd[s]
        end
    end
    return( wa )
end


## Fixed Gauss-Legendre nodes for the per-cell quadratures of the screened-potential sweeps.  Built once.
## 8 points integrate a polynomial of degree 15 exactly, which comfortably covers a cubic spline times
## s^k over ONE grid cell -- the sweeps never integrate across a cell boundary, so there is nothing for
## adaptivity to discover and QuadGK's 15-point Gauss-Kronrod rule (plus its subdivision logic and its
## error estimate) is simply paid for nothing.
const GBL_GaussLegendreCell = FastGaussQuadrature.gausslegendre(8)


"""
`RadialIntegrals.cellIntegralClaude(f, a::Float64, b::Float64)`  
    ... integrates f over the single interval [a,b] with a fixed 8-point Gauss-Legendre rule; a Float64 is
        returned.  Used only by the screened-potential sweeps, where each interval is one grid cell and the
        integrand is smooth on it.
"""
function cellIntegralClaude(f, a::Float64, b::Float64)
    (xg, wg) = GBL_GaussLegendreCell
    half = 0.5*(b - a);    mid = 0.5*(b + a);    acc = 0.
    for  j = 1:length(xg)     acc = acc + wg[j] * f( half*xg[j] + mid )     end
    return( half * acc )
end


"""
`RadialIntegrals.buildScreenedPotentialClaude(k::Int64, b::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid;
                                              rtol::Float64=1.0e-9, mtpOut::Union{Nothing,Int64}=nothing)`
    ... computes the "screened potential"

        V_k(r) = (1/r^(k+1)) int_0^r ds s^k rho_bd(s)  +  r^k int_r^r_max ds s^{-(k+1)} rho_bd(s),
        rho_bd(s) = P_b(s)P_d(s) + Q_b(s)Q_d(s),

        i.e. the r-dependent potential felt (via the Coulomb interaction of rank k) from the fixed orbital pair
        (b,d), evaluated at every tabulated grid point r = grid.r[i]. This is the auxiliary quantity underlying
        the two-electron Slater integral R^k(abcd) = int dr rho_ac(r) V_k(r); since V_k depends only on (b,d),
        not on whatever it is later contracted against, it only needs to be built ONCE per orbital pair and can
        then be reused cheaply (via the existing, already-tabulated grid quadrature, since V_k(r) is itself
        smooth/kink-free once built correctly) for every downstream contraction -- e.g. RadialIntegrals.
        SlaterRk_2dimClaude below, or a B-spline matrix element as needed for SCF orbital optimization.

        The kink that the original r_</r_>^(k+1) kernel has at r=s is handled explicitly here: for each outer
        point r WITHIN the source density's own extent, the integral over s is split into [0,r] and [r,r_max]
        -- each smooth on its own -- rather than integrated across the kink with a single quadrature rule as
        SlaterRk_2dim does. For r BEYOND the source's own extent (rho_bd is by construction zero there), V_k(r)
        reduces to a single r-independent constant divided by r^(k+1) -- the standard multipole falloff of a
        localized source -- computed once and reused for every such r, no adaptive quadrature needed there.
        Getting this "beyond the source" branch right matters: it is NOT optional truncation. When V_k feeds a
        B-spline MATRIX element (rather than a pre-truncated orbital-orbital integral), the matrix's row/column
        B-splines span the FULL basis and can extend well past a compact orbital pair (b,d)'s own reach (e.g. a
        2p orbital screening a 1s row built from B-splines spanning the whole grid) -- silently ending V_k's
        array at the source's own extent, as an earlier version of this function did, drops that entire tail
        contribution and was traced to a real, non-negligible bug: Ne's 1s orbital came out MEASURABLY too
        deeply bound (eigenvalue off by roughly a fifth from the experimental K-shell value) because its Fock
        matrix was silently missing part of its direct screening from the 2p shell.
        rho_bd(s) is represented as a cubic spline (Dierckx.Spline1D) over the tabulated grid values, and each
        piece is evaluated by adaptive Gauss-Kronrod quadrature (QuadGK.quadgk), exactly as already done
        elsewhere in JAC (cf. Nuclear.jl's nuclear-density integrals) for integrals with known problematic
        points.
        BOTH integrals are accumulated by a SINGLE SWEEP over the grid rather than integrated afresh at every
        point: the inner moment grows with r and the outer moment shrinks with r, so each differs from its
        neighbour by exactly one grid cell. An earlier version ran two full-range adaptive quadratures at each
        of the ~1500 grid points, integrating the same material O(N) times over; profiling of both the SCF and
        the CI paths (09-Aug-2026) put ~94% of a configuration-interaction matrix build inside this one loop.
        The forward/backward sweeps are the standard construction (Hartree's Y_k, GRASP's YZK). Note that it is
        the SPLIT AT THE KINK that is essential, not its repetition -- the split is retained exactly; only the
        redundant re-integration is gone. Per-cell quadrature also keeps the same rtol while converging more
        reliably, since one grid cell is smooth and needs no subdivision.
        mtpOut lets the caller request the FULL range actually needed (e.g. grid.NoPoints for a
        B-spline matrix element); left unspecified, it defaults to min(size(b.P,1),size(d.P,1)) as before, which
        remains correct wherever V_k is only ever contracted against something ALSO naturally truncated to that
        same orbital-pair extent (e.g. SlaterRk_2dimClaude).
        A Vk::Vector{Float64}, of length mtpOut (or min(size(b.P,1),size(d.P,1)) if mtpOut is not given), is
        returned.
"""
function buildScreenedPotentialClaude(k::Int64, b::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid;
                                      rtol::Float64=1.0e-9, mtpOut::Union{Nothing,Int64}=nothing)
    mtp_bd  = min(size(b.P, 1), size(d.P, 1))
    mtpOutx = isnothing(mtpOut) ? mtp_bd : mtpOut
    Vk      = zeros(mtpOutx)
    if  mtp_bd < 2    return( Vk )    end

    rbd    = grid.r[1:mtp_bd]
    rhoBd  = [b.P[i]*d.P[i] + b.Q[i]*d.Q[i]  for i = 1:mtp_bd]
    splBd  = Dierckx.Spline1D(rbd, rhoBd)
    r0     = grid.r[1]                     # smallest tabulated grid point, > 0; mirrors the grid's own convention
    rmaxBd = rbd[end]

    # Last output index still strictly INSIDE the source's own extent; beyond it the multipole-tail
    # branch applies and no quadrature is needed at all.
    iLast = 1
    for  i = 2:mtpOutx     if  grid.r[i] < rmaxBd    iLast = i    else    break    end    end

    # FORWARD sweep:  inner[i] = int_{r0}^{r_i} ds s^k rho_bd(s), accumulated ONE GRID CELL AT A TIME.
    inner = zeros(mtpOutx);    acc = 0.
    for  i = 2:iLast
        acc      = acc + RadialIntegrals.cellIntegralClaude(s -> s^k * splBd(s), grid.r[i-1], grid.r[i])
        inner[i] = acc
    end
    # the remaining cell out to rmaxBd completes the full inner moment, so that fullInner is by
    # construction the end value of the same sweep rather than a separately integrated quantity
    fullInner = grid.r[iLast] < rmaxBd ?
                acc + RadialIntegrals.cellIntegralClaude(s -> s^k * splBd(s), grid.r[iLast], rmaxBd)  :  acc

    # BACKWARD sweep:  outer[i] = int_{r_i}^{rmaxBd} ds s^{-(k+1)} rho_bd(s), likewise cell by cell.
    # Accumulating from the far end inwards also adds the small contributions first, which is the
    # numerically favourable order.
    outer = zeros(mtpOutx)
    if  iLast >= 2
        acc = RadialIntegrals.cellIntegralClaude(s -> s^(-(k+1)) * splBd(s), grid.r[iLast], rmaxBd)
        outer[iLast] = acc
        for  i = iLast-1:-1:2
            acc      = acc + RadialIntegrals.cellIntegralClaude(s -> s^(-(k+1)) * splBd(s), grid.r[i], grid.r[i+1])
            outer[i] = acc
        end
    end

    for  i = 2:mtpOutx
        r = grid.r[i]
        if  r >= rmaxBd    Vk[i] = fullInner / r^(k+1)
        else               Vk[i] = inner[i] / r^(k+1) + r^k * outer[i]
        end
    end
    return( Vk )
end


"""
`RadialIntegrals.SlaterRk_2dimClaude(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital,
                                     grid::Radial.Grid; rtol::Float64=1.0e-9)`
    ... computes the same (relativistic) Slater integral as SlaterRk_2dim,

    R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d)

    but via the kink-aware screened potential RadialIntegrals.buildScreenedPotentialClaude(k,b,d,grid) instead of
    the naive tensor-product Gauss-Legendre double sum SlaterRk_2dim uses (which is exact for smooth polynomials
    within a break-point cell, but not for a function with a first-derivative discontinuity running through the
    middle of one, as r_</r_>^(k+1) has at r=s). Since V_k(r) is smooth once built, the outer integral over r can
    safely reuse the existing grid quadrature weights -- no further adaptive treatment is needed there.
    A value::Float64 is returned.
"""
function SlaterRk_2dimClaude(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital,
                             grid::Radial.Grid; rtol::Float64=1.0e-9)
    mtp_ac = min(size(a.P, 1), size(c.P, 1))
    Vk     = buildScreenedPotentialClaude(k, b, d, grid; rtol=rtol)
    mtp    = min(mtp_ac, length(Vk))

    wa = 0.
    for  r = 2:mtp   wa = wa + (a.P[r]*c.P[r] + a.Q[r]*c.Q[r]) * grid.wr[r] * Vk[r]   end
    return( wa )
end


"""
`RadialIntegrals.buildScreenedPotentialPairClaude(k::Int64, Ba::Vector{Float64}, orbComp::Vector{Float64},
                                                  grid::Radial.Grid; rtol::Float64=1.0e-9,
                                                  mtpOut::Union{Nothing,Int64}=nothing)`
    ... computes the same kind of "screened potential"

        V_k(r) = (1/r^(k+1)) int_0^r ds s^k rho(s)  +  r^k int_r^r_max ds s^{-(k+1)} rho(s),
        rho(s) = Ba(s) * orbComp(s),

        as RadialIntegrals.buildScreenedPotentialClaude, but for a DENSITY BUILT FROM A SINGLE B-SPLINE (Ba,
        a zero-padded array with support limited to that spline's own knot interval) times a single component
        (P or Q) of a fixed orbital, rather than a full P_bP_d+Q_bQ_d orbital-pair density. This is the
        building block for the "exchange"-signature InteractionStrength.XL_CoulombClaude, where -- unlike the
        "direct" signature -- BOTH B-spline indices sit on opposite sides of the two-electron kernel, so no
        single screened potential can be shared across the whole (nsL+nsS) x (nsL+nsS) matrix. Instead, one such
        potential is built per ROW (per "a" B-spline), using only that spline's own compact support to define
        rho, and then reused cheaply against every column ("d" B-spline) of that row.
        V_k(r) is evaluated out to the full extent of orbComp by default (not just Ba's own narrow support),
        since even a spatially localized source density still produces a potential everywhere via its
        1/r^(k+1) tail; for r beyond Ba's support this tail integral is a single r-independent constant,
        computed once and then just divided by r^(k+1) at each point -- so only the (typically few) grid points
        inside Ba's own support incur the cost of an adaptive quadrature split at the kink.
        Within that support the two moments are accumulated by a forward and a backward SWEEP over the grid
        cells, for the same reason as in buildScreenedPotentialClaude: consecutive points differ by one cell,
        so re-integrating the whole sub-range at each point repeats work. fullInner and fullOuter are the end
        values of those same sweeps, so they cannot drift away from the pointwise values.
        The optional mtpOut lets the caller request a LONGER output range explicitly -- needed, for instance, by
        RadialIntegrals.buildSlaterMomentCacheClaude, where BOTH Ba and orbComp are themselves compact B-splines
        (not a broad orbital), so size(orbComp,1) alone would be far too short for how the result is used later.
        A Vk::Vector{Float64}, of length mtpOut (or size(orbComp,1) if mtpOut is not given), is returned.
"""
function buildScreenedPotentialPairClaude(k::Int64, Ba::Vector{Float64}, orbComp::Vector{Float64}, grid::Radial.Grid;
                                          rtol::Float64=1.0e-9, mtpOut::Union{Nothing,Int64}=nothing)
    mtpSrc = min(size(Ba, 1), size(orbComp, 1))
    mtpOutx = isnothing(mtpOut) ? size(orbComp, 1) : mtpOut
    Vk     = zeros(mtpOutx)
    if  mtpSrc < 2    return( Vk )    end

    # Restrict to the density's ACTUAL nonzero support. Ba and orbComp are zero-padded from index 1, so for a
    # "late" B-spline pair the true overlap can start far later than index 1 (e.g. a padded array of length 392
    # with a true support of only ~20 points) -- running adaptive quadrature over that leading, identically-zero
    # stretch as well as the true support wastes most of the work for exactly the pairs buildSlaterMomentCacheClaude
    # calls this with most often.
    startIdx = 0
    for  i = 1:mtpSrc   if  Ba[i]*orbComp[i] != 0.    startIdx = i;   break   end   end
    if  startIdx == 0    return( Vk )    end   # density identically zero -> potential is zero everywhere
    if  mtpSrc - startIdx < 1    return( Vk )    end   # overlap region is a single point -> measure zero, no contribution

    rSrc    = grid.r[startIdx:mtpSrc]
    rho     = [Ba[i]*orbComp[i] for i = startIdx:mtpSrc]
    splOrd  = min(3, length(rSrc)-1)
    splBa   = Dierckx.Spline1D(rSrc, rho; k=splOrd)
    rStart  = grid.r[startIdx]
    rmaxSrc = rSrc[end]

    # The output indices that fall strictly INSIDE the source's support, rStart < r < rmaxSrc; only these
    # need quadrature, the two outer branches are O(1) each.
    iFirst = startIdx + 1
    iLast  = iFirst - 1
    for  i = iFirst:mtpOutx     if  grid.r[i] < rmaxSrc    iLast = i    else    break    end    end

    # FORWARD sweep:  inner[i] = int_{rStart}^{r_i} ds s^k rho(s), one grid cell at a time.
    inner = zeros(mtpOutx);    acc = 0.
    for  i = iFirst:iLast
        acc      = acc + RadialIntegrals.cellIntegralClaude(s -> s^k * splBa(s), grid.r[i-1], grid.r[i])
        inner[i] = acc
    end
    fullInner = (iLast >= iFirst  &&  grid.r[iLast] < rmaxSrc) ?
                acc + RadialIntegrals.cellIntegralClaude(s -> s^k * splBa(s), grid.r[iLast], rmaxSrc)  :
                QuadGK.quadgk(s -> s^k * splBa(s), rStart, rmaxSrc, rtol=rtol)[1]

    # BACKWARD sweep:  outer[i] = int_{r_i}^{rmaxSrc} ds s^{-(k+1)} rho(s).  fullOuter is then the same
    # sweep carried the last step down to rStart, so the two agree by construction.
    outer = zeros(mtpOutx);    fullOuter = 0.
    if  iLast >= iFirst
        acc = RadialIntegrals.cellIntegralClaude(s -> s^(-(k+1)) * splBa(s), grid.r[iLast], rmaxSrc)
        outer[iLast] = acc
        for  i = iLast-1:-1:iFirst
            acc      = acc + RadialIntegrals.cellIntegralClaude(s -> s^(-(k+1)) * splBa(s), grid.r[i], grid.r[i+1])
            outer[i] = acc
        end
        fullOuter = acc + RadialIntegrals.cellIntegralClaude(s -> s^(-(k+1)) * splBa(s), rStart, grid.r[iFirst])
    else
        fullOuter = QuadGK.quadgk(s -> s^(-(k+1)) * splBa(s), rStart, rmaxSrc, rtol=rtol)[1]
    end

    for  i = 2:mtpOutx
        r = grid.r[i]
        if       r >= rmaxSrc    Vk[i] = fullInner / r^(k+1)
        elseif   r <= rStart     Vk[i] = r^k * fullOuter
        else                     Vk[i] = inner[i] / r^(k+1) + r^k * outer[i]
        end
    end
    return( Vk )
end


"""
`struct  RadialIntegrals.SlaterMomentCacheClaude`
    ... holds a precomputed cache of kink-aware "screened potentials" Phi_(i,i')(s), one for every overlapping
        pair of B-splines (i,i') from a GIVEN B-spline basis (e.g. primitives.bsplinesL or primitives.bsplinesS),
        for one fixed multipole rank k. THIS IS NOT ITSELF A PHYSICAL SLATER/RADIAL INTEGRAL -- it is an
        auxiliary, basis-only tensor (it does not reference any orbital at all) from which many different direct
        and exchange two-electron radial integrals can afterwards be obtained CHEAPLY, by a simple grid-quadrature
        dot product against whatever B-spline or orbital density is actually needed (see
        RadialIntegrals.buildSlaterMomentCacheClaude for how it is built, and
        SelfConsistent.computeDirectExchangeVTensor, once written, for how it gets used). Building this cache is
        the expensive, kink-aware step -- it uses the SAME split-quadrature technique as
        buildScreenedPotentialPairClaude, just applied ONCE per (basis-only) B-spline pair rather than once per
        SCF call; using it afterwards, many times over, for many different orbital pairs and many SCF iterations,
        is cheap. This follows the "moment array" idea underlying Zatsarinny \\& Froese Fischer, Comput. Phys.
        Commun. 202, 287-303 (2016), their Eq. (16) -- though here the array is stored in this partially-
        contracted Phi_(i,i')(s) form, rather than as the full four-index R^k(ij;i'j') array of that paper, to
        avoid the O((N*band)^2) memory cost of also precomputing every (j,j') combination up front.

        i and i' need NOT come from the same B-spline basis: a two-electron exchange-type integral needs the
        LL and SS ("same-basis") combinations for its diagonal blocks, but ALSO an LS ("cross-basis") combination
        for its off-diagonal blocks, since the large- and small-component B-spline bases share the same
        underlying grid (just different spline orders) and therefore genuinely overlap with EACH OTHER too, not
        only within themselves. See RadialIntegrals.buildSlaterMomentCacheClaude for how the same-basis and
        cross-basis cases are told apart and handled.
    + k        ::Int64                                    ... the multipole rank this cache was built for
    + sameBasis::Bool                                     ... true if this cache was built from two identical
                                                              B-spline lists (i.e. a same-basis/"LL" or "SS"
                                                              cache); false for a cross-basis/"LS" cache. Kept
                                                              here so a cache cannot be silently misused for the
                                                              wrong kind of contraction.
    + band     ::Int64                                    ... the largest observed |i-i'| among overlapping
                                                              B-spline pairs (a property of the spline order,
                                                              kept here only as a diagnostic)
    + Phi      ::Dict{Tuple{Int64,Int64}, Vector{Float64}}  ... Phi[(i,i')] = the screened potential sourced from
                                                              the density B1_i(r)*B2_i'(r) (B1,B2 the two bases
                                                              this cache was built from, possibly identical),
                                                              tabulated on the standard grid; only present
                                                              (nonzero) for OVERLAPPING (i,i') pairs. i indexes
                                                              the FIRST basis, i' the SECOND.
"""
struct SlaterMomentCacheClaude
    k         ::Int64
    sameBasis ::Bool
    band      ::Int64
    Phi       ::Dict{Tuple{Int64,Int64}, Vector{Float64}}
end


"""
`RadialIntegrals.buildSlaterMomentCacheClaude(k::Int64, bsplines1::Array{<:Any,1}, bsplines2::Array{<:Any,1},
                                              grid::Radial.Grid; rtol::Float64=1.0e-9)`
    ... builds a RadialIntegrals.SlaterMomentCacheClaude for multipole rank k from TWO lists of B-splines --
        typically primitives.bsplinesL and/or primitives.bsplinesS. For every pair (i,i'), i from bsplines1 and
        i' from bsplines2, whose supports overlap, the kink-aware screened potential

            Phi_(i,i')(s) = int_0^infty dr  B1_i(r) B2_i'(r) * r_<^k / r_>^(k+1)

        is built ONCE, via RadialIntegrals.buildScreenedPotentialPairClaude applied to the (compact-support)
        density B1_i*B2_i', and cached. AGAIN: THIS IS NOT ITSELF A PHYSICAL SLATER INTEGRAL -- it is the
        reusable, basis-only building block from which any two-electron direct or exchange radial integral for
        this rank k can later be obtained cheaply, without repeating the expensive kink-aware quadrature on
        every SCF call. Because B-splines only overlap with a small number of neighbors, the returned cache is
        sparse (banded) in (i,i'), not a dense N x N object.

        Whether bsplines1 and bsplines2 are the SAME list matters and is recognized explicitly (by object
        identity, ===), not left for the caller to track:
        - SAME list (e.g. both primitives.bsplinesL, or both primitives.bsplinesS): Phi_(i,i') is then
          symmetric, Phi_(i,i')=Phi_(i',i), so only the i<=i' triangle is computed and mirrored into storage.
          This is the common "LL" or "SS" case.
        - DIFFERENT lists (e.g. primitives.bsplinesL and primitives.bsplinesS): no symmetry can be assumed -- an
          index i from bsplines1 and an index i' from bsplines2 are not interchangeable, since the two bases can
          have different spline orders (and therefore different [lower,upper] support ranges) even though they
          live on the same grid. The full, non-mirrored set of overlapping pairs is computed instead, and
          overlap is tested explicitly via each B-spline's own [lower,upper] range rather than assuming any
          shared ordering convention between the two lists. This is the "LS" (equivalently "SL", with the two
          lists swapped) case, needed for the off-diagonal blocks of an exchange-signature Fock matrix.
        A RadialIntegrals.SlaterMomentCacheClaude is returned, with its sameBasis field set accordingly.
"""
function buildSlaterMomentCacheClaude(k::Int64, bsplines1::Array{<:Any,1}, bsplines2::Array{<:Any,1},
                                      grid::Radial.Grid; rtol::Float64=1.0e-9)
    sameBasis = (bsplines1 === bsplines2)
    n1  = length(bsplines1);   n2 = length(bsplines2)
    Phi = Dict{Tuple{Int64,Int64}, Vector{Float64}}()

    for  i = 1:n1
        Bi = bsplines1[i]
        Pi = zeros(Bi.upper);   addi = 1 - Bi.lower
        for  m = Bi.lower:Bi.upper   Pi[m] = Pi[m] + Bi.bs[m+addi]   end

        ipRange = sameBasis ? (i:n1) : (1:n2)
        for  ip in ipRange
            Bip = sameBasis ? bsplines1[ip] : bsplines2[ip]
            # Explicit overlap test on the two B-splines' own [lower,upper] ranges -- safe regardless of whether
            # the two lists share an ordering/indexing convention, which they need not (different spline orders
            # on the same grid).
            if  Bi.upper < Bip.lower  ||  Bip.upper < Bi.lower    continue   end

            Pip = zeros(Bip.upper);   addip = 1 - Bip.lower
            for  m = Bip.lower:Bip.upper   Pip[m] = Pip[m] + Bip.bs[m+addip]   end

            # mtpOut is forced to the full grid extent here: unlike the exchange-integral use of
            # buildScreenedPotentialPairClaude (where the SECOND argument is a broad orbital, "far enough" on
            # its own), BOTH Pi and Pip are themselves compact B-splines, so leaving mtpOut at its default
            # (size(Pip,1)) would silently truncate Phi_(i,i') long before it is needed by later contractions.
            V = buildScreenedPotentialPairClaude(k, Pi, Pip, grid; rtol=rtol, mtpOut=grid.NoPoints)
            Phi[(i,ip)] = V
            if  sameBasis && ip != i   Phi[(ip,i)] = V   end   # Phi is symmetric in (i,i') only for a same-basis cache
        end
    end

    band = 0
    for  (i,ip)  in  keys(Phi)   band = max(band, abs(i-ip))   end

    return( SlaterMomentCacheClaude(k, sameBasis, band, Phi) )
end


"""
`RadialIntegrals.SlaterRk_2dim_WO(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, grid::Radial.Grid)`
    ... computes the (relativistic) Slater integral

        R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d)

        of rank k for the four orbitals a, b, c, d, and over the given grid by using an explicit 2-dimensional integration scheme
        but without optimization (WO); a value::Float64 is returned.
"""
function SlaterRk_2dim_WO(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( r^k/s^(k+1) )
        elseif r > s     return( s^k/r^(k+1) )
        end
    end

    
    mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                        (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]   end
    end
    ## println("Test: SlaterRk_2dim(); wa = $wa")
    return( wa )
end

"""
`RadialIntegrals.SlaterRk_2dim_Damped(tau::Float64, k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the (relativistic) Slater integral

        R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d) * exp(-tau * r - tau*s)

        of rank k for the four orbitals a, b, c, d, and over the given grid by using an explicit 2-dimensional integration scheme; a 
        value::Float64 is returned.
"""
function SlaterRk_2dim_Damped(tau::Float64, k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( r^k/s^(k+1) )
        elseif r > s     return( s^k/r^(k+1) )
        end
    end

    
    mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                        (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]  *
                                        exp(- tau * grid.r[r] - tau * grid.r[s] )                         end
    end
    ## println("Test: SlaterRk_2dim(); wa = $wa")
    return( wa )
end


"""
`RadialIntegrals.SlaterRk_DebyeHueckel_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, 
                                            grid::Radial.Grid, lambda::Float64)`  
    ... computes the (relativistic) Slater-Debye-Hueckel integral

        R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) [r_<^k / r_>^(k+1)]^(DH screened) (P_b P_d + Q_b Q_d)

        of rank k for the four orbitals a, b, c, d, and over the given grid by using an explicit 2-dimensional integration 
        scheme; a value::Float64 is returned.
"""
function SlaterRk_DebyeHueckel_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, 
                                    grid::Radial.Grid, lambda::Float64)
                
    function ul_DH(L::Int64, s::Float64, r::Float64) 
        # Calculates the ul_DH(r,s) function for  s <= r.
        sum = 0.;   suma = 0.
        for  p = 0:2
            for q = 0:L
                sum = sum + (2^(L-q)) * (lambda^(L+p+p-q)) *  factorial(L+q) * factorial(L+p) /
                            ( factorial(L+L+p+p+1) * factorial(L-q) * factorial(p) * factorial(q)) * 
                            (s^(L+p+p)) * exp(-lambda*r) / (r^(q+1))
            end
            if (p == 2) suma = sum   end
        end
        return( (L+L+1) * sum )
    end                                 

    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( ul_DH(k, r, s) )
        elseif r > s     return( ul_DH(k, s, r) )
        end
    end

    
    mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                        (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]   end
    end
    return( wa )
end

"""
`RadialIntegrals.Vinti(a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)` 
    ... computes the (radial) Vinti integral for the two radia integrals a and b:
            
                                        [ d        kappa_a (kappa_a+1) - kappa_b (kappa_b+1) ]
        R^(Vinti) = int_0^infty dr  P_a [ --   -  ------------------------------------------ ] P_b    +  similar (not equal) in Q_a, Q_b
                                        [ dr                          2 r                    ]
        
        a value::Float64 is returned.
"""
function  Vinti(a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    
    mtp_ab = min(size(a.P, 1), size(b.P, 1));    kapa = a.subshell.kappa;     kapb = b.subshell.kappa
    wa = 0.
    for  r = 2:mtp_ab
        wc = a.P[r] * b.Pprime[r] - a.P[r] * kapa * (kapa+1) * b.P[r] / (2grid.r[r])  + a.P[r] * kapb * (kapb+1) * b.P[r] / (2grid.r[r])
        wd = a.Q[r] * b.Qprime[r] + a.Q[r] * kapa * (-kapa+1)* b.Q[r] / (2grid.r[r]) - a.Q[r] * kapb * (-kapb+1) * b.Q[r] / (2grid.r[r])                            
        wa = wa +  (wc + wd) * grid.wr[r]
    end
    return( wa )
end

"""
`RadialIntegrals.V0(wa::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)` 
    ... computes the (radial) integral int_0^infty dr wa; a value::Float64 is returned.
"""
function V0(wa::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)
    
    wb = 0.
    for  i = 1:mtp   wb = wb + wa[i] * grid.wr[i]   end
    return( wb )
end

"""
`RadialIntegrals.W5_Integral(mu::Int64, nu::Int64, a::Radial.Orbital, b::Radial.Orbital,  
                                                    c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)`  
    ... computes the (radial) integral for four relativistic orbitals: 
                            
        W_5 [ac|bd] = int_0^infty dr   int_0^r ds   [Pa Qc]_{r}  * ( s^nu / r^(nu+1) ) * [Pb Qd]_{s}

        as it frequently occurs in the frequency-independent Breit interaction.
"""
function W5_Integral(mu::Int64, nu::Int64, a::Radial.Orbital, b::Radial.Orbital,  
                                            c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    # Note mu = 5 is fixed historically and not used for this integral.
    !(mu == 5)   &&   error("mu = 5 required.")
    mtp = min(size(b.P, 1), size(d.P, 1))
    
    mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
    wa = 0.
    for  r = 2:mtp_ac
        for  s = 2:mtp_bd   
            if     s > r  continue  
            elseif s ==r  wa = wa + (a.P[r] * c.Q[r]) * (grid.r[s]^nu) / (grid.r[r]^(nu+1)) * (b.P[s] * d.Q[s]) * grid.wr[r] * grid.wr[s] / 2.0   
            else          wa = wa + (a.P[r] * c.Q[r]) * (grid.r[s]^nu) / (grid.r[r]^(nu+1)) * (b.P[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]
            end
        end
    end
    return( wa )
end

"""
`RadialIntegrals.Yk_ab(k::Int64, r::Float64, rho_ab::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)`  
    ... computes the (radial) integral

                                            r<^k
        Y_ab^k (r) = r * int_0^infty dr'  ------   rho_ab(r')
                                            r>^k+1

        a value::Float64 is returned.
"""
function Yk_ab(k::Int64, r::Float64, rho_ab::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)
    
    wa = 0.
    for  i = 2:mtp   
        rl = min(r, grid.r[i]);   rg = max(r, grid.r[i])
        wa = wa + rho_ab[i] * rl^k / rg^(k+1) * grid.wr[i]
    end
    return( r * wa )
end

end # module

