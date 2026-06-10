
"""
`module JAC.RadialIntegrals`  
... a submodel of JAC that contains methods for calculating radial one- and two-particle matrix elements. These integrals occur 
    frequently in atomic structure and collision theory, and their fast computations often appears essential.
"""
module  RadialIntegrals

using  GSL, QuadGK
using  ..AngularMomentum, ..Basics, ..BsplinesN, ..Defaults,  ..Radial, ..Math, ..ManyElectron, ..Nuclear


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
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f1(i::Int64)   dPb(j) = Math.derivative(b.P, j)
                                dQb(j) = Math.derivative(b.Q, j)
                                return( a.Q[i] * dPb(i) - a.P[i] * dQb(i) )                          end
        function f2(i::Int64)   return( (a.Q[i] * b.P[i] + a.P[i] * b.Q[i]) / grid.r[i] )            end
        function f3(i::Int64)   return( a.Q[i] * b.Q[i] )                                            end
        function f4(i::Int64)   return( -Zr[i] * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) / grid.r[i] )   end

        I1 = Math.integrateFit(f1, min(size(a.P, 1), size(b.P, 1)), grid) / grid.h
        I2 = Math.integrateFitTransform(f2, min(size(a.P, 1), size(b.P, 1)), grid)
        I3 = Math.integrateFitTransform(f3, min(size(a.P, 1), size(b.P, 1)), grid)
        I4 = Math.integrateFitTransform(f4, min(size(a.P, 1), size(b.P, 1)), grid)

        return( Defaults.INVERSE_FINE_STRUCTURE_CONSTANT * I1 + Defaults.INVERSE_FINE_STRUCTURE_CONSTANT * kappa * I2 -
                2 * Defaults.INVERSE_FINE_STRUCTURE_CONSTANT^2 * I3 + I4 )
                
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wa = wa + grid.wr[i] * (  wc * a.Q[i] * (b.Pprime[i] + kappa/grid.r[i] * b.P[i])  
                                    - wc * a.P[i] * (b.Qprime[i] - kappa/grid.r[i] * b.Q[i]) 
                                    - 2wc^2 * a.Q[i] * b.Q[i]
                                    - Zr[i] * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) / grid.r[i]  ) 
        end
        return( wa )
    else
        error("stop a")
    end
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
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wa = wa + grid.wr[i] * (  wc * a.Q[i] * (b.Pprime[i] + kappa/grid.r[i] * b.P[i])  
                                    - wc * a.P[i] * (b.Qprime[i] - kappa/grid.r[i] * b.Q[i]) 
                                    - 2wc^2 * a.Q[i] * b.Q[i]
                                    - Zr[i] * (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) / grid.r[i]  ) * exp(-tau * grid.r[i])
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.GrantILminus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^- (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b - Q_a P_b ] .
"""
function GrantILminus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) )       end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i] - a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.GrantILplus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^+ (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b + Q_a P_b ] .
"""
function GrantILplus(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) )       end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.GrantIL0(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        I_L^0 (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a Q_b ] .
"""
function GrantIL0(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) )       end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.GrantJL(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        J_L (q; a,b) = int_0^\\infty dr j_L (qr) [ P_a P_b + Q_a Q_b ] .
"""
function GrantJL(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) )       end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * GSL.sf_bessel_jl(L, q * grid.r[i]) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.GrantJL_cp(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, cp::CorePolarization)`  
    ... computes Grant's (radial) integral for two relativistic orbitals:  
        J_L (q; a,b) = int_0^\\infty dr   [ P_a P_b + Q_a Q_b ] * r * [1 - coreAlpha / (r^2 + coreRadius^2)^3/2 ].
"""
function GrantJL_cp(L::Int64, q::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid, cp::CorePolarization)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * grid.wr[i] * grid.r[i] * 
                                    (1.0 - cp.coreAlpha / (grid.r[i]^2 + cp.coreRadius^2)^(3/2) )   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.isotope_boson(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the boson-field shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_boson(a::Orbital, b::Orbital, potential::Array{Float64,1}, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1), size(potential, 1));   
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wa = wa + (a.P[i] * b.P[i]  +  a.Q[i] * b.Q[i]) * potential[i] * grid.wr[i]
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.isotope_field(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the field-shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_field(a::Orbital, b::Orbital, deltaPotential::Array{Float64,1}, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1), size(deltaPotential, 1));   
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wa = wa - (a.P[i] * b.P[i]  +  a.Q[i] * b.Q[i]) * deltaPotential[i] / grid.r[i] * grid.wr[i]
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.isotope_nms(a::Orbital, b::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the normal mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_nms(a::Orbital, b::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1));   lb = Basics.subshell_l(b.subshell);    jb2   = Basics.subshell_2j(b.subshell)
    alphaZ = Defaults.getDefaults("alpha") * Z
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wb = (a.Pprime[i] * b.Pprime[i]  +  a.Qprime[i] * b.Qprime[i])  + 
                    (lb*(lb+1) * a.P[i] * b.P[i]  +  (jb2-1)*jb2 * a.Q[i] * b.Q[i]) / (grid.r[i]^2)
            wc = - 2 * alphaZ * (a.Q[i] * b.Pprime[i]  +  b.Q[i] * a.Pprime[i]) / grid.r[i]
            wd = - alphaZ * (b.subshell.kappa - 1) * (a.Q[i] * b.P[i]  +  b.Q[i] * a.P[i]) / (grid.r[i]^2)
            wa = wa + (wb + wc + wd) * grid.wr[i]   
        end
        return( wa / 2. )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.isotope_smsB(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the specific mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_smsB(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(c.P, 1));   alphaZ = Defaults.getDefaults("alpha") * Z
    kapa = a.subshell.kappa;   mkapa = -kapa;    kapc = c.subshell.kappa;   mkapc = -kapc 
    
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wb = (- a.Q[i] * c.P[i] * AngularMomentum.sigma_reduced_me_ma(mkapa, kapc)  +
                    c.Q[i] * a.P[i] * AngularMomentum.sigma_reduced_me_mb(kapa,  mkapc)  )
            wa = wa - alphaZ / grid.r[i] * wb * grid.wr[i]   
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.isotope_smsC(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)`  
    ... computes the specific mass shift radial integral int_o^infty ... A value::Float64 is returned.
"""
function isotope_smsC(a::Orbital, c::Orbital, Z::Float64, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(c.P, 1));   alphaZ = Defaults.getDefaults("alpha") * Z
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wa = wa - alphaZ / grid.r[i] * (a.Q[i] * c.P[i] - c.Q[i] * a.P[i]) * grid.wr[i]   
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.nondiagonalD(pm::Int64, kappa::Int64, bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)`  
    ... computes the (radial and non-diagonal) D_kappa^+/- integral two the bsplines, all defined on grid
        <bspline1| +/- d/dr + kappa/r | bspline2>. -- pm = +1/-1 provides the phase for taking the derivative.
"""
function nondiagonalD(pm::Int64, kappa::Int64, bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid) 
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    mtp = min( bspline1.upper, bspline2.upper)
    n0  = max( bspline1.lower, bspline2.lower)
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f1(i::Int64)
            wa = pm * bspline1.bs[i] * bspline2.bp[i]
            return( wa )                            
        end
        function f2(i :: Int64)
            if  i == 1  wa = bspline1.bs[i] * kappa * bspline2.bs[i] / (0.3 * grid.r[2]) 
            else        wa = bspline1.bs[i] * kappa * bspline2.bs[i] / grid.r[i]   end
            return( wa )
        end

        I1 = Math.integrateTransform(f1, n0, mtp, grid)
        I2 = Math.integrateTransform(f2, n0, mtp, grid)
        return( I1+I2 )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = n0:mtp  
            wa = wa + pm * bspline1.bs[i] * bspline2.bp[i] * grid.wr[i] 
            if  i == 1  wa = wa + bspline1.bs[i] * kappa * bspline2.bs[i] / (0.3 * grid.r[2]) * grid.wr[i] 
            else        wa = wa + bspline1.bs[i] * kappa * bspline2.bs[i] / grid.r[i] * grid.wr[i]   end
        end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.overlap()`

+ (orbital1::Radial.Orbital, orbital2::Radial.Orbital, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral <orbital_a|orbital_b>  for two relativistic orbitals of the same 
        symmetry (kappa).
"""
function overlap(orbital1::Radial.Orbital, orbital2::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(orbital1.P, 1), size(orbital2.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f(i :: Int64)
            return( orbital1.P[i] * orbital2.P[i] + orbital1.Q[i] * orbital2.Q[i] )
        end

        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 1:grid.NoPoints 
            if i > mtp   break   end
            wa = wa + ( orbital1.P[i] * orbital2.P[i] + orbital1.Q[i] * orbital2.Q[i] ) * grid.wr[i]   
        end
        return( wa )
    else
        error("stop a")
    end
end


"""
+ (bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral <bspline1|bsplines>  for two bpslines as defined on grid.
"""
function overlap(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
        mtp = min( bspline1.upper, bspline2.upper)
        n0  = max( bspline1.lower, bspline2.lower)

        function f(i :: Int64)
            return( bspline1.bs[i] * bspline2.bs[i] )
        end

        return( Math.integrateTransform(f, n0, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 1:grid.NoPoints   wa = wa + bspline1.bs[i] * bspline2.bs[i] * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
+ (p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral of two (non-relativistic) radial orbital functions <p1|p2>  as defined on grid.
"""
function overlap(p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)
    
    mtp = min( length(p1List), length(p2List))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f(i :: Int64)
            return( p1List[i] * p2List[i] )
        end

        return( Math.integrateTransform(f, 1, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 1:grid.NoPoints 
            if i > mtp   break   end
            wa = wa + p1List[i] * p2List[i] * grid.wr[i]   
        end
        return( wa )
    else
        error("stop b")
    end
end


function overlap_old2022(p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)
    
    mtp = min( length(p1List), length(p2List))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f(i :: Int64)
            return( p1List[i] * p2List[i] )
        end

        return( Math.integrateTransform(f, 1, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        error("stop a")
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.qedDampedOverlap(lambda::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)` 
    ... computes the damped (radial) integral  int_0^infty (P_a P_b  +  Q_a Q_b) * e^{r/lambda} for the radial 
        orbitals a, b on the given grid. A value::Float64 is returned.
"""
function qedDampedOverlap(lambda::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wb = Base.MathConstants.e^(- grid.r[i]/lambda);     wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.qedLowFrequency(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)` 
    ... computes the (radial) integral for the low-frequency QED potential for the radial orbitals a, b on the given grid. 
        A value::Float64 is returned.
"""
function qedLowFrequency(a::Radial.Orbital, b::Radial.Orbital, nm::Nuclear.Model, grid::Radial.Grid, qgrid::Radial.GridGL)
    alpha = Defaults.getDefaults("alpha");    BZ = 0.074 + 0.035 * nm.Z * alpha
    mtp = min(size(a.P, 1), size(b.P, 1))
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 1:mtp   wb = Base.MathConstants.e^(-nm.Z * grid.r[i]) ;     wa = wa + (a.P[i]*wb*b.P[i] + a.Q[i]*wb*b.Q[i]) * grid.wr[i]   end
        wa = -BZ * nm.Z^4 * alpha^3 * wa
        else
        error("stop a")
    end
    
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
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGL()
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
    else
        error("stop a")
    end
    
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
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            wb = tIntegral(grid.r[i]) 
            wc = (-pot.Zr[i] / grid.r[i])
            wa = wa + (a.P[i]*wb*wc*b.P[i] + a.Q[i]*wb*wc*b.Q[i]) * grid.wr[i]   
        end
        wa = 2. * Defaults.getDefaults("alpha") / (3pi) * wa
    else
        error("stop a")
    end
    
    println("QED single-electron strength <$(a.subshell)| h^(simplified Uehling) | $(b.subshell)> = $wa ")
    return( wa )
end


"""
`RadialIntegrals.rkDiagonal()`   ... computes the (radial and diagonal) integral of r^k for two radial orbital functions.

+ (k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`  
    ... computes this integral for two relativistic orbitals:   < r^k >_ab = int_0^\\infty  dr  [P_a P_b + Q_a Q_b]  r^k
"""
function rkDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * (grid.r[i]^k) )    end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        if  k > -3   m0 = 2   else   m0 = 6   end    # Don't allow too small r-values
        for  i = m0:mtp   wa = wa + (a.P[i] * b.P[i] + a.Q[i] * b.Q[i]) * (grid.r[i]^k) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
+ (k::Int64, p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)`  
    ... computes this integral for two non-relativistic orbitals:   < r^k >_ab = int_0^\\infty  dr  [P_a P_b]  r^k
"""
function rkDiagonal(k::Int64, p1List::Array{Float64,1}, p2List::Array{Float64,1}, grid::Radial.Grid)
    mtp = min( length(p1List), length(p2List))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f(i :: Int64)
            return( p1List[i] * p2List[i] * (grid.r[i]^k) )
        end

        return( Math.integrateTransform(f, 2, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + p1List[i] * p2List[i] * (grid.r[i]^k) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.rkNonDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)` 
    ... computes the (radial and non-diagonal) integral of r^k for two relativistic orbitals:
        [ r^k ]_ab = int_0^\\infty dr [P_a Q_b + Q_a P_b] r^k
"""
function rkNonDiagonal(k::Int64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P, 1), size(b.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * (grid.r[i]^k) )    end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   wa = wa + (a.P[i] * b.Q[i] + a.Q[i] * b.P[i]) * (grid.r[i]^k) * grid.wr[i]   end
        return( wa )
    else
        error("stop a")
    end
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

    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(Ba, 1), size(Bc, 1));    mtp_bd = min(size(Bb, 1), size(Bd, 1))
        wac = zeros(mtp_ac);   wbd = zeros(mtp_bd)
        for  r = 2:mtp_ac   wac[r] = (Ba[r] * Bc[r]) * grid.wr[r]  end
        for  s = 2:mtp_bd   wbd[s] = (Bb[s] * Bd[s]) * grid.wr[s]  end
        wa = 0.
        for  r = 2:mtp_ac
            for  s = 2:mtp_bd   wa = wa + wac[r] * ul(grid.r[r], grid.r[s]) * wbd[s]   end
        end
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.SlaterRk_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the (relativistic) Slater integral

    R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d)

    of rank k for the four orbitals a, b, c, d, and over the given grid by using an explicit 2-dimensional integration scheme; a 
    value::Float64 is returned.
"""
function SlaterRk_2dim(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( r^k/s^(k+1) )
        elseif r > s     return( s^k/r^(k+1) )
        end
    end

    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function fs(r :: Int64, s :: Int64) :: Float64
            return( ul(grid.r[r], grid.r[s]) * ( b.P[s] * d.P[s] + b.Q[s] * d.Q[s] ) )
        end

        function f(r :: Int64) :: Float64
        function ff(i :: Int64) :: Float64    return( fs(r, i) )    end
        return( (a.P[r] * c.P[r] + a.Q[r] * c.Q[r] ) * Math.integrateFitTransform(ff, min(size(b.P, 1), size(d.P, 1)), grid) )
        end

        return Math.integrateFitTransform(f, min(size(a.P, 1), size(c.P, 1)), grid)
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
        wac = zeros(mtp_ac);   wbd = zeros(mtp_bd)
        for  r = 2:mtp_ac   wac[r] = (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * grid.wr[r]  end
        for  s = 2:mtp_bd   wbd[s] = (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[s]  end
        wa = 0.
        for  r = 2:mtp_ac
            for  s = 2:mtp_bd   wa = wa + wac[r] * ul(grid.r[r], grid.r[s]) * wbd[s]   end
        end
        return( wa )
    else
        error("stop a")
    end
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

    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function fs(r :: Int64, s :: Int64) :: Float64
            return( ul(grid.r[r], grid.r[s]) * ( b.P[s] * d.P[s] + b.Q[s] * d.Q[s] ) )
        end

        function f(r :: Int64) :: Float64
        function ff(i :: Int64) :: Float64    return( fs(r, i) )    end
        return( (a.P[r] * c.P[r] + a.Q[r] * c.Q[r] ) * Math.integrateFitTransform(ff, min(size(b.P, 1), size(d.P, 1)), grid) )
        end

        return Math.integrateFitTransform(f, min(size(a.P, 1), size(c.P, 1)), grid)
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
        wa = 0.
        for  r = 2:mtp_ac
            for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                            (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]   end
        end
        ## println("Test: SlaterRk_2dim(); wa = $wa")
        return( wa )
    else
        error("stop a")
    end
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

    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
        wa = 0.
        for  r = 2:mtp_ac
            for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                            (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]  *
                                            exp(- tau * grid.r[r] - tau * grid.r[s] )                         end
        end
        ## println("Test: SlaterRk_2dim(); wa = $wa")
        return( wa )
    else
        error("stop b")
    end
end


"""
`RadialIntegrals.SlaterRk_new(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Orbital, grid::Radial.Grid)`  
    ... computes the (relativistic) Slater integral

        R^k (abcd) = int_0^infty dr int_0^infty ds (P_a P_c + Q_a Q_c) r_<^k / r_>^(k+1) (P_b P_d + Q_b Q_d)

        of rank k for the four orbitals a, b, c, d, and over the given grid by using a factorized integration scheme; a 
        value::Float64 is returned.
"""
function SlaterRk_new(k::Int64, a::Radial.Orbital, b::Radial.Orbital, c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    function ul(r :: Float64, s :: Float64) :: Float64
        if     r <= s    return( r^k/s^(k+1) )
        elseif r > s     return( s^k/r^(k+1) )
        end
    end
    
    function rLowUp()  
        rtuple = Tuple{Int64,Int64}[]
        for  n = 2:grid.orderGL:100000000
            if   n > grid.NoPoints   break   end
            push!( rtuple, (n, n+grid.orderGL-1))
        end
        return( rtuple )
    end
            
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1));  rtuple = rLowUp()
        wa = 0.
        for  (rlow, rup)  in  rtuple
            for  (slow, sup)  in  rtuple
                war = 0.;    was = 0.
                if  sup < rlow
                    for  r = rlow:rup   if r > mtp_ac   break   end
                        war = war + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) / grid.r[r]^(k+1) * grid.wr[r]    end
                    for  s = slow:sup   if s > mtp_bd   break   end
                        was = was + (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.r[s]^k     * grid.wr[s]    end
                    wa = wa + war * was
                elseif  rup < slow
                    for  r = rlow:rup   if r > mtp_ac   break   end   
                        war = war + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * grid.r[r]^k     * grid.wr[r]    end
                    for  s = slow:sup   if s > mtp_bd   break   end   
                        was = was + (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) / grid.r[s]^(k+1) * grid.wr[s]    end
                    wa = wa + war * was
                elseif  rlow == slow  &&  rup == sup
                    for  r = rlow:rup       if r > mtp_ac   break   end   
                        for  s = slow:sup   if s > mtp_bd   break   end   
                            wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                        (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]       end
                    end
                else  error("stop b")
                end
            end
        end
        return( wa )
    else
        error("stop a")
    end
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

    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function fs(r :: Int64, s :: Int64) :: Float64
            return( ul(grid.r[r], grid.r[s]) * ( b.P[s] * d.P[s] + b.Q[s] * d.Q[s] ) )
        end

        function f(r :: Int64) :: Float64
        function ff(i :: Int64) :: Float64    return( fs(r, i) )    end
        return( (a.P[r] * c.P[r] + a.Q[r] * c.Q[r] ) * Math.integrateFitTransform(ff, min(size(b.P, 1), size(d.P, 1)), grid) )
        end

        return Math.integrateFitTransform(f, min(size(a.P, 1), size(c.P, 1)), grid)
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ac = min(size(a.P, 1), size(c.P, 1));    mtp_bd = min(size(b.P, 1), size(d.P, 1))
        wa = 0.
        for  r = 2:mtp_ac
            for  s = 2:mtp_bd   wa = wa + (a.P[r] * c.P[r] + a.Q[r] * c.Q[r]) * ul(grid.r[r], grid.r[s]) * 
                                            (b.P[s] * d.P[s] + b.Q[s] * d.Q[s]) * grid.wr[r] * grid.wr[s]   end
        end
        return( wa )
    else
        error("stop a")
    end
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
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        error("stop a")
    elseif  grid.meshType == Radial.MeshGL()
        mtp_ab = min(size(a.P, 1), size(b.P, 1));    kapa = a.subshell.kappa;     kapb = b.subshell.kappa
        wa = 0.
        for  r = 2:mtp_ab
            wc = a.P[r] * b.Pprime[r] - a.P[r] * kapa * (kapa+1) * b.P[r] / (2grid.r[r])  + a.P[r] * kapb * (kapb+1) * b.P[r] / (2grid.r[r])
            wd = a.Q[r] * b.Qprime[r] + a.Q[r] * kapa * (-kapa+1)* b.Q[r] / (2grid.r[r]) - a.Q[r] * kapb * (-kapb+1) * b.Q[r] / (2grid.r[r])                            
            wa = wa +  (wc + wd) * grid.wr[r]
        end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.Vlocal(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, potential::Radial.Potential, grid::Radial.Grid)`  
    ... computes the (radial) integral for the local potential and for two bpslines, all defined on grid 
        <bspline1| potential.Zr | bspline2>. -- Here, potential. V must provide the effective charge zz(r) = - V * r.
"""
function Vlocal(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, potential::Radial.Potential, grid::Radial.Grid)
    ## if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    mtp = min( bspline1.upper, bspline2.upper)
    n0  = max( bspline1.lower, bspline2.lower)
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()

        function f(i :: Int64)
            if  i == 1  wa = - bspline1.bs[i] * potential.Zr[i] * bspline2.bs[i] / (0.3 * grid.r[2]) 
            else        wa = - bspline1.bs[i] * potential.Zr[i] * bspline2.bs[i] / grid.r[i]   end
            return( wa )
        end

        return( Math.integrateTransform(f, n0, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = n0:mtp  
            if  i == 1  wa = wa - bspline1.bs[i] * potential.Zr[i] * bspline2.bs[i] / (0.3 * grid.r[2]) * grid.wr[i] 
            else        wa = wa - bspline1.bs[i] * potential.Zr[i] * bspline2.bs[i] / grid.r[i] * grid.wr[i]   end
        end
        return( wa )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.V0(wa::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)` 
    ... computes the (radial) integral int_0^infty dr wa; a value::Float64 is returned.
"""
function V0(wa::Array{Float64,1}, mtp::Int64, grid::Radial.Grid)
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    return( wa[i] )     end
        return( Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wb = 0.
        for  i = 1:mtp   wb = wb + wa[i] * grid.wr[i]   end
        return( wb )
    else
        error("stop a")
    end
end


"""
`RadialIntegrals.W5_Integral(mu::Int64, nu::Int64, orbitala::Radial.Orbital, orbitalc::Radial.Orbital,  
                                                    orbitalb::Radial.Orbital, orbitald::Radial.Orbital, grid::Radial.Grid)`  
    ... computes the (radial) integral for four relativistic orbitals: 
                            
        W_5 [ac|bd] = int_0^infty dr   int_0^r ds   [Pa Qc]_{r}  * ( s^nu / r^(nu+1) ) * [Pb Qd]_{s}

        as it frequently occurs in the frequency-independent Breit interaction.
"""
function W5_Integral(mu::Int64, nu::Int64, a::Radial.Orbital, b::Radial.Orbital,  
                                            c::Radial.Orbital, d::Radial.Orbital, grid::Radial.Grid)
    # Note mu = 5 is fixed historically and not used for this integral.
    !(mu == 5)   &&   error("mu = 5 required.")
    mtp = min(size(b.P, 1), size(d.P, 1))
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function fs(s :: Int64) :: Float64
            return( grid.r[s]^nu * ( b.P[s] * d.Q[s] ) )
        end

        function f(r :: Int64) :: Float64
        function ff(i :: Int64) :: Float64  return( Math.integrateFitTransform(fs, i, grid) )  end
            if  r > mtp   return( 0. )
            else          return(  a.P[r] * c.Q[r] * ff(r) / grid.r[r]^(nu+1) )
            end
        end

        return Math.integrateFitTransform(f, min(size(a.P, 1), size(c.P, 1)), grid)
    elseif  grid.meshType == Radial.MeshGL()
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
    else
        error("stop a")
    end
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
    
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGrasp()
        function f(i :: Int64)    rl = min(r, grid.r[i]);   rg = max(r, grid.r[i]);   return( rho_ab[i] * rl^k / rg^(k+1) )     end
        return( r * Math.integrateFitTransform(f, mtp, grid) )
    elseif  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp   
            rl = min(r, grid.r[i]);   rg = max(r, grid.r[i])
            wa = wa + rho_ab[i] * rl^k / rg^(k+1) * grid.wr[i]
        end
        return( r * wa )
    else
        error("stop a")
    end
end

end # module

