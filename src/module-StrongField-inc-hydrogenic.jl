
"""
`StrongField.hydrogenPnl(epsiloni::Float64, n::Int64, l::Int64, rGrid::Array{Float64,1})`  
    ... returns the non-relativistic hydrogen-like radial wave function with modified binding energy 
        epsiloni at the radial grid points rGrid.
"""
function  hydrogenPnl(epsiloni::Float64, n::Int64, l::Int64, rGrid::Array{Float64,1})
    Pnl = Float64[];    Z = n*sqrt(-2*epsiloni)
    
    for j = 1:length(rGrid)
        r = rGrid[j]
        p = r * sqrt( (2*Z/n)^3 * factorial(n-l-1)/(2*n*factorial(n+l)) ) * exp(-Z*r/n) * (2*Z*r/n)^l * 
            GSL.sf_laguerre_n(n-l-1,2*l+1,2*Z*r/n)
        push!( Pnl, p )
    end
    
    return( Pnl )
end


"""
`StrongField.hydrogenDPnlDr(epsiloni::Float64, n::Int64, l::Int64, rGrid::Array{Float64,1})`  
    ... returns the r-derivative of the non-relativistic hydrogen-like radial wave function with modified binding 
        nergy epsiloni at the radial grid points rGrid
"""
function  hydrogenDPnlDr(epsiloni::Float64, n::Int64, l::Int64, rGrid::Array{Float64,1})
    DPnl = Float64[];   Z = n*sqrt(-2*epsiloni)
    
    for j = 1:length(rGrid)
        r = rGrid[j]
        if  n-l-2 >= 0
            p = 1/n * 2^(l+1) * exp(-r*Z/n) * (r*Z/n)^l * sqrt( Z^3 * factorial(n-l-1) / (n^4*factorial(l+n)) ) *
                ( (n+l*n-r*Z) * GSL.sf_laguerre_n(n-l-1,2*l+1,2*r*Z/n) - 2*r*Z * GSL.sf_laguerre_n(n-l-2,2+2*l,2*r*Z/n) )
        elseif  n == 1  &&  l == 0
            p = -2 * exp(-r*Z) * sqrt(Z^3) * (r*Z - 1)
        else # not implemented
            p = 0
        end
        push!( DPnl, p )
    end
    
    return( DPnl )
end


"""
`StrongField.computeScalarProdBoundContHydrogenicUncoupled(epsilonp::Float64, lp::Int64, n::Int64, l::Int64, 
                                                            epsiloni::Float64, volkov::AbstractVolkovState, grid::Radial.Grid)`  
    ... computes the scalar product of the bound (hydrogenic) and continuum (Volkov) states in the one-particle picture.
"""
function  computeScalarProdBoundContHydrogenicUncoupled(epsilonp::Float64, lp::Int64, n::Int64, l::Int64, m::Int64, 
                                                        epsiloni::Float64, volkov::AbstractVolkovState, grid::Radial.Grid)
    # Compute <epsilonp lp||p||n l>        
    rvalues = grid.r;      rweights = grid.wr

    Pnl = hydrogenPnl( epsiloni, n, l, rvalues )
    if      typeof(volkov) == FreeVolkov        Pepsplp = VolkovP( epsilonp, lp, rvalues )
    elseif  typeof(volkov) == CoulombVolkov     Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rvalues )
    end

    noPoints = min(size(Pepsplp, 1), size(Pnl, 1))
    integral = ComplexF64(0.)

    #   Sum over grid and compute Gauss-Legendre sum
    for    k = 1:noPoints
        r = rvalues[k]
        integrand = conj( Pepsplp[k] ) * Pnl[k]
        # Gauss-Legendre sum
        integral  = integral + rweights[k] * integrand
    end

    return((-im)^l * integral)
end
    

"""
`StrongField.computeScalarProdBoundContHydrogenic(Pepsplp::Array{ComplexF64,1}, epsiloni::Float64, 
                                                    n::Int64, l::Int64, rGrid::Radial.Grid)`  
    ... computes the scalar product of the bound and continuum states in the one-particle picture for hydrogenic initial 
        states, and where the radial wave function Pepsplp of |epsilonp lp jp> needs to be provided as an argument 
        on the grid rGrid
"""
function  computeScalarProdBoundContHydrogenic(Pepsplp::Array{ComplexF64,1}, epsiloni::Float64, 
                                                n::Int64, l::Int64, rGrid::Radial.Grid)
    
    rgrid = rGrid.r;        orderGL = size(rgrid)[1];      weights = rGrid.wr
    
    #Prepare radial wave functions
    Pnl  = hydrogenPnl( epsiloni, n, l, rgrid )
    DPnl = hydrogenDPnlDr( epsiloni, n, l, rgrid )       
    
    integral = ComplexF64(0.)
    
    # Sum over grid and compute Gauss-Legendre sum
    for    k = 1:orderGL
        r = rgrid[k]
        integrand = conj( Pepsplp[k] ) * Pnl[k]
        # Gauss-Legendre sum
        integral = integral + weights[k] * integrand
    end
    
    return((-im)^l * integral)
end


"""
`StrongField.pReducedMEHydrogenic(Pepsplp::Array{ComplexF64,1} , lp::Int64, jp::Float64, epsiloni::Float64, 
                                    n::Int64, l::Int64, j::Float64, rGrid::Radial.Grid )`  
    ... computes the reduced matrix elements of the momentum operator <epsilonp lp jp||p||n l j> in the one-particle 
        picture for hydrogenic initial states, where the radial wave function Pepsplp of |epsilonp lp jp> needs to be 
        provided as an argument on the grid rGrid.
        
        FORMULATION in j-l-mj-BASIS:
"""
function  pReducedMEHydrogenic(Pepsplp::Array{ComplexF64,1}, lp::Int64, jp::Float64, epsiloni::Float64, 
                                n::Int64, l::Int64, j::Float64, rGrid::Radial.Grid )

    # Compute the factor in <epsilonp lp jp||p||n l j> = fac * <epsilonp lp||p||n l>
    fac = 0.
    
    for    kms = -1:0   # ms = 1/2 + kms
        ms = 1/2 + kms
        for    ml = -l:l
            mj = ml + ms 
            for    mlp = -lp:lp
                mjp = mlp + ms 
                fac = fac + AngularMomentum.ClebschGordan_old(  AngularJ64(lp), AngularM64(mlp),  AngularJ64(1//2),           
                                            AngularM64(Rational(ms)), AngularJ64(Rational(jp)),   AngularM64(Rational(mjp)) ) * 
                            AngularMomentum.ClebschGordan_old(  AngularJ64(l),  AngularM64(ml),   AngularJ64(1//2),           
                                            AngularM64(Rational(ms)), AngularJ64(Rational(j)),    AngularM64(Rational(mj)) )    
            end
        end
    end

    
    # Compute <epsilonp lp||p||n l>
    rgrid = rGrid.r;      orderGL = size(rgrid)[1];     weights = rGrid.wr
    
    Pnl  = hydrogenPnl( epsiloni, n, l, rgrid )
    DPnl = hydrogenDPnlDr( epsiloni, n, l, rgrid )
    
    integral = ComplexF64(0.)
    
    # Sum over grid and compute Gauss-Legendre sum
    for    k = 1:orderGL
        r = rgrid[k]
        integrand = conj( Pepsplp[k] )/r * ( r*DPnl[k] - ((lp-l)*(lp+l+1))/2 * Pnl[k] )
        # Gauss-Legendre sum
        integral = integral + weights[k] * integrand
    end

    # Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
    integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )
    
    return( fac * integral )
end


"""
`StrongField.pReducedMEHydrogenicUncoupled(epsilonp::Float64, lp::Int64, n::Int64, l::Int64, epsiloni::Float64, 
                                            volkov::AbstractVolkovState, grid::Radial.Grid)`  
    ... computes the reduced matrix elements of the momentum operator <epsilonp lp ||p||n l> in the 
        one-particle picture in the uncoupled basis.
"""
function  pReducedMEHydrogenicUncoupled(epsilonp::Float64, lp::Int64, n::Int64, l::Int64, epsiloni::Float64, 
                                        volkov::AbstractVolkovState, grid::Radial.Grid)
    # Compute <epsilonp lp||p||n l>     
    rvalues  = grid.r;      rweights = grid.wr

    Pnl = hydrogenPnl( epsiloni, n, l, rvalues )
    if      typeof(volkov) == FreeVolkov        Pepsplp = VolkovP( epsilonp, lp, rvalues )
    elseif  typeof(volkov) == CoulombVolkov     Pepsplp = CoulombVolkovP( epsilonp, lp, volkov.Z, rvalues )
    end
    DPnl = HydrogenDPnlDr( epsiloni, n, l, rvalues )

    noPoints = min(size(Pepsplp, 1), size(Pnl, 1))
    integral = ComplexF64(0.)

    # Sum over grid and compute Gauss-Legendre sum
    for    k = 1:noPoints
        r = rvalues[k]
        integrand = conj( Pepsplp[k] )/r * ( r*DPnl[k] - ((lp-l)*(lp+l+1))/2 * Pnl[k] )
        # Gauss-Legendre sum
        integral = integral + rweights[k] * integrand
    end

    # Note that GSL.sf_coupling_3j takes takes the input (2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
    integral = integral * (-im)^(lp+1) * (-1)^lp * GSL.sf_coupling_3j( 2*lp, 2*1, 2*l, 0, 0, 0 )

    return( integral )
end

