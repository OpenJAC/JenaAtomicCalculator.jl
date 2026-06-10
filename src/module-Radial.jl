
"""
`module  JAC.Radial`  
... a submodel of JAC that contains all structs and methods to deal with the radial grid, radial orbitals and potentials.
"""
module Radial


using  QuadGK, Printf, JenaAtomicCalculator, ..Basics 

export Grid, Potential, Orbital, SingleElecSpectrum


"""
`abstract type Radial.AbstractMesh` 
    ... defines an abstract mesh type for the use of the radial mesh
    
    + struct MeshNone       ... if an empty instance of Radial.Grid() need to be used.
    + struct MeshGrasp      ... to use a GRASP-like grid with Newton-Cortes integration rules.
    + struct MeshGL         ... to use a Gauss-Legendre type grid where the mesh points between the break points of the t-grid
                                are chosen to support a Gauss-Legendre integration in this potential.
"""
abstract type  AbstractMesh                     end
struct   MeshNone   <:  Radial.AbstractMesh     end
struct   MeshGrasp  <:  Radial.AbstractMesh     end
struct   MeshGL     <:  Radial.AbstractMesh     end


"""
`struct  Radial.Grid`  ... defines a type for the radial grid which contains all information about the grid parameters, the genration 
                            of the B-spline basis as well as for performing radial integrations.

    ** Physical grid parameter **
    + rnt        ::Float64           ... smalles grid point > 0.
    + h          ::Float64           ... stepsize in the construction of the exponential grid.
    + hp         ::Float64           ... asymptotic stepsize of the log-lin grid.
    + NoPoints   ::Int64             ... No. of grid points so that r[NoPoints] coincides also 
                                            with the largest break point of the B-spline knot.
    ** B-spline grid parameters and break points **
    + tL         ::Array{Float64,1}  ... radial break points for the B-splines of the large c.
    + tS         ::Array{Float64,1}  ... radial break points for the B-splines of the small c.
    + ntL        ::Int64             ... number of break points in the t-grid of the large c.
    + ntS        ::Int64             ... number of break points in the t-grid of the small c.
    + orderL     ::Int64             ... B-spline order of large components.
    + orderS     ::Int64             ... B-spline order of small components.
    + nsL        ::Int64             ... number of B-splines for large components.
    + nsS        ::Int64             ... number of B-splines for small components.
    + orderGL    ::Int64             
        ... order of the Gauss-Lengedre integration if mesh == Radial.MeshGL(); this order also determines
            the (number of) break points by taking the orderGL-th point from the physical grid points.
    ** Radial mesh points **
    + meshType   ::Radial.AbstractMesh
    + r          ::Array{Float64,1}  ... radial grid points
    + rp         ::Array{Float64,1}  ... derivative of the radial grid at the grid points
    + rpor       ::Array{Float64,1}  ... rp over r
    + wr         ::Array{Float64,1}  
        ... integration weights for all grid points, for instance, GL weights.
"""
struct  Grid
    rnt          ::Float64
    h            ::Float64
    hp           ::Float64
    NoPoints     ::Int64
    #
    tL           ::Array{Float64,1}
    tS           ::Array{Float64,1}
    ntL          ::Int64 
    ntS          ::Int64 
    orderL       ::Int64
    orderS       ::Int64
    nsL          ::Int64
    nsS          ::Int64
    orderGL      ::Int64
    #
    meshType     ::Radial.AbstractMesh
    r            ::Array{Float64,1}
    rp           ::Array{Float64,1}
    rpor         ::Array{Float64,1}
    wr           ::Array{Float64,1}
end


"""
`Radial.Grid()` ... constructor to define an 'empty' grid.  
"""
function Grid()
    Radial.Grid(0., 0., 0., 0,   Float64[], Float64[], 0, 0, 0, 0, 0, 0, 0,   Radial.MeshNone(), Float64[], Float64[], Float64[], Float64[] )
end


"""
`Radial.Grid(exponential::Bool; printout::Bool=false)` 
    ... constructor to define either a standard 'exponential' mesh (true) or a 'log-lin' mesh (false). 
        The standard mesh type is Radial.MeshGL().
"""
function Grid(exponential::Bool; printout::Bool=false)
    # The standard mesh type is Radial.MeshGL() with some proper parameters
    nth = orderL = 7;    orderS = 8;   orderGL = 7;    meshType = MeshGL()
    
    if  exponential
        NoPoints = 392;    NoPoints = NoPoints - rem(NoPoints,orderGL)
        grid = Radial.Grid(2.0e-6, 5.0e-2,     0., NoPoints,   Float64[], Float64[], 0, 0, orderL, orderS, 0, 0, orderGL, meshType,
                            Float64[], Float64[], Float64[], Float64[] )
    else
        NoPoints = 600;    NoPoints = NoPoints - rem(NoPoints,orderGL)
        grid = Radial.Grid(2.0e-6, 5.0e-2, 2.0e-2, NoPoints,   Float64[], Float64[], 0, 0, orderL, orderS, 0, 0, orderGL, meshType,
                            Float64[], Float64[], Float64[], Float64[] )
    end
    return( Radial.determineGrid(grid, printout=printout) )
end


"""
`Radial.Grid(gr::Radial.Grid;`

        rnt=..,     h=..,       hp=..,      rbox=..,    orderL=..,  orderS=..,  orderGL=..,    
        meshType..,  printout=..)
    ... constructor for modifying the given Radial.Grid by 'overwriting' the previously selected parameters.
"""
function Grid(gr::Radial.Grid;    
    rnt::Union{Nothing,Float64}=nothing,        h::Union{Nothing,Float64}=nothing,      hp::Union{Nothing,Float64}=nothing, 
    rbox::Union{Nothing,Float64}=nothing,       nth::Union{Nothing,Int64}=nothing,      orderL::Union{Nothing,Int64}=nothing,
    orderS::Union{Nothing,Int64}=nothing,       orderGL::Union{Nothing,Int64}=nothing,  meshType::Union{Nothing,Radial.AbstractMesh}=nothing, 
    printout::Bool=false)
    
    if  rnt      == nothing   rntx      = gr.rnt        else    rntx      = rnt       end 
    if  h        == nothing   hx        = gr.h          else    hx        = h         end 
    if  hp       == nothing   hpx       = gr.hp         else    hpx       = hp        end 
    if  rbox     == nothing   rboxx     = nothing       else    rboxx     = rbox      end 
    if  orderL   == nothing   orderLx   = gr.orderL     else    orderLx   = orderL    end 
    if  orderS   == nothing   orderSx   = gr.orderS     else    orderSx   = orderS    end 
    if  orderGL  == nothing   orderGLx  = gr.orderGL    else    orderGLx  = orderGL   end 
    if  meshType == nothing   meshTypex = gr.meshType   else    meshTypex = meshType  end 
    
    if      rboxx == nothing    NoPointsx = gr.NoPoints - rem(gr.NoPoints, orderGLx)
    elseif  rboxx  > 0.         NoPointsx = Radial.determineNoPoints(rntx, hx, hpx, rboxx, orderGLx)
    else    error("stop a")
    end
    
    grid = Radial.Grid(rntx, hx, hpx, NoPointsx, Float64[], Float64[], 0, 0, orderLx, orderSx, 0, 0, orderGLx, meshTypex,
                        Float64[], Float64[], Float64[], Float64[] )
    return( Radial.determineGrid(grid, printout=printout) )
end


# `Base.show(io::IO, grid::Radial.Grid)`  ... prepares a proper printout of the variable grid::Radial.Grid.
function Base.show(io::IO, grid::Radial.Grid) 

    sa = Base.string(grid::Radial.Grid);    print(io, sa * "\n")

    nr = grid.NoPoints;    ntL = length(grid.tL);    ntS = length(grid.tS);    if  nr < 6   return( nothing )    end
    print(io, "r:    ", grid.r[1:3],    "  ...  ", grid.r[nr-2:nr],         "\n") 
    if  grid.meshType == Radial.MeshGrasp()
        print(io, "rp:   ", grid.rp[1:3],   "  ...  ", grid.rp[nr-2:nr],    "\n") 
        print(io, "rpor: ", grid.rpor[1:3], "  ...  ", grid.rpor[nr-2:nr],  "\n") 
    end
    print(io, "wr:   ", grid.wr[1:3],   "  ...  ", grid.wr[nr-2:nr],        "\n") 
    ## print(io, "tL:   ", grid.tL[1:3],   "  ...  ", grid.tL[ntL-2:ntL],   "\n") 
    print(io, "tS:   ", grid.tS[1:3],   "  ...  ", grid.tS[ntS-2:ntS]           ) 
end


# `Base.string(grid::Radial.Grid)`  ... provides a String notation for the variable grid::Radial.Grid.
function Base.string(grid::Radial.Grid) 
    if   grid.NoPoints == 0
        sa = "Radial grid not defined;  NoPoints = $(grid.NoPoints) ..."
    else
        sa = "Radial grid:  rnt = $(grid.rnt),  h = $(grid.h),  hp = $(grid.hp),  NoPoints = $(grid.NoPoints),  "
        sa = sa * "ntL = $(grid.ntL),  ntS = $(grid.ntS), "
        sa = sa * "orderL = $(grid.orderL),  orderS = $(grid.orderS),  nsL = $(grid.nsL),  nsS = $(grid.nsS),  mesh = $(grid.meshType), ...  "
    end 
end


"""
`Radial.determineGrid(grid::Radial.Grid; printout::Bool=false)`  
    ... determines the detailed radial grid due to the given gridType and parameters; a gr::Radial.Grid is returned
"""
function determineGrid(grid::Radial.Grid; printout::Bool=false) 
    
    # Read the meshType and general parameters from the given grid
    meshType = grid.meshType
    rnt    = grid.rnt;       h      = grid.h;         hp      = grid.hp;        NoPoints = grid.NoPoints
    nth    = grid.orderGL;   orderL = grid.orderL;    orderS  = grid.orderS;    orderGL  = grid.orderGL
    nr     = grid.NoPoints;  r      = zeros(nr);      rp      = zeros(nr);      rpor     = zeros(nr);      nsL = 0;   nsS = 0
        
    # Ensure that the largest grid points is always consistent with the largest break point of the B-splines
    if  NoPoints != NoPoints - rem(NoPoints,orderGL)    error("stop a")     end

    # Now define the physical grid due to the given parameters: either an exponential or exponential-linear grid
    if  hp == 0.
        r[1] = rnt;   rp[1] = rnt;   eph = exp(h);    ett = 1.0
        
        for  i  in 2:nr
            ett     = eph * ett
            ettm1   = ett - 1.0
            r[i]    = rnt * ettm1
            rp[i]   = rnt * ett
            rpor[i] = ett/ettm1
        end 

    elseif  hp != 0.
        function rprime(r :: Float64)
            return( 1/ (1/(r + rnt) + h/hp) )
        end
    
        function f(r :: Float64, i :: Int)
            return( log( r/rnt + 1) + h/hp * r - (i - 1) * h )
        end
    
        function fprime(r :: Float64)
            return( 1. / (r + rnt) + h/hp )
        end
    
        r[1]    = 0.
        rp[1]   = rprime(0.) 
        rpor[1] = 0.
    
        rc = 0.
        rn = 0.
    
        for i = 2:nr
            rn = rn + rnt
            while ((abs((rc - rn)/rn) > 100 * eps(Float64))) 
                rc = rn
                rn = rc - f(rc, i) / fprime(rc)
            end
            r[i]    = rn
            rp[i]   = rprime(rn)
            rpor[i] = rp[i] / r[i]
        end
    else
        error("stop a")
    end
    
    # Define the radial break points and the number of B-splines
    ntL = orderL;    ntS = orderS;    tL = zeros( orderL );    tS = zeros( orderS )
    for  i = nth:nth:NoPoints
        ntL = ntL + 1;    push!( tL, r[i]);    ntS = ntS + 1;    push!( tS, r[i])
    end
    for   i = ntL+1:ntL+orderL-1   push!( tL, tL[ntL] )   end;    nsL = ntL - 1;   ntL = length(tL)
    for   i = ntS+1:ntS+orderS-1   push!( tS, tS[ntS] )   end;    nsS = ntS - 1;   ntS = length(tS)
    
    # Define the radial grid due to the definition
    wr = zeros(nr)
    if  meshType == Radial.MeshGL()
        nr = 0;   r = Float64[];   rp = Float64[];   rpor = Float64[];   wr = Float64[];    rlow = 0.
        wax = gauss(orderGL);   ra = wax[1];   wa = wax[2] 
        for  i = 1:length(tL)
            rstep = tL[i] - rlow
            if  rstep > 0
                bma = tL[i] - rlow;   bpa = tL[i] + rlow 
                for  j = 1:orderGL
                    nr = nr + 1;    push!(r, ra[j] * bma/2. + bpa/2. );    push!(wr, wa[j] * bma/2. )
                end
            end
            rlow = tL[i]
        end
        #
        if  printout
            println("Define a radial grid of type $meshType with $nr grid points")
            println(" [rnt=" * @sprintf("%.3e",rnt) * ", h=" * @sprintf("%.3e",h) * 
                    ", hp=" * @sprintf("%.3e",hp) * ", NoPoints=$NoPoints, r_max=" * @sprintf("%.3e",r[nr]) * ";")
            println("  B-splines with break points at every $(nth)th point, nsL=$nsL, nsS=$nsS, orderL=$orderL, orderS=$orderS, orderGL=$orderGL,  " *
                    "ntL=$ntL, ntS=$ntS] ")
        end
    end

    Grid(rnt, h, hp, NoPoints,  tL, tS, ntL, ntS, orderL, orderS, nsL, nsS, orderGL,   meshType, r, rp, rpor, wr)
end


"""
`struct  Radial.GridGL`  ... defines a type for Gauss-Legendre grid of given order but where the interval is hard-coded due to
                                the given keystring

    + nt           :Int64                 ... number of mesh points in the grid.
    + t            ::Array{Float64,1}     ... mesh points in the variable t.
    + wt           ::Array{Float64,1}     ... weights to the mesh points in t.
"""
struct  GridGL
    nt             ::Int64  
    t              ::Array{Float64,1}
    wt             ::Array{Float64,1}
end



"""
`Radial.GridGL()`  
    ... specified a default version of a Gauss-Legendre grid with 6 points in the interval [0.,1.].
"""
function GridGL()
    GridGL("Finite", 0., 1., 6; printout=false)
end


"""
`Radial.GridGL("QED", orderGL::Int64; printout::Bool=false)`  
    ... constructor to define Gauss-Legendre grid for the typical QED computation in the interval [1.0, infinity].
"""
function GridGL(sa::String, orderGL::Int64; printout::Bool=false)
    !(sa == "QED")  && error("Unrecognized keystring; sa = $sa")
    txlow = 1.;    t = Float64[];    wt = Float64[];    nt = 0
    for i = 1:100000
        # Define the exponential increase (1.5) and the maximum size (infinity=150.)
        txup = txlow * 1.3;    if  txup > 100.   break   end
        #
        wax = QuadGK.gauss(orderGL);   tx = wax[1];   wtx = wax[2] 
        bma  = txup - txlow;   bpa = txup + txlow
        for  j = 1:orderGL
            nt = nt + 1;    push!(t, tx[j] * bma/2. + bpa/2. );    push!(wt, wtx[j] * bma/2. )
        end
        txlow = txup
    end
    if  printout   println("Gauss-Legendre grid with $nt mesh points from t = 1 ... $txlow.")   end
        
    GridGL(nt, t, wt)
end


"""
`Radial.GridGL("Finite", tmin::Float64, tmax::Float64, orderGL::Int64; printout::Bool=false)`  
    ... constructor to define Gauss-Legendre grid in the interval [tmin, tmax].
"""
function GridGL(sa::String, tmin::Float64, tmax::Float64, orderGL::Int64; printout::Bool=false)
    !(sa == "Finite")  && error("Unrecognized keystring; sa = $sa")
    t = Float64[];    wt = Float64[]
    
    wax = QuadGK.gauss(orderGL);    t = wax[1];     wt = wax[2]        
    fac1 = 0.5 * ( tmax - tmin );   fac2 = 0.5 * ( tmax + tmin )
    
    for  j = 1:orderGL
            t[j]    = fac1 * t[j] + fac2
            wt[j]   = fac1 * wt[j]
    end
    
    if  printout   println("Gauss-Legendre grid with $orderGL mesh points from t = $tmin ... $tmax.")   end
        
    GridGL(orderGL, t, wt)
end


# `Base.show(io::IO, grid::GridGL)`  ... prepares a proper printout of the variable grid::GridGL.
function Base.show(io::IO, grid::GridGL) 
    print(io, "Gauss-Legendre grid with $(grid.nt) mesh points: \n")
    print(io, "t:   ", grid.t[1:5],    "  ...  ", grid.t[grid.nt-4:grid.nt],     "\n") 
    print(io, "w:   ", grid.wt[1:5],   "  ...  ", grid.wt[grid.nt-4:grid.nt],    "\n") 
end


"""
`struct  Radial.GridGH`  ... defines a type for Gauss-Hermite grid of given order

    + nt           ::Int64                ... number of mesh points in the grid.
    + t            ::Array{Float64,1}     ... mesh points in the variable t.
    + wt           ::Array{Float64,1}     ... weights to the mesh points in t.
"""
struct  GridGH
    nt             ::Int64  
    t              ::Array{Float64,1}
    wt             ::Array{Float64,1}
end


"""
`Radial.GridGH(orderGH::Int64; printout::Bool=false)`  
    ... constructor to define Gauss-Hermite grid.
"""
function GridGH(orderGH::Int64; printout::Bool=false)
    t = Float64[];    wt = Float64[]
    
    wax = gausshermite(orderGH) #this function is included in the package FastGaussQuadrature
    
    t = wax[1]
    wt = wax[2]
    
    if  printout   println("Gauss-Hermite grid with $orderGH mesh points.")   end
        
    GridGH(orderGH, t, wt)
end


# `Base.show(io::IO, grid::GridGH)`  ... prepares a proper printout of the variable grid::GridGH.
function Base.show(io::IO, grid::GridGH) 
    print(io, "Gauss-Hermite grid with $(grid.nt) mesh points: \n")
    print(io, "t:   ", grid.t[1:5],    "  ...  ", grid.t[grid.nt-4:grid.nt],    "\n") 
    print(io, "w:   ", grid.wt[1:5],   "  ...  ", grid.wt[grid.nt-4:grid.nt],    "\n") 
end



"""
`struct  Radial.Density`  ... defines a struct for a radial density distribution.

    + name           ::String            ... A name for the radial density.
    + Dr             ::Array{Float64,1}  ... radial density function D(r).
    + grid           ::RadialGrid        ... radial grid on which the density is defined.
"""
struct  Density
    name             ::String
    Dr               ::Array{Float64,1}
    grid             ::Radial.Grid
end


"""
`Radial.Density()`  ... constructor to define an 'empty' instance of the radial density.
"""
function Density()
    Density("", Float64[], Radial.Grid())
end


# `Base.show(io::IO, density::Radial.Density)`  ... prepares a proper printout of the variable density::Radial.Density.
function Base.show(io::IO, density::Radial.Density) 

    sa = Base.string(density);    print(io, sa * "\n")

    n = length(density.Dr);                    if  n < 6   return( nothing )    end
    print(io, "Dr:    ", density.Dr[1:23],    "  ...  ", density.Zr[n-22:n],    "\n") 
    print(io, density.grid)
end


# `Base.string(density::Radial.Density)`  ... provides a String notation for the variable density::Radial.Density.
function Base.string(density::Radial.Density) 
    if  length(density.Dr)  == 0
        sa = "Radial density not yet defined; kind = $(density.name) ..."
    else
        sa = "$(density.name) (radial) density ... defined on $(length(density.Zr)) grid points ..."
    end 
end



"""
`struct  Radial.Potential`  ... defines a struct for a local radial potential.

    + name           ::String            ... A name for the potential.
    + Zr             ::Array{Float64,1}  ... radial potential function Z(r) = - r * V(r) as usual in atomic structure theory.
    + grid           ::RadialGrid        ... radial grid on which the potential is defined.
"""
struct  Potential
    name             ::String
    Zr               ::Array{Float64,1}
    grid             ::Radial.Grid
end


"""
`Radial.Potential()`  ... constructor to define an 'empty' instance of the radial potential.
"""
function Potential()
    Potential("", Float64[], Radial.Grid())
end


# `Base.show(io::IO, potential::Radial.Potential)`  ... prepares a proper printout of the variable potential::Radial.Potential.
function Base.show(io::IO, potential::Radial.Potential) 

    sa = Base.string(potential);    print(io, sa * "\n")

    n = length(potential.Zr);                    if  n < 6   return( nothing )    end
    print(io, "Zr:    ", potential.Zr[1:23],    "  ...  ", potential.Zr[n-22:n],    "\n") 
    print(io, potential.grid)
end


# `Base.string(potential::Radial.Potential)`  ... provides a String notation for the variable potential::Radial.Potential.
function Base.string(potential::Radial.Potential) 
    if  length(potential.Zr)  == 0
        sa = "Radial potential not yet defined; kind = $(potential.name) ..."
    else
        sa = "$(potential.name) (radial) potential ... defined on $(length(potential.Zr)) grid points ..."
    end 
end


"""
`struct  Radial.Orbital`  
    ... defines a type for a single-electron radial orbital function with a large and small component, and which can refer to
        either the standard or an explicitly given grid due to the logical flag useStandardGrid. Bound-state orbitals with energy < 0 are 
        distinguished from free-electron orbitals by the flag isBound.
        

    + subshell        ::Subshell          ... Relativistic subshell.
    + isBound         ::Bool              ... Logical flag to distinguish between bound (true) and free-electron orbitals (false).
    + useStandardGrid ::Bool              ... Logical flag for using the standard grid (true) or an explicitly given grid (false).
    + energy          ::Float64           ... Single-electron energies of bound orbitals are always negative.
    + P               ::Array{Float64,1}  ... Large and ..
    + Q               ::Array{Float64,1}  ... small component of the radial orbital.
    + Pprime          ::Array{Float64,1}  ... dP/dr.
    + Qprime          ::Array{Float64,1}  ... dQ/dr.
    + grid            ::Array{Float64,1}  ... explic. defined radial grid array for P, Q, if StandardGrid = false.
"""
mutable struct Orbital
    subshell          ::Subshell
    isBound           ::Bool             
    useStandardGrid   ::Bool
    energy            ::Float64 
    P                 ::Array{Float64,1} 
    Q                 ::Array{Float64,1}
    Pprime            ::Array{Float64,1}
    Qprime            ::Array{Float64,1}
    grid              ::Radial.Grid
end


"""
`Radial.Orbital(subshell::Subshell, energy::Float64)`  
    ... constructor for given subshell and energy, and where useStandardGrid is set to true; the grid must be defined 
        explicitly and neither the large and small components nor their derivatives are yet defined in this case.
"""
function Orbital(subshell::Subshell, energy::Float64)
    if energy < 0    isBound = true    else    isBound = false    end
    useStandardGrid = true
    P = Array{Float64,1}[];    Q = Array{Float64,1}[];    grid = Radial.Grid()
    Pprime = Array{Float64,1}[];    Qprime = Array{Float64,1}[]

    Orbital(subshell, isBound, useStandardGrid, energy, P, Q, Pprime, Qprime, grid)
end


"""
`Radial.Orbital(label::String, energy::Float64)`  
    ... constructor for given string identifier and energy, and where useStandardGrid is set to true; the grid must be 
        defined explicitly and neither the large and small components nor their derivatives are yet defined in this case.
"""
function Orbital(label::String, energy::Float64)
    if energy < 0    isBound = true    else    isBound = false    end

    subshell = Subshell(label);    useStandardGrid = true
    P = Array{Float64,1}[];        Q = Array{Float64,1}[];    grid = Radial.Grid()
    Pprime = Array{Float64,1}[];    Qprime = Array{Float64,1}[]

    Orbital(subshell, isBound, useStandardGrid, energy, P, Q, Pprime, Qprime, grid)
end


# `Base.show(io::IO, orbital::Orbital)`  ... prepares a proper printout of the variable orbital::Orbital.
function Base.show(io::IO, orbital::Orbital) 
    n = length(orbital.P)

    if   orbital.useStandardGrid
        stdgrid = Defaults.getDefaults("standard grid")
        stdgrid.NoPoints == 0                     &&   return( print("Standard grid has not yet been defined.") )

        n = min(length(orbital.P), stdgrid.NoPoints)
        n > stdgrid.NoPoints                      &&   error("length of P does not match to standard grid; n=$n  NoPoints=$(stdgrid.NoPoints) ")    
        length(orbital.P) != length(orbital.Q)    &&   error("P and Q have different length")  
    else  
        !(n == length(orbital.P) == length(orbital.Q))   &&    error("P, Q, grid have different length")
    end

    sa = Base.string(orbital::Orbital);    print(io, sa * "\n")

    if n <= 6
        print(io, "Large component P: ", orbital.P[1:end], "\n") 
        print(io, "Small component Q: ", orbital.Q[1:end], "\n") 
        print(io, "Pprime:            ", orbital.Pprime[1:end], "\n") 
        print(io, "Qprime:            ", orbital.Qprime[1:end], "\n") 
    else 
        n = length(orbital.P)
        print(io, "Large component P: ", orbital.P[1:25], "  ...  ", orbital.P[n-10:n], "\n") 
        print(io, "Small component Q: ", orbital.Q[1:25], "  ...  ", orbital.Q[n-10:n], "\n")
        print(io, "Pprime:            ", orbital.Pprime[1:25], "  ...  ", orbital.Pprime[n-10:n], "\n") 
        print(io, "Qprime:            ", orbital.Qprime[1:25], "  ...  ", orbital.Qprime[n-10:n], "\n")
    end

    if orbital.useStandardGrid 
        n <= 6   &&   print(io, "Defined on Grid:   ", stdgrid.r[1:n], "\n") 
        n >  6   &&   print(io, "Defined on Grid:   ", stdgrid.r[1:3], "  ...  ", stdgrid.r[n-3:n], "\n") 
    else
        n <= 6   &&   print(io, "Defined on Grid:   ", orbital.grid.r[1:n], "\n") 
        n > 6    &&   print(io, "Defined on Grid:   ", orbital.grid.r[1:5], "  ...  ", orbital.grid.r[n-4:n], "\n") 
    end
end


# `Base.string(orbital::Orbital)`  ... provides a String notation for the variable orbital::Orbital with just basic information.
function Base.string(orbital::Orbital) 
    if  orbital.isBound           sa = "Bound-state orbital "     else    sa = "Free-electron orbital "         end
    if  orbital.useStandardGrid   sb = "the standard grid: "      else    sb = "an explicitly-defined grid: "   end
    energy = orbital.energy
    sc = string(orbital.subshell) * " with energy $energy a.u. ";   n = length(orbital.P)

    return( sa * sc * "is defined with $n (grid) points on " * sb )
end


"""
`struct  Radial.SingleSymOrbitals`  
    ... defines a type for a (relativistic and quasi-complete) single-electron spectrum for symmetry kappa with N positive 
        and/or N negative states. All these states are defined with regard to the same grid.

    + kappa        ::Int64                 ... symmetry of the one-electron spectrum
    + NoStates     ::Int64                 ... Number of positive and negative states (if onlyPositive = true)
    + onlyPositive ::Bool                  ... True if only the positive part is kept.
    + pOrbitals    ::Array{Orbital,1}      ... Positive-energy orbitals states, in increasing order.
    + nOrbitals    ::Array{Orbital,1}      ... Negative-energy orbitals states, in increasing order.
    + grid         ::RadialGrid            ... radial grid on which the states are represented.
"""
struct SingleSymOrbitals
    kappa          ::Int64
    NoStates       ::Int64
    onlyPositive   ::Bool
    pOrbitals      ::Array{Orbital,1} 
    nOrbitals      ::Array{Orbital,1} 
    grid           ::Radial.Grid
end


"""
    `Radial.SingleSymOrbitals()`  ... constructor for providing an 'empty' instance of this struct.
"""
function SingleSymOrbitals()
    SingleSymOrbitals(0, 0, true, Orbital[], Orbital[], Radial.Grid() )
end


# `Base.show(io::IO, symOrbitals::SingleSymOrbitals)`  ... prepares a proper printout of the variable symOrbitals::SingleSymOrbitals.
function Base.show(io::IO, symOrbitals::SingleSymOrbitals) 
    sa = "Single-symmetry orbital set for kappa=$(symOrbitals.kappa) with $(symOrbitals.NoStates) states "
    if  symOrbitals.onlyPositive   sa = sa * "from just the positive part of the spectrum."
    else                           sa = sa * "from the positive and negative part of the spectrum."
    end
    println(io, sa) 
end


"""
`struct  Radial.SingleElecSpectrum`  
    ... defines a type for a (relativistic and quasi-complete) single-electron spectrum for different symmetries kappa and either
        (individually) N_kappa positive  or  N_kappa positive and negative states. All these orbitals are defined with regard to 
        the same grid.

    + name         ::String                      ... A name for this spectrum, may contain information about the original 
                                                        type of basis functions.
    + symOrbitals  ::Array{SingleSymOrbitals,1}  ... Set of orbitals of the same kappa-symmetry.
"""
struct  SingleElecSpectrum
    name           ::String
    symOrbitals    ::Array{SingleSymOrbitals,1}
end


"""
    `Radial.SingleElecSpectrum()`  ... constructor for providing an 'empty' instance of this struct.
"""
function SingleElecSpectrum()
    SingleElecSpectrum("", SingleSymOrbitals[])
end


# `Base.show(io::IO, spectrum::SingleElecSpectrum)`  ... prepares a proper printout of the variable spectrum::SingleElecSpectrum.
function Base.show(io::IO, spectrum::SingleElecSpectrum) 
    sa = "Single-electron spectrum $(spectrum.name) for different kappa-symmetries and (individual) number of orbitals " *
            "positive and/or positive & negative energy."
    println(io, sa) 
end


"""
`Radial.OrbitalBunge1993(subshell::Subshell, Z::Int64)`  ... to calculate the radial orbital on the standard grid.
"""
function OrbitalBunge1993(subshell::Subshell, Z::Int64, grid::Radial.Grid)
    #x grid = Defaults.getDefaults("standard grid");    
    #x qn = SubshellQuantumNumbers( string(subshell) );    n = qn[1];    kappa = qn[2];   l = qn[3]
    n = subshell.n;    kappa = subshell.kappa;   l = Basics.subshell_l(subshell)

    wa = Basics.store("orbital functions: NR, Bunge (1993)", Z)
    if      l == 0  &&  n > wa[5][1]+0     error("orbitals with s symmetry for Z = $Z only available up to $(wa[5][1]+0) only.")
    elseif  l == 1  &&  n > wa[6][1]+1     error("orbitals with p symmetry for Z = $Z only available up to $(wa[6][1]+1) only.")
    elseif  l == 2  &&  n > wa[7][1]+2     error("orbitals with d symmetry for Z = $Z only available up to $(wa[7][1]+2) only.")
    elseif  l == 3  &&  n > wa[8][1]+3     error("orbitals with f symmetry for Z = $Z only available up to $(wa[8][1]+3) only.")
    elseif  l == 4  &&  n > wa[9][1]+4     error("orbitals with g symmetry for Z = $Z only available up to $(wa[9][1]+4) only.")
    end

    P = Float64[]
    for  k in 1:grid.NoPoints-5
        Rnl = 0.; r = grid.r[k]
        for i in 4:length(wa[5+l])
            njl = wa[5+l][i][1] 
            Zjl = wa[5+l][i][2]
            Cjl = wa[5+l][i][2+n-l]
            Njl = (2*Zjl)^(njl+0.5) / sqrt( factorial(2*njl) ) 
            Rnl = Rnl + Cjl*Njl*r^(njl-1)*exp(-Zjl*r);
        end
        push!(P, Rnl*r )
    end

    # Now define small component by apllying kinetic-balance condition
    Q = zeros(size(P, 1))
    #x println("small component not yet defined properly; P = $P ")

    dP(i) = JenaAtomicCalculator.Math.derivative(P, i)
    for i = 2:size(Q, 1)
        Q[i] = -1/(2 * JenaAtomicCalculator.Defaults.INVERSE_FINE_STRUCTURE_CONSTANT) * (dP(i) / grid.h / grid.rp[i] + kappa/grid.r[i]) * P[i]
    end

    o     = Orbital( subshell, true, true, -0.5, P, Q, grid)
    norma = RadialIntegrals.overlap(o, o, grid)
    o.P = o.P/sqrt(norma)
    o.Q = o.Q/sqrt(norma)

    normb = RadialIntegrals.overlap(o, o, grid)

    return( o )
end


"""
`Radial.OrbitalMcLean1981(subshell::Subshell,Z::Int64)`  ... to calculate the radial orbital on the standard grid.  
"""
function compute_McLean1981(subshell::Subshell,Z::Int64)
    grid = Defaults.getDefaults("standard grid");    
    #x qn = Basics.SubshellQuantumNumbers( string(subshell) );    n = qn[1];   l = qn[3]
    n = subshell.n;   l = Basics.subshell_l( subshell )

    wa = Basics.store("orbital functions: NR, McLean (1981)", Z)

    !(55 <= Z <= 92)   &&                     error("Nuclear charge must be 55 <= (Z=$Z) <= 92.")            
    wa == Any[]        &&                     error("No data available for neutral Z = $Z.")

    if      l == 0     &&  n > wa[5][1]+0     error("Orbitals with s symmetry for Z = $Z only available up to $(wa[5][1]+0) only.")
    elseif  l == 1     &&  n > wa[6][1]+1     error("orbitals with p symmetry for Z = $Z only available up to $(wa[6][1]+1) only.")
    elseif  l == 2     &&  n > wa[7][1]+2     error("orbitals with d symmetry for Z = $Z only available up to $(wa[7][1]+2) only.")
    elseif  l == 3     &&  n > wa[8][1]+3     error("orbitals with f symmetry for Z = $Z only available up to $(wa[8][1]+3) only.")
    elseif  l == 4     &&  n > wa[9][1]+4     error("orbitals with g symmetry for Z = $Z only available up to $(wa[9][1]+4) only.")
    end

    P = Float64[]
    for  k in 1:grid.NoPoints-5
        Rnl = 0.; r = grid.r[k]
        for i in 4:length(wa[5+l])
            njl = wa[5+l][i][1] 
            Zjl = wa[5+l][i][2]
            Cjl = wa[5+l][i][2+n-l]
            Njl = (2*Zjl)^(njl+0.5) / sqrt( factorial(2*njl) ) 
            Rnl = Rnl + Cjl*Njl*r^(njl-1)*exp(-Zjl*r);
        end
        println("k, r, Rnl = $k  $r  $Rnl ")
        push!(P, Rnl*r )
    end

    println("small component not yet defined properly; P = $P ")

    return( Orbital( subshell, -0.5) )
end


"""
`Radial.OrbitalPrimitiveSlater(subshell::Subshell, N::Int64, alpha0::Float64, beta0::Float64, grid::Radial.Grid)`  
    ... to calculate a list of (non-relativistic) radial Slater primitives::Array{Orbital,1} for the subshell-symmetry on 
        the given grid. All orbitals have the same subshell, isBound = true, useStandardGrid = true and energy = 0. 
        The small components are constructed by applying kinetic balance to the large component, and all the primitives are 
        properly normalized.
"""
function OrbitalPrimitiveSlater(subshell::Subshell, N::Int64, alpha0::Float64, beta0::Float64, grid::Radial.Grid) 
    function P(i::Int64, l::Int64, eta::Float64)   return( grid.r[i]^(l+1) * exp(-eta*grid.r[i]) )   end
    dP(i) = JenaAtomicCalculator.Math.derivative(Px, i)

    orbitals = Orbital[];    ll = Basics.subshell_l(subshell)
    for  i = 1:N
        eta = alpha0 * beta0^(i-1)
        # Px  = [for i = 1:grid.NoPoints P(i, ll, eta)]
        Qx  = zeros(size(Px, 1))
        for j = 2:size(Q, 1)
            Qx[j] = -1/(2 * JenaAtomicCalculator.Defaults.INVERSE_FINE_STRUCTURE_CONSTANT) * (dP(j) / grid.h / grid.rp[j] + kappa/grid.r[j]) * Px[j]
        end

        o     = Orbital( subshell, true, true, 0., Px, Qx, grid() )
        norma = JenaAtomicCalculator.RadialIntegrals.overlap(o, o, grid)
        o.P = o.P/sqrt(norma)
        o.Q = o.Q/sqrt(norma)

        normb = JenaAtomicCalculator.RadialIntegrals.overlap(o, o, grid)
        println("OrbitalPrimitiveSlater-aa:  for i = $i : norm-before = $norma, norm-after = $normb")
        
        push!( orbitals, o)
    end

    return( SingleElecSpectrum() )
end


###################################################################################################################
###################################################################################################################
###################################################################################################################


"""
`Radial.checkGrid(grid::Radial.Grid, shell::Shell, Z::Float64)`  
    ... checks that the radial grid can represent the hydrogenic orbitals for the given subshell 
        of nuclear charge Z; it simply tests that the boxSize > 4 * n^2 / Z [a_o].
        (More accurate one could use the mean expectation values of r of this shell).
        The program reports about this comparison and terminates if this is not the case. Nothing is returned.
"""
function checkGrid(grid::Radial.Grid, shell::Shell, Z::Float64)
    wa = 4.0 * shell.n^2 / Z
    println(">> Check radial grid for shell $shell of hydrogenic ions with Z = $Z;" *
            "  r[NoPoints] = $(grid.r[grid.NoPoints]) ! > ! 4 * n^2 / Z = $wa")
        
    if  grid.r[grid.NoPoints] < wa    error("Grid with inappropriate extent; enlarge the boxSize.")   end
    
    return( nothing )
end


"""
`Radial.checkGrid(grid::Radial.Grid, subshell::Subshell, relAccuracy::Float64)`  
    ... checks that the radial grid can represent the hydrogenic orbitals for the given subshell with
        a relative accuracy <= relAccuracy. The program reports about the relative accuracy and terminates with an 
        error message if this is not the case. Nothing is returned.
"""
function checkGrid(grid::Radial.Grid, subshell::Subshell, relAccuracy::Float64)
    # Generate a hydrogenic spectrum for the symmetry of the given subshell
    error("Not yet implemented.")
    
    return( nothing )
end


"""
`Radial.determineZbar(pot::Radial.Potential)`  
    ... determines the effective charge that is asymptotically seen by the electron in the potential pot. 
        A Zbar::Float64 is returned.
"""
function determineZbar(pot::Radial.Potential) 
    nx = 5   # Number of grid points for determining the effective charge Zbar
    mtp = size( pot.Zr, 1);    grid = pot.grid;    meanZbar = zeros(nx);    devsZbar = zeros(nx);   ny = 0
        for  i = mtp-nx+1:mtp
        ny  = ny + 1;    meanZbar[ny] = pot.Zr[mtp]
    end
    mZbar   = sum(meanZbar) / nx;   for  i = 1:nx    devsZbar[i] = (meanZbar[i] - mZbar)^2    end
    stdZbar = sqrt( sum(devsZbar) / nx )
    if  true  println(">> Radial potential with effective charge Zbar=" * @sprintf("%.4e",mZbar) *
                        " (Delta-Zbar=" * @sprintf("%.4e",stdZbar) * ") at r=" * @sprintf("%.4e",pot.grid.r[mtp]) * " a.u." )    end
    
    return( mZbar )
end


"""
`Radial.determineNoPoints(rnt::Float64, h::Float64, hp::Float64, rbox::Float64, orderGL::Int64)`  
    ... determines the number of points of the physical size if the grid parameters and the box-size is given; moreover,
        it is ensured that this NoPoints is consistent with the largest break point for the B-spline grid.
        An NoPoints::Int64 is returned.
"""
function determineNoPoints(rnt::Float64, h::Float64, hp::Float64, rbox::Float64, orderGL::Int64) 

    if  hp == 0.
        NoPoints = 1;  eph = exp(h);    ett = 1.0
        while true
            NoPoints = NoPoints + 1
            ett      = eph * ett
            ettm1    = ett - 1.0
            rmax     = rnt * ettm1
            if  rmax > rbox  &&  rem(NoPoints, orderGL) == 0    break   end
        end 

    elseif  hp != 0.
        NoPoints = 1;  rc = 0.;    rn = 0.

        function f(r :: Float64, i :: Int)
            return( log( r/rnt + 1) + h/hp * r - (i - 1) * h )
        end
    
        function fprime(r :: Float64)
            return( 1. / (r + rnt) + h/hp )
        end
        while true
            NoPoints = NoPoints + 1
            rn = rn + rnt
            while ((abs((rc - rn)/rn) > 100 * eps(Float64))) 
                rc = rn
                rn = rc - f(rc, NoPoints) / fprime(rc)
            end
            rmax    = rn
            if  rmax > rbox  &&  rem(NoPoints, orderGL) == 0    break   end
        end
    end
    
    return( NoPoints )
end


"""
`Radial.generateGrid(grid::Radial.Grid; boxSize::Union{Nothing,Float64}=nothing, 
                                        maximumFreeElectronEnergy::Union{Nothing,Float64}=nothing, 
                                        maximumPrincipalQN::Union{Nothing,Int64}=nothing, 
                                        NoPointsInsideNucleus::Union{Nothing,Int64}=nothing, 
                                        NoPointsInsideFirstBohrRadius::Union{Nothing,Int64}=nothing)
                                        
    ... to generate a grid that fulfills special requirements; th following schemes are supported:
    
        boxSize 
            ... apply a fixed box size (in atomic units) as needed in an average-atom model and elsewhere.
        boxSizeWithZeroWeights 
            ... apply a fixed box size (in atomic units) as needed in an average-atom model and elsewhere.
                Here, a standard grid is used but with all weights w [r] = 0 for r > boxSize.
        maximumFreeElectronEnergy
            ... provide a maximum free-electron energy (Hartree) and determine the linearized stepsize such,
                that 20 points per wavelength are used asymptotically.
        maximumPrincipalQN
            ... generate an grid with a rbox-size that is suitable to represent subshell orbitals
                with the given n; the value of rbox = 5 * <r_n> is taken, i.e. 5 times the mean hydrogenic 
                value (not yet).
        NoPointsInsideNucleus
            ... generate a grid (of given type) with the given number inside the nucleus; this affects the values
                of rnt and h (not yet).
        NoPointsInsideFirstBohrRadius
            ... generate a grid (of given type) with the given number inside the first Bohr radius; this affects 
                the values of rnt and h (not yet).
    
    Only on of these optional parameters can be selected at a given time. A proper grid::Radial.Grid is returned,
    along with a short reasoning of what has been selected.
"""
function generateGrid(grid::Radial.Grid; boxSize::Union{Nothing,Float64}=nothing, 
                                            boxSizeWithZeroWeights::Union{Nothing,Float64}=nothing,
                                            maximumFreeElectronEnergy::Union{Nothing,Float64}=nothing, 
                                            maximumPrincipalQN::Union{Nothing,Int64}=nothing, 
                                            NoPointsInsideNucleus::Union{Nothing,Int64}=nothing, 
                                            NoPointsInsideFirstBohrRadius::Union{Nothing,Int64}=nothing) 
    
    if  typeof(boxSize)                            != Nothing
        # Apply a fixed box size (in atomic units) as needed in an average-atom model and elsewhere.
        newGrid = Radial.Grid(grid; rbox=boxSize);    rnt = newGrid.rnt * boxSize / newGrid.tL[end]
        newGrid = Radial.Grid(newGrid; rnt = rnt)
        println(">> Generate a new grid with boxSize = $boxSize a.u. and the break points " *
                "\n   tL = $(newGrid.tL)   \n   tS = $(newGrid.tS)"  *
                "\n>> Note that the last grid point r[end] is always slightly smaller as it refers " *
                "to the last GL (zero of) integration along r.")
        #
    elseif  typeof(boxSizeWithZeroWeights)         != Nothing
        # Apply a fixed box size (in atomic units) as needed in an average-atom model and elsewhere.
        # All (integration) weight are set to zero for r > boxSizeWithZeroWeights.
        # First determine true boxsize
        boxSize = 0.;   for  t in grid.tL  if  t < boxSizeWithZeroWeights   boxSize = t   end   end
        println(">>  Generate a new grid with zero weights for boxSize = $boxSize a.u.; " * 
                "\n    modify the (grid) parameters of the given grid to bring boxSize closer to boxSizeWithZeroWeights.")
        wr = Float64[];   R = 0.;   iR = 0    
        for  (ir, r)  in  enumerate(grid.r)   
            if  r > boxSize      push!(wr, 0.)
            else                 push!(wr, grid.wr[ir]);   R = r;   iR = ir
            end 
        end
        ## println(">>> wr = $wr   ... just to see the number of zero weights for the new grid. \n")
        newGrid = Radial.Grid(grid.rnt, grid.h, grid.hp, grid.NoPoints, grid.tL, grid.tS, grid.ntL, grid.ntS, 
                                grid.orderL, grid.orderS, grid.nsL, grid.nsS, grid.orderGL, grid.meshType,
                                grid.r, grid.rp, grid.rpor, wr)
        println(">> Generate a new grid with boxSize = $boxSize a.u. and with all weight grid.wr[j] = 0 " * 
                "for r_j > boxSizeWithZeroWeights (j > $iR)" *
                "\n>> Note that the last non-zero weight at R = $R grid point r[end] is always slightly smaller than " *
                "\n   boxSize = $boxSize as it refers to the last GL (zero of) integration point.")
        #
    elseif  typeof(maximumFreeElectronEnergy)      != Nothing
        wavenb      = sqrt( 2maximumFreeElectronEnergy + maximumFreeElectronEnergy * Defaults.getDefaults("alpha")^2 )
        wavelgth    = 2pi / wavenb;     hp = wavelgth / 20
        newGrid     = Radial.Grid(grid; hp=hp)
        println(">> Generate a new grid for the minium wavelength = $wavelgth a.u. of free electrons and  hp = $hp")
        #
    elseif  typeof(maximumPrincipalQN)             != Nothing
        error("stop a")
        #
    elseif  typeof(NoPointsInsideNucleus)          != Nothing
        error("stop b")
    elseif  typeof(NoPointsInsideFirstBohrRadius)  != Nothing
        error("stop c")
    end
        
    return( newGrid )
end


end # module
