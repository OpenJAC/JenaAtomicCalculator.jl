
"""
`module  JAC.Radial`  
... a submodel of JAC that contains all structs and methods to deal with the radial grid, radial orbitals and potentials.
"""
module Radial


using  QuadGK, Printf, FastGaussQuadrature, JenaAtomicCalculator, ..Basics

export Grid, Potential, Orbital, SingleElecSpectrum


## Radial.AbstractMesh, with MeshNone, MeshGrasp and MeshGL, was REMOVED on 12-Aug-2026.
##
## Nothing in JAC could ever produce anything but a Gauss-Legendre mesh: meshType was assigned in exactly one
## place, Radial.Grid(exponential::Bool), always to MeshGL(), and no file in src/, examples/ or test/ ever
## passed the meshType= keyword.  All 26 MeshGrasp comparisons were therefore permanently false, and what
## hung off them was substantial: 25 dead branches in module-RadialIntegrals.jl, the grid fields rp and rpor
## (empty arrays on every grid that has ever existed), Math.integrateFit/integrateFitTransform, and the
## Basics.integrate(::NewtonCotes/::SimpsonRule/::TrapezRule) family with its four exported type names.
##
## That region had been unreachable long enough for a one-character typo to sit undetected in its central
## routine -- Math.derivative read `JenaAtomicCalculatorDefaults.weights`, so its interior-point branch, the
## common case, raised UndefVarError.  All of it is now gone, and Radial.Grid describes one thing: a
## Gauss-Legendre mesh laid over B-spline break points, which is what every JAC calculation has actually used.


"""
`abstract type Radial.AbstractGridGaussLegendreScheme`
    ... labels the kind of Gauss-Legendre grid to be constructed; it is used for dispatch and to avoid
        string comparisons.
    Concrete subtypes:
    + GridGaussLegendreQED    ... generate a GL grid for QED computations over the interval [1, infinity).
    + GridGaussLegendreFinite ... generate a GL grid over a finite interval [tmin, tmax].
"""
abstract type  AbstractGridGaussLegendreScheme                                              end
struct         GridGaussLegendreQED     <:  AbstractGridGaussLegendreScheme                 end
struct         GridGaussLegendreFinite  <:  AbstractGridGaussLegendreScheme                 end



"""
`struct  Radial.GridParameters`
    ... the ANALYTIC DEFINITION of the radial grid: the handful of numbers a user chooses, from which
        everything else follows.

    + rnt        ::Float64  ... smallest grid point > 0.
    + h          ::Float64  ... stepsize in the construction of the exponential grid.
    + hp         ::Float64  ... asymptotic stepsize of the log-lin grid; hp = 0 gives a purely exponential one.
    + NoPoints   ::Int64    ... No. of points of the physical grid, chosen so that its last point coincides
                                with the largest B-spline break point.  NOTE that this counts the PHYSICAL
                                grid of the definition, which is a different quantity from the number of
                                mesh points length(mesh.r), even where the two happen to coincide.
    + orderGL    ::Int64    ... order of the Gauss-Legendre integration; it also fixes the break points, by
                                taking every orderGL-th point of the physical grid.
"""
struct  GridParameters
    rnt          ::Float64
    h            ::Float64
    hp           ::Float64
    NoPoints     ::Int64
    orderGL      ::Int64
end


"""
`struct  Radial.KnotSequence`
    ... the B-SPLINE BREAK POINTS, separately for the large and the small component.

    + tL         ::Array{Float64,1}  ... radial break points for the B-splines of the large component.
    + tS         ::Array{Float64,1}  ... radial break points for the B-splines of the small component.
    + orderL     ::Int64             ... B-spline order of the large components.
    + orderS     ::Int64             ... B-spline order of the small components.

        The counts ntL, ntS, nsL and nsS are NOT stored: they are functions of the above and are provided as
        properties of Radial.Grid, so that they cannot fall out of step with the knots they describe.
"""
struct  KnotSequence
    tL           ::Array{Float64,1}
    tS           ::Array{Float64,1}
    orderL       ::Int64
    orderS       ::Int64
end


"""
`struct  Radial.RadialMesh`
    ... the MESH on which radial functions are tabulated and integrated: a composite Gauss-Legendre rule
        laid over the intervals between the break points of the knot sequence.

    + r          ::Array{Float64,1}  ... mesh points, the Gauss-Legendre nodes of each interval of tL.
    + wr         ::Array{Float64,1}  ... the Gauss-Legendre weights belonging to r.
"""
struct  RadialMesh
    r            ::Array{Float64,1}
    wr           ::Array{Float64,1}
end


"""
`struct  Radial.Grid`
    ... the radial grid: its analytic definition, the B-spline break points that follow from it, and the
        mesh on which radial integrations are performed.

    + parameters ::Radial.GridParameters  ... what was asked for.
    + knots      ::Radial.KnotSequence    ... the B-spline break points.
    + mesh       ::Radial.RadialMesh      ... the Gauss-Legendre mesh points and weights.

        Until 12-Aug-2026 these seventeen numbers and arrays sat side by side in one flat struct, so that
        nothing distinguished the definition from what was derived from it, and four of the fields (ntL, ntS,
        nsL, nsS) duplicated information already present in the knots.

        FOR COMPATIBILITY, AND BECAUSE THEY READ WELL, all the former field names remain available directly
        on the Grid: grid.r, grid.wr, grid.rnt, grid.tL, grid.nsL and the rest all work exactly as before
        (see Base.getproperty below).  ntL, ntS, nsL and nsS are now COMPUTED from the knot sequence rather
        than stored, which is the one behavioural difference: they can no longer be stale.
"""
struct  Grid
    parameters   ::GridParameters
    knots        ::KnotSequence
    mesh         ::RadialMesh
end


"""
`Base.getproperty(grid::Radial.Grid, s::Symbol)`
    ... forwards the field names of the former flat Radial.Grid to the three structs that now hold them, and
        computes ntL, ntS, nsL and nsS from the knot sequence.  About 700 accesses across JAC read these
        names; forwarding them is deliberate rather than transitional, since `grid.r` reads better than
        `grid.mesh.r` at the point of use.  A literal `grid.r` is constant-folded and costs nothing.
"""
function Base.getproperty(grid::Radial.Grid, s::Symbol)
    ## the analytic definition
    s === :rnt       &&  return( getfield(grid, :parameters).rnt      )
    s === :h         &&  return( getfield(grid, :parameters).h        )
    s === :hp        &&  return( getfield(grid, :parameters).hp       )
    s === :NoPoints  &&  return( getfield(grid, :parameters).NoPoints )
    s === :orderGL   &&  return( getfield(grid, :parameters).orderGL  )
    ## the knot sequence
    s === :tL        &&  return( getfield(grid, :knots).tL      )
    s === :tS        &&  return( getfield(grid, :knots).tS      )
    s === :orderL    &&  return( getfield(grid, :knots).orderL  )
    s === :orderS    &&  return( getfield(grid, :knots).orderS  )
    ## ... and the four counts that used to be stored alongside it
    s === :ntL       &&  return( length(getfield(grid, :knots).tL) )
    s === :ntS       &&  return( length(getfield(grid, :knots).tS) )
    s === :nsL       &&  return( length(getfield(grid, :knots).tL) - getfield(grid, :knots).orderL )
    s === :nsS       &&  return( length(getfield(grid, :knots).tS) - getfield(grid, :knots).orderS )
    ## the mesh
    s === :r         &&  return( getfield(grid, :mesh).r  )
    s === :wr        &&  return( getfield(grid, :mesh).wr )
    return( getfield(grid, s) )
end


# `Base.propertynames(grid::Radial.Grid)`  ... so that the forwarded names are visible to the REPL and to tab completion.
Base.propertynames(grid::Radial.Grid) =
    (:parameters, :knots, :mesh, :rnt, :h, :hp, :NoPoints, :orderGL, :tL, :tS, :orderL, :orderS,
     :ntL, :ntS, :nsL, :nsS, :r, :wr)


"""
`Radial.Grid()` ... constructor to define an 'empty' grid.
"""
function Grid()
    Radial.Grid( GridParameters(0., 0., 0., 0, 0), KnotSequence(Float64[], Float64[], 0, 0),
                 RadialMesh(Float64[], Float64[]) )
end


"""
`Radial.Grid(exponential::Bool; printout::Bool=false)` 
    ... constructor to define either a standard 'exponential' mesh (true) or a 'log-lin' mesh (false).
"""
function Grid(exponential::Bool; printout::Bool=false)
    orderL = 7;    orderS = 8;   orderGL = 7
    
    if  exponential
        NoPoints = 392;    NoPoints = NoPoints - rem(NoPoints,orderGL)
        grid = Radial.Grid( GridParameters(2.0e-6, 5.0e-2,     0., NoPoints, orderGL),
                            KnotSequence(Float64[], Float64[], orderL, orderS), RadialMesh(Float64[], Float64[]) )
    else
        NoPoints = 600;    NoPoints = NoPoints - rem(NoPoints,orderGL)
        grid = Radial.Grid( GridParameters(2.0e-6, 5.0e-2, 2.0e-2, NoPoints, orderGL),
                            KnotSequence(Float64[], Float64[], orderL, orderS), RadialMesh(Float64[], Float64[]) )
    end
    return( Radial.determineGrid(grid, printout=printout) )
end


"""
`Radial.Grid(gr::Radial.Grid;`

        rnt=..,     h=..,       hp=..,      rbox=..,    orderL=..,  orderS=..,  orderGL=..,  printout=..)
    ... constructor for modifying the given Radial.Grid by 'overwriting' the previously selected parameters.

        The keyword nth= was removed on 12-Aug-2026.  It was accepted and then silently dropped: the spacing
        of the B-spline break points is set by orderGL in Radial.determineGrid, and Grid(gr; nth=3) left the
        grid unchanged.  Nothing in src/, examples/ or test/ ever passed it.  A second knob for the same
        quantity is not restored here; use orderGL.
"""
function Grid(gr::Radial.Grid;
    rnt::Union{Nothing,Float64}=nothing,        h::Union{Nothing,Float64}=nothing,      hp::Union{Nothing,Float64}=nothing,
    rbox::Union{Nothing,Float64}=nothing,       orderL::Union{Nothing,Int64}=nothing,
    orderS::Union{Nothing,Int64}=nothing,       orderGL::Union{Nothing,Int64}=nothing,  printout::Bool=false)
    
    if  isnothing(rnt)        rntx      = gr.rnt        else    rntx      = rnt       end 
    if  isnothing(h)          hx        = gr.h          else    hx        = h         end 
    if  isnothing(hp)         hpx       = gr.hp         else    hpx       = hp        end 
    if  isnothing(rbox)       rboxx     = nothing       else    rboxx     = rbox      end 
    if  isnothing(orderL)     orderLx   = gr.orderL     else    orderLx   = orderL    end 
    if  isnothing(orderS)     orderSx   = gr.orderS     else    orderSx   = orderS    end 
    if  isnothing(orderGL)    orderGLx  = gr.orderGL    else    orderGLx  = orderGL   end 
    
    if      isnothing(rboxx)    NoPointsx = gr.NoPoints - rem(gr.NoPoints, orderGLx)
    elseif  rboxx  > 0.         NoPointsx = Radial.determineNoPoints(rntx, hx, hpx, rboxx, orderGLx)
    else    error("Radial.Grid(): rbox = $rboxx must be positive; it is the size of the radial box in " *
                 "atomic units.  Omit the keyword to keep the number of points of the given grid.")
    end
    
    grid = Radial.Grid( GridParameters(rntx, hx, hpx, NoPointsx, orderGLx),
                        KnotSequence(Float64[], Float64[], orderLx, orderSx), RadialMesh(Float64[], Float64[]) )
    return( Radial.determineGrid(grid, printout=printout) )
end


# `Base.show(io::IO, grid::Radial.Grid)`  ... prepares a proper printout of the variable grid::Radial.Grid.
function Base.show(io::IO, grid::Radial.Grid) 

    sa = Base.string(grid::Radial.Grid);    print(io, sa * "\n")

    ## Indexed by the length of the MESH, not by NoPoints, which counts the physical grid of the definition
    ## and is a different quantity even where the two happen to coincide.
    nr = length(grid.r);   ntS = length(grid.tS);    if  nr < 6   return( nothing )    end
    print(io, "r:    ", grid.r[1:3],    "  ...  ", grid.r[nr-2:nr],         "\n")
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
        sa = sa * "orderL = $(grid.orderL),  orderS = $(grid.orderS),  nsL = $(grid.nsL),  nsS = $(grid.nsS), ...  "
    end 
end


"""
`Radial.determineGrid(grid::Radial.Grid; printout::Bool=false)`  
    ... determines the detailed radial grid due to the given gridType and parameters; a gr::Radial.Grid is returned
"""
function determineGrid(grid::Radial.Grid; printout::Bool=false) 
    
    # Read the general parameters from the given grid
    rnt    = grid.rnt;       h      = grid.h;         hp      = grid.hp;        NoPoints = grid.NoPoints
    nth    = grid.orderGL;   orderL = grid.orderL;    orderS  = grid.orderS;    orderGL  = grid.orderGL
    ## rphys is the PHYSICAL grid of the analytic definition below.  It exists only to place the B-spline
    ## break points; the mesh that the returned Grid carries is the Gauss-Legendre one built further down.
    ## Its derivatives rp and rpor used to be computed here and stored in the Grid, where they were then
    ## overwritten with empty arrays for every grid JAC has ever built (12-Aug-2026).
    nr     = grid.NoPoints;  rphys  = zeros(nr);      nsL = 0;   nsS = 0
        
    # Ensure that the largest grid points is always consistent with the largest break point of the B-splines
    if  NoPoints != NoPoints - rem(NoPoints,orderGL)
        error("Radial.determineGrid(): NoPoints = $NoPoints is not a multiple of orderGL = $orderGL, so the " *
              "last grid point cannot coincide with the largest B-spline break point.  Use " *
              "NoPoints = $(NoPoints - rem(NoPoints,orderGL)) instead.")
    end

    # Now define the physical grid due to the given parameters: either an exponential or exponential-linear grid
    if  hp == 0.
        rphys[1] = rnt;   eph = exp(h);    ett = 1.0

        for  i  in 2:nr
            ett      = eph * ett
            ettm1    = ett - 1.0
            rphys[i] = rnt * ettm1
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
    
        rphys[1] = 0.

        rc = 0.
        rn = 0.

        for i = 2:nr
            rn = rn + rnt
            while ((abs((rc - rn)/rn) > 100 * eps(Float64)))
                rc = rn
                rn = rc - f(rc, i) / fprime(rc)
            end
            rphys[i] = rn
        end
    else
        ## Reachable only for hp = NaN, since every other Float64 satisfies one of the two branches above.
        error("Radial.determineGrid(): hp = $hp is not a usable asymptotic step size; use hp = 0 for a " *
              "purely exponential grid or hp > 0 for a log-linear one.")
    end
    
    # Define the radial break points and the number of B-splines
    ntL = orderL;    ntS = orderS;    tL = zeros( orderL );    tS = zeros( orderS )
    for  i = nth:nth:NoPoints
        ntL = ntL + 1;    push!( tL, rphys[i]);    ntS = ntS + 1;    push!( tS, rphys[i])
    end
    for   i = ntL+1:ntL+orderL-1   push!( tL, tL[ntL] )   end;    nsL = ntL - 1;   ntL = length(tL)
    for   i = ntS+1:ntS+orderS-1   push!( tS, tS[ntS] )   end;    nsS = ntS - 1;   ntS = length(tS)
    
    # Lay a Gauss-Legendre mesh over each interval between two break points; this is the mesh the Grid carries
    nr = 0;   r = Float64[];   wr = Float64[];    rlow = 0.
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
        println("Define a radial grid with $nr Gauss-Legendre mesh points")
        println(" [rnt=" * @sprintf("%.3e",rnt) * ", h=" * @sprintf("%.3e",h) *
                ", hp=" * @sprintf("%.3e",hp) * ", NoPoints=$NoPoints, r_max=" * @sprintf("%.3e",r[nr]) * ";")
        println("  B-splines with break points at every $(nth)th point, nsL=$nsL, nsS=$nsS, orderL=$orderL, orderS=$orderS, orderGL=$orderGL,  " *
                "ntL=$ntL, ntS=$ntS] ")
    end

    ## ntL, ntS, nsL and nsS are no longer stored; they follow from the knots.  They are still computed
    ## above because the printout below reports them, and they are asserted against the properties that
    ## replace them, so that the two definitions cannot drift apart unnoticed.
    newGrid = Grid( GridParameters(rnt, h, hp, NoPoints, orderGL), KnotSequence(tL, tS, orderL, orderS),
                    RadialMesh(r, wr) )
    @assert newGrid.ntL == ntL && newGrid.ntS == ntS && newGrid.nsL == nsL && newGrid.nsS == nsS
    return( newGrid )
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
    GridGL(GridGaussLegendreFinite(), 0., 1., 6; printout=false)
end


"""
`Radial.GridGL(::GridGaussLegendreQED, orderGL::Int64; printout::Bool=false)`
    ... constructor to define Gauss-Legendre grid for the typical QED computation in the interval [1.0, infinity).
"""
function GridGL(::GridGaussLegendreQED, orderGL::Int64; printout::Bool=false)
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
`Radial.GridGL(::GridGaussLegendreFinite, tmin::Float64, tmax::Float64, orderGL::Int64; printout::Bool=false)`
    ... constructor to define Gauss-Legendre grid in the finite interval [tmin, tmax].
"""
function GridGL(::GridGaussLegendreFinite, tmin::Float64, tmax::Float64, orderGL::Int64; printout::Bool=false)
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
    
    ## Qualified, and FastGaussQuadrature added to the module's `using` list (12-Aug-2026).  The call read
    ## `gausshermite(orderGH)` with the package named only in a comment, so this constructor raised an
    ## UndefVarError every time it was reached -- it cannot ever have worked.  Three call sites in
    ## module-Pulse.jl reach it.
    wax = FastGaussQuadrature.gausshermite(orderGH)
    
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

    ## Both the head/tail slices and the field name were wrong here (12-Aug-2026).  This method was copied
    ## from Radial.Potential and kept its `Zr`, a field Radial.Density does not have, so printing any density
    ## raised a FieldError; and the guard `n < 6` did not cover the [1:23] slices below it, so 6 <= n < 46
    ## gave a BoundsError.  The window is now taken from n itself.
    n = length(density.Dr);                    if  n < 6   return( nothing )    end
    nw = min(23, div(n, 2))
    print(io, "Dr:    ", density.Dr[1:nw],    "  ...  ", density.Dr[n-nw+1:n],    "\n")
    print(io, density.grid)
end


# `Base.string(density::Radial.Density)`  ... provides a String notation for the variable density::Radial.Density.
function Base.string(density::Radial.Density)
    if  length(density.Dr)  == 0
        sa = "Radial density not yet defined; kind = $(density.name) ..."
    else
        sa = "$(density.name) (radial) density ... defined on $(length(density.Dr)) grid points ..."
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

    ## The guard `n < 6` did not cover the [1:23] slices, so 6 <= n < 46 gave a BoundsError (12-Aug-2026);
    ## the window is now taken from n itself.
    n = length(potential.Zr);                    if  n < 6   return( nothing )    end
    nw = min(23, div(n, 2))
    print(io, "Zr:    ", potential.Zr[1:nw],    "  ...  ", potential.Zr[n-nw+1:n],    "\n")
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
    + grid            ::Radial.Grid       ... explic. defined radial grid for P, Q, if useStandardGrid = false.
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


## Radial.OrbitalBunge1993, Radial.compute_McLean1981 and Radial.OrbitalPrimitiveSlater were REMOVED on
## 12-Aug-2026.  None of the three could run, and none had a caller:
##
##   * the first two called Basics.store("orbital functions: NR, Bunge (1993)" / "... McLean (1981)", Z),
##     and Basics.store is DEFINED NOWHERE IN JAC.  It was meant to hold the Roothaan-Hartree-Fock
##     coefficient tables of Bunge et al., ADNDT 53 (1993) 113 and McLean & McLean, ADNDT 26 (1981) 197;
##     that data has never been part of this repository, so the functions were never more than a sketch;
##   * OrbitalBunge1993 and OrbitalPrimitiveSlater additionally read grid.rp -- an array that was empty on
##     every grid JAC ever built and no longer exists at all -- and constructed Radial.Orbital with SEVEN
##     arguments where it has nine fields.
##
## DECIDED 13-Aug-2026 by the maintainer: the Bunge and McLean routes are RETIRED, not deferred, and the
## compute themes RadialOrbitalBunge1993 and RadialOrbitalMcLean1981 have been removed with them.  Three
## reasons, none of them the missing table itself: both sets are NON-RELATIVISTIC while JAC is a Dirac code,
## and McLean's half covers Z = 55-92, exactly where that matters most; neither reaches beyond Z = 92 or past
## neutral ground configurations, which is most of what JAC actually computes; and Anderson acceleration has
## since removed much of the convergence cost that a better start orbital would have bought.  Should analytic
## start orbitals be wanted, use ManyElectron.StartFromThomasFermi with the screened potential
## Basics.ThomasFermiField (added 13-Aug-2026) -- it needs no external data and works at any Z, including the
## superheavy region these tables cannot reach.
##
## The working route to hydrogenic start orbitals remains Bsplines.generateOrbitalsHydrogenic, which is what
## AsfSettings(..., StartFromHydrogenic(), ...) already uses throughout JAC.







###################################################################################################################
###################################################################################################################
###################################################################################################################


"""
`Radial.determineZbar(pot::Radial.Potential)`  
    ... determines the effective charge that is asymptotically seen by the electron in the potential pot. 
        A Zbar::Float64 is returned.
"""
function determineZbar(pot::Radial.Potential) 
    ## The sampling loop read pot.Zr[mtp] rather than pot.Zr[i], so all nx samples were the SAME number
    ## (12-Aug-2026).  The mean was therefore just Zr[mtp] and the printed spread was identically zero --
    ## a quality indicator that could never fire, whatever the potential did near the outer boundary.
    nx = 5   # Number of grid points for determining the effective charge Zbar
    mtp = size( pot.Zr, 1);    meanZbar = zeros(nx);    devsZbar = zeros(nx);   ny = 0
    for  i = mtp-nx+1:mtp
        ny  = ny + 1;    meanZbar[ny] = pot.Zr[i]
    end
    mZbar   = sum(meanZbar) / nx;   for  i = 1:nx    devsZbar[i] = (meanZbar[i] - mZbar)^2    end
    stdZbar = sqrt( sum(devsZbar) / nx )
    println(">> Radial potential with effective charge Zbar=" * @sprintf("%.4e",mZbar) *
            " (Delta-Zbar=" * @sprintf("%.4e",stdZbar) * ") at r=" * @sprintf("%.4e",pot.grid.r[mtp]) * " a.u." )

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
    
    if  !isnothing(boxSize)
        # Apply a fixed box size (in atomic units) as needed in an average-atom model and elsewhere.
        newGrid = Radial.Grid(grid; rbox=boxSize);    rnt = newGrid.rnt * boxSize / newGrid.tL[end]
        newGrid = Radial.Grid(newGrid; rnt = rnt)
        println(">> Generate a new grid with boxSize = $boxSize a.u. and the break points " *
                "\n   tL = $(newGrid.tL)   \n   tS = $(newGrid.tS)"  *
                "\n>> Note that the last grid point r[end] is always slightly smaller as it refers " *
                "to the last GL (zero of) integration along r.")
        #
    elseif  !isnothing(boxSizeWithZeroWeights)
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
        newGrid = Radial.Grid( getfield(grid, :parameters), getfield(grid, :knots), RadialMesh(grid.r, wr) )
        println(">> Generate a new grid with boxSize = $boxSize a.u. and with all weight grid.wr[j] = 0 " * 
                "for r_j > boxSizeWithZeroWeights (j > $iR)" *
                "\n>> Note that the last non-zero weight at R = $R grid point r[end] is always slightly smaller than " *
                "\n   boxSize = $boxSize as it refers to the last GL (zero of) integration point.")
        #
    elseif  !isnothing(maximumFreeElectronEnergy)
        wavenb      = sqrt( 2maximumFreeElectronEnergy + maximumFreeElectronEnergy * Defaults.getDefaults("alpha")^2 )
        wavelgth    = 2pi / wavenb;     hp = wavelgth / 20
        newGrid     = Radial.Grid(grid; hp=hp)
        println(">> Generate a new grid for the minium wavelength = $wavelgth a.u. of free electrons and  hp = $hp")
        #
    ## The three schemes below are DOCUMENTED BUT NOT IMPLEMENTED; the docstring marks two of them
    ## "(not yet)".  They used to raise error("stop a"/"b"/"c"), which told the caller nothing.  Following
    ## the corePolarization.doApply precedent in module-PhotoEmission.jl, an unimplemented option now
    ## explains itself and names what to use instead (12-Aug-2026).
    elseif  !isnothing(maximumPrincipalQN)
        error("Radial.generateGrid(): the maximumPrincipalQN scheme is not implemented.  It is meant to " *
              "choose rbox = 5 <r_n> for the given n; until it exists, pass that box size directly with " *
              "Radial.Grid(grid; rbox = ...).")
    elseif  !isnothing(NoPointsInsideNucleus)
        error("Radial.generateGrid(): the NoPointsInsideNucleus scheme is not implemented.  It is meant to " *
              "set rnt and h so that the given number of points falls inside the nucleus; until it exists, " *
              "choose rnt and h directly with Radial.Grid(grid; rnt = ..., h = ...).")
    elseif  !isnothing(NoPointsInsideFirstBohrRadius)
        error("Radial.generateGrid(): the NoPointsInsideFirstBohrRadius scheme is not implemented.  It is " *
              "meant to set rnt and h from the number of points inside a_0; until it exists, choose rnt " *
              "and h directly with Radial.Grid(grid; rnt = ..., h = ...).")
    end
        
    return( newGrid )
end


end # module
