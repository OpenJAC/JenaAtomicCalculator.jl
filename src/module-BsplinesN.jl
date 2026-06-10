
"""
`module  JAC.BsplinesN`  
	... a submodel of JAC that contains all structs and methods to generate the B-spline basis and to 
	    solve the single-electron Dirac equation in a local potential. It also provides the major function
        calls to generate self-consistent fields; cf. JAC.SelfConsistent.
    
    !!  Remove in RadialIntegrals:       ... all methods that depend on Bspline
    !!  Remove in InteractionStrength:   ... all methods that depend on Bspline
    !!  Change in Hydrogenic             ... sequence of arguments
    !!  Change in module-Plasma-inc-average-atom  ... sequence of arguments
"""
module BsplinesN


using  BSplineKit, Printf, ..Basics, ..Defaults, ..Nuclear, ..Radial, JenaAtomicCalculator


"""
`struct  BsplinesN.Bspline`  
    ... defines a type for a (single) B-spline that is defined on a given radial grid from r[lower:upper].
        Note that only the non-zero values are specified for the B-spline function and its derivative.

    + lower        ::Int64               ... lower radial index (on the radial grid.r) from where the functions is nonzero.
    + upper        ::Int64               ... upper radial index up to which the functions is nonzero.
    + bs           ::Array{Float64,1}    ... radial B-spline functions as defined on the predefined grid.r[lower:upper]
    + bp           ::Array{Float64,1}    ... derivative of bs on the predefined grid grid.r[lower:upper]
"""
struct Bspline
    lower          ::Int64 
    upper          ::Int64 
    bs             ::Array{Float64,1}   
    bp             ::Array{Float64,1}   
end


"""
`struct  BsplinesN.Primitives`  ... defines a type for a set of primitive functions which typically belongs to a well-defined grid.

    + grid         ::Radial.Grid         ... radial grid on which the states are represented.
    + bsplinesL    ::Array{Bspline,1}    ... set of B-splines for the large components on the given radial grid.
    + bsplinesS    ::Array{Bspline,1}    ... set of B-splines for the small components on the given radial grid.
"""
struct Primitives
    grid           ::Radial.Grid
    bsplinesL      ::Array{Bspline,1}
    bsplinesS      ::Array{Bspline,1}
end


# `Base.show(io::IO, primitives::BsplinesN.Primitives)`  ... prepares a proper printout of the variable BsplinesN.Primitives.
function Base.show(io::IO, primitives::BsplinesN.Primitives) 
    println(io, "grid:               $(primitives.grid)  ")
    println(io, "bsplinesL:           (primitives.bsplinesL)  ")
    println(io, "bsplinesS:           (primitives.bsplinesS)  ")
end


"""
`BsplinesN.computeOverlap(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral <bspline1|bsplines>  for two bpslines as defined on grid.
"""
function computeOverlap(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower
    
    wa = 0.            
    for  i = lower:upper   wa = wa + bspline1.bs[i+add1] * bspline2.bs[i+add2] * grid.wr[i]   end
    return( wa )
end


"""
`BsplinesN.computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)`  
    ... computes the (radial and non-diagonal) D_kappa^+/- integral two the bsplines, all defined on grid
        <bspline1| +/- d/dr + kappa/r | bspline2>. -- pm = +1/-1 provides the phase for taking the derivative.
"""
function computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower
    
    wa = 0.
    for  i = lower:upper  
        wa = wa + pm * bspline1.bs[i+add1] * bspline2.bp[i+add2] * grid.wr[i] 
        if  i == 1  wa = wa + bspline1.bs[i+add1] * kappa * bspline2.bs[i+add2] / (0.3 * grid.r[2]) * grid.wr[i] 
        else        wa = wa + bspline1.bs[i+add1] * kappa * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]   end
    end
    return( wa )
end


"""
`BsplinesN.computeVlocal(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, pot::Radial.Potential, grid::Radial.Grid)`  
    ... computes the (radial) integral <bspline1| V_pot |bsplines>  for two bpslines and the given radial potential 
        as defined on grid.
"""
function computeVlocal(bspline1::BsplinesN.Bspline, bspline2::BsplinesN.Bspline, pot::Radial.Potential, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower
    
    wa = 0.
    for  i = lower:upper  
        if  i == 1  wa = wa - bspline1.bs[i+add1] * pot.Zr[i] * bspline2.bs[i+add2] / (0.3 * grid.r[2]) * grid.wr[i] 
        else        wa = wa - bspline1.bs[i+add1] * pot.Zr[i] * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]   
        end
    end
    return( wa )
end


"""
`BsplinesN.extractBsplineCoefficients(sh::Subshell, wc::Basics.Eigen, grid::Radial.Grid)`  
    ... Here, it is assumed that the matrix wc contains the (column) eigenvectors as associated with a single-electron Dirac 
        Hamiltonian matrix for symmetry kappa in Subshell(n, kappa). The procedure then extracts the (full) vector of B-spline
        coefficients for the radial orbital of subshell sh by applying the standard rules of atomic physics for the principal
        quantum number n.  A  vector::Array{Float64,1}  is returned, whose length is nsL+nsS in the original basis of B-spline 
        functions/primitives.
"""
function extractBsplineCoefficients(sh::Subshell, wc::Basics.Eigen, grid::Radial.Grid)
    nsL = grid.nsL - 2;    nsS = grid.nsS - 2
    l   = Basics.subshell_l(sh);   ni = nsS + sh.n - l;          
    en  = wc.values[ni];        
    ev  = wc.vectors[ni];       if  length(ev) != nsL + nsS    error("stop a")  end
    
    return(ev)
end


"""
`BsplinesN.generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::BsplinesN.Primitives)`  
    ... generates the Galerkin-A matrix for the given potential and B-spline primitives; a matrix::Array{Float64,2} is returned.
"""
function generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::BsplinesN.Primitives)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS

    # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function
    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")
    storage  = Dict{String,Array{Float64,2}}()
    # Set-up the overlap matrix
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = BsplinesN.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = BsplinesN.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    # Set-up the local Hamiltonian matrix
    wa = BsplinesN.setupLocalMatrix(sh.kappa, primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    wa[1:end,1:end] = wa[1:end,1:end] - energy * wb[1:end,1:end]

    return( wa )
end


"""
`BsplinesN.generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::BsplinesN.Primitives; printout::Bool=true)`  
    ... generates all single-electron orbitals from subshell list for the nuclear potential as specified by nm.
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::BsplinesN.Primitives; printout::Bool=true)
    # Extract the requested radial potential from nm
    if       nm.model == "point"    pot = Nuclear.pointNucleus(nm.Z, primitives.grid)
    elseif   nm.model == "Fermi"    pot = Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, primitives.grid) 
    elseif   nm.model == "uniform"  pot = Nuclear.uniformNucleus(nm.radius, nm.Z, primitives.grid)
    else                            error("stop a")
    end
    
    orbitals = BsplinesN.generateOrbitals(subshells, pot, nm, primitives; printout=printout)
    return( orbitals )
end


"""
`BsplinesN.generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model, 
                            primitives::BsplinesN.Primitives; printout::Bool=true)`  
    ... generates all single-electron orbitals from subshell list for the radial potential pot. 
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model, 
                          primitives::BsplinesN.Primitives; printout::Bool=true)
    orbitals = Dict{Subshell, Orbital}()
    kappas   = Int64[];   for sh in subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS
    
    # Define the storage for the calculations of matrices; this is necessary to use the BsplinesN.generateTTpMatrix!() function.
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    for kappa  in  kappas
        # Set-up the overlap matrix
        wb = zeros( nsL+nsS, nsL+nsS )
        
        # (1) Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = BsplinesN.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = BsplinesN.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
        
        # (2) Compute the local Hamiltonian matrix and diagonalize it
        wa = BsplinesN.setupLocalMatrix(kappa, primitives, pot, storage)
        w2 = Basics.diagonalize("generalized eigenvalues: LinearAlgebra", wa, wb)
        nsi = nsS    
        if  printout  Basics.tabulateKappaSymmetryEnergiesDirac(kappa, w2.values, nsi, nm)    end
        
        # (3) Collect all the requested single-electron orbitals
        for  sh in subshells
            if  sh.kappa == kappa    orbitals[sh] = BsplinesN.generateOrbitalFromPrimitives(sh, w2, primitives)    end
        end
    end
    
    return( orbitals )
end


"""
`BsplinesN.extractVectorFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::BsplinesN.Primitives)`  
    ... extracts the B-spline coefficient of the sh orbital from eigenvalues & eigenvectors. 
        A vector::Array{Float64,1} is returned.
"""
function extractVectorFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::BsplinesN.Primitives)
    nsL = primitives.grid.nsL;     nsS = primitives.grid.nsS
    l   = Basics.subshell_l(sh);   ni = nsS + sh.n - l;          if   sh.kappa > 0   ni = ni + 1 - 1  end
    vector = wc.vectors[ni];       if  length(vector) != nsL + nsS    error("stop a")                 end
    
    return( vector )   
end


"""
`BsplinesN.generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::BsplinesN.Primitives)`  
    ... generates the large and small components for the subshell sh from the primitives and their eigenvalues & eigenvectors. 
        A (normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::BsplinesN.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    l  = Basics.subshell_l(sh);   ni = nsS + sh.n - l;          if   sh.kappa > 0   ni = ni + 1 - 1  end
    en = wc.values[ni];        if  en < 0.    isBound = true  else   isBound = false                 end
    ev = wc.vectors[ni];       if  length(ev) != nsL + nsS    error("stop a")                        end
    
    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j+add]      end
    end 
    for  i = 1:nsS   
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end 
    
    # Determine the maximum number of grid points for this orbital and normalized it propery
    mtp = 0;   for j = primitives.grid.NoPoints:-1:1    if  abs(P[j])^2 + abs(Q[j])^2 > 1.0e-13   mtp = j;   break   end     end
    
    Px = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Pprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp]  
    Qx = zeros(mtp);    Qx[1:mtp] = Q[1:mtp];    Qprimex = zeros(mtp);    Qprimex[1:mtp] = Qprime[1:mtp]    
    for  j = 1:mtp      if  abs(Px[j])      < 1.0e-10    Px[j] = 0.       end
                        if  abs(Qx[j])      < 1.0e-10    Qx[j] = 0.       end 
                        if  abs(Pprimex[j]) < 1.0e-10    Pprimex[j] = 0.  end
                        if  abs(Qprimex[j]) < 1.0e-10    Qprimex[j] = 0.  end      end
                        
    # Ensure that the large component of all orbitals start 'positive'
    wSign     = sum( Px[1:30] )
    if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Pprimex[1:mtp] = -Pprimex[1:mtp] 
                     Qx[1:mtp] = -Qx[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end
    
    orbital   = Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid())
    
    # Renormalize the radial orbital   
    wN        = sqrt( JenaAtomicCalculator.RadialIntegrals.overlap(orbital, orbital, primitives.grid) )
    Px[1:mtp] = Px[1:mtp] / wN;    Pprimex[1:mtp] = Pprimex[1:mtp] / wN
    Qx[1:mtp] = Qx[1:mtp] / wN;    Qprimex[1:mtp] = Qprimex[1:mtp] / wN 
    
    orb = Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid())
    
    return( orb )   
end


"""
`BsplinesN.generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::BsplinesN.Primitives)`  
    ... generates the large and small components of a (relativistic) orbital for the subshell sh from the given primitives and the 
        eigenvector ev. A (non-normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::BsplinesN.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL   
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j+add]      end
    end 
    for  i = 1:nsS   
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end 
    
    Px      = zeros(mtp);    Qx      = zeros(mtp);    Px[1:mtp]      = P[1:mtp];         Qx[1:mtp]      = Q[1:mtp]    
    Pprimex = zeros(mtp);    Qprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp];    Qprimex[1:mtp] = Qprime[1:mtp]    
    
    return( Orbital(sh, false, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )   
end


"""
`BsplinesN.generatePrimitives(grid::Radial.Grid)`  
    ... generates the breaks, knots and the B-spline primitives of order k, both for the large and small components. 
        The function applies the given grid parameters; no primitive is defined beyond grid[n_max]. The definition of 
        the primitives follow the work of Zatsarinny and Froese Fischer, CPC 202 (2016) 287. --- A (set of) 
        primitives::BsplinesN.Primitives is returned.
"""
function generatePrimitives(grid::Radial.Grid)
    !(1 <= grid.orderL <= 11)   &&   error("Order should be 2 <= grid.orderL <= 11; obtained order = $(grid.orderL)")
    !(1 <= grid.orderS <= 11)   &&   error("Order should be 2 <= grid.orderS <= 11; obtained order = $(grid.orderS)")
    
    # Now determined the B-splines on the grid for the large and small components; initialize values
    primitivesL = BsplinesN.Bspline[];   primitivesS = BsplinesN.Bspline[];   lower = 0;   upper = 0
    
    # Generate B-spline basis for large component
    breaks = deepcopy( grid.tL[grid.orderL:end-grid.orderL+1] )
    BL = BSplineKit.BSplineBasis(BSplineOrder(grid.orderL), breaks)
    #
    for  (ib, bL)  in  enumerate(BL)
        bs = Float64[];   bp = Float64[];   needlower = true
        for  (ir,r)  in  enumerate(grid.r)  
            if  bL(r) > 0.    push!(bs, bL(r));  push!(bp, bL(r, Derivative(1)) )   
                              upper = ir;        if needlower   lower = ir;   needlower = false   end
            end
        end
        push!(primitivesL, Bspline(lower, upper, bs, bp) )
    end
    
    # Generate B-spline basis for large component
    breaks = deepcopy( grid.tS[grid.orderS:end-grid.orderS+1] )
    BL = BSplineKit.BSplineBasis(BSplineOrder(grid.orderS), breaks)
    #
    for  (ib, bL)  in  enumerate(BL)
        bs = Float64[];   bp = Float64[];   needlower = true
        for  (ir,r)  in  enumerate(grid.r)  
            if  bL(r) > 0.    push!(bs, bL(r));  push!(bp, bL(r, Derivative(1)) )   
                              upper = ir;        if needlower   lower = ir;   needlower = false   end
            end
        end
        push!(primitivesS, Bspline(lower, upper, bs, bp) )
    end
    
    return( BsplinesN.Primitives(grid, primitivesL, primitivesS) )
end


"""
`BsplinesN.generateTTpMatrix!(TTp::String, kappa::Int64, primitives::BsplinesN.Primitives, storage::Dict{String,Array{Float64,2}})`  
    ... returns the TTp block of the (single-electron) Dirac Hamiltonian matrix for an electron with symmetry kappa
        without any potential. The following TTp strings are allowed: ["LL-overlap", "SS-overlap", "LS-D_kappa^-", "LS-D_kappa^+"].
        
        Two modes are distinguished owing to the values that are available in the storage (Dict).
            * The TTp matrix block from the storage is returned, if an entry is known; it is assumed that this matrix
              block belong to the given set of primitives.
            * The TTp matrix is computed and set to the storage otherwise; from the TTp string, the key string
              key = string(kappa) * ":" * TTp is generated an applied in the storage dictionary.
              
        All B-splines are supposed to be defined for the same (radial) grid; a  matrix::Array{Float64,2}  is returned which 
        is quadratic for 'LL-overlap' and 'SS-overlap' and whose dimension depends on the number of B-splines for the large 
        and small component, otherwise.  
"""
function generateTTpMatrix!(TTp::String, kappa::Int64, primitives::BsplinesN.Primitives, storage::Dict{String,Array{Float64,2}})
    # Look up the dictionary of whether the requested matrix has been calculated before
    key = string(kappa) * ":" * TTp;      nsL = primitives.grid.nsL;   nsS = primitives.grid.nsS;
    wc  = Defaults.getDefaults("speed of light: c")
    
    wa  = get( storage, key, zeros(1,1) )
    if  wa != zeros(1,1)  
        ## println(">>>> Re-used $TTp matrix for kappa = $kappa ...")
        return( wa )    
    end
    
    # Now calculate and store the requested matrix
    if      TTp == "LL-overlap"
        wa = zeros( nsL, nsL ) 
        for  i = 1:nsL,  j = 1:nsL
            wa[i,j] = BsplinesN.computeOverlap(primitives.bsplinesL[i], primitives.bsplinesL[j], primitives.grid)
        end
    elseif  TTp == "SS-overlap"
        wa = zeros( nsS, nsS ) 
        for  i = 1:nsS,  j = 1:nsS
            wa[i,j] = BsplinesN.computeOverlap(primitives.bsplinesS[i], primitives.bsplinesS[j], primitives.grid)
         end
    elseif  TTp == "LS-D_kappa^-"
        wa = zeros( nsL, nsS ) 
        for  i = 1:nsL,  j = 1:nsS
            wa[i,j] = wc * BsplinesN.computeNondiagonalD(-1, kappa, primitives.bsplinesL[i], primitives.bsplinesS[j], primitives.grid)
        end
    elseif  TTp == "SL-D_kappa^+"
        wa = zeros( nsS, nsL ) 
        for  i = 1:nsS,  j = 1:nsL
            wa[i,j] = wc * BsplinesN.computeNondiagonalD( 1, kappa, primitives.bsplinesS[i], primitives.bsplinesL[j], primitives.grid)
        end
    else   println("TTp = $TTp ");    error("stop a")
    end
    
    storage[key] = copy(wa)
    return( wa )
end


"""
`BsplinesN.setupLocalMatrix(kappa::Int64, primitives::BsplinesN.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}}) 
        ...set-up the local parts of the generalized eigenvalue problem for the symmetry block kappa and the given (local) potential pot. 
        The B-spline (basis) functions are defined by primitivesL for the large component and primitivesS for the small one, respectively.
"""
function setupLocalMatrix(kappa::Int64, primitives::BsplinesN.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wa  = zeros( nsL+nsS, nsL+nsS );   wb  = zeros( nsL+nsS, nsL+nsS )
    
    # (1) Compute or fetch the diagonal 'overlap' blocks
    wb[1:nsL,1:nsL]                 = BsplinesN.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = BsplinesN.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # (2) Re-compute the diagonal blocks for the local potential
    for  i = 1:nsL,  j = 1:nsL   
        wa[i,j] = BsplinesN.computeVlocal(primitives.bsplinesL[i], primitives.bsplinesL[j], pot, primitives.grid)
    end
    for  i = 1:nsS,  j = 1:nsS    
        wa[nsL+i,nsL+j] = BsplinesN.computeVlocal(primitives.bsplinesS[i], primitives.bsplinesS[j], pot, primitives.grid)
    end
    
    # (3) Substract the rest mass from the 'SS' block
    wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] - 
                                      2 * Defaults.getDefaults("speed of light: c")^2 * wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS]
    
    # (4) Compute or fetch the diagonal 'D_kappa' blocks
    wa[1:nsL,nsL+1:nsL+nsS] = wa[1:nsL,nsL+1:nsL+nsS] + BsplinesN.generateTTpMatrix!("LS-D_kappa^-", kappa, primitives, storage)
    wa[nsL+1:nsL+nsS,1:nsL] = wa[nsL+1:nsL+nsS,1:nsL] + BsplinesN.generateTTpMatrix!("SL-D_kappa^+", kappa, primitives, storage)
    
    #=====
    # Test for 'real-symmetric matrix' ... this is not fullfilled if the last B-spline is included !!
    nx = 0
    for  i = 1:nsL+nsS    
        for  j = i+1:nsL+nsS    
            if  abs(  (wa[i,j] - wa[j,i])/(wa[i,j] + wa[j,i]) ) > 1.0e-7   nx = nx + 1    
                @show "setupLocalMatrix", i, j, wa[i,j], wa[j,i] 
            end
        end
    end
    ny = (nsL+nsS)^2/2 - (nsL+nsS)
    if  nx > 0    
        println(">>> setupLocalMatrix:: $nx (from $(ny)) non-symmetric H-matrix integrals for kappa = $kappa with relative deviation > 1.0e-7.")  end
    =====#
    
    return( wa )
end   


end # module
