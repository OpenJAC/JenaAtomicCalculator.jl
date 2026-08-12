
"""
`module  JAC.Bsplines`
	... a submodel of JAC that contains all structs and methods to generate the B-spline basis and to
	    solve the single-electron Dirac equation in a local potential. It also provides the major function
        calls to generate self-consistent fields; cf. JAC.SelfConsistent.
"""
module Bsplines


using  BSplineKit, Printf, ..Basics, ..Defaults, ..Nuclear, ..Radial, JenaAtomicCalculator


"""
`struct  Bsplines.Bspline`  
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
`struct  Bsplines.Primitives`  ... defines a type for a set of primitive functions which typically belongs to a well-defined grid.

    + grid         ::Radial.Grid         ... radial grid on which the states are represented.
    + bsplinesL    ::Array{Bspline,1}    ... set of B-splines for the large components on the given radial grid.
    + bsplinesS    ::Array{Bspline,1}    ... set of B-splines for the small components on the given radial grid.
"""
struct Primitives
    grid           ::Radial.Grid
    bsplinesL      ::Array{Bspline,1}
    bsplinesS      ::Array{Bspline,1}
end


# `Base.show(io::IO, primitives::Bsplines.Primitives)`  ... prepares a proper printout of the variable Bsplines.Primitives.
function Base.show(io::IO, primitives::Bsplines.Primitives) 
    println(io, "grid:               $(primitives.grid)  ")
    println(io, "bsplinesL:           (primitives.bsplinesL)  ")
    println(io, "bsplinesS:           (primitives.bsplinesS)  ")
end


"""
`Bsplines.computeOverlap(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)`  
    ... computes the (radial) overlap integral <bspline1|bsplines>  for two bpslines as defined on grid.
"""
function computeOverlap(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower
    
    wa = 0.            
    for  i = lower:upper   wa = wa + bspline1.bs[i+add1] * bspline2.bs[i+add2] * grid.wr[i]   end
    return( wa )
end


"""
`Bsplines.computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)`
    ... computes the (radial and non-diagonal) D_kappa^+/- integral two the bsplines, all defined on grid
        <bspline1| +/- d/dr + kappa/r | bspline2>. -- pm = +1/-1 provides the phase for taking the derivative.

        Note (30-Jul-2026): the earlier i==1 special case (substituting the ad hoc 0.3*grid.r[2] for the
        genuine, tiny but nonzero grid.r[1]) had no clear derivation and differed from Zatsarinny &
        Froese Fischer's own DBSR_HF reference code (dbsr_lib_dbs.f90's ZINTYM), which applies a single,
        uniform composite Gauss-Legendre quadrature with no special-casing at any grid point, since every
        grid.r[i] under MeshGL() is already a genuine (nonzero) Gauss-Legendre node, never literally r=0.
        Removed to match; see project_zeeman_hfs_bugs.md (30-Jul-2026).
"""
function computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower

    wa = 0.
    for  i = lower:upper
        wa = wa + pm * bspline1.bs[i+add1] * bspline2.bp[i+add2] * grid.wr[i]
        wa = wa + bspline1.bs[i+add1] * kappa * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]
    end
    return( wa )
end


"""
`Bsplines.computeVlocal(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, pot::Radial.Potential, grid::Radial.Grid)`
    ... computes the (radial) integral <bspline1| V_pot |bsplines>  for two bpslines and the given radial potential
        as defined on grid.

        Note (30-Jul-2026): the earlier i==1 special case (ad hoc 0.3*grid.r[2] substitute) removed for the
        same reason as computeNondiagonalD above -- see that docstring and project_zeeman_hfs_bugs.md.
"""
function computeVlocal(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, pot::Radial.Potential, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower

    wa = 0.
    for  i = lower:upper
        wa = wa - bspline1.bs[i+add1] * pot.Zr[i] * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]
    end
    return( wa )
end





"""
`Bsplines.generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::Bsplines.Primitives)`  
    ... generates the Galerkin-A matrix for the given potential and B-spline primitives; a matrix::Array{Float64,2} is returned.
"""
function generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::Bsplines.Primitives)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS

    # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function
    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")
    storage  = Dict{String,Array{Float64,2}}()
    # Set-up the overlap matrix
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    # Set-up the local Hamiltonian matrix
    wa = Bsplines.setupLocalMatrix(sh.kappa, primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    wa[1:end,1:end] = wa[1:end,1:end] - energy * wb[1:end,1:end]

    return( wa )
end


"""
`Bsplines.generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::Bsplines.Primitives; printout::Bool=true)`  
    ... generates all single-electron orbitals from subshell list for the nuclear potential as specified by nm.
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::Bsplines.Primitives; printout::Bool=true)
    # Extract the requested radial potential from nm
    if       nm.model == "point"    pot = Nuclear.pointNucleus(nm.Z, primitives.grid)
    elseif   nm.model == "Fermi"    pot = Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, primitives.grid) 
    elseif   nm.model == "uniform"  pot = Nuclear.uniformNucleus(nm.radius, nm.Z, primitives.grid)
    else                            error("stop a")
    end
    
    orbitals = Bsplines.generateOrbitals(subshells, pot, nm, primitives; printout=printout)
    return( orbitals )
end


"""
`Bsplines.checkGridRepresentation(subshells::Array{Subshell,1}, Z::Float64, primitives::Bsplines.Primitives;
                                  accuracy::Float64=1.0e-3, stopper::Bool=true)`
    ... checks whether the given radial grid can represent every subshell of the list, by solving the
        single-electron Dirac equation for a POINT nucleus of charge Z on this grid and comparing each level
        with the closed-form point-nucleus energy Basics.computeDiracEnergy(sh, Z). Every subshell whose
        energy deviates by more than `accuracy` is listed, and an error is raised; with stopper = false a
        loud warning is printed instead and the computation proceeds.

        WHY A POINT NUCLEUS IS THE RIGHT YARDSTICK, whatever the computation itself uses. The question asked
        here is not "is this orbital physically accurate" but "can this GRID resolve an orbital of this size
        and shape at all", and only the point-nucleus spectrum has closed-form energies to test against. An
        extended nucleus, or a Dirac-Fock rather than a bare nuclear field, changes each individual orbital
        by a modest factor -- it does not change the ORDER of the spectrum, nor the radial scale that the box
        has to accommodate. So a grid that fails this test will fail the real computation too.

        THE USUAL CAUSE IS A BOX THAT IS TOO LARGE, not one that is too small. The number of B-splines is
        fixed, so a box much wider than the orbitals spends them on empty space and starves the physical
        region: at Z = 10 the default box of 614 a.u. leaves 5f_7/2 wrong by 31%, while a box of 11 a.u.
        -- matched to the orbital -- gives it to 6e-5. A hydrogenic orbital (n,l) has its outer turning
        point at r_plus = (n^2/Z) (1 + sqrt(1 - l(l+1)/n^2)), and a box of roughly 2.5 r_plus is a good
        choice; cf. Radial.Grid(grid; rbox=..).
        A tuple  (isRepresentable::Bool, recommendedRbox::Float64)  is returned.
"""
function checkGridRepresentation(subshells::Array{Subshell,1}, Z::Float64, primitives::Bsplines.Primitives;
                                 accuracy::Float64=1.0e-3, stopper::Bool=true)
    grid    = primitives.grid;      nsL = grid.nsL;     nsS = grid.nsS
    pot     = Nuclear.pointNucleus(Z, grid)
    storage = Dict{String,Array{Float64,2}}()
    wb      = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    #
    offenders = Tuple{Subshell,Float64,Float64,Float64}[];       rbox = 0.
    for  kappa  in  unique( [sh.kappa  for sh in subshells] )
        wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
        w2 = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
        mm = Bsplines.findPositiveBranchStart(w2.values)
        for  sh  in  subshells
            if  sh.kappa != kappa                                        continue    end
            l  = Basics.subshell_l(sh);      ni = mm + sh.n - l - 1
            ex = Basics.computeDiracEnergy(sh, Z)
            if  ni < 1  ||  ni > length(w2.values)   en = NaN;   dev = Inf
            else                                     en = w2.values[ni];   dev = abs(en/ex - 1)    end
            if  dev > accuracy      push!(offenders, (sh, en, ex, dev))                            end
            wr   = (sh.n^2/Z) * (1 + sqrt( max(0., 1 - l*(l+1)/sh.n^2) ));      rbox = max(rbox, 2.5*wr)
        end
    end
    #
    if  length(offenders) > 0
        printstyled("\n>>> GRID CHECK FAILED: on this grid the following subshells are not represented to the requested " *
                    "accuracy of $accuracy\n>>> in the pure (point-nucleus, Z = $Z) Dirac spectrum:\n", color=:light_red)
        printstyled("      subshell      E(grid) [a.u.]      E(Dirac) [a.u.]     rel. deviation\n", color=:light_red)
        for  (sh, en, ex, dev)  in  offenders
            printstyled(@sprintf("    %10s    %+.8e     %+.8e      %.2e\n", string(sh), en, ex, dev), color=:light_red)
        end
        printstyled(@sprintf(">>> The present box is r_max = %.1f a.u.;  a box of about %.1f a.u. suits these subshells.\n",
                             grid.r[end], rbox), color=:light_red)
        printstyled(">>> Note that a box which is much TOO LARGE starves the basis just as badly as one that is too\n" *
                    ">>> small, since the number of B-splines is fixed;  use Radial.Grid(grid; rbox=..) to match it.\n",
                    color=:light_red)
        if  stopper   error("Bsplines.checkGridRepresentation(): the grid fails to represent " *
                            "$(length(offenders)) of $(length(subshells)) subshells to accuracy $accuracy.")   end
        return( (false, rbox) )
    end

    return( (true, rbox) )
end


"""
`Bsplines.generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model, 
                            primitives::Bsplines.Primitives; printout::Bool=true)`  
    ... generates all single-electron orbitals from subshell list for the radial potential pot. 
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model, 
                          primitives::Bsplines.Primitives; printout::Bool=true)
    orbitals = Dict{Subshell, Orbital}()
    kappas   = Int64[];   for sh in subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS
    
    # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateTTpMatrix!() function.
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    for kappa  in  kappas
        # Set-up the overlap matrix
        wb = zeros( nsL+nsS, nsL+nsS )
        
        # (1) Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
        
        # (2) Compute the local Hamiltonian matrix and diagonalize it
        wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
        w2 = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
        ## The offset handed to the tabulation is where the POSITIVE-energy branch begins, and that is not nsS
        ## (corrected 10-Aug-2026). diagonalizeLocalMatrix eliminates dropP+dropQ B-splines for this symmetry,
        ## so the branch starts at findPositiveBranchStart -- 46..61 on the standard grid against nsS = 63.
        ## With nsS the table printed a slice several states too high while labelling it 4f, 5f, ...: at Z = 10
        ## it showed 4f_5/2 = -6.03e-01 against the exact -3.13e+00, i.e. a "Delta-E/|E|" of +4.19, although the
        ## orbital actually extracted was correct to 6e-7. Purely a display defect, but a badly misleading one.
        nsi = Bsplines.findPositiveBranchStart(w2.values) - 1
        if  printout  Basics.tabulateKappaSymmetryEnergiesDirac(kappa, w2.values, nsi, nm)    end
        
        # (3) Collect all the requested single-electron orbitals
        for  sh in subshells
            if  sh.kappa == kappa    orbitals[sh] = Bsplines.generateOrbitalFromPrimitives(sh, w2, primitives)    end
        end
    end
    
    return( orbitals )
end





"""
`Bsplines.fitVectorToPrimitivesClaude(orb::Radial.Orbital, primitives::Bsplines.Primitives, matrixB::Array{Float64,2})`
    ... projects the (already CLEANED/truncated) tabulated orbital orb onto the B-spline primitives basis via
        the standard Galerkin/least-squares projection matrixB * p = rhs, rhs[i] = <B_i|orb.P-or-Q>, using the
        existing grid quadrature. Unlike pulling the RAW diagonalization eigenvector -- which reproduces the
        UNCLEANED tabulated function from BEFORE generateOrbitalFromPrimitives' own truncation-at-mtp and
        small-value cleanup are applied (the route the removed Bsplines.extractVectorFromPrimitives took) -- this
        function GUARANTEES the returned coefficient vector is fully self-consistent with orb's OWN (already
        cleaned) tabulated P, Q arrays, inheriting orb's own well-defined truncation instead of carrying whatever
        small numerical noise the raw eigenvector's tail coefficients happen to have.
        This matters wherever B-spline expansion coefficients are themselves summed/weighted directly, rather
        than only ever used to reconstruct one smooth tabulated function -- such sums do not automatically
        benefit from the cancellation that evaluating a single, already-cleaned tabulated function enjoys. See
        InteractionStrength.XL_CoulombTensorClaude and SelfConsistent.solveAverageLevelField, where using
        the raw eigenvector instead of this projection was traced to a real, non-negligible SCF discrepancy.
        The returned vector is explicitly re-normalized so that v'*matrixB*v = 1 EXACTLY (to floating-point
        precision), rather than trusting the least-squares fit to land there on its own -- Hamiltonian.
        projectHamiltonian's projection operator (I - S*bb') is only truly idempotent for an exactly
        S-normalized b; feeding it a vector off by even a small residual leaves the "projection" not quite a
        projection, which compounds under repeated application across SCF iterations.
        A vector::Vector{Float64}, of length nsL+nsS, is returned.
"""
function fitVectorToPrimitivesClaude(orb::Radial.Orbital, primitives::Bsplines.Primitives, matrixB::Array{Float64,2})
    nsL = primitives.grid.nsL;   nsS = primitives.grid.nsS;   grid = primitives.grid
    rhs = zeros(nsL+nsS)

    for  i = 1:nsL
        Bi = primitives.bsplinesL[i]
        Pi = zeros(Bi.upper);   add = 1 - Bi.lower
        for  j = Bi.lower:Bi.upper   Pi[j] = Pi[j] + Bi.bs[j+add]   end
        mtp = min(length(Pi), length(orb.P))
        s   = 0.;   for  r = 2:mtp   s = s + Pi[r] * orb.P[r] * grid.wr[r]   end
        rhs[i] = s
    end
    for  i = 1:nsS
        Bi = primitives.bsplinesS[i]
        Qi = zeros(Bi.upper);   add = 1 - Bi.lower
        for  j = Bi.lower:Bi.upper   Qi[j] = Qi[j] + Bi.bs[j+add]   end
        mtp = min(length(Qi), length(orb.Q))
        s   = 0.;   for  r = 2:mtp   s = s + Qi[r] * orb.Q[r] * grid.wr[r]   end
        rhs[nsL+i] = s
    end

    vector = matrixB \ rhs
    norm2  = transpose(vector) * matrixB * vector
    return( vector / sqrt(norm2) )
end


"""
`Bsplines.generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives)`  
    ... generates the large and small components for the subshell sh from the primitives and their eigenvalues & eigenvectors. 
        A (normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    l  = Basics.subshell_l(sh);   mm = Bsplines.findPositiveBranchStart(wc.values);   ni = mm + sh.n - l - 1
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
`Bsplines.generateOrbitalFromVectorClaude(sh::Subshell, energy::Float64, vector::Vector{Float64},
                                          primitives::Bsplines.Primitives)`
    ... generates a (normalized, cleaned) tabulated orbital directly from a given B-spline expansion coefficient
        vector, rather than from an eigenvector INDEXED out of a Basics.Eigen (as
        Bsplines.generateOrbitalFromPrimitives(sh,wc,primitives) requires). This is needed whenever the vector in
        hand is not literally an eigenvector of anything -- e.g. after SelfConsistent.
        orthonormalizeSameKappaClaude's Loewdin symmetric orthogonalization, which produces a LINEAR COMBINATION
        of eigenvectors that is itself not an eigenvector of the original problem.
        Reproduces generateOrbitalFromPrimitives(sh,wc,primitives)'s own reconstruction exactly (auto-detects
        mtp from where the density drops below 1e-13, zeros values below 1e-10, fixes the sign convention so the
        large component starts positive, and renormalizes to unit norm in the grid quadrature) -- just without
        the wc/ni lookup indirection.
        A (normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromVectorClaude(sh::Subshell, energy::Float64, vector::Vector{Float64}, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    if  length(vector) != nsL + nsS    error("stop a")    end

    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + vector[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + vector[i] * primitives.bsplinesL[i].bp[j+add]      end
    end
    for  i = 1:nsS
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + vector[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + vector[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end

    mtp = 0;   for j = primitives.grid.NoPoints:-1:1    if  abs(P[j])^2 + abs(Q[j])^2 > 1.0e-13   mtp = j;   break   end     end

    Px = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Pprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp]
    Qx = zeros(mtp);    Qx[1:mtp] = Q[1:mtp];    Qprimex = zeros(mtp);    Qprimex[1:mtp] = Qprime[1:mtp]
    for  j = 1:mtp      if  abs(Px[j])      < 1.0e-10    Px[j] = 0.       end
                        if  abs(Qx[j])      < 1.0e-10    Qx[j] = 0.       end
                        if  abs(Pprimex[j]) < 1.0e-10    Pprimex[j] = 0.  end
                        if  abs(Qprimex[j]) < 1.0e-10    Qprimex[j] = 0.  end      end

    wSign = sum( Px[1:min(30,mtp)] )
    if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Pprimex[1:mtp] = -Pprimex[1:mtp]
                     Qx[1:mtp] = -Qx[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end

    orbital = Orbital(sh, true, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid())

    wN        = sqrt( JenaAtomicCalculator.RadialIntegrals.overlap(orbital, orbital, primitives.grid) )
    Px[1:mtp] = Px[1:mtp] / wN;    Pprimex[1:mtp] = Pprimex[1:mtp] / wN
    Qx[1:mtp] = Qx[1:mtp] / wN;    Qprimex[1:mtp] = Qprimex[1:mtp] / wN

    return( Orbital(sh, true, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )
end


"""
`Bsplines.generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)`
    ... generates the large and small components of a (relativistic) orbital for the subshell sh from the given primitives and the
        eigenvector ev. A (non-normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)
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
`Bsplines.generatePrimitives(grid::Radial.Grid)`
    ... generates the breaks, knots and the B-spline primitives of order k, both for the large and small components.
        The function applies the given grid parameters; no primitive is defined beyond grid[n_max]. The definition of
        the primitives follow the work of Zatsarinny and Froese Fischer, CPC 202 (2016) 287. --- A (set of)
        primitives::Bsplines.Primitives is returned.

        Note (30-Jul-2026): the raw, full B-spline count (including the corner spline that is nonzero at
        r=0) is kept here on purpose -- the r=0 boundary condition is instead enforced per symmetry kappa,
        with the correct (l- and kappa-sign-dependent) number of leading splines, inside
        Bsplines.setupLocalMatrix -- see Bsplines.boundaryDropCounts and project_zeeman_hfs_bugs.md
        (30-Jul-2026), following the same approach as Zatsarinny & Froese Fischer's DBSR_HF reference code
        (hf_boundary.f90's Boundary_conditions, for a finite/non-point nucleus).
"""
function generatePrimitives(grid::Radial.Grid)
    !(1 <= grid.orderL <= 11)   &&   error("Order should be 2 <= grid.orderL <= 11; obtained order = $(grid.orderL)")
    !(1 <= grid.orderS <= 11)   &&   error("Order should be 2 <= grid.orderS <= 11; obtained order = $(grid.orderS)")

    # Now determined the B-splines on the grid for the large and small components; initialize values
    primitivesL = Bsplines.Bspline[];   primitivesS = Bsplines.Bspline[];   lower = 0;   upper = 0

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

    return( Bsplines.Primitives(grid, primitivesL, primitivesS) )
end


"""
`Bsplines.boundaryDropCounts(kappa::Int64, grid::Radial.Grid)`
    ... returns (dropP, dropQ, trailP, trailQ), the number of leading (near r=0) and trailing (near r=R_max)
        large- and small-component B-splines that must be excluded from the generalized eigenvalue problem
        for the given symmetry kappa, following Zatsarinny & Froese Fischer's DBSR_HF reference code
        (hf_boundary.f90's Boundary_conditions, the finite/non-point-nucleus branch), cross-checked directly
        against a live DBSR_HF run on Sc (30-Jul-2026, temporary diagnostic dump of its iprm mask): P needs
        l+1 excluded at r=0 (kappa-sign independent, since the large component's leading near-origin power
        depends only on l); Q needs l+2 excluded for kappa<0 (j=l+1/2, where the naive leading power
        l+1+kappa vanishes identically since kappa=-(l+1), pushing the true leading power one order higher)
        but only l for kappa>0 (j=l-1/2, kappa=l, no such cancellation) -- both capped at (order-1),
        matching DBSR_HF's own j>ksp-1/j>ksq-1 caps. The trailing (outer-boundary) counts are
        kappa/l-INDEPENDENT in DBSR_HF (confirmed: identical trailP=3, trailQ=2 for every orbital in the
        live Sc run) and tied to a small fixed parameter DBSR_HF calls ibzero (=2 by default): trailP =
        ibzero+1, trailQ = ibzero -- omitting this outer truncation (an earlier, incomplete version of this
        fix) left the DIAGONALIZED spectrum subtly perturbed even for kappa>0 states that were already
        correct before any boundary condition was applied at all. See project_zeeman_hfs_bugs.md
        (30-Jul-2026).
"""
function boundaryDropCounts(kappa::Int64, grid::Radial.Grid)
    l      = kappa > 0  ?  kappa  :  -kappa - 1
    dropP  = min(l + 1,                     grid.orderL - 1)
    dropQ  = kappa < 0  ?  min(l + 2, grid.orderS - 1)  :  min(l, grid.orderS - 1)
    ibzero = 2
    trailP = ibzero + 1
    trailQ = ibzero
    return( dropP, dropQ, trailP, trailQ )
end


"""
`Bsplines.findPositiveBranchStart(values::Array{Float64,1})`
    ... returns the 1-based index of the first (ascending-sorted) eigenvalue above -1.999*c^2, the standard
        threshold marking the boundary between the unphysical negative-energy continuum ("Dirac sea")
        branch and the physical positive branch, following Zatsarinny & Froese Fischer's DBSR_HF reference
        code (hf_solve_HF.f90's hf_eiv). The (n,l) bound state then sits at (this index) + (n-l) - 1 within
        the SAME values array -- robust to however many leading B-splines were eliminated for a given
        symmetry kappa (Bsplines.boundaryDropCounts/diagonalizeLocalMatrix), unlike a fixed index counted
        from nsL/nsS. See project_zeeman_hfs_bugs.md (30-Jul-2026).

    WARNING (3-Aug-2026): this threshold-based separation silently returns garbage if `values` was
    diagonalized against a potential with no (or a much-too-weak) attractive nuclear well -- e.g. a
    caller that passes only `Basics.computePotential(Basics.DFSField(1.0), grid, basis)` (the ELECTRONIC
    mean-field potential alone) without adding `Nuclear.nuclearPotential(nm, grid)` first. Without a real
    potential well, the Dirac equation has no clean energetic gap between the unphysical Dirac-sea branch
    and genuine atomic bound states, so this function can return an index that is STILL within (or
    immediately adjacent to) the spurious negative-continuum branch -- e.g. eigenvalues clustering right at
    -1.999*c^2 instead of the expected atomic scale (roughly -1 to -2000 Hartree, not ~-37500 for a typical
    ion). This silently produces orbitals that are numerical garbage, not "slightly wrong" -- found via
    `module-InternalRecombination.jl`, which was missing exactly this nuclear term (fixed there); always
    add the nuclear potential to any potential passed into `Bsplines.generateOrbitals`/this function. See
    project_bsplines_spurious_dirac_sea_bug.md for the full diagnostic.
"""
function findPositiveBranchStart(values::Array{Float64,1})
    c  = Defaults.getDefaults("speed of light: c")
    zz = -1.999 * c^2
    for  (i,v)  in  enumerate(values)
        if  v > zz   return( i )   end
    end
    error("Bsplines.findPositiveBranchStart(): no eigenvalue found above the negative-continuum threshold.")
end


"""
`Bsplines.diagonalizeLocalMatrix(kappa::Int64, matrixA::Array{Float64,2}, matrixB::Array{Float64,2}, primitives::Bsplines.Primitives)`
    ... enforces the r=0 boundary condition for symmetry kappa (Bsplines.boundaryDropCounts) by eliminating
        the corresponding leading large- and small-component B-splines from the generalized eigenvalue
        problem (matrixA, matrixB) before diagonalizing, following Zatsarinny & Froese Fischer's DBSR_HF
        reference code (hf_boundary.f90's Boundary_conditions, finite-nucleus branch). Each returned
        eigenvector is re-embedded into the full (nsL+nsS)-length vector space, with exact zeros at the
        eliminated positions, so that all downstream code (Bsplines.generateOrbitalFromPrimitives, the
        bVector-native AL/EOL machinery, etc.) is unaffected in how it reconstructs P(r)/Q(r) -- only the
        NUMBER of returned eigenpairs shrinks by (dropP+dropQ). Callers that pick an eigenpair by index must
        use Bsplines.findPositiveBranchStart on the returned values, NOT a fixed index counted from
        nsL/nsS -- see project_zeeman_hfs_bugs.md (30-Jul-2026). An eigen::Basics.Eigen is returned.
"""
function diagonalizeLocalMatrix(kappa::Int64, matrixA::Array{Float64,2}, matrixB::Array{Float64,2}, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    dropP, dropQ, trailP, trailQ = Bsplines.boundaryDropCounts(kappa, primitives.grid)
    keep = vcat( collect(dropP+1:nsL-trailP), collect(nsL+dropQ+1:nsL+nsS-trailQ) )
    #
    wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrixA[keep,keep], matrixB[keep,keep])
    #
    vectors = Vector{Float64}[]
    for  v  in  wc.vectors
        full = zeros(nsL+nsS);   full[keep] = v;   push!(vectors, full)
    end
    return( Basics.Eigen(wc.values, vectors) )
end


"""
`Bsplines.generateTTpMatrix!(TTp::String, kappa::Int64, primitives::Bsplines.Primitives, storage::Dict{String,Array{Float64,2}})`  
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
function generateTTpMatrix!(TTp::String, kappa::Int64, primitives::Bsplines.Primitives, storage::Dict{String,Array{Float64,2}})
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
            wa[i,j] = Bsplines.computeOverlap(primitives.bsplinesL[i], primitives.bsplinesL[j], primitives.grid)
        end
    elseif  TTp == "SS-overlap"
        wa = zeros( nsS, nsS ) 
        for  i = 1:nsS,  j = 1:nsS
            wa[i,j] = Bsplines.computeOverlap(primitives.bsplinesS[i], primitives.bsplinesS[j], primitives.grid)
         end
    elseif  TTp == "LS-D_kappa^-"
        wa = zeros( nsL, nsS ) 
        for  i = 1:nsL,  j = 1:nsS
            wa[i,j] = wc * Bsplines.computeNondiagonalD(-1, kappa, primitives.bsplinesL[i], primitives.bsplinesS[j], primitives.grid)
        end
    elseif  TTp == "SL-D_kappa^+"
        wa = zeros( nsS, nsL ) 
        for  i = 1:nsS,  j = 1:nsL
            wa[i,j] = wc * Bsplines.computeNondiagonalD( 1, kappa, primitives.bsplinesS[i], primitives.bsplinesL[j], primitives.grid)
        end
    else   println("TTp = $TTp ");    error("stop a")
    end
    
    storage[key] = copy(wa)
    return( wa )
end


"""
`Bsplines.setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}}) 
        ...set-up the local parts of the generalized eigenvalue problem for the symmetry block kappa and the given (local) potential pot. 
        The B-spline (basis) functions are defined by primitivesL for the large component and primitivesS for the small one, respectively.
"""
function setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wa  = zeros( nsL+nsS, nsL+nsS );   wb  = zeros( nsL+nsS, nsL+nsS )
    
    # (1) Compute or fetch the diagonal 'overlap' blocks
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # (2) Re-compute the diagonal blocks for the local potential
    for  i = 1:nsL,  j = 1:nsL   
        wa[i,j] = Bsplines.computeVlocal(primitives.bsplinesL[i], primitives.bsplinesL[j], pot, primitives.grid)
    end
    for  i = 1:nsS,  j = 1:nsS    
        wa[nsL+i,nsL+j] = Bsplines.computeVlocal(primitives.bsplinesS[i], primitives.bsplinesS[j], pot, primitives.grid)
    end
    
    # (3) Substract the rest mass from the 'SS' block
    wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = wa[nsL+1:nsL+nsS,nsL+1:nsL+nsS] - 
                                      2 * Defaults.getDefaults("speed of light: c")^2 * wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS]
    
    # (4) Compute or fetch the diagonal 'D_kappa' blocks
    wa[1:nsL,nsL+1:nsL+nsS] = wa[1:nsL,nsL+1:nsL+nsS] + Bsplines.generateTTpMatrix!("LS-D_kappa^-", kappa, primitives, storage)
    wa[nsL+1:nsL+nsS,1:nsL] = wa[nsL+1:nsL+nsS,1:nsL] + Bsplines.generateTTpMatrix!("SL-D_kappa^+", kappa, primitives, storage)
    
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


"""
`Bsplines.checkOrbitalConsistency(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                                  rTolerance::Float64=1.5, eTolerance::Float64=2.0, stopper::Bool=true)`
    ... checks the generated orbitals for the one failure mode that Bsplines.checkGridRepresentation cannot
        see: an SCF that has converged onto the WRONG STATE for some symmetry. It compares the two
        spin-orbit partners of every subshell -- same n, same l, kappa of either sign -- which must describe
        the same shell and therefore must be close in both mean radius and binding energy.

        WHY THIS IS NEEDED IN ADDITION TO checkGridRepresentation. That function tests hydrogenic orbitals at
        the FULL nuclear charge, so it only sees whether the grid can resolve a COMPACT orbital. It passes
        happily on a grid that cannot represent a diffuse, screened outer orbital: Ge II [Ar] 3d^10 4s^2 4f
        on r_max = 614 a.u. returned E(4f_7/2) = -1.5758 with <r> = 5.20 against E(4f_5/2) = -0.0619 with
        <r> = 10.87 -- a different state entirely -- and produced a Lande factor of -2.264 against the exact
        8/7, while the grid check reported no problem at all.

        THE TOLERANCES ARE CALIBRATED, not guessed. Genuine fine structure does separate the partners, and
        the more so the heavier the ion, so a tight criterion would fire on correct results. Measured for
        hydrogen-like ions on matched boxes, |E_1/E_2| and <r>_1/<r>_2 are
             Z =  10   2p  1.0013 / 0.9982      3d  1.0003 / 0.9996
             Z =  26   2p  1.0092 / 0.9877      3d  1.0020 / 0.9969
             Z =  54   2p  1.0426 / 0.9448      3d  1.0088 / 0.9867
             Z =  92   2p  1.1540 / 0.8195      3d  1.0268 / 0.9601      4f  1.0097 / 0.9844
        i.e. at worst 15% in energy and 18% in radius for a legitimate pair, against a factor 26 in energy
        and 2.1 in radius for the broken Ge II case. The defaults of 2.0 and 1.5 sit in that gap with room
        on both sides.
        A value::Bool is returned -- true if every partner pair is consistent.
"""
function checkOrbitalConsistency(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                                 rTolerance::Float64=1.5, eTolerance::Float64=2.0,
                                 eInversion::Float64=0.05, stopper::Bool=true)
    ## group the orbitals by (n, l); only pairs with both kappa signs can be compared
    groups = Dict{Tuple{Int64,Int64}, Array{Subshell,1}}()
    for  sh  in  keys(orbitals)
        key = (sh.n, Basics.subshell_l(sh));    groups[key] = push!( get(groups, key, Subshell[]), sh )
    end
    #
    offenders = Tuple{Subshell,Subshell,Float64,Float64}[]
    for  (key, shs)  in  groups
        length(shs) == 2   ||   continue
        a = orbitals[shs[1]];     b = orbitals[shs[2]]
        ra = JenaAtomicCalculator.RadialIntegrals.rkDiagonal(1, a, a, grid)
        rb = JenaAtomicCalculator.RadialIntegrals.rkDiagonal(1, b, b, grid)
        rRatio = (ra > 0. && rb > 0.)  ?  max(ra,rb)/min(ra,rb)  :  Inf
        eRatio = (a.energy * b.energy > 0.)  ?  max(abs(a.energy),abs(b.energy))/min(abs(a.energy),abs(b.energy))  :  Inf
        ## The MEAN RADIUS is one criterion; the energy RATIO is reported for diagnosis but does not trigger
        ## on its own. Energies are fragile here: a nearly-unbound subshell in a highly-ionised configuration
        ## can have its two partners straddle zero, which makes eRatio infinite while the orbitals are in
        ## fact identical -- exactly what the Cascade stepwise-decay test does with 3p_1/2 / 3p_3/2, whose
        ## radii agree to 0.1%. The radius alone already separates the calibration cases by a wide margin.
        ##
        ## THE ORDERING, however, IS a criterion, and a sharp one (added 12-Aug-2026). In a central field the
        ## two members of a spin-orbit pair are ordered by physics, not by magnitude: j = l-1/2 (kappa > 0)
        ## lies BELOW j = l+1/2 (kappa < 0). Hydrogen's 2p_1/2 below 2p_3/2 is the familiar case. An inverted
        ## pair therefore means the SCF has produced something that is not a spin-orbit doublet at all, no
        ## matter how similar the two radii happen to be.
        ## THIS IS THE CASE THE RADIUS TEST MISSED: neutral Pr I [Xe] 4f^3 6s^2 gave E(4f_7/2) = -0.1854 below
        ## E(4f_5/2) = -0.1154 -- inverted, and 20x too large -- while the two mean radii agreed to 6%, well
        ## inside rTolerance. The level structure that followed had the ^4I multiplet upside down and its fine
        ## structure collapsed by four orders of magnitude, and nothing warned.
        ## The test is applied only when BOTH partners are bound, which is what keeps the straddling-zero
        ## Cascade case out of it; there, one energy is positive and the ordering carries no meaning.
        inverted = 0.
        if  a.energy < 0.  &&  b.energy < 0.
            shLo = shs[1].kappa > 0  ?  shs[1]  :  shs[2]      ## j = l - 1/2, must be the more bound
            shHi = shs[1].kappa > 0  ?  shs[2]  :  shs[1]      ## j = l + 1/2
            eLo  = orbitals[shLo].energy;      eHi = orbitals[shHi].energy
            inverted = (eLo - eHi) / max(abs(eLo), abs(eHi))   ## > 0 means the pair is upside down
        end
        if  rRatio > rTolerance   ||   inverted > eInversion
            push!(offenders, (shs[1], shs[2], eRatio, rRatio))
        end
    end
    #
    if  length(offenders) > 0
        printstyled("\n>>> ORBITAL CHECK FAILED: these spin-orbit partners do not describe the same shell, which means\n" *
                    ">>> the SCF has converged onto the wrong state for one of them (usually an ill-matched radial box):\n",
                    color=:light_red)
        printstyled("      partners                    E ratio     <r> ratio     (limits $eTolerance / $rTolerance;\n" *
                    "      an inverted spin-orbit pair also triggers, at $eInversion)\n", color=:light_red)
        for  (sa, sb, er, rr)  in  offenders
            printstyled(@sprintf("    %-10s %-10s  %11.4g %13.4g\n", string(sa), string(sb), er, rr), color=:light_red)
        end
        printstyled(">>> Compare the two orbital energies and mean radii directly, and match the radial box to the\n" *
                    ">>> orbitals; note that a box much TOO LARGE starves the basis just as badly as one too small.\n",
                    color=:light_red)
        if  stopper   error("Bsplines.checkOrbitalConsistency(): $(length(offenders)) spin-orbit partner pair(s) " *
                            "describe different states; the orbitals are not trustworthy.")   end
        return( false )
    end

    return( true )
end


end # module
