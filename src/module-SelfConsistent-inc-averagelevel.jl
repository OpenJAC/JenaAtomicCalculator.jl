

# The AL (average-level) path: its spin-angular coefficients and its solver.

"""
`SelfConsistent.computeAngularCoefficients(scField::Basics.ALField, basis::Basis)` 
    ... computes all spin-angular coefficients for the average-level (AL) functional as obtained for the given
        basis. In the AL functional, only the diagonal matrix elements (from the trace of the Hamiltonian matrix)
        are considered, and the energy is averaged with regard to the number of CSF in the basis. This makes the total 
        energy of the system comparable and independent of the size of the basis. A Tuple of two lists with one- and 
        two-particle coefficients  tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}} 
        is returned.
"""
function computeAngularCoefficients(scField::Basics.ALField, basis::Basis)
    ncsf = length(basis.csfs);    coeffs1p = Coefficient1p[];     coeffs2p = Coefficient2p[] 
    
    # Compute angular coefficients in turn for all diagonal ME
    for  csf  in  basis.csfs
        coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus), 
                                                       csf, csf, basis.subshells)
        # Add to the existing list
        for  cf in coeffs   push!(coeffs1p, Coefficient1p(cf.nu, cf.a, cf.b, cf.T / ncsf) )   end
        
        coeffs = SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus), 
                                                       csf, csf, basis.subshells)
        # Add to the existing lists
        for  cf in coeffs   push!(coeffs2p, Coefficient2p(cf.nu, cf.a, cf.b, cf.c, cf.d, cf.V / ncsf) )   end
    end 

    # Condense angular coefficients if they refer to the same set of orbital; 
    # include symmetry <ab||cd> == <ba||dc> for symmetric interactions
    coeffs1px = Coefficient1p[];     coeffs2px = Coefficient2p[]
    
    hasConsidered = falses( length(coeffs1p) );   T = 0.
    for  (ic, cf) in enumerate(coeffs1p)
        if    hasConsidered[ic]   
        else  nu = cf.nu;   a = cf.a;   b = cf.b
            T = T + cf.T;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs1p)
                if    hasConsidered[icx]   
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cfx.b    T = T + cfx.T;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs1px, Coefficient1p(nu, a, b, T) );  T = 0.
        end 
    end
    
    hasConsidered = falses( length(coeffs2p) );   V = 0.
    for  (ic, cf) in enumerate(coeffs2p)
        if    hasConsidered[ic]   
        else  
            nu = cf.nu;       a = cf.a;   b = cf.b;   c = cf.c;   d = cf.d
            V  = V + cf.V;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs2p)
                if    hasConsidered[icx]   
                elseif  nu == cfx.nu &&    a == cfx.a  &&  b == cfx.b  &&  c == cfx.c &&  d == cfx.d    
                        V  = V + cfx.V;    hasConsidered[icx] = true
                # elseif  nu == cfx.nu &&  a == cfx.b && b == cf.a &&  c == cfx.d &&  d == cf.c    
                #         V = V + cfx.V;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs2px, Coefficient2p(nu, a, b, c, d, V) );   V = 0.
        end
    end

    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.sameKappaRotationPairs(coeffs1p::Array{Coefficient1p,1}, subshells::Array{Subshell,1},
                                       frozenSubshells::Array{Subshell,1})`
    ... selects the pairs of subshells between which the average-level functional has a genuine rotational
        degree of freedom:  equal kappa, and UNEQUAL generalised occupation.  An `Array{Tuple{Subshell,Subshell},1}`
        is returned, empty where the configuration has no such pair -- which is the case for every closed shell.

        THE OCCUPATION IS READ OFF THE RANK-0 DIAGONAL COEFFICIENT, which in JAC's convention IS the generalised
        occupation (see Rule 18:  the sqrt(2j+1) sits inside the coefficient, so a rank-0 diagonal coefficient is
        literally the occupation number).  Verified on Li-like Z = 26 over two configurations, where the
        coefficients give 1s = 2.000000, 2s = 0.333333, 2p_1/2 = 0.333333, 2p_3/2 = 0.333333 and the occupation
        averaged over the CSFs gives the same to every digit.  No separate machinery is needed, and `coeffs1p` is
        already in the caller's hand.

        WHY UNEQUAL OCCUPATIONS ARE THE TEST.  With EQUAL occupations the functional is invariant under the
        rotation -- the standard Hartree-Fock invariance within a symmetry block -- so the angle is not a degree
        of freedom at all and optimizing it is both pointless and ill-conditioned, the second derivative being
        zero.  Measured:  Be-like 1s^2 2s^2 gives dF/dtheta = 5.0e-07 Ha/rad against Li-like 1s^2 2s + 1s^2 2p's
        -0.17, i.e. 340,000 times flatter.  A FROZEN subshell is excluded for the reason the caller excludes it
        from Anderson acceleration:  the rotation would move an orbital the caller has pinned.
"""
function sameKappaRotationPairs(coeffs1p::Array{Coefficient1p,1}, subshells::Array{Subshell,1},
                                frozenSubshells::Array{Subshell,1})
    occupation = Dict{Subshell,Float64}()
    for  cf  in  coeffs1p
        if  cf.nu == 0  &&  cf.a == cf.b     occupation[cf.a] = cf.T     end
    end

    pairs = Tuple{Subshell,Subshell}[]
    for  i = 1:length(subshells),  j = i+1:length(subshells)
        sa = subshells[i];    sb = subshells[j]
        sa.kappa == sb.kappa                             ||  continue
        (sa in frozenSubshells  ||  sb in frozenSubshells)  &&  continue
        (haskey(occupation, sa)  &&  haskey(occupation, sb))  ||  continue
        abs( occupation[sa] - occupation[sb] ) > 1.0e-8   ||  continue
        push!(pairs, (sa, sb))
    end

    return( pairs )
end


"""
`SelfConsistent.rotateOrbitalPair(orbitals::Dict{Subshell, Orbital}, sa::Subshell, sb::Subshell, theta::Float64)`
    ... returns a copy of `orbitals` in which the pair `(sa, sb)` has been rotated into each other by `theta`.
        Orthonormality is preserved EXACTLY rather than approximately:  the Gram matrix of two orthonormal
        orbitals is the identity, and a rotation leaves it so, which is what makes every trial angle feasible
        and the search safe.  A `Dict{Subshell, Orbital}` is returned.
"""
function rotateOrbitalPair(orbitals::Dict{Subshell, Orbital}, sa::Subshell, sb::Subshell, theta::Float64)
    a = orbitals[sa];    b = orbitals[sb];    c = cos(theta);    sn = sin(theta)
    n = min(length(a.P), length(a.Q), length(b.P), length(b.Q))
    newA = Radial.Orbital(a.subshell, a.isBound, a.useStandardGrid, a.energy,
                          c*a.P[1:n]      .+ sn*b.P[1:n],      c*a.Q[1:n]      .+ sn*b.Q[1:n],
                          c*a.Pprime[1:n] .+ sn*b.Pprime[1:n], c*a.Qprime[1:n] .+ sn*b.Qprime[1:n], a.grid)
    newB = Radial.Orbital(b.subshell, b.isBound, b.useStandardGrid, b.energy,
                          -sn*a.P[1:n]      .+ c*b.P[1:n],      -sn*a.Q[1:n]      .+ c*b.Q[1:n],
                          -sn*a.Pprime[1:n] .+ c*b.Pprime[1:n], -sn*a.Qprime[1:n] .+ c*b.Qprime[1:n], b.grid)
    out = Dict{Subshell,Orbital}()
    for  (sh, orb)  in  orbitals     out[sh] = orb     end
    out[sa] = newA;    out[sb] = newB

    return( out )
end


"""
`SelfConsistent.optimizeSameKappaRotation(orbitals::Dict{Subshell, Orbital}, basis::Basis,
                                          coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
                                          grid::Radial.Grid, nucPot::Radial.Potential,
                                          frozenSubshells::Array{Subshell,1}; printout::Bool=false)`
    ... rotates each exposed same-kappa pair of the CONVERGED orbitals to the angle that minimizes the
        average-level functional, so that the result is STATIONARY with respect to that mixing and not merely
        orthonormal.  A `Dict{Subshell, Orbital}` is returned, unchanged where no pair is exposed -- so a closed
        shell passes through untouched.

        **IT RUNS ONCE, AFTER THE SCF HAS CONVERGED, AND THAT PLACEMENT IS THE WHOLE POINT.**  Applying it inside
        the iteration accomplishes nothing at all, which was measured before this was written the right way round:
        each iteration re-derives every orbital by diagonalizing its own Fock matrix, so the angle is simply reset
        on the next pass and the functional comes back to where it was.  The angle is not a quantity this
        iteration carries -- which is also WHY its fixed point is not a minimum.  Imposing orthogonality through
        the eigenproblem fixes the metric of a same-kappa pair and leaves the mixing angle wherever it falls;
        making the angle stationary inside the iteration needs the off-diagonal Lagrange multipliers that GRASP
        carries in `lagcon.f90`, which is a far larger piece of work.  This is the variational remedy that is
        available without it:  the returned orbitals are orthonormal and LOWER in the functional, which is what
        the average-level field exists to minimize, and they are not eigenvectors of a Fock matrix -- a property
        the stationarity conditions never required of them.

        MEASURED on Li-like Z = 26 over 1s^2 2s + 1s^2 2p, where 1s carries 2.0 and 2s carries 1/3:  the
        functional falls by 1.75e-05 Ha at theta = 2.1e-04, and at that same angle the E1 velocity-gauge rate
        moves from 1.2e-02 away from an independent code to 5e-04 while the length-gauge rate hardly moves.

        Up to six safeguarded Newton steps per pair, the derivatives from three evaluations at `-h, 0, +h`.  A
        step is taken only where the curvature is positive, is capped at 0.1 rad so a flat direction cannot throw
        the orbitals across the space, and is kept only where it actually lowers the functional;  the search
        stops as soon as a step fails to.
"""
function optimizeSameKappaRotation(orbitals::Dict{Subshell, Orbital}, basis::Basis,
                                   coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
                                   grid::Radial.Grid, nucPot::Radial.Potential,
                                   frozenSubshells::Array{Subshell,1}; printout::Bool=false)
    pairs = SelfConsistent.sameKappaRotationPairs(coeffs1p, basis.subshells, frozenSubshells)
    length(pairs) == 0  &&  return( orbitals )

    newOrbitals = orbitals
    for  (sa, sb)  in  pairs
        total = 0.
        for  step = 1:6
            h  = 1.0e-3
            f0 = SelfConsistent.computeFunctional(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
            fp = SelfConsistent.computeFunctional(coeffs1p, coeffs2p,
                                SelfConsistent.rotateOrbitalPair(newOrbitals, sa, sb,  h), grid, nucPot)
            fm = SelfConsistent.computeFunctional(coeffs1p, coeffs2p,
                                SelfConsistent.rotateOrbitalPair(newOrbitals, sa, sb, -h), grid, nucPot)
            #  NOTE the explicit 2.0*f0:  `2f0` is a Float32 LITERAL in Julia, not two times f0, and writing
            #  it that way made d2 come out negative and huge, so the guard below rejected every step and the
            #  rotation silently never happened.
            d1 = (fp - fm) / (2h);    d2 = (fp - 2.0*f0 + fm) / h^2
            d2 > 0.               ||  break
            theta = -d1 / d2
            abs(theta) > 0.1      &&  break
            abs(theta) < 1.0e-12  &&  break
            trial = SelfConsistent.rotateOrbitalPair(newOrbitals, sa, sb, theta)
            SelfConsistent.computeFunctional(coeffs1p, coeffs2p, trial, grid, nucPot) < f0  ||  break
            newOrbitals = trial;    total = total + theta
        end
        if  printout  &&  abs(total) > 0.
            println(">> [AL] same-kappa rotation $(string(sa))/$(string(sb)) by " *
                    @sprintf("%.3e", total) * " rad to make the functional stationary.")
        end
    end

    return( newOrbitals )
end


"""
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                       settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=GBL_AL_ANDERSON_DEPTH)`
    ... solves the average-level (AL) self-consistent field, validated to 5+ significant figures against
        literature for He/Be/Ne/Ar (see project memory), with a bVector-native architecture modeled directly
        on DBSR_HF (Zatsarinny & Froese Fischer, CPC 202, 287 (2016)): B-spline expansion coefficients
        (Dict{Subshell,Vector{Float64}}) are the SOLE, canonical, persistent per-iteration orbital state --
        no Dict{Subshell,Orbital} is maintained across iterations at all, and there is no
        Bsplines.fitVectorToPrimitives "fit back" step, since bVectors are never derived from a
        tabulated form to begin with; they come directly from Basics.diagonalize's eigenvector for each
        subshell, in turn. Orthogonality between same-kappa subshells (e.g. Ne's 1s/2s) is enforced by
        projecting the Fock matrix (Hamiltonian.projectHamiltonian, reused unchanged) against each
        ALREADY-PROCESSED lower same-kappa subshell's bVector directly inside the generalized eigenvalue
        problem, before diagonalizing -- matching DBSR_HF's hf_solve_HF.f90/hf_eiv sequential approach --
        rather than a post-hoc Löwdin symmetric orthogonalization of the whole
        same-kappa group after the per-orbital loop. The target eigenvalue index is shifted down by one for
        each such projection applied, mirroring hf_eiv's `mm = m + (orthogonalized-count) - 1`. A tabulated
        Orbital is only ever built as a disposable, read-only byproduct: once per unique "partner" subshell
        inside computeTwoElectronV (for the two-electron potential contraction), and once per
        subshell per iteration purely for reporting/energy-functional evaluation (reusing
        computeFunctional unchanged) -- never stored as competing state, never refit. A final,
        single export pass (Bsplines.generateOrbitalFromVector) produces a standard
        Dict{Subshell,Orbital} for the returned basis::Basis, so every downstream consumer (properties,
        processes, CI/DCB Hamiltonian construction) is unaffected by this being a bVector-native SCF.
        Reached via performSCF's scField = Basics.ALField() dispatch. A (new) basis::Basis is returned.

        WHAT accuracyScf MEANS HERE, and what it does not.  The iteration stops when the overlap defect
        `1 - |<b_old|b_new>|` of the worst subshell falls below settings.accuracyScf.  That defect is QUADRATIC
        in the orbital change: the change itself is `|| b_new - b_old ||_B = sqrt(2 * defect)`, so the default
        accuracyScf = 1e-6 accepts orbitals that are still moving by 1.4e-3, and leaves the argon energy 6.2e-5
        Ha short of the converged value.  A user quoting an AL energy at the default is therefore quoting ~1e-4
        and not ~1e-6.  Both numbers are printed each iteration, and the driver states at the end whether the
        field CONVERGED or merely STOPPED at maxIterationsScf -- the latter used to end in silence, so that a
        field which never converged printed exactly what a converged one printed.  Nothing about the criterion
        itself was changed, since that would move every number JAC has ever produced with it.
"""
function solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=GBL_AL_ANDERSON_DEPTH)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    # Anderson history over the CONCATENATED b-vectors. One full Gauss-Seidel sweep is the fixed-point map
    # g(x); the 0.5 damping already inside it stays, so depth 0 reproduces the previous behaviour exactly.
    xHistAL = Vector{Vector{Float64}}();    fHistAL = Vector{Vector{Float64}}()

    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once
    if  printout    println(">> [AL] (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    meanOcc = Basics.extractMeanOccupation(basis)

    # bVectors is the SOLE canonical orbital state from here on; initialized once from the starting
    # (hydrogenic) guess, then updated ONLY from diagonalization, never re-fit from a tabulated form.
    bVectors = Dict{Subshell, Vector{Float64}}()
    energies = Dict{Subshell, Float64}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitives(basis.orbitals[sh], primitives, matrixB)
        energies[sh] = basis.orbitals[sh].energy
    end

    # (2) Generate angular coefficients (unaffected by the kink fix / bVector-native rebuild)
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(Basics.ALField(), basis)

    # (3) Precompute kink-aware Slater-moment tensor caches for every rank that occurs; only the exchange
    # branches of computeTwoElectronV use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [AL] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    orbitals = Dict{Subshell, Orbital}()    # only ever a disposable, per-iteration reporting byproduct

    # FROZEN SUBSHELLS.  settings.frozenSubshells names orbitals that must NOT be varied, and until 31-Aug-2026
    # this driver ignored the field entirely: only the DFS mean-field driver honoured it, so an AL field refined
    # every subshell whatever the caller asked for.  That is not a cosmetic difference.  performSCF runs an AL
    # pass as the STARTING POINT of every EOL computation (see the :optimizedLevel branch), and the EOL solver
    # then pins its frozen orbitals to the values it finds ON ENTRY -- i.e. to whatever this driver had just done
    # to them.  A RAS step that freezes its spectroscopic shells, or a caller that hands in a converged set
    # through StartFromPrevious, therefore had those orbitals silently re-optimized anyway, and re-optimized on
    # the FULL layer basis: Basics.extractMeanOccupation averages flatly over every CSF present, so for
    # Cf^17+ with an SD layer into {7s,7p} the mean configuration that defined the "spectroscopic" 5f and 6p
    # carried the layer's 7s/7p occupation as well.  Frozen now means frozen at the value handed in.  The
    # orbitals still enter the Fock matrix and the same-kappa orthogonality projection of every shell refined
    # after them, exactly as an occupied orbital must; they are simply never themselves diagonalized.
    frozenSubshells = [ sh  for sh in basis.subshells  if  sh in settings.frozenSubshells ]
    pinnedB = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in frozenSubshells )
    pinnedE = Dict{Subshell, Float64}(         sh => energies[sh]        for sh in frozenSubshells )
    if  printout  &&  !isempty(frozenSubshells)
        println(">> [AL] frozen and NOT refined: " * join(string.(frozenSubshells), ", ") *
                ";  $(length(basis.subshells) - length(frozenSubshells)) of $(length(basis.subshells)) subshells are varied.")
    end

    isConverged = false;    NoIterations = 0;    lastDefect = 1.0;    lastStep = 1.0;    lastShell = basis.subshells[1]

    for  iter = 1:Basics.maxIterations(settings.scfRoute)
        println("\n> SCF interation $(iter) [AL]: ")
        newBVectors = Dict{Subshell, Vector{Float64}}();    newEnergies = Dict{Subshell, Float64}()
        processedBVectors = Dict{Subshell, Vector{Float64}}()
        dpm = Dict{Subshell, Float64}()

        # The subshell-independent part of every two-electron matrix is reused across this sweep and ONLY
        # across this sweep: bVectors is reassigned at the END of the iteration, so the partner orbitals are
        # fixed here but not between iterations, and a cache carried over would serve stale matrices.  Built
        # fresh each iteration for exactly that reason.  Measured redundancy: 3.5x for argon, 9.1x for Th+.
        directKernels   = Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}()
        exchangeKernels = Dict{Tuple{Int64,Subshell},Array{Float64,2}}()

        for  subshell  in  basis.subshells
            if  subshell in frozenSubshells
                # Kept at the pinned value, but still registered so that the shells refined after it are
                # projected orthogonal to it.
                newBVectors[subshell]       = copy(pinnedB[subshell])
                newEnergies[subshell]       = pinnedE[subshell]
                processedBVectors[subshell] = copy(pinnedB[subshell])
                dpm[subshell]               = 0.0
                continue
            end
            occ = meanOcc[subshell]
            # A SUBSHELL OF THE BASIS CAN CARRY ZERO MEAN OCCUPATION, and it was long assumed here that it
            # could not -- `Basics.extractMeanOccupation` averages over the CSFs OF THIS BASIS, not over every
            # CSF the configurations could produce, so a basis restricted to ONE SYMMETRY empties any subshell
            # that no CSF of that symmetry can occupy.  The clean example is Sc+ [Ar] 3d 4s at J = 1: coupling
            # 3d_5/2 with 4s_1/2 gives J = 2 or 3 only, so no J = 1 CSF occupies 3d_5/2 at all and its mean
            # occupation is exactly 0.  A correlation layer does the same to any virtual its CSFs leave empty.
            # Refining such a subshell asked `computeFockMatrix` for a mean field that does not exist, and it
            # refused (correctly) -- which aborted the whole SCF, and with it every layer of a RAS ladder that
            # had already succeeded.  There is nothing to refine: no mean field is defined, so the orbital is
            # carried forward unchanged, exactly as a frozen one is, and still registered so that the shells
            # refined after it are projected orthogonal to it.
            if  abs(occ) < 1.0e-12
                if  printout
                    println(">> $subshell carries zero mean occupation in this basis; no mean field is defined " *
                            "for it, so it is kept unchanged.")
                end
                newBVectors[subshell]       = copy(bVectors[subshell])
                newEnergies[subshell]       = energies[subshell]
                processedBVectors[subshell] = copy(bVectors[subshell])
                dpm[subshell]               = 0.0
                continue
            end
            print(">> Refine $subshell orbital with mean occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrix(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches;
                                                              directKernels=directKernels,
                                                              exchangeKernels=exchangeKernels)

            # Orthogonality: project against every ALREADY-PROCESSED lower same-kappa subshell's bVector,
            # directly inside the eigenvalue problem (DBSR_HF hf_eiv style), not post-hoc.
            count = Base.count( sh2 -> sh2.kappa == subshell.kappa  &&  sh2.shape == subshell.shape, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Bsplines.diagonalizeLocalMatrix(subshell.kappa, matrix, matrixB, primitives)
            l  = Basics.subshell_l(subshell)
            mm = Bsplines.findPositiveBranchStart(wc.values)
            ni = mm + subshell.n - l - count - 1
            # WHAT THIS INDEX MEANS, and why the offset is the interesting part.  The Hamiltonian has already
            # been projected against every processed same-kappa orbital, so wc.vectors[mm] is the LOWEST state
            # orthogonal to them, and  n - l - count - 1  is therefore the number of radial nodes this orbital
            # is being given ABOVE the minimum that orthogonality alone forces.  For a contiguous sequence it
            # is exactly zero -- 3s has n-l-count-1 = 3-0-2-1 = 0 -- so a spectroscopic computation always sits
            # at the minimum and nothing here concerns it.  The offset is nonzero only where the kappa sequence
            # has a GAP, i.e. only for a correlation orbital named above its first free slot: with 3d, 4d and 5d
            # occupied, a "7d" is given 4 nodes where orthogonality demands 3.
            #   THAT OFFSET IS A CONSTRAINT NOBODY ASKED FOR.  A correlation orbital is a basis function, not an
            # observable; the functional sees it only through the space the CSFs span, and the only genuine
            # requirements are square-integrability and orthonormality within each kappa.  Orthogonality to k
            # occupied same-kappa orbitals already forces at least k sign changes (a function orthogonal to the
            # first k eigenfunctions of a Sturm-Liouville operator cannot have fewer), so node counting is not an
            # extra imposition -- but the MINIMUM is what the constraint gives, and anything beyond it is the
            # label speaking rather than the physics.  GRASP varies under orthogonality alone and so lands on the
            # minimum whatever the orbital is called, which is exactly why its block naming (7s ... 7g, with 5g
            # and 6d absent) is sound there and does not transfer here.
            #   THE FIRST AVAILABLE ORTHOGONAL STATE IS THEREFORE WHAT IS TAKEN, making the label pure
            # bookkeeping as it is in GRASP.  This was measured before it was made the default (1-Sep-2026,
            # probe-firstfree.jl on a three-valence-electron ion at Z = 98).  Three AL runs on one space:
            #   7g as this routine used to build it :  2 nodes, <r> = 3.0579 a0
            #   7g as first available               :  0 nodes, <r> = 1.1721 a0
            #   5g named explicitly                 :  0 nodes, <r> = 1.1721 a0
            #   |<7g first-available | 5g explicit>| = 1.000000     -- the SAME orbital, to six digits
            #   |<7g as-built-before | 5g explicit>| = 0.000047     -- essentially ORTHOGONAL to it
            # So the old index did not produce a poorer version of the right orbital, it produced a different
            # basis function: with the correlated 5f shell at <r> = 0.85 a0, a g orbital at 3.06 a0 with two
            # radial nodes through the valence region cannot correlate it.  The cost is not marginal.  On the
            # same CI space (3 367 CSFs, references + S + D), a layer named {7s,7p,7d,7f,7g} was worth 387
            # cm^-1 of correlation energy on the ground level and one named {7s,7p,6d,6f,5g} was worth 5 409 --
            # a factor of fourteen, from nothing but the label -- and the first-free layer was better on every
            # observable, moving the 5f fine structure, the 6p_3/2 level and the clock transition all closer to
            # the published values.
            #   NOTHING SPECTROSCOPIC CHANGES.  The offset removed is n - l - count - 1, which is identically
            # zero for a contiguous kappa sequence (3s gives 3-0-2-1 = 0), so it can only ever act on a
            # correlation slot named above its first free label.  JAC_AL_NODECOUNT restores the old selection
            # for anyone who needs to reproduce an earlier number.
            if  ni > mm  &&  !haskey(ENV, "JAC_AL_NODECOUNT")
                if  printout   print("[first-free: node offset $(ni-mm) dropped] ")   end
                ni = mm
            end
            rawVector  = wc.vectors[ni];    newEnergy = wc.values[ni]

            # Damping (26-Jul-2026): sequential (Jacobi-style) per-orbital refinement combined with the
            # in-matrix orthogonality projection above reproduces the period-2 SCF oscillation already
            # documented for this projector when applied one-subshell-at-a-time (see memory
            # project_df_al_kink_bug.md) -- the projection formula implicitly assumes same-kappa orbitals
            # are varied simultaneously. Standard fix: linear mixing of the new and previous bVector before
            # acceptance (a common SCF damping technique), aligning sign in the B-metric first since a
            # generalized eigensolver may return either sign for an eigenvector.
            oldVector = bVectors[subshell]
            if  transpose(oldVector) * matrixB * rawVector < 0    rawVector = -rawVector    end
            damping = 0.5
            mixed     = damping * oldVector + (1.0 - damping) * rawVector
            newVector = mixed / sqrt( transpose(mixed) * matrixB * mixed )

            newBVectors[subshell]       = newVector
            newEnergies[subshell]       = newEnergy
            processedBVectors[subshell] = newVector

            oldVector = bVectors[subshell]
            ovlap     = abs( transpose(oldVector) * matrixB * newVector )
            dpm[subshell] = 1.0 - ovlap
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end

        # Disposable, read-only tabulation for reporting/energy purposes only -- never refit, never stored
        # as persistent state; rebuilt fresh from newBVectors/newEnergies every iteration.
        newOrbitals = Dict{Subshell, Orbital}()
        for  sh  in  basis.subshells
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, newEnergies[sh], newBVectors[sh], primitives)
        end
        # ANDERSON ACCELERATION on the orbitals themselves.  The mean-field driver already accelerates its
        # screening potential this way, with 1.6-2.6x fewer iterations and the same solution to ~1e-7; here
        # the iterate is the concatenated b-vector set and g(x) is one full sweep.  The sweep sign-aligns
        # each eigenvector against the previous one, so g is sign-consistent with x and the residual means
        # what it should.  depth <= 0 leaves the plain damped iteration untouched.
        if  andersonDepth > 0
            xNow = vcat( [ bVectors[sh]     for sh in basis.subshells ]... )
            gNow = vcat( [ newBVectors[sh]  for sh in basis.subshells ]... )
            fNow = gNow - xNow
            push!(xHistAL, copy(xNow));    push!(fHistAL, fNow)
            if  length(xHistAL) > andersonDepth + 1    popfirst!(xHistAL);   popfirst!(fHistAL)    end
            m = length(xHistAL) - 1
            if  m >= 1
                dF = zeros( length(fNow), m );    dX = zeros( length(fNow), m )
                for  j = 1:m
                    dF[:,j] = fHistAL[j+1] - fHistAL[j];    dX[:,j] = xHistAL[j+1] - xHistAL[j]
                end
                local gamma
                try
                    gamma = dF \ fNow
                catch
                    gamma = zeros(m)              # a rank-deficient window: fall back to the plain step
                end
                if  any(!isfinite, gamma)   gamma = zeros(m)   end
                xNew = xNow + fNow - (dX + dF) * gamma
                if  all(isfinite, xNew)
                    i0 = 0
                    for  sh  in  basis.subshells
                        v = xNew[i0+1 : i0+nsL+nsS];    i0 = i0 + nsL + nsS
                        nrm = sqrt( abs(transpose(v) * matrixB * v) )
                        if  nrm > 1.0e-12    newBVectors[sh] = v / nrm    end
                    end
                    for  sh  in  basis.subshells
                        newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, newEnergies[sh],
                                                                            newBVectors[sh], primitives)
                    end
                end
            end
        end
        # Under test: restore the same-kappa orthonormality that the damping step destroys.
        if  SelfConsistent.GBL_SCF_REORTHONORMALIZE
            (newOrbitals, newBVectors) = SelfConsistent.orthonormalizeSameKappa(newOrbitals, newBVectors,
                                                                basis.subshells, primitives, matrixB)
        end
        # Anderson acceleration works on the CONCATENATED b-vector of every subshell, and
        # orthonormalizeSameKappa rotates within each kappa block; both would move a frozen orbital.  Put the
        # pinned values back before anything is measured or accepted.
        for  sh  in  frozenSubshells
            newBVectors[sh] = copy(pinnedB[sh]);    newEnergies[sh] = pinnedE[sh]
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, pinnedE[sh], pinnedB[sh], primitives)
        end
        # Convergence is measured against what is ACTUALLY accepted, which Anderson may have moved. Only
        # when it is active, so that depth 0 remains a bit-for-bit control on the previous behaviour.
        if  andersonDepth > 0
            for  sh  in  basis.subshells
                dpm[sh] = 1.0 - abs( transpose(bVectors[sh]) * matrixB * newBVectors[sh] )
            end
        end
        eFunctional = SelfConsistent.computeFunctional(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        # The overlap defect 1 - |<b_old|b_new>| is QUADRATIC in the orbital change; the change itself is
        # || b_new - b_old ||_B = sqrt(2 * defect), and that is what a user means by "the orbitals still move
        # by".  Both are reported so that the tolerance cannot be misread by a factor of its own square root.
        overlapDefect = 0.;    worstShell = basis.subshells[1]
        for  (sh, d)  in dpm    if  d > overlapDefect   overlapDefect = d;   worstShell = sh   end    end
        orbitalConv = overlapDefect < 1.0 ? 1.0 - overlapDefect : 0.0
        orbitalStep = sqrt( 2.0 * max(0., overlapDefect) )
        lastStep    = orbitalStep;    lastDefect = overlapDefect;    lastShell = worstShell;    NoIterations = iter

        println(">> Total energy = $(eFunctional*1)   orbital-conv = $orbitalConv   orbital-acc = $(1.0 - orbitalConv)" *
                "   orbital-step = " * @sprintf("%.3e", orbitalStep) * " ($worstShell)")

        bVectors = newBVectors;    energies = newEnergies;    orbitals = newOrbitals
        if  overlapDefect < settings.accuracyScf    isConverged = true;    break    end
    end

    # Say which of the two happened.  An iteration that merely runs out of steps used to end in silence, so
    # that a stopped field and a converged one printed the same thing and were quoted the same way.
    if      isConverged
        println(">> [AL] converged after $NoIterations iterations: overlap defect " *
                @sprintf("%.2e", lastDefect) * " < accuracyScf = " * @sprintf("%.2e", settings.accuracyScf) *
                ", with the orbitals still moving by " * @sprintf("%.2e", lastStep) * " ($lastShell).")
    else
        println(">> [AL] STOPPED, NOT CONVERGED, after $NoIterations iterations (maxIterationsScf): overlap defect " *
                @sprintf("%.2e", lastDefect) * " has not reached accuracyScf = " *
                @sprintf("%.2e", settings.accuracyScf) * "; the orbitals are still moving by " *
                @sprintf("%.2e", lastStep) * " ($lastShell).  The energies below are NOT self-consistent.")
        # COLLECTED AS WELL, not merely printed.  solveMeanFieldBasis and solveAverageAtomField have
        # raised a warning here since long ago; AL and EOL did not, so the two fields recommended for
        # ordinary work -- AL for a reference layer, EOL for an added one -- were the only ones whose
        # non-convergence reached nothing but the console.  A driver that writes its multiplet to a
        # file then records unconverged energies with no trace that they are unconverged.
        Defaults.warn(AddWarning(), "SelfConsistent.solveAverageLevelField(): the AL field did NOT converge for " *
                      string(basis.subshells) * " -- overlap defect " * @sprintf("%.1e", lastDefect) *
                      " after $NoIterations iterations (maxIterationsScf).  The energies are NOT self-consistent.")
    end

    # THE SAME-KAPPA ROTATION, once, on the CONVERGED orbitals.  OFF by default;  see the note at
    # SelfConsistent.GBL_AL_SAME_KAPPA_ROTATION for the measurement and for why switching it on is an editorial
    # act, and the docstring of optimizeSameKappaRotation for why it cannot be done inside the iteration.
    if  SelfConsistent.GBL_AL_SAME_KAPPA_ROTATION
        orbitals = SelfConsistent.optimizeSameKappaRotation(orbitals, basis, coeffs1p, coeffs2p, grid, nucPot,
                                                            frozenSubshells; printout=printout)
    end

    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end
