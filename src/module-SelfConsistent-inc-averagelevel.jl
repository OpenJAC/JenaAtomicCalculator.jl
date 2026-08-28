

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
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                       settings::AsfSettings; printout::Bool=true)`
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

    isConverged = false;    NoIterations = 0;    lastDefect = 1.0;    lastStep = 1.0;    lastShell = basis.subshells[1]

    for  iter = 1:settings.maxIterationsScf
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
            occ = meanOcc[subshell]
            print(">> Refine $subshell orbital with mean occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrix(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches;
                                                              directKernels=directKernels,
                                                              exchangeKernels=exchangeKernels)

            # Orthogonality: project against every ALREADY-PROCESSED lower same-kappa subshell's bVector,
            # directly inside the eigenvalue problem (DBSR_HF hf_eiv style), not post-hoc.
            count = Base.count( sh2 -> sh2.kappa == subshell.kappa, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Bsplines.diagonalizeLocalMatrix(subshell.kappa, matrix, matrixB, primitives)
            l  = Basics.subshell_l(subshell)
            mm = Bsplines.findPositiveBranchStart(wc.values)
            ni = mm + subshell.n - l - count - 1
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
                    gamma = zeros(m)              ## a rank-deficient window: fall back to the plain step
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
    end

    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end
