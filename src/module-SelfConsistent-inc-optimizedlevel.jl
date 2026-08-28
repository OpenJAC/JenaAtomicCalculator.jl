

# The EOL (optimized-level) path: the CSF-pair coefficient cache, the CI matrix and its diagonalization,
# the generalized occupation, and the two EOL solvers.

"""
`SelfConsistent.cacheCsfPairCoefficientsEOL(sym::LevelSymmetry, basis::Basis)`
    ... caches, for every CSF pair (r,s) with symmetry sym in the given basis, the pure spin-angular
        coefficients (independent of the orbitals/radial functions) as returned by
        SpinAngular.computeCoefficientsScalar -- the same call Hamiltonian.setupMatrix/setupMatrixKinkAware
        make internally to build the CI Hamiltonian matrix. Here the intermediate coefficient lists are
        retained instead of being discarded, so they can be reused across the whole EOL outer SCF+CI loop
        (they depend only on the fixed CSF list, never on the current orbitals). A tuple
        (idxCsf::Array{Int64,1}, cache1p, cache2p) is returned, with (r,s) keys running over the LOCAL
        index (1:length(idxCsf)) into the symmetry block.
"""
function cacheCsfPairCoefficientsEOL(sym::LevelSymmetry, basis::Basis)
    idxCsf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == sym.J   &&   basis.csfs[idx].parity == sym.parity    push!(idxCsf, idx)    end
    end
    n = length(idxCsf)

    cache1p = Dict{Tuple{Int64,Int64}, Array{Coefficient1p,1}}()
    cache2p = Dict{Tuple{Int64,Int64}, Array{Coefficient2p,1}}()
    for  r = 1:n
        for  s = 1:n
            csfR = basis.csfs[idxCsf[r]];   csfS = basis.csfs[idxCsf[s]]
            cache1p[(r,s)] = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                                    csfR, csfS, basis.subshells)
            cache2p[(r,s)] = SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus),
                                                                    csfR, csfS, basis.subshells)
        end
    end

    return( (idxCsf, cache1p, cache2p) )
end


"""
`SelfConsistent.buildCIMatrixEOL(idxCsf::Array{Int64,1}, cache1p, cache2p, orbitals::Dict{Subshell, Orbital},
                                 grid::Radial.Grid, potential::Radial.Potential)`
    ... (re-) builds the CI Hamiltonian matrix for one symmetry block from CACHED, orbital-independent
        angular coefficients (see cacheCsfPairCoefficientsEOL) and the CURRENT radial functions in
        orbitals -- algebraically identical to Hamiltonian.setupMatrixKinkAware's pure-Coulomb contribution
        (kink-aware InteractionStrength.XL_CoulombKinkAware, matching the rest of the average-level (ALField) line; Breit
        and QED are added only once, at the final Hamiltonian.performCIKinkAware call, exactly as for the AL
        scheme), but without repeating the (unchanged) angular-coefficient computation on every outer
        SCF+CI iteration. The trailing radial1pCache/radial2pCache arguments memoize each radial integral by
        its bare subshell labels (never by CSF-pair index), so a radial integral shared by many different
        CSF pairs -- e.g. the same "1s-1s" self-interaction appearing in every CSF's diagonal term -- is
        evaluated once per outer iteration rather than once per (r,s) occurrence; pass the SAME two Dicts
        into every block diagonalized within one outer iteration to also share the cache across blocks
        (found by profiling to be the dominant EOL cost for multi-CSF cases -- see
        project_eol_implementation.md). A  matrix::Array{Float64,2}  is returned.
"""
function buildCIMatrixEOL(idxCsf::Array{Int64,1}, cache1p, cache2p, orbitals::Dict{Subshell, Orbital},
                          grid::Radial.Grid, potential::Radial.Potential,
                          radial1pCache::Dict{Tuple{Subshell,Subshell},Float64}          = Dict{Tuple{Subshell,Subshell},Float64}(),
                          radial2pCache::Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64} =
                                        Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}())
    # radial1pCache/radial2pCache: keyed purely by subshell labels (never by CSF-pair index), so a radial
    # integral shared by MANY different CSF pairs -- e.g. the same "1s-1s" self-interaction appearing in
    # every CSF's diagonal term -- is evaluated once per outer SCF+CI iteration and reused, instead of once
    # per (r,s) occurrence. Pass the SAME cache Dicts in across every block diagonalized within one outer
    # iteration (they also depend only on the current orbitals, not on which block/CSF-pair references them);
    # a fresh empty cache per call (the default) is still correct, just without the cross-block reuse.
    # Hermitian-symmetry shortcut (28-Jul-2026): only the UPPER triangle (r<=s) is computed -- exact, not
    # just safe, since this matrix feeds diagonalizeBlockEOL -> Basics.diagonalize(MatrixWithLinearAlgebra(),
    # ...), whose Symmetric(matrix) wrapper (default uplo=:U) already discards the lower triangle. See
    # Hamiltonian.setupMatrix's identical note for the confirming test.
    n = length(idxCsf);   matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = r:n
            me = 0.
            for  cf in cache1p[(r,s)]
                I_ab = get!(radial1pCache, (cf.a,cf.b)) do
                    RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
                end
                me = me + cf.T * I_ab
            end
            for  cf in cache2p[(r,s)]
                R_abcd = get!(radial2pCache, (cf.nu,cf.a,cf.b,cf.c,cf.d)) do
                    InteractionStrength.XL_CoulombKinkAware(cf.nu, orbitals[cf.a], orbitals[cf.b], orbitals[cf.c], orbitals[cf.d], grid)
                end
                me = me + cf.V * R_abcd
            end
            matrix[r,s] = me
        end
    end

    return( matrix )
end


"""
`SelfConsistent.diagonalizeBlockEOL(sym::LevelSymmetry, idxCsf::Array{Int64,1}, matrix::Array{Float64,2}, basis::Basis)`
    ... diagonalizes the (already-built) CI matrix of one symmetry block and reassigns the eigenvectors to
        Level instances w.r.t. the full basis (zero-padded outside this block), exactly as
        Hamiltonian.performCI/performCIKinkAware do internally per symmetry block -- factored out here so it
        can be called every EOL outer iteration without their Multiplet-merge/tabulate overhead.
        An  Array{Level,1}  is returned.
"""
function diagonalizeBlockEOL(sym::LevelSymmetry, idxCsf::Array{Int64,1}, matrix::Array{Float64,2}, basis::Basis)
    eigen  = Basics.diagonalize(MatrixWithLinearAlgebra(), matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        vector = zeros( length(basis.csfs) )
        for  (r, idx)  in  enumerate(idxCsf)    vector[idx] = eigen.vectors[ev][r]    end
        newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector )
        push!( levels, newlevel)
    end

    return( levels )
end


"""
`SelfConsistent.selectTargetLevelsEOL(mp::Multiplet, levelSelectionCI::LevelSelection)`
    ... determines the target level(s) for the EOL functional from an (energy-sorted) multiplet mp,
        using levelSelectionCI EXCLUSIVELY: either indices or symmetries may be given, never both.
        If symmetries is given, the target set is the LOWEST level of each listed symmetry (the classic
        EOL use case, e.g. the lowest of J=1/2^+ together with the lowest of J=3/2^+). If indices is
        given, the target set is those exact levels by their (global, energy-sorted) index. If
        levelSelectionCI is inactive or both arrays are empty, the target set defaults to the single
        lowest level overall -- a genuine OL (one-level) computation.
        An  Array{Level,1}  is returned.
"""
function selectTargetLevelsEOL(mp::Multiplet, levelSelectionCI::LevelSelection)
    if  !levelSelectionCI.active  ||  ( isempty(levelSelectionCI.indices) && isempty(levelSelectionCI.symmetries) )
        return( [ mp.levels[1] ] )
    elseif  !isempty(levelSelectionCI.indices)  &&  !isempty(levelSelectionCI.symmetries)
        error("stop a; levelSelectionCI must specify EITHER indices OR symmetries for the EOL scheme, not both.")
    elseif  !isempty(levelSelectionCI.symmetries)
        targetLevels = Level[]
        for  sym  in  levelSelectionCI.symmetries
            for  level  in  mp.levels                                   # mp.levels is energy-sorted; first match = lowest
                if  LevelSymmetry(level.J, level.parity) == sym    push!(targetLevels, level);   break    end
            end
        end
        return( targetLevels )
    else
        return( [ mp.levels[i]  for i in levelSelectionCI.indices ] )
    end
end


"""
`SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels::Array{Level,1}, basis::Basis)`
    ... computes the EOL generalized occupation number per subshell,
        q(nlj) = Σᵣ d²_r · q_r(nlj), with d²_r = Σᵢ WT_i · (c_r⁽ⁱ⁾)² the DIAGONAL (r=r) case of the same
        statistical-(2J+1)-weighted generalized weight used in combineAngularCoefficientsEOL, and q_r(nlj)
        the plain occupation of nlj in CSF r (basis.csfs[r].occupation). This REPLACES
        Basics.extractMeanOccupation(basis) -- which averages FLATLY over every CSF in the whole basis,
        appropriate only for AL's single-average-CSF philosophy -- with the occupation actually implied by
        the current target level(s)' own CI mixing, recomputed every outer iteration since the mixing
        coefficients change. Without this, computeFockMatrix's (1.0/occ) two-electron scaling uses a
        basis-wide average that can be wildly wrong for a multi-configuration target level (e.g. an
        essentially-pure 2s² level in a 2s²/2p² basis would otherwise see occ(2s) diluted by the unrelated
        2p² CSFs, inflating its two-electron potential many-fold). A  Dict{Subshell,Float64}  is returned.
"""
function computeGeneralizedOccupationEOL(blockCaches, targetLevels::Array{Level,1}, basis::Basis)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumWeights  = sum( twiceJp1(level.J)  for level in targetLevels )
    weights     = [ twiceJp1(level.J) / sumWeights  for level in targetLevels ]

    occs = Dict{Subshell, Float64}();   for  sh in basis.subshells   occs[sh] = 0.   end
    for  (_, (idxCsf, _, _))  in  blockCaches
        for  r  in  idxCsf
            drr = 0.
            for  (i, level)  in  enumerate(targetLevels)    drr = drr + weights[i] * level.mc[r]^2    end
            if  drr == 0.    continue    end
            for  (is, sh)  in  enumerate(basis.subshells)   occs[sh] = occs[sh] + drr * basis.csfs[r].occupation[is]   end
        end
    end

    return( occs )
end


"""
`SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels::Array{Level,1})`
    ... generalizes SelfConsistent.computeAngularCoefficients (AL's single-CSF-average analog: loop CSFs,
        weight 1/ncsf) to CSF PAIRS, weighted by the EOL generalized weight
        d²_rs = Σᵢ WT_i · c_r⁽ⁱ⁾ · c_s⁽ⁱ⁾, with WT_i the normalized statistical (2Jᵢ+1) weight of each
        target level i in targetLevels (its own mixing vector c⁽ⁱ⁾ = level.mc) -- never a user-supplied
        weight. Reuses the cached, orbital-independent per-CSF-pair coefficients from
        cacheCsfPairCoefficientsEOL for every relevant symmetry block in blockCaches; a level's mc vector
        is zero outside its own block, so multiple symmetry blocks combine correctly without
        special-casing. The dedup/condensation logic is identical to computeAngularCoefficients.
        A Tuple  (coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1})  is returned.
"""
function combineAngularCoefficientsEOL(blockCaches, targetLevels::Array{Level,1}; pairs::Symbol=:all)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumWeights  = sum( twiceJp1(level.J)  for level in targetLevels )
    weights     = [ twiceJp1(level.J) / sumWeights  for level in targetLevels ]

    coeffs1p = Coefficient1p[];   coeffs2p = Coefficient2p[]
    for  (_, (idxCsf, cache1p, cache2p))  in  blockCaches
        n = length(idxCsf)
        for  r = 1:n
            for  s = 1:n
                # pairs = :all (default, unchanged) | :diagonal (r == s only) | :offdiagonal (r != s only).
                # The split exists to test whether the off-diagonal CSF-pair contributions should be scaled
                # by the generalized occupation at all -- see the DA-term note in solveOptimizedLevelField.
                if      pairs == :diagonal      &&  r != s     continue
                elseif  pairs == :offdiagonal   &&  r == s     continue
                end
                drs = 0.
                for  (i, level)  in  enumerate(targetLevels)    drs = drs + weights[i] * level.mc[idxCsf[r]] * level.mc[idxCsf[s]]    end
                if  drs == 0.    continue    end
                for  cf in cache1p[(r,s)]   push!(coeffs1p, Coefficient1p(cf.nu, cf.a, cf.b, cf.T * drs) )   end
                for  cf in cache2p[(r,s)]   push!(coeffs2p, Coefficient2p(cf.nu, cf.a, cf.b, cf.c, cf.d, cf.V * drs) )   end
            end
        end
    end

    # Condense angular coefficients if they refer to the same set of orbitals -- identical to
    # SelfConsistent.computeAngularCoefficients' own condensation tail.
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
                end
            end
            push!(coeffs2px, Coefficient2p(nu, a, b, c, d, V) );   V = 0.
        end
    end

    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, grid::Radial.Grid,
        nucPot::Radial.Potential)`  
    ... the EOL energy as a plain scalar function of the orbital B-spline coefficient vectors, with the
        angular coefficients held fixed. This is exactly the functional solveOptimizedLevelField reports,
        just expressed in the variables that are actually varied, so that it can be differentiated. Note
        that it includes coeffs1p -- which the Fock matrix of the present scheme never receives, and which
        is the inconsistency the rotation-based path exists to remove. A value::Float64 is returned.
"""
function energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
                                   coeffs1p::Array{Coefficient1p,1},
                                   coeffs2p::Array{Coefficient2p,1},
                                   subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                   grid::Radial.Grid, nucPot::Radial.Potential)
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
    end
    return( SelfConsistent.computeFunctional(coeffs1p, coeffs2p, orbitals, grid, nucPot) )
end


"""
`SelfConsistent.positiveBranchSpectrum(subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
        nucPot::Radial.Potential, matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}})`  
    ... the positive-energy eigenvectors of the ONE-ELECTRON Dirac matrix of every kappa present, together
        with their eigenvalues. It depends only on the nuclear potential and the B-spline basis, NOT on the
        orbitals, so it is constant for a whole run and must be built once and passed in -- both
        virtualDirections and projectOntoPositiveBranch used to rebuild it on every call,
        i.e. twice per iteration. Measured on Be-like C^2+: an iteration cost 1.77 s of which the line
        search was only 0.15 s, and this was the bulk of the remainder.
        A Dict{Int64, Tuple{Array{Vector{Float64},1}, Vector{Float64}}} is returned, keyed by kappa.
"""
function positiveBranchSpectrum(subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, matrixB::Array{Float64,2},
                                       storage::Dict{String,Array{Float64,2}})
    spectrum = Dict{Int64, Tuple{Array{Vector{Float64},1}, Vector{Float64}}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        oneEl = Bsplines.setupLocalMatrix(kappa, primitives, nucPot, storage)
        wc    = Bsplines.diagonalizeLocalMatrix(kappa, oneEl, matrixB, primitives)
        mm    = Bsplines.findPositiveBranchStart(wc.values)
        spectrum[kappa] = ( [ wc.vectors[i]  for i = mm:length(wc.values) ],
                            [ wc.values[i]   for i = mm:length(wc.values) ] )
    end
    return( spectrum )
end


"""
`SelfConsistent.virtualDirections(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, nucPot::Radial.Potential, matrixB::Array{Float64,2},
        storage::Dict{String,Array{Float64,2}}; nVirtual::Int64=20)`  
    ... builds, for each subshell, an S-orthonormal set of ALLOWED rotation directions: positive-branch
        eigenvectors of the one-electron Dirac matrix of that kappa, projected free of every occupied
        orbital of the same kappa. Rotations among the occupied orbitals themselves are deliberately
        excluded -- they leave the CSF space invariant, are redundant, and would make the Hessian singular.
        A Dict{Subshell, Array{Vector{Float64},1}} is returned.
"""
function virtualDirections(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
                                  primitives::Bsplines.Primitives, nucPot::Radial.Potential,
                                  matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}};
                                  nVirtual::Int64=20, spectrum=nothing)
    posSpec = isnothing(spectrum) ?
              SelfConsistent.positiveBranchSpectrum(subshells, primitives, nucPot, matrixB, storage) : spectrum
    virtuals = Dict{Subshell, Array{Vector{Float64},1}}()
    for  sh  in  subshells
        # the occupied orbitals of this kappa, S-orthonormalized so that the projection below is exact
        occSame = Vector{Float64}[]
        for  s2  in  subshells
            if  s2.kappa != sh.kappa    continue    end
            v = copy(bVectors[s2])
            for  u in occSame    v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            if  nrm > 1.0e-10    push!(occSame, v / nrm)    end
        end
        # a one-electron reference spectrum for this kappa; orbital-independent, hence a fixed frame
        (posVecs, _) = posSpec[sh.kappa]
        dirs  = Vector{Float64}[]
        for  i = 1:length(posVecs)
            v = copy(posVecs[i])
            for  u in occSame    v = v - (transpose(u) * matrixB * v) * u    end
            for  u in dirs       v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            if  nrm > 1.0e-8    push!(dirs, v / nrm)    end
            if  length(dirs) >= nVirtual    break    end
        end
        virtuals[sh] = dirs
    end
    return( virtuals )
end


"""
`SelfConsistent.expandBVector(vec::Vector{Float64}, primitives::Bsplines.Primitives)`  
    ... expands a B-spline coefficient vector into its large and small radial components on the grid,
        (P, Q), exactly as the matrix contraction implies -- P(r) = sum_k vec[k] B_k^L(r) and likewise for Q.
        Deliberately NOT via Bsplines.generateOrbitalFromVector, which truncates at mtp and cleans
        small values: that round trip is lossy, and the whole point here is to reproduce the matrix product
        exactly. A tuple (P, Q) of Vector{Float64} over the full grid is returned.
"""
function expandBVector(vec::Vector{Float64}, primitives::Bsplines.Primitives)
    grid = primitives.grid;    nsL = grid.nsL;    nsS = grid.nsS
    P = zeros( grid.NoPoints );    Q = zeros( grid.NoPoints )
    for  k = 1:nsL
        bs = primitives.bsplinesL[k];   add = 1 - bs.lower
        for  r = bs.lower:min(bs.upper, grid.NoPoints)    P[r] = P[r] + vec[k] * bs.bs[r+add]    end
    end
    for  k = 1:nsS
        bs = primitives.bsplinesS[k];   add = 1 - bs.lower
        for  r = bs.lower:min(bs.upper, grid.NoPoints)    Q[r] = Q[r] + vec[nsL+k] * bs.bs[r+add]    end
    end
    return( (P, Q) )
end


"""
`SelfConsistent.screenedProduct(Vk::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
        primitives::Bsplines.Primitives)`  
    ... forms the vector whose i-th entry is  INT B_i(r) f(r) w_r V_k(r) dr,  with f = P on the large block
        and Q on the small one -- i.e. the product of the screened-potential matrix with a coefficient
        vector, WITHOUT ever building that matrix.

        The matrix that InteractionStrength.XL_CoulombKinkAware assembles has entries
        wm[i,k] = INT B_i B_k w_r V_k, a full double loop over roughly 110000 (i,k) pairs which re-expands
        its B-spline arrays inside the inner loop -- and it is then contracted with a single vector.
        Measured 10-Aug-2026 on Be-like C^2+: that assembly was 0.99 s of a 1.32 s gradient, and the
        gradient was ~80% of a rotation iteration.  Since the matrix is symmetric within each block, both
        products a caller needs use the same potential and differ only in the vector.
        A Vector{Float64} of length nsL+nsS is returned.
"""
function screenedProduct(Vk::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64},
                                primitives::Bsplines.Primitives)
    grid = primitives.grid;    nsL = grid.nsL;    nsS = grid.nsS
    out  = zeros( nsL + nsS )
    for  i = 1:nsL
        bs  = primitives.bsplinesL[i];    add = 1 - bs.lower
        mtp = min( bs.upper, length(Vk), length(P) );    wa = 0.
        for  r = max(2, bs.lower):mtp    wa = wa + bs.bs[r+add] * P[r] * grid.wr[r] * Vk[r]    end
        out[i] = wa
    end
    for  i = 1:nsS
        bs  = primitives.bsplinesS[i];    add = 1 - bs.lower
        mtp = min( bs.upper, length(Vk), length(Q) );    wa = 0.
        for  r = max(2, bs.lower):mtp    wa = wa + bs.bs[r+add] * Q[r] * grid.wr[r] * Vk[r]    end
        out[nsL+i] = wa
    end
    return( out )
end


"""
`SelfConsistent.computeOrbitalGradient(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, nucPot::Radial.Potential,
        storage::Dict{String,Array{Float64,2}})`  
    ... the ANALYTIC gradient dE/db_a of the same energy that energyFromBVectors evaluates, for every
        subshell, as a full B-spline coefficient vector. Nothing here is new machinery: the one-electron
        integral is I(a,b) = b_a^T H1 b_b with H1 = Bsplines.setupLocalMatrix, and the Slater integral is
        R^k(abcd) = b_a^T M b_c with M = InteractionStrength.XL_CoulombKinkAware(k, a, orb_b, c, orb_d,
        primitives) -- the matrix-valued overload that already exists for the Fock build. Each slot in which
        a subshell occurs contributes once, using the symmetry R^k(abcd) = R^k(badc) for the second pair.
        Must be checked against gradientByFiniteDifference before being trusted.
        A Dict{Subshell, Vector{Float64}} is returned.
"""
function computeOrbitalGradient(bVectors::Dict{Subshell, Vector{Float64}},
                                       coeffs1p::Array{Coefficient1p,1},
                                       coeffs2p::Array{Coefficient2p,1},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
    end
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells    grad[sh] = zeros(nsL+nsS)    end

    expanded = Dict{Subshell, Tuple{Vector{Float64},Vector{Float64}}}()
    scale    = Dict{Subshell, Float64}()
    for  sh  in  subshells
        (pR, qR) = SelfConsistent.expandBVector(bVectors[sh], primitives)
        np = min(length(pR), length(orbitals[sh].P));    nq = min(length(qR), length(orbitals[sh].Q))
        num = sum( orbitals[sh].P[1:np] .* pR[1:np] ) + sum( orbitals[sh].Q[1:nq] .* qR[1:nq] )
        den = sum( pR[1:np] .* pR[1:np] )             + sum( qR[1:nq] .* qR[1:nq] )
        scale[sh]    = den > 0. ? num/den : 1.0
        expanded[sh] = (orbitals[sh].P, orbitals[sh].Q)
    end

    # The screened potential depends only on the rank and the ORBITAL PAIR (b,d), not on the slot being
    # differentiated, and many angular coefficients share the same triple.  The orbitals are fixed for the
    # whole of one gradient evaluation, so the memo is local to this call and cannot go stale.  This is the
    # same redundancy that dominates the average-level Fock build (3.5x for argon, 9.1x for Th+), and after
    # the average-level field was memoised the rotation became the larger part of EOL: 66% of it for W+.
    vkCache = Dict{Tuple{Int64,Subshell,Subshell}, Vector{Float64}}()

    # one-electron:  d/db_a [ w * b_a^T H1 b_b ]  contributes to BOTH slots
    h1 = Dict{Int64, Array{Float64,2}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        h1[kappa] = Bsplines.setupLocalMatrix(kappa, primitives, nucPot, storage)
    end
    for  cf  in  coeffs1p
        if  cf.a.kappa != cf.b.kappa    continue    end          ## no one-electron element between kappas
        # The SAME chain-rule factor the two-electron part applies, and for the same reason.  The energy's
        # one-electron term is built from the ORBITALS, and an orbital is scale*expand(b) -- so
        # GrantIab(orb_a, orb_b) = scale_a * scale_b * (b_a^T H1 b_b) and the derivative carries that product.
        # Omitting it is harmless whenever a and b share a sign, which is why it went unnoticed: a DIAGONAL
        # term has scale^2 = +1 always, and a kappa whose subshells happen to be canonicalized alike gives +1
        # too.  It bites only on an OFF-DIAGONAL element between two subshells of OPPOSITE sign, and then it
        # flips the sign of that contribution outright.  Measured on a three-layer Be RAS, validating the
        # gradient along guaranteed-tangent directions: every kappa came out exact to 1.0000 except kappa = -2,
        # the single kappa holding an off-diagonal pair -- 2p_3/2 with scale +1 and 3p_3/2 with scale -1 --
        # where the ratio of finite difference to prediction ran -0.10, -0.56, -0.024.
        w  = cf.T * scale[cf.a] * scale[cf.b]
        hh = h1[cf.a.kappa]
        grad[cf.a] = grad[cf.a] + w * (hh * bVectors[cf.b])
        grad[cf.b] = grad[cf.b] + w * (transpose(hh) * bVectors[cf.a])
    end

    # two-electron:  R^k(abcd) = b_a^T M(a,orb_b;c,orb_d) b_c  and, by R^k(abcd) = R^k(badc),
    #                          = b_b^T M(b,orb_a;d,orb_c) b_d
    #
    # MATRIX-FREE (10-Aug-2026).  M is never formed.  Its entries are INT B_i B_k w_r V_L, so
    # (M v)_i = INT B_i(r) f(r) w_r V_L(r) with f the expansion of v -- and M is symmetric within each
    # block, so the two products a coefficient needs share one potential and differ only in the vector.
    # Each subshell's expansion is built ONCE per gradient rather than per coefficient.
    # The expansions must carry the SAME sign convention as the orbitals the energy is built from.
    # Bsplines.generateOrbitalFromVector canonicalizes each orbital so that P > 0 at small r, and does NOT
    # feed that back into bVectors -- so a raw expansion and its own orbital can be oppositely signed. The
    # screened potential Vk below comes from the ORBITALS (b,d) while the contracted vector came from the
    # RAW b-vector (a,c), and an off-diagonal CSF pair carries an ODD power of a correlating orbital's sign,
    # so the mismatch survives instead of cancelling. Measured on the Be RAS step-2 case, where 2p_1/2 and
    # 2p_3/2 are canonicalized against their raw vectors in every iteration: the analytic gradient agreed
    # with a finite difference to five digits while only s-orbitals were involved, then went 2.5x, 19x and
    # finally SIGN-WRONG as the 2p weight grew -- which is what stalled the line search. Same defect, and
    # the same remedy, as the cVector sign-matching in computeTwoElectronV.
    for  cf  in  coeffs2p
        for  (sA, sB, sC, sD)  in  [ (cf.a, cf.b, cf.c, cf.d), (cf.b, cf.a, cf.d, cf.c) ]
            # the same triangular-delta and parity guard XL_CoulombKinkAware applies before doing any work
            lA = Basics.subshell_l(sA);   jA = Basics.subshell_2j(sA)
            lB = Basics.subshell_l(sB);   jB = Basics.subshell_2j(sB)
            lC = Basics.subshell_l(sC);   jC = Basics.subshell_2j(sC)
            lD = Basics.subshell_l(sD);   jD = Basics.subshell_2j(sD)
            if  AngularMomentum.triangularDelta(jA+1, jC+1, cf.nu+cf.nu+1) *
                AngularMomentum.triangularDelta(jB+1, jD+1, cf.nu+cf.nu+1) == 0   ||
                rem(lA+lC+cf.nu, 2) == 1   ||   rem(lB+lD+cf.nu, 2) == 1     continue
            end
            xc = AngularMomentum.CL_reduced_me(sA, cf.nu, sC) * AngularMomentum.CL_reduced_me(sB, cf.nu, sD)
            if  rem(cf.nu, 2) == 1    xc = - xc    end
            Vk = get!(vkCache, (cf.nu, sB, sD)) do
                     RadialIntegrals.buildScreenedPotential(cf.nu, orbitals[sB], orbitals[sD], primitives.grid;
                                                                  mtpOut=primitives.grid.NoPoints)
                 end
            (pC, qC) = expanded[sC];     (pA, qA) = expanded[sA]
            grad[sA] = grad[sA] + (cf.V * xc * scale[sA]) * SelfConsistent.screenedProduct(Vk, pC, qC, primitives)
            grad[sC] = grad[sC] + (cf.V * xc * scale[sC]) * SelfConsistent.screenedProduct(Vk, pA, qA, primitives)
        end
    end
    return( grad )
end


"""
`SelfConsistent.gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
        virtuals::Dict{Subshell, Array{Vector{Float64},1}}, coeffs1p::Array{Coefficient1p,1},
        coeffs2p::Array{Coefficient2p,1}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)`  
    ... the gradient of the EOL energy with respect to the allowed orbital rotations, by central finite
        differences. Slow but free of any derivation, so it is the reference against which an analytic
        gradient must later be checked, and on its own it answers the question "is this converged solution
        actually stationary?". Each direction is S-orthogonal to the orbital itself, so normalization
        contributes only at second order and is omitted. A Dict{Subshell, Vector{Float64}} is returned.
"""
function gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
                                           virtuals::Dict{Subshell, Array{Vector{Float64},1}},
                                           coeffs1p::Array{Coefficient1p,1},
                                           coeffs2p::Array{Coefficient2p,1},
                                           subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                           grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells
        dirs = virtuals[sh];    g = zeros( length(dirs) )
        for  (iv, phi)  in  enumerate(dirs)
            bPlus  = copy(bVectors);    bPlus[sh]  = bVectors[sh] + hStep * phi
            bMinus = copy(bVectors);    bMinus[sh] = bVectors[sh] - hStep * phi
            ePlus  = SelfConsistent.energyFromBVectors(bPlus,  coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            eMinus = SelfConsistent.energyFromBVectors(bMinus, coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            g[iv]  = (ePlus - eMinus) / (2 * hStep)
        end
        grad[sh] = g
    end
    return( grad )
end


"""
`SelfConsistent.projectOntoPositiveBranch(bVectors::Dict{Subshell, Vector{Float64}},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, nucPot::Radial.Potential,
        matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}})`  
    ... projects every orbital onto the POSITIVE-ENERGY branch of its own kappa block, then S-orthonormalizes
        the orbitals of each kappa among themselves, entirely in B-spline coefficient space. Returns the
        projected (bVectors, maximum negative-branch weight removed).

        This is what keeps the rotation from collapsing into the Dirac sea. The Dirac operator is unbounded
        below, so the energy has no minimum once an orbital acquires negative-continuum character. The
        rotation DIRECTIONS are all positive-branch by construction, but the round trip
        b -> Bsplines.generateOrbitalFromVector -> b truncates at mtp and cleans small values, and that
        residue carries negative-branch content. Measured for Be with 48 directions: the 2s negative-branch
        weight jumps from 5e-12 to 1.2e-04 in one step and to 3.3e-03 in the next, at which point <h1> for 2s
        has gone from -1.6 to -43.5 and the energy to -69.9 Ha. Projecting after every step removes the
        residue before it can be amplified. The orthonormalization is done here, in coefficient space, rather
        than via orthonormalizeSameKappa, precisely to avoid that lossy round trip.
"""
function projectOntoPositiveBranch(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
                                          primitives::Bsplines.Primitives, nucPot::Radial.Potential,
                                          matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}};
                                          spectrum=nothing)
    out    = Dict{Subshell, Vector{Float64}}()
    worst  = 0.
    posSpec = isnothing(spectrum) ?
              SelfConsistent.positiveBranchSpectrum(subshells, primitives, nucPot, matrixB, storage) : spectrum
    posSet  = Dict{Int64, Array{Vector{Float64},1}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )    posSet[kappa] = posSpec[kappa][1]    end
    # (1) project each orbital on the positive branch of its kappa
    for  sh  in  subshells
        b = bVectors[sh];    v = zeros( length(b) )
        for  phi  in  posSet[sh.kappa]    v = v + (transpose(phi) * matrixB * b) * phi    end
        nrm2Full = abs( transpose(b) * matrixB * b );    nrm2Pos = abs( transpose(v) * matrixB * v )
        worst    = max( worst, 1.0 - nrm2Pos/max(nrm2Full,1.0e-30) )
        out[sh]  = v / sqrt( max(nrm2Pos, 1.0e-30) )
    end
    # (2) S-orthonormalize within each kappa, in coefficient space, so the positive span is preserved --
    # but ONLY where it is actually needed.  Gram-Schmidt is sequential and asymmetric: it leaves the first
    # orbital of a kappa untouched and pushes the whole correction onto the later ones, so applying it to
    # an already-orthonormal set rotates the orbitals for nothing.  Measured on Li, doing it unconditionally
    # cost 8.5e-07 Ha while removing a negative-branch weight of only 3.4e-14 -- more than the entire
    # discrepancy that sent us looking.  Skip it when the block is orthonormal to tolerance.
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        shk = [ sh for sh in subshells if sh.kappa == kappa ]
        dev = 0.
        for  (i, sha) in enumerate(shk),  (j, shb) in enumerate(shk)
            ov  = transpose(out[sha]) * matrixB * out[shb]
            dev = max( dev, abs( ov - (i == j ? 1.0 : 0.0) ) )
        end
        if  dev < 1.0e-9    continue    end
        # SYMMETRIC (Loewdin) orthogonalisation, S^(-1/2), in place of Gram-Schmidt.  Gram-Schmidt is
        # sequential and asymmetric: it leaves the FIRST orbital of a kappa untouched and pushes the whole
        # correction onto the later ones, so perturbing an earlier orbital of a kappa silently moves the later
        # ones as well.  Loewdin treats every orbital of the block alike and is the orthonormal set CLOSEST to
        # the input in a least-squares sense, so it introduces no ordering of its own.
        nk  = length(shk)
        ovl = zeros(nk, nk)
        for  (i, sha) in enumerate(shk),  (j, shb) in enumerate(shk)
            ovl[i,j] = transpose(out[sha]) * matrixB * out[shb]
        end
        ovl = 0.5 * (ovl + transpose(ovl))          ## exact symmetry before the eigendecomposition
        wa  = LinearAlgebra.eigen(ovl)
        sinv = zeros(nk, nk)
        for  k = 1:nk
            if  wa.values[k] < 1.0e-12    continue    end     ## a linearly dependent block keeps its input
            sinv = sinv + (wa.vectors[:,k] * transpose(wa.vectors[:,k])) / sqrt(wa.values[k])
        end
        newv = [ zeros(length(out[shk[1]]))  for i = 1:nk ]
        for  i = 1:nk,  j = 1:nk    newv[i] = newv[i] + sinv[j,i] * out[shk[j]]    end
        for  (i, sh) in enumerate(shk)
            nrm = sqrt( abs(transpose(newv[i]) * matrixB * newv[i]) )
            if  nrm > 1.0e-10    out[sh] = newv[i] / nrm    end
        end
    end
    return( out, worst )
end


"""
`SelfConsistent.solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model,
        primitives::Bsplines.Primitives, settings::AsfSettings; printout::Bool=true)`  
    ... EOL by direct minimization of the energy over orbital ROTATIONS, instead of by solving a per-subshell
        generalized eigenvalue problem. It sits beside the older solveOptimizedLevelField, which is left
        untouched; this one is what SelfConsistent.performSCF dispatches to for Basics.EOLField, and hence
        what every layer of a RasExpansion runs.

        Why: the eigenvalue formulation is ill-posed exactly when a correlating CSF's weight becomes small.
        Measured (09-Aug-2026) for Be 1s^2 2s^2 + 1s^2 2p^2, the present scheme converges to a DEGENERATE
        stationary point -- the 2p^2 weight collapses to -4e-5, so the 2p orbital no longer enters the energy
        at all and its gradient vanishes for a trivial reason (|g| = 9e-7). Its energy is 0.0195 Ha above the
        true minimum. Minimizing over rotations cannot fall into that trap: the orbital is driven by the
        energy gradient rather than selected from a spectrum, and it is free to change character -- which a
        correlation orbital must be, and which maximum-overlap selection (tested, refuted) forbids.

        Each outer step re-diagonalizes the CI, rebuilds the angular coefficients from the target level(s),
        and takes one backtracking-line-search step along the negative gradient projected on the
        active-virtual rotations. Occupied-occupied rotations are excluded as redundant.
        A multiplet::Multiplet of the target block(s) is returned.

        settings.frozenSubshells IS HONOURED: the orbitals it names are held exactly fixed and take no part in
        the gradient, the search direction, the curvature history or the line search, while still entering the
        energy and the CI matrix and still being orthogonalized against. That distinction is what a RAS layer
        means by freezing an orbital, and Basics.generate sets the field for every step of a RasExpansion.
        Freezing everything is allowed and is reported as such: the returned multiplet is then the CI result on
        the orbitals as given.
"""
function solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true, nVirtual::Int64=16,
                                         method::Symbol=:lbfgs, cgRestart::Int64=10, lbfgsMemory::Int64=12)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)

    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitives(basis.orbitals[sh], primitives, matrixB)
    end

    if  settings.levelSelectionCI.active  &&  !isempty(settings.levelSelectionCI.symmetries)
        relevantSyms = unique( settings.levelSelectionCI.symmetries )
    else
        relevantSyms = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    end
    blockCaches = Dict()
    for  sym  in  relevantSyms    blockCaches[sym] = SelfConsistent.cacheCsfPairCoefficientsEOL(sym, basis)   end

    # Built ONCE: it depends only on the nuclear potential and the B-spline basis, never on the orbitals.
    posSpectrum = SelfConsistent.positiveBranchSpectrum(basis.subshells, primitives, nucPot, matrixB, storage)
    (bVectors, _) = SelfConsistent.projectOntoPositiveBranch(bVectors, basis.subshells, primitives,
                                                                     nucPot, matrixB, storage; spectrum=posSpectrum)
    # FROZEN SUBSHELLS.  settings.frozenSubshells names orbitals that must NOT be varied, which is what a RAS
    # layer means by "frozen" -- Basics.generate sets it for every step of a RasExpansion, so this is the
    # ordinary case for this driver and not an exotic one.  The optimizer below therefore runs over
    # activeSubshells alone, while the ENERGY, the CI matrix, the virtual directions and the positive-branch
    # projection continue to see EVERY subshell: a frozen orbital still contributes to the energy and still
    # has to be orthogonalized against.  Nothing has to be re-imposed afterwards, because virtualDirections
    # builds each subshell's directions orthogonal to all occupied orbitals of the same kappa, frozen ones
    # included -- so a step taken by an active orbital leaves the kappa block orthonormal by construction.
    # The pinned vectors are taken AFTER the initial projection, so that they sit on the positive branch like
    # everything else, and are written back after every later projection: projectOntoPositiveBranch
    # renormalizes, and its Gram-Schmidt would rotate them on the occasions it fires.
    frozenSubshells = [ sh  for sh in basis.subshells  if    sh in settings.frozenSubshells ]
    activeSubshells = [ sh  for sh in basis.subshells  if  !(sh in settings.frozenSubshells) ]
    pinnedB         = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in frozenSubshells )
    restoreFrozen!  = function(bs)
        for  sh  in  frozenSubshells    bs[sh] = copy(pinnedB[sh])    end
        return( bs )
    end
    if  printout  &&  !isempty(frozenSubshells)
        println(">> [EOL-C3] frozen and NOT optimized: " * join(string.(frozenSubshells), ", ") *
                ";  $(length(activeSubshells)) of $(length(basis.subshells)) subshells are varied.")
    end
    ePrevious = 0.;   tStep = 1.0;   multiplet = Multiplet("EOL-ByRotation", Level[])
    # Set by every exit below.  A loop that simply runs out of iterations used to end in silence, which was the
    # fifth of five ways this driver can stop and the only one left unreported.
    stopReason = "";   gNorm = 0.;   iterDone = 0
    # Direction state, all held in b-space: the virtual directions are rebuilt every iteration, so anything
    # stored in THAT basis would be meaningless one step later.
    dirPrev = Dict{Subshell, Vector{Float64}}();   gPrev = Dict{Subshell, Vector{Float64}}()
    sgPrev  = 0.;    iterSinceRestart = 0
    bPrev   = Dict{Subshell, Vector{Float64}}()
    sHist   = Vector{Dict{Subshell, Vector{Float64}}}();   yHist = Vector{Dict{Subshell, Vector{Float64}}}()
    rhoHist = Float64[]
    for  iter = 1:settings.maxIterationsScf
        orbitals = Dict{Subshell, Orbital}()
        for  sh  in  basis.subshells
            orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
        end
        tempBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
        radial1p  = Dict{Tuple{Subshell,Subshell},Float64}()
        radial2p  = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
        levels    = Level[]
        for  sym  in  relevantSyms
            (idxCsf, cache1p, cache2p) = blockCaches[sym]
            mtx = SelfConsistent.buildCIMatrixEOL(idxCsf, cache1p, cache2p, orbitals, grid, nucPot, radial1p, radial2p)
            append!( levels, SelfConsistent.diagonalizeBlockEOL(sym, idxCsf, mtx, tempBasis) )
        end
        multiplet    = Basics.sortByEnergy( Multiplet("EOL-ByRotation", levels) )
        targetLevels = SelfConsistent.selectTargetLevelsEOL(multiplet, settings.levelSelectionCI)
        (coeffs1p, coeffs2p) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)

        # e0 MUST be the energy of the point the line search starts from, measured the way the line search
        # measures its trials -- on the POSITIVE-BRANCH-PROJECTED vectors.  Taking it from the raw bVectors
        # instead made the comparison `eTrial < e0` a comparison between two different functions, and once
        # the two drifted apart no step could ever win: traced on a three-layer Be RAS, the driver reported
        # "no descent found" at iteration 61 while e0 sat 5.42 mHa BELOW the energy of its own starting
        # point, and the trial energy was flat in t across five decades because every t was being measured
        # against that offset.  The direction was never at fault -- dg = -0.047 there, and its tangential
        # part -0.046.  Re-projecting here is a no-op whenever the iterate is already on the manifold, which
        # is what makes it safe: it costs one projection per iteration and removes a whole class of stall.
        (bVectors, _) = SelfConsistent.projectOntoPositiveBranch(bVectors, basis.subshells, primitives,
                                                        nucPot, matrixB, storage; spectrum=posSpectrum)
        restoreFrozen!(bVectors)
        e0   = SelfConsistent.energyFromBVectors(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                         primitives, grid, nucPot)
        grad = SelfConsistent.computeOrbitalGradient(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                             primitives, nucPot, storage)
        virt = SelfConsistent.virtualDirections(bVectors, basis.subshells, primitives, nucPot,
                                                        matrixB, storage; nVirtual=nVirtual, spectrum=posSpectrum)
        # Project the gradient on the allowed rotations, and PRECONDITION each component by the
        # orbital-energy denominator (eps_v - eps_a) -- the diagonal of the rotation Hessian, i.e. the
        # standard first-order estimate  kappa_av = -g_av / (eps_v - eps_a).  Plain steepest descent
        # converges far too slowly here: it left the Li control 1.8e-6 Ha short after 60 steps.
        # The denominator is floored, since a near-degenerate pair would otherwise produce a huge step.
        gProj = Dict{Subshell, Vector{Float64}}();   step = Dict{Subshell, Vector{Float64}}()
        denom = Dict{Subshell, Vector{Float64}}()
        gNorm = 0.;    sNorm = 0.
        for  sh  in  activeSubshells
            h1k  = Bsplines.setupLocalMatrix(sh.kappa, primitives, nucPot, storage)
            epsA = transpose(bVectors[sh]) * h1k * bVectors[sh]
            gv   = [ transpose(phi) * grad[sh]  for phi in virt[sh] ]
            sv   = zeros( length(gv) );    dv = zeros( length(gv) )
            for  (iv, phi)  in  enumerate(virt[sh])
                epsV    = transpose(phi) * h1k * phi
                dv[iv]  = max( epsV - epsA, 0.05 )
                sv[iv]  = - gv[iv] / dv[iv]
            end
            gProj[sh] = gv;    step[sh] = sv;    denom[sh] = dv
            gNorm = gNorm + sum( gv.^2 );    sNorm = sNorm + sum( sv.^2 )
        end
        gNorm = sqrt(gNorm);    sNorm = sqrt(sNorm);    iterDone = iter
        # Every subshell frozen is a legitimate request, not an error: the multiplet built at the top of this
        # iteration IS the answer, being the CI result on the given orbitals.  Reported separately because the
        # collapsed-direction exit below would otherwise fire and call it a failure.
        if  isempty(activeSubshells)
            stopReason = "all subshells frozen";   println(">> [EOL-C3] every subshell of this basis is listed in " *
                    "settings.frozenSubshells, so there is nothing to optimize;  the multiplet is the CI result " *
                    "on the orbitals as given.")
            break
        end

        # Assemble the search direction in b-space.  Plain preconditioned steepest descent zigzags here:
        # the energy falls steadily while |grad| merely oscillates, each step undoing part of the previous
        # one.  Polak-Ribiere conjugacy reuses the previous direction to cancel that, at the cost of two
        # dot products and one stored vector per subshell.
        sVec = Dict{Subshell, Vector{Float64}}();   gVec = Dict{Subshell, Vector{Float64}}()
        for  sh  in  activeSubshells
            sv = zeros(nsL+nsS);    gv = zeros(nsL+nsS)
            for  (iv, phi)  in  enumerate(virt[sh])
                sv = sv + step[sh][iv]  * phi
                gv = gv + gProj[sh][iv] * phi
            end
            sVec[sh] = sv;    gVec[sh] = gv
        end
        # method = :conjugate is the DEFAULT and the one to use.  :lbfgs is kept because it is measurably
        # better where the basis is stable -- Be 1s^2 2s^2 + 1s^2 2p^2 at 12 iterations reaches -14.618710
        # against conjugacy's -14.616507, i.e. it does change the RATE and not merely the constant -- but it
        # LOSES on the harder Be RAS step-2 case, -14.617374 against -14.619313, and there it discards its
        # own curvature history at iterations 8, 15, 21 and 23.  WHY IT FAILS IS NOT KNOWN.  Two candidates
        # were proposed and BOTH MEASURED SMALL, so neither should be repeated as an explanation:
        #   * the virtual space is NOT churning -- successive frames overlap to |1-<new,old>| = 1e-4..3e-3
        #     with no change in the number of directions, one 0.67 rotation excepted, and that one does not
        #     coincide with any of the four history discards;
        #   * the objective DOES drift, since coeffs1p/coeffs2p are rebuilt from the re-diagonalized mixing
        #     vector every iteration, but only by max|dV| ~ 1e-4..3e-3 against max|V| ~ 1.8.
        # Conjugacy tolerates whatever this is because it carries ONE previous direction and restarts every
        # ten; L-BFGS accumulates five pairs and does not.  Anyone picking this up should measure first.
        dotAll(u, v) = sum( sum(u[sh] .* v[sh])  for sh in activeSubshells )
        # H_0 is the DIAGONAL PRECONDITIONER, not the usual gamma*I: the (eps_v - eps_a) denominators carry
        # real physics and throwing them away for a scalar would be a step backwards.  applyPrecond is what
        # turns gVec into -sVec, so it is exactly the operator already in use.
        applyPrecond = function(q)
            r = Dict{Subshell, Vector{Float64}}()
            for  sh  in  activeSubshells
                rv = zeros(nsL+nsS)
                for  (iv, phi)  in  enumerate(virt[sh])
                    rv = rv + ( (transpose(phi) * q[sh]) / denom[sh][iv] ) * phi
                end
                r[sh] = rv
            end
            return( r )
        end
        # Record the curvature pair (s,y) of the step just taken.  Only pairs with <y,s> > 0 are kept, which
        # is what keeps the implicit inverse Hessian positive definite -- and hence the direction a descent
        # direction -- on a surface that is not convex.
        if  method == :lbfgs  &&  !isempty(bPrev)
            sPair = Dict{Subshell, Vector{Float64}}();   yPair = Dict{Subshell, Vector{Float64}}()
            for  sh  in  activeSubshells
                sPair[sh] = bVectors[sh] - bPrev[sh];    yPair[sh] = gVec[sh] - gPrev[sh]
            end
            ys = dotAll(yPair, sPair)
            if  ys > 1.0e-14
                push!(sHist, sPair);    push!(yHist, yPair);    push!(rhoHist, 1.0/ys)
                while  length(sHist) > lbfgsMemory
                    popfirst!(sHist);   popfirst!(yHist);   popfirst!(rhoHist)
                end
            end
        end
        beta = 0.
        dir  = Dict{Subshell, Vector{Float64}}()
        if      method == :lbfgs  &&  !isempty(sHist)
            # two-loop recursion, giving d = -H grad with H built from the stored pairs around H_0
            q = Dict{Subshell, Vector{Float64}}( sh => copy(gVec[sh])  for sh in activeSubshells )
            alphas = zeros( length(sHist) )
            for  i = length(sHist):-1:1
                alphas[i] = rhoHist[i] * dotAll(sHist[i], q)
                for  sh  in  activeSubshells    q[sh] = q[sh] - alphas[i] * yHist[i][sh]    end
            end
            r = applyPrecond(q)
            for  i = 1:length(sHist)
                bb = rhoHist[i] * dotAll(yHist[i], r)
                for  sh  in  activeSubshells    r[sh] = r[sh] + (alphas[i] - bb) * sHist[i][sh]    end
            end
            for  sh  in  activeSubshells    dir[sh] = -r[sh]    end
        elseif  method == :conjugate
            if  !isempty(dirPrev)  &&  iterSinceRestart < cgRestart  &&  abs(sgPrev) > 1.0e-30
                num = 0.
                for  sh  in  activeSubshells    num = num + sum( (gVec[sh] - gPrev[sh]) .* sVec[sh] )    end
                beta = max( 0., num / sgPrev )                   ## Polak-Ribiere+, i.e. restart on beta < 0
            end
            for  sh  in  activeSubshells
                dir[sh] = beta == 0. ? sVec[sh] : sVec[sh] + beta * dirPrev[sh]
            end
        else
            for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
        end
        # Guard: a conjugate direction must still descend.  If it does not, fall back to steepest descent.
        dg = 0.;   for sh in activeSubshells   dg = dg + sum( gVec[sh] .* dir[sh] )   end
        if  dg >= 0.
            for  sh  in  activeSubshells    dir[sh] = sVec[sh]    end
            beta = 0.
        end
        iterSinceRestart = beta == 0. ? 0 : iterSinceRestart + 1
        sgPrev = 0.;   for sh in activeSubshells   sgPrev = sgPrev + sum( gVec[sh] .* sVec[sh] )   end
        gPrev  = gVec;    dirPrev = dir
        bPrev  = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in activeSubshells )
        if  sNorm < 1.0e-14
            # Until 18-Aug-2026 this was the one exit of four that said NOTHING, so a run could end here and
            # be read as a completed optimisation.  It means the PRECONDITIONED direction has collapsed, which
            # is not the same as a converged gradient and must not be reported as one.
            stopReason = "direction collapsed";   println(">> [EOL-C3] STOPPED at iteration $iter: the preconditioned direction has collapsed " *
                    "(|s| = $sNorm), with |grad| = $gNorm.  This is NOT convergence.")
            break
        end
        if  printout
            println(">> [EOL-C3] iter $iter:  E = $(multiplet.levels[1].energy)   |grad| = $gNorm   step = $tStep")
        end
        # Converged when the GRADIENT is small. settings.accuracyScf is its tolerance, which is the honest
        # test: the old sole criterion, |E - E_prev| < accuracyScf, halts when PROGRESS is slow, which is a
        # statement about the optimizer and not about the solution -- and it is why every EOL value quoted
        # before 16-Aug-2026 was an upper bound. It cannot stand ALONE, though: |grad| PLATEAUS at a floor
        # set by the basis and the projection, so on Li a pure gradient test ran 48 further iterations after
        # the energy had stopped moving, for nothing. Both tests are kept, and the driver says which fired.
        if  gNorm < settings.accuracyScf
            stopReason = "converged";   println(">> [EOL-C3] CONVERGED at iteration $iter: |grad| = $gNorm < accuracyScf = " *
                    "$(settings.accuracyScf), tStep = $tStep.")
            break
        end

        # backtracking line search along -gProj; halve until the energy actually falls
        accepted = false
        for  trial = 1:24
            newB = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
            for  sh  in  activeSubshells
                v = bVectors[sh] + tStep * dir[sh]
                newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
            end
            # The projection must happen BEFORE the acceptance test, not after it.  Projecting and
            # re-orthonormalizing moves the orbitals, so accepting `newB` on the strength of its own energy
            # and then storing the PROJECTED vector stores something that was never tested -- and the
            # projection gives a little of the gain back each step.  That broke monotonicity: on Li the
            # driver reduced the gradient 19-fold while the energy ROSE by 3.7e-7 Ha, which a descent
            # method cannot do.  Testing the projected vector restores  E(new) < E(old)  by construction,
            # and with it the guarantee that the CI eigenvalue falls too (it is bounded above by this
            # fixed-coefficient functional, and equals it at the previous orbitals).
            (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
            eTrial = SelfConsistent.energyFromBVectors(projB, coeffs1p, coeffs2p, basis.subshells,
                                                               primitives, grid, nucPot)
            if  eTrial < e0
                bVectors = restoreFrozen!(projB)
                if  printout  &&  negW > 1.0e-8
                    println(">> [EOL-C3] removed negative-branch weight $negW from the step.")
                end
                accepted = true;    tStep = min(1.0, 1.3*tStep);    break
            end
            tStep = tStep / 2
        end
        if  printout
        end
        if  !accepted  &&  method == :lbfgs  &&  !isempty(sHist)
            # A line search that finds no descent along an L-BFGS direction means the stored curvature is
            # no longer describing this surface -- unsurprising, since virtualDirections is rebuilt every
            # iteration and the pairs then mix vectors from different subspaces. Discard the history and
            # retry from the preconditioned gradient rather than giving up.
            if  printout    println(">> [EOL-C3] L-BFGS history discarded at iteration $iter; retrying.")   end
            empty!(sHist);   empty!(yHist);   empty!(rhoHist);   tStep = 1.0
            for  trial = 1:24
                newB = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in basis.subshells )
                for  sh  in  activeSubshells
                    v = bVectors[sh] + tStep * sVec[sh]
                    newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
                end
                (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
                eTrial = SelfConsistent.energyFromBVectors(projB, coeffs1p, coeffs2p, basis.subshells,
                                                                   primitives, grid, nucPot)
                if  eTrial < e0    bVectors = restoreFrozen!(projB);   accepted = true;   break    end
                tStep = tStep / 2
            end
        end
        if  !accepted
            stopReason = "no descent";   println(">> [EOL-C3] STOPPED at iteration $iter: no descent found along the search direction " *
                    "(tStep fell to $tStep), with |grad| = $gNorm.  This is NOT convergence.")
            break
        end
        # A STAGNANT ENERGY IS ONLY CONVERGENCE IF THE STEP IS STILL HEALTHY.  The test compares successive
        # CI eigenvalues, and a variational energy is STATIONARY at a minimum: |dE| falls as the SQUARE of the
        # orbital error while |grad| falls linearly, so |dE| < 1e-11 is reached long before convergence
        # whenever the steps have become small -- and then it reports a collapsed line search as a converged
        # calculation.  Measured at the moment it used to fire: tStep = 4.1e-9 on a Be RAS correlation layer
        # and 3.5e-11 on a carbon one, i.e. down eight to eleven orders from unity, with |grad| still 0.037
        # and 0.998.  The collapse is TEMPORARY -- allowed to continue, carbon recovers a step of 0.125 and
        # converges at |grad| = 4.7e-6, sixty-five mHa BELOW where it used to stop.  So the remedy is not a
        # smaller threshold but the extra condition: stagnation ends the iteration only when the step that
        # produced it was of usable size.
        stepFloor = 1.0e-6
        if  iter > 1  &&  abs(multiplet.levels[1].energy - ePrevious) < 1.0e-11  &&  tStep > stepFloor
            stopReason = "stationary energy";   println(">> [EOL-C3] stopped at iteration $iter on a stationary energy: |dE| = " *
                    "$(abs(multiplet.levels[1].energy - ePrevious)) < 1.0e-11 with a healthy step " *
                    "tStep = $tStep and |grad| = $gNorm.  A converging UPPER BOUND, not a converged gradient.")
            break
        end
        ePrevious = multiplet.levels[1].energy
    end
    if  stopReason == ""
        println(">> [EOL-C3] STOPPED after $iterDone iterations: the limit maxIterationsScf = " *
                "$(settings.maxIterationsScf) was reached with |grad| = $gNorm and tStep = $tStep.  " *
                "This is NOT convergence; raise maxIterationsScf to see where it goes.")
    end

    return( multiplet )
end


"""
`SelfConsistent.solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true)`
    ... solves the self-consistent field for the extended-optimal-level (EOL) functional: orbitals are
        optimized against a statistically-(2J+1)-weighted combination of one or more target ASF levels,
        with the combination weights coming from the levels' own CI mixing coefficients -- never a
        user-supplied weight. The target level(s) are selected EXCLUSIVELY via settings.levelSelectionCI
        (either indices or symmetries, never both; if inactive/empty, defaults to indices=[1], a genuine
        OL/single-level computation). This nests a CI diagonalization inside the AL SCF loop
        (GRASP rmcdhf90's scf.f90, algorithm 5.1 of Froese Fischer, Comput. Phys. Rep. 3 (1986) 290):
        diagonalize -> build generalized coefficients from the target levels' mixing vectors -> refine
        every orbital (reusing computeFockMatrix/computeTwoElectronV unchanged) ->
        re-diagonalize -> repeat, converging on both orbital self-consistency and the weighted-average
        energy. Per-CSF-pair angular coefficients (cacheCsfPairCoefficientsEOL) are computed once and
        reused every outer iteration, since they depend only on the fixed CSF list, never on the current
        orbitals or mixing coefficients. A (new) multiplet::Multiplet is returned.

        KNOWN LIMITATION (confirmed 28-Jul-2026, NOT yet fixed): when two or more CSFs of the SAME symmetry
        block compete for the same correlation channel (e.g. Be's 1s^2 2p_1/2^2 vs 1s^2 2p_3/2^2, both
        correlating with 1s^2 2s^2), this implementation can converge to a spurious, winner-take-all fixed
        point where one competing CSF's mixing coefficient is driven to ~0 while a comparably-important
        partner is not -- confirmed against a DFS-Field reference occupying the same CSF space, which lands
        both more bound AND with both CSFs contributing substantially. Root cause: the off-diagonal
        (CSF-pair) coefficients folded into computeTwoElectronV's two-electron potential scale
        LINEARLY in a shrinking CSF's own mixing coefficient, while computeGeneralizedOccupationEOL's
        occupation (the (1.0/occ) divisor in computeFockMatrix) scales QUADRATICALLY in it -- so the
        ratio diverges as that coefficient shrinks, rather than settling. This is the concrete manifestation
        of the "DA/inhomogeneous-term mechanism" gap vs. GRASP's setcof.f90 (which treats within-level
        off-diagonal coupling as a separate inhomogeneous/source term, not folded into the same per-orbital
        homogeneous eigenvalue division) -- see project_eol_implementation.md. Flooring `occ` before the
        division was tried and REJECTED as a fix (non-monotonic in the floor constant). Safe for
        single-CSF-per-block cases (validated: He, Li) and multi-CSF cases where every competing CSF's own
        weight stays comfortably bounded away from zero; NOT yet safe/reliable for genuine near-degenerate
        competing correlation (Be's 2p^2 case, and by extension most 3+ layer RAS scenarios). A real fix
        needs the actual inhomogeneous-term mechanism -- deferred, substantial future work.
"""
function solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                  settings::AsfSettings; printout::Bool=true)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid

    if  !settings.levelSelectionCI.active  ||  ( isempty(settings.levelSelectionCI.indices) && isempty(settings.levelSelectionCI.symmetries) )
        println(">> [EOL] No levelSelectionCI given; defaulting to indices = [1] -- a genuine OL (single lowest level) computation.")
    elseif  !isempty(settings.levelSelectionCI.indices)  &&  !isempty(settings.levelSelectionCI.symmetries)
        error("stop a; settings.levelSelectionCI must specify EITHER indices OR symmetries for the EOL scheme, not both.")
    end

    # Determine which symmetry block(s) need to be (re-) diagonalized every outer iteration: an explicit
    # symmetries-only selection only ever needs those blocks; index-based (or default) selection needs
    # every block present in the basis to resolve the global, energy-sorted ordering.
    if  settings.levelSelectionCI.active  &&  !isempty(settings.levelSelectionCI.symmetries)
        relevantSyms = unique( settings.levelSelectionCI.symmetries )
    else
        relevantSyms = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    end

    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once --
    # identical to solveAverageLevelField
    if  printout    println(">> [EOL] (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)

    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitives(basis.orbitals[sh], primitives, matrixB)
    end
    orbitals = basis.orbitals

    # (2) Cache the (orbital-independent) per-CSF-pair angular coefficients once for every relevant block
    if  printout    println(">> [EOL] Caching per-CSF-pair angular coefficients for symmetries $(relevantSyms) ...")    end
    blockCaches = Dict{LevelSymmetry, Tuple{Array{Int64,1}, Dict{Tuple{Int64,Int64},Array{Coefficient1p,1}},
                                             Dict{Tuple{Int64,Int64},Array{Coefficient2p,1}}}}()
    for  sym  in  relevantSyms
        blockCaches[sym] = SelfConsistent.cacheCsfPairCoefficientsEOL(sym, basis)
    end

    # (3) Initial CI diagonalization (starting orbitals) to get the first mixing vectors for the target level(s)
    function diagonalizeAllBlocks(currentOrbitals::Dict{Subshell, Orbital})
        tempBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, currentOrbitals)
        # One shared radial-integral cache per outer iteration, reused across every block: the SAME
        # subshell-labeled integral (e.g. "1s-1s") often recurs across many CSF pairs and even across
        # different symmetry blocks, but depends only on currentOrbitals, not on which block/pair asked for it.
        radial1pCache = Dict{Tuple{Subshell,Subshell},Float64}()
        radial2pCache = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
        levels = Level[]
        for  sym  in  relevantSyms
            (idxCsf, cache1p, cache2p) = blockCaches[sym]
            matrix = SelfConsistent.buildCIMatrixEOL(idxCsf, cache1p, cache2p, currentOrbitals, grid, nucPot,
                                                      radial1pCache, radial2pCache)
            append!( levels, SelfConsistent.diagonalizeBlockEOL(sym, idxCsf, matrix, tempBasis) )
        end
        mp = Basics.sortByEnergy( Multiplet("EOL", levels) )
        return( mp )
    end

    mp           = diagonalizeAllBlocks(orbitals)
    targetLevels = SelfConsistent.selectTargetLevelsEOL(mp, settings.levelSelectionCI)
    previousMc   = [ copy(level.mc)  for level in targetLevels ]
    if  SelfConsistent.GBL_EOL_UNSCALED_OFFDIAGONAL
        (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels; pairs=:diagonal)
        (coeffs1pOff, coeffs2pOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels; pairs=:offdiagonal)
    else
        (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)
        coeffs2pOff                = Coefficient2p[]
    end
    genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels, basis)

    weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
        sumW = sum( twiceJp1(level.J)  for level in targetLevels )
        [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    end
    weightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

    # (4) Precompute kink-aware Slater-moment tensor caches for every rank that occurs, exactly as
    # solveAverageLevelField does; only the exchange branches of computeTwoElectronV use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [EOL] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildScreenedPotentialCache(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    for  iter = 1:settings.maxIterationsScf
        println("\n> SCF+CI iteration $(iter) [EOL]: ")
        for  level  in  targetLevels
            println("   target level  J=$(level.J)  parity=$(level.parity)  energy=$(level.energy)")
        end

        newBVectors = Dict{Subshell, Vector{Float64}}()
        # Pre-seed with every FROZEN subshell's (fixed) bVector, before any active subshell is refined: this
        # makes Hamiltonian.projectHamiltonian orthogonalize active subshells against frozen same-kappa ones
        # (e.g. a frozen 1s vs. an active, higher-n correlation orbital of the same kappa) and keeps the
        # `count`-based target-eigenvalue-index shift correct, exactly as if the frozen orbitals had already
        # been "processed" this iteration -- which, since they never change, they effectively have.
        processedBVectors = Dict{Subshell, Vector{Float64}}( sh => bVectors[sh]  for sh in settings.frozenSubshells  if  sh in basis.subshells )
        dpm = Dict{Subshell, Float64}()
        # fresh each sweep, for the reason given in solveAverageLevelField: the partner orbitals are fixed
        # within a sweep but not between sweeps, so a cache carried over would serve stale matrices
        directKernels   = Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}()
        exchangeKernels = Dict{Tuple{Int64,Subshell},Array{Float64,2}}()

        for  subshell  in  basis.subshells
            if  subshell in settings.frozenSubshells
                # Carry the frozen bVector forward unchanged into newBVectors -- required, since newBVectors
                # is a fresh Dict every iteration and the final newOrbitals tabulation loop below reads
                # bVectors[sh] for EVERY sh in basis.subshells, frozen or not (a KeyError otherwise).
                newBVectors[subshell] = bVectors[subshell]
                continue
            end
            occ = genOcc[subshell]
            # A subshell can legitimately carry ZERO generalized occupation: the EOL functional is built from
            # the target level(s) alone, so a subshell that no target level has any weight on (e.g. 3d when the
            # only target is the J=0+ level of 1s^2 2s^2 + 1s^2 2p^2 + 1s^2 2s3s + 1s^2 2s3d) simply does not
            # enter the energy at all.  The functional is then stationary with respect to that orbital and
            # there is nothing to optimize -- so it is carried forward unchanged, exactly like a frozen one.
            # Refining it regardless used to produce a Fock matrix of NaN (computeFockMatrix divides by
            # occ, and Inf * 0 = NaN), which surfaced far downstream as the thoroughly misleading
            # "Bsplines.findPositiveBranchStart(): no eigenvalue found above the negative-continuum threshold"
            # -- the signature of a missing nuclear well, which was not the problem at all.  AL never meets
            # this because its MEAN occupation averages over every CSF and so is never zero.
            if  abs(occ) < 1.0e-12
                println(">> Subshell $subshell carries zero generalized occupation in the target level(s); " *
                        "the EOL functional does not depend on it, so it is kept unchanged.")
                newBVectors[subshell] = bVectors[subshell]
                continue
            end
            print(">> Refine $subshell orbital with generalized occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrix(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches; coeffs2pUnscaled=coeffs2pOff,
                                                              directKernels=directKernels,
                                                              exchangeKernels=exchangeKernels)

            count = Base.count( sh2 -> sh2.kappa == subshell.kappa, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Bsplines.diagonalizeLocalMatrix(subshell.kappa, matrix, matrixB, primitives)
            l  = Basics.subshell_l(subshell)
            mm = Bsplines.findPositiveBranchStart(wc.values)
            oldVector = bVectors[subshell]
            # TESTED AND REFUTED (09-Aug-2026): selecting the eigenvector of maximum OVERLAP with the
            # previous orbital instead of the counted index changes nothing here (identical to eight
            # decimals), so the counted index was never the problem -- and it actively harms a correlating
            # orbital, which must be free to change character.  Do not re-propose it.
            ni = mm + subshell.n - l - count - 1
            rawVector = wc.vectors[ni]

            if  transpose(oldVector) * matrixB * rawVector < 0    rawVector = -rawVector    end
            damping = 0.5
            mixed     = damping * oldVector + (1.0 - damping) * rawVector
            newVector = mixed / sqrt( transpose(mixed) * matrixB * mixed )

            newBVectors[subshell]       = newVector
            processedBVectors[subshell] = newVector

            ovlap = abs( transpose(oldVector) * matrixB * newVector )
            dpm[subshell] = 1.0 - ovlap
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end

        bVectors = newBVectors
        newOrbitals = Dict{Subshell, Orbital}()
        for  sh  in  basis.subshells
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
        end
        # The EOL driver damps exactly as the AL one does, so it loses same-kappa orthogonality in exactly
        # the same way and needs the same repair.  Measured before this was added: Si^2+ [Ne] 3s^2 + 3p^2
        # gave a worst same-kappa overlap of 2.4e-06 under EOL against 9.3e-10 under AL, and Si^+ reached
        # 6.9e-05.  Wiring the switch into solveAverageLevelField alone was an oversight.
        if  SelfConsistent.GBL_SCF_REORTHONORMALIZE
            (newOrbitals, bVectors) = SelfConsistent.orthonormalizeSameKappa(newOrbitals, bVectors,
                                                        basis.subshells, primitives, matrixB)
        end
        orbitals = newOrbitals

        # Re-diagonalize CI with the refined orbitals; refresh the target level(s)' mixing vectors and energies
        mp           = diagonalizeAllBlocks(orbitals)
        targetLevels = SelfConsistent.selectTargetLevelsEOL(mp, settings.levelSelectionCI)

        # Damping (27-Jul-2026): the outer CI-mixing refresh, feeding straight back into coeffs2p/genOcc every
        # iteration, reproduces the same kind of period-2 oscillation the inner bVector update needed damping
        # for (see solveAverageLevelField/project_df_al_kink_bug.md) -- here in the mixing-vector <->
        # orbital <-> generalized-occupation three-way loop instead of just orbitals <-> orthogonality. Same
        # standard fix: linear mixing of the new and previous mixing vector per target level before use,
        # aligning sign first (a CI eigensolver may return either sign) and renormalizing. The RAW (undamped)
        # targetLevels/energies are still used for reporting and the convergence check -- only the vectors
        # feeding coeffs2p/genOcc are damped.
        dampedLevels = Level[]
        for  (i, level)  in  enumerate(targetLevels)
            newMc = level.mc
            if  transpose(previousMc[i]) * newMc < 0    newMc = -newMc    end
            damping = 0.5
            mixed   = damping * previousMc[i] + (1.0 - damping) * newMc
            mixed   = mixed / sqrt( transpose(mixed) * mixed )
            push!( dampedLevels, Level(level.J, level.M, level.parity, level.index, level.energy,
                                        level.relativeOcc, level.hasStateRep, level.basis, mixed) )
            previousMc[i] = mixed
        end

        if  SelfConsistent.GBL_EOL_UNSCALED_OFFDIAGONAL
            (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels; pairs=:diagonal)
            (coeffs1pOff, coeffs2pOff) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels; pairs=:offdiagonal)
        else
            (coeffs1p, coeffs2p)       = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels)
            coeffs2pOff                = Coefficient2p[]
        end
        genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, dampedLevels, basis)

        weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
            sumW = sum( twiceJp1(level.J)  for level in targetLevels )
            [ twiceJp1(level.J) / sumW  for level in targetLevels ]
        end
        newWeightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

        orbitalConv = maximum( values(dpm) ) < 1.0 ? 1.0 - maximum( values(dpm) ) : 0.0
        energyDiff  = abs( newWeightedEnergy - weightedEnergy )
        println(">> Weighted-average energy = $newWeightedEnergy   orbital-conv = $orbitalConv   " *
                "orbital-acc = $(1.0 - orbitalConv)   energy-diff = $energyDiff")

        weightedEnergy = newWeightedEnergy
        if  abs(1.0 - orbitalConv) < settings.accuracyScf  &&  energyDiff < settings.accuracyScf    break   end
    end

    finalBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    multiplet  = Hamiltonian.performCIKinkAware(finalBasis, nuclearModel, grid, settings; printout=printout)
    return( multiplet )
end
