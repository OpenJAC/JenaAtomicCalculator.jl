
"""
`module  JAC.SelfConsistent`  
	... a submodel of JAC that contains all structs and methods to generate self-consistent fields of different 
	    kind and complexity.
"""
module SelfConsistent

using  Printf, ..Basics, ..Bsplines, ..Defaults, ..Hamiltonian, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial,
       ..RadialIntegrals, ..SpinAngular


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
    ncsf = length(basis.csfs);    coeffs1p = SpinAngular.Coefficient1p[];     coeffs2p = SpinAngular.Coefficient2p[] 
    
    # Compute angular coefficients in turn for all diagonal ME
    for  csf  in  basis.csfs
        coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus, true), 
                                                       csf, csf, basis.subshells)
        # Add to the existing list
        for  cf in coeffs   push!(coeffs1p, SpinAngular.Coefficient1p(cf.nu, cf.a, cf.b, cf.T / ncsf) )   end
        
        coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.TwoParticleOperator(0, Basics.plus, true), 
                                                       csf, csf, basis.subshells)
        # Add to the existing lists
        for  cf in coeffs   push!(coeffs2p, SpinAngular.Coefficient2p(cf.nu, cf.a, cf.b, cf.c, cf.d, cf.V / ncsf) )   end
    end 

    # Condense angular coefficients if they refer to the same set of orbital; 
    # include symmetry <ab||cd> == <ba||dc> for symmetric interactions
    coeffs1px = SpinAngular.Coefficient1p[];     coeffs2px = SpinAngular.Coefficient2p[]
    
    hasConsidered = falses( length(coeffs1p) );   T = 0.
    for  (ic, cf) in enumerate(coeffs1p)
        if    hasConsidered[ic]   
        else  nu = cf.nu;   a = cf.a;   b = cf.b
            T = T + cf.T;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs1p)
                if    hasConsidered[icx]   
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cf.b     T = T + cfx.T;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs1px, SpinAngular.Coefficient1p(nu, a, b, T) );  T = 0.
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
                ## elseif  nu == cfx.nu &&  a == cfx.b && b == cf.a &&  c == cfx.d &&  d == cf.c    
                ##         V = V + cfx.V;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs2px, SpinAngular.Coefficient2p(nu, a, b, c, d, V) );   V = 0.
        end
    end

    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.cacheCsfPairCoefficientsEOL(sym::LevelSymmetry, basis::Basis)`
    ... caches, for every CSF pair (r,s) with symmetry sym in the given basis, the pure spin-angular
        coefficients (independent of the orbitals/radial functions) as returned by
        SpinAngular.computeCoefficientsScalar -- the same call Hamiltonian.setupMatrix/setupMatrixClaude
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

    cache1p = Dict{Tuple{Int64,Int64}, Array{SpinAngular.Coefficient1p,1}}()
    cache2p = Dict{Tuple{Int64,Int64}, Array{SpinAngular.Coefficient2p,1}}()
    for  r = 1:n
        for  s = 1:n
            csfR = basis.csfs[idxCsf[r]];   csfS = basis.csfs[idxCsf[s]]
            cache1p[(r,s)] = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus, true),
                                                                    csfR, csfS, basis.subshells)
            cache2p[(r,s)] = SpinAngular.computeCoefficientsScalar(SpinAngular.TwoParticleOperator(0, Basics.plus, true),
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
        orbitals -- algebraically identical to Hamiltonian.setupMatrixClaude's pure-Coulomb contribution
        (kink-aware InteractionStrength.XL_CoulombClaude, matching the rest of the Claude2 SCF line; Breit
        and QED are added only once, at the final Hamiltonian.performCIClaude call, exactly as for the AL
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
    #
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
                jj = Basics.subshell_2j(cf.a)
                me = me + cf.T * sqrt( jj + 1) * I_ab
            end
            for  cf in cache2p[(r,s)]
                R_abcd = get!(radial2pCache, (cf.nu,cf.a,cf.b,cf.c,cf.d)) do
                    InteractionStrength.XL_CoulombClaude(cf.nu, orbitals[cf.a], orbitals[cf.b], orbitals[cf.c], orbitals[cf.d], grid)
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
        Hamiltonian.performCI/performCIClaude do internally per symmetry block -- factored out here so it
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
        coefficients change. Without this, computeFockMatrixClaude2's (1.0/occ) two-electron scaling uses a
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


## EXPERIMENTAL SWITCH (09-Aug-2026), default false = the behaviour that has always been in place.
## When true, solveOptimizedLevelField builds the Fock matrix as
##      h  +  (1/occ) * V(diagonal CSF pairs)  +  V(off-diagonal CSF pairs)
## instead of  h + (1/occ) * V(all pairs).  This tests whether the off-diagonal (CI-coupling) part should
## be scaled by the generalized occupation at all: it carries weight ~ c_r c_s while occ carries ~ c_r^2,
## so dividing it by occ introduces a c_1/c_2 factor that grows without bound as a correlating CSF's own
## coefficient shrinks -- the suspected mechanism behind the winner-take-all collapse documented in
## solveOptimizedLevelField's KNOWN LIMITATION.  Set with SelfConsistent.setEolUnscaledOffDiagonal(true).
GBL_EOL_UNSCALED_OFFDIAGONAL = false


"""
`SelfConsistent.setEolUnscaledOffDiagonal(flag::Bool)`  
    ... sets the experimental switch that keeps the off-diagonal CSF-pair contributions OUT of the
        (1/occ) scaling in the EOL Fock matrix; nothing is returned.
"""
function setEolUnscaledOffDiagonal(flag::Bool)
    global GBL_EOL_UNSCALED_OFFDIAGONAL = flag
    return( nothing )
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

    coeffs1p = SpinAngular.Coefficient1p[];   coeffs2p = SpinAngular.Coefficient2p[]
    for  (_, (idxCsf, cache1p, cache2p))  in  blockCaches
        n = length(idxCsf)
        for  r = 1:n
            for  s = 1:n
                ## pairs = :all (default, unchanged) | :diagonal (r == s only) | :offdiagonal (r != s only).
                ## The split exists to test whether the off-diagonal CSF-pair contributions should be scaled
                ## by the generalized occupation at all -- see the DA-term note in solveOptimizedLevelField.
                if      pairs == :diagonal      &&  r != s     continue
                elseif  pairs == :offdiagonal   &&  r == s     continue
                end
                drs = 0.
                for  (i, level)  in  enumerate(targetLevels)    drs = drs + weights[i] * level.mc[idxCsf[r]] * level.mc[idxCsf[s]]    end
                if  drs == 0.    continue    end
                for  cf in cache1p[(r,s)]   push!(coeffs1p, SpinAngular.Coefficient1p(cf.nu, cf.a, cf.b, cf.T * drs) )   end
                for  cf in cache2p[(r,s)]   push!(coeffs2p, SpinAngular.Coefficient2p(cf.nu, cf.a, cf.b, cf.c, cf.d, cf.V * drs) )   end
            end
        end
    end

    # Condense angular coefficients if they refer to the same set of orbitals -- identical to
    # SelfConsistent.computeAngularCoefficients' own condensation tail.
    coeffs1px = SpinAngular.Coefficient1p[];     coeffs2px = SpinAngular.Coefficient2p[]

    hasConsidered = falses( length(coeffs1p) );   T = 0.
    for  (ic, cf) in enumerate(coeffs1p)
        if    hasConsidered[ic]
        else  nu = cf.nu;   a = cf.a;   b = cf.b
            T = T + cf.T;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs1p)
                if    hasConsidered[icx]
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cf.b     T = T + cfx.T;    hasConsidered[icx] = true
                end
            end
            push!(coeffs1px, SpinAngular.Coefficient1p(nu, a, b, T) );  T = 0.
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
            push!(coeffs2px, SpinAngular.Coefficient2p(nu, a, b, c, d, V) );   V = 0.
        end
    end

    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.orthonormalizeSameKappaClaude(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
                                              subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                              matrixB::Array{Float64,2})`
    ... applies a Loewdin SYMMETRIC orthogonalization to every group of subshells sharing the same kappa, treating
        them as a SIMULTANEOUSLY-valid set rather than the sequence in which they happen to have been refined --
        this is the explicit "orthonormalize the orbital set" step shown as its OWN, separate box (after, not
        instead of, the per-orbital "apply projection" step) in the block diagram of Zatsarinny \\& Froese
        Fischer, Comput. Phys. Commun. 202, 287-303 (2016), Fig. 1. The per-orbital
        Hamiltonian.projectHamiltonian step alone does not guarantee a simultaneously-consistent same-kappa set
        when subshells are refined one at a time with lagged (Jacobi-style) data from the previous outer
        iteration -- this step corrects that afterwards, once per outer iteration, independent of processing
        order.

        For a kappa-group of n>1 subshells with B-spline expansion coefficient vectors b_1...b_n (columns of a
        (nsL+nsS) x n matrix B) and S-metric overlap matrix M (M_ij = b_i' * matrixB * b_j), the orthonormalized
        set is B_new = B * M^(-1/2), computed via an eigendecomposition of the small (n x n) matrix M. Unlike
        Gram-Schmidt, this treats every member of the group symmetrically -- no member is privileged by
        processing order -- matching the paper's "simultaneous" variation of same-symmetry orbitals directly.
        Groups of size 1 (no same-kappa partner) are left untouched, already trivially orthonormal.
        A (orbitals::Dict{Subshell,Orbital}, bVectors::Dict{Subshell,Vector{Float64}}) tuple is returned.
"""
function orthonormalizeSameKappaClaude(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       matrixB::Array{Float64,2})
    orbitalsOut = deepcopy(newOrbitals)
    bVectorsOut = deepcopy(newbVectors)

    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        group = [sh for sh in subshells if sh.kappa == kappa]
        if  length(group) < 2    continue    end

        B = hcat( [newbVectors[sh] for sh in group]... )        # (nsL+nsS) x n
        M = transpose(B) * matrixB * B                          # n x n, symmetric overlap matrix of the group

        eig  = Basics.diagonalize(MatrixWithLinearAlgebra(), M)
        n    = length(group)
        Mm12 = zeros(n, n)
        for  k = 1:n
            vk    = eig.vectors[k]
            Mm12 .+= (1.0/sqrt(eig.values[k])) .* (vk * transpose(vk))
        end

        Bnew = B * Mm12
        for  (idx, sh)  in  enumerate(group)
            vec             = Bnew[:,idx]
            en              = newOrbitals[sh].energy
            bVectorsOut[sh] = vec
            orbitalsOut[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, en, vec, primitives)
        end
    end

    return( orbitalsOut, bVectorsOut )
end


"""
`SelfConsistent.computeFunctionalClaude(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                        orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)`
    ... computes the MCDHF energy functional using InteractionStrength.XL_CoulombClaude (kink-aware
        two-electron Slater integral) for the two-electron term. Used by the ALField/EOLField code line, cf.
        solveAverageLevelField. An energy::Float64 is returned.
"""
function computeFunctionalClaude(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                 orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)
    energy = 0.

    # Collect one-electron contributions -- unchanged, no kink in this integral
    for  cf  in  coeffs1p
        jj     = Basics.subshell_2j(cf.a)
        energy = energy + cf.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
    end

    # Collect two-electron contributions via the kink-aware integral
    for  cf  in  coeffs2p
        energy = energy + cf.V * InteractionStrength.XL_CoulombClaude(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                                              orbitals[cf.c], orbitals[cf.d], grid)
    end

    return( energy )
end


"""
`SelfConsistent.energyFromBVectorsClaude3(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, grid::Radial.Grid,
        nucPot::Radial.Potential)`  
    ... the EOL energy as a plain scalar function of the orbital B-spline coefficient vectors, with the
        angular coefficients held fixed. This is exactly the functional solveOptimizedLevelField reports,
        just expressed in the variables that are actually varied, so that it can be differentiated. Note
        that it includes coeffs1p -- which the Fock matrix of the present scheme never receives, and which
        is the inconsistency the Claude3 path exists to remove. A value::Float64 is returned.
"""
function energyFromBVectorsClaude3(bVectors::Dict{Subshell, Vector{Float64}},
                                   coeffs1p::Array{SpinAngular.Coefficient1p,1},
                                   coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                   subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                   grid::Radial.Grid, nucPot::Radial.Potential)
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, 0.0, bVectors[sh], primitives)
    end
    return( SelfConsistent.computeFunctionalClaude(coeffs1p, coeffs2p, orbitals, grid, nucPot) )
end


"""
`SelfConsistent.virtualDirectionsClaude3(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, nucPot::Radial.Potential, matrixB::Array{Float64,2},
        storage::Dict{String,Array{Float64,2}}; nVirtual::Int64=20)`  
    ... builds, for each subshell, an S-orthonormal set of ALLOWED rotation directions: positive-branch
        eigenvectors of the one-electron Dirac matrix of that kappa, projected free of every occupied
        orbital of the same kappa. Rotations among the occupied orbitals themselves are deliberately
        excluded -- they leave the CSF space invariant, are redundant, and would make the Hessian singular.
        A Dict{Subshell, Array{Vector{Float64},1}} is returned.
"""
function virtualDirectionsClaude3(bVectors::Dict{Subshell, Vector{Float64}}, subshells::Array{Subshell,1},
                                  primitives::Bsplines.Primitives, nucPot::Radial.Potential,
                                  matrixB::Array{Float64,2}, storage::Dict{String,Array{Float64,2}};
                                  nVirtual::Int64=20)
    virtuals = Dict{Subshell, Array{Vector{Float64},1}}()
    for  sh  in  subshells
        ## the occupied orbitals of this kappa, S-orthonormalized so that the projection below is exact
        occSame = Vector{Float64}[]
        for  s2  in  subshells
            if  s2.kappa != sh.kappa    continue    end
            v = copy(bVectors[s2])
            for  u in occSame    v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            if  nrm > 1.0e-10    push!(occSame, v / nrm)    end
        end
        ## a one-electron reference spectrum for this kappa; orbital-independent, hence a fixed frame
        oneEl = Bsplines.setupLocalMatrix(sh.kappa, primitives, nucPot, storage)
        wc    = Bsplines.diagonalizeLocalMatrix(sh.kappa, oneEl, matrixB, primitives)
        mm    = Bsplines.findPositiveBranchStart(wc.values)
        dirs  = Vector{Float64}[]
        for  i = mm:length(wc.values)
            v = copy(wc.vectors[i])
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
`SelfConsistent.computeOrbitalGradientClaude3(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, nucPot::Radial.Potential,
        storage::Dict{String,Array{Float64,2}})`  
    ... the ANALYTIC gradient dE/db_a of the same energy that energyFromBVectorsClaude3 evaluates, for every
        subshell, as a full B-spline coefficient vector. Nothing here is new machinery: the one-electron
        integral is I(a,b) = b_a^T H1 b_b with H1 = Bsplines.setupLocalMatrix, and the Slater integral is
        R^k(abcd) = b_a^T M b_c with M = InteractionStrength.XL_CoulombClaude(k, a, orb_b, c, orb_d,
        primitives) -- the matrix-valued overload that already exists for the Fock build. Each slot in which
        a subshell occurs contributes once, using the symmetry R^k(abcd) = R^k(badc) for the second pair.
        Must be checked against gradientByFiniteDifferenceClaude3 before being trusted.
        A Dict{Subshell, Vector{Float64}} is returned.
"""
function computeOrbitalGradientClaude3(bVectors::Dict{Subshell, Vector{Float64}},
                                       coeffs1p::Array{SpinAngular.Coefficient1p,1},
                                       coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, 0.0, bVectors[sh], primitives)
    end
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells    grad[sh] = zeros(nsL+nsS)    end

    ## one-electron:  d/db_a [ w * b_a^T H1 b_b ]  contributes to BOTH slots
    h1 = Dict{Int64, Array{Float64,2}}()
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        h1[kappa] = Bsplines.setupLocalMatrix(kappa, primitives, nucPot, storage)
    end
    for  cf  in  coeffs1p
        if  cf.a.kappa != cf.b.kappa    continue    end          ## no one-electron element between kappas
        w = cf.T * sqrt( Basics.subshell_2j(cf.a) + 1 )
        hh = h1[cf.a.kappa]
        grad[cf.a] = grad[cf.a] + w * (hh * bVectors[cf.b])
        grad[cf.b] = grad[cf.b] + w * (transpose(hh) * bVectors[cf.a])
    end

    ## two-electron:  R^k(abcd) = b_a^T M(a,orb_b;c,orb_d) b_c  and, by R^k(abcd) = R^k(badc),
    ##                          = b_b^T M(b,orb_a;d,orb_c) b_d
    for  cf  in  coeffs2p
        mac = InteractionStrength.XL_CoulombClaude(cf.nu, cf.a, orbitals[cf.b], cf.c, orbitals[cf.d], primitives)
        grad[cf.a] = grad[cf.a] + cf.V * (mac * bVectors[cf.c])
        grad[cf.c] = grad[cf.c] + cf.V * (transpose(mac) * bVectors[cf.a])
        mbd = InteractionStrength.XL_CoulombClaude(cf.nu, cf.b, orbitals[cf.a], cf.d, orbitals[cf.c], primitives)
        grad[cf.b] = grad[cf.b] + cf.V * (mbd * bVectors[cf.d])
        grad[cf.d] = grad[cf.d] + cf.V * (transpose(mbd) * bVectors[cf.b])
    end
    return( grad )
end


"""
`SelfConsistent.gradientByFiniteDifferenceClaude3(bVectors::Dict{Subshell, Vector{Float64}},
        virtuals::Dict{Subshell, Array{Vector{Float64},1}}, coeffs1p::Array{SpinAngular.Coefficient1p,1},
        coeffs2p::Array{SpinAngular.Coefficient2p,1}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)`  
    ... the gradient of the EOL energy with respect to the allowed orbital rotations, by central finite
        differences. Slow but free of any derivation, so it is the reference against which an analytic
        gradient must later be checked, and on its own it answers the question "is this converged solution
        actually stationary?". Each direction is S-orthogonal to the orbital itself, so normalization
        contributes only at second order and is omitted. A Dict{Subshell, Vector{Float64}} is returned.
"""
function gradientByFiniteDifferenceClaude3(bVectors::Dict{Subshell, Vector{Float64}},
                                           virtuals::Dict{Subshell, Array{Vector{Float64},1}},
                                           coeffs1p::Array{SpinAngular.Coefficient1p,1},
                                           coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                           subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                           grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)
    grad = Dict{Subshell, Vector{Float64}}()
    for  sh  in  subshells
        dirs = virtuals[sh];    g = zeros( length(dirs) )
        for  (iv, phi)  in  enumerate(dirs)
            bPlus  = copy(bVectors);    bPlus[sh]  = bVectors[sh] + hStep * phi
            bMinus = copy(bVectors);    bMinus[sh] = bVectors[sh] - hStep * phi
            ePlus  = SelfConsistent.energyFromBVectorsClaude3(bPlus,  coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            eMinus = SelfConsistent.energyFromBVectorsClaude3(bMinus, coeffs1p, coeffs2p, subshells, primitives, grid, nucPot)
            g[iv]  = (ePlus - eMinus) / (2 * hStep)
        end
        grad[sh] = g
    end
    return( grad )
end


"""
`SelfConsistent.solveOptimizedLevelFieldClaude3(basis::Basis, nuclearModel::Nuclear.Model,
        primitives::Bsplines.Primitives, settings::AsfSettings; printout::Bool=true)`  
    ... EOL by direct minimization of the energy over orbital ROTATIONS, instead of by solving a per-subshell
        generalized eigenvalue problem. Built beside solveOptimizedLevelField, which is left untouched and is
        still what Basics.EOLField dispatches to; this one is called explicitly.

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
"""
function solveOptimizedLevelFieldClaude3(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true, nVirtual::Int64=16)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    nucPot  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)

    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitivesClaude(basis.orbitals[sh], primitives, matrixB)
    end

    if  settings.levelSelectionCI.active  &&  !isempty(settings.levelSelectionCI.symmetries)
        relevantSyms = unique( settings.levelSelectionCI.symmetries )
    else
        relevantSyms = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    end
    blockCaches = Dict()
    for  sym  in  relevantSyms    blockCaches[sym] = SelfConsistent.cacheCsfPairCoefficientsEOL(sym, basis)   end

    ePrevious = 0.;   tStep = 1.0;   multiplet = Multiplet("EOL-Claude3", Level[])
    for  iter = 1:settings.maxIterationsScf
        orbitals = Dict{Subshell, Orbital}()
        for  sh  in  basis.subshells
            orbitals[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, 0.0, bVectors[sh], primitives)
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
        multiplet    = Basics.sortByEnergy( Multiplet("EOL-Claude3", levels) )
        targetLevels = SelfConsistent.selectTargetLevelsEOL(multiplet, settings.levelSelectionCI)
        (coeffs1p, coeffs2p) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)

        e0   = SelfConsistent.energyFromBVectorsClaude3(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                         primitives, grid, nucPot)
        grad = SelfConsistent.computeOrbitalGradientClaude3(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                             primitives, nucPot, storage)
        virt = SelfConsistent.virtualDirectionsClaude3(bVectors, basis.subshells, primitives, nucPot,
                                                        matrixB, storage; nVirtual=nVirtual)
        ## Project the gradient on the allowed rotations, and PRECONDITION each component by the
        ## orbital-energy denominator (eps_v - eps_a) -- the diagonal of the rotation Hessian, i.e. the
        ## standard first-order estimate  kappa_av = -g_av / (eps_v - eps_a).  Plain steepest descent
        ## converges far too slowly here: it left the Li control 1.8e-6 Ha short after 60 steps.
        ## The denominator is floored, since a near-degenerate pair would otherwise produce a huge step.
        gProj = Dict{Subshell, Vector{Float64}}();   step = Dict{Subshell, Vector{Float64}}()
        gNorm = 0.;    sNorm = 0.
        for  sh  in  basis.subshells
            h1k  = Bsplines.setupLocalMatrix(sh.kappa, primitives, nucPot, storage)
            epsA = transpose(bVectors[sh]) * h1k * bVectors[sh]
            gv   = [ transpose(phi) * grad[sh]  for phi in virt[sh] ]
            sv   = zeros( length(gv) )
            for  (iv, phi)  in  enumerate(virt[sh])
                epsV    = transpose(phi) * h1k * phi
                sv[iv]  = - gv[iv] / max( epsV - epsA, 0.05 )
            end
            gProj[sh] = gv;    step[sh] = sv
            gNorm = gNorm + sum( gv.^2 );    sNorm = sNorm + sum( sv.^2 )
        end
        gNorm = sqrt(gNorm);    sNorm = sqrt(sNorm)
        if  sNorm < 1.0e-14    break    end
        if  printout
            println(">> [EOL-C3] iter $iter:  E = $(multiplet.levels[1].energy)   |grad| = $gNorm   step = $tStep")
        end
        if  gNorm < 1.0e-8    break    end

        ## backtracking line search along -gProj; halve until the energy actually falls
        accepted = false
        for  trial = 1:24
            newB = Dict{Subshell, Vector{Float64}}()
            for  sh  in  basis.subshells
                v = copy(bVectors[sh])
                for  (iv, phi)  in  enumerate(virt[sh])    v = v + tStep * step[sh][iv] * phi    end
                newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
            end
            eTrial = SelfConsistent.energyFromBVectorsClaude3(newB, coeffs1p, coeffs2p, basis.subshells,
                                                               primitives, grid, nucPot)
            if  eTrial < e0
                newOrbs = Dict{Subshell, Orbital}()
                for  sh in basis.subshells
                    newOrbs[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, 0.0, newB[sh], primitives)
                end
                (_, bVectors) = SelfConsistent.orthonormalizeSameKappaClaude(newOrbs, newB, basis.subshells,
                                                                             primitives, matrixB)
                accepted = true;    tStep = min(1.0, 1.3*tStep);    break
            end
            tStep = tStep / 2
        end
        if  !accepted
            if  printout    println(">> [EOL-C3] no descent found; stopping at iteration $iter.")    end
            break
        end
        if  iter > 1  &&  abs(multiplet.levels[1].energy - ePrevious) < settings.accuracyScf    break    end
        ePrevious = multiplet.levels[1].energy
    end
    return( multiplet )
end


"""
`SelfConsistent.initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Initialized a many-electron basis from the given list of configurations, the nuclear model as well as ASF settings.
        It assumes that a proper set of primitives::Primitives has been initialized before. The initial set of orbitals in this
        basis is determined by the settings::AsfSettings.  A basis::Basis is returned.
"""
function initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                         settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    NoElectrons = configs[1].NoElectrons;   subshells = Subshell[];   coreSubshells = Subshell[];     csfs = CsfR[] 
    orbitals    = Dict{Subshell, Orbital}()
    
    # Perform some simple tests: Number of electrons must be equal in all configurations
    for  conf in configs   if  conf.NoElectrons != NoElectrons    error("stop a")   end     end
    
    # Generate a full set of relativistic CSF from the given configurations and collect the associated level symmetries
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">>> include ", relconfList[i])    end   end
    subshells = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshells; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshells)
        append!( csfList, newCsfs)
    end
    
    # Select CSF with requested symmetry if needed
    if  length(levelSymmetries) == 0
        csfs = csfList          # Take all relativistic CSF into account
    else
        for  csf in csfList
            if  LevelSymmetry(csf.J, csf.parity)  in  levelSymmetries   push!(csfs, csf)    end
        end
    end

    # Determine the number of electrons and the list of coreSubshells
    for  (k,sh)  in  enumerate(subshells)
        mocc = Basics.subshell_2j(sh) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc     is_filled = false;           break   end
        end
        if   is_filled    push!( coreSubshells, subshells[k])      else    break   end
    end
        
    # Check that the radial grid is able to represent these subshells at all, BEFORE any orbital is generated.
    # This sits here rather than inside Bsplines.generateOrbitalsHydrogenic so that it applies whichever way the
    # orbitals are seeded -- StartFromPrevious inherits a grid just as much as StartFromHydrogenic does.
    Bsplines.checkGridRepresentation(subshells, nuclearModel.Z, primitives;
                                     accuracy = settings.gridAccuracy, stopper = settings.gridStopper)

    # Initialize the orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("> Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = Bsplines.generateOrbitalsHydrogenic(subshells, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("> Start SCF process from given list of orbitals.energy")    end
        # Taking starting orbitals for the given dictionary; non-relativistic orbitals with a proper nuclear charge
        # are adapted if no orbital is found
        orbitals = Dict{Subshell, Orbital}()
        for  sh in subshells
            if  haskey(settings.startScfFrom.orbitals, sh)  
                orbitals[sh] = settings.startScfFrom.orbitals[sh]
            else
                println("Start orbitals do not contain an Orbital for subshell $sh ")
                orb          = HydrogenicIon.radialOrbital(subsh, nuclearModel.Z, grid)
                orb          = Orbital(orb.subshell, orb.isBound, true, orb.energy, orb.P, orb.Q, orb.Pprime, orb.Qprime, Radial.Grid())
                orbitals[sh] = orb 
            end
        end
    else  error("stop b")
    end
    
    basis = Basis(true, NoElectrons, subshells, csfs, coreSubshells, orbitals)
    return( basis )
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    
    # Generate primitives and initialize the many-electron basis
    Defaults.setDefaults("standard grid", grid)
    primitives = Bsplines.generatePrimitives(grid)    
    basis      = SelfConsistent.initializeBasis(configs, nm, primitives, settings; levelSymmetries, printout)
    
    # Solve a self-consistent field for this basis
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()]
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   typeof(settings.scField) == Basics.EOLField
        # solveOptimizedLevelField already returns a complete, correctly-diagonalized multiplet (built via
        # its own internal, kink-aware Hamiltonian.performCIClaude call on the converged orbitals) -- return
        # it directly. Re-diagonalizing below would be redundant AND wrong: EOLField isn't in the
        # ALField check just below, so it would fall through to the bare,
        # non-kink-aware Hamiltonian.performCI, silently discarding the kink-aware result.
        return( SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   settings.scField in [Basics.ALField()]
        mp = Hamiltonian.performCIClaude(basis, nm, grid, settings, printout=printout)
    else
        mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)
    end

    return( mp )
end


"""
`SelfConsistent.performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)`
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)
    
    # Generate primitives
    primitives = Bsplines.generatePrimitives(grid)    
    
    # Solve a self-consistent field for this basis
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()]
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   typeof(settings.scField) == Basics.EOLField
        # See the identical note in the other performSCF overload just above: solveOptimizedLevelField
        # already returns a complete, correctly (kink-aware) diagonalized multiplet -- return it directly
        # rather than redundantly re-diagonalizing with the wrong, non-kink-aware Hamiltonian.performCI.
        return( SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   settings.scField in [Basics.ALField()]
        mp = Hamiltonian.performCIClaude(basis, nm, grid, settings, printout=printout)
    else
        mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)
    end

    return( mp )
end


"""
`SelfConsistent.rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid,
                               settings::AsfSettings)`
    ... rotates pairwise the orbitals to enhance the covergence of the SCF procedures; it follows that 
        rotation schemes of C.F. Fischer ... . The pairwise rotation of the orbitals of the same symmetry 
        is done from the inner-to-outer subshells. A rotation is not necessary if two subshells are filled
        of if both are fixed. For each rotation, a rotation parameter beta is determined.
        A pair (newOrbitals::Dict{Subshell, Orbital}, betaMax::Float64) us returned which contains the rotated
        orbitals as well as the maximum beta parameter that occurred during the pairwise rotation of the orbital set.
"""
function rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, 
                        settings::AsfSettings)
    betaMax = 0.;    newOrbitals = orbitals
    
    @warn("Not yet properly implemented.")
    
    return( newOrbitals, betaMax)
end


"""
`SelfConsistent.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64,
                                           nm::Nuclear.Model, grid::Radial.Grid)`
    ... determines the chemical potential so that Sum_i f(epsilon_i, mu, temp) = Z.
        The Newton-Raphson methods is used to iterate to the chemical potential; a chemMu::Float64 is returned.

        Note: this general finite-temperature Fermi-Dirac root-finding utility was moved here from module Plasma
              (where it originated as Plasma.determineChemicalPotential), since Plasma.perform(::AverageAtomScheme,
              ...) needs SelfConsistent.solveAverageAtomField below, and solveAverageAtomField itself needs this
              function internally at every SCF iteration; keeping it in Plasma would have made the two modules
              depend on each other circularly. Nothing here is Plasma-scheme-specific.
"""
function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)
    function g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = - nm.Z
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            wa  = wa + occ / (exp(wb) + 1)
        end
        return( wa )
    end
    #
    function gprime(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = 0.
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            wc  = exp( wb )
            wa  = wa + occ * wc^2 / temp / (wc+1)^2
        end
        return( wa )
    end
    # Iterate for the chemical potential
    chemMu = -0.1;     newMu = 0.;     nx = 0
    while true
        nx = nx + 1
        newMu = chemMu - g(chemMu, orbitals, temp, nm) / gprime(chemMu, orbitals, temp, nm)
        if  abs(newMu - chemMu) < 1.0e-4  break
        else    chemMu = newMu
        end
    end

    chemMu = chemMu - 0.0011  ## Seems to bring better stability in the SCF computations

    println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
    return ( chemMu )
end


"""
`SelfConsistent.solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                                      temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)`
    ... solves the self-consistent field for a given local average-atom potential as specified by scField
        A (new) set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                               temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)
    # Determine the chemical potential
    chemMu    = determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, primitives.grid);
    
    # Extract the kappa's from orbitals
    kappas = Int64[];     for (k,v)  in  orbitals     push!(kappas, k.kappa)    end;    kappas = unique(kappas);

    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, grid)
            
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  32
                println(">> Maximum number of SCF iterations = 32 is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) Re-compute the local potential
            wp  = Basics.computePotential(scField, grid, previousOrbitals, chemMu, temp)
            pot = Basics.add(nuclearPotential, wp)
            
            # (2) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (3) Solve the generalized eigenvalue problem
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
            
            # (4) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.0e-6   go_on = true   end     ## accuracyScf
            for  (k,v)  in  orbitals
                if      k.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(k, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[k], newOrbital)
                    if  wcOrbital > 1.0e-6   accuracyScf = wcOrbital;   go_on = true   end     ## accuracyScf
                        sa = "  $k::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    ## println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
                    previousOrbitals[k] = newOrbital
                end
            end
            # (5) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        chemMu              = determineChemicalPotential(previousOrbitals, temp, radiusWS, nuclearModel, primitives.grid)
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)    
    return( analyzedOrbitals )
end


"""
`SelfConsistent.computeTwoElectronVClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                           bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                           tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})`
    ... computes the direct and exchange two-electron potential matrix for the given subshell, taking
        bVectors::Dict{Subshell,Vector{Float64}} -- B-spline expansion coefficients -- as the SOLE, canonical
        per-subshell state, rather than a persistent Dict{Subshell,Orbital}. Every coeffs2p entry relevant to
        subshell involves exactly one "partner" subshell (the diagonal self-term is its own partner); a
        partner's tabulated form is built only as a disposable, read-only byproduct
        (Bsplines.generateOrbitalFromVectorClaude), cached within this call so repeated coefficients sharing
        the same partner do not re-evaluate it, and discarded once this function returns -- never re-fit back
        into bVectors, since bVectors are never derived FROM a tabulated form here in the first place; they
        come directly from diagonalization (solveAverageLevelField/solveOptimizedLevelField). Shared
        by the ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is
        returned.
"""
function computeTwoElectronVClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                    bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                    tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})
    nsL     = primitives.grid.nsL;        nsS = primitives.grid.nsS
    matrixV = zeros( nsL+nsS, nsL+nsS )
    partnerOrbitals = Dict{Subshell, Orbital}()
    partnerCVectors = Dict{Subshell, Vector{Float64}}()

    for  cf  in  coeffs2p
        if      subshell == cf.a  &&  subshell == cf.b  &&  subshell == cf.c  &&  subshell == cf.d
            partner = subshell;   pattern = :diagonal
        elseif  subshell == cf.a  &&  subshell == cf.c  &&  cf.b == cf.d
            partner = cf.b;       pattern = :direct
        elseif  subshell == cf.a  &&  subshell == cf.d  &&  cf.b == cf.c
            partner = cf.b;       pattern = :exchange
        elseif  subshell == cf.b  &&  subshell == cf.d  &&  cf.a == cf.c
            partner = cf.a;       pattern = :direct
        elseif  subshell == cf.b  &&  subshell == cf.c  &&  cf.a == cf.d
            partner = cf.a;       pattern = :exchange
        elseif  subshell == cf.a  &&  subshell == cf.b  &&  cf.c == cf.d  &&  cf.c != subshell
            # "Pair-excitation" pattern (a=b=subshell, c=d=partner, subshell!=partner): arises only from
            # genuine cross-CSF (off-diagonal) two-particle matrix elements, e.g. an EOL coefficient
            # combination mixing two different configurations (never produced by AL's single-CSF-average
            # scheme, hence never needed before). Algebraically IDENTICAL to the :exchange pattern above --
            # R^k(a,b,c,d) with (a,b,c,d)=(X,X,Y,Y) equals R^k(X,Y,Y,X) exactly, since the second
            # integration variable's density is just a product of two scalar functions (Y(r)*X(r) times
            # X(r)*Y(r) are the same product, whichever argument slot each factor nominally occupies) --
            # so it is handled by the identical XL_CoulombTensorClaude call, unweighted, as :exchange.
            partner = cf.c;       pattern = :exchange
        elseif  subshell == cf.c  &&  subshell == cf.d  &&  cf.a == cf.b  &&  cf.a != subshell
            partner = cf.a;       pattern = :exchange
        else
            continue
        end
        # NOTE (26-Jul-2026): a previous `iseven(subshell_l + partner_l + cf.nu)` filter here, intended to
        # drop Breit-only coefficients, was REMOVED -- it used the parity of two unrelated orbitals' l-values
        # plus nu, which is not a real Breit/Coulomb selection rule and was silently discarding legitimate
        # direct (nu=0) Coulomb terms between different-l shells (e.g. 1s-2p, 2s-2p). This is why Ne/Ar were
        # wrong while He/Be (pure s-shells, where the bogus filter never triggered) validated fine. coeffs2p
        # here comes from Basics.ALField()'s own (Coulomb-only) angular coefficients, so no Breit terms are
        # present to filter in the first place. Confirmed: summing ALL entries reproduces Ne's literature
        # total energy to ~3e-3 Ha (residual consistent with import precision), vs -157.6 Ha when filtered.

        partnerOrb = get!(partnerOrbitals, partner) do
            Bsplines.generateOrbitalFromVectorClaude(partner, 0.0, bVectors[partner], primitives)
        end

        if      pattern == :diagonal
            matrixV = matrixV + 2 * cf.V *
                      InteractionStrength.XL_CoulombClaude(cf.nu, subshell, partnerOrb, subshell, partnerOrb, primitives)
        elseif  pattern == :direct
            matrixV = matrixV + cf.V *
                      InteractionStrength.XL_CoulombClaude(cf.nu, subshell, partnerOrb, subshell, partnerOrb, primitives)
        else    # :exchange
            # partnerOrb was built by Bsplines.generateOrbitalFromVectorClaude, which silently canonicalizes its
            # sign so that P is positive at small r (see that function's `wSign` step) -- a convention applied
            # to the TABULATED reconstruction only, never fed back into bVectors[partner] itself. Whenever that
            # canonicalization actually flips the sign, partnerOrb and the raw bVectors[partner] end up
            # oppositely signed, even though XL_CoulombTensorClaude's "b"/"c" arguments (from partnerOrb) and its
            # "cVector" argument (from bVectors[partner]) are meant to represent the identical orbital -- an odd
            # total power of partner's sign then leaks into the exchange matrix element (confirmed empirically
            # this session: exchange terms flip sign in exact lockstep with this mismatch, while diagonal/direct
            # terms -- which only ever use partnerOrb, an even number of times -- stay sign-invariant as they
            # should). Fix: mirror generateOrbitalFromVectorClaude's own small-r wSign test on the raw bVector,
            # and pass a sign-matched cVector so it always agrees with partnerOrb's canonical sign.
            cVector = get!(partnerCVectors, partner) do
                bV      = bVectors[partner]
                wSign   = 0.0
                for  i = 1:nsL
                    Bi    = primitives.bsplinesL[i];   add = 1 - Bi.lower
                    for  j = Bi.lower:min(Bi.upper,30)   wSign = wSign + bV[i] * Bi.bs[j+add]   end
                end
                wSign < 0.  ?  -bV  :  bV
            end
            (cacheLL, cacheLS, cacheSS) = tensorCaches[cf.nu]
            matrixV = matrixV + cf.V *
                      InteractionStrength.XL_CoulombTensorClaude(cf.nu, subshell, partnerOrb, partnerOrb,
                                                                  cVector, subshell,
                                                                  cacheLL, cacheLS, cacheSS, primitives)
        end
    end

    return( matrixV )
end


"""
`SelfConsistent.computeFockMatrixClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                         bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                         nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                         occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})`
    ... computes the full per-subshell Hamiltonian matrix (one-electron + direct + exchange), matching the
        role of DBSR_HF's hf_matrix.f90: Bsplines.setupLocalMatrix provides the one-electron (kinetic +
        nuclear + rest-mass) block natively in B-spline basis (reused unchanged), and
        computeTwoElectronVClaude2 provides the two-electron block, divided by the subshell's own mean
        occupation -- algebraically identical to DBSR_HF's own hfm=qsum(i)*dhl (one-electron scaled up)
        followed by hfm=hfm/qsum(i) (whole matrix scaled down) once expanded, given the two-electron
        coefficients already carry their full (not per-electron) occupation weighting (verified this
        session against the DBSR_HF av_energy_coef reference for all of Ne's coefficients). Shared by the
        ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrix::Array{Float64,2} is returned.
"""
function computeFockMatrixClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                  bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                  nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                  occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}};
                                  coeffs2pUnscaled::Array{SpinAngular.Coefficient2p,1}=SpinAngular.Coefficient2p[])
    # occ == 0 would give 1/occ = Inf and, against the zero two-electron matrix such a subshell has,
    # Inf * 0 = NaN in every element -- a whole Fock matrix of NaN that only surfaces much later, and
    # then as a completely unrelated-looking complaint about the negative-energy continuum.  A subshell
    # with no occupation has no mean field to be refined in; the caller must decide what to do with it
    # (solveOptimizedLevelField carries it forward unchanged), so refuse here rather than return NaN.
    if  abs(occ) < 1.0e-12
        error("SelfConsistent.computeFockMatrixClaude2(): subshell $subshell has occupation $occ. There is no " *
              "mean field to define for an unoccupied subshell; the caller must skip it instead.")
    end
    matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
    matrixV = computeTwoElectronVClaude2(subshell, coeffs2p, bVectors, primitives, tensorCaches)
    ## coeffs2pUnscaled, when given, contributes WITHOUT the 1/occ scaling. Empty by default, so the
    ## returned matrix is bit-for-bit what it always was unless a caller opts in.
    if  length(coeffs2pUnscaled) > 0
        matrixU = computeTwoElectronVClaude2(subshell, coeffs2pUnscaled, bVectors, primitives, tensorCaches)
        return( matrix + (1.0/occ) * matrixV + matrixU )
    end
    return( matrix + (1.0/occ) * matrixV )
end


"""
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                       settings::AsfSettings; printout::Bool=true)`
    ... solves the average-level (AL) self-consistent field, validated to 5+ significant figures against
        literature for He/Be/Ne/Ar (see project memory), with a bVector-native architecture modeled directly
        on DBSR_HF (Zatsarinny & Froese Fischer, CPC 202, 287 (2016)): B-spline expansion coefficients
        (Dict{Subshell,Vector{Float64}}) are the SOLE, canonical, persistent per-iteration orbital state --
        no Dict{Subshell,Orbital} is maintained across iterations at all, and there is no
        Bsplines.fitVectorToPrimitivesClaude "fit back" step, since bVectors are never derived from a
        tabulated form to begin with; they come directly from Basics.diagonalize's eigenvector for each
        subshell, in turn. Orthogonality between same-kappa subshells (e.g. Ne's 1s/2s) is enforced by
        projecting the Fock matrix (Hamiltonian.projectHamiltonian, reused unchanged) against each
        ALREADY-PROCESSED lower same-kappa subshell's bVector directly inside the generalized eigenvalue
        problem, before diagonalizing -- matching DBSR_HF's hf_solve_HF.f90/hf_eiv sequential approach --
        rather than a post-hoc Löwdin symmetric orthogonalization of the whole
        same-kappa group after the per-orbital loop. The target eigenvalue index is shifted down by one for
        each such projection applied, mirroring hf_eiv's `mm = m + (orthogonalized-count) - 1`. A tabulated
        Orbital is only ever built as a disposable, read-only byproduct: once per unique "partner" subshell
        inside computeTwoElectronVClaude2 (for the two-electron potential contraction), and once per
        subshell per iteration purely for reporting/energy-functional evaluation (reusing
        computeFunctionalClaude unchanged) -- never stored as competing state, never refit. A final,
        single export pass (Bsplines.generateOrbitalFromVectorClaude) produces a standard
        Dict{Subshell,Orbital} for the returned basis::Basis, so every downstream consumer (properties,
        processes, CI/DCB Hamiltonian construction) is unaffected by this being a bVector-native SCF.
        Reached via performSCF's scField = Basics.ALField() dispatch. A (new) basis::Basis is returned.
"""
function solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                settings::AsfSettings; printout::Bool=true)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid

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
        bVectors[sh] = Bsplines.fitVectorToPrimitivesClaude(basis.orbitals[sh], primitives, matrixB)
        energies[sh] = basis.orbitals[sh].energy
    end

    # (2) Generate angular coefficients (unaffected by the kink fix / bVector-native rebuild)
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(Basics.ALField(), basis)

    # (3) Precompute kink-aware Slater-moment tensor caches for every rank that occurs; only the exchange
    # branches of computeTwoElectronVClaude2 use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [AL] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    orbitals = Dict{Subshell, Orbital}()    # only ever a disposable, per-iteration reporting byproduct

    for  iter = 1:settings.maxIterationsScf
        println("\n> SCF interation $(iter) [AL]: ")
        newBVectors = Dict{Subshell, Vector{Float64}}();    newEnergies = Dict{Subshell, Float64}()
        processedBVectors = Dict{Subshell, Vector{Float64}}()
        dpm = Dict{Subshell, Float64}()

        for  subshell  in  basis.subshells
            occ = meanOcc[subshell]
            print(">> Refine $subshell orbital with mean occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrixClaude2(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches)

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
            newOrbitals[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, newEnergies[sh], newBVectors[sh], primitives)
        end
        eFunctional = SelfConsistent.computeFunctionalClaude(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        orbitalConv = maximum( values(dpm) ) < 1.0 ? 1.0 - maximum( values(dpm) ) : 0.0

        println(">> Total energy = $(eFunctional*1)   orbital-conv = $orbitalConv   orbital-acc = $(1.0 - orbitalConv)")

        bVectors = newBVectors;    energies = newEnergies;    orbitals = newOrbitals
        if  abs(1.0 - orbitalConv) < settings.accuracyScf    break   end
    end

    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end


function normX(sa::String, matrix::Array{Float64,2})
    absD = absND = 0.
    for  i = 1:size(matrix,1)                    absD  = absD  + abs(matrix[i,i])   
        for  j = 1:size(matrix,2)   if i != j    absND = absND + abs(matrix[i,j])  end   end
    end 
    return(nothing)
end


"""
`SelfConsistent.solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                    settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        the (local) mean-field potential as specified by the settings::AsfSettings. A (new) basis::Basis is returned.
"""
function solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                             settings::AsfSettings; printout::Bool=true)
    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)

    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) First re-define an (arbitrary) 'level' that represents the mean occupation for the local potential
            wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
            NoCsf  = length(wBasis.csfs)
            wmc    = zeros( NoCsf );   wN = 0.
            for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
            for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
            wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
            # (2) Re-compute the local potential
            wp  = Basics.computePotential(settings.scField, primitives.grid, wLevel)
            pot = Basics.add(nuclearPotential, wp)
            # (3) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (4) Solve the generalized eigenvalue problem
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
            # (5) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
            for  sh in basis.subshells
                if      sh in settings.frozenSubshells   ## do nothing
                elseif  sh.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(sh, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                    if  wcOrbital > settings.accuracyScf   accuracyScf = wcOrbital;   go_on = true   end
                        sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    previousOrbitals[sh] = newOrbital
                end
            end
            # (6) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, analyzedOrbitals)
    return( newBasis )
end
    

"""
`SelfConsistent.solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true)`
    ... solves the self-consistent field for the extended-optimal-level (EOL) functional: orbitals are
        optimized against a statistically-(2J+1)-weighted combination of one or more target ASF levels,
        with the combination weights coming from the levels' own CI mixing coefficients -- never a
        user-supplied weight. The target level(s) are selected EXCLUSIVELY via settings.levelSelectionCI
        (either indices or symmetries, never both; if inactive/empty, defaults to indices=[1], a genuine
        OL/single-level computation). This nests a CI diagonalization inside the AL/Claude2 SCF loop
        (GRASP rmcdhf90's scf.f90, algorithm 5.1 of Froese Fischer, Comput. Phys. Rep. 3 (1986) 290):
        diagonalize -> build generalized coefficients from the target levels' mixing vectors -> refine
        every orbital (reusing computeFockMatrixClaude2/computeTwoElectronVClaude2 unchanged) ->
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
        (CSF-pair) coefficients folded into computeTwoElectronVClaude2's two-electron potential scale
        LINEARLY in a shrinking CSF's own mixing coefficient, while computeGeneralizedOccupationEOL's
        occupation (the (1.0/occ) divisor in computeFockMatrixClaude2) scales QUADRATICALLY in it -- so the
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
        bVectors[sh] = Bsplines.fitVectorToPrimitivesClaude(basis.orbitals[sh], primitives, matrixB)
    end
    orbitals = basis.orbitals

    # (2) Cache the (orbital-independent) per-CSF-pair angular coefficients once for every relevant block
    if  printout    println(">> [EOL] Caching per-CSF-pair angular coefficients for symmetries $(relevantSyms) ...")    end
    blockCaches = Dict{LevelSymmetry, Tuple{Array{Int64,1}, Dict{Tuple{Int64,Int64},Array{SpinAngular.Coefficient1p,1}},
                                             Dict{Tuple{Int64,Int64},Array{SpinAngular.Coefficient2p,1}}}}()
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
        coeffs2pOff                = SpinAngular.Coefficient2p[]
    end
    genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels, basis)

    weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
        sumW = sum( twiceJp1(level.J)  for level in targetLevels )
        [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    end
    weightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

    # (4) Precompute kink-aware Slater-moment tensor caches for every rank that occurs, exactly as
    # solveAverageLevelField does; only the exchange branches of computeTwoElectronVClaude2 use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [EOL] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
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
            # Refining it regardless used to produce a Fock matrix of NaN (computeFockMatrixClaude2 divides by
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

            matrix = SelfConsistent.computeFockMatrixClaude2(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches; coeffs2pUnscaled=coeffs2pOff)

            count = Base.count( sh2 -> sh2.kappa == subshell.kappa, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Bsplines.diagonalizeLocalMatrix(subshell.kappa, matrix, matrixB, primitives)
            l  = Basics.subshell_l(subshell)
            mm = Bsplines.findPositiveBranchStart(wc.values)
            oldVector = bVectors[subshell]
            ## TESTED AND REFUTED (09-Aug-2026): selecting the eigenvector of maximum OVERLAP with the
            ## previous orbital instead of the counted index changes nothing here (identical to eight
            ## decimals), so the counted index was never the problem -- and it actively harms a correlating
            ## orbital, which must be free to change character.  Do not re-propose it.
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
            newOrbitals[sh] = Bsplines.generateOrbitalFromVectorClaude(sh, 0.0, bVectors[sh], primitives)
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
            coeffs2pOff                = SpinAngular.Coefficient2p[]
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
    multiplet  = Hamiltonian.performCIClaude(finalBasis, nuclearModel, grid, settings; printout=printout)
    return( multiplet )
end

end # module
