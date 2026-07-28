
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
    n = length(idxCsf);   matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = 1:n
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
function combineAngularCoefficientsEOL(blockCaches, targetLevels::Array{Level,1})
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumWeights  = sum( twiceJp1(level.J)  for level in targetLevels )
    weights     = [ twiceJp1(level.J) / sumWeights  for level in targetLevels ]

    coeffs1p = SpinAngular.Coefficient1p[];   coeffs2p = SpinAngular.Coefficient2p[]
    for  (_, (idxCsf, cache1p, cache2p))  in  blockCaches
        n = length(idxCsf)
        for  r = 1:n
            for  s = 1:n
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
`SelfConsistent.computeConvergence(aOrbitals::Dict{Subshell, Orbital}, bOrbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)` 
    ... computes a measure for the orbital convergence, i.e. the convergence of all generated orbitals.
        Such a measure is given by the modulus of all products of overlap integrals 
        vonv = | Prod_a   <a|a>  =!= 1, if all orbitals are the same. The procedure assumes that all orbitals
        in aOrbitals and bOrbitals are normalized. The procedure terminates with an error message if any
        orbital in aOrbitals is not found also in aOrbitals. A value conv::Float64 is returned.
"""
function computeConvergence(aOrbitals::Dict{Subshell, Orbital}, bOrbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)
    conv = 1.
    
    for (asubsh, aorb)  in  aOrbitals
        borb = bOrbitals[asubsh]   # Deal with situations that no orbital is found
        conv = conv * RadialIntegrals.overlap(aorb, borb, grid)
    end
    
    return( abs(conv) )
end


"""
`SelfConsistent.computeConvergence(aVectors::Dict{Subshell, Vector{Float64}}, bVectors::Dict{Subshell, Vector{Float64}},
                                   matrixB::Array{Float64,2})` 
    ... computes a measure for the orbital convergence due to the B-spline expansion coefficients, i.e. the convergence 
        of all generated orbitals in terms of their expansion coefficients. Such a measure is given by the modulus 
        of all products of overlap integrals 
        
            conv = | Prod_a   aVector[a] * B * bVector[a] | =!= 1, if all orbitals are the same, and where 
        
        matrixB contains the overlap matrix of the primitives. Formally, this should be always equivalent to
        convergence criterium that is directly based on the orbital representation and, therefore, enables an
        independent test of the generation of orbitals. The procedure assumes that all coefficient vectors 
        are properly normalized. The procedure terminates with an error message if any vector for an orbital 
        in aVectors is not found also in bVectors. A value conv::Float64 is returned.
"""
function computeConvergence(aVectors::Dict{Subshell, Vector{Float64}}, bVectors::Dict{Subshell, Vector{Float64}},
                            matrixB::Array{Float64,2})
    conv = 1.
    
    for (asubsh, aVector)  in  aVectors
        bVector = bVectors[asubsh]   # Deal with situations that no orbital is found
        conv    = conv * (transpose(aVector) * matrixB * bVector)
    end

    return( abs(conv) )
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
`SelfConsistent.computeDirectExchangeV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                       primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital})`
    ... computes the direct and exchange contributions to the one-electron Hamiltonian matrix of the given subshelll.
        These contributions and their position in the Hamiltonian matrix can be derived from the position of subshell
        in the individual coefficient. The coefficient need to be "doubled" if the interaction refer to the same
        subshell due to the variation of the functional. Likely, moreover, the coefficient need to be multiplied
        by sqrt(2j + 1) where j refers to subshell (??)
        A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function computeDirectExchangeV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                                primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital})
    nsL     = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    matrixV = zeros( nsL+nsS, nsL+nsS )
    
    for  cf  in  coeffs2p
        ## @show cf
        if      subshell == cf.a  &&  subshell == cf.b  &&  subshell == cf.c  &&  subshell == cf.d  
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + 2 * cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.c  &&  cf.b == cf.d                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
             else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.d  &&  cf.b == cf.c                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu ) 
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], orbitals[cf.c], subshell, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.d  &&  cf.a == cf.c                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.a], subshell, orbitals[cf.c], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.c  &&  cf.a == cf.d                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu ) 
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.a], orbitals[cf.d], subshell, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell != cf.a  &&  subshell != cf.b  &&  subshell != cf.c  &&  subshell != cf.d  continue
        else    error("No proper interpretation of coeff = $cf")
        end
    end

    return( matrixV )
end


"""
`SelfConsistent.computeDirectExchangeVClaude(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                             primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital},
                                             bVectors::Dict{Subshell, Vector{Float64}},
                                             tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})`
    ... computes the same direct and exchange contributions to the one-electron Hamiltonian matrix as
        computeDirectExchangeV, but using InteractionStrength.XL_CoulombClaude (kink-aware, per-call) wherever
        the refined subshell occupies the (a,c) "direct" position of the coefficient, and
        InteractionStrength.XL_CoulombTensorClaude (kink-aware, precomputed-tensor-based, no per-call
        quadrature) wherever it occupies the (a,d)/(b,c) "exchange" position. tensorCaches must hold, for every
        rank cf.nu that occurs in coeffs2p, a (cacheLL,cacheLS,cacheSS) triple as built by
        RadialIntegrals.buildSlaterMomentCacheClaude -- see solveAverageLevelFieldClaude, which builds this once
        per SCF run, before iterating. bVectors must hold the current B-spline expansion coefficients for every
        subshell (as obtained from Bsplines.extractVectorFromPrimitives), needed by the exchange branches (see
        InteractionStrength.XL_CoulombTensorClaude's docstring for why it is orbital "c", not "b", whose
        expansion coefficients are required there). Isolated from computeDirectExchangeV; only used by the
        ALFieldClaude code line, cf. solveAverageLevelFieldClaude.
        A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function computeDirectExchangeVClaude(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                      primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital},
                                      bVectors::Dict{Subshell, Vector{Float64}},
                                      tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})
    nsL     = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    matrixV = zeros( nsL+nsS, nsL+nsS )

    for  cf  in  coeffs2p
        if      subshell == cf.a  &&  subshell == cf.b  &&  subshell == cf.c  &&  subshell == cf.d
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + 2 * cf.V *
                          InteractionStrength.XL_CoulombClaude(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.c  &&  cf.b == cf.d
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_CoulombClaude(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
             else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.d  &&  cf.b == cf.c
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                (cacheLL, cacheLS, cacheSS) = tensorCaches[cf.nu]
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_CoulombTensorClaude(cf.nu, subshell, orbitals[cf.b], orbitals[cf.c],
                                                                      bVectors[cf.c], subshell,
                                                                      cacheLL, cacheLS, cacheSS, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.d  &&  cf.a == cf.c
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_CoulombClaude(cf.nu, subshell, orbitals[cf.a], subshell, orbitals[cf.c], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.c  &&  cf.a == cf.d
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu )
                (cacheLL, cacheLS, cacheSS) = tensorCaches[cf.nu]
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_CoulombTensorClaude(cf.nu, subshell, orbitals[cf.a], orbitals[cf.d],
                                                                      bVectors[cf.d], subshell,
                                                                      cacheLL, cacheLS, cacheSS, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell != cf.a  &&  subshell != cf.b  &&  subshell != cf.c  &&  subshell != cf.d  continue
        else    error("No proper interpretation of coeff = $cf")
        end
    end

    return( matrixV )
end


"""
`SelfConsistent.computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                  orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)` 
    ... computes the value of the energy functional upon which this MCDHF optimization is based on.
        The procedure assumes that this functional is simply given by the form 
        
             EF = Sum_t  coeff1p(a_t,b_t) * I(a_t,b_t)  +   
                  Sum_t  coeff2p(k_t; a_t,b_t,c_t,d_t) * <a_t,b_t | R^(k_t) | c_t,d_t>
                  
        where  <a_t,b_t | R^(k_t) | c_t,d_t>  denotes the rank_k two-electron interaction strengths of the orbitals
        a, ,b , c, and d. An energy::Float64 is returned. 
"""
function computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                           orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)
    energy = 0.

    # Collect one-electron contributions
    for  cf  in  coeffs1p
        jj     = Basics.subshell_2j(cf.a)
        energy = energy + cf.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
    end
    
    # Collect two-electron contributions
    for  cf  in  coeffs2p
        energy = energy + cf.V * InteractionStrength.XL_Coulomb(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                                       orbitals[cf.c], orbitals[cf.d], grid)
    end

    return( energy )
end


"""
`SelfConsistent.computeFunctionalClaude(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                        orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)`
    ... computes the same MCDHF energy functional as computeFunctional, but using InteractionStrength.
        XL_CoulombClaude (kink-aware two-electron Slater integral) instead of InteractionStrength.XL_Coulomb for
        the two-electron term. Isolated from computeFunctional; only used by the ALFieldClaude code line, cf.
        solveAverageLevelFieldClaude. An energy::Float64 is returned.
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
    elseif   settings.scField in [Basics.ALFieldClaude()]
        basis     = SelfConsistent.solveAverageLevelFieldClaude(basis, nm, primitives, settings; printout=printout)
    elseif   settings.scField in [Basics.ALFieldClaude2()]
        basis     = SelfConsistent.solveAverageLevelFieldClaude2(basis, nm, primitives, settings; printout=printout)
    elseif   typeof(settings.scField) == Basics.EOLField
        # solveOptimizedLevelField already returns a complete, correctly-diagonalized multiplet (built via
        # its own internal, kink-aware Hamiltonian.performCIClaude call on the converged orbitals) -- return
        # it directly. Re-diagonalizing below would be redundant AND wrong: EOLField isn't in the
        # ALFieldClaude/ALFieldClaude2 check just below, so it would fall through to the bare,
        # non-kink-aware Hamiltonian.performCI, silently discarding the kink-aware result.
        return( SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   settings.scField in [Basics.ALFieldClaude(), Basics.ALFieldClaude2()]
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
    elseif   settings.scField in [Basics.ALFieldClaude()]
        basis     = SelfConsistent.solveAverageLevelFieldClaude(basis, nm, primitives, settings; printout=printout)
    elseif   settings.scField in [Basics.ALFieldClaude2()]
        basis     = SelfConsistent.solveAverageLevelFieldClaude2(basis, nm, primitives, settings; printout=printout)
    elseif   typeof(settings.scField) == Basics.EOLField
        # See the identical note in the other performSCF overload just above: solveOptimizedLevelField
        # already returns a complete, correctly (kink-aware) diagonalized multiplet -- return it directly
        # rather than redundantly re-diagonalizing with the wrong, non-kink-aware Hamiltonian.performCI.
        return( SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   settings.scField in [Basics.ALFieldClaude(), Basics.ALFieldClaude2()]
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
    chemMu    = determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, primitives.grid);        @show chemMu
    
    # Extract the kappa's from orbitals
    kappas = Int64[];     for (k,v)  in  orbitals     push!(kappas, k.kappa)    end;    kappas = unique(kappas);   @show kappas

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
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), wa, wb)
            
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
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                       settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        for the average-level (AL) functional. In addition, the settings::AsfSettings are taken into account.
        A (new) basis::Basis is returned.
"""
function solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                settings::AsfSettings; printout::Bool=true)
    nsL    = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    rotate = false
    
    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once at the beginning   
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot   = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    meanOcc  = Basics.extractMeanOccupation(basis)
    orbitals = deepcopy( basis.orbitals )
    bVectors = Dict{Subshell, Vector{Float64}}()
    v = zeros(nsL+nsS);     for  (k, orb)  in  orbitals   bVectors[k] = v   end
    
    # (2) Generate angular coefficients for all diagonal matrix elements including the division by nCSF 
    # to obtain mean energy functional
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(settings.scField, basis::Basis)
    for  cf in  coeffs1p   wx = Basics.twice(Basics.subshell_j(cf.a)) + 1;   @show cf, cf.T * sqrt(wx)   end
    
    # (3) Iterate the orbital set of the given basis until convergence is achieved; initialize bVectors in order to keep 
    # the B-spline coefficients of each orbital for the orbital projections
 
    for iter = 1:settings.maxIterationsScf
        println("\n> SCF interation $(iter): ")
        newOrbitals = Dict{Subshell, Orbital}();   newbVectors = Dict{Subshell, Vector{Float64}}()
        
        # (4) Rotate orbitals pairwise if required; such rotations are not needed if both subshells are filled
        #     determine a maximum rotation parameter
        if rotate
            orbitals = SelfConsistent.rotateOrbitals(basis.subshells, orbitals, grid, settings)
        end 
        
        # (5) Cycle through all subshells from the inner-to-outer subshells; 
        for subshell  in  basis.subshells
            occ      = meanOcc[subshell];  occm1 = 1.0 / occ   
            orb      = orbitals[subshell] 
            print(">> Refine $subshell orbital with mean occ = $occ ... ")
            
            # (6) Set-up the Hamiltonian matrix for this shell, including the Dirac Hamiltonian + direct + exchange potentials
            matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
            matrixV = SelfConsistent.computeDirectExchangeV(subshell, coeffs2p, primitives, orbitals)
            matrix  = matrix + occm1 * matrixV

            # (7) Apply projections of matrix upon the other orbitals of the same symmetry (accounts for Langrange multipliers)
            matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, bVectors)
            
            # (8) Diagonalize the Hamiltonian matrix and refine the orbital and bVectors
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrix, matrixB)
            newOrb                = Bsplines.generateOrbitalFromPrimitives(subshell, wc, primitives)
            newOrbitals[subshell] = newOrb
            newbVectors[subshell] = Bsplines.extractVectorFromPrimitives(subshell,   wc, primitives)
            
            # (9) Report about the new orbital
            ovlap = abs( RadialIntegrals.overlap(orb, newOrb, grid) )
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end
        
        # (10) Compute energy functional and analyze the convergence of orbitals
        eFunctional = SelfConsistent.computeFunctional(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        orbitalConv = SelfConsistent.computeConvergence(orbitals, newOrbitals, grid)
        bVectorConv = SelfConsistent.computeConvergence(bVectors, newbVectors, matrixB)
        
        println(">> Total energy = $(eFunctional*1)   orbital-conv = $orbitalConv   orbital-acc = $(1.0 - orbitalConv)   " *
                "B-conv = $bVectorConv")
        
        if  abs(1.0 - orbitalConv) < settings.accuracyScf    break   end
        orbitals = deepcopy(newOrbitals);    bVectors = newbVectors
    end
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end


"""
`SelfConsistent.solveAverageLevelFieldClaude(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                             settings::AsfSettings; printout::Bool=true)`
    ... solves the same average-level (AL) self-consistent field as solveAverageLevelField, but using the
        kink-aware two-electron Slater-integral treatment (computeDirectExchangeVClaude, computeFunctionalClaude)
        throughout, both for the SCF orbital-optimization matrix and for the reported energy functional. Isolated
        from solveAverageLevelField; only reached via performSCF's scField = Basics.ALFieldClaude() dispatch.
        A (new) basis::Basis is returned.
"""
function solveAverageLevelFieldClaude(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                      settings::AsfSettings; printout::Bool=true)
    nsL    = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    rotate = false

    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once at the beginning
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot   = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    meanOcc  = Basics.extractMeanOccupation(basis)
    orbitals = deepcopy( basis.orbitals )

    # Claude fix (19-Jul-2026): bVectors must hold each subshell's ACTUAL initial B-spline expansion
    # coefficients, matching orbitals -- NOT all-zero vectors. The original solveAverageLevelField zero-inits
    # bVectors because it only ever uses them for Hamiltonian.projectHamiltonian (where an all-zero "other
    # orbital" on iteration 1 is a harmless no-op, since there is nothing yet to project against). But this
    # Claude line ALSO uses bVectors as the "cVector" input to InteractionStrength.XL_CoulombTensorClaude
    # (the exchange-integral tensor contraction), where an all-zero vector on iteration 1 silently makes the
    # ENTIRE exchange contribution zero for that iteration, not a harmless no-op.
    #
    # Second Claude fix (19-Jul-2026): it is not enough to just be nonzero -- bVectors must also be obtained by
    # Bsplines.fitVectorToPrimitivesClaude (projecting the CLEANED, already-truncated tabulated orbital back onto
    # the B-spline basis), NOT by Bsplines.extractVectorFromPrimitives (the RAW diagonalization eigenvector).
    # generateOrbitalFromPrimitives truncates the tabulated P,Q arrays at their own mtp and zeros small values,
    # but the raw eigenvector carries no such cleanup in its tail coefficients. Reconstructing one smooth
    # tabulated function from those (as the original row-wise XL_CoulombClaude does, via orbitals[.].P/.Q)
    # tolerates that noise just fine; summing the SAME coefficients directly, weighted against many separate
    # cached Phi_(i,m) entries (as XL_CoulombTensorClaude does), does not enjoy the same cancellation, and was
    # traced to a real SCF divergence starting from iteration 2 (initially tiny, ~1e-9-level orbital differences
    # were amplified into an ~0.4 relative matrixV disagreement within a couple of iterations). Fitting bVectors
    # from the already-cleaned orbitals ties their tail behaviour to the SAME, already-trusted truncation
    # criterion used everywhere else, instead of inheriting whatever noise the diagonalization leaves behind.
    bVectors = Dict{Subshell, Vector{Float64}}()
    for  sh  in  basis.subshells
        bVectors[sh] = Bsplines.fitVectorToPrimitivesClaude(orbitals[sh], primitives, matrixB)
    end

    # (2) Generate angular coefficients for all diagonal matrix elements including the division by nCSF
    # to obtain mean energy functional; the angular-coefficient logic itself is unaffected by the kink fix,
    # so the ALField() method is called directly rather than dispatching on settings.scField (ALFieldClaude)
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(Basics.ALField(), basis::Basis)
    for  cf in  coeffs1p   wx = Basics.twice(Basics.subshell_j(cf.a)) + 1;   @show cf, cf.T * sqrt(wx)   end

    # (2a) Precompute the kink-aware Slater-moment tensor caches (RadialIntegrals.SlaterMomentCacheClaude), once
    # per multipole rank appearing anywhere in coeffs2p -- this is the one-time, basis-only cost that lets every
    # exchange-type Fock-matrix contraction below (InteractionStrength.XL_CoulombTensorClaude) run without any
    # further per-call adaptive quadrature, for every subshell and every SCF iteration. Direct-type contractions
    # (InteractionStrength.XL_CoulombClaude) do not need these caches -- a single screened potential per call is
    # already cheap for that case -- so building a cache per rank here even though only the exchange branches use
    # it is a modest, deliberate over-build in exchange for simplicity.
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    # (3) Iterate the orbital set of the given basis until convergence is achieved; initialize bVectors in order to keep
    # the B-spline coefficients of each orbital for the orbital projections

    for iter = 1:settings.maxIterationsScf
        println("\n> SCF interation $(iter) [Claude]: ")
        newOrbitals = Dict{Subshell, Orbital}();   newbVectors = Dict{Subshell, Vector{Float64}}()

        # (4) Rotate orbitals pairwise if required; such rotations are not needed if both subshells are filled
        #     determine a maximum rotation parameter
        if rotate
            orbitals = SelfConsistent.rotateOrbitals(basis.subshells, orbitals, grid, settings)
        end

        # (5) Cycle through all subshells from the inner-to-outer subshells;
        for subshell  in  basis.subshells
            occ      = meanOcc[subshell];  occm1 = 1.0 / occ
            orb      = orbitals[subshell]
            print(">> Refine $subshell orbital with mean occ = $occ ... ")

            # (6) Set-up the Hamiltonian matrix for this shell, including the Dirac Hamiltonian + direct + exchange potentials
            matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
            matrixV = SelfConsistent.computeDirectExchangeVClaude(subshell, coeffs2p, primitives, orbitals, bVectors, tensorCaches)
            matrix  = matrix + occm1 * matrixV

            # (7) Apply projections of matrix upon the other orbitals of the same symmetry (accounts for Langrange multipliers)
            # DISABLED again, paused for reassessment (19-Jul-2026). Progress made: bVectors not being exactly
            # S-normalized (b'*matrixB*b == 1) -- required for the projection operator (I - S*bb') to be truly
            # idempotent -- was a real, confirmed bug (Bsplines.fitVectorToPrimitivesClaude now re-normalizes its
            # output, kept regardless of this flag's state). That fix reduced Be's period-2 oscillation amplitude
            # by >100x (from an ~1 Hartree swing to ~0.004 Hartree), but did not eliminate it, and the residual
            # oscillation still converges to the wrong energy. Leading hypothesis: the paper's projection formula
            # assumes same-kappa orbitals are varied SIMULTANEOUSLY, which may not be numerically compatible, as
            # implemented here, with JAC's one-subshell-at-a-time, lagged (Gauss-Seidel) iteration order. Left
            # disabled so the tensor-exchange work (validated correct on He/Be/Ne) is not left in a regressed
            # state; the last main() message before pausing this thread proposed a post-diagonalization
            # Gram-Schmidt orthogonalization as a more robust alternative, not yet tried.
            ## matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, bVectors)

            # (8) Diagonalize the Hamiltonian matrix and refine the orbital and bVectors
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrix, matrixB)
            newOrb                = Bsplines.generateOrbitalFromPrimitives(subshell, wc, primitives)
            newOrbitals[subshell] = newOrb
            # Claude fix: fit from the CLEANED newOrb (Bsplines.fitVectorToPrimitivesClaude), not the raw
            # eigenvector (Bsplines.extractVectorFromPrimitives) -- see the note above where bVectors is first
            # initialized for why this matters.
            newbVectors[subshell] = Bsplines.fitVectorToPrimitivesClaude(newOrb, primitives, matrixB)

            # (9) Report about the new orbital
            ovlap = abs( RadialIntegrals.overlap(orb, newOrb, grid) )
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end

        # (9a) Orthonormalize the just-refined orbital SET within each kappa group, treating same-kappa
        # subshells as simultaneously valid rather than privileging whichever was processed first/last in the
        # loop above -- the "orthonormalize the orbital set" step of Zatsarinny & Froese Fischer, Comput. Phys.
        # Commun. 202, 287-303 (2016), Fig. 1, kept as its own step (not folded into the per-orbital
        # Hamiltonian.projectHamiltonian call above, which remains disabled pending further study).
        (newOrbitals, newbVectors) = SelfConsistent.orthonormalizeSameKappaClaude(newOrbitals, newbVectors, basis.subshells,
                                                                                   primitives, matrixB)

        # (10) Compute energy functional and analyze the convergence of orbitals
        eFunctional = SelfConsistent.computeFunctionalClaude(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        orbitalConv = SelfConsistent.computeConvergence(orbitals, newOrbitals, grid)
        bVectorConv = SelfConsistent.computeConvergence(bVectors, newbVectors, matrixB)

        println(">> Total energy = $(eFunctional*1)   orbital-conv = $orbitalConv   orbital-acc = $(1.0 - orbitalConv)   " *
                "B-conv = $bVectorConv")

        # Claude fix: update orbitals/bVectors to the just-computed newOrbitals/newbVectors BEFORE testing for
        # convergence, not after -- the original solveAverageLevelField skips this update when the break fires,
        # so it returns orbitals that are one iteration stale relative to the eFunctional/orbitalConv just
        # printed. Invisible whenever consecutive iterations are already nearly identical (as for He), but not
        # for a system that is still changing significantly between iterations at the point of convergence.
        orbitals = deepcopy(newOrbitals);    bVectors = newbVectors
        if  abs(1.0 - orbitalConv) < settings.accuracyScf    break   end
    end

    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end


"""
`SelfConsistent.computeTwoElectronVClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                           bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                           tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})`
    ... computes the direct and exchange two-electron potential matrix for the given subshell, exactly as
        computeDirectExchangeVClaude does, but taking bVectors::Dict{Subshell,Vector{Float64}} -- B-spline
        expansion coefficients -- as the SOLE, canonical per-subshell state, rather than a persistent
        Dict{Subshell,Orbital}. Every coeffs2p entry relevant to subshell involves exactly one "partner"
        subshell (the diagonal self-term is its own partner); a partner's tabulated form is built only as a
        disposable, read-only byproduct (Bsplines.generateOrbitalFromVectorClaude), cached within this call
        so repeated coefficients sharing the same partner do not re-evaluate it, and discarded once this
        function returns -- never re-fit back into bVectors, since bVectors are never derived FROM a
        tabulated form here in the first place; they come directly from diagonalization
        (solveAverageLevelFieldClaude2). Isolated from computeDirectExchangeVClaude/computeDirectExchangeV;
        only used by the ALFieldClaude2 code line. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is
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
        session against the DBSR_HF av_energy_coef reference for all of Ne's coefficients). Isolated from
        computeDirectExchangeVClaude; only used by the ALFieldClaude2 code line. A (nsL+nsS) x (nsL+nsS)
        matrix::Array{Float64,2} is returned.
"""
function computeFockMatrixClaude2(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                  bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                  nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                  occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}})
    matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
    matrixV = computeTwoElectronVClaude2(subshell, coeffs2p, bVectors, primitives, tensorCaches)
    return( matrix + (1.0/occ) * matrixV )
end


"""
`SelfConsistent.solveAverageLevelFieldClaude2(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                              settings::AsfSettings; printout::Bool=true)`
    ... solves the same average-level (AL) self-consistent field as solveAverageLevelField/
        solveAverageLevelFieldClaude, but with a bVector-native architecture modeled directly on DBSR_HF
        (Zatsarinny & Froese Fischer, CPC 202, 287 (2016)): B-spline expansion coefficients
        (Dict{Subshell,Vector{Float64}}) are the SOLE, canonical, persistent per-iteration orbital state --
        no Dict{Subshell,Orbital} is maintained across iterations at all, and there is no
        Bsplines.fitVectorToPrimitivesClaude "fit back" step, since bVectors are never derived from a
        tabulated form to begin with; they come directly from Basics.diagonalize's eigenvector for each
        subshell, in turn. Orthogonality between same-kappa subshells (e.g. Ne's 1s/2s) is enforced by
        projecting the Fock matrix (Hamiltonian.projectHamiltonian, reused unchanged) against each
        ALREADY-PROCESSED lower same-kappa subshell's bVector directly inside the generalized eigenvalue
        problem, before diagonalizing -- matching DBSR_HF's hf_solve_HF.f90/hf_eiv sequential approach --
        rather than solveAverageLevelFieldClaude's post-hoc Löwdin symmetric orthogonalization of the whole
        same-kappa group after the per-orbital loop. The target eigenvalue index is shifted down by one for
        each such projection applied, mirroring hf_eiv's `mm = m + (orthogonalized-count) - 1`. A tabulated
        Orbital is only ever built as a disposable, read-only byproduct: once per unique "partner" subshell
        inside computeTwoElectronVClaude2 (for the two-electron potential contraction), and once per
        subshell per iteration purely for reporting/energy-functional evaluation (reusing
        computeFunctionalClaude unchanged) -- never stored as competing state, never refit. A final,
        single export pass (Bsplines.generateOrbitalFromVectorClaude) produces a standard
        Dict{Subshell,Orbital} for the returned basis::Basis, so every downstream consumer (properties,
        processes, CI/DCB Hamiltonian construction) is unaffected by this being a bVector-native SCF.
        Isolated from solveAverageLevelField/solveAverageLevelFieldClaude; only reached via performSCF's
        scField = Basics.ALFieldClaude2() dispatch. A (new) basis::Basis is returned.
"""
function solveAverageLevelFieldClaude2(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                       settings::AsfSettings; printout::Bool=true)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid

    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once
    if  printout    println(">> [Claude2] (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
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

    # (2) Generate angular coefficients (unaffected by the kink fix / bVector-native rebuild; ALField()
    # is called directly, exactly as solveAverageLevelFieldClaude already does)
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(Basics.ALField(), basis)

    # (3) Precompute kink-aware Slater-moment tensor caches for every rank that occurs, exactly as
    # solveAverageLevelFieldClaude does; only the exchange branches of computeTwoElectronVClaude2 use them
    neededRanks = unique( [ cf.nu for cf in coeffs2p ] )
    if  printout    println(">> [Claude2] Precompute kink-aware Slater-moment tensor caches for ranks $(neededRanks) ...")    end
    tensorCaches = Dict{Int64, NTuple{3,RadialIntegrals.SlaterMomentCacheClaude}}()
    for  L  in  neededRanks
        cacheLL = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesL, grid; rtol=1.0e-6)
        cacheLS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesL, primitives.bsplinesS, grid; rtol=1.0e-6)
        cacheSS = RadialIntegrals.buildSlaterMomentCacheClaude(L, primitives.bsplinesS, primitives.bsplinesS, grid; rtol=1.0e-6)
        tensorCaches[L] = (cacheLL, cacheLS, cacheSS)
    end

    orbitals = Dict{Subshell, Orbital}()    # only ever a disposable, per-iteration reporting byproduct

    for  iter = 1:settings.maxIterationsScf
        println("\n> SCF interation $(iter) [Claude2]: ")
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

            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrix, matrixB)
            l  = Basics.subshell_l(subshell)
            ni = nsS + subshell.n - l - count
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
    @show  sa, size(matrix), absD, absND
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
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), wa, wb)
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
    # identical to solveAverageLevelFieldClaude2
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
    (coeffs1p, coeffs2p) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, targetLevels)
    genOcc = SelfConsistent.computeGeneralizedOccupationEOL(blockCaches, targetLevels, basis)

    weights = let  twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
        sumW = sum( twiceJp1(level.J)  for level in targetLevels )
        [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    end
    weightedEnergy = sum( weights[i] * targetLevels[i].energy  for i = 1:length(targetLevels) )

    # (4) Precompute kink-aware Slater-moment tensor caches for every rank that occurs, exactly as
    # solveAverageLevelFieldClaude2 does; only the exchange branches of computeTwoElectronVClaude2 use them
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
            print(">> Refine $subshell orbital with generalized occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrixClaude2(subshell, coeffs2p, bVectors, primitives, nucPot,
                                                              storage, occ, tensorCaches)

            count = Base.count( sh2 -> sh2.kappa == subshell.kappa, keys(processedBVectors) )
            if  count > 0
                matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, processedBVectors)
            end

            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrix, matrixB)
            l  = Basics.subshell_l(subshell)
            ni = nsS + subshell.n - l - count
            rawVector = wc.vectors[ni]

            oldVector = bVectors[subshell]
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
        # for (see solveAverageLevelFieldClaude2/project_df_al_kink_bug.md) -- here in the mixing-vector <->
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

        (coeffs1p, coeffs2p) = SelfConsistent.combineAngularCoefficientsEOL(blockCaches, dampedLevels)
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
