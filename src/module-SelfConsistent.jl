
"""
`module  JAC.SelfConsistent`  
	... a submodel of JAC that contains all structs and methods to generate self-consistent fields of different 
	    kind and complexity.
"""
module SelfConsistent

using  Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..Hamiltonian, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial,
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
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cfx.b    T = T + cfx.T;    hasConsidered[icx] = true
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


## EXPERIMENTAL SWITCH (09-Aug-2026), default false = the behaviour that has always been in place.
## When true, solveOptimizedLevelField builds the Fock matrix as
##      h  +  (1/occ) * V(diagonal CSF pairs)  +  V(off-diagonal CSF pairs)
## instead of  h + (1/occ) * V(all pairs).  This tests whether the off-diagonal (CI-coupling) part should
## be scaled by the generalized occupation at all: it carries weight ~ c_r c_s while occ carries ~ c_r^2,
## so dividing it by occ introduces a c_1/c_2 factor that grows without bound as a correlating CSF's own
## coefficient shrinks -- the suspected mechanism behind the winner-take-all collapse documented in
## solveOptimizedLevelField's KNOWN LIMITATION.  Set with SelfConsistent.setEolUnscaledOffDiagonal(true).
GBL_EOL_UNSCALED_OFFDIAGONAL = false

## Re-orthonormalize the same-kappa orbitals after the damping step. DEFAULT TRUE since 10-Aug-2026;
## the switch is kept only so the two behaviours can still be compared.
##
## Hamiltonian.projectHamiltonian makes each raw eigenvector orthogonal to the already-processed same-kappa
## orbitals, but the damping that follows, mixed = 0.5*old + 0.5*raw, mixes it back with the PREVIOUS
## iteration's vector, which is NOT orthogonal to them -- and nothing restored it. The CSF expansion assumes
## an orthonormal orbital set, so the resulting energies were not legitimate variational numbers.
## Measured on Li 1s^2 2s + 1s^2 3s + 1s^2 3d (three s-orbitals in kappa = -1), converged:
##
##                     <2s|3s>      <1s|2s>      E
##     as it was      -1.128e-03   -5.558e-05   -7.4335291982
##     re-orthonorm.  -1.373e-11   -5.358e-13   -7.4335284248
##
## The energy RISES by 7.7e-07 Ha, which is the honest direction: the non-orthogonal set was giving a
## slightly-too-low number. The whole approved test suite is blind to this (44/44 either way), which is why
## it went unnoticed -- TestFrames.testMethod_OrbitalOrthonormality now asserts it directly.
## The fix needs no new code: orthonormalizeSameKappa was written for exactly this and had never
## been called from anywhere.
GBL_SCF_REORTHONORMALIZE = true

## Anderson depth for the AVERAGE-LEVEL field, separate from the mean-field one above because the iterate is
## different: there it is the screening potential, here the orbitals themselves.  0 = the plain damped
## iteration exactly as before.
##
## ON since 17-Aug-2026, at the same depth 2 the mean-field driver uses.  It reaches the SAME solution --
## Ar 3s^2 3p^6 agrees to 4.3e-9 once accuracyScf is tight enough to converge at all -- and is 1.4x to 1.9x
## faster, the gain GROWING with the accuracy demanded (1.47x at 1e-6, 1.79x at 1e-9, 1.89x at 1e-12).
## On Be 4-config it is also the more STABLE of the two: across accuracyScf = 1e-6, 1e-9, 1e-12 it drifts by
## 1e-5 where the plain iteration swings 6.5e-5 NON-MONOTONICALLY, and it gets there in 211 s against 1093 s.
##
## It could not be switched on until 18aaf5f.  Anderson perturbs the orbital tails just enough to move the
## old mtp cut, which made TestFrames.testMethod_OrbitalOrthonormality report 6.5e-08 where the plain
## iteration gave 6.8e-10 -- an artefact of the truncated integral, not of the orbitals, which were
## orthonormal to 1e-17 throughout.  With the tails kept, both give ~1e-16.
##
## STILL OPEN, and NOT fixed by any of this: plain AL does not converge for a multi-configuration basis with
## near-degenerate CSFs (Be 4-config), and the default accuracyScf = 1e-6 hides it by stopping early.
## Tightening the tolerance is therefore not a general cure -- right for Ar-like cases, worse for Be.
GBL_AL_ANDERSON_DEPTH = 2

## Anderson depth for the mean-field (DFS/HS) SCF.  0 = the plain iteration exactly as before; a positive
## value routes performSCF to SelfConsistent.solveMeanFieldBasisAnderson, which reaches the SAME self-consistent
## solution in fewer iterations.  DEFAULT 0 so that nothing changes unless it is asked for.
##
## The plain iteration converges linearly, the residual shrinking by a constant factor r per step; measured
## 12-Aug-2026, r = 0.44 (Ar 1s^2..3p^6), 0.57 (Ne 1s^2 2s^2 2p^6), 0.69 (Fe [Ar] 3d^6 4s^2).  Anderson
## mixing builds the next screening potential from a least-squares combination of the last few iterates and
## their residuals, cancelling the slowest-decaying error rather than waiting for it to decay:
##
##                        plain        depth 2      agreement of the converged orbital energies
##     Ne  2s^2 2p^6      28 it        13 it        7.9e-07
##     Ar  3s^2 3p^6      23 it        14 it        6.9e-09
##     Fe  3d^6 4s^2      45 it        18 it        6.5e-07
##     Ne+ 1s hole        17 it        11 it        1.2e-07
##
## Depth 2 is the measured optimum; 3 is nearly equal, and LARGER IS WORSE (Ne: 24 it at depth 5, 36 at 12),
## the usual ill-conditioning of a long Anderson history.  Note that depth 0 in the Anderson driver itself is a
## JACOBI sweep and does NOT converge in 60 iterations -- the Gauss-Seidel ordering of the original driver is
## what makes the plain iteration viable at all.
##
## STANDARD SINCE 12-Aug-2026 (was 0 = the plain iteration when this was first added).  Setting it to 0
## restores the old path exactly, which is how the two were compared.  Making it the default is a deliberate
## editorial act: both iterations reach the same self-consistent solution, but only to within accuracyScf,
## so results can move by ~1e-6.  What that costs was measured by regenerating every approved reference with
## it on: 27 of 29 came out BITWISE IDENTICAL, one moved by 3.0e-09, and the only large apparent change --
## test-Cascade-StepwiseDecay -- was not numerical at all, but two DEGENERATE levels swapping index labels
## (all 606 transition rows the same set).
GBL_SCF_ANDERSON_DEPTH = 2


"""
`SelfConsistent.setScfAndersonDepth(depth::Int64)`
    ... sets the Anderson-mixing depth of the mean-field SCF; 0 restores the plain iteration. Nothing is returned.
"""
function setScfAndersonDepth(depth::Int64)
    global GBL_SCF_ANDERSON_DEPTH = depth
    return( nothing )
end


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
`SelfConsistent.setScfReorthonormalize(flag::Bool)`  
    ... sets the experimental switch that re-orthonormalizes same-kappa orbitals after the SCF damping
        step; nothing is returned.
"""
function setScfReorthonormalize(flag::Bool)
    global GBL_SCF_REORTHONORMALIZE = flag
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
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cfx.b    T = T + cfx.T;    hasConsidered[icx] = true
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
`SelfConsistent.orthonormalizeSameKappa(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
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
function orthonormalizeSameKappa(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
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
            orbitalsOut[sh] = Bsplines.generateOrbitalFromVector(sh, en, vec, primitives)
        end
    end

    return( orbitalsOut, bVectorsOut )
end


"""
`SelfConsistent.computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                        orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)`
    ... computes the MCDHF energy functional using InteractionStrength.XL_CoulombKinkAware (kink-aware
        two-electron Slater integral) for the two-electron term. Used by the ALField/EOLField code line, cf.
        solveAverageLevelField. An energy::Float64 is returned.
"""
function computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                 orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)
    energy = 0.

    # Collect one-electron contributions -- unchanged, no kink in this integral
    for  cf  in  coeffs1p
        jj     = Basics.subshell_2j(cf.a)
        energy = energy + cf.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
    end

    # Collect two-electron contributions via the kink-aware integral
    for  cf  in  coeffs2p
        energy = energy + cf.V * InteractionStrength.XL_CoulombKinkAware(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                                              orbitals[cf.c], orbitals[cf.d], grid)
    end

    return( energy )
end


"""
`SelfConsistent.energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
        coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
        subshells::Array{Subshell,1}, primitives::Bsplines.Primitives, grid::Radial.Grid,
        nucPot::Radial.Potential)`  
    ... the EOL energy as a plain scalar function of the orbital B-spline coefficient vectors, with the
        angular coefficients held fixed. This is exactly the functional solveOptimizedLevelField reports,
        just expressed in the variables that are actually varied, so that it can be differentiated. Note
        that it includes coeffs1p -- which the Fock matrix of the present scheme never receives, and which
        is the inconsistency the rotation-based path exists to remove. A value::Float64 is returned.
"""
function energyFromBVectors(bVectors::Dict{Subshell, Vector{Float64}},
                                   coeffs1p::Array{SpinAngular.Coefficient1p,1},
                                   coeffs2p::Array{SpinAngular.Coefficient2p,1},
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
        coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1},
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
                                       coeffs1p::Array{SpinAngular.Coefficient1p,1},
                                       coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}})
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    orbitals = Dict{Subshell, Orbital}()
    for  sh  in  subshells
        orbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
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
    ##
    ## MATRIX-FREE (10-Aug-2026).  M is never formed.  Its entries are INT B_i B_k w_r V_L, so
    ## (M v)_i = INT B_i(r) f(r) w_r V_L(r) with f the expansion of v -- and M is symmetric within each
    ## block, so the two products a coefficient needs share one potential and differ only in the vector.
    ## Each subshell's expansion is built ONCE per gradient rather than per coefficient.
    ## The expansions must carry the SAME sign convention as the orbitals the energy is built from.
    ## Bsplines.generateOrbitalFromVector canonicalizes each orbital so that P > 0 at small r, and does NOT
    ## feed that back into bVectors -- so a raw expansion and its own orbital can be oppositely signed. The
    ## screened potential Vk below comes from the ORBITALS (b,d) while the contracted vector came from the
    ## RAW b-vector (a,c), and an off-diagonal CSF pair carries an ODD power of a correlating orbital's sign,
    ## so the mismatch survives instead of cancelling. Measured on the Be RAS step-2 case, where 2p_1/2 and
    ## 2p_3/2 are canonicalized against their raw vectors in every iteration: the analytic gradient agreed
    ## with a finite difference to five digits while only s-orbitals were involved, then went 2.5x, 19x and
    ## finally SIGN-WRONG as the 2p weight grew -- which is what stalled the line search. Same defect, and
    ## the same remedy, as the cVector sign-matching in computeTwoElectronV.
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

    for  cf  in  coeffs2p
        for  (sA, sB, sC, sD)  in  [ (cf.a, cf.b, cf.c, cf.d), (cf.b, cf.a, cf.d, cf.c) ]
            ## the same triangular-delta and parity guard XL_CoulombKinkAware applies before doing any work
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
            Vk = RadialIntegrals.buildScreenedPotential(cf.nu, orbitals[sB], orbitals[sD], primitives.grid;
                                                              mtpOut=primitives.grid.NoPoints)
            (pC, qC) = expanded[sC];     (pA, qA) = expanded[sA]
            grad[sA] = grad[sA] + (cf.V * xc * scale[sA]) * SelfConsistent.screenedProduct(Vk, pC, qC, primitives)
            grad[sC] = grad[sC] + (cf.V * xc * scale[sC]) * SelfConsistent.screenedProduct(Vk, pA, qA, primitives)
        end
    end
    return( grad )
end


"""
`SelfConsistent.gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
        virtuals::Dict{Subshell, Array{Vector{Float64},1}}, coeffs1p::Array{SpinAngular.Coefficient1p,1},
        coeffs2p::Array{SpinAngular.Coefficient2p,1}, subshells::Array{Subshell,1},
        primitives::Bsplines.Primitives, grid::Radial.Grid, nucPot::Radial.Potential; hStep::Float64=1.0e-4)`  
    ... the gradient of the EOL energy with respect to the allowed orbital rotations, by central finite
        differences. Slow but free of any derivation, so it is the reference against which an analytic
        gradient must later be checked, and on its own it answers the question "is this converged solution
        actually stationary?". Each direction is S-orthogonal to the orbital itself, so normalization
        contributes only at second order and is omitted. A Dict{Subshell, Vector{Float64}} is returned.
"""
function gradientByFiniteDifference(bVectors::Dict{Subshell, Vector{Float64}},
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
    ## (1) project each orbital on the positive branch of its kappa
    for  sh  in  subshells
        b = bVectors[sh];    v = zeros( length(b) )
        for  phi  in  posSet[sh.kappa]    v = v + (transpose(phi) * matrixB * b) * phi    end
        nrm2Full = abs( transpose(b) * matrixB * b );    nrm2Pos = abs( transpose(v) * matrixB * v )
        worst    = max( worst, 1.0 - nrm2Pos/max(nrm2Full,1.0e-30) )
        out[sh]  = v / sqrt( max(nrm2Pos, 1.0e-30) )
    end
    ## (2) S-orthonormalize within each kappa, in coefficient space, so the positive span is preserved --
    ## but ONLY where it is actually needed.  Gram-Schmidt is sequential and asymmetric: it leaves the first
    ## orbital of a kappa untouched and pushes the whole correction onto the later ones, so applying it to
    ## an already-orthonormal set rotates the orbitals for nothing.  Measured on Li, doing it unconditionally
    ## cost 8.5e-07 Ha while removing a negative-branch weight of only 3.4e-14 -- more than the entire
    ## discrepancy that sent us looking.  Skip it when the block is orthonormal to tolerance.
    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        shk = [ sh for sh in subshells if sh.kappa == kappa ]
        dev = 0.
        for  (i, sha) in enumerate(shk),  (j, shb) in enumerate(shk)
            ov  = transpose(out[sha]) * matrixB * out[shb]
            dev = max( dev, abs( ov - (i == j ? 1.0 : 0.0) ) )
        end
        if  dev < 1.0e-9    continue    end
        done = Vector{Float64}[]
        for  sh  in  shk
            v = out[sh]
            for  u  in  done    v = v - (transpose(u) * matrixB * v) * u    end
            nrm = sqrt( abs(transpose(v) * matrixB * v) )
            if  nrm > 1.0e-10    v = v / nrm;    push!(done, v);    out[sh] = v    end
        end
    end
    return( out, worst )
end


"""
`SelfConsistent.solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model,
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
function solveOptimizedLevelFieldByRotation(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                         settings::AsfSettings; printout::Bool=true, nVirtual::Int64=16,
                                         method::Symbol=:conjugate, cgRestart::Int64=10, lbfgsMemory::Int64=5)
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

    ## Built ONCE: it depends only on the nuclear potential and the B-spline basis, never on the orbitals.
    posSpectrum = SelfConsistent.positiveBranchSpectrum(basis.subshells, primitives, nucPot, matrixB, storage)
    (bVectors, _) = SelfConsistent.projectOntoPositiveBranch(bVectors, basis.subshells, primitives,
                                                                     nucPot, matrixB, storage; spectrum=posSpectrum)
    ePrevious = 0.;   tStep = 1.0;   multiplet = Multiplet("EOL-ByRotation", Level[])
    ## Direction state, all held in b-space: the virtual directions are rebuilt every iteration, so anything
    ## stored in THAT basis would be meaningless one step later.
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

        e0   = SelfConsistent.energyFromBVectors(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                         primitives, grid, nucPot)
        grad = SelfConsistent.computeOrbitalGradient(bVectors, coeffs1p, coeffs2p, basis.subshells,
                                                             primitives, nucPot, storage)
        virt = SelfConsistent.virtualDirections(bVectors, basis.subshells, primitives, nucPot,
                                                        matrixB, storage; nVirtual=nVirtual, spectrum=posSpectrum)
        ## Project the gradient on the allowed rotations, and PRECONDITION each component by the
        ## orbital-energy denominator (eps_v - eps_a) -- the diagonal of the rotation Hessian, i.e. the
        ## standard first-order estimate  kappa_av = -g_av / (eps_v - eps_a).  Plain steepest descent
        ## converges far too slowly here: it left the Li control 1.8e-6 Ha short after 60 steps.
        ## The denominator is floored, since a near-degenerate pair would otherwise produce a huge step.
        gProj = Dict{Subshell, Vector{Float64}}();   step = Dict{Subshell, Vector{Float64}}()
        denom = Dict{Subshell, Vector{Float64}}()
        gNorm = 0.;    sNorm = 0.
        for  sh  in  basis.subshells
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
        gNorm = sqrt(gNorm);    sNorm = sqrt(sNorm)

        ## Assemble the search direction in b-space.  Plain preconditioned steepest descent zigzags here:
        ## the energy falls steadily while |grad| merely oscillates, each step undoing part of the previous
        ## one.  Polak-Ribiere conjugacy reuses the previous direction to cancel that, at the cost of two
        ## dot products and one stored vector per subshell.
        sVec = Dict{Subshell, Vector{Float64}}();   gVec = Dict{Subshell, Vector{Float64}}()
        for  sh  in  basis.subshells
            sv = zeros(nsL+nsS);    gv = zeros(nsL+nsS)
            for  (iv, phi)  in  enumerate(virt[sh])
                sv = sv + step[sh][iv]  * phi
                gv = gv + gProj[sh][iv] * phi
            end
            sVec[sh] = sv;    gVec[sh] = gv
        end
        ## method = :conjugate is the DEFAULT and the one to use.  :lbfgs is kept because it is measurably
        ## better where the basis is stable -- Be 1s^2 2s^2 + 1s^2 2p^2 at 12 iterations reaches -14.618710
        ## against conjugacy's -14.616507, i.e. it does change the RATE and not merely the constant -- but it
        ## LOSES on the harder Be RAS step-2 case, -14.617374 against -14.619313, and there it discards its
        ## own curvature history at iterations 8, 15, 21 and 23.  WHY IT FAILS IS NOT KNOWN.  Two candidates
        ## were proposed and BOTH MEASURED SMALL, so neither should be repeated as an explanation:
        ##   * the virtual space is NOT churning -- successive frames overlap to |1-<new,old>| = 1e-4..3e-3
        ##     with no change in the number of directions, one 0.67 rotation excepted, and that one does not
        ##     coincide with any of the four history discards;
        ##   * the objective DOES drift, since coeffs1p/coeffs2p are rebuilt from the re-diagonalized mixing
        ##     vector every iteration, but only by max|dV| ~ 1e-4..3e-3 against max|V| ~ 1.8.
        ## Conjugacy tolerates whatever this is because it carries ONE previous direction and restarts every
        ## ten; L-BFGS accumulates five pairs and does not.  Anyone picking this up should measure first.
        dotAll(u, v) = sum( sum(u[sh] .* v[sh])  for sh in basis.subshells )
        ## H_0 is the DIAGONAL PRECONDITIONER, not the usual gamma*I: the (eps_v - eps_a) denominators carry
        ## real physics and throwing them away for a scalar would be a step backwards.  applyPrecond is what
        ## turns gVec into -sVec, so it is exactly the operator already in use.
        applyPrecond = function(q)
            r = Dict{Subshell, Vector{Float64}}()
            for  sh  in  basis.subshells
                rv = zeros(nsL+nsS)
                for  (iv, phi)  in  enumerate(virt[sh])
                    rv = rv + ( (transpose(phi) * q[sh]) / denom[sh][iv] ) * phi
                end
                r[sh] = rv
            end
            return( r )
        end
        ## Record the curvature pair (s,y) of the step just taken.  Only pairs with <y,s> > 0 are kept, which
        ## is what keeps the implicit inverse Hessian positive definite -- and hence the direction a descent
        ## direction -- on a surface that is not convex.
        if  method == :lbfgs  &&  !isempty(bPrev)
            sPair = Dict{Subshell, Vector{Float64}}();   yPair = Dict{Subshell, Vector{Float64}}()
            for  sh  in  basis.subshells
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
            ## two-loop recursion, giving d = -H grad with H built from the stored pairs around H_0
            q = Dict{Subshell, Vector{Float64}}( sh => copy(gVec[sh])  for sh in basis.subshells )
            alphas = zeros( length(sHist) )
            for  i = length(sHist):-1:1
                alphas[i] = rhoHist[i] * dotAll(sHist[i], q)
                for  sh  in  basis.subshells    q[sh] = q[sh] - alphas[i] * yHist[i][sh]    end
            end
            r = applyPrecond(q)
            for  i = 1:length(sHist)
                bb = rhoHist[i] * dotAll(yHist[i], r)
                for  sh  in  basis.subshells    r[sh] = r[sh] + (alphas[i] - bb) * sHist[i][sh]    end
            end
            for  sh  in  basis.subshells    dir[sh] = -r[sh]    end
        elseif  method == :conjugate
            if  !isempty(dirPrev)  &&  iterSinceRestart < cgRestart  &&  abs(sgPrev) > 1.0e-30
                num = 0.
                for  sh  in  basis.subshells    num = num + sum( (gVec[sh] - gPrev[sh]) .* sVec[sh] )    end
                beta = max( 0., num / sgPrev )                   ## Polak-Ribiere+, i.e. restart on beta < 0
            end
            for  sh  in  basis.subshells
                dir[sh] = beta == 0. ? sVec[sh] : sVec[sh] + beta * dirPrev[sh]
            end
        else
            for  sh  in  basis.subshells    dir[sh] = sVec[sh]    end
        end
        ## Guard: a conjugate direction must still descend.  If it does not, fall back to steepest descent.
        dg = 0.;   for sh in basis.subshells   dg = dg + sum( gVec[sh] .* dir[sh] )   end
        if  dg >= 0.
            for  sh  in  basis.subshells    dir[sh] = sVec[sh]    end
            beta = 0.
        end
        iterSinceRestart = beta == 0. ? 0 : iterSinceRestart + 1
        sgPrev = 0.;   for sh in basis.subshells   sgPrev = sgPrev + sum( gVec[sh] .* sVec[sh] )   end
        gPrev  = gVec;    dirPrev = dir
        bPrev  = Dict{Subshell, Vector{Float64}}( sh => copy(bVectors[sh])  for sh in basis.subshells )
        if  sNorm < 1.0e-14    break    end
        if  printout
            println(">> [EOL-C3] iter $iter:  E = $(multiplet.levels[1].energy)   |grad| = $gNorm   step = $tStep")
        end
        ## Converged when the GRADIENT is small. settings.accuracyScf is its tolerance, which is the honest
        ## test: the old sole criterion, |E - E_prev| < accuracyScf, halts when PROGRESS is slow, which is a
        ## statement about the optimizer and not about the solution -- and it is why every EOL value quoted
        ## before 16-Aug-2026 was an upper bound. It cannot stand ALONE, though: |grad| PLATEAUS at a floor
        ## set by the basis and the projection, so on Li a pure gradient test ran 48 further iterations after
        ## the energy had stopped moving, for nothing. Both tests are kept, and the driver says which fired.
        if  gNorm < settings.accuracyScf
            if  printout   println(">> [EOL-C3] converged: |grad| = $gNorm < $(settings.accuracyScf).")   end
            break
        end

        ## backtracking line search along -gProj; halve until the energy actually falls
        accepted = false
        for  trial = 1:24
            newB = Dict{Subshell, Vector{Float64}}()
            for  sh  in  basis.subshells
                v = bVectors[sh] + tStep * dir[sh]
                newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
            end
            ## The projection must happen BEFORE the acceptance test, not after it.  Projecting and
            ## re-orthonormalizing moves the orbitals, so accepting `newB` on the strength of its own energy
            ## and then storing the PROJECTED vector stores something that was never tested -- and the
            ## projection gives a little of the gain back each step.  That broke monotonicity: on Li the
            ## driver reduced the gradient 19-fold while the energy ROSE by 3.7e-7 Ha, which a descent
            ## method cannot do.  Testing the projected vector restores  E(new) < E(old)  by construction,
            ## and with it the guarantee that the CI eigenvalue falls too (it is bounded above by this
            ## fixed-coefficient functional, and equals it at the previous orbitals).
            (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
            eTrial = SelfConsistent.energyFromBVectors(projB, coeffs1p, coeffs2p, basis.subshells,
                                                               primitives, grid, nucPot)
            if  eTrial < e0
                bVectors = projB
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
            ## A line search that finds no descent along an L-BFGS direction means the stored curvature is
            ## no longer describing this surface -- unsurprising, since virtualDirections is rebuilt every
            ## iteration and the pairs then mix vectors from different subspaces. Discard the history and
            ## retry from the preconditioned gradient rather than giving up.
            if  printout    println(">> [EOL-C3] L-BFGS history discarded at iteration $iter; retrying.")   end
            empty!(sHist);   empty!(yHist);   empty!(rhoHist);   tStep = 1.0
            for  trial = 1:24
                newB = Dict{Subshell, Vector{Float64}}()
                for  sh  in  basis.subshells
                    v = bVectors[sh] + tStep * sVec[sh]
                    newB[sh] = v / sqrt( abs(transpose(v) * matrixB * v) )
                end
                (projB, negW) = SelfConsistent.projectOntoPositiveBranch(newB, basis.subshells,
                                                    primitives, nucPot, matrixB, storage; spectrum=posSpectrum)
                eTrial = SelfConsistent.energyFromBVectors(projB, coeffs1p, coeffs2p, basis.subshells,
                                                                   primitives, grid, nucPot)
                if  eTrial < e0    bVectors = projB;   accepted = true;   break    end
                tStep = tStep / 2
            end
        end
        if  !accepted
            if  printout    println(">> [EOL-C3] no descent found; stopping at iteration $iter.")    end
            break
        end
        if  iter > 1  &&  abs(multiplet.levels[1].energy - ePrevious) < 1.0e-11
            if  printout
                println(">> [EOL-C3] stopped on a stagnant energy (|dE| < 1.0e-11) with |grad| = $gNorm; " *
                        "this is a converging UPPER BOUND, not a converged gradient.")
            end
            break
        end
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
    # orbitals are seeded -- StartFromPrevious inherits a grid just as much as StartFromHydrogenic does. The
    # occupations are handed over so that each subshell is tested at the charge it actually feels; at the bare
    # charge the check rejects the valence orbital of any heavy near-neutral system, whose box must be matched
    # to a screened orbital some thirty times more extended than the bare-Z one.
    occupations = Dict{Shell,Int64}()
    for  conf in configs
        for  (sh, occ)  in conf.shells
            if  occ > 0     occupations[sh] = max( get(occupations, sh, 0), occ )    end
        end
    end
    Bsplines.checkGridRepresentation(subshells, nuclearModel.Z, primitives; occupations = occupations,
                                     accuracy = settings.gridAccuracy, stopper = settings.gridStopper)

    # Initialize the orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("> Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = Bsplines.generateOrbitalsHydrogenic(subshells, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromThomasFermi
        if  printout   println("> Start SCF process with orbitals in a Thomas-Fermi potential.")   end
        ## The nucleus screened by a statistical model of the electron cloud.  Unlike every self-consistent
        ## field this needs no density, so it is available before any orbital exists -- which is the point of
        ## a start potential.  Bsplines.generateOrbitals then does what it does for any other potential.
        tfPot     = Basics.add( Nuclear.nuclearPotential(nuclearModel, primitives.grid),
                                Basics.computePotential(Basics.ThomasFermiField(), primitives.grid,
                                                        nuclearModel.Z, NoElectrons) )
        orbitals  = Bsplines.generateOrbitals(subshells, tfPot, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("> Start SCF process from given list of orbitals.energy")    end
        ## Take what the given set provides and fall back to a hydrogenic orbital for anything it does not.
        ## THE FALLBACK NEVER RAN BEFORE (fixed 13-Aug-2026): it called HydrogenicIon.radialOrbital(subsh, ...,
        ## grid) where neither `subsh` nor `grid` exists in this method -- the loop variable is `sh` and only
        ## `primitives` is passed -- and radialOrbital takes a Shell rather than a Subshell in any case, so the
        ## branch could only ever have thrown.  It went unnoticed because StartFromPrevious had exactly one
        ## caller, which always supplied a complete set.  Warm-starting one cascade block from another does
        ## not: consecutive blocks differ in which subshells are occupied.  The missing ones are now generated
        ## by the same B-spline routine the StartFromHydrogenic branch above uses, in ONE call.
        orbitals = Dict{Subshell, Orbital}()
        missingSubshells = Subshell[]
        for  sh in subshells
            if  haskey(settings.startScfFrom.orbitals, sh)   orbitals[sh] = settings.startScfFrom.orbitals[sh]
            else                                             push!(missingSubshells, sh)
            end
        end
        if  length(missingSubshells) > 0
            if  printout   println("> Start orbitals do not cover $missingSubshells; these are taken hydrogenic.")   end
            hydrogenic = Bsplines.generateOrbitalsHydrogenic(missingSubshells, nuclearModel, primitives; printout=false)
            for  sh in missingSubshells   orbitals[sh] = hydrogenic[sh]   end
        end
    else  error("stop b")
    end
    
    basis = Basis(true, NoElectrons, subshells, csfs, coreSubshells, orbitals)
    return( basis )
end


"""
`SelfConsistent.checkScFieldIsSupported(scField::Basics.AbstractScField)`
    ... verifies that the given field is one that performSCF can actually iterate, and raises an explanatory
        error if it is not. Several members of Basics.AbstractScField are POTENTIALS rather than fields: they
        answer Basics.providesPotential but not Basics.providesScfDriver, and one of them is not
        self-consistent even in principle. Checked HERE, before a grid, a set of primitives and a many-electron
        basis have been built, so that an unsupported choice costs nothing and says what to do instead.
        Nothing is returned.
"""
function checkScFieldIsSupported(scField::Basics.AbstractScField)
    Basics.providesScfDriver(scField)   &&   return( nothing )
    sa = "\n\nSelfConsistent.performSCF(): $(nameof(typeof(scField))) is not a self-consistent field that " *
         "this driver can iterate.\n"
    if      typeof(scField) == Basics.ThomasFermiField
        sa = sa * ">>> It is not self-consistent even in principle: it needs the nuclear charge and the "     *
                  "electron number only, and no\n    density at all. That is exactly what makes it useful "   *
                  "as a STARTING potential -- use ManyElectron.StartFromThomasFermi,\n    or ask for the "    *
                  "potential itself with Basics.computePotential.\n"
    elseif  typeof(scField) in [Basics.AaDFSField, Basics.AaHSField]
        sa = sa * ">>> It is an AVERAGE-ATOM potential for finite temperature: its Basics.computePotential "  *
                  "takes a chemical potential\n    and a temperature, so it belongs to the plasma line "      *
                  "rather than to this bound-state SCF. Use Plasma.AverageAtomScheme\n    (see "              *
                  "examples/example-Ja.jl).\n"
    else
        sa = sa * ">>> It is a screened POTENTIAL, not a field: it owns a Basics.computePotential method but " *
                  "no SCF driver of its own.\n    Use Basics.computePotential(...) to obtain it, or choose "  *
                  "one of the fields listed below.\n"
    end
    sa = sa * ">>> The fields performSCF iterates are:  " * join(Basics.scfDriverFields(), ", ") *
              "  (NuclearField with a hydrogenic start only).\n"
    error(sa)
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model,
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)`
    ... performs a SCF computation for which NO grid is given, so that the radial box is derived from the
        configurations themselves by Basics.recommendedGrid; a multiplet::Multiplet is returned.

        This exists because choosing the box is the step most likely to be got wrong, and getting it wrong
        does not look like a grid problem: the record attributes four separate "bugs" -- an E3 rate 1000x too
        small, a Zeeman kappa <= -3 failure, a MultipolePolarizibility defect and a Breit sign flip -- to a
        box that did not match the orbitals, and each was first blamed on the angular machinery.  The derived
        box beats JAC's hand-chosen default grid for every system it has been measured on (see
        Basics.recommendedGrid), by 2.9e-3 Ha for argon and 2.3e-2 Ha for Ti+.

        ONE CASE NEEDS THE GRID GIVEN BY HAND.  The estimate cannot tell a spectroscopic Rydberg shell from a
        CORRELATION shell of the same n and l, and reads both as diffuse: a beryllium basis carrying 3s and 3d
        for correlation is given 67 a.u., which is right for a real 1s^2 2s 3s state and far too generous for
        a correlation orbital that contracts onto the valence region.  Where the high-n shells are there to
        correlate rather than to be occupied, pass a grid, or pass `rbox` to Basics.recommendedGrid.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model,
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    grid = Basics.recommendedGrid(configs, nm, printout=printout)

    return( SelfConsistent.performSCF(configs, nm, grid, settings,
                                      levelSymmetries=levelSymmetries, printout=printout) )
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
    
    SelfConsistent.checkScFieldIsSupported(settings.scField)

    # Generate primitives and initialize the many-electron basis
    Defaults.setDefaults("standard grid", grid)
    primitives = Bsplines.generatePrimitives(grid)    
    basis      = SelfConsistent.initializeBasis(configs, nm, primitives, settings; levelSymmetries, printout)
    
    # Solve a self-consistent field for this basis
    scfProc = Basics.scfProcedure(settings.scField)
    if   scfProc == :meanFieldIteration
        ## GBL_SCF_ANDERSON_DEPTH = 0 keeps the plain iteration; a positive depth reaches the SAME
        ## self-consistent solution in fewer iterations (see the note at the switch).
        if  GBL_SCF_ANDERSON_DEPTH > 0
            basis = SelfConsistent.solveMeanFieldBasisAnderson(basis, nm, primitives, settings; printout=printout,
                                                            andersonDepth=GBL_SCF_ANDERSON_DEPTH)
        else
            basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout)
        end 
    elseif   scfProc == :hydrogenicStartOnly  &&  settings.startScfFrom == StartFromHydrogenic()
        # Return the basis as already generated.
    elseif   scfProc == :averageLevel
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   scfProc == :optimizedLevel
        # EOL is done by ORBITAL ROTATION.  The older solveOptimizedLevelField, which folds the off-diagonal
        # CSF-pair terms into the same (1/occ)-scaled Fock matrix, converges to a DEGENERATE stationary point
        # whenever two near-degenerate CSFs compete for one correlation channel: the correlating weight runs
        # to zero, the correlation orbital then no longer enters the energy, and its gradient vanishes for a
        # trivial reason.  Measured on Be 1s^2 2s^2 + 1s^2 2p^2 it lands 19.4 mHa ABOVE the average-level
        # field on the very level it is asked to optimize, with the 2p weight collapsed from 0.25 to 0.0001.
        # Rotating the orbitals instead escapes that point and reaches 5.3 mHa BELOW AL.  See example-Ao.jl.
        #
        # The rotation is a LOCAL optimizer, so it starts from an average-level basis rather than from the
        # initial guess -- that is how it was validated, and a hydrogenic start has no reason to lie in its
        # basin.  Both solvers return a complete, correctly (kink-aware) diagonalized multiplet, so return it
        # directly; falling through would re-diagonalize with the bare, non-kink-aware Hamiltonian.performCI.
        alSettings = AsfSettings(settings; scField = Basics.ALField())
        basis      = SelfConsistent.solveAverageLevelField(basis, nm, primitives, alSettings; printout=printout)
        return( SelfConsistent.solveOptimizedLevelFieldByRotation(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Now that the orbitals are final, check that no symmetry has converged onto the wrong state. This
    # catches what Bsplines.checkGridRepresentation cannot: that check tests HYDROGENIC orbitals -- since
    # 17-Aug-2026 at the SCREENED charge rather than the bare one, which is what lets a heavy near-neutral
    # system through at all -- and a hydrogenic test can never see a symmetry that has converged onto a
    # different state, Ge II 4f on a 614 a.u. box being the case that motivated it. Note the EOLField branch
    # above returns early and is therefore not covered here.
    Bsplines.checkOrbitalConsistency(basis.orbitals, grid; stopper = settings.gridStopper)

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   scfProc == :averageLevel
        mp = Hamiltonian.performCIKinkAware(basis, nm, grid, settings, printout=printout)
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
    
    SelfConsistent.checkScFieldIsSupported(settings.scField)

    # Generate primitives
    primitives = Bsplines.generatePrimitives(grid)    
    
    # Solve a self-consistent field for this basis
    scfProc = Basics.scfProcedure(settings.scField)
    if   scfProc == :meanFieldIteration
        ## GBL_SCF_ANDERSON_DEPTH = 0 keeps the plain iteration; a positive depth reaches the SAME
        ## self-consistent solution in fewer iterations (see the note at the switch).
        if  GBL_SCF_ANDERSON_DEPTH > 0
            basis = SelfConsistent.solveMeanFieldBasisAnderson(basis, nm, primitives, settings; printout=printout,
                                                            andersonDepth=GBL_SCF_ANDERSON_DEPTH)
        else
            basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout)
        end 
    elseif   scfProc == :hydrogenicStartOnly  &&  settings.startScfFrom == StartFromHydrogenic()
        # Return the basis as already generated.
    elseif   scfProc == :averageLevel
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   scfProc == :optimizedLevel
        # See the note in the other performSCF overload just above: EOL is done by orbital rotation, started
        # from an average-level basis, and returns a complete, correctly (kink-aware) diagonalized multiplet.
        alSettings = AsfSettings(settings; scField = Basics.ALField())
        basis      = SelfConsistent.solveAverageLevelField(basis, nm, primitives, alSettings; printout=printout)
        return( SelfConsistent.solveOptimizedLevelFieldByRotation(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Now that the orbitals are final, check that no symmetry has converged onto the wrong state. This
    # catches what Bsplines.checkGridRepresentation cannot: that check tests HYDROGENIC orbitals -- since
    # 17-Aug-2026 at the SCREENED charge rather than the bare one, which is what lets a heavy near-neutral
    # system through at all -- and a hydrogenic test can never see a symmetry that has converged onto a
    # different state, Ge II 4f on a 614 a.u. box being the case that motivated it. Note the EOLField branch
    # above returns early and is therefore not covered here.
    Bsplines.checkOrbitalConsistency(basis.orbitals, grid; stopper = settings.gridStopper)

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   scfProc == :averageLevel
        mp = Hamiltonian.performCIKinkAware(basis, nm, grid, settings, printout=printout)
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
              (where it originated as `determineChemicalPotential`), since Plasma.perform(::AverageAtomScheme,
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
                ## Collected as well: in a long run this line scrolls away, and nobody learns afterwards that a
                ## field never converged.  The accuracy is rounded so that repeated identical failures collapse
                ## into one counted entry; see Defaults.warn.
                Defaults.warn(AddWarning(), "SelfConsistent.solveAverageAtomField(): the SCF did NOT converge -- " *
                              "stopped at accuracy " * @sprintf("%.1e", accuracyScf) * " after 32 iterations.")
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
`SelfConsistent.computeTwoElectronV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                           bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                           tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}})`
    ... computes the direct and exchange two-electron potential matrix for the given subshell, taking
        bVectors::Dict{Subshell,Vector{Float64}} -- B-spline expansion coefficients -- as the SOLE, canonical
        per-subshell state, rather than a persistent Dict{Subshell,Orbital}. Every coeffs2p entry relevant to
        subshell involves exactly one "partner" subshell (the diagonal self-term is its own partner); a
        partner's tabulated form is built only as a disposable, read-only byproduct
        (Bsplines.generateOrbitalFromVector), cached within this call so repeated coefficients sharing
        the same partner do not re-evaluate it, and discarded once this function returns -- never re-fit back
        into bVectors, since bVectors are never derived FROM a tabulated form here in the first place; they
        come directly from diagonalization (solveAverageLevelField/solveOptimizedLevelField). Shared
        by the ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is
        returned.
"""
function computeTwoElectronV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                    bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                    tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}})
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
            # so it is handled by the identical XL_CoulombTensor call, unweighted, as :exchange.
            partner = cf.c;       pattern = :exchange
        elseif  subshell == cf.c  &&  subshell == cf.d  &&  cf.a == cf.b  &&  cf.a != subshell
            partner = cf.a;       pattern = :exchange
        else
            continue
        end
        # SpinAngular does not pre-filter its coefficients on the C^L parity rule, so combinations whose
        # reduced matrix element <a||C^nu||c> vanishes arrive here carrying a NONZERO V -- seven of them even
        # for Be 1s^2 2s^2 + 1s^2 2p^2. The radial kernels annihilate them, so skipping is numerically inert;
        # it only avoids building a matrix already known to be zero, and stops the kernels' own guards from
        # firing thousands of times per SCF. The rule depends on the pattern, which is exactly what the
        # 26-Jul note below got wrong: the direct call passes a = c and b = d, so it collapses to `nu even`,
        # whereas the exchange call pairs (subshell, partner) and needs l_sh + l_partner + nu even.
        if      pattern in [:diagonal, :direct]   &&   isodd(cf.nu)                                     continue
        elseif  pattern == :exchange              &&
                isodd( Basics.subshell_l(subshell) + Basics.subshell_l(partner) + cf.nu )               continue
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
            Bsplines.generateOrbitalFromVector(partner, 0.0, bVectors[partner], primitives)
        end

        if      pattern == :diagonal
            matrixV = matrixV + 2 * cf.V *
                      InteractionStrength.XL_CoulombKinkAware(cf.nu, subshell, partnerOrb, subshell, partnerOrb, primitives)
        elseif  pattern == :direct
            matrixV = matrixV + cf.V *
                      InteractionStrength.XL_CoulombKinkAware(cf.nu, subshell, partnerOrb, subshell, partnerOrb, primitives)
        else    # :exchange
            # partnerOrb was built by Bsplines.generateOrbitalFromVector, which silently canonicalizes its
            # sign so that P is positive at small r (see that function's `wSign` step) -- a convention applied
            # to the TABULATED reconstruction only, never fed back into bVectors[partner] itself. Whenever that
            # canonicalization actually flips the sign, partnerOrb and the raw bVectors[partner] end up
            # oppositely signed, even though XL_CoulombTensor's "b"/"c" arguments (from partnerOrb) and its
            # "cVector" argument (from bVectors[partner]) are meant to represent the identical orbital -- an odd
            # total power of partner's sign then leaks into the exchange matrix element (confirmed empirically
            # this session: exchange terms flip sign in exact lockstep with this mismatch, while diagonal/direct
            # terms -- which only ever use partnerOrb, an even number of times -- stay sign-invariant as they
            # should). Fix: mirror generateOrbitalFromVector's own small-r wSign test on the raw bVector,
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
                      InteractionStrength.XL_CoulombTensor(cf.nu, subshell, partnerOrb, partnerOrb,
                                                                  cVector, subshell,
                                                                  cacheLL, cacheLS, cacheSS, primitives)
        end
    end

    return( matrixV )
end


"""
`SelfConsistent.computeFockMatrix(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                         bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                         nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                         occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}})`
    ... computes the full per-subshell Hamiltonian matrix (one-electron + direct + exchange), matching the
        role of DBSR_HF's hf_matrix.f90: Bsplines.setupLocalMatrix provides the one-electron (kinetic +
        nuclear + rest-mass) block natively in B-spline basis (reused unchanged), and
        computeTwoElectronV provides the two-electron block, divided by the subshell's own mean
        occupation -- algebraically identical to DBSR_HF's own hfm=qsum(i)*dhl (one-electron scaled up)
        followed by hfm=hfm/qsum(i) (whole matrix scaled down) once expanded, given the two-electron
        coefficients already carry their full (not per-electron) occupation weighting (verified this
        session against the DBSR_HF av_energy_coef reference for all of Ne's coefficients). Shared by the
        ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrix::Array{Float64,2} is returned.
"""
function computeFockMatrix(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1},
                                  bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                  nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                  occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}};
                                  coeffs2pUnscaled::Array{SpinAngular.Coefficient2p,1}=SpinAngular.Coefficient2p[])
    # occ == 0 would give 1/occ = Inf and, against the zero two-electron matrix such a subshell has,
    # Inf * 0 = NaN in every element -- a whole Fock matrix of NaN that only surfaces much later, and
    # then as a completely unrelated-looking complaint about the negative-energy continuum.  A subshell
    # with no occupation has no mean field to be refined in; the caller must decide what to do with it
    # (solveOptimizedLevelField carries it forward unchanged), so refuse here rather than return NaN.
    if  abs(occ) < 1.0e-12
        error("SelfConsistent.computeFockMatrix(): subshell $subshell has occupation $occ. There is no " *
              "mean field to define for an unoccupied subshell; the caller must skip it instead.")
    end
    matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
    matrixV = computeTwoElectronV(subshell, coeffs2p, bVectors, primitives, tensorCaches)
    ## coeffs2pUnscaled, when given, contributes WITHOUT the 1/occ scaling. Empty by default, so the
    ## returned matrix is bit-for-bit what it always was unless a caller opts in.
    if  length(coeffs2pUnscaled) > 0
        matrixU = computeTwoElectronV(subshell, coeffs2pUnscaled, bVectors, primitives, tensorCaches)
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
    ## Anderson history over the CONCATENATED b-vectors. One full Gauss-Seidel sweep is the fixed-point map
    ## g(x); the 0.5 damping already inside it stays, so depth 0 reproduces the previous behaviour exactly.
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

        for  subshell  in  basis.subshells
            occ = meanOcc[subshell]
            print(">> Refine $subshell orbital with mean occ = $occ ... ")

            matrix = SelfConsistent.computeFockMatrix(subshell, coeffs2p, bVectors, primitives, nucPot,
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
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, newEnergies[sh], newBVectors[sh], primitives)
        end
        ## ANDERSON ACCELERATION on the orbitals themselves.  The mean-field driver already accelerates its
        ## screening potential this way, with 1.6-2.6x fewer iterations and the same solution to ~1e-7; here
        ## the iterate is the concatenated b-vector set and g(x) is one full sweep.  The sweep sign-aligns
        ## each eigenvector against the previous one, so g is sign-consistent with x and the residual means
        ## what it should.  depth <= 0 leaves the plain damped iteration untouched.
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
        ## Under test: restore the same-kappa orthonormality that the damping step destroys.
        if  SelfConsistent.GBL_SCF_REORTHONORMALIZE
            (newOrbitals, newBVectors) = SelfConsistent.orthonormalizeSameKappa(newOrbitals, newBVectors,
                                                                basis.subshells, primitives, matrixB)
        end
        ## Convergence is measured against what is ACTUALLY accepted, which Anderson may have moved. Only
        ## when it is active, so that depth 0 remains a bit-for-bit control on the previous behaviour.
        if  andersonDepth > 0
            for  sh  in  basis.subshells
                dpm[sh] = 1.0 - abs( transpose(bVectors[sh]) * matrixB * newBVectors[sh] )
            end
        end
        eFunctional = SelfConsistent.computeFunctional(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        ## The overlap defect 1 - |<b_old|b_new>| is QUADRATIC in the orbital change; the change itself is
        ## || b_new - b_old ||_B = sqrt(2 * defect), and that is what a user means by "the orbitals still move
        ## by".  Both are reported so that the tolerance cannot be misread by a factor of its own square root.
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

    ## Say which of the two happened.  An iteration that merely runs out of steps used to end in silence, so
    ## that a stopped field and a converged one printed the same thing and were quoted the same way.
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
                ## Collected as well; see the note at solveAverageAtomField above.
                Defaults.warn(AddWarning(), "SelfConsistent.solveMeanFieldBasis(): the SCF did NOT converge for " *
                              string(basis.subshells) * " -- stopped at accuracy " * @sprintf("%.1e", accuracyScf) *
                              " after $(settings.maxIterationsScf) iterations.")
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
`SelfConsistent.solveMeanFieldBasisAnderson(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                          settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=2)`
    ... solves the same mean-field (DFS/HS) self-consistent problem as SelfConsistent.solveMeanFieldBasis, but
        accelerates the iteration with Anderson mixing. The physics is identical -- the same functional, the same
        fixed point; only the route to it differs, and the returned basis::Basis is the SAME self-consistent
        solution, reached in fewer iterations.

        WHY. The plain iteration is an undamped fixed-point (Picard) map: build the screening potential from the
        current orbitals, solve the one-electron Dirac equation in it, take the resulting orbitals as the new
        ones. Such a map converges LINEARLY, the residual shrinking by a roughly constant factor r each step, so
        the cost to reach a given accuracy scales as log(tol)/log(r). Measured on this code (12-Aug-2026), r is
        0.44 for Ar 1s^2..3p^6, 0.57 for Ne 1s^2 2s^2 2p^6 and 0.69 for Fe [Ar] 3d^6 4s^2 -- 17, 25 and 37
        iterations respectively to reach 1e-6, and the open-3d case additionally sits at residual ~1 for seven
        iterations before descending at all. Anderson mixing forms the next screening potential from a
        least-squares combination of the last few iterates and their residuals, which cancels the slowest-decaying
        error components instead of waiting for them to decay geometrically.

        WHAT IS MIXED. The fixed-point variable is the ELECTRONIC screening potential on the radial grid, not the
        orbitals: it is one vector, the nuclear part is constant and drops out, and the map
        potential -> orbitals -> potential is exactly the self-consistency condition. andersonDepth = 0 recovers
        the plain iteration (and is then Jacobi rather than Gauss-Seidel, see below); the default 2 is the measured optimum here (3 is nearly equal, larger is WORSE: 24 iterations at 5 and 36 at 12 for Ne), and it
        choice and rarely needs changing.

        ONE DELIBERATE DIFFERENCE FROM solveMeanFieldBasis. There the potential is rebuilt inside the kappa loop,
        so each symmetry is solved in a field the earlier symmetries have already improved (a Gauss-Seidel
        sweep); here it is built once per iteration and all symmetries see the same field (a Jacobi sweep), which
        is what makes a single well-defined residual per iteration -- and hence Anderson -- possible. Both have
        the SAME fixed point, so a converged result is unaffected; only an iteration stopped short of
        convergence could differ.

        A basis::Basis with the self-consistent orbitals is returned.
"""
function solveMeanFieldBasisAnderson(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                   settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=2)
    storage  = Dict{String,Array{Float64,2}}()
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)
    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end

    ## The screening potential built from a given set of orbitals -- the right-hand side of the self-consistency
    ## condition, and the vector Anderson works on.
    function screeningZr(orbitals)
        wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
        NoCsf  = length(wBasis.csfs)
        wmc    = zeros( NoCsf );   wN = 0.
        for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
        for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
        wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
        return( Basics.computePotential(settings.scField, primitives.grid, wLevel).Zr )
    end

    xUsed    = screeningZr(previousOrbitals)          ## the screening potential the next solve will use
    xHistory = Vector{Float64}[];   fHistory = Vector{Float64}[]
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " *
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
                ## Collected as well; see the note at solveAverageAtomField above.
                Defaults.warn(AddWarning(), "SelfConsistent.solveMeanFieldBasisAnderson(): the SCF did NOT converge for " *
                              string(basis.subshells) * " -- stopped at accuracy " * @sprintf("%.1e", accuracyScf) *
                              " after $(settings.maxIterationsScf) iterations.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        pot = Basics.add(nuclearPotential, Radial.Potential("mean field", xUsed, primitives.grid))
        #
        for kappa in kappas
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
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
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
        #
        ## The self-consistency residual: what the orbitals just obtained say the screening potential should be,
        ## minus what was actually used to obtain them. It vanishes exactly at the self-consistent solution.
        gNew = screeningZr(previousOrbitals);    fNew = gNew - xUsed
        if  andersonDepth <= 0
            xUsed = gNew                                                        ## plain (Jacobi) iteration
        else
            push!(xHistory, copy(xUsed));   push!(fHistory, fNew)
            if  length(xHistory) > andersonDepth + 1
                popfirst!(xHistory);   popfirst!(fHistory)
            end
            m = length(xHistory) - 1
            if  m < 1     xUsed = gNew                                          ## no history yet
            else
                ## Least squares: choose the combination of the last m residual DIFFERENCES that best cancels
                ## the current residual, then apply it to the iterates as well (standard Anderson, beta = 1).
                dF = zeros( length(fNew), m );   dX = zeros( length(fNew), m )
                for j = 1:m   dF[:,j] = fHistory[j+1] - fHistory[j];   dX[:,j] = xHistory[j+1] - xHistory[j]   end
                local gamma
                try
                    gamma = dF \ fNew
                catch
                    gamma = zeros(m)      ## a rank-deficient window: fall back to the plain step
                end
                if  any(!isfinite, gamma)   gamma = zeros(m)   end
                xUsed = xUsed + fNew - (dX + dF) * gamma
            end
        end
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
            newOrbitals[sh] = Bsplines.generateOrbitalFromVector(sh, 0.0, bVectors[sh], primitives)
        end
        ## The EOL driver damps exactly as the AL one does, so it loses same-kappa orthogonality in exactly
        ## the same way and needs the same repair.  Measured before this was added: Si^2+ [Ne] 3s^2 + 3p^2
        ## gave a worst same-kappa overlap of 2.4e-06 under EOL against 9.3e-10 under AL, and Si^+ reached
        ## 6.9e-05.  Wiring the switch into solveAverageLevelField alone was an oversight.
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
    multiplet  = Hamiltonian.performCIKinkAware(finalBasis, nuclearModel, grid, settings; printout=printout)
    return( multiplet )
end

end # module
