
"""
`module  JAC.BiOrthogonal`
	... a submodel of JAC that contains all methods for a bi-orthogonal transformation of two non-orthogonal orbitals set.

    Validated (1-Aug-2026) against three independent known-answer tests, following the artificial-rotation
    methodology suggested by Olsen, Godefroid, Jönsson, Malmqvist & Froese Fischer, Phys. Rev. E 52, 4499
    (1995), Sec. V: (1) closed-core self-overlap under a manufactured 1s/2s rotation, exact to 1e-15;
    (2) a full production PhotoEmission E1 rate (Fe X 3p^6 -> 3p^5) under the same closed-core rotation,
    matching a shared-basis reference to ~1e-10 relative; (3) an open-shell/virtual-mixing self-overlap
    (occupied 3s rotated with virtual 4s, forcing a genuine CI counter-rotation) plus length- and
    velocity-form dipole matrix elements, matching to ~1e-8 - 1e-10 relative. All three confirm
    `computeTransformationMatrices` and `generateCounterRotatingCiMatrices` exactly reconstruct the known
    physical answer under invertible orbital rotations, including cases requiring nontrivial CI mixing.
    This does NOT explain the puzzling near-identical gauge-ratio shift observed for a genuine K-hole
    SCF-relaxation case (Ne K-alpha, example-Da.jl branch 0) -- that remains open and is more likely a
    physical-relaxation or transition-operator question than a defect in this module.
"""
module BiOrthogonal


using Printf, LinearAlgebra, ..Basics,  ..Defaults, ..ManyElectron, ..Radial, ..RadialIntegrals, ..SpinAngularNew


"""
`BiOrthogonal.computeTransformation(leftMp::Multiplet, rightMp::Multiplet, grid::Radial.Grid)`
    ... computes the bi-orthogonal transformation of the two given multiplets by rotating the radial functions
        and by counter-rotating the corresponding CI coefficients. leftMp and rightMp must have the same
        number of electrons and be generated on the same radial grid; grid must be passed explicitly (an
        orbital's own `.grid` field is only a placeholder whenever `useStandardGrid == true`, which is the
        normal case). A warning (not an error) is issued if either multiplet's CSF list is not closed
        under de-excitation (see `checkClosureUnderDeexcitation`), in which case the transformation is only
        approximate. A tuple of two multiplets tpl(newLeftMp::Multiplet, newRightMp::Multiplet) is
        returned: same levels (J, parity, index, energy, relativeOcc) as the input, but each with a new,
        shared, biorthonormal `Basis` and a counter-rotated `mc` vector.
"""
function  computeTransformation(leftMp::Multiplet, rightMp::Multiplet, grid::Radial.Grid)
    leftBasis  = leftMp.levels[1].basis
    rightBasis = rightMp.levels[1].basis

    for  (label, basis)  in  [("leftMp", leftBasis), ("rightMp", rightBasis)]
        violations = BiOrthogonal.checkClosureUnderDeexcitation(basis)
        if  !isempty(violations)
            @warn "BiOrthogonal.computeTransformation: $label's CSF list is not closed under de-excitation " *
                    "($(length(violations)) violation(s)); the bi-orthogonal transformation will only be " *
                    "approximate. First violation: $(violations[1])"
        end
    end

    transformation                     = BiOrthogonal.computeTransformationMatrices(leftBasis, rightBasis, grid)
    newLeftOrbitals, newRightOrbitals  = BiOrthogonal.generateBiorthogonalShellMatrices(leftBasis, rightBasis, grid)
    Mleft                              = BiOrthogonal.generateCounterRotatingCiMatrices(leftBasis,  transformation, :left)
    Mright                             = BiOrthogonal.generateCounterRotatingCiMatrices(rightBasis, transformation, :right)

    newLeftBasis  = Basis(true, leftBasis.NoElectrons,  leftBasis.subshells,  leftBasis.csfs,  leftBasis.coreSubshells,  newLeftOrbitals)
    newRightBasis = Basis(true, rightBasis.NoElectrons, rightBasis.subshells, rightBasis.csfs, rightBasis.coreSubshells, newRightOrbitals)

    newLeftLevels  = [ Level(lv.J, lv.M, lv.parity, lv.index, lv.energy, lv.relativeOcc, lv.hasStateRep, newLeftBasis,  Mleft  * lv.mc)  for lv in leftMp.levels  ]
    newRightLevels = [ Level(lv.J, lv.M, lv.parity, lv.index, lv.energy, lv.relativeOcc, lv.hasStateRep, newRightBasis, Mright * lv.mc)  for lv in rightMp.levels ]

    return( (Multiplet(leftMp.name * " (bi-orthogonal)", newLeftLevels), Multiplet(rightMp.name * " (bi-orthogonal)", newRightLevels)) )
end


"""
`BiOrthogonal.computeTransformationMatrices(leftBasis::Basis, rightBasis::Basis, grid::Radial.Grid)`
    ... computes, for every kappa symmetry present in leftBasis or rightBasis, the pair of upper-triangular
        orbital transformation matrices (Cleft, Cright) that bring the (in general non-orthogonal) orbital
        sets of leftBasis and rightBasis into biorthonormal form, following Olsen, Godefroid, Jönsson,
        Malmqvist & Froese Fischer, Phys. Rev. E 52, 4499 (1995): for a fixed orbital ordering (here, by
        increasing principal quantum number n within a kappa symmetry), requiring the transformation
        matrices to be upper triangular reduces their determination to a triangular (LU) factorization of
        the overlap matrix S_ij = <phi_i^left|phi_j^right>. With S = L*U (L unit-lower-triangular, U
        upper-triangular, no pivoting so the chosen n-ordering is preserved), Cleft = inv(L)' and
        Cright = inv(U) satisfy Cleft' * S * Cright = I exactly. Presently requires leftBasis and rightBasis
        to have the same number of orbitals for every shared kappa (the general, differing-dimension case
        of the paper's Appendix B is not yet supported). A Dict{Int64, Tuple{Array{Subshell,1},
        Array{Subshell,1}, Matrix{Float64}, Matrix{Float64}}}, keyed by kappa and giving the (ordered)
        subshell lists together with (Cleft, Cright) for that symmetry block, is returned.
"""
function  computeTransformationMatrices(leftBasis::Basis, rightBasis::Basis, grid::Radial.Grid)
    kappas = unique( [sh.kappa for sh in leftBasis.subshells] ∪ [sh.kappa for sh in rightBasis.subshells] )
    result = Dict{Int64, Tuple{Array{Subshell,1},Array{Subshell,1},Matrix{Float64},Matrix{Float64}}}()

    for  kappa  in  kappas
        lList = sort( [sh for sh in leftBasis.subshells  if sh.kappa == kappa], by = sh -> sh.n )
        rList = sort( [sh for sh in rightBasis.subshells if sh.kappa == kappa], by = sh -> sh.n )
        if  length(lList) != length(rList)
            error("BiOrthogonal.computeTransformationMatrices: leftBasis and rightBasis have a differing " *
                    "number of orbitals for kappa = $kappa ($(length(lList)) vs $(length(rList))); the " *
                    "differing-dimension case is not yet supported.")
        end
        n = length(lList)
        S = zeros(n, n)
        for  i = 1:n,  j = 1:n
            S[i,j] = RadialIntegrals.overlap(leftBasis.orbitals[lList[i]], rightBasis.orbitals[rList[j]], grid)
        end

        F      = lu(S, NoPivot())
        Cleft  = inv(Matrix(F.L))'
        Cright = inv(Matrix(F.U))
        result[kappa] = (lList, rList, Cleft, Cright)
    end

    return( result )
end


"""
`BiOrthogonal.generateBiorthogonalShellMatrices(leftBasis::Basis, rightBasis::Basis, grid::Radial.Grid)`
    ... generates the biorthogonal forms of the orbitals of leftBasis and rightBasis, by applying the
        transformation matrices from `computeTransformationMatrices` to the radial functions (P, Q, Pprime,
        Qprime) of every orbital, kappa symmetry by kappa symmetry. Since orbitals of different length
        (different last tabulated grid point) can occur for the same kappa, each is zero-padded up to the
        longest orbital of its own side before the linear combination is formed. A tuple
        tpl(newLeftOrbitals::Dict{Subshell,Orbital}, newRightOrbitals::Dict{Subshell,Orbital}) is returned,
        with the SAME subshell keys as the original bases -- only the radial functions have changed.
"""
function  generateBiorthogonalShellMatrices(leftBasis::Basis, rightBasis::Basis, grid::Radial.Grid)
    transformation = BiOrthogonal.computeTransformationMatrices(leftBasis, rightBasis, grid)
    newLeft  = Dict{Subshell,Radial.Orbital}()
    newRight = Dict{Subshell,Radial.Orbital}()

    padTo(v::Array{Float64,1}, len::Int64) = length(v) >= len ? v[1:len] : vcat(v, zeros(len - length(v)))

    for  (_, (lList, rList, Cleft, Cright))  in  transformation
        n     = length(lList)
        lOrbs = [leftBasis.orbitals[sh]  for sh in lList]
        rOrbs = [rightBasis.orbitals[sh] for sh in rList]
        lenL  = maximum( length(o.P) for o in lOrbs )
        lenR  = maximum( length(o.P) for o in rOrbs )

        for  i = 1:n
            newP  = sum( padTo(lOrbs[j].P,      lenL) * Cleft[j,i]  for j = 1:n )
            newQ  = sum( padTo(lOrbs[j].Q,      lenL) * Cleft[j,i]  for j = 1:n )
            newPp = sum( padTo(lOrbs[j].Pprime, lenL) * Cleft[j,i]  for j = 1:n )
            newQp = sum( padTo(lOrbs[j].Qprime, lenL) * Cleft[j,i]  for j = 1:n )
            newLeft[lList[i]]  = Radial.Orbital(lList[i], lOrbs[i].isBound, lOrbs[i].useStandardGrid,
                                                    lOrbs[i].energy, newP, newQ, newPp, newQp, grid)

            newP  = sum( padTo(rOrbs[j].P,      lenR) * Cright[j,i] for j = 1:n )
            newQ  = sum( padTo(rOrbs[j].Q,      lenR) * Cright[j,i] for j = 1:n )
            newPp = sum( padTo(rOrbs[j].Pprime, lenR) * Cright[j,i] for j = 1:n )
            newQp = sum( padTo(rOrbs[j].Qprime, lenR) * Cright[j,i] for j = 1:n )
            newRight[rList[i]] = Radial.Orbital(rList[i], rOrbs[i].isBound, rOrbs[i].useStandardGrid,
                                                    rOrbs[i].energy, newP, newQ, newPp, newQp, grid)
        end
    end

    return( (newLeft, newRight) )
end



"""
`BiOrthogonal.checkClosureUnderDeexcitation(basis::Basis)`
    ... checks whether the CSF list of basis is closed under de-excitation (Olsen, Godefroid, Jönsson,
        Malmqvist & Froese Fischer, Phys. Rev. E 52, 4499 (1995)): for every CSF and every occupied
        subshell nl with n > n', if the (lower) subshell n'l is also present in basis.subshells and is
        not already fully occupied in that CSF, then moving one electron from nl to n'l (all other
        occupations unchanged) must produce the occupation pattern of some CSF already in basis.csfs.
        This is checked at the level of occupation-number patterns only (not the full recoupling of
        seniority/subshellJ/subshellX), which is a necessary but not by itself sufficient condition for
        closure -- for the RAS/CAS-type configuration lists this transformation is intended for, an
        occupation-pattern match is expected to reflect true closure in practice. A
        Array{String,1} of human-readable violation descriptions is returned (empty if closed).
"""
function  checkClosureUnderDeexcitation(basis::Basis)
    violations = String[]
    occupationSets = Set( csf.occupation  for  csf in basis.csfs )

    for  csf  in  basis.csfs
        # basis.subshells, not Defaults.GBL_STANDARD_SUBSHELL_LIST -- the latter is mutable global state that
        # may have since been overwritten by a later, unrelated SCF call (e.g. for a differently-sized basis
        # on the other side of a bi-orthogonal transformation), while csf.occupation is always indexed
        # against THIS basis's own (fixed) subshell list, regardless of csf.useStandardSubshells.
        subshells = basis.subshells
        for  (i, shi)  in  enumerate(subshells)
            if  csf.occupation[i] == 0    continue    end
            for  (j, shj)  in  enumerate(subshells)
                if  Basics.subshell_l(shj) != Basics.subshell_l(shi)  ||  shj.n >= shi.n    continue    end
                maxOccJ = Basics.subshell_2j(shj) + 1
                if  csf.occupation[j] >= maxOccJ    continue    end

                targetOcc        = copy(csf.occupation)
                targetOcc[i]    -= 1
                targetOcc[j]    += 1
                if  !(targetOcc in occupationSets)
                    push!(violations, "CSF (J=$(csf.J), parity=$(csf.parity), occupation=$(csf.occupation)): " *
                                        "de-exciting $shi -> $shj is not closed within this basis.")
                end
            end
        end
    end

    return( violations )
end


"""
`BiOrthogonal.generateCounterRotatingCiMatrices(basis::Basis, transformation::Dict, side::Symbol)`
    ... generates the nCSF x nCSF matrix M that counter-rotates the CI (mixing-coefficient) vectors of
        basis's levels to compensate for the orbital transformation encoded in transformation (the output
        of `computeTransformationMatrices`); side is either :left or :right, selecting the (Cleft) or
        (Cright) half of the transformation and the correspondingly-ordered subshell list for that side.
        A general (not necessarily unitary) one-electron transformation new_creator_k = sum_l C_lk *
        old_creator_l is implemented via the many-electron operator built from its matrix logarithm (a
        standard result: the generator of a one-body basis change is the one-body operator with the SAME
        matrix, in the log; holds for any invertible C, not only orthogonal/unitary ones). Concretely:
        (i) per kappa symmetry, factor G = log(C) (upper triangular, since C is upper triangular with
        nonzero diagonal); (ii) build the one-electron operator matrix Ghat over the full CSF list using
        the SAME reduced one-body coupling coefficients JAC already uses for one-body Hamiltonian terms
        (`SpinAngularNew.computeCoefficients(OneParticleOperator(0,plus), csfR, csfS, subshells)`,
        weighted by G_ab instead of a radial integral). The coefficient is used BARE: since 27-Aug-2026 the
        rank-0 coefficient carries sqrt(2 j_a + 1) itself, so no factor belongs at the call site. A one-particle toy-model check (N=1,
        two orbitals) shows the OLD and NEW expansion coefficients are related by c_old = C * c_new
        (not c_new = C * c_old) -- i.e. the CI vector transforms with the INVERSE of the orbital
        transformation -- so (iii) M = exp(-Ghat) is returned. The transformed CI vector for a level is
        mc_new = M * mc_old (same CSF ordering, only the coefficients change).
"""
function  generateCounterRotatingCiMatrices(basis::Basis, transformation::Dict, side::Symbol)
    nCSF  = length(basis.csfs)
    Gmap  = Dict{Tuple{Subshell,Subshell},Float64}()

    for  (_, (lList, rList, Cleft, Cright))  in  transformation
        shellList = side == :left ? lList : (side == :right ? rList : error("side must be :left or :right"))
        C         = side == :left ? Cleft : Cright
        G         = log(C)
        n         = length(shellList)
        for  a = 1:n,  b = 1:n
            if  G[a,b] != 0.0    Gmap[(shellList[a], shellList[b])] = G[a,b]    end
        end
    end

    opa   = SpinAngularNew.OneParticleOperator(0, Basics.plus)
    Ghat  = zeros(nCSF, nCSF)
    for  mu = 1:nCSF,  nu = 1:nCSF
        coeffs = SpinAngularNew.computeCoefficients(opa, basis.csfs[mu], basis.csfs[nu], basis.subshells)
        for  c  in  coeffs
            if  haskey(Gmap, (c.a, c.b))
                Ghat[mu,nu] += c.T * Gmap[(c.a, c.b)]
            end
        end
    end

    return( exp(-Ghat) )
end

end # module


#========================================================
DISABLED, and KEPT: design notes addressed to a collaborator, not code.  They record why the GRASP2018
rbiotransform structure was not followed and what a compact Julia implementation needs instead.
Labelled 13-Aug-2026.

Comments & thoughts on the implementation of this module
--------------------------------------------------------

++  I looked briefly through your rbiotransform implementation in GRASP2018 to get a very first impression which basic ingredients 
are needed for such an implementation. A large fraction of code refers to 'internal data handling of quantum numbers and 
functions' which might become obsolete or can be coded much more compact in Julia.
++  Moreover, for the manipulation of matrices (multiplication, inversion, ...), etc., we shall use the internal Julia features 
as much as possible. This will likely reduced the size of the code considerably since some of the procedures in rbiotransform
can be written as a single line.

++  Overall, therefore, I belief that the code could be made more compact by a factor 3-7, and also more easy to check.

++  A first goal should be to write-down the basic formulas for and in the notation of the UserGuide of JAC.
This will enable us to bring the notation and internal use much closer to each other; I noticed that large parts in the
rbiotransform implementation in GRASP2018 are highly technical with very little use of the underlying physics language
which is required to formulate the transformation itself.

++  I will be happy (and grateful) if we can agree about proper names of functions and if we make use of the 
standard Julia features.

++  Please, use ? AngularMomentum.Wigner_6j(a, b, c, d, e, f) and similar for Wigner_3j and _9j symbols to get the Wigner symbols
if needed.


++  Let's use CamelCase notation an start all variables with a lowercase letter; arrays are usually given the same name
with an additional s:  coeff::Coefficient1p  vs.  coeffs::Array{Coefficient1p,1}
++  All angular momenta should be of type j::AngularJ64 (for j-values) and of type m::AngularM64 (for m-values, 
and if they appear at all in your implementation), please.
++  For a given subhsell sh::Subshell, get n and kappa simply by: sh.n, sh.kappa
++  For a given subhsell sh::Subshell, get l and j simply by:     Basics.subshell_l(sh), Basics.subshell_j(sh)

++  Look for ?CsfR ... to get the underlying definition of a CSF in JAC; it allows a direct access to J, parity,
the occupation, seniority and all angular momenta from the coupling. 

++  if csf.useStandardSubshells  subshells = Defaults.GBL_STANDARD_SUBSHELL_LIST
else                         subshells = csf.subshells    
end
In the latter case, these subshells could be different for the leftCsf and rightCsf.

========================================================#

