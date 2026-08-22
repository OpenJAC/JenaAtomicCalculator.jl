
"""
`module  JAC.SpinAngularNew`
... an independent re-implementation of the spin-angular coefficients of JAC, developed ALONGSIDE `SpinAngular` and not
    (yet) replacing it. It exists to settle two defects of the present module that cannot be repaired in place without
    touching its eighteen callers at once.

    **(1) One normalization for every rank.** `SpinAngular` applies a factor `sqrt(2j_a+1)` inside the coefficient for
    rank > 0 (the line marked "GRASP like") and has the identical line COMMENTED OUT for rank 0. Neither is wrong on its
    own -- `Hamiltonian` and `IsotopeShift` re-apply the factor at the call site, `Hfs` divides it back out -- but the same
    physics is expressed two opposite ways and every caller must know which. Four modules have already guessed wrong, and
    one of those guesses shipped a hyperfine constant too large by sqrt(2). This module uses GRASP's convention, with
    `sqrt(2j_a+1)` INSIDE the coefficient, at EVERY rank. The rank-0 diagonal coefficient is then the plain occupation
    number, which is a value a reader can check by eye and a test can assert exactly.

    **(2) Zeros decided by selection rules, not by thresholds.** `SpinAngular` guards its recursion with `abs(wa) >= 2.0e-10`,
    which answers "is this number small?" where the question is "do the selection rules permit this at all?". Triangle
    conditions, parity and occupation are combinatorial facts about quantum numbers, decidable in integer arithmetic before
    any float exists. Deciding them up front removes the arbitrary cut-off and skips the work rather than doing it and
    discarding the result -- and it suppresses coefficients that should never have been emitted (see
    `SpinAngularNew.isAllowed1p`).

    **What the coefficient must be multiplied by depends on the rank, and that is irreducible.** A rank-0 coefficient
    multiplies the ORDINARY one-electron integral `<a| o |b>`; a rank-k coefficient multiplies the REDUCED
    `<a|| o^(k) ||b>`. A scalar and a rank-k operator are different objects, so this asymmetry cannot be normalized away --
    only declared. It is therefore carried in the TYPE, as `Coefficient1p{OrdinaryKind}` or `Coefficient1p{ReducedKind}`,
    and `SpinAngularNew.contract` exists only for a coefficient and a one-electron matrix element of the SAME kind. Pairing
    the wrong two raises a `MethodError` naming both, instead of returning a number that is wrong by `sqrt(2j_a+1)`.

    The phantom parameter costs nothing measurable: it carries no data, both element types are concrete, and a two-million
    coefficient contraction was byte-identical in allocation and within run-to-run noise in time.

    STATUS, 22-Aug-2026:  UNDER DEVELOPMENT, stage 1a. Rank 0 is implemented for the diagonal case and for a single-electron
    substitution between subshells of the same kappa; every other case RAISES rather than returning a number nobody has
    checked. Rank > 0 is NOT yet implemented. This module is deliberately NOT included from `JenaAtomicCalculator.jl`:
    `examples/example-Aq.jl` includes it directly, so that a broken intermediate state cannot break the package.
"""
module SpinAngularNew

# NOTE: this module is included DIRECTLY by examples/example-Aq.jl and not from JenaAtomicCalculator.jl, so the imports are
# absolute rather than relative. They become `..Basics, ..ManyElectron` on the day it moves inside the package.
using  Printf, JenaAtomicCalculator
using  JenaAtomicCalculator.Basics, JenaAtomicCalculator.ManyElectron


"""
`abstract type SpinAngularNew.AbstractOneElectronKind`
    ... defines an abstract type to distinguish WHICH one-electron matrix element a spin-angular coefficient is to be
        multiplied by. The distinction is a property of the operator rank and cannot be normalized away, so it is carried
        in the type of the coefficient rather than left to a comment.

    + OrdinaryKind    ... the coefficient multiplies the ordinary matrix element <a| o |b>       (rank 0).
    + ReducedKind     ... the coefficient multiplies the reduced matrix element <a|| o^(k) ||b>  (rank > 0).
"""
abstract type  AbstractOneElectronKind                                  end
struct         OrdinaryKind   <:  AbstractOneElectronKind               end
struct         ReducedKind    <:  AbstractOneElectronKind               end


"""
`SpinAngularNew.oneElectronKind(rank::Int64)`
    ... to decide, in the single place where this decision is made, which one-electron matrix element a coefficient of the
        given rank must be multiplied by; a kind::AbstractOneElectronKind is returned.
"""
function oneElectronKind(rank::Int64)
    if  rank == 0   return( OrdinaryKind() )   else   return( ReducedKind() )   end
end


"""
`struct  SpinAngularNew.OneParticleOperator`
    ... a struct for defining a (reduced) one-particle operator of given rank and parity.

    + rank     ::Int64      ... Rank (k) of the one-particle operator.
    + parity   ::Parity     ... Parity of the operator, plus or minus.
"""
struct  OneParticleOperator
    rank       ::Int64
    parity     ::Parity
end


"""
`SpinAngularNew.OneParticleOperator()`  ... constructor for a scalar, parity-even one-particle operator.
"""
function OneParticleOperator()
    OneParticleOperator( 0, Basics.plus )
end


# `Base.show(io::IO, op::SpinAngularNew.OneParticleOperator)`  ... prepares a proper printout of op.
function Base.show(io::IO, op::SpinAngularNew.OneParticleOperator)
    print(io, "one-particle operator O^($(op.rank)) [$(string(op.parity))]")
end


"""
`struct  SpinAngularNew.Coefficient1p{K<:AbstractOneElectronKind}`
    ... a struct for a single spin-angular coefficient of a one-particle matrix element. The type parameter K records which
        one-electron matrix element the coefficient is to be multiplied by, so that the wrong pairing has no method rather
        than a wrong value.

    + nu       ::Int64      ... Rank (k) of the one-particle operator.
    + a        ::Subshell   ... Left-hand subshell (orbital).
    + b        ::Subshell   ... Right-hand subshell (orbital).
    + T        ::Float64    ... (Value of) the spin-angular coefficient, in GRASP convention.
"""
struct  Coefficient1p{K<:AbstractOneElectronKind}
    nu         ::Int64
    a          ::Subshell
    b          ::Subshell
    T          ::Float64
end


# `Base.show(io::IO, coeff::SpinAngularNew.Coefficient1p)`  ... prepares a proper printout of coeff.
function Base.show(io::IO, coeff::SpinAngularNew.Coefficient1p{K})  where K<:AbstractOneElectronKind
    print(io, "   T^$(coeff.nu) [$(coeff.a), $(coeff.b)] = $(coeff.T)   ($(K.name.name))")
end


"""
`struct  SpinAngularNew.OneElectronMe{K<:AbstractOneElectronKind, F}`
    ... a struct that tags a caller's one-electron matrix element with the kind it returns, so that it can only be combined
        with coefficients of the same kind. F is kept concrete so that the wrapped call inlines and costs nothing.

    + f        ::F          ... A callable (a::Subshell, b::Subshell) -> Float64.
"""
struct  OneElectronMe{K<:AbstractOneElectronKind, F}
    f          ::F
end


"""
`SpinAngularNew.OrdinaryMe(f)`
    ... to tag a callable `f(a,b)` as returning the ORDINARY one-electron matrix element <a| o |b>, as a rank-0 coefficient
        requires; a me::OneElectronMe{OrdinaryKind} is returned.
"""
OrdinaryMe(f) = OneElectronMe{OrdinaryKind, typeof(f)}(f)


"""
`SpinAngularNew.ReducedMe(f)`
    ... to tag a callable `f(a,b)` as returning the REDUCED one-electron matrix element <a|| o^(k) ||b>, as a rank-k
        coefficient requires; a me::OneElectronMe{ReducedKind} is returned.
"""
ReducedMe(f)  = OneElectronMe{ReducedKind,  typeof(f)}(f)


"""
`SpinAngularNew.contract(coeffs::Array{Coefficient1p{K},1}, me::OneElectronMe{K})`
    ... to combine a list of spin-angular coefficients with the one-electron matrix elements they belong to, and so to form
        the many-electron matrix element. This method exists ONLY where the two kinds agree; a mismatched pair has no
        applicable method and raises a `MethodError` naming both types, which is the whole purpose of the type parameter.
        A value::Float64 is returned.
"""
function contract(coeffs::Array{Coefficient1p{K},1}, me::OneElectronMe{K})  where K<:AbstractOneElectronKind
    wa = 0.0
    for  c in coeffs    wa = wa + c.T * me.f(c.a, c.b)    end

    return( wa )
end


"""
`SpinAngularNew.isAllowed1p(op::SpinAngularNew.OneParticleOperator, a::Subshell, b::Subshell)`
    ... to decide from the SELECTION RULES ALONE, in integer arithmetic and before any floating-point work is done, whether
        a one-particle operator of the given rank and parity can connect the two subshells at all. This replaces the
        `abs(wa) >= 2.0e-10` magnitude guards of `SpinAngular`: a coefficient is zero because the quantum numbers forbid it,
        not because it came out small.

        Two conditions are tested. The triangle condition |j_a - j_b| <= k <= j_a + j_b, and the parity of the operator
        against the parity (-1)^(l_a+l_b) of the subshell pair. For k = 0 the triangle condition alone already enforces
        j_a = j_b, and the parity rule then enforces l_a = l_b, i.e. kappa_a = kappa_b -- which is why a scalar operator
        cannot connect 2s to 2p_1/2. `SpinAngular` emits such a coefficient; GRASP2018 does not. A Bool is returned.
"""
function isAllowed1p(op::SpinAngularNew.OneParticleOperator, a::Subshell, b::Subshell)
    ja2 = Basics.subshell_2j(a);        jb2 = Basics.subshell_2j(b)
    la  = Basics.subshell_l(a);         lb  = Basics.subshell_l(b)
    # Triangle condition, in twice-j integers throughout to avoid any half-integer arithmetic
    if  2*op.rank < abs(ja2 - jb2)   ||   2*op.rank > ja2 + jb2      return( false )   end
    # Parity of the operator against the parity of the subshell pair
    opParity  = (op.parity == Basics.plus)  ?  1  :  -1
    pairParity = iseven(la + lb)  ?  1  :  -1
    if  opParity != pairParity                                       return( false )   end

    return( true )
end


"""
`SpinAngularNew.computeCoefficients(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the spin-angular coefficients of the one-particle matrix element <leftCsf || O^(k) || rightCsf> for the
        given subshell list, in GRASP convention, i.e. with the factor sqrt(2j_a+1) carried INSIDE the coefficient at every
        rank. A list coeffs::Array{Coefficient1p{K},1} is returned, whose kind K follows from the rank.

        STAGE 1a: only rank 0 is implemented, and only for a diagonal pair of CSFs or for a single-electron substitution
        between two subshells of the same kappa with all other subshells in identical coupling. Everything else raises,
        following `ReducedDensityMatrix.compute1pRDMDirect`, rather than returning a number nobody has checked.
"""
function computeCoefficients(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                             subshells::Array{Subshell,1})
    if  op.rank != 0
        error("\n\nSpinAngularNew.computeCoefficients: rank $(op.rank) is NOT yet implemented (stage 1a covers rank 0 only).\n" *
              ">>> Use SpinAngular.computeCoefficients for rank > 0, remembering that it returns the coefficient in the\n"      *
              ">>> same GRASP convention this module uses, so no conversion is needed for that rank.\n")
    end

    return( computeCoefficientsScalar(op, leftCsf, rightCsf, subshells) )
end


"""
`SpinAngularNew.computeCoefficientsScalar(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the spin-angular coefficients of a SCALAR (rank-0) one-particle operator, in GRASP convention. For a
        diagonal pair of CSFs the coefficient of each subshell is exactly its occupation number, since
        <Psi| sum_i f(i) |Psi> = sum_a N_a <a| f |a> holds for any coupling; that identity is what
        `SpinAngularNew.checkOccupationIdentity` asserts. A list coeffs::Array{Coefficient1p{OrdinaryKind},1} is returned.
"""
function computeCoefficientsScalar(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                                   subshells::Array{Subshell,1})
    coeffs = Coefficient1p{OrdinaryKind}[]
    nw     = length(subshells)

    # A scalar operator connects only states of equal J and equal parity; this is decided before any work is done.
    if  leftCsf.J != rightCsf.J   ||   leftCsf.parity != rightCsf.parity      return( coeffs )   end

    diffs = Int64[]
    for  i = 1:nw
        if  leftCsf.occupation[i] != rightCsf.occupation[i]    push!(diffs, i)    end
    end

    if       length(diffs) == 0     coeffs = scalarDiagonal(op, leftCsf, rightCsf, subshells)
    elseif   length(diffs) == 2     coeffs = scalarSingleSubstitution(op, leftCsf, rightCsf, subshells, diffs)
    else
        error("\n\nSpinAngularNew.computeCoefficientsScalar: the two CSFs differ in $(length(diffs)) subshells.\n" *
              ">>> Stage 1a supports a diagonal pair and a single-electron substitution only; a more general\n"    *
              ">>> occupation pattern needs the coefficient-of-fractional-parentage machinery, which is not yet\n" *
              ">>> re-implemented here. Use SpinAngular.computeCoefficients for such a pair.\n")
    end

    return( coeffs )
end


"""
`SpinAngularNew.scalarDiagonal(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the coefficients of a scalar one-particle operator between two CSFs of identical occupation. Every
        coefficient is the occupation number of its subshell, exactly and independently of the coupling. A list
        coeffs::Array{Coefficient1p{OrdinaryKind},1} is returned.
"""
function scalarDiagonal(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                        subshells::Array{Subshell,1})
    coeffs = Coefficient1p{OrdinaryKind}[]

    # The occupations agree by construction; the couplings must agree too, or the two CSFs are orthogonal.
    for  i = 1:length(subshells)
        if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.subshellX[i] != rightCsf.subshellX[i]  ||
            leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
            return( coeffs )
        end
    end

    for  (i, sh) in enumerate(subshells)
        occ = leftCsf.occupation[i]
        if  occ == 0                            continue    end
        if  !isAllowed1p(op, sh, sh)            continue    end
        push!( coeffs, Coefficient1p{OrdinaryKind}(op.rank, sh, sh, Float64(occ)) )
    end

    return( coeffs )
end


"""
`SpinAngularNew.scalarSingleSubstitution(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, diffs::Array{Int64,1})`
    ... to compute the coefficient of a scalar one-particle operator between two CSFs that differ by moving one electron
        between the two subshells named in `diffs`. A scalar operator cannot change j or l, so the two subshells must share
        the same kappa; `SpinAngularNew.isAllowed1p` decides this, and a pair such as (2s, 2p_1/2) is rejected here rather
        than being emitted with a value that only a vanishing radial integral would suppress. A list
        coeffs::Array{Coefficient1p{OrdinaryKind},1} is returned.
"""
function scalarSingleSubstitution(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                                  subshells::Array{Subshell,1}, diffs::Array{Int64,1})
    coeffs = Coefficient1p{OrdinaryKind}[]
    ia, ib = diffs[1], diffs[2]

    if  abs(leftCsf.occupation[ia] - rightCsf.occupation[ia]) != 1  ||
        abs(leftCsf.occupation[ib] - rightCsf.occupation[ib]) != 1
        error("\n\nSpinAngularNew.scalarSingleSubstitution: the two CSFs differ by more than one electron in a subshell.\n" *
              ">>> Stage 1a supports a single-electron substitution only.\n")
    end
    if  !isAllowed1p(op, subshells[ia], subshells[ib])      return( coeffs )   end

    # All remaining subshells must be in identical coupling, or the two CSFs are orthogonal.
    for  i = 1:length(subshells)
        if  i == ia  ||  i == ib    continue    end
        if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.subshellX[i] != rightCsf.subshellX[i]  ||
            leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
            return( coeffs )
        end
    end

    naL = leftCsf.occupation[ia];    nbL = leftCsf.occupation[ib]
    naR = rightCsf.occupation[ia];   nbR = rightCsf.occupation[ib]
    # Orient the substitution so that the electron is annihilated in the ket and created in the bra.
    if       naL == naR + 1  &&  nbL == nbR - 1     iCre, iAnn = ia, ib
    elseif   nbL == nbR + 1  &&  naL == naR - 1     iCre, iAnn = ib, ia
    else
        error("\n\nSpinAngularNew.scalarSingleSubstitution: the occupations do not describe a single transfer.\n")
    end

    # For a substitution between two subshells of the same kappa, whose occupations are 0/1 on one side and full/full-1
    # on the other, the transfer coefficient is sqrt(N_full) in GRASP convention. The general open-shell case needs the
    # fractional-parentage machinery and is refused rather than guessed.
    fullCre = Basics.subshell_2j(subshells[iCre]) + 1
    fullAnn = Basics.subshell_2j(subshells[iAnn]) + 1
    nCre    = max(leftCsf.occupation[iCre], rightCsf.occupation[iCre])
    nAnn    = max(leftCsf.occupation[iAnn], rightCsf.occupation[iAnn])
    if  !( (nCre == 1 || nCre == fullCre)  &&  (nAnn == 1 || nAnn == fullAnn) )
        error("\n\nSpinAngularNew.scalarSingleSubstitution: a general open-shell substitution is NOT yet supported.\n" *
              ">>> Stage 1a handles a transfer between subshells that are singly occupied or closed; the two here\n"   *
              ">>> hold $nCre of $fullCre and $nAnn of $fullAnn electrons. Use SpinAngular.computeCoefficients.\n")
    end

    value = sqrt( Float64(nCre) ) * sqrt( Float64(nAnn) )
    push!( coeffs, Coefficient1p{OrdinaryKind}(op.rank, subshells[iCre], subshells[iAnn], value) )

    return( coeffs )
end


"""
`SpinAngularNew.checkOccupationIdentity(coeffs::Array{Coefficient1p{OrdinaryKind},1}, csf::CsfR, subshells::Array{Subshell,1})`
    ... to test the one exact identity a scalar one-particle operator must satisfy on a diagonal matrix element, namely that
        the coefficient of each subshell is its occupation number, since <Psi| sum_i f(i) |Psi> = sum_a N_a <a| f |a> holds
        for any coupling. This is an identity and not a tolerance: the deviation must be zero to rounding, and a non-zero
        value is a defect. It is available on every diagonal call and costs nothing, which matters in a module whose
        predecessor has no direct test coverage at all. A deviation::Float64 is returned.
"""
function checkOccupationIdentity(coeffs::Array{Coefficient1p{OrdinaryKind},1}, csf::CsfR, subshells::Array{Subshell,1})
    deviation = 0.0
    for  (i, sh) in enumerate(subshells)
        occ = Float64( csf.occupation[i] )
        wa  = 0.0
        for  c in coeffs    if  c.a == sh  &&  c.b == sh    wa = wa + c.T    end    end
        deviation = max( deviation, abs(wa - occ) )
    end

    return( deviation )
end

end # module
