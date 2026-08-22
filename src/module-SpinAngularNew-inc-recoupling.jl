
#
# Recoupling for one-particle operators of rank k > 0.
#
# This file is included from module-SpinAngularNew.jl and shares its namespace.
#

"""
`SpinAngularNew.openShells(csf::CsfR, subshells::Array{Subshell,1})`
    ... to list the indices of the subshells of `csf` that are neither empty nor closed. A closed subshell couples to
        J = 0 and contributes nothing to a tensor of rank k > 0, and an empty one contributes nothing at all, so only the
        open subshells carry the recoupling. A list indices::Array{Int64,1} is returned.
"""
function openShells(csf::CsfR, subshells::Array{Subshell,1})
    indices = Int64[]
    for  (i, sh) in enumerate(subshells)
        occ = csf.occupation[i];    full = Basics.subshell_2j(sh) + 1
        if  occ != 0  &&  occ != full    push!(indices, i)    end
    end

    return( indices )
end


"""
`SpinAngularNew.computeCoefficientsNonScalar(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the spin-angular coefficients of a one-particle operator of rank k > 0, in GRASP convention. A list
        coeffs::Array{Coefficient1p{ReducedKind},1} is returned.

        STAGE 1b: implemented for CSFs whose open subshells hold exactly ONE electron each, and at most two of them --
        the case in which the coupling tree is a product of single-electron angular momenta and no coefficient of
        fractional parentage is needed. A subshell holding two or more electrons requires the CFP machinery and RAISES
        rather than returning a number nobody has checked.
"""
function computeCoefficientsNonScalar(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                                      subshells::Array{Subshell,1})
    coeffs = Coefficient1p{ReducedKind}[]

    # A one-body operator changes the occupation of at most TWO subshells by exactly ONE electron each; any other pattern
    # is an EXACT ZERO. Everything else that is not yet implemented must RAISE -- returning an empty list where a real
    # coefficient exists would be a silent wrong answer, which is the failure mode this module exists to prevent.
    diffs = Int64[]
    for  i = 1:length(subshells)
        if  leftCsf.occupation[i] != rightCsf.occupation[i]    push!(diffs, i)    end
    end
    if  length(diffs) >  2                                                return( coeffs )   end
    if  length(diffs) == 1                                                return( coeffs )   end
    if  length(diffs) == 2
        if  abs(leftCsf.occupation[diffs[1]] - rightCsf.occupation[diffs[1]]) != 1  ||
            abs(leftCsf.occupation[diffs[2]] - rightCsf.occupation[diffs[2]]) != 1  return( coeffs )   end
        error("\n\nSpinAngularNew.computeCoefficientsNonScalar: a single-electron SUBSTITUTION at rank $(op.rank).\n" *
              ">>> Stage 1b computes rank > 0 for CSF pairs of equal occupation only. This pair moves one electron\n"   *
              ">>> between $(subshells[diffs[1]]) and $(subshells[diffs[2]]), for which GRASP does return a coefficient,\n" *
              ">>> so an empty list would be a silent wrong answer. Use SpinAngular.computeCoefficients for this pair.\n")
    end

    lOpen = openShells(leftCsf,  subshells);     rOpen = openShells(rightCsf, subshells)
    # A closed-shell-only CSF carries no rank-k > 0 one-particle coefficient at all.
    if  length(lOpen) == 0  ||  length(rOpen) == 0                        return( coeffs )   end
    if  length(lOpen) > 2
        error("\n\nSpinAngularNew.computeCoefficientsNonScalar: $(length(lOpen)) open subshells.\n" *
              ">>> Stage 1b handles at most two. Use SpinAngular.computeCoefficients.\n")
    end
    for  i in lOpen
        if  leftCsf.occupation[i] != 1  ||  rightCsf.occupation[i] != 1
            error("\n\nSpinAngularNew.computeCoefficientsNonScalar: a subshell holds more than one electron.\n" *
                  ">>> Stage 1b handles singly-occupied open subshells only; two or more electrons in one\n"     *
                  ">>> subshell need the coefficients of fractional parentage, which are not yet re-implemented\n" *
                  ">>> here. Use SpinAngular.computeCoefficients for such a CSF.\n")
        end
    end

    if       length(lOpen) == 1   coeffs = nonScalarOneOpenShell( op, leftCsf, rightCsf, subshells, lOpen[1])
    elseif   length(lOpen) == 2   coeffs = nonScalarTwoOpenShells(op, leftCsf, rightCsf, subshells, lOpen)
    end

    return( coeffs )
end


"""
`SpinAngularNew.nonScalarOneOpenShell(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64)`
    ... to compute the rank-k coefficient of a CSF with a single, singly-occupied open subshell. The total angular momentum
        is then that of the one electron, so the coupling tree is trivial. In GRASP's convention the coefficient is
        sqrt(2j_a+1)/sqrt(2J+1), and since J = j_a here the two roots cancel and the value is EXACTLY 1 -- which is what
        GRASP returns for every such pair. A list coeffs::Array{Coefficient1p{ReducedKind},1} is returned.
"""
function nonScalarOneOpenShell(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                               subshells::Array{Subshell,1}, ia::Int64)
    coeffs = Coefficient1p{ReducedKind}[]
    sh     = subshells[ia]
    if  !isAllowed1p(op, sh, sh)                                          return( coeffs )   end
    if  leftCsf.J != rightCsf.J                                           return( coeffs )   end

    # GRASP normalises by sqrt(2J_bra+1), not by sqrt(2k+1) -- MEASURED, see the header note. With a single open
    # subshell J is j_a itself, so the two roots cancel and the coefficient is exactly 1.
    ja2   = Basics.subshell_2j(sh)
    value = sqrt(ja2 + 1.0) / sqrt(Basics.twice(leftCsf.J) + 1.0)
    push!( coeffs, Coefficient1p{ReducedKind}(op.rank, sh, sh, value) )

    return( coeffs )
end


"""
`SpinAngularNew.nonScalarTwoOpenShells(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, open::Array{Int64,1})`
    ... to compute the rank-k coefficients of a CSF with two singly-occupied open subshells, coupled as |(j_a j_b) J>. The
        operator acts on one electron at a time, and the standard two-subsystem reduction of Edmonds gives, for the tensor
        acting on the FIRST subsystem,

            <(j_a j_b) J || T^(k)(1) || (j_a j_b) J'> = (-1)^(j_a+j_b+J'+k) sqrt((2J+1)(2J'+1)) {j_a J j_b; J' j_a k}
                                                        * <j_a || t^(k) || j_a>

        and the corresponding expression with j_a and j_b exchanged, and the phase carrying J rather than J', for the
        second. The coefficient returned is everything except the one-electron reduced matrix element, multiplied by the
        GRASP normalization sqrt(2j_a+1)/sqrt(2J_bra+1).

        THAT NORMALIZATION WAS MEASURED, NOT READ OFF. `oneparticlejj1.f90:61` divides by sqrt(2k+1), but GRASP's own
        recoupling factor carries a further J-dependence, so the NET convention is the one above. Calibrating the textbook
        two-subsystem formula against GRASP gave a ratio sqrt((2J_bra+1)/(2k+1)) that fitted 15 of 15 coefficients across
        ranks 1, 2 and 3 and two values of J -- and the corrected form then predicted, out of sample, the exact 1.0 that
        GRASP returns for every single-open-shell pair. A convention confirmed on one case is a coincidence; this one is
        confirmed on fifteen and then tested on a set it was not fitted to. A list coeffs::Array{Coefficient1p{ReducedKind},1}
        is returned.
"""
function nonScalarTwoOpenShells(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                                subshells::Array{Subshell,1}, open::Array{Int64,1})
    coeffs = Coefficient1p{ReducedKind}[]
    ia, ib = open[1], open[2]
    sha, shb = subshells[ia], subshells[ib]
    ja  = AngularJ64( Basics.subshell_2j(sha)//2 );    jb = AngularJ64( Basics.subshell_2j(shb)//2 )
    Jf  = leftCsf.J;                                   Ji = rightCsf.J
    k   = op.rank
    wJf = Basics.twice(Jf) + 1.0;                      wJi = Basics.twice(Ji) + 1.0

    # ... the operator acting on the electron in subshell ia
    if  isAllowed1p(op, sha, sha)
        phase = (-1)^Int64( (Basics.twice(ja) + Basics.twice(jb) + Basics.twice(Ji))//2 + k )
        wa    = phase * sqrt(wJf * wJi) * AngularMomentum.Wigner_6j(ja, Jf, jb, Ji, ja, AngularJ64(k))
        wa    = wa * sqrt(Basics.twice(ja) + 1.0) / sqrt(wJf)
        if  abs(wa) > 0.0    push!( coeffs, Coefficient1p{ReducedKind}(k, sha, sha, wa) )    end
    end

    # ... and on the electron in subshell ib
    if  isAllowed1p(op, shb, shb)
        phase = (-1)^Int64( (Basics.twice(ja) + Basics.twice(jb) + Basics.twice(Jf))//2 + k )
        wa    = phase * sqrt(wJf * wJi) * AngularMomentum.Wigner_6j(jb, Jf, ja, Ji, jb, AngularJ64(k))
        wa    = wa * sqrt(Basics.twice(jb) + 1.0) / sqrt(wJf)
        if  abs(wa) > 0.0    push!( coeffs, Coefficient1p{ReducedKind}(k, shb, shb, wa) )    end
    end

    return( coeffs )
end
