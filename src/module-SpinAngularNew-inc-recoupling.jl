
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
              ">>> Rank > 0 is computed for CSF pairs of EQUAL OCCUPATION only. This pair moves one electron\n"   *
              ">>> between $(subshells[diffs[1]]) and $(subshells[diffs[2]]), for which GRASP does return a coefficient,\n" *
              ">>> so an empty list would be a silent wrong answer. Use SpinAngular.computeCoefficients for this pair.\n")
    end

    # The occupations are equal from here on, so both CSFs have the same open subshells. A CSF of closed subshells only
    # carries no rank-k > 0 one-particle coefficient at all.
    lOpen = openShells(leftCsf, subshells)
    if  length(lOpen) == 0                                                return( coeffs )   end

    coeffs = nonScalarGeneral(op, leftCsf, rightCsf, subshells, lOpen)

    return( coeffs )
end


"""
`SpinAngularNew.shellReducedW(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64, Jket::AngularJ64, kj::Int64)`
    ... to compute the reduced matrix element of the shell operator W^(kj) = (a^+ x a~)^(kj) within a single subshell j^N,

            <j^N v J || W^(kj) || j^N v' J'>

        assembled from the quasispin representation. The coefficients of fractional parentage themselves are NOT
        re-derived here: `SpinAngular.completelyReducedWkk` holds G. Gaigalas's completely reduced (j Q J ||| W^(kq kj) |||
        j Q' J') as exact data -- stored as [sign, num, den] and returned as sign*sqrt(num/den) -- and re-typing a correct
        table would add risk and nothing else. What is re-implemented is the ASSEMBLY: the quasispin Wigner-Eckart step
        that turns the completely reduced element into the one for a shell of N electrons.

        The quasispin rank follows from the angular rank, kq = 1 for even kj and kq = 0 for odd kj, and the projection is
        M_Q = (N - (2j+1)/2)/2 on both sides since the operator conserves particle number. For kj = 0 the result is the
        closed form -N sqrt((2J+1)/(2j+1)).

        VERIFIED against `SpinAngular.irreducibleTensor(SchemeEta_W(), ...)` to ratio 1.000000 on every case tested, which
        isolates this step from the outer normalization. A value::Float64 is returned.
"""
function shellReducedW(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64, Jket::AngularJ64,
                       kj::Int64)
    SA = JenaAtomicCalculator.SpinAngular
    Qb = SA.qshellTermQ(j, senBra);           Qk = SA.qshellTermQ(j, senKet)
    MQ = SA.qshellTermM(j, N)
    ib = SA.getTermNumber(j, N, Qb, Jbra);    ik = SA.getTermNumber(j, N, Qk, Jket)

    if  kj == 0
        if  ib != ik    return( 0.0 )    end
        return( -N * sqrt( (Basics.twice(Jbra) + 1.0) / (Basics.twice(j) + 1.0) ) )
    end

    kq = iseven(kj) ? 1 : 0
    if  AngularMomentum.triangularDelta(Qb, AngularJ64(kq), Qk) == 0     return( 0.0 )   end
    wa = AngularMomentum.ClebschGordan(Qk, MQ, AngularJ64(kq), AngularM64(0), Qb, MQ)
    wa = wa * SA.completelyReducedWkk(ib, ik, kq, kj)
    wa = wa / sqrt( (Basics.twice(Qb) + 1.0) * 2.0 )

    return( wa )
end


"""
`SpinAngularNew.chainRecoupling(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, k::Int64)`
    ... to compute the recoupling factor for a one-particle tensor of rank k acting on the subshell `ip` of a CSF whose
        subshells are coupled as a chain X_1 = J_1, X_q = X_{q-1} x J_q, X_n = J.

        The tensor is peeled outwards, one subshell at a time. For every q > ip the operator sits in the FIRST subsystem
        with J_q as spectator,

            (-1)^(X_{q-1}+J_q+X'_q+k) sqrt((2X_q+1)(2X'_q+1)) { X_{q-1} X_q J_q ; X'_q X'_{q-1} k }

        and at q = ip it sits in the SECOND subsystem with X_{ip-1} as spectator,

            (-1)^(X_{ip-1}+J'_ip+X_ip+k) sqrt((2X_ip+1)(2X'_ip+1)) { J_ip X_ip X_{ip-1} ; X'_ip J'_ip k }

        with X_0 = 0. Both limits that were already verified fall out of this, which is why it replaces them rather than
        sitting beside them: with every other subshell closed each factor collapses to 1, giving the single-open-subshell
        result; and with two singly-occupied subshells the two expressions above reduce term for term to the Edmonds
        two-subsystem formulae. A value::Float64 is returned.
"""
function chainRecoupling(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, k::Int64)
    nw = length(leftCsf.occupation)
    wa = 1.0
    kJ = AngularJ64(k)

    # ... the outer subshells, q > ip, operator in the first subsystem
    for  q = ip+1:nw
        Xqm = (q == 1) ? AngularJ64(0) : leftCsf.subshellX[q-1]
        Ypm = (q == 1) ? AngularJ64(0) : rightCsf.subshellX[q-1]
        Xq  = leftCsf.subshellX[q];              Yq = rightCsf.subshellX[q]
        Jq  = leftCsf.subshellJ[q]
        if  Jq != rightCsf.subshellJ[q]          return( 0.0 )   end
        ph  = Int64( (Basics.twice(Xqm) + Basics.twice(Jq) + Basics.twice(Yq))//2 ) + k
        wa  = wa * (-1)^ph * sqrt((Basics.twice(Xq)+1.0)*(Basics.twice(Yq)+1.0)) *
                   AngularMomentum.Wigner_6j(Xqm, Xq, Jq, Yq, Ypm, kJ)
        if  wa == 0.0    return( 0.0 )   end
    end

    # ... and the acting subshell itself, operator in the second subsystem
    Xpm = (ip == 1) ? AngularJ64(0) : leftCsf.subshellX[ip-1]
    Xp  = leftCsf.subshellX[ip];             Yp  = rightCsf.subshellX[ip]
    Jp  = leftCsf.subshellJ[ip];             Jpp = rightCsf.subshellJ[ip]
    ph  = Int64( (Basics.twice(Xpm) + Basics.twice(Jpp) + Basics.twice(Xp))//2 ) + k
    wa  = wa * (-1)^ph * sqrt((Basics.twice(Xp)+1.0)*(Basics.twice(Yp)+1.0)) *
               AngularMomentum.Wigner_6j(Jp, Xp, Xpm, Yp, Jpp, kJ)

    return( wa )
end


"""
`SpinAngularNew.nonScalarGeneral(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, openList::Array{Int64,1})`
    ... to compute the rank-k coefficients of a CSF pair of equal occupation with ANY number of open subshells, each
        holding any number of electrons. Each open subshell contributes in turn, its own shell matrix element from
        `SpinAngularNew.shellReducedW` and its place in the coupling tree from `SpinAngularNew.chainRecoupling`:

            T^(k)(a,a)  =  - R_chain * <j^N v J_a || W^(k) || j^N v' J'_a> * sqrt(2j_a+1) / ( sqrt(2k+1) sqrt(2J_bra+1) )

        This is the same expression as for a single open subshell, with the recoupling factor no longer equal to one. A
        list coeffs::Array{Coefficient1p{ReducedKind},1} is returned.
"""
function nonScalarGeneral(op::SpinAngularNew.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                          subshells::Array{Subshell,1}, openList::Array{Int64,1})
    coeffs = Coefficient1p{ReducedKind}[]
    k      = op.rank

    for  ip in openList
        sh = subshells[ip]
        if  !isAllowed1p(op, sh, sh)    continue    end
        # every OTHER subshell must be unchanged in its own coupling, or the two CSFs are orthogonal
        ok = true
        for  i = 1:length(subshells)
            if  i == ip    continue    end
            if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
                ok = false;    break
            end
        end
        if  !ok    continue    end

        j  = AngularJ64( Basics.subshell_2j(sh)//2 )
        N  = leftCsf.occupation[ip]
        wS = shellReducedW(j, N, leftCsf.seniorityNr[ip], leftCsf.subshellJ[ip],
                                 rightCsf.seniorityNr[ip], rightCsf.subshellJ[ip], k)
        if  wS == 0.0    continue    end
        wR = chainRecoupling(leftCsf, rightCsf, ip, k)
        if  wR == 0.0    continue    end

        value = -wR * wS * sqrt(Basics.twice(j) + 1.0) /
                     ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) )
        if  value != 0.0    push!( coeffs, Coefficient1p{ReducedKind}(k, sh, sh, value) )    end
    end

    return( coeffs )
end
