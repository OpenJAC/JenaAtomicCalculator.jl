
#
# Recoupling for one-particle operators of rank k > 0.
#
# This file is included from module-SpinAngular.jl and shares its namespace.
#

"""
`SpinAngular.openShells(csf::CsfR, subshells::Array{Subshell,1})`
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
`SpinAngular.computeCoefficientsNonScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the spin-angular coefficients of a one-particle operator of rank k > 0, in GRASP convention. A list
        coeffs::Array{Coefficient1p{ReducedKind},1} is returned.

        STAGE 1b: implemented for CSFs whose open subshells hold exactly ONE electron each, and at most two of them --
        the case in which the coupling tree is a product of single-electron angular momenta and no coefficient of
        fractional parentage is needed. A subshell holding two or more electrons requires the CFP machinery and RAISES
        rather than returning a number nobody has checked.
"""
function computeCoefficientsNonScalar(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
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
        iCre = leftCsf.occupation[diffs[1]] > rightCsf.occupation[diffs[1]] ? diffs[1] : diffs[2]
        iAnn = (iCre == diffs[1]) ? diffs[2] : diffs[1]
        return( nonScalarSubstitution(op, leftCsf, rightCsf, subshells, iCre, iAnn) )
    end

    # The occupations are equal from here on, so both CSFs have the same open subshells. A CSF of closed subshells only
    # carries no rank-k > 0 one-particle coefficient at all.
    lOpen = openShells(leftCsf, subshells)
    if  length(lOpen) == 0                                                return( coeffs )   end

    coeffs = nonScalarGeneral(op, leftCsf, rightCsf, subshells, lOpen)

    return( coeffs )
end


"""
`SpinAngular.shellReducedW(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64, Jket::AngularJ64, kj::Int64)`
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

        VERIFIED against `SpinAngularGaigalas.irreducibleTensor(SchemeEta_W(), ...)` to ratio 1.000000 on every case tested, which
        isolates this step from the outer normalization. A value::Float64 is returned.
"""
#
# The memo for `shellReducedW`. A MODULE-LEVEL cache is used here deliberately, and the reasoning is worth stating
# because JAC's own precedent, `InteractionStrength.XLCache`, argues the opposite -- that a cache should be a
# parameter and not a global.
#
# THAT ARGUMENT TURNS ON A HAZARD THAT DOES NOT EXIST HERE. The XL key holds subshell LABELS, and a label identifies
# an orbital only within one basis: two bases can both contain a "2s_1/2" whose radial functions differ entirely, so
# a global store can outlive the basis whose labels it uses and hand back numbers for the wrong orbitals. The key
# below holds no labels at all -- only j, the occupation, the two seniorities, the two shell J values and the rank.
# Those are pure angular-momentum quantum numbers, and the value is a pure angular quantity: the reduced W of a
# 3d_5/2 shell with N = 3, v = 3, J = 9/2 is the same number in every basis, every atom, and every calculation.
# There is nothing for it to outlive.
#
# The key space is also bounded -- j <= 9/2, N <= 10, and the ranks that survive the triangle conditions -- so the
# store cannot grow without limit. MEASURED before being written: over three configurations the same key recurs 9.7
# times on average (882 invocations, 91 distinct), whereas caching the RECOUPLING instead would gain nothing at all,
# its repeat factor being 1.04.
#
const SHELL_W_CACHE = Dict{NTuple{7,Int64}, Float64}()


# The single-operator matrix elements are memoised for the same reason as the rank-k ones, and it matters more: a
# profile of the one-electron-move sweep is dominated not by the recoupling but by EXACT RATIONAL ARITHMETIC -- BigInt
# allocation inside the Wigner-symbol package, reached through the Clebsch-Gordan in `shellReducedA`. The arguments are
# a handful of small quantum numbers and the same ones recur across every CSF pair of a calculation.
const SHELL_A_CACHE = Dict{NTuple{8,Int64}, Float64}()


"""
`SpinAngular.clearCaches()`
    ... to empty the memo of `SpinAngular.shellReducedW`. Not needed for correctness -- the cached quantities are
        basis-independent -- but useful for timing a cold run. Returns the number of entries discarded.
"""
function clearCaches()
    n = length(SHELL_W_CACHE) + length(PARTNER_CACHE) + length(SHELL_A_CACHE)
    empty!(SHELL_W_CACHE);    empty!(PARTNER_CACHE);    empty!(SHELL_A_CACHE)
    empty!(PARTNER_CACHE)

    return( n )
end


"""
`SpinAngular.shellReducedW(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64,`
                          `Jket::AngularJ64, kj::Int64)`
    ... gives the reduced matrix element of the W^(kj) double tensor within one subshell j^N, between the states of
        seniority senBra/Jbra and senKet/Jket. **This is the MEMOISED wrapper**: it returns a cached value where one
        exists and otherwise calls `shellReducedWUncached` and stores the result under the key of all seven
        arguments in doubled form. A value::Float64 is returned. The cache `SHELL_W_CACHE` is session-lived and
        unbounded; see Rule 19 on whether it should be scoped.
"""
function shellReducedW(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64, Jket::AngularJ64,
                       kj::Int64)
    key = (Basics.twice(j), N, senBra, Basics.twice(Jbra), senKet, Basics.twice(Jket), kj)
    haskey(SHELL_W_CACHE, key)  &&  return( SHELL_W_CACHE[key] )
    wa  = shellReducedWUncached(j, N, senBra, Jbra, senKet, Jket, kj)
    SHELL_W_CACHE[key] = wa

    return( wa )
end


"""
`SpinAngular.shellReducedWUncached(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64, Jket::AngularJ64, kj::Int64)`
    ... the body of `SpinAngular.shellReducedW`, without the memo. Kept separate so that the cache can be tested
        against it directly rather than trusted. A value::Float64 is returned.
"""
function shellReducedWUncached(j::AngularJ64, N::Int64, senBra::Int64, Jbra::AngularJ64, senKet::Int64,
                               Jket::AngularJ64, kj::Int64)
    SA = SpinAngularTables
    # AN EMPTY OR SINGLY OCCUPIED SUBSHELL NEEDS NO COEFFICIENT OF FRACTIONAL PARENTAGE, and taking those two cases
    # here is what lets this routine work above j = 9/2, where the quasispin tables it otherwise consults stop.
    # A singly occupied subshell has exactly one state, seniority 1 with J = j, and measured against the tables
    # wherever they do reach -- j = 1/2 ... 9/2, kj = 0 ... 3 -- the element is -sqrt(2 kj + 1), INDEPENDENT of j:
    # -1, -sqrt(3), -sqrt(5), -sqrt(7). See the companion closed form in shellReducedA.
    if  N == 0                                                                          return( 0.0 )   end
    if  N == 1
        if  senBra != 1  ||  senKet != 1  ||  Jbra != j  ||  Jket != j                   return( 0.0 )   end
        return( -sqrt(2.0*kj + 1.0) )
    end
    # Beyond those, the tables are needed and their reach must be stated rather than silently returned as a zero
    # or thrown as a BoundsError from inside SpinAngular.
    if  Basics.twice(j) > 9
        error("\n\nSpinAngular.shellReducedW: the quasispin/CFP tables reach j <= 9/2 and j = $(j) was asked "   *
              "for\n>>> with occupation N = $N, which is beyond the closed forms for an empty or singly occupied "  *
              "subshell.\n")
    end
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
`SpinAngular.chainRecoupling(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, kJ::AngularJ64; top::Int64 = 0)`
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
function chainRecoupling(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, kJ::AngularJ64; top::Int64 = 0)
    nw = (top == 0) ? length(leftCsf.occupation) : top
    wa = peelRange(leftCsf, rightCsf, ip+1, nw, kJ)
    if  wa == 0.0                                                         return( 0.0 )   end

    return( wa * actingFactor(leftCsf, rightCsf, ip, kJ) )
end


"""
`SpinAngular.nonScalarGeneral(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, openList::Array{Int64,1})`
    ... to compute the rank-k coefficients of a CSF pair of equal occupation with ANY number of open subshells, each
        holding any number of electrons. Each open subshell contributes in turn, its own shell matrix element from
        `SpinAngular.shellReducedW` and its place in the coupling tree from `SpinAngular.chainRecoupling`:

            T^(k)(a,a)  =  - R_chain * <j^N v J_a || W^(k) || j^N v' J'_a> * sqrt(2j_a+1) / ( sqrt(2k+1) sqrt(2J_bra+1) )

        This is the same expression as for a single open subshell, with the recoupling factor no longer equal to one. A
        list coeffs::Array{Coefficient1p{ReducedKind},1} is returned.
"""
function nonScalarGeneral(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
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
        wR = chainRecoupling(leftCsf, rightCsf, ip, AngularJ64(k))
        if  wR == 0.0    continue    end

        value = -wR * wS * sqrt(Basics.twice(j) + 1.0) /
                     ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) )
        if  value != 0.0    push!( coeffs, Coefficient1p{ReducedKind}(k, sh, sh, value) )    end
    end

    return( coeffs )
end


"""
`SpinAngular.shellReducedA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64, Jket::AngularJ64, mq::AngularM64)`
    ... to compute the reduced matrix element of a single creation or annihilation operator within one subshell,
        <j^Nbra v J || a^(+/-) || j^Nket v' J'>, from the quasispin representation. As with `shellReducedW`, the
        coefficients of fractional parentage are reused as DATA -- `SpinAngular.completlyReducedCfpByIndices` -- and only
        the assembly is re-implemented. `mq` is +1/2 for creation and -1/2 for annihilation. A value::Float64 is returned.
"""
function shellReducedA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64,
                       Nket::Int64, senKet::Int64, Jket::AngularJ64, mq::AngularM64)
    key = (Basics.twice(j), Nbra, senBra, Basics.twice(Jbra), Nket, senKet, Basics.twice(Jket), Basics.twice(mq))
    haskey(SHELL_A_CACHE, key)  &&  return( SHELL_A_CACHE[key] )
    wa = shellReducedAUncached(j, Nbra, senBra, Jbra, Nket, senKet, Jket, mq)
    SHELL_A_CACHE[key] = wa

    return( wa )
end


"""
`SpinAngular.shellReducedAUncached(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64, Jket::AngularJ64, mq::AngularM64)`
    ... the body of `SpinAngular.shellReducedA`, evaluated without consulting the cache. Kept separate so that the
        cache is a wrapper and the physics is in one place. A value::Float64 is returned.
"""
function shellReducedAUncached(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64,
                               Nket::Int64, senKet::Int64, Jket::AngularJ64, mq::AngularM64)
    SA = SpinAngularTables
    # AN EMPTY SUBSHELL AGAINST A SINGLY OCCUPIED ONE NEEDS NO COEFFICIENT OF FRACTIONAL PARENTAGE: there is only
    # one state on each side, and the reduced matrix element is -/+ sqrt(2j+1) in closed form, + for the
    # annihilation j^1 -> j^0 and - for the creation j^0 -> j^1.
    #
    # Taking it here is not an optimisation. The quasispin tables consulted below -- qshellTermQ, qshellTermM,
    # getTermNumber, completlyReducedCfpByIndices -- cover j <= 9/2 ONLY, while a continuum partial wave reaches
    # j = 21/2. Until 27-Aug-2026 the lookup ran off the end of those tables and this routine returned 0, so every
    # two-particle coefficient moving an electron into or out of a high-j subshell was silently dropped: the
    # electron-impact excitation cross section of H-like C5+ came out 27 % low with nothing failing.
    if  (Nbra == 0  &&  Nket == 1)  ||  (Nbra == 1  &&  Nket == 0)
        nilBra = Nbra == 0
        if  nilBra  &&  (senBra != 0  ||  Basics.twice(Jbra) != 0  ||  senKet != 1  ||  Jket != j)   return( 0.0 )   end
        if !nilBra  &&  (senKet != 0  ||  Basics.twice(Jket) != 0  ||  senBra != 1  ||  Jbra != j)   return( 0.0 )   end
        return( (Basics.twice(mq) > 0  ?  -1.0  :  1.0) * sqrt(Basics.twice(j) + 1.0) )
    end
    # Beyond that closed form the tables are needed, and outside their reach the honest answer is to say so rather
    # than to return a zero that cannot be told from a selection rule.
    if  Basics.twice(j) > 9
        error("\n\nSpinAngular.shellReducedA: the quasispin/CFP tables reach j <= 9/2 and j = $(j) was asked "  *
              "for,\n>>> with occupations $Nket -> $Nbra, which is beyond the closed form for an empty against a "  *
              "singly\n>>> occupied subshell. Returning 0 here would be indistinguishable from a selection rule.\n")
    end
    Qb = SA.qshellTermQ(j, senBra);            Qk = SA.qshellTermQ(j, senKet)
    if  AngularMomentum.triangularDelta(Qb, AngularJ64(1//2), Qk) == 0        return( 0.0 )   end
    if  AngularMomentum.triangularDelta(Jbra, j, Jket) == 0                   return( 0.0 )   end
    bMQ = SA.qshellTermM(j, Nket);             aMQ = SA.qshellTermM(j, Nbra)
    ib  = SA.getTermNumber(j, Nbra, Qb, Jbra); ik = SA.getTermNumber(j, Nket, Qk, Jket)

    wa = - AngularMomentum.ClebschGordan(Qk, bMQ, AngularJ64(1//2), mq, Qb, aMQ)
    wa = wa * SA.completlyReducedCfpByIndices(ib, ik) / sqrt(Basics.twice(Qb) + 1.0)

    return( wa )
end


"""
`SpinAngular.substitutionRecoupling(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ja::AngularJ64, jb::AngularJ64, k::Int64)`
    ... to compute the recoupling factor for the two-subshell operator (A^(ja)(ia) x B^(jb)(ib))^(k) with ia < ib, i.e.
        for a CSF pair that differs by moving one electron between two subshells.

        The chain is cut at `ib`. Beyond it the operator of total rank k is peeled outwards exactly as for a
        single-subshell tensor. AT `ib` the two ranks join, which is where a NINE-j appears rather than a six-j:

            sqrt((2X_ib+1)(2X'_ib+1)(2k+1)) * { X_{ib-1}  X'_{ib-1}  ja ;  J_ib  J'_ib  jb ;  X_ib  X'_ib  k }

        and below it the operator of rank ja acting on subshell ia is reduced through the SAME `chainRecoupling` used for
        the equal-occupation case, restricted to the sub-chain 1 ... ib-1. Reusing it there rather than writing a second
        peeling loop is the reason this stays short. A value::Float64 is returned.
"""
function substitutionRecoupling(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ja::AngularJ64, jb::AngularJ64,
                                k::Int64)
    nw = length(leftCsf.occupation)
    wa = 1.0
    kJ = AngularJ64(k)

    # ... beyond ib: the same first-subsystem peeling as for one subshell
    for  q = ib+1:nw
        Xqm = leftCsf.subshellX[q-1];        Ypm = rightCsf.subshellX[q-1]
        Xq  = leftCsf.subshellX[q];          Yq  = rightCsf.subshellX[q]
        Jq  = leftCsf.subshellJ[q]
        if  Jq != rightCsf.subshellJ[q]      return( 0.0 )   end
        ph  = Int64( (Basics.twice(Xqm) + Basics.twice(Jq) + Basics.twice(Yq))//2 ) + k
        wa  = wa * (-1)^ph * sqrt((Basics.twice(Xq)+1.0)*(Basics.twice(Yq)+1.0)) *
                   AngularMomentum.Wigner_6j(Xqm, Xq, Jq, Yq, Ypm, kJ)
        if  wa == 0.0    return( 0.0 )   end
    end

    # ... at ib the two ranks couple, so a nine-j
    Xbm = (ib == 1) ? AngularJ64(0) : leftCsf.subshellX[ib-1]
    Ybm = (ib == 1) ? AngularJ64(0) : rightCsf.subshellX[ib-1]
    Xb  = leftCsf.subshellX[ib];         Yb  = rightCsf.subshellX[ib]
    Jb  = leftCsf.subshellJ[ib];         Jbp = rightCsf.subshellJ[ib]
    # A NINE-J WITH A ZERO ARGUMENT IS NOT A NINE-J. When the total rank is zero -- which is every call from the
    # two-particle direct term, where the two rank-k tensors couple to a scalar -- the symbol collapses exactly to a
    # six-j. That is worth doing rather than caching: a 9j costs 27.8 us here against 0.98 us for a 6j, so the
    # degenerate case is 28 times cheaper, and the identity was checked on 23328 argument combinations (1475 of them
    # non-zero) to a worst deviation of 5.6e-17 before being used.
    #
    #   {a b c; d e f; g h 0} = delta(c,f) delta(g,h) (-1)^(b+c+d+g) / sqrt((2c+1)(2g+1)) * {a b c; e d g}
    #
    if  Basics.twice(kJ) == 0  &&  Basics.twice(ja) == Basics.twice(jb)  &&  Basics.twice(Xb) == Basics.twice(Yb)
        ph9 = Basics.twice(Ybm) + Basics.twice(ja) + Basics.twice(Jb) + Basics.twice(Xb)
        if  iseven(ph9)
            nine = (-1)^Int64(ph9//2) * AngularMomentum.Wigner_6j(Xbm, Ybm, ja, Jbp, Jb, Xb) /
                   sqrt((Basics.twice(ja) + 1.0) * (Basics.twice(Xb) + 1.0))
        else
            nine = AngularMomentum.Wigner_9j(Xbm, Ybm, ja,  Jb, Jbp, jb,  Xb, Yb, kJ)
        end
    else
        nine = AngularMomentum.Wigner_9j(Xbm, Ybm, ja,  Jb, Jbp, jb,  Xb, Yb, kJ)
    end
    wa  = wa * sqrt((Basics.twice(Xb)+1.0)*(Basics.twice(Yb)+1.0)*(2.0*k+1.0)) * nine
    if  wa == 0.0                                                         return( 0.0 )   end

    # ... and below it, the rank-ja operator on subshell ia through the sub-chain 1 ... ib-1
    if  ib > 1    wa = wa * chainRecoupling(leftCsf, rightCsf, ia, ja; top = ib-1)    end

    return( wa )
end


"""
`SpinAngular.nonScalarSubstitution(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64)`
    ... to compute the rank-k coefficient of a CSF pair that differs by moving ONE electron: created in subshell `iCre`
        of the bra, annihilated in subshell `iAnn` of the ket. Three ingredients beyond the recoupling, each of which is
        easy to drop and each of which changes the answer:

        (1) the two single-subshell matrix elements, `shellReducedA` with mq = +1/2 and -1/2;
        (2) an ORDERING phase, since the recoupling is set up with the lower subshell index first, so a creation on the
            HIGHER index costs (-1)^(j_a + j_b - k + 1);
        (3) the JORDAN-WIGNER phase (-1)^(occupation between the two subshells, + 1) -- the sign from anticommuting the
            operator past the electrons that sit between them in the subshell ordering. It depends on the OTHER
            subshells' occupations, not on the two taking part, which is what makes it easy to forget.

        A list coeffs::Array{Coefficient1p{ReducedKind},1} is returned.
"""
function nonScalarSubstitution(op::SpinAngular.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                               subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64)
    coeffs = Coefficient1p{ReducedKind}[]
    shC    = subshells[iCre];                 shA = subshells[iAnn]
    if  !isAllowed1p(op, shC, shA)                                        return( coeffs )   end

    k    = op.rank
    jC   = AngularJ64( Basics.subshell_2j(shC)//2 );    jA = AngularJ64( Basics.subshell_2j(shA)//2 )
    ia   = min(iCre, iAnn);                   ib = max(iCre, iAnn)
    ja   = AngularJ64( Basics.subshell_2j(subshells[ia])//2 )
    jb   = AngularJ64( Basics.subshell_2j(subshells[ib])//2 )

    # every subshell other than the two taking part must be unchanged, or the two CSFs are orthogonal
    for  i = 1:length(subshells)
        if  i == ia  ||  i == ib    continue    end
        if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
            return( coeffs )
        end
    end

    wR = substitutionRecoupling(leftCsf, rightCsf, ia, ib, ja, jb, k)
    if  wR == 0.0                                                         return( coeffs )   end

    wC = shellReducedA(jC, leftCsf.occupation[iCre],  leftCsf.seniorityNr[iCre],  leftCsf.subshellJ[iCre],
                           rightCsf.occupation[iCre], rightCsf.seniorityNr[iCre], rightCsf.subshellJ[iCre],
                           AngularM64(1//2))
    if  wC == 0.0                                                         return( coeffs )   end
    wA = shellReducedA(jA, leftCsf.occupation[iAnn],  leftCsf.seniorityNr[iAnn],  leftCsf.subshellJ[iAnn],
                           rightCsf.occupation[iAnn], rightCsf.seniorityNr[iAnn], rightCsf.subshellJ[iAnn],
                           AngularM64(-1//2))
    if  wA == 0.0                                                         return( coeffs )   end

    wa = wR * wC * wA
    # ... (2) ordering, when the creation sits on the higher subshell index
    if  iCre == ib
        wa = wa * (-1)^Int64( (Basics.twice(ja) + Basics.twice(jb) - 2*k + 2)//2 )
    end
    # ... (3) Jordan-Wigner string over the subshells strictly between the two
    occup = 0
    for  i = ia:ib-1    occup = occup + leftCsf.occupation[i]    end
    wa = (-1)^(occup + 1) * wa

    # ... and the SAME outer normalization as the equal-occupation case, which is a result rather than an assumption:
    # with the phases right, the residual against GRASP came out 1.000, sqrt(3), sqrt(5) on sixteen coefficients, i.e.
    # exactly sqrt(2J_bra+1) at J_bra = 0, 1, 2.
    value = -wa * sqrt(Basics.twice(jC) + 1.0) /
                  ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) )
    if  value != 0.0    push!( coeffs, Coefficient1p{ReducedKind}(k, shC, shA, value) )    end

    return( coeffs )
end


"""
`SpinAngular.peelRange(leftCsf::CsfR, rightCsf::CsfR, from::Int64, to::Int64, kJ::AngularJ64)`
    ... to peel a tensor of total rank `kJ` outwards through the subshells `from` ... `to`, each of which carries NO
        operator and therefore enters as a spectator in the first subsystem,

            (-1)^(X_{q-1}+J_q+X'_q+k) sqrt((2X_q+1)(2X'_q+1)) { X_{q-1} X_q J_q ; X'_q X'_{q-1} k } .

        This is the outer loop of `SpinAngular.chainRecoupling` lifted out unchanged, so that the same peeling serves
        an operator acting on one subshell and one acting on several. An empty range returns 1. A value::Float64 is
        returned.
"""
function peelRange(leftCsf::CsfR, rightCsf::CsfR, from::Int64, to::Int64, kJ::AngularJ64)
    wa = 1.0
    k2 = Basics.twice(kJ)

    for  q = from:to
        Xqm = (q == 1) ? AngularJ64(0) : leftCsf.subshellX[q-1]
        Ypm = (q == 1) ? AngularJ64(0) : rightCsf.subshellX[q-1]
        Xq  = leftCsf.subshellX[q];              Yq = rightCsf.subshellX[q]
        Jq  = leftCsf.subshellJ[q]
        if  Jq != rightCsf.subshellJ[q]          return( 0.0 )   end
        ph  = Int64( (Basics.twice(Xqm) + Basics.twice(Jq) + Basics.twice(Yq) + k2)//2 )
        wa  = wa * (-1)^ph * sqrt((Basics.twice(Xq)+1.0)*(Basics.twice(Yq)+1.0)) *
                   AngularMomentum.Wigner_6j(Xqm, Xq, Jq, Yq, Ypm, kJ)
        if  wa == 0.0    return( 0.0 )   end
    end

    return( wa )
end


"""
`SpinAngular.actingFactor(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, kJ::AngularJ64)`
    ... to compute the factor contributed by the LOWEST subshell that carries an operator, where the tensor sits in the
        second subsystem with X_{ip-1} as spectator,

            (-1)^(X_{ip-1}+J'_ip+X_ip+k) sqrt((2X_ip+1)(2X'_ip+1)) { J_ip X_ip X_{ip-1} ; X'_ip J'_ip k } .

        Lifted out of `SpinAngular.chainRecoupling` unchanged, for the same reason as `SpinAngular.peelRange`. A
        value::Float64 is returned.
"""
function actingFactor(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, kJ::AngularJ64)
    # BELOW THE LOWEST ACTING SUBSHELL THE OPERATOR IS THE IDENTITY, so bra and ket must agree there -- in their
    # intermediate couplings X_q as well as in J_q. Nothing in the algebra below enforces it: the expression takes
    # X_{ip-1} from the BRA alone, and silently returns the diagonal answer for a pair whose couplings differ further
    # down. That produced a spurious 4s_1/2 coefficient for two CSFs of 1s^2 3d^2 4s that differ only in X_3, and the
    # give-away was that it was not symmetric under exchanging bra and ket (+0.745 against -0.447 at equal J).
    for  q = 1:ip-1
        if  leftCsf.subshellX[q] != rightCsf.subshellX[q]  ||  leftCsf.subshellJ[q] != rightCsf.subshellJ[q]
            return( 0.0 )
        end
    end
    Xpm = (ip == 1) ? AngularJ64(0) : leftCsf.subshellX[ip-1]
    Xp  = leftCsf.subshellX[ip];             Yp  = rightCsf.subshellX[ip]
    Jp  = leftCsf.subshellJ[ip];             Jpp = rightCsf.subshellJ[ip]
    ph  = Int64( (Basics.twice(Xpm) + Basics.twice(Jpp) + Basics.twice(Xp) + Basics.twice(kJ))//2 )

    return( (-1)^ph * sqrt((Basics.twice(Xp)+1.0)*(Basics.twice(Yp)+1.0)) *
                      AngularMomentum.Wigner_6j(Jp, Xp, Xpm, Yp, Jpp, kJ) )
end


"""
`SpinAngular.junctionFactor(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, below::AngularJ64, here::AngularJ64, above::AngularJ64)`
    ... to compute the factor contributed by a subshell at which a rank JOINS the chain: the accumulated rank `below`
        coming up from the subshells 1 ... ip-1 couples with the rank `here` of the operator acting on subshell ip to the
        accumulated rank `above`. This is where a NINE-j appears rather than a six-j,

            sqrt((2X_ip+1)(2X'_ip+1)(2above+1)) * { X_{ip-1} X'_{ip-1} below ; J_ip J'_ip here ; X_ip X'_ip above } ,

        with the degenerate collapse to a six-j taken whenever `above` is zero -- which is every junction of a two-body
        operator at its outermost acting subshell, and 28 times cheaper. A value::Float64 is returned.
"""
function junctionFactor(leftCsf::CsfR, rightCsf::CsfR, ip::Int64, below::AngularJ64, here::AngularJ64,
                        above::AngularJ64)
    Xbm = (ip == 1) ? AngularJ64(0) : leftCsf.subshellX[ip-1]
    Ybm = (ip == 1) ? AngularJ64(0) : rightCsf.subshellX[ip-1]
    Xb  = leftCsf.subshellX[ip];         Yb  = rightCsf.subshellX[ip]
    Jb  = leftCsf.subshellJ[ip];         Jbp = rightCsf.subshellJ[ip]

    if  Basics.twice(above) == 0  &&  Basics.twice(below) == Basics.twice(here)  &&
        Basics.twice(Xb)    == Basics.twice(Yb)
        ph9 = Basics.twice(Ybm) + Basics.twice(below) + Basics.twice(Jb) + Basics.twice(Xb)
        if  iseven(ph9)
            nine = (-1)^Int64(ph9//2) * AngularMomentum.Wigner_6j(Xbm, Ybm, below, Jbp, Jb, Xb) /
                   sqrt((Basics.twice(below) + 1.0) * (Basics.twice(Xb) + 1.0))
        else
            nine = AngularMomentum.Wigner_9j(Xbm, Ybm, below,  Jb, Jbp, here,  Xb, Yb, above)
        end
    else
        nine = AngularMomentum.Wigner_9j(Xbm, Ybm, below,  Jb, Jbp, here,  Xb, Yb, above)
    end

    return( sqrt((Basics.twice(Xb)+1.0)*(Basics.twice(Yb)+1.0)*(Basics.twice(above)+1.0)) * nine )
end


"""
`SpinAngular.treeRecoupling(leftCsf::CsfR, rightCsf::CsfR, sites::Array{Int64,1}, ranks::Array{AngularJ64,1}, inter::Array{AngularJ64,1})`
    ... to compute the recoupling factor for an operator that acts on ANY number of subshells, whose shell tensors are
        coupled ALONG THE SUBSHELL CHAIN itself: `sites` holds the acting subshell indices in increasing order, `ranks[m]`
        the rank of the tensor on `sites[m]`, and `inter[m]` the accumulated rank after that site, so that `inter[1]`
        equals `ranks[1]` and `inter[end]` is the total rank of the operator.

        The whole chain is then three kinds of factor and nothing else: `SpinAngular.peelRange` between acting
        subshells, `SpinAngular.junctionFactor` at each acting subshell above the first, and
        `SpinAngular.actingFactor` at the lowest one. That is what makes this ONE routine instead of a case tree:
        `chainRecoupling(ip, k)` is the single site `[ip]`, and `substitutionRecoupling(ia, ib, ja, jb, k)` is the two
        sites `[ia, ib]` -- both are checked against it rather than assumed. A value::Float64 is returned.
"""
function treeRecoupling(leftCsf::CsfR, rightCsf::CsfR, sites::Array{Int64,1}, ranks::Array{AngularJ64,1},
                        inter::Array{AngularJ64,1})
    nw = length(leftCsf.occupation)
    ns = length(sites)

    wa = peelRange(leftCsf, rightCsf, sites[ns]+1, nw, inter[ns])
    if  wa == 0.0                                                         return( 0.0 )   end

    for  m = ns:-1:2
        wa = wa * junctionFactor(leftCsf, rightCsf, sites[m], inter[m-1], ranks[m], inter[m])
        if  wa == 0.0                                                     return( 0.0 )   end
        wa = wa * peelRange(leftCsf, rightCsf, sites[m-1]+1, sites[m]-1, inter[m-1])
        if  wa == 0.0                                                     return( 0.0 )   end
    end

    return( wa * actingFactor(leftCsf, rightCsf, sites[1], ranks[1]) )
end
