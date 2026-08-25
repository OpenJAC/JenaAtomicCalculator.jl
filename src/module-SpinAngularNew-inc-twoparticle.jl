
#
# Two-particle (electron-electron) coefficients.
#
# This file is included from module-SpinAngularNew.jl and shares its namespace.
#

"""
`abstract type SpinAngularNew.AbstractTwoElectronKind`
    ... defines an abstract type to distinguish WHICH two-electron radial quantity a two-particle coefficient is to be
        multiplied by. This is not a hypothetical hazard: JAC and GRASP2018 decompose the electron-electron interaction
        onto DIFFERENT quantities, so their coefficients are not comparable term by term, and a coefficient carried from
        one into the other is wrong by the reduced C^L factors rather than obviously broken.

    + EffectiveStrengthKind   ... multiplies the effective strength X^L(abcd). This is JAC's convention and the GENERAL
                                  one: the SAME coefficient serves the Coulomb, Breit and Gaunt operators, because the
                                  tensorial structure that distinguishes them lives in X^L rather than in the
                                  coefficient. `module-Hamiltonian.jl:282-292` uses one `coeff.V` with both
                                  `XL_Coulomb` and `XL_Breit`, which is the whole argument for the design.
    + GraspCoulombKind        ... multiplies the plain Slater integral R^k(abcd). This is GRASP's COULOMB convention and
                                  ONLY that -- it does not generalise. GRASP's Breit path uses a different callback
                                  (BREID), carries a sixth label ITYPE in 1..6, and multiplies one of SIX integral
                                  routines BRINT1..BRINT6 chosen by that tag, with the sign of the tag routing the
                                  contribution. There is no single "plain radial integral" for Breit, so this kind
                                  exists ONLY to express a comparison with the GRASP Coulomb oracle, and is not a way
                                  anyone should compute anything.
"""
abstract type  AbstractTwoElectronKind                                        end
struct         EffectiveStrengthKind  <:  AbstractTwoElectronKind             end
struct         GraspCoulombKind      <:  AbstractTwoElectronKind             end


"""
`struct  SpinAngularNew.Coefficient2p{K<:AbstractTwoElectronKind}`
    ... a struct for a single spin-angular coefficient of a two-particle matrix element. As for the one-particle case the
        type parameter records which radial quantity the coefficient belongs to, so that mixing the two conventions has
        no method rather than a plausible value.

    + nu       ::Int64      ... Rank (k) of the two-particle interaction strength.
    + a        ::Subshell   ... First left-hand subshell.
    + b        ::Subshell   ... Second left-hand subshell.
    + c        ::Subshell   ... First right-hand subshell.
    + d        ::Subshell   ... Second right-hand subshell.
    + V        ::Float64    ... (Value of) the spin-angular coefficient.
"""
struct  Coefficient2p{K<:AbstractTwoElectronKind}
    nu         ::Int64
    a          ::Subshell
    b          ::Subshell
    c          ::Subshell
    d          ::Subshell
    V          ::Float64
end


# `Base.show(io::IO, coeff::SpinAngularNew.Coefficient2p)`  ... prepares a proper printout of coeff.
function Base.show(io::IO, coeff::SpinAngularNew.Coefficient2p{K})  where K<:AbstractTwoElectronKind
    print(io, "   V^$(coeff.nu) [$(coeff.a), $(coeff.b); $(coeff.c), $(coeff.d)] = $(coeff.V)   ($(K.name.name))")
end


"""
`SpinAngularNew.toGraspCoulomb(coeff::SpinAngularNew.Coefficient2p{EffectiveStrengthKind})`
    ... to convert a coefficient of the effective strength X^L into the corresponding coefficient of GRASP's COULOMB
        convention, i.e. of the plain Slater integral R^k, by restoring the reduced C^k factors that X^L carries:

            V_grasp  =  V_effective * (-1)^k * <kappa_a||C^k||kappa_c> * <kappa_b||C^k||kappa_d>

        THE (-1)^k IS NOT A FITTED PHASE. It is read from `InteractionStrength.XL_CoulombReference`, which forms
        `xc = CL_reduced_me(a,L,c) * CL_reduced_me(b,L,d)` and then `if rem(L,2) == 1   xc = -xc   end` before
        multiplying by the Slater integral. Omitting it left exactly the 14 odd-k EXCHANGE terms of 1s^2 2s^2 2p^2
        disagreeing with GRASP by a sign, and nothing else; restoring it made all 61 agree.

        THIS IS NOT A RELATION EITHER CODE STATES. GRASP has no effective strength and never writes this equation; it is
        a numerical bridge, built here so that a comparison against the GRASP oracle can be expressed at all, and
        verified on ALL 61 surviving coefficients of 1s^2 2s^2 2p^2, direct and exchange, every rank.

        It applies to the COULOMB interaction ONLY. GRASP's Breit coefficients multiply one of six different integral
        kinds selected by a type tag, so no analogous single conversion exists there -- which is precisely why the
        effective strength is the better internal representation.

        A coeff::Coefficient2p{GraspCoulombKind} is returned.
"""
function toGraspCoulomb(coeff::SpinAngularNew.Coefficient2p{EffectiveStrengthKind})
    wa = coeff.V * AngularMomentum.CL_reduced_me(coeff.a, coeff.nu, coeff.c) *
                   AngularMomentum.CL_reduced_me(coeff.b, coeff.nu, coeff.d)
    if  isodd(coeff.nu)    wa = -wa    end

    return( Coefficient2p{GraspCoulombKind}(coeff.nu, coeff.a, coeff.b, coeff.c, coeff.d, wa) )
end


"""
`SpinAngularNew.computeCoefficients(op::SpinAngularNew.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... to compute the spin-angular coefficients of the two-particle (electron-electron) interaction for a DIAGONAL CSF
        pair, in JAC's convention, i.e. as coefficients of the effective strength X^L. A list
        coeffs::Array{Coefficient2p{EffectiveStrengthKind},1} is returned.

        Three terms make up the answer and each was settled separately, each verified against `SpinAngular` on
        configurations it was not fitted to:

        * the DIRECT term between two distinct subshells, `twoParticleDirect`, general at any occupations;
        * the EXCHANGE term between two distinct subshells, `twoParticleExchange`, a Racah sum over the direct channel;
        * the SAME-SUBSHELL term, `twoParticleSameShell`, the two-body quasispin object less its normal-ordering
          correction.

        OFF-DIAGONAL CSF PAIRS ARE NOT COVERED and raise. A pair differing in coupling or occupation needs recoupling
        this method does not do, and returning the diagonal answer for it would be a wrong number rather than a
        missing one.
"""
coeffs2pEmpty() = Coefficient2p{EffectiveStrengthKind}[]


function computeCoefficients(op::SpinAngularNew.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                             subshells::Array{Subshell,1})
    if  op.rank != 0
        error("\n\nSpinAngularNew.computeCoefficients: only the scalar (rank-0) two-particle operator is defined.\n")
    end
    # EQUAL OCCUPATIONS are required; the couplings may differ. Each term is bra/ket aware -- the direct and
    # exchange terms through `substitutionRecoupling`, the same-subshell term through its own orthogonality guard --
    # so a pair differing only in coupling is handled rather than refused, and gives zero where it should.
    if  leftCsf.J != rightCsf.J  ||  leftCsf.parity != rightCsf.parity     return( coeffs2pEmpty() )   end
    # ... a pair that MOVES ONE electron between two subshells has its own assembly; anything further apart is
    #     orthogonal for a two-body operator, and returning the equal-occupation answer for it would be a wrong
    #     number rather than a missing one.
    if  leftCsf.occupation != rightCsf.occupation
        diff = leftCsf.occupation - rightCsf.occupation
        if  count(!=(0), diff) == 2  &&  sum(abs, diff) == 2
            mv   = findall(!=(0), diff)
            iCre = diff[mv[1]] > 0 ? mv[1] : mv[2];      iAnn = diff[mv[1]] > 0 ? mv[2] : mv[1]
            return( twoParticleMoveOne(leftCsf, rightCsf, subshells, iCre, iAnn) )
        elseif  count(!=(0), diff) == 4  &&  sum(abs, diff) == 4
            cre = sort([i for i in findall(!=(0), diff) if diff[i] > 0])
            ann = sort([i for i in findall(!=(0), diff) if diff[i] < 0])
            return( twoParticleMoveTwo(leftCsf, rightCsf, subshells, cre, ann) )
        elseif  count(!=(0), diff) == 2  &&  sum(abs, diff) == 4
            mv   = findall(!=(0), diff)
            iCre = diff[mv[1]] > 0 ? mv[1] : mv[2];      iAnn = diff[mv[1]] > 0 ? mv[2] : mv[1]
            return( twoParticleMoveTwoBothDoubled(leftCsf, rightCsf, subshells, iCre, iAnn) )
        elseif  count(!=(0), diff) == 3  &&  sum(abs, diff) == 4
            dbl = findfirst(i -> abs(diff[i]) == 2, 1:length(diff))
            sgl = sort([i for i in findall(!=(0), diff) if abs(diff[i]) == 1])
            if  !isnothing(dbl)  &&  length(sgl) == 2
                return( twoParticleMoveTwoSameShell(leftCsf, rightCsf, subshells, dbl, sgl, diff) )
            end
        elseif  count(!=(0), diff) > 4  ||  sum(abs, diff) > 4
            return( coeffs2pEmpty() )
        end
        # ... every occupation pattern a two-body operator can connect is handled above; anything reaching here is a
        #     pattern that should not exist, so it raises rather than returning an empty list that would look like a
        #     legitimate zero.
        error("\n\nSpinAngularNew.computeCoefficients (two-particle): unreachable occupation pattern, "  *
              "count = $(count(!=(0), diff)), sum = $(sum(abs, diff)).\n")
    end

    coeffs = Coefficient2p{EffectiveStrengthKind}[]
    nw     = length(subshells)
    for  ia = 1:nw
        leftCsf.occupation[ia] == 0  &&  continue
        sha = subshells[ia]
        # ... the same-subshell term, which needs two electrons in the one shell
        if  leftCsf.occupation[ia] >= 2
            for  k = 0:Basics.subshell_2j(sha)
                v = twoParticleSameShell(leftCsf, rightCsf, subshells, ia, k)
                abs(v) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, sha, sha, sha, sha, v))
            end
        end
        # ... and the direct and exchange terms with every higher subshell. The direct vector is formed ONCE per
        #     subshell pair and both families are read off it: every exchange coefficient is a Racah sum over the
        #     whole vector, so computing them rank by rank would re-evaluate it kMax + 2 times over.
        for  ib = ia+1:nw
            leftCsf.occupation[ib] == 0  &&  continue
            shb  = subshells[ib]
            ja   = AngularJ64( Basics.subshell_2j(sha)//2 );   jb = AngularJ64( Basics.subshell_2j(shb)//2 )
            kMax = Int64( (Basics.subshell_2j(sha) + Basics.subshell_2j(shb))//2 )
            dvec = twoParticleDirectVector(leftCsf, rightCsf, subshells, ia, ib, kMax)
            for  k = 0:kMax
                vd = dvec[k+1]
                abs(vd) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, sha, shb, sha, shb, vd))
                ve = exchangeFromDirect(dvec, ja, jb, k)
                abs(ve) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, sha, shb, shb, sha, ve))
            end
        end
    end

    return( coeffs )
end

"""
`SpinAngularNew.twoParticleDirect(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64, k::Int64)`
    ... to compute the DIRECT two-particle coefficient V^k(a,b;a,b) for two distinct open subshells `ia` and `ib`, in
        the effective-strength convention. The direct part of the electron-electron interaction is the scalar product
        of two rank-k tensors, one on each subshell, so it is built from the same pieces as the one-particle case:

            V^k(a,b;a,b)  =  <j_a^N||W^(k)||j_a^N> <j_b^N||W^(k)||j_b^N> R_chain / ( sqrt(2J_bra+1) sqrt(2k+1) )

        with R_chain = `substitutionRecoupling(leftCsf, rightCsf, ia, ib, k, k, 0)` -- the SAME routine the
        single-electron substitution uses, called with both ranks equal to k and coupled to zero, which is what makes
        this short.

        THE NORMALIZATION WAS CALIBRATED AND THEN TESTED OUT OF SAMPLE. Four values of 1s 2s fixed the two factors;
        twenty further coefficients of 1s 2p, 2s 3d and 1s 3d, spanning J = 0, 1, 2, 3 and j = 1/2, 3/2, 5/2, then
        reproduced SpinAngular exactly, every ratio 1.0000000.

        ITS REACH WAS MEASURED, NOT ASSUMED. Applied to EVERY pair of distinct subshells, classified by occupation,
        it is exact in all four classes: closed/closed (19), closed/open (9), closed/single (72), single/single (68) --
        168 coefficients, every ratio 1.0. So it is general for two distinct subshells at any occupations, which is
        more than was expected of it; what it does NOT cover is the same-subshell term (a,a,a,a).

        THIS MATTERS BECAUSE THE FIRST ATTEMPT AT THIS TERM WAS WRONG. A closed form fitted to closed-shell data gave
        the k = 0 direct term correctly and everything else wrongly, because a closed shell forces J = 0 and hides the
        J-dependence entirely. The expression above carries J through the recoupling and reproduces, for instance,
        +0.5 at J = 0 and -0.1666667 at J = 1 for the same k = 1 term of 1s 2s. A value::Float64 is returned.
"""
function twoParticleDirect(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64,
                           k::Int64)
    sha = subshells[ia];    shb = subshells[ib]
    ja  = AngularJ64( Basics.subshell_2j(sha)//2 );    jb = AngularJ64( Basics.subshell_2j(shb)//2 )

    wa = shellReducedW(ja, leftCsf.occupation[ia], leftCsf.seniorityNr[ia],  leftCsf.subshellJ[ia],
                                                   rightCsf.seniorityNr[ia], rightCsf.subshellJ[ia], k)
    if  wa == 0.0                                                         return( 0.0 )   end
    wb = shellReducedW(jb, leftCsf.occupation[ib], leftCsf.seniorityNr[ib],  leftCsf.subshellJ[ib],
                                                   rightCsf.seniorityNr[ib], rightCsf.subshellJ[ib], k)
    if  wb == 0.0                                                         return( 0.0 )   end
    wR = substitutionRecoupling(leftCsf, rightCsf, ia, ib, AngularJ64(k), AngularJ64(k), 0)

    return( wa * wb * wR / ( sqrt(Basics.twice(leftCsf.J) + 1.0) * sqrt(2.0*k + 1.0) ) )
end

"""
`SpinAngularNew.twoParticleExchange(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64, k::Int64)`
    ... to compute the EXCHANGE two-particle coefficient V^k(a,b;b,a) for two distinct subshells at ANY occupations,
        in the effective-strength convention. The exchange channel is not the direct one relabelled: it expands over
        the direct channel at intermediate ranks, by the Racah transformation between the two coupling schemes,

            V^k(a,b;b,a)  =  SUM_K  (2K+1) { j_a  j_b  k ;  j_b  j_a  K }  V^K(a,b;a,b)

        so it needs no machinery of its own -- `SpinAngularNew.twoParticleDirect` supplies every term of the sum.

        HOW THE FORM WAS ARRIVED AT, since a fitted 6j would be worth little. For singly occupied subshells the sum
        collapses to one symbol, and that collapsed case was settled first BY ELIMINATION: of five candidate 6j
        orderings, three VANISH at ranks where the coefficient does not, which refutes them outright, and of the two
        survivors only the phase (-1)^(j_a+j_b+k) holds across 68 coefficients. The general sum above was then written
        and required to REPRODUCE that collapsed case before replacing it, which it does.

        A value::Float64 is returned.
"""
function twoParticleExchange(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64,
                             k::Int64)
    ja = AngularJ64( Basics.subshell_2j(subshells[ia])//2 )
    jb = AngularJ64( Basics.subshell_2j(subshells[ib])//2 )
    kJ = AngularJ64(k)
    wa = 0.0
    # ... the intermediate rank cannot exceed either subshell's own triangle
    kMax = Int64( (Basics.twice(ja) + Basics.twice(jb))//2 )
    for  K = 0:kMax
        d = twoParticleDirect(leftCsf, rightCsf, subshells, ia, ib, K)
        if  d == 0.0    continue    end
        wa = wa + (2.0*K + 1.0) * AngularMomentum.Wigner_6j(ja, jb, kJ, jb, ja, AngularJ64(K)) * d
    end

    return( wa )
end

"""
`SpinAngularNew.shellReducedWByIndex(j::AngularJ64, Nbra::Int64, ibra::Int64, Nket::Int64, iket::Int64, kj::Int64)`
    ... as `SpinAngularNew.shellReducedW`, but addressing the two shell terms by their quasispin INDEX rather than by
        seniority and J. Needed because the two-body assembly sums over INTERMEDIATE shell terms, which are naturally
        enumerated by index. A value::Float64 is returned.
"""
function shellReducedWByIndex(j::AngularJ64, Nbra::Int64, ibra::Int64, Nket::Int64, iket::Int64, kj::Int64)
    SA = JenaAtomicCalculator.SpinAngular
    tb = SA.qspaceTerms(ibra);    tk = SA.qspaceTerms(iket)
    if  kj == 0
        return( ibra == iket ? -Nbra * sqrt((Basics.twice(tb.J)+1.0)/(Basics.twice(j)+1.0)) : 0.0 )
    end
    kq = iseven(kj) ? 1 : 0
    if  AngularMomentum.triangularDelta(tb.Q, AngularJ64(kq), tk.Q) == 0    return( 0.0 )   end
    wa = AngularMomentum.ClebschGordan(tk.Q, SA.qshellTermM(j, Nket), AngularJ64(kq), AngularM64(0),
                                       tb.Q, SA.qshellTermM(j, Nbra))
    wa = wa * SA.completelyReducedWkk(ibra, iket, kq, kj) / sqrt((Basics.twice(tb.Q) + 1.0) * 2.0)

    return( wa )
end


"""
`SpinAngularNew.shellWW(j::AngularJ64, N::Int64, ibra::Int64, iket::Int64, kj::Int64)`
    ... to compute the TWO-body quasispin object <j^N || (W^(k) W^(k))^(0) || j^N> within one subshell, as a closure
        sum over the shell's own intermediate terms:

            sum_r  (-1)^((2k - 2J_bra + 2J_r)/2)  <bra||W^(k)||r> <r||W^(k)||ket>   / sqrt((2k+1)(2J_bra+1))

        restricted to the even-occupation terms r, since W conserves particle number, and filtered by the six-j
        {k k 0; J_ket J_bra J_r} which decides whether the intermediate can contribute at all.

        THIS IS NOT THE SAME-SUBSHELL COEFFICIENT, and mistaking it for one cost two failed attempts. Compared
        directly against the coefficient it disagrees everywhere -- for j = 3/2, N = 2, J = 0 the coefficient is 0.25
        at every rank while this object gives 1.0, 0, 2.236, 0. The object was right; what was missing was the factor
        of one half and the normal-ordering subtraction that `SpinAngularNew.twoParticleSameShell` supplies.

        A value::Float64 is returned.
"""
function shellWW(j::AngularJ64, N::Int64, ibra::Int64, iket::Int64, kj::Int64)
    SA = JenaAtomicCalculator.SpinAngular
    tb = SA.qspaceTerms(ibra);    tk = SA.qspaceTerms(iket);    wa = 0.0
    for  r = tb.min_even:tb.max_even
        tr = SA.qspaceTerms(r)
        if  SA.qspacedelta(tr.Q, SA.qshellTermM(tr.j, N)) == 0    continue    end
        if  AngularMomentum.Wigner_6j(AngularJ64(kj), AngularJ64(kj), AngularJ64(0), tk.J, tb.J, tr.J) == 0.0
            continue
        end
        ph = (-1)^Int64( (2*kj - Basics.twice(tb.J) + Basics.twice(tr.J))//2 )
        wa = wa + ph * shellReducedWByIndex(j, N, ibra, N, r, kj) * shellReducedWByIndex(j, N, r, N, iket, kj)
    end

    return( wa / sqrt((2.0*kj + 1.0) * (Basics.twice(tb.J) + 1.0)) )
end


"""
`SpinAngularNew.twoParticleSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, k::Int64)`
    ... to compute the SAME-SUBSHELL two-particle coefficient V^k(a,a;a,a), in the effective-strength convention:

            V^k(a,a;a,a)  =  0.5 * ( WW(k)/sqrt(2k+1)  -  (-1)^(2j+k) W(0)/sqrt(2j+1) ) / sqrt(2J+1)

        The second term is a NORMAL-ORDERING correction and is the whole reason two earlier attempts failed. Both
        operators act on one shell, so the two-body object counts each electron with itself; that self-interaction has
        to be removed, and it is removed by the rank-0 one-body element, not by anything rank-dependent. Two attempts
        compared the two-body object directly against the coefficient and concluded the closure route was dead. It was
        not: the object was correct and the wrapper was missing.

        Verified against `SpinAngular` on 254 coefficients over four configurations -- 2p^2, 2p^4, 3d^2 and 3d^3,
        spanning j = 1/2, 3/2, 5/2, seniorities 0, 2 and 3, and every allowed rank -- worst ratio 1.000000000000.

        A value::Float64 is returned.
"""
function twoParticleSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, k::Int64)
    SA = JenaAtomicCalculator.SpinAngular
    sh = subshells[ia];    N = leftCsf.occupation[ia]
    if  N < 2                                                             return( 0.0 )   end
    # (W^(k) x W^(k))^(0) is a SCALAR on the subshell it acts in, so it cannot change any coupling: every other
    # subshell must be in an identical state, and the running couplings X must agree throughout, or the two CSFs are
    # orthogonal and the element vanishes. Without this the routine returns the DIAGONAL value on a pair that differs
    # only in another subshell's coupling -- which is exactly the defect this module documents in GRASP2018 at
    # example-Aq.jl branch k, and it was present here until measured against SpinAngular on off-diagonal pairs.
    if  leftCsf.subshellX != rightCsf.subshellX  ||  leftCsf.subshellJ != rightCsf.subshellJ
        return( 0.0 )
    end
    for  i = 1:length(subshells)
        if  i != ia  &&  leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]     return( 0.0 )   end
    end
    j  = AngularJ64( Basics.subshell_2j(sh)//2 )
    Jb = leftCsf.subshellJ[ia];    Jk = rightCsf.subshellJ[ia]
    ib = SA.getTermNumber(j, N, SA.qshellTermQ(j, leftCsf.seniorityNr[ia]),  Jb)
    ik = SA.getTermNumber(j, N, SA.qshellTermQ(j, rightCsf.seniorityNr[ia]), Jk)
    w0 = shellReducedWByIndex(j, N, ib, N, ik, 0) / sqrt(Basics.twice(j) + 1.0)
    wa = 0.5 * ( shellWW(j, N, ib, ik, k)/sqrt(2.0*k + 1.0) -
                 (-1)^Int64(Basics.twice(j) + k) * w0 ) / sqrt(Basics.twice(Jb) + 1.0)

    return( wa )
end

"""
`SpinAngularNew.twoParticleDirectVector(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64, kMax::Int64)`
    ... to compute the direct coefficients V^K(a,b;a,b) for every rank K = 0 ... kMax at once. A
        vector::Array{Float64,1} indexed as K+1 is returned.

        WHY THIS EXISTS. Every exchange coefficient is a Racah sum over the WHOLE direct vector, so computing the
        exchange rank by rank re-evaluates the same direct terms once per rank: for one subshell pair with kMax = 3
        that is 20 evaluations of which 4 are distinct, and for kMax = 4 it is 30 of which 5 are distinct. The
        redundancy is a factor kMax + 2 and grows with j. Forming the vector once removes it.
"""
function twoParticleDirectVector(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, ia::Int64, ib::Int64,
                                 kMax::Int64)
    return( [twoParticleDirect(leftCsf, rightCsf, subshells, ia, ib, K)  for K = 0:kMax] )
end


"""
`SpinAngularNew.exchangeFromDirect(direct::Array{Float64,1}, ja::AngularJ64, jb::AngularJ64, k::Int64)`
    ... to form the exchange coefficient at rank k from a direct vector already computed, by the Racah transformation

            V^k(a,b;b,a)  =  SUM_K  (2K+1) { j_a  j_b  k ;  j_b  j_a  K }  V^K(a,b;a,b)

        This is the same expression `SpinAngularNew.twoParticleExchange` evaluates; the difference is only that the
        direct vector is supplied rather than recomputed. A value::Float64 is returned.
"""
function exchangeFromDirect(direct::Array{Float64,1}, ja::AngularJ64, jb::AngularJ64, k::Int64)
    wa = 0.0;    kJ = AngularJ64(k)
    for  K = 0:length(direct)-1
        direct[K+1] == 0.0  &&  continue
        wa = wa + (2.0*K + 1.0) * AngularMomentum.Wigner_6j(ja, jb, kJ, jb, ja, AngularJ64(K)) * direct[K+1]
    end

    return( wa )
end


"""
`SpinAngularNew.moveOnePrimaryVector(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64, kMax::Int64, wCre::Float64, wAnn::Float64)`
    ... to compute the primary one-electron-move coefficients for ALL ranks 0 ... kMax of one spectator subshell at once,
        returning them as a vector. Only the spectator's own rank changes with k: the ordering of the three acting
        subshells along the chain, and their places in it, do not -- so they are formed once and the rank vectors are
        REUSED rather than rebuilt for every rank. This is the same reason `SpinAngularNew.twoParticleDirectVector`
        exists for the equal-occupation case, and it is worth doing for the same reason: the work is unchanged and the
        allocation is not. A vector prim::Array{Float64,1} is returned.
"""
function moveOnePrimaryVector(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64,
                              iSpec::Int64, kMax::Int64, wCre::Float64, wAnn::Float64)
    prim = zeros(Float64, kMax+1)
    if  wCre == 0.0  ||  wAnn == 0.0                                      return( prim )   end
    jA  = AngularJ64( Basics.subshell_2j(subshells[iCre])//2 )
    jD  = AngularJ64( Basics.subshell_2j(subshells[iAnn])//2 )
    jS  = AngularJ64( Basics.subshell_2j(subshells[iSpec])//2 )

    ord   = sortperm( [iCre, iAnn, iSpec] )
    sites = [iCre, iAnn, iSpec][ord]
    base  = AngularJ64[ jA, jD, AngularJ64(0) ]
    ranks = Vector{AngularJ64}(undef, 3);      inter = Vector{AngularJ64}(undef, 3)
    between = min(iCre,iAnn) < iSpec < max(iCre,iAnn)
    occup   = 0
    for  i = min(iCre,iAnn):max(iCre,iAnn)-1   occup = occup + leftCsf.occupation[i]   end

    for  k = 0:kMax
        base[3] = AngularJ64(k)
        for  m = 1:3    ranks[m] = base[ord[m]]    end
        inter[1] = ranks[1];    inter[2] = ranks[3];    inter[3] = AngularJ64(0)
        if  AngularMomentum.triangularDelta(inter[1], ranks[2], inter[2]) == 0    continue    end
        wS = shellReducedW(jS, leftCsf.occupation[iSpec], leftCsf.seniorityNr[iSpec], leftCsf.subshellJ[iSpec],
                               rightCsf.seniorityNr[iSpec], rightCsf.subshellJ[iSpec], k)
        if  wS == 0.0                                                             continue    end
        wR = treeRecoupling(leftCsf, rightCsf, sites, ranks, inter)
        if  wR == 0.0                                                             continue    end

        ph = (-1)^(occup + 1)
        if  iAnn < iCre
            ph = ph * (-1)^Int64( (Basics.twice(jA) + Basics.twice(jD) - 2*k + 2)//2 )
        end
        if  between
            ph = ph * (-1)^Int64( (Basics.twice(jA) + Basics.twice(jD) + 2*k + 2)//2 )
        end
        prim[k+1] = ph * wR * wCre * wAnn * wS / ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) )
    end

    return( prim )
end


"""
`SpinAngularNew.moveOneQuads(subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64)`
    ... to determine the two index quadruples (a,b,c,d) under which the one-electron-move terms of a given spectator
        subshell are emitted, together with their four subshell angular momenta and whether the spectator lies BETWEEN
        the two subshells that change occupation.

        The rank of a two-particle coefficient couples (a,c) and (b,d), so the quadruple says which pairing a value
        belongs to, and there are two of them for each spectator: the DIRECT pairing, spectator with spectator and
        acceptor with donor, and the CROSSED one, acceptor with spectator and spectator with donor. Which of the two the
        assembly of `SpinAngularNew.moveOnePrimary` yields is decided by the arrangement, and that is the whole reason
        this small routine exists rather than being inlined: with the spectator outside it yields the direct pairing,
        with the spectator between it yields the crossed one. A tuple (qPrimary, qPartner, js, between) is returned.
"""
function moveOneQuads(subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64)
    shA = subshells[iCre];    shD = subshells[iAnn];    shS = subshells[iSpec]
    between = min(iCre,iAnn) < iSpec < max(iCre,iAnn)
    aIsS    = iSpec < iCre
    (a, b)  = aIsS ? (shS, shA) : (shA, shS)
    firstSD = between ? !aIsS : aIsS
    (c, d)  = firstSD ? (shS, shD) : (shD, shS)
    js      = ( AngularJ64(Basics.subshell_2j(a)//2), AngularJ64(Basics.subshell_2j(b)//2),
                AngularJ64(Basics.subshell_2j(c)//2), AngularJ64(Basics.subshell_2j(d)//2) )

    return( ((a,b,c,d), (a,b,d,c), js, between) )
end


"""
`SpinAngularNew.moveOnePrimary(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64, k::Int64)`
    ... to compute the primary two-particle coefficient, at rank `k`, of a CSF pair that differs by ONE electron -- created
        in subshell `iCre` of the bra, annihilated in `iAnn` of the ket -- with the two-body operator's second index pair
        acting on the distinct spectator subshell `iSpec`. Three subshells act at once, and the operator is assembled from
        pieces that already exist rather than from a new case tree:

            a creation `shellReducedA(+1/2)` on iCre, an annihilation `shellReducedA(-1/2)` on iAnn, a spectator tensor
            `shellReducedW` of rank k on iSpec, and ONE call to `treeRecoupling` over the three sites,

        divided by sqrt(2k+1) sqrt(2J_bra+1) as in the equal-occupation case. Three phases decide the sign, and each was
        established on data rather than assumed:

        (1) the JORDAN-WIGNER string (-1)^(occupation between the two changing subshells, + 1), taken unconditionally --
            the spectator's own electrons included when it lies between them;
        (2) an ORDERING phase (-1)^(j_a + j_d - k + 1) when the creation sits on the HIGHER subshell index, the same
            phase the one-particle substitution carries;
        (3) a RE-PAIRING phase (-1)^(j_a + j_d + k + 1) when the spectator lies BETWEEN the two, because the coupling
            tree then joins the acceptor with the spectator instead of with the donor. It is a function of the three
            ranks alone, as a re-pairing of three tensors coupled to zero must be; it was found as a parity rule exact
            on 518 coefficients and then held on 8820 further ones from six configuration sets not used to find it.

        A value::Float64 is returned.
"""
function moveOnePrimary(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64,
                        iSpec::Int64, k::Int64)
    jA = AngularJ64( Basics.subshell_2j(subshells[iCre])//2 )
    jD = AngularJ64( Basics.subshell_2j(subshells[iAnn])//2 )
    wCre = shellReducedA(jA, leftCsf.occupation[iCre],  leftCsf.seniorityNr[iCre],  leftCsf.subshellJ[iCre],
                             rightCsf.occupation[iCre], rightCsf.seniorityNr[iCre], rightCsf.subshellJ[iCre],
                             AngularM64(1//2))
    wAnn = shellReducedA(jD, leftCsf.occupation[iAnn],  leftCsf.seniorityNr[iAnn],  leftCsf.subshellJ[iAnn],
                             rightCsf.occupation[iAnn], rightCsf.seniorityNr[iAnn], rightCsf.subshellJ[iAnn],
                             AngularM64(-1//2))

    return( moveOnePrimaryFrom(leftCsf, rightCsf, subshells, iCre, iAnn, iSpec, k, wCre, wAnn) )
end


"""
`SpinAngularNew.moveOnePrimaryFrom(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64, k::Int64, wCre::Float64, wAnn::Float64)`
    ... as `SpinAngularNew.moveOnePrimary`, but taking the creation and annihilation matrix elements ALREADY FORMED.
        They depend on neither the spectator subshell nor the rank, so a caller sweeping over both computes them once
        and passes them in; `moveOnePrimary` itself is the convenience entry that forms them and delegates here. A
        value::Float64 is returned.
"""
function moveOnePrimaryFrom(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64,
                            iSpec::Int64, k::Int64, wCre::Float64, wAnn::Float64)
    jA = AngularJ64( Basics.subshell_2j(subshells[iCre])//2 )
    jD = AngularJ64( Basics.subshell_2j(subshells[iAnn])//2 )
    jS = AngularJ64( Basics.subshell_2j(subshells[iSpec])//2 )

    ord   = sortperm( [iCre, iAnn, iSpec] )
    sites = [iCre, iAnn, iSpec][ord];        ranks = [jA, jD, AngularJ64(k)][ord]
    inter = [ranks[1], ranks[3], AngularJ64(0)]
    if  AngularMomentum.triangularDelta(inter[1], ranks[2], inter[2]) == 0            return( 0.0 )   end

    wR = treeRecoupling(leftCsf, rightCsf, sites, ranks, inter);   if  wR == 0.0      return( 0.0 )   end
    wA = wCre;                                                     if  wA == 0.0      return( 0.0 )   end
    wD = wAnn;                                                     if  wD == 0.0      return( 0.0 )   end
    wS = shellReducedW(jS, leftCsf.occupation[iSpec], leftCsf.seniorityNr[iSpec], leftCsf.subshellJ[iSpec],
                           rightCsf.seniorityNr[iSpec], rightCsf.subshellJ[iSpec], k)
    if  wS == 0.0                                                                    return( 0.0 )   end

    occup = 0
    for  i = min(iCre,iAnn):max(iCre,iAnn)-1    occup = occup + leftCsf.occupation[i]    end
    ph = (-1)^(occup + 1)
    if  iAnn < iCre
        ph = ph * (-1)^Int64( (Basics.twice(jA) + Basics.twice(jD) - 2*k + 2)//2 )
    end
    if  min(iCre,iAnn) < iSpec < max(iCre,iAnn)
        ph = ph * (-1)^Int64( (Basics.twice(jA) + Basics.twice(jD) + 2*k + 2)//2 )
    end

    return( ph * wR * wA * wD * wS / ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) ) )
end


const PARTNER_CACHE = Dict{NTuple{6,Int64}, Matrix{Float64}}()


"""
`SpinAngularNew.partnerMatrix(js::NTuple{4,AngularJ64}, nk::Int64, between::Bool)`
    ... to return the matrix M with M[k+1, K+1] = (2K+1) { j_a j_b k ; ... K } that turns a whole family of coefficients
        into its partner, so that the partner follows by ONE matrix-vector product instead of a six-j per (k, K).

        The point of caching it is that the matrix depends only on the FOUR SUBSHELL ANGULAR MOMENTA and not on the CSF
        pair, so the same one serves every pair of a calculation; a six-j costs about a microsecond and this family of
        them was being re-evaluated for each of thousands of pairs. A matrix::Matrix{Float64} is returned.
"""
function partnerMatrix(js::NTuple{4,AngularJ64}, nk::Int64, between::Bool)
    ja, jb, jc, jd = js
    key = (Basics.twice(ja), Basics.twice(jb), Basics.twice(jc), Basics.twice(jd), nk, between ? 1 : 0)
    haskey(PARTNER_CACHE, key)  &&  return( PARTNER_CACHE[key] )

    M = zeros(Float64, nk, nk)
    for  k = 0:nk-1,  K = 0:nk-1
        w6 = between ? AngularMomentum.Wigner_6j(ja, jb, AngularJ64(k), jc, jd, AngularJ64(K)) :
                       AngularMomentum.Wigner_6j(ja, jb, AngularJ64(k), jd, jc, AngularJ64(K))
        M[k+1, K+1] = (2.0*K + 1.0) * w6
    end
    PARTNER_CACHE[key] = M

    return( M )
end


"""
`SpinAngularNew.moveOnePartner(prim::Array{Float64,1}, js::NTuple{4,AngularJ64}, k::Int64, between::Bool)`
    ... to compute, at rank `k`, the coefficient of the OTHER of the two pairings of a one-electron-move term from the whole
        vector `prim` of primary coefficients, by the Racah sum that relates the two,

            V^k(a,b;d,c)  =  sum_K (2K+1) { j_a j_b k ; j_d j_c K } V^K(a,b;c,d) ,

        the same relation the equal-occupation exchange term uses. The bottom row is read on the PRIMARY quadruple's own
        c and d, so that the sum runs in the correct direction: with the spectator outside the primary is the direct
        pairing and the transform produces the crossed one, with the spectator between it is the other way round, and
        taking the two directions to be the same sum leaves 214 of 872 coefficients wrong. A value::Float64 is returned.
"""
function moveOnePartner(prim::Array{Float64,1}, js::NTuple{4,AngularJ64}, k::Int64, between::Bool)
    nk = length(prim)
    if  k >= nk
        M = partnerMatrix(js, k+1, between)
        wa = 0.0
        for  K = 0:nk-1    prim[K+1] != 0.0  &&  (wa = wa + M[k+1, K+1] * prim[K+1])    end
        return( wa )
    end
    M  = partnerMatrix(js, nk, between)
    wa = 0.0
    for  K = 0:nk-1    prim[K+1] != 0.0  &&  (wa = wa + M[k+1, K+1] * prim[K+1])    end

    return( wa )
end


"""
`SpinAngularNew.twoParticleMoveOne(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64)`
    ... to compute all two-particle coefficients of a CSF pair that differs by ONE electron, moved from subshell `iAnn`
        of the ket to subshell `iCre` of the bra. Every subshell that keeps its occupation and holds an electron acts in
        turn as the spectator on which the operator's second index pair sits, contributing two coefficients per rank --
        the primary from `SpinAngularNew.moveOnePrimary` and its partner from `SpinAngularNew.moveOnePartner`.

        A spectator that COINCIDES with the acceptor or the donor is included too, through
        `SpinAngularNew.moveOneSameShell`: three of the four one-electron operators then fall on a single subshell and
        the shell matrix element becomes the coupled `shellReducedAW` or `shellReducedWA`. So this method covers the
        one-electron-move case COMPLETELY, and nothing about it raises.

        A list coeffs::Array{Coefficient2p{EffectiveStrengthKind},1} is returned.
"""
function twoParticleMoveOne(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64)
    coeffs = Coefficient2p{EffectiveStrengthKind}[]
    # ... the creation on the acceptor and the annihilation on the donor depend on NEITHER the spectator NOR the rank,
    #     so they are formed once here rather than once per spectator and rank.
    jCre = AngularJ64( Basics.subshell_2j(subshells[iCre])//2 )
    jAnn = AngularJ64( Basics.subshell_2j(subshells[iAnn])//2 )
    wCre = shellReducedA(jCre, leftCsf.occupation[iCre],  leftCsf.seniorityNr[iCre],  leftCsf.subshellJ[iCre],
                               rightCsf.occupation[iCre], rightCsf.seniorityNr[iCre], rightCsf.subshellJ[iCre],
                               AngularM64(1//2))
    wAnn = shellReducedA(jAnn, leftCsf.occupation[iAnn],  leftCsf.seniorityNr[iAnn],  leftCsf.subshellJ[iAnn],
                               rightCsf.occupation[iAnn], rightCsf.seniorityNr[iAnn], rightCsf.subshellJ[iAnn],
                               AngularM64(-1//2))

    # ... the two SAME-SUBSHELL terms, where the spectator coincides with the acceptor or with the donor and three of
    #     the four one-electron operators fall on one subshell
    for  iSpec in [iCre, iAnn]
        if  iSpec == iCre  &&  rightCsf.occupation[iCre] < 1                continue    end
        if  iSpec == iAnn  &&  rightCsf.occupation[iAnn] < 2                continue    end
        ok = true
        for  i = 1:length(subshells)
            if  i == iCre  ||  i == iAnn    continue    end
            if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
                ok = false;    break
            end
        end
        if  !ok    continue    end
        q = moveOneSameShellQuad(subshells, iCre, iAnn, iSpec)
        for  k = 0:Basics.subshell_2j(subshells[iSpec])
            v = moveOneSameShell(leftCsf, rightCsf, subshells, iCre, iAnn, iSpec, k)
            abs(v) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, q[1], q[2], q[3], q[4], v))
        end
    end

    for  iSpec = 1:length(subshells)
        if  iSpec == iCre  ||  iSpec == iAnn                                        continue    end
        if  rightCsf.occupation[iSpec] < 1                                          continue    end
        if  leftCsf.occupation[iSpec] != rightCsf.occupation[iSpec]                  continue    end
        # ... every subshell that neither changes occupation nor is the spectator must keep its own coupling
        ok = true
        for  i = 1:length(subshells)
            if  i == iCre  ||  i == iAnn  ||  i == iSpec    continue    end
            if  leftCsf.subshellJ[i] != rightCsf.subshellJ[i]  ||  leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]
                ok = false;    break
            end
        end
        if  !ok    continue    end

        qP, qX, js, between = moveOneQuads(subshells, iCre, iAnn, iSpec)
        # ... the spectator tensor of rank k needs k <= 2 j_spec, and the pairing it joins bounds it by j_spec + j of
        #     the further subshell; the larger of the two is a safe upper limit and the triangle tests do the rest.
        kMax = Int64( ( Basics.subshell_2j(subshells[iSpec]) +
                        max(Basics.subshell_2j(subshells[iCre]), Basics.subshell_2j(subshells[iAnn])) )//2 )
        prim = moveOnePrimaryVector(leftCsf, rightCsf, subshells, iCre, iAnn, iSpec, kMax, wCre, wAnn)
        if  all(iszero, prim)    continue    end

        for  k = 0:kMax
            abs(prim[k+1]) > 1.0e-14  &&
                push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, qP[1], qP[2], qP[3], qP[4], prim[k+1]))
            vx = moveOnePartner(prim, js, k, between)
            abs(vx) > 1.0e-14  &&
                push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, qX[1], qX[2], qX[3], qX[4], vx))
        end
    end

    return( coeffs )
end


"""
`SpinAngularNew.shellReducedAW(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64, Jket::AngularJ64, kW::Int64, K::AngularJ64, mq::AngularM64)`
    ... to compute the reduced matrix element of the COUPLED product of a single creation or annihilation operator with
        a rank-kW tensor acting on the SAME subshell, <j^Nbra v J || (a^(j) x W^(kW))^(K) || j^Nket v' J'>. This is the
        object needed whenever three of a two-body operator's four one-electron operators fall on one subshell, i.e.
        when the spectator subshell coincides with the acceptor or the donor.

        It is built by CLOSURE over the intermediate subshell terms rather than from a new table: `a` carries the
        subshell from Nbra to Nket electrons, `W^(kW)` preserves that number, and the intermediate terms are exactly the
        ones the Gaigalas tables already enumerate -- `SpinAngular.qspaceTerms` gives each term its own index range
        `min_odd ... max_odd` for precisely this purpose. So the coefficients of fractional parentage are reused as DATA
        once more, and only the assembly is re-implemented:

            <bra||(a x W^(kW))^(K)||ket>  =  (-1)^(J_bra + J_ket + K) sqrt(2K+1)
                                             sum_r { j kW K ; J_ket J_bra J_r } <bra||a||r> <r||W^(kW)||ket>

        `mq` is +1/2 for creation and -1/2 for annihilation. A value::Float64 is returned.
"""
function shellReducedAW(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64,
                        Jket::AngularJ64, kW::Int64, K::AngularJ64, mq::AngularM64)
    SA   = JenaAtomicCalculator.SpinAngular
    Qbra = SA.qshellTermQ(j, senBra)
    if  AngularMomentum.triangularDelta(Jbra, K, Jket) == 0                return( 0.0 )   end
    # ... the triangle INEQUALITIES are not the whole rule: the three ranks must also sum to an integer, which
    #     `AngularMomentum.triangularDelta` does not test, so j = 5/2 with kW = 2 and K = 1 passes it and is
    #     nonetheless impossible.
    if  AngularMomentum.triangularDelta(j, AngularJ64(kW), K) == 0         return( 0.0 )   end
    if  isodd(Basics.twice(j) + 2*kW + Basics.twice(K))                    return( 0.0 )   end
    if  isodd(Basics.twice(Jbra) + Basics.twice(K) + Basics.twice(Jket))   return( 0.0 )   end
    ibra = SA.getTermNumber(j, Nbra, Qbra, Jbra)
    if  ibra >= 64                                                        return( 0.0 )   end
    bT   = SA.qspaceTerms(ibra)
    nu0  = Int64((Basics.twice(j) + 1)//2)
    # ... a term is only reachable at a given electron number when its quasispin PROJECTION fits inside its quasispin,
    #     |M_Q| <= Q; without that test the enumeration offers states that do not exist and the 3-j throws.
    if  SA.qspacedelta(Qbra, SA.qshellTermM(j, Nbra)) == 0                return( 0.0 )   end

    wa = 0.0
    for  ir = bT.min_odd:bT.max_odd
        rT  = SA.qspaceTerms(ir)
        if  Basics.twice(rT.j) != Basics.twice(j)                         continue        end
        senR = nu0 - Basics.twice(rT.Q)
        if  senR < 0  ||  senR > Nket                                     continue        end
        if  SA.qspacedelta(rT.Q, SA.qshellTermM(j, Nket)) == 0            continue        end
        w6 = AngularMomentum.Wigner_6j(j, AngularJ64(kW), K, Jket, Jbra, rT.J)
        if  w6 == 0.0                                                     continue        end
        wA = shellReducedA(j, Nbra, senBra, Jbra, Nket, senR, rT.J, mq)
        if  wA == 0.0                                                     continue        end
        wW = shellReducedW(j, Nket, senR, rT.J, senKet, Jket, kW)
        if  wW == 0.0                                                     continue        end
        wa = wa + w6 * wA * wW
    end
    # ... the phase is GLOBAL, (-1)^(J_bra + J_ket + K), and not one per intermediate term. Putting it inside the sum
    #     instead reproduced 4879 of 9043 values and corrupted the rest, since a wrong sign per term changes the sum by
    #     an amount that is not a constant factor -- which is exactly how it showed up: a tail of ratios like -2.14
    #     and -5.72 beside the correct ones, rather than a uniform -1.
    ph = Int64( (Basics.twice(Jbra) + Basics.twice(Jket) + Basics.twice(K))//2 )

    return( (-1)^ph * wa * sqrt(Basics.twice(K) + 1.0) )
end


"""
`SpinAngularNew.shellReducedWA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64, Jket::AngularJ64, kW::Int64, K::AngularJ64, mq::AngularM64)`
    ... to compute the reduced matrix element of the coupled product taken in the OTHER order,
        <j^Nbra v J || (W^(kW) x a^(j))^(K) || j^Nket v' J'>, where the single operator acts FIRST and the rank-kW tensor
        second, so that the intermediate subshell terms carry `Nbra` electrons rather than `Nket`.

        THIS IS NOT THE SAME OPERATOR AS `SpinAngularNew.shellReducedAW` WITH A PHASE. The reordering identity
        (A x B)^K = (-1)^(k_A + k_B - K) (B x A)^K holds for tensors that commute, and `a` and `W = (a^+ a)` on the SAME
        subshell do not. The two are genuinely different objects, which is why the predecessor carries two schemes for
        them, and the difference is invisible until the subshell is FULL on one side: with a filled j = 1/2 donor shell
        in the ket, W^(1) on it vanishes because a closed shell is a scalar, so the ket-side intermediate kills a rank-1
        term that physically exists. That is exactly how the need for this routine showed up.

        A value::Float64 is returned.
"""
function shellReducedWA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64,
                        Jket::AngularJ64, kW::Int64, K::AngularJ64, mq::AngularM64)
    SA   = JenaAtomicCalculator.SpinAngular
    Qket = SA.qshellTermQ(j, senKet)
    if  AngularMomentum.triangularDelta(Jbra, K, Jket) == 0                return( 0.0 )   end
    if  AngularMomentum.triangularDelta(j, AngularJ64(kW), K) == 0         return( 0.0 )   end
    if  isodd(Basics.twice(j) + 2*kW + Basics.twice(K))                    return( 0.0 )   end
    if  isodd(Basics.twice(Jbra) + Basics.twice(K) + Basics.twice(Jket))   return( 0.0 )   end
    iket = SA.getTermNumber(j, Nket, Qket, Jket)
    if  iket >= 64                                                        return( 0.0 )   end
    kT   = SA.qspaceTerms(iket)
    nu0  = Int64((Basics.twice(j) + 1)//2)
    if  SA.qspacedelta(Qket, SA.qshellTermM(j, Nket)) == 0                return( 0.0 )   end

    wa = 0.0
    for  ir = kT.min_odd:kT.max_odd
        rT  = SA.qspaceTerms(ir)
        if  Basics.twice(rT.j) != Basics.twice(j)                         continue        end
        senR = nu0 - Basics.twice(rT.Q)
        if  senR < 0  ||  senR > Nbra                                     continue        end
        if  SA.qspacedelta(rT.Q, SA.qshellTermM(j, Nbra)) == 0            continue        end
        w6 = AngularMomentum.Wigner_6j(AngularJ64(kW), j, K, Jket, Jbra, rT.J)
        if  w6 == 0.0                                                     continue        end
        wW = shellReducedW(j, Nbra, senBra, Jbra, senR, rT.J, kW)
        if  wW == 0.0                                                     continue        end
        wA = shellReducedA(j, Nbra, senR, rT.J, Nket, senKet, Jket, mq)
        if  wA == 0.0                                                     continue        end
        wa = wa + w6 * wW * wA
    end
    ph = Int64( (Basics.twice(Jbra) + Basics.twice(Jket) + Basics.twice(K))//2 )

    return( (-1)^ph * wa * sqrt(Basics.twice(K) + 1.0) )
end


"""
`SpinAngularNew.moveOneSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64, k::Int64)`
    ... to compute, at rank `k`, the two-particle coefficient of a CSF pair differing by ONE electron in the case where
        the spectator subshell COINCIDES with the acceptor `iCre` or the donor `iAnn`. Three of the operator's four
        one-electron operators then fall on that one subshell -- two creations and an annihilation when it is the
        acceptor, a creation and two annihilations when it is the donor -- so the shell matrix element is the coupled
        `SpinAngularNew.shellReducedAW` rather than a plain creation, and only TWO subshells act in total.

        Because the two site ranks must couple to zero, the coupled rank K is fixed to the OTHER subshell's j, and there
        is no sum. Two phases beyond the Jordan-Wigner string, both established on data:

        (1) the Jordan-Wigner string over the occupations between the two subshells;
        (2) the reordering phase (-1)^(j_spec + k - K), carried in BOTH cases.

        A value::Float64 is returned.
"""
function moveOneSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64,
                          iSpec::Int64, k::Int64)
    iOth = (iSpec == iCre) ? iAnn : iCre
    jS   = AngularJ64( Basics.subshell_2j(subshells[iSpec])//2 )
    jO   = AngularJ64( Basics.subshell_2j(subshells[iOth])//2 )
    K    = jO
    mqS  = (iSpec == iCre) ? AngularM64(1//2)  : AngularM64(-1//2)
    mqO  = (iOth  == iCre) ? AngularM64(1//2)  : AngularM64(-1//2)

    # ... W GOES ON THE SIDE THAT HOLDS FEWER ELECTRONS, which decides the ordering: (a x W) when the spectator is the
    #     acceptor, (W x a) when it is the donor. The two are different operators, not one with a phase, and taking
    #     (a x W) for the donor is right until the donor subshell is CLOSED in the ket -- W^(k>0) on a closed shell
    #     vanishes, so a rank-1 term that physically exists is silently dropped.
    wS = (iSpec == iCre) ?
         shellReducedAW(jS, leftCsf.occupation[iSpec],  leftCsf.seniorityNr[iSpec],  leftCsf.subshellJ[iSpec],
                            rightCsf.occupation[iSpec], rightCsf.seniorityNr[iSpec], rightCsf.subshellJ[iSpec],
                            k, K, mqS) :
         shellReducedWA(jS, leftCsf.occupation[iSpec],  leftCsf.seniorityNr[iSpec],  leftCsf.subshellJ[iSpec],
                            rightCsf.occupation[iSpec], rightCsf.seniorityNr[iSpec], rightCsf.subshellJ[iSpec],
                            k, K, mqS)
    if  wS == 0.0                                                                        return( 0.0 )   end
    wO = shellReducedA(jO, leftCsf.occupation[iOth],  leftCsf.seniorityNr[iOth],  leftCsf.subshellJ[iOth],
                           rightCsf.occupation[iOth], rightCsf.seniorityNr[iOth], rightCsf.subshellJ[iOth],
                           mqO);                                           if  wO == 0.0   return( 0.0 )   end

    lo, hi = minmax(iSpec, iOth)
    ranks  = (lo == iSpec) ? [K, jO] : [jO, K]
    wR     = treeRecoupling(leftCsf, rightCsf, [lo,hi], ranks, [ranks[1], AngularJ64(0)])
    if  wR == 0.0                                                                        return( 0.0 )   end

    occup = 0
    for  i = lo:hi-1    occup = occup + leftCsf.occupation[i]    end
    ph = (-1)^(occup + 1) * (-1)^Int64( (Basics.twice(jS) + 2*k - Basics.twice(K))//2 )

    return( ph * wR * wS * wO / ( sqrt(2.0*k + 1.0) * sqrt(Basics.twice(leftCsf.J) + 1.0) ) )
end


"""
`SpinAngularNew.moveOneSameShellQuad(subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64)`
    ... to determine the index quadruple (a,b,c,d) under which the same-subshell one-electron-move term is emitted. The
        rule is the one every other topology in this module follows: `a` and `b` are the two CREATIONS in subshell-index
        order and `c` and `d` the two ANNIHILATIONS in index order. With the spectator on the acceptor the creations are
        (S,S) and the annihilations (S,O); with it on the donor the creations are (O,S) and the annihilations (S,S). The
        two orders of `c` and `d` are the SAME integral here, since X^k(abcd) = X^k(badc) and a equals b, so there is one
        family and not two. A tuple (a,b,c,d)::NTuple{4,Subshell} is returned.
"""
function moveOneSameShellQuad(subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64, iSpec::Int64)
    iOth = (iSpec == iCre) ? iAnn : iCre
    sS   = subshells[iSpec];    sO = subshells[iOth]
    if  iSpec == iCre    return( iSpec < iOth ? (sS, sS, sS, sO) : (sS, sS, sO, sS) )
    else                 return( iOth  < iSpec ? (sO, sS, sS, sS) : (sS, sO, sS, sS) )
    end
end


"""
`SpinAngularNew.moveTwoWeight(jc::NTuple{2,AngularJ64}, ja::NTuple{2,AngularJ64}, jsite::NTuple{4,AngularJ64}, swapped::NTuple{2,Bool}, crossPair::Bool, R::AngularJ64, k::Int64)`
    ... to compute the weight with which one value `R` of the chain's free intermediate rank enters a four-subshell
        two-electron-move coefficient of rank `k`. Four subshells act, each with one operator, and the physical operator
        pairs each creation with its own annihilation -- a pairing that in general CROSSES the subshell chain, which is
        why a sum over `R` appears here where the one-electron cases needed only a phase.

        Only TWO patterns can occur for a given family, decided by where the two physical pairs sit among the four
        subshells in index order:

        + `:chain`   ... the pairs are (1,2) and (3,4), already the chain's own pairing, so R is forced to k and the
                         weight is 1/sqrt(2k+1);
        + `:cross13` ... the pairs are (1,3) and (2,4); the nine-j with a zero column collapses to
                         (-1)^(j_2 + j_3 + R + k + 1) sqrt(2R+1) { j_1 j_2 R ; j_4 j_3 k };
        + `:cross14` ... the pairs are (1,4) and (2,3), giving (-1)^(j_2 + j_3 + k + 1) sqrt(2R+1) { j_1 j_2 R ; j_3 j_4 k },
                         whose phase carries NO R, because exchanging j_3 and j_4 in the previous line cancels it.

        On top of that each physical pair whose ANNIHILATION sits on the lower subshell index contributes the ordering
        phase (-1)^(j_cre + j_ann + k + 1) -- the same phase the one-particle substitution and the one-electron move
        already carry, so it is one rule across the module rather than three. A value::Float64 is returned.
"""
function moveTwoWeight(jc::NTuple{2,AngularJ64}, ja::NTuple{2,AngularJ64}, jsite::NTuple{4,AngularJ64},
                       swapped::NTuple{2,Bool}, pattern::Symbol, R::AngularJ64, k::Int64)
    wa = if      pattern == :chain
                 Basics.twice(R) == 2*k ? 1.0/sqrt(2.0*k + 1.0) : 0.0
         elseif  pattern == :cross13
                 (-1)^Int64( (Basics.twice(jsite[2]) + Basics.twice(jsite[3]) + Basics.twice(R) + 2*k + 2)//2 ) *
                 sqrt(Basics.twice(R) + 1.0) *
                 AngularMomentum.Wigner_6j(jsite[1], jsite[2], R, jsite[4], jsite[3], AngularJ64(k))
         else
                 # ... exchanging j_3 and j_4 in the :cross13 form brings a factor (-1)^(j_3 + j_4 - R), and the R in it
                 #     CANCELS the R in that form's phase -- so this one alone has no R in its sign. Every candidate I
                 #     tried with an R there failed; the algebra says why.
                 (-1)^Int64( (Basics.twice(jsite[2]) + Basics.twice(jsite[3]) + 2*k + 2)//2 ) *
                 sqrt(Basics.twice(R) + 1.0) *
                 AngularMomentum.Wigner_6j(jsite[1], jsite[2], R, jsite[3], jsite[4], AngularJ64(k))
         end
    if  wa == 0.0                                                         return( 0.0 )   end

    for  m = 1:2
        if  swapped[m]
            wa = wa * (-1)^Int64( (Basics.twice(jc[m]) + Basics.twice(ja[m]) + 2*k + 2)//2 )
        end
    end

    return( wa )
end


"""
`SpinAngularNew.moveTwoRecoupling(leftCsf::CsfR, rightCsf::CsfR, sites::Array{Int64,1}, jsite::NTuple{4,AngularJ64}, R::AngularJ64)`
    ... to compute the RECOUPLING alone of a four-subshell two-electron move at one value `R` of the free intermediate
        rank. Only this factor depends on R; the four one-electron matrix elements and the Jordan-Wigner string do not,
        and `SpinAngularNew.twoParticleMoveTwo` forms them once rather than once per R. For four operators that string
        is a PREFIX count -- the occupations lying before each acting subshell -- and not a count of the gaps between
        them; the two differ as soon as the creations and annihilations interleave, which they do in half of the
        arrangements. A value::Float64 is returned.
"""
function moveTwoRecoupling(leftCsf::CsfR, rightCsf::CsfR, sites::Array{Int64,1}, jsite::NTuple{4,AngularJ64},
                          R::AngularJ64)
    if  AngularMomentum.triangularDelta(jsite[1], jsite[2], R) == 0        return( 0.0 )   end
    if  AngularMomentum.triangularDelta(R, jsite[3], jsite[4]) == 0        return( 0.0 )   end
    ranks = [ jsite[1], jsite[2], jsite[3], jsite[4] ]
    inter = [ jsite[1], R, jsite[4], AngularJ64(0) ]

    return( treeRecoupling(leftCsf, rightCsf, sites, ranks, inter) )
end


"""
`SpinAngularNew.twoParticleMoveTwo(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, cre::Array{Int64,1}, ann::Array{Int64,1})`
    ... to compute all two-particle coefficients of a CSF pair that differs by TWO electrons moved between FOUR distinct
        subshells, `cre` holding the two that gain an electron and `ann` the two that lose one, each in increasing
        subshell order. Every subshell carries exactly one operator, so no coupled shell tensor is needed and the whole
        assembly is `SpinAngularNew.moveTwoChain` summed over the chain's free intermediate rank with the weights of
        `SpinAngularNew.moveTwoWeight`.

        The DIRECT family pairs each creation with the annihilation of the same rank order and is built directly; its
        partner, with the two annihilations exchanged, follows by the same Racah sum the equal-occupation exchange term
        and the one-electron move already use, so the third possible pairing pattern never has to be derived.

        A list coeffs::Array{Coefficient2p{EffectiveStrengthKind},1} is returned.
"""
function twoParticleMoveTwo(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, cre::Array{Int64,1},
                            ann::Array{Int64,1})
    coeffs = Coefficient2p{EffectiveStrengthKind}[]
    sites  = sort( vcat(cre, ann) )
    mqs    = [ (i in cre) ? AngularM64(1//2) : AngularM64(-1//2)  for i in sites ]
    posOf  = Dict(s => n for (n,s) in enumerate(sites))
    jsite  = Tuple( AngularJ64( Basics.subshell_2j(subshells[i])//2 )  for i in sites )

    Rs = AngularJ64[]
    for  t = 0:Basics.twice(jsite[1]) + Basics.twice(jsite[2])
        R = AngularJ64(t//2)
        if  AngularMomentum.triangularDelta(jsite[1], jsite[2], R) != 0  &&
            AngularMomentum.triangularDelta(R, jsite[3], jsite[4]) != 0     push!(Rs, R)   end
    end
    if  isempty(Rs)                                                       return( coeffs )   end
    # ... the four one-electron matrix elements and the Jordan-Wigner string do NOT depend on R, so they are formed
    #     ONCE here instead of inside the loop over the intermediate rank, where they were being recomputed nR times.
    shell = 1.0
    for  (n,i) in enumerate(sites)
        shell = shell * shellReducedA(jsite[n], leftCsf.occupation[i],  leftCsf.seniorityNr[i],  leftCsf.subshellJ[i],
                                                rightCsf.occupation[i], rightCsf.seniorityNr[i], rightCsf.subshellJ[i],
                                                mqs[n])
        if  shell == 0.0                                                  return( coeffs )   end
    end
    jw = 0
    for  i in sites    jw = jw + sum(leftCsf.occupation[1:i-1]; init = 0)    end
    shell = (-1)^jw * shell / sqrt(Basics.twice(leftCsf.J) + 1.0)

    chain = [ moveTwoRecoupling(leftCsf, rightCsf, sites, jsite, R) * shell  for R in Rs ]
    if  all(iszero, chain)                                                return( coeffs )   end

    # ... the two families are INDEPENDENT contractions, not an exchange-related pair: with four distinct subshells
    #     R^k(c1,c2,a1,a2) and R^k(c1,c2,a2,a1) are different integrals. Each is assembled with its own pairing of
    #     creations to annihilations, and its own pattern against the subshell chain.
    for  (aOrder, qd) in [ ([ann[1], ann[2]], (subshells[cre[1]], subshells[cre[2]], subshells[ann[1]], subshells[ann[2]])),
                           ([ann[2], ann[1]], (subshells[cre[1]], subshells[cre[2]], subshells[ann[2]], subshells[ann[1]])) ]
        jc = ( AngularJ64(Basics.subshell_2j(subshells[cre[1]])//2), AngularJ64(Basics.subshell_2j(subshells[cre[2]])//2) )
        ja = ( AngularJ64(Basics.subshell_2j(subshells[aOrder[1]])//2), AngularJ64(Basics.subshell_2j(subshells[aOrder[2]])//2) )
        p1 = Set([posOf[cre[1]], posOf[aOrder[1]]]);   p2 = Set([posOf[cre[2]], posOf[aOrder[2]]])
        pattern = Set([p1,p2]) == Set([Set([1,2]),Set([3,4])]) ? :chain :
                  Set([p1,p2]) == Set([Set([1,3]),Set([2,4])]) ? :cross13 : :cross14
        swapped = ( posOf[aOrder[1]] < posOf[cre[1]], posOf[aOrder[2]] < posOf[cre[2]] )
        kMax = Int64( min(Basics.subshell_2j(subshells[cre[1]]) + Basics.subshell_2j(subshells[aOrder[1]]),
                          Basics.subshell_2j(subshells[cre[2]]) + Basics.subshell_2j(subshells[aOrder[2]]))//2 )
        for  k = 0:kMax
            v = 0.0
            for  (n,R) in enumerate(Rs)
                chain[n] == 0.0  &&  continue
                v = v + moveTwoWeight(jc, ja, jsite, swapped, pattern, R, k) * chain[n]
            end
            abs(v) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, qd[1], qd[2], qd[3], qd[4], v))
        end
    end

    return( coeffs )
end


"""
`SpinAngularNew.shellReducedAA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64, Jket::AngularJ64, K::AngularJ64, mq::AngularM64)`
    ... to compute the reduced matrix element of TWO creation or TWO annihilation operators coupled on one subshell,
        <j^Nbra v J || (a^(j) x a^(j))^(K) || j^Nket v' J'>, with `mq` = +1/2 for two creations and -1/2 for two
        annihilations, so that Nbra and Nket differ by two. This is the object a two-electron move needs whenever BOTH
        electrons leave, or both arrive in, the same subshell.

        Built by the same closure as `SpinAngularNew.shellReducedAW`, over the intermediate subshell terms of j^(Nket+2mq)
        that the Gaigalas tables already enumerate:

            <bra||(a x a)^(K)||ket>  =  (-1)^(J_bra + J_ket + K) sqrt(2K+1)
                                        sum_r { j j K ; J_ket J_bra J_r } <bra||a||r> <r||a||ket>

        A value::Float64 is returned.
"""
function shellReducedAA(j::AngularJ64, Nbra::Int64, senBra::Int64, Jbra::AngularJ64, Nket::Int64, senKet::Int64,
                        Jket::AngularJ64, K::AngularJ64, mq::AngularM64)
    SA   = JenaAtomicCalculator.SpinAngular
    if  AngularMomentum.triangularDelta(Jbra, K, Jket) == 0                return( 0.0 )   end
    if  AngularMomentum.triangularDelta(j, j, K) == 0                      return( 0.0 )   end
    if  isodd(Basics.twice(Jbra) + Basics.twice(K) + Basics.twice(Jket))   return( 0.0 )   end
    Qbra = SA.qshellTermQ(j, senBra)
    if  SA.qspacedelta(Qbra, SA.qshellTermM(j, Nbra)) == 0                 return( 0.0 )   end
    ibra = SA.getTermNumber(j, Nbra, Qbra, Jbra)
    if  ibra >= 64                                                         return( 0.0 )   end
    bT   = SA.qspaceTerms(ibra)
    nu0  = Int64((Basics.twice(j) + 1)//2)
    Nmid = Nket + Basics.twice(mq)

    wa = 0.0
    for  ir = bT.min_odd:bT.max_odd
        rT = SA.qspaceTerms(ir)
        if  Basics.twice(rT.j) != Basics.twice(j)                          continue        end
        senR = nu0 - Basics.twice(rT.Q)
        if  senR < 0  ||  senR > Nmid                                      continue        end
        if  SA.qspacedelta(rT.Q, SA.qshellTermM(j, Nmid)) == 0             continue        end
        w6 = AngularMomentum.Wigner_6j(j, j, K, Jket, Jbra, rT.J)
        if  w6 == 0.0                                                      continue        end
        w1 = shellReducedA(j, Nbra, senBra, Jbra, Nmid, senR, rT.J, mq)
        if  w1 == 0.0                                                      continue        end
        w2 = shellReducedA(j, Nmid, senR, rT.J, Nket, senKet, Jket, mq)
        if  w2 == 0.0                                                      continue        end
        wa = wa + w6 * w1 * w2
    end
    ph = Int64( (Basics.twice(Jbra) + Basics.twice(Jket) + Basics.twice(K))//2 )

    return( (-1)^ph * wa * sqrt(Basics.twice(K) + 1.0) )
end


"""
`SpinAngularNew.twoParticleMoveTwoSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, dbl::Int64, sgl::Array{Int64,1}, diff::Array{Int64,1})`
    ... to compute all two-particle coefficients of a CSF pair that differs by TWO electrons where BOTH of them leave, or
        both arrive in, the SAME subshell `dbl`, the other two operators sitting on the single subshells `sgl`. The
        doubled subshell then carries the coupled two-creation or two-annihilation tensor `shellReducedAA`, and only
        three subshells act.

        THE WEIGHT IS THE SAME ONE AS THE FOUR-SUBSHELL cross13 CASE, and that is a derivation rather than a
        coincidence: the physical pairing (S,x1),(S,x2) occupies slots (1,3),(2,4) of the quadruple (S,S,x1,x2), which
        is exactly the pattern that transformation was derived for. It gives the exact MAGNITUDE of every coefficient
        immediately.

        Three signs complete it, and each was found by reading what the predecessor does rather than by fitting:

        (1) the JORDAN-WIGNER string runs BETWEEN THE TWO SINGLE SUBSHELLS only -- not over all three acting subshells,
            which is what a prefix count would give and is what the four-subshell case needs;
        (2) which single subshell enters the weight's phase depends on WHICH KIND of operator is doubled: the first when
            the doubled subshell creates, the second when it annihilates. This is the `irez` distinction of
            `SpinAngular.twoParticle15to18order`, and without it half of the doubled-annihilation coefficients come out
            with the wrong sign while every doubled-creation one is right;
        (3) a RE-PAIRING phase (-1)^(j_x1 + j_x2 + 1) when the doubled subshell lies BETWEEN the two single ones,
            because the coupling chain then joins it to a single subshell first where the operator couples the two
            singles together. It depends on the two single subshells' j alone, as a re-pairing of three ranks coupled to
            zero must, and it was verified on 1210 coefficients.

        A list coeffs::Array{Coefficient2p{EffectiveStrengthKind},1} is returned.
"""
function twoParticleMoveTwoSameShell(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, dbl::Int64,
                                     sgl::Array{Int64,1}, diff::Array{Int64,1})
    coeffs = Coefficient2p{EffectiveStrengthKind}[]
    jS     = AngularJ64( Basics.subshell_2j(subshells[dbl])//2 )
    jx     = [ AngularJ64( Basics.subshell_2j(subshells[i])//2 )  for i in sgl ]
    isCre  = diff[dbl] > 0
    mqD    = isCre ? AngularM64(1//2) : AngularM64(-1//2)
    sites  = sort( vcat([dbl], sgl) )

    # ... the one-electron matrix elements on the two single subshells do not depend on the rank, so they come first
    wSingle = 1.0
    for  i in sgl
        wSingle = wSingle *
            shellReducedA(AngularJ64(Basics.subshell_2j(subshells[i])//2),
                          leftCsf.occupation[i],  leftCsf.seniorityNr[i],  leftCsf.subshellJ[i],
                          rightCsf.occupation[i], rightCsf.seniorityNr[i], rightCsf.subshellJ[i],
                          diff[i] > 0 ? AngularM64(1//2) : AngularM64(-1//2))
        if  wSingle == 0.0                                                return( coeffs )   end
    end

    jw = sum(leftCsf.occupation[min(sgl[1],sgl[2]) : max(sgl[1],sgl[2])-1]; init = 0)
    ph = (-1)^(jw + 1)
    if  min(sgl[1],sgl[2]) < dbl < max(sgl[1],sgl[2])
        ph = ph * (-1)^Int64( (Basics.twice(jx[1]) + Basics.twice(jx[2]) + 2)//2 )
    end
    jPh  = isCre ? jx[1] : jx[2]
    cre  = isCre ? (subshells[dbl], subshells[dbl]) : (subshells[sgl[1]], subshells[sgl[2]])
    ann  = isCre ? (subshells[sgl[1]], subshells[sgl[2]]) : (subshells[dbl], subshells[dbl])
    kMax = Int64( min(Basics.subshell_2j(cre[1]) + Basics.subshell_2j(ann[1]),
                      Basics.subshell_2j(cre[2]) + Basics.subshell_2j(ann[2]))//2 )

    for  k = 0:kMax
        tot = 0.0
        for  t = 0:2*Basics.twice(jS)
            if  isodd(t)                                                  continue   end
            K = AngularJ64(t//2)
            if  AngularMomentum.triangularDelta(jS, jS, K) == 0           continue   end
            ranks = [ i == dbl ? K : AngularJ64(Basics.subshell_2j(subshells[i])//2)  for i in sites ]
            inter = [ ranks[1], ranks[3], AngularJ64(0) ]
            if  AngularMomentum.triangularDelta(inter[1], ranks[2], inter[2]) == 0    continue   end
            wR = treeRecoupling(leftCsf, rightCsf, sites, ranks, inter);   if  wR == 0.0   continue   end
            wD = shellReducedAA(jS, leftCsf.occupation[dbl],  leftCsf.seniorityNr[dbl],  leftCsf.subshellJ[dbl],
                                    rightCsf.occupation[dbl], rightCsf.seniorityNr[dbl], rightCsf.subshellJ[dbl],
                                    K, mqD);                               if  wD == 0.0   continue   end
            w6 = AngularMomentum.Wigner_6j(jS, jS, K, jx[2], jx[1], AngularJ64(k))
            if  w6 == 0.0                                                             continue   end
            phK = Int64( (Basics.twice(jS) + Basics.twice(jPh) + Basics.twice(K) + 2*k + 2)//2 )
            tot = tot + (-1)^phK * sqrt(Basics.twice(K) + 1.0) * w6 * wR * wD
        end
        v = ph * wSingle * tot / sqrt(Basics.twice(leftCsf.J) + 1.0)
        abs(v) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, cre[1], cre[2], ann[1], ann[2], v))
    end

    return( coeffs )
end


"""
`SpinAngularNew.twoParticleMoveTwoBothDoubled(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64, iAnn::Int64)`
    ... to compute all two-particle coefficients of a CSF pair in which TWO electrons leave one subshell `iAnn` and TWO
        arrive in another `iCre`, so that only TWO subshells act and each carries a coupled pair tensor: `shellReducedAA`
        with two creations on one and two annihilations on the other. Their ranks must be equal, since the operator is a
        scalar, so there is one sum over that common rank K and no free intermediate.

        The weight is once more the transformation already derived for the four-subshell cross13 case -- the physical
        pairing (S1,S2),(S1,S2) occupies slots (1,3),(2,4) of the quadruple (S1,S1,S2,S2) -- so nothing new is fitted.

        THE ONE THING THIS CASE NEEDS THAT NO OTHER DOES IS A FACTOR OF ONE HALF, and it is not a fudge. The two-body
        operator carries a 1/2, which in every other topology is cancelled because the two Wick contractions are DISTINCT
        and contribute equally. Here both physical pairs are (S1,S2): the two contractions COINCIDE, nothing cancels the
        1/2, and it survives into the coefficient. Leaving it out makes every coefficient of this class exactly twice too
        large with the sign already right, which is how it was found.

        A list coeffs::Array{Coefficient2p{EffectiveStrengthKind},1} is returned.
"""
function twoParticleMoveTwoBothDoubled(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1}, iCre::Int64,
                                       iAnn::Int64)
    coeffs = Coefficient2p{EffectiveStrengthKind}[]
    jC     = AngularJ64( Basics.subshell_2j(subshells[iCre])//2 )
    jA     = AngularJ64( Basics.subshell_2j(subshells[iAnn])//2 )
    lo, hi = minmax(iCre, iAnn);        sites = [lo, hi]
    kMax   = Int64( (Basics.subshell_2j(subshells[iCre]) + Basics.subshell_2j(subshells[iAnn]))//2 )

    for  k = 0:kMax
        tot = 0.0
        for  t = 0:2*min(Basics.twice(jC), Basics.twice(jA))
            if  isodd(t)                                                  continue   end
            K = AngularJ64(t//2)
            if  AngularMomentum.triangularDelta(jC, jC, K) == 0           continue   end
            if  AngularMomentum.triangularDelta(jA, jA, K) == 0           continue   end
            wR = treeRecoupling(leftCsf, rightCsf, sites, [K, K], [K, AngularJ64(0)])
            if  wR == 0.0                                                 continue   end
            wC = shellReducedAA(jC, leftCsf.occupation[iCre],  leftCsf.seniorityNr[iCre],  leftCsf.subshellJ[iCre],
                                    rightCsf.occupation[iCre], rightCsf.seniorityNr[iCre], rightCsf.subshellJ[iCre],
                                    K, AngularM64(1//2));                  if  wC == 0.0   continue   end
            wA = shellReducedAA(jA, leftCsf.occupation[iAnn],  leftCsf.seniorityNr[iAnn],  leftCsf.subshellJ[iAnn],
                                    rightCsf.occupation[iAnn], rightCsf.seniorityNr[iAnn], rightCsf.subshellJ[iAnn],
                                    K, AngularM64(-1//2));                 if  wA == 0.0   continue   end
            w6 = AngularMomentum.Wigner_6j(jC, jC, K, jA, jA, AngularJ64(k))
            if  w6 == 0.0                                                 continue   end
            ph  = Int64( (Basics.twice(jC) + Basics.twice(jA) + Basics.twice(K) + 2*k + 2)//2 )
            tot = tot + (-1)^ph * sqrt(Basics.twice(K) + 1.0) * w6 * wR * wC * wA
        end
        v = 0.5 * tot / sqrt(Basics.twice(leftCsf.J) + 1.0)
        abs(v) > 1.0e-14  &&  push!(coeffs, Coefficient2p{EffectiveStrengthKind}(k, subshells[iCre], subshells[iCre],
                                                                                subshells[iAnn], subshells[iAnn], v))
    end

    return( coeffs )
end
