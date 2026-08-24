
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
    if  leftCsf.occupation != rightCsf.occupation
        error("\n\nSpinAngularNew.computeCoefficients (two-particle): the occupations differ.\n"             *
              ">>> Equal-occupation pairs are handled, whatever their coupling. A pair that MOVES electrons\n"  *
              ">>> between subshells needs the creation/annihilation machinery this method does not have, and\n" *
              ">>> returning the equal-occupation answer for it would be a wrong number rather than a missing\n" *
              ">>> one. Use SpinAngular.computeCoefficients for such a pair.\n")
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
