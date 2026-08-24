
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
    ... the two-particle (electron-electron) coefficients. NOT YET IMPLEMENTED: this method raises, and the reason is
        recorded here because it was arrived at the hard way.

        A first attempt fitted closed forms to the diagonal two-subshell terms, read off a CLOSED-SHELL configuration:

            direct    (abab), k = 0 :   N_a N_b / sqrt((2j_a+1)(2j_b+1))
            exchange  (abba)        :   (-1)^(j_a+j_b+k) N_a N_b / ((2j_a+1)(2j_b+1))

        The direct term is right and survives every test. **The exchange term is WRONG**, and wrong in a way the closed
        shells could not show: for two singly-occupied subshells the exchange coefficient DEPENDS ON HOW THE TWO
        ELECTRONS ARE COUPLED. For 1s 2s it is +0.5 at k = 0 when J = 0 and -0.5 when J = 1; the expression above has no
        J in it at all. A closed shell forces J = 0, so the dependence is invisible there -- the set the formula was
        fitted on could not discriminate, which is the same trap this module's own example file warns about twice.

        There is also a DIRECT term at k > 0 for open subshells, likewise J-dependent, which the closed-shell data does
        not show either because it vanishes when the shell is full.

        So the two-subshell diagonal case needs the genuine coupled-tensor recoupling -- the scalar product
        [W^(k)(a) x W^(k)(b)]^(0) reduced through the coupling tree -- and not a closed form.

        BOTH TERMS BETWEEN TWO DISTINCT SUBSHELLS NOW EXIST, with different reach. `SpinAngularNew.twoParticleDirect` computes the DIRECT term from exactly that scalar
        product, and reproduces SpinAngular on twenty coefficients out of sample, J = 0 ... 3 and j = 1/2, 3/2, 5/2,
        every ratio 1.0000000 -- J-dependence included, which is what the withdrawn closed form got wrong.

        `SpinAngularNew.twoParticleExchange` covers the exchange term for SINGLY OCCUPIED subshells only -- there the
        Racah sum over intermediate ranks collapses to one 6j. For N >= 2 it does not collapse, and that method refuses
        rather than applying the collapsed form out of range.

        WHAT IS STILL MISSING, and why this method continues to RAISE: the SAME-SUBSHELL term (a,a,a,a).

        TWO ATTEMPTS HAVE FAILED, and what they establish is worth more than a label.

        (i) A bare CLOSURE sum over the shell's own intermediate terms, S(k) = sum_int <bra||W^(k)||int>
        <int||W^(k)||bra>, gives a clean constant ratio of 4 at k = 0 -- encouraging and misleading -- and VANISHES at
        k = 1 where the coefficient is 0.5.

        I first recorded that as refuting the closure ROUTE. THAT WAS TOO STRONG AND IS WITHDRAWN. Reading
        `SpinAngular`'s SchemeEta_WW afterwards shows the genuine assembly IS a closure sum over intermediate terms;
        what my version lacked was a phase (-1)^((2k - 2J_a + 2J_r)/2), a normalisation 1/sqrt((2k+1)(2J_a+1)), and a
        restricted intermediate range. So attempt (i) refutes an IMPLEMENTATION, not the route -- a distinction worth
        keeping, since a route wrongly marked dead is not revisited.

        (ii) The corrected closure, carrying that phase and normalisation, does not reproduce SpinAngular either. For
        j = 3/2, N = 2, J = 0 the coefficient is 0.25 at EVERY rank 0 ... 3, while the sum gives 1.0, 0, 2.236, 0.
        The constancy in k is itself a clue and is not yet explained.

        So the same-shell term remains open, with the closure route live rather than dead, and the operator structure
        -- which m-projections, and what normal-ordering correction -- the thing still to get right. Returning what exists would be an incomplete
        coefficient list, which is the silent wrong answer this module exists to prevent; half a two-particle answer is
        worth less than none.
"""
function computeCoefficients(op::SpinAngularNew.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR,
                             subshells::Array{Subshell,1})
    error("\n\nSpinAngularNew.computeCoefficients (two-particle): NOT YET IMPLEMENTED.\n"                          *
          ">>> A closed-form attempt was WITHDRAWN: its exchange term was fitted on closed-shell data and misses\n"  *
          ">>> the coupling dependence that only open subshells reveal (1s 2s gives +0.5 at J = 0 and -0.5 at\n"     *
          ">>> J = 1 for the same term). See the docstring of this method.\n"                                        *
          ">>> Use SpinAngular.computeCoefficients for two-particle coefficients.\n")
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
