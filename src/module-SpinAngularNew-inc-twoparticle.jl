
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

        HALF OF THAT IS NOW DONE. `SpinAngularNew.twoParticleDirect` computes the DIRECT term from exactly that scalar
        product, and reproduces SpinAngular on twenty coefficients out of sample, J = 0 ... 3 and j = 1/2, 3/2, 5/2,
        every ratio 1.0000000 -- J-dependence included, which is what the withdrawn closed form got wrong.

        THE EXCHANGE TERM IS STILL OPEN. A single 6j {j_a j_b J; j_b j_a k} was tried and REFUTED: it vanishes at ranks
        where the coefficient does not, so the exchange channel needs the full Racah sum over intermediate ranks rather
        than one symbol. Until that exists this method RAISES, because returning direct terms alone would be an
        incomplete coefficient list -- the silent wrong answer this module exists to prevent.
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
