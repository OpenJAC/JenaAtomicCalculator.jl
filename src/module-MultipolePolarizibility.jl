
"""
`module  JAC.MultipolePolarizibility`
... a submodel of JAC that contains all methods for computing static multipole polarizibilities.

    ONE NOTATION, WRITTEN `alpha^(Mp,K)`, AND NO SUBSCRIPTS -- the maintainer's instruction, 19-Oct-2026, in the form
    JAC already uses for `O^(Mp, emission)`, `T^(Mp, absorption)` and `H^(sms,C)`: a single parenthesised superscript
    carrying a comma-separated label list, multipole first and rank second. The literature splits these two labels
    between a superscript and a subscript and then disagrees about which is which, and the collision is not
    hypothetical: `alpha_2` means the RANK-2 DIPOLE polarizability in the Stark-shift literature and the SCALAR
    QUADRUPOLE polarizability in the dispersion-coefficient literature, which differ by orders of magnitude. Here they
    are `alpha^(E1,2)` and `alpha^(E2,0)` and cannot be confused.

    THE CONVENTION IS FIXED BY AN OBSERVABLE, NOT BY A CONVENIENT FACTOR (Rule 18). `alpha^(Mp,K)(J)` is defined as the
    coefficient of the angular shape A_K(J,M) in the polarizability of the individual magnetic substate,

        alpha(J,M)  =  sum_K  A_K(J,M) alpha^(Mp,K)(J),       A_K(J,M) = C^(JM)_(JM,K0) / C^(JJ)_(JJ,K0),   A_0 = 1

    so that alpha^(Mp,0) is the usual scalar polarizability and every alpha^(Mp,K>0) is the anisotropy carried by the
    stretched M = J substate. For K = 2 the shape is the familiar [3M^2 - J(J+1)]/[J(2J-1)], so alpha^(E1,0) and
    alpha^(E1,2) ARE the alpha_0 and alpha_2 of Angel & Sandars (1968) below, unchanged.

    THE GENERAL RANK, for any multipole and any rank the level can carry, is computed by
    MultipolePolarizibility.computeMultipolePolarizabilities from

        alpha^(Mp,K)(J)  =  sum_J'  c(L,K,J,J') |<J' || O^(Mp) || J>|^2 / (W_J' - W_J)

        c(L,K,J,J')      =  2 (-1)^(J+J') C^(JJ)_(JJ,K0) (2K+1)/sqrt(2J+1) * (L L K; 0 0 0) * {J J' L; L K J}

    with a Wigner 3-j and a Wigner 6-j symbol, and the three selection rules are carried BY the formula rather than
    imposed on it: the 3-j vanishes for odd K and for K > 2L, the Clebsch-Gordan for K > 2J. Setting K = 0 collapses it
    analytically to 2/((2J+1)(2L+1)), which at L = 1 is the Angel & Sandars scalar prefactor. Odd K -- the vector or
    axial polarizability, K = 1 -- exists only for a circularly polarized field and is not computed here.

    HOW MUCH OF THIS IS PROVEN, AND HOW. The angular machinery is exact and was settled numerically rather than by
    hand, since a factor of this shape is exactly what gets written with a wrong phase and compensated elsewhere:
      * the substate polarizability was built directly from Clebsch-Gordan coefficients and decomposed in the A_K
        shapes over 70 (J,L) combinations -- residual <= 1.7e-16, which proves the rank set complete;
      * the ratio of the extracted coefficient to (-1)^(J+J') {J J' L; L K J} is INDEPENDENT of J' to <= 9.1e-15,
        which proves the 6-j structure, and the normalisation and phase then matched with ZERO mismatches over 421
        samples;
      * with REAL matrix elements, alpha(J,M) rebuilt as sum_K A_K alpha^(E1,K) reproduces the direct substate sum to
        5.8e-16 on H 2p_3/2 with nd perturbers -- the first NON-ZERO alpha^(E1,2) this module has produced;
      * against the Dirac B-spline values of Tang, Zhang, Zhang, Jiang & Mitroy, arXiv:1205.6883, on H 1s and on a box
        MATCHED to the perturbers: alpha^(E1,0) 1.00000, alpha^(E2,0) 0.99999, alpha^(E3,0) 0.99999, alpha^(E4,0)
        1.0031 of the reference. That the K = 0 normalisation holds at FOUR multipole orders rules out a (2L+1)-type
        error, which would have shown as a factor of 5 at E2.
    `tools/probe-generalRankPolarizability.jl` is the harness and records its own results.

    AND WHAT IS NOT THE FORMULA BUT THE PERTURBER SET, which is where every remaining error lives. The discrete series
    stands in for the continuum, and it saturates far more slowly the higher the multipole: a few bound perturbers give
    0.81, 0.34, 0.08 and 0.01 of the reference at E1, E2, E3 and E4, and reaching 1.0 needs the series run to n ~ 16
    (E1, E2, E4) or n ~ 22 (E3). WORSE, THE RADIAL BOX IS AMPLIFIED BY THE MULTIPOLE ORDER, because the matrix elements
    weight r^L: the very box that is harmless for E1 (614 a.u. with the grid guard overridden, which is what
    examples/example-Cj.jl branch b does deliberately and correctly) gives alpha^(E3,0) 4.9 times and alpha^(E4,0)
    47500 TIMES the reference. Rule 12 with a multipole lever on it. Three guards warn about this -- see
    MultipolePolarizibility.checkPerturberSet, which also records the two candidate detectors that were measured and
    REJECTED for not discriminating.

    Stage 1 (implemented here): the static electric-dipole (E1) scalar and tensor polarizabilities
    alpha_0(J), alpha_2(J), following the classic Angel & Sandars (1968) rank-0/rank-2 decomposition
    (formula and prefactors confirmed against van Leeuwen & Hogervorst, Z. Phys. A 316, 149 (1984),
    Eq. (5)):

        alpha_0(J) = -2/(3(2J+1)) * sum_J' |<J||P||J'>|^2 / (W_J - W_J')

        alpha_2(J) = -2*sqrt(10/3)*sqrt(J(2J-1)/((2J+3)(J+1)(2J+1))) * (-1)^(2J)
                     * sum_J' (-1)^(J-J') {J J' 1; 1 2 J} |<J||P||J'>|^2 / (W_J - W_J')

    (alpha_2 is defined to be exactly 0 for J<1, where the tensor polarizability does not exist).

    The sum over J' runs, in principle, over the complete (opposite-parity) many-electron spectrum,
    including the continuum -- literally enumerating this is out of reach. Following the same
    approach already used by module-MultiPhotonTransition.jl for its own (structurally identical)
    sum-over-virtual-intermediate-states problem, the sum is instead carried out over a finite,
    user-supplied set of intermediate ASF -- a *local, approximate Green multiplet*
    (MultipolePolarizibility.Settings.gMultiplet) that brings in the physically relevant perturber
    levels (2p_1/2, 2p_3/2, ... for an ns ground level, say). This multiplet is generated
    INDEPENDENTLY, ahead of time, by the caller -- simply as a plain Atomic.Computation of energies
    & ASF (no properties) for the relevant opposite-parity configurations, e.g. 2p, 3p, 4p, ... (see
    examples/example-Cj.jl) -- and is simply consumed here; this module performs no SCF/CI generation
    of its own (mirroring MultiPhotonTransition.Settings.gMultiplet exactly). A high-n tail of such
    perturber configurations (e.g. up to 40p, 60p, ...) may also be included deliberately: these
    orbitals are no longer physically accurate bound Rydberg states at the default grid/basis size
    (see project notes), but remain properly normalized, mutually orthogonal ASF and so provide a
    simple, explicit (if crude) discrete approximation ("pseudo-continuum") to the true continuum's
    contribution -- their individual contributions can be tabulated and inspected exactly like any
    other perturber, without any change to the formula or code below.

    Every contribution (and the total alpha_0/alpha_2) is computed from the frequency-INDEPENDENT static
    multipole-moment operator (MultipoleMoment.emmStaticAmplitude, i.e. Johnson's r^k C_k), which is the
    operator a static polarizability -- and a C_6/C_8/C_10 dispersion coefficient -- actually requires. There
    is no gauge freedom in a static quantity, so the two components of the EmProperty are identical; they are
    retained only because Outcome carries EmProperty fields throughout JAC.

    Until 10-Aug-2026 the Coulomb component was instead taken from PhotoEmission.amplitude at omega = 0,
    i.e. from InteractionStrength.MabEmission -- the frequency-DEPENDENT transition operator -- and the two
    columns were presented side by side as a gauge-agreement diagnostic. That was a category error rather
    than a calibration one. The velocity and length forms are connected by <f|p|i> = i (E_f - E_i) <f|r|i>,
    so evaluating the velocity form at omega = 0 drops exactly the 1/omega that relates them; the result was
    alpha_0(Coulomb) ~ 1.5e-6 against 1.164 in Babushkin, six orders apart and completely insensitive to the
    radial box -- which is what identified it as the wrong operator rather than a numerical problem. There is
    no gauge diagnostic to be had from a static quantity, and none is offered here.

    Stage 2 IS DONE, 19-Oct-2026: general multipoles beyond E1, the general rank-K decomposition, and the DYNAMIC
    polarizability at both real and imaginary frequency. What remains deferred is the MAGNETIC multipoles, for which
    `MultipoleMoment.emmStaticAmplitude` has no operator and the literature is thin, and the dispersion coefficients
    themselves, which are a separate, later consumer of this module.

    THE FREQUENCY ENTERS ONLY THROUGH THE DENOMINATOR, so the dynamic case costs nothing beyond the static one. With
    dE = W_J' - W_J, MultipolePolarizibility.frequencyFactor replaces 1/dE by dE/(dE^2 - omega^2) at a real frequency
    and by dE/(dE^2 + omega^2) at the imaginary frequency i*omega, and nothing else changes: the same reduced matrix
    elements, the same angular coefficients, the same perturber set. `Settings.omegas` lists the REAL frequencies
    wanted (omega = 0 is always computed as well, so alpha0/alpha2 keep their meaning); the imaginary axis is reached
    by calling computeMultipolePolarizabilities directly, since i*omega has no physical meaning of its own and exists
    for dispersion coefficients.

    THE TWO AXES BEHAVE QUITE DIFFERENTLY, AND THE IMAGINARY ONE IS THE EASY ONE. On the real axis the factor diverges
    whenever omega meets an excitation energy, and ABOVE the lowest excitation energy the true alpha(omega) is COMPLEX
    -- the atom absorbs -- while a sum over a discrete perturber set can only return a principal-value-like real part
    with no width. On the imaginary axis dE^2 + omega^2 cannot vanish, so there are no poles at all; that is why a
    dispersion coefficient is built there. MultipolePolarizibility.checkFrequency warns about both real-axis traps.

    AND alpha(i*omega) MAKES A DISPERSION COEFFICIENT A DOUBLE SUM RATHER THAN A QUADRATURE. On a discrete perturber
    spectrum the Casimir-Polder integral is closed form, int_0^inf [a/(a^2+w^2)][b/(b^2+w^2)] dw = pi/(2(a+b)), so
    with A_a = c(L,0,J,J_a) |<J_a||O||J>|^2 one gets C_6 = (3/2) sum_ab A_a A_b / (dE_a + dE_b), which is Dalgarno's
    oscillator-strength form. Measured on H(1s)-H(1s) with np to n = 16, this reproduces the published C_6 -- see
    `tools/probe-dynamicPolarizability.jl`. The coefficients themselves belong in a module of their own.

    UNPARKED 13-Oct-2026 ON THE MAINTAINER'S INSTRUCTION, AND THE MODULE REPRODUCES EXACT HYDROGEN.
    It had been parked since 10-Aug-2026, and the whole of the reason was a stale number. The story is kept
    here because it cost two months and the mechanism is general.

      * **alpha_0(H 1s) = 4.499752 a.u. against the exact 4.500000** -- a relative error of 6e-05 -- from
        np perturbers n = 2...16 with the grid guard switched off (see below). The approach to it is what
        shows the pseudo-continuum is real rather than fortuitous: 3.4918 at 2p..4p, 3.6258 at 2p..8p --
        saturating towards the BOUND-series limit of about 3.663 -- then 4.2486 at 2p..12p as the high-n
        states cross threshold and begin to act as a continuum, 4.4997 at 2p..16p, and 4.4998 at 2p..20p,
        i.e. converged. **So the question the old block left open -- whether such a tail converges to the
        right number -- is answered, and it does.**

      * THE BOUND PART IS EXACT TERM BY TERM. A 2p+3p+4p sum gives 3.491829 against the closed-form 3.4921,
        and perturber by perturber 2.9595 / 0.4004 / 0.1320 against hydrogen's 2.9596 / 0.4005 / 0.1320, with
        each spin-orbit pair splitting 1 : 2 by its statistical weight. Both gauge columns are identical by
        construction, as a static quantity requires.

      * **THE 1.164 THAT GOT IT PARKED WAS THE RIGHT ANSWER DIVIDED BY THREE, AND THE DIVISION HAD ALREADY
        BEEN REPAIRED WHEN THE NOTE WAS WRITTEN.** Dividing today's alpha_0 by 3 reproduces 1.1639429, i.e.
        the recorded value to seven figures. The repair was `0bdff9f6` (the Mab/Mba consolidation) of
        09-Aug-2026 -- the day BEFORE the parking note -- and nobody re-ran the module afterwards. A working
        module therefore sat parked for two months behind a number that was obsolete the day it was recorded.
        `testModule_MultipolePolarizibility` now asserts the closed-form hydrogen values and has been shown to
        fail on exactly this error; the old assertion, 0 < alpha_0 < 4.5, passed 1.164 as happily as 3.4918.

    WHAT IT IS GOOD FOR NOW: static scalar and tensor polarizabilities as absolute values, not merely as lower
    bounds, PROVIDED the perturber set includes the high-n tail. What it still does not do is choose that tail
    for the caller -- see the grid note below.

    THE GRID, AND WHY THE GUARD MUST BE SWITCHED OFF FOR THE TAIL. Beyond about 5p on a 129 a.u. box,
    `Bsplines.checkGridRepresentation` refuses the perturbers, and it is RIGHT to: they are not accurate
    Rydberg states there. But that is exactly what a pseudo-continuum is -- box-confined, properly normalised,
    mutually orthogonal ASF standing in for the true continuum. So a caller who wants the converged value
    passes `AsfSettings(...; gridStopper = false)` when GENERATING THE PERTURBER MULTIPLET, deliberately and
    with this paragraph as the reason. The test of whether the basis is adequate is then not the guard but the
    CONVERGENCE OF THE SUM, which for hydrogen is checked against an exact answer above.

    THE UNITS, because comparing against a measured Stark shift without them has already produced one false alarm.
    `alpha_0` is returned and printed in ATOMIC UNITS. The Stark-shift literature quotes MHz/(kV/cm)^2, and the two
    differ by 1 a.u. = 1.244160e-04 MHz/(kV/cm)^2. Dividing the printed a.u. value by a published MHz/(kV/cm)^2 one
    yields a meaningless ratio near 2000-3000, which was recorded for eight days as "the module is too large by a
    factor of ~2600". CONVERTED, the same two Ba 6s8s cases give 1.63 against 5.68 and 1.88 against 4.68 -- a factor
    ~3 LOW, which is the lower bound this block describes, and consistent with H 1s at 0.46 of the exact value.
    So: convert before forming any ratio.

    **AND THE "NEVER ABOVE AN EXPERIMENT" THAT USED TO END THIS PARAGRAPH IS WITHDRAWN -- 14-Oct-2026, with a
    measured counterexample.** The lower-bound property is real but much narrower than that sentence claimed,
    and the two limits matter in practice:

      * it holds only for the LOWEST level of its parity, where every energy denominator carries the same sign
        so that each omitted perturber can only ADD. For an excited level such as Ba 6s8s the perturbers lying
        BELOW it contribute with the opposite sign, and the sum is then bounded in neither direction;
      * and even for a ground state it bounds only what THE GIVEN MEAN FIELD would yield with a complete
        perturber set -- it says nothing whatever about agreement with EXPERIMENT, because the field itself can
        overshoot. Measured on Li 2s against Puchalski et al.'s 164.0740 a.u., same perturber set throughout:
        DFS gives 183.7863 (12 % ABOVE), HS gives 227.9872 (39 % ABOVE), AL gives 167.6745 (2.2 % above).
        DFS underbinds the valence 2s to -0.079 Ha against the true -0.198, shrinking the 2s-2p denominator.

    **THAT IS ALSO THE USEFUL DIAGNOSTIC, so it is worth stating the other way round:** a ground-state value
    coming out ABOVE a trusted reference CANNOT be blamed on an incomplete perturber set, and points at the
    mean field instead. That is how the Li case was solved -- see `examples/example-Cj.jl` branch (C), which
    now runs on `ALField` for exactly this reason.
"""
module MultipolePolarizibility


using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..MultipoleMoment,
      ..Nuclear, ..PhotoEmission, ..Radial, ..TableStrings


"""
`MultipolePolarizibility.AU_ALPHA_IN_MHZ_PER_KVCM2`
    ... the factor that converts a static E1 polarizability from atomic units into MHz/(kV/cm)^2, the unit in which
        the Stark-shift literature quotes alpha_0 and alpha_2; a `factor::Float64` is returned. It is built as
        (1/2) * (E_h/h in Hz) * (1 kV/cm / E_au)^2 * 1e-6, with E_au = 5.142206751e9 V/cm the atomic unit of electric
        field strength, and the leading 1/2 is the customary half of DeltaE = -(1/2) alpha E^2 that experimenters
        absorb into the quoted number. Its value is 1.244160e-04, so a polarizability of 13098 a.u. is 1.63 in the
        literature's unit. Dividing a printed a.u. value by a published MHz/(kV/cm)^2 one -- rather than converting
        first -- yields a spurious ratio of some thousands, which is exactly how this module was once believed to be
        too large by a factor of 2600 when it is in fact a factor ~3 LOW; see the module docstring.
"""
const AU_ALPHA_IN_MHZ_PER_KVCM2 = 1.0e-6 * Defaults.convertUnits("energy: from atomic to Hz", 1.0) *
                                  0.5 * (1000.0 / 5.142206751e9)^2


"""
`struct  MultipolePolarizibility.Contribution`
    ... defines a type to keep the individual contribution of ONE intermediate (perturber) ASF from gMultiplet to the
        static multipole polarizability alpha^(Mp,K) of some other (fixed) level. One such term is built per perturber
        AND per multipole, and it carries the multipole it belongs to in its own `multipole` field, so that the
        contributions of several requested multipoles can be told apart in one Outcome.

    + nuLevel        ::Level        ... The intermediate (perturber) level nuLevel = J' from gMultiplet.
    + deltaE         ::Float64      ... Energy denominator W_J - W_J' (atomic units) entering this term.
    + contribAlpha0  ::EmProperty   ... This perturber's own additive contribution to alpha_0(J), in a.u.,
                                        in both the Coulomb and Babushkin gauge (see module docstring).
    + contribAlpha2  ::EmProperty   ... This perturber's own additive contribution to alpha_2(J), in a.u.,
                                        in both gauges; identically 0. for J<1, where the tensor
                                        polarizability does not exist.
    + multipole      ::EmMultipole  ... The multipole Mp whose polarizability this term contributes to; E1 for
                                        every term built by computeScalarTensorPolarizability.
    + reducedMe      ::Float64      ... The static reduced matrix element <J' || O^(Mp) || J> ITSELF, kept
                                        because the whole frequency dependence of a polarizability sits in the
                                        denominator: with reducedMe and deltaE stored, alpha^(Mp,K) at any real
                                        or imaginary frequency -- and hence any dispersion coefficient -- is a
                                        re-summation of these same terms and needs no new structure calculation.
    + contribAlphas  ::Array{Tuple{Int64,EmProperty},1}
                                    ... (K, contribution) for EVERY rank K this perturber feeds, in ascending K;
                                        contribAlpha0/contribAlpha2 are the K = 0 and K = 2 members of it.
    + omega          ::Float64      ... The frequency at which this term was evaluated; 0. for a static one. A
                                        Contribution belongs to ONE (multipole, frequency) evaluation, so a caller
                                        that asked for several frequencies can tell the terms apart by this field.
    + imaginaryOmega ::Bool         ... True if `omega` is to be read as the IMAGINARY frequency i*omega, which is
                                        what a dispersion coefficient needs and where the sum has no poles at all.
"""
struct Contribution
    nuLevel          ::Level
    deltaE           ::Float64
    contribAlpha0    ::EmProperty
    contribAlpha2    ::EmProperty
    multipole        ::EmMultipole
    reducedMe        ::Float64
    contribAlphas    ::Array{Tuple{Int64,EmProperty},1}
    omega            ::Float64
    imaginaryOmega   ::Bool
end


"""
`struct  MultipolePolarizibility.Outcome`  ... defines a type to keep the outcome of a static-polarizability computation.

    + level          ::Level                                 ... Atomic level to which the outcome refers to.
    + alpha0         ::EmProperty                             ... static scalar (rank 0) electric-dipole
                                                                    polarizability alpha_0(J), in atomic units,
                                                                    in both the Coulomb and Babushkin gauge.
    + alpha2         ::EmProperty                             ... static tensor (rank 2) electric-dipole
                                                                    polarizability alpha_2(J), in atomic units,
                                                                    in both gauges; identically 0. for J<1,
                                                                    where the tensor polarizability does not exist.
    + contributions  ::Array{MultipolePolarizibility.Contribution,1}
                                                                ... individual, additive contribution of every
                                                                    intermediate ASF in gMultiplet that the requested
                                                                    multipoles can reach from level; each term names
                                                                    its own multipole (see module docstring).
    + alphas         ::Array{Tuple{EmMultipole,Int64,Float64,EmProperty},1}
                                                                ... (Mp, K, omega, alpha^(Mp,K)(omega)) for every
                                                                    multipole requested, every rank K this level can
                                                                    carry and every frequency in Settings.omegas,
                                                                    in the order computed; omega = 0. is ALWAYS
                                                                    present, and alpha0/alpha2 above are its (E1,0)
                                                                    and (E1,2) members. These frequencies are REAL;
                                                                    the imaginary axis is reached by calling
                                                                    computeMultipolePolarizabilities directly, since
                                                                    i*omega is an internal device for dispersion
                                                                    coefficients rather than a property anyone asks
                                                                    for on its own.
"""
struct Outcome
    level          ::Level
    alpha0         ::EmProperty
    alpha2         ::EmProperty
    contributions  ::Array{MultipolePolarizibility.Contribution,1}
    alphas         ::Array{Tuple{EmMultipole,Int64,Float64,EmProperty},1}
end


"""
`MultipolePolarizibility.Outcome()`
    ... constructor for an `empty` instance of MultipolePolarizibility.Outcome for the computation of static
        polarizibilities.
"""
function Outcome()
    Outcome(Level(), EmProperty(0.), EmProperty(0.), MultipolePolarizibility.Contribution[],
            Tuple{EmMultipole,Int64,Float64,EmProperty}[])
end


# `Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)`
#		 ... prepares a proper printout of the variable outcome::MultipolePolarizibility.Outcome.
function Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)
    println(io, "level:                     $(outcome.level)  ")
    println(io, "alpha0:                    $(outcome.alpha0)  ")
    println(io, "alpha2:                    $(outcome.alpha2)  ")
    println(io, "contributions:             $(length(outcome.contributions)) perturber level(s)  ")
end


"""
`struct  MultipolePolarizibility.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of computing multipolar polarizibilities.

    + multipoles                ::Array{EmMultipole,1}  ... Multipoles to be considered for the polarity. Stage 1
                                                              only supports E1 (the classic scalar/tensor
                                                              polarizability); any other multipole present in the
                                                              list is skipped with a warning.
    + gMultiplet                ::Multiplet             ... A local, approximate Green multiplet of intermediate
                                                              (opposite-parity) ASF that the sum-over-states runs
                                                              over -- generated INDEPENDENTLY, ahead of time, by the
                                                              caller (e.g. via AtomicState.GreenExpansion/generate(),
                                                              see examples/example-Aj.jl) and simply consumed here;
                                                              mirrors MultiPhotonTransition.Settings.gMultiplet.
    + omegas                    ::Array{Float64,1}      ... List of REAL frequencies (in atomic units) at which the
                                                              dynamic polarizability is wanted, beside the static
                                                              case. omega = 0. is always computed whether or not it
                                                              appears here, so alpha0/alpha2 of the Outcome keep
                                                              their meaning; duplicates are dropped. A frequency at
                                                              or above the lowest excitation energy is WARNED about
                                                              rather than refused -- see
                                                              MultipolePolarizibility.checkFrequency.
    + printBefore               ::Bool                  ... True if a list of selected levels is printed before the
                                                            actual computations start.
    + levelSelection            ::LevelSelection        ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    multipoles                  ::Array{EmMultipole,1}
    gMultiplet                  ::Multiplet
    omegas                      ::Array{Float64,1}
    printBefore                 ::Bool
    levelSelection              ::LevelSelection
end


"""
`MultipolePolarizibility.Settings()`
    ... constructor for an `empty` instance of MultipolePolarizibility.Settings for the computation dynamic polarizibilities.
"""
function Settings()
    Settings(EmMultipole[], Multiplet(), Float64[], false, LevelSelection() )
end


"""
`MultipolePolarizibility.Settings(set::MultipolePolarizibility.Settings;`

        multipoles=.., gMultiplet=.., omegas=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::MultipolePolarizibility.Settings.
"""
function Settings(set::MultipolePolarizibility.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
        gMultiplet::Union{Nothing,Multiplet}=nothing,
        omegas::Union{Nothing,Array{Float64,1}}=nothing,   printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(multipoles)       multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(gMultiplet)       gMultipletx     = set.gMultiplet     else   gMultipletx     = gMultiplet     end
    if  isnothing(omegas)           omegasx         = set.omegas         else   omegasx         = omegas         end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( multipolesx, gMultipletx, omegasx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings:MultipolePolarizibility.Settings)`
#		... prepares a proper printout of the variable settings::MultipolePolarizibility.Settings.
function Base.show(io::IO, settings::MultipolePolarizibility.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
    println(io, "omegas:                   $(settings.omegas)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`MultipolePolarizibility.allowedRanks(multipole::EmMultipole, J::AngularJ64)`  
    ... returns the ranks K for which a level of angular momentum J can carry a non-vanishing static multipole
        polarizability alpha^(Mp,K) in the field of the given 2^L-pole `multipole`. Three rules act together: K <= 2L,
        because the rank-K component is built by coupling the two 2^L-pole operators of the second-order expression;
        K <= 2J, by the Wigner-Eckart theorem applied between the substates of one level; and K EVEN, because the two
        coupled operators are identical. Odd K does exist physically -- K = 1 is the vector (axial) polarizability --
        but only for a circularly polarized field, and it is deliberately not returned here. A ranks::Array{Int64,1}
        is returned, in ascending order; it always contains K = 0.
"""
function allowedRanks(multipole::EmMultipole, J::AngularJ64)
    Jd = Basics.twice(J) / 2.0

    return( [K   for K = 0:min(2*multipole.L, round(Int64, 2*Jd))   if iseven(K)] )
end


"""
`MultipolePolarizibility.checkFrequency(multipole::EmMultipole, level::Level, omega::Float64, imaginaryOmega::Bool,
                                        contributions::Array{MultipolePolarizibility.Contribution,1})`  
    ... inspects the frequency a dynamic polarizability was asked for and WARNS LOUDLY, to stdout, about the two ways
        it goes wrong while still returning a finite, plausible-looking number. Both concern the REAL axis only; on the
        imaginary axis the denominator is dE^2 + omega^2, which cannot vanish, and that is precisely why a dispersion
        coefficient is evaluated there. Nothing is returned.

        THE FIRST IS THAT THE TRUE POLARIZABILITY IS COMPLEX AND THIS SUM IS REAL. Above the lowest excitation energy
        the atom absorbs, alpha(omega) acquires an imaginary part, and a sum over a DISCRETE perturber set cannot
        represent it -- what comes back is a principal-value-like real part with no width. The number is not a
        polarizability there, however reasonable it looks.

        THE SECOND IS A POLE. As omega approaches a perturber's excitation energy the factor dE/(dE^2 - omega^2)
        diverges, so the result is carried by ONE term and is only as accurate as that ONE energy. For a genuine
        low-lying perturber that is a real resonance; for a high-n pseudo-state standing in for the continuum -- which
        is most of the set, and whose energy means nothing individually -- it is an artefact of the discretisation.
"""
function checkFrequency(multipole::EmMultipole, level::Level, omega::Float64, imaginaryOmega::Bool,
                        contributions::Array{MultipolePolarizibility.Contribution,1})
    if  omega == 0.0  ||  imaginaryOmega  ||  isempty(contributions)    return( nothing )    end

    # dE = W_J' - W_J > 0 for every perturber above the level;  only those can be resonant with a positive omega
    dEs = [ -c.deltaE  for c in contributions  if -c.deltaE > 0.0 ]
    if  isempty(dEs)    return( nothing )    end
    dEmin = minimum(dEs)

    if  omega >= dEmin
        printstyled("\n>>> omega = " * @sprintf("%.6f", omega) * " a.u. IS AT OR ABOVE THE LOWEST EXCITATION ENERGY " *
                    @sprintf("%.6f", dEmin) * " a.u. of this perturber set,\n" *
                    "    so the true alpha^($multipole, K) is COMPLEX there and this real sum over a DISCRETE set " *
                    "cannot represent its\n    imaginary part. What is returned is a principal-value-like real part " *
                    "with no width -- not a polarizability.\n", color=:light_red)
    end

    # the nearest perturber in frequency, as a fraction of its own energy: that is what sets the sensitivity
    wa    = argmin( [abs(dE - omega) for dE in dEs] )
    dNear = dEs[wa]
    if  abs(dNear - omega) < 2.0e-2 * dNear
        printstyled("\n>>> omega = " * @sprintf("%.6f", omega) * " a.u. SITS ON A POLE: a perturber lies at " *
                    @sprintf("%.6f", dNear) * " a.u., i.e. " *
                    @sprintf("%.2f", 100*abs(dNear-omega)/dNear) * " % away, so\n" *
                    "    alpha^($multipole, K) is carried by that ONE term and is only as accurate as that ONE " *
                    "energy. For a genuine\n    low-lying perturber that is a resonance; for a high-n pseudo-state " *
                    "standing in for the continuum it is an\n    artefact of the discretisation, since such a " *
                    "state's individual energy means nothing.\n", color=:light_red)
    end

    return( nothing )
end


"""
`MultipolePolarizibility.frequencyFactor(dE::Float64, omega::Float64, imaginaryOmega::Bool)`  
    ... returns the energy-denominator factor which carries the WHOLE frequency dependence of a multipole
        polarizability, so that alpha^(Mp,K)(omega) = sum_J' c(L,K,J,J') |<J'||O^(Mp)||J>|^2 * this factor, with
        dE = W_J' - W_J:

            omega = 0           1 / dE
            real omega          dE / (dE^2 - omega^2)
            imaginary i*omega   dE / (dE^2 + omega^2)

        All three agree at omega = 0. THE IMAGINARY AXIS HAS NO POLES -- dE^2 + omega^2 cannot vanish -- which is why
        a dispersion coefficient is built there and why alpha(i*omega) is numerically the EASIER quantity, despite
        having no direct physical meaning of its own. On the real axis the factor diverges whenever omega meets an
        excitation energy; MultipolePolarizibility.checkFrequency warns about that. A value::Float64 is returned.

        For an EXCITED level of interest some perturbers lie BELOW it, dE < 0, and the same expressions hold: such a
        term contributes with the opposite sign, which is the de-excitation channel and not an error.
"""
function frequencyFactor(dE::Float64, omega::Float64, imaginaryOmega::Bool)
    if      omega == 0.0        return( 1.0 / dE )
    elseif  imaginaryOmega      return( dE / (dE*dE + omega*omega) )
    else                        return( dE / (dE*dE - omega*omega) )
    end
end


"""
`MultipolePolarizibility.checkPerturberSet(multipole::EmMultipole, level::Level, gMultiplet::Multiplet,
                                            grid::Radial.Grid,
                                            contributions::Array{MultipolePolarizibility.Contribution,1},
                                            alphaScalar::Float64)`  
    ... inspects the perturber set that was actually used and WARNS LOUDLY, to stdout, about the two conditions under
        which a multipole polarizability comes out wrong while looking perfectly reasonable. The first is that NO level
        of `gMultiplet` can be reached at all: a 2^L-pole operator needs perturbers of its own parity and within its own
        triangle, so a set of np levels -- right for E1 -- reaches NOTHING for E2, and the sum then returns a clean zero
        rather than an error. The second is that the sum HAS NOT SATURATED: the discrete perturber series stands in for
        the continuum, so the highest-lying member must already contribute negligibly, and where it does not the result
        is a lower bound of unknown quality rather than a polarizability. The third, for L >= 2 only, re-runs
        `Bsplines.checkOrbitalConsistency` on the perturber orbitals, because an ill-matched radial box is far more
        dangerous for a high multipole than for E1 -- see the measured figures in the module docstring. Nothing is
        returned.

        NO FOURTH GUARD IS OFFERED, AND THAT IS DELIBERATE. Two candidate detectors for a box-corrupted sum were
        built and MEASURED on one healthy and two catastrophic E4 cases, and NEITHER discriminates: where the largest
        contribution sits in the energy-ordered series put the healthy case FURTHER up (78 % of the way) than either
        corrupted one (43 %, 37 %), and box/<r> of the dominant perturber ran 1.67/1.45 for the healthy cases against
        1.33/1.27 for the corrupted ones -- overlapping, and in the wrong direction. An uncalibrated guard here would
        be worse than none, since it would be read as a clean bill of health.
"""
function checkPerturberSet(multipole::EmMultipole, level::Level, gMultiplet::Multiplet, grid::Radial.Grid,
                           contributions::Array{MultipolePolarizibility.Contribution,1}, alphaScalar::Float64)
    # the parity a perturber must carry: an electric 2^L-pole has parity (-1)^L, a magnetic one (-1)^(L+1)
    nReach = 0
    for  nuLevel in gMultiplet.levels
        if  MultipolePolarizibility.isReachable(multipole, level, nuLevel)   nReach = nReach + 1    end
    end

    if  nReach == 0
        wantsSame = multipole.electric ? iseven(multipole.L) : isodd(multipole.L)
        printstyled("\n>>> NO PERTURBER REACHED for alpha^($multipole, K) of level $(level.index): " *
                    "the given gMultiplet holds $(length(gMultiplet.levels)) level(s), and not one of them is\n" *
                    "    reachable from J^P = $(level.J)$(string(level.parity)) by a $multipole operator. " *
                    "A $multipole needs perturbers of the " * (wantsSame ? "SAME" : "OPPOSITE") * " parity with\n" *
                    "    |J - $(multipole.L)| <= J' <= J + $(multipole.L), so the polarizability below is a " *
                    "STRUCTURAL ZERO and not a physical one.\n", color=:light_red)
        return( nothing )
    end

    if  isempty(contributions)  ||  alphaScalar == 0.0   return( nothing )   end

    # the highest-lying perturber stands at the top of the discrete stand-in for the continuum; its own share of
    # alpha^(Mp,0) is the cheapest honest measure of whether the series has run out
    wa       = sortperm( [c.nuLevel.energy  for c in contributions] )
    topContr = contributions[wa[end]]
    topShare = 0.0
    for  (K, cAlpha) in topContr.contribAlphas
        if  K == 0   topShare = abs( cAlpha.Babushkin / alphaScalar )    end
    end

    if  topShare > 1.0e-2
        printstyled("\n>>> PERTURBER SUM NOT SATURATED for alpha^($multipole, 0) of level $(level.index): " *
                    "the HIGHEST-lying perturber still carries\n" *
                    @sprintf("    %.2f %% of the total (%d perturbers, top one at %.3f a.u. above the level). ",
                             100*topShare, length(contributions), topContr.nuLevel.energy - level.energy) *
                    "The discrete series stands in for\n" *
                    "    the continuum, so this is a LOWER BOUND of unknown quality rather than a converged " *
                    "polarizability; extend gMultiplet\n    to higher n and watch alpha^($multipole, 0) settle.\n",
                    color=:light_red)
    end

    ## For a rank-L polarizability the matrix elements weight r^L, so the radial box matters far more than it does
    ## for E1 -- and overriding the grid guard, which is standard and right for an E1 pseudo-continuum (see
    ## examples/example-Cj.jl branch b), is ruinous above it. Rather than invent a detector, JAC's own calibrated
    ## check is re-run here on the perturber orbitals.
    if  multipole.L >= 2
        orbitals = Dict{Subshell,Orbital}()
        for  nuLevel in gMultiplet.levels
            for  (sub, orb) in nuLevel.basis.orbitals    orbitals[sub] = orb    end
        end
        if  !isempty(orbitals)
            wa = Bsplines.checkOrbitalConsistency(orbitals, grid; stopper=false)
            if  !wa.ok
                printstyled(">>> AND THE ORBITAL CHECK ABOVE MATTERS MORE FOR alpha^($multipole, K) THAN IT WOULD FOR E1: " *
                            "the matrix elements\n    weight r^$(multipole.L), so an ill-matched box is amplified by the " *
                            "multipole order. Measured on H 1s, 19-Oct-2026, against\n    the Dirac values of " *
                            "arXiv:1205.6883: on a MATCHED 129 a.u. box alpha^(E3,0) and alpha^(E4,0) come out at 1.0000 " *
                            "and 1.0031\n    of the reference, while on a 614 a.u. box they come out at 4.9 and 47500 " *
                            "times it. Match the box.\n", color=:light_red)
            end
        end
    end

    return( nothing )
end


"""
`MultipolePolarizibility.computeMultipolePolarizabilities(multipole::EmMultipole, level::Level, gMultiplet::Multiplet,
                                                           grid::Radial.Grid; omega::Float64=0.0,
                                                           imaginaryOmega::Bool=false)`  
    ... computes the static multipole polarizabilities alpha^(Mp,K)(level) of `level` for the 2^L-pole `multipole` and
        for EVERY rank K that the level can carry (MultipolePolarizibility.allowedRanks), by summing the static reduced
        matrix elements <J' || O^(Mp) || J> = MultipoleMoment.emmStaticAmplitude(L, ...) over every level of
        `gMultiplet`. The convention is the one fixed in the module docstring: alpha^(Mp,K) is the coefficient of the
        angular shape A_K(J,M) = C^(JM)_(JM,K0) / C^(JJ)_(JJ,K0) in the substate polarizability alpha(J,M), so that
        alpha^(Mp,0) is the scalar polarizability and alpha^(Mp,K>0) is the anisotropy carried by the M = J substate.
        A static polarizability has no gauge freedom, so both components of every EmProperty hold the same number. The
        perturber set is checked and warned about (MultipolePolarizibility.checkPerturberSet). A tuple
        (alphas::Array{Tuple{Int64,EmProperty},1}, contributions::Array{MultipolePolarizibility.Contribution,1}) is
        returned, the latter giving each perturber's own additive share of every rank.

        THE FREQUENCY ENTERS ONLY THROUGH THE DENOMINATOR, which is why the dynamic case costs nothing beyond the
        static one: `omega` (real, atomic units) or, with `imaginaryOmega = true`, the imaginary frequency i*omega
        replaces 1/dE by MultipolePolarizibility.frequencyFactor and NOTHING else changes -- the same reduced matrix
        elements, the same angular coefficients, the same perturber set. A caller sweeping many frequencies therefore
        pays for one structure calculation, and `Contribution.reducedMe` is stored so that even that can be reused.
        The frequency is checked and warned about by MultipolePolarizibility.checkFrequency, which matters on the real
        axis only.
"""
function computeMultipolePolarizabilities(multipole::EmMultipole, level::Level, gMultiplet::Multiplet,
                                          grid::Radial.Grid; omega::Float64=0.0, imaginaryOmega::Bool=false)
    Ks            = MultipolePolarizibility.allowedRanks(multipole, level.J)
    alphaK        = zeros( Float64, length(Ks) )
    contributions = MultipolePolarizibility.Contribution[]

    for  nuLevel in gMultiplet.levels
        deltaE = level.energy - nuLevel.energy
        if  deltaE == 0.0   continue    end

        # Both emmStaticAmplitude's and PhotoEmission.amplitude's SpinAngular/CSF machinery need BOTH levels' CSFs
        # indexed against the SAME subshell list; nuLevel (from gMultiplet's typically much larger pseudo-basis) and
        # level (the level of interest's own, typically much smaller basis) do not share one -- reduce both onto their
        # union first (same idiom as AutoIonization/InternalConversion's amplitude construction).
        subshellList = Basics.generate(OrderedSubshellList(), nuLevel.basis, level.basis)
        newNuLevel   = Basics.generateLevelWithSymmetryReducedBasis(nuLevel, subshellList)
        newLevel     = Basics.generateLevelWithSymmetryReducedBasis(level,   subshellList)
        dStatic      = MultipoleMoment.emmStaticAmplitude(multipole.L, newNuLevel, newLevel, grid)
        if  abs2(dStatic) == 0.0   continue    end

        ## The coefficient below multiplies |<J'||O||J>|^2 / (E_J' - E_J), i.e. the denominator with the SIGN OPPOSITE
        ## to deltaE, which is kept as E_J - E_J' because that is what Contribution has always stored and printed.
        ## THE STATIC BRANCH IS WRITTEN EXACTLY AS IT WAS, not folded into the general one, so that a static result is
        ## BIT-identical and not merely equal to rounding: c*|d|^2*(1/dE) and -c*|d|^2/dE differ in the last bit, and
        ## test/approved is not to move for a change that adds a frequency.
        contribAlphas = Tuple{Int64,EmProperty}[]
        for  (i, K) in enumerate(Ks)
            coeff = MultipolePolarizibility.rankCoefficient(multipole.L, K, level.J, nuLevel.J)
            if  omega == 0.0
                wa = - coeff * abs2(dStatic) / deltaE
            else
                wa =   coeff * abs2(dStatic) *
                       MultipolePolarizibility.frequencyFactor(-deltaE, omega, imaginaryOmega)
            end
            alphaK[i] = alphaK[i] + wa
            push!( contribAlphas, (K, EmProperty(wa, wa)) )
        end

        ## alpha0/alpha2 are retained on Contribution so that every existing caller and printout is unaffected; they
        ## are the K = 0 and K = 2 members of contribAlphas and are not computed a second time.
        c0 = EmProperty(0.);   c2 = EmProperty(0.)
        for  (K, cAlpha) in contribAlphas
            if      K == 0   c0 = cAlpha
            elseif  K == 2   c2 = cAlpha
            end
        end
        push!( contributions, MultipolePolarizibility.Contribution(nuLevel, deltaE, c0, c2, multipole, dStatic,
                                  contribAlphas, omega, imaginaryOmega) )
    end

    alphas = Tuple{Int64,EmProperty}[]
    for  (i, K) in enumerate(Ks)   push!( alphas, (K, EmProperty(alphaK[i], alphaK[i])) )    end
    ## THE PERTURBER SET DOES NOT DEPEND ON THE FREQUENCY, so it is checked at omega = 0 only. A caller sweeping
    ## thirty frequencies would otherwise get thirty copies of the same warning, and a guard that repeats is a guard
    ## that gets scrolled past.
    scalar = isempty(alphas) ? 0.0 : alphas[1][2].Babushkin
    if  omega == 0.0
        MultipolePolarizibility.checkPerturberSet(multipole, level, gMultiplet, grid, contributions, scalar)
    end
    MultipolePolarizibility.checkFrequency(multipole, level, omega, imaginaryOmega, contributions)

    return( (alphas, contributions) )
end


"""
`MultipolePolarizibility.computeScalarTensorPolarizability(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)`
    ... computes the static E1 scalar and tensor polarizabilities of `level`, i.e. alpha^(E1,0) and alpha^(E1,2) in the
        notation of the module docstring, which the Stark-shift literature writes alpha_0 and alpha_2 (Angel & Sandars
        1968). This is the E1 rank-0/rank-2 FACE of MultipolePolarizibility.computeMultipolePolarizabilities and adds
        nothing of its own: it is kept because it is the entry point that module-StarkShift.jl and several examples
        call, and because those two ranks are what a Stark shift of a linearly polarized field needs. A static
        polarizability carries no gauge freedom, so both components of each EmProperty hold the same number. A tuple
        (alpha0::EmProperty, alpha2::EmProperty, contributions::Array{MultipolePolarizibility.Contribution,1}) is
        returned, the latter giving each perturber level's own additive share of alpha0/alpha2 in both gauges.

        UNTIL 19-Oct-2026 THIS FUNCTION CARRIED ITS OWN LOOP, with the two Angel & Sandars prefactors written out. The
        general function reproduces that loop EXACTLY -- 0.0e+00 relative difference on the total and 3.2e-16 on the
        worst of its eight perturber terms, on the H 1s case of test/approved -- which is what made it safe to have
        one implementation instead of two that can drift apart.
"""
function computeScalarTensorPolarizability(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)
    alphas, contributions = MultipolePolarizibility.computeMultipolePolarizabilities(E1, level, gMultiplet, grid)
    alpha0 = EmProperty(0.);     alpha2 = EmProperty(0.)
    for  (K, alphaK) in alphas
        if      K == 0   alpha0 = alphaK
        elseif  K == 2   alpha2 = alphaK
        end
    end

    return  (alpha0, alpha2, contributions)
end


"""
`MultipolePolarizibility.computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, gMultiplet::Multiplet,
                                                        grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
    ... to compute the static scalar/tensor E1 polarizabilities for a given level; an outcome::MultipolePolarizibility.Outcome is
        returned for which alpha0/alpha2/contributions are now evaluated explicitly. If E1 is not among
        settings.multipoles (or gMultiplet is empty), alpha0=alpha2=0. is returned without further computation.
"""
function  computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, gMultiplet::Multiplet,
                                       grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
    if  isempty(settings.multipoles)  ||  isempty(gMultiplet.levels)   return  outcome   end

    alphas = Tuple{EmMultipole,Int64,Float64,EmProperty}[];   contributions = MultipolePolarizibility.Contribution[]
    alpha0 = EmProperty(0.);                                  alpha2        = EmProperty(0.)

    ## omega = 0. IS ALWAYS COMPUTED, whether or not the caller listed it, because alpha0/alpha2 of the Outcome are
    ## defined as the STATIC scalar and tensor polarizabilities and several callers read them. It is also computed
    ## FIRST, so that with an empty settings.omegas -- the default -- the loop below is term for term what it was
    ## before a frequency existed here.
    frequencies = unique( vcat(0.0, settings.omegas) )

    for  multipole in settings.multipoles,  omega in frequencies
        mpAlphas, mpContributions = MultipolePolarizibility.computeMultipolePolarizabilities(multipole, outcome.level,
                                                                 gMultiplet, grid; omega=omega)
        for  (K, alphaK) in mpAlphas             push!( alphas, (multipole, K, omega, alphaK) )    end
        for  contribution in mpContributions     push!( contributions, contribution )               end

        ## alpha0/alpha2 remain the STATIC E1 rank-0 and rank-2 members, so that every existing caller, printout and
        ## approved reference sees exactly what it saw before.
        if  multipole == E1  &&  omega == 0.0
            for  (K, alphaK) in mpAlphas
                if      K == 0   alpha0 = alphaK
                elseif  K == 2   alpha2 = alphaK
                end
            end
        end
    end

    return  MultipolePolarizibility.Outcome(outcome.level, alpha0, alpha2, contributions, alphas)
end


"""
`MultipolePolarizibility.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                                settings::MultipolePolarizibility.Settings; output=true)`
    ... to compute (as selected) the static E1 scalar/tensor polarizabilities for the levels of the given multiplet and
        as specified by the given settings. The intermediate (perturber) ASF are taken directly from
        settings.gMultiplet, which the caller must generate independently, ahead of time (see the module docstring
        and examples/example-Aj.jl/example-Cj.jl). The results are printed in neat tables to screen but nothing is
        returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
            settings::MultipolePolarizibility.Settings; output=true)
    println("")
    printstyled("MultipolePolarizibility.computeOutcomes(): The computation of multipole polarizibilities starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    outcomes = MultipolePolarizibility.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    MultipolePolarizibility.displayOutcomes(outcomes)    end
    #
    # Calculate all amplitudes and requested properties
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    newOutcomes = MultipolePolarizibility.Outcome[]
    for  outcome in outcomes
        newOutcome = MultipolePolarizibility.computeAmplitudesProperties(outcome, settings.gMultiplet, grid, settings)
        push!( newOutcomes, newOutcome)
        MultipolePolarizibility.displayContributions(stdout, newOutcome)
        if  printSummary    MultipolePolarizibility.displayContributions(iostream, newOutcome)   end
    end
    # Print all results to screen
    MultipolePolarizibility.displayResults(stdout, newOutcomes)
    if  printSummary    MultipolePolarizibility.displayResults(iostream, newOutcomes)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`MultipolePolarizibility.determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings)`
    ... to determine a list of Outcomes's for the computation of the multipole polarizibilities for the given multiplet. It takes
        into account the particular selections and settings. An Array{MultipolePolarizibility.Outcome,1} is returned. Apart from the
        level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings)
    outcomes = MultipolePolarizibility.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, MultipolePolarizibility.Outcome(level, EmProperty(0.), EmProperty(0.),
                                 MultipolePolarizibility.Contribution[], Tuple{EmMultipole,Int64,Float64,EmProperty}[]) )
        end
    end
    return( outcomes )
end


"""
`MultipolePolarizibility.displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})`
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected MultipolePolarizibility levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    #
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`MultipolePolarizibility.displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})`
    ... to display the energies and static scalar/tensor E1 polarizabilities for the selected levels. A neat table is
        printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})
    nx = 141
    println(stream, " ")
    println(stream, "  Static electric-dipole scalar/tensor polarizibilities (Coulomb and Babushkin gauge, not calibrated " *
                     "onto a common scale -- see module docstring):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(15, "alpha_0 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_2 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_2 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[MHz/(kV/cm)^2]"; na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))    * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Coulomb) )                  * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin) )                * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Coulomb) )                  * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Babushkin) )                * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin * MultipolePolarizibility.AU_ALPHA_IN_MHZ_PER_KVCM2) ) * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    # THE DOMINANT ENERGY DENOMINATOR, PRINTED BECAUSE IT IS THE ONE NUMBER THAT DIAGNOSES THE WHOLE FAMILY OF
    # FAILURES HERE (item 70, 14-Oct-2026).  alpha_0 ~ |<a||P||b>|^2 / (W_b - W_a), so it is INVERSELY
    # proportional to this gap, and a mean field that misplaces the lowest perturber misplaces alpha_0 by the
    # same factor -- while every internal check in this module still passes, because the formula is right and
    # only its input is wrong.  A reader who can compare ONE number against a measured level energy can catch
    # that in seconds;  a reader given only alpha_0 cannot catch it at all.
    #
    # MEASURED, Li [He]2s against Puchalski et al.'s 164.0740 a.u. with ONE perturber set throughout:
    #     DFSField  183.7863 (+12 %, gap 0.056168 Ha)    HSField  227.9872 (+39 %)    ALField  167.6745 (+2.2 %,
    #     gap 0.067688 Ha against the experimental 0.067700).  DFS puts the Li 2s at -0.079 Ha where the true
    #     value is -0.198, i.e. underbound by a factor 2.5, and the gap and alpha_0 move together exactly as
    #     the formula says they must.  `AsfSettings()` DEFAULTS to DFSField, so this is the common case, not
    #     an exotic one.
    for  outcome in outcomes
        # the STATIC terms only: with several frequencies requested the list holds one set per frequency, and the
        # note below is about alpha_0, which is the static quantity
        statics = [c  for c in outcome.contributions  if c.omega == 0.0  &&  c.multipole == E1]
        isempty(statics)   &&   continue
        dominant = statics[1]
        for  c in statics
            abs(c.contribAlpha0.Babushkin) > abs(dominant.contribAlpha0.Babushkin)   &&   (dominant = c)
        end
        share = abs(outcome.alpha0.Babushkin) > 0.0 ?
                abs(dominant.contribAlpha0.Babushkin / outcome.alpha0.Babushkin) : 0.0
        println(stream, "  Level $(outcome.level.index):  the DOMINANT perturber is " *
                        "$(LevelSymmetry(dominant.nuLevel.J, dominant.nuLevel.parity)) at " *
                        @sprintf("Delta-E = %.6f", abs(dominant.deltaE)) * " a.u. " *
                        "(" * @sprintf("%.6f", Defaults.convertUnits("energy: from atomic",
                                                                      abs(dominant.deltaE))) * " " *
                        TableStrings.inUnits("energy") * "), carrying " *
                        @sprintf("%.1f", 100*share) * " % of alpha_0.")
    end
    println(stream, "  ")
    println(stream, "  >>> alpha_0 is INVERSELY proportional to that gap, so CHECK IT against a measured level")
    println(stream, "  >>> energy before trusting alpha_0: a mean field that underbinds the valence electron")
    println(stream, "  >>> shrinks the gap and inflates alpha_0 by the same factor, with every internal check")
    println(stream, "  >>> in this module still passing.  Measured on Li [He]2s against a literature 164.0740:")
    println(stream, "  >>>     DFSField  +12 %      HSField  +39 %      ALField  +2.2 %")
    println(stream, "  >>> and `AsfSettings()` DEFAULTS to DFSField.  Prefer ALField or EOLField for a")
    println(stream, "  >>> polarizability;  see examples/example-Cj.jl branch (C).")
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    # THE GENERAL-RANK TABLE IS PRINTED ONLY WHEN THERE IS SOMETHING IN IT BEYOND (E1,0) AND (E1,2), i.e. only when
    # a multipole other than E1 or a rank other than 0 and 2 was actually asked for. The table above then stays the
    # whole of the output for every caller that wants a Stark shift, and no stored reference moves.
    general = false
    for  outcome in outcomes
        for  (mp, K, omega, alphaK) in outcome.alphas
            if  mp != E1  ||  !(K in [0, 2])  ||  omega != 0.0    general = true    end
        end
    end
    if  general
        println(stream, " ")
        println(stream, "  MULTIPOLE polarizibilities alpha^(Mp,K)(omega), in the convention of the module docstring: " *
                         "alpha^(Mp,K) is the")
        println(stream, "  coefficient of the angular shape A_K(J,M) in the substate polarizability, so alpha^(Mp,0) is " *
                         "the scalar one and")
        println(stream, "  alpha^(Mp,K>0) is the anisotropy carried by the stretched M = J substate. The frequencies " *
                         "are REAL, and omega = 0 is")
        println(stream, "  the static case:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                         sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                         sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center( 8, "Mp";    na=2);                         sb = sb * TableStrings.hBlank(10)
        sa = sa * TableStrings.center( 5, "K";     na=3);                         sb = sb * TableStrings.hBlank( 8)
        sa = sa * TableStrings.center(14, "omega";  na=3)
        sb = sb * TableStrings.center(14, "[a.u.]"; na=3)
        sa = sa * TableStrings.center(18, "alpha^(Mp,K)"; na=3)
        sb = sb * TableStrings.center(18, "[a.u.]"; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        for  outcome in outcomes
            sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            for  (mp, K, omega, alphaK) in outcome.alphas
                sa = "  "
                sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * TableStrings.center( 8, string(mp);  na=2)
                sa = sa * TableStrings.center( 5, string(K);   na=3)
                sa = sa * TableStrings.flushright(14, @sprintf("%.6f", omega) )              * "    "
                sa = sa * TableStrings.flushright(18, @sprintf("%.6e", alphaK.Babushkin) )   * "    "
                println(stream, sa)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  >>> A perturber set is specific to ONE multipole: an np set serves E1 and reaches NOTHING")
        println(stream, "  >>> for E2, which needs same-parity nd levels.  And the discrete series stands in for the")
        println(stream, "  >>> continuum, so EVERY alpha^(Mp,0) here is a lower bound until it is seen to settle;")
        println(stream, "  >>> the higher the multipole, the deeper the series has to run.  Measured on H 1s against")
        println(stream, "  >>> the Dirac values of arXiv:1205.6883, on a box MATCHED to the perturbers:  E1 needs")
        println(stream, "  >>> np to n = 16 for 1.00000 of the reference, E2 nd to 16 for 0.99999, E3 nf to 22 for")
        println(stream, "  >>> 0.99999, E4 ng to 16 for 1.0031 -- while a few bound perturbers give only 0.81, 0.34,")
        println(stream, "  >>> 0.08 and 0.01 of it.  See MultipolePolarizibility.checkPerturberSet for the guards.")
        println(stream, "  >>> AND A REAL FREQUENCY HAS TWO TRAPS OF ITS OWN: above the lowest excitation energy the")
        println(stream, "  >>> true alpha is COMPLEX and this discrete real sum cannot represent its imaginary part,")
        println(stream, "  >>> and near an excitation energy the sum is carried by ONE term and is only as good as")
        println(stream, "  >>> that one energy.  MultipolePolarizibility.checkFrequency warns about both.")
        println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    end
    #
    return( nothing )
end


"""
`MultipolePolarizibility.isReachable(multipole::EmMultipole, level::Level, nuLevel::Level)`  
    ... returns true if the static 2^L-pole operator of `multipole` can connect `level` to the perturber `nuLevel` at
        all, i.e. if the two angular momenta satisfy the triangle |J - L| <= J' <= J + L AND the two parities match the
        parity of the operator, which is (-1)^L for an electric and (-1)^(L+1) for a magnetic multipole. This is the
        selection rule that decides whether a given gMultiplet is of any use for a given multipole -- an np set serves
        E1 and reaches nothing at all for E2. A value::Bool is returned.
"""
function isReachable(multipole::EmMultipole, level::Level, nuLevel::Level)
    Jd = Basics.twice(level.J) / 2.0;      Jnud = Basics.twice(nuLevel.J) / 2.0
    if  Jnud < abs(Jd - multipole.L) - 1.0e-10   ||   Jnud > Jd + multipole.L + 1.0e-10    return( false )    end
    wantsSameParity = multipole.electric ? iseven(multipole.L) : isodd(multipole.L)

    return( wantsSameParity == (level.parity == nuLevel.parity) )
end


"""
`MultipolePolarizibility.rankCoefficient(L::Int64, K::Int64, J::AngularJ64, Jnu::AngularJ64)`  
    ... returns the angular coefficient c(L,K,J,J') which multiplies |<J' || O^(Mp) || J>|^2 / (E_J' - E_J) in the
        rank-K static multipole polarizability alpha^(Mp,K)(J) of the module docstring, i.e.

            c(L,K,J,J')  =  2 (-1)^(J+J') C^(JJ)_(JJ,K0) (2K+1) / sqrt(2J+1) * (L L K; 0 0 0) * {J J' L; L K J}

        with a Wigner 3-j and a Wigner 6-j symbol. The three selection rules are carried by the formula itself rather
        than imposed on it: the 3-j symbol vanishes for odd K and for K > 2L, and the Clebsch-Gordan coefficient for
        K > 2J. For K = 0 the expression collapses analytically to 2 / ((2J+1)(2L+1)), which for L = 1 is the scalar
        prefactor of Angel & Sandars (1968) that this module has always used. A value::Float64 is returned.

        THE FORMULA WAS SETTLED NUMERICALLY RATHER THAN BY HAND, 19-Oct-2026, because a factor of this shape is exactly
        what gets written with the wrong phase and then compensated elsewhere (Rule 18). The substate polarizability
        alpha(J,M) was built directly from the Clebsch-Gordan coefficients, decomposed in the A_K shapes, and the
        coefficient read off: the decomposition residual is <= 1.7e-16 over 70 (J,L) combinations, which proves the
        rank set complete, and the ratio of the extracted coefficient to (-1)^(J+J') {J J' L; L K J} is INDEPENDENT of
        J' to <= 9.1e-15, which proves the 6-j structure. The remaining normalisation came out as 2 sqrt(2J+1) times
        the Clebsch-Gordan/3-j product above, and the phase as (-1)^(J+J') with ZERO mismatches over 421 samples.
"""
function rankCoefficient(L::Int64, K::Int64, J::AngularJ64, Jnu::AngularJ64)
    Jd = Basics.twice(J) / 2.0;      Jnud = Basics.twice(Jnu) / 2.0
    if  Jnud < abs(Jd - L) - 1.0e-10   ||   Jnud > Jd + L + 1.0e-10    return( 0.0 )    end
    if  isodd(K)   ||   K > 2*L   ||   K > round(Int64, 2*Jd)          return( 0.0 )    end

    gK  = AngularMomentum.ClebschGordan(Jd, Jd, 1.0*K, 0.0, Jd, Jd)
    c3j = AngularMomentum.Wigner_3j(1.0*L, 1.0*L, 1.0*K, 0.0, 0.0, 0.0)
    c6j = AngularMomentum.Wigner_6j(Jd, Jnud, 1.0*L, 1.0*L, 1.0*K, Jd)
    ph  = (-1.0)^round(Int64, Jd + Jnud)

    return( 2.0 * ph * gK * (2*K + 1.0) / sqrt(2*Jd + 1.0) * c3j * c6j )
end


"""
`MultipolePolarizibility.rankShape(J::AngularJ64, K::Int64, M::AngularM64)`  
    ... returns the angular shape A_K(J,M) = C^(JM)_(JM,K0) / C^(JJ)_(JJ,K0) in which the substate polarizability of a
        level is expanded, alpha(J,M) = sum_K A_K(J,M) alpha^(Mp,K)(J). The normalisation is fixed by A_0 = 1 and
        A_K(J,J) = 1, so that alpha^(Mp,0) is the usual scalar polarizability and every alpha^(Mp,K>0) is the anisotropy
        carried by the stretched M = J substate -- which is what makes the number directly readable. For K = 2 the
        shape is the familiar [3M^2 - J(J+1)] / [J(2J-1)] of the Stark-shift literature. A value::Float64 is returned.
"""
function rankShape(J::AngularJ64, K::Int64, M::AngularM64)
    if  K == 0   return( 1.0 )    end
    Jd = Basics.twice(J) / 2.0;      Md = Basics.twice(M) / 2.0
    if  K > round(Int64, 2*Jd)   return( 0.0 )    end

    return( AngularMomentum.ClebschGordan(Jd, Md, 1.0*K, 0.0, Jd, Md) /
            AngularMomentum.ClebschGordan(Jd, Jd, 1.0*K, 0.0, Jd, Jd) )
end


"""
`MultipolePolarizibility.displayContributions(stream::IO, outcome::MultipolePolarizibility.Outcome)`
    ... to display, for the single level outcome.level, the individual (additive) contribution of every
        intermediate ASF in gMultiplet that couples to it via E1, sorted by descending |contribAlpha0| so that
        the dominant perturber levels (e.g. 2p_1/2, 2p_3/2 for an ns ground level) appear first. A neat table is
        printed but nothing is returned otherwise. Nothing is printed if outcome.contributions is empty.
"""
function  displayContributions(stream::IO, outcome::MultipolePolarizibility.Outcome)
    if  isempty(outcome.contributions)   return( nothing )   end
    nx = 138
    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
    println(stream, " ")
    println(stream, "  Individual perturber contributions to alpha_0/alpha_2 of level $(TableStrings.level(outcome.level.index)) " *
                     "($(string(sym))), Coulomb and Babushkin gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Perturber"; na=2);                          sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J'^P'";     na=4);                          sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy";    na=4);
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Delta E";   na=4)
    sb = sb * TableStrings.center(14, "[a.u.]"; na=4)
    sa = sa * TableStrings.center(15, "d-alpha_0 [Coul]"; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_0 [Bab]" ; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_2 [Coul]"; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_2 [Bab]" ; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    ordered = sort(outcome.contributions, by = c -> -abs(c.contribAlpha0.Babushkin))
    for  c in ordered
        sa  = "  ";    nuSym = LevelSymmetry( c.nuLevel.J, c.nuLevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(c.nuLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(nuSym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", c.nuLevel.energy)) * "    "
        sa = sa * TableStrings.flushright(14, @sprintf("%.6e", c.deltaE) )                          * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha0.Coulomb) )           * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha0.Babushkin) )         * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha2.Coulomb) )           * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha2.Babushkin) )         * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module
