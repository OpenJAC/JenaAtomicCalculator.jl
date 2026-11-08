
"""
`module  JAC.PhotoDoubleIonization`  
... a submodel of JAC that computes photo-double ionization, i.e. the absorption of ONE photon with TWO electrons emitted into the
    continuum. The final state is built as a pair of outgoing partial waves: the first electron couples to the final ion to give an
    intermediate symmetry, the second couples to that to give the total symmetry of the scattering state, and a multipole must reach that
    total symmetry from the initial level. The amplitude is second order, with both time orderings of the knock-out (TS1) mechanism, and
    the intermediate states are levels of the (N-1)-electron ion with ONE continuum partial wave added.

    STATUS, 15-Sep-2026: THE PARKING IS LIFTED BY THE MAINTAINER -- see the closing paragraph of this block for
    what that does and does not mean. The 15-Aug text below is kept because it is the diagnosis, and it is still
    the best statement of what is wrong.

    STATUS, 15-Aug-2026 (SUPERSEDED, KEPT AS THE DIAGNOSIS): POSTPONED BY THE MAINTAINER. Do NOT pick this module up as a task without
    saying so
    first. It was rebuilt on this date from the parked quasi-shell version, and the rebuild is UNFINISHED:
    the structure works and passes several internal tests, but the ABSOLUTE SCALE IS WRONG BY THREE TO FOUR
    ORDERS OF MAGNITUDE and no number it produces means anything yet.

    DO NOT MISREAD THE STATE. That the module now runs end to end, returns finite numbers and passes its
    symmetry tests is NOT evidence that it is nearly right. The three self-tests below are necessary and
    cheap, and they were all passing while the result was still wrong by 10^3-10^4. In particular, a
    beautifully symmetric energy-sharing distribution says only that the pair set is exchange-closed.

    WHAT WORKS, and was verified on He (Z=2, 1s^2 -> He^2+ + 2e) at 200/400/800 eV:
      * The five-layer structure Line -> Sharing -> PartialWavePair -> Channel -> MultipoleAmplitude.
        For He at E1 every pair correctly serves the single total symmetry 1-, and the eight pairs at
        maxKappa=2 are the relativistic decomposition of the textbook eps-s eps-p and eps-p eps-d channels.
      * The double-ionization threshold comes out at 77.77 eV against 79.0 eV measured.
      * The energy-sharing distribution is mirror-symmetric to 1e-10 .. 1e-14, and EACH TIME ORDERING is
        separately symmetric, which is the stronger test: a defect in one can otherwise be masked by the other.
      * The distribution is U-shaped and deepens with photon energy, as it must far above threshold.
      * The total cross section reproduces the weighted sum of the differential ones to 0.00%.

    WHAT IS WRONG. sigma(2+)/sigma(+) comes out 0.001% at 200 eV rising to 0.006% at 800 eV, against a
    measured 3-4% near 200 eV falling to the asymptotic 1.68%. So the result is ~3600x too small AND trends
    the wrong way. The single-ionization cross section computed from the SAME orbitals, normalization and
    prefactor lands within 15% of the literature, so the deficit is specific to the second-order
    double-ionization amplitude, not to the machinery underneath it.

    THE REFERENCE THAT SETTLES THE NORMALIZATION is Kornberg & Miraglia, Phys. Rev. A 48, 3714 (1993),
    "Double photoionization of helium: use of a correlated two-electron continuum wave function". Their Eq. (2) reads

        d^5 sigma(2+) / (d eps1 d Omega1 d Omega2)  =  4 pi^2 alpha a0^2 k1 k2 C^(G) |e.T_fi^(G)(k1,k2)|^2

    with C^(L) = E_gamma, C^(V) = 1/E_gamma, C^(A) = E_gamma^(-3). This was read on 15-Aug-2026 as THREE missing factors.
    SETTLED 03-Nov-2026: only ONE of the three was real, and the other two must NOT be applied. Re-deriving it cost
    nothing and removed a factor of two from every number this module had produced, so the derivation is kept here.

      1. k1 k2 IS NOT MISSING -- it is the momentum-to-energy conversion that JAC has already performed. K&M normalize
         their final state to delta^3(k - k'); JAC normalizes its continuum orbitals to the amplitude sqrt(2/(pi q))
         (see Continuum.jl), which is per unit ENERGY. Since dE = q dq the two differ by |eps,Omega> = sqrt(k) |k>,
         hence |T_energy|^2 = k1 k2 |T_momentum|^2. Their k1 k2 IS that conversion. Applying it would convert twice.
         The "flattens the U by a factor 2.4" shape prediction recorded with it is therefore VOID -- do not test it.
      2. C^(G) IS NOT MISSING either -- JAC carries the gauge factor INSIDE the multipole operator, which is why ONE
         gauge-independent prefactor correctly serves both Coulomb and Babushkin here. A missing E_gamma^2 would be a
         factor 54 at 200 eV and could not hide behind a single-ionization cross section within 15% of the literature.
      3. THE FACTOR 1/2 WAS REAL, and it was supplied on 03-Nov-2026: "the correctly normalized total cross section is
         one half of the integral in the shake-off region". dsigma/d eps1 is defined over the whole range
         [0, omega - I_double], so integrating it counts every physical final state twice, once as (eps1, eps2) and
         once as (eps2, eps1). It is applied where the total is formed; see the comment there.

    THE SIGN OF THE DISCREPANCY IS WHAT IDENTIFIED IT, with no computation at all, and that is the method worth reusing.
    This amplitude carries knock-out (TS1) ONLY -- no shake-off, no ground-state correlation -- and an uncorrelated
    initial state underestimates the high-energy shake-off asymptote by a factor of 3. A model missing an entire
    mechanism cannot OVERSHOOT on physics. Ours overshot by 1.4-2.7x, which points at the prefactor and not at the
    physics. Where a discrepancy has the wrong SIGN for the known deficiency, look for a factor before looking further.

    THE OTHER SUSPECT is the intermediate space, which is severely truncated: the numbers above used four
    Green levels, maxKappa=1 and a 5-point Gauss-Legendre integral over the WHOLE excess-energy range, for an
    integrand peaked near threshold. The rise of the ratio with photon energy is the signature of exactly
    that. A convergence scan was started and did not finish; it is the first thing to run on resuming.

    WHAT IS NOT A PROBLEM. The Babushkin/Coulomb ratio of ~2.7 is NOT evidence of a defect: Kornberg &
    Miraglia report a large length-velocity discrepancy at all energies, with the two forms separated by
    about an order of magnitude over much of the range and velocity the trustworthy one at high energy.
    Do not spend time hunting it.

    THE PAPER WAS READ ON 04-Nov-2026, and three things in it change how it may be used. Mb/eV is the unit of
    their Fig. 2, and the column printed by `displayResults` was barn per HARTREE under a [barn]/[eV] header
    until that date.

    1. THEIR TABLE II HAS TWO DOUBLE-IONIZATION COLUMNS AND THE SMALLER ONE IS NOT THEIR MODEL. At 200 eV,
       sigma(+) = 5.44e-02 Mb, sigma(2+)_C2 = 1.69e-03 Mb and sigma(2+)_C3 = 1.3e-02 Mb. C3 is the correlated
       three-Coulomb-wave final state the paper is ABOUT; C2 is the independent-electron product of two Coulomb
       waves. They differ by 7.7x here. Anything quoted as "the K&M value" must say WHICH, and earlier notes in
       this repository used 1.69e-03 without saying it was C2.
    2. NEITHER COLUMN IS A GOOD BENCHMARK AT 200 eV -- the paper says so itself. "At intermediate energies, the
       calculated ratio is up to five times larger than the data", and C3 at 200 eV gives 23.9% against a
       measured 3-4%. At the other end, "at E_gamma = 84 eV the theory gives a ratio of 0.004% while the data is
       0.6%", a factor 140 low. Their velocity form is trustworthy only ABOVE about 1 keV, where it tends to
       1.75% against the 1.68% asymptote (C2 gives 1.60%). SO COMPARE AGAINST THE DATA, NOT AGAINST K&M, and use
       their Table II for sigma(+) and for the collected measurements rather than as a target.
    3. AND THE SHAPE AT 120 eV IS THE ONE PLACE THIS MODULE CAN BE TESTED AGAINST EXPERIMENT DIRECTLY. Their
       Fig. 2 plots d sigma(2+)/d eps1 at 120, 150 and 200 eV with the measured curve of Wehlitz et al. overlaid
       at 120 eV. The paper's own verdict: "At lower energies, the theoretical shape does not agree with the
       experimental U-shaped form, indicating the inadequacy of the present approximation in this range. However,
       as the energy increases the expected U-shaped form is reached; at 1.0 keV the curve is parabolic." So the
       measured distribution at 120 eV IS U-shaped and C3 does not reproduce it. This module gives an
       outermost/midpoint ratio of 1.81 at 120 eV, i.e. a U of the right sense where C3 has none -- a qualitative
       success, and not yet a quantitative one, on two counts. The 1.81 is itself a maxKappa = 2 number, and
       completing the partial-wave sum FLATTENS the U -- 5.55 -> 4.42 at 200 eV -- so the 120 eV run has to be
       repeated at maxKappa = 3 before its depth may be compared with anything. Fig. 2 was digitized on
       04-Nov-2026 and the measured depth is in hand; whether read-off values belong in this repository at all
       is the maintainer's call, so none are quoted here.

    ---------------------------------------------------------------------------------------------------------
    THE PARKING IS LIFTED, 15-Sep-2026, ON THE MAINTAINER'S DECISION. What changed, and what did not:

    FIXED SINCE 15-Aug: the photon was subtracted instead of added in one time ordering; the intermediate
    integral stopped at the wrong energy; the residue term was missing, so the amplitude was exactly real;
    the cross-section constant was wrong by 4 pi/(alpha^2 omega^2). And on 15-Sep the OUTPUT PATH: this
    module's `displayResults` was called by `computeLines` and had never been written, while `displayLines`
    read a `channels` field that the 15-Aug rebuild had moved onto the PartialWavePair. Either would raise.
    So until today the public entry point could not run to completion whatever else was right.

    THE ABSOLUTE SCALE IS NO LONGER 3-4 ORDERS OUT; THE 15-Aug TEXT ABOVE IS STALE ON THIS POINT. That figure
    predates the four fixes of 06-Sep, one of which was a cross-section constant wrong by 4 pi/(alpha^2 omega^2).
    THE RATIO, MEASURED SELF-CONSISTENTLY ON 04-Nov-2026 -- sigma(+) and sigma(2+) from the SAME orbitals, grid,
    normalization and prefactor, He at 200 eV, box 40 a.u., 9 sharings, maxKappa = 2, NoIntermediateEnergies = 9
    (tools/comparePdiHeliumKornberg.jl), against Kornberg & Miraglia Table II at 200 eV:

                               JAC Coulomb     JAC Babushkin       K&M C2       K&M C3
        sigma(+)   [Mb]         5.9499e-02        5.1565e-02     5.4400e-02   5.4400e-02
        sigma(2+)  [Mb]         2.3015e-03        4.2347e-03     1.6900e-03   1.3000e-02
        sigma(2+)/sigma(+)          3.868 %           8.212 %        3.107 %      23.90 %

    THE COMPARISON THAT COUNTS IS WITH THE MEASUREMENT, not with either K&M column -- AND IT DOES NOT AGREE ONCE
    THE PARTIAL-WAVE SUM IS CONVERGED. The maxKappa = 2 row above gives 3.868 % in velocity, inside the measured
    3-4 % band. That reading was written here on 04-Nov-2026 and WITHDRAWN the same day: converging maxKappa
    (2 -> 3 -> 4 moves the total +21.6 %, then +0.8 %, with a tail under 0.1 %) gives

        sigma(2+)                2.8208e-03 Mb velocity      6.7563e-03 Mb length
        sigma(2+)/sigma(+)           4.741 % velocity            13.10 % length      (5.185 % on K&M's sigma(+))

    so velocity sits ABOVE the measured band by about a fifth rather than inside it. What looked like agreement
    was a 22 % truncation error cancelling part of a real overshoot -- and the lesson outlasts the number: A
    RESULT LANDING INSIDE A MEASURED BAND IS NOT EVIDENCE UNTIL THE CONVERGENCE ERROR IS SMALLER THAN THE BAND.
    Landing between C2 and the data is still the right place to land -- K&M's own C3 gives 23.9 % there, and
    their paper says the ratio is "up to five times larger than the data" at intermediate energies, so agreeing
    with C3 would still have been a failure -- but this is an ORDER-AND-SIGN result and must not be quoted as
    agreement. IT IS NO LONGER AN OPEN DEFECT ON THE PRIORITY LIST: it was carried as item 6 (Nov-2026
    numbering) and CLOSED on 05-Nov-2026, parked as challenge S35 with the question OPEN, after BOTH of its
    named suspects had been tested and exonerated -- the bound Green space and the shared orbital set, below.
    It returns if an explicit application needs an absolute helium double-ionization cross section, or any
    absolute second-order amplitude from this junction. Three cautions:
      * EVERY NUMERICAL TRUNCATION IS NOW MEASURED AND EVERY ONE IS UNDER 1 %: box 0.05 % over 24 -> 40 a.u.,
        sharings +0.58 % over 5 -> 21, intermediate mesh -0.87 % then +0.24 % over 9 -> 17 -> 25 (non-monotone,
        hence converged), partial waves +0.8 % over maxKappa 3 -> 4. Against a 20-50 % discrepancy with the
        measurement, so THE REMAINING ERROR IS IN THE MODEL AND NOT IN THE ARITHMETIC. The module is knock-out
        only and a model missing a mechanism cannot overshoot on physics, WHICH LEAVES THE SECOND-ORDER
        CONSTRUCTION ITSELF. The four-level (N-1) Green space was the other candidate and was TESTED late
        on 04-Nov-2026: opened to nine levels (He+ n = 3) on an orbital span widened to n = 3 and held
        COMMON to both runs, so that the Green space was the only thing varied, it moved sigma(2+) by
        -0.01 % and Babushkin/Coulomb from 2.377 to 2.374, and no energy-sharing node by more than
        0.13 %. It accounts for NEITHER the overshoot nor the gauge gap. One caveat on that null: the
        n = 3 orbitals the configuration-averaged mean field generates are over-screened -- Z_eff = 0.67,
        0.56, 0.44 for 3s, 3p, 3d against the 2 a He+ state should see -- so what is exonerated is the
        Green space AS THIS SHARED ORBITAL SET CAN REPRESENT IT. Widening the span itself did move the
        total, +3.7 % Coulomb and +5.1 % Babushkin, which is a basis effect and not a Green-space one.
      * AND THE SHARED ORBITAL SET IS NOT THE LEVER EITHER. The V^(N-1) test that the caveat above asked for
        was built and run on 05-Nov-2026, together with four further sets, each converged at maxKappa = 3:

            set        mean field / span       eps(1s)   1st-order gauge   2nd-order gap   Bab/Cou   sigma(2+)
            V^(N-1)    AL on 1s,2s,2p,n=3      -2.0001        4.19 %          138.56 %      2.386    4582 barn
            AL n=3     AL on 1s^2 + n=2,n=3    -1.7144        6.99 %          121.29 %      2.213    3933
            DFS n=3    DFS on 1s^2 + n=2,n=3   -1.2103       13.13 %          137.70 %      2.377    2902
            DFS n=2    DFS on 1s^2 + n=2       -1.1008       13.33 %          132.70 %      2.346    2798
            HS n=2     HS on 1s^2 + n=2        -1.0033       17.47 %          147.85 %      2.479    2268

        The first-order gauge deviation falls by a FACTOR 4.2 and the second-order gap does not follow it:
        121-148 %, a factor 1.11, and not even monotone -- the best first-order set gives 138.56 % while
        DFS n=2, three times worse at first order, gives 132.70 %. A straight-line fit gives r = 0.551,
        which on five points is indistinguishable from noise, and a slope of 0.999 points per point. Take
        that trend entirely at face value -- the reading most generous to the orbital hypothesis -- and the
        gap at ZERO first-order error is still 124.6 %, so the conclusion survives whether or not the trend
        is believed: THE ORBITAL AXIS CANNOT CLOSE THIS GAP.
      * WHAT THE ORBITAL SET DOES CONTROL IS THE MAGNITUDE, and that is the more useful half: sigma(2+) runs
        2268 -> 4582 barn, a factor 2.02, STRICTLY MONOTONE in eps(1s), while the gauge gap wanders in a
        10 % band. The two symptoms therefore SEPARATE -- whatever is wrong multiplies both gauges by an
        orbital-independent ~2.4 while the absolute scale rides on the 1s binding -- so a candidate
        mechanism must be nearly gauge-symmetric and nearly orbital-independent. That is a sharper
        constraint on the defect than the overshoot alone.
      * AND THE TWO ANCHORS CANNOT BE SATISFIED TOGETHER, which anyone retuning the orbitals must know
        first. Gauge consistency of <eps|r|1s> optimizes at the He+ limit, eps(1s) = -2.0; the observable
        sigma(+) optimizes near eps(1s) = -1.0, where HS reproduces K&M to -0.0 % and gives the closest
        ratio of the five, 4.172 %; and the exact He ground state (AL on 1s^2 alone, eps(1s) = -0.917942)
        is at the wrong end of BOTH. No single shared set satisfies both, and the proper treatment -- two
        sets with a biorthogonal transformation -- collides head-on with the one-shared-orbital-set
        precondition this construction rests on. The worked record is examples/example-Dj.jl branch (d).
      * The LENGTH form is far outside the band and CONVERGING THE SUM MAKES IT WORSE: 8.2 % at maxKappa = 2
        against 13.10 % at maxKappa = 4, where velocity moved 3.868 -> 4.741 %. That is the direction K&M's own
        large length-velocity discrepancy predicts. Only velocity is being compared -- and Babushkin/Coulomb
        was itself STILL MOVING at the edge of the scan, 2.346 -> 2.395 over maxKappa 3 -> 4, so the gauge
        ratio is not converged even where the TOTAL is.
      * THE EARLIER "2.70 % at 200 eV" WAS NOT MERELY A FACTOR 2 OUT. Halving it to 1.35 % is correct arithmetic
        and still wrong, because its sigma(+) was never checked and was about 3x the literature; the two errors
        partly cancelled. A ratio is only worth quoting when BOTH members come from the same calculation, and the
        04-Nov table above is the first time that has been true here. The 400 eV and 800 eV entries of that old
        series are withdrawn rather than halved: they have no companion sigma(+).
    "Three orders of magnitude too small" is long dead and must not be repeated.

    WHAT IS STILL NOT TRUSTWORTHY, and the caller is warned of it at run time rather than stopped:
      * THE SHARING INTEGRAL CONVERGES; THE "+29 %" ABOVE WAS AN ARTEFACT OF TWO CONFOUNDS AND IS WITHDRAWN.
        The 15-Sep scan varied two things at once. (a) `NoIntermediateEnergies` did not exist until 30-Sep-2026,
        so asking for more sharings silently made the second-order sum more complete too -- worth a factor 2.7 at
        FIXED sharing. (b) The asymptotic amplitude was read as a maximum over the last 300 grid points until
        03-Nov-2026, a fixed POINT count that spans less than half a wavelength below ~14 eV, and the added
        Gauss-Legendre nodes march towards exactly that slow-electron band. With the intermediate mesh held fixed
        and the amplitude FITTED, the total moves +0.58 % over 5 -> 21 points. THE LEADING ERROR IS THE
        PARTIAL-WAVE SUM, identified 04-Nov-2026 and converged at maxKappa = 3; everything else is under 1 %.
      * THE LENGTH (Babushkin) FORM, AND WHICH TRUNCATION ACTUALLY MOVES IT. The 15-Sep reading -- that
        quadrupling the intermediate set leaves velocity under 1 % at 200 eV while length nearly DOUBLES --
        attributed the near-doubling to the wrong mesh, because maxKappa was never varied beside it. Measured
        04-Nov-2026 at 200 eV, box 40 a.u., 9 sharings, one variable per run: 9 -> 17 -> 25 intermediate
        energies move Coulomb -0.87 %, +0.24 % and Babushkin -1.09 %, +0.56 %, leaving the gauge ratio at
        1.840, 1.836, 1.842 -- so THAT mesh is converged at 9, and it is blind to the gauge gap. maxKappa
        2 -> 3 -> 4 moves Coulomb +21.6 %, +0.8 % and Babushkin +55.0 %, +2.9 %, and walks the gauge ratio to
        2.346 then 2.395. The length form IS the unstable one and the gauge RATIO is still not a convergence
        indicator, but what destabilizes it is the PARTIAL-WAVE sum and not the intermediate mesh. This agrees
        with Kornberg & Miraglia, who find velocity the trustworthy form at high energy.
      * THE ENERGY TREND is not monotonic in either gauge; see the numbers above.

    WHY THE TWO GAUGES ARE NOT EXPECTED TO AGREE, which is a property of second-order amplitudes and not a
    defect of this module. The length and velocity forms of a one-photon matrix element are connected by the
    commutator identity <a|p|b> = i (E_a - E_b) <a|r|b>, i.e. by the LEVEL DIFFERENCE. In a second-order
    amplitude the photon operator carries the PHYSICAL frequency omega, while the intermediate state is
    off shell: E_n - E_i is not omega. The difference between the two forms is then governed by the factor
    (E_i - E_n - omega), which vanishes only ON SHELL, and the two gauges are equivalent only through
    CLOSURE over a COMPLETE intermediate set. Putting omega in the operator is the correct S-matrix
    prescription and is what this module does; it is ALSO the reason the two forms separate here. See
    Fortun et al., arXiv:1610.03571, for the master identity and the statement that extended gauge
    invariance holds only at exact resonance. Do not "fix" this by replacing omega with a level difference:
    that would restore gauge agreement matrix element by matrix element while making the amplitude wrong.

    THE SUITE NOW GUARDS THIS MODULE. `TestFrames.testModule_PhotoDoubleIonization` asserts the output path,
    the mirror symmetry of the sharing distribution, and that the quadrature adds up -- three invariants that
    hold whatever the scale turns out to be. It deliberately asserts NOTHING about the absolute value and
    nothing about gauge equality, for the reason given above.
"""
module PhotoDoubleIonization


using Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..Radial, ..RadialIntegrals, ..Nuclear,
                ..ManyElectron, ..PhotoEmission, ..TableStrings

"""
`struct  PhotoDoubleIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing photo-double ionization lines.

    + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
    + photonEnergies          ::Array{Float64,1}        ... List of photon energies [in user-selected units].  
    + electronEnergyShift     ::Float64                 ... An overall energy shift of the two emitted electrons [in user-selected units];
        it corrects the total electron energy for a shift of the i-f transition energy and thereby separates the energies at which the
        continuum orbitals are generated from the sharing coordinates.
    + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
    + intermediateEnergyFactor ::Float64               ... The top of the INTERMEDIATE energy integral, as a multiple of the on-shell
                                                           energy E_i + omega - E_green(min).  It must EXCEED 1, or the pole of the
                                                           photon-first ordering sits exactly on the upper endpoint, where no principal
                                                           value exists; and it is a convergence parameter, because the intermediate
                                                           state is VIRTUAL and its electron energy is not bounded by energy
                                                           conservation. The grid must support a continuum orbital at the resulting top
                                                           energy.  Default 3.
    + NoIntermediateEnergies  ::Int64                   ... Number of points in the integral over the INTERMEDIATE electron energy of the
        second-order amplitude. This is a different convergence question from `NoEnergySharings` and must be set separately: the sharing
        mesh resolves the OUTER integral over how the excess energy is divided between the two emitted electrons, whereas this mesh decides
        how completely the second-order sum over intermediate states is carried out. **Until 30-Sep-2026 both were taken from
        `NoEnergySharings`**, so asking for more sharings silently made the intermediate sum more complete as well, and the differential
        cross section at a FIXED sharing then grew with the number of sharings -- measured at equal sharing, a factor 2.7 between 3 and 41
        sharings, with the partial-wave set unchanged. 0 means fall back to `NoEnergySharings`, i.e. the old, entangled behaviour.
    + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
    + calcDifferentialCs      ::Bool                    ... True, if the energy-differential cs are to be calculated and false otherwise.  
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + eeInteraction           ::AbstractEeInteraction   ... Type of the electron-electron interaction in the second-order treatment.
    + gMultiplet              ::Multiplet               ... Mean-field multiplet of the (N-1)-ELECTRON IONIC CORES of the
        intermediate levels, NOT of the intermediate levels themselves.  `PhotoDoubleIonization.generateIntermediateLevels` builds each
        intermediate by generating a continuum orbital IN THE FIELD OF one of these levels and attaching it, so a multiplet of
        N-electron intermediates handed in here yields (N+1)-electron intermediates and every junction of the second-order sum then
        pairs unequal electron numbers.  This field said "intermediate levels" until 01-Nov-2026 and `examples/example-Dj.jl` followed
        it: for the two-electron 1s^2 initial state it supplied the two-electron configurations 1s4s and 1s4p, and the code made them
        three-electron -- measured, and the reason that branch has never carried a `Last successful` date.  The same role is played, and
        was always correctly documented, by `AutoIonization.Settings.gCoreMultiplet`.
    + continuumExchange       ::Union{Nothing,Bool}     ... True/false to force exchange between the two ejected electrons
        and the bound electrons on or off; `nothing` -- the default -- lets `Basics.perform` decide automatically from the
        final-ion `AsfSettings.scField` (`Continuum.resolveExchangeDefault`, item 99).
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    photonEnergies            ::Array{Float64,1}
    electronEnergyShift       ::Float64
    NoEnergySharings          ::Int64
    NoIntermediateEnergies    ::Int64
    intermediateEnergyFactor  ::Float64
    maxKappa                  ::Int64
    calcDifferentialCs        ::Bool
    printBefore               ::Bool
    lineSelection             ::LineSelection
    eeInteraction             ::AbstractEeInteraction
    gMultiplet                ::Multiplet
    continuumExchange         ::Union{Nothing,Bool}
end


"""
`PhotoDoubleIonization.Settings()`  ... constructor for the default values of PhotoDoubleIonization line computations
        A `settings::PhotoDoubleIonization.Settings` is returned.
"""
function Settings()
    Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], 0., 0, 0, 3., 0, false, false,
                LineSelection(), CoulombInteraction(), Multiplet(), nothing)
end


"""
`PhotoDoubleIonization.Settings(set::PhotoDoubleIonization.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergyShift::Union{Nothing,Float64}=nothing,
    NoEnergySharings::Union{Nothing,Int64}=nothing,                         NoIntermediateEnergies::Union{Nothing,Int64}=nothing,
    intermediateEnergyFactor::Union{Nothing,Float64}=nothing,
    maxKappa::Union{Nothing,Int64}=nothing,                                 calcDifferentialCs::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                               lineSelection::Union{Nothing,LineSelection}=nothing,
    eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,            gMultiplet::Union{Nothing,Multiplet}=nothing,
    continuumExchange::Union{Nothing,Bool}=nothing)`

    ... constructor for modifying the given PhotoDoubleIonization.Settings by 'overwriting' the previously selected parameters.
        A `settings::PhotoDoubleIonization.Settings` is returned.
"""
function Settings(set::PhotoDoubleIonization.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergyShift::Union{Nothing,Float64}=nothing,
    NoEnergySharings::Union{Nothing,Int64}=nothing,                         NoIntermediateEnergies::Union{Nothing,Int64}=nothing,
    intermediateEnergyFactor::Union{Nothing,Float64}=nothing,
    maxKappa::Union{Nothing,Int64}=nothing,                                 calcDifferentialCs::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                               lineSelection::Union{Nothing,LineSelection}=nothing,
    eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,            gMultiplet::Union{Nothing,Multiplet}=nothing,
    continuumExchange::Union{Nothing,Bool}=nothing)

    if  isnothing(multipoles)           multipolesx         = set.multipoles         else  multipolesx         = multipoles          end
    if  isnothing(gauges)               gaugesx             = set.gauges             else  gaugesx             = gauges              end
    if  isnothing(photonEnergies)       photonEnergiesx     = set.photonEnergies     else  photonEnergiesx     = photonEnergies      end
    if  isnothing(electronEnergyShift)  electronEnergyShiftx= set.electronEnergyShift else electronEnergyShiftx= electronEnergyShift end
    if  isnothing(NoEnergySharings)     NoEnergySharingsx   = set.NoEnergySharings   else  NoEnergySharingsx   = NoEnergySharings    end
    if  isnothing(NoIntermediateEnergies)  NoIntermediateEnergiesx = set.NoIntermediateEnergies  else  NoIntermediateEnergiesx = NoIntermediateEnergies  end
    if  isnothing(intermediateEnergyFactor)  intermediateEnergyFactorx = set.intermediateEnergyFactor  else  intermediateEnergyFactorx = intermediateEnergyFactor  end
    if  isnothing(maxKappa)             maxKappax           = set.maxKappa           else  maxKappax           = maxKappa            end
    if  isnothing(calcDifferentialCs)   calcDifferentialCsx = set.calcDifferentialCs else  calcDifferentialCsx = calcDifferentialCs  end
    if  isnothing(printBefore)          printBeforex        = set.printBefore        else  printBeforex        = printBefore         end
    if  isnothing(lineSelection)        lineSelectionx      = set.lineSelection      else  lineSelectionx      = lineSelection       end
    if  isnothing(eeInteraction)        eeInteractionx      = set.eeInteraction      else  eeInteractionx      = eeInteraction       end
    if  isnothing(gMultiplet)           gMultipletx         = set.gMultiplet         else  gMultipletx         = gMultiplet          end
    if  isnothing(continuumExchange)    continuumExchangex  = set.continuumExchange  else  continuumExchangex  = continuumExchange   end

    Settings( multipolesx, gaugesx, photonEnergiesx, electronEnergyShiftx, NoEnergySharingsx, NoIntermediateEnergiesx,
                intermediateEnergyFactorx, maxKappax,
                calcDifferentialCsx, printBeforex, lineSelectionx, eeInteractionx, gMultipletx, continuumExchangex)
end


function Base.show(io::IO, settings::PhotoDoubleIonization.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "electronEnergyShift:      $(settings.electronEnergyShift)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "NoEnergySharings:         $(settings.NoEnergySharings)  ")
    println(io, "NoIntermediateEnergies:   $(settings.NoIntermediateEnergies)  ")
    println(io, "intermediateEnergyFactor: $(settings.intermediateEnergyFactor)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "calcDifferentialCs:       $(settings.calcDifferentialCs)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "eeInteraction:            $(settings.eeInteraction)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
    println(io, "continuumExchange:        $(settings.continuumExchange)  ")
end


"""
`struct  PhotoDoubleIonization.Channel`
    ... defines a type for one total symmetry of the complete final state, i.e. of the final ion together with both emitted electrons, and
        collects the multipole amplitudes that reach this symmetry from the initial level.

    + symmetry       ::LevelSymmetry                    ... total angular momentum and parity of the complete scattering state.
    + amplitudes     ::Array{MultipoleAmplitude,1}      ... one entry per contributing multipole, each holding BOTH gauges.
"""
struct  Channel
    symmetry         ::LevelSymmetry
    amplitudes       ::Array{MultipoleAmplitude,1}
end


function Base.show(io::IO, channel::PhotoDoubleIonization.Channel)
    println(io, "symmetry:               $(channel.symmetry)  ")
    println(io, "amplitudes:             $(channel.amplitudes)  ")
end


"""
`struct  PhotoDoubleIonization.PartialWavePair`
    ... defines a type for one pair of outgoing partial waves, i.e. for the two electrons that leave the ion together. The two electrons are
        coupled in sequence: the first couples to the final ion to give xSymmetry, and the second couples to that intermediate symmetry to
        give the total symmetry of a channel. One pair generally serves several total symmetries, which is why the channels are held here
        rather than the pair being repeated for each of them; the two continuum orbitals are therefore generated once per pair.

    + kappa1         ::Int64                ... partial wave of the first emitted electron.
    + energy1        ::Float64              ... energy at which the orbital of the first electron is generated.
    + phase1         ::Float64              ... scattering phase of the first partial wave.
    + xSymmetry      ::LevelSymmetry        ... intermediate symmetry of (final ion + electron 1).
    + kappa2         ::Int64                ... partial wave of the second emitted electron.
    + energy2        ::Float64              ... energy at which the orbital of the second electron is generated.
    + phase2         ::Float64              ... scattering phase of the second partial wave.
    + channels       ::Array{PhotoDoubleIonization.Channel,1}   ... total symmetries that this pair of partial waves serves.
"""
struct  PartialWavePair
    kappa1           ::Int64
    energy1          ::Float64
    phase1           ::Float64
    xSymmetry        ::LevelSymmetry
    kappa2           ::Int64
    energy2          ::Float64
    phase2           ::Float64
    channels         ::Array{PhotoDoubleIonization.Channel,1}
end


function Base.show(io::IO, pair::PhotoDoubleIonization.PartialWavePair)
    println(io, "kappa1:                 $(pair.kappa1)  ")
    println(io, "energy1:                $(pair.energy1)  ")
    println(io, "phase1:                 $(pair.phase1)  ")
    println(io, "xSymmetry:              $(pair.xSymmetry)  ")
    println(io, "kappa2:                 $(pair.kappa2)  ")
    println(io, "energy2:                $(pair.energy2)  ")
    println(io, "phase2:                 $(pair.phase2)  ")
    println(io, "channels:               $(pair.channels)  ")
end


"""
`struct  PhotoDoubleIonization.Sharing`
    ... defines a type for one division of the excess energy between the two emitted electrons. The list of sharings of a line IS the
        energy-sharing distribution, which is the primary observable of this process, and each sharing owns two quantities that belong
        nowhere else: the quadrature weight of this sharing point, and the differential cross section summed over all pairs and channels
        at this point.

        The sharing coordinates epsilon1, epsilon2 are what Basics.determineEnergySharings produces and define the quadrature point. They
        are NOT the same as the energies at which the orbitals are generated, which additionally carry settings.electronEnergyShift and
        are held by the partial-wave pair.

    + omega          ::Float64         ... energy of the incident photon.
    + epsilon1       ::Float64         ... sharing coordinate of (free) electron 1.
    + epsilon2       ::Float64         ... sharing coordinate of (free) electron 2.
    + weight         ::Float64         ... Gauss-Legendre weight of this sharing for energy-integrated quantities.
    + differentialCs ::EmProperty      ... energy-differential cross section dsigma/depsilon1 at this sharing.
    + partialWavePairs ::Array{PhotoDoubleIonization.PartialWavePair,1}   ... pairs of partial waves at this sharing.
"""
struct  Sharing
    omega            ::Float64
    epsilon1         ::Float64
    epsilon2         ::Float64
    weight           ::Float64
    differentialCs   ::EmProperty
    partialWavePairs ::Array{PhotoDoubleIonization.PartialWavePair,1}
end


function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)
    println(io, "omega:                  $(sharing.omega)  ")
    println(io, "epsilon1:               $(sharing.epsilon1)  ")
    println(io, "epsilon2:               $(sharing.epsilon2)  ")
    println(io, "weight:                 $(sharing.weight)  ")
    println(io, "differentialCs:         $(sharing.differentialCs)  ")
    println(io, "partialWavePairs:       $(sharing.partialWavePairs)  ")
end


"""
`struct  PhotoDoubleIonization.Line`
    ... defines a type for a photo-double ionization line between an initial and a final level, for one photon energy, in which a single
        photon is absorbed and TWO electrons are emitted into the continuum. The cross section is the sharing-integrated total.

    + initialLevel   ::Level                  ... initial-(state) level.
    + finalLevel     ::Level                  ... final-(state) level of the doubly-ionized ion.
    + photonEnergy   ::Float64                ... energy of the absorbed photon.
    + crossSection   ::EmProperty             ... total cross section, integrated over all energy sharings.
    + sharings       ::Array{PhotoDoubleIonization.Sharing,1}  ... the energy-sharing distribution of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    photonEnergy     ::Float64
    crossSection     ::EmProperty
    sharings         ::Array{PhotoDoubleIonization.Sharing,1}
end


"""
`PhotoDoubleIonization.Line()`  ... constructor for an empty PhotoDoubleIonization line; a line::PhotoDoubleIonization.Line is returned.
"""
function Line()
    Line(Level(), Level(), 0., EmProperty(0., 0.), PhotoDoubleIonization.Sharing[])
end


function Base.show(io::IO, line::PhotoDoubleIonization.Line)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "sharings:          $(line.sharings)  ")
end


"""
`PhotoDoubleIonization.amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                                 nLevels::Array{Level,1}, nWeights::Array{Float64,1}, grid::Radial.Grid; display::Bool=false,
                                 printout::Bool=false,
                                 nOnShell::Union{Nothing,Array{Bool,1}}=nothing, nBlock::Union{Nothing,Array{Int64,1}}=nothing,
                                 nEnergy::Union{Nothing,Array{Float64,1}}=nothing,
                                 energyRange::Union{Nothing,Tuple{Float64,Float64}}=nothing)`
    ... computes the second-order photo-double ionization amplitude

                <alpha_f J_f || O^(Mp, absorption) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>  / (E_i - E_n)
            +   <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, absorption) || alpha_i J_i>  / (E_i - omega - E_n)

        for the interaction with a photon of frequency omega, multipolarity Mp and the given gauge; a value::ComplexF64 is returned.
        The first term is the ordering in which the two electrons interact first and the photon is absorbed afterwards, the second the
        ordering in which the photon is absorbed first. Together they are the knock-out (TS1) mechanism.

        The intermediate levels are handed in rather than taken from a multiplet, because they are not bound states: each of them is a
        level of the (N-1)-electron ion with ONE continuum partial wave added, and such a state has one continuum electron exactly as the
        final state has two. A Green multiplet of bound levels cannot contribute at all, since it differs from the final state in two
        orbitals and the one-body photon operator connects only one.

        Because those partial waves are normalized PER ENERGY INTERVAL, so are the intermediate levels, and the sum over them is really an
        integral over the intermediate electron energy. nWeights carries the quadrature weight of each intermediate level and must have the
        same length as nLevels; passing ones would silently turn the integral back into a plain sum.

        The intermediate level is in general multi-configurational, so the two matrix elements of each ordering run over two INDEPENDENT
        CSF indices t and tp; a single index with mc[t]^2 would keep only the diagonal terms. Each ordering has its own symmetry condition
        and the two are mutually exclusive whenever symi != symf, which is why neither may skip the loop over intermediate levels.
"""
function amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                   nLevels::Array{Level,1}, nWeights::Array{Float64,1}, grid::Radial.Grid; display::Bool=false, printout::Bool=false,
                   nOnShell::Union{Nothing,Array{Bool,1}}=nothing, nBlock::Union{Nothing,Array{Int64,1}}=nothing,
                   nEnergy::Union{Nothing,Array{Float64,1}}=nothing,
                   energyRange::Union{Nothing,Tuple{Float64,Float64}}=nothing)
    length(nLevels) == length(nWeights)  ||  error("nLevels and nWeights must have the same length.")
    onShell = isnothing(nOnShell) ? falses(length(nLevels)) : nOnShell
    length(nLevels) == length(onShell)   ||  error("nLevels and nOnShell must have the same length.")
    if  length(nLevels) == 0    return( ComplexF64(0.) )    end
    # THE INTERMEDIATE MUST CARRY THE SAME NUMBER OF ELECTRONS AS THE INITIAL STATE, and this is checked rather than assumed because
    # the one thing that must not happen is a plausible number from a sum whose junctions pair different electron counts.  It is a
    # caller error and not a physical limitation: `settings.gMultiplet` must hold the (N-1)-electron ionic CORES, to which
    # generateIntermediateLevels attaches the continuum electron -- see that field's docstring.
    for  nLevel in nLevels
        if  nLevel.basis.NoElectrons != initialLevel.basis.NoElectrons
            error("PhotoDoubleIonization.amplitude(): an intermediate level carries $(nLevel.basis.NoElectrons) electrons where the "  *
                  "initial level carries $(initialLevel.basis.NoElectrons).  Every junction of the second-order sum would then pair "   *
                  "unequal electron numbers.  `settings.gMultiplet` must be the (N-1)-ELECTRON IONIC CORES, not the intermediate "     *
                  "levels themselves:  generateIntermediateLevels generates a continuum orbital in the field of each of its levels "   *
                  "and attaches it.  See the docstring of PhotoDoubleIonization.Settings' gMultiplet field.")
        end
    end
    # The principal value needs the levels grouped by their pole and the integration interval they span; without
    # all three the old, non-convergent assembly is used, so an un-updated caller changes no number silently.
    doPV = !isnothing(nBlock)  &&  !isnothing(nEnergy)  &&  !isnothing(energyRange)
    doPV  &&  ( length(nBlock) == length(nLevels) && length(nEnergy) == length(nLevels)  ||
                error("nBlock and nEnergy must have the same length as nLevels.") )

    # Always ensure the same subshell list for all initial, intermediate and final levels.  Each intermediate
    # level carries a continuum subshell of its own kappa, so every one of them has to enter the merge.
    subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
    for  nLevel in nLevels   subshells = Basics.merge(subshells, nLevel.basis.subshells)   end
    iLevel    = Level(initialLevel, subshells)
    fLevel    = Level(finalLevel,   subshells)

    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity);    eni = iLevel.energy

    if  printout   printstyled("Compute photo-double $(Mp) ionization amplitude for the transition " *
                               "[$(iLevel.index)-$(fLevel.index)] ... ", color=:light_green)    end
    amplitude = ComplexF64(0.)
    # The photon-first ordering is the one with the pole, so its numerators are accumulated PER INTERMEDIATE
    # LEVEL and assembled afterwards; the electron-first ordering has no pole and is summed in place as before.
    numer     = zeros(ComplexF64, length(nLevels))

    for  r = 1:nf
        symr = LevelSymmetry(fLevel.basis.csfs[r].J, fLevel.basis.csfs[r].parity);      if  symr != symf    continue    end
        for  s = 1:ni
            syms = LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity);  if  syms != symi    continue    end
            for  (k, nLevel) in enumerate(nLevels)
                nLevel = Level(nLevel, subshells)
                symn   = LevelSymmetry(nLevel.J, nLevel.parity);    enn = nLevel.energy;    wn = nWeights[k]
                nn     = length(nLevel.basis.csfs)

                #   Compute <alpha_f J_f || O^(Mp, kind) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
                #   THIS JUNCTION REACHES CSF PAIRS DIFFERING IN ONE ORBITAL, WHERE V^(e-e) IS NOT THE WHOLE
                #   OPERATOR -- and it nevertheless STAYS, measured rather than assumed, 02-Nov-2026.  The
                #   operator a one-orbital junction requires is H - E, and on example-Dj.jl's Green core 1s at
                #   omega = 400 eV and equal sharing the two read V_ee = +2.470622e-01 against H - E =
                #   +1.560546e-01, a factor 1.583, the one-body part contributing -9.100758e-02 with the right
                #   sign.  But ManyElectron.matrixElement_HminusE drops the overlap term -E <r|s>, and here that
                #   term is the LARGEST piece: the continuum orbital is generated in the field of the GREEN
                #   level, so <101s_1/2|1s_initial> = -0.0235 rather than 0, and |E_i <r|s>| ~ 3.2e-01, i.e.
                #   2.05x the whole H - E element.  Swapping the operator would add 0.091 and omit 0.32.  The
                #   precondition is therefore CHECKED instead, in checkSecondOrderConditions, and the swap waits
                #   for a case whose Green space shares the initial state's orbitals.  The junctions differing in
                #   TWO orbitals are unaffected either way -- the two operators agree there to 1.000 exactly.
                if  symn == symi
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        OMp = ManyElectron.matrixElement_Mab(Mp, gauge, omega, fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            Vee       = ManyElectron.matrixElement_Vee(CoulombInteraction(), nLevel.basis, tp, iLevel.basis, s, grid)
                            amplitude = amplitude + wn * fLevel.mc[r] * OMp * nLevel.mc[t] * nLevel.mc[tp] * Vee *
                                                    iLevel.mc[s] / (eni - enn)
                        end
                    end
                end

                #   Compute <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, kind) || alpha_i J_i>
                #   THE SIGN OF omega HERE IS THE ONE THAT DECIDES WHETHER THIS MODULE HAS A POLE.  In this time
                #   ordering the PHOTON ACTS FIRST, so the intermediate state has already absorbed it and the
                #   denominator is E_i + omega - E_n; the other ordering above absorbs no photon before reaching
                #   |n> and correctly reads E_i - E_n.  Until 6-Sep-2026 this line carried -omega, which is the
                #   EMISSION denominator, and it moved the whole integrand off shell: measured at omega = 200 eV
                #   it swept -226 to -344 eV and never approached zero, so the intermediate integral was smooth,
                #   converged to six figures in eleven points, and described no resonance.  With +omega the
                #   denominator passes through zero at exactly one electron photoionized and the other still
                #   bound -- which IS the TS1 mechanism this module is built on -- and that on-shell region
                #   needs a contour prescription rather than plain quadrature.
                if  symn == symf
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        Vee = ManyElectron.matrixElement_Vee(CoulombInteraction(), fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            OMp       = ManyElectron.matrixElement_Mab(Mp, gauge, omega, nLevel.basis, tp, iLevel.basis, s, grid)
                            # THE PROPAGATOR IS CHOSEN BY THE LEVEL'S OWN KIND, not by testing its energy.  A
                            # retarded Green function is P[1/(E_i+omega-E_n)] - i*pi*delta(E_i+omega-E_n); the
                            # quadrature levels carry the principal value, the on-shell level carries the residue
                            # and takes no denominator at all.  Deciding this by how close enn happens to lie to
                            # the pole would misclassify any quadrature node that drifted near it.
                            # THE DELTA-RESIDUE IS HERE AND THE PRINCIPAL VALUE IS NOT -- MEASURED 30-Sep-2026,
                            # and it is why the second-order sum has no converged value.  A retarded propagator is
                            # 1/(x + i eta) = PV(1/x) - i pi delta(x).  The -i pi below IS the delta piece, taken at
                            # the one designated on-shell node; every OTHER node then divides by the bare x, and a
                            # plain Gauss-Legendre sum of 1/x through a pole is not a principal value.  It has no
                            # limit: refining the mesh puts nodes ever nearer the pole.  Measured on He at 200 eV,
                            # maxKappa = 2, with the SHARING mesh held at 9 and only this one refined, the total runs
                            # 4.35e-05 -> 4.18 -> 4.44 -> 5.15 -> 6.35 -> 7.99 -> 9.64e-05 for 3, 5, 9, 15, 25, 41,
                            # 61 points -- more than doubling, with the increments NOT shrinking (+16, +23, +26,
                            # +21 %).  The OUTER sharing quadrature, by contrast, converges to five digits by 15
                            # points once the two meshes are separated, so it was never the problem.  A correct PV
                            # needs nodes placed symmetrically about the pole, a subtraction of the singular part, or
                            # the interval split AT the pole -- none of which a plain rule does.  Priority item 33.
                            numer[k] = numer[k] + fLevel.mc[r] * Vee * nLevel.mc[t] * nLevel.mc[tp] * OMp *
                                                  iLevel.mc[s]
                        end
                    end
                end
            end
        end
    end
    amplitude = amplitude + PhotoDoubleIonization.assemblePrincipalValue(numer, nLevels, nWeights, onShell,
                                    nBlock, nEnergy, energyRange, eni + omega, doPV)

    if  printout   printstyled("done. \n", color=:light_green)    end

    if  display
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " PhotoDouble^($Mp, absorption) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end

    return( amplitude )
end


"""
`PhotoDoubleIonization.assemblePrincipalValue(numer::Array{ComplexF64,1}, nLevels::Array{Level,1}, nWeights::Array{Float64,1},
        onShell::Array{Bool,1}, nBlock::Union{Nothing,Array{Int64,1}}, nEnergy::Union{Nothing,Array{Float64,1}},
        energyRange::Union{Nothing,Tuple{Float64,Float64}}, ePole::Float64, doPV::Bool)`  
    ... assembles the photon-first second-order sum from the per-level numerators, as a CAUCHY PRINCIPAL VALUE plus the delta residue;
        a value::ComplexF64 is returned.

        WHY A PLAIN SUM CANNOT DO IT.  The retarded propagator is  1/(x + i eta) = PV(1/x) - i pi delta(x), and in the photon-first
        ordering x = E_i + omega - E_n passes through zero at exactly one electron photoionized and the other still bound -- the TS1
        mechanism this module is built on.  Until 05-Oct-2026 the delta residue was taken at one designated node while every other node
        divided by the BARE x, and a Gauss-Legendre sum of 1/x through a pole is not a principal value and has no limit: each refinement
        places nodes nearer the singularity.  Measured on He at 200 eV with maxKappa = 2 and the SHARING mesh held at 9, the total ran
        4.35e-05, 4.18, 4.44, 5.15, 6.35, 7.99, 9.64e-05 for 3, 5, 9, 15, 25, 41, 61 intermediate points -- more than doubling, with the
        increments NOT shrinking (+16, +23, +26, +21 %).

        WHAT IS DONE INSTEAD -- SUBTRACTION OF THE SINGULAR PART, which needs no new mesh and no interval splitting:

            PV int_a^b  N(e) / (e* - e)  de   =   int_a^b  [N(e) - N(e*)] / (e* - e)  de   +   N(e*) ln|(e* - a)/(b - e*)|

        The first integrand is regular at the pole -- it tends to -N'(e*) -- so the existing Gauss-Legendre rule integrates it correctly,
        and the second term is analytic.  The one thing this needs is N(e*), the numerator AT the pole, and that is already computed: it
        is the flagged on-shell level, which the module generates precisely there.  Nothing else about the quadrature changes.

        WHY THE LEVELS MUST BE GROUPED.  There is one pole PER (Green level, kappa) group, at e* = E_i + omega - E_gLevel, while all
        groups share one node set.  Subtracting a group's own N(e*) from another group's nodes would be meaningless, so `nBlock` labels
        the groups and the subtraction is done inside each.

        AND A GROUP WHOSE POLE LIES OUTSIDE THE INTERVAL IS LEFT ALONE.  `generateIntermediateLevels` adds an on-shell level whenever
        e* > 0, which does not mean e* < b: a pole above the top of the intermediate mesh is not enclosed by the integration range, there
        is nothing to regularize, and the residue must NOT be taken either -- that channel is simply not reached by this quadrature.  Such
        a group therefore keeps the plain weighted sum, and its on-shell level is skipped.

        + numer        ::Array{ComplexF64,1}   ... the numerator of each intermediate level, already summed over all CSF indices.
        + nWeights     ::Array{Float64,1}      ... quadrature weight of each level; the on-shell level carries 1 and is not integrated.
        + onShell      ::Array{Bool,1}         ... true for the one level per group that sits AT the pole.
        + nBlock       ::Array{Int64,1}        ... the (Green level, kappa) group of each level, or nothing to fall back.
        + nEnergy      ::Array{Float64,1}      ... continuum electron energy of each level, or nothing to fall back.
        + energyRange  ::Tuple{Float64,Float64}... the interval (a, b) the intermediate quadrature spans, or nothing to fall back.
        + ePole        ::Float64               ... E_i + omega, i.e. the total energy whose crossing defines the pole.
        + doPV         ::Bool                  ... false reproduces the pre-05-Oct-2026 assembly exactly, for an un-updated caller.
"""
function assemblePrincipalValue(numer::Array{ComplexF64,1}, nLevels::Array{Level,1}, nWeights::Array{Float64,1},
                                onShell::Array{Bool,1}, nBlock::Union{Nothing,Array{Int64,1}},
                                nEnergy::Union{Nothing,Array{Float64,1}},
                                energyRange::Union{Nothing,Tuple{Float64,Float64}}, ePole::Float64, doPV::Bool)
    value = ComplexF64(0.)
    if  !doPV
        for  k = 1:length(numer)
            value = value + nWeights[k] * numer[k] *
                    ( onShell[k] ? ComplexF64(0., -pi) : ComplexF64(1. / (ePole - nLevels[k].energy)) )
        end
        return( value )
    end

    aRange = energyRange[1];    bRange = energyRange[2]
    for  block  in  sort(unique(nBlock))
        ks     = findall(isequal(block), nBlock)
        kStars = filter(k -> onShell[k], ks)
        kStar  = isempty(kStars) ? 0 : kStars[1]
        # Is the pole actually enclosed by the integration range?  Only then is there anything to regularize.
        enclosed = kStar > 0  &&  aRange < nEnergy[kStar] < bRange
        if  !enclosed
            for  k  in  ks
                onShell[k]  &&  continue
                value = value + nWeights[k] * numer[k] / (ePole - nLevels[k].energy)
            end
            continue
        end
        nStar = numer[kStar];    eStar = nEnergy[kStar]
        for  k  in  ks
            k == kStar  &&  continue
            value = value + nWeights[k] * (numer[k] - nStar) / (ePole - nLevels[k].energy)
        end
        value = value + nStar * log( (eStar - aRange) / (bRange - eStar) )
        value = value + nStar * ComplexF64(0., -pi)
    end

    return( value )
end


"""
`PhotoDoubleIonization.generateIntermediateLevels(symn::LevelSymmetry, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                            energyGrid::Radial.GridGL, contSettings::Continuum.Settings, maxKappa::Int64;
                                            nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                            primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                            spectators::Union{Nothing,Array{Tuple{Subshell,Float64},1}}=nothing,
                                            onShellTotalEnergy::Union{Nothing,Float64}=nothing)`
    ... builds the intermediate levels of the second-order sum for ONE required total symmetry symn, by adding a continuum partial wave to
        every level of the (N-1)-electron Green multiplet. Every kappa with abs(kappa) <= maxKappa that couples a Green level to symn is
        taken, at every energy of energyGrid.

        These states are not bound: they carry one continuum electron, which is what lets the one-body photon operator connect them to a
        final state carrying two. Because the partial waves are normalized per energy interval, the corresponding quadrature weight of
        energyGrid is returned with each level, and the caller must use it -- the sum over intermediate states is an integral over the
        intermediate electron energy.

        If onShellTotalEnergy is given, ONE further level is added per Green level and kappa, at the energy e* that makes the intermediate
        state degenerate with the initial level plus the photon. That level is the pole of the photon-first ordering and is flagged rather
        than distinguished by its energy, so that the caller multiplies it by -i*pi instead of dividing it by a denominator: a retarded
        Green function is the principal value MINUS i*pi times this residue, and without it the amplitude of an above-threshold process
        comes out exactly real, i.e. with the open channel closed. Levels are added only where e* > 0; a closed channel has no residue.

        A tuple (nLevels, nWeights, nOnShell, nBlock, nEnergy)::Tuple{Array{Level,1}, Array{Float64,1}, Array{Bool,1},
        Array{Int64,1}, Array{Float64,1}} is returned. nBlock labels the (Green level, kappa) group that each intermediate level belongs
        to, and nEnergy carries its continuum electron energy; together they let the caller assemble ONE principal-value integral per
        group, which is what the pole of the photon-first ordering requires -- a group is exactly the set of levels sharing one pole.
"""
function generateIntermediateLevels(symn::LevelSymmetry, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                          energyGrid::Radial.GridGL, contSettings::Continuum.Settings, maxKappa::Int64;
                                          nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                          primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                          spectators::Union{Nothing,Array{Tuple{Subshell,Float64},1}}=nothing,
                                          onShellTotalEnergy::Union{Nothing,Float64}=nothing)
    nLevels = Level[];    nWeights = Float64[];    nOnShell = Bool[]
    nBlock  = Int64[];    nEnergy  = Float64[];    block    = 0

    # A ONE-BODY operator connects this intermediate level to the final state, so the continuum electron that
    # the photon does NOT act upon must be the very same orbital in both -- same kappa, same energy AND the same
    # subshell label, since two continuum subshells that differ only in their label count as different orbitals
    # and make the one-body matrix element vanish.  The caller therefore passes the pair's own subshells.
    if  !isnothing(spectators)
        for  gLevel in gMultiplet.levels
            gLevel = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
            symg   = LevelSymmetry(gLevel.J, gLevel.parity)
            for  (subsh, en) in spectators
                if  !(subsh.kappa in AngularMomentum.allowedKappaSymmetries(symg, symn))    continue    end
                nOrbital, nPhase = Continuum.generateOrbitalForLevel(en, subsh, gLevel, nm, grid, contSettings;
                                                                    nuclearPot=nuclearPot, primitives=primitives)
                block = block + 1
                push!(nLevels,  Basics.generateLevelWithExtraElectron(nOrbital, symn, gLevel))
                push!(nWeights, 1.0);    push!(nOnShell, false)
                push!(nBlock,   block);  push!(nEnergy,  en)
            end
        end
        return( nLevels, nWeights, nOnShell, nBlock, nEnergy )
    end

    for  gLevel in gMultiplet.levels
        # The Green multiplet's basis holds CSFs of every symmetry it was built from; only the CSFs of this
        # level's own symmetry may be coupled to a partial wave, so the basis is reduced first.
        gLevel = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
        symg   = LevelSymmetry(gLevel.J, gLevel.parity)
        for  kappa in AngularMomentum.allowedKappaSymmetries(symg, symn)
            if  abs(kappa) > maxKappa    continue    end
            shn = Subshell(103, kappa);    block = block + 1
            for  (ie, en)  in  enumerate(energyGrid.t)
                nOrbital, nPhase = Continuum.generateOrbitalForLevel(en, shn, gLevel, nm, grid, contSettings;
                                                                    nuclearPot=nuclearPot, primitives=primitives)
                push!(nLevels,  Basics.generateLevelWithExtraElectron(nOrbital, symn, gLevel))
                push!(nWeights, energyGrid.wt[ie]);    push!(nOnShell, false)
                push!(nBlock,   block);                push!(nEnergy,  en)
            end
            # The residue of the retarded Green function: the SAME numerator, with the partial wave generated at the
            # energy that puts the intermediate state on shell, carrying weight 1 because the orbitals are normalized
            # per unit energy and the density of states is therefore unity.
            if  !isnothing(onShellTotalEnergy)
                enStar = onShellTotalEnergy - gLevel.energy
                if  enStar > 0.
                    nOrbital, nPhase = Continuum.generateOrbitalForLevel(enStar, shn, gLevel, nm, grid, contSettings;
                                                                        nuclearPot=nuclearPot, primitives=primitives)
                    push!(nLevels,  Basics.generateLevelWithExtraElectron(nOrbital, symn, gLevel))
                    push!(nWeights, 1.0);    push!(nOnShell, true)
                    push!(nBlock,   block);  push!(nEnergy,  enStar)
                end
            end
        end
    end

    return( nLevels, nWeights, nOnShell, nBlock, nEnergy )
end


"""
`PhotoDoubleIonization.checkGreenSpaceOnGrid(gMultiplet::Multiplet, grid::Radial.Grid)`
    ... tests that the caller-supplied Green space actually lives on the grid this computation integrates on, and returns a
        list::Array{String,1} of the problems found -- empty when there is none.  `settings.gMultiplet` is built by the CALLER,
        normally from its own `Representation`, and a Representation handed a DIFFERENT grid returns orbitals carrying
        `useStandardGrid = true` and an EMPTY stored grid.  Their radial arrays are then read index by index against a mesh they
        do not belong to, so every matrix element through them is nonsense rather than merely inaccurate, and nothing says so.
        MEASURED on `examples/example-Dj.jl`, whose branch 1 passed `Radial.Grid(true)` (392 points to 614 a.u.) while branch 2
        computed on 1302 points to 10 a.u.:  <a|a> came back 0.3581 for 1s_1/2, 0.1149 for 2s_1/2 and 0.1242 / 0.1241 for the two
        2p levels, where it must be 1 -- norms wrong by factors of 3 to 9.  NO choice of transition operator repairs that, which
        is why this check belongs in front of the operator question of item 98 rather than after it.
        Only clearly BOUND orbitals are tested and the tolerance is loose, because a continuum orbital follows a different
        normalisation convention entirely (see `Continuum.checkContinuumAmplitude`) and a bound orbital may legitimately be a
        box-confined pseudostate.  Calibrated by mutation: on the two-grid case every one of the four bound shells is reported,
        and on the same case with one grid nothing is.
"""
function checkGreenSpaceOnGrid(gMultiplet::Multiplet, grid::Radial.Grid)
    problems = String[]
    if  length(gMultiplet.levels) == 0      return( problems )      end
    for  (sh, orb)  in  gMultiplet.levels[1].basis.orbitals
        orb.isBound  ||  continue
        ovl = try   RadialIntegrals.overlap(orb, orb, grid)   catch   NaN   end
        if  isnan(ovl)  ||  abs(ovl - 1.) > 0.05
            push!( problems, "the gMultiplet orbital $(string(sh)) has <a|a> = $(round(ovl, digits=6)) on the grid this " *
                             "computation runs on, where it must be 1.  The Green space was almost certainly built on a " *
                             "DIFFERENT grid:  give the `Representation` that builds `settings.gMultiplet` the SAME grid as " *
                             "the `Atomic.Computation`, or call setDefaults(\"standard grid\", grid) before building it.  " *
                             "Every matrix element through this orbital is otherwise read against the wrong radii." )
        end
    end

    return( problems )
end


"""
`PhotoDoubleIonization.checkSecondOrderConditions(initialMultiplet::Multiplet, gMultiplet::Multiplet, grid::Radial.Grid)`
    ... tests whether the Green space and the initial state are built on ONE orbital set, and returns a list::Array{String,1} of the
        problems found -- empty where they are.  This is the precondition under which the one-orbital junctions of the second-order sum
        in `PhotoDoubleIonization.amplitude` mean what they say, and it is a property of the CALLER's input rather than of this module.

        WHY IT IS CHECKED HERE AND WHY THE SUM NEVERTHELESS USES `ManyElectron.matrixElement_Vee`.  The operator a one-orbital junction
        requires is `H - E`, and `ManyElectron.matrixElement_HminusE` computes it under a stated orthonormality assumption: it drops the
        overlap term `-E <r|s>`, which vanishes identically on one orbital set.  Where the Green space carries its own SCF orbitals that
        term does not vanish, and on this module's own example it is the LARGEST piece of all.  Measured 02-Nov-2026 on
        `examples/example-Dj.jl` -- Green space [1s, 2s, 2p] of Be^3+, initial state 1s^2, each from its own SCF, both on the one grid
        that branch now uses -- at the junction whose Green core is `1s`, omega = 400 eV and equal sharing:

            V_ee  = +2.470622e-01        H - E = +1.560546e-01        one-body part = -9.100758e-02        V_ee/(H-E) = 1.583
            <1s_green|1s_initial> = +0.994787      <2s_green|1s_initial> = +0.081996      worst departure 8.2e-02
            <101s_1/2(eps = 0.5076)|1s_initial> = -0.023524, so  |E_i <r|s>| ~ 3.2e-01 = 2.05 x the whole H - E element

        The one-body part has the right sign and a sensible size, so this is NOT the sign flip measured for
        `TwoElectronOnePhoton`;  it fails the other way.  The continuum orbital is generated in the field of the GREEN level, hence
        orthogonal to the GREEN bound set and not to the initial state's, and that single overlap of -0.0235 times `E_i` is twice the
        matrix element it would correct.  Migrating this module's junctions to `H - E` would therefore ADD a 0.091 correction while
        SILENTLY OMITTING a 0.32 term of the same origin -- a known factor 1.583 traded for an unquantified larger one.  So the sum
        keeps `V_ee`, the ratio is recorded at the junction itself, and the migration waits for a case that satisfies this condition.
        Building one needs the three multiplets on one frozen orbital set, which `Atomic.Computation` + `perform` cannot express today.

        Only clearly BOUND orbitals are compared, since a continuum orbital follows a different normalisation convention entirely;  the
        continuum-bound overlap above is the dominant term but is not reachable once per computation, as those orbitals are generated per
        energy sharing and partial-wave pair.  Calibrated by mutation: the example-Dj case is reported, and the same case with the Green
        multiplet replaced by the initial state's own basis is not.
"""
function checkSecondOrderConditions(initialMultiplet::Multiplet, gMultiplet::Multiplet, grid::Radial.Grid)
    problems = String[]
    if  length(gMultiplet.levels) == 0  ||  length(initialMultiplet.levels) == 0    return( problems )      end
    gBasis = gMultiplet.levels[1].basis;        iBasis = initialMultiplet.levels[1].basis
    worst  = 0.;                                worstName = ""
    for  (sha, orba)  in  gBasis.orbitals,  (shb, orbb)  in  iBasis.orbitals
        if  sha.kappa != shb.kappa                  continue     end
        if  !(orba.isBound  &&  orbb.isBound)       continue     end
        ovl = try   RadialIntegrals.overlap(orba, orbb, grid)   catch   NaN   end
        dev = sha == shb  ?  abs(abs(ovl) - 1.)  :  abs(ovl)
        if  !isnan(dev)  &&  dev > worst
            worst = dev;    worstName = "<$(string(sha))_green | $(string(shb))_initial> = $(round(ovl, digits=6))"
        end
    end
    if  worst > 1.0e-4
        push!( problems, "the Green space and the initial state are not built on one orbital set -- worst departure "                  *
                         "$(round(worst, digits=6)) at $worstName.  The one-orbital junctions of the second-order sum are then "        *
                         "inexact by a factor of order unity, and the overlap term that would correct them is not available:  see "    *
                         "PhotoDoubleIonization.checkSecondOrderConditions, where this is measured.  Build the Green space on the "    *
                         "initial state's own orbitals with Hamiltonian.performCIwithFrozenOrbitals to remove this." )
    end

    return( problems )
end


"""
`PhotoDoubleIonization.computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
        nrContinuum::Int64, settings::PhotoDoubleIonization.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing,
        primitives::Union{Nothing,Bsplines.Primitives}=nothing, printout::Bool=false)`
    ... computes the amplitudes of every partial-wave pair of the given line, at every energy sharing, and from them the energy-differential
        cross section of each sharing. TWO continuum orbitals are generated per pair -- one for each emitted electron, at the energies the
        pair carries -- and they are shared by all channels of that pair, which is what the pair layer is for.

        The two electrons are coupled in sequence: the first to the final ion, giving the intermediate xSymmetry, and the second to that,
        giving the total symmetry of the channel. An electric multipole is evaluated once per gauge into a single EmPropertyC, a magnetic
        one once into an EmPropertyC with equal components, so that no gauge is iterated over here.

        The energy-differential cross section dsigma/depsilon1 of each sharing is the INCOHERENT sum of abs2 over pairs, channels and
        multipoles, times the same prefactor 4 pi^2 alpha omega / (2 (2J_i+1)) that photoionization uses -- the operator and the per-energy
        normalization of the continuum orbitals are the same, and here BOTH outgoing electrons carry it. The total cross section is the
        sharing-integrated differential one TIMES 1/2, because the two electrons are identical and the sharing integral over the full
        range therefore counts every final state twice; `Basics.determineEnergySharings` deliberately supplies no statistics factor and
        says so. A newLine::PhotoDoubleIonization.Line with all amplitudes and both cross sections evaluated is returned.
"""
function  computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                            nrContinuum::Int64, settings::PhotoDoubleIonization.Settings;
                                            nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                            primitives::Union{Nothing,Bsplines.Primitives}=nothing, printout::Bool=false)
    newSharings  = PhotoDoubleIonization.Sharing[]
    contSettings = Continuum.Settings(isnothing(settings.continuumExchange) ? false : settings.continuumExchange, nrContinuum)
    symi         = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    # The sum over intermediate states is an integral over the intermediate electron energy, because those states
    # are normalized per energy interval.  ITS NUMBER OF POINTS IS ITS OWN SETTING -- until 30-Sep-2026 it was taken
    # from `NoEnergySharings`, which tied two unrelated convergence questions together: the outer integral over the
    # energy SHARING, and the completeness of the second-order sum.  Measured consequence, He at 200 eV, maxKappa = 2:
    # dsigma/deps_1 AT EQUAL SHARING grew 7.30e-06 -> 2.00e-05, a factor 2.7, from 3 to 41 sharings, with the eight
    # partial-wave pairs unchanged -- so the differential cross section was not a function of the sharing at all, and
    # a convergence scan in `NoEnergySharings` was really watching the intermediate sum move.  See priority item 33.
    Ji2          = Basics.twice(line.initialLevel.J)
    # THE SAME WRONG PREFACTOR THAT WAS ALREADY FOUND AND CORRECTED IN PhotoIonization.  Until
    # 6-Sep-2026 this read  4 pi^2 alpha omega / (2 (2J_i+1)), which is character for character the
    # expression measured to be wrong at module-PhotoIonization.jl:466-478: there the same amplitudes
    # gave 0.001311 Mb with it and 0.2291 Mb with 8 pi^3/(alpha omega), and only the second is the
    # right order for a K-shell edge.  The two differ by 4 pi / (alpha^2 omega^2), i.e. by two powers
    # of omega -- which is the energy trend this module could not explain -- and by 1/alpha^2 = 18779.
    # The dimensions agree between the two modules although they compute different things: PDI forms a
    # cross section PER UNIT ENERGY, and the extra inverse energy comes from the final state carrying
    # TWO continuum orbitals normalized per unit energy where photoionization carries one.
    csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy / (2*(Ji2 + 1))
    # THE INTERMEDIATE ELECTRON IS ONE ELECTRON, NOT TWO.  Until 6-Sep-2026 this limit was
    # maximum(sh.epsilon1 + sh.epsilon2), i.e. the energy shared by the two FINAL electrons,
    # omega - I_double.  But the intermediate state carries a SINGLE continuum electron, whose
    # energy runs to omega - I_single: for helium the two differ by |E(He+)| = 54.4 eV, and the
    # on-shell point of the second-order sum -- one electron photoionized, the other still bound,
    # which IS the TS1 mechanism this module describes -- lay ABOVE the top of the integration
    # range at every photon energy.  Measured 6-Sep-2026: 176.7 eV against a grid ending at
    # 122.2 eV at omega = 200 eV, and 776.7 against 722.2 at 800 eV, so 45 % of the range was
    # missing at 200 eV and 7.5 % at 800 eV.  The correction is required on its own terms.  It is
    # however NOT numerically large as the code stands -- measured against the committed tree it
    # moves sigma(2+) by +3.2 % at 200 eV, +1.1 % at 400 eV and -0.2 % at 800 eV -- and it does
    # NOT explain the energy trend of the STATUS block, which is a factor of six.  The added
    # region only carries weight once the second-order denominator has the sign that puts the
    # on-shell point inside it; see the note at that denominator.

    # AND THE LIMIT MUST LIE ABOVE THE ON-SHELL POINT, NOT AT IT -- found 05-Oct-2026, and it is why priority
    # item 33's three proposed remedies could not work.  `onShellEnergy` below is exactly the limit this code used
    # until today, so for the LOWEST Green level the pole of the photon-first ordering sat precisely ON the upper
    # endpoint.  An endpoint pole is not a principal value and has no value at all: the integral diverges
    # logarithmically from one side, and no symmetric node placement, no subtraction and no splitting AT the pole
    # can regularize what has no other side.  Measured: subtracting the singular part on the old range left the
    # divergence untouched -- the totals merely shifted by a constant -2.37e-05 (the dropped residue) while the
    # increments stayed at +26, +30, +29, +23 % for 15, 25, 41, 61 points.

    # THE PHYSICAL POINT IS THAT AN INTERMEDIATE STATE IS VIRTUAL.  The second-order sum runs over a COMPLETE set,
    # and the intermediate electron's energy is not restricted by energy conservation -- only the FINAL state's is.
    # Cutting the integral off at the energetically allowed maximum is therefore wrong on its own terms, quite
    # apart from where it puts the pole.  `intermediateEnergyFactor` extends the range beyond it; it must exceed 1,
    # and because it truncates a complete set it is a convergence parameter like `NoIntermediateEnergies`.
    # PRIORITY ITEM 33 WAS CLOSED HERE, 05-Oct-2026.  Should the intermediate sum again fail to converge, the
    # quadrature RULE is not the place to look: subtraction of the singular part was implemented first, on the old
    # range, and left a divergent sequence shifted by a constant.  What matters is that the pole stays INSIDE the
    # interval, i.e. that `intermediateEnergyFactor` exceeds 1 -- and that factor is itself unconverged at its
    # default of 3, which is 1.7 % below its factor-8 value.
    eGreenMin      = minimum(lv.energy for lv in settings.gMultiplet.levels)
    onShellEnergy  = line.initialLevel.energy + line.photonEnergy - eGreenMin
    intFactor      = settings.intermediateEnergyFactor > 1. ? settings.intermediateEnergyFactor : 3.
    maxIntEnergy   = intFactor * onShellEnergy
    noIntEnergies    = settings.NoIntermediateEnergies > 0 ? settings.NoIntermediateEnergies : settings.NoEnergySharings
    intermediateGrid = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.01, maxIntEnergy, noIntEnergies; printout=false)

    for  sharing in line.sharings
        newPairs = PhotoDoubleIonization.PartialWavePair[]
        for  pw in sharing.partialWavePairs
            sh1 = Subshell(101, pw.kappa1);    sh2 = Subshell(102, pw.kappa2)
            # A fresh subshell list: appending to the level's own list would mutate the level that was passed in.
            subshells = copy(line.initialLevel.basis.subshells);    push!(subshells, sh1);    push!(subshells, sh2)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshells)
            redFLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshells[1:end-2])

            # Electron 1 couples to the final ion and gives the intermediate symmetry of the pair.
            cOrbital1, phase1 = Continuum.generateOrbitalForLevel(pw.energy1, sh1, line.finalLevel, nm, grid, contSettings;
                                                                 nuclearPot=nuclearPot, primitives=primitives)
            xLevel            = Basics.generateLevelWithExtraElectron(cOrbital1, pw.xSymmetry, redFLevel)
            # Electron 2 couples to that intermediate symmetry and gives the total symmetry of each channel.
            cOrbital2, phase2 = Continuum.generateOrbitalForLevel(pw.energy2, sh2, line.finalLevel, nm, grid, contSettings;
                                                                 nuclearPot=nuclearPot, primitives=primitives)

            newChannels = PhotoDoubleIonization.Channel[]
            for  ch in pw.channels
                cLevel  = Basics.generateLevelWithExtraElectron(cOrbital2, ch.symmetry, xLevel)
                # The intermediate states of the second-order sum carry ONE continuum electron and are built here,
                # for the two symmetries the two time orderings require: symi for the ordering in which the photon
                # is absorbed last, and the total symmetry of this channel for the ordering in which it is absorbed
                # first.  They come with the quadrature weights of their per-energy normalization.
                # Only the photon-first ordering has a pole: its denominator E_i + omega - E_n crosses zero at the
                # sequential point, one electron ejected and the other still bound.  The other ordering's denominator
                # is E_i - E_n, which stays negative and far from zero, so it is given no on-shell energy.
                nLevelsI, nWeightsI, nOnShellI, nBlockI, nEnergyI = PhotoDoubleIonization.generateIntermediateLevels(symi, settings.gMultiplet, nm, grid,
                                            intermediateGrid, contSettings, settings.maxKappa;
                                            nuclearPot=nuclearPot, primitives=primitives,
                                            spectators=[(sh1, pw.energy1), (sh2, pw.energy2)])
                nLevelsF, nWeightsF, nOnShellF, nBlockF, nEnergyF = PhotoDoubleIonization.generateIntermediateLevels(ch.symmetry, settings.gMultiplet, nm, grid,
                                            intermediateGrid, contSettings, settings.maxKappa;
                                            nuclearPot=nuclearPot, primitives=primitives,
                                            onShellTotalEnergy = line.initialLevel.energy + line.photonEnergy)
                nLevels  = vcat(nLevelsI,  nLevelsF)
                nWeights = vcat(nWeightsI, nWeightsF)
                nOnShell = vcat(nOnShellI, nOnShellF)
                # The two sets are assembled as ONE sum, so the second set's group labels are offset: a group is
                # the set of levels sharing one pole, and a label reused across the two sets would subtract one
                # group's on-shell numerator from the other group's nodes.
                nOffset  = isempty(nBlockI) ? 0 : maximum(nBlockI)
                nBlock   = vcat(nBlockI,  nBlockF .+ nOffset)
                nEnergy  = vcat(nEnergyI, nEnergyF)
                newAmps  = MultipoleAmplitude[]
                for  ma in ch.amplitudes
                    mp = ma.multipole
                    if  string(mp)[1] == 'E'
                        ampC = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Coulomb,   line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout, nOnShell=nOnShell,
                                                               nBlock=nBlock, nEnergy=nEnergy, energyRange=(0.01, maxIntEnergy))
                        ampB = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Babushkin, line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout, nOnShell=nOnShell,
                                                               nBlock=nBlock, nEnergy=nEnergy, energyRange=(0.01, maxIntEnergy))
                        push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampC, ampB)))
                    else
                        ampM = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Magnetic,  line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout, nOnShell=nOnShell,
                                                               nBlock=nBlock, nEnergy=nEnergy, energyRange=(0.01, maxIntEnergy))
                        push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampM, ampM)))
                    end
                end
                push!(newChannels, PhotoDoubleIonization.Channel(ch.symmetry, newAmps))
            end
            push!(newPairs, PhotoDoubleIonization.PartialWavePair(pw.kappa1, pw.energy1, phase1, pw.xSymmetry,
                                                                 pw.kappa2, pw.energy2, phase2, newChannels))
        end

        # Incoherent over pairs, channels and multipoles; abs2 of an EmPropertyC keeps the two gauges apart by itself.
        dcs = EmProperty(0., 0.)
        for  pw in newPairs,  ch in pw.channels,  ma in ch.amplitudes    dcs = dcs + abs2(ma.amplitude)    end
        # The same prefactor as for photoionization: the electron-photon operator and the per-energy
        # normalization of the continuum orbitals are the same, and BOTH outgoing electrons already carry it.
        dcs = csFactor * dcs
        push!(newSharings, PhotoDoubleIonization.Sharing(sharing.omega, sharing.epsilon1, sharing.epsilon2,
                                                         sharing.weight, dcs, newPairs))
    end

    # The total is the sharing-integrated differential cross section, TIMES 1/2 BECAUSE THE TWO ELECTRONS ARE IDENTICAL.
    # dsigma/d eps1 is defined over the whole range [0, omega - I_double], so that integral counts each physical final state
    # twice -- once as (eps1, eps2) and once as (eps2, eps1) -- and the total carries a further factor 1/2.

    # THESE TWO LINES CLAIMED THE OPPOSITE UNTIL 03-Nov-2026, namely that the double counting was "already handled inside
    # Basics.determineEnergySharings and must NOT be applied again".  That routine's own docstring denies it in as many words:
    # "NO PARTICLE-STATISTICS FACTOR IS INCLUDED HERE, DELIBERATELY ... that factor belongs to the PHYSICS of the calling
    # process, not to a quadrature rule, and it must be written at the point where the total is formed so that it can be read
    # there."  The `/2` carried in its weights is the Jacobian of the map [-1,1] -> [0, energy] and would be there for
    # distinguishable shares just the same.  Its 08-Aug-2026 audit of consumers recorded that PhotoDoubleIonization "store[s]
    # and display[s] the weight but never integrate[s] with it", which was true on that date -- and then the 15-Aug rebuild
    # added the integral below, making this module a consumer, while the statistics factor was never supplied.

    # Kornberg & Miraglia state the same 1/2 for their own total.  Measured consequence at omega = 200 eV:
    # sigma(2+)/sigma(+) goes 2.70 % -> 1.35 %, i.e. BELOW the measured 3-4 %, which is the direction a knock-out-only
    # amplitude must err in because it carries no shake-off at all; 2.70 % overshot a measurement this model cannot reach.
    crossSection = EmProperty(0., 0.)
    for  sh in newSharings   crossSection = crossSection + sh.weight * sh.differentialCs   end
    crossSection = 0.5 * crossSection

    newLine = PhotoDoubleIonization.Line(line.initialLevel, line.finalLevel, line.photonEnergy, crossSection, newSharings)

    return( newLine )
end


"""
`PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::PhotoDoubleIonization.Settings; output::Bool=true)`  
    ... to compute the photo-double ionization transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{PhotoDoubleIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::PhotoDoubleIonization.Settings; output::Bool=true)
    println("")
    printstyled("PhotoDoubleIonization.computeLines(): The computation of photo-double ionization properties starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    # The 15-Aug-2026 barrier was removed on 15-Sep-2026 when the maintainer lifted the parking. It is replaced
    # by a WARNING and not by silence: the module is usable for shapes, ratios and trends, and is not usable for
    # a bare cross section. Saying so on every run is cheaper than a reader discovering it from a wrong number.
    printstyled(">>> PhotoDoubleIonization: THIS AMPLITUDE IS KNOCK-OUT (TS1) ONLY -- it carries no shake-off and no\n" *
                ">>> ground-state correlation, so its DOMAIN IS NEAR THRESHOLD, where knock-out dominates. Knock-out\n" *
                ">>> gives a nearly FLAT energy-sharing distribution; shake-off prevails at high energy and alone\n" *
                ">>> sets the asymptotic 1.65-1.67 % ratio, which this module therefore cannot reach. Measured for\n" *
                ">>> He at a 40 a.u. box, 9 sharings, maxKappa = 2: the outermost/midpoint ratio of dsigma/d eps1\n" *
                ">>> is 1.30 at 100 eV, 1.81 at 120 eV, 3.32 at 150 eV and 5.55 at 200 eV, i.e. flat near threshold\n" *
                ">>> and deeply U-shaped above it. THOSE FOUR ARE maxKappa = 2 NUMBERS AND THE U IS TOO DEEP AT IT:\n" *
                ">>> completing the partial-wave sum raises the EQUAL-sharing node by +36 % against +10 % at the\n" *
                ">>> extreme one, so at 200 eV the ratio falls 5.55 -> 4.48 -> 4.42 over maxKappa 2 -> 3 -> 4. Read\n" *
                ">>> shapes and trends here; quote an absolute number only with the box, the sharing count, maxKappa\n" *
                ">>> and NoIntermediateEnergies stated beside it.\n" *
                ">>> The COULOMB (velocity) form is the stable one, and WHICH TRUNCATION DESTABILIZES THE OTHER\n" *
                ">>> was misattributed until 04-Nov-2026. Measured with one variable per run: 9 -> 17 -> 25 inter-\n" *
                ">>> mediate energies move Coulomb -0.87 %, +0.24 % and Babushkin -1.09 %, +0.56 %, and leave the\n" *
                ">>> gauge ratio at 1.840, 1.836, 1.842 -- that mesh is CONVERGED at 9. The near-doubling belongs\n" *
                ">>> to maxKappa: 2 -> 3 -> 4 moves Coulomb +21.6 %, +0.8 % and Babushkin +55.0 %, +2.9 %, walking\n" *
                ">>> the gauge ratio to 2.346 then 2.395. So prefer velocity, converge maxKappa to 3, and do not\n" *
                ">>> read the gauge ratio as a convergence indicator -- Kornberg & Miraglia report a large length-\n" *
                ">>> velocity discrepancy at ALL energies. See the STATUS block in the module docstring.\n",
                color=:yellow)

    # The nuclear potential and the B-spline primitives are constant for the whole computation and are built once.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)

    # ONCE PER COMPUTATION, not once per amplitude.  The first version of this check sat inside `amplitude` and fired 6656
    # times on one run of `examples/example-Dj.jl`, which is noise rather than a guard.
    for  problem in PhotoDoubleIonization.checkGreenSpaceOnGrid(settings.gMultiplet, grid)
        @warn("PhotoDoubleIonization.computeLines(): " * problem)
    end
    for  problem in PhotoDoubleIonization.checkSecondOrderConditions(initialMultiplet, settings.gMultiplet, grid)
        @warn("PhotoDoubleIonization.computeLines(): " * problem)
    end
    lines = PhotoDoubleIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoDoubleIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.
    for  line in lines,  sh in line.sharings   maxEnergy = max(maxEnergy, sh.epsilon1, sh.epsilon2)   end
    minEnergy = maxEnergy
    for  line in lines,  sh in line.sharings   minEnergy = min(minEnergy, sh.epsilon1, sh.epsilon2)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid; minEnergy=minEnergy)
    # Calculate all amplitudes and requested properties
    newLines = PhotoDoubleIonization.Line[]
    for  line in lines
        println("\n>> Calculate photo-double ionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
        newLine = PhotoDoubleIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings;
                                                                    nuclearPot=nuclearPot, primitives=primitives) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoDoubleIonization.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoDoubleIonization.displayResults(iostream, newLines, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoDoubleIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)`  
    ... to determine a list of PhotoDoubleIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{PhotoDoubleIonization.Line,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)
    lines = PhotoDoubleIonization.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                # Add lines for all photon energies
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    excess   = omega_au - (fLevel.energy - iLevel.energy)
                    if  excess < 0.    continue   end  
                    sharings = PhotoDoubleIonization.determineSharingsAndChannels(fLevel, iLevel, omega_au, settings) 
                    push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, omega_au, EmProperty(0., 0.), sharings) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoDoubleIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64,
                                                          settings::PhotoDoubleIonization.Settings)`
    ... determines the energy sharings of the two emitted electrons and, for each of them, the pairs of outgoing partial waves together
        with the total symmetries they serve. The coupling is read off the types: the first electron couples to the final ion to give
        xSymmetry, the second couples to that to give the total symmetry of a channel, and that total symmetry is what a multipole must
        reach from the initial level. Every kappa with abs(kappa) > settings.maxKappa is discarded.

        The distinct pairs (kappa1, xSymmetry, kappa2) are collected FIRST and every multipole that reaches one of their total symmetries
        is then attached, so that a pair which serves several symmetries appears once and its two continuum orbitals are generated once.
        The gauge is not iterated over at all: an electric multipole yields one EmPropertyC holding both gauges, a magnetic one an
        EmPropertyC with equal components.

        An Array{PhotoDoubleIonization.Sharing,1} is returned, with all amplitudes and cross sections still zero.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64,
                                            settings::PhotoDoubleIonization.Settings)
    symi     = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf     = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    shift_au = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    sharings = PhotoDoubleIonization.Sharing[]

    # The excess energy is what remains of the photon energy after the two electrons have been removed.
    eSharings = Basics.determineEnergySharings(omega - (finalLevel.energy - initialLevel.energy), settings.NoEnergySharings)

    for  es in eSharings
        epsilon1 = es[1];    epsilon2 = es[2];    weight = es[3]

        # Collect, for every distinct pair of partial waves, which total symmetry is reached by which multipole.
        pairs = Dict{Tuple{Int64,LevelSymmetry,Int64}, Dict{LevelSymmetry,Array{EmMultipole,1}}}()
        for  mp in settings.multipoles
            for  symt in AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                for  kappa1 = -settings.maxKappa-1:settings.maxKappa
                    if  kappa1 == 0  ||  abs(kappa1) > settings.maxKappa    continue    end
                    for  symx in AngularMomentum.allowedTotalSymmetries(symf, kappa1)
                        for  kappa2 in AngularMomentum.allowedKappaSymmetries(symt, symx)
                            if  abs(kappa2) > settings.maxKappa    continue    end
                            key = (kappa1, symx, kappa2)
                            haskey(pairs, key)         ||  (pairs[key] = Dict{LevelSymmetry,Array{EmMultipole,1}}())
                            haskey(pairs[key], symt)   ||  (pairs[key][symt] = EmMultipole[])
                            mp in pairs[key][symt]     ||  push!(pairs[key][symt], mp)
                        end
                    end
                end
            end
        end

        newPairs = PhotoDoubleIonization.PartialWavePair[]
        for  (key, symtDict) in pairs
            kappa1, symx, kappa2 = key
            newChannels = PhotoDoubleIonization.Channel[]
            for  (symt, mpList) in symtDict
                amplitudes = MultipoleAmplitude[]
                for  mp in mpList   push!(amplitudes, MultipoleAmplitude(mp, EmPropertyC(0.)))    end
                push!(newChannels, PhotoDoubleIonization.Channel(symt, amplitudes))
            end
            push!(newPairs, PhotoDoubleIonization.PartialWavePair(kappa1, epsilon1 + shift_au, 0., symx,
                                                                 kappa2, epsilon2 + shift_au, 0., newChannels))
        end
        push!(sharings, PhotoDoubleIonization.Sharing(omega, epsilon1, epsilon2, weight, EmProperty(0., 0.), newPairs))
    end

    return( sharings )
end


"""
`PhotoDoubleIonization.displayResults(stream::IO, lines::Array{PhotoDoubleIonization.Line,1},`
        `settings::PhotoDoubleIonization.Settings)`
    ... to display the total photo-double ionization cross sections and, for each line, the energy-sharing distribution from which they
        were obtained.  A neat table is printed to stream; nothing::Nothing is returned.

        TWO TABLES, because the two quantities are not equally trustworthy.  The total cross section is the sharing integral
        sigma = SUM_k w_k (dsigma/d eps_1)_k over the Gauss-Legendre points, and it is a single number per line.  The DISTRIBUTION is the
        primary observable of this process and is what a reader must see in order to judge the integral: an integrand that is not smooth
        across the sharing coordinate, or that fails to fall towards the two ends, says that the quadrature is not converged whatever the
        total happens to be.  The second table therefore prints the integrand and its weights rather than only their sum.

        The symmetry of the distribution about equal sharing is a property of the pair set rather than of the physics, and is the internal
        check that costs nothing to look at here: eps_1 <-> eps_2 must give the same differential cross section when the pair set is
        exchange-closed.
"""
function  displayResults(stream::IO, lines::Array{PhotoDoubleIonization.Line,1}, settings::PhotoDoubleIonization.Settings)

    # First table: the sharing-integrated total cross section, one row per line
    nx = 108
    println(stream, " ")
    println(stream, "  Total photo-double ionization cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(12, "omega";      na=4)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(10, "sharings";   na=3);                        sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " *
                                      TableStrings.inUnits("cross section"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))

    for  line in lines
        sa   = "";     isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))          * "    "
        sa = sa * TableStrings.center(10, string(length(line.sharings)); na=3)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb))   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    # Second table: the energy-sharing distribution, i.e. the integrand of the number above
    nx = 118
    println(stream, " ")
    println(stream, "  Energy-sharing distributions   ... the integrand of the totals above; sigma = SUM_k weight_k * dsigma/d eps_1")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(12, "epsilon_1"; na=2)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "epsilon_2"; na=2)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "weight";    na=3);                         sb = sb * TableStrings.hBlank(15)
    sa = sa * TableStrings.center(34, "Cou -- dsigma/d eps_1 -- Bab"; na=2)
    sb = sb * TableStrings.center(34, TableStrings.inUnits("differential cross section") *
                                      "        " * TableStrings.inUnits("differential cross section"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))

    for  line in lines
        sa = TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "  "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "  "
            sb = sb * @sprintf("%.5e", sharing.weight)                                                 * "   "
            # an ENERGY-DIFFERENTIAL cross section needs BOTH conversions.  Until 04-Nov-2026 only the cross-section one
            # was applied while the header read [barn]/[eV], so this column carried barn per HARTREE -- a factor 27.2114,
            # and the epsilon columns beside it were already in eV.  Caught by the quadrature identity, which is unit-blind:
            # the weights of Basics.determineEnergySharings sum to the excess in a.u., so 0.5*sum(w*d) reproduces the total
            # only if d is per Hartree.  `energy-diff. cross section` and TableStrings.inUnits("differential cross section")
            # existed for this from the start and were wired up by no module at all.
            sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy-diff. cross section: from atomic",
                                                            sharing.differentialCs.Coulomb))     * "   "
            sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy-diff. cross section: from atomic",
                                                            sharing.differentialCs.Babushkin))
            println(stream, sb)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotoDoubleIonization.displayLines(stream::IO, lines::Array{PhotoDoubleIonization.Line,1})`
    ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{PhotoDoubleIonization.Line,1})

    # First, print lines and sharings
    nx = 94
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & sharings:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
    sb = sb * TableStrings.flushleft(54, "  i -- f         omega        epsilon_1    epsilon_2"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "

        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))    * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "   "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))


    # Second, print lines and channles
    nx = 130
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & channels:   ... channels are shown only for the first sharing")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(100, "Partial-wave pairs and the total symmetries they serve"; na=5);              
    sb = sb * TableStrings.flushleft(100, "kappa_1 + kappa_2    J^P_x        -->   J^P_total   multipoles"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        if  length(line.sharings) == 0      println(stream, sa * "  (no sharings)");     continue    end
        ic = 0
        for  pair  in  line.sharings[1].partialWavePairs
            for  ch  in  pair.channels
                ic = ic + 1
                if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
                sb = sb * string(Subshell(1, pair.kappa1))[end-4:end] * " + " *
                          string(Subshell(1, pair.kappa2))[end-4:end] * "      "
                sb = sb * TableStrings.center(10, string(pair.xSymmetry); na=3)
                sb = sb * "-->   " * TableStrings.center(10, string(ch.symmetry); na=3)
                sb = sb * join([string(ma.multipole) for ma in ch.amplitudes], ", ")
                println(stream,  sb )
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end

end # module
