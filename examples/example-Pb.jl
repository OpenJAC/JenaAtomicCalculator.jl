
println("Pb) RAYLEIGH scattering of photons at atoms and ions: gamma + |i> --> |i> + gamma, second order in the radiation field.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... the omega^4 law IS reproduced in both gauges and the basis is exonerated, but the
    #                              ABSOLUTE cross sections rest on an underived prefactor and nothing here constrains a
    #                              PHASE.  A date on a cross-section branch would read as "these cross sections are
    #                              right", which is more than has been shown.  See the report below.
    #
    # Branch a (the omega^4 law): elastic Rayleigh scattering on the ground level of BERYLLIUM-LIKE NEON,
    #   1s^2 2s^2 ^1S_0, at three photon energies well below the 2s -> 2p excitation.
    #
    #   WHY THIS TARGET.  It is the smallest case that exercises the whole second-order machinery honestly.  The
    #   initial level has J = 0 and even parity, so a single E1 photon can reach only J = 1 ODD intermediates --
    #   ONE intermediate symmetry, 1^-.  The configuration 1s^2 2s 2p supplies exactly two levels of that
    #   symmetry, ^3P_1 and ^1P_1.  So gMultiplet holds TWO intermediate levels, both of them namable, and the
    #   sum over intermediate states can be checked by hand if it ever has to be.  That is the whole point of
    #   working with a short, explicit gMultiplet rather than a Green expansion of twenty levels per symmetry.
    #
    #   WHAT IS BEING TESTED, and why it does not depend on the unverified prefactor.  For ELASTIC scattering
    #   omega_out = omega_in, so PhotonScattering.rayleighCrossSection's factor omega_in * omega_out^3 becomes
    #   omega^4 -- the Rayleigh law, the reason the sky is blue.  Far below any resonance the amplitude itself is
    #   nearly energy-independent, so the cross section must rise as omega^4 and a RATIO of two cross sections is
    #   completely insensitive to the constant in front of the formula.  Doubling the photon energy must
    #   therefore multiply the cross section by very nearly 16, and doubling again by 16 once more.
    #
    #   This is the same style of test as the Z^5 scan for one-photon annihilation (example-Of.jl branch b) and
    #   the detailed-balance test for bound-free pair creation (example-Pa.jl branch b): it interrogates the
    #   AMPLITUDE while the normalization is still open, which no absolute number can do.
    #
    #   The deviation from exactly 16 is itself physics rather than error: the amplitude grows slowly as the
    #   photon energy approaches the 2s-2p resonance, so the measured exponent should sit slightly ABOVE 4 and
    #   should grow with energy.  An exponent BELOW 4, or one that moves the wrong way, would point at the
    #   amplitude.
    #
    #   NOT THE OLD MODULE.  This is a fresh implementation in PhotonScattering, not JAC.RayleighCompton, which
    #   was found on 21-Aug-2026 not to compute at all -- its sum over intermediate states runs for the first
    #   level only.  See examples/example-Dg.jl branch a for that diagnosis.
    #
    #   REPORT (21-Aug-2026, first run).  The basis check passed at once -- "the intermediate basis spans all 1
    #   required symmetries" -- confirming the single 1^- expected, and the computation ran clean.
    #
    #       omega [eV]     Coulomb          Babushkin        exponent(Cou)   exponent(Bab)   Cou/Bab
    #          1.0         2.592621e-25     2.231595e-30           -               -         1.162e+05
    #          2.0         4.188359e-24     5.768218e-28         4.014           8.014       7.261e+03
    #          4.0         6.967999e-23     1.535441e-25         4.056           8.056       4.538e+02
    #
    #   WHAT PASSES -- and it is the amplitude, which is what this branch was built to test.  THE COULOMB GAUGE
    #   REPRODUCES THE omega^4 RAYLEIGH LAW, with exponents 4.014 and 4.056.  It passes in the particular way
    #   predicted above BEFORE the run: the exponent sits slightly ABOVE 4 and INCREASES with energy, which is
    #   the amplitude growing as omega approaches the 2s-2p resonance.  The sum over the two intermediate levels,
    #   both time orderings, both energy denominators and the angular coupling therefore all work.
    #
    #   WHAT FAILS.  The Babushkin gauge is exactly FOUR POWERS steeper, 8.014 and 8.056.  This is not a fitted
    #   number: the Coulomb/Babushkin ratio falls by 16.0 per doubling at BOTH steps, i.e. by 2^4 to three
    #   significant figures.  So the two gauge forms differ by exactly omega^4 in the cross section and hence by
    #   exactly omega^2 in the amplitude.
    #
    #   WHERE IT LIVES -- MEASURED, not inferred.  A single vertex was computed for the same levels over
    #   omega = 0.5, 1, 2, 4 eV:
    #
    #       |A_Coulomb|     7.513941e-07  7.514237e-07  7.515424e-07  7.520171e-07     flat to four digits
    #       |A_Babushkin|   1.825198e-07  3.650593e-07  7.301977e-07  1.460712e-06     exactly ~ omega
    #       ratio Cou/Bab   4.116782      2.058361      1.029231      5.148292e-01     falls x0.5000 per doubling
    #
    #   So ONE vertex differs by omega between the gauges, TWO vertices by omega^2, and the cross section by
    #   omega^4 -- the entire discrepancy, with nothing left over.  The structural reason is in
    #   InteractionStrength.multipoleTransition: the length form carries j_L(qr), the velocity form j_L(qr)/(qr)
    #   and j'_L(qr), and for small argument j_1(x) ~ x/3 against j_1(x)/x ~ 1/3.
    #
    #   TWO EXPLANATIONS WERE PROPOSED AND BOTH REFUTED, which is recorded because the discipline is the point.
    #   First: that the factor sat in the Coulomb amplitude -- refuted, the ratio FALLS with omega, so it is
    #   Babushkin that carries it.  Second: that the gauges cross at the on-shell point -- refuted, they cross at
    #   2.058 eV while this system's 1^- transition energies are 10.816 eV and 25.003 eV, matching neither.
    #
    #   THE QUESTION IS ANSWERED, and the answer puts the fault here.  Length-velocity equivalence is an ON-SHELL
    #   IDENTITY: the forms are related by [H,r] = -i p/m, and turning <f|p|i> into omega <f|r|i> uses
    #   (E_f - E_i) = omega.  Off shell that step is FALSE, so the amplitudes have no reason to agree, and the
    #   Bessel asymmetry is exactly how it shows.  Nothing is wrong in InteractionStrength or PhotoEmission.
    #   Two on-shell facts fix the direction: helium 1s^2 -> 1s2p at its true energy gives f_Cou = 0.378 against
    #   f_Bab = 0.421, ratio 0.897, nowhere near omega or 1/omega; and GeneralizedOscillatorStrength, which is a
    #   LENGTH-form calculation by construction, reproduces PhotoExcitation's BABUSHKIN f to 1.000035 as K -> 0.
    #
    #   AND MEASURING THE ON-SHELL RATIO IN *THIS* SYSTEM EXPOSED A SECOND, INDEPENDENT FAULT.  The identity says
    #   the gauges should agree at each level's own transition energy.  They do not:
    #
    #       1^- #1   dE = 10.816 eV   |A_Cou| = 7.560e-07   |A_Bab| = 3.953e-06   ratio 0.191
    #       1^- #2   dE = 25.003 eV   |A_Cou| = 1.059e-03   |A_Bab| = 1.434e-03   ratio 0.738
    #
    #   Neither is near unity, against 0.897 for helium.  The initial 1s^2 2s^2 and intermediate 1s^2 2s 2p
    #   multiplets come from two SEPARATELY CONVERGED SCF runs, so the orbitals are non-orthogonal and no
    #   biorthogonal transformation is applied -- which is what PhotoEmission's calcBiorthogonal exists for.  Note
    #   also that the two levels differ by THREE ORDERS OF MAGNITUDE in amplitude and that it is the WEAK one whose
    #   gauges disagree worst: #1 is an intercombination line, its E1 amplitude surviving only through
    #   singlet-triplet mixing, i.e. a small difference of large numbers.  Strong cancellation is where gauge
    #   agreement fails first, which is why JAC computes a cancellation factor at all.
    #
    #   SO TWO FAULTS, conflated until they were measured apart:
    #     (1) the off-shell omega weighting -- structural, understood, belongs to rayleighAmplitude;
    #     (2) basis quality -- non-orthogonal orbitals, compounded by an intercombination level.
    #
    #   CONSEQUENCE, stated plainly: NO ABSOLUTE NUMBER IN THIS BRANCH SHOULD BE TRUSTED YET.  The omega^4
    #   EXPONENT survives, a scaling exponent being insensitive to overall wave-function quality -- which is
    #   exactly why a ratio test was the right first instrument -- but the magnitudes are not evidence of anything.
    #
    #   NOT YET DONE, and the sharpest test still missing: a HERMITICITY check.  Everything used today --
    #   omega^4 here, Z^5 for annihilation, detailed balance for pair creation -- constrains a MAGNITUDE, and every
    #   one of them would pass with a sign error in one time ordering.  For RAMAN, exchanging initial and final
    #   levels while swapping omega_in and omega_out must relate the amplitude to its own conjugate in a definite
    #   way, with the two orderings swapping into each other.  That constrains a PHASE.  A parallel session lost
    #   an afternoon to an anapole that passed selection rules, reality and plausible magnitudes and was caught
    #   only by exchanging two levels and finding -3 where -1 was required.
    #
    #   FIXED 21-Aug-2026, and the fix is confirmed four ways.  PhotonScattering.offShellFactor rescales each
    #   BABUSHKIN vertex by dE/omega, returning it to the point where the length-velocity identity holds, while
    #   the energy denominator keeps the true photon energy.  Re-running the same scan:
    #
    #       omega [eV]     Coulomb          Babushkin        exp(Cou)   exp(Bab)   Bab/Cou
    #          1.0         2.592621e-25     8.720013e-25         -          -       3.3634
    #          2.0         4.188359e-24     1.408711e-23       4.014      4.014     3.3634
    #          4.0         6.967999e-23     2.343618e-22       4.056      4.056     3.3634
    #
    #     (1) BOTH gauges now obey the omega^4 law, with identical exponents.
    #     (2) Coulomb is BIT-IDENTICAL to the pre-fix run -- checked, not assumed -- so the correction touched
    #         only Babushkin, which is what offShellFactor was written to do.
    #     (3) The gauge ratio is CONSTANT to five digits, varying by a factor 1.0000 across the scan where it
    #         previously varied by 256.  The power-law difference is gone entirely.
    #     (4) The check that was NOT designed in, and the most convincing: the residual 3.3634 sits BETWEEN the
    #         two independently measured on-shell ratios, 1/0.191 = 5.228 for 1^- #1 and 1/0.738 = 1.354 for
    #         1^- #2 -- exactly where a weighted combination of the two levels' contributions belongs.  Those
    #         numbers came from a separate calculation made for another purpose, so the agreement is not
    #         something the fix could have manufactured.
    #
    #   The gauge comparison has thereby become what a gauge check is meant to be: a measure of WAVE-FUNCTION
    #   QUALITY.  The residual 3.36 is correlation-limited -- single-configuration Dirac-Fock, a 2x2 CI, and one
    #   of the two intermediate levels an intercombination line -- and is NOT a defect.
    #
    #   THE BASIS WAS EXONERATED ALONG THE WAY, by a three-point comparison worth keeping:
    #
    #       separate DFS         ratios 0.191 / 0.738      good states, DIFFERENT one-body Hamiltonians
    #       common mean field    ratios 0.168 / 0.726      good states, SHARED Hamiltonian
    #       common nuclear field ratios 0.0155 / 0.0293    BAD states,  shared Hamiltonian
    #
    #   Sharing the Hamiltonian changes essentially nothing; destroying state quality changes everything.  So the
    #   on-shell disagreement is correlation, not the non-orthogonality it was first blamed on, and no
    #   biorthogonal transformation is called for.  Both runs used scField = DFSField() -- checked, since a
    #   different default would have varied two things at once.  Note that Basics.NuclearField() is a DIAGNOSTIC
    #   ONLY and a CONFOUNDED one: it shares the Hamiltonian but also removes all screening, so it varies state
    #   quality too, which is why it made the ratio worse rather than better.
    #
    #   STILL NOT DATED, and the reason is now narrow.  The SHAPE is verified in both gauges; the ABSOLUTE
    #   MAGNITUDE is not, rayleighCrossSection's prefactor never having been derived.  And no test used here
    #   constrains a PHASE: a sign error in one time ordering would survive every check above.  The Raman
    #   Hermiticity check -- exchange initial and final levels, swap omega_in and omega_out, require the amplitude
    #   to relate to its own conjugate with the two orderings swapping -- is the missing instrument.
    #
    grid       = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm         = Nuclear.Model(10.)

    # The intermediate levels: 1s^2 2s 2p supplies ^3P_1 and ^1P_1, the only two levels of symmetry 1^- that a
    # single E1 photon can reach from 1s^2 2s^2 ^1S_0.
    gMultiplet = SelfConsistent.performSCF([Configuration("1s^2 2s 2p")], nm, grid, AsfSettings())

    settings   = PhotonScattering.Settings(PhotonScattering.Settings();
                        process        = PhotonScattering.RayleighScattering(),
                        approximation  = PhotonScattering.SecondOrderGreen(),
                        photonEnergies = [1.0, 2.0, 4.0],           # eV, far below the 2s -> 2p excitation
                        multipoles     = [E1],
                        gMultiplet     = gMultiplet,
                        selfTolerance  = 1.0e-6,
                        printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(), name="Rayleigh scattering on Be-like neon", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... not yet run.
    #
    # Branch b (the guard, deliberately provoked): the same computation with an intermediate basis that does NOT
    #   span the required symmetry, to check that PhotonScattering.checkIntermediateBasis says so BEFORE any
    #   amplitude is evaluated.
    #
    #   WHY THIS DESERVES A BRANCH OF ITS OWN.  A second-order amplitude whose intermediate basis lacks the
    #   needed symmetry evaluates to exactly zero, and a column of exact zeros reads like a selection rule.  Two
    #   separate pieces of work were caught by precisely this on 21-Aug-2026: the old RayleighCompton dies with
    #   `error("stop a: Green channel not found for symmetry ...")` from deep inside its amplitude loop, after
    #   the whole channel table has already printed; and in a parallel session an anapole amplitude returned
    #   identically 0.0 for sixteen level pairs, which was an inverted parity template rather than a rule.
    #
    #   Here the gMultiplet is built from 1s^2 2s 3s, which carries NO level of symmetry 1^-, so every channel of
    #   the computation is starved.  The expected behaviour is a @warn naming the missing symmetry and stating
    #   that the resulting zeros are an artefact of the basis and not a selection rule -- and then a table of
    #   zeros that the reader has been told how to interpret.
    #
    grid       = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm         = Nuclear.Model(10.)
    gMultiplet = SelfConsistent.performSCF([Configuration("1s^2 2s 3s")], nm, grid, AsfSettings())

    settings   = PhotonScattering.Settings(PhotonScattering.Settings();
                        process        = PhotonScattering.RayleighScattering(),
                        approximation  = PhotonScattering.SecondOrderGreen(),
                        photonEnergies = [2.0], multipoles = [E1],
                        gMultiplet     = gMultiplet, selfTolerance = 1.0e-6, printBefore = true )

    wa = Atomic.Computation(Atomic.Computation(), name="Rayleigh with a starved intermediate basis", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      22-Aug-2026
    # Last successful:  unknown ... see the report below.
    #
    # Branch c (HERMITICITY -- a check that constrains a PHASE): the Raman amplitude for i --> f is compared with
    #   the amplitude for the REVERSE process f --> i with the two photon energies exchanged.
    #
    #   WHY THIS BRANCH EXISTS.  Every other test in this enterprise constrains a MAGNITUDE -- the omega^4 law
    #   above, the Z^5 scan of example-Of.jl, the detailed-balance test of example-Pa.jl, the Lorentzian of
    #   example-Pc.jl.  All four would pass unchanged with a sign error in one of the two time orderings.  A
    #   parallel session lost an afternoon to exactly that: an anapole amplitude that passed selection rules,
    #   reality and plausible magnitudes, and was caught only by exchanging two levels and finding -3 where -1
    #   was required.  A check that constrains a phase is worth more than one that constrains a magnitude.
    #
    #   THE RELATION, derived rather than asserted.  Write the reverse process with incoming omega_out and
    #   outgoing omega_in; energy conservation makes that consistent, since omega_in - omega_out = E_f - E_i.
    #   Both denominators are then UNCHANGED,
    #
    #       E_f + omega_out - E_nu  =  E_i + omega_in  - E_nu       (AbsorbThenEmit)
    #       E_f - omega_in  - E_nu  =  E_i - omega_out - E_nu       (EmitThenAbsorb)
    #
    #   while the numerators conjugate, JAC defining Absorption as the conjugate of Emission with the levels
    #   exchanged.  Hence, term by term,
    #
    #       A(f --> i; omega_out, omega_in)  =  conj( A(i --> f; omega_in, omega_out) )
    #
    #   and EACH ORDERING MAPS TO ITSELF.  Earlier notes in this repository said the two orderings "swap into
    #   each other" under the exchange; that was wrong and is corrected here.
    #
    #   WHAT IT CATCHES, and what it does not -- stated because a test whose reach is unknown is worth little.
    #     CATCHES: a denominator built from the wrong photon energy (omega_in where omega_out belongs, or the
    #       reverse); an inconsistent operator ordering between the two time-ordering branches; a conjugation
    #       applied in one place and not the other.  These are the errors that produce a plausible magnitude.
    #     DOES NOT CATCH: a global sign on ONE ordering applied consistently in both directions -- it conjugates
    #       along with everything else.  For that, the second check below is the instrument.
    #
    #   THE SECOND CHECK, nearly free and sharper for exactly the case the first one misses.  As omega --> 0 both
    #   denominators tend to the SAME value E_i - E_nu, so the two orderings must ADD.  If one carried the wrong
    #   sign they would cancel instead, and the low-frequency amplitude would nearly vanish rather than tending
    #   to twice a single ordering.  So: |T1 + T2| must approach |T1| + |T2|, not |T1| - |T2|.
    #
    #   REPORT (22-Aug-2026).
    #
    #   TEST 1 PASSED, AND THE PASS IS WORTHLESS.  Every channel returned |rev - conj(fwd)|/|fwd| = 0.0 exactly.
    #   It also returned |Im|/|A| = 0.0: the amplitudes are PURELY REAL.  Bound-bound E1 matrix elements built
    #   from real orbitals, divided by real non-resonant denominators, give a real amplitude -- so conj is the
    #   IDENTITY and the relation A(f->i) = conj(A(i->f)) degenerates into A = A.  A perfect score on a tautology.
    #   Recorded rather than deleted because seeing 0.0 and writing "Hermiticity verified" was one step away, and
    #   a test that cannot fail is worse than no test: it manufactures confidence.
    #
    #   THE OBVIOUS REPAIR DOES NOT WORK EITHER, which is why this branch keeps two tests instead of one.  Running
    #   the exchange resonantly WOULD make the amplitude complex -- but the width breaks the relation by
    #   construction: the reverse denominator is delta + i*Gamma/2, not its conjugate, because a decaying state
    #   genuinely violates time reversal.  So the exchange test is VACUOUS where it applies and INAPPLICABLE where
    #   it would bite.  It is kept only to record that, and to print |Im|/|A| so the emptiness is visible.
    #
    #   TEST 2 PASSED, AND THIS ONE HAS TEETH.  At omega = 0.20 eV, elastic:
    #
    #       Coulomb     T1 = -1.2941e-06   T2 = -1.2758e-06   |T1+T2| = 2.5700e-06   ratio to constructive 1.0
    #       Babushkin   T1 = -3.1623e-06   T2 = -3.1176e-06   |T1+T2| = 6.2798e-06   ratio to constructive 1.0
    #
    #   The two orderings carry the SAME sign and add.  A wrong relative sign would have given 1.83e-08 instead of
    #   2.57e-06 -- a factor of 140, so this discriminates by two orders of magnitude rather than by a percent.
    #
    #   A SECOND CONFIRMATION NOT DESIGNED IN: T1 and T2 agree to 1.4 %.  At low omega both denominators tend to
    #   the same E_i - E_nu, so the orderings must become nearly EQUAL, not merely same-signed.  They do.  That is
    #   the low-frequency limit behaving correctly for a reason independent of the additivity itself.
    #
    #   WHAT IS THEREFORE NOW ESTABLISHED, and it is the one thing every magnitude test missed: the RELATIVE SIGN
    #   between the two time orderings is right.  What is still NOT established is any overall phase convention,
    #   which for a real amplitude is not a meaningful question, and the absolute normalization, which remains
    #   underived here as everywhere else in this module.
    #
    grid  = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm    = Nuclear.Model(10.)
    # ONE multiplet for initial and final alike, so both are on a single energy scale.  Mixing two separately
    # converged SCF runs corrupts E_f - E_i and with it omega_out; that error cost a run in example-Pc.jl.
    both  = [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]
    mult  = SelfConsistent.performSCF(both, nm, grid, AsfSettings())
    gMult = SelfConsistent.performSCF([Configuration("1s^2 2s 2p")], nm, grid, AsfSettings())

    settings = PhotonScattering.Settings(PhotonScattering.Settings();
                    process = PhotonScattering.RamanScattering(), approximation = PhotonScattering.SecondOrderGreen(),
                    multipoles = [E1], gMultiplet = gMult, selfTolerance = 1.0e-6 )

    iLevel = mult.levels[1]
    fLevel = mult.levels[ findfirst(l -> l.parity == Basics.plus && l.index != 1, mult.levels) ]
    # omega_in is set FROM the level splitting so that omega_out is a fixed 10 eV and the channel is always OPEN.
    # Choosing omega_in blindly closes it: the first even-parity excited level here lies 36 eV up, so omega_in = 5 eV
    # gave omega_out = -31 eV and the run died in a Bessel function. A Raman channel needs omega_in > E_f - E_i.
    outOm  = Defaults.convertUnits("energy: to atomic", 10.0)
    inOm   = outOm + (fLevel.energy - iLevel.energy)
    println("\n  Raman  $(iLevel.index) --> $(fLevel.index),  omega_in = " *
            string(Defaults.convertUnits("energy: from atomic", inOm)) * " eV,  omega_out = 10.0 eV")

    println("\n  TEST 1 -- exchange:   A(f->i; w_out,w_in)  =?  conj( A(i->f; w_in,w_out) )")
    for ch in PhotonScattering.determineSecondOrderChannels(fLevel, iLevel, settings)
        fwd = PhotonScattering.rayleighAmplitude(ch, fLevel, iLevel, inOm, outOm, grid, settings; printout=false)
        rev = PhotonScattering.rayleighAmplitude(ch, iLevel, fLevel, outOm, inOm, grid, settings; printout=false)
        if  abs(fwd) > 1.0e-30
            println("    $(ch.gauge) $(ch.timeOrdering):  fwd = $fwd   |Im|/|A| = " * string(abs(imag(fwd))/abs(fwd)) *
                    "   |rev - conj(fwd)|/|fwd| = " * string(abs(rev - conj(fwd))/abs(fwd)))
        end
    end

    # TEST 2 -- the one with teeth.  As omega -> 0 both denominators tend to the SAME value E_i - E_nu, so the two
    # time orderings must ADD.  A wrong relative sign makes them cancel instead.  Run at a SMALL elastic omega, where
    # the two denominators are equal in sign; at the Raman energies above they legitimately differ in sign.
    smallOm = Defaults.convertUnits("energy: to atomic", 0.20)
    println("\n  TEST 2 -- low-frequency additivity at omega = 0.20 eV, elastic (level 1 -> level 1):")
    for gauge in [Basics.Coulomb, Basics.Babushkin]
        t1 = ComplexF64(0.);   t2 = ComplexF64(0.)
        for ch in PhotonScattering.determineSecondOrderChannels(iLevel, iLevel, settings)
            if  ch.gauge != gauge    continue    end
            a = PhotonScattering.rayleighAmplitude(ch, iLevel, iLevel, smallOm, smallOm, grid, settings; printout=false)
            if  ch.timeOrdering == PhotonScattering.AbsorbThenEmit()   t1 = t1 + a   else   t2 = t2 + a   end
        end
        constructive = abs(t1) + abs(t2);   destructive = abs(abs(t1) - abs(t2))
        println("    $gauge:  T1 = $t1   T2 = $t2")
        println("        |T1+T2| = " * string(abs(t1+t2)) * "   constructive would be " * string(constructive) *
                "   destructive would be " * string(destructive))
        println("        ratio to constructive = " * string(abs(t1+t2)/constructive) * "   (must be ~1)")
    end
    #
end
#
setDefaults("print summary: close", "")
