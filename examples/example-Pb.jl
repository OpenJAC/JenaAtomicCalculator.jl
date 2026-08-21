
println("Pb) RAYLEIGH scattering of photons at atoms and ions: gamma + |i> --> |i> + gamma, second order in the radiation field.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... first run of a new implementation; see the report below.
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
    #   THE OPEN QUESTION, deliberately left open.  Are JAC's two gauge amplitudes meant to be EQUAL at the same
    #   omega?  PhotoEmission.computeAmplitudesProperties forms its rate as 8 pi alpha omega |amp|^2 / (2J_i+1),
    #   with ONE gauge-independent omega, which would require them to be.  If they are, the fault is local: this
    #   file uses a real-photon amplitude off shell and must supply the missing weighting.  If they are not, the
    #   question belongs to PhotoEmission and to every module forming a product of two such amplitudes --
    #   MultiPhotonTransition among them, which carries six dated branches.  For a real transition omega is pinned
    #   at the transition energy, so PhotoEmission never varies it and this cannot show there; a second-order sum
    #   varies omega freely, which is why it surfaced here first.  ONE RUN WOULD SETTLE IT: compute a real E1 rate
    #   for a known transition in both gauges and see whether they agree.
    #
    #   The branch is NOT dated.  One gauge reproducing a known law is a genuine result and is recorded as such,
    #   but an implementation in which the two gauges obey different power laws is not verified.
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
end
#
setDefaults("print summary: close", "")
