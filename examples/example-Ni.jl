#
println("Ni) Tests of the empirical Belyaev-Yakovleva (2017) simplified model for inelastic ion/atom-hydrogen collisions.")
setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "cm^2")
setDefaults("unit: rate",          "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  true
    #
    # Last successful:  23-Jul-2026
    # Branch 1: low-level building blocks (Empirical.nonadiabaticRadius, landauZenerProbability) and
    #   Empirical.neutralizationReducedRate, validated against the worked example of Belyaev & Yakovleva,
    #   A&A 608, A33 (2017): Ba2+ + H^- -> Ba+(j) + H neutralization, Table 1 (channels) and Fig. 4/Table 2
    #   (reduced rate coefficients) at T = 6000 K.
    # System: Ba2+ + H^- (ionic) vs. Ba+(nl) + H (19 covalent channels from Table 1 of the paper).
    # Checks:
    #   - N_if(T=6000K, Ef=-4.29573 eV) [[Ba+(6d), the paper's largest-rate channel]] = 7.588e-8 cm^3/s here vs.
    #     the paper's quoted 7.59e-8 cm^3/s -- agreement to 3 significant figures.
    #   - The full reduced-rate curve peaks near Ef ~ -4 to -4.4 eV at ~7.6e-8 cm^3/s and falls off by several
    #     orders of magnitude toward both shallower and deeper bound energies, exactly the shape of the paper's
    #     Fig. 4 ("most optimal window" [-2.8,-8.3] eV, paper's Sec. 2.2).
    #
    println("\n  Empirical.nonadiabaticRadius / landauZenerProbability / neutralizationReducedRate for Ba2+ + H^-:\n")
    Z  = 1.0                                            # Ba+ covalent charge (Ba2+ ionic charge = Z+1 = 2)
    MBa = 137.327 / Defaults.ELECTRON_MASS_U;   MH = 1.00794 / Defaults.ELECTRON_MASS_U
    mu  = MBa * MH / (MBa + MH)
    EHminus = Empirical.hydrogenAnionEnergy()
    T   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)

    # Table 1 of Belyaev & Yakovleva (2017): (name, bound energy [eV], statistical probability)
    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    channels = [ Empirical.InelasticHChannel("Ba+(6s)",  eV(-10.0080), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(5d)",  eV(-9.34416), 1/20),
                 Empirical.InelasticHChannel("Ba+(6p)",  eV(-7.35615), 1/12),
                 Empirical.InelasticHChannel("Ba+(7s)",  eV(-4.75665), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(6d)",  eV(-4.29573), 1/20),
                 Empirical.InelasticHChannel("Ba+(4f)",  eV(-4.00879), 1/28),
                 Empirical.InelasticHChannel("Ba+(7p)",  eV(-3.83309), 1/12),
                 Empirical.InelasticHChannel("Ba+(5f)",  eV(-2.87539), 1/28),
                 Empirical.InelasticHChannel("Ba+(8s)",  eV(-2.81381), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(7d)",  eV(-2.58669), 1/20),
                 Empirical.InelasticHChannel("Ba+(8p)",  eV(-2.37789), 1/12),
                 Empirical.InelasticHChannel("Ba+(5g)",  eV(-2.19370), 1/36),
                 Empirical.InelasticHChannel("Ba+(6f)",  eV(-1.99196), 1/28),
                 Empirical.InelasticHChannel("Ba+(9s)",  eV(-1.86429), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(8d)",  eV(-1.73768), 1/20),
                 Empirical.InelasticHChannel("Ba+(9p)",  eV(-1.62363), 1/12),
                 Empirical.InelasticHChannel("Ba+(6g)",  eV(-1.52427), 1/36),
                 Empirical.InelasticHChannel("Ba+(7f)",  eV(-1.42341), 1/28),
                 Empirical.InelasticHChannel("Ba+(10s)", eV(-1.32732), 1/4 ) ]

    N6d = Empirical.neutralizationReducedRate(T, channels[5].E, EHminus, Z, mu; printout=true)
    factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    println("  N_if(Ba+(6d)) = $(round(factor*N6d, sigdigits=4)) cm^3/s   [paper: 7.59e-8 cm^3/s]")

    println("\n  Full neutralization reduced-rate curve at T = 6000 K:")
    for  ch  in  channels
        N = Empirical.neutralizationReducedRate(T, ch.E, EHminus, Z, mu)
        println("    $(ch.name):  Ef = $(round(Defaults.convertUnits("energy: from atomic to eV", ch.E), digits=3)) eV" *
                "   N_if = $(round(factor*N, sigdigits=4)) cm^3/s")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  unknown -- deliberately left undated, see below.
    # Branch 2: Empirical.inelasticHCollisionRateMatrix -- the full K_if(T) rate-coefficient matrix (neutralization,
    #   ion-pair formation, excitation, de-excitation) for the same Ba2+/Ba+ + H(H^-) system, at several temperatures,
    #   compared to Table 2 of Belyaev & Yakovleva (2017).
    # Checks:
    #   - K(ionic -> Ba+(6d)) matches N_if from branch 1 exactly (pstat_ionic = 1 here, cf. paper's Step 4).
    #   - K(ionic -> Ba+(6d)) agrees with the paper's Table 2 value (7.59e-08 cm^3/s at T=6000K) to 3 significant figures.
    #   - OPEN, PAUSED ISSUE: de-excitation/excitation entries (e.g. Ba+(7p) -> Ba+(6d)) are only order-of-magnitude/
    #     qualitative-trend estimates -- checked against the paper's Table 2 value of 8.94e-10 cm^3/s (T=6000K), this
    #     gives ~4.6e-10 cm^3/s, a factor ~1.9 too small, growing to a factor ~2.6 for larger energy defects (checked
    #     across all 6 de-excitation pairs from the Ba+(7p) row of Table 2). The likely root cause was identified via
    #     Belyaev, Phys. Rev. A 48, 4299 (1993) [examples/papers/a93.pra-belyaev-charge-exchange.pdf]: this module's
    #     de-excitation formula treats the process as an isolated 3-state (i, ionic, f) problem, ignoring that the
    #     ionic curve also crosses the OTHER ~17 covalent channels along the way -- exactly the omission that
    #     paper's Eq. (3.8)/multichannel discussion warns can cause a factor-2 (or larger) error. A first attempt to
    #     apply that paper's general multichannel formula directly (treating the ionic curve as entrance, excluding
    #     state i, applying Eq. 3.8 to the other 18 channels) made the result *worse* (2 orders of magnitude too
    #     small), so the correct mapping between the two papers' formalisms is not yet understood -- deliberately
    #     PAUSED here rather than guessed at further; see project_inelastic_h_collisions.md for the full trail.
    #     Do not rely on the de-excitation/excitation numbers quantitatively; neutralization/ion-pair-formation are
    #     well validated and safe to use.
    #
    println("\n  Empirical.inelasticHCollisionRateMatrix for Ba2+ + H^- / Ba+ + H:\n")
    Z  = 1.0
    MBa = 137.327 / Defaults.ELECTRON_MASS_U;   MH = 1.00794 / Defaults.ELECTRON_MASS_U
    mu  = MBa * MH / (MBa + MH)
    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    channels = [ Empirical.InelasticHChannel("Ba+(6s)",  eV(-10.0080), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(5d)",  eV(-9.34416), 1/20),
                 Empirical.InelasticHChannel("Ba+(6p)",  eV(-7.35615), 1/12),
                 Empirical.InelasticHChannel("Ba+(7s)",  eV(-4.75665), 1/4 ),
                 Empirical.InelasticHChannel("Ba+(6d)",  eV(-4.29573), 1/20),
                 Empirical.InelasticHChannel("Ba+(4f)",  eV(-4.00879), 1/28),
                 Empirical.InelasticHChannel("Ba+(7p)",  eV(-3.83309), 1/12) ]     # a 7-channel subset for speed

    for  TK  in  [2000., 6000., 10000.]
        local T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", TK)
        result = Empirical.inelasticHCollisionRateMatrix(T, channels, 1.0, Z, mu; printout=true)
    end
    #
    setDefaults("print summary: close", "")
    #
end
