#
println("Nh) Tests of empirical single-electron-capture (charge-exchange) cross sections and plasma rate coefficients.")
setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "cm^2")
setDefaults("unit: rate",          "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  17-Jul-2026
    # Branch 1: Empirical.chargeExchangeCrossSection(q, Ip, approx) -- the direct (q, Ip) method, comparing the two
    #   slow-collision approximations for the classic He2+ + H(1s) -> He+ + H+ benchmark system.
    # Systems: He2+ (q=2) capturing from atomic hydrogen (Ip = 13.6 eV).
    # Checks:
    #   - OverBarrierModel1980 and NiehausScaling1986 agree to within a factor of a few (7.9e-16 vs 2.8e-15 cm^2),
    #     both squarely in the ~1e-16 ... 1e-14 cm^2 range of measured slow-collision CX cross sections for this
    #     and similar systems (e.g. Fritsch & Lin, Shah & Gilbody-type experiments); *not* an exact benchmark, since
    #     both formulas are order-of-magnitude/global-fit estimates, good to a factor of a few at best.
    #   - A collision velocity of v = 1.5 a.u. (well above the adiabatic regime) triggers the expected warning.
    #
    println("\n  Empirical.chargeExchangeCrossSection(q, Ip, approx) for He2+ + H(1s):\n")
    setDefaults("nuclear: charge", 1.)
    IpH = Empirical.ionizationPotential(1, Configuration("1s^1"))
    sigOBM  = Empirical.chargeExchangeCrossSection(2.0, IpH, Empirical.OverBarrierModel1980();  printout=true)
    sigNieh = Empirical.chargeExchangeCrossSection(2.0, IpH, Empirical.NiehausScaling1986();     printout=true)
    println("  ratio Niehaus/OBM = $(round(sigNieh/sigOBM, digits=2)) (both estimates only, not an exact benchmark)")
    println("\n  Velocity regime warning (v = 1.5 a.u., well above the adiabatic limit):")
    Empirical.chargeExchangeCrossSection(2.0, IpH, Empirical.OverBarrierModel1980(); velocity=1.5)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 2: Empirical.chargeExchangeCrossSection(q, targetConf, targetZ, approx) -- the Configuration-based
    #   convenience wrapper, demonstrating that atom-ion and ion-ion charge exchange use the *same* function, just
    #   with the donor's own (possibly already ionized) Configuration and Z.
    # Systems: (a) atom-ion: C6+ (q=6) capturing from neutral He (Ip = 24.6 eV);
    #          (b) ion-ion: C6+ (q=6) capturing from He+ (Ip = 54.4 eV) -- a much more tightly bound donor.
    # Checks:
    #   - The wrapper reproduces the direct (q, Ip) method to machine precision (verified in testModule_Empirical).
    #   - The ion-ion cross section (He+ donor) is markedly smaller than the atom-ion one (He donor), since the
    #     more tightly bound donor electron requires the projectile to approach much closer before the barrier
    #     drops -- a smaller critical radius R_c and, for NiehausScaling1986, the explicit 1/Ip^2 suppression.
    #
    println("\n  Empirical.chargeExchangeCrossSection(q, targetConf, targetZ, approx): atom-ion vs ion-ion CX:\n")
    setDefaults("nuclear: charge", 2.)
    sigAtomIon = Empirical.chargeExchangeCrossSection(6.0, Configuration("1s^2"), 2.0, Empirical.NiehausScaling1986(); printout=true)
    sigIonIon  = Empirical.chargeExchangeCrossSection(6.0, Configuration("1s^1"), 2.0, Empirical.NiehausScaling1986(); printout=true)
    sigAtomIonCm2 = Defaults.convertUnits("cross section: from atomic to cm^2", sigAtomIon)
    sigIonIonCm2  = Defaults.convertUnits("cross section: from atomic to cm^2", sigIonIon)
    println("  C6+ + He  (atom-ion) sigma = $(round(sigAtomIonCm2, sigdigits=4)) cm^2")
    println("  C6+ + He+ (ion-ion)  sigma = $(round(sigIonIonCm2,  sigdigits=4)) cm^2")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 3: Empirical.chargeExchangePlasmaAlpha(T, q, Ip, Mproj, Mtarget, approx) -- the Maxwellian relative-
    #   velocity-folded plasma rate coefficient, analytic (alpha = <v> sigma) since sigma is velocity-independent
    #   within the slow-collision regime -- no Gauss-Legendre quadrature is needed, unlike the electron-impact
    #   plasma rate coefficients elsewhere in this module.
    # System: alpha particles (He2+, q=2, M ~ 4 amu) charge-exchanging with atomic hydrogen (M ~ 1 amu) at several
    #   plasma temperatures -- a fusion-edge-plasma-relevant system.
    # Checks:
    #   - alpha^(CX) grows as sqrt(T) (via the Maxwellian mean speed), since sigma itself does not depend on T here.
    #   - At T = 10 keV the thermal relative speed approaches the adiabatic limit (v ~ 0.5-1 a.u.) and the expected
    #     regime warning fires, flagging that the constant-sigma approximation is no longer trustworthy there.
    #
    println("\n  Empirical.chargeExchangePlasmaAlpha(T, q, Ip, Mproj, Mtarget, approx) for He2+ + H:\n")
    setDefaults("nuclear: charge", 1.)
    IpH   = Empirical.ionizationPotential(1, Configuration("1s^1"))
    Mproj = Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U * 4.0    ##  alpha particle, M ~ 4 amu
    Mtarg = Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U * 1.0    ##  atomic hydrogen, M ~ 1 amu
    println("    T_e [eV]     alpha^(CX) [cm^3/s]")
    for  TeV  in  [0.1, 1.0, 100.0, 1000.0]
        Tau = Defaults.convertUnits("energy: from eV to atomic", TeV)
        alp = Empirical.chargeExchangePlasmaAlpha(Tau, 2.0, IpH, Mproj, Mtarg, Empirical.NiehausScaling1986())
        alpx = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0) * alp
        println("    $TeV          $(round(alpx, sigdigits=4))")
    end
    println("\n  Regime warning at T = 10 keV (thermal speed approaches the adiabatic limit):")
    T10k = Defaults.convertUnits("energy: from eV to atomic", 1.0e4)
    Empirical.chargeExchangePlasmaAlpha(T10k, 2.0, IpH, Mproj, Mtarg, Empirical.NiehausScaling1986(); printout=true)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 4: Empirical.chargeExchangeCaptureShell/StateSelectiveCrossSection -- the state-selective (n,l)
    #   extension of the total capture cross section, via the standard level-matching estimate n_c = q/sqrt(2 Ip)
    #   and a statistical (2l+1)/n^2 split within the single dominant shell.
    # System: C6+ + H(1s) -- a classic example where q and Ip are matched so that capture goes predominantly into
    #   a well-known, textbook shell (n = 6).
    # Checks:
    #   - n_c = 6.0 to 3 digits (q=6, Ip=13.6 eV), and the rounded shell n = 6 accordingly.
    #   - The l-resolved fractions sum to 1 and increase monotonically with l (the (2l+1)/n^2 statistical weight).
    #   - The Configuration-based wrapper reproduces the direct (q, Ip) call.
    #   - A weakly-matched case (q=1, Ip=500 eV, i.e. a tightly bound donor far outmatching the projectile charge)
    #     floors at the physical minimum n = 1.
    #
    println("\n  Empirical.chargeExchangeCaptureShell/StateSelectiveCrossSection for C6+ + H(1s):\n")
    setDefaults("nuclear: charge", 1.)
    IpH = Empirical.ionizationPotential(1, Configuration("1s^1"))
    capShell = Empirical.chargeExchangeCaptureShell(6.0, IpH)
    println("  q=6, Ip(H): nc = $(round(capShell.nc, digits=3)), dominant shell n = $(capShell.n)")
    stateSel = Empirical.chargeExchangeStateSelectiveCrossSection(6.0, IpH, Empirical.NiehausScaling1986(); printout=true)
    ##  Configuration-based wrapper, and a weakly-matched case
    stateSelWrap = Empirical.chargeExchangeStateSelectiveCrossSection(6.0, Configuration("1s^1"), 1.0,
                                                                       Empirical.NiehausScaling1986(); printout=false)
    println("  wrapper n = $(stateSelWrap.n), sigmaTotal matches direct: $(stateSelWrap.sigmaTotal == stateSel.sigmaTotal)")
    weakShell = Empirical.chargeExchangeCaptureShell(1.0, Defaults.convertUnits("energy: from eV to atomic", 500.0))
    println("  weakly-matched (q=1, Ip=500 eV): nc = $(round(weakShell.nc, digits=4)), floored n = $(weakShell.n)")
    #
    setDefaults("print summary: close", "")
    #
end
