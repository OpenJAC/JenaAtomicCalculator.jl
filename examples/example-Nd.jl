#
println("Nd) Tests of empirical electron-impact excitation (Van Regemorter) and ionization (Lotz) cross sections and rate coefficients.")
setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "barn")
setDefaults("unit: rate",          "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  16-Jul-2026
    # Branch 1: Empirical.impactIonizationCrossSection(energies, iConf, fConf, Lotz1967()) --
    #   simplified Lotz (1967) formula sigma = A xi ln(eps/P)/(eps P) with A = 4.5e-14 cm^2 eV^2.
    # Systems: He (Z=2) 1s^2 and Ne (Z=10) 2p^6; several electron energies spanning the thresholds.
    # Note: the simplified constant A = 4.5e-14 applies best to ions and overestimates light neutral
    #   atoms by ~40% (Lotz's own accuracy claim is +40/-30%); the full Lotz formula carries per-subshell
    #   constants a_i = 2.9 ... 4.5e-14 and two near-threshold parameters (b_i, c_i) that are not included.
    # Checks:
    #   - Cross section = 0 below the threshold P (He: 24.59 eV, Ne 2p: 21.6 eV).
    #   - The Lotz form peaks at eps = e * P (He: ~67 eV, Ne: ~59 eV) and falls off ~ln(eps)/eps.
    #   - He sigma(100 eV) = 5.13e-17 cm^2 = 5.13e7 barn, vs. the experimental ~3.6e-17 cm^2 (~1.4x high).
    #   - Neutral Ne 2p^6 (xi = 6) peaks near 2.1e-16 cm^2, ~2.5x above the experimental ~0.8e-16 cm^2:
    #     the simplified formula degrades for many-electron neutral shells, where Lotz's full per-subshell
    #     constants (a_i < 4.5e-14, b_i, c_i) would be needed.
    #
    println("\n  Empirical.impactIonizationCrossSection(energies, iConf, fConf, Lotz1967()):\n")
    epsEV = [15.0, 30.0, 67.0, 100.0, 500.0, 2000.0]
    eps   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in epsEV]
    ##  (a) He 1s^2 -> 1s^1
    setDefaults("nuclear: charge", 2.)
    cssHe = Empirical.impactIonizationCrossSection(eps, Configuration("1s^2"), Configuration("1s^1"),
                                                   Empirical.Lotz1967(); printout=true)
    ##  (b) Ne 2p^6 -> 2p^5
    setDefaults("nuclear: charge", 10.)
    cssNe = Empirical.impactIonizationCrossSection(eps, Configuration("1s^2 2s^2 2p^6"), Configuration("1s^2 2s^2 2p^5"),
                                                   Empirical.Lotz1967(); printout=false)
    println("    eps [eV]      He 1s^2 [barn]      Ne 2p^6 [barn]")
    for  (eEV, csa, csb)  in  zip(epsEV, cssHe, cssNe)
        csax = round(Defaults.convertUnits("cross section: from atomic", csa), sigdigits=4)
        csbx = round(Defaults.convertUnits("cross section: from atomic", csb), sigdigits=4)
        println("    $eEV          $csax            $csbx")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 2: Empirical.impactExcitationCrossSection(energies, iConf, fConf, VanRegemorter1962(); aSource=..) --
    #   Van Regemorter (1962) formula for optically allowed (E1) transitions; aSource selects how the underlying
    #   Einstein-A value is estimated (ScaledHydrogenic() or GivenEinsteinA(..)).
    # Systems:
    #   (a) H 1s -> 2p with the *exact* Einstein-A = 6.268e8 1/s via GivenEinsteinA; the implied
    #       oscillator strength f = (g_f/g_i) A/(2 alpha^3 omega^2) = 0.4162 (literature value).
    #   (b) Ne^+ 1s -> 2p (the inverse of K-alpha) with ScaledHydrogenic; deltaE = 848.6 eV.
    # Note: for the neutral H target the effective Gaunt factor gbar = 0.074 x (1 + x^2) vanishes at
    #   threshold, while for the Ne^+ ion gbar = 0.2 keeps the cross section finite there.
    # Checks:
    #   - Cross section = 0 below threshold (10.2 eV resp. 848.6 eV).
    #   - H: sigma(2 deltaE) = 0.794 a_o^2 = 2.22e7 barn; accuracy ~factor 2 (Van Regemorter's claim).
    #   - Ne^+: finite cross section just above threshold (ion), falling ~1/eps at high energies.
    #
    println("\n  Empirical.impactExcitationCrossSection(energies, iConf, fConf, VanRegemorter1962(); aSource=..):\n")
    ##  (a) H 1s -> 2p with the exact Einstein-A
    setDefaults("nuclear: charge", 1.)
    Aau   = 6.268e8 / (1.0 / Defaults.convertUnits("time: from atomic to sec", 1.0))
    given = Empirical.GivenEinsteinA(Basics.E1, 0.375, Aau)
    epsEV = [8.0, 15.0, 20.4, 40.8, 102.0, 408.0]
    eps   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in epsEV]
    cssH  = Empirical.impactExcitationCrossSection(eps, Configuration("1s^1"), Configuration("2p^1"), Empirical.VanRegemorter1962();
                                                   aSource=given, printout=true)
    ##  (b) Ne^+ 1s -> 2p with ScaledHydrogenic
    setDefaults("nuclear: charge", 10.)
    epsEV2 = [500.0, 900.0, 1700.0, 4243.0, 8486.0]
    eps2   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in epsEV2]
    cssNe  = Empirical.impactExcitationCrossSection(eps2, Configuration("1s^2 2s^2 2p^5"), Configuration("1s^1 2s^2 2p^6"),
                                                    Empirical.VanRegemorter1962(); aSource=Empirical.ScaledHydrogenic(), printout=true)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 3: Plasma rate coefficients alpha^(EIE) and alpha^(EII) for a Maxwellian electron gas --
    #   Gauss-Legendre integration of the cross sections over the electron distribution; the mesh starts
    #   at the threshold (an ion's EIE cross section is finite there) and follows the temperature.
    # Systems: H 1s -> 2p excitation (exact Einstein-A) and He 1s^2 ionization, for several T_e.
    # Checks:
    #   - Both coefficients grow steeply with T_e while exp(-threshold/T_e) is small (threshold-dominated).
    #   - alpha^(EII)(He; 5 eV) = 1.38e-10 cm^3/s;  alpha^(EIE)(H Ly-alpha; 2 eV) = 3.32e-11 cm^3/s.
    #   - Quantity: rate coefficients [cm^3/s]; multiply by n_e for the rate per atom/ion [1/s].
    #
    println("\n  Empirical plasma rate coefficients alpha^(EIE) and alpha^(EII) for Maxwellian electrons:\n")
    fac   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    ##  (a) H 1s -> 2p excitation
    setDefaults("nuclear: charge", 1.)
    Aau   = 6.268e8 / (1.0 / Defaults.convertUnits("time: from atomic to sec", 1.0))
    given = Empirical.GivenEinsteinA(Basics.E1, 0.375, Aau)
    println("    T_e [eV]     alpha^(EIE) H 1s->2p [cm^3/s]")
    for  TeV  in  [1.0, 2.0, 5.0, 10.0]
        eD  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", TeV))
        alp = Empirical.impactExcitationPlasmaAlpha(eD, Configuration("1s^1"), Configuration("2p^1"), aSource=given)
        println("    $TeV          $(round(fac*alp, sigdigits=4))")
    end
    ##  (b) He 1s^2 ionization
    setDefaults("nuclear: charge", 2.)
    println("\n    T_e [eV]     alpha^(EII) He 1s^2 [cm^3/s]")
    for  TeV  in  [2.0, 5.0, 10.0, 20.0]
        eD  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", TeV))
        alp = Empirical.impactIonizationPlasmaAlpha(eD, Configuration("1s^2"), Configuration("1s^1"))
        println("    $TeV          $(round(fac*alp, sigdigits=4))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 4: Empirical.forbiddenExcitationCrossSection/PlasmaAlpha(..., ConstantCollisionStrength()) --
    #   electron-impact excitation of an optically forbidden (M1 or E2) transition via a constant, order-unity
    #   collision strength Omega [M. J. Seaton, Proc. Phys. Soc. 79, 1105 (1962)].
    # Systems: H 1s -> 2s (M1, Delta-l = 0, the classic "forbidden" 2s-1s transition) and H 1s -> 3d (E2,
    #   Delta-l = 2); both with the default Omega = 1.
    # Checks:
    #   - Cross section = 0 below threshold (10.2 eV resp. 12.1 eV).
    #   - 1s -> 2p (E1, Delta-l = 1) is rejected, pointing to VanRegemorter1962(); 1s -> 4f (Delta-l = 3) is
    #     rejected outright (no useful collision-strength estimate attempted beyond M1/E2).
    #   - The Maxwellian-folded plasma rate coefficient reproduces the textbook formula
    #     C_ij = 8.629e-6/(g_i sqrt(T_e[K])) exp(-deltaE/T_e) [cm^3/s] (Osterbrock & Ferland) to <0.01%,
    #     an independent, literature-based cross check of the Gauss-Legendre folding used here.
    #
    println("\n  Empirical.forbiddenExcitationCrossSection(energies, iConf, fConf, ConstantCollisionStrength()):\n")
    setDefaults("nuclear: charge", 1.)
    h1s = Configuration("1s^1");   h2s = Configuration("2s^1");   h3d = Configuration("3d^1")
    epsEV3 = [5.0, 10.2, 12.1, 15.0, 40.8]
    eps3   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in epsEV3]
    ##  (a) H 1s -> 2s (M1)
    cssM1  = Empirical.forbiddenExcitationCrossSection(eps3, h1s, h2s, Empirical.ConstantCollisionStrength(); printout=true)
    ##  (b) H 1s -> 3d (E2)
    cssE2  = Empirical.forbiddenExcitationCrossSection(eps3, h1s, h3d, Empirical.ConstantCollisionStrength(); printout=true)
    ##  (c) error paths: E1 (rejected, use VanRegemorter1962) and Delta-l = 3 (unsupported)
    try     Empirical.forbiddenExcitationCrossSection(eps3, h1s, Configuration("2p^1"), Empirical.ConstantCollisionStrength())
    catch e
        println("  1s -> 2p (E1) correctly rejected: ", sprint(showerror, e))
    end
    try     Empirical.forbiddenExcitationCrossSection(eps3, h1s, Configuration("4f^1"), Empirical.ConstantCollisionStrength())
    catch e
        println("  1s -> 4f (Delta-l = 3) correctly rejected: ", sprint(showerror, e))
    end
    ##  (d) Maxwellian plasma rate coefficient vs. the textbook formula
    fac3   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    deltaE = Empirical.scaledBindingEnergy(1.0, Shell("1s"), h1s, PeriodicTable.Williams2000()) -
             Empirical.scaledBindingEnergy(1.0, Shell("2s"), h2s, PeriodicTable.Williams2000())
    deltaEeV = Defaults.convertUnits("energy: from atomic to eV", deltaE)
    gi     = Basics.extractFromConfiguration(Basics.Multiplicity(), h1s)
    println("\n    T_e [eV]     alpha^(EIE) H 1s->2s [cm^3/s]   textbook [cm^3/s]   ratio")
    for  TeV  in  [1.0, 5.0, 20.0, 50.0]
        eD   = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", TeV))
        alp  = fac3 * Empirical.forbiddenExcitationPlasmaAlpha(eD, h1s, h2s, printout=false)
        TeK  = TeV * 11604.518
        text = 8.629e-6 / (gi * sqrt(TeK)) * exp(-deltaEeV*11604.518/TeK)
        println("    $TeV          $(round(alp, sigdigits=4))              $(round(text, sigdigits=4))         $(round(alp/text, digits=4))")
    end
    #
    setDefaults("print summary: close", "")
    #
end
