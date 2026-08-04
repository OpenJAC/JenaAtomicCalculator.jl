#
println("Nc) Tests of empirical photoionization/photorecombination cross sections and plasma rate coefficients.")
setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "barn")
setDefaults("unit: rate",          "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  16-Jul-2026
    # Branch 1: Empirical.photoionizationCrossSection(omegas, iConf, fConf, ScaledHydrogenic()) --
    #   Kramers (1923) empirical formula for the PI cross section, scaled by the true binding energy and
    #   corrected by the bound-free Gaunt factor.
    # System: Ne (Z=10) 2p photoionization; several photon energies spanning the threshold.
    # Note: ScaledHydrogenic takes the threshold from the tabulated (Williams2000) binding energy, i.e. 21.6 eV
    #   for Ne 2p, and applies sigma = 32 pi alpha N_e /(3 sqrt(3) n) * bE^2/omega^3 * g^(bf), where N_e = 6 is
    #   the number of equivalent electrons in the 2p shell and g^(bf) the bound-free Gaunt factor of
    #   Empirical.boundFreeGauntFactor (exact for n = 1; Menzel-Pekeris n^(-2/3) scaling otherwise).
    # Checks:
    #   - Cross section = 0 below the 21.6 eV threshold.
    #   - Doubling omega (600 -> 1200 eV) reduces the cross section by 8 * g(27.8)/g(55.6) = 8.91; the Gaunt
    #     factor steepens Kramers' omega^-3 tail towards the exact omega^-(7/2).
    #   - At threshold 13.0 Mb against an experimental ~6.3 Mb (2.1x; l-averaged Kramers limit for p subshells).
    #     For He 1s^2, where the formula is clean, it gives 6.97 Mb against an experimental ~7.4 Mb (0.94x),
    #     and for H 1s the exact Stobbe threshold 6.30 Mb.
    #
    println("\n  Empirical.photoionizationCrossSection(omegas, iConf, fConf, ScaledHydrogenic()) -- Ne 2p:\n")
    Z     = 10.;  setDefaults("nuclear: charge", Z)
    iConf = Configuration("1s^2 2s^2 2p^6")
    fConf = Configuration("1s^2 2s^2 2p^5")
    ##  Photon energies spanning the true 2p threshold (21.6 eV); one point below, the others above.
    omegaEV  = [15.0, 25.0, 30.0, 600.0, 1200.0]
    omegas   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in omegaEV]
    css      = Empirical.photoionizationCrossSection(omegas, iConf, fConf, Empirical.ScaledHydrogenic(); printout=true)
    println("    omega [eV]    cross section [barn]")
    for  (omEV, cs)  in  zip(omegaEV, css)
        csx = round(Defaults.convertUnits("cross section: from atomic", cs), digits=3)
        println("    $omEV            $csx")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 2: Empirical.photorecombinationCrossSection(energies, iConf, fConf, ScaledHydrogenic()) --
    #   Einstein-Milne relation applied to the Kramers PI cross section to obtain PR cross sections.
    # System: Ne (Z=10) 2p recombination; several electron kinetic energies.
    # Note: the captured electron binds with the tabulated 2p binding energy of 21.6 eV, the same threshold
    #   as seen by the inverse (photoionization) process in Branch 1.
    # Note: the N_e = 6 of the PI cross section (Branch 1) cancels here against the statistical ratio
    #   g_f/g_i = g(2p^6)/g(2p^5) = 1/6 of the Einstein-Milne relation, since the 2p^5 ion offers exactly one
    #   vacancy to be filled. N_e counts the electrons that could be ejected, g_f/g_i the holes to be filled.
    # Checks:
    #   - PR cross section diverges as E_kin -> 0 (sigma_PR ~ 1/E_kin for radiative recombination).
    #   - sigma_PR * E_kin * (E_kin + 21.6 eV) / g^(bf) = const = 24556.3 barn eV^2 for all five energies;
    #     without the Gaunt factor g^(bf)(omega/21.6 eV) the product itself would be constant.
    #   - Cross section decreases with increasing electron energy.
    #
    println("\n  Empirical.photorecombinationCrossSection(energies, iConf, fConf, ScaledHydrogenic()) -- Ne 2p:\n")
    Z     = 10.;  setDefaults("nuclear: charge", Z)
    iConf = Configuration("1s^2 2s^2 2p^5")
    fConf = Configuration("1s^2 2s^2 2p^6")
    ekinEV   = [0.1, 1.0, 5.0, 20.0, 100.0]
    ekins    = [Defaults.convertUnits("energy: from eV to atomic", e) for e in ekinEV]
    prCss    = Empirical.photorecombinationCrossSection(ekins, iConf, fConf, Empirical.ScaledHydrogenic(); printout=false)
    println("    E_kin [eV]    PR cross section [barn]")
    for  (eEV, cs)  in  zip(ekinEV, prCss)
        csx = round(Defaults.convertUnits("cross section: from atomic", cs), digits=3)
        println("    $eEV            $csx")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 3: PI plasma rate per ion R^(PI) and PR plasma rate coefficients alpha^(PR: spontaneous/total).
    # System: Ne (Z=10) 2p; electron temperature T_e = 10 eV (low-temperature plasma).
    # Note: R^(PI) is a *rate* [1/s], not a rate coefficient. The convolution
    #   R^(PI: per ion) = int d(omega) n(omega;T) c sigma^(PI)(omega) already contains the photon number density
    #   of the radiation field, so that nothing has to be multiplied to it. The PR coefficient alpha^(PR) [cm^3/s],
    #   in contrast, still needs to be multiplied by the electron density n_e to obtain a rate.
    # Note: the total PR coefficient includes the stimulated capture,
    #   alpha^(PR: total) = int d(eps) f_e(eps;T_e) v(eps) sigma^(PR: spont)(eps) [1 + nbar(eps + bE; T)],
    #   where [1 + nbar] must stay *inside* the electron integral: PR is a bound-free process, whose emitted photon
    #   energy omega = eps + bE varies across the whole electron distribution, so that no single nbar can be
    #   factored out. The electron temperature T_e and the radiation temperature T are independent (non-LTE).
    # Checks:
    #   - R^(PI) is small at T = 10 eV, since only the tail of the Planck field passes the 21.6 eV threshold.
    #   - A vacuum photon field must reproduce alpha^(PR: spontaneous) exactly (nbar = 0).
    #   - The stimulated part grows with the radiation temperature: x 1.07 at T = 10 eV, x 36.2 at T = 1 keV.
    #
    println("\n  Empirical plasma rates and rate coefficients for Ne 2p at T_e = 10 eV:\n")
    Z     = 10.;  setDefaults("nuclear: charge", Z)
    iConf = Configuration("1s^2 2s^2 2p^6")
    fConf = Configuration("1s^2 2s^2 2p^5")
    T     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.16e5)    ##  T = T_e = 10 eV
    distP = Distribution.PhotonPlanck(T)
    distE = Distribution.ElectronMaxwell(T)
    ratPI = Empirical.photoionizationPlasmaRatePerIon(   distP, iConf, fConf, approx=Empirical.ScaledHydrogenic(), printout=true)
    alpPR = Empirical.photorecombinationPlasmaAlpha(distE, fConf, iConf, approx=Empirical.ScaledHydrogenic(), printout=true)
    ##  The PI rate is expressed in the user-defined rate unit, the PR coefficient in cm^3/s
    fac    = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    ratPIx = round(Defaults.convertUnits("rate: from atomic", ratPI), sigdigits=3)
    alpPRx = round(fac*alpPR, sigdigits=3)
    println("  R^(PI: per ion) [1/s]            = $ratPIx")
    println("  alpha^(PR: spontaneous) [cm^3/s] = $alpPRx")
    ##  Total PR coefficients for photon fields of different radiation temperature
    println("\n    photon field         alpha^(PR: total) [cm^3/s]    ratio to spontaneous")
    for  (label, dist)  in  [("vacuum        ", Distribution.PhotonVacuumField(0.)),
                             ("Planck,  10 eV", Distribution.PhotonPlanck(T)),
                             ("Planck,   1 keV", Distribution.PhotonPlanck(100*T))]
        alp  = Empirical.photorecombinationPlasmaAlpha(distE, dist, fConf, iConf, approx=Empirical.ScaledHydrogenic())
        alpx = round(fac*alp, sigdigits=4);    ratio = round(alp/alpPR, digits=4)
        println("    $label       $alpx                   $ratio")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 4: Empirical.boundFreeGauntFactor(n, x) -- the quantum correction to Kramers' formula.
    # Basis: for n = 1 the factor is exact, the ratio of Stobbe's (1930) closed-form 1s cross section to the
    #   Kramers form; for n >= 2 the deviation from unity is scaled as n^(-2/3) (Menzel & Pekeris 1935).
    # Checks:
    #   - g(1, x -> 1) = 0.7973: at threshold, Kramers' 7.907 Mb becomes the exact Stobbe 6.304 Mb.
    #   - g crosses unity near x = 4 and falls off as x^(-1/2) at large x (exact tail omega^-(7/2) vs
    #     Kramers omega^-3): g * sqrt(x) approaches 2 pi = 6.28 for x -> infinity.
    #   - g(n, threshold) rises towards 1 with n (semiclassical limit): 0.797, 0.872, 0.903, ...
    #
    println("\n  Empirical.boundFreeGauntFactor(n, x):\n")
    println("    x = omega/omega_thr   g(n=1)      g(n=2)      g(n=5)      g(n=1)*sqrt(x)")
    for  x  in  [1.000001, 1.5, 2.0, 4.0, 10.0, 40.0, 100.0, 1000.0]
        g1 = Empirical.boundFreeGauntFactor(1, x);   g2 = Empirical.boundFreeGauntFactor(2, x)
        g5 = Empirical.boundFreeGauntFactor(5, x)
        println("    $(rpad(x,12))          $(rpad(round(g1,digits=4),8))    $(rpad(round(g2,digits=4),8))    " *
                "$(rpad(round(g5,digits=4),8))    $(round(g1*sqrt(x),digits=3))")
    end
    println("\n    At threshold (x -> 1): g(n) = ", [round(Empirical.boundFreeGauntFactor(n, 1.000001), digits=4) for n = 1:6])
    #
    setDefaults("print summary: close", "")
    #
end
