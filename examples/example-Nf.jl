#
println("Nf) Tests of empirical stopping powers: projectile x material x approximation.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  17-Jul-2026
    # Branch 1: Empirical.stoppingPower(energies, projectile, material, approx) -- the four approximations
    #   compared for an electron in a free-electron gas. All share
    #   - dE/dx = n_e (4 pi z^2 e^4/(m_e v^2)) ln Lambda and differ only in the Coulomb logarithm:
    #   Bohr1913 (classical, heavy projectile), Bethe1931 (quantal, with the -ln 2 electron-projectile
    #   correction), KozmaFranson1992 (piecewise: quantal above / classical below 14 eV) and Axelrod1980
    #   (relativistic, hbar omega_p as effective ionization potential).
    # System: ElectronProjectile() in FreeElectronGas(n_e = 1e20 cm^-3), hbar omega_p = 0.37 eV.
    # Checks:
    #   - KF92 equals Bethe1931 above 14 eV; below, it follows the classical form.
    #   - The Bohr and KF92 logarithms differ by ln 2 (electron-projectile correction: 1.123 x 1.7811 = 2).
    #   - Axelrod1980/Bethe1931 -> 1 in the nonrelativistic limit and rises to ~2.5 at 500 keV,
    #     where the nonrelativistic 1/E fall-off is no longer valid.
    #   - All stopping powers decrease with energy in the nonrelativistic regime.
    #
    println("\n  Electron in a free-electron gas with n_e = 1e20 cm^-3:\n")
    neAu  = 1.0e20 * Defaults.convertUnits("length: from atomic to cm", 1.0)^3
    ele   = Empirical.ElectronProjectile();    feg = Empirical.FreeElectronGas(neAu)
    epsEV = [5.0, 20.0, 100.0, 1000.0, 1.0e4, 1.0e5, 5.0e5]
    eps   = [Defaults.convertUnits("energy: from eV to atomic", e) for e in epsEV]
    spBoh = Empirical.stoppingPower(eps, ele, feg, Empirical.Bohr1913())
    spBet = Empirical.stoppingPower(eps, ele, feg, Empirical.Bethe1931())
    spKF  = Empirical.stoppingPower(eps, ele, feg, Empirical.KozmaFranson1992(), printout=true)
    spAx  = Empirical.stoppingPower(eps, ele, feg, Empirical.Axelrod1980())
    aoCm  = Defaults.convertUnits("length: from atomic to cm", 1.0)
    println("    E [eV]      Bohr1913 [eV/cm]   Bethe1931 [eV/cm]  KozmaFranson1992   Axelrod1980")
    for  (i, eEV)  in  enumerate(epsEV)
        f(x) = round(Defaults.convertUnits("energy: from atomic to eV", x)/aoCm, sigdigits=4)
        println("    $(rpad(eEV,10))  $(rpad(f(spBoh[i]),18)) $(rpad(f(spBet[i]),18)) $(rpad(f(spKF[i]),18)) $(f(spAx[i]))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 2: the Quantity relation -- stopping power vs. stopping cross section.
    #   The stopping power - dE/dx [eV/cm] contains the electron density n_e; the stopping cross
    #   section S = (-dE/dx)/n_e [eV cm^2] is (nearly) a property of the projectile alone and depends
    #   on the density only logarithmically, through hbar omega_p ~ sqrt(n_e) inside ln Lambda.
    # Checks:
    #   - Over 12 orders of magnitude in n_e (1e14 ... 1e26 cm^-3), - dE/dx changes by ~11 orders,
    #     while S changes by less than one order (a factor ~9, purely the Coulomb logarithm, which
    #     shrinks towards solid density as hbar omega_p approaches the projectile energy).
    #
    println("\n  Stopping power vs. stopping cross section for a 1 keV electron (KozmaFranson1992):\n")
    aoCm  = Defaults.convertUnits("length: from atomic to cm", 1.0)
    ele   = Empirical.ElectronProjectile()
    eps   = [Defaults.convertUnits("energy: from eV to atomic", 1000.0)]
    println("    n_e [1/cm^3]    - dE/dx [eV/cm]     S = (-dE/dx)/n_e [eV cm^2]")
    for  neCm  in  [1.0e14, 1.0e17, 1.0e20, 1.0e23, 1.0e26]
        neAu2 = neCm * aoCm^3
        sp   = Empirical.stoppingPower(eps, ele, Empirical.FreeElectronGas(neAu2), Empirical.KozmaFranson1992())[1]
        spx  = Defaults.convertUnits("energy: from atomic to eV", sp) / aoCm
        println("    $(rpad(neCm,15)) $(rpad(round(spx, sigdigits=4),19)) $(round(spx/neCm, sigdigits=4))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 3: CSDA range of an electron in a fully ionized plasma -- an application.
    #   In the continuous slowing-down approximation, the path length until thermalization is
    #       R = int_Emin^E0  dE / (-dE/dx),
    #   here integrated with the KozmaFranson1992 loss function down to the density-dependent floor
    #   E_min at which ln Lambda = 1; below this floor the loss function vanishes (the clamped Coulomb
    #   logarithm) and the electron is quasi-thermal, so that the CSDA picture no longer applies.
    # System: electrons of 0.1, 1 and 10 keV in plasmas of n_e = 1e18 and 1e20 cm^-3.
    # Checks:
    #   - The range grows roughly as E0^2/ln(E0) (from the 1/E fall-off of the loss function).
    #   - The range scales as ~1/n_e (up to the Coulomb logarithm): 100x lower density, ~100x longer range.
    #
    println("\n  CSDA range of an electron (KozmaFranson1992), R = int dE/(-dE/dx) down to ln Lambda = 1:\n")
    aoCm  = Defaults.convertUnits("length: from atomic to cm", 1.0)
    ele   = Empirical.ElectronProjectile()
    println("    E0 [keV]    n_e [1/cm^3]    range R [cm]")
    for  neCm  in  [1.0e18, 1.0e20]
        neAu3  = neCm * aoCm^3
        feg3   = Empirical.FreeElectronGas(neAu3)
        omegaP = sqrt(4pi * neAu3)
        for  E0keV  in  [0.1, 1.0, 10.0]
            E0   = Defaults.convertUnits("energy: from eV to atomic", 1000.0*E0keV)
            ##  Density-dependent floor: ln Lambda = 1 at v^3 = e * 1.7811 * omega_p
            Emin = (exp(1.0) * 1.7811 * omegaP)^(2/3) / 2
            ##  Simple trapezoidal quadrature on a logarithmic energy grid
            nPts = 400;    R = 0.
            es   = [Emin * (E0/Emin)^(k/nPts)  for k = 0:nPts]
            sps  = Empirical.stoppingPower(es, ele, feg3, Empirical.KozmaFranson1992())
            for  k = 1:nPts
                R = R + (es[k+1] - es[k]) * 0.5 * (1.0/sps[k+1] + 1.0/sps[k])
            end
            println("    $(rpad(E0keV,11)) $(rpad(neCm,15)) $(round(R*aoCm, sigdigits=4))")
        end
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 4: protons in solid aluminum and in a partially ionized gas -- the projectile and material axes.
    #   (a) IonProjectile(1, M_p) in NeutralAtomGas(13, 26.98, n_atom of solid Al): the classic Bethe
    #       stopping with Ibar = 159 eV from the Segre/Roy-Reed estimate. The bare formula omits the shell,
    #       Barkas and relativistic corrections; at 10 MeV it yields 34.2 MeV cm^2/g against ~40 of
    #       PSTAR-type tabulations (~15% low), and it degrades quickly below ~1 MeV, where such
    #       corrections dominate (cf. ICRU-49).
    #   (b) The equal-velocity identity: a proton with E_p = (M_p/m_e) x 100 eV = 184 keV has the same
    #       velocity as a 100 eV electron; the Bethe logarithms then differ by exactly ln 2 (kappa = 2
    #       for the distinguishable proton vs. kappa = 1 for the electron).
    #   (c) An electron in partially ionized Al vapor (chi_e = 3): the loss decomposes exactly into the
    #       bound (Ibar) and free (hbar omega_p) populations, cf. Milne et al. (1999, Eq. 3).
    #
    println("\n  Protons in solid Al; the equal-velocity identity; partially ionized material:\n")
    aoCm  = Defaults.convertUnits("length: from atomic to cm", 1.0)
    Mp    = Defaults.PROTON_MASS_U / Defaults.ELECTRON_MASS_U
    prot  = Empirical.IonProjectile(1.0, Mp)
    ##  (a) solid Al: rho = 2.699 g/cm^3
    gPerU = Defaults.ELECTRON_MASS_IN_G / Defaults.ELECTRON_MASS_U
    natAl = 2.699 / (26.98 * gPerU) * aoCm^3
    alSol = Empirical.NeutralAtomGas(13, 26.98, natAl)
    epsMeV = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    eps    = [Defaults.convertUnits("energy: from eV to atomic", 1.0e6*e) for e in epsMeV]
    sps    = Empirical.stoppingPower(eps, prot, alSol, Empirical.Bethe1931(), printout=true)
    rho    = natAl / aoCm^3 * 26.98 * gPerU
    println("    E_p [MeV]    - dE/dx [MeV/cm]    - (dE/dx)/rho [MeV cm^2/g]")
    for  (i, eMeV)  in  enumerate(epsMeV)
        spx = Defaults.convertUnits("energy: from atomic to eV", sps[i]) / aoCm * 1.0e-6
        println("    $(rpad(eMeV,12)) $(rpad(round(spx, sigdigits=4),19)) $(round(spx/rho, sigdigits=4))")
    end
    ##  (b) equal-velocity identity
    neAu4 = 1.0e20 * aoCm^3;    feg4 = Empirical.FreeElectronGas(neAu4)
    ele4  = Empirical.ElectronProjectile()
    e100  = Defaults.convertUnits("energy: from eV to atomic", 100.0)
    lE    = Empirical.stoppingPower([e100],      ele4, feg4, Empirical.Bethe1931())[1]
    lP    = Empirical.stoppingPower([e100 * Mp], prot, feg4, Empirical.Bethe1931())[1]
    println("\n    Equal velocity: lnLambda(p) - lnLambda(e) = $(round((lP - lE) * 2*e100/(4pi*neAu4), digits=6))   [ln 2 = 0.693147]")
    ##  (c) electron in partially ionized Al vapor
    natV  = 1.0e19 * aoCm^3
    lPig  = Empirical.stoppingPower([Defaults.convertUnits("energy: from eV to atomic", 1.0e4)], ele4,
                                    Empirical.PartiallyIonizedGas(13, 26.98, 3.0, natV), Empirical.Axelrod1980(), printout=true)
    #
    setDefaults("print summary: close", "")
    #
end
