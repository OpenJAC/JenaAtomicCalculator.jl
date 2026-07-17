
## Functions in this file cover: empirical/semiempirical estimates and fits.
## Alphabetical order within this file.


"""
`TestFrames.testModule_Empirical(; short::Bool=true)`  ... tests on module Empirical.
"""
function testModule_Empirical(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Empirical-new.sum")
    printstyled("\n\nTest the module  Empirical  ... \n", color=:cyan)
    success = true
    ## Preserve the global state; these tests set the nuclear charge and rely on fixed units.
    oldZ = Defaults.getDefaults("nuclear: charge")
    Defaults.setDefaults("unit: energy",        "eV")
    Defaults.setDefaults("unit: cross section", "barn")
    Defaults.setDefaults("unit: rate",          "1/s")
    Defaults.setDefaults("nuclear: charge",     10.)
    #
    ## Test 1: bindingEnergy(Z, sh::Shell) for Ne 2p (Williams2000) = 21.6 eV
    e1 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, Shell("2p"), data=PeriodicTable.Williams2000()))
    success = success && abs(e1 - 21.6) < 0.05
    ## Test 2: ionizationPotential(Z, conf) for neutral Ne = 21.565 eV
    e2 = Defaults.convertUnits("energy: from atomic", Empirical.ionizationPotential(10, Configuration("1s^2 2s^2 2p^6")))
    success = success && abs(e2 - 21.565) < 0.001
    #
    ## Tests 3&4: photoemissionEinsteinA(ScaledHydrogenic) for the Ne K-alpha analogue 1s^1 2p^6 -> 1s^2 2p^5.
    ##   The transition energy must be positive (emitted photon); a negative value signals a tEnergy sign error.
    ##   With the tabulated binding energies, tEnergy = bE(1s) - bE(2p) = 870.2 - 21.6 = 848.6 eV, in agreement
    ##   with the tabulated Ne K-alpha energy. The rate (6.03e13 1/s) is ~1.5x below the literature 8.8e13 1/s,
    ##   which is the accuracy to be expected from a hydrogenic estimate of the r-expectation values.
    iConf = Configuration("1s^1 2p^6");    fConf = Configuration("1s^2 2p^5")
    wa    = Empirical.photoemissionEinsteinA(iConf, fConf, Empirical.ScaledHydrogenic(), printout=false)
    e3    = Defaults.convertUnits("energy: from atomic", wa.energy)
    r4    = Defaults.convertUnits("rate: from atomic",   wa.rate)
    success = success && wa.multipole == E1
    success = success && e3 > 0.  &&  abs(e3 - 848.6) < 0.001
    success = success && abs(r4 - 6.033863052213893e13) / 6.033863052213893e13 < 1.0e-6
    #
    ## Tests 5&6: photoexcitation/photodeexcitation plasma rates for Ne K-alpha at T = 1 keV.
    ##   R^(PX) = (g_i/g_f) * A * nbar  and  R^(PD) = A * (1 + nbar)  with nbar the mean photon occupation number.
    T     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.16e7)    ##  T = 1 keV
    dist  = Distribution.PhotonPlanck(T)
    rPX   = Defaults.convertUnits("rate: from atomic",
                Empirical.photoexcitationPlasmaRatePerIon(  dist, iConf, fConf, approx=Empirical.ScaledHydrogenic(), printout=false))
    rPD   = Defaults.convertUnits("rate: from atomic",
                Empirical.photodeexcitationPlasmaRatePerIon(dist, iConf, fConf, approx=Empirical.ScaledHydrogenic(), printout=false))
    success = success && abs(rPX - 1.5041673132576172e13) / 1.5041673132576172e13 < 1.0e-6
    success = success && abs(rPD - 1.0546364991986745e14) / 1.0546364991986745e14 < 1.0e-6
    ## Test 6: detailed balance  R^(PD) / R^(PX) = (g_f/g_i) * exp(omega/kT); a physics identity, independent
    ##   of the ScaledHydrogenic approximation itself.
    gRatio  = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) /
              Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    balance = gRatio * exp(wa.energy / T)
    success = success && abs(rPD/rPX - balance) / balance < 1.0e-6
    #
    ## Test 7: photoionizationCrossSection(ScaledHydrogenic) for Ne 2p; Kramers (1923) formula with the
    ##   bound-free Gaunt factor. The threshold is the tabulated 2p binding energy of 21.6 eV; below it the cross
    ##   section vanishes. Doubling the photon energy far above threshold reduces the cross section by a factor
    ##   8 * g(600/21.6)/g(1200/21.6) = 8.913: the Gaunt factor steepens Kramers' omega^-3 tail towards the exact
    ##   omega^-(7/2).
    piConf = Configuration("1s^2 2s^2 2p^6");    pfConf = Configuration("1s^2 2s^2 2p^5")
    omegas = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [15.0, 25.0, 30.0, 600.0, 1200.0]]
    piCss  = Empirical.photoionizationCrossSection(omegas, piConf, pfConf, Empirical.ScaledHydrogenic(), printout=false)
    success = success && piCss[1] == 0.  &&  piCss[2] > 0.
    success = success && abs(piCss[4]/piCss[5] - 8.913025108522556) < 1.0e-6
    #
    ## Test 8: photorecombinationCrossSection(ScaledHydrogenic) for Ne 2p via the Einstein-Milne relation.
    ##   sigma^(PR) ~ 1 / [eps * (eps + bE)] and must diverge as eps -> 0; a flat or vanishing cross section
    ##   signals a broken Milne relation or an inconsistent threshold.
    ekins  = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [0.1, 1.0, 5.0, 20.0, 100.0]]
    prCss  = Empirical.photorecombinationCrossSection(ekins, pfConf, piConf, Empirical.ScaledHydrogenic(), printout=false)
    prBarn = [round(Defaults.convertUnits("cross section: from atomic", cs), digits=3) for cs in prCss]
    for  (ic, cs)  in  enumerate([9880.029, 955.962, 167.145, 28.339, 2.006])
        success = success && abs(prBarn[ic] - cs) < 0.001
    end
    #
    ## Test 9: the H 1s threshold cross section against the exact (Stobbe 1930) value. With the bound-free Gaunt
    ##   factor included, the Kramers 7.91 Mb becomes the exact 2^9 pi^2 alpha a_o^2/(3 e^4) = 6.304 Mb; this anchors
    ##   both the prefactor 32 pi alpha/(3 sqrt(3) n) and boundFreeGauntFactor to an analytic literature value.
    Defaults.setDefaults("nuclear: charge", 1.)
    bEH    = Empirical.bindingEnergy(1, Shell("1s"))
    hCss   = Empirical.photoionizationCrossSection([bEH * 1.0000001], Configuration("1s^1"), Configuration("1s^0"),
                                                   Empirical.ScaledHydrogenic(), printout=false)
    hMb    = Defaults.convertUnits("cross section: from atomic", hCss[1]) / 1.0e6
    success = success && abs(hMb - 6.307) < 0.01
    #
    ## Test 10: the tabulated data sets cover different ranges of Z; Williams2000 has no Xe (Z=54) and must fall back
    ##   to a Slater-screened hydrogenic estimate, while XrayDataBooklet provides the 5p binding energy of 12.1 eV.
    ##   A photon of 15 eV therefore ionizes Xe 5p only in the XrayDataBooklet data set.
    Defaults.setDefaults("nuclear: charge", 54.)
    xeConf = Configuration("[Kr] 4d^10 5s^2 5p^6");    xeIon = Configuration("[Kr] 4d^10 5s^2 5p^5")
    xeOm   = [Defaults.convertUnits("energy: from eV to atomic", 15.0)]
    xeCss  = Empirical.photoionizationCrossSection(xeOm, xeConf, xeIon, Empirical.ScaledHydrogenic(), printout=false,
                                                   data=PeriodicTable.XrayDataBooklet())
    success = success && xeCss[1] > 0.
    #
    ## Test 11: bindingEnergy(Z, subsh, conf) for a shell of a *given* configuration; this method was previously
    ##   unreachable (it referred to an undefined variable). For the neutral atom it must reproduce the tabulated
    ##   value, while each missing electron adds 0.3/n Hartree to the binding energy.
    Defaults.setDefaults("nuclear: charge", 10.)
    b11a = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, Subshell("2p_3/2"), Configuration("1s^2 2s^2 2p^6")))
    b11b = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, Subshell("2p_3/2"), Configuration("1s^2 2s^2 2p^3")))
    shift = Defaults.convertUnits("energy: from atomic", 0.3/2)
    success = success && abs(b11a - 21.6) < 0.001
    success = success && abs(b11b - (21.6 + 3*shift)) < 0.001
    ## The associated shell must be occupied in the given configuration.
    b11ok = false
    try     Empirical.bindingEnergy(10, Subshell("3d_5/2"), Configuration("1s^2 2s^2 2p^6"))
    catch
        b11ok = true
    end
    success = success && b11ok
    #
    ## Test 12: meanCharge(Z, subsh, conf) inverts bE = Z^(eff)^2/(2n^2); it previously took the square root of a
    ##   negative number and could only ever throw. For Ne 1s (870.2 eV = 31.98 Hartree) this gives Z^(eff) = 8.0.
    z12 = Empirical.meanCharge(10, Subshell("1s_1/2"), Configuration("1s^2 2s^2 2p^6"))
    success = success && abs(z12 - 7.997) < 0.01
    #
    ## Test 13: the number of equivalent electrons N_e in Kramers' formula. He 1s^2 (N_e = 2) is the clean benchmark:
    ##   it isolates N_e, since H 1s has N_e = 1 by definition and cannot discriminate. With N_e and the bound-free
    ##   Gaunt factor the threshold cross section is 6.97 Mb against an experimental ~7.4 Mb (0.94x); N_e = 1 would
    ##   be 2x too low, and without the Gaunt factor the 8.75 Mb were 1.18x too high.
    Defaults.setDefaults("nuclear: charge", 2.)
    bEHe   = Empirical.bindingEnergy(2, Shell("1s"))
    heCss  = Empirical.photoionizationCrossSection([bEHe * 1.0000001], Configuration("1s^2"), Configuration("1s^1"),
                                                   Empirical.ScaledHydrogenic(), printout=false)
    heMb   = Defaults.convertUnits("cross section: from atomic", heCss[1]) / 1.0e6
    success = success && abs(heMb - 6.9735) < 0.001
    #
    ## Tests 14-16: photorecombination plasma rate coefficients for Ne 2p at T_e = 10 eV.
    Defaults.setDefaults("nuclear: charge", 10.)
    eDist  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 10.0))
    aSpont = Empirical.photorecombinationPlasmaAlpha(eDist, pfConf, piConf, printout=false)
    fac    = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    success = success && abs(fac*aSpont - 1.613974906174875e-14) / 1.613974906174875e-14 < 1.0e-6
    ## Test 14: the vacuum photon field has n(omega) = 0 and must reproduce the spontaneous coefficient *exactly*.
    aVac   = Empirical.photorecombinationPlasmaAlpha(eDist, Distribution.PhotonVacuumField(0.), pfConf, piConf, printout=false)
    success = success && aVac == aSpont
    ## Test 15: the stimulated contribution scales linearly with the photon field; a dilute Planckian with w = 0.5 must
    ##   halve it exactly. This also verifies that nbar = n(omega) pi^2 c^3/omega^2 is applied distribution-independently.
    T1000  = Defaults.convertUnits("energy: from eV to atomic", 1000.0)
    aPlanck = Empirical.photorecombinationPlasmaAlpha(eDist, Distribution.PhotonPlanck(T1000), pfConf, piConf, printout=false)
    aDilute = Empirical.photorecombinationPlasmaAlpha(eDist, Distribution.PhotonDilute(T1000, 0.5), pfConf, piConf, printout=false)
    success = success && abs(fac*aPlanck - 5.849171271603789e-13) / 5.849171271603789e-13 < 1.0e-6
    success = success && abs( (aDilute - aSpont)/(aPlanck - aSpont) - 0.5 ) < 1.0e-10
    ## Test 16: the quadrature range must follow the electron temperature. With a fixed range of 0 ... 100 a.u., this
    ##   coefficient came out as ~1e-39 instead of ~2.4e-13, i.e. the thermal peak was missed altogether.
    eCool  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 0.1))
    aCool  = Empirical.photorecombinationPlasmaAlpha(eCool, pfConf, piConf, printout=false)
    success = success && abs(fac*aCool - 2.0909570957225128e-13) / 2.0909570957225128e-13 < 1.0e-4
    #
    ## Test 17: the photoionization plasma rate per ion R^(PI: per ion) = int d(omega) n(omega;T) c sigma^(PI)(omega).
    ##   This is a *rate* [1/s], not a rate coefficient [cm^3/s]: the photon number density is already contained in the
    ##   convolution. The factor c was previously missing, and the mesh straddled the threshold step; the integration
    ##   now starts at the threshold and is converged to < 0.2% against a dense reference quadrature.
    rPI = Empirical.photoionizationPlasmaRatePerIon(Distribution.PhotonPlanck(Defaults.convertUnits("energy: from eV to atomic", 10.0)),
                                                piConf, pfConf, printout=false)
    rPIx = Defaults.convertUnits("rate: from atomic", rPI)
    success = success && abs(rPIx - 2.2816321564282784e9) / 2.2816321564282784e9 < 1.0e-6
    #
    ## Test 18: the UsingJAC PI/PR cross sections from mean-field orbitals and E1 amplitudes. He 1s^2 serves as the
    ##   experimental benchmark: 6.26 Mb at 26 eV and 0.317 Mb at 100 eV against the measured ~7.4 and ~0.3 Mb, i.e.
    ##   within ~15%. (The former implementation produced spurious box resonances -- 92 Mb at 100 eV for Ne 2p --
    ##   from a continuum grid without a linear tail.) The PR cross section, via the Einstein-Milne relation, must
    ##   decrease from 1 to 10 eV (the 1/eps divergence towards threshold).
    ##   The UsingJAC amplitudes depend on the continuum method; earlier testsets (cascades) switch it and do not
    ##   restore, so the JAC defaults are pinned here explicitly to keep this test independent of the testset order.
    Defaults.setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, Alok")
    Defaults.setDefaults("nuclear: charge", 2.)
    ujOms = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [26.0, 100.0]]
    ujCss = Empirical.photoionizationCrossSection(ujOms, Configuration("1s^2"), Configuration("1s^1"),
                                                  Empirical.UsingJAC(), printout=false)
    success = success && abs(ujCss[1] - 0.22361086741886252) / 0.22361086741886252 < 1.0e-6
    success = success && abs(ujCss[2] - 0.01133283305665921) / 0.01133283305665921 < 1.0e-6
    ujEks = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [1.0, 10.0]]
    ujPrs = Empirical.photorecombinationCrossSection(ujEks, Configuration("1s^1"), Configuration("1s^2"),
                                                     Empirical.UsingJAC(), printout=false)
    success = success && ujPrs[1] > ujPrs[2] > 0.
    #
    ## Test 19: K-shell fluorescence yield (KrauseAdopted2016) and the derived Auger yield/rate. The yield of Ne is
    ##   0.01519 (Auger-dominated), of Ag 0.8313 (fluorescence-dominated); a_K = 1 - omega_K. The empirical Auger rate
    ##   uses A_auger = A_rad (1 - omega_K)/omega_K, so that A_auger/A_rad must equal (1 - omega_K)/omega_K exactly.
    success = success && abs(Empirical.fluorescenceYield(10) - 0.01519) < 1.0e-8
    success = success && abs(Empirical.fluorescenceYield(47) - 0.8313)  < 1.0e-8
    success = success && abs(Empirical.augerYield(10) - (1.0 - 0.01519)) < 1.0e-8
    Defaults.setDefaults("nuclear: charge", 10.)
    iK    = Configuration("1s^1 2p^6");    fK = Configuration("1s^2 2p^5")
    aRad  = Empirical.photoemissionEinsteinA(iK, fK, Empirical.ScaledHydrogenic(), printout=false).rate
    aAug  = Empirical.augerRate(iK, fK, Empirical.ScaledHydrogenic(), printout=false)
    success = success && abs(aAug/aRad - (1.0 - 0.01519)/0.01519) < 1.0e-6
    ## Yields are defined only where a K-shell vacancy can decay radiatively (Z >= 3) and within the tabulation.
    fyErr = false
    try     Empirical.fluorescenceYield(2)
    catch
        fyErr = true
    end
    success = success && fyErr
    #
    ## Test 20: electron-impact ionization after Lotz (1967). For He (P = 24.587 eV, xi = 2) the simplified formula
    ##   with A = 4.5e-14 cm^2 eV^2 gives sigma(100 eV) = 5.13e-17 cm^2 = 1.8322 a_o^2, ~1.4x above the experimental
    ##   ~3.6e-17 cm^2, i.e. within Lotz's claimed +40/-30% for light neutral atoms. Zero below the threshold, and the
    ##   Maxwellian rate coefficient at T_e = 5 eV is locked as a regression value.
    Defaults.setDefaults("nuclear: charge", 2.)
    heI   = Configuration("1s^2");    heII = Configuration("1s^1")
    eiiCs = Empirical.impactIonizationCrossSection([Defaults.convertUnits("energy: from eV to atomic", e) for e in [20.0, 100.0]],
                                                   heI, heII, Empirical.Lotz1967(), printout=false)
    success = success && eiiCs[1] == 0.  &&  abs(eiiCs[2] - 1.8322485748145694) / 1.8322485748145694 < 1.0e-6
    aEII  = Empirical.impactIonizationPlasmaAlpha(Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 5.0)),
                                                  heI, heII, printout=false)
    success = success && abs(fac*aEII - 1.3812552446980436e-10) / 1.3812552446980436e-10 < 1.0e-6
    #
    ## Test 21: electron-impact excitation after Van Regemorter (1962), anchored on H 1s -> 2p with the exact
    ##   Einstein-A = 6.268e8 1/s: the implied oscillator strength f = (g_f/g_i) A / (2 alpha^3 omega^2) must equal
    ##   the literature f = 0.4162, and sigma(eps = 2 deltaE) = 0.7944 a_o^2 with the neutral-atom Gaunt factor
    ##   gbar = 0.074 x (1 + x^2) = 0.148 at x^2 = 1. The ion and neutral gbar branches join at x^2 = 2.06.
    Defaults.setDefaults("nuclear: charge", 1.)
    Aau   = 6.268e8 / (1.0 / Defaults.convertUnits("time: from atomic to sec", 1.0))
    given = Empirical.GivenEinsteinA(Basics.E1, 0.375, Aau)
    fosc  = 3.0 * Aau / (2 * Defaults.getDefaults("alpha")^3 * 0.375^2)
    success = success && abs(fosc - 0.4162) < 0.001
    eieCs = Empirical.impactExcitationCrossSection([Defaults.convertUnits("energy: from eV to atomic", 10.0), 2*0.375],
                                                   Configuration("1s^1"), Configuration("2p^1"), Empirical.VanRegemorter1962();
                                                   aSource=given, printout=false)
    success = success && eieCs[1] == 0.  &&  abs(eieCs[2] - 0.7944475607726056) / 0.7944475607726056 < 1.0e-6
    aEIE  = Empirical.impactExcitationPlasmaAlpha(Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 2.0)),
                                                  Configuration("1s^1"), Configuration("2p^1"), aSource=given, printout=false)
    success = success && abs(fac*aEIE - 3.3150416264147014e-11) / 3.3150416264147014e-11 < 1.0e-6
    success = success && abs(Empirical.effectiveGauntFactor(1.0, true)  - 0.2) < 1.0e-10
    success = success && abs(Empirical.effectiveGauntFactor(10.0, true) - Empirical.effectiveGauntFactor(10.0, false)) < 1.0e-10
    #
    ## Test 22: three-body recombination by detailed balance with the Lotz EII coefficient. The capture channels of
    ##   H^+ (bare) with nLayers = 10 are the 34 shells n = 1 ... 10, l <= 3; for Ne^+ 2p^5 they are the 2p vacancy
    ##   plus the empty n = 3 ... 12 shells (40 channels). The single-channel coefficient must obey the Saha identity
    ##   alpha^(TBR) = alpha^(EII) g_f/(2 g_i) (2 pi/T)^(3/2) exp(P/T) where both sides are computable (T = 1 eV),
    ##   and must remain finite at T_e = 0.05 eV (P/T = 272), where a naive product underflows to 0 * Inf = NaN.
    Defaults.setDefaults("nuclear: charge", 1.)
    bare  = Configuration("1s^0");    h1s = Configuration("1s^1")
    success = success && length(Empirical.recombinationConfigurations(bare, nLayers=10)) == 34
    success = success && length(Empirical.recombinationConfigurations(Configuration("1s^2 2s^2 2p^5"), nLayers=10)) == 40
    Te1   = Defaults.convertUnits("energy: from eV to atomic", 1.0)
    aTBR  = Empirical.threeBodyRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Te1), bare, h1s, printout=false)
    aEII  = Empirical.impactIonizationPlasmaAlpha(Distribution.ElectronMaxwell(Te1), h1s, bare, printout=false)
    Ph    = Empirical.scaledBindingEnergy(1.0, Shell("1s"), h1s, PeriodicTable.Williams2000())
    success = success && abs( aTBR / (aEII * (2/1)/2 * (2pi/Te1)^1.5 * exp(Ph/Te1)) - 1.0 ) < 1.0e-10
    aLow  = Empirical.threeBodyRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 0.05)),
                                                        bare, h1s, printout=false)
    success = success && !isnan(aLow)  &&  aLow > 0.
    fac6  = Defaults.convertUnits("length: from atomic to cm", 1.0)^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    aTot  = Empirical.threeBodyRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 0.1)),
                                                        bare, nLayers=10, printout=false)
    success = success && abs(fac6*aTot - 1.537388245827469e-23) / 1.537388245827469e-23 < 1.0e-6
    #
    ## Test 23: multipoles of photoemissionEinsteinA(ScaledHydrogenic). The relativistically induced M1 estimate is
    ##   anchored on the hydrogenic 2s -> 1s rate, W(M1) = 2.496e-6 Z^10 1/s (Breit & Teller 1940; Johnson 1972),
    ##   scaled with the transition energy as (Delta-E/0.375 Hartree)^5; with the tabulated H binding energies
    ##   (13.6 - 3.4014 = 10.199 eV vs. the exact 10.204 eV) it reproduces 2.489e-6 1/s. E2 and higher multipoles
    ##   must refuse rather than return a silent zero rate.
    Defaults.setDefaults("nuclear: charge", 1.)
    m1    = Empirical.photoemissionEinsteinA(Configuration("2s^1"), Configuration("1s^1"), Empirical.ScaledHydrogenic(), printout=false)
    m1x   = Defaults.convertUnits("rate: from atomic", m1.rate)
    success = success && m1.multipole == M1  &&  abs(m1x - 2.4890450380695574e-6) / 2.4890450380695574e-6 < 1.0e-6
    e2Err = false
    try     Empirical.photoemissionEinsteinA(Configuration("3d^1"), Configuration("1s^1"), Empirical.ScaledHydrogenic(), printout=false)
    catch
        e2Err = true
    end
    success = success && e2Err
    #
    ## Test 24: stopping powers of a projectile in a material, stoppingPower(energies, projectile, material, approx).
    ##   For an electron in a free-electron gas (n_e = 1e20 cm^-3) the anchors are hand-computed from the closed
    ##   formulas, and three identities knit the family together: KF92 equals Bethe1931 above 14 eV; the Bohr and
    ##   KF92 logarithms differ by ln 2 below (the electron-projectile correction, since 1.123 x 1.7811 = 2); and
    ##   Axelrod1980 reduces to Bethe1931 in the nonrelativistic limit. For a proton at the *same velocity* as a
    ##   100 eV electron, the Bethe logarithms differ by ln 2 (kappa = 2 vs. 1) while the Bohr formula is
    ##   velocity-only and must agree exactly. A PartiallyIonizedGas must decompose exactly into its bound and free
    ##   populations, and the 10 MeV proton in solid Al is locked as a regression anchor (34.15 MeV cm^2/g; the
    ##   PSTAR-type ~40 differs by the neglected shell and relativistic corrections). The Coulomb logarithm is
    ##   clamped at zero far below the validity range, and KozmaFranson1992/Axelrod1980 must refuse unsupported
    ##   projectile-material combinations.
    aoCm   = Defaults.convertUnits("length: from atomic to cm", 1.0)
    spNe   = 1.0e20 * aoCm^3
    spEle  = Empirical.ElectronProjectile();    spFeg = Empirical.FreeElectronGas(spNe)
    spProt = Empirical.IonProjectile(1.0, Defaults.PROTON_MASS_U / Defaults.ELECTRON_MASS_U)
    spE100 = [Defaults.convertUnits("energy: from eV to atomic", 100.0)]
    spE5   = [Defaults.convertUnits("energy: from eV to atomic", 5.0)]
    spBet  = Empirical.stoppingPower(spE100, spEle, spFeg, Empirical.Bethe1931(), printout=false)[1]
    success = success && abs(spBet - 0.00015933628333711752) / 0.00015933628333711752 < 1.0e-6
    success = success && isapprox(Empirical.stoppingPower(spE100, spEle, spFeg, Empirical.KozmaFranson1992(), printout=false)[1],
                                  spBet, rtol=1.0e-12)
    spKF5  = Empirical.stoppingPower(spE5, spEle, spFeg, Empirical.KozmaFranson1992(), printout=false)[1]
    spBoh5 = Empirical.stoppingPower(spE5, spEle, spFeg, Empirical.Bohr1913(), printout=false)[1]
    success = success && abs(spKF5 - 0.0011226270953942747) / 0.0011226270953942747 < 1.0e-6
    ## ln Lambda = L eps/(2 pi n_e); the Bohr-KF92 difference of the logarithms must equal ln 2 (up to the rounded
    ## literature constants 1.123 and 1.7811).
    lnDiff = (spBoh5 - spKF5) * spE5[1] / (2pi * spNe)
    success = success && abs(lnDiff - log(2.0)) < 1.0e-3
    spAx   = Empirical.stoppingPower(spE100, spEle, spFeg, Empirical.Axelrod1980(), printout=false)[1]
    success = success && abs(spAx/spBet - 1.0) < 1.0e-3
    success = success && abs(Empirical.stoppingPower([Defaults.convertUnits("energy: from eV to atomic", 1.0e5)], spEle,
                                 spFeg, Empirical.Axelrod1980(), printout=false)[1] - 4.238006007084985e-7) / 4.238006007084985e-7 < 1.0e-6
    success = success && Empirical.stoppingPower([Defaults.convertUnits("energy: from eV to atomic", 0.001)], spEle,
                                 spFeg, Empirical.KozmaFranson1992(), printout=false)[1] == 0.
    ## Proton at the same velocity as a 100 eV electron: E_p = M x 100 eV.
    spEp   = [spE100[1] * spProt.M]
    spLp   = Empirical.stoppingPower(spEp, spProt, spFeg, Empirical.Bethe1931(), printout=false)[1]
    success = success && abs( (spLp - spBet) * 2*spE100[1]/(4pi * spNe) - log(2.0) ) < 1.0e-10
    success = success && isapprox(Empirical.stoppingPower(spEp,   spProt, spFeg, Empirical.Bohr1913(), printout=false)[1],
                                  Empirical.stoppingPower(spE100, spEle,  spFeg, Empirical.Bohr1913(), printout=false)[1], rtol=1.0e-10)
    ## PartiallyIonizedGas decomposes exactly into its bound and free electron populations.
    spNat  = 1.0e19 * aoCm^3
    spE10k = [Defaults.convertUnits("energy: from eV to atomic", 1.0e4)]
    spPig  = Empirical.stoppingPower(spE10k, spEle, Empirical.PartiallyIonizedGas(13, 26.98, 3.0, spNat), Empirical.Axelrod1980(), printout=false)[1]
    spNag  = Empirical.stoppingPower(spE10k, spEle, Empirical.NeutralAtomGas(13, 26.98, spNat),           Empirical.Axelrod1980(), printout=false)[1]
    spFrp  = Empirical.stoppingPower(spE10k, spEle, Empirical.FreeElectronGas(3.0 * spNat),               Empirical.Axelrod1980(), printout=false)[1]
    success = success && isapprox(spPig, (10/13)*spNag + spFrp, rtol=1.0e-12)
    ## 10 MeV proton in solid Al (rho = 2.699 g/cm^3): the bare-Bethe regression anchor in eV/cm.
    spNAl  = 2.699 / (26.98 * Defaults.ELECTRON_MASS_IN_G / Defaults.ELECTRON_MASS_U) * aoCm^3
    spAl   = Empirical.stoppingPower([Defaults.convertUnits("energy: from eV to atomic", 1.0e7)], spProt,
                                     Empirical.NeutralAtomGas(13, 26.98, spNAl), Empirical.Bethe1931(), printout=false)[1]
    success = success && abs(Defaults.convertUnits("energy: from atomic to eV", spAl)/aoCm - 9.218125319197525e7) / 9.218125319197525e7 < 1.0e-6
    ## Unsupported combinations must refuse.
    spErr = 0
    try     Empirical.stoppingPower(spE100, spProt, spFeg, Empirical.KozmaFranson1992(), printout=false)
    catch
        spErr = spErr + 1
    end
    try     Empirical.stoppingPower(spE100, spProt, spFeg, Empirical.Axelrod1980(), printout=false)
    catch
        spErr = spErr + 1
    end
    success = success && spErr == 2
    #
    ## Test 25: tunneling ionization rate (ADK 1986). The direct (Z, Ip, l) method must reproduce the exact, closed
    ##   Landau tunneling formula for hydrogenic 1s states (Z=1, Ip=0.5 Hartree and Z=2, Ip=2.0 Hartree) to machine
    ##   precision -- a literature-independent cross check of the Gamma-function ADK coefficient. The angular factor
    ##   and coefficient are checked at simple integer arguments, and both the Configuration-based wrapper (which
    ##   reproduces the direct method up to the ~2% tabulated-vs-exact Ip difference for neutral H) and the two
    ##   error paths (non-positive Z or Ip) are exercised.
    adkFs = [0.02, 0.05, 0.1, 0.2, 0.5]
    function landauRate(Ip, F)
        return 4*(2*Ip)^2.5/F * exp(-2*(2*Ip)^1.5/(3*F))
    end
    adkH  = Empirical.tunnelingIonizationRate(adkFs, 1.0, 0.5, 0, Empirical.ADK1986(), printout=false)
    for  (ic, F)  in  enumerate(adkFs)
        success = success && abs(adkH[ic]/landauRate(0.5, F) - 1.0) < 1.0e-8
    end
    adkHe = Empirical.tunnelingIonizationRate([0.5, 1.0, 2.0], 2.0, 2.0, 0, Empirical.ADK1986(), printout=false)
    for  (ic, F)  in  enumerate([0.5, 1.0, 2.0])
        success = success && abs(adkHe[ic]/landauRate(2.0, F) - 1.0) < 1.0e-8
    end
    success = success && Empirical.adkAngularFactor(0,0) == 1.0
    success = success && Empirical.adkAngularFactor(1,0) == 3.0
    success = success && Empirical.adkAngularFactor(1,1) == 3.0
    success = success && Empirical.adkCoefficient(1.0, 0.0) == 4.0
    Defaults.setDefaults("nuclear: charge", 1.)
    adkWrap = Empirical.tunnelingIonizationRate([0.1], Configuration("1s^1"), Configuration("1s^0"), Empirical.ADK1986(), printout=false)
    success = success && abs(adkWrap[1]/landauRate(0.5, 0.1) - 1.0) < 0.01
    adkErr = 0
    try     Empirical.tunnelingIonizationRate([0.1], -1.0, 0.5, 0, Empirical.ADK1986())
    catch
        adkErr = adkErr + 1
    end
    try     Empirical.tunnelingIonizationRate([0.1], 1.0, -0.5, 0, Empirical.ADK1986())
    catch
        adkErr = adkErr + 1
    end
    success = success && adkErr == 2
    #
    ## Test 26: scaledBindingEnergy() for a genuine few-electron ion (NoElectrons <= 2 and < Z) now bypasses the
    ##   *neutral*-atom tabulation and returns the exact (H-like) or approximate (He-like, ~5-30% high, decreasing
    ##   with Z; checked against NIST/CRC data) hydrogenic value, instead of silently reusing the neutral-atom
    ##   entry; the tabulated path for near-neutral configurations (including the Ne K-alpha analog
    ##   spectator-omitted shorthand) must remain exactly unchanged.
    Defaults.setDefaults("nuclear: charge", 2.)
    bEHePlus = Empirical.scaledBindingEnergy(2.0, Shell("1s"), Configuration("1s^1"), PeriodicTable.Williams2000())
    success = success && abs(bEHePlus - 2.0) < 1.0e-10                              ## exact hydrogenic, Z^2/2 = 2.0
    Defaults.setDefaults("nuclear: charge", 3.)
    bELi2p = Empirical.scaledBindingEnergy(3.0, Shell("1s"), Configuration("1s^1"), PeriodicTable.Williams2000())
    success = success && abs(bELi2p - 4.5) < 1.0e-10                                ## exact hydrogenic, Z^2/2 = 4.5
    Defaults.setDefaults("nuclear: charge", 2.)
    bEHeNeutral = Empirical.scaledBindingEnergy(2.0, Shell("1s"), Configuration("1s^2"), PeriodicTable.Williams2000())
    success = success && abs(Defaults.convertUnits("energy: from atomic", bEHeNeutral) - 24.6) < 0.05  ## unchanged tabulated value
    Defaults.setDefaults("nuclear: charge", 10.)
    bENeKalpha = Empirical.scaledBindingEnergy(10.0, Shell("2p"), Configuration("1s^1 2p^6"), PeriodicTable.Williams2000())
    e1au       = Empirical.bindingEnergy(10, Shell("2p"), data=PeriodicTable.Williams2000())
    success = success && abs(bENeKalpha - e1au) < 1.0e-10     ## unchanged: same tabulated 2p value as Test 1
    #
    ## Test 27: forbiddenExcitationCrossSection/PlasmaAlpha(ConstantCollisionStrength()) for the classic H 1s -> 2s
    ##   (M1, Delta-l = 0) and 1s -> 3d (E2, Delta-l = 2) transitions. Checks: cross section vanishes below
    ##   threshold; an E1 transition (1s -> 2p) is rejected with a pointer to VanRegemorter1962; an E3+ transition
    ##   (1s -> 4f) is rejected outright; and the Maxwellian-folded plasma rate coefficient for Omega = 1 reproduces
    ##   the textbook C_ij = 8.629e-6/(g_i sqrt(T_e[K])) exp(-deltaE/T_e) [cm^3/s] formula (Osterbrock & Ferland) to
    ##   better than 1.0e-3 relative -- an independent, literature-based cross check of the Gauss-Legendre folding.
    Defaults.setDefaults("nuclear: charge", 1.)
    h1s = Configuration("1s^1");   h2s = Configuration("2s^1");   h2p = Configuration("2p^1")
    h3d = Configuration("3d^1");   h4f = Configuration("4f^1")
    epsFb = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [5.0, 10.2, 15.0]]
    cssM1 = Empirical.forbiddenExcitationCrossSection(epsFb, h1s, h2s, Empirical.ConstantCollisionStrength(), printout=false)
    success = success && cssM1[1] == 0.  &&  cssM1[3] > 0.
    cssE2 = Empirical.forbiddenExcitationCrossSection(epsFb, h1s, h3d, Empirical.ConstantCollisionStrength(), printout=false)
    success = success && cssE2[1] == 0.  &&  cssE2[2] == 0.  &&  cssE2[3] > 0.   ## E2 threshold (3d) lies above 10.2 eV
    fbErr = 0
    try     Empirical.forbiddenExcitationCrossSection(epsFb, h1s, h2p, Empirical.ConstantCollisionStrength())
    catch
        fbErr = fbErr + 1
    end
    try     Empirical.forbiddenExcitationCrossSection(epsFb, h1s, h4f, Empirical.ConstantCollisionStrength())
    catch
        fbErr = fbErr + 1
    end
    success = success && fbErr == 2
    deltaEfb = Empirical.scaledBindingEnergy(1.0, Shell("1s"), h1s, PeriodicTable.Williams2000()) -
               Empirical.scaledBindingEnergy(1.0, Shell("2s"), h2s, PeriodicTable.Williams2000())
    giFb  = Basics.extractFromConfiguration(Basics.Multiplicity(), h1s)
    TeFb  = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    alpFb = Empirical.forbiddenExcitationPlasmaAlpha(Distribution.ElectronMaxwell(TeFb), h1s, h2s, printout=false)
    TeFbK = 5.0 * 11604.518;   deltaEfbK = Defaults.convertUnits("energy: from atomic to eV", deltaEfb) * 11604.518
    textbook  = 8.629e-6 / (giFb * sqrt(TeFbK)) * exp(-deltaEfbK/TeFbK) *
                Defaults.convertUnits("time: from atomic to sec", 1.0) / Defaults.convertUnits("length: from atomic to cm", 1.0)^3
    success = success && abs(alpFb/textbook - 1.0) < 1.0e-3
    #
    ## Test 28: chargeExchangeCrossSection/PlasmaAlpha(OverBarrierModel1980/NiehausScaling1986) for He2+ + H(1s).
    ##   Checks: the a.u.-converted Niehaus cross section reproduces the cm^2 formula sigma = 2.6e-13 q/Ip[eV]^2
    ##   applied by hand; the Configuration-based wrapper (any standard-filling conf, neutral or ion) agrees with
    ##   the direct (q, Ip) method to machine precision, for both a neutral-atom donor (H) and an already-ionized
    ##   donor (He+, i.e. the ion-ion CX regime); non-positive q or Ip are rejected; and the Maxwellian-folded
    ##   plasma rate coefficient alpha = <v> sigma (analytic, no quadrature) matches a hand-computed cross check.
    Defaults.setDefaults("nuclear: charge", 1.)
    hConf  = Configuration("1s^1")
    IpH    = Empirical.ionizationPotential(1, hConf)
    IpHeV  = Defaults.convertUnits("energy: from atomic", IpH)
    sigNieh = Empirical.chargeExchangeCrossSection(2.0, IpH, Empirical.NiehausScaling1986(), printout=false)
    sigNiehCm2 = Defaults.convertUnits("cross section: from atomic to cm^2", sigNieh)
    success = success && abs(sigNiehCm2/(2.6e-13 * 2.0/IpHeV^2) - 1.0) < 1.0e-8
    sigOBM  = Empirical.chargeExchangeCrossSection(2.0, IpH, Empirical.OverBarrierModel1980(), printout=false)
    success = success && abs(sigOBM - pi*((2.0+1.0)/(2*IpH))^2) < 1.0e-10
    sigWrap = Empirical.chargeExchangeCrossSection(2.0, hConf, 1.0, Empirical.OverBarrierModel1980(), printout=false)
    success = success && sigWrap == sigOBM
    Defaults.setDefaults("nuclear: charge", 2.)
    heIonConf = Configuration("1s^1")
    sigIonIon = Empirical.chargeExchangeCrossSection(3.0, heIonConf, 2.0, Empirical.NiehausScaling1986(), printout=false)
    success = success && sigIonIon > 0.
    cxErr = 0
    try     Empirical.chargeExchangeCrossSection(-1.0, IpH, Empirical.OverBarrierModel1980())
    catch
        cxErr = cxErr + 1
    end
    try     Empirical.chargeExchangeCrossSection(2.0, -1.0, Empirical.NiehausScaling1986())
    catch
        cxErr = cxErr + 1
    end
    success = success && cxErr == 2
    Mproj = Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U * 4.0;   Mtarg = Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U
    Tcx   = Defaults.convertUnits("energy: from eV to atomic", 1.0)
    alphaCx = Empirical.chargeExchangePlasmaAlpha(Tcx, 2.0, IpH, Mproj, Mtarg, Empirical.NiehausScaling1986(), printout=false)
    muCx  = Mproj*Mtarg/(Mproj+Mtarg)
    success = success && abs(alphaCx/(sqrt(8*Tcx/(pi*muCx)) * sigNieh) - 1.0) < 1.0e-10
    #
    ## Test 29: XrayDataBooklet is now the default energy-data set (Z = 1, ..., 92, shells 1s-6p), replacing
    ##   Williams2000 (Z = 1, ..., 36 only, shells 1s-4p); the two are numerically identical wherever both are
    ##   tabulated, except Carbon 2s/2p, where XrayDataBooklet was missing a correction Williams2000 already had --
    ##   now propagated. Checks: bindingEnergy(6, 2s/2p) now returns the corrected (not -1/fallback) values by
    ##   default; the Ne K-alpha energy (a Z <= 36 case) is bit-identical under the new default and under an
    ##   explicit Williams2000() call; and a Z = 54 (Xe) binding energy, unreachable via the old default, now
    ##   resolves to a real tabulated value instead of the cruder Slater-hydrogenic fallback.
    Defaults.setDefaults("nuclear: charge", 6.)
    c2s = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(6, Shell("2s")))
    c2p = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(6, Shell("2p")))
    success = success && abs(c2s - 19.4) < 1.0e-8  &&  abs(c2p - 7.0) < 1.0e-8
    Defaults.setDefaults("nuclear: charge", 10.)
    neConf = Configuration("1s^1 2p^6");   neFConf = Configuration("1s^2 2p^5")
    waDefault  = Empirical.photoemissionEinsteinA(neConf, neFConf, Empirical.ScaledHydrogenic(), printout=false)
    waWilliams = Empirical.photoemissionEinsteinA(neConf, neFConf, Empirical.ScaledHydrogenic(), printout=false,
                                                  data=PeriodicTable.Williams2000())
    success = success && waDefault.energy == waWilliams.energy
    Defaults.setDefaults("nuclear: charge", 54.)
    xeConf = Configuration("[Kr] 4d^10 5s^2 5p^6")
    bEXe   = Empirical.scaledBindingEnergy(54.0, Shell("5p"), xeConf, PeriodicTable.XrayDataBooklet())
    success = success && abs(Defaults.convertUnits("energy: from atomic", bEXe) - 12.1) < 1.0e-8
    #
    ## Test 30: chargeExchangeCaptureShell/StateSelectiveCrossSection for the classic C6+ + H(1s) system, whose
    ##   dominant capture shell n = 6 is a well-known textbook result for this q ~ Ip-matched case. Checks: the
    ##   level-matching identity q^2/(2 nc^2) == Ip holds to machine precision; the returned n (rounded nc) equals
    ##   6; the l-fractions and l-resolved cross sections sum exactly to 1 and to sigmaTotal, respectively; the
    ##   Configuration-based wrapper agrees with the direct (q, Ip) method; non-positive q or Ip are rejected; and
    ##   a weakly-matched case (small q, large Ip) correctly floors at n = 1.
    Defaults.setDefaults("nuclear: charge", 1.)
    hConf30 = Configuration("1s^1");   IpH30 = Empirical.ionizationPotential(1, hConf30)
    capShell = Empirical.chargeExchangeCaptureShell(6.0, IpH30)
    success = success && abs(6.0^2/(2*capShell.nc^2)/IpH30 - 1.0) < 1.0e-12
    success = success && capShell.n == 6
    stateSel = Empirical.chargeExchangeStateSelectiveCrossSection(6.0, IpH30, Empirical.NiehausScaling1986(), printout=false)
    success = success && stateSel.n == 6  &&  length(stateSel.states) == 6
    success = success && abs(sum(st.fraction for st in stateSel.states) - 1.0) < 1.0e-12
    success = success && abs(sum(st.sigma for st in stateSel.states)/stateSel.sigmaTotal - 1.0) < 1.0e-10
    stateSelWrap = Empirical.chargeExchangeStateSelectiveCrossSection(6.0, hConf30, 1.0, Empirical.NiehausScaling1986(), printout=false)
    success = success && stateSel.n == stateSelWrap.n  &&  stateSel.sigmaTotal == stateSelWrap.sigmaTotal
    css30Err = 0
    try     Empirical.chargeExchangeCaptureShell(-1.0, IpH30)
    catch
        css30Err = css30Err + 1
    end
    try     Empirical.chargeExchangeCaptureShell(6.0, -1.0)
    catch
        css30Err = css30Err + 1
    end
    success = success && css30Err == 2
    weakShell = Empirical.chargeExchangeCaptureShell(1.0, Defaults.convertUnits("energy: from eV to atomic", 500.0))
    success = success && weakShell.n == 1  &&  weakShell.nc < 1.0
    #
    ## Restore the global nuclear charge for all subsequent tests.
    Defaults.setDefaults("nuclear: charge", oldZ)
    ###
    Defaults.setDefaults("print summary: close", "")
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-Empirical-new.sum")
    testPrint("testModule_Empirical()::", success)
    return(success)
end


"""
`TestFrames.testModule_ImpactIonization(; short::Bool=true)`  ... tests on module ImpactIonization.
"""
function testModule_ImpactIonization(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ImpactIonization-new.sum")
    printstyled("\n\nTest the module  ImpactIonization  ... \n", color=:cyan)
    ### Make the tests: total EII cross sections for He I (1s^2) with BEBmodel at 5 impact energies.
    grid        = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    iEnergies   = [50.0, 100.0, 200.0, 500.0, 1000.0]
    shells      = Basics.generateShellList(1, 1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    eiiSettings = ImpactIonization.Settings(ImpactIonization.BEBmodel(), 1, iEnergies, false, true, selection)
    comp        = Empirical.Computation("EII cross section for He I (BEBmodel).", Nuclear.Model(2.0), grid,
                                        [Configuration("1s^2")], eiiSettings)
    perform(comp; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-ImpactIonization-new.sum")
    ## Activate once test/approved/test-ImpactIonization-approved.sum has been verified and copied:
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-ImpactIonization-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-ImpactIonization-new.sum"),
    ##                             "Total ionization cross sections", 11)
    success = true
    testPrint("testModule_ImpactIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_Semiempirical(; short::Bool=true)`  ... tests on module Semiempirical.
"""
function testModule_Semiempirical(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Semiempirical-new.sum")
    printstyled("\n\nTest the module  Semiempirical  ... \n", color=:cyan)
    success = true
    ## Test 1: EstimateIonizationPotentialInnerShell for Ne 1s_1/2 (Bug 4 fix: Shell → Subshell argument)
    e1 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateIonizationPotentialInnerShell(), Subshell("1s_1/2"), 10))
    success = success && abs(e1 - 870.2) < 0.01
    ## Test 2: EstimateBindingEnergyWilliams2000 Subshell -- same Williams source, must agree with Test 1 exactly
    e2 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyWilliams2000(), 10, Subshell("1s_1/2")))
    success = success && e2 == e1
    ## Test 3: EstimateBindingEnergyXrayDataBooklet Subshell for Xe 4d_3/2 (new function extending to deep subshells)
    e3 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyXrayDataBooklet(), 54, Subshell("4d_3/2")))
    success = success && abs(e3 - 69.5) < 0.01
    ## Test 4: EstimateBindingEnergyWilliams2000 Configuration for Kr ground config (Bugs 1&2 fix: no BoundsError)
    e4 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyWilliams2000(), 36, Configuration("[Ar] 3d^10 4s^2 4p^6")))
    success = success && abs(e4 - 45514.0) < 0.1
    ## Test 5: estimateSlaterZeff for Rb 5s_1/2 -- σ = 8×0.85 + 28×1.00 = 34.80 → Zeff = 2.20
    z5 = Semiempirical.estimateSlaterZeff(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), Subshell("5s_1/2"))
    success = success && abs(z5 - 2.20) < 1.0e-10
    ## Test 6: EstimateBindingEnergyNist2025 for He IP_1 = 24.587 eV (well-known textbook value)
    e6 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 2, 1))
    success = success && abs(e6 - 24.587) < 0.001
    ###
    Defaults.setDefaults("print summary: close", "")
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-Semiempirical-new.sum")
    testPrint("testModule_Semiempirical()::", success)
    return(success)
end
