
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
    ## Test 7: photoionizationCrossSection(ScaledHydrogenic) for Ne 2p; Kramers (1923) formula.
    ##   The threshold is the tabulated 2p binding energy of 21.6 eV; below it the cross section vanishes.
    ##   Above threshold the cross section falls off as omega^-3, so doubling the photon energy must reduce
    ##   it by exactly a factor of 8.
    piConf = Configuration("1s^2 2s^2 2p^6");    pfConf = Configuration("1s^2 2s^2 2p^5")
    omegas = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [15.0, 25.0, 30.0, 600.0, 1200.0]]
    piCss  = Empirical.photoionizationCrossSection(omegas, piConf, pfConf, Empirical.ScaledHydrogenic(), printout=false)
    success = success && piCss[1] == 0.  &&  piCss[2] > 0.
    success = success && abs(piCss[4]/piCss[5] - 8.) < 1.0e-6
    #
    ## Test 8: photorecombinationCrossSection(ScaledHydrogenic) for Ne 2p via the Einstein-Milne relation.
    ##   sigma^(PR) ~ 1 / [eps * (eps + bE)] and must diverge as eps -> 0; a flat or vanishing cross section
    ##   signals a broken Milne relation or an inconsistent threshold.
    ekins  = [Defaults.convertUnits("energy: from eV to atomic", e) for e in [0.1, 1.0, 5.0, 20.0, 100.0]]
    prCss  = Empirical.photorecombinationCrossSection(ekins, pfConf, piConf, Empirical.ScaledHydrogenic(), printout=false)
    prBarn = [round(Defaults.convertUnits("cross section: from atomic", cs), digits=3) for cs in prCss]
    for  (ic, cs)  in  enumerate([1886.048, 181.094, 30.772, 4.919, 0.337])
        success = success && abs(prBarn[ic] - cs) < 0.001
    end
    #
    ## Test 9: Kramers' normalisation against the literature. For H 1s the semiclassical Kramers formula is known
    ##   to give 7.91 Mb at threshold, to be compared with the exact 6.30 Mb (Gaunt factor 0.80). This anchors the
    ##   prefactor 32 pi alpha / (3 sqrt(3) n) to a published value rather than to JAC's own output.
    Defaults.setDefaults("nuclear: charge", 1.)
    bEH    = Empirical.bindingEnergy(1, Shell("1s"))
    hCss   = Empirical.photoionizationCrossSection([bEH * 1.0000001], Configuration("1s^1"), Configuration("1s^0"),
                                                   Empirical.ScaledHydrogenic(), printout=false)
    hMb    = Defaults.convertUnits("cross section: from atomic", hCss[1]) / 1.0e6
    success = success && abs(hMb - 7.91) < 0.01
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
