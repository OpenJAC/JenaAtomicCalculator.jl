#
println("Nb) Tests of empirical Einstein-A values and photoemission/excitation/deexcitation plasma rates.")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate",   "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  15-Jul-2026
    # Branch 1: Empirical.photoemissionEinsteinA(iConf, fConf, ScaledHydrogenic()) --
    #   transition energy from the tabulated binding energies; hydrogenic scaling of the E1 rate.
    # Systems: Ne (Z=10) and Mg (Z=12), K-alpha analogue: 1s-hole -> 2p-hole.
    # Note: the energy is now the difference of two tabulated binding energies and is therefore accurate,
    #   while the rate still rests on hydrogenic r-expectation values and is only good to a factor ~1.5.
    # Checks:
    #   - multipole should be E1 (|Delta l| = 1 between 1s and 2p).
    #   - energy > 0 (emitted photon); Ne = 870.2 - 21.6 = 848.6 eV, Mg = 1303.0 - 49.5 = 1253.5 eV,
    #     to be compared with the literature K-alpha energies of 848.6 eV and 1253.6 eV.
    #   - rate > 0 for E1 transitions; Ne ~ 6.0e13 1/s vs. the literature 8.8e13 1/s.
    #
    println("\n  Empirical.photoemissionEinsteinA(iConf, fConf, ScaledHydrogenic()): \n")
    ##  Ne (Z=10): K-alpha analogue  1s^1 2p^6  ->  1s^2 2p^5
    setDefaults("nuclear: charge", 10.)
    wa = Empirical.photoemissionEinsteinA(Configuration("1s^1 2p^6"), Configuration("1s^2 2p^5"), Empirical.ScaledHydrogenic(), printout=true)
    wb = round(Defaults.convertUnits("energy: from atomic", wa.energy), digits=1)
    wc = round(Defaults.convertUnits("rate: from atomic",   wa.rate),   digits=3)
    println("  Ne K-alpha (Z=10):  multipole=$(wa.multipole),  energy=$wb eV,  rate=$wc 1/s")
    ##  Mg (Z=12): K-alpha analogue  1s^1 2p^6 3s^2  ->  1s^2 2p^5 3s^2
    setDefaults("nuclear: charge", 12.)
    wa = Empirical.photoemissionEinsteinA(Configuration("1s^1 2p^6 3s^2"), Configuration("1s^2 2p^5 3s^2"), Empirical.ScaledHydrogenic(), printout=true)
    wb = round(Defaults.convertUnits("energy: from atomic", wa.energy), digits=1)
    wc = round(Defaults.convertUnits("rate: from atomic",   wa.rate),   digits=3)
    println("  Mg K-alpha (Z=12):  multipole=$(wa.multipole),  energy=$wb eV,  rate=$wc 1/s")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 2: Empirical.photoemissionEinsteinA(iConf, fConf, UsingJAC()) --
    #   mean-field computation for the transition energy and E1/M1 rate.
    # System: Ne (Z=10), K-alpha: 1s^1 2p^6 -> 1s^2 2p^5.
    # Checks:
    #   - Energy should agree with Branch 1 to within ~30% (both are approximations).
    #   - Rate from UsingJAC should be closer to literature (Ne K-alpha A ~ 8.8e13 s^-1).
    #
    println("\n  Empirical.photoemissionEinsteinA(iConf, fConf, UsingJAC()): \n")
    Z     = 10.;  setDefaults("nuclear: charge", Z)
    iConf = Configuration("1s^1 2p^6");   fConf = Configuration("1s^2 2p^5")
    res   = Empirical.photoemissionEinsteinA(iConf, fConf, Empirical.UsingJAC(), printout=true)
    eVx   = round(Defaults.convertUnits("energy: from atomic", res.energy), digits=1)
    ratx  = round(Defaults.convertUnits("rate: from atomic",   res.rate),   digits=3)
    println("  Ne K-alpha (Z=$Z):  multipole=$(res.multipole),  energy=$eVx eV,  rate=$ratx 1/s")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 3: Empirical.photoexcitationPlasmaRatePerIon and photodeexcitationPlasmaRatePerIon --
    #   plasma rates at a given temperature using a Planck photon distribution.
    # System: Ne K-alpha transition (848.6 eV); temperature T = 1 keV (typical coronal plasma).
    # Note: with omega/kT = 0.849 the mean occupation number nbar = 1/(exp(0.849)-1) = 0.748 is of order
    #   unity, so that the stimulated term is *not* negligible here: R^PD = A (1 + nbar) = 1.75 A.
    # Checks:
    #   - R^PD / A = 1 + nbar = 1.748  (A = 6.034e13 1/s, R^PD = 1.055e14 1/s).
    #   - R^PX = (g_i/g_f) A nbar = (2/6) x 6.034e13 x 0.748 = 1.504e13 1/s.
    #   - Detailed balance: R^PD / R^PX = (g_f/g_i) exp(omega/kT) = 3 x 2.336 = 7.011.
    #
    println("\n  Empirical.photoexcitation/deexcitation plasma rates for Ne K-alpha:\n")
    Z     = 10.;  setDefaults("nuclear: charge", Z)
    iConf = Configuration("1s^1 2p^6");   fConf = Configuration("1s^2 2p^5")
    T     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.16e7)    ##  T = 1 keV
    dist  = Distribution.PhotonPlanck(T)
    rPX   = Empirical.photoexcitationPlasmaRatePerIon(  dist, iConf, fConf, approx=Empirical.ScaledHydrogenic(), printout=true)
    rPD   = Empirical.photodeexcitationPlasmaRatePerIon(dist, iConf, fConf, approx=Empirical.ScaledHydrogenic(), printout=true)
    rPXx  = round(Defaults.convertUnits("rate: from atomic", rPX), digits=3)
    rPDx  = round(Defaults.convertUnits("rate: from atomic", rPD), digits=3)
    println("  R^PX (per ion) = $rPXx 1/s")
    println("  R^PD (total)   = $rPDx 1/s")
    #
    setDefaults("print summary: close", "")
    #
end
