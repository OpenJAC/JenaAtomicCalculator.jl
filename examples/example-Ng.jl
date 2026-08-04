#
println("Ng) Tests of the (quasiclassical) tunneling ionization rate (ADK1986), module Empirical.")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate",   "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  17-Jul-2026
    # Branch 1: Empirical.tunnelingIonizationRate(fields, Z, Ip, l, ADK1986()) -- the direct (Z, Ip, l) method,
    #   validated against the exact, closed Landau tunneling formula for hydrogenic 1s states.
    # Systems: H (Z=1, Ip=0.5 Hartree) and He^+ (Z=2, Ip=2.0 Hartree), both l=0, m=0.
    # Note: for a pure hydrogenic 1s state (l*=0, n*=1 by construction) the Gamma-function ADK coefficient
    #   C^2_(1,0) = 4 exactly, and the ADK rate then reduces *exactly* to Landau's closed formula -- a
    #   literature-independent cross check of the normalization used here (verified to machine precision).
    # Checks:
    #   - ADK/Landau = 1 to 8+ digits for both H and He^+ across the tabulated field range.
    #   - The rate rises by many orders of magnitude as F approaches and passes the critical field F_c = kappa^3/16
    #     (H: F_c = 0.0625 a.u.; He^+: F_c = 0.5 a.u.): the tunneling barrier progressively thins.
    #
    println("\n  Empirical.tunnelingIonizationRate(fields, Z, Ip, l, ADK1986()) vs. the exact Landau formula:\n")
    landauRate(Ip, F) = 4*(2*Ip)^2.5/F * exp(-2*(2*Ip)^1.5/(3*F))
    fieldsEV = [0.02, 0.05, 0.1, 0.2, 0.5]
    ratesH   = Empirical.tunnelingIonizationRate(fieldsEV, 1.0, 0.5, 0, Empirical.ADK1986(); printout=true)
    ratesHe  = Empirical.tunnelingIonizationRate(fieldsEV, 2.0, 2.0, 0, Empirical.ADK1986())
    println("    F [a.u.]      ADK(H) [1/s]        Landau(H) [1/s]     ratio      ADK(He+)/Landau(He+)")
    for  (i, F)  in  enumerate(fieldsEV)
        lH  = landauRate(0.5, F);   lHe = landauRate(2.0, F)
        ##  ratesH/lH is a ratio of two atomic-unit rates and needs no conversion; the magnitudes shown for context
        ##  are converted to the user-defined rate unit (1/s), exactly as the library's own printout above does.
        adkHx = Defaults.convertUnits("rate: from atomic", ratesH[i]);   lHx = Defaults.convertUnits("rate: from atomic", lH)
        println("    $(rpad(F,12))  $(rpad(round(adkHx, sigdigits=5),19)) $(rpad(round(lHx, sigdigits=5),19)) " *
                "$(rpad(round(ratesH[i]/lH, digits=8),10)) $(round(ratesHe[i]/lHe, digits=8))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 2: Empirical.tunnelingIonizationRate(fields, iConf, fConf, ADK1986()) -- the Configuration-based
    #   convenience wrapper, now fixed for H-like/He-like ions, and the narrower residual gap that remains.
    # Systems: (a) neutral H (near-exact: the neutral-atom table entry IS the hydrogenic value here);
    #   (b) neutral Ne 2p (l=1) at a field derived from a typical strong-field laser intensity;
    #   (c) He^+ -> He^2+: Empirical.scaledBindingEnergy() now recognizes iConf = "1s^1" at Z=2 as a genuine
    #   H-like ion (NoElectrons=1 < Z, <=2) and returns the exact hydrogenic 54.4 eV automatically, so the
    #   wrapper and the direct (Z, Ip, l) method now agree; a residual gap remains only for stripped ions with
    #   3 or more electrons (e.g. Li-like O^5+), which are indistinguishable from a spectator-omitted shorthand
    #   by electron count alone and still fall back to the *neutral*-element tabulated value.
    # Checks:
    #   - (a) at F=0.1 a.u. matches the direct-method (exact Landau) result of Branch 1 to 0.45% (the small
    #     difference between the tabulated 13.6 eV and the exact 13.6057 eV, amplified through the exponential;
    #     the mismatch grows to ~2% at the smaller F=0.02 a.u. tested in Branch 1).
    #   - (c) the wrapper and the direct method now agree to machine precision at F=1 a.u. (both use Ip=54.4 eV).
    #
    println("\n  Empirical.tunnelingIonizationRate(fields, iConf, fConf, ADK1986()):\n")
    ##  (a) neutral H
    setDefaults("nuclear: charge", 1.)
    rH = Empirical.tunnelingIonizationRate([0.1], Configuration("1s^1"), Configuration("1s^0"), Empirical.ADK1986(); printout=true)
    ##  (b) neutral Ne 2p at an intensity-derived field
    setDefaults("nuclear: charge", 10.)
    F1 = Empirical.electricFieldFromIntensity(1.0e15)    ##  a typical strong-field/HHG intensity [W/cm^2]
    rNe = Empirical.tunnelingIonizationRate([F1], Configuration("1s^2 2s^2 2p^6"), Configuration("1s^2 2s^2 2p^5"),
                                            Empirical.ADK1986())
    println("  Ne 2p at I = 1e15 W/cm^2 (F = $(round(F1, digits=4)) a.u.): rate = $(round(rNe[1], sigdigits=4)) [1/s]")
    ##  (c) He+ -> He2+: now fixed -- wrapper and direct method agree
    setDefaults("nuclear: charge", 2.)
    rHeWrap = Empirical.tunnelingIonizationRate([1.0], Configuration("1s^1"), Configuration("1s^0"), Empirical.ADK1986(); printout=true)
    rHeDirect = Empirical.tunnelingIonizationRate([1.0], 2.0, 2.0, 0, Empirical.ADK1986())
    println("  He+ -> He2+ at F=1 a.u.:  wrapper (fixed, Ip=54.4 eV) = $(round(rHeWrap[1], sigdigits=4))," *
            "  direct (Ip=54.4 eV) = $(round(rHeDirect[1], sigdigits=4))  [1/s]  (now identical)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  17-Jul-2026
    # Branch 3: the tunneling regime and its boundary -- rate vs. field strength across F_c.
    #   System: neutral Ar (Z=1 residual, Ip = tabulated 3p binding energy), scanning the field from well inside
    #   the tunneling regime (F << F_c) to well above the critical (barrier-suppression) field F_c = kappa^3/16.
    # Checks:
    #   - The rate rises approximately as exp(-2 kappa^3/(3F)): a log-rate plot is close to linear in 1/F for
    #     F << F_c.
    #   - No error is raised for F > F_c (over-barrier regime), but the printout flags F/F_c > 1: ADK is known to
    #     substantially *overestimate* the true rate there (the tunneling picture itself breaks down).
    #
    println("\n  Tunneling ionization rate of Ar 3p across the critical field F_c:\n")
    setDefaults("nuclear: charge", 18.)
    iC = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6");   fC = Configuration("1s^2 2s^2 2p^6 3s^2 3p^5")
    Fs = [0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2]
    rates = Empirical.tunnelingIonizationRate(Fs, iC, fC, Empirical.ADK1986(); printout=true)
    println("    F [a.u.]    rate [1/s]")
    for  (i, F)  in  enumerate(Fs)
        ratex = Defaults.convertUnits("rate: from atomic", rates[i])
        println("    $(rpad(F,11)) $(round(ratex, sigdigits=4))")
    end
    #
    setDefaults("print summary: close", "")
    #
end
