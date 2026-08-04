#
println("Hd) Tests of the K-shell fluorescence and Auger yields (KrauseAdopted2016) and the empirical Auger rates.")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate",   "1/s")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  16-Jul-2026
    # Branch 1: Empirical.fluorescenceYield(Z) and Empirical.augerYield(Z) -- the KrauseAdopted2016 data set.
    # Data: K-shell fluorescence yields omega_K for Z = 3 ... 100, adopted from the xraylib project
    #   (K-shell revision of 2016); base data ORNL-5399 (Krause et al. 1978). Auger yield a_K = 1 - omega_K.
    # Checks:
    #   - omega_K rises monotonically with Z: light elements decay by Auger (Ne: a_K = 0.985), heavy
    #     elements by fluorescence (Au: a_K = 0.040).
    #   - The crossing omega_K = a_K = 0.5 falls between Z = 30 (Zn, 0.486) and Z = 31 (Ga, 0.517).
    #   - Z = 1, 2 raise an error: a K-shell vacancy of H/He^+ cannot decay radiatively (no L electron).
    #
    println("\n  K-shell fluorescence and Auger yields (KrauseAdopted2016):\n")
    println("    Z    element    omega_K      a_K = 1 - omega_K")
    for  (Z, el)  in  [(5,"B"), (10,"Ne"), (14,"Si"), (18,"Ar"), (26,"Fe"), (30,"Zn"), (31,"Ga"),
                       (36,"Kr"), (47,"Ag"), (54,"Xe"), (74,"W"), (79,"Au"), (92,"U")]
        wK = Empirical.fluorescenceYield(Z);    aK = Empirical.augerYield(Z)
        println("    $Z    $el         $(rpad(wK, 10))   $(round(aK, digits=4))")
    end
    ##  The error path for Z = 2
    try
        Empirical.fluorescenceYield(2)
    catch
        println("\n  Z = 2 correctly raises an error: no radiative K-vacancy decay without an L electron.")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 2: Empirical.augerRate(iConf, fConf, ScaledHydrogenic()) -- empirical K-shell Auger rates
    #   from A_auger = A_rad (1 - omega_K)/omega_K, with A_rad the empirical K-alpha rate.
    # Physics: K-shell Auger rates are nearly independent of Z along an isoelectronic sequence, whereas
    #   radiative rates grow ~Z^4; the whole Z-dependence of the branching is carried by omega_K.
    # Note: the branching ratio A_auger/A_rad = (1 - omega_K)/omega_K is exact by construction, but the
    #   *absolute* Auger rate inherits the ~1.5x uncertainty of the ScaledHydrogenic A_rad, amplified by
    #   the large (1 - omega_K)/omega_K at low Z; e.g., the implied Ne K width is ~2.6 eV vs. a measured
    #   ~0.27 eV. Read the rates as order-of-magnitude estimates, most reliable near omega_K ~ 0.5.
    # Checks:
    #   - A_auger stays within one order of magnitude (~1e15 ... 1e16 1/s) from Ne to Fe, while A_rad
    #     grows by ~2 orders -- the expected near-Z-independence of the Auger rate.
    #   - A_auger/A_rad equals (1 - omega_K)/omega_K for every Z.
    #
    println("\n  Empirical.augerRate(iConf, fConf, ScaledHydrogenic()) for K-shell holes:\n")
    println("    Z    element    A_rad [1/s]      A_auger [1/s]    A_auger/A_rad    (1-w)/w")
    for  (Z, el, iC, fC)  in  [
            (10., "Ne", Configuration("1s^1 2p^6"),                Configuration("1s^2 2p^5")),
            (12., "Mg", Configuration("1s^1 2p^6 3s^2"),           Configuration("1s^2 2p^5 3s^2")),
            (14., "Si", Configuration("1s^1 2p^6 3s^2 3p^2"),      Configuration("1s^2 2p^5 3s^2 3p^2")),
            (18., "Ar", Configuration("1s^1 2s^2 2p^6 3s^2 3p^6"), Configuration("1s^2 2s^2 2p^5 3s^2 3p^6")),
            (26., "Fe", Configuration("1s^1 2s^2 2p^6 3s^2 3p^6 3d^6 4s^2"),
                        Configuration("1s^2 2s^2 2p^5 3s^2 3p^6 3d^6 4s^2")) ]
        setDefaults("nuclear: charge", Z)
        Arad = Empirical.photoemissionEinsteinA(iC, fC, Empirical.ScaledHydrogenic()).rate
        Aaug = Empirical.augerRate(iC, fC, Empirical.ScaledHydrogenic())
        wK   = Empirical.fluorescenceYield(round(Int64, Z))
        Ax   = round(Defaults.convertUnits("rate: from atomic", Arad), sigdigits=3)
        Bx   = round(Defaults.convertUnits("rate: from atomic", Aaug), sigdigits=3)
        println("    $Z    $el       $Ax        $Bx        $(round(Aaug/Arad, digits=3))        $(round((1-wK)/wK, digits=3))")
    end
    #
    setDefaults("print summary: close", "")
    #
end
