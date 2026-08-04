#
println("Nj) Tests of empirical dielectronic recombination (DR), Arnaud & Rothenflug (1985).")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  24-Jul-2026
    # Branch 1: H-like sequence DR, Empirical.dielectronicRecombinationPlasmaAlpha(eDist, iConf) via the
    #   4 tabulated ions of Table IIB (Arnaud & Rothenflug 1985): O VIII, Mg XII, Ca XX, Fe XXVI.
    # Checks:
    #   - alpha^(DR) is positive and lies in the physically expected ~1e-13 - 1e-10 cm^3/s range for all 4 ions,
    #     evaluated at their own peak temperature T = 2*T0/3 (the analytic maximum of T^(-3/2) exp(-T0/T)).
    #   - The peak really is a local maximum: alpha^(DR) at T = 2*T0/3 exceeds its value at T shifted +/-20%.
    #   Verified: O VIII 1.316e-12, Mg XII 6.89e-13, Ca XX 4.803e-13, Fe XXVI 3.255e-13 cm^3/s, all local maxima.
    #
    println("\n  H-like sequence DR, Table IIB ions, each at its own peak temperature T = 2 T0/3:\n")
    fac3 = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    ions = [ ("O VIII",  8), ("Mg XII", 12), ("Ca XX", 20), ("Fe XXVI", 26) ]
    for  (name, Z)  in  ions
        setDefaults("nuclear: charge", Float64(Z))
        hLike  = Configuration("1s^1")
        ADI, T0 = PeriodicTable.dielectronicRecombinationParameters_Arnaud1985(Z)
        local Tpeak = 2*T0/3
        local eD = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", Tpeak))
        aPeak  = Empirical.dielectronicRecombinationPlasmaAlpha(eD, hLike)
        eDlo   = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 0.8*Tpeak))
        eDhi   = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.2*Tpeak))
        aLo    = Empirical.dielectronicRecombinationPlasmaAlpha(eDlo, hLike)
        aHi    = Empirical.dielectronicRecombinationPlasmaAlpha(eDhi, hLike)
        isMax  = aPeak > aLo  &&  aPeak > aHi
        println("    $name (Z=$Z):  T_peak = $(round(Tpeak, sigdigits=4)) K   " *
                "alpha^(DR) = $(round(fac3*aPeak, sigdigits=4)) cm^3/s   local max: $isMax")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 2: He-like sequence DR, general Z-scaling formula (Sect. 3.2.2), for Fe XXV -- a classic
    #   astrophysical benchmark ion. Unlike branch 1 (H-like), this sequence's peak is *not* pre-tabulated, so
    #   the peak location chi_peak = (2/3) D(z) is derived analytically here (same calculus as T = 2 T0/3 in
    #   branch 1, applied to chi^(-3/2) exp(-D/chi) instead of T^(-3/2) exp(-T0/T)) and converted to T_peak via
    #   the actual I_Li from PeriodicTable.ionizationPotentials_Nist2025, rather than guessing a temperature or
    #   magnitude ahead of time.
    #   Verified: I_Li(Fe) = 2045.8 eV, T_peak = 5.084e7 K, alpha^(DR) = 3.716e-13 / 3.869e-13 / 3.779e-13 cm^3/s
    #   at 0.8/1.0/1.2 x T_peak -- a genuine local max, and close to branch 1's H-like Fe XXVI value (3.255e-13).
    # Checks:
    #   - alpha^(DR) is positive at T_peak and exceeds its value at T shifted +/-20% (local max).
    #   - Magnitude is broadly consistent with the (independently validated) H-like Fe XXVI result from branch 1
    #     (3.255e-13 cm^3/s) -- both are similarly highly-charged ions, so a similar order of magnitude is the
    #     relevant sanity check, not an absolute number guessed without the code.
    #   - Z = 2 (He itself, z = Z-2 = 0) is correctly rejected (division by z^2 in the A(z), D(z) formulas).
    #
    println("\n  He-like sequence DR for Fe XXV (Z=26), at its analytic peak temperature:\n")
    fac3 = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    setDefaults("nuclear: charge", 26.)
    heLike = Configuration("1s^2")
    z      = 26 - 2
    Dz     = 3.0*(z+1)^2/z^2 / (1 + 0.015*z^3/(z+1)^3)
    chiPeak = 2*Dz/3
    ILi_eV  = PeriodicTable.ionizationPotentials_Nist2025(26)[26-2]
    kTpeak  = chiPeak * Defaults.convertUnits("energy: from eV to atomic", ILi_eV)
    Tpeak   = Defaults.convertUnits("temperature: from atomic to Kelvin", kTpeak)
    println("    I_Li(Fe) = $(round(ILi_eV, digits=1)) eV,   D(z=24) = $(round(Dz, digits=4)),   " *
            "chi_peak = $(round(chiPeak, digits=4)),   T_peak = $(round(Tpeak, sigdigits=4)) K")
    for  factor  in  [0.8, 1.0, 1.2]
        local eD = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", factor*Tpeak))
        aDR = Empirical.dielectronicRecombinationPlasmaAlpha(eD, heLike, printout=(factor == 1.0))
        println("    T_e = $(round(factor*Tpeak, sigdigits=4)) K (x$factor):   alpha^(DR)(Fe XXV) = " *
                "$(round(fac3*aDR, sigdigits=4)) cm^3/s")
    end
    #
    println("\n  Rejection of the unphysical z = Z-2 = 0 case (He itself):")
    setDefaults("nuclear: charge", 2.)
    try
        Empirical.dielectronicRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(1.0), Configuration("1s^2"))
        println("    ERROR: should have thrown!")
    catch err
        println("    Correctly rejected (Z=2, z=0):  ", sprint(showerror, err))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 3: DR vs. RR -- the physical motivation for implementing DR at all. Same capture channel,
    #   Fe XXV (He-like) + e -> Fe XXIV (Li-like, ground state), computed two independent ways:
    #   Empirical.dielectronicRecombinationPlasmaAlpha (Arnaud 1985 fit) vs.
    #   Empirical.photorecombinationPlasmaAlpha (Milne relation from the photoionization cross section).
    # Checks:
    #   - Both coefficients positive.
    #   - DR/RR grows strongly with T and crosses 1 within the tested range -- the textbook reason DR cannot be
    #     neglected in astrophysical ionization-balance modeling once the plasma is hot enough.
    #   Verified: T_e=3e6K: DR/RR=4.13e-9 (RR dominates); T_e=1e7K: DR/RR=0.0867; T_e=3e7K: DR/RR=7.08
    #   (DR now dominates by 7x) -- a clean, physically meaningful crossover between 1e7 and 3e7 K.
    #
    println("\n  DR vs. RR for the same channel, Fe XXV + e -> Fe XXIV(2s):\n")
    fac3 = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    setDefaults("nuclear: charge", 26.)
    heLike = Configuration("1s^2");   liLike = Configuration("1s^2 2s^1")
    for  TK  in  [3.0e6, 1.0e7, 3.0e7]
        local eD  = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", TK))
        aDR = Empirical.dielectronicRecombinationPlasmaAlpha(eD, heLike)
        aRR = Empirical.photorecombinationPlasmaAlpha(eD, heLike, liLike)
        println("    T_e = $TK K:   alpha^(DR) = $(round(fac3*aDR, sigdigits=4)) cm^3/s   " *
                "alpha^(RR) = $(round(fac3*aRR, sigdigits=4)) cm^3/s   " *
                "DR/RR = $(round(aDR/aRR, sigdigits=3))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  true
    #
    # Last successful:  24-Jul-2026
    # Branch 4: scope check -- an isoelectronic sequence beyond H-like/He-like must be rejected cleanly, not
    #   silently mishandled, since Arnaud & Rothenflug (1985) does not tabulate self-contained DR parameters
    #   for Li-like and beyond (see project notes).
    # System: Li-like carbon, C^3+ (1s^2 2s^1), Z = 6.
    # Check: an informative error is raised, not a wrong number.
    #   Verified: raises "...supported only for the H-like...and He-like...sequences; got iConf = ... with
    #   NoElectrons = 3", correctly identifying and naming the unsupported case.
    #
    println("\n  Scope check: Li-like sequence (unsupported) must raise an informative error, not a wrong number:\n")
    setDefaults("nuclear: charge", 6.)
    liLikeC = Configuration("1s^2 2s^1")
    eD      = Distribution.ElectronMaxwell(Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.0e6))
    try
        Empirical.dielectronicRecombinationPlasmaAlpha(eD, liLikeC)
        println("    ERROR: should have thrown!")
    catch err
        println("    Correctly rejected (Li-like, NoElectrons=3):  ", sprint(showerror, err))
    end
    #
    setDefaults("print summary: close", "")
    #
end
