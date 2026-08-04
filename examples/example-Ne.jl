#
println("Ne) Tests of empirical three-body recombination and the capture channels of an ion.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  16-Jul-2026
    # Branch 1: Empirical.recombinationConfigurations(iConf; nLayers) and the single-channel
    #   Empirical.threeBodyRecombinationPlasmaAlpha(eDist, iConf, fConf) by detailed balance (Saha)
    #   with the Lotz EII coefficient of the inverse process.
    # Systems: bare H^+ and Ne^+ 2p^5; the channel H^+ + 2e -> H(1s) + e in detail.
    # Checks:
    #   - H^+ has 34 capture channels for nLayers = 10 (n = 1 ... 10, l <= 3), Ne^+ 2p^5 has 40
    #     (the 2p vacancy plus the empty n = 3 ... 12 shells).
    #   - Saha identity: alpha^(TBR) = alpha^(EII) g_f/(2 g_i) (2 pi/T_e)^(3/2) exp(P/T_e), exact at T_e = 1 eV.
    #   - Numerical stability: the coefficient stays finite at T_e = 0.05 eV, where P/T_e = 272 and a naive
    #     product alpha^(EII) exp(P/T_e) underflows to 0 * Inf = NaN.
    #
    println("\n  Capture channels and the single-channel TBR coefficient:\n")
    fac6  = Defaults.convertUnits("length: from atomic to cm", 1.0)^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    bare  = Configuration("1s^0");    h1s = Configuration("1s^1")
    setDefaults("nuclear: charge", 1.)
    println("    H^+ (bare)   : ", length(Empirical.recombinationConfigurations(bare, nLayers=10)), " capture channels")
    println("    Ne^+ (2p^5)  : ", length(Empirical.recombinationConfigurations(Configuration("1s^2 2s^2 2p^5"), nLayers=10)),
            " capture channels")
    ##  Saha identity at T_e = 1 eV
    Te    = Defaults.convertUnits("energy: from eV to atomic", 1.0)
    eD    = Distribution.ElectronMaxwell(Te)
    aTBR  = Empirical.threeBodyRecombinationPlasmaAlpha(eD, bare, h1s, printout=true)
    aEII  = Empirical.impactIonizationPlasmaAlpha(eD, h1s, bare)
    P     = Empirical.scaledBindingEnergy(1.0, Shell("1s"), h1s, PeriodicTable.Williams2000())
    saha  = (2/1) / 2 * (2pi/Te)^1.5 * exp(P/Te)
    println("    Saha identity: alpha^(TBR) / [alpha^(EII) * Saha] = ", round(aTBR/(aEII*saha), digits=8), "   [must be 1]")
    ##  Low-temperature stability
    aLow  = Empirical.threeBodyRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 0.05)),
                                                        bare, h1s)
    println("    alpha^(TBR)(1s channel; T_e = 0.05 eV) = ", round(fac6*aLow, sigdigits=4), " cm^6/s   [finite, no NaN]")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 2: the *total* TBR coefficient, summed over all capture channels with the truncation policy
    #   n_max = n_valence + nLayers (default 10).
    # System: bare H^+ for several (low) electron temperatures.
    # Note: this truncated direct-capture sum is *not* the classical collisional-radiative coefficient of
    #   Mansbach & Keck (1969): their scaling alpha ~ T_e^(-9/2) rests on the cascade through levels with
    #   binding energy ~ T_e, which lie above the truncation at low T_e (H: bE ~ T_e = 0.05 eV would need
    #   n ~ 16 > n_max = 10). Capture into such barely bound shells is undone by the next collision, which
    #   is the physical reason for truncating.
    # Checks:
    #   - alpha^(TBR: total)(H^+; 0.1 eV) = 1.537e-23 cm^6/s.
    #   - The measured local slope d(ln alpha)/d(ln T) is only ~ -1.2 ... -1.4 between 0.05 and 0.4 eV,
    #     i.e. much weaker than the classical -4.5, and the magnitude stays ~20x below the classical
    #     estimate -- both consequences of the (physical) truncation.
    #
    println("\n  Total TBR coefficient of bare H^+, summed over the capture channels:\n")
    fac6  = Defaults.convertUnits("length: from atomic to cm", 1.0)^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    setDefaults("nuclear: charge", 1.)
    bare  = Configuration("1s^0")
    Ts    = [0.05, 0.1, 0.2, 0.4];    alphas = Float64[]
    for  TeV  in  Ts
        eD  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", TeV))
        alp = Empirical.threeBodyRecombinationPlasmaAlpha(eD, bare, nLayers=10, printout=(TeV == 0.1))
        push!(alphas, fac6*alp)
        println("    T_e = $TeV eV:   alpha^(TBR: total) = $(round(fac6*alp, sigdigits=4)) cm^6/s")
    end
    for  i = 1:length(Ts)-1
        slope = log(alphas[i+1]/alphas[i]) / log(Ts[i+1]/Ts[i])
        println("    local slope between $(Ts[i]) and $(Ts[i+1]) eV:  ", round(slope, digits=2), "   [Mansbach-Keck cascade: -4.5]")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  16-Jul-2026
    # Branch 3: three-body vs. radiative recombination -- the crossover density.
    # System: H^+ + e -> H(1s), at T_e = 1 eV. Radiative capture (per ion) grows ~ n_e, three-body
    #   capture grows ~ n_e^2; the two balance at the crossover density
    #       n_e^* = alpha^(PR) / alpha^(TBR)   [1/cm^3],
    #   above which the plasma recombines predominantly by three-body capture into this channel.
    # Checks:
    #   - Both coefficients positive; n_e^* of order 1e16 ... 1e17 cm^-3 at T_e = 1 eV, the textbook
    #     ballpark for the transition from a radiation-dominated to a collision-dominated plasma.
    #   - Quantity check: [cm^3/s] / [cm^6/s] = [1/cm^3], a density -- the units tell the physics.
    #
    println("\n  Crossover density between radiative and three-body recombination, H(1s) at T_e = 1 eV:\n")
    fac3  = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    fac6  = Defaults.convertUnits("length: from atomic to cm", 1.0)^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    setDefaults("nuclear: charge", 1.)
    bare  = Configuration("1s^0");    h1s = Configuration("1s^1")
    eD    = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", 1.0))
    aPR   = Empirical.photorecombinationPlasmaAlpha(eD, bare, h1s, printout=true)
    aTBR  = Empirical.threeBodyRecombinationPlasmaAlpha(eD, bare, h1s, printout=true)
    println("    alpha^(PR)  [cm^3/s] = ", round(fac3*aPR,  sigdigits=4))
    println("    alpha^(TBR) [cm^6/s] = ", round(fac6*aTBR, sigdigits=4))
    println("    crossover density n_e^* = alpha^(PR)/alpha^(TBR) = ", round(fac3*aPR/(fac6*aTBR), sigdigits=3), " 1/cm^3")
    #
    setDefaults("print summary: close", "")
    #
end
