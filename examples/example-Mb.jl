
println("Mb) Test  computeTransitionRates(),  computeBranchingFractions()  and  displaySpectrum()  for ForPhotoEmission and ForAutoIonization.")

setDefaults("unit: energy", "eV")
setDefaults("unit: rate",   "1/s")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # computeTransitionRates(ForPhotoEmission): K-alpha radiative rates for the 1s-hole of neutral neon.
    # Initial: 1s vacancy;  final: 2s or 2p vacancy (same number of electrons).
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
    computeTransitionRates(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # computeTransitionRates(ForPhotoEmission): K-alpha radiative rates for the 1s-hole of neutral argon.
    setDefaults("nuclear: charge", 18.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6 3s^2 3p^6")]
    finalConfigs   = [Configuration("1s^2 2s 2p^6 3s^2 3p^6"), Configuration("1s^2 2s^2 2p^5 3s^2 3p^6")]
    computeTransitionRates(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # computeTransitionRates(ForPhotoEmission): radiative rates among the low-lying levels of neutral lithium.
    setDefaults("nuclear: charge", 3.0)
    initialConfigs = [Configuration("1s^2 2p"), Configuration("1s^2 3s"), Configuration("1s^2 3p")]
    finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]
    computeTransitionRates(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # computeTransitionRates(ForAutoIonization): KLL Auger rates for the 1s-hole of neutral neon.
    # Initial: 1s vacancy (9 electrons);  final: Ne^2+ with two L-shell vacancies (8 electrons).
    setDefaults("method: continuum, Galerkin")
    setDefaults("method: normalization, pure sine")
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"), Configuration("[He] 2s^2 2p^4")]
    computeTransitionRates(Basics.ForAutoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # computeTransitionRates(ForAutoIonization): KLL Auger rates for the 1s-hole of neutral argon.
    # Initial: 1s vacancy (17 electrons);  final: Ar^2+ with two L-shell vacancies (16 electrons).
    setDefaults("nuclear: charge", 18.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6 3s^2 3p^6")]
    finalConfigs   = [Configuration("1s^2 2p^6 3s^2 3p^6"), Configuration("1s^2 2s 2p^5 3s^2 3p^6"),
                      Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")]
    computeTransitionRates(Basics.ForAutoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  27Jun2026
    # computeBranchingFractions(ForPhotoEmission): K-alpha branching fractions for the 1s-hole of neutral neon.
    # BF shows how the K-alpha_1 (2p_3/2 -> 1s) and K-alpha_2 (2p_1/2 -> 1s) channels share the radiative decay.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
    computeBranchingFractions(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  27Jun2026
    # computeBranchingFractions(ForAutoIonization): KLL Auger branching fractions for the 1s-hole of neutral neon.
    # BF shows how the 10 KLL channels share the total Auger decay; L23L23 1D2 dominates (~61%).
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"), Configuration("[He] 2s^2 2p^4")]
    computeBranchingFractions(Basics.ForAutoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  28Jun2026
    # displaySpectrum(ForPhotoEmission): photon emission spectrum for low-lying levels of neutral lithium.
    # Three initial configs (2p, 3s, 3p) -> two final configs (2s, 2p): several E1 lines visible.
    using Plots
    setDefaults("nuclear: charge", 3.0)
    initialConfigs = [Configuration("1s^2 2p"), Configuration("1s^2 3s"), Configuration("1s^2 3p")]
    finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]
    displaySpectrum(Basics.ForPhotoEmission(), initialConfigs, finalConfigs,
                    plotfile="spectrum-li-photon.pdf")
    #
elseif  true
    # Last successful:  28Jun2026
    # displaySpectrum(ForAutoIonization): Ne KLL Auger electron spectrum.
    # Three final configurations; L23L23 1D2 line dominates at ~810 eV.
    using Plots
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"), Configuration("[He] 2s^2 2p^4")]
    displaySpectrum(Basics.ForAutoIonization(), initialConfigs, finalConfigs,
                    plotfile="spectrum-ne-auger.pdf")
    #
end


setDefaults("print summary: close", "")
