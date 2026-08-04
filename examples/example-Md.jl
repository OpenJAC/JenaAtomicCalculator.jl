
println("Md) Test  computeLifetimes(ForPhotoEmission(), ...)  and  computeLifetimes(ForAutoIonization(), ...)  for several atoms.")

setDefaults("unit: energy", "eV")
setDefaults("unit: rate",   "1/s")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # computeLifetimes(ForPhotoEmission): radiative (K-alpha) lifetime of the Ne 1s-hole state.
    # Final configs are generated automatically from the initial K-hole configuration.
    setDefaults("nuclear: charge", 10.0)
    configs = [Configuration("1s 2s^2 2p^6")]
    computeLifetimes(Basics.ForPhotoEmission(), configs)
    #
elseif  false
    # Last successful:  26Jun2026
    # computeLifetimes(ForAutoIonization): KLL Auger lifetime of the Ne 1s-hole state.
    # Final configs (Ne^2+ two-L-vacancy states) are generated automatically.
    setDefaults("nuclear: charge", 10.0)
    configs = [Configuration("1s 2s^2 2p^6")]
    computeLifetimes(Basics.ForAutoIonization(), configs)
    #
elseif  true
    # Last successful:  26Jun2026
    # Ar 1s-hole: both radiative and KLL Auger lifetimes, to estimate the K-shell fluorescence yield.
    # ω_K = Γ_rad / (Γ_rad + Γ_Auger)  ... read Γ_rad and Γ_Auger [1/s] from the two output tables.
    # Experimental ω_K(Ar) ≈ 12%.
    setDefaults("nuclear: charge", 18.0)
    configs = [Configuration("1s 2s^2 2p^6 3s^2 3p^6")]
    computeLifetimes(Basics.ForPhotoEmission(),    configs)    ## → Γ_rad  [1/s]
    computeLifetimes(Basics.ForAutoIonization(),   configs)    ## → Γ_Auger [1/s]
    #
end


setDefaults("print summary: close", "")
