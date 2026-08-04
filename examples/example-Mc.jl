
println("Mc) Test  computeCrossSections(ForPhotoIonization(), ...)  and  computeCrossSections(ForPhotoRecombination(), ...)  for several atoms.")

setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "barn")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # Photoionization cross sections for neutral neon; 2s and 2p outer-shell channels.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("[Ne]")]
    finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
    computeCrossSections(Basics.ForPhotoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # Photoionization cross sections for neutral argon; 3s and 3p outer-shell channels.
    # Corresponds to Figure 3(b) of the pedestrian review article.
    setDefaults("nuclear: charge", 18.0)
    initialConfigs = [Configuration("[Ar]")]
    finalConfigs   = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")]
    computeCrossSections(Basics.ForPhotoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # Photoionization cross sections for neutral neon; all subshell channels including K-shell.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("[Ne]")]
    finalConfigs   = [Configuration("1s 2s^2 2p^6"), Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
    computeCrossSections(Basics.ForPhotoIonization(), initialConfigs, finalConfigs)
    #
elseif  false
    # Last successful:  26Jun2026
    # Photorecombination cross sections: electron capture into He-like neon (Z=10).
    # Initial: H-like Ne^9+ (1s^1);  capture into 1s --> He-like Ne^8+ (1s^2).
    # Time-reverse of the Ne 1s photoionization channel.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s")]
    intoShells     = [Shell("1s")]
    computeCrossSections(Basics.ForPhotoRecombination(intoShells), initialConfigs)
    #
elseif  true
    # Last successful:  26Jun2026
    # Photorecombination cross sections: electron capture into Li-like neon (Z=10).
    # Initial: He-like Ne^8+ (1s^2);  capture into 2s or 2p --> Li-like Ne^7+.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s^2")]
    intoShells     = [Shell("2s"), Shell("2p")]
    computeCrossSections(Basics.ForPhotoRecombination(intoShells), initialConfigs)
    #
end


setDefaults("print summary: close", "")
