
println("Mf) Test  estimateCrossSections(ForImpactIonization(), ...)  for several atoms.")

setDefaults("unit: energy",        "eV")
setDefaults("unit: cross section", "barn")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # BEB electron impact-ionization cross sections for all subshells of neutral neon (Z=10).
    # Energy grid covers all shells: 2p (~21 eV), 2s (~48 eV), 1s (~870 eV).
    setDefaults("nuclear: charge", 10.0)
    initialConfigs   = [Configuration("[Ne]")]
    electronEnergies = [25., 50., 100., 200., 500., 1000., 2000., 5000., 10000.]
    estimateCrossSections(Basics.ForImpactIonization(), initialConfigs, electronEnergies=electronEnergies)
    #
elseif  true
    # Last successful:  26Jun2026
    # BEB electron impact-ionization cross sections for all subshells of neutral argon (Z=18).
    # Energy grid covers all shells: 3p (~16 eV), 3s (~33 eV), 2p/2s (~250-330 eV), 1s (~3200 eV).
    setDefaults("nuclear: charge", 18.0)
    initialConfigs   = [Configuration("[Ar]")]
    electronEnergies = [20., 50., 100., 300., 500., 1000., 3000., 5000., 10000., 50000.]
    estimateCrossSections(Basics.ForImpactIonization(), initialConfigs, electronEnergies=electronEnergies)
    #
end


setDefaults("print summary: close", "")
