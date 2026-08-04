
println("Ma) Test  computeLevelEnergies()  for several atoms, configurations and isoelectronic sequences.")

setDefaults("unit: energy", "eV")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # Level energies for neon ground configuration and the first excited 2p^5 3s configuration.
    # Simple demo; also used in demos/A4-Pluto-for-pedestrians.jl.
    setDefaults("nuclear: charge", 10.0)
    configs = [Configuration("[Ne]"), Configuration("[He] 2s^2 2p^5 3s")]
    computeLevelEnergies(Basics.ForGivenConfigs(), configs)
    #
elseif  false
    # Last successful:  26Jun2026
    # Level energies for the argon ground configuration and the K-hole configuration;
    # corresponds to Figure 3(a) of the pedestrian review article.
    setDefaults("nuclear: charge", 18.0)
    configs = [Configuration("[Ar]"), Configuration("1s 2s^2 2p^6 3s^2 3p^6")]
    computeLevelEnergies(Basics.ForGivenConfigs(), configs)
    #
elseif  false
    # Last successful:  26Jun2026
    # Level energies for several low-lying configurations of neutral oxygen (open-shell case).
    setDefaults("nuclear: charge", 8.0)
    configs = [Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")]
    computeLevelEnergies(Basics.ForGivenConfigs(), configs)
    #
elseif  false
    # Last successful:  27Jun2026
    # He-like isoelectronic sequence Z=4..10: configuration-averaged energies of 1s^2 (ground)
    # and 1s 2p (first excited).  Shows Z^2 scaling and the 1s->2p gap across Be^2+ ... Ne^8+.
    using Plots
    Zvalues = collect(4.0:1.0:10.0)
    configs  = [Configuration("1s^2"), Configuration("1s 2p")]
    computeLevelEnergies(Basics.ForIsoelectronicSequence(), Zvalues, configs,
                         plotfile="isoelectronic-he-like.pdf")
    #
elseif  true
    # Last successful:  27Jun2026
    # Li-like isoelectronic sequence Z=6..14: configuration-averaged energies of 1s^2 2s and 1s 2p^2.
    # Shows the 2s-2p energy splitting and its Z-dependence across C^3+ ... Si^11+.
    using Plots
    Zvalues = collect(6.0:2.0:14.0)
    configs  = [Configuration("1s^2 2s"), Configuration("1s 2p^2")]
    computeLevelEnergies(Basics.ForIsoelectronicSequence(), Zvalues, configs,
                         plotfile="isoelectronic-li-like.pdf")
    #
end


setDefaults("print summary: close", "")
