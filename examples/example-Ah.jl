#
println("Ah) Apply & test a mean-field basis, mean-field multiplets as well as a configuration-interaction (CI) expansions.")

if  true
    # Last successful:  unknown ...
    # Compute 
    name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
    refConfigs  = [Configuration("[He] 2s^2 2p^4")]
    mfSettings  = MeanFieldSettings()
    #
    wa1         = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldBasis(mfSettings) )
    println("wa1 = $wa1")

    wb1 = generate(wa1, output=true)
    #
    wa2         = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldMultiplet(mfSettings) )
    println("wa2 = $wa2")

    wb2 = generate(wa2, output=true)

    orbitals    = wb1["mean-field basis"].orbitals
    ciSettings  = CiSettings(CoulombInteraction(), LevelSelection())
    from        = [Shell("2s")]
    #
    frozen      = [Shell("1s")]
    to          = [Shell("2s"), Shell("2p")]
    excitations = RasStep()
    #             RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    wc          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                                CiExpansion(orbitals, excitations, ciSettings) )
    println("wc = $wc")

    wd = generate(wc, output=true)
    #
end

