
#
println("Aj) Apply & test for a Green (-function) expansion.")

if  false
    # Last successful:  unknown ...
    # Compute 
    name          = "Lithium 1s^2 2s ground configuration"
    refConfigs    = [Configuration("[He] 2s")]
    greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, LevelSelection())
    #
    wa          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                                ## GreenExpansion( Atomic.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                ## GreenExpansion( Atomic.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                    GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                                    [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)], 3, greenSettings) )
    println(wa)

    wb = generate(wa, output=true)
    #
elseif true
    # Last successful:  unknown ...
    # Test example for MS "Approximate atomic Green functions" in molecules (2021)
    name            = "Approximate Green function for neon-like ions."
    refConfigs      = [Configuration("1s^2 2s^2 2p^6")]
    levelSymmetries = [LevelSymmetry(1, Basics.minus)]  ## , LevelSymmetry(2, Basics.minus)]
    greenSettings   = GreenSettings( 50, [0, 1, 2], 0.01, true, LevelSelection())
    greenExpansion  = GreenExpansion(AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
    # greenExpansion  = GreenExpansion(AtomicState.DampedSpaceCI(), Basics.DeExciteTwoElectrons(), 
    # greenExpansion  = GreenExpansion(AtomicState.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                     levelSymmetries, 10, greenSettings)
    #
    rep = Representation(name, Nuclear.Model(18.), Radial.Grid(true), refConfigs, greenExpansion)
    #
    gExp = generate(rep, output=true)
    gChannels = gExp["Green channels"]
    cha1      = gChannels[1]
    @show cha1.gMultiplet.levels[1].J, cha1.gMultiplet.levels[1].energy
    @show cha1.gMultiplet.levels[2].J, cha1.gMultiplet.levels[2].energy
    @show cha1.gMultiplet.levels[3].J, cha1.gMultiplet.levels[3].energy
    Basics.displayLevels(stdout, gChannels)
    #
end



