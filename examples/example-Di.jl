
println("Di) Test of the DoubleAutoIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-DoubleAutoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   
setDefaults("unit: rate", "1/s")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)


if  true
    # Last successful:  unknown ...
    # Green function for double Auger decay of 1s photoionized neon
    name             = "Double Auger decay of 1s photoionized neon"
    refConfigs       = [Configuration("1s 2s^2 2p^6")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)]
    greenSettings    = GreenSettings(7, [0, 1, 2], 0.01, true, LevelSelection())
    greenexpansion   = GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 9, greenSettings)
    greenRep         = Representation(name, Nuclear.Model(10.), grid, refConfigs, greenexpansion ) 
    greenOut         = generate(greenRep, output=true)
    daGreen          = greenOut["Green channels"]
    #
elseif true
    # Last successful:  unknown ...
    # Double Auger decay of 1s photoionized neon
    daSettings = DoubleAutoIonization.Settings(daGreen, 2, true, 0.1, 1000., 2, CoulombInteraction(), LineSelection(true, indexPairs=[(1,4)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  = [Configuration("1s 2s^2 2p^6")],
                            ## finalConfigs    = [Configuration("1s^2 2p^5")], 
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^3"), Configuration("1s^2 2s 2p^4"), Configuration("1s^2 2p^5")], 
                            processSettings = daSettings )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")


