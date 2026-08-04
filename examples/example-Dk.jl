
println("Dk) Test of the RadiativeAuger module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-RadiativeAuger.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  true
    # Last successful:  unknown ...
    # Green function for radiative Auger decay of 1s photoionized neon
    name             = "radiative Auger decay of 1s photoionized neon"
    refConfigs       = [Configuration("1s 2s^2 2p^6")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(1//2, Basics.plus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation(name, Nuclear.Model(10.), Radial.Grid(true), refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) ) 
    greenOut         = generate(greenRep, output=true)
    raGreen          = greenOut["Green channels"]
    #
elseif true
    # Last successful:  unknown ...
    # Radiative Auger decay of 1s photoionized neon
    raSettings = RadiativeAuger.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], raGreen, 2, 3, true, CoulombInteraction(), 
                                         LineSelection(true, indexPairs=[(1,1), (1,2)]))

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(10.), 
                            initialConfigs=[Configuration("1s 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6") ], 
                            process = RAuger(), processSettings = raSettings )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

