
println("Dg) Apply & test the RayleighCompton module with ASF from an internally generated initial- and final-state multiplets.")

setDefaults("print summary: open", "zzz-RayleighCompton.sum")
setDefaults("unit: rate", "1/s")  


if  false
    # Last successful:  13May2024
    # Green function for Rayleigh scattering on the ground-configuration levels of B-like neon
    grid             = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    name             = "Rayleigh scattering on Be-like 1s^2 2s^2 2p neon"
    refConfigs       = [Configuration("[He] 2s^2 2p")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.minus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation(name, Nuclear.Model(10.), Radial.Grid(true), refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) ) 
    greenOut         = generate(greenRep, output=true)
    green            = greenOut["Green channels"]
    #
elseif true
    # Last successful:  13May2024
    # Green function for Rayleigh scattering on the ground-configuration levels of B-like neon
    grid             = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    asfSettings      = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    rayleighSettings = RayleighCompton.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 100.], green, true, true, true, 
                                                ExpStokes(), [SolidAngle(1.0, 0.0)], LineSelection())
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  = [Configuration("1s^2 2s^2 2p")],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p")], 
                            processSettings = rayleighSettings )
    @show wa
    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
