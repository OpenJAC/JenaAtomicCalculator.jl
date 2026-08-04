
println("Dh) Apply & test the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-MultiPhotonDeExcitation.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: energy", "A")   
setDefaults("unit: cross section", "barn")   
setDefaults("unit: rate", "1/s")   

## grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
name = "two-photon absorption or emission of xenon"
    

if  true
    # Last successful:  unknown ...
    # Two-photon absorption of Xe: [Xe] --> [Xe] (5p^-1) 6p: Generation of MeanFieldMultiplet
    name        = "[Xe] (5p^-1) 6s"
    refConfigs  = [Configuration("[Kr] 4d^10 5s^2 5p^5 6s"), Configuration("[Kr] 4d^10 5s^2 5p^5 5d")]
    mfSettings  = MeanFieldSettings()
    #
    wa          = Representation(name, Nuclear.Model(54.01), Radial.Grid(true), refConfigs, MeanFieldMultiplet(mfSettings) )
    wb          = generate(wa, output=true)
    #
elseif  true
    # Last successful:  unknown ...
    # Two-photon absorption of Xe: [Xe] --> [Xe] (5p^-1) 6p: Calculation of absorption cross sections
    nMultiplet  = wb["mean-field multiplet"]
    tpaSettings = MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic(), 
                                                   [E1], [UseCoulomb,JAC.UseBabushkin], nMultiplet, 4,  true, 0.015,
                                                   LineSelection(true, [(1,2), (1,5), (1,6)], Tuple{LevelSymmetry,LevelSymmetry}[]) )
    wc = Atomic.Computation(Atomic.Computation(), name="two-photon absorption", grid=grid, nuclearModel=Nuclear.Model(54.01), 
                            initialConfigs = [Configuration("[Kr] 4d^10 5s^2 5p^6")],
                            finalConfigs   = [Configuration("[Kr] 4d^10 5s^2 5p^5 6p")],  
                            processSettings= tpaSettings )
    ## @show wc
    wd = perform(wc)
    #
elseif  false
    # Last successful:  unknown ...
    # Green function for two-photon absorption or emission of hydrogen-like xenon; this is no longer used (August 2023)
    ## refConfigs       = [Configuration("1s^2 2p")]
    refConfigs       = [Configuration("2p")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(3//2, Basics.minus)]
    greenSettings    = GreenSettings(6, [0, 1], 0.01, true, LevelSelection())

    greenRep         = Representation(name, Nuclear.Model(4.0), grid, refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) )  
                                   ## GreenExpansion( AtomicState.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                   ## GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron()
    greenOut         = generate(greenRep, output=true)
    mpGreen          = greenOut["Green channels"]
    #
elseif true
    # Last successful:  unknown ...
    #
    println("Calculate two-photon absorption cross sections by monochromatic and equally-polarized photons from the same beam")
    mpSettings       = MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic(), 
                                                        [E1], [UseCoulomb, UseBabushkin], mpGreen, 0, true, LineSelection())
    wa = Atomic.Computation(Atomic.Computation(), name=name, grid=grid, nuclearModel=Nuclear.Model(4.0), 
                            initialConfigs= [Configuration("1s")],
                            finalConfigs  = [Configuration("2s")], 
                            processSettings= mpSettings  )
    wb = perform(wa)
    #
elseif  true
    # Last successful:  unknown ...
    #
    println("Calculate two-photon emission cross sections: energy-differential + total cross sections")
    mpSettings       = MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonEmission(), 
                                                        [E1], [UseCoulomb, UseBabushkin], mpGreen, 6, true, LineSelection())
    
    wc = Atomic.Computation(Atomic.Computation(), name=name, grid=grid, nuclearModel=Nuclear.Model(53.), 
                            initialConfigs= [Configuration("1s^2 3s")],
                            finalConfigs  = [Configuration("1s^2 2s")], 
                            processSettings= mpSettings    )
    wd = perform(wc)
    #
end
#
setDefaults("print summary: close", "")


