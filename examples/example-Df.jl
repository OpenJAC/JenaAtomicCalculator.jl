
println("Df) Apply & test the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")

setDefaults("print summary: open", "zzz-Dielectronic.sum")
setDefaults("unit: rate", "1/s")
setDefaults("unit: strength", "cm^2 eV")   
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, Galerkin") setDefaults("method: continuum, asymptotic Coulomb") 
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


if  false
    # Last successful:  13May2024
    # K-LL DR into initially He-like tungsten; comparison with Tu et al. (PP, 2016)
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    drSettings    = Dielectronic.Settings(Dielectronic.Settings(), multipoles = [E1], gauges = [JAC.UseCoulomb, JAC.UseBabushkin],
                                          printBefore = true, electronEnergyShift = -230. )
                                          
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.), 
                            initialConfigs      = [Configuration("1s^2")],
                            intermediateConfigs = [Configuration("1s 2s^2"),  Configuration("1s 2s 2p"),  Configuration("1s 2p^2")],  
                            finalConfigs        = [Configuration("1s^2 2s"),  Configuration("1s^2 2p")],  
                            processSettings     = drSettings )

    wb = perform(wa)
    #
elseif  true
    # Last successful:  13May2024
    # K-LL DR into initially Be-like tungsten; comparison with Tu et al. (PP, 2016)
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    drSettings    = Dielectronic.Settings(Dielectronic.Settings(), multipoles = [E1], gauges = [JAC.UseCoulomb, JAC.UseBabushkin],
                                          printBefore = true, electronEnergyShift = -230. )

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.), 
                            initialConfigs      = [Configuration("1s^2 2s^2")],
                            intermediateConfigs = [Configuration("1s 2s^2 2p^2")],  
                            finalConfigs        = [Configuration("1s^2 2s^2 2p")],  
                            processSettings     = drSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown ...
    # K-LL DR into initially O-like tungsten; comparison with Tu et al. (PP, 2016)
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    drSettings    = Dielectronic.Settings([E1], [JAC.UseCoulomb, JAC.UseBabushkin], true, PathwaySelection(true, indexTriples=[(1,1,0)]), 
                                          -300., 0., 0., CoulombInteraction())
    grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s");   setDefaults("unit: strength", "cm^2 eV")    

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.), 
                            initialConfigs      = [Configuration("1s^2 2s^2 2p^4")],
                            intermediateConfigs = [Configuration("1s 2s^2 2p^6")],  
                            finalConfigs        = [Configuration("1s^2 2s^2 2p^5")],  
                            process = Dierec(), processSettings = drSettings )

    wb = perform(wa)
    #
elseif  true
    # Last successful:  unknown ...
    # K-LL DR into initially He-like carbon; comparison with Xu et al. (PRA, 2016)
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    drSettings    = Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                          PathwaySelection(true, indexTriples=[(1,1,0), (1,5,0), (1,6,0), (1,10,0), (1,11,0)]), 
                                          0., 0., 0., CoulombInteraction())
    grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 12.0)
    setDefaults("unit: rate", "1/s");   setDefaults("unit: strength", "cm^2 eV")    

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), 
                            initialConfigs      = [Configuration("1s^2")],
                            intermediateConfigs = [Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2"), Configuration("1s 2s 3s"),
                                                   Configuration("1s 2s 3p"), Configuration("1s 2s 3d"), Configuration("1s 2p 3s"), Configuration("1s 2p 3p"),
                                                   Configuration("1s 2p 3d")],  
                            finalConfigs        = [Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s^2 3s"), Configuration("1s^2 3p"),
                                                   Configuration("1s^2 3d"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],  
                            initialAsfSettings  = asfSettings, finalAsfSettings=asfSettings,
                            process = Dierec(), processSettings = drSettings)

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown ...
    # DR if initially Li-like Beryllium; comparison with Mohamed et al. (PRA, 2002) ... however difficult and the identification of resonances and
    # rates does not work at present
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    drSettings    = Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, PathwaySelection(), 0., 0., 0., CoulombInteraction())
    grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")   

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(4.), 
                            initialConfigs=[Configuration("1s^2 2s")],
                            intermediateConfigs=[Configuration("1s^2 2s^2"),  Configuration("1s^2 2s 2p"),  Configuration("1s^2 2p^2"),
                                                 Configuration("1s^2 2p 3s"), Configuration("1s^2 2p 3p"),  Configuration("1s^2 2p 3d") ],  
                            finalConfigs       =[Configuration("1s^2 2s^2"),  Configuration("1s^2 2s 2p"),  Configuration("1s^2 2p^2"),
                                                 Configuration("1s^2 2s 3s"), Configuration("1s^2 2s 3p"),  Configuration("1s^2 2s 3d")],  
                            process = Dierec(), processSettings = drSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown ...
    # DR strength for sodium-like silicon: Comparison with Schmidt et al. (PRA, 2007) ... need to be  adapted
    drSettings    = Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                          PathwaySelection(false, indexTriples=[(1,5,0), (1,6,0)]), 0., 0., 0., CoulombInteraction())

    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 1500)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(14.0), 
                        initialConfigs=[Configuration("[Ne] 3s")],
                        intermediateConfigs=[Configuration("[Ne] 3s^2"), Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3s 3d"),
                                             Configuration("[Ne] 3p 4d"), Configuration("[Ne] 3p 4f") ],  
                        finalConfigs  =[Configuration("[Ne] 3s^2"),
                                        Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3s 3d"), Configuration("[Ne] 3s 4d"), Configuration("[Ne] 3s 4f") ],
                        process = Dierec(), 
                        processSettings=Dielectronic.Settings([E1, M1, M2, E2], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              false, Tuple{Int64,Int64,Int64}[(1,10,1), (1,10,7)], 0., 0., 0., CoulombInteraction())  )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown ...
    # Need to be adapted !!!
    refConfigs = [Configuration("[Ne] 3s")]
    fromShells = [Shell("3s")]
    toShells   = [Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f"), Shell("5s")]
    intermediateConfigs = Basics.generateConfigurations(refConfigs, fromShells, toShells)
    grid=JAC.Radial.Grid(true)
    #== grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 2000)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.0, "Fermi"), 
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p^2"), Configuration("1s^2 2s 3s"), 
                                    Configuration("1s^2 2s 3p"), Configuration("1s^2 2s 3d"), Configuration("1s^2 2p 3s"), Configuration("1s^2 2p 3p"), 
                                    Configuration("1s^2 2p 3d")],
                            ## configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                            ##          Configuration("1s^2 7s"), Configuration("1s^2 7p"), Configuration("1s^2 7d"), Configuration("1s^2 7f")], 
                            asfSettings=AsfSettings(AsfSettings(), scField=Basics.DFSField())  ) ==#

    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 5.0e-2, NoPoints = 1500)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(74.0), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), 
                                            Configuration("1s^2 7s"), Configuration("1s^2 7p"), Configuration("1s^2 7d"), Configuration("1s^2 7f")],
                            intermediateConfigs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2p 7s"), Configuration("1s^2 2p 7p"),  
                                                Configuration("1s^2 2p 7d"), Configuration("1s^2 2p 7f")],  
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2s 7s"), Configuration("1s^2 2s 7p"),  
                                            Configuration("1s^2 2s 7d"), Configuration("1s^2 2s 7f")],
                            process = Dierec(), 
                            processSettings=Dielectronic.Settings([E1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                                false, Tuple{Int64,Int64,Int64}[(1,10,1), (1,10,7)], 0., 0., 0., CoulombInteraction())  )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")


