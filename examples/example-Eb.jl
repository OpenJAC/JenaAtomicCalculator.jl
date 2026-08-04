
println("Eb) Test of the PhotoExcitationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.")

setDefaults("print summary: open", "zzz-PhotoExcitationFluores.sum")

if  true
    # Last successful:  unknown ...
    # 1s photoexcitation of Li-like carbon
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    excSettings   = PhotoExcitationFluores.Settings([E1], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, true, true, ExpStokes(), 
                                                    [SolidAngle(1.0, 0.0)], 0., PathwaySelection())
    grid          = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")  

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), 
                            initialConfigs      = [Configuration("1s^2 2s")],
                            intermediateConfigs = [Configuration("1s 2s 2p")],
                            finalConfigs        = [Configuration("1s^2 2s")], 
                            processSettings = excSettings  )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

## [E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], false, false, false, true, true, Tuple{Int64,Int64,Int64}[(1,1,0), (1,2,0), (1,2,0)], 0., 0., 0., CoulombInteraction()
