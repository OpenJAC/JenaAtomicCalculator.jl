

println("Ec) Test of the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")

setDefaults("print summary: open", "zzz-PhotoExcitationAutoIon.sum")


if  true
    # Last successful:  unknown ...
    # 1s photoexcitation and subsequent autoionization of Li-like carbon
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())  # not used
    excSettings   = PhotoExcitationAutoion.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, true, true, ExpStokes(), 
                                                    [SolidAngle(1.0, 0.0)], 0., 4, PathwaySelection())
    grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")  

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), 
                            initialConfigs      = [Configuration("1s^2 2s")],
                            intermediateConfigs = [Configuration("1s 2s 2p")],
                            finalConfigs        = [Configuration("1s^2")], 
                            processSettings = excSettings  )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
