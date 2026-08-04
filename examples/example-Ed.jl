
println("Ed) Test of the Photoionization and PlasmaShift module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoIonization-Plasma.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  true
    # Last successful:  unknown ...
    # Compute 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(false), nuclearModel=Nuclear.Model(36.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                            processSettings=PlasmaShift.PhotoSettings(PlasmaShift.DebyeHueckel(), 0.25, 0., 0, [E1], [JAC.UseCoulomb], 
                                                                    [3000., 3200.], true, LineSelection()) )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

