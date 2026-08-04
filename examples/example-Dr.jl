#
println("Du) Test of the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-CoulombIonization.sum")


if  true
    # Last successful:  unknown ...
    # Compute 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                            process = Coulion(), 
                            processSettings=CoulombIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 4000.], false, true, 
                                            false, Tuple{Int64,Int64}[]) )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

