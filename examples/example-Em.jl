

println("Em) Test of the PhotoIonizationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")

if  true
    # Last successful:  unknown ...
    # Compute 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                            process = PhotoIonAuto(), 
                            processSettings=PhotoIonizationAutoion.Settings()  )

    wb = perform(wa)
    #
end

