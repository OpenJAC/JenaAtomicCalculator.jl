#
println("Dr) Test of the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-CoulombIonization.sum")


if  true
    # Last successful:  unknown ...
    # Compute 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.), 
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                            # WAS a POSITIONAL CoulombIonization.Settings of the retired shape, led by multipoles and
                            # gauges -- fields the struct no longer has at all, so the first value landed on
                            # ionEnergies::Array{Float64,1} and raised "Cannot convert EmMultipole to Float64". The
                            # struct is now (ionEnergies, electronEnergies, calcAlignment, printBefore, lineSelection,
                            # zerosGL, lValues). [3000., 4000.] is carried over as the ELECTRON energies, which is how
                            # examples/example-Dp.jl uses this same settings object; ionEnergies is left at its default
                            # rather than guessed. `process = Coulion()` is gone too -- the settings object names the
                            # process now.
                            processSettings=CoulombIonization.Settings(CoulombIonization.Settings(),
                                                electronEnergies=[3000., 4000.], calcAlignment=false,
                                                printBefore=false) )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

