
println("Ch) Apply & test the RadiativeOpacity module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-DecayYield.sum")
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.), 
                            configs=[Configuration("1s 2s^2 2p^6")],
                            propertySettings=[ DecayYield.Settings("SCA", true, LevelSelection()) ] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


