
#
println("Bb) Apply & test the parity non-conservation, Schiff moment and anapole moment amplitudes.")

if  true
    # Last successful:  unknown ...
    # Compute 
    grid = JAC.Radial.Grid(true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")] )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    nModel = Nuclear.Model(26.);    flow = 6;    fup = 8;   ilow = 1;   iup = 3

    println("\n\nWeak-charge and Schiff-moment amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                ##x println("Levels with indices f = $(finalLevel.index) and i = $(initialLevel.index)")
                ParityNonConservation.weakChargeAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
                ParityNonConservation.schiffMomentAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
            end
        end
    end
    #
end



