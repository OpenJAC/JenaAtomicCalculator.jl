#
println("Ba) Apply & test the dipole, em multipole and momentum-transfer amplitudes.")

if  true
    # Last successful:  unknown ...
    # Compute 
    grid = Radial.Grid(true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")] )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    flow = 6;    fup = 8;   ilow = 1;   iup = 3
    println("\n\nDipole amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                ##x println("Levels with indices f = $(finalLevel.index) and i = $(initialLevel.index)")
                MultipoleMoment.dipoleAmplitude(finalLevel, initialLevel, grid; display=true)
            end
        end
    end


    println("\n\nElectromagnetic multipole-transition (absorption) amplitudes from PhotoEmission.amplitude(..):\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                PhotoEmission.amplitude("absorption", E1, Basics.Coulomb,  1.0, finalLevel, initialLevel, grid; display=true)
                PhotoEmission.amplitude("absorption", M1, Basics.Magnetic, 1.0, finalLevel, initialLevel, grid; display=true)
            end
        end
    end


    println("\n\nElectromagnetic multipole-transition (absorption) amplitudes from MultipoleMoment.transitionAmplitude(..):\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                MultipoleMoment.transitionAmplitude(E1, Basics.Velocity, 1.0, finalLevel, initialLevel, grid; display=true)
                MultipoleMoment.transitionAmplitude(M1, Basics.Magnetic, 1.0, finalLevel, initialLevel, grid; display=true)
            end
        end
    end


    println("\n\nElectromagnetic multipole-moment amplitudes from MultipoleMoment.amplitude(..):\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                MultipoleMoment.amplitude(E1, Basics.Velocity, 1.0, finalLevel, initialLevel, grid; display=true)
                MultipoleMoment.amplitude(M1, Basics.Magnetic, 1.0, finalLevel, initialLevel, grid; display=true)
            end
        end
    end


    println("\n\nMomentum transfer amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                FormFactor.amplitude(1.0, finalLevel, initialLevel, grid; display=true)
            end
        end
    end
    #
end

