#
println("Ba) Apply & test the dipole, em multipole and momentum-transfer amplitudes.")

if  true
    # Last successful:  28-Aug-2026 for the PhotoEmission section, which is what this file is used to check.
    #   REPAIRED TODAY: `PhotoEmission.amplitude` dispatches on a TYPE, not a string, and this file still passed
    #   "absorption". It had therefore been dead with `MethodError: no method matching amplitude(::String, ...)`.
    #   Now `Basics.Absorption()`; note the qualification is needed, the type is not exported into Main.
    #   WHAT IT SHOWS, and it is the reason this file is a good PhotoEmission check: the E1 and M1 amplitudes obey
    #   the parity selection rule EXACTLY, pair by pair. E1 is non-zero only between OPPOSITE parities and M1 only
    #   between the SAME parity, with the other identically 0.0 -- e.g. level 6 (3/2-) <- 1 (1/2+) gives
    #   E1 = -4.2323e-4 with M1 = 0, while 6 (3/2-) <- 2 (1/2-) gives M1 = -1.1860e-6 with E1 = 0.
    #   THE FILE STILL STOPS LATER, and deliberately: the MultipoleMoment.transitionAmplitude section below is
    #   PARKED under Rule 13 (fce090d, priority item 61) and raises with an explanation rather than returning an
    #   unvalidated number. That section will run again when item 61 is settled; nothing here needs fixing for it.
    # Previously: unknown ...
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
                PhotoEmission.amplitude(Basics.Absorption(), E1, Basics.Coulomb,  1.0, finalLevel, initialLevel, grid; display=true)
                PhotoEmission.amplitude(Basics.Absorption(), M1, Basics.Magnetic, 1.0, finalLevel, initialLevel, grid; display=true)
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

