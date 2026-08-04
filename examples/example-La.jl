#
println("Ja) Apply & test the Liouville (time-evolution) computations for the stimulated Raman scattering.")

if  true
    # Last successful:  unknown ...
    nm          = Nuclear.Model(10.0)
    grid        = Radial.Grid(true)
    refConfigs  = [Configuration("[Ne]"), Configuration("1s 2s^2 2p^6 3p"), Configuration("1s^2 2s^2 2p^5 3p")]
    pulse1      = Pulse.FelPulse("GaussianSimplified", 890.0, 3.5e16, 8.5, 50.,  8.0)
    pulse2      = Pulse.FelPulse("GaussianSimplified", 880.0, 3.0e16, 8.5, 50., 10.0)
    pulses      = JenaAtomicCalculator.Pulse.AbstractPulse[pulse1, pulse2]
    notations   = ["     1S_0 ground       == |1>", "1P1 exited        == |2>", "J=1 stimulated    == |3>", "J=1 spontaneous   == |4>", 
                   "               autoionizing/loss == |5>", ]
    settings    = Liouville.Settings(true)
    scheme      = Liouville.StimulatedRamanScheme(LevelSelection(true, indices=[1,2,3, 3]), notations, 0., 0., true)
    
    wa          = Liouville.Computation(Liouville.Computation(), scheme=scheme, pulses=pulses, freeTime=30.,
                                        approach=Liouville.FirstOrderTimeEvolution(), 
                                        nuclearModel=nm, grid=grid, refConfigs=refConfigs,
                                        settings=settings)
    @show wa
    wb          = perform(wa, output=true)
end
