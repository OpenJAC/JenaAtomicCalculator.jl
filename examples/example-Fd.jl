
println("Fd) Photoabsorption of Ne^+: AverageSCA model.")

using JLD


if  false
    # Last successful:  unknown ...
    # Calculation of the Ne^+  2s^2 2p^5; 2s --> 2p, 3p, 4p, 5p excitation cross sections
    setDefaults("print summary: open", "Ne-plus-photoexcitation.sum")
    setDefaults("method: continuum, asymptotic Coulomb")                                                      ## "method: continuum, Galerkin"
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, NoPoints=1110)
    name = "Photoabsorption calculations for Ne^+ for energies = [1, 4] a.u."
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotonExcitationScheme([PhotoExc()], [E1], 0.5, 4.0, 1, [Shell("2s"), Shell("2p")], 
                                                                     [Shell("2s"), Shell("2p"), Shell("3p"), Shell("4p"), Shell("5p")]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Last successful:  unknown ...
    # Calculation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoionization.sum")

    name = "Photoionization of Si- "
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, NoPoints=1110)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotonIonizationScheme([Photo()], 1, [1.0, 3.0, 5.0]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Last successful:  unknown ...
    # Simulation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoabsorption-simulation.sum")

    data = [JLD.load("zzz-cascade-ionizing-computations-2020-07-10T19.jld"), JLD.load("zzz-cascade-excitation-computations-2020-07-10T19.jld")]
    ##x JLD2.@load "zzz-cascade-ionizing-computations-2020-06-23T09.jld"
    ##x resIon  = results
    ##x JLD2.@load "zzz-cascade-excitation-computations-2020-06-23T09.jld"
    ##x resExc = results
    ##x data = [resExc]

    name = "Simulation of photoabsorption for Ne^+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.PhotoAbsorption(), 
                              settings=Cascade.SimulationSettings(0., 0., 1., 4., 3., [(1, 1.0)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Last successful:  unknown ...
    # Example for Atomic cascade computation (symmetry, 2021)
    # Calculation of the Mg  1s 2s^2 2p^6 3s^2 decay cascade after 1s photoionization
    setDefaults("print summary: open", "Mg-1s-decay.sum")
    
    decayScheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 3, Dict{Int64,Float64}(), 0., Shell[], Shell[])
    
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    name = "Computation of Mg  1s 2s^2 2p^6 3s^2 decay cascade after 1s photoionization"
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid, 
                               approach=Cascade.AverageSCA(), scheme=decayScheme,
                               ## initialConfigs=[Configuration("1s 2s^2 2p^6 3s^2")] )
                               initialConfigs=[Configuration("1s 2s^2 2p^6 3s^2 3p")] )
    @show wa
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  unknown ...
    # Simulation of the Mg  1s 2s^2 2p^6 3s^2 decay cascade after 1s photoionization
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")
    
    ## simulationSettings = Cascade.SimulationSettings(0., 0., 0., 0., 0., [(1, 1.0)])
    ## simulationSettings = Cascade.SimulationSettings(0., 0., 0., 0., 0., [(4, 1.0)])
    simulationSettings = Cascade.SimulationSettings(0., 0., 0., 0., 0., [(1, 0.25), (2, 0.25), (3, 0.25), (4, 0.25)])

    data = [JLD.load("zzz-cascade-decay-computations-2021-03-14T19.jld")]
    name = "Simulation of Mg  1s 2s^2 2p^6 3s^2 decay cascade after 1s photoionization"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.IonDistribution(), 
                              settings=simulationSettings, computationData=data )
    @show wc
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
