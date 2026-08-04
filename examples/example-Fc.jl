
println("Fc) Three-step cascade computations and simulations for the decay of the neon 1s^-1 3p hole states: AverageSCA model.")

using JLD
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")


if false 
    # Last successful:  unknown ...
    # Compute 
    #
elseif  false    ## Stepwise decay cascade
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-old-Cascade-computation-following-decay.sum")

    name = "Cascade after neon 1s --> 3p excitation"
    grid = Radial.Grid(false)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                            scheme=Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0, Shell[], Shell[]),
                            initialConfigs=[Configuration("1s^1 2s^2 2p^6 3p")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
else
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

    data = [JLD.load("zzz-cascade-decay-computations-2020-07-13T08.jld")]
    ## data = [wb]
    name = "Simulation of the neon 1s^-1 3p decay"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.IonDistribution(), 
                            settings=Cascade.SimulationSettings(0., 0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
