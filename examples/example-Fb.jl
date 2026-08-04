

println("Fb) Three-step cascade computation and simulation for the photo-ionization of Si^- and its subsequent decay: AverageSCA model.")

using JLD2
#
println("aa")
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")


if true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-Cascade-computation-photoionization.sum")

    name = "Photoionization of Si- "
    grid = Radial.Grid(false)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                            scheme=Cascade.PhotonIonizationScheme([Photo], 1, [30.0, 80.0]),
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")] )
    println(wa)
    @show name
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false    ## Stepwise decay cascade
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-Cascade-computation-following-decay.sum")

    using JLD2
    JLD2.@load "zzz-cascade-ionizing-computations-2020-02-02T20.jld"
    iniMultiplets = results["generated multiplets:"]

    name = "Si- (1s^-1) decay cascade"
    grid = Radial.Grid(false)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                            scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 1, 0, Shell[], Shell[]),
                            initialMultiplets=iniMultiplets )
    println(wa)
    @show name
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
else
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

    JLD2.@load "zzz-cascade-ionizing-computations-2020-02-02T20.jld"
    resIon  = results
    JLD2.@load "zzz-cascade-decay-computations-2020-02-02T20.jld"
    resDecay = results

    data = [resIon, resDecay]
    name = "Simulation after Si- 1s and 2s ionization and subsequent decay"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution()], ## , Cascade.FinalLevelDistribution()], 
                            settings=Cascade.SimulationSettings(0., 0., 0., 0., 30., [(78, 2.0), (79, 1.0), (80, 0.5)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
