
#
println("Fg) Expansion opacity calculations.")

using JLD
#
setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)


if  false
    # Last successful:  unknown ...
    # Calculation of the expansion opacity cascade for P^2+
    setDefaults("print summary: open", "opacity-cascade.sum")
    name = "Expansion opacities for P^2+"
    # asfSettings = AsfSettings(AsfSettings(), maxIterationsScf = 36)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(15.), grid=grid, approach=Cascade.SCA(),
                               scheme=Cascade.ExpansionOpacityScheme([E1], 0., 3.0, 0., 1, [Shell("3s"), Shell("3p")], 
                                                                                           [Shell("3s"), Shell("3p"), Shell("3d")], true),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  unknown ...
    # Simulation of the expansion opacity cascade for P^2+
    setDefaults("print summary: open", "opacity-simulation.sum")
    setDefaults("unit: energy", "eV")
    depValues = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3 ];    binning = 0.1
    ## setDefaults("unit: energy", "A")
    ## depValues  = convertUnits("energy: from predefined to atomic unit", [1.0e2, 3.0e2, 1.0e3, 1.0e4, 1.0e5])
    ## binning    = convertUnits("energy: from predefined to atomic unit", 10.0)  ## 1 nm, 10 Kayser


    data = [JLD.load("zzz-cascade-expansion-opacity-computations-2022-08-18T11.jld")]
    name = "Simulation of expansion opacities for P^2+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(), 
                                                                        Cascade.WavelengthOpacityDependence(binning), 1.0e-5, 1000., 1.0e5, 0., depValues), 
                              settings=Cascade.SimulationSettings(true, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last successful:  unknown ...
    # Simulation of the Rosseland opacity cascade for P^2+
    setDefaults("print summary: open", "opacity-simulation.sum")
    setDefaults("unit: energy", "eV")
    binning      = 0.01
    ionDensities = [1.0e-5, 3.0e-5, 5.0e-5];   temperatures = [1.0e3, 1.0e4, 1.0e5]


    data = [JLD.load("zzz-cascade-expansion-opacity-computations-2022-08-16T20.jld")]
    name = "Simulation of Rosseland opacities for P^2+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.RosselandOpacities(Basics.BoltzmannLevelPopulation(),
                                                               Cascade.TemperatureOpacityDependence(binning), ionDensities, temperatures, 1.0, 0.), 
                              settings=Cascade.SimulationSettings(true, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
