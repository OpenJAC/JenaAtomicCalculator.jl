
println("Fh) Photoabsorption cascade computations.")

using JLD
#
setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
setDefaults("unit: energy", "eV");                     setDefaults("unit: cross section", "barn")
grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)


if  false
    # Last successful:  unknown ...
    # Calculation of a photoabsorption cascade that is purely based on the photoionization of the target;
    # this includes the direct and excitation+autoionization scheme
    setDefaults("print summary: open", "photoabsorption-cascade.sum")
    name    = "Photoabsorption of neon-like Mg"
    pScheme = Cascade.PhotoAbsorptionScheme([E1, M1], [en for en = 50.0:50.0:300.0], Float64[], 
                                            [Shell("2s"), Shell("2p")], Shell[Shell("3s"), Shell("3p")], 
                                            LevelSelection(false),  [0, 1], true, true, 0., 1000. )
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid, scheme=pScheme,
                               approach=Cascade.AverageSCA(),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last successful:  unknown ...
    # Simulation of the (total) photoabsorption cross sections
    setDefaults("print summary: open", "photoabsorption-simulation.sum")

    data     = [JLD.load("zzz-cascade-photoabsorption-computations-2023-08-18T20.jld")]
    name     = "Simulation of total photoabsorption CS for neon-like Mg"
    energies = [en for en = 30.0:1.0:350.0]
    property = Cascade.PhotoAbsorptionSpectrum(false, true, 0.1, 0.01, energies,  Shell[], [(1,1.0)], Configuration[])
    settings = Cascade.SimulationSettings(false, false, 0.)
    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, 
                              property=property, settings=settings, computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
