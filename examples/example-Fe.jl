
println("Fe) DR rate coefficients for neon-like Ar: AverageSCA model.")

using JLD2
#
setDefaults("unit: energy", "Hartree")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, Galerkin") setDefaults("method: continuum, asymptotic Coulomb") 
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
grid   = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
temps  = [1.0e+3,   1.0e+4,   1.0e+5,   1.0e+6,   1.0e+7,  1.0e+8]
# temps_au    = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.0) .* temps;     @show temps_au



if  false
    # Last successful:  unknown ...
    # Calculation of the Ne^+  2s^2 2p^6; 2s, 2p --> 3s, 3p, 3d, ... 5d electron capture and subsequent autoionization and stabilization
    setDefaults("print summary: open", "Ne-like-DR-rate.sum")
    name = "DR rate coefficient for neon-like argon and for temperatures < 500 K"
    asfSettings = AsfSettings(AsfSettings(), maxIterationsScf = 36)
    decayShells = Basics.generateShellList(3,3,2)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, asfSettings=asfSettings,
                               approach=Cascade.AverageSCA(),
                               scheme=Cascade.DielectronicRecombinationScheme( [E1], false, Shell("5d"), 5.0, 0.0, 0.0, 1, 
                                                                               [Shell("2s"), Shell("2p")], [Shell("3s"), Shell("3p")],
                                                                               [Shell("5s"), Shell("5p")], decayShells),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  unknown ...
    # Calculation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoionization.sum")

    name = "Photoionization of Si- "
    scheme=Cascade.PhotoIonizationScheme([E1], [0.5], [4.0], [Shell("2s"), Shell("2p")], 
                                         LevelSelection(), [0,1], 0., 0.)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=scheme,
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
else
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
end
