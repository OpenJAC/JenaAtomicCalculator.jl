
println("Fj) Three-step cascade computation and simulation for the photo-ionization of Si^- and its subsequent decay: AverageSCA model.")

# PROVENANCE (05-Aug-2026): this is the former contents of example-Fb.jl, moved here unchanged so that Fb.jl
# could be rebuilt around the cascade SIMULATIONS. The body below is a byte-exact copy taken from git
# (commit 077d6ca, "Track the example scripts"), not a retyping. It is kept because the PhotoIonizationScheme
# scenario it sketches is to be picked up later; nothing here has been repaired.
#
# KNOWN BREAKAGES, none of them fixed:
#   * Cascade.PhotonIonizationScheme does not exist; the type is Cascade.PhotoIonizationScheme, and it takes
#     nine fields (multipoles, photonEnergies, electronEnergies, excitationFromShells, excitationToShells,
#     initialLevelSelection, lValues, electronEnergyShift, minCrossSection), not the three given below.
#   * Cascade.StepwiseDecayScheme takes seven fields (processes, maxElectronLoss, chargeStateShifts,
#     NoShakeDisplacements, decayShells, shakeFromShells, shakeToShells); branch b passes five, and passes the
#     processes as types (Auger, Radiative) rather than as instances (Auger(), Radiative()).
#   * The .jld files named below are hard-coded artefacts from February 2020 and no longer exist.
#   * Cascade.Simulation now takes a single `property=`, not a `properties=` array, and
#     Cascade.SimulationSettings carries three fields (printTree, printLongTree, initialPhotonEnergy), not the
#     six given in branch c; the initial level occupations moved into the property itself.
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
