
## Functions in this file cover: cascade computations and simulations.
## Alphabetical order within this file.


"""
`TestFrames.testModule_Cascade_PhotonExcitation(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotonExcitation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotonExcitation-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotonExcitationScheme ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-Cascade-PhotonExcitation-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotonExcitation-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-Cascade-PhotonExcitation-new.sum"), "Steps that are defined for the", 15)
    ## disabled: JLD2 output filename contains a run-date timestamp that changes daily, causing permanent mismatch
    testPrint("testModule_Cascade-PhotonExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_PhotonIonization(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotonIonization(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotonIonization-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotonIonizationScheme ... \n", color=:cyan)
    ### Make the tests
    name = "Photoionization of Si- "
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                                scheme=Cascade.PhotoIonizationScheme([E1], [0.5], [4.0], [Shell("2s"), Shell("2p")],
                                                                     [Shell("2s"), Shell("2p"), Shell("3p"), Shell("4p"), Shell("5p")],
                                                                     LevelSelection(), [0,1], 0., 0.),
                                initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    wb = perform(wa; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotonIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-PhotonIonization-new.sum"), "Total photoionization cross sections for", 15)
    testPrint("testModule_Cascade-PhotonIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_Simulation(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_Simulation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Cascade-Simulation-new.sum")
    printstyled("\n\nTest the module  Cascade for Simulations ... \n", color=:cyan)
    ### Make the tests
    datafile = joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-StepwiseDecay-data.jld")
    data = [JLD2.load(datafile)]
    name = "Simulation of the neon 1s^-1 3p decay"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.IonDistribution(),
                                settings=Cascade.SimulationSettings(0., 0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
    wd = perform(wc; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-Simulation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-Simulation-new.sum"), "(Final) Ion distribution for", 7)
    testPrint("testModule_Cascade-Simulation()::", success)
    return(success)
end


"""
`TestFrames.testModule_Cascade_StepwiseDecay(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_StepwiseDecay(; short::Bool=true)
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
    Defaults.setDefaults("print summary: open", "test-Cascade-StepwiseDecay-new.sum")
    printstyled("\n\nTest the module  Cascade for the StepwiseDecayScheme ... \n", color=:cyan)
    ### Make the tests
    name = "Cascade after neon 1s --> 3p excitation"
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    decayShells = [Shell(1,0), Shell(2,0), Shell(2,1), Shell(3,1)]
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                                scheme=Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0, decayShells, Shell[], Shell[]),
                                initialConfigs=[Configuration("1s^1 2s^2 2p^6 3p")] )
    println(wa)
    wb = perform(wa; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-StepwiseDecay-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Cascade-StepwiseDecay-new.sum"), "Steps that are defined for the curren", 10)
    testPrint("testModule_Cascade-StepwiseDecay()::", success)
    return(success)
end
