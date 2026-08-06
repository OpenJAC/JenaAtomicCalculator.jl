
## Functions in this file cover: atomic processes.
## Alphabetical order within this file.


"""
`TestFrames.testModule_AutoIonization(; short::Bool=true)`  ... tests on module AutoIonization.
"""
function testModule_AutoIonization(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-AutoIonization-new.sum")
    printstyled("\n\nTest the module  AutoIonization  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(36.),
                            initialConfigs=[Configuration("1s^2 2s^2 2p"), Configuration("1s 2s^2 2p^2")],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")],
                            processSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy=true, printBefore=true,
                                                                      lineSelection=LineSelection(true, indexPairs=[(3,1), (4,1), (5,1), (6,1)]),
                                                                      maxKappa=2) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-AutoIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-AutoIonization-new.sum"), "Auger rates and intrinsic", 5)
    testPrint("testModule_AutoIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_CoulombExcitation(; short::Bool=true)`  ... tests on module CoulombExcitation.
"""
function testModule_CoulombExcitation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-CoulombExcitation-new.sum")
    printstyled("\n\nTest the module  CoulombExcitation  ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-CoulombExcitation-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-AutoIonization-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-AutoIonization-new.sum"), "AutoIonization rates and intrinsic angular parameters:", 25)
    testPrint("testModule_CoulombExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_DielectronicRecombination(; short::Bool=true)`  ... tests on module DielectronicRecombination.
"""
function testModule_DielectronicRecombination(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DielectronicRecombination-new.sum")
    printstyled("\n\nTest the module  DielectronicRecombination  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 7.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid,
                            nuclearModel=Nuclear.Model(26.),
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ],
                            processSettings=DielectronicRecombination.Settings(DielectronicRecombination.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin],
                                                                  pathwaySelection=PathwaySelection(true, indexTriples=[(1,1,0)]) )
)
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-DielectronicRecombination-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-DielectronicRecombination-new.sum"),
                            "Total Auger rates", 7)
    testPrint("testModule_DielectronicRecombination()::", success)
    return(success)
end


"""
`TestFrames.testModule_MultiPhotonDeExcitation(; short::Bool=true)`  ... tests on module MultiPhotonDeExcitation.
"""
function testModule_MultiPhotonDeExcitation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-MultiPhotonDeExcitation-new.sum")
    printstyled("\n\nTest the module  MultiPhotonDeExcitation  ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-MultiPhotonDeExcitation-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultiPhotonDeExcitation-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-MultiPhotonDeExcitation-new.sum"), "xxx", 100)
    testPrint("testModule_MultiPhotonDeExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_PhotoEmission(; short::Bool=true)`  ... tests on module PhotoEmission.
"""
function testModule_PhotoEmission(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PhotoEmission-new.sum")
    printstyled("\n\nTest the module  PhotoEmission  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(36.),
                            initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            processSettings=PhotoEmission.Settings([E1, M1, E2, M2], [UseCoulomb, UseBabushkin], true, true, CorePolarization(),
                                LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)]), 0., 0., 10000., false ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoEmission-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoEmission-new.sum"), "Einstein coefficients, t", 100)
    testPrint("testModule_PhotoEmission()::", success)
    return(success)
end


"""
`TestFrames.testModule_PhotoExcitation(; short::Bool=true)`  ... tests on module PhotoExcitation.
"""
function testModule_PhotoExcitation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PhotoExcitation-new.sum")
    printstyled("\n\nTest the module  PhotoExcitation  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(36.),
                            initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, false, true,
                                                                        LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes(0., 0., 0.) ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitation-new.sum"),
                            "Photoexcitation integrated cross sections", 200)
    testPrint("testModule_PhotoExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_PhotoIonization(; short::Bool=true)`  ... tests on module PhotoIonization.
"""
function testModule_PhotoIonization(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PhotoIonization-new.sum")
    printstyled("\n\nTest the module  PhotoIonization  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(36.),
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                            processSettings=PhotoIonization.Settings(PhotoIonization.Settings(), multipoles=[E1, M1], photonEnergies=[3000., 4000.],
                                                                     calcAnisotropy=true, printBefore=true,
                                                                     lineSelection=LineSelection(true, indexPairs=[(1,1), (1,2)])) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"), "Total photoionization c", 3)
    ## Check the summed (grand-total) cross sections separately: they must stay consistent with the
    ## line-resolved table above, of which they are the sum over all final levels.
    success = success  &&
              testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"),
                                "Total photoionization cross sections, summed", 15)
    testPrint("testModule_PhotoIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_PhotoRecombination(; short::Bool=true)`  ... tests on module PhotoRecombination.
"""
function testModule_PhotoRecombination(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PhotoRecombination-new.sum")
    printstyled("\n\nTest the module  PhotoRecombination  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 5.0e-2, hp = 1.0e-2, rbox = 6.5)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.),
                            initialConfigs=[Configuration("1s^2")],
                            finalConfigs  =[Configuration("1s^2 2s"), Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d")],
                            processSettings=PhotoRecombination.Settings([E1, M1], [UseCoulomb, UseBabushkin], [10.],
                                                    [2.18, 21.8, 218.0], false, false, false, false, true, 2, LineSelection() ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoRecombination-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoRecombination-new.sum"), "Photorecombination cross sections", 10)
    testPrint("testModule_PhotoRecombination()::", success)
    return(success)
end


"""
`TestFrames.testModule_RayleighCompton(; short::Bool=true)`  ... tests on module RayleighCompton.
"""
function testModule_RayleighCompton(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-RayleighCompton-new.sum")
    printstyled("\n\nTest the module  RayleighCompton  ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-RayleighCompton-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-RayleighCompton-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-RayleighCompton-new.sum"), "xxx", 100)
    testPrint("testModule_RayleighCompton()::", success)
    return(success)
end


"""
`TestFrames.testModule_HyperfineInduced(; short::Bool=true)`  ... tests on module HyperfineInduced.
"""
function testModule_HyperfineInduced(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-HyperfineInduced-new.sum")
    printstyled("\n\nTest the module  HyperfineInduced  ... \n", color=:cyan)
    ### Make the tests
    ## NUCLEAR HYPERFINE MIXING in H-like 229Th89+ -- deliberately the smallest system that exercises the whole
    ## chain: two nuclear states in one hyperfine basis, the mixing that follows, both terms of the amplitude
    ## (nuclear radiation and the borrowed electronic one), and the level lifetimes. One electron, one electronic
    ## level, four hyperfine levels, five lines -- it runs in seconds.
    ##
    ## Basics.NuclearField() rather than the default DFS: a DFS potential self-interacts badly on a one-electron
    ## system (H 1s comes out at -0.194 instead of -0.5 a.u.) and would corrupt exactly the hyperfine matrix
    ## elements under test.
    elemM = Nuclear.reducedTransitionAmplitude(M1, 0.008, 229, AngularJ64(3//2))
    gsTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.plus, energy = 0.0,
                           mu =  0.360, multipoleM = [M1], elementM = [elemM])
    isTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus, energy = 8.356,
                           mu = -0.378, multipoleM = [M1], elementM = [elemM])
    asfTh = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid  = Radial.Grid(Radial.Grid(false), rnt = 1.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid,
                            nuclearModel   = Nuclear.Model(90.0, "Fermi"),
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("1s")],
                            initialAsfSettings = asfTh, finalAsfSettings = asfTh,
                            processSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                                multipoles = [M1], hfMultipoles = [M1], gauges = [UseCoulomb],
                                isomers = Nuclear.Isomer[gsTh, isTh], calcOverview = false,
                                lineSelection = LineSelection(), printBefore = false, calcLifetimes = true ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-HyperfineInduced-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-HyperfineInduced-new.sum"),
                                "Hyperfine-induced transition rates", 11)
    testPrint("testModule_HyperfineInduced()::", success)
    return(success)
end
