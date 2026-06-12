
## Functions in this file cover: atomic properties.
## Alphabetical order within this file.


"""
`TestFrames.testModule_AlphaVariation(; short::Bool=true)`  ... tests on module AlphaVariation.
"""
function testModule_AlphaVariation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-AlphaVariation-new.sum")
    printstyled("\n\nTest the module  AlphaVariation  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ AlphaVariation.Settings(true, true, LevelSelection() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-AlphaVariation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-AlphaVariation-new.sum"), "Alpha variation parameters:", 1)
    testPrint("testModule_AlphaVariation()::", success)
    return(success)
end


"""
`TestFrames.testModule_DecayYield(; short::Bool=true)`  ... tests on module DecayYield.
"""
function testModule_DecayYield(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DecayYield-new.sum")
    printstyled("\n\nTest the module  DecayYield  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox=10.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(12.),
                            configs=[Configuration("1s 2s^2 2p^6")],
                            propertySettings = [ DecayYield.Settings("SCA", true, false, LevelSelection() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-DecayYield-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-DecayYield-new.sum"), "Fluorescence and Auger", 4)
    testPrint("testModule_DecayYield()::", success)
    return(success)
end


"""
`TestFrames.testModule_Einstein(; short::Bool=true)`  ... tests on module Einstein.
"""
function testModule_Einstein(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Einstein-new.sum")
    printstyled("\n\nTest the module  Einstein  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(36.),
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            propertySettings = [ Einstein.Settings([E1, M1, E2, M2], true,
                                                    LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0)]), 0., 0., 10000. )] )

    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Einstein-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Einstein-new.sum"), "Einstein coefficients, t", 80)
    testPrint("testModule_Einstein()::", success)
    return(success)
end


"""
`TestFrames.testModule_FormFactor(; short::Bool=true)`  ... tests on module FormFactor.
"""
function testModule_FormFactor(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-FormFactor-new.sum")
    printstyled("\n\nTest the module  FormFactor  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ FormFactor.Settings([0.1], true, LevelSelection() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-FormFactor-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-FormFactor-new.sum"), "Standard and modifi", 6)
    testPrint("testModule_FormFactor()::", success)
    return(success)
end


"""
`TestFrames.testModule_Hfs(; short::Bool=true)`  ... tests on module Hfs.
"""
function testModule_Hfs(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Hfs-b-new.sum")
    printstyled("\n\nTest the module  Hfs  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0, 0.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ Hfs.Settings(true, true, true, true, true, false, LevelSelection() )] )

    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    println("aaa  ")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Hfs-b-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Hfs-b-new.sum"), "Level  J Parity          Hartrees", 20)
    println("bbb  success = $success")
    testPrint("testModule_Hfs()::", success)
    return(success)
end


"""
`TestFrames.testModule_IsotopeShift(; short::Bool=true)`  ... tests on module IsotopeShift.
"""
function testModule_IsotopeShift(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-IsotopeShift-new.sum")
    printstyled("\n\nTest the module  IsotopeShift  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ IsotopeShift.Settings(true, true, true, false, true, 0.0, LevelSelection())] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-IsotopeShift-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-IsotopeShift-new.sum"), "IsotopeShift parameters and amplitudes:", 15)
    testPrint("testModule_IsotopeShift()::", success)
    return(success)
end


"""
`TestFrames.testModule_LandeZeeman(; short::Bool=true)`  ... tests on module LandeZeeman.
"""
function testModule_LandeZeeman(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-LandeZeeman-new.sum")
    printstyled("\n\nTest the module  LandeZeeman  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26., "Fermi", 58., 3.75, AngularJ64(5//2), 1.0, 2.0, 0.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ LandeZeeman.Settings(true, true, true, false, true, true, 0.,
                                                                      LevelSelection(), Multiplet() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-LandeZeeman-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-LandeZeeman-new.sum"), "Lande g_J factors and Zeeman amplitudes:", 30)
    testPrint("testModule_LandeZeeman()::", success)
    return(success)
end


"""
`TestFrames.testModule_MultipolePolarizibility(; short::Bool=true)`  ... tests on module MultipolePolarizibility.
"""
function testModule_MultipolePolarizibility(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-MultipolePolarizibility-new.sum")
    printstyled("\n\nTest the module  MultipolePolarizibility  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ MultipolePolarizibility.Settings(EmMultipole[], 0, 0, Float64[], false, LevelSelection() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultipolePolarizibility-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-MultipolePolarizibility-new.sum"),
                            "Multipole polarizibilities and amplitudes:", 5)
    testPrint("testModule_MultipolePolarizibility()::", success)
    return(success)
end
