
## Functions in this file cover: plasma shift and plasma environment effects.
## Alphabetical order within this file.


#==  disabled: PlasmaShift is no longer an Atomic property; it is a Plasma.Computation scheme (like Cascade),
##   not yet fully implemented. Postponed until the Plasma module is reworked.
"""
`TestFrames.testModule_PlasmaShift(; short::Bool=true)`  ... tests on module PlasmaShift.
"""
function testModule_PlasmaShift(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PlasmaShift-new.sum")
    printstyled("\n\nTest the module  PlasmaShift  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings = [ Plasma.Settings()] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PlasmaShift-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PlasmaShift-new.sum"), "Plasma screening", 3)
    testPrint("testModule_PlasmaShift()::", success)
    return(success)
end  ==#
