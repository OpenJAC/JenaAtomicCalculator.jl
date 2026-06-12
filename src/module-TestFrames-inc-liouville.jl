
## Functions in this file cover: open quantum systems and Liouville-space computations.
## Alphabetical order within this file.


"""
`TestFrames.testModule_Liouville(; short::Bool=true)`  ... tests on module Liouville.
"""
function testModule_Liouville(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Liouville-new.sum")
    printstyled("\n\nTest the module  Liouville  ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-Liouville-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Liouville-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-Liouville-new.sum"), "xxx", 10)
    testPrint("testModule_Liouville()::", success)
    return( success )
end
