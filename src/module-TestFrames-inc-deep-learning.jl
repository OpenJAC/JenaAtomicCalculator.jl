
## Functions in this file cover: deep-learning and neural-network based estimates.
## Alphabetical order within this file.


"""
`TestFrames.testModule_DeepLearning(; short::Bool=true)`  ... tests on module DeepLearning.
"""
function testModule_DeepLearning(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DeepLearning-new.sum")
    printstyled("\n\nTest the module  DeepLearning  ... \n", color=:cyan)
    ### Make the tests
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-DeepLearning-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-DeepLearning-approved.sum"),
    ##                             joinpath(@__DIR__, "..", "test", "test-DeepLearning-new.sum"), "xxx", 10)
    testPrint("testModule_DeepLearning()::", success)
    return( success )
end
