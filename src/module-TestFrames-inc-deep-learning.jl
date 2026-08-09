
## Functions in this file cover: deep-learning and neural-network based estimates.
## Alphabetical order within this file.


"""
`TestFrames.testModule_DeepLearning(; short::Bool=true)`  ... tests on module DeepLearning.
"""
function testModule_DeepLearning(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DeepLearning-new.sum")
    printstyled("\n\nTest the module  DeepLearning  ... \n", color=:cyan)
    ### Make the tests
    ## REWRITTEN 08-Aug-2026; it previously asserted nothing at all. See the note in
    ## testModule_CoulombExcitation for what these scaffold checks cover.
    ##
    ## NOT COVERED: any network, any training, any level estimate. Stage 2 of this work produced an honest
    ## NEGATIVE result (the computed-energy baseline beat both the MLP and the GBT), so there is no accuracy
    ## claim here to protect -- but the plumbing should at least be known to construct.
    success = true
    applic = DeepLearning.Application()
    if  length(sprint(show, applic))   == 0     success = false;   println("** empty show for DeepLearning.Application")   end
    if  length(Base.string(applic))    == 0     success = false;   println("** empty string for DeepLearning.Application") end
    ## the copy-constructor must carry a changed value through
    newApplic = DeepLearning.Application(applic; name="renamed application")
    if  newApplic.name != "renamed application"
                                                success = false;   println("** name not carried by the copy-constructor")  end
    if  typeof(newApplic.request) != typeof(applic.request)
                                                success = false;   println("** request altered by an unrelated keyword")    end
    ###
    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Scaffold checks only for DeepLearning; NO physics is compared against approved data.")
    testPrint("testModule_DeepLearning()::", success)
    return( success )
end
