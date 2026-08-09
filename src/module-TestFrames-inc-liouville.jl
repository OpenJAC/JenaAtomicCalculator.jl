
## Functions in this file cover: open quantum systems and Liouville-space computations.
## Alphabetical order within this file.


"""
`TestFrames.testModule_Liouville(; short::Bool=true)`  ... tests on module Liouville.
"""
function testModule_Liouville(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Liouville-new.sum")
    printstyled("\n\nTest the module  Liouville  ... \n", color=:cyan)
    ### Make the tests
    ## REWRITTEN 08-Aug-2026; it previously asserted nothing at all -- an empty body and an unconditional
    ## success = true, so the suite counted a pass for a module it never touched. See the note in
    ## testModule_CoulombExcitation for what these scaffold checks cover.
    ##
    ## NOT COVERED: any density-matrix evolution, Raman amplitude or level population.
    success = true
    settings = Liouville.Settings()
    if  length(sprint(show, settings)) == 0     success = false;   println("** empty show for Liouville.Settings")   end
    ## every declared scheme must construct and print; a scheme that exists as a type but cannot be shown is the
    ## defect that hid for months in MultiPhotonTransition
    scheme = Liouville.StimulatedRamanScheme()
    if  length(sprint(show, scheme))   == 0     success = false;   println("** empty show for StimulatedRamanScheme")   end
    if  length(Base.string(scheme))    == 0     success = false;   println("** empty string for StimulatedRamanScheme")  end
    ###
    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Scaffold checks only for Liouville; NO physics is compared against approved data.")
    testPrint("testModule_Liouville()::", success)
    return( success )
end
