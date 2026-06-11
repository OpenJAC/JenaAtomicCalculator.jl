
"""
`module  JAC.TestFrames`  
... a submodel of JAC that contains all methods for testing individual methods, modules, ....
"""
module TestFrames


using  Printf, SymEngine, JLD2, JenaAtomicCalculator, 
       ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

export testDummy

    
"""
`TestFrames.testCompareLines(oldLine, newLine; rtol=1.0e-6)`
    ... compares two output lines token by token: tokens that are parseable as Float64 are compared
        with relative tolerance rtol; all other tokens are compared exactly as strings. Returns true
        if every token agrees within tolerance, false otherwise.
"""
function testCompareLines(oldLine::String, newLine::String; rtol::Float64=1.0e-6)
    oldTokens = split(strip(oldLine))
    newTokens = split(strip(newLine))
    length(oldTokens) != length(newTokens)  &&  return( false )
    for (ot, nt) in zip(oldTokens, newTokens)
        ov = tryparse(Float64, ot)
        nv = tryparse(Float64, nt)
        if  ov !== nothing  &&  nv !== nothing
            denom = max(abs(ov), abs(nv), 1.0e-100)
            abs(ov - nv) / denom > rtol  &&  return( false )
        else
            ot != nt  &&  return( false )
        end
    end
    return( true )
end


function testCompareFiles(fold::String, fnew::String, sa::String, noLines::Int64; rtol::Float64=1.0e-6)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # Compare the test computations with previous results
    oldLines = readlines(fold)
    newLines = readlines(fnew)
    iold = 0;  for i=1:length(oldLines)   line = oldLines[i];   if  occursin(sa, line)  iold = i;   break   end   end
    inew = 0;  for i=1:length(newLines)   line = newLines[i];   if  occursin(sa, line)  inew = i;   break   end   end
    if  iold == 0   ||   inew == 0    success = false
        println(stdout, "Tries to compare two inappropriate files fold = $(fold); fnew =$(fnew) on string $sa  ($iold, $inew)")
        if printTest   println(iostream, "Tries to compare two inappropriate files fold = $(fold); fnew =$(fnew) on string $sa  ($iold, $inew)")  end
        return( success )
    end
    #
    ii = inew + 1
    for  i = iold+2:iold+noLines
        ii = ii + 1
        if   length(oldLines[i]) < 5    continue    end
        if   !testCompareLines(oldLines[i], newLines[ii]; rtol=rtol)    success = false
            if  printTest  println(iostream, "    *** Old::  " * oldLines[i])
                           println(iostream, "    *** New::  " * newLines[ii])
            end
        end
    end

    return( success )
end


function testPrint(sa::String, success::Bool)
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    ok(succ) =  succ ? "[OK]" : "[Fail]"
    sb = sa * TableStrings.hBlank(110);   sb = sb[1:100] * ok(success);    println(iostream, sb)
    return( nothing )
end


"""
`TestFrames.testStructConstructors(; short::Bool=true)`
    ... tests that all major struct constructors and Base.show methods work without error.
        This fast structural test catches breakage caused by adding or removing struct fields
        without updating all constructor call sites. No physics computations are performed.
"""
function testStructConstructors(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    items = Tuple{String,Function}[
        # Core structs: Nuclear.Model with all simple constructors
        ("Nuclear.Model(Z)",                        () -> Nuclear.Model(26.)                            ),
        ("Nuclear.Model(Z,M)",                      () -> Nuclear.Model(26., 55.845)                    ),
        ("Nuclear.Model(Z,model)",                  () -> Nuclear.Model(26., "Fermi")                   ),
        ("Nuclear.Isomer()",                        () -> Nuclear.Isomer()                              ),
        # Radial grid
        ("Radial.Grid(false)",                      () -> Radial.Grid(false)                            ),
        # Configurations, shells
        ("Configuration(string)",                   () -> Configuration("1s^2 2s^2 2p^6")              ),
        ("Shell(string)",                           () -> Shell("2p")                                   ),
        ("Subshell(string)",                        () -> Subshell("2p_1/2")                            ),
        # General settings structs
        ("AsfSettings()",                           () -> AsfSettings()                                 ),
        ("LevelSelection()",                        () -> LevelSelection()                              ),
        ("LineSelection()",                         () -> LineSelection()                               ),
        # Atomic property settings
        ("Einstein.Settings()",                     () -> Einstein.Settings()                           ),
        ("Hfs.Settings()",                          () -> Hfs.Settings()                               ),
        ("IsotopeShift.Settings()",                 () -> IsotopeShift.Settings()                      ),
        ("LandeZeeman.Settings()",                  () -> LandeZeeman.Settings()                       ),
        ("StarkShift.Settings()",                   () -> StarkShift.Settings()                        ),
        ("AlphaVariation.Settings()",               () -> AlphaVariation.Settings()                    ),
        ("FormFactor.Settings()",                   () -> FormFactor.Settings()                        ),
        ("DecayYield.Settings()",                   () -> DecayYield.Settings()                        ),
        ("MultipolePolarizibility.Settings()",      () -> MultipolePolarizibility.Settings()           ),
        ("ReducedDensityMatrix.Settings()",         () -> ReducedDensityMatrix.Settings()              ),
        ("RadiativeOpacity.Settings()",             () -> RadiativeOpacity.Settings()                  ),
        # Basic process settings
        ("PhotoEmission.Settings()",                () -> PhotoEmission.Settings()                     ),
        ("PhotoExcitation.Settings()",              () -> PhotoExcitation.Settings()                   ),
        ("PhotoIonization.Settings()",              () -> PhotoIonization.Settings()                   ),
        ("PhotoRecombination.Settings()",           () -> PhotoRecombination.Settings()                ),
        ("AutoIonization.Settings()",               () -> AutoIonization.Settings()                    ),
        ("ElectronCapture.Settings()",              () -> ElectronCapture.Settings()                   ),
        ("DielectronicRecombination.Settings()",    () -> DielectronicRecombination.Settings()         ),
        ("PhotoExcitationFluores.Settings()",       () -> PhotoExcitationFluores.Settings()            ),
        ("PhotoExcitationAutoion.Settings()",       () -> PhotoExcitationAutoion.Settings()            ),
        ("RayleighCompton.Settings()",              () -> RayleighCompton.Settings()                   ),
        ("ParticleScattering.Settings()",           () -> ParticleScattering.Settings()                ),
        ("BeamPhotoExcitation.Settings()",          () -> BeamPhotoExcitation.Settings()               ),
        ("HyperfineInduced.Settings()",             () -> HyperfineInduced.Settings()                  ),
        ("ResonantInelastic.Settings()",            () -> ResonantInelastic.Settings()                 ),
        ("ImpactExcitation.Settings()",             () -> ImpactExcitation.Settings()                  ),
        ("CoulombExcitation.Settings()",            () -> CoulombExcitation.Settings()                 ),
        # Advanced process settings
        ("MultiPhotonDeExcitation.Settings()",      () -> MultiPhotonDeExcitation.Settings()           ),
        ("MultiPhotonIonization.Settings()",        () -> MultiPhotonIonization.Settings()             ),
        ("MultiPhotonDoubleIon.Settings()",         () -> MultiPhotonDoubleIon.Settings()              ),
        ("PhotoDoubleIonization.Settings()",        () -> PhotoDoubleIonization.Settings()             ),
        ("PhotoIonizationFluores.Settings()",       () -> PhotoIonizationFluores.Settings()            ),
        ("PhotoIonizationAutoion.Settings()",       () -> PhotoIonizationAutoion.Settings()            ),
        ("ImpactExcitationAutoion.Settings()",      () -> ImpactExcitationAutoion.Settings()           ),
        ("RadiativeAuger.Settings()",               () -> RadiativeAuger.Settings()                    ),
        ("InternalConversion.Settings()",           () -> InternalConversion.Settings()                ),
        ("InternalRecombination.Settings()",        () -> InternalRecombination.Settings()             ),
        ("TwoElectronOnePhoton.Settings()",         () -> TwoElectronOnePhoton.Settings()              ),
        ("DoubleAutoIonization.Settings()",         () -> DoubleAutoIonization.Settings()              ),
        # Further module settings
        ("Plasma.Settings()",                       () -> Plasma.Settings()                            ),
        ("Liouville.Settings()",                    () -> Liouville.Settings()                         ),
        ("StrongField.Settings()",                  () -> StrongField.Settings()                       ),
    ]

    for (name, fn) in items
        try
            obj = fn()
            show(devnull, obj)
        catch e
            success = false
            println(stdout,   "    *** FAIL $name: $(sprint(showerror, e))")
            if printTest   println(iostream, "    *** FAIL $name: $(sprint(showerror, e))")   end
        end
    end

    testPrint("testStructConstructors()::", success)
    return( success )
end


"""
`TestFrames.testEvaluation_Wigner_3j_specialValues(; short::Bool=true)`
    ... tests on special values for the Wigner 3j symbols.
"""
function testEvaluation_Wigner_3j_specialValues(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    w3j = RacahAlgebra.selectW3j(2);                    println(">> w3j-original     = $w3j")
    wa  = RacahAlgebra.symmetricForms(w3j)
    wb  = RacahAlgebra.evaluate(wa[1], special=true);   println(">> wb-special value  = $wb")
    wc  = RacahAlgebra.evaluate(wa[2], special=true);   println(">> wc-special value  = $wc")
    if  wb != wc
        success = false
        if printTest   info(iostream, "$w3j:   $wb != $wc")   end
    end

    testPrint("testEvaluation_Wigner_3j_specialValues()::", success)
    return(success)  
end


"""
`TestFrames.testEvaluation_Wigner_6j_specialValues(; short::Bool=true)`  
    ... tests on special values for the Wigner 6j symbols.
"""
function testEvaluation_Wigner_6j_specialValues(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    w6j = RacahAlgebra.selectW6j(2);                    println(">> w6j-original      = $w6j")
    wa  = RacahAlgebra.symmetricForms(w6j)
    wb  = RacahAlgebra.evaluate(wa[1], special=true);   println(">> wb-special value  = $wb")
    wc  = RacahAlgebra.evaluate(wa[2], special=true);   println(">> wc-special value  = $wc")
    if  wb != wc
        success = false
        if printTest   info(iostream, "$w6j:   $wb != $wc")   end
    end

    testPrint("testEvaluation_Wigner_6j_specialValues()::", success)
    return(success)  
end


"""
`TestFrames.testEvaluation_Wigner_9j_specialValues(; short::Bool=true)`  
    ... tests on special values for the Wigner 9j symbols.
"""
function testEvaluation_Wigner_9j_specialValues(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    w9j = RacahAlgebra.selectW9j(2);                     println(">> w9j-original     = $w9j")
    wa  = RacahAlgebra.symmetricForms(w9j)
    wb  = RacahAlgebra.evaluate(wa[5],  special=true);   println(">> wb-special value = $wb")
    wc  = RacahAlgebra.evaluate(wa[11], special=true);   println(">> wc-special value = $wc")
    if  wb != wc
        success = false
        if printTest   info(iostream, "$w9j:   $wb != $wc")   end
    end

    testPrint("testEvaluation_Wigner_9j_specialValues()::", success)
    return(success)  
end


"""
`TestFrames.testEvaluation_sumRulesForOneWnj(; short::Bool=true)`  
    ... tests on special values for the Wigner 3-j symbols.
"""
function testEvaluation_sumRulesForOneWnj(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    rex = RacahAlgebra.selectRacahExpression(1);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(2);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(3);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(4);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(5);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end
    

    testPrint("testEvaluation_sumRulesForOneWnj()::", success)
    return(success)  
end


"""
`TestFrames.testEvaluation_sumRulesForTwoWnj(; short::Bool=true)`  
    ... tests on special values for the Wigner 3-j symbols.
"""
function testEvaluation_sumRulesForTwoWnj(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    rex = RacahAlgebra.selectRacahExpression(6);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(7);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end

    rex = RacahAlgebra.selectRacahExpression(8);         println(">> rex-original  = $rex")
    wa  = RacahAlgebra.evaluate(rex);                    println(">> rex-evaluated = $wa")
    if  isnothing(wa)
        success = false
        if printTest   info(iostream, "No simplification found for $rex")   end
    end
    

    testPrint("testEvaluation_sumRulesForTwoWnj()::", success)
    return(success)  
end



"""
`Basics.testMethod_integrate_ongrid(; short::Bool=true)`  ... tests the integration on grid.
"""
function testMethod_integrate_ongrid(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # Test the integration on the grid for an analytical function
    grid = Radial.Grid(Radial.Grid(false), rnt=2.0e-10, h=5.0e-3, NoPoints=9000)
    function f1(r)    
        A = 10^5;    gamma = 3;    a = 300
        return( A * r^gamma * exp(-a * r^2) )
    end
    exact1 = 5/9

    integrand = Vector{Float64}(size(grid.r, 1))
    for i = 1:size(grid.r, 1)
        integrand[i] = f1(grid.r[i])
    end
    integrand = integrand .* grid.rp[1:size(integrand, 1)]   # adopt to the form of Grasp92

    integral  = Basic.integrate(NewtonCotes(), integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... Newton-Cotes:  I = $integral,  Err = $err")  end
    end

    integral  = Basic.integrate(SimpsonRule(), integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... Simpson rule:  I = $integral,  Err = $err")  end
    end

    integral  = Basic.integrate(TrapezRule(), integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... trapez rule:  I = $integral,  Err = $err")  end
    end
    
    println(iostream, "Warning(testMethod_integrate_ongrid): test of integration with Grasp orbitals has been set silent.")
    #=
    # Test the integration with orbital functions from Grasp92
    grid = Radial.Grid(true)

    orbitals1 = Basics.readOrbitalFileGrasp92("../test/approved/Ne-0+-scf.exp.out", grid)
    orbitals2 = Basics.readOrbitalFileGrasp92("../test/approved/Ne-1+-scf.exp.out", grid)

    for i = 1:size(orbitals1, 1)
        for j = 1:size(orbitals2, 1)
            orb1      = orbitals1[i];   orb2 = orbitals2[j]
            if   orb1.subshell.kappa != orb2.subshell.kappa    break   end

            mtp       = min(size(orb1.P, 1), size(orb2.P, 1))
            integrand = ( orb1.P[1:mtp] .* orb2.P[1:mtp] + orb1.Q[1:mtp] .* orb2.Q[1:mtp] ) .* grid.rp[1:mtp]
    
            integrala = Basic.integrate(NewtonCotes(), integrand, grid)^2
            integralb = Basic.integrate(SimpsonRule(), integrand, grid)^2
            integralc = Basic.integrate(TrapezRule(),  integrand, grid)^2
    
            info(iostream, "<$(string(orb1.subshell)) | $(string(orb2.subshell))> = $integrala, $integralb, $integralc")
        end
    end  =#

    testPrint("testMethod_integrate_ongrid()::", success)
    return(success)  
end


"""
`TestFrames.testMethod_Wigner_3j(; short::Bool=true)`  ... tests on Wigner 3j symbols.
"""
function testMethod_Wigner_3j(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # Wigner_3j(1,2,1,0,0,0) = 0.36514837167011074230
    a = c = AngularJ64(1);    b = AngularJ64(2);    ma = mb = mc = AngularM64(0)
    wa = AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.36514837167011074230) > 1.0e-12
        success = false
        if printTest   info(iostream, "Wigner_3j(1,2,1,0,0,0) = 0.36514837167011074230 ... but obtains value = $wa")   end
    end

    # Wigner_3j(3,6,3,0,0,0) = 0.182482967150452976281
    a = c = AngularJ64(3);    b = AngularJ64(6);    ma = mb = mc = AngularM64(0)
    wa = AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.18248296715045297628) > 1.0e-12
        success = false
        if printTest   info(iostream, "Wigner_3j(3,6,3,0,0,0) = 0.182482967150452976281 ... but obtains value = $wa")   end
    end

    # Wigner_3j(1/2,1,1/2,1/2,0,-1/2) = 0.40824829046386301637
    a = c = AngularJ64(1//2);    b = AngularJ64(1);    ma =  AngularM64(1//2);   mb =  AngularM64(0);    mc = AngularM64(-1//2)
    wa = AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.40824829046386301637) > 1.0e-12    
        success = false
        if printTest   info(iostream, "Wigner_3j(1/2,1,1/2,1/2,0,-1/2) = 0.40824829046386301637 ... but obtains value = $wa")  end
    end

    testPrint("testMethod_Wigner_3j()::", success)
    return(success)  
end



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
`TestFrames.testModule_Cascade_PhotonExcitation(; short::Bool=true)`  ... tests on module Cascade.
"""
function testModule_Cascade_PhotonExcitation(; short::Bool=true) 
    Defaults.setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
    Defaults.setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb") 
    Defaults.setDefaults("print summary: open", "test-Cascade-PhotonExcitation-new.sum")
    printstyled("\n\nTest the module  Cascade for the PhotonExcitationScheme ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)
    name = "Photoabsorption calculations for Ne^+ for energies = [1, 4] a.u."
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotoExcitationScheme([E1], 0.5, 4.0, 1, [Shell("2s"), Shell("2p")], 
                                                                     [Shell("2s"), Shell("2p"), Shell("3p"), Shell("4p"), Shell("5p")],
                                                                     LevelSelection(), [0,1], 0., 0.),
                                                                     initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    wb = perform(wa; output=true)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Cascade-PhotonExcitation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Cascade-PhotonExcitation-new.sum"), "Steps that are defined for the", 15) 
    testPrint("testModule_Cascade-PhotonExcitation()::", success)
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
`TestFrames.testModule_DielectronicRecombination(; short::Bool=true)`  ... tests on module DielectronicRecombination.
"""
function testModule_DielectronicRecombination(; short::Bool=true)     Defaults.setDefaults("print summary: open", "test-DielectronicRecombination-new.sum")
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
                            "Partial (Auger) capture", 10) 
    testPrint("testModule_DielectronicRecombination()::", success)
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
`TestFrames.testModule_MultipoleMoment(; short::Bool=true)`  ... tests on module MultipoleMoment.
"""
function testModule_MultipoleMoment(; short::Bool=true) 
    Defaults.setDefaults("print summary: open", "test-MultipoleMoment-new.sum")
    ### Make the tests
    printstyled("\n\nTest the module  MultipoleMoment  ... \n", color=:cyan)
    grid = Radial.Grid(true)
    wa = Atomic.Computation(Atomic.Computation(), 
                            name="xx",  nuclearModel=Nuclear.Model(26.), grid=grid, 
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")], printout=true )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    flow = 6;    fup = 8;   ilow = 1;   iup = 3
    println("\n\nDipole amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                MultipoleMoment.dipoleAmplitude(finalLevel, initialLevel, grid; display=true)
            end
        end
    end
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultipoleMoment-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-MultipoleMoment-new.sum"), "Dipole amplitude", 2) 
    testPrint("testModule_MultipoleMoment()::", success)
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



"""
`TestFrames.testModule_ParityNonConservation(; short::Bool=true)`  ... tests on module ParityNonconservation.
"""
function testModule_ParityNonConservation(; short::Bool=true) 
    Defaults.setDefaults("print summary: open", "test-ParityNonConservation-new.sum")
    printstyled("\n\nTest the module  ParityNonConservation  ... \n", color=:cyan)
    ### Make the tests
    grid = Defaults.getDefaults("standard grid")
    wa = Atomic.Computation(Atomic.Computation(), name="xx",  nuclearModel=Nuclear.Model(26.), 
                            grid=grid,
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")] )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    nModel = Nuclear.Model(26.);    flow = 6;    fup = 8;   ilow = 1;   iup = 3

    println("\n\nWeak-charge and Schiff-moment amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                ParityNonConservation.weakChargeAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
                ParityNonConservation.schiffMomentAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
            end
        end
    end
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-ParityNonConservation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-ParityNonConservation-new.sum"), "weak-charge amplitude", 12) 
    testPrint("testModule_ParityNonConservation()::", success)
    return(success)  
end


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
                            processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, true, 
                                                                        LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes(0., 0., 0.) ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitation-new.sum"),  
                            "Photoexcitation cross sections for (completely)", 200) 
    testPrint("testModule_PhotoExcitation()::", success)
    return(success)  
end



"""
`TestFrames.testModule_PhotoExcitationAutoion(; short::Bool=true)`  ... tests on module PhotoExcitationAutoion.
"""
function testModule_PhotoExcitationAutoion(; short::Bool=true) 
    Defaults.setDefaults("print summary: open", "test-PhotoExcitationAutoion-new.sum")
    printstyled("\n\nTest the module  PhotoExcitationAutoion  ... \n", color=:cyan)
    ### Make the tests
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=grid, 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2"), Configuration("1s 2p^2")],
                            finalConfigs=[Configuration("1s^2")], 
                            processSettings=PhotoExcitationAutoion.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, 
                                                                            PathwaySelection(true, indexTriples=[(1,1,1)]), 2 ) )

wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-PhotoExcitationAutoion-new.sum")
    success = true
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitationAutoion-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitationAutoion-new.sum"), "xxx", 10) 
    testPrint("testModule_PhotoExcitationAutoion()::", success)
    return(success)  
end



"""
`TestFrames.testModule_PhotoExcitationFluores(; short::Bool=true)`  ... tests on module PhotoExcitationFluores.
"""
function testModule_PhotoExcitationFluores(; short::Bool=true) 
    Defaults.setDefaults("print summary: open", "test-PhotoExcitationFluores-new.sum")
    printstyled("\n\nTest the module  PhotoExcitationFluores  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                            processSettings=PhotoExcitationFluores.Settings([E1, M1], [UseCoulomb, UseBabushkin],
                                                    false, false, false, true, ExpStokes(), SolidAngle[], 0.,
                                                    PathwaySelection(true, indexTriples=[(1,1,1)]) )  )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-PhotoExcitationFluores-new.sum")
    success = true
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitationFluores-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitationFluores-new.sum"),  "xxx", 10) 
    testPrint("testModule_PhotoExcitationFluores()::", success)
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
                                LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)]), 0., 0., 10000. ) )
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
`TestFrames.testRepresentation_MeanFieldBasis_CiExpansion(; short::Bool=true)`  ... tests on the representation .
"""
function testRepresentation_MeanFieldBasis_CiExpansion(; short::Bool=true) 
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
    refConfigs  = [Configuration("[He] 2s^2 2p^4")]
    mfSettings  = MeanFieldSettings()
    #
    wa          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, MeanFieldBasis(mfSettings) )
    wb = generate(wa, output=true)
    #
    orbitals    = wb["mean-field basis"].orbitals
    ciSettings  = CiSettings(CoulombInteraction(), LevelSelection() )
    from        = [Shell("2s")]
    #
    frozen      = [Shell("1s")]
    to          = [Shell("2s"), Shell("2p")]
    excitations = RasStep()
    #             RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    wc          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                                    CiExpansion(orbitals, excitations, ciSettings) )
    println("wc = $wc")
    wd = generate(wc, output=true)
    
    if  abs(orbitals[Subshell("1s_1/2")].energy + 18.705283) > 1.0e-3
        success = false
        if printTest   @info(iostream, "orbital energy $(orbitals[Subshell("1s_1/2")].energy) != -18.705283")     end
        @info(iostream, "orbital energy $(orbitals[Subshell("1s_1/2")].energy) != -18.705283")
    end
    if  abs(wd["CI multiplet"].levels[1].energy + 74.840309)  > 1.0e-2
        success = false
        if printTest   @info(iostream, "levels[1].energy $(wd["CI multiplet"].levels[1].energy) != -74.840309")   end
        @info(iostream, "levels[1].energy $(wd["CI multiplet"].levels[1].energy) != -74.840309")
    end

    testPrint("testRepresentation_MeanFieldBasis_CiExpansion()::", success)
    return(success)  
end



"""
`TestFrames.testRepresentation_RasExpansion(; short::Bool=true)`  ... tests on the representation .
"""
function testRepresentation_RasExpansion(; short::Bool=true) 
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1,2,3]) )
    from        = [Shell("2s")]
    #
    frozen      = [Shell("1s")]
    to          = [Shell("2s"), Shell("2p")]
    step1       = RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    append!(frozen, [Shell("2s"), Shell("2p")])
    append!(to,     [Shell("3s"), Shell("3p"), Shell("3d")])
    step2       = RasStep(step1; seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    append!(frozen, [Shell("3s"), Shell("3p"), Shell("3d")])
    append!(to,     [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
    step3       = RasStep(step2, seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    wa          = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs, 
                                 RasExpansion([LevelSymmetry(0, Basics.plus)], 4, [step1, step2, step3], rasSettings) )
    wb = generate(wa, output=true)
    if  abs(wb["step3"].levels[1].energy + 14.7540801622)  > 1.0e-3
        success = false
        if printTest   info(iostream, "levels[1].energy $(wd["CI multiplet"].levels[1].energy) != 14.754080162")   end
    end

    testPrint("testRepresentation_RasExpansion()::", success)
    return(success)  
end



"""
`TestFrames.testRepresentation_GreenExpansion(; short::Bool=true)`  ... tests on the representation .
"""
function testRepresentation_GreenExpansion(; short::Bool=true) 
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    
    name          = "Lithium 1s^2 2s ground configuration"
    refConfigs    = [Configuration("[He] 2s")]
    greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, LevelSelection() )
    #
    wa          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                                    ## GreenExpansion( AtomicState.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                    ## GreenExpansion( AtomicState.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                    GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                                    [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)], 3, greenSettings) )
    wb = generate(wa, output=true)

    if  abs(wb["Green channels"][1].gMultiplet.levels[1].energy + 64.080705)  > 1.0e-3
        success = false
        if printTest   info(iostream, "gMultiplet.levels[1].energy $(wb["Green channels"][1].gMultiplet.levels[1].energy) != -64.080705")   end
    end

    testPrint("testRepresentation_GreenExpansion()::", success)
    return(success)  
end

end # module
