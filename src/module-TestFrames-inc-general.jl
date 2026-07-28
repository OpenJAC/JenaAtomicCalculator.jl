
## Functions in this file cover: methods, structs, evaluations, representations, amplitudes.
## Alphabetical order within this file.


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

    # Kept deliberately small (2 layers, single EOL target level): a real per-layer EOL SCF optimization
    # runs at every layer, and a 3-layer, multi-target-level scope was found to run for many minutes -- far
    # too slow for a routine regression test. The 3-layer, multi-level scope is instead exercised (and
    # independently verified) by examples/example-Ai.jl Scenario B.
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1]) )
    coreShells  = [Shell("1s")]
    fromShells  = [Shell("2s")]
    layers      = [ RasLayer(Shell[]; se=false, de=false),   # layer 1: reference SCF only, no correlation
                     RasLayer([Shell("2p")]) ]               # layer 2: add 2p as a correlating shell
    #
    wa          = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs,
                                 RasExpansion([LevelSymmetry(0, Basics.plus)], 4, coreShells, fromShells, layers, rasSettings) )
    wb = generate(wa, output=true)
    if  abs(wb["step2"].levels[1].energy + 14.589901130961)  > 1.0e-3
        success = false
        if printTest   info(iostream, "levels[1].energy $(wb["step2"].levels[1].energy) != -14.589901130961")   end
    end

    testPrint("testRepresentation_RasExpansion()::", success)
    return(success)
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
