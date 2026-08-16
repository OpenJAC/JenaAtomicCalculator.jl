
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
`TestFrames.testMethod_HydrogenicRates(; short::Bool=true)`
    ... known-answer test of DielectronicRecombination.computeHydrogenicRate(ni, li, nf, lf, Zeff), the
        Infeld-Hull recursion that supplies the whole physical content of the HydrogenicCorrections and of the
        Rydberg-tail extrapolation: the radiative rate by which a captured Rydberg spectator stabilizes a DR
        resonance through decay channels not represented in the final multiplet.

        The reference values are NOT stored. They are recomputed here from an INDEPENDENT implementation --
        the hydrogenic radial dipole integral evaluated by direct Simpson integration over the analytic
        Schroedinger radial functions (generalized Laguerre polynomials by recursion) -- which shares nothing
        with the Infeld-Hull scheme except the physics. That is why this test needs no approved .sum file: it
        compares a number against a number computed on the spot, not against stored output.

        Three checks: the independent sweep over all E1-allowed transitions with n <= 6; the exact Z^4 scaling
        of a hydrogenic rate; and the guards, which must return exactly zero for an upward step, for
        Delta l != 1 and for ni == nf rather than a small wrong number.
"""
function testMethod_HydrogenicRates(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    alpha = Defaults.getDefaults("alpha")

    ## --- the independent implementation -------------------------------------------------------------------
    logfact(k::Int64) = k <= 1 ? 0.0 : sum(log, 1:k)
    ## generalized Laguerre L_k^a(x) by the standard three-term recursion
    function laguerreGen(k::Int64, a::Int64, x::Float64)
        if      k == 0    return( 1.0 )
        elseif  k == 1    return( 1.0 + a - x )
        end
        lm1 = 1.0;    l0 = 1.0 + a - x
        for  j = 1:k-1
            lp1 = ( (2j + 1 + a - x) * l0 - (j + a) * lm1 ) / (j + 1);    lm1 = l0;    l0 = lp1
        end
        return( l0 )
    end
    ## normalized hydrogenic radial function R_nl(r), nuclear charge Z, atomic units
    function radialR(n::Int64, l::Int64, Z::Float64, r::Float64)
        rho  = 2Z * r / n
        logN = 3*log(2Z/n) + logfact(n-l-1) - log(2n) - logfact(n+l)
        return( exp(0.5*logN) * exp(-rho/2) * rho^l * laguerreGen(n-l-1, 2l+1, rho) )
    end
    ## A(ni li --> nf lf) = 4/3 alpha^3 omega^3 max(li,lf)/(2li+1) |<r>|^2   [atomic units]
    function rateIndependent(ni::Int64, li::Int64, nf::Int64, lf::Int64, Z::Float64)
        if  abs(li-lf) != 1  ||  ni <= nf  ||  li >= ni  ||  lf >= nf  ||  li < 0  ||  lf < 0    return( 0.0 )   end
        rmax = 40.0 * ni^2 / Z;    npts = 200001;    h = rmax / (npts - 1);    s = 0.0
        for  k = 0:npts-1
            r = k * h
            f = radialR(ni, li, Z, r) * radialR(nf, lf, Z, r) * r^3
            w = (k == 0 || k == npts-1) ? 1.0 : (isodd(k) ? 4.0 : 2.0)
            s = s + w * f
        end
        rint  = s * h / 3
        omega = Z^2 / 2 * (1/nf^2 - 1/ni^2)
        return( 4/3 * alpha^3 * omega^3 * max(li, lf) / (2li + 1) * rint^2 )
    end

    ## --- (1) sweep against the independent implementation --------------------------------------------------
    nBad = 0;   worst = 0.0;   worstLabel = ""
    for  ni = 2:6,  li = 0:ni-1,  nf = 1:ni-1,  lf = 0:nf-1
        if  abs(li - lf) != 1    continue    end
        aJAC = DielectronicRecombination.computeHydrogenicRate(ni, li, nf, lf, 1.0)
        aInd = rateIndependent(ni, li, nf, lf, 1.0)
        if  aInd <= 0.0    continue    end
        dev = abs(aJAC/aInd - 1.0)
        if  dev > worst    worst = dev;    worstLabel = "($ni,$li -> $nf,$lf)"    end
        if  dev > 1.0e-3   nBad = nBad + 1    end
    end
    if  nBad > 0
        success = false
        if printTest   info(iostream, "computeHydrogenicRate: $nBad transitions deviate by more than 0.1 % from an " *
                                      "independent radial-dipole integration; worst $worstLabel at $(worst*100) %")   end
    end

    ## --- (2) exact Z^4 scaling ------------------------------------------------------------------------------
    a1 = DielectronicRecombination.computeHydrogenicRate(2, 1, 1, 0, 1.0)
    for  Z in [2.0, 6.0, 18.0, 54.0, 75.0]
        ratio = DielectronicRecombination.computeHydrogenicRate(2, 1, 1, 0, Z) / (a1 * Z^4)
        if  abs(ratio - 1.0) > 1.0e-10
            success = false
            if printTest   info(iostream, "computeHydrogenicRate: A(2p->1s) at Z = $Z is not Z^4 times its Z = 1 " *
                                          "value; ratio = $ratio")   end
        end
    end

    ## --- (3) the guards must return exactly zero ------------------------------------------------------------
    for  (ni, li, nf, lf, why) in [ (2, 1, 3, 0, "upward step ni < nf"), (3, 0, 2, 0, "Delta l = 0"),
                                    (3, 2, 2, 0, "Delta l = 2"),         (2, 1, 2, 0, "ni == nf") ]
        wa = DielectronicRecombination.computeHydrogenicRate(ni, li, nf, lf, 1.0)
        if  wa != 0.0
            success = false
            if printTest   info(iostream, "computeHydrogenicRate: the guard for '$why' returned $wa instead of 0.0")  end
        end
    end

    testPrint("testMethod_HydrogenicRates()::", success)
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
                                joinpath(@__DIR__, "..", "test", "test-MultipoleMoment-new.sum"), "Dipole amplitude", 8)
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
    # The expected value changed TWICE on 16-Aug-2026. First when Basics.EOLField was wired to the
    # orbital-rotation solver, and again when computeOrbitalGradient was corrected: it contracted a screened
    # potential built from the sign-canonicalized ORBITALS with a RAW b-vector, so an off-diagonal CSF pair
    # carried an odd power of a correlating orbital's sign and the gradient drifted from the functional --
    # 2.5x, then 19x, then SIGN-WRONG as the 2p weight grew, which stalled the line search at iteration 3 of
    # 24.  With the full s/wN factor on every slot the analytic slope matches a finite difference to five
    # digits, the driver runs to its iteration limit, and this case reaches -14.617216 instead of -14.612567.  The former -14.589901130961 was the COLLAPSED answer: with the old solver this very case
    # converged cleanly, to its own 1e-8, onto a degenerate stationary point whose mixing vector was
    # [0.971, -0.0000482, 0.238] -- the 2p_3/2 correlation channel eliminated outright.  An independent
    # DFS-Field run of the identical layer structure gave -14.605300 Ha with both channels contributing.
    # The rotation solver reaches -14.612567, i.e. BELOW that reference, so this test previously asserted
    # the defect.  Do not restore the old number.
    if  abs(wb["step2"].levels[1].energy + 14.617216447357)  > 1.0e-3
        success = false
        if printTest   info(iostream, "levels[1].energy $(wb["step2"].levels[1].energy) != -14.617216447357")   end
    end

    testPrint("testRepresentation_RasExpansion()::", success)
    return(success)
end


"""
`TestFrames.testMethod_SettingsCopyConstructors(; short::Bool=true)`
    ... exercises the keyword copy-constructor of EVERY submodule Settings that has one; returns success::Bool.
"""
function testMethod_SettingsCopyConstructors(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-SettingsCopyConstructors-new.sum")
    printstyled("\n\nTest all Settings copy-constructors  ... \n", color=:cyan)
    ### Make the tests
    ## WHY THIS TEST EXISTS (added 09-Aug-2026). JAC has more than 400 keyword copy-constructor guards written in
    ## the aligned form
    ##      if  isnothing(a)   ax = set.a   else   ax = a   end
    ## and that alignment, which makes them readable, also makes a one-character slip invisible. Four such slips
    ## were live when this test was written, each broken in a DIFFERENT way and none caught by anything:
    ##   * CoulombExcitation   -- the if-branch assigned the parameter `ionEnergies` rather than the local
    ##                            `ionEnergiesx` the constructor reads: raised whenever the keyword was OMITTED.
    ##   * ParticleScattering  -- the guard tested `processTypey`, a name defined nowhere: raised ALWAYS.
    ##   * PhotoDoubleIonization -- one guard read `set.electronEnergies`, a field the struct does not have; a
    ##                            second assigned `maxKappasx` where the constructor reads `maxKappax`, so it
    ##                            raised whenever that keyword WAS given.
    ##   * BeamPhotoExcitation -- tested `x == missing`, which evaluates to `missing` and not to a Bool, so the
    ##                            `if` itself raised; two keyword types were copied from their neighbours as well.
    ## A line-by-line scan found the first three; the last was found only by CALLING every constructor, which is
    ## what this test does. Reflection is used deliberately: a hand-written list would go stale the moment a new
    ## module is added, and this class of defect is precisely the one nobody thinks to add a test for.
    success = true
    JAC     = JenaAtomicCalculator
    ## KNOWN FAILURE, EXCLUDED BY DECISION rather than skipped silently: Liouville's copy-constructor does not
    ## convert, but that module was postponed by the user on 09-Aug-2026 and is not to be touched. If it is ever
    ## revived, delete it from this list FIRST and let the test tell you what is wrong with it.
    excluded = [:Liouville]
    modules  = Module[]
    for  nm in names(JAC, all=true, imported=false)
        isdefined(JAC, nm) || continue
        val = getfield(JAC, nm)
        if  val isa Module  &&  val !== JAC  &&  !(nm in excluded)    push!(modules, val)    end
    end
    nTested = 0
    for  M in sort(modules, by=string)
        isdefined(M, :Settings) || continue
        S = getfield(M, :Settings)
        S isa Type || continue
        ## only those that offer both a zero-argument and a copy-constructor
        (hasmethod(S, Tuple{}) && hasmethod(S, Tuple{S})) || continue
        try
            set0 = Base.invokelatest(S)
            Base.invokelatest(S, set0)
            nTested = nTested + 1
        catch e
            success = false
            println("** $(M).Settings copy-constructor raises: " * first(split(sprint(showerror, e), "\n")))
        end
    end
    ## a floor, so that a refactor which silently stops finding the modules cannot pass unnoticed
    if  nTested < 30
        success = false;   println("** only $nTested Settings copy-constructors were found; expected at least 30")
    end
    println("    ... $nTested Settings copy-constructors exercised, $(length(excluded)) excluded by decision.")
    ###
    Defaults.setDefaults("print summary: close", "")
    testPrint("testMethod_SettingsCopyConstructors()::", success)
    return( success )
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
        ("Nuclear.Model(Z,model)",                  () -> Nuclear.Model(26., FermiNucleus())            ),
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
        ("MultiPhotonTransition.Settings()",      () -> MultiPhotonTransition.Settings()           ),
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
`TestFrames.testMethod_OrbitalOrthonormality(; short::Bool=true)`  
    ... asserts that a converged self-consistent field returns an ORTHONORMAL orbital set, which the CSF
        expansion assumes and which nothing else in the suite checks. Needs no reference data: the
        requirement is exact, so the test is against zero rather than against a tabulated number.

        This is the check that would have caught the defect found on 10-Aug-2026, where the damping step
        `mixed = 0.5*old + 0.5*raw` destroyed the orthogonality that Hamiltonian.projectHamiltonian had
        just enforced, and nothing restored it. Converged Li 1s^2 2s + 1s^2 3s + 1s^2 3d gave
        <2s|3s> = -1.128e-03 while every one of the 44 approved tests passed. A same-kappa block with
        THREE or more orbitals is what exposes it -- two orbitals only reach ~5e-05 -- so the case below
        is chosen for its three s-orbitals, not for its physics.
        Returns true if the worst same-kappa overlap deviates from delta_ab by less than 1.0e-08.
"""
function testMethod_OrbitalOrthonormality(; short::Bool=true)
    success = true
    printstyled("\n\nTest the orthonormality of a converged SCF orbital set: \n", color=:light_green)
    printstyled(  "-------------------------------------------------------- \n", color=:light_green)

    grid    = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 30.0)
    configs = [Configuration("1s^2 2s"), Configuration("1s^2 3s"), Configuration("1s^2 3d")]
    multiplet = SelfConsistent.performSCF(configs, Nuclear.Model(3.), grid,
                    AsfSettings(AsfSettings(); scField=Basics.ALField(), maxIterationsScf=24); printout=false)
    basis    = multiplet.levels[1].basis;    orbitals = basis.orbitals
    worst    = 0.;    worstPair = ""
    for  (i, sha)  in  enumerate(basis.subshells),  (j, shb)  in  enumerate(basis.subshells)
        if  sha.kappa != shb.kappa   ||   j < i    continue    end
        oa  = orbitals[sha];    ob = orbitals[shb]
        mtp = min( size(oa.P,1), size(ob.P,1), length(grid.wr) )
        ov  = sum( grid.wr[k] * (oa.P[k]*ob.P[k] + oa.Q[k]*ob.Q[k])  for k = 1:mtp )
        dev = abs( ov - (i == j ? 1.0 : 0.0) )
        if  dev > worst    worst = dev;   worstPair = "<$(string(sha)) | $(string(shb))>"    end
    end
    println("  Worst same-kappa deviation from orthonormality:  $worstPair = $worst ")
    if  worst > 1.0e-8
        success = false
        println("  *** The converged orbitals are NOT orthonormal; the CSF expansion assumes they are, so " *
                "the energies of this basis are not legitimate variational numbers.")
    end

    testPrint("testMethod_OrbitalOrthonormality()::", success)
    return( success )
end


"""
`TestFrames.testMethod_BreitInteraction(; short::Bool=true)`
    ... asserts that the Breit interaction is computed with all of its tensorial parts present. Nothing else
        in the suite does: no approved test exercises InteractionStrength.XL_Breit at all, so before this
        method the Breit interaction could be broken outright without a single test noticing. Needs no
        reference data -- every requirement below is either exact or a selection rule.

        This is the check that would have caught the defect found on 14-Aug-2026. Commit 8f0930b gave
        AngularMomentum.CL_reduced_me the parity selection rule it had genuinely been missing -- correct in
        itself, and needed by the Coulomb, SMS and electric-multipole call sites -- but XL_Breit_coefficients
        handed ONE angular prefactor to all four of its blocks, and they do not all want the same parity.
        The 'T' block at nu = L is guarded by l_a+l_c+L ODD while the shared prefactor was now non-zero only
        for EVEN; the two are mutually exclusive, so that block became unreachable. It carries the DOMINANT
        magnetic term, and Gaunt came out a factor 3.3 too small for four days while the whole suite stayed
        at 45/45. Four things are asserted:
        (1) the two parities are COMPLEMENTARY -- of <kappa_a||C^L||kappa_c> and <-kappa_a||C^L||kappa_c>
            exactly one is non-zero, and the non-zero one carries the full magnitude. Both branches are
            needed, because the nu = L integrand is a.P * c.Q, LARGE component against SMALL;
        (2) besselPhiPsi, on which both radial kernels rest, reproduces the elementary closed forms
            phi_0 = sin(x)/x, psi_0 = cos(x), phi_1 = 3(sin(x)/x - cos(x))/x^2, psi_1 = cos(x) + x sin(x);
        (3) for four 2p_3/2 orbitals at L = 1 every other block is vetoed by its own parity guard, so
            X^1_Breit consists of exactly eight ('T', nu = 1) coefficients and must NOT vanish. This case
            isolates the dominant magnetic term, and it was EXACTLY ZERO during the defect above;
        (4) the frequency correction is quadratic in the factor of CoulombBreit(factor), as an O(omega^2)
            retardation must be. A wrong kernel, a wrong power of omega or a mis-paired coefficient each
            break the quadratic law.
        Returns true if all four hold.
"""
function testMethod_BreitInteraction(; short::Bool=true)
    success = true
    printstyled("\n\nTest that the Breit interaction retains all of its tensorial parts: \n", color=:light_green)
    printstyled(  "------------------------------------------------------------------- \n", color=:light_green)

    ## (1) The even- and odd-parity reduced matrix elements must be complementary; XL_Breit_coefficients
    ##     needs BOTH, since its 'T' block at nu = L wants l_a+l_c+L odd and every other block wants it even.
    for  (sa, sc, L)  in  [ (Subshell("2p_3/2"), Subshell("2p_3/2"), 1), (Subshell("2p_3/2"), Subshell("2p_1/2"), 1),
                            (Subshell("3d_5/2"), Subshell("3p_3/2"), 1), (Subshell("3p_3/2"), Subshell("3p_1/2"), 2),
                            (Subshell("2p_1/2"), Subshell("2s_1/2"), 1) ]
        even  = AngularMomentum.CL_reduced_me(sa, L, sc)
        odd   = AngularMomentum.CL_reduced_me(Subshell(sa.n, -sa.kappa), L, sc)
        isOdd = isodd( Basics.subshell_l(sa) + Basics.subshell_l(sc) + L )
        wanted, other = isOdd ? (odd, even) : (even, odd)
        if  abs(other) > 1.0e-10
            success = false
            println("  *** <$(string(sa))||C^$L||$(string(sc))>: the branch that the parity rule must " *
                    "annihilate is non-zero ($other).")
        end
        if  abs(wanted) < 1.0e-10
            success = false
            println("  *** <$(string(sa))||C^$L||$(string(sc))>: BOTH parity branches vanish, so the Breit " *
                    "operator loses this term entirely.")
        end
    end
    println("  (1) The even- and odd-parity reduced matrix elements are complementary.")

    ## (2) besselPhiPsi carries both radial kernels; check it against the elementary closed forms.
    worst = 0.
    for  x  in  [0.05, 0.5, 1.0, 2.0]
        p0, q0 = InteractionStrength.besselPhiPsi(0, x)
        p1, q1 = InteractionStrength.besselPhiPsi(1, x)
        worst  = max( worst, abs(p0 - sin(x)/x), abs(q0 - cos(x)),
                             abs(p1 - 3*(sin(x)/x - cos(x))/x^2), abs(q1 - (cos(x) + x*sin(x))) )
    end
    println("  (2) besselPhiPsi against phi_0, psi_0, phi_1, psi_1:  worst deviation = $worst ")
    if  worst > 1.0e-12
        success = false
        println("  *** The normalised Bessel functions are wrong; both Breit radial kernels rest on them.")
    end

    ## (3) and (4) need one converged orbital set; Ne-like Fe is small and gives a Breit term of usable size.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    mp   = SelfConsistent.performSCF([Configuration("1s^2 2s^2 2p^6")], Nuclear.Model(26.), grid,
                                     AsfSettings(); printout=false)
    orbs = mp.levels[1].basis.orbitals
    p32  = orbs[Subshell("2p_3/2")];    s2 = orbs[Subshell("2s_1/2")]

    ## (3) Four 2p_3/2 orbitals at L = 1: l_a+l_c+L = 3 is odd, so the nu = L-1, nu = L+1 and 'S' blocks are
    ##     all vetoed and the dominant magnetic term stands alone.
    xcList = InteractionStrength.XL_Breit_coefficients(1, p32, p32, p32, p32)
    kinds  = unique([ (xc.kind, xc.nu) for xc in xcList ])
    xValue = InteractionStrength.XL_Breit(1, p32, p32, p32, p32, grid, Basics.CoulombGaunt())
    println("  (3) (2p_3/2)^4 at L = 1:  $(length(xcList)) coefficients, kinds = $kinds,  X^1_Breit = $xValue ")
    if  kinds != [('T',1)]  ||  length(xcList) != 8
        success = false
        println("  *** This case no longer isolates the nu = L block, so it no longer tests what it is for.")
    end
    if  abs(xValue) < 1.0e-6
        success = false
        println("  *** The DOMINANT magnetic term of the Breit interaction has vanished. This is exactly " *
                "the 10-Aug-2026 defect: a parity rule applied to an angular prefactor that the nu = L " *
                "block needs in its OTHER parity.")
    end

    ## (4) The retardation enters at O(omega^2) and CoulombBreit(factor) scales omega, so the correction
    ##     divided by factor^2 must be constant up to an O(factor^2) remainder.
    v0     = InteractionStrength.XL_Breit(1, p32, p32, s2, s2, grid, Basics.CoulombBreit(0.))
    ratios = Float64[]
    for  f  in  [0.25, 0.5, 1.0]
        vf = InteractionStrength.XL_Breit(1, p32, p32, s2, s2, grid, Basics.CoulombBreit(f))
        push!( ratios, (vf - v0) / f^2 )
    end
    spread = (maximum(ratios) - minimum(ratios)) / abs( sum(ratios)/length(ratios) )
    println("  (4) Frequency correction / factor^2 = $(round.(ratios, sigdigits=7)),  relative spread = $spread ")
    if  spread > 1.0e-2
        success = false
        println("  *** The frequency correction is not quadratic in the factor of CoulombBreit(factor), so " *
                "the retardation kernel, its power of omega or its coefficients are wrong.")
    end

    testPrint("testMethod_BreitInteraction()::", success)
    return( success )
end


"""
`TestFrames.testMethod_Opacities(; short::Bool=true)`
    ... asserts that the two mean opacities are the functionals they claim to be. Needs no reference data and
        no cascade: every requirement below is exact, or a textbook constant.

        This is the check that would have caught the defect found on 14-Aug-2026, where
        `simulateRosselandOpacities`, as it was then called, formed  sum_i w_i kappa_i  -- an ARITHMETIC, i.e. Planck, mean --
        while carrying the Rosseland weight and calling the result Rosseland. The Rosseland mean is the
        HARMONIC one, 1/kappa_R = sum_i w_i/kappa_i, and the distinction is the whole physics: a harmonic
        mean is dominated by the most TRANSPARENT frequencies, an arithmetic one by the most opaque, so on a
        line spectrum they differ by orders of magnitude and move in opposite directions. The property had
        never been run, which is why it survived.
        Five things are asserted:
        (1) the Rosseland weights are normalised -- sum_i w_i = 1 to better than 1e-5 on eight nodes;
        (2) the GREY limit: for a constant kappa_nu, BOTH means return that constant exactly. This is what
            fixes the normalisation of each functional independently of the other;
        (3) the DIRECTION: given one opaque bin among transparent ones, the Rosseland mean must lie near the
            TRANSPARENT value and the Planck mean near the opaque one, so that kappa_R << kappa_P. An
            arithmetic-for-harmonic swap makes the two identical and fails here;
        (4) an EMPTY bin sends the Rosseland mean to exactly zero and returns a non-empty explanation, while
            the Planck mean stays finite;
        (5) Thomson scattering off a fully ionised hydrogen plasma gives kappa_es = sigma_T/m_H, the textbook
            0.20 (1 + X) = 0.40 cm^2/g -- a literature value, not an internal one.
        Returns true if all five hold.
"""
function testMethod_Opacities(; short::Bool=true)
    success = true
    printstyled("\n\nTest the mean opacities against their defining limits: \n", color=:light_green)
    printstyled(  "------------------------------------------------------ \n", color=:light_green)

    ## (1) the Rosseland weights must be normalised; this residue is the accuracy floor of everything below.
    ulist, weights = Cascade.rosselandWeights(8)
    devWeights     = abs( sum(weights) - 1.0 )
    println("  (1) Rosseland weights on 8 nodes:  |sum w_i - 1| = $devWeights ")
    if  devWeights > 1.0e-5
        success = false
        println("  *** The Rosseland weight is not normalised; every mean built on it is scaled wrongly.")
    end

    ## (2) GREY limit: a constant kappa_nu must come back unchanged from BOTH means.
    worstGrey = 0.
    for  kappa  in  [0.4, 7.3]
        kappas   = Basics.EmProperty[ Basics.EmProperty(kappa, kappa)  for i = 1:8 ]
        rMean, _ = Cascade.applyOpacityMean(Cascade.RosselandMean(), kappas, weights)
        pMean, _ = Cascade.applyOpacityMean(Cascade.PlanckMean(),    kappas, weights)
        worstGrey = max( worstGrey, abs(rMean.Coulomb - kappa)/kappa, abs(pMean.Coulomb - kappa)/kappa )
    end
    println("  (2) Grey limit, both means against a constant kappa:  worst relative deviation = $worstGrey ")
    if  worstGrey > 1.0e-5
        success = false
        println("  *** A constant spectral opacity is not returned unchanged; the means are mis-normalised.")
    end

    ## (3) DIRECTION: one opaque bin among seven transparent ones must separate the two means widely, with
    ##     the harmonic mean sitting near the TRANSPARENT value. This is what an arithmetic-for-harmonic
    ##     substitution cannot reproduce: it would make the two means equal.
    kappas   = Basics.EmProperty[ Basics.EmProperty(i == 4 ? 100.0 : 0.01, i == 4 ? 100.0 : 0.01)  for i = 1:8 ]
    rMean, _ = Cascade.applyOpacityMean(Cascade.RosselandMean(), kappas, weights)
    pMean, _ = Cascade.applyOpacityMean(Cascade.PlanckMean(),    kappas, weights)
    println("  (3) One opaque bin (100.0) among seven transparent (0.01):  Rosseland = $(rMean.Coulomb), " *
            "Planck = $(pMean.Coulomb) ")
    if  !(rMean.Coulomb < 0.1)   ||   !(pMean.Coulomb > 1.0)   ||   !(pMean.Coulomb > 100*rMean.Coulomb)
        success = false
        println("  *** The two means no longer point in opposite directions. The Rosseland mean must be " *
                "dominated by the TRANSPARENT bins and the Planck mean by the opaque one; if they agree, " *
                "an arithmetic mean is being returned in place of the harmonic one.")
    end

    ## (4) an empty bin annihilates the harmonic mean, and must say so rather than return a usable-looking number.
    kappas[6]      = Basics.EmProperty(0., 0.)
    rMean, rNote   = Cascade.applyOpacityMean(Cascade.RosselandMean(), kappas, weights)
    pMean, _       = Cascade.applyOpacityMean(Cascade.PlanckMean(),    kappas, weights)
    println("  (4) With one EMPTY bin:  Rosseland = $(rMean.Coulomb),  Planck = $(pMean.Coulomb) ")
    if  rMean.Coulomb != 0.  ||  rNote == ""  ||  pMean.Coulomb <= 0.
        success = false
        println("  *** An empty bin must send the Rosseland mean to exactly zero AND return an explanation, " *
                "while leaving the Planck mean finite.")
    end

    ## (5) Thomson scattering of a fully ionised hydrogen plasma: n_e = rho/m_H, hence kappa_es = sigma_T/m_H.
    mHydrogen = 1.6726219e-24;      sigmaThomson = 6.6524587321e-25            ## [g], [cm^2]
    rho       = 1.0e-7;             nElectron    = rho / mHydrogen             ## [g/cm^3], [1/cm^3]
    property  = Cascade.MeanOpacities(Cascade.RosselandMean(),
                                      Cascade.AbstractOpacityContribution[ Cascade.ScatteringOpacity(nElectron) ],
                                      Cascade.TemperatureOpacityDependence(0.05), [nElectron], [rho], [1.0e6], 1., 0.)
    kappas    = Cascade.spectralOpacityContribution(property.contributions[1], Cascade.Data[], property,
                                                    nElectron, rho, 1.0e6, ulist)
    rMean, _  = Cascade.applyOpacityMean(Cascade.RosselandMean(), kappas, weights)
    devThomson = abs( rMean.Coulomb - sigmaThomson/mHydrogen ) / (sigmaThomson/mHydrogen)
    println("  (5) Fully ionised hydrogen, Thomson only:  kappa = $(rMean.Coulomb) cm^2/g  against the " *
            "textbook sigma_T/m_H = $(sigmaThomson/mHydrogen);  relative deviation = $devThomson ")
    if  devThomson > 1.0e-5
        success = false
        println("  *** The scattering opacity does not reproduce kappa_es = 0.20 (1 + X) = 0.40 cm^2/g for " *
                "X = 1, which is the one literature value this module can be held to without a line list.")
    end

    ## (6) Bound-free: hydrogen 1s at threshold. The Gaunt-corrected Kramers cross section is 6.30 Mb, and
    ##     with n_ion = rho/m_H the mass opacity must be exactly sigma/m_H -- which checks the whole unit
    ##     chain a_0^2 -> cm^2 -> cm^2/g, the piece most likely to be wrong by a power of a_0.
    rhoH      = 1.0e-7;    nH = rhoH / mHydrogen
    boundFree = Cascade.BoundFreeOpacity(1.0, [ (Configuration("1s"), Configuration("1s^0")) ])
    property  = Cascade.MeanOpacities(Cascade.RosselandMean(),
                                      Cascade.AbstractOpacityContribution[ boundFree ],
                                      Cascade.TemperatureOpacityDependence(0.05), [nH], [rhoH], [1.0e6], 1., 0.)
    ## u kT = 0.5 Ha exactly at the first node fixes the photon energy on the H 1s threshold.
    kTthr     = 0.5 / ulist[1] * 1.0000001
    Tthr      = kTthr / Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.0)
    kappas    = Cascade.spectralOpacityContribution(boundFree, Cascade.Data[], property, nH, rhoH, Tthr, ulist)
    ## kappa = n sigma / rho, so sigma [cm^2] = kappa rho / n, and 1 Mb = 1.0e-18 cm^2.
    sigmaMb    = kappas[1].Coulomb * rhoH / nH * 1.0e18
    devKramers = abs(sigmaMb - 6.30) / 6.30
    println("  (6) Bound-free, H 1s at threshold:  sigma = $sigmaMb Mb against the Gaunt-corrected " *
            "Kramers value 6.30 Mb;  relative deviation = $devKramers ")
    if  devKramers > 2.0e-3
        success = false
        println("  *** The bound-free contribution does not reproduce the 6.30 Mb threshold cross section " *
                "of hydrogen 1s, so either the cross section or the a_0^2 -> cm^2/g unit chain is wrong.")
    end
    ## and it must vanish BELOW threshold, which is what gives a bound-free opacity its edges
    kappasLow = Cascade.spectralOpacityContribution(boundFree, Cascade.Data[], property, nH, rhoH, 0.2*Tthr, ulist)
    if  kappasLow[1].Coulomb != 0.
        success = false
        println("  *** The bound-free opacity does not vanish below the ionization threshold; it must have " *
                "an edge there.")
    end

    testPrint("testMethod_Opacities()::", success)
    return( success )
end


"""
`TestFrames.testMethod_DocstringPointers(; short::Bool=true)`
    ... asserts that every `Module.name` written in a docstring under src/ actually resolves. Needs no
        reference data: the requirement is exact, so the test is against zero rather than a tabulated number.

        A docstring that names a function which does not exist is worse than one that says nothing, because
        it is read as a promise. Dangling pointers turned up by hand three times in two days in Aug-2026 --
        `XL_CoulombTensor` naming two functions that had never existed, `ScreenedPotentialCache` pointing at
        "computeDirectExchangeVTensor, once written", and `computeCImatrix` calling an `XL_Breit_WO` that was
        never defined, which meant half a reference cross-check could not have run. Each was found by
        accident. This method finds them on purpose.

        WHAT IT DOES NOT FLAG, and why:
        * file names. "see module-Cascade.jl" matches the same pattern and is not a pointer, so a trailing
          `.jl` is skipped;
        * a name whose only fault is a missing `!`. `propagateProbability` does not exist while
          `propagateProbability!` does; that is a real defect but a mechanical one, so it is counted and
          reported SEPARATELY, to keep it from hiding the pointers that need judgement;
        * anything inside a `#= ... =#` block. JAC uses these to disable whole regions of code, and such a
          region carries its own docstrings whose functions are legitimately absent -- the block is off. The
          first version of this method flagged `computeStepAugerAverageSCA`, which sits inside the
          disabled block of module-Cascade-inc-stepwise-decay.jl, and that would have been a false alarm.

        THE CONVENTION FOR HISTORY, which this method quietly enforces: a docstring explaining what a
        function used to be called must write the old name WITHOUT its module prefix -- "simulateRosseland-
        Opacities, as it was then called" -- so that a deliberate historical reference is not written in the
        same form as a live pointer. Prefixing it would make the sentence a lie the moment it is read as an
        instruction, and this method would rightly complain.
        Returns true if no unresolved pointer remains.
"""
function testMethod_DocstringPointers(; short::Bool=true)
    success = true
    printstyled("\n\nTest that every Module.name named in a docstring resolves: \n", color=:light_green)
    printstyled(  "---------------------------------------------------------- \n", color=:light_green)

    docPattern = r"\"\"\"(.*?)\"\"\""s
    ## NOTE THE ABSENCE OF A TRAILING \b, and do not put one back.  Julia's `!` is not a word character, so
    ## `name!\b` never matches: the boundary would have to fall between `!` and whatever follows, which is
    ## usually `(`, and two non-word characters make no boundary.  The pattern then backtracks to the bangless
    ## name and reports every mutating function as unresolved.  The first version of this method did exactly
    ## that and produced TEN false positives -- Cascade.propagateProbability!, pushLevels!, walkPathways! and
    ## the rest -- every one of which is written correctly in its docstring.  A negative lookahead does the
    ## job that \b cannot.
    refPattern = r"\b([A-Z][A-Za-z0-9]*)\.([a-z][A-Za-z0-9_]*!?)(?![A-Za-z0-9_!])"
    dangling   = Dict{String,Vector{String}}();    missingBang = Dict{String,Vector{String}}()
    noTokens   = 0

    for  file  in  filter(x -> endswith(x, ".jl"), readdir(@__DIR__, join=true))
        ## Strip the disabled `#= ... =#` regions first; the docstrings inside them describe code that is
        ## deliberately switched off, so their functions are absent on purpose.
        source = replace( read(file, String), r"#=.*?=#"s => "" )
        for  doc  in  eachmatch(docPattern, source)
            for  ref  in  eachmatch(refPattern, doc.captures[1])
                modName, fnName = ref.captures[1], ref.captures[2];    noTokens = noTokens + 1
                fnName == "jl"                                        &&  continue
                isdefined(JenaAtomicCalculator, Symbol(modName))      ||  continue
                theModule = getfield(JenaAtomicCalculator, Symbol(modName))
                theModule isa Module                                  ||  continue
                isdefined(theModule, Symbol(fnName))                  &&  continue
                if  isdefined(theModule, Symbol(fnName * "!"))
                    push!( get!(missingBang, basename(file), String[]), "$modName.$fnName" )
                else
                    push!( get!(dangling,    basename(file), String[]), "$modName.$fnName" )
                end
            end
        end
    end

    noBang     = sum( length(unique(v)) for v in values(missingBang); init=0 )
    noDangling = sum( length(unique(v)) for v in values(dangling);    init=0 )
    println("  $noTokens Module.name tokens scanned in the docstrings of src/. ")
    println("  Unresolved:  $noDangling dangling,  $noBang missing a `!`. ")
    if  noBang > 0
        success = false
        println("  *** These name a function that exists only WITH a trailing `!`; the name as written does " *
                "not resolve:")
        for  (file, refs)  in  sort(collect(missingBang), by=x->x[1])
            println("        $file:  " * join(sort(unique(refs)), ", "))
        end
    end
    if  noDangling > 0
        success = false
        println("  *** These name a function that does not exist under any name. A docstring pointer is read " *
                "as a promise, so each must be repaired, renamed, or rewritten as prose:")
        for  (file, refs)  in  sort(collect(dangling), by=x->x[1])
            println("        $file:  " * join(sort(unique(refs)), ", "))
        end
    end

    testPrint("testMethod_DocstringPointers()::", success)
    return( success )
end
