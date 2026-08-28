
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
`TestFrames.testModule_WeakInteractionMoment(; short::Bool=true)`
    ... tests on module WeakInteractionMoment. A success::Bool is returned.

        This test asserts SELECTION RULES and STRUCTURE, not stored numbers, and deliberately so. Its predecessor compared against an
        approved file whose every entry had im = 2*re exactly -- the fingerprint of the stub `(SUM coeff.T) * (1.0 + 2.0im)` that the
        module used to call -- and in which the weak-charge and Schiff-moment amplitudes were bit-identical to each other. Numbers of
        that kind can be re-blessed indefinitely without ever becoming right. What is checked here instead cannot be satisfied by any
        stub and needs no external reference:

          (i)   all three operators vanish EXACTLY between levels of the same parity, since all three are P-odd;
          (ii)  the rank-0 weak charge vanishes EXACTLY unless J_f = J_i, as a pseudoscalar must;
          (iii) the weak-charge amplitude is purely imaginary and the Schiff-moment amplitude purely real, which follows from gamma_5
                being antisymmetric in the large and small components while the Schiff operator is not;
          (iv)  H_W is Hermitian, so exchanging the two levels conjugates the amplitude -- which, for a purely imaginary rank-0
                quantity, means exact antisymmetry;
          (v)   the weak-charge and Schiff amplitudes are no longer equal to one another.
"""
function testModule_WeakInteractionMoment(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-WeakInteractionMoment-new.sum")
    printstyled("\n\nTest the module  WeakInteractionMoment  ... \n", color=:cyan)
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    nModel = Nuclear.Model(55.)
    asf    = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    mp     = SelfConsistent.performSCF([Configuration("1s"), Configuration("2s"), Configuration("2p")], nModel, grid, asf; printout=false)
    success = true;    tol = 1.0e-14

    for  a  in  mp.levels,  b  in  mp.levels
        wc = WeakInteractionMoment.weakChargeAmplitude(a, b, nModel, grid)
        sm = WeakInteractionMoment.schiffMomentAmplitude(a, b, nModel, grid)
        an = WeakInteractionMoment.anapoleAmplitude(a, b, nModel, grid)
        # (i) P-odd operators vanish between equal parities
        if  a.parity == b.parity  &&  (wc != 0.  ||  sm != 0.  ||  an != 0.)
            success = false;    println(">> same-parity pair $(a.index)<-$(b.index) does not vanish: $wc $sm $an")
        end
        # (ii) the rank-0 weak charge needs J_f = J_i
        if  a.J != b.J  &&  wc != 0.
            success = false;    println(">> rank-0 weak charge nonzero for J_f != J_i at $(a.index)<-$(b.index): $wc")
        end
        # (iii) the two operators have opposite reality
        if  abs(wc) > 0.  &&  abs(real(wc)) > tol*abs(wc)
            success = false;    println(">> weak-charge amplitude not purely imaginary at $(a.index)<-$(b.index): $wc")
        end
        if  abs(sm) > 0.  &&  abs(imag(sm)) > tol*abs(sm)
            success = false;    println(">> Schiff amplitude not purely real at $(a.index)<-$(b.index): $sm")
        end
        # (iv-bis) the anapole is purely imaginary too, and Hermitian
        if  abs(an) > 0.
            if  abs(real(an)) > tol*abs(an)
                success = false;    println(">> anapole amplitude not purely imaginary at $(a.index)<-$(b.index): $an")
            end
            anR = WeakInteractionMoment.anapoleAmplitude(b, a, nModel, grid)
            # Hermiticity: antisymmetric for J_f = J_i, symmetric otherwise -- both are exact
            dev = a.J == b.J ? abs(an + anR) : abs(an - anR)
            if  dev > tol*abs(an)
                success = false;    println(">> anapole fails Hermiticity under exchange at $(a.index)<-$(b.index)")
            end
        end
        # (iv) Hermiticity, and (v) the two operators must differ
        if  abs(wc) > 0.
            if  abs(wc + WeakInteractionMoment.weakChargeAmplitude(b, a, nModel, grid)) > tol*abs(wc)
                success = false;    println(">> weak charge not antisymmetric under exchange at $(a.index)<-$(b.index)")
            end
            if  abs(wc - sm) <= tol*abs(wc)
                success = false;    println(">> weak-charge and Schiff amplitudes coincide at $(a.index)<-$(b.index)")
            end
        end
    end
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_WeakInteractionMoment()::", success)
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
    # NO grid is given: the box is derived from refConfigs by Basics.recommendedGrid, which is the point.  Until
    # 17-Aug-2026 this ran on Radial.Grid(true) -- 614 a.u. and 62 splines for a FOUR-ELECTRON atom -- and on
    # that grid a correlation layer RAISES the energy, which is not physics.  The reference box is 20.9 a.u.
    # with 96 splines; the uncorrelated first layer agrees to eight digits between the two, so nothing about
    # the reference SCF is traded away.
    wa          = Representation(name, Nuclear.Model(4.), refConfigs,
                                 RasExpansion([LevelSymmetry(0, Basics.plus)], 4, coreShells, fromShells, layers, rasSettings) )
    wb = generate(wa, output=true)
    # The expected value changed TWICE on 16-Aug-2026. First when Basics.EOLField was wired to the
    # orbital-rotation solver, and again when computeOrbitalGradient was corrected: it contracted a screened
    # potential built from the sign-canonicalized ORBITALS with a RAW b-vector, so an off-diagonal CSF pair
    # carried an odd power of a correlating orbital's sign and the gradient drifted from the functional --
    # 2.5x, then 19x, then SIGN-WRONG as the 2p weight grew, which stalled the line search at iteration 3 of
    # 24.  With the full s/wN factor on every slot the analytic slope matches a finite difference to five
    # digits, the driver runs to its iteration limit, and this case reached -14.617216 instead of -14.612567.
    #   THIRD change, same day: the search direction became conjugate (Polak-Ribiere+).  Plain preconditioned
    # steepest descent zigzagged -- |grad| stayed flat near 0.057 over 40 iterations while the energy crept
    # down -- and conjugacy drives it to 0.013 for the same cost, reaching -14.619313 here.  That value was a
    # converging upper bound rather than a converged number, and moved with the iteration limit; see the FIFTH
    # change below, after which it no longer does.  The former -14.589901130961 was the COLLAPSED answer: with the old solver this very case
    # converged cleanly, to its own 1e-8, onto a degenerate stationary point whose mixing vector was
    # [0.971, -0.0000482, 0.238] -- the 2p_3/2 correlation channel eliminated outright.  An independent
    # DFS-Field run of the identical layer structure gave -14.605300 Ha with both channels contributing.
    # The rotation solver reaches -14.612567, i.e. BELOW that reference, so this test previously asserted
    # the defect.  Do not restore the old number.
    # 17-Aug-2026, the FOURTH change and the reason for the two above: the box now comes from refConfigs and the
    # EOL default became L-BFGS.  Those are connected.  L-BFGS had been rejected for "losing to conjugate
    # gradients" in every earlier RAS test, and every earlier RAS test ran on the 614 a.u. grid, where it does
    # lose -- by 1.46 mHa.  On the reference-sized box it WINS by 0.26 mHa and runs three times faster (8 s
    # against 23 s), and with the iteration budget varied it wins at every budget: 24 iterations reach what
    # conjugacy needs about 70 to reach.  A quasi-Newton method builds its curvature model from gradient
    # DIFFERENCES, so it was punished twice over, first by the gradient defect fixed the day before and then by
    # a basis too starved to be worth modelling.  Benchmarking an optimiser on a basis that cannot represent
    # the problem measures the basis.
    #
    # 27-Aug-2026, the FIFTH change, and the first that alters WHICH PROBLEM is being solved rather than how
    # well it is solved: solveOptimizedLevelFieldByRotation now honours AsfSettings.frozenSubshells, which it
    # had ignored entirely.  Every layer of this expansion sets it -- step 1 freezes 1s, step 2 freezes 1s and
    # 2s -- so until now no layer froze anything and each re-optimized every orbital, which is not what a RAS
    # layer means.  The energy therefore RISES, as a constrained minimum must: -14.619380954693 -> -14.618871269,
    # i.e. +5.10e-04 Ha, while |grad| at the iteration limit falls from 2.99e-02 to 4.45e-05, a factor of 670.
    #   THE NEW VALUE IS CONVERGED, which the old one was not, and that is the reason for preferring it rather
    # than mere novelty.  Run at 24, 60 and 120 iterations it gives -14.618871268972836, -14.618871275673818
    # and -14.618871275673818 -- identical at the last two, so the residual iteration-limit dependence is
    # 6.7e-09 Ha, five orders below the shift itself.  Step 1 is unmoved across all three budgets and now exits
    # on a stationary energy at iteration 15 instead of running out at 24.  The tolerance below is left at
    # 1.0e-3, which is now very loose for a number stable to 6.7e-09; tightening it is a separate decision.
    #
    # Do not restore any earlier number: each was obtained on a grid, a gradient or an unconstrained layer that
    # no longer exists here.
    if  abs(wb["step2"].levels[1].energy + 14.618871268972836)  > 1.0e-3
        success = false
        if printTest   info(iostream, "levels[1].energy $(wb["step2"].levels[1].energy) != -14.618871268972836")   end
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
        ("PhotoDoubleIonization.Settings()",        () -> PhotoDoubleIonization.Settings()             ),
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
`TestFrames.testMethod_SpinAngular(; short::Bool=true)`
    ... tests `SpinAngular` against statements that hold independently of any other implementation: exact identities
        that a correct spin-angular module must satisfy, closed forms whose values are known analytically, and
        selection rules whose violation is a defect rather than a disagreement. A success::Bool is returned.

        WHY IT DOES NOT COMPARE AGAINST `SpinAngularGaigalas`. That comparison exists, in
        `tools/diag-spinangular-compare.jl`, and it is valuable -- it found three real defects between 27 and
        28-Aug-2026. But it is a COMPARISON and not a REFERENCE: it answers "do the two agree", which says something
        about this module only while the other one is present and right. It is also blind to any error the two
        share, and the parity fault repaired on 27-Aug was exactly that, an assumption their CALLERS held in common
        for years. This test therefore asserts things that are true of the physics, so that it keeps its meaning if
        `SpinAngularGaigalas` is ever removed.

        THE SEVEN CHECKS, and each is chosen because a known defect would have broken it:
        1. the occupation identity -- a diagonal rank-0 coefficient IS the occupation number, for any coupling;
        2. the single-creation closed form, <j^0||a||j^1> = sqrt(2j+1), UP TO j = 21/2, which the quasispin tables
           cannot reach and which returned 0 there until 28-Aug (priority item 60);
        3. the single-subshell W closed form, <j^1||W^(k)||j^1> = -sqrt(2k+1), independent of j, same reach;
        4. a one-electron substitution out of a closed subshell, which must give sqrt(2 j_a + 1) -- silently dropped
           for l >= 1 until 28-Aug (item 59);
        5. the parity and angular-momentum selection rules, which a scalar operator cannot violate;
        6. hermiticity of the rank-0 one-particle coefficients under exchanging the two CSFs;
        7. a GRASP2018-traceable absolute value: for closed shells plus ONE electron every rank > 0 coefficient is
           exactly 1, which is what GRASP returns for that case.
"""
function testMethod_SpinAngular(; short::Bool=true)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    csfsOf = function(confs)
        relconfs = ConfigurationR[]
        for  c in confs
            append!(relconfs, Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(c)))
        end
        subshells = Basics.generateSubshellList(relconfs);    csfs = CsfR[]
        for  rc in relconfs    append!(csfs, Basics.generateCsfRs(rc, subshells))    end
        return( (csfs, subshells) )
    end

    # (1) THE OCCUPATION IDENTITY. <Psi| sum_i f(i) |Psi> = sum_a N_a <a|f|a> holds for ANY coupling, so a diagonal
    #     rank-0 coefficient must be exactly its subshell's occupation number. This is an identity, not a tolerance.
    for  confs in ( ["1s^2 2s^2", "1s^2 2p^2"], ["1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"], ["1s^2 2p^6 3d^1 4d^1"] )
        csfs, subshells = csfsOf(confs)
        for  csf in csfs
            coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                           csf, csf, subshells)
            dev    = SpinAngular.checkOccupationIdentity(coeffs, csf, subshells)
            if  dev > 1.0e-12
                success = false
                if printTest   info(iostream, "SpinAngular: the occupation identity fails by $dev for $confs")   end
            end
        end
    end

    # (2) and (3) THE CLOSED FORMS, TAKEN ABOVE j = 9/2 ON PURPOSE. The quasispin tables reach 9/2 and no further;
    #     a continuum partial wave reaches 21/2. Both routines returned 0 up there until 28-Aug-2026, which cost the
    #     electron-impact excitation cross section 27 %.
    for  j2 = 1:2:21
        j  = AngularJ64(j2//2)
        wa = SpinAngular.shellReducedA(j, 0, 0, AngularJ64(0), 1, 1, j, AngularM64(-1//2))
        wb = SpinAngular.shellReducedA(j, 1, 1, j, 0, 0, AngularJ64(0), AngularM64(1//2))
        if  abs(wa - sqrt(j2 + 1.0)) > 1.0e-12  ||  abs(wb + sqrt(j2 + 1.0)) > 1.0e-12
            success = false
            if printTest   info(iostream, "SpinAngular: <j^0||a||j^1> at j = $j gives $wa and $wb, not " *
                                          "+-$(sqrt(j2+1.0))")   end
        end
        for  kj = 0:min(4, j2)
            wc = SpinAngular.shellReducedW(j, 1, 1, j, 1, j, kj)
            if  abs(wc + sqrt(2kj + 1.0)) > 1.0e-12
                success = false
                if printTest   info(iostream, "SpinAngular: <j^1||W^($kj)||j^1> at j = $j gives $wc, not " *
                                              "$(-sqrt(2kj+1.0))")   end
            end
        end
    end

    # (4) A ONE-ELECTRON SUBSTITUTION OUT OF A CLOSED SUBSHELL into an empty one of the same kappa. The coefficient
    #     is sqrt(2 j_a + 1). Dropped silently for l >= 1 until 28-Aug-2026, and invisible for l = 0.
    let  (csfs, subshells) = csfsOf(["1s^2 2s^2 2p^6", "1s^2 2s^2 2p^5 3p^1"]);   seen = 0
        for  r in csfs,  s in csfs
            if  r.J != s.J  ||  r.parity != s.parity    continue    end
            for  coeff in SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                          r, s, subshells)
                if  coeff.a == coeff.b  ||  abs(coeff.T) < 1.0e-12    continue    end
                seen = seen + 1
                if  abs(abs(coeff.T) - sqrt(Basics.subshell_2j(coeff.a) + 1.0)) > 1.0e-12
                    success = false
                    if printTest   info(iostream, "SpinAngular: substitution $(coeff.a) x $(coeff.b) gives " *
                                                  "$(coeff.T), not +-sqrt(2j+1)")   end
                end
            end
        end
        if  seen == 0
            success = false
            if printTest   info(iostream, "SpinAngular: no substitution coefficient found at all for " *
                                          "2p^6 <-> 2p^5 3p; they are being dropped")   end
        end
    end

    # (5) THE SELECTION RULES. A scalar operator connects neither different J nor different parity, and no one-body
    #     operator can change a subshell occupation by two.
    let  (csfs, subshells) = csfsOf(["1s^2 2s^2 2p^2", "1s^2 2s^1 2p^3"])
        for  r in csfs,  s in csfs
            if  r.J == s.J  &&  r.parity == s.parity    continue    end
            if  length(SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                       r, s, subshells)) != 0
                success = false
                if printTest   info(iostream, "SpinAngular: a scalar operator connects two CSFs of different " *
                                              "symmetry, which it cannot")   end
            end
        end
    end

    # (6) HERMITICITY. Exchanging the two CSFs must exchange the two subshells of every coefficient and leave the
    #     value alone. This is what caught the 4th appearance of the X-coupling defect on 24-Aug-2026.
    let  (csfs, subshells) = csfsOf(["1s^2 2s^2", "1s^2 2p^2"])
        for  r in csfs,  s in csfs
            if  r.J != s.J  ||  r.parity != s.parity    continue    end
            fw = Dict( (c.a, c.b) => c.T  for c in
                       SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(0, Basics.plus), r, s, subshells) )
            bw = Dict( (c.b, c.a) => c.T  for c in
                       SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(0, Basics.plus), s, r, subshells) )
            for  key in union(keys(fw), keys(bw))
                if  abs( get(fw, key, 0.0) - get(bw, key, 0.0) ) > 1.0e-12
                    success = false
                    if printTest   info(iostream, "SpinAngular: hermiticity fails at $key, " *
                                                  "$(get(fw,key,0.0)) against $(get(bw,key,0.0))")   end
                end
            end
        end
    end

    # (7) A VALUE TRACEABLE TO GRASP2018. For closed shells plus ONE electron every rank > 0 coefficient GRASP
    #     returns is exactly 1. Note WHY this is only a weak anchor and is not relied on alone: an implementation
    #     that returned 1 for every allowed pair would score perfectly here, which is the trap recorded on
    #     23-Aug-2026 when this very case was used as the whole test set.
    let  (csfs, subshells) = csfsOf(["1s^2 2s^1", "1s^2 2p^1"]);   seen = 0
        for  r in csfs,  s in csfs,  k = 1:2
            for  coeff in SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(k, Basics.plus),
                                                          r, s, subshells)
                if  abs(coeff.T) < 1.0e-12    continue    end
                seen = seen + 1
                if  abs(abs(coeff.T) - 1.0) > 1.0e-12
                    success = false
                    if printTest   info(iostream, "SpinAngular: rank-$k coefficient $(coeff.T) for a " *
                                                  "closed-shell-plus-one case, where GRASP2018 gives 1")   end
                end
            end
        end
        if  seen == 0
            success = false
            if printTest   info(iostream, "SpinAngular: no rank > 0 coefficient found for 1s^2 2s + 1s^2 2p")   end
        end
    end

    testPrint("testMethod_SpinAngular()::", success)
    return(success)
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
