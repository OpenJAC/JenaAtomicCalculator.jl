
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
    # 60 iterations, not 24, because this test now asserts CONVERGENCE and not merely a number: with the
    # gradient corrected (items 121, 122 and the L-BFGS curvature pairs) step 1 converges at 12 and step 2 at
    # 45, so a budget of 24 would assert a failure.  The cost is a few seconds.
    rasSettings = RasSettings([1], 60, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1]) )
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
    # 31-Aug-2026, the SIXTH change, and the same KIND of change as the fifth: solveAverageLevelField now
    # honours AsfSettings.frozenSubshells, which it too had ignored entirely -- only the DFS mean-field driver
    # observed the field.  That mattered here because performSCF runs an AL pass as the STARTING POINT of every
    # EOL computation (the :optimizedLevel branch), and the EOL solver pins its frozen orbitals to what it finds
    # ON ENTRY.  So although Basics.generate hands each RAS step the previous step's converged orbitals for the
    # shells it freezes, that AL pass then re-optimized them anyway -- on THIS layer's basis, whose flat CSF
    # average includes the 2p correlation configurations.  Step 2's "frozen" 1s and 2s were therefore neither
    # step 1's orbitals nor free ones, and the minimum was less constrained than a RAS layer means.
    #   The energy RISES accordingly, as a constrained minimum must: -14.618871268972836 -> -14.614058864452650,
    # i.e. +4.81e-03 Ha.  It is a CONVERGED number and a better-behaved one: both layers now reach the gradient
    # tolerance outright, step 1 at iteration 12 and step 2 at 44, with |grad| = 8.5e-07 and 8.7e-07, where the
    # fifth change left step 2 at 4.45e-05 on the iteration limit.  Run at 24, 60 and 120 iterations it gives
    # -14.614058793496298, -14.614058864452650 and -14.614058864452650 -- identical at the last two, so the
    # residual budget dependence is 7.1e-08 Ha, five orders below the shift itself; step 1 is bit-identical at
    # all three.
    #
    # Do not restore any earlier number: each was obtained on a grid, a gradient or an unconstrained layer that
    # no longer exists here.
    if  abs(wb["step2"].levels[1].energy + 14.614058864452650)  > 1.0e-3
        success = false
        if printTest   info(iostream, "levels[1].energy $(wb["step2"].levels[1].energy) != -14.614058864452650")   end
    end

    # AND THE FIELD MUST CONVERGE, WHICH IS A DIFFERENT ASSERTION FROM THE ENERGY ABOVE.  Every EOL exit
    # other than convergence now raises a collected warning, so a converged run leaves the report clean and a
    # stalled one does not.  Reading the report is also what a user has: `generate` ends with PrintWarnings()
    # followed by ResetWarnings(), so by the time it returns the file is the only record.
    #   THIS IS THE CHECK THAT WAS MISSING.  The energy assertion above passed unchanged through every change
    # made to this solver between 29 and 31-Aug-2026, INCLUDING the two that made it worse -- a unit-step
    # reset that killed the line search at iteration 11, and a 4x growth factor that lost a convergence -- and
    # through the three gradient defects that preceded them.  A tolerance of 1e-3 on an energy cannot see any
    # of that; the stop reason can.
    #   NEGATIVE CONTROL, measured 31-Aug-2026 rather than assumed: forcing the field to :steepest leaves both
    # layers unconverged at the 60-iteration limit, step 2 with |grad| = 5.3e-03 -- and its energy is
    # -14.618760171, which deviates by 1.1e-04 and so PASSES the energy assertion above.  The check below is
    # what fails, which is the whole reason it exists.
    let  report = "jac-warn.report"
        text = isfile(report) ? read(report, String) : ""
        if  occursin("solveOptimizedLevelFieldByRotation", text)
            success = false
            if printTest   info(iostream, "the EOL field did not converge in this RAS expansion; " *
                                          "jac-warn.report holds: $text")   end
        end
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
    ## Nothing is excluded. Liouville stood here from 09-Aug-2026 as a "known failure ... its copy-constructor
    ## does not convert", which was not quite what was wrong with it: the module had NO copy-constructor at all,
    ## so the hasmethod guard below would have passed over it in any case and the exclusion decided nothing. One
    ## was written on 28-Aug-2026 and it converts, so the entry is removed rather than corrected. Postponing a
    ## module is a decision about its PHYSICS; its constructor is ordinary code and is tested like any other.
    excluded = Symbol[]
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

    ## (5) The two QUADRATURES of the same integral must agree.  CoulombBreit(factor, :direct) sums the full
    ##     double loop over the grid; CoulombBreit(factor, :swept) uses the fact that every Breit kernel
    ##     factorises into a function of r_< times one of r_>, and accumulates the inner integral in one pass.
    ##     They differ only in the ORDER of the summation, so they cannot be bitwise equal -- and must still
    ##     agree to the accuracy of the quadrature.  Both frequency regimes are covered: factor = 0. takes the
    ##     analytic omega -> 0 kernels, factor = 1. the Bessel ones, and those are DIFFERENT code in both routes.
    ##     KEEP THIS TEST.  An unused route rots: XL_BreitDamped sat with an `error("stop a")` body behind a
    ##     signature one argument short of its only caller, unnoticed, because no example ever selected it.
    p12   = orbs[Subshell("2p_1/2")];       s1 = orbs[Subshell("1s_1/2")]
    quads = [ (1, p32, p32, s2, s2), (1, p32, p32, p32, p32), (1, p12, p32, s2, s2),
              (2, p32, p32, p32, p32), (1, s1, p32, s2, p32),  (2, p32, s2, p32, s2) ]
    worst = 0.;     nUsed = 0
    for  f  in  [0., 1.]
        for  q  in  quads
            vd = InteractionStrength.XL_Breit(q..., grid, Basics.CoulombBreit(f, :direct))
            vs = InteractionStrength.XL_Breit(q..., grid, Basics.CoulombBreit(f, :swept))
            abs(vd) < 1.0e-10   &&   continue
            worst = max(worst, abs(vs - vd) / abs(vd));     nUsed = nUsed + 1
        end
    end
    println("  (5) :swept against :direct on $nUsed non-vanishing strengths, both frequencies:  " *
            "worst relative deviation = $worst ")
    if  nUsed < 6
        success = false
        println("  *** Too few strengths survived the 1e-10 floor for this comparison to test anything.")
    end
    if  worst > 1.0e-8
        success = false
        println("  *** The swept and the direct Breit quadrature disagree by more than the quadrature's own " *
                "accuracy. They are the same integral summed in a different order, so this is not " *
                "reassociation: one of the separable kernels in XL_Breit_densitiesSwept is wrong.")
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
`TestFrames.testMethod_DensityAtNucleus(; short::Bool=true)`
    ... tests Basics.densityAtNucleus, the electron density of one orbital averaged over the nuclear volume;
        a success::Bool is returned.  Every check below is an EXACT statement that needs no reference data.
"""
function testMethod_DensityAtNucleus(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DensityAtNucleus-new.sum")
    printstyled("\n\nTest the method  Basics.densityAtNucleus()  ... \n", color=:cyan)
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")

    ## The hydrogen-like 1s density at the nucleus is Z^3/pi exactly in the non-relativistic point-nucleus
    ## limit, so a light H-like ion is an absolute anchor: no table, no fitted number.  Relativistic and
    ## finite-size corrections are of order (Z alpha)^2 and grow with Z.  Measured 04-Sep-2026: the ratio to
    ## Z^3/pi is 1.0005, 1.0019 and 1.0042 at Z = 1, 2, 3 -- a MONOTONE rise, which is the relativistic
    ## enhancement rather than numerical noise, so 1 % is a tolerance the method comfortably meets and a
    ## failure would mean something real.
    d = Dict{Float64,Float64}()
    for  Z  in (1.0, 2.0, 3.0)
        grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6/Z, h = 5.0e-2, hp = 2.0e-2, rbox = 30.0/Z)
        nm   = Nuclear.Model(Z, Nuclear.UniformNucleus())
        asf  = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
        wa   = Atomic.Computation(Atomic.Computation(); name = "H-like", grid = grid, nuclearModel = nm,
                                  configs = [Configuration("1s")], asfSettings = asf)
        bas  = redirect_stdout(devnull) do;  perform(wa; output = true)["multiplet:"].levels[1].basis  end
        d[Z] = Basics.densityAtNucleus(bas.orbitals[Subshell("1s_1/2")], grid, nm)
        exact = Z^3/pi
        if  abs(d[Z]/exact - 1.0) > 0.01
            success = false
            if printTest   info(iostream, "the H-like 1s density at Z = $Z is $(d[Z]) against the exact " *
                                          "Z^3/pi = $exact, a deviation of $(abs(d[Z]/exact - 1))")   end
        end
    end
    ## and the Z^3 scaling itself, which divides the absolute normalisation out
    for  (Z, factor)  in  ((2.0, 8.0), (3.0, 27.0))
        if  abs( d[Z]/d[1.0]/factor - 1.0 ) > 0.01
            success = false
            if printTest   info(iostream, "the density does not scale as Z^3: d($Z)/d(1) = $(d[Z]/d[1.0]) " *
                                          "against $factor")   end
        end
    end

    ## A p_1/2 orbital has a NON-ZERO density at the nucleus, carried entirely by the small component -- a
    ## purely relativistic effect, and the reason p_1/2 orbitals take part in electron capture at all.  It
    ## must be non-zero and it must be far smaller than the s density of the same shell.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(54.0, Nuclear.UniformNucleus())
    asf  = AsfSettings(AsfSettings(); scField = Basics.DFSField())
    wa   = Atomic.Computation(Atomic.Computation(); name = "Xe", grid = grid, nuclearModel = nm,
               configs = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6")],
               asfSettings = asf)
    bas  = redirect_stdout(devnull) do;  perform(wa; output = true)["multiplet:"].levels[1].basis  end
    dS   = Basics.densityAtNucleus(bas.orbitals[Subshell("3s_1/2")],  grid, nm)
    dP   = Basics.densityAtNucleus(bas.orbitals[Subshell("3p_1/2")],  grid, nm)
    if  !(0. < dP < 0.5*dS)
        success = false
        if printTest   info(iostream, "the 3p_1/2 density at the nucleus is $dP against $dS for 3s; it must " *
                                      "be positive (a relativistic small-component effect) and much smaller")  end
    end

    ## The guard: a POINT nucleus has no volume to average over and must refuse rather than return a number.
    raised = false
    try     Basics.densityAtNucleus(bas.orbitals[Subshell("3s_1/2")], grid, Nuclear.Model(54.0, PointNucleus()))
    catch
        raised = true          ## NOTE: the assignment must NOT sit on the `catch` line -- Julia would read
    end                        ##       `raised` as the name of the exception variable and then fail on `=`.
    if  !raised
        success = false
        if printTest   info(iostream, "a point nucleus did not raise, although it has no volume to average over")  end
    end

    println(iostream, "Basics.densityAtNucleus: the H-like 1s density against the exact Z^3/pi at Z = 1, 2, 3 " *
                      "and its Z^3 scaling; a non-zero but much smaller p_1/2 density, which exists only " *
                      "relativistically; and the refusal for a point nucleus. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testMethod_DensityAtNucleus()::", success)
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


"""
`TestFrames.testModule_AngularMomentum(; short::Bool=true)`  ... tests on module AngularMomentum; a success::Bool is returned.
"""
function testModule_AngularMomentum(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-AngularMomentum-new.sum")
    printstyled("\n\nTest the module  AngularMomentum  ... \n", color=:cyan)
    # testMethod_Wigner_3j already checks three TABULATED values. This test checks IDENTITIES instead, and that is
    # a different thing: a tabulated value pins one point, whereas an orthogonality sum or a closed-form special
    # value constrains a whole family at once and cannot be satisfied by a lookup that happens to hold three
    # entries. Every check below is exact algebra, so the tolerances are at machine level rather than physical.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    jj(k2) = iseven(k2) ? AngularJ64(k2 ÷ 2) : AngularJ64(k2 // 2)      # from DOUBLED angular momenta, so that
    mm(k2) = iseven(k2) ? AngularM64(k2 ÷ 2) : AngularM64(k2 // 2)      # integer and half-integer share one loop

    # (1) 3j orthogonality:  sum_{ma,mb} (2c+1) * 3j(a,b,c; ma,mb,-M)^2 = 1  for every allowed c and M.
    for  (a2, b2)  in  [(2,2), (3,2), (4,6), (5,3), (1,1)]
        for  c2 = abs(a2-b2):2:(a2+b2)
            for  mc2 = -c2:2:c2
                wa = 0.
                for  ma2 = -a2:2:a2,  mb2 = -b2:2:b2
                    if  ma2 + mb2 != mc2    continue    end
                    wa = wa + (c2+1) * AngularMomentum.Wigner_3j(jj(a2), jj(b2), jj(c2), mm(ma2), mm(mb2), mm(-mc2))^2
                end
                if  abs(wa - 1.0) > 1.0e-10
                    success = false
                    if printTest   info(iostream, "3j orthogonality for 2a,2b,2c,2M = $a2,$b2,$c2,$mc2 gives $wa, must be 1")   end
                end
            end
        end
    end

    # (2) the closed-form special value  3j(j,0,j; m,0,-m) = (-1)^(j+m) / sqrt(2j+1).
    #     The phase is (j+m) and NOT (j-m). Both give the same answer for integer j, since 2m is then even, so the
    #     wrong one passes every integer case and fails every half-integer one -- which is what it did when this
    #     test was first written. The phase that holds for both follows from the Clebsch-Gordan relation
    #     (j1 j2 j3; m1 m2 m3) = (-1)^(j1-j2-m3)/sqrt(2j3+1) * <j1 m1 j2 m2 | j3 -m3>, which at j2 = 0 gives
    #     (-1)^(j+m); e.g. 3j(1/2, 0, 1/2; 1/2, 0, -1/2) = -1/sqrt(2), not +1/sqrt(2). Half-integer j is in the
    #     loop for exactly this reason.
    for  j2 = 0:1:8
        for  m2 = -j2:2:j2
            wa = AngularMomentum.Wigner_3j(jj(j2), jj(0), jj(j2), mm(m2), mm(0), mm(-m2))
            wb = (-1.)^((j2+m2)/2) / sqrt(j2 + 1.)
            if  abs(wa - wb) > 1.0e-12
                success = false
                if printTest   info(iostream, "3j(j,0,j; m,0,-m) for 2j,2m = $j2,$m2 gives $wa, must be $wb")   end
            end
        end
    end

    # (3) the closed-form 6j special value  {a b c; 0 c b} = (-1)^(a+b+c) / sqrt( (2b+1)(2c+1) )
    for  (a2, b2)  in  [(2,4), (3,3), (4,4), (6,2), (1,5)]
        for  c2 = abs(a2-b2):2:(a2+b2)
            wa = AngularMomentum.Wigner_6j(jj(a2), jj(b2), jj(c2), jj(0), jj(c2), jj(b2))
            wb = (-1.)^((a2+b2+c2)/2) / sqrt((b2+1.) * (c2+1.))
            if  abs(wa - wb) > 1.0e-12
                success = false
                if printTest   info(iostream, "6j{a b c; 0 c b} for 2a,2b,2c = $a2,$b2,$c2 gives $wa, must be $wb")   end
            end
        end
    end

    # (4) a violated triangle must give EXACTLY zero, not something small. Only INTEGER j appear here: m = 0 is
    #     not a legal projection of a half-integer j, and WignerSymbols rejects the pair outright rather than
    #     returning the zero the test is looking for.
    for  (a2, b2, c2)  in  [(2,2,10), (2,4,14), (4,2,14), (6,6,26)]
        wa = AngularMomentum.Wigner_3j(jj(a2), jj(b2), jj(c2), mm(0), mm(0), mm(0))
        if  wa != 0.
            success = false
            if printTest   info(iostream, "3j with a violated triangle 2a,2b,2c = $a2,$b2,$c2 gives $wa, must be exactly 0")   end
        end
    end

    println(iostream, "AngularMomentum: 3j orthogonality, the (j,0,j) and {a b c; 0 c b} closed forms, and exact " *
                      "zeros for violated triangles. Identities rather than tabulated points; no approved data.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_AngularMomentum()::", success)
    return( success )
end


"""
`TestFrames.testModule_HydrogenicIon(; short::Bool=true)`  ... tests on module HydrogenicIon; a success::Bool is returned.
"""
function testModule_HydrogenicIon(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-HydrogenicIon-new.sum")
    printstyled("\n\nTest the module  HydrogenicIon  ... \n", color=:cyan)
    # The Dirac energy is checked by its STRUCTURE rather than against a re-typed copy of the same formula, which
    # would be circular: the l-degeneracy at equal j, the non-relativistic limit, the (alpha Z)^2 scaling of the
    # departure from it, and the ordering of the fine-structure pair. A wrong formula fails at least one of those.
    # The r^k expectation values are then checked against a SECOND ROUTE -- numerical quadrature over the module's
    # own radial orbital -- so the closed form and the tabulated orbital must agree with each other.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    alpha   = Defaults.getDefaults("alpha")
    en(sh, Z) = redirect_stdout(devnull) do;  HydrogenicIon.energy(sh, Z)  end

    # (1) Dirac l-degeneracy: 2s_1/2 and 2p_1/2 share n and j, so their energies must agree to machine precision
    e2s = en(Subshell("2s_1/2"), 10.);    e2p = en(Subshell("2p_1/2"), 10.)
    if  abs(e2s - e2p) > 1.0e-12 * abs(e2s)
        success = false
        if printTest   info(iostream, "2s_1/2 and 2p_1/2 must be degenerate; got $e2s and $e2p")   end
    end
    # (2) the fine-structure pair is ordered: 2p_3/2 lies ABOVE 2p_1/2
    e2p3 = en(Subshell("2p_3/2"), 10.)
    if  !(e2p3 > e2p)
        success = false
        if printTest   info(iostream, "2p_3/2 = $e2p3 must lie above 2p_1/2 = $e2p")   end
    end
    # (3) the non-relativistic limit, and the (alpha Z)^2 law for the departure from it. The RATIO of the relative
    #     departures at Z = 2 and Z = 1 must be 4; that is implementation-independent, unlike the value itself.
    rel(Z) = (en(Subshell("1s_1/2"), Z) - (-Z^2/2)) / (-Z^2/2)
    r1 = rel(1.);   r2 = rel(2.)
    if  abs(r1) > 1.0e-4
        success = false
        if printTest   info(iostream, "1s_1/2 at Z=1 departs from -Z^2/2 by $r1, far more than (alpha Z)^2")   end
    end
    if  abs(r2/r1 - 4.0) > 1.0e-2
        success = false
        if printTest   info(iostream, "the relativistic shift scales as $(r2/r1) between Z=1 and Z=2, must be 4")   end
    end
    # (4) rkExpectation against numerical quadrature over the module's own radial orbital
    # rbox = 60 rather than 40, by Rule 12: hydrogenic 3s has its outer turning point at r+ = 18 a.u. and wants
    # a box near 45. At 40 the truncation shows up directly as <3s|3s> = 0.9999896 instead of 1.
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 60.0)
    for  sh  in  [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3d")]
        P = HydrogenicIon.radialOrbital(sh, 1., grid)
        for  (srk, k)  in  [("r", 1), ("r^2", 2), ("1/r", -1)]
            wa = HydrogenicIon.rkExpectation(srk, sh, 1.)
            wb = RadialIntegrals.rkDiagonal(k, P, P, grid)
            if  abs(wa - wb) > 1.0e-3 * abs(wa)
                success = false
                if printTest   info(iostream, "<$srk> for $sh: closed form $wa against quadrature $wb")   end
            end
        end
    end

    println(iostream, "HydrogenicIon: Dirac l-degeneracy at equal j, fine-structure ordering, the (alpha Z)^2 " *
                      "scaling of the non-relativistic departure, and <r^k> by two independent routes.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_HydrogenicIon()::", success)
    return( success )
end


"""
`TestFrames.testModule_Nuclear(; short::Bool=true)`  ... tests on module Nuclear; a success::Bool is returned.
"""
function testModule_Nuclear(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Nuclear-new.sum")
    printstyled("\n\nTest the module  Nuclear  ... \n", color=:cyan)
    # Nuclear is named 48 times inside other test bodies, but almost always as Nuclear.Model(...) setting up
    # somebody else's computation -- the constructor is exercised and nothing else is. The deformed Fermi shapes
    # added later have never been checked at all. They are tested here by INVERSES and LIMITS, both of which hold
    # exactly and neither of which needs a reference value: a round trip must return what it started from, and a
    # deformed shape at zero deformation must reproduce the spherical one.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # (1) axisRatio and deformationFromAxisRatio are inverses
    for  beta2  in  [-0.3, -0.1, 0.0, 0.05, 0.2, 0.4]
        q  = Nuclear.axisRatio(beta2)
        b2 = Nuclear.deformationFromAxisRatio(q)
        if  abs(b2 - beta2) > 1.0e-10
            success = false
            if printTest   info(iostream, "axisRatio round trip: beta2 = $beta2 -> q = $q -> $b2")   end
        end
    end
    # (2) the spherical limit: at beta2 = 0 the volume factor is exactly 1 and the shape is r-independent
    model0 = Nuclear.DeformedFermiNucleus(beta2 = 0.0)
    vf     = Nuclear.deformedVolumeFactor(model0)
    if  abs(vf - 1.0) > 1.0e-10
        success = false
        if printTest   info(iostream, "deformedVolumeFactor at beta2 = 0 is $vf, must be 1")   end
    end
    for  x  in  [-1.0, -0.5, 0.0, 0.3, 1.0]
        sh = Nuclear.deformedShape(x, model0)
        if  abs(sh - 1.0) > 1.0e-10
            success = false
            if printTest   info(iostream, "deformedShape($x) at beta2 = 0 is $sh, must be 1")   end
        end
    end
    # (3) the two-parameter Fermi radius is invertible: R -> b -> R
    for  R  in  [2.0, 3.5, 5.5, 7.0]
        b  = Nuclear.computeFermiBParameter(R)
        Rx = Nuclear.fermiRrms(b)
        if  abs(Rx - R) > 1.0e-6 * R
            success = false
            if printTest   info(iostream, "Fermi radius round trip: R = $R -> b = $b -> $Rx")   end
        end
    end
    # (4) the DEFORMED radius is invertible too, at zero and at finite deformation
    for  beta2  in  [0.0, 0.25]
        model = Nuclear.DeformedFermiNucleus(beta2 = beta2)
        for  R  in  [3.5, 5.5]
            c  = Nuclear.computeDeformedFermiC(R, 0.524, model)
            Rx = Nuclear.deformedFermiRrms(c, 0.524, model)
            if  abs(Rx - R) > 1.0e-4 * R
                success = false
                if printTest   info(iostream, "deformed radius round trip at beta2 = $beta2: R = $R -> c = $c -> $Rx")   end
            end
        end
    end

    println(iostream, "Nuclear: axisRatio and the Fermi/deformed-Fermi radii are checked as INVERSES, and the " *
                      "deformed shape at beta2 = 0 against the spherical limit. No reference value is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_Nuclear()::", success)
    return( success )
end


"""
`TestFrames.testModule_RadialIntegrals(; short::Bool=true)`  ... tests on module RadialIntegrals; a success::Bool is returned.
"""
function testModule_RadialIntegrals(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-RadialIntegrals-new.sum")
    printstyled("\n\nTest the module  RadialIntegrals  ... \n", color=:cyan)
    # The quadrature is checked against the CLOSED FORMS of the non-relativistic hydrogen atom, where <r>, <r^2>
    # and <1/r> are exact. The non-relativistic route is used on purpose: it removes the O((alpha Z)^2) difference
    # that a Dirac orbital would carry, so any discrepancy here is the quadrature and nothing else. The 1/r case
    # is the one worth having -- rkDiagonal skips the innermost points for negative k, which is a deliberate and
    # documented compromise, and this pins how much it costs.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # rbox = 60 rather than 40, by Rule 12: hydrogenic 3s has its outer turning point at r+ = 18 a.u. and wants a
    # box near 45. Measured at rbox = 40 the truncation alone gives <3s|3s> = 0.9999896, i.e. 1.0e-5 short, which
    # would be read here as a quadrature error when it is nothing of the kind. The box is part of the test.
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 60.0)
    Z    = 1.

    # (1) normalization: <nl|nl> = 1
    for  sh  in  [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3d")]
        P  = HydrogenicIon.radialOrbital(sh, Z, grid)
        wa = RadialIntegrals.overlap(P, P, grid)
        if  abs(wa - 1.0) > 1.0e-6
            success = false
            if printTest   info(iostream, "<$sh|$sh> = $wa, must be 1")   end
        end
    end
    # (2) orthogonality of different n at the same l
    for  (sha, shb)  in  [(Shell("1s"), Shell("2s")), (Shell("1s"), Shell("3s")), (Shell("2s"), Shell("3s"))]
        wa = RadialIntegrals.overlap(HydrogenicIon.radialOrbital(sha, Z, grid),
                                     HydrogenicIon.radialOrbital(shb, Z, grid), grid)
        if  abs(wa) > 1.0e-6
            success = false
            if printTest   info(iostream, "<$sha|$shb> = $wa, must be 0")   end
        end
    end
    # (3) r^k against the closed forms  <r> = (3n^2 - l(l+1))/2Z,  <r^2> = n^2(5n^2 + 1 - 3l(l+1))/2Z^2,
    #     <1/r> = Z/n^2
    for  sh  in  [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3d")]
        n = sh.n;    l = sh.l;    P = HydrogenicIon.radialOrbital(sh, Z, grid)
        for  (k, exact)  in  [( 1, (3n^2 - l*(l+1)) / (2Z)),
                              ( 2, (5n^2 + 1 - 3*l*(l+1)) * n^2 / (2*Z^2)),
                              (-1, Z / n^2)]
            wa = RadialIntegrals.rkDiagonal(k, P, P, grid)
            if  abs(wa - exact) > 1.0e-3 * abs(exact)
                success = false
                if printTest   info(iostream, "<r^$k> for $sh = $wa, closed form $exact")   end
            end
        end
    end

    println(iostream, "RadialIntegrals: overlap normalization and orthogonality, and r^k for k = 1, 2, -1 against " *
                      "the exact non-relativistic hydrogenic closed forms. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_RadialIntegrals()::", success)
    return( success )
end


"""
`TestFrames.testModule_Bsplines(; short::Bool=true)`  ... tests on module Bsplines; a success::Bool is returned.
"""
function testModule_Bsplines(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Bsplines-new.sum")
    printstyled("\n\nTest the module  Bsplines  ... \n", color=:cyan)
    # Bsplines carries the two Rule 12 guards, and a guard that cannot fail is worse than none -- it reads as a
    # clean bill of health. So the guards are tested from BOTH SIDES here: they must pass on a box matched to the
    # orbitals and REFUSE on one deliberately too small. That is the whole point of this test, and it is the part
    # that a single-sided check would have missed for as long as the guards have existed.
    #
    # The basis itself is checked against the exact point-nucleus Dirac energies. That is a genuine closed-form
    # comparison and not a second opinion from the same machinery: Basics.computeDiracEnergy evaluates the
    # analytic formula, while the eigenvalue comes from diagonalizing the B-spline Galerkin matrix.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    subshells = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
    Z         = 20.
    nm        = Nuclear.Model(Z, PointNucleus())
    # Z = 20, so the 2p turning point sits near r+ = 4/20 * 2 = 0.4 a.u.; rbox = 8 is comfortable for all four.
    goodGrid  = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 8.0)
    goodPrim  = Bsplines.generatePrimitives(goodGrid)
    orbitals  = redirect_stdout(devnull) do;  Bsplines.generateOrbitalsHydrogenic(subshells, nm, goodPrim, printout=false)  end

    # (1) the B-spline eigenvalues against the analytic Dirac formula
    for  sh  in  subshells
        exact = Basics.computeDiracEnergy(sh, Z)
        wa    = orbitals[sh].energy
        if  abs(wa - exact) > 1.0e-6 * abs(exact)
            success = false
            if printTest   info(iostream, "$sh: B-spline energy $wa against the Dirac formula $exact")   end
        end
    end
    # (2) the generated orbitals are orthonormal
    for  sh  in  subshells
        wa = RadialIntegrals.overlap(orbitals[sh], orbitals[sh], goodGrid)
        if  abs(wa - 1.0) > 1.0e-6
            success = false
            if printTest   info(iostream, "<$sh|$sh> = $wa, must be 1")   end
        end
    end
    wa = RadialIntegrals.overlap(orbitals[Subshell("1s_1/2")], orbitals[Subshell("2s_1/2")], goodGrid)
    if  abs(wa) > 1.0e-6
        success = false
        if printTest   info(iostream, "<1s_1/2|2s_1/2> = $wa, must be 0")   end
    end
    # (3) BOTH GUARDS PASS on the matched box
    # NOTE THE TWO RETURN TYPES, which are not the same: checkGridRepresentation gives back (ok::Bool, rbox), the
    # second element being the box it would recommend, while checkOrbitalConsistency gives a bare Bool.
    ok1, rboxWanted = redirect_stdout(devnull) do;  Bsplines.checkGridRepresentation(subshells, Z, goodPrim, stopper=false)  end
    ok2 = redirect_stdout(devnull) do;  Bsplines.checkOrbitalConsistency(orbitals, goodGrid, stopper=false)      end
    if  !ok1
        success = false
        if printTest   info(iostream, "checkGridRepresentation refuses a box that is matched to the orbitals")   end
    end
    if  !ok2
        success = false
        if printTest   info(iostream, "checkOrbitalConsistency refuses orbitals generated on a matched box")   end
    end
    # (4) AND THE GRID GUARD REFUSES a box far too small for the same orbitals. Without this the guard could be
    #     returning true unconditionally and everything above would still pass.
    tinyGrid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 0.05)
    tinyPrim = Bsplines.generatePrimitives(tinyGrid)
    ok3, _ = redirect_stdout(devnull) do;  Bsplines.checkGridRepresentation(subshells, Z, tinyPrim, stopper=false)  end
    if  ok3
        success = false
        if printTest   info(iostream, "checkGridRepresentation accepts a box of 0.05 a.u. for Z = 20 orbitals; " *
                                      "the guard is not guarding")   end
    end

    println(iostream, "Bsplines: the Galerkin eigenvalues against the analytic Dirac energies, orbital " *
                      "orthonormality, and BOTH Rule 12 guards exercised from both sides -- passing on a "  *
                      "matched box and refusing a 0.05 a.u. one. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_Bsplines()::", success)
    return( success )
end


"""
`TestFrames.testModule_StarkZeeman(; short::Bool=true)`  ... tests on module StarkZeeman; a success::Bool is returned.
"""
function testModule_StarkZeeman(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-StarkZeeman-new.sum")
    printstyled("\n\nTest the module  StarkZeeman  ... \n", color=:cyan)
    # Three exact statements about a field-dressed level, none of which needs a reference value:
    #
    #   * with both fields off, the eigenvalues must reproduce the zero-field level energies WITH their (2J+1)
    #     degeneracies -- the perturbation matrix is then diagonal and the diagonalization must not disturb it;
    #   * the CENTRE OF GRAVITY cannot move. Diagonalization preserves the trace, and the Zeeman diagonal is
    #     proportional to M, which sums to zero over a level, so the sum of all eigenvalues equals the sum of the
    #     unperturbed sublevel energies EXACTLY, at any field strength. This is the check with teeth: it holds for
    #     any correct implementation and is broken by a wrong sign, a wrong g-factor placement or a double count;
    #   * at fields far below the fine-structure scale the splitting is LINEAR in B, so doubling the field doubles
    #     the spread. That is a statement about the physics rather than about a number.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="StarkZeeman test: B-like C 2p", grid=grid,
                              nuclearModel = Nuclear.Model(6.), configs = [Configuration("1s^2 2s^2 2p")] )
    multiplet = redirect_stdout(devnull) do;  perform(wa; output=true)["multiplet:"]  end
    levels    = multiplet.levels
    nSub      = sum([Basics.twice(lev.J) + 1  for lev in levels])
    sum0      = sum([(Basics.twice(lev.J) + 1) * lev.energy  for lev in levels])

    cache = Dict{Float64, Array{Float64,1}}()
    run(bF) = redirect_stdout(devnull) do
        set = StarkZeeman.Settings(StarkZeeman.Settings(); includeEField = false, includeBField = (bF > 0.),
                                   bField = bF, printBefore = false)
        StarkZeeman.computeOutcomes(multiplet, Nuclear.Model(6.), grid, set; output=true)
    end
    # each field value is computed ONCE; the checks below ask for 0, 1, 2 and 100 T and two of them twice
    energies(bF) = get!(cache, bF) do;  sort([o.energy  for o in run(bF)])  end

    # (1) zero field: the sublevel energies reproduce the level energies with their degeneracies
    out0 = run(0.)
    if  length(out0) != nSub
        success = false
        if printTest   info(iostream, "zero field gives $(length(out0)) sublevels, must be $nSub")   end
    else
        for  o  in  out0
            if  minimum([abs(o.energy - lev.energy)  for lev in levels]) > 1.0e-12
                success = false
                if printTest   info(iostream, "zero-field eigenvalue $(o.energy) matches no level energy")   end
            end
        end
    end
    # (2) the centre of gravity is invariant, at every field strength
    for  bF  in  [0., 1.0, 100.0]
        wb = sum(energies(bF))
        if  abs(wb - sum0) > 1.0e-10 * max(abs(sum0), 1.0)
            success = false
            if printTest   info(iostream, "at B = $bF T the eigenvalues sum to $wb, must be $sum0; the centre " *
                                          "of gravity cannot move")   end
        end
    end
    # (3) the SHIFTS are linear in B. It is the shifts and not the raw spread: 1s^2 2s^2 2p carries a 2p_1/2 /
    #     2p_3/2 fine-structure splitting some four orders larger than a one-tesla Zeeman effect, so
    #     max(E) - min(E) measures the fine structure and barely moves with the field -- measured 1.0133 instead
    #     of 2 when this test was first written. Subtracting the zero-field eigenvalues, level by level in sorted
    #     order, leaves the pure Zeeman shift, and that must double when B does.
    e0 = energies(0.)
    d1 = energies(1.0) .- e0;    d2 = energies(2.0) .- e0
    s1 = maximum(d1) - minimum(d1);    s2 = maximum(d2) - minimum(d2)
    if  s1 <= 0.
        success = false
        if printTest   info(iostream, "a 1 T field produces no splitting at all")   end
    elseif  abs(s2/s1 - 2.0) > 1.0e-3
        success = false
        if printTest   info(iostream, "the Zeeman shifts grow by $(s2/s1) when B doubles, must be 2 in the " *
                                      "linear regime")   end
    end

    println(iostream, "StarkZeeman: the zero-field limit with its degeneracies, the invariance of the centre of " *
                      "gravity at 0, 1 and 100 T, and linearity of the splitting in B. No reference value.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_StarkZeeman()::", success)
    return( success )
end


"""
`TestFrames.testModule_Hamiltonian(; short::Bool=true)`  ... tests on module Hamiltonian; a success::Bool is returned.
"""
function testModule_Hamiltonian(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Hamiltonian-new.sum")
    printstyled("\n\nTest the module  Hamiltonian  ... \n", color=:cyan)
    # WHAT IS DELIBERATELY NOT ASSERTED HERE: that the CI matrix is symmetric. It is not, and that is by design.
    # setupMatrix computes only the UPPER triangle and leaves the lower one at zero, because Basics.diagonalize
    # wraps the matrix in LinearAlgebra.Symmetric(..., :U) and never reads below the diagonal. An `H == transpose(H)`
    # check would therefore FAIL on correct code -- the obvious test for a CI matrix is the wrong one here, and the
    # convention itself is what is pinned instead, so that a later refactor filling the lower triangle inconsistently
    # cannot pass unnoticed.
    #
    # The rest are exact: the trace is invariant under diagonalization, a one-dimensional block must return its own
    # diagonal element, enlarging the CSF space cannot RAISE the ground state (the variational principle), and the
    # kink-aware matrix must agree with the plain one, being the same operator through a different radial quadrature.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(8.)

    multiplet = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="Hamiltonian test", grid=grid, nuclearModel=nm,
                configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]); output=true)["multiplet:"]
    end
    basis    = multiplet.levels[1].basis
    settings = AsfSettings()
    sym      = LevelSymmetry(multiplet.levels[1].J, multiplet.levels[1].parity)
    cache    = InteractionStrength.XLCache()
    matrix   = redirect_stdout(devnull) do;  Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, cache)  end
    n        = size(matrix, 1)

    # (1) the upper-triangle convention: everything strictly below the diagonal is untouched
    below = maximum([abs(matrix[r,s])  for r = 1:n, s = 1:n  if r > s];  init = 0.0)
    if  below != 0.
        success = false
        if printTest   info(iostream, "the lower triangle of the CI matrix is not zero (largest $below); either " *
                                      "the upper-triangle convention changed or it is now filled inconsistently")   end
    end
    # (2) the trace is invariant under diagonalization
    eig    = Hamiltonian.diagonalizeCiMatrix(matrix, LevelSelection())
    trace  = sum([matrix[r,r]  for r = 1:n])
    sumEig = sum(eig.values)
    if  abs(sumEig - trace) > 1.0e-8 * max(abs(trace), 1.0)
        success = false
        if printTest   info(iostream, "the eigenvalues sum to $sumEig but the trace is $trace")   end
    end
    # (3) a 1x1 block must hand back its own diagonal element
    eig1 = Hamiltonian.diagonalizeCiMatrix(reshape([-3.75], 1, 1), LevelSelection())
    if  abs(eig1.values[1] + 3.75) > 1.0e-12
        success = false
        if printTest   info(iostream, "a 1x1 block gives $(eig1.values[1]) instead of -3.75")   end
    end
    # (4) the variational principle: adding CSFs cannot RAISE the ground state
    small = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="Hamiltonian test: small space", grid=grid,
                nuclearModel=nm, configs=[Configuration("1s^2 2s^2")]); output=true)["multiplet:"]
    end
    eSmall = minimum([lev.energy  for lev in small.levels])
    eLarge = minimum([lev.energy  for lev in multiplet.levels])
    if  eLarge > eSmall + 1.0e-10
        success = false
        if printTest   info(iostream, "enlarging the CSF space RAISED the ground state, $eSmall -> $eLarge; " *
                                      "the variational principle forbids it")   end
    end
    # (5) two routes to the same matrix: the kink-aware Slater integral must not change the operator
    matrixK = redirect_stdout(devnull) do
        Hamiltonian.setupMatrixKinkAware(sym, basis, nm, grid, settings, InteractionStrength.XLCache())
    end
    if  size(matrixK) != size(matrix)
        success = false
        if printTest   info(iostream, "the kink-aware matrix has size $(size(matrixK)), the plain one $(size(matrix))")   end
    else
        scale = maximum(abs.(matrix));    dev = maximum(abs.(matrixK .- matrix)) / max(scale, 1.0)
        if  dev > 1.0e-4
            success = false
            if printTest   info(iostream, "the kink-aware and plain CI matrices differ by $dev relative; they " *
                                          "are the same operator through different quadratures")   end
        end
    end

    println(iostream, "Hamiltonian: the upper-triangle convention, trace invariance, the 1x1 block, the "     *
                      "variational bound under CSF enlargement, and the kink-aware matrix against the plain " *
                      "one. Symmetry is NOT asserted -- the lower triangle is zero by design.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_Hamiltonian()::", success)
    return( success )
end


"""
`TestFrames.testModule_InteractionStrength(; short::Bool=true)`  ... tests on module InteractionStrength;
    a success::Bool is returned.
"""
function testModule_InteractionStrength(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-InteractionStrength-new.sum")
    printstyled("\n\nTest the module  InteractionStrength  ... \n", color=:cyan)
    # A TWO-ROUTE test that the file itself sets up: XL_CoulombKinkAware's own docstring says it "computes the same
    # effective Coulomb interaction strength as XL_Coulomb(L, a, b, c, d, grid), including the same triangular-delta
    # veto and angular prefactor", differing only in using RadialIntegrals.SlaterRkKinkAware for the radial integral
    # instead of SlaterRk. Two independent quadratures of one quantity must agree, and nothing outside the module is
    # needed to say so.
    #
    # This matters beyond tidiness: the kink-aware route is what the ALField and EOLField code lines use, through
    # Hamiltonian.setupMatrixKinkAware and SelfConsistent.computeTwoElectronV, so a divergence between the two would
    # move SCF energies while leaving the ordinary CI path untouched -- exactly the kind of split that is hard to
    # attribute afterwards.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(10.)
    orb  = Dict{String,Orbital}()
    for  s  in  ["1s_1/2", "2s_1/2", "2p_1/2", "2p_3/2"]
        orb[s] = HydrogenicIon.orbital(Subshell(s), nm, grid)
    end

    # (1) the two routes agree, over ranks and orbital combinations that include both diagonal and exchange-like
    #     orderings and both spin-orbit partners
    nCompared = 0
    for  (L, sa, sb, sc, sd)  in  [(0, "1s_1/2","1s_1/2","1s_1/2","1s_1/2"), (0, "1s_1/2","2s_1/2","1s_1/2","2s_1/2"),
                                   (0, "2p_1/2","2p_1/2","2p_1/2","2p_1/2"), (1, "1s_1/2","2p_1/2","2p_1/2","1s_1/2"),
                                   (1, "1s_1/2","2p_3/2","2p_3/2","1s_1/2"), (2, "2p_3/2","2p_3/2","2p_3/2","2p_3/2"),
                                   (1, "2s_1/2","2p_3/2","2p_3/2","2s_1/2"), (0, "2p_1/2","2p_3/2","2p_1/2","2p_3/2")]
        wa = InteractionStrength.XL_Coulomb(         L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
        wb = InteractionStrength.XL_CoulombKinkAware(L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
        if  abs(wa) < 1.0e-14  &&  abs(wb) < 1.0e-14    continue    end     # both vetoed; nothing to compare
        nCompared = nCompared + 1
        if  abs(wa - wb) > 1.0e-4 * max(abs(wa), abs(wb))
            success = false
            if printTest   info(iostream, "X^$L($sa,$sb,$sc,$sd): plain $wa against kink-aware $wb")   end
        end
    end
    # A comparison that compared nothing would pass silently, so the count is asserted too.
    if  nCompared < 6
        success = false
        if printTest   info(iostream, "only $nCompared of the eight combinations were non-vanishing; the " *
                                      "comparison is not exercising the routes")   end
    end
    # (2) the selection rules give EXACTLY zero, in both routes. Parity: l_a + l_c + L must be even.
    for  (L, sa, sb, sc, sd)  in  [(1, "1s_1/2","1s_1/2","1s_1/2","1s_1/2"), (1, "2s_1/2","2s_1/2","2s_1/2","2s_1/2"),
                                   (5, "1s_1/2","2p_1/2","2p_1/2","1s_1/2")]
        wa = InteractionStrength.XL_Coulomb(         L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
        wb = InteractionStrength.XL_CoulombKinkAware(L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
        if  wa != 0.  ||  wb != 0.
            success = false
            if printTest   info(iostream, "X^$L($sa,$sb,$sc,$sd) must be exactly 0 by the selection rules; " *
                                          "got $wa and $wb")   end
        end
    end
    # (3) the memoised method must return exactly what the uncached one does -- a cache that drifts from its own
    #     function is silent, and this family is used inside the SCF where a wrong cached value would simply become
    #     the answer.
    cache = InteractionStrength.XLCache()
    for  (L, sa, sb, sc, sd)  in  [(0, "1s_1/2","2s_1/2","1s_1/2","2s_1/2"), (1, "1s_1/2","2p_3/2","2p_3/2","1s_1/2")]
        wa = InteractionStrength.XL_Coulomb(L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
        wb = InteractionStrength.XL_Coulomb(L, orb[sa], orb[sb], orb[sc], orb[sd], grid, cache)
        wc = InteractionStrength.XL_Coulomb(L, orb[sa], orb[sb], orb[sc], orb[sd], grid, cache)   # now from the cache
        if  wa != wb  ||  wb != wc
            success = false
            if printTest   info(iostream, "X^$L($sa,$sb,$sc,$sd) uncached $wa, first cached $wb, second $wc")   end
        end
    end

    println(iostream, "InteractionStrength: XL_Coulomb against XL_CoulombKinkAware over eight (L, abcd) "     *
                      "combinations, exact zeros where the selection rules forbid, and the memoised method "  *
                      "against the plain one. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_InteractionStrength()::", success)
    return( success )
end


"""
`TestFrames.testModule_SelfConsistent(; short::Bool=true)`  ... tests on module SelfConsistent; a success::Bool is returned.
"""
function testModule_SelfConsistent(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-SelfConsistent-new.sum")
    printstyled("\n\nTest the module  SelfConsistent  ... \n", color=:cyan)
    # An SCF has no closed-form answer to be checked against, but it does have a DEFINING PROPERTY: the converged
    # solution is a FIXED POINT. Restarting the cycle from its own output must change nothing. That is exactly what
    # "converged" claims, it needs no reference value, and it is broken by any of the ways an SCF can stop early --
    # a stagnant-step test that fires too soon, an iteration limit reached quietly, a mixing scheme that oscillates
    # about the solution rather than settling on it.
    #
    # The second check is narrower and has a history: settings.frozenSubshells must be honoured EXACTLY. A frozen
    # orbital is not "changed very little", it is not changed at all, so the assertion is bit-identity of P and Q
    # and not a tolerance. The EOL rotation solver ignored this setting until it was fixed earlier today, and
    # nothing else stops that returning; here the DFS path is pinned, and MEASURED to honour it -- the frozen 1s
    # comes back with dE/E = 0.0 and an identical P-vector while the free 2s moves.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(8.);    cfg = [Configuration("1s^2 2s^2")]
    scf(set) = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="SCF test", grid=grid, nuclearModel=nm,
                configs=cfg, asfSettings=set); output=true)["multiplet:"]
    end

    orb1 = scf(AsfSettings()).levels[1].basis.orbitals
    shs  = sort( collect(keys(orb1)) )

    # (1) the converged orbitals are orthonormal
    for  sh  in  shs
        wa = RadialIntegrals.overlap(orb1[sh], orb1[sh], grid)
        if  abs(wa - 1.0) > 1.0e-8
            success = false
            if printTest   info(iostream, "<$sh|$sh> = $wa after the SCF, must be 1")   end
        end
    end
    wa = RadialIntegrals.overlap(orb1[Subshell("1s_1/2")], orb1[Subshell("2s_1/2")], grid)
    if  abs(wa) > 1.0e-8
        success = false
        if printTest   info(iostream, "<1s_1/2|2s_1/2> = $wa after the SCF, must be 0")   end
    end
    # (2) THE FIXED POINT. Restarted from its own output, the cycle must return the same orbitals. Measured
    #     28-Aug-2026 at dE/E = 2e-9 and an overlap of 1.0 to twelve digits; the tolerances below leave a factor
    #     of about fifty. Bit-identity is NOT asserted -- the cycle still performs an iteration and the last
    #     digits move -- so what is pinned is that the physics is stationary, which is the actual claim.
    # (3) A NON-CONVERGED FIELD MUST SAY SO WHERE A SCRIPT CAN SEE IT.  solveMeanFieldBasis and
    #     solveAverageAtomField have raised a COLLECTED warning on non-convergence for a long time; the AL and
    #     EOL fields only printed to the console, so a driver that writes its multiplet to a file recorded
    #     unconverged energies with nothing to mark them.  That is the failure mode this checks: run AL with
    #     maxIterationsScf = 1, which cannot converge, and require that the warning reaches the report.
    #     It is read from jac-warn.report and NOT from Defaults.GBL_WARNINGS, because `perform` ends with
    #     PrintWarnings() followed by ResetWarnings() -- so by the time it returns the array is empty and the
    #     file is the only record.  Reading the file is also what a user actually has.
    let  report = "jac-warn.report"
        rm(report, force = true)
        redirect_stdout(devnull) do
            scf(AsfSettings(AsfSettings(); scField = Basics.ALField(), maxIterationsScf = 1))
        end
        text = isfile(report) ? read(report, String) : ""
        if  !( occursin("solveAverageLevelField", text)  &&  occursin("did NOT converge", text) )
            success = false
            if printTest   info(iostream, "an AL field stopped at maxIterationsScf = 1 left no non-convergence " *
                                          "warning in $report; the file holds: $text")   end
        end
    end

    orb2 = scf(AsfSettings(AsfSettings(); startScfFrom = StartFromPrevious(orb1))).levels[1].basis.orbitals
    for  sh  in  shs
        de = abs(orb2[sh].energy - orb1[sh].energy) / abs(orb1[sh].energy)
        ov = RadialIntegrals.overlap(orb1[sh], orb2[sh], grid)
        if  de > 1.0e-7
            success = false
            if printTest   info(iostream, "restarting from the converged orbitals moved $sh by dE/E = $de; a " *
                                          "converged SCF is a fixed point")   end
        end
        if  abs(abs(ov) - 1.0) > 1.0e-8
            success = false
            if printTest   info(iostream, "the restarted $sh overlaps its own input by $ov, must be 1")   end
        end
    end
    # (3) A FROZEN SUBSHELL DOES NOT MOVE AT ALL -- bit-identical P and Q, and the same energy, not a tolerance.
    #     The free subshell of the same run must still be allowed to move, or "frozen" would be indistinguishable
    #     from "nothing happened".
    frozen = Subshell("1s_1/2");    free = Subshell("2s_1/2")
    orb3   = scf(AsfSettings(AsfSettings(); startScfFrom = StartFromPrevious(orb1),
                             frozenSubshells = [frozen])).levels[1].basis.orbitals
    if  orb3[frozen].P != orb1[frozen].P  ||  orb3[frozen].Q != orb1[frozen].Q  ||
        orb3[frozen].energy != orb1[frozen].energy
        success = false
        if printTest   info(iostream, "the frozen $frozen changed: energy $(orb1[frozen].energy) -> " *
                                      "$(orb3[frozen].energy); frozenSubshells is not being honoured")   end
    end
    if  orb3[free].P == orb1[free].P
        success = false
        if printTest   info(iostream, "the free $free did not move either, so this run proves nothing about " *
                                      "freezing")   end
    end

    println(iostream, "SelfConsistent: orthonormality of the converged orbitals, the FIXED-POINT property under " *
                      "restart from its own output, and frozenSubshells honoured bit-for-bit while a free "       *
                      "subshell still moves. No reference value is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_SelfConsistent()::", success)
    return( success )
end


"""
`TestFrames.testModule_Radial(; short::Bool=true)`  ... tests on module Radial; a success::Bool is returned.
"""
function testModule_Radial(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Radial-new.sum")
    printstyled("\n\nTest the module  Radial  ... \n", color=:cyan)
    # The grid is the one object every other module depends on, and both of its constructions carry a recipe that is
    # written out in the docstring.  A recipe that is only described can drift away from the code silently, so it is
    # re-derived here from the documented formula and compared, rather than compared against a stored number.
    #
    # Note which quantity witnesses the box.  `rbox` is a REQUEST: it enters only through determineNoPoints, which
    # rounds the point count UP to a multiple of orderGL, so the last break point tL[end] overshoots the request by
    # up to a percent or two and is NOT an exact witness.  `hp = rbox/300` is exact, and so is the MINIMALITY of the
    # point count, which is what is tested below.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    #
    # (1) Basics.recommendedGrid: rbox = max over shells of  r_+ + tailFactor*n/Zeff  with
    #     r_+ = (n^2/Zeff)(1 + sqrt(1 - l(l+1)/n^2))  and  Zeff = max(Z - N + 1, Z - slaterScreening, 1),
    #     and hp = rbox/300.  Both are re-derived here from the documented formula.
    occupations = Dict( Shell("1s") => 2, Shell("2s") => 2, Shell("2p") => 6, Shell("3s") => 1 )
    Z = 11.;    NoElectrons = sum( values(occupations) );    tailFactor = 16.;    rMax = 0.
    for  (sh, occ)  in occupations
        if  occ <= 0    continue    end
        n = sh.n;    l = sh.l
        Zeff = max( Z - NoElectrons + 1., Z - Basics.slaterScreening(sh, occupations), 1.0 )
        wa   = (n*n/Zeff) * (1.0 + sqrt( max(0., 1.0 - l*(l+1)/(n*n)) ))  +  tailFactor * n / Zeff
        if  wa > rMax   rMax = wa    end
    end
    grid = Basics.recommendedGrid(occupations, Z)
    if  abs(grid.hp - rMax/300.) > 1.0e-12
        success = false
        if printTest   info(iostream, "recommendedGrid: hp = $(grid.hp), the recipe gives rbox/300 = $(rMax/300.)")   end
    end
    #
    # (2) The two recommendedGrid methods must agree: the configuration route only collects the largest occupation
    #     of every shell and then calls the occupation route.
    grid2 = Basics.recommendedGrid([Configuration("1s^2 2s^2 2p^6 3s")], Nuclear.Model(Z, PointNucleus()))
    if  grid2.hp != grid.hp  ||  grid2.NoPoints != grid.NoPoints
        success = false
        if printTest   info(iostream, "recommendedGrid: the configuration and occupation routes disagree, " *
                                      "hp $(grid2.hp) vs $(grid.hp), NoPoints $(grid2.NoPoints) vs $(grid.NoPoints)")   end
    end
    #
    # (3) An explicit rbox or hp overrides the estimate; that is the user's way in and must be honoured exactly.
    grid3 = Basics.recommendedGrid(occupations, Z; rbox = 40., hp = 0.05)
    if  grid3.hp != 0.05
        success = false
        if printTest   info(iostream, "recommendedGrid: hp = $(grid3.hp) although 0.05 was given explicitly")   end
    end
    #
    # (4) Radial.determineNoPoints returns the SMALLEST multiple of orderGL whose mesh reaches beyond rbox.  The mesh
    #     point r(i) solves  log(r/rnt + 1) + (h/hp) r = (i-1) h,  which is solved here by BISECTION -- deliberately a
    #     different method from the Newton iteration the function itself uses, so the two cannot fail together.
    rOfI = function(rnt::Float64, h::Float64, hp::Float64, i::Int64)
        f  = r -> log(r/rnt + 1.) + (h/hp)*r - (i-1)*h
        lo = 0.;    hi = 1.
        while  f(hi) < 0.    hi = 2hi    end
        for  k = 1:200
            mid = 0.5*(lo+hi)
            if  f(mid) < 0.    lo = mid    else    hi = mid    end
        end
        return( 0.5*(lo+hi) )
    end
    rnt = 2.0e-6;    h = 5.0e-2;    orderGL = 7
    for  (rbox, hp)  in  [(30., 0.1), (5.95, 0.02), (249., 0.83), (12., 0.04), (80., 0.2666)]
        NoPoints = Radial.determineNoPoints(rnt, h, hp, rbox, orderGL)
        rHere    = rOfI(rnt, h, hp, NoPoints);    rBefore = rOfI(rnt, h, hp, NoPoints - orderGL)
        if  rem(NoPoints, orderGL) != 0  ||  rHere <= rbox  ||  rBefore > rbox
            success = false
            if printTest   info(iostream, "determineNoPoints($rbox, $hp) = $NoPoints is not the smallest multiple of " *
                                          "$orderGL reaching beyond the box: r = $rHere, r before = $rBefore")   end
        end
        # the docstring's own count, log(rbox/rnt)/h + rbox/hp, must be reproduced to within one orderGL
        predicted = log(rbox/rnt)/h + rbox/hp + 1.
        if  abs(NoPoints - predicted) > orderGL
            success = false
            if printTest   info(iostream, "determineNoPoints($rbox, $hp) = $NoPoints, the documented formula gives " *
                                          "$predicted")   end
        end
    end
    #
    # (5) The documented reason for scaling hp with the box: hp = rbox/300 holds the point count near 600 whatever the
    #     box.  The docstring quotes 602 points at rbox = 5.1 and 679 at rbox = 249, and both are reproduced exactly.
    for  (rbox, quoted)  in  [(5.1, 602), (249., 679)]
        NoPoints = Radial.determineNoPoints(rnt, h, rbox/300., rbox, orderGL)
        if  NoPoints != quoted
            success = false
            if printTest   info(iostream, "hp = rbox/300 at rbox = $rbox gives $NoPoints points, the docstring " *
                                          "quotes $quoted")   end
        end
    end
    #
    # (6) Radial.generateGrid(maximumPrincipalQN = n) is Rule 12 written as code: rbox = 2.5 r_+.  The resulting grid
    #     must carry exactly the point count that box asks for.
    for  (n, lValue, Zeff)  in  [(4, 0, 1.), (3, 2, 2.5), (5, 1, 3.)]
        rPlus   = (n^2 / Zeff) * (1.0 + sqrt(1.0 - lValue*(lValue+1)/n^2))
        newGrid = redirect_stdout(devnull) do
            Radial.generateGrid(Radial.Grid(false); maximumPrincipalQN = n, lValue = lValue, Zeff = Zeff)
        end
        wanted = Radial.determineNoPoints(newGrid.rnt, newGrid.h, newGrid.hp, 2.5*rPlus, newGrid.orderGL)
        if  newGrid.NoPoints != wanted
            success = false
            if printTest   info(iostream, "generateGrid(n=$n, l=$lValue, Zeff=$Zeff): NoPoints = $(newGrid.NoPoints), " *
                                          "rbox = 2.5 r_+ = $(2.5*rPlus) asks for $wanted")   end
        end
    end
    #
    # (7) Radial.determineZbar samples the last five points of Zr, so a potential with a constant Zr must return that
    #     constant exactly.  Until 12-Aug-2026 the loop read Zr[mtp] rather than Zr[i] and the spread could not fire.
    gridZ = Radial.Grid(true)
    Zbar  = redirect_stdout(devnull) do
        Radial.determineZbar( Radial.Potential("test", 7.0 * ones(length(gridZ.r)), gridZ) )
    end
    if  abs(Zbar - 7.0) > 1.0e-12
        success = false
        if printTest   info(iostream, "determineZbar of a constant Zr = 7 potential returned $Zbar")   end
    end
    #
    # (8) The guards.  A grid that cannot be built must say so rather than return a silently unusable one.
    guards = [ ("no shell",            () -> Basics.recommendedGrid(Dict{Shell,Int64}(), Z)),
               ("Z <= 0",              () -> Basics.recommendedGrid(occupations, -1.)),
               ("rbox <= 0",           () -> Basics.recommendedGrid(occupations, Z; rbox = -3.)),
               ("no configuration",    () -> Basics.recommendedGrid(Configuration[], Nuclear.Model(Z, PointNucleus()))),
               ("no scheme",           () -> Radial.generateGrid(Radial.Grid(false))),
               ("two schemes",         () -> Radial.generateGrid(Radial.Grid(false); boxSize = 10., maximumPrincipalQN = 3,
                                                                 Zeff = 1.)),
               ("n < 1",               () -> Radial.generateGrid(Radial.Grid(false); maximumPrincipalQN = 0, Zeff = 1.)),
               ("Zeff not given",      () -> Radial.generateGrid(Radial.Grid(false); maximumPrincipalQN = 3)),
               ("Zeff <= 0",           () -> Radial.generateGrid(Radial.Grid(false); maximumPrincipalQN = 3, Zeff = 0.)),
               ("l > n-1",             () -> Radial.generateGrid(Radial.Grid(false); maximumPrincipalQN = 2, lValue = 5,
                                                                 Zeff = 1.)) ]
    for  (label, thunk)  in  guards
        raised = false
        try     redirect_stdout(devnull) do
                    thunk()
                end
        catch
            raised = true
        end
        if  !raised
            success = false
            if printTest   info(iostream, "the guard against '$label' did not raise")   end
        end
    end
    #
    # (9) The mesh itself: strictly increasing, positive, and with no negative integration weight.
    if  !all(diff(grid.r) .> 0.)  ||  grid.r[1] <= 0.  ||  !all(grid.wr .>= 0.)
        success = false
        if printTest   info(iostream, "the mesh of the recommended grid is not strictly increasing and positive")   end
    end

    println(iostream, "Radial: the recommendedGrid recipe and the generateGrid box re-derived from their own " *
                      "documented formulae, determineNoPoints checked for minimality against an independent " *
                      "bisection solve of the mesh equation, determineZbar on a constant potential, and ten guards. " *
                      "No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_Radial()::", success)
    return( success )
end


"""
`TestFrames.testModule_LSjj(; short::Bool=true)`  ... tests on module LSjj; a success::Bool is returned.
"""
function testModule_LSjj(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-LSjj-new.sum")
    printstyled("\n\nTest the module  LSjj  ... \n", color=:cyan)
    # For a shell l^N at a fixed J, the jj-LS transformation is a change between two complete orthonormal couplings
    # of the same space, so its matrix is ORTHOGONAL: M'M = M M' = 1.  That is exact, holds table by table, and needs
    # nothing external -- which is what makes it the right test here.  It is also a genuinely different check from
    # the GRASP comparison the tables were validated against, since it tests the tables AGAINST EACH OTHER.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # LinearAlgebra is not among TestFrames' imports, so the distance of a matrix from the unit matrix is measured
    # here directly rather than through I.
    maxOffUnit = function(A)
        devs = 0.
        for  i = 1:size(A,1),  j = 1:size(A,2)
            devs = max( devs, abs( A[i,j] - (i == j ? 1. : 0.) ) )
        end
        return( devs )
    end
    #
    # The matrix of one (l, N, J) block is assembled from the tabulated entries themselves: the distinct LS labels
    # (w, Q, L, S) give the rows, the distinct jj labels (Nm, Qm, Jm, Qp, Jp) the columns, and every element is asked
    # for through getLSjjCoefficient -- so a label pair the table does not carry contributes the 0 it should.
    blockMatrix = function(l::Int64, N::Int64, table, JJ::Int64, mapNm::Bool)
        sub = [me for me in table if me.qn.JJ == JJ]
        LSs = sort(unique( [(me.qn.w, me.qn.QQ, me.qn.LL, me.qn.SS)                                   for me in sub] ))
        jjs = sort(unique( [(mapNm ? 2l - me.qn.Nm : me.qn.Nm, me.qn.Qm, me.qn.Jm, me.qn.Qp, me.qn.Jp) for me in sub] ))
        M   = zeros( length(LSs), length(jjs) )
        for  (i, ls)  in enumerate(LSs),  (j, jj)  in enumerate(jjs)
            qn     = LSjj.LS_jj_qn(ls[1], ls[2], ls[3], ls[4], JJ, jj[1], jj[2], jj[3], jj[4], jj[5])
            M[i,j] = LSjj.getLSjjCoefficient(l, N, qn)
        end
        return( M )
    end
    checkBlocks = function(l::Int64, N::Int64, table, mapNm::Bool)
        for  JJ  in  sort(unique( [me.qn.JJ for me in table] ))
            M = blockMatrix(l, N, table, JJ, mapNm)
            if  size(M,1) != size(M,2)
                success = false
                if printTest   info(iostream, "the l = $l, N = $N, 2J = $JJ block is $(size(M,1)) x $(size(M,2)) and " *
                                              "so cannot be orthogonal")   end
                continue
            end
            devs = max( maxOffUnit(M' * M), maxOffUnit(M * M') )
            if  devs > 1.0e-10
                success = false
                if printTest   info(iostream, "the l = $l, N = $N, 2J = $JJ block is not orthogonal, worst " *
                                              "deviation $devs")   end
            end
        end
    end
    #
    # (1) Every tabulated shell: p^3..p^6, d^3..d^10 and f^3..f^7, 108 blocks in all.
    for  (l, N, table)  in  [ (1, 3, LSjj.LS_jj_p_3),  (1, 4, LSjj.LS_jj_p_4),  (1, 5, LSjj.LS_jj_p_5),
                              (1, 6, LSjj.LS_jj_p_6),  (2, 3, LSjj.LS_jj_d_3),  (2, 4, LSjj.LS_jj_d_4),
                              (2, 5, LSjj.LS_jj_d_5),  (2, 6, LSjj.LS_jj_d_6),  (2, 7, LSjj.LS_jj_d_7),
                              (2, 8, LSjj.LS_jj_d_8),  (2, 9, LSjj.LS_jj_d_9),  (2,10, LSjj.LS_jj_d_10),
                              (3, 3, LSjj.LS_jj_f_3),  (3, 4, LSjj.LS_jj_f_4),  (3, 5, LSjj.LS_jj_f_5),
                              (3, 6, LSjj.LS_jj_f_6),  (3, 7, LSjj.LS_jj_f_7) ]
        checkBlocks(l, N, table, false)
    end
    #
    # (2) f^8..f^13 are NOT tabulated: getLSjjCoefficient reaches them through the Dyall-Grant electron-hole
    #     conjugation, Nm -> 2l - Nm with the phase (-1)^((Qm + Qp - Q)/2).  Their labels are therefore taken from the
    #     conjugate table with Nm mapped back.
    #
    #     WHAT THIS DOES NOT COVER, and it is worth knowing.  The phase's exponent splits as (Qm + Qp)/2 - Q/2, i.e.
    #     into a part fixed by the jj label alone and a part fixed by the LS label alone, so multiplying it in is a
    #     diagonal similarity D1 * M * D2 with D1, D2 = +-1 -- and M'M is then unchanged.  Measured over all four
    #     conjugated shells: every entry factorizes that way, none is an exception.  So NO orthogonality test can see
    #     that phase, and deleting it from the module leaves this test green (checked).  What is tested here is that
    #     the conjugated MAGNITUDES form an orthogonal matrix; the sign convention needs a comparison against the
    #     tables, which is where it was established.
    for  (N, table)  in  [ (8, LSjj.LS_jj_f_6), (9, LSjj.LS_jj_f_5), (10, LSjj.LS_jj_f_4), (11, LSjj.LS_jj_f_3) ]
        checkBlocks(3, N, table, true)
    end
    #
    # (3) The guards: an occupation or an l with no table must say so rather than return the 0. of a missed lookup,
    #     which is indistinguishable from a coefficient that is legitimately zero.
    guards = [ ("p^7",           () -> LSjj.getLSjjCoefficient(1,  7, LSjj.LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0))),
               ("d^11",          () -> LSjj.getLSjjCoefficient(2, 11, LSjj.LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0))),
               ("l = 4 (g)",     () -> LSjj.getLSjjCoefficient(4,  3, LSjj.LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0))),
               ("N = 1, Jm = Jp",() -> LSjj.getLSjjCoefficient(1,  1, LSjj.LS_jj_qn(0, 0, 0, 1, 1, 0, 0, 1, 0, 1))) ]
    for  (label, thunk)  in  guards
        raised = false
        try     thunk()
        catch
            raised = true
        end
        if  !raised
            success = false
            if printTest   info(iostream, "the guard against '$label' did not raise")   end
        end
    end

    println(iostream, "LSjj: the jj-LS transformation matrix of every tabulated block of p^3..p^6, d^3..d^10 and " *
                      "f^3..f^7, together with the f^8..f^13 electron-hole conjugation branch, is orthogonal to " *
                      "machine precision -- 153 blocks. Four guards are checked. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_LSjj()::", success)
    return( success )
end


"""
`TestFrames.testModule_BiOrthogonal(; short::Bool=true)`  ... tests on module BiOrthogonal; a success::Bool is returned.
"""
function testModule_BiOrthogonal(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-BiOrthogonal-new.sum")
    printstyled("\n\nTest the module  BiOrthogonal  ... \n", color=:cyan)
    # The transformation has a defining property, stated in its own docstring: Cleft' * S * Cright = 1 exactly, where
    # S is the overlap matrix of the two orbital sets.  It is an identity rather than a number, so it needs no
    # reference value and cannot drift.  Two bases from DIFFERENT nuclear charges are used, so that S is genuinely far
    # from the unit matrix -- against S = 1 the property would hold for any pair of inverse triangular matrices and
    # the test would be vacuous.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # LinearAlgebra is not among TestFrames' imports, so the distance of a matrix from the unit matrix is measured
    # here directly rather than through I.
    maxOffUnit = function(A)
        devs = 0.
        for  i = 1:size(A,1),  j = 1:size(A,2)
            devs = max( devs, abs( A[i,j] - (i == j ? 1. : 0.) ) )
        end
        return( devs )
    end
    # The configurations are chosen so that no kappa block is 1 x 1: 1s, 2s and 3s make kappa = -1 a 3 x 3 block, and
    # 2p together with 3p makes each of kappa = 1 and kappa = -2 a 2 x 2 one.  A 1 x 1 block has a trivial triangular
    # factorization and would satisfy the identity below for free, so the size of the blocks is part of the test.
    configs = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), Configuration("1s^2 2s 3s"),
               Configuration("1s^2 2s 3p")]
    multipletOf = function(Z::Float64)
        nm   = Nuclear.Model(Z, PointNucleus())
        grid = Basics.recommendedGrid(configs, nm)
        wb   = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = "Z=$Z", grid = grid, nuclearModel = nm,
                                        configs = configs); output = true )
        end
        return( (wb["multiplet:"], grid) )
    end
    leftMp,  grid  = multipletOf(4.)
    rightMp, _     = multipletOf(5.)
    leftBasis      = leftMp.levels[1].basis
    rightBasis     = rightMp.levels[1].basis
    overlapMatrix  = function(basA::Basis, lList::Array{Subshell,1}, basB::Basis, rList::Array{Subshell,1})
        n = length(lList);    S = zeros(n, n)
        for  i = 1:n,  j = 1:n
            S[i,j] = RadialIntegrals.overlap(basA.orbitals[lList[i]], basB.orbitals[rList[j]], grid)
        end
        return( S )
    end
    #
    # (1) The same basis on both sides.  Its orbitals are orthonormal, so S = 1, the triangular factorization is
    #     trivial and both transformation matrices must come back as the unit matrix.
    for  (kappa, (lList, rList, Cleft, Cright))  in  BiOrthogonal.computeTransformationMatrices(leftBasis, leftBasis, grid)
        S = overlapMatrix(leftBasis, lList, leftBasis, rList)
        if  maxOffUnit(S) > 1.0e-8
            success = false
            if printTest   info(iostream, "the orbitals of kappa = $kappa are not orthonormal, worst deviation " *
                                          "$(maxOffUnit(S))")   end
        end
        if  maxOffUnit(Cleft) > 1.0e-10  ||  maxOffUnit(Cright) > 1.0e-10
            success = false
            if printTest   info(iostream, "for one basis with itself the kappa = $kappa transformation is not the " *
                                          "unit matrix")   end
        end
    end
    #
    # (2) Two different bases: the defining property itself.  The overlap is also required to be genuinely non-trivial,
    #     so that a degenerate S cannot make the check pass for the wrong reason.
    biggestOverlapDeviation = 0.
    for  (kappa, (lList, rList, Cleft, Cright))  in  BiOrthogonal.computeTransformationMatrices(leftBasis, rightBasis, grid)
        S = overlapMatrix(leftBasis, lList, rightBasis, rList)
        biggestOverlapDeviation = max( biggestOverlapDeviation, maxOffUnit(S) )
        devs = maxOffUnit( Cleft' * S * Cright )
        if  devs > 1.0e-10
            success = false
            if printTest   info(iostream, "Cleft' * S * Cright is not the unit matrix for kappa = $kappa, worst " *
                                          "deviation $devs")   end
        end
    end
    if  biggestOverlapDeviation < 1.0e-2
        success = false
        if printTest   info(iostream, "the two bases overlap almost perfectly ($biggestOverlapDeviation), so the " *
                                      "biorthonormality check above is vacuous")   end
    end
    #
    # (3) Transforming a multiplet against ITSELF must leave every mixing coefficient where it was; the counter-
    #     rotation of the CI vectors is then the unit matrix as well.
    newLeftMp, newRightMp = redirect_stdout(devnull) do
        BiOrthogonal.computeTransformation(leftMp, leftMp, grid)
    end
    for  i = 1:length(leftMp.levels)
        devs = max( maximum(abs.(newLeftMp.levels[i].mc  - leftMp.levels[i].mc)),
                    maximum(abs.(newRightMp.levels[i].mc - leftMp.levels[i].mc)) )
        if  devs > 1.0e-10
            success = false
            if printTest   info(iostream, "the mixing coefficients of level $i moved by $devs under the identity " *
                                          "transformation")   end
        end
    end
    #
    # (4) THE DIFFERING-DIMENSION CASE, Appendix B of the same paper, implemented 04-Sep-2026.  Until then this
    #     check asserted the opposite -- that the code REFUSES a rectangular overlap -- and it is rewritten here
    #     rather than deleted, because the property it should now assert is the stronger one.  With one subshell
    #     dropped from the right-hand basis, that kappa has n_left > n_right, so this also exercises the MIRROR
    #     branch (the paper writes only the case where the right side is larger).  The condition is no longer
    #     Cleft' * S * Cright = 1 but = (I 0): a unit matrix of the order of the SMALLER side, and zero elsewhere.
    #     The trimmed basis is built by dropping one subshell from the one already computed, so no second SCF is
    #     needed.
    dropped     = leftBasis.subshells[end]
    keptShells  = [sh for sh in leftBasis.subshells if sh != dropped]
    keptOrbitals= Dict( sh => leftBasis.orbitals[sh] for sh in keptShells )
    trimmedBasis= Basis(true, leftBasis.NoElectrons, keptShells, leftBasis.csfs, leftBasis.coreSubshells, keptOrbitals)
    sawRectangular = false
    for  (kappa, (lList, rList, Cleft, Cright))  in  BiOrthogonal.computeTransformationMatrices(leftBasis, trimmedBasis, grid)
        nl = length(lList);    nr = length(rList)
        nl != nr   &&   (sawRectangular = true)
        S    = [RadialIntegrals.overlap(leftBasis.orbitals[lList[i]], trimmedBasis.orbitals[rList[j]], grid)
                for i = 1:nl, j = 1:nr]
        lhs  = Cleft' * S * Cright
        devs = 0.
        for  i = 1:nl,  j = 1:nr
            devs = max( devs, abs(lhs[i,j] - (i == j ? 1.0 : 0.0)) )
        end
        if  devs > 1.0e-10
            success = false
            if printTest   info(iostream, "with a right basis missing $dropped, Cleft' * S * Cright is not (I 0) " *
                                          "for kappa = $kappa (n = $nl vs $nr), worst deviation $devs")   end
        end
    end
    if  !sawRectangular
        success = false
        if printTest   info(iostream, "dropping the subshell $dropped left no kappa with differing dimensions, so " *
                                      "the Appendix B branch was never exercised")   end
    end

    println(iostream, "BiOrthogonal: Cleft' * S * Cright = 1 to machine precision for two bases whose overlap " *
                      "differs from the unit matrix by $(round(biggestOverlapDeviation, digits=3)); the " *
                      "transformation of a basis with itself is the identity, and leaves the mixing coefficients " *
                      "unchanged; and, for bases of DIFFERING dimension per kappa, Cleft' * S * Cright = (I 0), " *
                      "the Appendix B case. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_BiOrthogonal()::", success)
    return( success )
end


"""
`TestFrames.testModule_Continuum(; short::Bool=true)`  ... tests on module Continuum; a success::Bool is returned.
"""
function testModule_Continuum(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Continuum-new.sum")
    printstyled("\n\nTest the module  Continuum  ... \n", color=:cyan)
    # Continuum supplies the free electron to every process module in the package and had no test of its own until
    # 29-Aug-2026. Nothing here is compared against a stored number: a free wave has closed forms for its amplitude,
    # its phase and its derivative, and those are what is asserted. The energy-scale normalization sqrt(2/(pi q)) and
    # the centrifugal phase -l*pi/2 enter every cross section the package prints, so an error in either is an error in
    # all of them at once.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    wc   = Defaults.getDefaults("speed of light: c")
    # An asymptotic step of 0.02 a.u. against a de Broglie wavelength of 4.44 a.u. at 1 Hartree: 222 points per
    # oscillation, comfortably inside the 15 the module demands.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 40.0)
    mtp      = grid.NoPoints - 200
    settings = Continuum.Settings(false, mtp)

    # A guard that cannot fail is worse than no guard, since it reads as a clean bill of health; every refusal below
    # is therefore exercised from both sides.
    function isRefused(f::Function)
        try     redirect_stdout(devnull) do;   f()   end
        catch
            return( true )
        end
        return( false )
    end

    # Elementary spherical Bessel functions, written out in sines and cosines so that the check on the Bessel wave
    # below does not ask GSL to confirm GSL. Stable only away from the origin, which is where they are used.
    function jSmall(l::Int64, z::Float64)
        if      l == 0   return( sin(z)/z )
        elseif  l == 1   return( sin(z)/z^2 - cos(z)/z )
        elseif  l == 2   return( (3/z^3 - 1/z)*sin(z) - 3*cos(z)/z^2 )
        elseif  l == 3   return( (15/z^4 - 6/z^2)*sin(z) - (15/z^3 - 1/z)*cos(z) )
        else    error("no elementary form is kept here for l = $l")
        end
    end

    # (1) gridConsistency returns exactly the point its own docstring names, and refuses all four ways it says it does.
    nrCont = redirect_stdout(devnull) do;  Continuum.gridConsistency(1.0, grid)  end
    if  nrCont != grid.NoPoints - 200
        success = false
        if printTest   info(iostream, "gridConsistency returns $nrCont, but the documented recipe is " *
                                      "NoPoints-200 = $(grid.NoPoints - 200)")   end
    end
    zeroHp   = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 0.0,    rbox = 40.0)
    fewPts   = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 0.2)
    shortBox = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 1.0e-2, hp = 1.0e-3, rbox = 1.5)
    guards   = [ ("a grid with hp = 0 and so no asymptotic region",    () -> Continuum.gridConsistency(1.0,   zeroHp)),
                 ("a mesh far too coarse for the requested energy",    () -> Continuum.gridConsistency(400.0, grid)),
                 ("a grid of fewer than 600 points",                   () -> Continuum.gridConsistency(1.0,   fewPts)),
                 ("a box whose normalization point is inside 2 a.u.",  () -> Continuum.gridConsistency(1.0,   shortBox)) ]
    for  (what, f)  in  guards
        if  !isRefused(f)
            success = false
            if printTest   info(iostream, "gridConsistency accepts $what")   end
        end
    end
    # The threshold is stated as 15 points per oscillation, i.e. hp = wavelength/15 = 0.2962 a.u. at 1 Hartree. Testing
    # it from both sides fixes the guard AT its stated value rather than merely somewhere.
    wavelgth = 2pi / sqrt( 2.0 + Defaults.getDefaults("alpha")^2 )
    justFine = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 0.29, rbox = 400.0)
    justOver = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 0.30, rbox = 400.0)
    if  isRefused(() -> Continuum.gridConsistency(1.0, justFine))
        success = false
        if printTest   info(iostream, "gridConsistency refuses hp = 0.29 although wavelength/15 = $wavelgth/15")   end
    end
    if  !isRefused(() -> Continuum.gridConsistency(1.0, justOver))
        success = false
        if printTest   info(iostream, "gridConsistency accepts hp = 0.30 although wavelength/15 = $(wavelgth/15)")   end
    end

    # (2) THE PURE SINE WAVE AND ITS EXACT INVARIANT. P = N sin(qr - l*pi/2) with N = sqrt(2/(pi q)) obeys
    #     P^2 + (P'/q)^2 = N^2 at EVERY mesh point -- one identity tying the amplitude, the derivative and the phase
    #     together, and true point by point rather than on average.
    for  (energy, sh)  in  [(1.0, Subshell("1000s_1/2")), (1.0, Subshell("1000p_3/2")), (2.5, Subshell("1000d_5/2"))]
        cOrb = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(energy, sh, grid, settings)  end
        q    = sqrt(2energy);    wa = 0.
        for  i = 1:mtp   wa = max(wa, abs( cOrb.P[i]^2 + (cOrb.Pprime[i]/q)^2 - 2/(pi*q) ))   end
        if  wa > 1.0e-14 * 2/(pi*q)
            success = false
            if printTest   info(iostream, "$sh at E = $energy: P^2 + (P'/q)^2 departs from 2/(pi q) by $wa")   end
        end
    end
    # ... and that same invariant, multiplied by q, is 2/pi at every energy. This fixes the EXPONENT of q in the
    # energy-scale normalization, which a single energy cannot see.
    for  energy  in  [0.5, 1.0, 2.0, 5.0]
        cOrb = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(energy, Subshell("1000p_1/2"), grid, settings)  end
        q    = sqrt(2energy);    wa = (cOrb.P[mtp]^2 + (cOrb.Pprime[mtp]/q)^2) * q
        if  abs(wa - 2/pi) > 1.0e-13
            success = false
            if printTest   info(iostream, "at E = $energy the energy-normalization gives q*N^2 = $wa, not 2/pi")   end
        end
    end

    # (3) THE CENTRIFUGAL PHASE, and its SIGN, which no amplitude check can see. sin(x - (l+2)pi/2) = -sin(x - l pi/2),
    #     so the waves for l and l+2 are exact negatives and those for l and l+4 identical. The three subshells below
    #     also carry kappa = -1, +2 and -5, so Basics.subshell_l is exercised for both signs.
    cS = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(1.0, Subshell("1000s_1/2"), grid, settings)  end
    cD = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(1.0, Subshell("1000d_3/2"), grid, settings)  end
    cG = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(1.0, Subshell("1000g_7/2"), grid, settings)  end
    scale = maximum( abs.(cS.P) )
    if  maximum( abs.(cS.P .+ cD.P) ) > 1.0e-12 * scale
        success = false
        if printTest   info(iostream, "the l = 0 and l = 2 pure sine waves are not exact negatives of one another")   end
    end
    if  maximum( abs.(cS.P .- cG.P) ) > 1.0e-12 * scale
        success = false
        if printTest   info(iostream, "the l = 0 and l = 4 pure sine waves do not agree")   end
    end

    # (4) GENERATE THEN NORMALIZE IS THE IDENTITY. A wave that is already energy-normalized must come back with the
    #     phase 0 and the factor 1. The phase does so exactly; the FACTOR does not, and the reason is worth stating,
    #     because it is the difference between a derived bound and a tolerance chosen to make a test pass: the routine
    #     takes the DISCRETE maximum of |P|, and the nearest sample can sit half a step from the crest, so the maximum
    #     falls short by at most (q*hp)^2/8. The honest assertion is therefore 1 <= N <= 1/(1 - (q*hp)^2/8).
    for  (energy, sh)  in  [(1.0, Subshell("1000s_1/2")), (1.0, Subshell("1000p_3/2")),
                            (2.5, Subshell("1000d_5/2")), (0.5, Subshell("1000f_5/2"))]
        cOrb = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(energy, sh, grid, settings)  end
        newOrb, phi, normF = redirect_stdout(devnull) do;  Continuum.normalizeOrbitalPureSine(cOrb, grid, settings)  end
        q     = sqrt(2energy)
        # phi is returned modulo pi, so 0 and pi are the same answer and both must be accepted -- AND that reduction
        # is a REAL BLIND SPOT of this check, not an oversight in it. Flipping the sign of the centrifugal term in
        # normalizeOrbitalPureSine, +l*pi/2 -> -l*pi/2, moves the phase by exactly l*pi, an integer multiple of pi
        # for every l, so the reduced value is unchanged: measured at 4e-13 across l = 0...4, far below the 1e-10
        # asserted here. NO test of this function's returned phase can see that sign. The sign is nevertheless the
        # right one -- an orbital behaving as sin(qr - l*pi/2 + delta) gives at - kr + l*pi/2 = delta -- and it rests
        # on that derivation, not on anything below.
        if  min( abs(phi), abs(pi - phi) ) > 1.0e-10
            success = false
            if printTest   info(iostream, "$sh at E = $energy: renormalizing an already-normalized sine gives the " *
                                          "phase $phi, which is neither 0 nor pi")   end
        end
        wa = maximum( abs.(newOrb.P) ) / maximum( abs.(cOrb.P) )
        if  wa < 1.0 - 1.0e-12   ||   wa > 1/(1 - (q*grid.hp)^2/8)
            success = false
            if printTest   info(iostream, "$sh at E = $energy: renormalizing gives N = $wa, outside the sampling " *
                                          "bound [1, $(1/(1 - (q*grid.hp)^2/8))]")   end
        end
    end

    # (5) AND THE NORMALIZATION REFUSES what it cannot reach. |N| > 30 means the orbital never grew to its asymptotic
    #     amplitude inside the box, and the routine says so instead of returning a factor of a thousand.
    cOrb  = redirect_stdout(devnull) do;  Continuum.generateOrbitalPureSine(1.0, Subshell("1000s_1/2"), grid, settings)  end
    small = Orbital( cOrb.subshell, cOrb.isBound, cOrb.useStandardGrid, cOrb.energy,
                     cOrb.P/1000, cOrb.Q/1000, cOrb.Pprime/1000, cOrb.Qprime/1000, cOrb.grid )
    if  !isRefused(() -> Continuum.normalizeOrbitalPureSine(small, grid, settings))
        success = false
        if printTest   info(iostream, "normalizeOrbitalPureSine accepts an orbital needing N = 1000, although it " *
                                      "refuses above 30")   end
    end

    # (6) THE BESSEL WAVE AGAINST ELEMENTARY CLOSED FORMS: P = r j_l(qr) and dP/dr = (l+1) j_l(z) - z j_{l+1}(z), with
    #     j_0 ... j_3 written out above in sines and cosines. The derivative is the half that matters -- it is what the
    #     phase is read from downstream, and it carried a missing factor z until 29-Aug-2026.
    for  sh  in  [Subshell("1000s_1/2"), Subshell("1000p_3/2"), Subshell("1000d_5/2")]
        cOrb = redirect_stdout(devnull) do;  Continuum.generateOrbitalBessel(1.0, sh, grid, settings)  end
        q    = sqrt(2.0);    l = Basics.subshell_l(sh);    wa = 0.;   wb = 0.
        for  i = mtp-800:mtp-100
            z  = q * grid.r[i]
            wa = max(wa, abs( cOrb.P[i]      - grid.r[i]*jSmall(l, z) ))
            wb = max(wb, abs( cOrb.Pprime[i] - ((l+1)*jSmall(l, z) - z*jSmall(l+1, z)) ))
        end
        if  wa > 1.0e-12
            success = false
            if printTest   info(iostream, "$sh: the Bessel wave departs from r*j_l(qr) by $wa")   end
        end
        if  wb > 1.0e-12
            success = false
            if printTest   info(iostream, "$sh: d/dr of the Bessel wave departs from (l+1) j_l - z j_(l+1) by $wb")   end
        end
    end

    # (7) THE ASYMPTOTIC COULOMB WAVE: its four arrays must belong to ONE wave. The asymptotic ratio of the small to
    #     the large component of a free Dirac wave is R = sqrt(E/(E+2c^2)) -- a textbook result, not something read off
    #     this code -- so with P = A cos(theta) and Q = A R sin(theta) the arrays satisfy two identities in theta,
    #         P^2 + (Q/R)^2 = const        and        Q'*Q = -R^2 * P'*P,
    #     the second of which needs no knowledge of dtheta/dr and is exactly what a wrong amplitude on a derivative
    #     violates. Q' carried the LARGE-component amplitude until 29-Aug-2026, which is this identity out by 1/R^2.
    nm     = Nuclear.Model(10.0, PointNucleus())
    pot    = redirect_stdout(devnull) do;  Nuclear.nuclearPotential(nm, grid)  end
    energy = 1.0;      R = sqrt( energy / (energy + 2*wc^2) )
    for  sh  in  [Subshell("1000s_1/2"), Subshell("1000p_3/2"), Subshell("1000d_5/2")]
        cOrb = redirect_stdout(devnull) do;  Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)  end
        amp2 = cOrb.P[mtp]^2 + (cOrb.Q[mtp]/R)^2;     wa = 0.;   wb = 0.;   scale = 0.
        for  i = 2:mtp
            wa    = max(wa,    abs( cOrb.P[i]^2 + (cOrb.Q[i]/R)^2 - amp2 ))
            wb    = max(wb,    abs( cOrb.Qprime[i]*cOrb.Q[i] + R^2 * cOrb.Pprime[i]*cOrb.P[i] ))
            scale = max(scale, abs( R^2 * cOrb.Pprime[i]*cOrb.P[i] ))
        end
        if  wa > 1.0e-12 * amp2
            success = false
            if printTest   info(iostream, "$sh: P^2 + (Q/R)^2 is not constant along the asymptotic Coulomb wave; " *
                                          "it moves by $wa against $amp2")   end
        end
        if  wb > 1.0e-12 * scale
            success = false
            if printTest   info(iostream, "$sh: Q'*Q + R^2 P'*P = $wb against a scale of $scale; the small and large " *
                                          "components do not belong to one wave")   end
        end
    end

    # (8) twoFzero, the 2F0 asymptotic series translated from Salvat's RADIAL. It has one property that can be checked
    #     exactly: for a negative-integer first argument the series TERMINATES, so 2F0(-n, b; 1/z) is a polynomial of
    #     degree n in 1/z and can be summed here in three lines.
    for  (na, b, z)  in  [(1, 2.5+0.7im, 30.0+4.0im), (2, 1.3-0.4im, 25.0+0.0im), (3, 0.5+0.0im, 40.0+10.0im)]
        a  = ComplexF64(-na, 0.);      wa = Continuum.twoFzero(a, ComplexF64(b), ComplexF64(z))
        wb = 0.0+0.0im;                wx = 1.0+0.0im
        for  k = 0:na    wb = wb + wx;    wx = wx * (a+k)*(b+k)/((k+1)*z)    end
        if  abs(wa - wb) > 1.0e-14 * abs(wb)
            success = false
            if printTest   info(iostream, "twoFzero(-$na, $b; 1/$z) = $wa against the terminating series $wb")   end
        end
    end

    println(iostream, "Continuum: the energy-scale normalization and the centrifugal phase of the free wave as the " *
                      "exact invariant P^2 + (P'/q)^2 = 2/(pi q); generate-then-normalize as an identity, with the " *
                      "renormalization factor held inside a bound DERIVED from the mesh; the Bessel wave and its "  *
                      "derivative against elementary sines and cosines; the asymptotic Coulomb wave's small and "   *
                      "large components tied by the free-Dirac ratio sqrt(E/(E+2c^2)); twoFzero against a "         *
                      "terminating 2F0; and all five refusals -- four in gridConsistency and one in "               *
                      "normalizeOrbitalPureSine -- exercised from both sides. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_Continuum()::", success)
    return( success )
end


"""
`TestFrames.testMethod_ThomasReicheKuhn(; short::Bool=true)`  ... tests the Thomas-Reiche-Kuhn sum rule for the
    dipole oscillator strengths of hydrogen; a success::Bool is returned.
"""
function testMethod_ThomasReicheKuhn(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ThomasReicheKuhn-new.sum")
    printstyled("\n\nTest the Thomas-Reiche-Kuhn sum rule  ... \n", color=:cyan)
    # THE ONE CHECK IN THIS SUITE THAT PINS AN ABSOLUTE SCALE. Everything else asserted about a cross section or an
    # oscillator strength here is a RELATION -- a ratio, a symmetry, an invariance -- and three separate scale
    # errors survived for weeks in 2026 precisely because every check they met was of that kind: the plasma
    # photoionization normalisation (item 64, wrong by 4 pi/(alpha omega)^2), MultipolePolarizibility (challenge
    # 113, wrong by ~2600 with the two gauges agreeing EXACTLY with each other), and example-Ak.jl, which ran for
    # years while reporting only counts. A sum rule cannot be satisfied by a wrongly scaled spectrum.
    #
    # Thomas-Reiche-Kuhn: the dipole oscillator strengths from a given state sum to the number of electrons. For
    # hydrogen that is 1, and the BOUND states carry only 0.5650 of it -- the continuum holds the remaining
    # 0.4350. So the rule is only reachable if the B-spline pseudo-continuum is included, and it is that
    # requirement which makes this an absolute test rather than a comparison.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # NuclearField, not the default DFS: a ONE-electron system in a mean field built from its own density is
    # repelled by itself (see SelfConsistent.checkOneElectronSelfInteraction). gridStopper = false is the POINT
    # rather than a workaround -- the grid guard correctly reports that the high-n np states are not physical
    # bound states on this box, which is exactly what a pseudo-continuum is.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 120.0)
    nm   = Nuclear.Model(1., PointNucleus())
    asf  = AsfSettings(AsfSettings(); scField = Basics.NuclearField(), gridStopper = false)
    peS  = PhotoExcitation.Settings(PhotoExcitation.Settings(); multipoles=[E1],
                                    gauges=[UseCoulomb, UseBabushkin], printBefore=false)
    wa   = Atomic.Computation(Atomic.Computation(), name="Thomas-Reiche-Kuhn", grid=grid, nuclearModel=nm,
               initialConfigs = [Configuration("1s")],                        initialAsfSettings = asf,
               finalConfigs   = [Configuration("$(n)p")  for n = 2:20],       finalAsfSettings   = asf,
               processSettings = peS)
    lines = redirect_stdout(devnull) do;  perform(wa; output=true)["photo-excitation lines:"]  end

    # (1) the ground state is the exact hydrogenic one; without this the rest could be a sum rule of the wrong atom
    if  abs(lines[1].initialLevel.energy + 0.5) > 1.0e-4
        success = false
        if printTest   info(iostream, "the 1s level is at $(lines[1].initialLevel.energy), not -0.5")   end
    end

    # (2) THE INDIVIDUAL OSCILLATOR STRENGTHS AGAINST THE CLOSED FORM. Summed over the two fine-structure
    #     components, f(1s->np) is 0.4162, 0.0791, 0.0290, 0.0139, 0.0078 for n = 2..6 -- textbook hydrogen, and
    #     an exact statement rather than a stored number.
    # Rounded to THREE decimals, not more: the two fine-structure components of one n differ in omega only at the
    # sixth (0.375005 against 0.375006), so a finer key splits the doublet and compares one component against the
    # pair's total. Three decimals still separates n = 2..7 cleanly (0.375, 0.444, 0.469, 0.480, 0.486, 0.490).
    byOmega = Dict{Float64,Vector{EmProperty}}()
    for  l in lines   push!( get!(byOmega, round(l.omega, digits=3), EmProperty[]), l.oscStrength )   end
    fSorted = [ sum(v).Coulomb  for v in [byOmega[k] for k in sort(collect(keys(byOmega)))] ]
    for  (i, fExact)  in  enumerate([0.41620, 0.07910, 0.02899, 0.01394, 0.00780])
        if  abs(fSorted[i] - fExact) > 1.0e-3 * fExact
            success = false
            if printTest   info(iostream, "f(1s->$(i+1)p) = $(fSorted[i]) against the exact $fExact")   end
        end
    end

    # (3) THE FINE-STRUCTURE RATIO. f(np_3/2)/f(np_1/2) is 2 by statistical weight, exactly in the
    #     non-relativistic limit; the residue is the O((alpha Z)^2) = 5.3e-5 correction that item 46 measured, so
    #     1e-3 is the right tolerance and a factor error would be caught while that correction is not.
    for  k  in  keys(byOmega)
        v = byOmega[k];   length(v) == 2 || continue
        r = max(v[1].Coulomb, v[2].Coulomb) / min(v[1].Coulomb, v[2].Coulomb)
        if  abs(r - 2.0) > 1.0e-3
            success = false
            if printTest   info(iostream, "the doublet at omega = $k has f-ratio $r, not 2")   end
        end
    end

    # (4) THE PSEUDO-CONTINUUM MUST ACTUALLY BE THERE, or (5) below would be a statement about bound states only.
    #     Measured: 22 of the 38 final levels lie ABOVE threshold.
    nPositive = count(l -> l.finalLevel.energy > 0., lines)
    if  nPositive < 15
        success = false
        if printTest   info(iostream, "only $nPositive of $(length(lines)) final levels are above threshold; the " *
                                      "sum below would then be over bound states alone")   end
    end

    # (5) THE SUM RULE ITSELF, from both sides. It must EXCEED 0.5650 -- the sum over ALL bound states, which no
    #     bound spectrum can pass -- and it must not exceed 1, the electron number, which no spectrum may pass.
    #     Measured here: 0.688. Extending the basis walks it up to 0.958 at n <= 60, and it approaches 1 from
    #     below as it must.
    sumC = sum(l.oscStrength.Coulomb   for l in lines)
    sumB = sum(l.oscStrength.Babushkin for l in lines)
    for  (lab, s)  in  [("Coulomb", sumC), ("Babushkin", sumB)]
        if  s <= 0.5650
            success = false
            if printTest   info(iostream, "the $lab oscillator-strength sum is $s, at or below the bound-state " *
                                          "limit 0.5650, so the pseudo-continuum is contributing nothing")   end
        end
        if  s >= 1.0
            success = false
            if printTest   info(iostream, "the $lab oscillator-strength sum is $s, which EXCEEDS the electron " *
                                          "number: Thomas-Reiche-Kuhn is violated")   end
        end
    end

    # (6) and the two gauges agree, because with NuclearField the orbitals are the exact hydrogenic ones
    if  abs(sumC - sumB) > 1.0e-6 * sumC
        success = false
        if printTest   info(iostream, "the two gauges give $sumC and $sumB for the sum")   end
    end

    println(iostream, "Thomas-Reiche-Kuhn: hydrogen 1s -> np, n = 2..20, on a B-spline basis in which " *
                      "$nPositive of $(length(lines)) final levels lie above threshold. The oscillator strengths " *
                      "reproduce the closed-form 0.4162, 0.0791, 0.0290, 0.0139, 0.0078 for n = 2..6; the "      *
                      "fine-structure doublets carry them in the statistical ratio 2; and the sum comes to "     *
                      "$(round(sumC, digits=6)), which EXCEEDS the all-bound-states limit of 0.5650 -- so the "  *
                      "pseudo-continuum is contributing -- while staying below the electron number of 1, which " *
                      "the sum rule forbids it to pass. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testMethod_ThomasReicheKuhn()::", success)
    return( success )
end
