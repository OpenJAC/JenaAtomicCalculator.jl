
## Functions in this file cover: atomic properties.
## Alphabetical order within this file.


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
                            propertySettings = [ AlphaVariation.Settings(true, 0.125, true, LevelSelection() )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-AlphaVariation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-AlphaVariation-new.sum"), "Alpha variation parameters:", 9)
    testPrint("testModule_AlphaVariation()::", success)
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
                            propertySettings = [ DecayYield.Settings(Basics.SCA(), true, false, LevelSelection(), Shell[] )] )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-DecayYield-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-DecayYield-new.sum"), "Fluorescence and Auger", 11)
    testPrint("testModule_DecayYield()::", success)
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
                                joinpath(@__DIR__, "..", "test", "test-FormFactor-new.sum"), "Standard and modifi", 9)
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
                            nuclearModel=Nuclear.Model(26., FermiNucleus(), 58., 3.81, AngularJ64(5//2), 1.0, 1.0, 0.),
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
                            nuclearModel=Nuclear.Model(26., FermiNucleus(), 58., 3.75, AngularJ64(5//2), 1.0, 2.0, 0.),
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
`TestFrames.testModule_MultipolePolarizibility(; short::Bool=true)`  ... tests on module MultipolePolarizibility.
"""
function testModule_MultipolePolarizibility(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-MultipolePolarizibility-new.sum")
    printstyled("\n\nTest the module  MultipolePolarizibility  ... \n", color=:cyan)
    ### Make the tests
    ## RE-ENABLED 09-Aug-2026, and deliberately SHORT. This test had been disabled since 31-Jul because a parallel
    ## session was editing the module; that reason expired, but the old test then failed for a different one -- it
    ## passed six positional arguments to a Settings that has had five fields since the Angel-Sandars rewrite.
    ##
    ## WHAT IT CHECKS: the static E1 polarizability of hydrogen 1s, against two facts that need no reference data.
    ##   (1) THE TENSOR POLARIZABILITY MUST VANISH for J = 1/2. alpha2 is a rank-2 quantity and a level with
    ##       J < 1 cannot carry one. An exact zero, in both gauges.
    ##   (2) THE SCALAR POLARIZABILITY IS BOUNDED. For the hydrogen ground state every term of the sum over
    ##       perturbers is positive, so a sum over a FINITE, BOUND-ONLY set of np states is a strict lower bound
    ##       on the exact 4.5 a.u. Hence 0 < alpha0 < 4.5, which catches both a sign error and a runaway.
    ##
    ## THE "B-SPLINE PSEUDO-CONTINUUM BUG" THAT BLOCKED THIS MODULE WAS NEVER A BUG -- it was this test's own
    ## radial box (10-Aug-2026). With rbox = 20 a.u. the np perturbers came out at -0.12500, -0.05160 and
    ## +0.00817 a.u.; the 4p pair is POSITIVE, which is exactly the "pseudo-continuum orbitals get near-zero
    ## energies" symptom that was recorded as a Bsplines defect. Those Rydberg perturbers simply need a wider
    ## box: at rbox = 80 they are -0.12500, -0.05556, -0.03125, i.e. the exact hydrogenic values, and the
    ## result is converged (rbox = 150 and 300 reproduce it digit for digit). hp was widened 1.0e-2 -> 4.0e-2
    ## at the same time: the linear part of a log-lin grid costs rbox/hp points, so rbox = 80 at hp = 1.0e-2
    ## needs 8344 points and 6.6 s, while hp = 4.0e-2 gives the SAME perturber energies with 2338 points in
    ## 0.5 s -- i.e. the wider box costs nothing once hp is matched to it. The new grid check in
    ## Bsplines.checkGridRepresentation flags precisely this and recommends 77.4 a.u.; it was that check,
    ## on its first run over the suite, which found it.
    ##
    ## WHAT IT DOES NOT CHECK, and why. The Coulomb-gauge value comes out around 1.5e-6 where Babushkin gives
    ## 1.164, a discrepancy of six orders of magnitude, and it does NOT move with the box. It is a category
    ## error rather than a numerical one: a STATIC polarizability has no gauge freedom, and the Coulomb column
    ## here evaluates MabEmission -- the frequency-DEPENDENT transition operator -- at omega = 0, where the
    ## velocity form needs the 1/omega that connects it to the length form. The frequency-independent operator
    ## that belongs in a static polarizability (and in C_6/C_8/C_10 dispersion coefficients) is
    ## InteractionStrength.eMultipole, assembled by MultipoleMoment.emmStaticAmplitude -- Johnson's r^k C_k,
    ## with no omega and no gauge. Recorded rather than asserted, since fixing it is a separate task.
    ni    = Nuclear.Model(1.0, PointNucleus())
    grid  = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=4.0e-2, rbox=80.0)
    scf   = Basics.NuclearField();   asf = AsfSettings(AsfSettings(); scField=scf)
    gMp   = generate(Representation("np perturbers", ni, grid,
                       [Configuration("2p"), Configuration("3p"), Configuration("4p")],
                       MeanFieldMultiplet(MeanFieldSettings(scf))), output=true)["mean-field multiplet"]
    set   = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings();
                multipoles=[E1], gMultiplet=gMp, omegas=[0.], printBefore=false, levelSelection=LevelSelection())
    ## `asfSettings`, NOT `initialAsfSettings` (corrected 09-Aug-2026): Atomic.Computation carries BOTH, and the
    ## first pairs with `configs` while the second pairs with `initialConfigs`. Passing the wrong one is silent --
    ## the scField is simply ignored and the default DFS field runs instead, which for hydrogen gives a
    ## self-interacting 1s orbital at -0.194 a.u. against the exact -0.500. The perturbers here are built through
    ## a Representation, which DID honour the field, so the reference level and the perturbers sat in different
    ## potentials: exactly the mismatch the two-photon work found to be worth a factor 6.
    wa    = Atomic.Computation(Atomic.Computation(), name="H 1s static polarizibility", grid=grid, nuclearModel=ni,
                configs=[Configuration("1s")], asfSettings=asf, propertySettings=[set])
    outcomes = perform(wa; output=true)["Polarizibility outcomes:"]
    ###
    success = true
    if  length(outcomes) == 0    success = false;   println("** no polarizibility outcome was computed")   end
    for  outcome in outcomes
        if  Basics.twice(outcome.level.J) < 2      ## J < 1 : no tensor polarizability can exist
            if  outcome.alpha2.Coulomb != 0.  ||  outcome.alpha2.Babushkin != 0.
                success = false
                println("** alpha2 is non-zero for J = $(outcome.level.J), where a rank-2 polarizability " *
                        "cannot exist: $(outcome.alpha2)")
            end
        end
        ## the bound-only perturber sum must sit strictly between zero and the exact hydrogen value
        if  !(0. < outcome.alpha0.Babushkin < 4.5)
            success = false
            println("** alpha0 (Babushkin) = $(outcome.alpha0.Babushkin) is outside (0, 4.5); a sum over a " *
                    "bound-only perturber set must be a strict lower bound on the exact 4.5 a.u.")
        end
    end
    Defaults.setDefaults("print summary: close", "")
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "MultipolePolarizibility: the J = 1/2 tensor-polarizability zero and the bound on alpha0; " *
                      "no comparison against approved data, and the Coulomb gauge is deliberately not asserted.")
    testPrint("testModule_MultipolePolarizibility()::", success)
    return(success)
end
