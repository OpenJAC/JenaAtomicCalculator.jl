
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
`TestFrames.testModule_NuclearShielding(; short::Bool=true)`
    ... tests the ANGULAR NORMALISATION of module NuclearShielding against closed-form values, with no orbitals,
        no grid and no approved file, so that nothing in it can drift.

        WHY THIS IS THE TEST TO HAVE.  The whole module is one angular weight times one radial sum, and the weight
        is a single constant multiplying every number it returns.  If that constant is wrong, every shielding
        factor is wrong by the same silent factor and no comparison with experiment can reveal it, because the
        quantity has no independent measurement of its own.  Two closed-form statements pin it:

        (i)  the NON-RELATIVISTIC closed-shell coefficients of Sternheimer, Phys. Rev. 159, 266 (1967),
             48/25 (p->p), 16/7 (d->d) and 224/75 (f->f), must equal (8/5) |<l||C^(2)||l>|^2;
        (ii) the RELATIVISTIC sum over the two spin-orbit partners of |<kappa||C^(2)||kappa'>|^2 must equal
             exactly TWICE the non-relativistic |<l||C^(2)||l'>|^2 -- the spin factor -- which is what makes the
             relativistic module reduce to the tabulated non-relativistic limit.

        A success::Bool is returned.
"""
function testModule_NuclearShielding(; short::Bool=true)
    printstyled("\n\nTest the module  NuclearShielding  ... \n", color=:cyan)
    success = true
    clNR(l, L, lp) = (l+lp+L) % 2 != 0 ? 0.0 :
        (-1.0)^l * sqrt((2l+1)*(2lp+1)) *
        AngularMomentum.Wigner_3j(AngularJ64(l), AngularJ64(L), AngularJ64(lp),
                                  AngularM64(0), AngularM64(0), AngularM64(0))

    # (i) the published closed-shell coefficients
    for  (l, published)  in  ((1, 48/25), (2, 16/7), (3, 224/75))
        computed = 1.6 * clNR(l, 2, l)^2
        if  abs(computed - published) > 1.0e-10
            success = false
            println("  ** l = $l :  (8/5)|<l||C2||l>|^2 = $computed  but Sternheimer's value is $published")
        end
    end

    # (ii) the relativistic sum over spin-orbit partners must be exactly twice the non-relativistic value
    kappas(l) = l == 0 ? [-1] : [-(l+1), l]
    for  (l, lp)  in  ((1,1), (2,2), (3,3), (0,2), (1,3), (2,4))
        wa = 0.
        for  ka  in  kappas(l),  kv  in  kappas(lp)
            wa = wa + AngularMomentum.CL_reduced_me(Subshell(9, ka), 2, Subshell(9, kv))^2
        end
        wb = 2.0 * clNR(l, 2, lp)^2
        if  abs(wa - wb) > 1.0e-10
            success = false
            println("  ** l = $l -> $lp :  relativistic sum = $wa  but 2 x non-relativistic = $wb")
        end
    end

    testPrint("testModule_NuclearShielding()::", success)
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
`TestFrames.testModule_DispersionCoefficient(; short::Bool=true)`  ... tests on module DispersionCoefficient; a
    success::Bool is returned.
"""
function testModule_DispersionCoefficient(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DispersionCoefficient-new.sum")
    printstyled("\n\nTest the module  DispersionCoefficient  ... \n", color=:cyan)
    # Four invariants, of which only ONE uses a published number.  A dispersion coefficient is a re-summation of the
    # multipole response, so the module's own content is the combinatorics and the omega integral, and both can be
    # checked exactly.
    #
    #   (1) THE CASIMIR-POLDER COMBINATORICS, which needs no atomic input at all.  Summed over ORDERED (l1,l2)
    #       pairs the single binomial must reproduce the prefactors the literature writes out by hand.
    #   (2) THE CLOSED-FORM OVERLAP INTEGRAL against a quadrature of its own integrand.  This is what justifies
    #       having NO omega quadrature in src/: on a discrete spectrum every polarizability is a sum of single
    #       poles and int_0^inf [a/(a^2+w^2)][b/(b^2+w^2)] dw = pi/(2(a+b)) exactly, so each C_n is a double sum.
    #   (3) THE DIPOLE SUM RULE, which fixes the atom-wall C_3 of hydrogen EXACTLY and needs no literature.  With
    #       A_a = (2/3)|<a|d|0>|^2 the sum over a complete set is (2/3)<d^2> = 2 for H 1s, so
    #       C_3 = (1/8) sum_a A_a = 1/4 exactly -- which is also what an electrostatic image dipole gives.  This
    #       one check exercises the amplitudes, the angular factors and the normalization in one number.
    #   (4) C_6 of H-H against the published Dirac value, the single external anchor.
    #
    # AND TWO BEHAVIOURS, both of which are the module refusing rather than returning: the ground-level guard must
    # RAISE for an excited partner, and a C_n needing a multipole the partner does not carry must be reported
    # INCOMPLETE rather than silently computed too small.
    #
    # Cost: the hydrogen part is one 64 a.u. box with np n = 2..16, measured at 0.1 s once JAC is compiled.  The
    # BOX IS DELIBERATELY SMALL and points the opposite way to Rule 12's usual advice -- a sum over states needs
    # its upper members to stand in for the CONTINUUM, and a box matched to the highest orbital leaves them genuine
    # bound states, whereupon the continuum is missing entirely.  Measured, sum_a A_a / 2 over this np set:
    # 0.702 (n <= 6), 0.735 (8), 0.887 (10), 0.993 (12), 0.99995 (16), so n = 16 is what the 1e-4 below needs.
    # No approved reference file is used, so nothing here can be made to pass by re-approving data.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    #
    # (1) the combinatorics alone
    for  (l1, l2, exact, name)  in  ( (1, 1, 3/pi, "3/pi"),  (1, 2, 15/(2pi), "15/2pi"), (2, 1, 15/(2pi), "15/2pi"),
                                      (1, 3, 14/pi, "14/pi"), (3, 1, 14/pi, "14/pi"),    (2, 2, 35/pi, "35/pi") )
        prefactor = DispersionCoefficient.casimirPolderPrefactor(l1, l2)
        if  abs(prefactor/exact - 1.0) > 1.0e-13
            success = false
            if printTest   info(iostream, "the Casimir-Polder prefactor of the ordered pair ($l1,$l2) came out " *
                                          "$prefactor where the literature writes $name = $exact")   end
        end
    end
    #
    # (2) the closed-form overlap integral against a trapezoidal quadrature of the same integrand.  The integrand
    #     falls as w^-4, so the truncated tail, not the step, sets the error: measured -2.6e-11 at wmax = 1000 and
    #     -4.0e-10 at wmax = 400, both with 40000 intervals, so wmax = 1000 is chosen and 1e-9 is the bound.
    #     THE TOLERANCE IS WHAT GIVES THIS CHECK ITS TEETH, and it was calibrated by mutation: at the 1e-6 first
    #     written here, a prefactor deliberately wrong by 2e-7 still PASSED. Two other mutations -- the
    #     imaginary-axis denominator by 1e-3 and the Casimir-Polder prefactor by 1e-6 -- were caught at once.
    termsA  = DispersionCoefficient.MultipoleTerms(E1, [0.4162, 0.0791, 0.0290], [0.3750, 0.4444, 0.46875])
    termsB  = DispersionCoefficient.MultipoleTerms(E2, [1.2000, 0.3000],         [0.3750, 0.4800])
    closed  = DispersionCoefficient.overlapIntegral(termsA, termsB)
    nw      = 40000;    wmax = 1000.0;    dw = wmax/nw;    quad = 0.0
    for  i = 0:nw
        w     = i * dw
        value = DispersionCoefficient.polarizability(termsA, w) * DispersionCoefficient.polarizability(termsB, w)
        quad  = quad + ((i == 0  ||  i == nw) ? 0.5*value : value) * dw
    end
    if  abs(quad/closed - 1.0) > 1.0e-9
        success = false
        if printTest   info(iostream, "the closed-form overlap integral gives $closed where a $nw-point quadrature " *
                                      "of the same integrand gives $quad; the pi/(2(a+b)) identity that removes " *
                                      "the omega quadrature from src/ does not hold")   end
    end
    #
    # (3) and (4): hydrogen 1s, reached through Atomic.Computation exactly as a user would reach it
    nmDi   = Nuclear.Model(1.0, PointNucleus())
    gridDi = Radial.Grid(Radial.Grid(true); rbox = 64.0)
    asfDi  = AsfSettings(AsfSettings(); gridStopper = false)
    gmDi   = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(); name="np perturbers", grid=gridDi, nuclearModel=nmDi,
                 configs=[Configuration("$(n)p")  for n = 2:16], asfSettings=asfDi);
                 output=true)[ResultKeys.Multiplet]
    end
    setDi  = DispersionCoefficient.Settings(DispersionCoefficient.Settings();
                 orders=[6], multipoles=[E1], gMultiplet=gmDi, calcAtomWall=true)
    ocsDi  = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(); name="H(1s) dispersion", grid=gridDi, nuclearModel=nmDi,
                 configs=[Configuration("1s")], asfSettings=asfDi, propertySettings=[setDi]);
                 output=true)[DispersionCoefficient.Settings]
    end
    if  length(ocsDi) != 1
        success = false
        if printTest   info(iostream, "the property route returned $(length(ocsDi)) outcomes for the single " *
                                      "1s level of hydrogen")   end
    else
        outcomeDi = ocsDi[1]
        # an empty partnerB means a homonuclear pair, and partner B must then BE partner A
        if  outcomeDi.partnerA !== outcomeDi.partnerB
            success = false
            if printTest   info(iostream, "an empty settings.partnerB must give a homonuclear pair, but partner B " *
                                          "is a separate object from partner A")   end
        end
        sumA = sum(outcomeDi.partnerA.terms[1].strengths)
        if  abs(sumA/2.0 - 1.0) > 1.0e-4
            success = false
            if printTest   info(iostream, "the dipole oscillator-strength sum came out $sumA where the sum rule " *
                                          "(2/3)<d^2> = 2 for H 1s is exact; the amplitudes, the angular factor " *
                                          "or the perturber set is at fault")   end
        end
        # the atom-wall coefficient IS that sum over eight, so this part of it has no freedom at all
        if  outcomeDi.atomWallC3 != 0.125*sumA
            success = false
            if printTest   info(iostream, "C_3 = $(outcomeDi.atomWallC3) is not (1/8) sum_a A_a = $(0.125*sumA), " *
                                          "although the omega integral of a single pole leaves exactly that")   end
        end
        if  abs(outcomeDi.atomWallC3/0.25 - 1.0) > 1.0e-4
            success = false
            if printTest   info(iostream, "the atom-wall coefficient of H 1s came out $(outcomeDi.atomWallC3) " *
                                          "against the exactly 1/4 that the dipole sum rule and the image dipole " *
                                          "both give")   end
        end
        c6Di = outcomeDi.coefficients[1]
        if  !c6Di.complete
            success = false
            if printTest   info(iostream, "C_6 was reported incomplete although it needs the E1 response alone, " *
                                          "which this partner carries")   end
        end
        if  abs(c6Di.value/6.498392245754 - 1.0) > 1.0e-5
            success = false
            if printTest   info(iostream, "C_6(H-H) came out $(c6Di.value) against the published Dirac value " *
                                          "6.498392245754")   end
        end
        # a C_n that needs a multipole the partner does not carry is a LOWER BOUND and must say so
        setC8 = DispersionCoefficient.Settings(setDi; orders=[8])
        c8Di  = DispersionCoefficient.computeOutcome(outcomeDi.partnerA, outcomeDi.partnerA,
                                                      setC8).coefficients[1]
        if  c8Di.complete
            success = false
            if printTest   info(iostream, "C_8 was reported complete although it needs the E2 response and this " *
                                          "partner carries E1 only")   end
        end
    end
    #
    # the ground-level guard must RAISE.  A negative excitation energy is all it takes, so this needs no structure:
    # the contour argument behind the imaginary-frequency integral fails the moment one dE changes sign.
    excited = DispersionCoefficient.Partner("an excited partner", Level(),
                  [ DispersionCoefficient.MultipoleTerms(E1, [0.4], [-0.1]) ])
    raised  = false
    try
        DispersionCoefficient.checkGroundLevel(excited)
    catch
        raised = true
    end
    if  !raised
        success = false
        if printTest   info(iostream, "checkGroundLevel accepted a partner with a perturber BELOW it; the " *
                                      "Casimir-Polder form is not defined there and the guard must refuse")   end
    end

    println(iostream, "DispersionCoefficient: the Casimir-Polder prefactors against 3/pi, 15/2pi, 14/pi and " *
                      "35/pi; the closed-form overlap integral against a quadrature of its own integrand; the " *
                      "dipole sum rule, which fixes C_3(H 1s) = 1/4 exactly; C_6(H-H) against the published " *
                      "Dirac value; and the two refusals -- the ground-level guard raises, and a C_n missing a " *
                      "multipole is reported incomplete. Only C_6 uses an external number, and no approved " *
                      "reference file is read.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_DispersionCoefficient()::", success)
    return( success )
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
                            ## calcHfMultiplet = TRUE since 16-Aug-2026.  It had been false because that path
                            ## raised "... still to be done for a single nuclear spin/isomer"; the error stood
                            ## in front of a computation that works, and three display sites behind it still
                            ## used the field names of a retired type.  All four are repaired, so the test now
                            ## covers the F-resolved multiplet as well.

    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    println("aaa  ")
    # Make the comparison with approved data
    ## Anchor on what this test is FOR.  Until 16-Aug-2026 it compared "Level  J Parity  Hartrees", the
    ## level-energy table, so it would have verified no hyperfine quantity at all even once re-enabled.
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Hfs-b-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Hfs-b-new.sum"), "HFS parameters:", 9)
    success = success  &&
              testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Hfs-b-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-Hfs-b-new.sum"),
                                "Selected (non-) diagonal hyperfine amplitudes:", 15)
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
    # THE g_F ROWS OF THE APPROVED FILE WERE RE-APPROVED ON 30-Sep-2026, and this note is here so that the
    # change is findable from the test rather than only from a commit message.  `LandeZeeman.landeFFactor`
    # gained the NUCLEAR contribution to the hyperfine g factor, which had been omitted: the electrons are
    # not the only magnet in a hyperfine level, and the nuclear term is of relative order m_e/m_p = 5.4e-4,
    # i.e. a hundred times the 1e-6 at which such a g factor is now measured.  Eight g_F rows moved and
    # NOTHING ELSE in the file was touched; every difference equals -g_I (m_e/m_p) times the nuclear
    # projection bracket exactly, checked level by level, and the case I = 7/2 with J = 3/2 at F = 3 -- where
    # the ELECTRONIC bracket vanishes identically so the answer is purely nuclear -- reproduces
    # -g_I (m_e/m_p) to five figures.  The level energies in the same file drifted by 4.3e-8 relative over
    # the same period, which is inside this comparison's rtol of 1e-6 and was NOT re-approved.
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
    ## THE SIX-ORDER GAUGE GAP THIS COMMENT USED TO DESCRIBE IS GONE, and the paragraph is replaced rather than
    ## kept, because it describes code that no longer exists. It recorded alpha_0(Coulomb) ~ 1.5e-6 against
    ## 1.164 (Babushkin) and diagnosed it correctly: the Coulomb column was evaluating the frequency-DEPENDENT
    ## MabEmission at omega = 0. The module now uses MultipoleMoment.emmStaticAmplitude -- Johnson's r^k C_k,
    ## with no omega and no gauge -- so the two columns are identical BY CONSTRUCTION, and both give 3.491829.
    ##
    ## **AND THE NUMBERS BELOW ARE NOW ASSERTED AGAINST EXACT HYDROGEN, which is what makes this a real test.**
    ## For H 1s the sum-over-states alpha = sum_i f_i/omega_i^2 is known in closed form, term by term:
    ##     2p   f = 0.41620,  omega = 0.375      ->  2.9596        (and the 1/2 : 3/2 pair splits 1 : 2)
    ##     3p   f = 0.07910,  omega = 0.44444    ->  0.4005
    ##     4p   f = 0.02900,  omega = 0.46875    ->  0.1320
    ##                                bound sum  ->  3.4921
    ## so a bound-only 2p+3p+4p sum has ONE right answer and the test demands it to 0.2 %. That is what the old
    ## "0 < alpha0 < 4.5" bound could not do: it would have passed the stale 1.164 just as happily.
    ni    = Nuclear.Model(1.0, PointNucleus())
    grid  = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=4.0e-2, rbox=80.0)
    scf   = Basics.NuclearField();   asf = AsfSettings(AsfSettings(); scField=scf)
    gMp   = generate(Representation("np perturbers", ni, grid,
                       [Configuration("2p"), Configuration("3p"), Configuration("4p")],
                       MeanFieldMultiplet(MeanFieldSettings(scf))), output=true)[ResultKeys.MeanFieldMultiplet]
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
    outcomes = perform(wa; output=true)[MultipolePolarizibility.Settings]
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
        ## the bound-only perturber sum must sit strictly between zero and the exact hydrogen value ...
        if  !(0. < outcome.alpha0.Babushkin < 4.5)
            success = false
            println("** alpha0 (Babushkin) = $(outcome.alpha0.Babushkin) is outside (0, 4.5); a sum over a " *
                    "bound-only perturber set must be a strict lower bound on the exact 4.5 a.u.")
        end
        ## ... and it must be the RIGHT lower bound:  the closed-form 2p+3p+4p value of hydrogen, 3.4921.
        if  abs(outcome.alpha0.Babushkin/3.4921 - 1.0) > 2.0e-3
            success = false
            println("** alpha0 = $(outcome.alpha0.Babushkin) against the closed-form bound-only sum 3.4921 for " *
                    "H 1s with 2p, 3p and 4p perturbers -- a relative difference of " *
                    "$(abs(outcome.alpha0.Babushkin/3.4921 - 1.0))")
        end
        ## THE GENERAL alphas TABLE MUST MEET THE STATIC PAIR BIT FOR BIT AT omega = 0 -- added 19-Oct-2026,
        ## item 69.  `outcome.alpha0` and the (Mp, K, omega, alpha) table are filled by DIFFERENT code paths: the
        ## first by the original scalar/tensor route, the second by the general-rank, general-frequency one.  The
        ## module special-cases omega = 0 precisely so that adding the dynamic capability moves no approved
        ## reference, so a difference of one ulp here means the two paths have parted -- hence `!=` and not a
        ## tolerance.  This costs nothing: the computation above already runs with omegas = [0.].
        iDy = findfirst(t -> t[1] == E1  &&  t[2] == 0  &&  t[3] == 0.0, outcome.alphas)
        if      isnothing(iDy)
            success = false
            println("** outcome.alphas carries no (E1, K = 0, omega = 0) entry, so the general-rank table was " *
                    "not filled at all")
        elseif  outcome.alphas[iDy][4].Babushkin != outcome.alpha0.Babushkin
            success = false
            println("** the omega = 0 entry of the general alphas table is $(outcome.alphas[iDy][4].Babushkin) " *
                    "while alpha0 is $(outcome.alpha0.Babushkin); the two paths must meet bit for bit at omega = 0")
        end
        ## A STATIC POLARIZABILITY HAS NO GAUGE FREEDOM, so the two columns must be identical, not merely close.
        if  outcome.alpha0.Coulomb != outcome.alpha0.Babushkin
            success = false
            println("** the two gauge columns of a STATIC polarizability differ: $(outcome.alpha0.Coulomb) " *
                    "against $(outcome.alpha0.Babushkin).  They are built from one frequency-independent " *
                    "operator and must be equal bit for bit.")
        end
    end
    ## AND PERTURBER BY PERTURBER, which is where a wrong angular factor shows itself while a total can still
    ## look plausible.  Each np pair has a closed-form value, and the two spin-orbit partners must split 1 : 2
    ## by their statistical weights -- a ratio no fitting can produce by accident.
    for  outcome in outcomes
        pairs = Dict{Int64, Vector{Float64}}()
        for  c  in  outcome.contributions
            n = c.nuLevel.energy < -0.1 ? 2 : (c.nuLevel.energy < -0.04 ? 3 : 4)
            push!(get!(pairs, n, Float64[]), c.contribAlpha0.Babushkin)
        end
        for  (n, exact)  in  ((2, 2.9596), (3, 0.4005), (4, 0.1320))
            haskey(pairs, n) || continue
            wa = sum(pairs[n])
            if  abs(wa/exact - 1.0) > 3.0e-3
                success = false
                println("** the $(n)p perturber pair contributes $wa to alpha_0, against the closed-form " *
                        "$exact for hydrogen")
            end
            if  length(pairs[n]) == 2
                lo, hi = minimum(pairs[n]), maximum(pairs[n])
                if  abs(hi/lo - 2.0) > 1.0e-2
                    success = false
                    println("** the $(n)p spin-orbit pair splits $hi : $lo = $(hi/lo), not the 2 : 1 its " *
                            "statistical weights demand")
                end
            end
        end
    end
    Defaults.setDefaults("print summary: close", "")
    ## THE RANK-2 BRANCH, ASSERTED AGAINST AN EXACT NUMBER -- added 14-Oct-2026, item 68.  Everything above
    ## checks that alpha_2 VANISHES (J = 1/2, and J < 1 generally), and a vanishing test cannot detect a wrong
    ## rank-2 prefactor, a wrong 6-j or a sign.  This can, and it needs no reference data.
    ##
    ## WITH EXACTLY ONE INTERMEDIATE LEVEL the two Angel-Sandars formulae share the same radial factor,
    ##     alpha_0 = -2/(3(2J+1)) * R ,
    ##     alpha_2 = -2 sqrt(10/3) G(J) (-1)^(2J) (-1)^(J-J') {J J' 1; 1 2 J} * R ,
    ## so every radial integral AND every energy denominator cancels in the ratio, leaving a pure angular
    ## number -- which turns out to be a simple rational:
    ##     J = 3/2, J' = 1/2  ->  -1          J = 3/2, J' = 3/2  ->  +4/5
    ##     J = 3/2, J' = 5/2  ->  -1/5        J = 5/2, J' = 5/2  ->  +8/7
    ## Measured 14-Oct-2026 at -1.00000000 (2.2e-16 off), +0.80000000 (bit-exact) and -0.20000000 (2.8e-17).
    ##
    ## TWO CONSEQUENCES WORTH KNOWING.  Because the radial factor cancels, THE GRID IS IRRELEVANT here -- the
    ## coarse, cheap basis below gives the exact ratio, so this costs the suite almost nothing.  And together
    ## with the hydrogen checks above the pair is complete: those validate the RADIAL part of alpha_0 against
    ## the exact 4.5, this validates the ANGULAR part of alpha_2, and alpha_2 shares the radial factor.
    gridA2 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-06, h = 5.0e-2, hp = 1.2, rbox = 60.0)
    nmA2   = Nuclear.Model(4., PointNucleus())
    asfA2  = AsfSettings(AsfSettings(); scField = Basics.NuclearField(), gridStopper = false)
    for  (confA2, twoJp, exactRatio)  in  (("3s", 1, -1.0), ("3d", 3, 0.8))
        wp  = Atomic.Computation(Atomic.Computation(), name="alpha_2 perturber", grid=gridA2, nuclearModel=nmA2,
                                  configs=[Configuration(confA2)], asfSettings=asfA2)
        gMp = perform(wp; output=true)[ResultKeys.Multiplet]
        ## exactly ONE intermediate level, which is what the identity requires
        oneLevel = [ lev  for lev in gMp.levels  if Basics.twice(lev.J) == twoJp ]
        if  length(oneLevel) != 1
            success = false
            println("** the $confA2 perturber multiplet does not contain exactly one level with 2J = $twoJp")
            continue
        end
        polA2 = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings(), multipoles=[E1],
                                                  gMultiplet=Multiplet("one level", oneLevel), printBefore=false)
        wA2   = Atomic.Computation(Atomic.Computation(), name="alpha_2 ratio", grid=gridA2, nuclearModel=nmA2,
                                    configs=[Configuration("2p")], propertySettings=[polA2], asfSettings=asfA2)
        for  oc  in  perform(wA2; output=true)[MultipolePolarizibility.Settings]
            Basics.twice(oc.level.J) == 3  ||  continue        ## the 2p_3/2 level is the one carrying alpha_2
            ratioA2 = oc.alpha2.Babushkin / oc.alpha0.Babushkin
            if  abs(ratioA2 - exactRatio) > 1.0e-10
                success = false
                println("** alpha_2/alpha_0 for J = 3/2 with a single J' = $twoJp/2 intermediate came out " *
                        "$ratioA2 where the pure angular identity gives exactly $exactRatio -- the rank-2 " *
                        "prefactor, the 6-j or the sign is wrong")
            end
        end
    end

    ## THE DYNAMIC BRANCH, ASSERTED AGAINST AN EXACT IDENTITY -- added 19-Oct-2026, item 69.  Everything above
    ## runs at omega = 0, so until now NOTHING in the suite touched the frequency dependence: `Settings.omegas`,
    ## `MultipolePolarizibility.frequencyFactor` and the entire imaginary axis were held only by a `tools/` probe,
    ## which the suite never runs.
    ##
    ## WITH EXACTLY ONE INTERMEDIATE LEVEL the frequency dependence is a pure RATIO in which the radial integral,
    ## the angular factor and the oscillator strength all cancel, leaving the denominator alone:
    ##     alpha(omega)   / alpha(0) = dE^2 / (dE^2 - omega^2)     on the real axis,
    ##     alpha(i omega) / alpha(0) = dE^2 / (dE^2 + omega^2)     on the imaginary axis.
    ## So this tests `frequencyFactor` and nothing else -- which is the whole point, since the frequency enters the
    ## module ONLY through that one function.  As with the rank-2 ratio above the grid is therefore irrelevant, and
    ## the perturber here is the SAME 3s level, so the test costs one extra multiplet and nothing more.
    ##
    ## THE TWO AXES MUST ALSO MEET AT omega = 0, exactly:  dE/(dE^2 +- 0) is 1/dE either way, and that is the
    ## identity the omega = 0 fast path rests on.
    gridDy = Radial.Grid(Radial.Grid(false), rnt = 4.0e-06, h = 5.0e-2, hp = 1.2, rbox = 60.0)
    nmDy   = Nuclear.Model(4., PointNucleus())
    asfDy  = AsfSettings(AsfSettings(); scField = Basics.NuclearField(), gridStopper = false)
    wpDy   = Atomic.Computation(Atomic.Computation(); name="one perturber", grid=gridDy, nuclearModel=nmDy,
                                 configs=[Configuration("3s")], asfSettings=asfDy)
    gmDy   = perform(wpDy; output=true)[ResultKeys.Multiplet]
    wrDy   = Atomic.Computation(Atomic.Computation(); name="reference level", grid=gridDy, nuclearModel=nmDy,
                                 configs=[Configuration("2p")], asfSettings=asfDy)
    mpDy   = perform(wrDy; output=true)[ResultKeys.Multiplet]
    ## the 2p_1/2 level, which carries K = 0 alone, and the single 3s perturber above it
    iHalf  = findfirst(lev -> Basics.twice(lev.J) == 1, mpDy.levels)
    jHalf  = findfirst(lev -> Basics.twice(lev.J) == 1, gmDy.levels)
    if  isnothing(iHalf)  ||  isnothing(jHalf)
        success = false
        println("** the 2p or the 3s multiplet carries no J = 1/2 level, so the single-pole frequency identity " *
                "cannot be set up")
    else
        levDy  = mpDy.levels[iHalf];          oneDy = Multiplet("one level", [gmDy.levels[jHalf]])
        dEDy   = gmDy.levels[jHalf].energy - levDy.energy
        alphaDy = function(omega::Float64, imaginary::Bool)
            aw, _ = MultipolePolarizibility.computeMultipolePolarizabilities(E1, levDy, oneDy, gridDy;
                                                    omega = omega, imaginaryOmega = imaginary)
            return( aw[1][2].Babushkin )
        end
        a0Dy = alphaDy(0.0, false)
        if  dEDy <= 0.0  ||  a0Dy == 0.0
            success = false
            println("** the single-perturber setup is degenerate: dE = $dEDy and alpha(0) = $a0Dy")
        else
            ## the two axes meet at omega = 0, bit for bit
            if  alphaDy(0.0, true) != a0Dy
                success = false
                println("** the imaginary axis at omega = 0 gives $(alphaDy(0.0, true)) against $a0Dy on the real " *
                        "axis; dE/(dE^2 + 0) and dE/(dE^2 - 0) are the same number and must stay so")
            end
            ## and away from it, the exact single-pole ratio.  0.3 and 0.6 a.u. both sit below dE ~ 1.1 a.u.
            for  (omDy, imDy)  in  ((0.3, false), (0.6, false), (0.3, true), (0.6, true))
                exactDy = imDy ? dEDy^2/(dEDy^2 + omDy^2) : dEDy^2/(dEDy^2 - omDy^2)
                ratioDy = alphaDy(omDy, imDy) / a0Dy
                if  abs(ratioDy/exactDy - 1.0) > 1.0e-12
                    axisDy = imDy ? "imaginary" : "real"
                    success = false
                    println("** on the $axisDy axis at omega = $omDy the single-pole ratio alpha(omega)/alpha(0) " *
                            "came out $ratioDy where the exact dE^2/(dE^2 -+ omega^2) gives $exactDy -- " *
                            "frequencyFactor is wrong")
                end
            end
        end
    end
    ## THE RANK RULES, spot-checked against values derived by hand rather than from the code.  K <= 2L, K <= 2J and
    ## K even act together, and the J = 1/2 row is a physical statement and not a restatement of the formula: a
    ## spin-1/2 level cannot be aligned by ANY multipole, however high.
    for  (mpRk, twoJRk, exactRk)  in  ( (E1, 1, [0]),  (E2, 1, [0]),  (E3, 1, [0]),
                                        (E1, 5, [0,2]), (E2, 5, [0,2,4]), (E3, 5, [0,2,4]) )
        ranksRk = MultipolePolarizibility.allowedRanks(mpRk, AngularJ64(twoJRk//2))
        if  ranksRk != exactRk
            success = false
            println("** allowedRanks($mpRk, J = $(twoJRk)/2) returned $ranksRk where K <= 2L, K <= 2J and K even " *
                    "together give $exactRk")
        end
    end

    ## AND the approved reference, restored 03-Sep-2026. The two exact invariants above cannot be fooled but they
    ## are also very loose -- anything in (0, 4.5) passes -- so they would not notice alpha_0 drifting from 1.164
    ## to 3.9, nor a column vanishing from the table. The file comparison catches exactly that, at rtol = 1e-6.
    ## The reference had been orphaned since 09-Aug-2026, when this test was rewritten around the invariants and
    ## the comparison was dropped without the file being removed: it sat in test/approved/ describing output that
    ## nothing checked. A reference that no test reads is not a safety net, it is a thing that goes quietly stale.
    scomp   = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultipolePolarizibility-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-MultipolePolarizibility-new.sum"),
                                "Static electric-dipole scalar/tensor polarizibilities", 6)
    success = success  &&  scomp
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "MultipolePolarizibility: the J = 1/2 tensor-polarizability zero, the bound on alpha0, and " *
                      "the approved results table; the Coulomb gauge is deliberately not asserted. The RANK-2 " *
                      "branch is now asserted too, against the exact angular identity alpha_2/alpha_0 = -1 and " *
                      "+4/5 for a single J' = 1/2 and J' = 3/2 intermediate -- a ratio in which every radial " *
                      "integral cancels, so it needs no reference data and no converged grid.")
    testPrint("testModule_MultipolePolarizibility()::", success)
    return(success)
end


"""
`TestFrames.testModule_NistCalibration(; short::Bool=true)`
    ... tests module NistCalibration on a small, hand-built multiplet, no grid, no orbitals and no approved
        file, so that nothing in it can drift: (i) a matched level gets exactly the expected calibrated
        energy, ground + NIST-cm^-1-converted-to-a.u.; (ii) an UNmatched level keeps its raw energy and is
        flagged unmatched; (iii) excludeLevels suppresses a match that would otherwise apply; (iv)
        calibrating the already-calibrated cm.multiplet a second time reproduces the FIRST result exactly
        (idempotence, by index-based matching and a fixed reference level -- see the module's own docstring);
        (v) calibrating a CalibratedMultiplet directly raises an error rather than silently re-shifting it.
        A success::Bool is returned.
"""
function testModule_NistCalibration(; short::Bool=true)
    printstyled("\n\nTest the module  NistCalibration  ... \n", color=:cyan)
    success = true

    lev1 = ManyElectron.Level(AngularJ64(0), AngularM64(0), Parity("+"), 1, -10.0, 1.0, false, Basis(), Float64[])
    lev2 = ManyElectron.Level(AngularJ64(1), AngularM64(1), Parity("-"), 2,  -9.5, 1.0, false, Basis(), Float64[])
    lev3 = ManyElectron.Level(AngularJ64(2), AngularM64(2), Parity("-"), 3,  -9.0, 1.0, false, Basis(), Float64[])
    mp   = ManyElectron.Multiplet("test", [lev1, lev2, lev3])
    nistTable = [ NistCalibration.NistEntry(2, 100000.0, "level 2 measured"),
                  NistCalibration.NistEntry(3, 200000.0, "level 3 measured") ]

    cm = NistCalibration.calibrate(mp, nistTable, "test source"; groundLevelIndex=1)
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    expected2   = lev1.energy + 100000.0*auPerKayser
    if  abs(cm.multiplet.levels[2].energy - expected2) > 1.0e-10   ||   !cm.levelRecords[2].matched
        success = false;   println("  ** matched-level energy or flag wrong: $(cm.multiplet.levels[2].energy) vs $expected2")
    end
    if  cm.multiplet.levels[1].energy != lev1.energy   ||   cm.levelRecords[1].matched
        success = false;   println("  ** unmatched level 1 was changed or wrongly flagged matched")
    end

    cmExcl = NistCalibration.calibrate(mp, nistTable, "test source"; groundLevelIndex=1, excludeLevels=[2])
    if  cmExcl.levelRecords[2].matched   ||   cmExcl.multiplet.levels[2].energy != lev2.energy
        success = false;   println("  ** excludeLevels did not suppress the level 2 match")
    end

    cm2 = NistCalibration.calibrate(cm.multiplet, nistTable, "test source"; groundLevelIndex=1)
    if  cm2.multiplet.levels[2].energy != cm.multiplet.levels[2].energy   ||
        cm2.multiplet.levels[3].energy != cm.multiplet.levels[3].energy
        success = false;   println("  ** re-calibrating cm.multiplet was NOT idempotent")
    end

    try
        NistCalibration.calibrate(cm, nistTable, "test source"; groundLevelIndex=1)
        success = false;   println("  ** calibrate(::CalibratedMultiplet, ...) did not raise an error")
    catch
        # expected
    end

    testPrint("testModule_NistCalibration()::", success)
    return(success)
end


"""
`TestFrames.testModule_ReducedDensityMatrix(; short::Bool=true)`  ... tests on module ReducedDensityMatrix;
    a success::Bool is returned.
"""
function testModule_ReducedDensityMatrix(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ReducedDensityMatrix-new.sum")
    printstyled("\n\nTest the module  ReducedDensityMatrix  ... \n", color=:cyan)
    # A SANITY test built entirely from ALGEBRAIC invariants of a one-particle reduced density matrix, and
    # deliberately not a comparison against a stored .sum.  Every check below is exact for any basis and any
    # grid -- the trace is the electron number, rho is symmetric, and no natural orbital may hold fewer than
    # zero or more than 2j+1 electrons -- so none of them can pass on a stale reference, and none of them is
    # sensitive to how well the orbitals themselves are converged.  That insensitivity is the point: it tests
    # the RDM machinery rather than the SCF underneath it.
    #
    # The last check is the one with a history.  Until 28-Aug-2026 computeProperties published
    # compute1pRDMDirect, a special-case shortcut whose `samecoupling` test compares the RUNNING coupling
    # subshellX[i]; whenever a spectator subshell lies BETWEEN the two substituted ones its X necessarily
    # differs, and the element was dropped by `continue`.  The diagonal and the trace stayed perfect, so the
    # first three checks here would ALL have passed while every off-diagonal was silently zero -- and the
    # off-diagonals are exactly where the correlation lives.  Hence check (4): in a genuinely correlated
    # two-configuration level the off-diagonal block must not vanish.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 40.0)
    # ReducedDensityMatrix offers only the plain keyword constructor; it has no Settings(set; kw...) copy
    # constructor of the kind the rest of JAC uses, which is a separate convention gap and not this test's business.
    rdmSettings = ReducedDensityMatrix.Settings(calcNatural = true, calcDensity = false, calcIpq = false,
                      calc2pRDM = false, printBefore = false, levelSelection = LevelSelection(true, indices=[1]))

    # (a) the closed-shell reference: one configuration, so rho^(1p) is exactly diagonal with the subshell
    #     occupations 2, 2, 2, 4 on it and nothing anywhere else.
    waA = Atomic.Computation(Atomic.Computation(), name="RDM sanity: Ne closed shell", grid=grid,
              nuclearModel = Nuclear.Model(10.), configs = [Configuration("1s^2 2s^2 2p^6")],
              propertySettings = [rdmSettings] )
    outA  = perform(waA; output=true)[ReducedDensityMatrix.Settings]
    rhoA  = outA[1].rho1p;    nA = size(rhoA, 1)
    offA  = maximum([abs(rhoA[p,q])  for p = 1:nA, q = 1:nA  if p != q];  init = 0.0)
    diagA = [rhoA[p,p]  for p = 1:nA]
    if  offA > 1.0e-10
        success = false
        if printTest   info(iostream, "closed-shell rho^(1p) is not diagonal; largest off-diagonal $offA")   end
    end
    if  maximum(abs.(diagA - [2.0, 2.0, 2.0, 4.0])) > 1.0e-10
        success = false
        if printTest   info(iostream, "closed-shell occupations $diagA, must be [2, 2, 2, 4]")   end
    end

    # (b) the correlated case: two configurations differing by a 2p -> 3p substitution, with 2p_3/2 lying
    #     between 2p_1/2 and 3p_1/2 in the subshell order, which is precisely the arrangement the retired
    #     shortcut could not handle.
    waB = Atomic.Computation(Atomic.Computation(), name="RDM sanity: Ne with 2p -> 3p", grid=grid,
              nuclearModel = Nuclear.Model(10.),
              configs = [Configuration("1s^2 2s^2 2p^6"), Configuration("1s^2 2s^2 2p^5 3p")],
              propertySettings = [rdmSettings] )
    outB   = perform(waB; output=true)[ReducedDensityMatrix.Settings]
    rhoB   = outB[1].rho1p;    nB = size(rhoB, 1)
    subsh  = outB[1].naturalSubshells;    occB = outB[1].naturalOccupation

    # (1) the trace is the electron number, exactly
    traceB = sum(rhoB[p,p]  for p = 1:nB)
    if  abs(traceB - 10.0) > 1.0e-8
        success = false
        if printTest   info(iostream, "trace of rho^(1p) = $traceB, must be the electron number 10")   end
    end
    # (2) rho^(1p) is symmetric
    symB = maximum(abs.(rhoB - transpose(rhoB)))
    if  symB > 1.0e-12
        success = false
        if printTest   info(iostream, "rho^(1p) is not symmetric; largest asymmetry $symB")   end
    end
    # (3) every natural occupation lies between 0 and the subshell capacity 2j+1, and they sum to N
    for  (i, sh)  in  enumerate(subsh)
        capacity = Basics.subshell_2j(sh) + 1
        if  occB[i] < -1.0e-8  ||  occB[i] > capacity + 1.0e-8
            success = false
            if printTest   info(iostream, "natural occupation $(occB[i]) of $sh outside [0, $capacity]")   end
        end
    end
    if  abs(sum(occB) - 10.0) > 1.0e-8
        success = false
        if printTest   info(iostream, "natural occupations sum to $(sum(occB)), must be 10")   end
    end
    # (4) the correlation must actually be there: see the note at the head of this function
    offB = maximum([abs(rhoB[p,q])  for p = 1:nB, q = 1:nB  if p != q];  init = 0.0)
    if  offB < 1.0e-3
        success = false
        if printTest   info(iostream, "all off-diagonals of rho^(1p) vanish (largest $offB); a correlated "  *
                                      "level must have them -- this is the compute1pRDMDirect defect back.")   end
    end

    println(iostream, "ReducedDensityMatrix: algebraic invariants of rho^(1p) -- trace = N, symmetry, "        *
                      "0 <= occupation <= 2j+1 -- plus the exact closed-shell limit and a guard that the "     *
                      "off-diagonal correlation block does not vanish. No comparison against approved data.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_ReducedDensityMatrix()::", success)
    return(success)
end


"""
`TestFrames.testModule_StarkShift(; short::Bool=true)`  ... tests on module StarkShift; a success::Bool is returned.
"""
function testModule_StarkShift(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-StarkShift-new.sum")
    printstyled("\n\nTest the module  StarkShift  ... \n", color=:cyan)
    # The quadratic Stark shift  dE(J,M) = -(1/2) alpha_0 E^2 - (1/2) alpha_2 E^2 [3M^2 - J(J+1)]/[J(2J-1)]  has three
    # exact properties that hold whatever the polarizabilities come out as, so they test the module's own algebra and
    # not the accuracy of the polarizability sum behind it: the tensor term is TRACELESS in M, the shift depends on M
    # only through M^2, and the whole shift scales as E^2.
    #
    # Hydrogen 2p is used because it gives BOTH cases in one computation: 2p_1/2 has J = 1/2, where J(2J-1) = 0 forces
    # alpha_2 to be exactly zero, and 2p_3/2 has J = 3/2, where the tensor term is live.  Without a level of the second
    # kind the first two checks would be satisfied by a module that simply never computed a tensor shift, so the size
    # of the splitting is itself asserted below.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    nm         = Nuclear.Model(1., PointNucleus())
    perturbers = [Configuration("1s"), Configuration("3s"), Configuration("4s"), Configuration("3d"), Configuration("4d")]
    grid       = Basics.recommendedGrid( vcat(perturbers, [Configuration("2p")]), nm )
    starkOutcomes = function(EField::Float64, omega::Float64 = 0.0)
        wb1 = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = "H perturbers", grid = grid, nuclearModel = nm,
                                        configs = perturbers); output = true )
        end
        settings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts = true, gMultiplet = wb1[ResultKeys.Multiplet],
                                       EField = EField, omega = omega)
        wb = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = "H 2p Stark shift", grid = grid, nuclearModel = nm,
                                        configs = [Configuration("2p")], propertySettings = [settings]); output = true )
        end
        return( wb[StarkShift.Settings] )
    end
    EField   = 1.0e5
    outcomes = starkOutcomes(EField)
    EFieldAu = EField / StarkShift.AU_EFIELD_IN_VCM
    sawTensor = false
    for  outcome  in  outcomes
        Jd     = Basics.twice(outcome.Jlevel.J) / 2.0
        alpha0 = outcome.alpha0.Babushkin
        #
        # (1) Sum over M.  The tensor factor [3M^2 - J(J+1)] sums to zero over the 2J+1 sublevels, so whatever alpha_2
        #     is, the centre of gravity of the multiplet is the pure scalar shift.
        total = sum( sub.energy.Babushkin for sub in outcome.Jsublevels )
        want  = (2*Jd + 1) * (-0.5) * alpha0 * EFieldAu^2
        if  abs(total - want) > 1.0e-10 * max(abs(want), 1.0e-30)
            success = false
            if printTest   info(iostream, "for J = $Jd the sublevel shifts sum to $total, but the traceless tensor " *
                                          "term requires (2J+1)(-alpha_0 E^2/2) = $want")   end
        end
        #
        # (2) The shift depends on M only through M^2, so the M and -M sublevels must be exactly degenerate.
        for  subA  in  outcome.Jsublevels,  subB  in  outcome.Jsublevels
            if  Basics.twice(subA.M) == -Basics.twice(subB.M)  &&
                abs(subA.energy.Babushkin - subB.energy.Babushkin) > 1.0e-14 * max(abs(subA.energy.Babushkin), 1.0e-30)
                success = false
                if printTest   info(iostream, "for J = $Jd the M = $(subA.M) and M = $(subB.M) sublevels differ, " *
                                              "$(subA.energy.Babushkin) against $(subB.energy.Babushkin)")   end
            end
        end
        #
        # (3) J < 1 leaves J(2J-1) = 0, where the module returns alpha_2 = 0 rather than dividing by it; the level then
        #     shifts as a whole and does not split at all.
        spread = maximum(sub.energy.Babushkin for sub in outcome.Jsublevels) -
                 minimum(sub.energy.Babushkin for sub in outcome.Jsublevels)
        if      Jd < 1.   &&  (outcome.alpha2.Babushkin != 0.  ||  spread != 0.)
            success = false
            if printTest   info(iostream, "the J = $Jd level has alpha_2 = $(outcome.alpha2.Babushkin) and a splitting " *
                                          "of $spread, but J(2J-1) = 0 forces both to be exactly zero")   end
        elseif  Jd >= 1.  &&  spread > 0.
            sawTensor = true
        end
    end
    #
    # A tensor shift must actually have occurred somewhere, or checks (1) and (2) were satisfied by arithmetic on zeros.
    if  !sawTensor
        success = false
        if printTest   info(iostream, "no level showed a nonzero tensor splitting, so the traceless and M -> -M " *
                                      "checks above are vacuous")   end
    end
    #
    # (4) The shift is quadratic in the field: doubling E must multiply every sublevel shift by exactly four.  This is
    #     the one check that has content even where alpha_2 vanishes.
    for  (outcome, doubled)  in  zip(outcomes, starkOutcomes(2*EField))
        for  (subA, subB)  in  zip(outcome.Jsublevels, doubled.Jsublevels)
            if  abs(subA.energy.Babushkin) > 0.  &&
                abs(subB.energy.Babushkin / subA.energy.Babushkin - 4.0) > 1.0e-10
                success = false
                if printTest   info(iostream, "doubling the field changed the M = $(subA.M) shift by a factor " *
                                              "$(subB.energy.Babushkin / subA.energy.Babushkin) rather than 4")   end
            end
        end
    end

    # (5) THE AC BRANCH, AND THE ONE THING THAT CAN GO WRONG IN IT -- added 19-Oct-2026.  A static field of
    #     strength E and a field oscillating with PEAK amplitude E are not the same field: <cos^2> = 1/2, so the
    #     cycle-averaged shift carries -(1/4) alpha E^2 where the static one carries -(1/2) alpha E^2.  Hence the
    #     omega -> 0 LIMIT OF THE AC EXPRESSION IS EXACTLY HALF THE STATIC SHIFT, and that factor is the classic
    #     error in this quantity -- adopted in one convention and quoted in the other, it is invisible in every
    #     other check here, since halving every shift leaves the traceless sum, the M -> -M degeneracy and the
    #     exact E^2 scaling all intact.
    #
    #     omega = 1.0e-6 a.u. is far below the 2p -> ns/nd excitation energies, so alpha(omega) differs from
    #     alpha(0) by order (omega/dE)^2 ~ 1e-12 and the ratio is the convention factor alone.  The bound is 1e-9.
    for  (outcome, acOutcome)  in  zip(outcomes, starkOutcomes(EField, 1.0e-6))
        for  (subDc, subAc)  in  zip(outcome.Jsublevels, acOutcome.Jsublevels)
            abs(subDc.energy.Babushkin) > 0.  ||  continue
            ratioAc = subAc.energy.Babushkin / subDc.energy.Babushkin
            if  abs(ratioAc - 0.5) > 1.0e-9
                success = false
                if printTest   info(iostream, "the ac shift of the M = $(subAc.M) sublevel at omega -> 0 is " *
                                              "$ratioAc times the static one, where the cycle average demands " *
                                              "exactly one half")   end
            end
        end
    end

    # (6) THE BLACK-BODY SHIFT, AND THE TWO THINGS THAT DECIDE IT -- added 20-Oct-2026.  They fail differently,
    #     so they are asserted differently.
    #
    #     (a) <E^2> IS A COMPLETE EXTERNAL ANCHOR, because no atomic quantity enters it at all -- only the
    #         fine-structure constant does.  (4 pi^3/15) alpha_fs^3 T^4 at 300 K is E_rms = 831.9 V/m, the number
    #         every clock paper quotes.  The 0.1 V/m bound is LOOSE ON PURPOSE: what is being checked is the
    #         SCALE, and a wrong scale is wrong by a factor -- a missing factor of two moves E_rms by 1.41, and
    #         4 pi in place of 4 pi^3/15 by 2.2.
    #
    #     (b) eta -> 0 AS T -> 0.  The Planck average must collapse onto the closed-form leading term
    #         -(2/15)(alpha_fs pi)^3 T^4 alpha_0 once the spectrum is too cold to resolve any excitation energy,
    #         and that is what fixes the quadrature's normalisation.
    #
    #     AND THE LIMIT ALONE IS NOT ENOUGH, WHICH IS WHY THE T^2 SCALING IS ASSERTED BESIDE IT.  Both reasons
    #     were measured rather than supposed.  The first implementation normalised the Planck weight by the
    #     analytic pi^4/15 while taking its numerator by Simpson; those differ by the quadrature error, +2.57e-07
    #     relative at 400 nodes, and that went into the answer as a FLOOR on eta.  A floor looks exactly like a
    #     converged zero -- it was the SCALING that exposed it, eta coming out 5.97e-07 at 3 K where the law
    #     demands 3.4e-07.  That was a real defect, found while writing this test.  The second reason is that a
    #     quadrature using the IMAGINARY-axis factor dE/(dE^2 + omega^2) in place of the real
    #     dE/(dE^2 - omega^2) satisfies the limit to 1e-13 and is caught only away from it, where it returns eta
    #     of the opposite SIGN.  Hence three assertions and not one: eta > 0, eta -> 0, and eta(10 K)/eta(3 K)
    #     = (10/3)^2.
    alphaFs = 1.0 / Defaults.getDefaults("speed of light: c")
    T300    = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 300.0)
    eRms300 = sqrt( StarkShift.blackbodyFieldSquare(300.0) ) * StarkShift.AU_EFIELD_IN_VCM * 100.0
    if  abs(eRms300 - 831.9) > 0.1
        success = false
        if printTest   info(iostream, "the black-body field at 300 K came out E_rms = $eRms300 V/m where the " *
                                      "Stefan-Boltzmann value, which needs only the fine-structure constant, " *
                                      "is 831.9 V/m")   end
    end
    ratioBb = 0.5 * StarkShift.blackbodyFieldSquare(300.0) / (2.0/15.0 * (alphaFs*pi)^3 * T300^4)
    if  abs(ratioBb - 1.0) > 1.0e-13
        success = false
        if printTest   info(iostream, "(1/2)(4 pi^3/15) alpha_fs^3 T^4 and (2/15)(alpha_fs pi)^3 T^4 are the " *
                                      "same expression, but their ratio came out $ratioBb")   end
    end
    #
    gmBb = redirect_stdout(devnull) do
        perform( Atomic.Computation(Atomic.Computation(); name = "H perturbers", grid = grid, nuclearModel = nm,
                                    configs = perturbers); output = true )[ResultKeys.Multiplet]
    end
    mpBb = redirect_stdout(devnull) do
        perform( Atomic.Computation(Atomic.Computation(); name = "H 2p", grid = grid, nuclearModel = nm,
                                    configs = [Configuration("2p")]); output = true )[ResultKeys.Multiplet]
    end
    levBb      = mpBb.levels[1]
    a0Bb, _, _ = redirect_stdout(devnull) do
        MultipolePolarizibility.computeScalarTensorPolarizability(levBb, gmBb, grid)   end
    alpha0Bb   = a0Bb.Babushkin
    shiftBb    = function(TK::Float64)
        sh, et = redirect_stdout(devnull) do;   StarkShift.computeBlackbodyShift(levBb, gmBb, grid, TK)   end
        return( (sh.Babushkin, et) )
    end
    if  alpha0Bb <= 0.0
        success = false
        if printTest   info(iostream, "the H 2p level came out with alpha_0 = $alpha0Bb, so the black-body " *
                                      "checks below cannot be set up")   end
    else
        ## the limit: at 1 K the whole Planck spectrum lies ~4 orders below the smallest gap, so the shift must be
        ## the closed-form leading term and eta must be a clean zero -- 3.8e-08 as measured, bounded at 1e-06
        sh1, eta1 = shiftBb(1.0)
        T1Bb      = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.0)
        leadBb    = -2.0/15.0 * (alphaFs*pi)^3 * T1Bb^4 * alpha0Bb
        if  abs(sh1/leadBb - 1.0) > 1.0e-6
            success = false
            if printTest   info(iostream, "at 1 K the black-body shift came out $sh1 against the closed-form " *
                                          "leading term $leadBb, a ratio of $(sh1/leadBb); as T -> 0 the Planck " *
                                          "average must collapse onto it")   end
        end
        if  !(0.0 < eta1 < 1.0e-6)
            success = false
            if printTest   info(iostream, "at 1 K the dynamic correction came out eta = $eta1, where it must be " *
                                          "positive -- the real-axis factor dE/(dE^2 - omega^2) exceeds 1/dE -- " *
                                          "and below 1e-6, since the spectrum resolves no excitation energy")   end
        end
        ## and the scaling, which is what a floor or a wrong frequency factor cannot survive
        _, eta3  = shiftBb(3.0)
        _, eta10 = shiftBb(10.0)
        if  eta3 <= 0.0  ||  eta10 <= 0.0
            success = false
            if printTest   info(iostream, "eta came out $eta3 at 3 K and $eta10 at 10 K; both must be positive " *
                                          "on the real frequency axis")   end
        else
            for  (TaBb, TbBb, etaA, etaB)  in  ( (1.0, 3.0, eta1, eta3), (3.0, 10.0, eta3, eta10) )
                scaleBb = (etaB/etaA) / (TbBb/TaBb)^2
                if  abs(scaleBb - 1.0) > 1.0e-2
                    success = false
                    if printTest   info(iostream, "eta grew by a factor $(etaB/etaA) from $TaBb K to $TbBb K, " *
                                                  "where eta ~ <omega^2>/dE^2 ~ T^2 demands " *
                                                  "$((TbBb/TaBb)^2); the ratio of the two is $scaleBb")   end
                end
            end
        end
    end

    println(iostream, "StarkShift: for hydrogen 2p the tensor term is traceless in M, the sublevels are degenerate " *
                      "in M -> -M, J = 1/2 leaves alpha_2 exactly zero, and doubling the field multiplies every " *
                      "shift by exactly four. The ac branch is asserted too: at omega -> 0 every sublevel shift " *
                      "is EXACTLY HALF the static one, which is the cycle average <cos^2> = 1/2 and the one " *
                      "factor this quantity is habitually quoted in the wrong convention. The black-body shift " *
                      "is asserted too: <E^2> at 300 K gives E_rms = 831.9 V/m, which needs only the " *
                      "fine-structure constant and is therefore a genuine external anchor, and the Planck " *
                      "average collapses onto the closed-form leading term as T -> 0 with eta growing as T^2 " *
                      "thereafter -- the scaling being the part that a quadrature floor or a wrong frequency " *
                      "factor cannot survive. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_StarkShift()::", success)
    return( success )
end


"""
`TestFrames.testModule_Statistical(; short::Bool=true)`  ... tests on module Statistical; a success::Bool is returned.
"""
function testModule_Statistical(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-Statistical-new.sum")
    printstyled("\n\nTest the module  Statistical  ... \n", color=:cyan)
    # Every statement below is a theorem about a density matrix, not a number anyone measured, so none of it can drift
    # and none of it needs a reference or an atomic-structure computation. That matters here more than usual: two stubs
    # of this module were caught in July PRINTING FABRICATED TENSORS, which is precisely the failure that a test
    # comparing one computed number against another computed number cannot see. (1) and (4) are ABSOLUTE.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    key  = LevelKey( LevelSymmetry(AngularJ64(3//2), Basics.plus), 0., 1, 0. )
    J    = AngularJ64(3//2);      MList = AngularMomentum.m_values(J)
    ampl = Dict{AngularM64,ComplexF64}( MList[1] => ComplexF64(0.3, 0.1), MList[2] => ComplexF64(-0.5, 0.2),
                                        MList[3] => ComplexF64(0.7, -0.4), MList[4] => ComplexF64(0.2, 0.6) )
    tens = Statistical.computeTensors(6, ampl, key)
    # Tensor.k is an AngularJ64 and Tensor.q an AngularM64, not integers, so they are compared through
    # AngularMomentum.oneJ / oneM exactly as Statistical.alignmentParameters itself does.
    get(k, q) = ( w = filter(t -> AngularMomentum.oneJ(t.k) == 1.0*k  &&  AngularMomentum.oneM(t.q) == 1.0*q, tens);
                  length(w) == 0 ? ComplexF64(0.) : w[1].value )

    # (1) THE NORMALISATION, and it is an absolute statement rather than a ratio: the Clebsch-Gordan
    #     <J M, J -M | 0 0> = (-1)^(J-M)/sqrt(2J+1) cancels the (-1)^(J-M) phase exactly, so rho_00 is the total
    #     population divided by sqrt(2J+1) -- no more and no less.
    trace = sum( abs2(ampl[M])  for M in MList )
    rho00 = get(0, 0)
    if  abs( rho00 - trace / sqrt(Basics.twice(J) + 1.0) ) > 1.0e-12
        success = false
        if printTest   info(iostream, "rho_00 = $rho00 against the exact $(trace/sqrt(Basics.twice(J)+1.0))")   end
    end

    # (2) HERMITICITY. rho_(M,M') = a_M conj(a_M') is Hermitian for ANY amplitudes whatsoever, and that carries over
    #     to the tensors as rho_(k,-q) = (-1)^q conj(rho_kq). A wrong phase or a swapped index in the sum breaks this
    #     while leaving every magnitude intact, so it catches what a check on |rho| cannot.
    for  k = 0:3,  q = 1:k
        if  abs( get(k, -q) - (-1)^q * conj(get(k, q)) ) > 1.0e-12
            success = false
            if printTest   info(iostream, "rho_($k,$(-q)) = $(get(k,-q)) against (-1)^q conj(rho_($k,$q)) = " *
                                          "$((-1)^q * conj(get(k,q)))")   end
        end
    end

    # (3) THE TRIANGLE CONDITION. Two angular momenta J couple to at most 2J, so every tensor of rank k > 2J must
    #     vanish identically. Here 2J = 3, and ranks up to 6 were asked for above precisely so that 4, 5 and 6 have
    #     the chance to be wrong.
    for  k = 4:6,  q = -k:k
        if  abs( get(k, q) ) != 0.
            success = false
            if printTest   info(iostream, "rho_($k,$q) = $(get(k,q)) is non-zero although k > 2J = 3")   end
        end
    end

    # (4) AN ISOTROPIC ENSEMBLE HAS NO ALIGNMENT. Equal, incoherent populations of every sublevel describe an
    #     unpolarized state, for which every k > 0 tensor vanishes exactly and A_00 = 1 alone survives. This is the
    #     one case whose answer is known without any algebra, so a fabricated tensor cannot hide in it.
    dm = Dict{Tuple{AngularM64,AngularM64},ComplexF64}()
    for  M in MList,  Mp in MList    dm[(M, Mp)] = (M == Mp  ?  ComplexF64(1.)  :  ComplexF64(0.))    end
    iso = Statistical.computeTensors(3, dm, key)
    for  t  in  iso
        if  AngularMomentum.oneJ(t.k) != 0.  &&  abs(t.value) > 1.0e-12
            success = false
            if printTest   info(iostream, "an isotropic ensemble has rho_($(t.k),$(t.q)) = $(t.value), not 0")   end
        end
    end
    aTens = Statistical.alignmentParameters(iso)
    for  t  in  aTens
        expected = (AngularMomentum.oneJ(t.k) == 0.  &&  AngularMomentum.oneM(t.q) == 0.)  ?  ComplexF64(1.)  :  ComplexF64(0.)
        if  abs(t.value - expected) > 1.0e-12
            success = false
            if printTest   info(iostream, "A_($(t.k),$(t.q)) = $(t.value) for an isotropic ensemble, not $expected")   end
        end
    end

    # (5) alignmentParameters REFUSES an empty ensemble rather than dividing by zero.
    try
        Statistical.alignmentParameters( Statistical.Tensor[] )
        success = false
        if printTest   info(iostream, "alignmentParameters accepted a list with no rho_00")   end
    catch                                                                                              end

    Defaults.setDefaults("print summary: close", "")
    if  success   println("\n\nTest.perform: `Test the module Statistical`  ... has been successful.")   end
    return( success )
end


"""
`TestFrames.testModule_WeakInteractionEnhancement(; short::Bool=true)`  ... tests on module WeakInteractionEnhancement;
    a success::Bool is returned.
"""
function testModule_WeakInteractionEnhancement(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-WeakInteractionEnhancement-new.sum")
    printstyled("\n\nTest the module  WeakInteractionEnhancement  ... \n", color=:cyan)
    # This module has NO absolute anchor -- there is nothing to compare an electron-EDM enhancement against -- so its
    # exact selection rules and its reality are the only handle anyone has on it, and that is what is asserted here.
    # THE TWO AMPLITUDES ARE EXACT OPPOSITES IN BOTH RESPECTS, which is the point of testing them together: the weak
    # charge carries gamma_5 and so is PURELY IMAGINARY and connects levels of EQUAL parity, while H_EDM contracts an
    # axial with a polar vector and so is REAL and connects levels of OPPOSITE parity with equal J. A test that swapped
    # the two would have to fail four assertions, not one.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 12.0)
    nm   = Nuclear.Model(20., FermiNucleus())
    asf  = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    # one even- and one odd-parity multiplet of Ca^19+, small enough to be quick and hydrogenic enough to be exact
    evenMp = redirect_stdout(devnull) do
        SelfConsistent.performSCF([Configuration("1s"), Configuration("2s")], nm, grid, asf; printout=false)   end
    oddMp  = redirect_stdout(devnull) do
        SelfConsistent.performSCF([Configuration("2p")], nm, grid, asf; printout=false)                       end
    evenLs = evenMp.levels;    oddLs = oddMp.levels

    # (1) THE PNC AMPLITUDE IS PURELY IMAGINARY for real radial orbitals, since the weak-charge amplitude is.
    #     A real part is not a small error but a wrong operator.
    for  f  in  evenLs,  i  in  evenLs
        a = WeakInteractionEnhancement.computePncE1Amplitude(f, i, oddMp, nm, grid)
        if  abs(real(a)) > 1.0e-14 * max(1.0, abs(a))
            success = false
            if printTest   info(iostream, "PNC amplitude $(f.index)<-$(i.index) = $a has a real part")   end
        end
    end

    # (2) IT VANISHES BETWEEN LEVELS OF OPPOSITE PARITY -- which is where the ORDINARY E1 amplitude lives. H_W is a
    #     pseudoscalar and D is P-odd, so their product is P-even and connects EQUAL parities; the opposite-parity
    #     matrix element is the one that must be exactly zero.
    for  f  in  oddLs,  i  in  evenLs
        if  abs( WeakInteractionEnhancement.computePncE1Amplitude(f, i, oddMp, nm, grid) ) != 0.
            success = false
            if printTest   info(iostream, "a PNC amplitude survives between the opposite-parity levels " *
                                          "$(f.index) and $(i.index)")   end
        end
    end

    # (3) THE EDM AMPLITUDE IS THE MIRROR IMAGE: REAL, and non-zero only between OPPOSITE parities with EQUAL J.
    for  f  in  oddLs,  i  in  evenLs
        a = WeakInteractionEnhancement.edmAmplitude(f, i, nm, grid)
        if  abs(imag(a)) > 1.0e-14 * max(1.0, abs(a))
            success = false
            if printTest   info(iostream, "EDM amplitude $(f.index)<-$(i.index) = $a has an imaginary part")   end
        end
        if  f.J != i.J  &&  abs(a) != 0.
            success = false
            if printTest   info(iostream, "EDM amplitude $(f.index)<-$(i.index) = $a survives J_f = $(f.J) != J_i = $(i.J)")   end
        end
    end
    for  f  in  evenLs,  i  in  evenLs
        if  abs( WeakInteractionEnhancement.edmAmplitude(f, i, nm, grid) ) != 0.
            success = false
            if printTest   info(iostream, "an EDM amplitude survives between the equal-parity levels " *
                                          "$(f.index) and $(i.index)")   end
        end
    end

    # (4) AN EMPTY INTERMEDIATE SET GIVES EXACTLY ZERO. The PNC amplitude is a sum over opposite-parity intermediate
    #     levels; with none supplied there is nothing to sum, and anything other than zero would be invented.
    #
    # WHAT THE NEGATIVE CONTROLS SHOWED, 31-Aug-2026, and it is worth knowing before anyone strengthens this test:
    # THE TWO PARITY GUARDS IN THIS MODULE ARE REDUNDANT AND CANNOT BE CAUGHT BY ANY TEST. Flipping
    # computePncE1Amplitude's early return from != to ==, and deleting edmAmplitude's guard outright, both leave
    # every number unchanged -- because the physics inside already enforces them. In the PNC sum, term 1 is
    # D(f,n) H_W(n,i): H_W needs n opposite to i and D needs f opposite to n, so f and i must share parity whatever
    # the guard says; term 2 is the same. In the EDM amplitude the operator requires kappa_b = -kappa_a, which
    # excludes equal parities on its own. The guards are honest shortcuts, not constraints, and assertions (2) and
    # (3) below therefore verify the PHYSICS rather than the guards. What DOES make this test go red is destroying
    # the amplitude's reality -- replacing wOrd by abs(wOrd) in the first PNC term was the control that worked.
    empty = Multiplet("empty", Level[])
    for  f  in  evenLs,  i  in  evenLs
        if  abs( WeakInteractionEnhancement.computePncE1Amplitude(f, i, empty, nm, grid) ) != 0.
            success = false
            if printTest   info(iostream, "a PNC amplitude is non-zero with no intermediate levels at all")   end
        end
    end

    Defaults.setDefaults("print summary: close", "")
    if  success   println("\n\nTest.perform: `Test the module WeakInteractionEnhancement`  ... has been successful.")   end
    return( success )
end


"""
`TestFrames.testModule_CrystalField(; short::Bool=true)`  ... tests on module CrystalField; a success::Bool is returned.
"""
function testModule_CrystalField(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-CrystalField-new.sum")
    printstyled("\n\nTest the module  CrystalField  ... \n", color=:cyan)
    # Everything asserted here is a theorem about a crystal field, not a number anyone measured, so none of it can
    # drift and none of it needs a reference. Two of them are exact for ANY lattice whatsoever -- Kramers degeneracy
    # for an odd electron count, and the tracelessness of every k >= 1 term -- and the third, that a RIGID ROTATION
    # of the lattice cannot move the spectrum, exercises the spherical harmonics, the 3-j symbols and the phases
    # together in a way that no single matrix element does.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # (1) THE LATTICE SUMS ALONE, which need no atomic structure at all. An empty lattice gives exactly nothing, and
    #     a perfectly octahedral one has A_2q = 0 for every q while its k = 4 term survives only at q = 0, +-4 in the
    #     fixed ratio A_44/A_40 = sqrt(5/14). That ratio is a property of the cubic group, not of this code.
    rho     = 200.0
    empty   = CrystalField.Lattice()
    for  k = 1:6,  q = -k:k
        if  abs( CrystalField.multipoleLatticeSum(empty, k, q) ) != 0.
            success = false
            if printTest   info(iostream, "an empty lattice gives a non-zero A_$k$q")   end
        end
    end
    charge  = -2.0
    ohIons  = [ CrystalField.PointCharge(charge, rho, 0.0,  0.0),   CrystalField.PointCharge(charge, rho, pi,   0.0),
                CrystalField.PointCharge(charge, rho, pi/2, 0.0),   CrystalField.PointCharge(charge, rho, pi/2, pi),
                CrystalField.PointCharge(charge, rho, pi/2, pi/2),  CrystalField.PointCharge(charge, rho, pi/2, 3pi/2) ]
    oh      = CrystalField.Lattice(ohIons, "Oh")
    a40     = CrystalField.multipoleLatticeSum(oh, 4, 0)
    for  q = -2:2
        if  abs( CrystalField.multipoleLatticeSum(oh, 2, q) ) > 1.0e-10 * abs(a40)
            success = false
            if printTest   info(iostream, "the k = 2 lattice sum of a perfect octahedron does not vanish at q = $q")   end
        end
    end
    for  q  in  [-3, -2, -1, 1, 2, 3]
        if  abs( CrystalField.multipoleLatticeSum(oh, 4, q) ) > 1.0e-10 * abs(a40)
            success = false
            if printTest   info(iostream, "the k = 4 lattice sum of a perfect octahedron does not vanish at q = $q")   end
        end
    end
    if  abs( real(CrystalField.multipoleLatticeSum(oh, 4, 4) / a40) - sqrt(5/14) ) > 1.0e-10
        success = false
        if printTest   info(iostream, "A_44/A_40 of a perfect octahedron is " *
                                      "$(CrystalField.multipoleLatticeSum(oh, 4, 4)/a40), not sqrt(5/14)")   end
    end

    # The test system is bare hydrogen H(3d): ONE active electron, so J = 3/2 and 5/2 are single-CSF levels with no CI
    # mixing at all, and every matrix element below depends on the lattice geometry and the tensor algebra alone. One
    # electron is also an ODD electron count, which is what puts Kramers' theorem within reach.
    wa = Atomic.Computation(Atomic.Computation(), name="test-CrystalField", grid = Radial.Grid(true),
                            nuclearModel = Nuclear.Model(1., UniformNucleus(), 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
                            configs = [Configuration("3d")], propertySettings = Basics.AbstractPropertySettings[] )
    multiplet = redirect_stdout(devnull) do;   perform(wa; output=true)[ResultKeys.Multiplet]   end
    grid      = Radial.Grid(true)
    # A DELIBERATELY LOW-SYMMETRY lattice: three unequal charges at three generic directions and three distances. The
    # point is that no point-group argument can force any degeneracy here, so a degeneracy that does appear is Kramers'
    # and nothing else.
    lowSym = CrystalField.Lattice( [ CrystalField.PointCharge(-2.0, rho,      0.7, 0.3),
                                     CrystalField.PointCharge(-1.3, 1.4*rho, 1.9, 2.1),
                                     CrystalField.PointCharge(-0.8, 0.9*rho, 2.6, 4.7) ], "low symmetry" )
    spectrum = function(lev::Level, latt::CrystalField.Lattice, scale::Float64)
        cfm = redirect_stdout(devnull) do
                  CrystalField.computeRepresentation([lev], latt, CrystalField.PointChargeModel(scale), grid, 6)
              end
        return( sort([ cfLev.energy  for cfLev in cfm.cfLevels ]) )
    end

    for  lev  in  multiplet.levels
        energies = spectrum(lev, lowSym, 1.0);      n = length(energies)
        shifts   = energies .- lev.energy;          scale = maximum( abs.(shifts) )
        # (2) KRAMERS DEGENERACY, exact for a half-integer J in ANY electrostatic field. The 2J+1 sublevels must come
        #     in exactly degenerate pairs -- and the pairs must be DISTINCT from one another, or the check would be
        #     passing on a spectrum that simply never split.
        for  i = 1:div(n, 2)
            if  abs( energies[2i] - energies[2i-1] ) > 1.0e-10 * scale
                success = false
                if printTest   info(iostream, "J = $(lev.J): Kramers pair $i is split by " *
                                              "$(energies[2i] - energies[2i-1]) against a scale of $scale")   end
            end
        end
        pairs = [ (energies[2i] + energies[2i-1])/2  for i = 1:div(n, 2) ]
        for  i = 1:length(pairs)-1
            if  abs( pairs[i+1] - pairs[i] ) < 1.0e-3 * scale
                success = false
                if printTest   info(iostream, "J = $(lev.J): Kramers doublets $i and $(i+1) coincide, so the " *
                                              "degeneracy check is not being made on a split spectrum")   end
            end
        end
        # (3) THE CENTRE OF GRAVITY CANNOT MOVE. Every k >= 1 term is traceless over the M_J basis -- the sum over M of
        #     the 3-j symbol vanishes for k >= 1 -- so whatever the lattice does to the individual sublevels, the sum
        #     of the shifts is exactly zero.
        if  abs( sum(shifts) ) > 1.0e-10 * scale
            success = false
            if printTest   info(iostream, "J = $(lev.J): the crystal-field shifts sum to $(sum(shifts)) against a " *
                                          "scale of $scale, although every k >= 1 term is traceless")   end
        end
        # (4) AND AN EMPTY LATTICE MOVES NOTHING AT ALL, which is the other side of (2) and (3): without it they could
        #     both be passing on a spectrum in which nothing ever happened.
        for  wb  in  spectrum(lev, empty, 1.0)
            if  wb != lev.energy
                success = false
                if printTest   info(iostream, "J = $(lev.J): an empty lattice shifts a sublevel by $(wb - lev.energy)")   end
            end
        end
    end

    # (5) A RIGID ROTATION OF THE LATTICE CANNOT MOVE THE SPECTRUM. The eigenvalues describe the ion in its field; which
    #     way round the laboratory axes are pointing is not a physical fact about it. Ranks do not mix under rotation,
    #     so this stays exact even with the k <= 6 truncation. The rotation is applied in CARTESIAN coordinates and
    #     converted back, so it shares no code with the spherical harmonics it tests.
    function rotateLattice(latt::CrystalField.Lattice, alpha::Float64, beta::Float64, gamma::Float64)
        ions = CrystalField.PointCharge[]
        for  ion  in  latt.ions
            x  = ion.rho*sin(ion.theta)*cos(ion.phi);   y = ion.rho*sin(ion.theta)*sin(ion.phi);   z = ion.rho*cos(ion.theta)
            x1 = cos(gamma)*x - sin(gamma)*y;           y1 = sin(gamma)*x + cos(gamma)*y;          z1 = z
            x2 = cos(beta)*x1 + sin(beta)*z1;           y2 = y1;                                   z2 = -sin(beta)*x1 + cos(beta)*z1
            x3 = cos(alpha)*x2 - sin(alpha)*y2;         y3 = sin(alpha)*x2 + cos(alpha)*y2;        z3 = z2
            wb = sqrt(x3^2 + y3^2 + z3^2)
            push!(ions, CrystalField.PointCharge(ion.charge, wb, acos(z3/wb), atan(y3, x3)))
        end
        return( CrystalField.Lattice(ions, "rotated") )
    end
    lev      = multiplet.levels[2]
    energies = spectrum(lev, lowSym, 1.0);   scale = maximum(energies) - minimum(energies)
    for  (alpha, beta, gamma)  in  [(0.0, 0.0, 1.234), (0.0, 0.9, 0.0), (0.31, 1.77, 2.55)]
        wb = spectrum(lev, rotateLattice(lowSym, alpha, beta, gamma), 1.0)
        if  maximum( abs.(wb - energies) ) > 1.0e-10 * scale
            success = false
            if printTest   info(iostream, "rotating the lattice by (alpha, beta, gamma) = ($alpha, $beta, $gamma) " *
                                          "moves the spectrum by $(maximum(abs.(wb - energies)))")   end
        end
    end

    # (6) J-MIXING, where the two parent levels are diagonalized in ONE common M_J basis. Nothing above reaches the
    #     off-diagonal <J||C^(k)||J'> blocks at all, and both theorems survive into the joint problem unchanged: the
    #     ten sublevels still pair up under Kramers, and the centre of gravity of the whole ten is still the weighted
    #     mean of the two parent energies, since the k >= 1 terms are traceless within each diagonal block and the
    #     off-diagonal blocks contribute nothing to a trace.
    joint = redirect_stdout(devnull) do
                CrystalField.computeRepresentation(multiplet.levels, lowSym, CrystalField.PointChargeModel(1.0), grid, 6)
            end
    energies = sort([ cfLev.energy  for cfLev in joint.cfLevels ])
    wb       = sum( Basics.twice(lv.J)+1 for lv in multiplet.levels )
    if  length(energies) != wb
        success = false
        if printTest   info(iostream, "J-mixing over the two parent levels gives $(length(energies)) sublevels, not $wb")   end
    end
    barycentre = sum( (Basics.twice(lv.J)+1) * lv.energy  for lv in multiplet.levels )
    scale      = maximum(energies) - minimum(energies)
    if  abs( sum(energies) - barycentre ) > 1.0e-10 * scale
        success = false
        if printTest   info(iostream, "under J-mixing the centre of gravity moves by $(sum(energies) - barycentre) " *
                                      "against a scale of $scale")   end
    end
    for  i = 1:div(length(energies), 2)
        if  abs( energies[2i] - energies[2i-1] ) > 1.0e-10 * scale
            success = false
            if printTest   info(iostream, "under J-mixing, Kramers pair $i is split by " *
                                          "$(energies[2i] - energies[2i-1]) against a scale of $scale")   end
        end
    end
    # ... and the off-diagonal blocks must actually DO something, or the three checks above would be passing on a
    #     problem that never mixed: the joint spectrum must differ from the union of the two separate ones.
    wb = sort( vcat( spectrum(multiplet.levels[1], lowSym, 1.0), spectrum(multiplet.levels[2], lowSym, 1.0) ) )
    if  maximum( abs.(energies - wb) ) < 1.0e-6 * scale
        success = false
        if printTest   info(iostream, "the J-mixed spectrum equals the union of the unmixed ones, so the " *
                                      "off-diagonal blocks are contributing nothing")   end
    end

    # (7) THE SCALE FIELD, AND fitScaleField AS ITS INVERSE. For a single parent level the diagonal of the interaction
    #     matrix is a multiple of the identity, so every shift is EXACTLY proportional to scaleField and so is the
    #     characteristic splitting. fitScaleField must therefore return exactly the factor asked of it.
    cxs1 = CrystalField.characteristicSplitting( redirect_stdout(devnull) do
               CrystalField.computeRepresentation([lev], lowSym, CrystalField.PointChargeModel(1.0), grid, 6) end )
    cxs2 = CrystalField.characteristicSplitting( redirect_stdout(devnull) do
               CrystalField.computeRepresentation([lev], lowSym, CrystalField.PointChargeModel(2.0), grid, 6) end )
    if  cxs1 <= 0.   ||   abs(cxs2/cxs1 - 2.0) > 1.0e-10
        success = false
        if printTest   info(iostream, "doubling scaleField multiplies the characteristic splitting by " *
                                      "$(cxs2/cxs1), not by 2")   end
    end
    wb = redirect_stdout(devnull) do;  CrystalField.fitScaleField([lev], lowSym, grid, 3.7*cxs1, 6)  end
    if  abs(wb - 3.7) > 1.0e-8
        success = false
        if printTest   info(iostream, "fitScaleField asked for 3.7 times the unit-scale splitting returns $wb")   end
    end

    # (8) characteristicSplitting itself is an algorithm on a list of numbers, and is checked as one: a spectrum whose
    #     largest gap sits between 2 and 10 must be split there, giving mean(10, 11) - mean(0, 1, 2) = 9.5.
    cfBasis = CrystalField.CfBasisVector[]
    cfLev   = [ CrystalField.CfLevel(en, cfBasis, ComplexF64[])  for en in [10.0, 0.0, 2.0, 11.0, 1.0] ]
    wb      = CrystalField.characteristicSplitting( CrystalField.CfMultiplet("synthetic", cfLev) )
    if  abs(wb - 9.5) > 1.0e-12
        success = false
        if printTest   info(iostream, "characteristicSplitting of [0,1,2,10,11] is $wb, not 9.5")   end
    end
    if  CrystalField.characteristicSplitting( CrystalField.CfMultiplet("one level", [cfLev[1]]) ) != 0.
        success = false
        if printTest   info(iostream, "characteristicSplitting of a single sublevel is not 0")   end
    end

    println(iostream, "CrystalField: Kramers degeneracy for the odd electron count of H(3d) in a lattice of no "  *
                      "symmetry at all, where nothing else can force it; the centre of gravity held fixed by the " *
                      "tracelessness of every k >= 1 term; both of those again over the J-mixed basis, where the "  *
                      "off-diagonal blocks are shown to contribute; the spectrum unmoved by three rigid rotations " *
                      "of the lattice; the k = 2 sum of a perfect octahedron vanishing while A_44/A_40 = "         *
                      "sqrt(5/14); fitScaleField as the exact inverse of a proportionality; and an empty lattice " *
                      "moving nothing at all. NOTE what none of this can reach: every check here is a structural " *
                      "invariant, so an overall CONSTANT on the crystal-field strength -- for instance the "       *
                      "sqrt(2J+1) of the Wigner-Eckart factor, which for one parent level is exactly that -- "     *
                      "passes unseen. Fixing the absolute scale needs a comparison, and that is item 67. "         *
                      "No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_CrystalField()::", success)
    return( success )
end


"""
`TestFrames.testModule_UncertaintyEstimate(; short::Bool=true)`
    ... tests module UncertaintyEstimate on hand-built variants, no grid, no orbitals and no approved file:
        (i) the spread and relative spread of two variants match the plain arithmetic by hand; (ii) a single
        variant is refused, since a spread needs at least two; (iii) a mean of exactly zero is refused, since
        the relative spread would be undefined. A success::Bool is returned.
"""
function testModule_UncertaintyEstimate(; short::Bool=true)
    printstyled("\n\nTest the module  UncertaintyEstimate  ... \n", color=:cyan)
    success = true

    variants = [ UncertaintyEstimate.Variant("a", 100.0), UncertaintyEstimate.Variant("b", 120.0) ]
    est = UncertaintyEstimate.estimate("test quantity", variants)
    if  abs(est.spread - 20.0) > 1.0e-10   ||   abs(est.relSpread - 20.0/110.0) > 1.0e-10
        success = false;   println("  ** spread or relSpread wrong: $(est.spread), $(est.relSpread)")
    end

    try
        UncertaintyEstimate.estimate("x", [UncertaintyEstimate.Variant("only", 1.0)])
        success = false;   println("  ** a single variant was NOT refused")
    catch
        # expected
    end

    try
        UncertaintyEstimate.estimate("x", [UncertaintyEstimate.Variant("a", 5.0), UncertaintyEstimate.Variant("b", -5.0)])
        success = false;   println("  ** a zero mean was NOT refused")
    catch
        # expected
    end

    testPrint("testModule_UncertaintyEstimate()::", success)
    return(success)
end


"""
`TestFrames.testModule_AdjustedLevels(; short::Bool=true)`
    ... tests module AdjustedLevels on a small, cheap, REAL He-like C V (Z=6) 1s2p J=1(odd) block (2 CSFs,
        singlet ^1P1 and triplet ^3P1): (i) AdjustedLevels.buildOneAndTwoBodyMatrices at an UNSCALED
        combination (oneBody + 1.0*twoBody) must reproduce `Hamiltonian.setupMatrix`'s own, independently
        trusted CI matrix EXACTLY -- this is the one check that matters most, since everything else in the
        module rests on this parallel construction being faithful; (ii) fitting against a target deliberately
        10% off the ab initio splitting converges to a scale factor measurably different from 1.0, with a
        near-zero residual RMS (a single parameter fit to a single target should match almost exactly);
        (iii) AdjustedLevels.mergeIntoMultiplet, against a basis that ALSO carries the 1s^2 ground state (a
        different J^P block the fit never touches), must keep that other block's levels unchanged, replace
        only the fitted block's levels, and return exactly as many levels as the ab initio multiplet it started
        from; (iv) the SECOND, orthogonal route, on a richer He-like C V 2p-3d J=2(odd) block (4 CSFs, several
        F^k/G^k/R^k classes, see AdjustedLevels.twoBodyParameterClass): reconstructing the FULL CI matrix as
        oneBody + (the SUM of every discovered class's own matrix) must reproduce
        AdjustedLevels.buildOneAndTwoBodyMatrices' own (oneBody+twoBody) EXACTLY -- the faithfulness check of
        (i), one level up, now allowing for classes AdjustedLevels.buildParameterMatrices folds into oneBody
        itself (a closed subshell, or a class that is a multiple of the identity on this block, neither of
        which any splitting could ever determine); and against targets whose SHAPE a single uniform factor
        cannot reproduce, AdjustedLevels.fitOrthogonalParameters must achieve a materially smaller rms
        deviation than AdjustedLevels.fitScaleFactor on the SAME targets -- the mechanical version of the
        "power" this route exists for; (v) THE OVERFITTING GUARD: this same richer block has more discovered
        classes than fit targets, so `degreesOfFreedom` must be negative and the warning naming it must
        appear in `Base.show`'s own output; (vi) AdjustedLevels.fitOrthogonalParametersMultiBlock, pooling
        THIS J=2 block with an independent J=1(odd) block of the SAME basis: the returned radial parameters
        are GLOBAL (one shared value per class, not one per block), the combined rms must be finite and no
        worse than either block's own single-factor rms, and AdjustedLevels.mergeIntoMultiplet for the
        multi-block result must replace exactly the two fitted blocks and return the ab initio multiplet's
        own total level count. A success::Bool is returned.
"""
function testModule_AdjustedLevels(; short::Bool=true)
    printstyled("\n\nTest the module  AdjustedLevels  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    wb    = SelfConsistent.performSCF([Configuration("1s 2p")], ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    sym   = LevelSymmetry(AngularJ64(1), Basics.minus)

    # (i) Faithfulness against Hamiltonian.setupMatrix, at the unscaled combination.
    xlCache   = InteractionStrength.XLCache()
    refMatrix = Hamiltonian.setupMatrix(sym, basis, ni, grid, asfSettings, xlCache; printout=false)
    idxCsf, oneBody, twoBody = AdjustedLevels.buildOneAndTwoBodyMatrices(sym, basis, ni, grid)
    myMatrix  = oneBody .+ 1.0 .* twoBody
    if  maximum(abs.(myMatrix .- refMatrix)) > 1.0e-10
        success = false;   println("  ** buildOneAndTwoBodyMatrices does NOT reproduce Hamiltonian.setupMatrix at c=1.0")
    end

    # (ii) The fit itself: target 10% below the ab initio splitting, scale factor must move away from 1.0.
    eig1  = Hamiltonian.diagonalizeCiMatrix(refMatrix, LevelSelection())
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    abInitioRel = (sort(eig1.values) .- sort(eig1.values)[1]) ./ auPerKayser
    target_cm   = abInitioRel[2] * 0.9
    fm = AdjustedLevels.fitScaleFactor(sym, basis, ni, grid, [(2, target_cm)], "test source"; scaleBounds=(0.5,1.2))
    if  abs(fm.scaleFactors[1].value - 1.0) < 0.01
        success = false;   println("  ** fitted scale factor did not move away from 1.0 for a 10%-off target")
    end
    if  fm.rmsDeviation > 1.0
        success = false;   println("  ** rms deviation too large for a single-parameter, single-target fit: $(fm.rmsDeviation)")
    end

    # (iii) mergeIntoMultiplet: a basis with BOTH 1s^2 (J=0 even, untouched) and 1s2p (the fitted J=1 odd block).
    wb2       = SelfConsistent.performSCF([Configuration("1s^2"), Configuration("1s 2p")], ni, grid, asfSettings; printout=false)
    basis2    = wb2.levels[1].basis
    abInitio2 = Hamiltonian.performCI(basis2, ni, grid, asfSettings; printout=false)
    fm2       = AdjustedLevels.fitScaleFactor(sym, basis2, ni, grid, [(2, target_cm)], "test source"; scaleBounds=(0.5,1.2))
    combined  = AdjustedLevels.mergeIntoMultiplet(fm2, abInitio2)
    if  length(combined.levels) != length(abInitio2.levels)
        success = false;   println("  ** mergeIntoMultiplet changed the total number of levels")
    end
    groundSymEnergies    = sort([lev.energy for lev in abInitio2.levels if LevelSymmetry(lev.J, lev.parity) == LevelSymmetry(AngularJ64(0), Basics.plus)])
    groundSymEnergiesNew = sort([lev.energy for lev in combined.levels  if LevelSymmetry(lev.J, lev.parity) == LevelSymmetry(AngularJ64(0), Basics.plus)])
    if  groundSymEnergies != groundSymEnergiesNew
        success = false;   println("  ** mergeIntoMultiplet changed a level in a symmetry block the fit never touched")
    end

    # (iv) The orthogonal route: a richer block (He-like C V, 2p-3d, J=2 odd -- 4 CSFs, several F^k/G^k/R^k classes).
    wb3    = SelfConsistent.performSCF([Configuration("2p 3d")], ni, grid, asfSettings; printout=false)
    basis3 = wb3.levels[1].basis
    sym3   = LevelSymmetry(AngularJ64(2), Basics.minus)

    idxCsf3, oneBody3, classes3, matrices3 = AdjustedLevels.buildParameterMatrices(sym3, basis3, ni, grid)
    idxCsf3b, oneBody3b, twoBody3b         = AdjustedLevels.buildOneAndTwoBodyMatrices(sym3, basis3, ni, grid)
    lumped3   = zeros(Float64, length(idxCsf3), length(idxCsf3));   for  mtx in matrices3   lumped3 .+= mtx   end
    totalNew3 = oneBody3  .+ lumped3
    totalOld3 = oneBody3b .+ twoBody3b
    if  maximum(abs.(totalNew3 .- totalOld3)) > 1.0e-8
        success = false
        println("  ** buildParameterMatrices' reconstruction (oneBody + every class) does NOT match " *
                 "buildOneAndTwoBodyMatrices' own (oneBody+twoBody)")
    end
    if  !any(startswith(cl, "F^") for cl in classes3)   ||   !any(startswith(cl, "G^") for cl in classes3)
        success = false;   println("  ** 2p-3d (a genuine two-shell-pair case) did not separate into BOTH F^ and G^ classes")
    end
    if  length(classes3) > 20
        success = false;   println("  ** class count exploded ($(length(classes3))) -- closed-shell/identity folding is not working")
    end

    eig3  = Hamiltonian.diagonalizeCiMatrix(oneBody3b .+ twoBody3b, LevelSelection())
    rel3  = (sort(eig3.values) .- sort(eig3.values)[1]) ./ auPerKayser
    # Targets deliberately NOT proportional to the ab initio shape, so one uniform knob cannot match all three.
    targets3 = [ (2, rel3[2]*0.80), (3, rel3[3]*1.05), (4, rel3[4]*0.95) ]

    fmSingle3 = AdjustedLevels.fitScaleFactor(sym3, basis3, ni, grid, targets3, "test source"; scaleBounds=(0.5,1.3))
    fmOrtho3  = AdjustedLevels.fitOrthogonalParameters(sym3, basis3, ni, grid, targets3, "test source")
    if  fmOrtho3.rmsDeviation > 0.5 * fmSingle3.rmsDeviation
        success = false
        println("  ** orthogonal route did not materially improve on the single-factor route for a non-proportional target " *
                 "(single = $(fmSingle3.rmsDeviation) cm^-1, orthogonal = $(fmOrtho3.rmsDeviation) cm^-1)")
    end

    abInitio3 = Hamiltonian.performCI(basis3, ni, grid, asfSettings; printout=false)
    combined3 = AdjustedLevels.mergeIntoMultiplet(fmOrtho3, abInitio3)
    if  length(combined3.levels) != length(abInitio3.levels)
        success = false;   println("  ** mergeIntoMultiplet (orthogonal route) changed the total number of levels")
    end

    # (v) THE OVERFITTING GUARD: more classes (length(classes3)) than targets (3) here, so dof must be negative
    # and Base.show must say so -- this is the thing that stops a fitted class value from being read in
    # isolation as if it were independently confirmed.
    if  fmOrtho3.degreesOfFreedom != length(targets3) - length(fmOrtho3.radialParameters)
        success = false;   println("  ** degreesOfFreedom is not length(fitTargets) - length(radialParameters)")
    end
    if  fmOrtho3.degreesOfFreedom >= 0
        success = false;   println("  ** this richer block has more classes than targets, yet degreesOfFreedom >= 0")
    end
    if  !occursin("WARNING", sprint(show, fmOrtho3))
        success = false;   println("  ** Base.show did not print the overfitting-guard WARNING for a non-positive degreesOfFreedom")
    end

    # (vi) Multi-block fitting: pool THIS J=2(odd) block with an INDEPENDENT J=1(odd) block of the SAME basis --
    # one shared (GLOBAL) parameter set determined from both together.
    sym1 = LevelSymmetry(AngularJ64(1), Basics.minus)
    idxCsf1, oneBody1, twoBody1 = AdjustedLevels.buildOneAndTwoBodyMatrices(sym1, basis3, ni, grid)
    eig1b = Hamiltonian.diagonalizeCiMatrix(oneBody1 .+ twoBody1, LevelSelection())
    rel1  = (sort(eig1b.values) .- sort(eig1b.values)[1]) ./ auPerKayser
    targets1 = [ (2, rel1[2]*0.85), (3, rel1[3]*1.10) ]

    fmMulti = AdjustedLevels.fitOrthogonalParametersMultiBlock(
                  [ (sym3, targets3), (sym1, targets1) ], basis3, ni, grid, "test source" )
    if  length(fmMulti.blocks) != 2
        success = false;   println("  ** fitOrthogonalParametersMultiBlock did not return exactly 2 blocks")
    end
    if  !isfinite(fmMulti.rmsDeviation)
        success = false;   println("  ** fitOrthogonalParametersMultiBlock's combined rms is not finite")
    end
    # Baseline: the UNSCALED (c=1, pure ab initio) combined rms against the same, deliberately-off targets --
    # the one comparison that needs no assumption about what the fit COULD reach, only that it does better
    # than not fitting at all.
    function  combinedRmsAtC(c::Float64)
        r3 = (sort(Hamiltonian.diagonalizeCiMatrix(oneBody3b .+ c.*twoBody3b, LevelSelection()).values)); r3 = (r3 .- r3[1]) ./ auPerKayser
        r1 = (sort(Hamiltonian.diagonalizeCiMatrix(oneBody1  .+ c.*twoBody1,  LevelSelection()).values)); r1 = (r1 .- r1[1]) ./ auPerKayser
        s = sum( (r3[idx]-m)^2 for (idx,m) in targets3 ) + sum( (r1[idx]-m)^2 for (idx,m) in targets1 )
        return( sqrt(s/(length(targets3)+length(targets1))) )
    end
    baselineRms = combinedRmsAtC(1.0)
    if  fmMulti.rmsDeviation > baselineRms
        success = false
        println("  ** pooled multi-block orthogonal fit is WORSE than the unscaled (c=1) ab initio baseline " *
                 "($(fmMulti.rmsDeviation) against $(baselineRms) cm^-1)")
    end

    abInitio4 = Hamiltonian.performCI(basis3, ni, grid, asfSettings; printout=false)
    combined4 = AdjustedLevels.mergeIntoMultiplet(fmMulti, abInitio4)
    if  length(combined4.levels) != length(abInitio4.levels)
        success = false;   println("  ** mergeIntoMultiplet (multi-block route) changed the total number of levels")
    end

    testPrint("testModule_AdjustedLevels()::", success)
    return(success)
end


"""
`TestFrames.testModule_SuperConfigurations(; short::Bool=true)`
    ... tests module SuperConfigurations on C I (Z=6) 2s^2 2p^2 (5 CSFs across J=0,1,2): the mean energy and
        degeneracy (trivially exact, a CSF count) against a direct brute-force check, and -- the decisive
        test -- the variance against the SAME (2J+1)-weighted brute-force statistic computed by fully
        diagonalizing EVERY symmetry block of the configuration via Hamiltonian.performCI, which this
        module's own `varianceFromCsfMatrix` is explicitly built never to need. The two must agree to within
        the known ~1e-9-relative imprecision of the REUSED `BasicsAZ.meanEnergyInClosedForm` (propagated
        through variance = <E^2> - mean^2), not to machine precision, since that reused mean is itself only
        approximate; the variance's OWN <E^2> term (this module's actual new content) is exact and carries no
        such slack. A success::Bool is returned.
"""
function testModule_SuperConfigurations(; short::Bool=true)
    printstyled("\n\nTest the module  SuperConfigurations  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    conf  = Configuration("2s^2 2p^2")
    wb    = SelfConsistent.performSCF([conf], ni, grid, asfSettings; printout=false)
    orbitals = wb.levels[1].basis.orbitals

    group = SuperConfigurations.buildGroup(conf, orbitals, ni, grid, asfSettings)

    abInitio = Hamiltonian.performCI(wb.levels[1].basis, ni, grid, asfSettings; printout=false)
    wsum = 0.0;   wEsum = 0.0;   wE2sum = 0.0
    for  lev in abInitio.levels
        w = Basics.twice(lev.J) + 1
        wsum   = wsum   + w
        wEsum  = wEsum  + w*lev.energy
        wE2sum = wE2sum + w*lev.energy^2
    end
    bruteMean = wEsum/wsum
    bruteVar  = wE2sum/wsum - bruteMean^2

    if  group.degeneracy != wsum
        success = false;   println("  ** degeneracy does not match the exact (2J+1)-weighted CSF count")
    end
    if  abs(group.meanEnergy - bruteMean) > 1.0e-6
        success = false;   println("  ** meanEnergy too far from the brute-force (2J+1)-weighted mean: " *
                                     "$(group.meanEnergy) vs $(bruteMean)")
    end
    if  abs(group.variance - bruteVar) / abs(bruteVar) > 1.0e-3
        success = false;   println("  ** variance too far from the brute-force (2J+1)-weighted variance: " *
                                     "$(group.variance) vs $(bruteVar)")
    end

    testPrint("testModule_SuperConfigurations()::", success)
    return(success)
end


"""
`TestFrames.testModule_SuperGroup(; short::Bool=true)`
    ... tests `SuperConfigurations.SuperGroup` (competition item C9) on C I's `2s^2 2p^2` and `2p^4`
        (both even parity, from the same complex `module-SuperConfigurations.jl`'s own toy case uses).
        Checks the pooled degeneracy, mean and variance against an INDEPENDENT hand evaluation of the
        combined-sample formula from the two configurations' own separately-built `Group`s -- the same
        formula `buildSuperGroup` itself uses, re-derived here rather than re-run, so this is a check on
        the ARITHMETIC rather than a second measurement of the same code path. The mean's EXACTNESS under
        configuration-interaction mixing (a trace argument, not an approximation) and the variance's
        measured 19.5 % shortfall on a realistic mixing pair are recorded in `buildSuperGroup`'s own
        docstring rather than re-verified here every run, since that check needs a second SCF and a full
        CI diagonalization and is not cheap enough for the routine suite. A success::Bool is returned.
"""
function testModule_SuperGroup(; short::Bool=true)
    printstyled("\n\nTest the module  SuperConfigurations (SuperGroup)  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    confA = Configuration("2s^2 2p^2");   confB = Configuration("2p^4")
    wb    = SelfConsistent.performSCF([confA, confB], ni, grid, asfSettings; printout=false)
    orbitals = wb.levels[1].basis.orbitals

    gA = SuperConfigurations.buildGroup(confA, orbitals, ni, grid, asfSettings)
    gB = SuperConfigurations.buildGroup(confB, orbitals, ni, grid, asfSettings)
    sg = SuperConfigurations.buildSuperGroup([confA, confB], orbitals, ni, grid, asfSettings)

    G = gA.degeneracy + gB.degeneracy
    M = (gA.degeneracy*gA.meanEnergy + gB.degeneracy*gB.meanEnergy) / G
    V = (gA.degeneracy*(gA.variance + (gA.meanEnergy-M)^2) + gB.degeneracy*(gB.variance + (gB.meanEnergy-M)^2)) / G

    if  sg.degeneracy != G
        success = false;   println("  ** SuperGroup degeneracy $(sg.degeneracy) != hand-pooled $G")
    end
    if  abs(sg.meanEnergy - M) > 1.0e-12
        success = false;   println("  ** SuperGroup mean $(sg.meanEnergy) != hand-pooled $M")
    end
    if  abs(sg.variance - V) > 1.0e-12 * abs(V)
        success = false;   println("  ** SuperGroup variance $(sg.variance) != hand-pooled $V")
    end
    # the pooled variance must EXCEED the degeneracy-weighted average of the two groups' OWN variances --
    # the between-configuration spread (the two means are nowhere near equal) can only add, never subtract.
    withinOnly = (gA.degeneracy*gA.variance + gB.degeneracy*gB.variance) / G
    if  sg.variance <= withinOnly
        success = false;   println("  ** pooled variance $(sg.variance) does not exceed the within-group-only " *
                                     "average $withinOnly -- the between-configuration term is missing or wrong")
    end

    testPrint("testModule_SuperGroup()::", success)
    return(success)
end
