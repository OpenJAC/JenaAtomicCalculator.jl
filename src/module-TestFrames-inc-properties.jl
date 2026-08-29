
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
    outA  = perform(waA; output=true)["RDM outcomes:"]
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
    outB   = perform(waB; output=true)["RDM outcomes:"]
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
    starkOutcomes = function(EField::Float64)
        wb1 = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = "H perturbers", grid = grid, nuclearModel = nm,
                                        configs = perturbers); output = true )
        end
        settings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts = true, gMultiplet = wb1["multiplet:"],
                                       EField = EField)
        wb = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = "H 2p Stark shift", grid = grid, nuclearModel = nm,
                                        configs = [Configuration("2p")], propertySettings = [settings]); output = true )
        end
        return( wb["Stark-shift outcomes:"] )
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

    println(iostream, "StarkShift: for hydrogen 2p the tensor term is traceless in M, the sublevels are degenerate " *
                      "in M -> -M, J = 1/2 leaves alpha_2 exactly zero, and doubling the field multiplies every " *
                      "shift by exactly four. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_StarkShift()::", success)
    return( success )
end
