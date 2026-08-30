
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
    multiplet = redirect_stdout(devnull) do;   perform(wa; output=true)["multiplet:"]   end
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
