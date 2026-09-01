
## Functions in this file cover: atomic processes.
## Alphabetical order within this file.


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
                                joinpath(@__DIR__, "..", "test", "test-AutoIonization-new.sum"), "Auger rates and intrinsic", 10)
    testPrint("testModule_AutoIonization()::", success)
    return(success)
end


"""
`TestFrames.testModule_CoulombExcitation(; short::Bool=true)`  ... tests on module CoulombExcitation.
"""
function testModule_CoulombExcitation(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-CoulombExcitation-new.sum")
    printstyled("\n\nTest the module  CoulombExcitation  ... \n", color=:cyan)
    ### Make the tests
    ## THIS TEST USED TO ASSERT NOTHING (rewritten 08-Aug-2026). Its body between the two markers was EMPTY, its
    ## file comparison was commented out -- and pointed at the AutoIonization approved file, a copy-paste
    ## artifact -- and it returned success = true unconditionally. It therefore reported a pass for a module
    ## that was never exercised, in every run of the suite.
    ##
    ## WHAT IS CHECKED NOW is the SCAFFOLD, not the physics: that the settings construct, copy-construct and
    ## print. That is a deliberately modest bar, and it is not nothing -- three bug classes met this month would
    ## have been caught by exactly these lines: a `Base.show` declared without its `io` argument while using it
    ## (six methods at once in MultiPhotonTransition), a Settings field commented out while still referenced,
    ## and a copy-constructor that silently fails to carry a value through. The last of those was found HERE:
    ## the ionEnergies guard assigned the parameter instead of the local, so every copy-construction that did
    ## not pass ionEnergies explicitly raised an UndefVarError.
    ##
    ## NOT COVERED: any Coulomb-excitation amplitude, cross section or alignment parameter. When an approved
    ## reference exists, it belongs here.
    success = true
    settings = CoulombExcitation.Settings()
    if  length(sprint(show, settings)) == 0     success = false;   println("** empty show for CoulombExcitation.Settings")   end
    ## the copy-constructor must carry a changed value through AND leave the others intact
    newSettings = CoulombExcitation.Settings(settings; printBefore=true, zerosGL=7)
    if  !newSettings.printBefore                success = false;   println("** printBefore not carried by the copy-constructor")   end
    if  newSettings.zerosGL != 7                success = false;   println("** zerosGL not carried by the copy-constructor")       end
    if  newSettings.calcAlignment != settings.calcAlignment
                                                success = false;   println("** calcAlignment altered by an unrelated keyword")     end
    ##
    ## THE PHYSICS CHECK, added 09-Aug-2026: sigma(Mi,Mf) = sigma(-Mi,-Mf), EXACTLY.
    ##
    ## This is the very symmetry whose violation exposed the root-cause bug of the July work on this module --
    ## computeAmplitude's magnetic term was missing the -i prefactor that RATIP's coulex_pure_matrix() applies
    ## to every magnetic contribution but which Eq. (8) of Surzhykov et al., PRA 77, 042722 (2008) does not show.
    ## With the factor restored the equality holds to the last bit; without it, it fails visibly. It is
    ## parameter-free, needs no reference data, and is precisely the guard that would catch a reintroduction.
    ##
    ## COST: about 50 s, dominated by the Z = 92 self-consistent field, not by the Coulomb-excitation quadrature
    ## (verified -- shrinking the grid and halving zerosGL changed nothing). One impact energy is used.
    ceSettings = CoulombExcitation.Settings(CoulombExcitation.Settings(); ionEnergies=[100.], calcAlignment=true,
                                            printBefore=false, zerosGL=4)
    wa = Atomic.Computation(Atomic.Computation(), name="U90+ K-L Coulomb excitation",
                            grid = Radial.Grid(Radial.Grid(false), rnt=2.0e-5, h=5.0e-2, hp=1.5e-2, rbox=6.0),
                            nuclearModel = Nuclear.Model(92.),
                            initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s 2p")],
                            processSettings = ceSettings )
    ceLines = perform(wa; output=true)["Coulomb excitation lines:"]
    if  length(ceLines) == 0    success = false;   println("** no Coulomb-excitation lines were computed")   end
    for  line in ceLines
        partial = Dict( (mL.Mi, mL.Mf) => mL.partialCs   for mL in line.mLines )
        for  ((Mi, Mf), cs) in partial
            mMi = AngularM64(-Basics.twice(Mi)//2);    mMf = AngularM64(-Basics.twice(Mf)//2)
            if  haskey(partial, (mMi, mMf))  &&  cs != 0.
                if  abs(partial[(mMi, mMf)] - cs) / abs(cs) > 1.0e-10
                    success = false
                    println("** sigma(Mi,Mf) != sigma(-Mi,-Mf) for ($Mi,$Mf): " *
                            "$(cs) vs $(partial[(mMi, mMf)]) -- the magnetic-term phase is the first suspect")
                end
            end
        end
    end
    ###
    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "CoulombExcitation: scaffold checks plus the exact sigma(Mi,Mf) = sigma(-Mi,-Mf) symmetry; " *
                      "no cross section is compared against approved data.")
    testPrint("testModule_CoulombExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_DielectronicRecombination(; short::Bool=true)`  ... tests on module DielectronicRecombination.
"""
function testModule_DielectronicRecombination(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-DielectronicRecombination-new.sum")
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
                            "Total Auger rates", 7)
    testPrint("testModule_DielectronicRecombination()::", success)
    return(success)
end


"""
`TestFrames.testModule_MultiPhotonTransition(; short::Bool=true)`  ... tests on module MultiPhotonTransition.
"""
function testModule_MultiPhotonTransition(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-MultiPhotonTransition-new.sum")
    printstyled("\n\nTest the module  MultiPhotonTransition  ... \n", color=:cyan)
    ### Make the tests
    ## THIS TEST USED TO COMPUTE NOTHING. It set `success = true`, had its file comparison commented out, and
    ## contained no computation at all between the two setDefaults calls -- so the suite counted a pass for a
    ## module that had never been exercised. What follows are three real checks of the scaffold; the PHYSICS
    ## check (H 2s -> 1s = 8.2206 /s) belongs here once the amplitudes are verified.
    success = true
    ## (1) every scheme must construct, copy-construct and print.
    for  scheme in [MultiPhotonTransition.TwoPhotonEmissionScheme(), MultiPhotonTransition.TwoPhotonAbsorptionScheme(),
                    MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme(),
                    MultiPhotonTransition.ThreePhotonEmissionScheme(), MultiPhotonTransition.ThreePhotonAbsorptionScheme()]
        set = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings(); scheme=scheme)
        if  typeof(set.scheme) != typeof(scheme)    success = false;   println("** scheme not carried: $scheme")   end
        sprint(show, scheme)      ## must not raise
    end
    ## (2) every property must be printable. All six of these once declared Base.show WITHOUT an io argument
    ##     while using `io` in the body, so each raised an UndefVarError; that is what this guards against.
    for  property in [MultiPhotonTransition.EnergyDiffCs(), MultiPhotonTransition.TotalAlpha0(),
                      MultiPhotonTransition.TotalCsLinear(), MultiPhotonTransition.TotalCsRightCircular(),
                      MultiPhotonTransition.TotalCsUnpolarized(), MultiPhotonTransition.TotalCsDensityMatrix()]
        if  length(sprint(show, property)) == 0    success = false;   println("** empty show: $property")    end
    end
    ## (3) three-photon EMISSION must still FAIL, and fail informatively rather than in a bare MethodError.
    ##     ABSORPTION was in this list until 08-Aug-2026 and has been REMOVED because it now runs: the two halves
    ##     are deliberately in different states, since absorption fixes the three photon energies and needs no
    ##     sharings, while emission fixes only their sum and needs the simplex. A test asserting that a working
    ##     scheme raises would be a test of the wrong thing.
    for  scheme in [MultiPhotonTransition.ThreePhotonEmissionScheme()]
        set  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings(); scheme=scheme)
        threw = false
        try   MultiPhotonTransition.computeLines(scheme, Multiplet(), Multiplet(), Radial.Grid(true), set)
        catch e
            threw = true
            if  !occursin("not yet implemented", sprint(showerror, e))
                success = false;   println("** three-photon error is not informative: $e")
            end
        end
        if  !threw    success = false;   println("** three-photon scheme did not raise: $scheme")    end
    end
    ## (4) the three-photon SHARINGS must reproduce the exact simplex moments; this is the sum rule that makes
    ##     them a finished piece rather than plausible code. Four points per direction suffice for all four.
    for  (name, exact, quad, err) in MultiPhotonTransition.checkSharings_3p(0.5, 4)
        if  err > 1.0e-12    success = false;   println("** three-photon simplex moment '$name' off by $err")    end
    end
    ## (5) the three-photon ABSORPTION scheme must carry its two photon energies through Settings, and its
    ##     monochromatic sentinel must split the transition energy into three equal parts.
    scheme3 = MultiPhotonTransition.ThreePhotonAbsorptionScheme(
                  MultiPhotonTransition.AbstractMultiPhotonProperty[MultiPhotonTransition.TotalCsUnpolarized()], 0., 0.)
    omegas  = MultiPhotonTransition.determineOmegas_3pAbsorption(0.6, scheme3)
    if  !(omegas[1] ≈ 0.2  &&  omegas[2] ≈ 0.2  &&  omegas[3] ≈ 0.2)
        success = false;   println("** three-photon monochromatic sentinel gave $omegas, not three times 0.2")
    end
    ###
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_MultiPhotonTransition()::", success)
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
                                LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)]), 0., 0., 10000., false ) )
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
                            processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, false, true,
                                                                        LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes(0., 0., 0.) ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitation-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitation-new.sum"),
                            "Photoexcitation integrated cross sections", 200)
    testPrint("testModule_PhotoExcitation()::", success)
    return(success)
end


"""
`TestFrames.testModule_MultiPhotonIonization(; short::Bool=true)`  ... tests on module MultiPhotonIonization;
    a success::Bool is returned.
"""
function testModule_MultiPhotonIonization(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-MultiPhotonIonization-new.sum")
    printstyled("\n\nTest the module  MultiPhotonIonization  ... \n", color=:cyan)
    # WHAT IS **NOT** ASSERTED HERE, and the reason belongs in the file rather than in a list: GAUGE EQUALITY.
    # For a one-photon amplitude the length and velocity forms agree at the physical photon energy and that makes
    # a powerful test (see testModule_PhotoIonization). A SECOND-ORDER amplitude is different: every intermediate
    # step is evaluated OFF SHELL -- the operator carries the driving frequency omega, not the level difference --
    # and the two forms are then equivalent only through CLOSURE over a COMPLETE intermediate set. This module
    # sums bound states, so the set is incomplete and the gauges are EXPECTED to differ, by a lot. Asserting their
    # equality here would produce a test that fails for a correct code.
    #
    # The three statements below hold regardless of that, and none needs a stored reference.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")

    function setup(Z::Float64)
        # everything scales as 1/Z, so the grids are scaled copies and the Z-law below is exact, not approximate
        grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6/Z, h = 5.0e-2, hp = 1.0e-2/Z, rbox = 90.0/Z)
        nm   = Nuclear.Model(Z, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
        prim = Bsplines.generatePrimitives(grid)
        shs  = Subshell[]
        for  n = 1:8
            push!(shs, Subshell(n,-1))
            if n >= 2   push!(shs, Subshell(n,1));   push!(shs, Subshell(n,-2))   end
            if n >= 3   push!(shs, Subshell(n,2));   push!(shs, Subshell(n,-3))   end
        end
        orbs = redirect_stdout(devnull) do;  Bsplines.generateOrbitalsHydrogenic(shs, nm, prim; printout=false)  end
        return( (orbs, Nuclear.nuclearPotential(nm, grid)) )
    end

    # (1) THE HYDROGENIC Z^-6 LAW, which is an ABSOLUTE statement and the strongest thing available here.
    #     sigma^(2) is not an area: it has units of cm^4 s, because a two-photon rate goes as the flux SQUARED.
    #     In atomic units that is length^4 x time, and under Z-scaling length -> 1/Z and time -> 1/Z^2, so at a
    #     FIXED omega/I every hydrogen-like ion must give exactly Z^-6 times the hydrogen value. No wrong power of
    #     omega, of Z, or of the flux normalisation can satisfy that, and it needs no reference data at all.
    ref = 0.
    for  Z  in [1.0, 2.0]
        orbs, pot = setup(Z)
        iOrb  = orbs[Subshell("1s_1/2")]
        omega = 0.30 * Z^2                      # the same fraction of the threshold for every Z
        cs, _ = redirect_stdout(devnull) do
            MultiPhotonIonization.crossSection2pOneElectron(:linear, iOrb, omega, omega, [E1],
                                                            Basics.Babushkin, orbs, pot)   end
        if  Z == 1.0    ref = cs
        elseif  abs( cs * Z^6 / ref - 1.0 ) > 2.0e-2
            success = false
            if printTest   info(iostream, "MultiPhotonIonization, Z = $Z: Z^6 sigma2 = $(cs*Z^6) against " *
                                          "$(ref) at Z = 1, a departure of $(abs(cs*Z^6/ref - 1.0)) from Z^-6")  end
        end
    end

    # (2) THE MONOCHROMATIC AND BICHROMATIC CONVENTIONS MUST MEET as omega_1 -> omega_2. The single-beam factor is
    #     a CONVENTION, and the module pins it by requiring exactly this; a mismatched factor of 2 between the two
    #     is the classic error in two-photon work and would show up here as a factor of 2.
    orbs, pot = setup(1.0)
    iOrb = orbs[Subshell("1s_1/2")]
    csMono, _ = redirect_stdout(devnull) do
        MultiPhotonIonization.crossSection2pOneElectron(:linear, iOrb, 0.30, 0.30, [E1], Basics.Babushkin, orbs, pot)  end
    csBi, _   = redirect_stdout(devnull) do
        MultiPhotonIonization.crossSection2pOneElectron(:linear, iOrb, 0.30 - 1.0e-6, 0.30 + 1.0e-6, [E1],
                                                        Basics.Babushkin, orbs, pot)  end
    if  abs(csMono/csBi - 1.0) > 1.0e-3
        success = false
        if printTest   info(iostream, "MultiPhotonIonization: monochromatic $csMono and bichromatic $csBi do not " *
                                      "meet as omega_1 -> omega_2; ratio $(csMono/csBi)")   end
    end

    # (3) THE RESONANCE GUARD REFUSES RATHER THAN RETURNING A NUMBER, and it is RELATIVE for a reason worth
    #     keeping: an ABSOLUTE tolerance of 1e-6 a.u. did NOT trip on hydrogen's own 1s -> 2p resonance, whose
    #     computed position is 0.3750046, leaving a denominator of 4.6e-6 that looked "safe" while inflating that
    #     one term by 2e5. Dropping the term instead of refusing would remove most of the amplitude and hand back
    #     a small number in place of no number.
    try
        redirect_stdout(devnull) do
            MultiPhotonIonization.crossSection2pOneElectron(:linear, iOrb, 0.3750046, 0.3750046, [E1],
                                                            Basics.Babushkin, orbs, pot)   end
        success = false
        if printTest   info(iostream, "MultiPhotonIonization: the 1s -> 2p resonance at omega = 0.3750046 was " *
                                      "not refused")   end
    catch                                                                                              end

    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_MultiPhotonIonization()::", success)
    return( success )
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
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"), "Total photoionization c", 10)
    ## Check the summed (grand-total) cross sections separately: they must stay consistent with the
    ## line-resolved table above, of which they are the sum over all final levels.
    success = success  &&
              testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"),
                                "Total photoionization cross sections, summed", 15)
    # ----------------------------------------------------------------------------------------------------------
    # TWO EXACT INVARIANTS, added 01-Sep-2026, and they test something the comparison above CANNOT.
    #
    # `testCompareFiles` is a REGRESSION check: it asks whether today's numbers match a stored .sum file. That
    # file was produced by this same code, so a systematically wrong NORMALISATION agrees with it forever. This
    # module had exactly such a fault -- until 28-Aug-2026 one route used 4 pi^2 alpha omega |d|^2, the textbook
    # expression for a DIPOLE matrix element, while JAC's amplitude is built on j_L(alpha omega r) and already
    # carries a power of alpha omega, so the prefactor must go as 1/(alpha omega). A correct formula applied to
    # the wrong object, invisible to any stored reference.
    #
    # Both checks below are statements that must hold exactly, so neither needs a reference and neither can drift.
    # They test the TOTAL cross section and deliberately not the M_f-resolved partial one: the partials do not
    # sum to the total, they warn at the point of use, and nothing in the package consumes them (challenge 63).
    function hydrogenicCs(Z::Float64, x::Float64)
        # x = omega/omega_threshold, so every Z is the SAME physical point on the cross-section curve.
        # NuclearField is asked for explicitly: a one-electron system in a field built from its own density is
        # repelled by itself, which is what put these two gauges 24-30 % apart for weeks.
        omega = x * Z^2/2                                                # threshold = Z^2/2 Hartree for 1s
        # THE GRID IS A SCALED COPY, and it is NOT Basics.recommendedGrid. That function sizes the box from the
        # BOUND orbitals alone and returns rbox = 18.09 a.u. for H 1s -- too short to hold the outgoing
        # photoelectron, on which these same two gauges come out 1.090, 1.073 and 0.840 at x = 1.5, 2, 3 instead
        # of 1, and the cross section is 1.4 % low. Nothing warns: Continuum.gridConsistency checks the STEP
        # against the electron's wavelength and Bsplines.checkGridRepresentation checks BOUND orbitals, so a box
        # too short for the continuum falls between them. Measured 01-Sep-2026; carried as its own item.
        # Everything scales as 1/Z, so the three charges below are the same calculation in scaled units and the
        # 1/Z^2 law is exact rather than approximate.
        grid  = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6/Z, h = 5.0e-2, hp = 1.0e-2/Z, rbox = 30.0/Z)
        asf   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
        piS   = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles = [E1],
                        gauges = [UseCoulomb, UseBabushkin], printBefore = false,
                        photonEnergies = [Defaults.convertUnits("energy: from atomic to eV", omega)])
        wc = Atomic.Computation(Atomic.Computation(), name = "H-like", grid = grid,
                 nuclearModel   = Nuclear.Model(Z, PointNucleus()),
                 initialConfigs = [Configuration("1s")],   initialAsfSettings = asf,
                 finalConfigs   = [Configuration("1s^0")], finalAsfSettings   = asf, processSettings = piS)
        wd = redirect_stdout(devnull) do;  perform(wc; output = true)["photoionization lines:"]  end
        return( wd[1].crossSection )
    end

    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # (1) THE TWO GAUGES AGREE ON A ONE-ELECTRON SYSTEM. For an exact hydrogenic wavefunction the length and
    #     velocity forms of the dipole are equivalent, so Coulomb/Babushkin is 1 exactly. This ratio sat at
    #     1.24-1.30 for weeks and was blamed on the continuum machinery; the cause was the self-interaction above.
    for  Z  in [1.0, 2.0]
        cs = hydrogenicCs(Z, 2.0)
        if  abs(cs.Coulomb / cs.Babushkin - 1.0) > 1.0e-3
            success = false
            if printTest   info(iostream, "PhotoIonization, Z = $Z: gauges give Cou/Bab = " *
                                          "$(cs.Coulomb/cs.Babushkin), not 1")   end
        end
    end

    # (2) THE HYDROGENIC Z-SCALING, an absolute statement about the whole curve rather than one point. At fixed
    #     x = omega/omega_threshold every H-like ion is a scaled copy of hydrogen and the 1s cross section goes
    #     exactly as 1/Z^2; departures are the O((alpha Z)^2) relativistic terms. A wrong power of omega or of Z
    #     anywhere in the prefactor cannot satisfy this.
    csRef = hydrogenicCs(1.0, 2.0).Babushkin
    for  Z  in [2.0, 3.0]
        cs = hydrogenicCs(Z, 2.0).Babushkin
        if  abs( cs * Z^2 / csRef - 1.0 ) > 2.0e-2
            success = false
            if printTest   info(iostream, "PhotoIonization, Z = $Z: Z^2 sigma = $(cs*Z^2) against sigma(Z=1) = " *
                                          "$csRef, a departure of $(abs(cs*Z^2/csRef - 1.0)) from the 1/Z^2 law")  end
        end
    end

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
`TestFrames.testModule_RayleighCompton(; short::Bool=true)`  ... tests on module RayleighCompton.
"""
function testModule_RayleighCompton(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-RayleighCompton-new.sum")
    printstyled("\n\nTest the module  RayleighCompton  ... \n", color=:cyan)
    ### Make the tests
    ## REWRITTEN 08-Aug-2026; it previously asserted nothing at all. See the note in testModule_CoulombExcitation
    ## for what these scaffold checks do and do not cover.
    ##
    ## NOT COVERED: any Rayleigh or Compton scattering amplitude, cross section or Stokes parameter.
    success = true
    settings = RayleighCompton.Settings()
    if  length(sprint(show, settings)) == 0     success = false;   println("** empty show for RayleighCompton.Settings")   end
    newSettings = RayleighCompton.Settings(settings; multipoles=[E1, M1], calcStokes=true)
    if  newSettings.multipoles != [E1, M1]      success = false;   println("** multipoles not carried by the copy-constructor")   end
    if  !newSettings.calcStokes                 success = false;   println("** calcStokes not carried by the copy-constructor")   end
    if  newSettings.calcAngular != settings.calcAngular
                                                success = false;   println("** calcAngular altered by an unrelated keyword")      end
    ###
    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Scaffold checks only for RayleighCompton; NO physics is compared against approved data.")
    testPrint("testModule_RayleighCompton()::", success)
    return(success)
end


"""
`TestFrames.testModule_HyperfineInduced(; short::Bool=true)`  ... tests on module HyperfineInduced.
"""
function testModule_HyperfineInduced(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-HyperfineInduced-new.sum")
    printstyled("\n\nTest the module  HyperfineInduced  ... \n", color=:cyan)
    ### Make the tests
    ## NUCLEAR HYPERFINE MIXING in H-like 229Th89+ -- deliberately the smallest system that exercises the whole
    ## chain: two nuclear states in one hyperfine basis, the mixing that follows, both terms of the amplitude
    ## (nuclear radiation and the borrowed electronic one), and the level lifetimes. One electron, one electronic
    ## level, four hyperfine levels, five lines -- it runs in seconds.
    ##
    ## Basics.NuclearField() rather than the default DFS: a DFS potential self-interacts badly on a one-electron
    ## system (H 1s comes out at -0.194 instead of -0.5 a.u.) and would corrupt exactly the hyperfine matrix
    ## elements under test.
    elemM = Nuclear.reducedTransitionAmplitude(M1, 0.008, 229, AngularJ64(3//2))
    gsTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.plus, energy = 0.0,
                           mu =  0.360, multipoleM = [M1], elementM = [elemM])
    isTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus, energy = 8.356,
                           mu = -0.378, multipoleM = [M1], elementM = [elemM])
    asfTh = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid  = Radial.Grid(Radial.Grid(false), rnt = 1.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid,
                            nuclearModel   = Nuclear.Model(90.0, FermiNucleus()),
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("1s")],
                            initialAsfSettings = asfTh, finalAsfSettings = asfTh,
                            processSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                                multipoles = [M1], hfMultipoles = [M1], gauges = [UseCoulomb],
                                isomers = Nuclear.Isomer[gsTh, isTh], calcOverview = false,
                                lineSelection = LineSelection(), printBefore = false, calcLifetimes = true ) )
    wb = perform(wa)
    ###
    Defaults.setDefaults("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-HyperfineInduced-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-HyperfineInduced-new.sum"),
                                "Hyperfine-induced transition rates", 11)
    testPrint("testModule_HyperfineInduced()::", success)
    return(success)
end


"""
`TestFrames.testModule_ParticleScattering(; short::Bool=true)`  ... tests on module ParticleScattering.
"""
function testModule_ParticleScattering(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ParticleScattering-new.sum")
    printstyled("\n\nTest the module  ParticleScattering  ... \n", color=:cyan)
    ## The module was rebuilt onto Dirac partial waves on 17-Aug-2026; before that its only coverage was a
    ## ParticleScattering.Settings() entry in testStructConstructors, i.e. a test that could not fail. The checks below
    ## are cheap and each of them was verified to FAIL under a matching perturbation of the source.
    success = true
    oldEnergyUnit = Defaults.getDefaults("unit: energy")
    Defaults.setDefaults("unit: energy", "eV")
    #
    ## Test 1: kappa is the only quantum number stored in a PartialWave, and l and j must follow from it.
    for  kappa in [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
        sh = Subshell(101, kappa)
        l  = (kappa < 0) ? -kappa - 1 : kappa
        success = success && Basics.subshell_l(sh) == l
        success = success && Basics.subshell_j(sh) == AngularJ64( (2*abs(kappa) - 1)//2 )
    end
    #
    ## Test 2: e-He elastic scattering at 300 eV, plane wave, static field with local exchange.
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.5e-2, rbox = 10.0)
    thetas = [0.3, 1.0, 2.0, 3.0]
    ## 100 eV rather than 300: the assertions below are identities and bounds, none of which needs a high energy, and
    ## the partial-wave series is much shorter there. Keeping this test cheap matters -- it is one of ~50 in the suite.
    psSet  = ParticleScattering.Settings(ParticleScattering.Settings(), impactEnergies=[100.0], polarThetas=thetas,
                                         polarPhis=[0.0], printBefore=false, epsPartialWave=1.0e-6, maxL=60)
    wc     = Atomic.Computation(Atomic.Computation(), name="e-He elastic", grid=grid, nuclearModel=Nuclear.Model(2.0),
                                initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s^2")],
                                processSettings = psSet )
    wd     = redirect_stdout(devnull) do
        perform(wc, output=true)
    end
    event  = wd["particle-scattering events:"][1]
    pws    = event.partialWaves
    #
    ## Test 3: the series must have converged well inside maxL, i.e. the epsPartialWave criterion ended it, not the
    ##   backstop. This is the check the earlier hard-coded l = 0:30 loop made impossible.
    lMax    = maximum( Basics.subshell_l(Subshell(101, pw.kappa)) for pw in pws )
    success = success && 0 < lMax < 60
    ## ... and BOTH spin-orbit partners must be present for every l >= 1; a single-kappa treatment carries no spin.
    success = success && length(pws) == 2*lMax + 1
    #
    ## Test 4: the total elastic cross section computed two independent ways -- by Gauss-Legendre quadrature of
    ##   |f|^2 + |g|^2 over the scattering angle, and from the partial-wave sum of sin^2(delta) -- must agree.
    sigmaQuad  = event.integrated.sigmaElastic
    sigmaPhase = ParticleScattering.elasticCrossSectionFromPhases(pws, event.impactEnergy)
    success = success && sigmaQuad > 0.  &&  abs(sigmaQuad/sigmaPhase - 1.0) < 1.0e-8
    #
    ## Test 5: the transport cross sections are positive and ordered sigma_1 < sigma_el for forward-peaked scattering.
    success = success && event.integrated.sigmaMomentumTransfer > 0.
    success = success && event.integrated.sigmaViscosity        > 0.
    success = success && event.integrated.sigmaMomentumTransfer < event.integrated.sigmaElastic
    #
    ## Test 6: the Sherman function is bounded by unity at every angle, and it is NOT identically zero -- a non-zero S
    ##   is only possible because the two kappa branches carry different phase shifts.
    success = success && all( abs(obs.sherman) <= 1.0 + 1.0e-10  for obs in event.angular )
    success = success && any( abs(obs.sherman) >  1.0e-12        for obs in event.angular )
    ## ... and it vanishes in the forward direction, where the spin-flip amplitude has no P_l^1 to build on.
    success = success && abs( ParticleScattering.angularObservables(pws, event.impactEnergy, 0., 0.).sherman ) < 1.0e-10
    #
    ## Test 7: the f/g projection must REFUSE a target that it does not describe. It rests on a spinless target, so an
    ##   initial level with J /= 0 has to raise rather than return a plausible-looking number.
    psErr = 0
    badEvent = ParticleScattering.Event(ParticleScattering.Electron(), ParticleScattering.ElasticScattering(),
                                        ParticleScattering.StaticFieldFurnessMcCarthy(), Beam.PlaneWave(),
                                        Level(AngularJ64(1), AngularM64(0), Basics.plus, 0, 0., 0., false, Basis(), Float64[]),
                                        Level(), 1.0, pws, ParticleScattering.ScatteringChannel[],
                                        ParticleScattering.AngularObservables[], ParticleScattering.IntegratedObservables())
    try     ParticleScattering.assertSpinlessTarget(badEvent)
    catch;  psErr = psErr + 1
    end
    ## ... and a twisted beam must raise too, rather than fall back on the plane-wave superposition.
    try     ParticleScattering.beamObservables(Beam.BesselBeam(), pws, event.impactEnergy, 1.0, 0.)
    catch;  psErr = psErr + 1
    end
    success = success && psErr == 2
    #
    ## Test 8: an ABSOLUTE check against PUBLISHED values, and the only assertion here that is not internal.
    ##   Jablonski, Salvat & Powell, J. Phys. Chem. Ref. Data 33, 409 (2004), Table 3 give the transport cross section
    ##   of helium at 100 eV as 0.928 a0^2 from their Dirac-Hartree-Fock potential and 0.834 from Thomas-Fermi-Dirac,
    ##   against 0.754 derived from measurement. Their sigma_tr is our sigma_1 and is the same integral (their Eq. (15)
    ##   of NSRDS-64), so this is like-for-like. With the Furness-McCarthy exchange that ELSEPA also uses -- the default
    ##   model, hence the main run above -- we get 0.9432, i.e. 1.6 % from DHF. The window below is about 4 % wide and
    ##   pins the absolute normalisation of the amplitudes, the phase-shift convention, the transport integral and the
    ##   exchange term together; no internal identity can do that. It costs NOTHING, reusing the run already made.
    ##   (The high-energy Rutherford limit is the other absolute check and lives in examples/example-Ob.jl, where it
    ##   can afford the energies it needs.)
    success = success && 0.90 < event.integrated.sigmaMomentumTransfer < 1.00
    ## Test 9: the POSITRON path, which shares every piece of machinery with the electron one and differs only in the
    ##   potential -- the electrostatic interaction changes sign and no exchange term applies. Two consequences are
    ##   checked, and both are textbook rather than conventional: the s-wave phase shift must change SIGN, an attractive
    ##   potential advancing the phase and a repulsive core retarding it; and the positron, kept away from the nucleus,
    ##   must scatter LESS than an electron in the same static field, most of all backwards. Cheap: at 100 eV the
    ##   positron series is as short as the electron one.
    posSet = ParticleScattering.Settings(ParticleScattering.Settings(), projectile=ParticleScattering.Positron(),
                                         interaction=ParticleScattering.StaticField(), impactEnergies=[100.0],
                                         polarThetas=[Float64(pi)], polarPhis=[0.0], printBefore=false,
                                         epsPartialWave=1.0e-4, maxL=40)   ## signs and orderings only: precision not needed
    eSet   = ParticleScattering.Settings(posSet, projectile=ParticleScattering.Electron())
    posEv, eEv = map( (posSet, eSet) ) do  st
        wcx = Atomic.Computation(Atomic.Computation(), name="positron/electron", grid=grid,
                                 nuclearModel=Nuclear.Model(2.0), initialConfigs=[Configuration("1s^2")],
                                 finalConfigs=[Configuration("1s^2")], processSettings=st )
        (redirect_stdout(devnull) do;  perform(wcx, output=true);  end)["particle-scattering events:"][1]
    end
    dPos = ParticleScattering.phaseShift(posEv.partialWaves, -1)
    dEle = ParticleScattering.phaseShift(eEv.partialWaves,   -1)
    success = success && dPos < 0.  &&  dEle > 0.
    success = success && posEv.integrated.sigmaElastic < eEv.integrated.sigmaElastic
    success = success && posEv.angular[1].dcs          < eEv.angular[1].dcs
    #
    ## Test 10: the exchange model must MATTER and must be energy-dependent. The Furness-McCarthy term used by ELSEPA
    ##   deepens the well relative to the pure static field, so it must raise the elastic cross section; and it is not
    ##   the Slater term of JAC's DFS field, which is built for a bound electron and (measured against Jablonski,
    ##   Salvat & Powell 2004, Table 3) overshoots -- see examples/example-Od.jl.
    ##   No extra computation is needed: the main run above already uses the default Furness-McCarthy model, at the same
    ##   energy and on the same grid as eSet, so `event` IS the FM case and eEv the static one.
    success = success && event.integrated.sigmaElastic > eEv.integrated.sigmaElastic
    success = success && ParticleScattering.phaseShift(event.partialWaves, -1) > dEle

    Defaults.setDefaults("unit: energy", oldEnergyUnit)
    Defaults.setDefaults("print summary: close", "")
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "ParticleScattering: kappa bookkeeping, partial-wave convergence, the two routes to sigma_el, " *
                      "the transport cross sections, the Sherman function, and the two guards.")
    testPrint("testModule_ParticleScattering()::", success)
    return(success)
end


"""
`TestFrames.testModule_GeneralizedOscillatorStrength(; short::Bool=true)`  ... tests on module
    GeneralizedOscillatorStrength; a success::Bool is returned.
"""
function testModule_GeneralizedOscillatorStrength(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-GeneralizedOscillatorStrength-new.sum")
    printstyled("\n\nTest the module  GeneralizedOscillatorStrength  ... \n", color=:cyan)
    # A SANITY test, deliberately not a comparison against a stored .sum.  Both checks below are absolute --
    # one against a closed-form result and one against a statistical weight -- so neither can pass on a stale
    # reference, which is the failure mode a file comparison cannot rule out.  Atomic hydrogen is used because
    # Inokuti (1971) gives its GOS in closed form, Eq. (4.2):  f_2p(q) = 2^13 3^-9 [1 + (4/9) q]^-6.
    # Basics.NuclearField() is essential rather than a detail: with the default DFS field a ONE-electron system
    # acquires a spurious self-interaction, and the initial and final orbitals then solve different one-body
    # operators.  See examples/example-Dpnew.jl, branch a, for the full acceptance test this is cut down from.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    f2pExact(q) = 2.0^13 * 3.0^-9 * (1 + (4/9)*q)^-6
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 1.2e-1, hp = 3.0e-2, rbox = 30.0)
    Defaults.setDefaults("standard grid", grid)
    nModel = Nuclear.Model(1., UniformNucleus(), 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0)
    asfH   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    Ks     = [1.0]
    gosSettings = GeneralizedOscillatorStrength.Settings(GeneralizedOscillatorStrength.Settings();
                      qValues = Ks, calcOpticalLimit = false, printBefore = false)
    wa = Atomic.Computation(Atomic.Computation(), name="GOS sanity: H 1s -> 2p", grid=grid, nuclearModel=nModel,
             initialConfigs = [Configuration("1s")],   initialAsfSettings = asfH,
             finalConfigs   = [Configuration("2p")],   finalAsfSettings   = asfH,
             processSettings = gosSettings)
    lines = perform(wa; output=true)["generalized oscillator strengths:"]

    f2p = zeros(length(Ks));    fsSplit = Dict{String,Vector{Float64}}()
    for  line in lines
        f2p .+= line.gosValues
        fsSplit[string(LevelSymmetry(line.finalLevel.J, line.finalLevel.parity))] = line.gosValues
    end
    # (1) the summed GOS against the closed form; 2e-5 is reached, and 1e-3 leaves room for the grid
    ratio = f2p[1] / f2pExact(Ks[1]^2)
    if  abs(ratio - 1.0) > 1.0e-3
        success = false
        if printTest   info(iostream, "f_2p(q=1) = $(f2p[1]) against the exact $(f2pExact(Ks[1]^2)), ratio $ratio")   end
    end
    # (2) the two spin-orbit partners carry the GOS in the ratio of their statistical weights.  That ratio is
    #     2 only in the NON-RELATIVISTIC limit -- 2p_1/2 and 2p_3/2 are different radial functions -- so the
    #     measured 1.999967 is not an error but the O((alpha Z)^2) = 5e-5 correction at Z = 1.  The tolerance is
    #     set at 1e-3, which still catches what this check is for: a module that mishandled the fine-structure
    #     sum would be wrong by a FACTOR, not by two parts in 1e5, and would look perfectly healthy otherwise.
    if  haskey(fsSplit, "3/2 -")  &&  haskey(fsSplit, "1/2 -")
        wRatio = fsSplit["3/2 -"][1] / fsSplit["1/2 -"][1]
        if  abs(wRatio - 2.0) > 1.0e-3
            success = false
            if printTest   info(iostream, "f(2p_3/2)/f(2p_1/2) = $wRatio, must be 2 by statistical weight")   end
        end
    else
        success = false
        if printTest   info(iostream, "the two 2p levels were not both returned: $(keys(fsSplit))")   end
    end

    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_GeneralizedOscillatorStrength()::", success)
    return(success)
end


"""
`TestFrames.testModule_PhotoRecombinationInterference(; short::Bool=true)`  ... tests on module
    PhotoRecombinationInterference; a success::Bool is returned.
"""
function testModule_PhotoRecombinationInterference(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PhotoRecombinationInterference-new.sum")
    printstyled("\n\nTest the module  PhotoRecombinationInterference  ... \n", color=:cyan)
    # A SANITY test on the module's own internal consistency, needing no stored reference.  With includeDR =
    # false the coherent sum must collapse to plain radiative recombination, so the module has to reproduce
    # PhotoRecombination exactly -- not approximately, since it is then evaluating the same amplitude -- and
    # both the resonant and the interference term must be identically zero.  Li-like Fe23+ is used because it
    # is the case of examples/example-Dvnew.jl, branch a, and is cheap at maxKappa = 2 and one energy.
    success  = true;   printTest, iostream = Defaults.getDefaults("test flag/stream")
    grid     = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 4.0)
    Defaults.setDefaults("standard grid", grid)
    nModel   = Nuclear.Model(26.)
    energies = [2000.0]
    iConfigs = [Configuration("1s^2 2s")];   mConfigs = [Configuration("1s 2s^2 2p")]
    fConfigs = [Configuration("1s^2 2s^2")]

    prSet   = PhotoRecombination.Settings(PhotoRecombination.Settings(); multipoles=[E1], gauges=[UseCoulomb],
                  electronEnergies=energies, maxKappa=2, calcAnisotropy=false, printBefore=false)
    rrLines = perform( Atomic.Computation(Atomic.Computation(), name="RR reference", grid=grid, nuclearModel=nModel,
                  initialConfigs=iConfigs, finalConfigs=fConfigs, processSettings=prSet);
                  output=true )["photo recombination lines:"]
    priSet   = PhotoRecombinationInterference.Settings(PhotoRecombinationInterference.Settings();
                  multipoles=[E1], gauges=[UseCoulomb], electronEnergies=energies, maxKappa=2,
                  includeRR=true, includeDR=false, calcAnisotropy=false, calcPolarization=false, printBefore=false)
    priPaths = perform( Atomic.Computation(Atomic.Computation(), name="RR limit of the interference module",
                  grid=grid, nuclearModel=nModel, initialConfigs=iConfigs, intermediateConfigs=mConfigs,
                  finalConfigs=fConfigs, processSettings=priSet);
                  output=true )["photorecombination-interference pathways:"]

    if  length(priPaths) == 0
        success = false
        if printTest   info(iostream, "no interference pathways were returned")   end
    end
    for  p  in  priPaths
        # (1) the resonant and interference terms must vanish IDENTICALLY when DR is switched off
        if  p.crossSectionDR.Coulomb != 0.0   ||   p.interference.Coulomb != 0.0
            success = false
            if printTest   info(iostream, "with includeDR = false, sigma_DR = $(p.crossSectionDR.Coulomb) and " *
                                "interference = $(p.interference.Coulomb); both must be identically zero")   end
        end
        # (2) what remains must BE radiative recombination, evaluated by the other module
        for  l  in  rrLines
            if  abs(l.electronEnergy - p.electronEnergy) < 1.0e-10  &&  l.finalLevel.index == p.finalLevel.index
                ratio = p.crossSectionRR.Coulomb / l.crossSection.Coulomb
                if  abs(ratio - 1.0) > 1.0e-8
                    success = false
                    if printTest   info(iostream, "the RR limit gives $(p.crossSectionRR.Coulomb) against " *
                                        "PhotoRecombination's $(l.crossSection.Coulomb), ratio $ratio")   end
                end
            end
        end
    end

    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_PhotoRecombinationInterference()::", success)
    return(success)
end


"""
`TestFrames.testModule_ResonantImpactIonization(; short::Bool=true)`  ... tests on module ResonantImpactIonization;
    a success::Bool is returned.
"""
function testModule_ResonantImpactIonization(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ResonantImpactIonization-new.sum")
    printstyled("\n\nTest the module  ResonantImpactIonization  ... \n", color=:cyan)
    # The module's own example, example-Fl.jl, is a cascade that takes 72-88 s -- far too heavy for the suite, and
    # its headline check (the 158 branching sums coming out at exactly 1) is one its own comment warns about: it
    # tests the ARITHMETIC and never the completeness, so a decay route left out of the configuration lists is
    # missing from the denominator too and the sum is 1 regardless.
    #
    # What IS checked here are the SCALING LAWS and BOUNDS of the strength formulas, which are exact, which no
    # implementation can satisfy by accident, and which need no computation at all:
    #     S_0 = pi^2/k^2 * A_capture * (2J_d + 1)/(2J_i + 1)
    # is linear in the capture rate, proportional to the statistical-weight ratio, and falls as 1/k^2; while the
    # sequential and simultaneous strengths are S_0 times branchings that cannot exceed 1, so neither may ever
    # exceed S_0 -- more strength into one route than was captured in the first place would be a defect that no
    # comparison against a literature total would necessarily reveal.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    E = 1.0;    A = 1.0e13;    RII = ResonantImpactIonization

    # (1) exact linearity in the capture rate
    s1 = RII.resonanceStrength(E, A, 3, 1);    s2 = RII.resonanceStrength(E, 2A, 3, 1)
    if  abs(s2 - 2*s1) > 1.0e-12 * abs(s2)
        success = false
        if printTest   info(iostream, "resonanceStrength is not linear in the capture rate: $s1 and $s2")   end
    end
    # (2) the statistical-weight ratio (2J_d+1)/(2J_i+1), exactly
    for  (jd2, ji2)  in  [(1,1), (3,1), (5,3), (7,1), (0,4)]
        wa = RII.resonanceStrength(E, A, jd2, ji2);    wb = RII.resonanceStrength(E, A, ji2, ji2)
        if  abs(wa/wb - (jd2 + 1.)/(ji2 + 1.)) > 1.0e-12
            success = false
            if printTest   info(iostream, "weight ratio for 2Jd,2Ji = $jd2,$ji2 is $(wa/wb), must be " *
                                          "$((jd2+1.)/(ji2+1.))")   end
        end
    end
    # (3) the 1/k^2 fall-off. At these energies k^2 = 2E to within the relativistic correction, so the strength
    #     ratio between E and 2E must be 2 to better than a part in 1e3.
    r = RII.resonanceStrength(0.5, A, 3, 1) / RII.resonanceStrength(1.0, A, 3, 1)
    if  abs(r - 2.0) > 1.0e-3
        success = false
        if printTest   info(iostream, "the strength ratio between E and 2E is $r, must be 2 to O(alpha^2 E)")   end
    end
    # (4) a non-positive resonance energy and a vanishing width both give EXACTLY zero
    if  RII.resonanceStrength(0., A, 3, 1) != 0.  ||  RII.resonanceStrength(-1., A, 3, 1) != 0.
        success = false
        if printTest   info(iostream, "resonanceStrength must be exactly 0 at non-positive energy")   end
    end
    if  RII.sequentialStrength(E, A, 3, 1, A, 0., A, A) != 0.  ||  RII.simultaneousStrength(E, A, 3, 1, A, 0.) != 0.
        success = false
        if printTest   info(iostream, "a vanishing total width must give exactly 0")   end
    end
    # (5) NEITHER route may exceed the formation strength, and both saturate it when the branchings are 1
    s0 = RII.resonanceStrength(E, A, 3, 1)
    for  (aToN, gD, aN, gN)  in  [(0.3A, A, 0.5A, A), (A, A, A, A), (0.1A, A, 0.9A, A)]
        wa = RII.sequentialStrength(E, A, 3, 1, aToN, gD, aN, gN)
        if  wa > s0 * (1.0 + 1.0e-12)
            success = false
            if printTest   info(iostream, "sequentialStrength $wa exceeds the formation strength $s0")   end
        end
        if  aToN == gD  &&  aN == gN  &&  abs(wa - s0) > 1.0e-12 * s0
            success = false
            if printTest   info(iostream, "with both branchings 1, sequentialStrength $wa must equal $s0")   end
        end
    end
    for  (aDouble, gD)  in  [(0.25A, A), (A, A), (0.75A, A)]
        wa = RII.simultaneousStrength(E, A, 3, 1, aDouble, gD)
        if  wa > s0 * (1.0 + 1.0e-12)
            success = false
            if printTest   info(iostream, "simultaneousStrength $wa exceeds the formation strength $s0")   end
        end
        if  aDouble == gD  &&  abs(wa - s0) > 1.0e-12 * s0
            success = false
            if printTest   info(iostream, "with a branching of 1, simultaneousStrength $wa must equal $s0")   end
        end
    end
    # (6) shakeProbability is a probability, and it is EXACTLY zero when the electron's own orbital is among the
    #     final ones -- the projection is then complete and nothing can be shaken off.
    grid  = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    orb10 = HydrogenicIon.orbital(Subshell("1s_1/2"), Nuclear.Model(10.), grid)
    orb12 = HydrogenicIon.orbital(Subshell("1s_1/2"), Nuclear.Model(12.), grid)
    p0    = RII.shakeProbability(orb10, Radial.Orbital[orb10], grid)
    if  abs(p0) > 1.0e-6
        success = false
        if printTest   info(iostream, "shakeProbability against the electron's own orbital is $p0, must be 0")   end
    end
    p1 = RII.shakeProbability(orb10, Radial.Orbital[orb12], grid)
    if  !(0. <= p1 <= 1.)
        success = false
        if printTest   info(iostream, "shakeProbability = $p1 is outside [0,1]")   end
    end

    println(iostream, "ResonantImpactIonization: the strength formulas by their exact scaling laws -- linearity " *
                      "in the capture rate, the statistical-weight ratio, 1/k^2 -- and the bound that neither "   *
                      "route may exceed the formation strength. No cascade and no literature value.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_ResonantImpactIonization()::", success)
    return( success )
end


"""
`TestFrames.testModule_TwoElectronOnePhoton(; short::Bool=true)`  ... tests on module TwoElectronOnePhoton;
    a success::Bool is returned.
"""
function testModule_TwoElectronOnePhoton(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-TwoElectronOnePhoton-new.sum")
    printstyled("\n\nTest the module  TwoElectronOnePhoton  ... \n", color=:cyan)
    # A TWO-ROUTE test, cut down from example-Dn.jl branch a). The same physical rate for Be-like Ne is reached in
    # two independent ways, and to leading order in the configuration mixing they must agree:
    #
    #   (ii) EXPLICIT second-order perturbation theory -- this module. The initial state is the PURE 2s2p
    #        configuration and the coupling to the one-hole intermediates is carried by an explicit sum over a
    #        Green-function multiplet, with the energy denominators written out.
    #   (i)  IMPLICIT, through ordinary CI mixing -- a plain PhotoEmission computation whose initial multiplet
    #        already contains 1s2p beside 2s2p, so the one-body E1 operator connects the small 1s2p admixture to
    #        the ground state. The two multiplets are built on different, non-orthogonal orbital sets, hence
    #        calcBiorthogonal = true.
    #
    # The agreement is not a coincidence to be tuned: the CI coefficient of route (i) IS c_n = <n|V|i>/(E_i - E_n),
    # the very ratio route (ii) forms inside its second-order sum. So this is a known-answer test in which the
    # answer is supplied by a different part of the code, and neither route can be fixed to match the other by
    # adjusting a stored number, because there is none.
    #
    # The COULOMB gauge is compared. Babushkin is deliberately not asserted: these amplitudes are small and
    # near-cancelling, the two routes truncate the cancellation differently, and the length form is the more
    # sensitive of the two -- example-Dn.jl records the weak line's Babushkin rate moving 37 % under a change that
    # left Coulomb at 0.14 %. A tolerance loose enough to admit that would test nothing.
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # hp = 1.2e-2 where example-Dn.jl uses 0.6e-2, and the coarser grid costs nothing in accuracy: MEASURED at
    # hp = 0.6e-2 / 1.2e-2 / 2.0e-2 (1967 / 1134 / 798 points) the ratio below comes out 1.0516 / 1.0517 / 1.0517
    # and the TEOP rate 1.53476e10 / 1.53477e10 / 1.53478e10 1/s -- five figures unchanged over a factor 2.5 in
    # grid size, and the finest reproduces example-Dn.jl's recorded 1.534757e10 exactly. Grid-insensitivity on
    # that scale is itself evidence the comparison is converged rather than balanced on a resolution.
    #
    # NO TIMING CLAIM IS MADE FOR THE GRID CHANGE. The obvious experiment -- loop over the three grids in one
    # session -- reports 74.5 s, 3.6 s and 2.3 s, which looks like a twentyfold saving and is not: the FIRST
    # iteration pays for compiling the TEOP and biorthogonal paths and the later ones run warm. Measured cold, the
    # same test takes 69.7 s at 0.6e-2 and 56.6 s at 1.2e-2, so the grid is worth about 13 s and compilation is
    # the rest. Inside the suite most of that compilation is shared with other tests and is not paid again here.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.2e-2, rbox = 10.0)
    finalConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("2p^2")]

    gComp      = Atomic.Computation(Atomic.Computation(), name="TEOP test: Green-function multiplet", grid=grid,
                                    nuclearModel = Nuclear.Model(10.),
                                    configs = [Configuration("1s 2p"), Configuration("1s 2s")] )
    gMultiplet = perform(gComp; output=true)["multiplet:"]

    # route (ii): the module itself
    teopSettings = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,UseBabushkin], false,
                                                 LineSelection(true, indexPairs=[(i,1) for i = 1:4]), 0.,
                                                 CoulombInteraction(), gMultiplet)
    wb = Atomic.Computation(Atomic.Computation(), name="TEOP test: explicit second order", grid=grid,
                            nuclearModel   = Nuclear.Model(10.),
                            initialConfigs = [Configuration("2s 2p")],
                            finalConfigs   = finalConfigs,      processSettings = teopSettings )
    teopLines = perform(wb; output=true)["two-electron-one-photon lines:"]

    # route (i): CI mixing plus the biorthogonal transformation
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin],
                        printBefore=false, calcBiorthogonal=true,
                        lineSelection=LineSelection(true, indexPairs=[(i,1) for i = 1:8]) )
    wc = Atomic.Computation(Atomic.Computation(), name="TEOP test: CI mixing + biorthogonal", grid=grid,
                            nuclearModel   = Nuclear.Model(10.),
                            initialConfigs = [Configuration("2s 2p"), Configuration("1s 2p")],
                            finalConfigs   = finalConfigs,      processSettings = photoSettings )
    photoLines = perform(wc; output=true)["radiative lines:"]

    # THE TWO ROUTES ARE MATCHED BY PHOTON ENERGY, NOT BY SIZE, and that is the one subtlety of this test. Route
    # (i)'s multiplet contains the 1s2p levels as well as the 2s2p ones, so it also produces the ORDINARY Kalpha
    # lines -- at ~0.92 keV and A ~ 9.9e12 1/s, three orders above anything TEOP. Taking the largest rate finds
    # those, compares an allowed one-electron transition against a second-order one, and fails by a factor of 600.
    # The TEOP counterpart sits near 1.93 keV instead, and the two routes describe the SAME transition, so their
    # photon energies must very nearly coincide -- which is asserted here before the rates are compared at all.
    if  length(teopLines) == 0  ||  length(photoLines) == 0
        success = false
        if printTest   info(iostream, "one of the two routes returned no line at all: $(length(teopLines)) and " *
                                      "$(length(photoLines))")   end
    else
        teopLine = teopLines[argmax([line.teopRate.Coulomb  for line in teopLines])]
        toGround = [line  for line in photoLines  if line.finalLevel.index == 1]
        photoLine = toGround[argmin([abs(line.omega - teopLine.omega)  for line in toGround])]
        dOmega    = abs(photoLine.omega - teopLine.omega) / teopLine.omega
        if  dOmega > 0.01
            success = false
            if printTest   info(iostream, "no route-(i) line matches the TEOP photon energy $(teopLine.omega); " *
                                          "closest is $(photoLine.omega), off by $dOmega")   end
        end
        aTeop  = teopLine.teopRate.Coulomb
        aPhoto = photoLine.photonRate.Coulomb
        ratio  = aTeop / aPhoto
        # (1) the two routes agree. 15 % admits the 5 % of example-Dn.jl with room for the truncation to differ,
        #     and still excludes everything this test is for: a wrong second-order denominator, a missing time
        #     ordering or a lost factor would move this by an order of magnitude, not by a tenth.
        if  !(0.85 < ratio < 1.15)
            success = false
            if printTest   info(iostream, "second order gives $aTeop 1/s, CI + biorthogonal $aPhoto 1/s, " *
                                          "ratio $ratio; the two routes must agree to leading order")   end
        end
        # (2) and it really is a SECOND-order rate: an allowed one-electron E1 line at Z = 10 runs at 1e13-1e14
        #     per second, so a TEOP rate anywhere near that would mean the amplitude is not the one intended.
        #     NOTE THE UNITS. `teopRate` and `photonRate` are stored in ATOMIC UNITS; the 1/s figures in the
        #     printed tables and in example-Dn.jl are converted on the way out. 1 a.u. of rate is 4.1341e16 /s,
        #     so the 1.5348e10 /s of that example is 3.712e-7 here. The ratio in check (1) is unaffected, being
        #     a ratio; this bound is not, and was written in 1/s at first, which is what the restored diagnostic
        #     caught in one line.
        aTeopPerSec = Defaults.convertUnits("rate: from atomic to 1/s", aTeop)
        if  !(1.0e6 < aTeopPerSec < 1.0e12)
            success = false
            if printTest   info(iostream, "the TEOP rate $aTeopPerSec 1/s is outside the range of a " *
                                          "second-order process at Z = 10")   end
        end
    end

    println(iostream, "TwoElectronOnePhoton: the same Be-like Ne rate by explicit second-order perturbation " *
                      "theory and by CI mixing with a biorthogonal transformation; the two must agree to "    *
                      "leading order in the mixing. Coulomb gauge only. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_TwoElectronOnePhoton()::", success)
    return( success )
end
