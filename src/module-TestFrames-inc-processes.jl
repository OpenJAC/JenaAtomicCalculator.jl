
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
                                joinpath(@__DIR__, "..", "test", "test-AutoIonization-new.sum"), "Auger rates and intrinsic", 5)
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
    ## Check the summed (grand-total) cross sections separately: they must stay consistent with the
    ## line-resolved table above, of which they are the sum over all final levels.
    success = success  &&
              testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoIonization-approved.sum"),
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"),
                                "Total photoionization cross sections, summed", 15)
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
                            nuclearModel   = Nuclear.Model(90.0, "Fermi"),
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
