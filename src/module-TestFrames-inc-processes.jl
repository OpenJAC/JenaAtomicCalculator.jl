
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
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    thetas = [0.3, 1.0, 2.0, 3.0]
    psSet  = ParticleScattering.Settings(ParticleScattering.Settings(), impactEnergies=[300.0], polarThetas=thetas,
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
                                        ParticleScattering.StaticFieldExchange(), Beam.PlaneWave(),
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
    ## Test 8: an ABSOLUTE check against an exact analytic limit, and the only one here that is not internal. At a large
    ##   scattering angle and a high impact energy the projectile probes inside the electron cloud, where the potential is
    ##   the bare nuclear one, so the backward cross section must approach Rutherford for the full nuclear charge,
    ##   d sigma/d Omega -> Z^2 / (16 E^2 sin^4(theta/2)). Measured for e + He at 1 keV and theta = 180 deg, the ratio is
    ##   1.069: above unity because screening, exchange and the finite electron cloud all still raise it there, and by a
    ##   margin that shrinks with energy (1.355 at 200 eV, 1.003 at 2 keV; see examples/example-Ob.jl). Both bounds below
    ##   are physical, and together they pin the absolute normalisation of the amplitudes, which no internal check can do.
    ruthSet = ParticleScattering.Settings(ParticleScattering.Settings(), impactEnergies=[1000.0],
                                          polarThetas=[Float64(pi)], polarPhis=[0.0], printBefore=false,
                                          epsPartialWave=1.0e-7, maxL=150)
    wcr = Atomic.Computation(Atomic.Computation(), name="Rutherford limit", grid=grid, nuclearModel=Nuclear.Model(2.0),
                             initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s^2")],
                             processSettings = ruthSet )
    wdr = redirect_stdout(devnull) do
        perform(wcr, output=true)
    end
    evr   = wdr["particle-scattering events:"][1]
    ratio = evr.angular[1].dcs / ( 4.0 / (16 * evr.impactEnergy^2) )
    success = success && 1.0 < ratio < 1.15
    #
    Defaults.setDefaults("unit: energy", oldEnergyUnit)
    Defaults.setDefaults("print summary: close", "")
    _, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "ParticleScattering: kappa bookkeeping, partial-wave convergence, the two routes to sigma_el, " *
                      "the transport cross sections, the Sherman function, and the two guards.")
    testPrint("testModule_ParticleScattering()::", success)
    return(success)
end
