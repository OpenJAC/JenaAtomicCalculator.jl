using Test
using JenaAtomicCalculator, ..Defaults, ..TestFrames

## Clear the .sum files of the PREVIOUS run, before this one starts.  Each test that opens a summary writes a
## test-<Name>-new.sum here and never removes it, so they accumulate -- 51 had collected by 29-Aug-2026.
## Deliberately done at the START and not at the end, unlike the zzz-*.jld block below: after a FAILING run those
## files are the evidence of what the test actually produced, and deleting them on the way out would throw away
## the one thing worth reading.  Cleaning on the way in gives the same tidiness and keeps the evidence.
## test/approved/ is a subdirectory and readdir() does not descend, so the 33 approved references are out of reach.
let  dir = @__DIR__
    stale = filter(f -> startswith(f, "test-") && endswith(f, "-new.sum"), readdir(dir))
    for  f in stale   rm(joinpath(dir, f), force = true)   end
    if  length(stale) > 0
        printstyled("\nCleared $(length(stale)) .sum file(s) from the previous run in $(dir) \n", color = :cyan)
    end
end

@testset "Name" begin
    printstyled("\nPerform tests on the JAC program; this may take a while .... \n", color=:cyan)
    # To see WHY a test failed, uncomment the next line: it turns on the per-test diagnostics that every
    # `if printTest   info(iostream, "...")   end` in TestFrames writes.  They are off by default because a
    # passing suite has nothing to say.
    #     Defaults.setDefaults("print test: open", joinpath(@__DIR__, "runtests.report"))
    # The line that stood here until 29-Aug-2026 read `Defaults.Constants.define(...)`, which is not an API
    # this package has -- there is no Constants submodule and no `define` -- so uncommenting it would have
    # RAISED rather than enabled anything.  Nobody found out, because nobody uncommented it; and the 77
    # diagnostics behind it were themselves calling an `info` that was only defined on 28-Aug-2026.

    @testset "JAC methods" begin
        @test TestFrames.testMethod_Wigner_3j()
        @test TestFrames.testMethod_HydrogenicRates()
        @test TestFrames.testMethod_OrbitalOrthonormality()
        @test TestFrames.testMethod_BreitInteraction()
        @test TestFrames.testMethod_Opacities()
        @test TestFrames.testMethod_SpinAngular()
        @test TestFrames.testMethod_DocstringPointers()
        @test TestFrames.testModule_AngularMomentum()   ## the four added 28-Aug-2026 are closed-form or
        @test TestFrames.testModule_HydrogenicIon()     ## identity checks: inverses, limits, orthogonality
        @test TestFrames.testModule_Nuclear()           ## and exact hydrogenic values -- so none of them
        @test TestFrames.testModule_RadialIntegrals()   ## can pass on a stale stored reference
        @test TestFrames.testModule_Bsplines()      ## added 28-Aug-2026; exercises BOTH Rule 12 guards from
                                                    ## both sides -- they must refuse a 0.05 a.u. box
        @test TestFrames.testModule_StarkZeeman()   ## added 28-Aug-2026; the centre of gravity cannot move
        @test TestFrames.testModule_Hamiltonian()          ## added 28-Aug-2026; trace, the variational bound,
        @test TestFrames.testModule_InteractionStrength()  ## and two quadratures of one operator
        @test TestFrames.testModule_SelfConsistent()  ## added 28-Aug-2026; the converged SCF is a FIXED POINT,
                                                      ## and a frozen subshell does not move at all
    end

    @testset "JAC structs" begin
        @test TestFrames.testStructConstructors()
        @test TestFrames.testMethod_SettingsCopyConstructors()
    end

    @testset "JAC evaluations" begin
        @test TestFrames.testEvaluation_Wigner_3j_specialValues() 
        @test TestFrames.testEvaluation_Wigner_6j_specialValues() 
        @test TestFrames.testEvaluation_Wigner_9j_specialValues() 
        @test TestFrames.testEvaluation_sumRulesForOneWnj() 
        @test TestFrames.testEvaluation_sumRulesForTwoWnj() 
        @test RacahAlgebra.testSpecialValuesW3j()
        @test RacahAlgebra.testSpecialValuesW6j()
        @test RacahAlgebra.testSpecialValuesW9j()
        @test RacahAlgebra.testSumRules()
    end

    @testset "JAC representations" begin
        @test TestFrames.testRepresentation_MeanFieldBasis_CiExpansion() 
        @test TestFrames.testRepresentation_RasExpansion() 
        @test TestFrames.testRepresentation_GreenExpansion() 
    end

    @testset "JAC amplitudes" begin
        @test TestFrames.testModule_MultipoleMoment() 
        @test TestFrames.testModule_WeakInteractionMoment() 
    end

    @testset "JAC properties" begin
        @test TestFrames.testModule_Einstein()
        @test  TestFrames.testModule_Hfs()             ## re-enabled 16-Aug-2026; the blocker was calcHfMultiplet, not calcNondiagonal
        @test TestFrames.testModule_LandeZeeman() 
        @test TestFrames.testModule_IsotopeShift()   
        @test TestFrames.testModule_AlphaVariation() 
        @test TestFrames.testModule_FormFactor() 
        @test TestFrames.testModule_DecayYield()
        @test TestFrames.testModule_MultipolePolarizibility()
        @test TestFrames.testModule_ReducedDensityMatrix()  ## added 28-Aug-2026; algebraic invariants of
                                                            ## rho^(1p), no stored .sum
    end

    @testset "JAC processes" begin
        @test TestFrames.testModule_PhotoEmission()
        @test TestFrames.testModule_PhotoExcitation()
        @test TestFrames.testModule_PhotoIonization()
        @test TestFrames.testModule_PhotoRecombination()
        @test TestFrames.testModule_AutoIonization()  
        @test TestFrames.testModule_DielectronicRecombination()  
        @test TestFrames.testModule_HyperfineInduced()
        @test TestFrames.testModule_RayleighCompton()
        @test TestFrames.testModule_MultiPhotonTransition()
        @test TestFrames.testModule_CoulombExcitation()
        @test TestFrames.testModule_GeneralizedOscillatorStrength()   ## added 28-Aug-2026; sanity only --
        @test TestFrames.testModule_PhotoRecombinationInterference()  ## absolute checks, no stored .sum
        @test TestFrames.testModule_ParticleScattering()   ## added 17-Aug-2026 with the Dirac rebuild; replaces a
                                                           ## Settings() constructor entry that could not fail
        @test TestFrames.testModule_TwoElectronOnePhoton()      ## added 28-Aug-2026: a two-route agreement, and
        @test TestFrames.testModule_ResonantImpactIonization()  ## exact scaling laws -- no stored .sum in either
    end

    @testset "JAC cascades" begin
        @test TestFrames.testModule_Cascade_StepwiseDecay()
        @test TestFrames.testModule_Cascade_PhotonIonization()
        @test TestFrames.testModule_Cascade_PhotonExcitation()
        @test TestFrames.testModule_Cascade_PhotoAbsorption()
        @test TestFrames.testModule_Cascade_DielectronicCapture()
        @test TestFrames.testModule_Cascade_ResonantIonization()
        @test TestFrames.testModule_Cascade_EiiRateCoefficients()
        @test  TestFrames.testModule_Cascade_Simulation()   ## re-enabled 16-Aug-2026; data file regenerated and the Simulation API call brought up to date
    end

    @testset "JAC empirical" begin
        @test TestFrames.testModule_Empirical()
        @test TestFrames.testModule_ImpactIonization()
        @test TestFrames.testModule_Semiempirical()
    end

    @testset "JAC plasma" begin
    end

    @testset "JAC strongfield" begin
        ## no tests defined yet
    end

    @testset "JAC Liouville" begin
        @test TestFrames.testModule_Liouville()
    end

    @testset "JAC DeepLearning" begin
        @test TestFrames.testModule_DeepLearning()
    end

end

## Clean the run debris.  Each cascade computation dumps its results to a zzz-*.jld file whose name carries the date and the
## hour, so the suite leaves one file behind per scheme and per hour it is run in, without bound -- 109 files and 51 MB had
## accumulated since 5-Aug-2026.  The dumps exist so that a single cascade can be reloaded by hand (the schemes print the
## JLD2.load call for it); nothing reads them back, and least of all this suite.  The glob is anchored on the zzz- prefix and
## the .jld suffix and does not descend, so test/approved/ -- which holds a .jld reference -- is out of its reach.
let  dir = @__DIR__
    debris = filter(f -> startswith(f, "zzz-") && endswith(f, ".jld"), readdir(dir))
    bytes  = sum(f -> filesize(joinpath(dir, f)), debris; init = 0)
    for  f in debris   rm(joinpath(dir, f), force = true)   end
    printstyled("\nRemoved $(length(debris)) cascade dump file(s), $(round(bytes / 2^20, digits = 1)) MB, from $(dir) \n",
                color = :cyan)
end
