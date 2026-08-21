using Test
using JenaAtomicCalculator, ..Defaults, ..TestFrames

@testset "Name" begin
    printstyled("\nPerform tests on the JAC program; this may take a while .... \n", color=:cyan)
    ## Defaults.Constants.define("print test: open", pwd() * "/runtests.report")

    @testset "JAC methods" begin
        @test TestFrames.testMethod_Wigner_3j()
        @test TestFrames.testMethod_HydrogenicRates()
        @test TestFrames.testMethod_OrbitalOrthonormality()
        @test TestFrames.testMethod_BreitInteraction()
        @test TestFrames.testMethod_Opacities()
        @test TestFrames.testMethod_DocstringPointers()
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
        @test TestFrames.testModule_ParticleScattering()   ## added 17-Aug-2026 with the Dirac rebuild; replaces a
                                                           ## Settings() constructor entry that could not fail
    end

    @testset "JAC cascades" begin
        @test TestFrames.testModule_Cascade_StepwiseDecay()
        @test TestFrames.testModule_Cascade_PhotonIonization()
        @test TestFrames.testModule_Cascade_PhotonExcitation()
        @test TestFrames.testModule_Cascade_PhotoAbsorption()
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
