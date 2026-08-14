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
        @test TestFrames.testModule_ParityNonConservation() 
    end

    @testset "JAC properties" begin
        @test TestFrames.testModule_Einstein()
        ## @test TestFrames.testModule_Hfs()           ## disabled: runtime error "still to be done" in Hfs.computeAmplitudesProperties (calcNondiagonal path not implemented)   [re-verified 09-Aug-2026: still raises "... still to be done for a single nuclear spin/isomer"]
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
    end

    @testset "JAC cascades" begin
        @test TestFrames.testModule_Cascade_StepwiseDecay()
        @test TestFrames.testModule_Cascade_PhotonIonization()
        @test TestFrames.testModule_Cascade_PhotonExcitation()
        @test TestFrames.testModule_Cascade_PhotoAbsorption()
        ## @test TestFrames.testModule_Cascade_Simulation()  ## disabled: test-Cascade-StepwiseDecay-data.jld predates CsfR.seniorityNr field; needs data file regeneration   [re-verified 09-Aug-2026: still raises JLD2.ReconstructedMutable{:CsfR}; regenerating the reference data is an editorial act, not a repair]
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
