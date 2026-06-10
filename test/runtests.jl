using Test
using JenaAtomicCalculator, ..Defaults, ..TestFrames

@testset "Name" begin
    printstyled("\nPerform tests on the JAC program; this may take a while .... \n", color=:cyan)
    ## Defaults.Constants.define("print test: open", pwd() * "/runtests.report")

    @testset "JAC methods" begin
        @test TestFrames.testMethod_Wigner_3j()
    end

    @testset "JAC structs" begin
        @test TestFrames.testStructConstructors()
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
        ## @test TestFrames.testModule_PlasmaShift()  ## disabled: PlasmaShift is now a Plasma.Computation scheme, not an Atomic property; needs Plasma module rework
        ## @test TestFrames.testModule_Hfs()           ## disabled: runtime error "still to be done" in Hfs.computeAmplitudesProperties (calcNondiagonal path not implemented)
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
        ## @test TestFrames.testModule_PhotoExcitationFluores()  ## disabled: Atomic.Computation constructor signature changed; test call needs updating
        ## @test TestFrames.testModule_PhotoExcitationAutoion()  ## disabled: PhotoExcitationAutoion.Settings positional call outdated; test call needs updating
        @test TestFrames.testModule_RayleighCompton()
        @test TestFrames.testModule_MultiPhotonDeExcitation()
        @test TestFrames.testModule_CoulombExcitation()
    end

    @testset "JAC cascades" begin
        @test TestFrames.testModule_Cascade_StepwiseDecay()
        ## @test TestFrames.testModule_Cascade_PhotonIonization()
        ## @test TestFrames.testModule_Cascade_PhotonExcitation()
        ## @test TestFrames.testModule_Cascade_Simulation()
    end

end
