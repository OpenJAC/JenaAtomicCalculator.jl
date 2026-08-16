#
println("Ka) Train & apply a neural network for estimating missing energies in the NIST tables.")


if  true
    # Last successful:  unknown ...
    # Compute the atomic model for Ar^+, Ar^2+, Ar^3+ ions and with nMax = 4;
    # extract the NIST levels from externally prepared files
    filenames   = String[]
    nistLevels  = DeepLearning.extractNistLevels(filesnames)
    nistConfs   = DeepLearning.extractNistConfigurations(nistLevels)
    atomicModel = DeepLearning.generateAtomicModelForLE_Arn4()
    #
elseif  true
    # Last successful:  unknown ...
    # Generate the x- (and y-) feature vectors for the training of a neutral network (NN);
    # write them out in a format suitable for externals NN tools
    testProbability   = 0.1  # Select how many test feature vectors should be separated from the whole set
    (trainXy, testXy) = DeepLearning.generateFeatureVectors(nistConfs, atomicModel, testProbability, nistLevels)
    DeepLearning.writeFeatureVectors(trainXy, "nn.training-data-$(today()).txt")
    DeepLearning.writeFeatureVectors(testXy,  "nn.test-data-$(today()).txt")
    #
elseif  true
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ...
    # Run the network for one or several configurations and compared with predicted excitation energies
    # with the data available from the NIST tables.
    requestConfs = [Configuration("[Ne] 3s^2 3p^3 4f")]
    request      = DeepLearning.LevelEstimationRequest(requestConfs, nistLevels)
    application  = DeepLearning.Application(Application(), name="Ar-NN with nMax=4", atomicModel=atomicModel,
                                            request=request)
                                           
    Basics.run(application)
    #
end
    
    
