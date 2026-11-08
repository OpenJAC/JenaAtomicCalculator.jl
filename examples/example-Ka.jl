#
println("Ka) Train & apply a neural network for estimating missing energies in the NIST tables.")


if  true
    # Last visit:  07-Nov-2026
    # CAUTION, 08-Nov-2026 (priority list item 11): generateAtomicModelForLE_Arn4()'s grid (rbox=10.0 a.u.,
    # r_max ~ 13 a.u. once built) flags 4s/4p/4d/4f (Z_eff ~ 2) as BOX TOO SMALL -- ~96 a.u. would fit them
    # comfortably. ANSWERED DIRECTLY by the ONE-SCF check item 77/86 prescribes, not by the heuristic alone:
    # at the CURRENT box the mean field for Ar+ converges cleanly and every spin-orbit pair is physically
    # sound (4s=-0.3068, 4p=-0.2104/-0.2096, 4d=-0.09409/-0.09404, 4f=-0.03215/-0.03223 Ha -- correctly
    # ordered by n and l, correctly small splittings, no inversion). THE OBVIOUS "FIX" BREAKS IT: widening to
    # rbox=40 makes `SelfConsistent.performSCF`'s own `Bsplines.checkOrbitalConsistency` REFUSE outright --
    # 4f_5/2/4f_7/2 land on a 5.235x energy ratio, i.e. two DIFFERENT states, not two fine-structure partners
    # -- reproducing `example-Ac.jl`'s own "wrong turn (b)" (box too large starves this basis) in an
    # UNRELATED system, which is why that trap should be read as general rather than Ac.jl-specific. THE
    # WARNING IS A FALSE POSITIVE; the box stays as it is. (This is independent of why `Last successful`
    # below is "unknown" -- this branch needs real NIST input files it does not have, `filenames = String[]`.)
    # Last successful:  unknown ...
    # Compute the atomic model for Ar^+, Ar^2+, Ar^3+ ions and with nMax = 4;
    # extract the NIST levels from externally prepared files
    filenames   = String[]
    nistLevels  = DeepLearning.extractNistLevels(filenames)
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
    # Last visit:  07-Nov-2026   (NEEDS AN EARLIER BRANCH OF THIS FILE)
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ...
    # Run the network for one or several configurations and compared with predicted excitation energies
    # with the data available from the NIST tables.
    # `nistLevels` and `atomicModel` are built by branch a of this file (lines 9-12); without them this branch
    # died on an UndefVarError naming only the first of them.
    if  !@isdefined(nistLevels)   error("Run branch a of this file first; it builds nistLevels and atomicModel.")   end
    requestConfs = [Configuration("[Ne] 3s^2 3p^3 4f")]
    request      = DeepLearning.LevelEstimationRequest(requestConfs, nistLevels)
    application  = DeepLearning.Application(Application(), name="Ar-NN with nMax=4", atomicModel=atomicModel,
                                            request=request)
                                           
    Basics.run(application)
    #
end
    
    
