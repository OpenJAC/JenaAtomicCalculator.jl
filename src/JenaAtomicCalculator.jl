#
#  Use git:                    git status ....;   git add <filenames>;   git commit -m "..";   git push;   git rm  <filenames>
#  Use Jupyter notebooks:      using IJulia;   notebook()
#  Activation:                 ];   pkg> up;   pkg> activate
#  Working with JAC:           using Revise;   using JAC;   include("../src/jac.jl");   pkg> test
#  
#  Copy to desktop             scp -r JAC.jl/ fritzsch@10.140.119.236:~/fri/.
"""
`module JenaAtomicCalculator`  
    ... Jena Atomic Calculator (JAC) provides tools for performing atomic (structure) calculations at various degrees of complexity 
        and sophistication. It has been designed to not only calculate atomic level structures and properties [such as g-factors or
        hyperfine and isotope-shift parameters] but also transition amplitudes between bound-state levels [for the dipole 
        operator, etc.] and, in particular, (atomic) transition probabilities, Auger rates, photoionization cross sections, 
        radiative and dielectronic recombination rates as well as cross sections for several other -- elementary or composed --
        processes. 

"""
module JenaAtomicCalculator

const JAC = JenaAtomicCalculator


# Restrict the size and functionality of code by just including certain modules, while others are not taken into account.
# Of course, this "exclusion" may cause later errors if missing functionality is invoked.  The selective use of code
# is introduced mainly for development purposes and for keeping the storage requirement moderate (if one wishes to focus
# on certain classes of applications)
incProperties           = true
incBasicProcesses       = true
incAdvancedProcesses    = true
incCascades             = true  ## Requires: incBasicProcesses
incPlasma               = true  ## Requires: incProperties
incLiouville            = true  ## Requires: incProperties
incStrongField          = true
incAtomicCompass        = true
incRacahAlgebra         = true
incDeepLearning         = true

#==
incAdvancedProcesses    = false
incCascades             = false  ## Requires: incBasicProcesses
incPlasma               = false  ## Requires: incProperties
incLiouville            = false  ## Requires: incProperties
incStrongField          = false
incAtomicCompass        = false
incRacahAlgebra         = false  ==#



using  Dates,  Printf,  BSplineKit, LinearAlgebra, SpecialFunctions, QuadGK, Cubature, GSL, JLD2, SymEngine, 
       HypergeometricFunctions  ## , Interact, GaussQuadrature, IJulia, FortranFiles

export AbstractCImethod, AbstractComputeTheme, AbstractConfigurationRestriction, AbstractConfigurationTheme,
       AbstractDisplayTheme, AbstractEeInteraction, AbstractEmpiricalSettings, AbstractLineShiftSettings,
       AbstractNeutralNetwork, AbstractNeutralNetworkRequest, AbstractNuclearModel, AbstractPlasmaModel, AbstractPlotTheme,
       AbstractPotential, AbstractProcessSettings, AbstractPropertySettings, AbstractQedModel, AbstractStartOrbitals,
       add, AddElectrons, AllShells, AlphaVariation, analyze, AnapoleMoment, AngularCoeffs1pGrasp92,
       AngularCoeffs1pRatip2013, AngularCoeffsEeRatip2013, AngularJ, AngularJ64, AngularM64, AngularMomentum,
       Application, AsfSettings, Atomic, AtomicCompass, AtomicFeatures, AtomicModel, AtomicState, AtomicStructure,
       Auger, AugerInPlasma, AutoIonization, AverageAtom,
       Basics, Basis, Beam, BeamPhotoExcitation, BiOrthogonal, BreitInteraction, Bsplines, ByMultipoles, ByNumber, ByParity,
       Cartesian2DFieldVector, Cartesian3DFieldVector, CartesianVector, Cascade, checkConfigurations, CiExpansion,
       CImatrixWithSymmetryJP, CiSettings, ClebschGordan, CloseCoupling, ClosedCore, ClosedShells, ClosedSubshells,
       Compton, compute, computeBranchingFractions, computeChargeStateDistribution, computeCrossSections,
       computeForPedestrians, computeLevelEnergies, computeLifetimes, computeResonanceStrength,
       computeTransitionRates, Configuration, ConfigurationR, Continuum, ContractShells, convertUnits,
       CorePolarization, Coulex, Coulion, CoulombBreit, CoulombExcitation, CoulombGaunt, CoulombInteraction,
       CoulombIonization, CrystalField, CrystalFieldEmission, CsfR, CurrentSettings,
       DecayYield, DeepLearning, DefaultQuantizationAxis, Defaults, DiagonalCoulomb, diagonalize, DielectronicRecombination,
       Dierec, displayConfiguration, displayConfigurations, displayCouplings, displaySpectrum, Distribution, Djpq,
       DoubleAuger, DoubleAutoIonization,
       E1, E2, E3, E4, Einstein, ElecCapture, ElectricDipoleMoment, ElectronCapture, EmMultipole, Empirical, EmProperty,
       EmPropertyC, estimateCrossSections, evaluate, ExcitationLevel, ExciteElectrons, ExpandShells,
       ExpStokes, extractConfiguration, extractConfigurations, extractFromConfiguration, extractFromConfigurations,
       FermiNucleus, FineStructure, FineStructureLS, ForAutoIonization, ForDielectronicCapture, ForDielectronicRecombination,
       ForElectronCapture, ForGivenConfigs, ForHollowIons, ForImpactIonization, ForIsoelectronicSequence,
       FormFactor, ForPedestrians,
       ForPhotoEmission, ForPhotoIonization, ForPhotoRecombination, ForRasExcitations, ForStepwiseDecay, FromBasis,
       FromMultiplet, FullCIeigen,
       GeneralizedConfigurations, generate, generateConfigurations, getDefaults, GetParity, Green,
       GreenChannel, GreenExpansion, GreenSettings, GroundConfiguration, Gui,
       Hamiltonian, HarmonicQuantizationAxis, Hfs, HighHarmonic, HundsRules, HydrogenicIon, HyperfineInduced,
       HyperfineStructure,
       ImpactExc, ImpactExcAuto, ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, Integral, integrate,
       InteractionStrength, InternalConv, InternalConversion, InternalRecombination, interpolate, IsOccupied,
       IsotopeShift, IsotopicFraction,
       Kronecker,
       LandeF, LandeJ, LandeZeeman, LeadingConfiguration, LeadingConfigurationR, LeftCircular, Level, LevelSelection,
       LevelSymmetry, LineSelection, Liouville, LSjj, LSjjSettings,
       M1, M2, M3, M4, ManyElectron, MeanConfiguration, MeanFieldBasis, MeanFieldMultiplet, MeanFieldSettings,
       MeanOccupation, minus, Model, MultiPDI, MultiPhotonDE, MultiPhotonDoubleIon, MultiPhotonIonization,
       MultiPhotonTransition, MultiPI, Multiplet, Multiplicity, MultipoleMoment, MultipolePolarizibility,
       NoAmplitude, NoneQed, NonrelativisticBasis, NoProcess, NoProperty, Nuclear, NumberOfElectrons,
       OccupationDifference, OneElectronSettings, OneElectronSpectrum, OpenShellNumber, OpenShells, OpenSubshells, oplus,
       Orbital,
       PairA1P, PairAnnihilation1Photon, PairAnnihilation2Photon, PairProduction, Parity, ParityNonConservation,
       ParticleScattering, PathwaySelection, perform, PeriodicTable, Photo, PhotoDouble, PhotoDoubleIonization,
       PhotoEmission, PhotoExc, PhotoExcAuto, PhotoExcFluor, PhotoExcitation, PhotoExcitationAutoion,
       PhotoExcitationFluores, PhotoIonAuto, PhotoIonFluor, PhotoIonization, PhotoIonizationAutoion,
       PhotoIonizationFluores, PhotoRecombination, PhysicalConstants, Plasma, plus, PointNucleus, PrintWarnings, Pulse,
       QedPetersburg, QedSydney,
       RacahAlgebra, RacahExpression, Radial, RadialIntegrals, RadialOrbitalBunge1993, RadialOrbitalHydrogenic,
       RadialOrbitalMcLean1981, RadialOrbitalsBoth, RadialOrbitalsLarge, RadialOrbitalsSmall,
       RadialOrbitalThomasFermi, RadialPotentials, Radiative, RadiativeAuger, RadiativeOpacity, RasExpansion,
       RasLayer, RasSettings, RasStep, RAuger, RayleighCompton, READI, Rec, recast, REDA, ReducedDensityMatrix,
       RelativisticConfigurations, RemoveElectrons, Representation, RequestMaximumOccupation,
       RequestMinimumOccupation, ResonantInelastic, RestrictExcitations, RestrictMaximumDisplacements,
       RestrictNoElectronsTo, RestrictParity, RestrictToShellDoubles, run,
       SchiffMoment, SelfConsistent, Semiempirical, setDefaults, Shell, ShellSelection, SolidAngle, Spectroscopy,
       SphericalMesh, SphericalTensor, SpinAngular, SplitByEnergy, StarkShift, StarkZeeman, StartFromHydrogenic,
       StartFromPrevious, StaticField,
       StaticQuantizationAxis, StrongField, StrongField2, Subshell, SuperConfiguration,
       tabulate, TestFrames, TimeHarmonicField, tools, TotalAM, Triangle, TwoElectronOnePhoton,
       UniformNucleus, UseBabushkin, UseCoulomb, UseGauge,
       ValenceOccupation, ValenceShells,
       W3j, W6j, W9j, WeightedCartesian,
       Ylm
     
# Basic data and data structures
include("module-Basics.jl");            using ..Basics
include("module-Radial.jl");            using ..Radial
include("module-Math.jl");              using ..Math
include("module-Defaults.jl");          using ..Defaults
include("module-Distribution.jl");      using ..Defaults
include("module-ManyElectron.jl");      using ..ManyElectron
include("module-Nuclear.jl");           using ..Nuclear


# Specialized functions/methods to manipulate these data
include("module-AngularMomentum.jl")
## include("module-AngularCoefficients-Ratip2013.jl")  ## keep for internal test purposes only
include("module-SpinAngular.jl");       using ..SpinAngular
include("module-Bsplines.jl");          using ..Bsplines
include("module-Pulse.jl");             using ..Pulse
include("module-Beam.jl")
include("module-Continuum.jl")
include("module-RadialIntegrals.jl");   using ..RadialIntegrals
include("module-HydrogenicIon.jl")
include("module-InteractionStrength.jl")
include("module-InteractionStrengthQED.jl")
include("module-Hamiltonian.jl");       using ..Hamiltonian
include("module-SelfConsistent.jl");    using ..SelfConsistent
include("module-PeriodicTable.jl")
include("module-TableStrings.jl")
include("module-Tools.jl")
include("module-AtomicState.jl");       using ..AtomicState
include("module-LSjj.jl");              using ..LSjj
include("module-BiOrthogonal.jl");      using ..BiOrthogonal

include("module-PhotoEmission.jl")

if  incProperties
# Functions/methods for atomic amplitudes
include("module-MultipoleMoment.jl")
include("module-ParityNonConservation.jl")
# Functions/methods for atomic properties
include("module-Einstein.jl")    
include("module-Hfs.jl")
include("module-IsotopeShift.jl")
include("module-LandeZeeman.jl")
include("module-FormFactor.jl")
include("module-ReducedDensityMatrix.jl")
include("module-AlphaVariation.jl")
include("module-RadiativeOpacity.jl")
include("module-MultipolePolarizibility.jl")
include("module-StarkShift.jl")
include("module-StarkZeeman.jl")
include("module-CrystalField.jl")
include("module-CrystalFieldEmission.jl")
end

if  incBasicProcesses
# Functions/methods for atomic processes
include("module-PhotoExcitation.jl")
include("module-PhotoIonization.jl")
include("module-PhotoRecombination.jl")
include("module-AutoIonization.jl")
include("module-ElectronCapture.jl")
include("module-DielectronicRecombination.jl")
include("module-PhotoExcitationFluores.jl")
include("module-PhotoExcitationAutoion.jl")
include("module-RayleighCompton.jl")
include("module-ParticleScattering.jl")
include("module-BeamPhotoExcitation.jl") 
include("module-HyperfineInduced.jl") 
include("module-ResonantInelastic.jl") 
include("module-DecayYield.jl")
include("module-ImpactExcitation.jl")
include("module-CoulombExcitation.jl")
end

if incAdvancedProcesses
# Functions/methods for more advanced atomic processes
include("module-MultiPhotonTransition.jl")
include("module-CoulombIonization.jl")
include("module-PhotoDoubleIonization.jl")
include("module-PhotoIonizationFluores.jl")
include("module-PhotoIonizationAutoion.jl")
include("module-ImpactExcitationAutoion.jl")
include("module-RadiativeAuger.jl")
include("module-MultiPhotonIonization.jl")
include("module-MultiPhotonDoubleIon.jl")
include("module-InternalConversion.jl") 
include("module-InternalRecombination.jl") 
include("module-TwoElectronOnePhoton.jl") 
include("module-DoubleAutoIonization.jl")
#= Further processes, not yet included into the code
include("module-REDA.jl")
include("module-READI.jl")
include("module-PairProduction.jl")
include("module-PairAnnihilation1Photon.jl")
include("module-PairAnnihilation2Photon.jl")  =#
end

# Functions/methods for atomic responses and time evolutions
# include("module-Statistical.jl")

if incStrongField
# Functions/methods for the computation of atomic responses
## include("module-HighHarmonic.jl")
include("module-StrongField.jl") 
end

if incAtomicCompass
# Functions/methods for the computation of atomic-compass simulations
include("module-AtomicCompass.jl") 
end

# Functions/methods for semi-empirical estimations
include("module-ImpactIonization.jl")
include("module-Semiempirical.jl")
include("module-Empirical.jl");         using ..Empirical

# Functions/methods for atomic computations
include("module-Atomic.jl");            using ..Atomic

if  incPlasma
# Functions/methods for plasma computations
include("module-Plasma.jl");            using ..Plasma
end

if  incLiouville
# Functions/methods for plasma computations
include("module-Liouville.jl");         using ..Liouville
end


if  incCascades
# Functions/methods for cascade computations
include("module-Cascade.jl");           using ..Cascade
end

# ForPedestrians depends on Cascade and must be included after it
include("module-ForPedestrians.jl");    using ..ForPedestrians

# Functions/methods for symbolic computations
if  incRacahAlgebra
include("module-RacahAlgebra.jl");      using ..RacahAlgebra
include("module-SphericalTensor.jl");   using ..SphericalTensor
end

if  incDeepLearning
include("module-AtomicFeatures.jl");    using ..AtomicFeatures
include("module-DeepLearning.jl");      using ..DeepLearning
end

# Basic functions/methods to manipulate these data
include("module-BasicsAZ.jl")
include("module-ManyElectronAZ.jl")

# Specialized macros

# All test functions/methods stay with the JAC root module
include("module-TestFrames.jl");        using ..TestFrames
    
function __init__()
    # The following variables need to be initialized at runtime to enable precompilation
    global JAC_SUMMARY_IOSTREAM = stdout
    global JAC_TEST_IOSTREAM    = stdout
end

println("\nWelcome to JenaAtomicCalculator (JAC):  A community approach to the computation of atomic structures, " *
        "cascades and time evolutions [(C) Copyright by Stephan Fritzsche, Jena (2018-2026)].")
        

end


