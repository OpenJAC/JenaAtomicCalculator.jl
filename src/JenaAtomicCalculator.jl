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

export AbstractCImethod, AbstractConfigurationRestriction, AbstractEeInteraction, AbstractPotential, AbstractQedModel, AbstractStartOrbitals,
       AbstractProcessSettings, AbstractEmpiricalSettings, AbstractPlasmaModel, AbstractPropertySettings, AbstractLineShiftSettings, AbstractPlotTheme,
       AbstractNeutralNetwork, AbstractNeutralNetworkRequest, Application,
       add, AllShells, AlphaVariation, analyze, AnapoleMoment,
       AngularJ64, AngularM64, AngularJ, AngularMomentum, AddElectrons, 
       AsfSettings, Atomic, AtomicState, AtomicStructure, Auger, AugerInPlasma, AutoIonization, AverageAtom, AtomicCompass,
       AtomicModel, AtomicFeatures,
       Basics, Basis, Beam, BeamPhotoExcitation, BiOrthogonal, BreitInteraction, Bsplines, ByMultipoles, ByNumber, ByParity,
       CartesianVector, Cartesian2DFieldVector, Cartesian3DFieldVector, CiSettings, CiExpansion, ClebschGordan, CloseCoupling, 
       compute, convertUnits, Compton, Configuration, ConfigurationR, 
       Cascade, Continuum, CorePolarization, Coulex, CoulombExcitation, Coulion, CoulombBreit, CoulombGaunt, 
       CoulombInteraction, CoulombIonization, CrystalField, CrystalFieldEmission, CsfR, ClosedCore, ClosedShells, ClosedSubshells, ContractShells, checkConfigurations,
       computeBranchingFractions,  computeChargeStateDistribution,  computeCrossSections,  computeForPedestrians,  computeLevelEnergies,  computeLifetimes,  computeResonanceStrength, computeTransitionRates,  displaySpectrum,
       diagonalize, Defaults, DecayYield, DielectronicRecombination, Dierec, Djpq, DoubleAutoIonization, DoubleAuger, DeepLearning,
       DiagonalCoulomb, DefaultQuantizationAxis, displayCouplings, displayConfiguration,  displayConfigurations, Distribution,
       Eimex, ElectronCapture, ElecCapture, estimate, ElectricDipoleMoment, Einstein, EmMultipole, evaluate, ExpStokes,
       Empirical, ExciteElectrons, ExcitationLevel, ExpandShells, estimateCrossSections, extractConfiguration,  extractConfigurations,  
       extractFromConfiguration, extractFromConfigurations, EmProperty, EmPropertyC,
       E1, M1, E2, M2, E3, M3, E4, M4,
       FormFactor, FullCIeigen, ForPedestrians, ForGivenConfigs,
       ForAutoIonization, ForDielectronicCapture, ForDielectronicRecombination, ForElectronCapture, ForHollowIons, ForImpactIonization, 
       ForPhotoEmission, ForPhotoIonization, ForPhotoRecombination, ForRasExcitations, ForStepwiseDecay, FineStructure, FineStructureLS, FromBasis, 
       generate, GreenSettings, GreenChannel, GreenExpansion, getDefaults, Green, Gui, GroundConfiguration, GeneralizedConfigurations,
       GetParity, generateConfiguration, generateConfigurations,
       Hamiltonian, Hfs, HyperfineInduced, HighHarmonic, HydrogenicIon, HarmonicQuantizationAxis, HundsRules, HyperfineStructure,
       interpolate, integrate, Integral, ImpactExcAuto, ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, 
       InteractionStrength, InternalConv, InternalConversion, InternalRecombination, Isotope, IsotopeShift, IsotopicFraction, IsOccupied,
       Kronecker, 
       LandeF, LandeJ, LandeZeeman, Level, LevelSelection, LevelSymmetry, LineSelection, LSjj, LSjjSettings, LeftCircular, 
       LeadingConfiguration, LeadingConfigurationR, Liouville,
       ManyElectron, MeanFieldSettings, MeanFieldBasis, MeanFieldMultiplet, minus, Model, modify, 
       MultiPhotonDE, MultiPhotonDeExcitation, MultiPhotonDoubleIon, 
       MultiPI, MultiPDI, MultiPhotonIonization, MultipoleMoment, MultipolePolarizibility, Multiplet, MeanConfiguration, MeanOccupation, 
       Multiplicity,
       NoAmplitude, Nuclear, NoneQed, NoProcess, NoProperty, NonrelativisticBasis, NumberOfElectrons,
       OneElectronSettings, OneElectronSpectrum, Orbital, oplus, OccupationDifference, OpenShellNumber, OpenShells, OpenSubshells,
       PairA1P, PairAnnihilation1Photon, PairAnnihilation2Photon, PairProduction, Parity, ParityNonConservation, ParticleScattering, 
       PathwaySelection, PeriodicTable, perform, 
       Photo, PhotoDouble, PhotoDoubleIonization, PhotoEmission, PhotoExc, PhotoExcAuto, PhotoExcFluor, 
       PhotoExcitation, PhotoExcitationAutoion, PhotoExcitationFluores, PhotoIonAuto, PhotoIonFluor, PhotoIonization, 
       PhotoIonizationAutoion, PhotoIonizationFluores, PhotoRecombination,   
       Plasma, plus, PrintWarnings, provide, Pulse,
       QedPetersburg, QedSydney, 
       RacahAlgebra, RacahExpression, Radial, RadialIntegrals, Radiative, RadiativeAuger, RAuger, RasSettings, RasStep,
       RasLayer, RasExpansion, RayleighCompton, recast, Rec, REDA, READI, Representation, ReducedDensityMatrix, RadiativeOpacity,
       RadialPotentials, RadialOrbitalsLarge, RadialOrbitalsSmall, RadialOrbitalsBoth,
       RestrictMaximumDisplacements, RestrictNoElectronsTo, RestrictParity, RestrictToShellDoubles, RequestMinimumOccupation, RequestMaximumOccupation,
       ResonantInelastic, RemoveElectrons, RestrictExcitations, RelativisticConfigurations, run,
       SchiffMoment, Semiempirical, setDefaults, Shell, ShellSelection, SolidAngle, Spectroscopy, SphericalTensor, SpinAngular, StarkShift,
       StartFromHydrogenic, StartFromPrevious, StrongField, StrongField2, Subshell, StaticQuantizationAxis, StaticField, SelfConsistent, 
       SuperConfiguration,
       tabulate, TestFrames, tools, Triangle, TwoElectronOnePhoton, TimeHarmonicField, TotalAM,
       UseBabushkin, UseCoulomb, UseGauge,
       ValenceOccupation, ValenceShells,
       WeightedCartesian, W3j, W6j, W9j,
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
include("module-StarkShift.jl")
include("module-FormFactor.jl")
include("module-ReducedDensityMatrix.jl")
include("module-AlphaVariation.jl")
include("module-RadiativeOpacity.jl")
include("module-MultipolePolarizibility.jl")
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
include("module-MultiPhotonDeExcitation.jl")
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


