

"""
`struct  Cascade.Data{T}`  ... defines a type for communicating different line lists to the cascade simulations

    + lines    ::Array{T,1}       ... List of lines to type T.
"""  
struct  Data{T}
    lines      ::Array{T,1}
end 


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

#== 
# These data structures have become obsolet ... and should be replaced by Cascade.Data{T} or Cascade.Data[]
"""
abstract type  Cascade.AbstractData`  
    ... defines an abstract type to distinguish different output data of a cascade computation; see also:

    + struct DecayData       ... to comprise the amplitudes/rates from a stepswise decay cascade.
    + struct PhotoIonData    ... to comprise the amplitudes/rates from a photo-ionization part of a cascade.
    + struct ExcitationData  ... to comprise the amplitudes/rates from a photo-excitation part of a cascade.
    + struct CaptureData     ... to comprise the amplitudes/rates from a electron-capture part of a cascade.
"""
abstract type  AbstractData      end


"""
`struct  Cascade.DecayData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of radiative and Auger lines.

    + linesR         ::Array{PhotoEmission.Line,1}          ... List of radiative lines.
    + linesA         ::Array{AutoIonization.Line,1}         ... List of Auger lines.
"""  
struct  DecayData  <:  Cascade.AbstractData
    linesR           ::Array{PhotoEmission.Line,1}
    linesA           ::Array{AutoIonization.Line,1}
end 


"""
`Cascade.DecayData()`  ... (simple) constructor for cascade  decay data.
"""
function DecayData()
    DecayData(Array{PhotoEmission.Line,1}[], Array{AutoIonization.Line,1}[] )
end


# `Base.show(io::IO, data::Cascade.DecayData)`  ... prepares a proper printout of the variable data::Cascade.DecayData.
function Base.show(io::IO, data::Cascade.DecayData) 
    println(io, "linesR:                  $(data.linesR)  ")
    println(io, "linesA:                  $(data.linesA)  ")
end


"""
`struct  Cascade.PhotoIonData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of photoionization lines.

    + photonEnergy   ::Float64                              ... Photon energy for this (part of the) cascade.
    + linesP         ::Array{PhotoIonization.Line,1}        ... List of photoionization lines.
"""  
struct  PhotoIonData  <:  Cascade.AbstractData
    photonEnergy     ::Float64
    linesP           ::Array{PhotoIonization.Line,1}
end 


"""
`Cascade.PhotoIonData()`  ... (simple) constructor for cascade photo-ionization data.
"""
function PhotoIonData()
    PhotoIonData(0., Array{PhotoIonization.Line,1}[] )
end


# `Base.show(io::IO, data::Cascade.PhotoIonData)`  ... prepares a proper printout of the variable data::Cascade.PhotoIonData.
function Base.show(io::IO, data::Cascade.PhotoIonData) 
    println(io, "photonEnergy:            $(data.photonEnergy)  ")
    println(io, "linesP:                  $(data.linesP)  ")
end


"""
`struct  Cascade.PhotoIonizationData  <:  Cascade.AbstractData`  
    ... defines a type for an atomic cascade with photoionization lines. This has been implemented in addition to the 
        Cascade.PhotoIonData() since the photon energy is provided by the individual lines.

    + linesP         ::Array{PhotoIonization.Line,1}        ... List of photoionization lines.
"""  
struct  PhotoIonizationData  <:  Cascade.AbstractData
    linesP           ::Array{PhotoIonization.Line,1}
end 


"""
`Cascade.PhotoIonizationData()`  ... (simple) constructor for cascade photo-ionization data.
"""
function PhotoIonizationData()
    PhotoIonizationData(Array{PhotoIonization.Line,1}[] )
end


# `Base.show(io::IO, data::Cascade.PhotoIonizationData)`  ... prepares a proper printout of the variable data::Cascade.PhotoIonizationData.
function Base.show(io::IO, data::Cascade.PhotoIonizationData) 
    println(io, "linesP:                  $(data.linesP)  ")
end


"""
`struct  Cascade.ExcitationData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of excitation lines.

    + linesE         ::Array{PhotoExcitation.Line,1}        ... List of photoexcitation lines.
"""  
struct  ExcitationData  <:  Cascade.AbstractData
    linesE           ::Array{PhotoExcitation.Line,1}
end 


"""
`Cascade.ExcitationData()`  ... (simple) constructor for cascade excitation data.
"""
function ExcitationData()
    ExcitationData(Array{PhotoExcitation.Line,1}[] )
end


# `Base.show(io::IO, data::Cascade.ExcitationData)`  ... prepares a proper printout of the variable data::Cascade.ExcitationData.
function Base.show(io::IO, data::Cascade.ExcitationData) 
    println(io, "linesE:                  $(data.linesE)  ")
end


"""
`struct  Cascade.RecombinationData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of photorecombination lines.

    + linesR         ::Array{PhotoRecombination.Line,1}        ... List of photorecombination lines.
"""  
struct  RecombinationData  <:  Cascade.AbstractData
    linesR           ::Array{PhotoRecombination.Line,1}
end 


"""
`Cascade.RecombinationData()`  ... (simple) constructor for cascade recombination data.
"""
function RecombinationData()
    RecombinationData(Array{PhotoRecombination.Line,1}[] )
end


# `Base.show(io::IO, data::Cascade.RecombinationData)`  ... prepares a proper printout of the variable data::Cascade.RecombinationData.
function Base.show(io::IO, data::Cascade.RecombinationData) 
    println(io, "linesR:                  $(data.linesR)  ")
end
==#


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
`abstract type  Cascade.AbstractSimulationProperty`  
    ... defines an abstract and various singleton types for the different properties that can be obtained from the simulation of 
        cascade data.

    + struct DecayPathes                ... determine the major 'decay pathes' of the cascade.
    + struct DrRateCoefficients         ... simulate the DR (plasma) rate coefficients for given plasma temperatures. 
    + struct ElectronCoincidence        ... simulate electron-coincidence spectra (not yet).
    + struct FinalLevelDistribution     ... simulate the 'final-level distribution' as it is found after all cascade 
                                            processes are completed.
    + struct IonDistribution            ... simulate the 'ion distribution' as it is found after all cascade processes are completed.
    + struct MeanLineWidths             ... simulate the mean line widths of a line near to a given energy (not yet). 
    + struct MeanRelaxationTime         ... simulate the mean relaxation times in which 70%, 80%, of the occupied levels
                                            decay to the ground configuration.        
    + struct ElectronIntensities        ... simulate the electron-line intensities as function of electron energy.
    + struct PhotoAbsorptionCS          ... simulate the (total) photoabsorption cross sections for a given set of photo-excitation 
                                            and ionization processes.
    + struct PhotoResonances            ... simulate the position and strength of photo-resonances for a given set of photo-excitation 
                                            amplitudes & cross sections (not yet).
    + struct PhotonIntensities          ... simulate the photon-line intensities as function of electron energy. 
    + struct PiRateCoefficients         ... simulate the PI (plasma) rate coefficients for given plasma temperatures (not yet). 
    + struct TimeBinnedPhotonIntensity  ... simulate the photon-line intensities as function of electron energy 
                                            in a given time interval (not yet).
"""
abstract type  AbstractSimulationProperty                              end


"""
`struct  Cascade.DrRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the DR plasma rate coefficients as function of the (free) electron energy and
        plasma temperature.

    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + electronEnergyShift ::Float64     
        ... (total) energy shifts that apply to all resonances when alpha^(DR) is computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + nDetailed           ::Int64       
        ... principal quantum of the `last' shell for which Auger and radiatiative amplitudes have been calculated
            by the Cascade.Computation; all contributions of this shell are scaled for nDetailed < n <= nMax
            also for higher shells by a simple scaling rule.
    + nMax                ::Int64 
        ... Maximum n (principal quantum number), for which contributions are scaled; NO scaling is taken
            into account for nMax <= nDetailed.
    + resonanceSelection  ::DielectronicRecombination.ResonanceSelection
"""  
struct  DrRateCoefficients   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64 
    electronEnergyShift   ::Float64 
    temperatures          ::Array{Float64,1}
    nDetailed             ::Int64 
    nMax                  ::Int64 
    resonanceSelection    ::DielectronicRecombination.ResonanceSelection
end 


"""
`Cascade.DrRateCoefficients()`  ... (simple) constructor for cascade DrRateCoefficients.
"""
function DrRateCoefficients()
    DrRateCoefficients(1, 0.,  Float64[], 0, 0, DielectronicRecombination.ResonanceSelection() )
end


# `Base.show(io::IO, dist::Cascade.DrRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.DrRateCoefficients.
function Base.show(io::IO, dist::Cascade.DrRateCoefficients) 
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "electronEnergyShift:      $(dist.electronEnergyShift)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "nDetailed:                $(dist.nDetailed)  ")
    println(io, "nMax:                     $(dist.nMax)  ")
    println(io, "resonanceSelection:       $(dist.resonanceSelection)  ")
end


"""
`struct  Cascade.RrRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the RR plasma rate coefficients as function of the (free) electron energy and
        plasma temperature.

    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + multipoles          ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the excitation processes.
    + finalLevelSelection ::LevelSelection    
        ... Specifies the selected final levels of some given final-state configurations; these final level numbers/
            symmetries always refer to single configurations.
    + finalConfigurations ::Array{Configuration,1}
    

"""  
struct  RrRateCoefficients   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64 
    temperatures          ::Array{Float64,1}
    multipoles            ::Array{EmMultipole} 
    finalLevelSelection   ::LevelSelection
    finalConfigurations   ::Array{Configuration,1}
end 


"""
`Cascade.RrRateCoefficients()`  ... (simple) constructor for cascade RrRateCoefficients.
"""
function RrRateCoefficients()
    RrRateCoefficients(1, Float64[], EmMultipole[], LevelSelection(), Configuration[])
end


# `Base.show(io::IO, dist::Cascade.RrRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.RrRateCoefficients.
function Base.show(io::IO, dist::Cascade.RrRateCoefficients) 
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "multipoles:               $(dist.multipoles)  ")
    println(io, "finalLevelSelection:      $(dist.finalLevelSelection)  ")
    println(io, "finalConfigurations:      $(dist.finalConfigurations)  ")
end


"""
`struct  Cascade.PiRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the PI plasma rate coefficients as function of the (free) photon energy distribution and
        the plasma temperature. These coefficients included a convolution of the direct photoionization cross sections over the 
        Maxwell distribution of electron and (if requested) add the pre-evaluated resonant parts afterwards.
        
        For the implementation: (1) Perform GL integration over free-electron contributions analog to the photoabsorption;
        (2) for the resonant part, simply add convoluted summation terms for all final states (to be worked out in detail);
        (3) use a electronEnergies grid ... instead of the photonEnergies grid.
        (4) in the direct part, sum over the final states before convolution

    + includeResonantPart ::Bool        ... True, if the resonant contributions are to be included.
    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + multipoles          ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the photoionization processes.
    + finalConfigurations ::Array{Configuration,1}
    

"""  
struct  PiRateCoefficients   <:  Cascade.AbstractSimulationProperty
    includeResonantPart   ::Bool
    initialLevelNo        ::Int64 
    temperatures          ::Array{Float64,1}
    multipoles            ::Array{EmMultipole} 
    finalConfigurations   ::Array{Configuration,1}
end 


"""
`Cascade.PiRateCoefficients()`  ... (simple) constructor for cascade PiRateCoefficients.
"""
function PiRateCoefficients()
    PiRateCoefficients(false, 1, Float64[], EmMultipole[], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PiRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.PiRateCoefficients.
function Base.show(io::IO, dist::Cascade.PiRateCoefficients) 
    println(io, "includeResonantPart:      $(dist.includeResonantPart)  ")
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "multipoles:               $(dist.multipoles)  ")
    println(io, "finalConfigurations:      $(dist.finalConfigurations)  ")
end


"""
abstract type  Cascade.AbstractOpacityDependence`  
    ... defines an abstract type to distinguish different dependencies for the opacity; see also:

    + struct FrequencyOpacityDependence     ... to deal with omega-dependence opacities [omega].
    + struct WavelengthOpacityDependence    ... to deal with wavelength-dependence opacities [lambda].
    + struct TemperatureOpacityDependence   ... to deal with temperature-normalized dependent opacities [u].
"""
abstract type  AbstractOpacityDependence      end


"""
`struct  Cascade.FrequencyOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with omega-dependence opacities [omega] and a binning that need to be given in Hartree.
"""  
struct  FrequencyOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.WavelengthOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with wavelength-dependence opacities [lambda] and a binning that need always to be given in nm.
"""  
struct  WavelengthOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.TemperatureOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with temperature-normalized dependent opacities [u] and a binning that need to be given Hartree.
"""  
struct  TemperatureOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.ExpansionOpacities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the expansion opacity as function of the wavelength as well as (parametrically) the density
        and expansion time.

    + levelPopulation        ::Basics.AbstractLevelPopulation  
        ... to specify the kind of level population that is considered for the given opacity calculations.
    + opacityDependence      ::Cascade.AbstractOpacityDependence    
        ... to specify the dependence of the opacities [omega, lambda, (temperature-normalized) u]
    + ionDensity             ::Float64       ... ion density [in g/cm^3]
    + temperature            ::Float64       ... temperature [in K]
    + expansionTime          ::Float64       ... (expansion/observation) time [in sec]
    + transitionEnergyShift  ::Float64     
        ... (total) energy shifts that apply to all transition energies; the amplitudes are re-scaled accordingly.
    + dependencyValues       ::Array{Float64,1}
        ... values [in a.u.] for which the expansion opacity is to be calculated.
"""  
struct  ExpansionOpacities      <:  Cascade.AbstractSimulationProperty
    levelPopulation          ::Basics.AbstractLevelPopulation
    opacityDependence        ::Cascade.AbstractOpacityDependence
    ionDensity               ::Float64
    temperature              ::Float64
    expansionTime            ::Float64
    transitionEnergyShift    ::Float64
    dependencyValues         ::Array{Float64,1}
end 


"""
`Cascade.ExpansionOpacities()`  ... (simple) constructor for expansion opacity simulations.
"""
function ExpansionOpacities()
    ExpansionOpacities(Basics.BoltzmannLevelPopulation(), WavelengthOpacityDependence(0.01), 1, 1000., 1.,  0., Float64[])
end


# `Base.show(io::IO, opacities::Cascade.ExpansionOpacities)`  ... prepares a proper printout of the opacities::Cascade.ExpansionOpacities.
function Base.show(io::IO, opacities::Cascade.ExpansionOpacities) 
    println(io, "levelPopulation:            $(opacities.levelPopulation)  ")
    println(io, "opacityDependence:          $(opacities.opacityDependence)  ")
    println(io, "ionDensity:                 $(opacities.ionDensity)  ")
    println(io, "temperature:                $(opacities.temperature)  ")
    println(io, "expansionTime:              $(opacities.expansionTime)  ")
    println(io, "transitionEnergyShift:      $(opacities.transitionEnergyShift)  ")
    println(io, "dependencyValues:           $(opacities.dependencyValues)  ")
end


"""
`struct  Cascade.RosselandOpacities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the Rosseland opacity as function of the temperature or density as well as (parametrically) 
        the expansion time. Usually, different values are given for either the temperature or density.

    + levelPopulation        ::Basics.AbstractLevelPopulation  
        ... to specify the kind of level population that is considered for the given opacity calculations.
    + opacityDependence      ::Cascade.TemperatureOpacityDependence    
        ... to specify the dependence of the opacities in terms of (the temperature-normalized) u]
    + ionDensities           ::Array{Float64,1}      ... list of ion densities [in g/cm^3]
    + temperatures           ::Array{Float64,1}      ... list of temperatures [in K]
    + expansionTime          ::Float64               ... (expansion/observation) time [in sec]
    + transitionEnergyShift  ::Float64     
        ... (total) energy shifts that apply to all transition energies; the amplitudes are re-scaled accordingly.
"""  
struct  RosselandOpacities      <:  Cascade.AbstractSimulationProperty
    levelPopulation          ::Basics.AbstractLevelPopulation
    opacityDependence        ::Cascade.TemperatureOpacityDependence
    ionDensities             ::Array{Float64,1}
    temperatures             ::Array{Float64,1}
    expansionTime            ::Float64
    transitionEnergyShift    ::Float64
end 


"""
`Cascade.RosselandOpacities()`  ... (simple) constructor for RosselandOpacities opacity simulations.
"""
function RosselandOpacitiesOpacities()
    RosselandOpacitiesOpacities(Basics.BoltzmannLevelPopulation(), TemperatureOpacityDependence(0.01), [1.0], [1.0], 1.,  0.)
end


# `Base.show(io::IO, opacities::Cascade.RosselandOpacities)`  
#   ... prepares a proper printout of opacities::Cascade.RosselandOpacities.
function Base.show(io::IO, opacities::Cascade.RosselandOpacities) 
    println(io, "levelPopulation:            $(opacities.levelPopulation)  ")
    println(io, "opacityDependence:          $(opacities.opacityDependence)  ")
    println(io, "ionDensities:               $(opacities.ionDensities)  ")
    println(io, "temperatures:               $(opacities.temperatures)  ")
    println(io, "expansionTime:              $(opacities.expansionTime)  ")
    println(io, "transitionEnergyShift:      $(opacities.transitionEnergyShift)  ")
end


"""
`struct  Cascade.ElectronIntensities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the electron-line intensities as function of electron energy.

    + minElectronEnergy   ::Float64     ... Minimum electron energy for the simulation of electron spectra.
    + maxElectronEnergy   ::Float64     ... Maximum electron energy for the simulation of electron spectra.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
"""  
struct  ElectronIntensities   <:  Cascade.AbstractSimulationProperty
    minElectronEnergy     ::Float64
    maxElectronEnergy     ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
end 


"""
`Cascade.ElectronIntensities()`  ... (simple) constructor for cascade ElectronIntensities.
"""
function ElectronIntensities()
    ElectronIntensities(0., 1.0e6,  [(1, 1.0)])
end


# `Base.show(io::IO, data::Cascade.ElectronIntensities)`  ... prepares a proper printout of the variable data::Cascade.ElectronIntensities.
function Base.show(io::IO, dist::Cascade.ElectronIntensities) 
    println(io, "minElectronEnergy:        $(dist.minElectronEnergy)  ")
    println(io, "maxElectronEnergy:        $(dist.maxElectronEnergy)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
end


"""
`struct  Cascade.FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the 'final-level distribution' as it is found after all cascade processes are completed.

    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
    + finalConfigs        ::Array{Configuration,1}   
        ... List of final configurations whose level population are to be "listed" finally; these configuration
            only determine the printout but not the propagation of the probabilities. If this list is empty, all 
            the levels are shown.
"""  
struct  FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
    finalConfigs          ::Array{Configuration,1}
end 


"""
`Cascade.FinalLevelDistribution()`  ... (simple) constructor for cascade FinalLevelDistribution.
"""
function FinalLevelDistribution()
    FinalLevelDistribution([(1, 1.0)], Configuration[], Configuration[])
end


# `Base.show(io::IO, data::Cascade.FinalLevelDistribution)`  
#       ... prepares a proper printout of the variable data::Cascade.FinalLevelDistribution.
function Base.show(io::IO, dist::Cascade.FinalLevelDistribution) 
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    println(io, "finalConfigs:             $(dist.finalConfigs)  ")
end


"""
`struct  Cascade.IonDistribution   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the 'ion distribution' as it is found after all cascade processes are completed.

    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  IonDistribution   <:  Cascade.AbstractSimulationProperty
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.IonDistribution()`  ... (simple) constructor for cascade IonDistribution data.
"""
function IonDistribution()
    IonDistribution([(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, data::Cascade.IonDistribution)`  ... prepares a proper printout of the variable data::Cascade.IonDistribution.
function Base.show(io::IO, dist::Cascade.IonDistribution) 
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


"""
`struct  Cascade.MeanRelaxationTime   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the mean relaxation times in which 70%, 80%, of the occupied levels decay to the 
        ground configuration.

    + timeStep            ::Float64     ... Time-step for following the decay of levels [a.u.]
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
    + groundConfigs       ::Array{Configuration,1}   
        ... List of ground configurations into which the decay is considered.
"""  
struct  MeanRelaxationTime   <:  Cascade.AbstractSimulationProperty
    timeStep              ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
    groundConfigs         ::Array{Configuration,1}
end 


"""
`Cascade.MeanRelaxationTime()`  ... (simple) constructor for cascade MeanRelaxationTime.
"""
function MeanRelaxationTime()
    MeanRelaxationTime(10., [(1, 1.0)], Configuration[], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.MeanRelaxationTime)`  ... prepares a proper printout of the variable data::Cascade.MeanRelaxationTime.
function Base.show(io::IO, dist::Cascade.MeanRelaxationTime) 
    println(io, "timeStep:                 $(dist.timeStep)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    println(io, "groundConfigs:            $(dist.groundConfigs)  ")
end


"""
`struct  Cascade.PhotonIntensities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the photon-line intensities as function of photon energy.

    + minPhotonEnergy     ::Float64     ... Minimum photon energy for the simulation of photon spectra.
    + maxPhotonEnergy     ::Float64     ... Maximum photon energy for the simulation of photon spectra.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  PhotonIntensities   <:  Cascade.AbstractSimulationProperty
    minPhotonEnergy       ::Float64
    maxPhotonEnergy       ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.PhotonIntensities()`  ... (simple) constructor for cascade PhotonIntensities.
"""
function PhotonIntensities()
    PhotonIntensities(0., 1.0e6,  [(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PhotonIntensities)`  ... prepares a proper printout of the variable data::Cascade.PhotonIntensities.
function Base.show(io::IO, dist::Cascade.PhotonIntensities) 
    println(io, "minPhotonEnergy:          $(dist.minPhotonEnergy)  ")
    println(io, "maxPhotonEnergy:          $(dist.maxPhotonEnergy)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


"""
`struct  Cascade.PhotoAbsorptionSpectrum   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the total photo-absorption cross sections in a given interval of photon energies
        as well as for a given set of photo-ionization and photo-excitation cross sections

    + includeIonization   ::Bool             ... True, if photo-ionization cross sections are to be considered.
    + includeExcitation   ::Bool             ... True, if photo-excitation lines are to be considered.
    + resonanceWidth      ::Float64          ... Widths of the resonances (user-defined units)
    + csScaling           ::Float64          ... Scaling factor do enhance the strengths of resonances (default=1.0)
    + photonEnergies      ::Array{Float64,1} ... Photon energies (in user-selected units) for the simulation of photon spectra
                                                    to describe the interval and resolution of the absorption cross sections.
    + shells              ::Array{Shell,1}   
        ... Shells that should be included for a partial absorption cross sections; a cross section contribution is considered,
            if the occupation of one of these shells is lowered in the leading configurations.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population;
            at least one of these tuples must be given.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  PhotoAbsorptionSpectrum   <:  Cascade.AbstractSimulationProperty
    includeIonization     ::Bool
    includeExcitation     ::Bool
    resonanceWidth        ::Float64
    csScaling             ::Float64
    photonEnergies        ::Array{Float64,1}
    shells                ::Array{Shell,1}
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.PhotoAbsorptionSpectrum()`  ... (simple) constructor for cascade PhotoAbsorptionSpectrum.
"""
function PhotoAbsorptionSpectrum()
    PhotoAbsorptionSpectrum(true, false, 0.1, 1., [1.0],  Shell[], [(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PhotoAbsorptionSpectrum)`  ... prepares a proper printout of the variable data::Cascade.PhotoAbsorptionSpectrum.
function Base.show(io::IO, dist::Cascade.Cascade.PhotoAbsorptionSpectrum) 
    println(io, "includeIonization:        $(dist.includeIonization)  ")
    println(io, "includeExcitation:        $(dist.includeExcitation)  ")
    println(io, "resonanceWidth:           $(dist.resonanceWidth)  ")
    println(io, "csScaling:                $(dist.csScaling)  ")
    println(io, "photonEnergies:           $(dist.photonEnergies)  ")
    println(io, "shells:                   $(dist.shells)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


struct   DecayPathes                  <:  AbstractSimulationProperty   end
struct   ElectronCoincidence          <:  AbstractSimulationProperty   end
struct   TimeBinnedPhotonIntensity    <:  AbstractSimulationProperty   end
struct   MeanLineWidth                <:  AbstractSimulationProperty   end
struct   PhotoResonances              <:  AbstractSimulationProperty   end



"""
abstract type  Cascade.AbstractSimulationMethod`  
    ... defines a abstract and a list of singleton data types for the properties that can be 'simulated' from a given
        list of lines.

    + struct ProbPropagation     ... to propagate the (occupation) probabilites of the levels until no further changes occur.
    + struct MonteCarlo          ... to simulate the cascade decay by a Monte-Carlo approach of possible pathes (not yet considered).
    + struct RateEquations       ... to solve the cascade by a set of rate equations (not yet considered).
"""
abstract type  AbstractSimulationMethod                  end
struct   ProbPropagation  <:  AbstractSimulationMethod   end
struct   MonteCarlo       <:  AbstractSimulationMethod   end
struct   RateEquations    <:  AbstractSimulationMethod   end


"""
`struct  Cascade.SimulationSettings`  ... defines settings for performing the simulation of some cascade (data).

    + printTree           ::Bool        ... Print the cascade tree in a short form
    + printLongTree       ::Bool        ... Print the cascade tree in a long form
    + initialPhotonEnergy ::Float64     ... Photon energy for which photoionization data are considered. 
"""
struct  SimulationSettings
    printTree             ::Bool
    printLongTree         ::Bool 
    initialPhotonEnergy   ::Float64
end 


"""
`Cascade.SimulationSettings()`  ... constructor for an 'empty' instance of a Cascade.Block.
"""
function SimulationSettings()
    SimulationSettings(false, false, 0.)
end


# `Base.show(io::IO, settings::SimulationSettings)`  ... prepares a proper printout of the variable settings::SimulationSettings.
function Base.show(io::IO, settings::SimulationSettings) 
    println(io, "printTree:                $(settings.printTree)  ")
    println(io, "printLongTree:            $(settings.printLongTree)  ")
    println(io, "initialPhotonEnergy:      $(settings.initialPhotonEnergy)  ")
end


"""
`struct  Cascade.Simulation`  ... defines a simulation on some given cascade (data).

    + name            ::String                              ... Name of the simulation
    + property        ::Cascade.AbstractSimulationProperty 
        ... Property that is to be considered in this simulation of the cascade (data).
    + method          ::Cascade.AbstractSimulationMethod    
        ... Method that is used in the cascade simulation; cf. Cascade.SimulationMethod.
    + settings        ::Cascade.SimulationSettings          ... Settings for performing these simulations.
    + computationData ::Array{Dict{String,Any},1}           ... Date on which the simulations are performed
"""
struct  Simulation
    name              ::String
    property          ::Cascade.AbstractSimulationProperty
    method            ::Cascade.AbstractSimulationMethod
    settings          ::Cascade.SimulationSettings 
    computationData   ::Array{Dict{String,Any},1}
end 


"""
`Cascade.Simulation()`  ... constructor for an 'default' instance of a Cascade.Simulation.
"""
function Simulation()
    Simulation("Default cascade simulation", Cascade.PhotoAbsorptionSpectrum(), Cascade.ProbPropagation(), 
                Cascade.SimulationSettings(), Array{Dict{String,Any},1}[] )
end


"""
`Cascade.Simulation(sim::Cascade.Simulation;`
    
            name=..,               property=..,             method=..,              settings=..,     computationData=.. )
            
    ... constructor for re-defining the computation::Cascade.Simulation.
"""
function Simulation(sim::Cascade.Simulation;                              
    name::Union{Nothing,String}=nothing,                                  property::Union{Nothing,Cascade.AbstractSimulationProperty}=nothing,
    method::Union{Nothing,Cascade.AbstractSimulationMethod}=nothing,      settings::Union{Nothing,Cascade.SimulationSettings}=nothing,    
    computationData::Union{Nothing,Array{Dict{String,Any},1}}=nothing )

    if  name            == nothing   namex            = sim.name              else  namex            = name                end 
    if  property        == nothing   propertyx        = sim.property          else  propertyx        = property            end 
    if  method          == nothing   methodx          = sim.method            else  methodx          = method              end 
    if  settings        == nothing   settingsx        = sim.settings          else  settingsx        = settings            end 
    if  computationData == nothing   computationDatax = sim.computationData   else  computationDatax = computationData     end 
    
    Simulation(namex, propertyx, methodx, settingsx, computationDatax)
end


# `Base.show(io::IO, simulation::Cascade.Simulation)`  ... prepares a proper printout of the variable simulation::Cascade.Simulation.
function Base.show(io::IO, simulation::Cascade.Simulation) 
    println(io, "name:              $(simulation.name)  ")
    println(io, "property:          $(simulation.property)  ")
    println(io, "method:            $(simulation.method)  ")
    println(io, "computationData:   ... based on $(length(simulation.computationData)) cascade data sets ")
    println(io, "> settings:      \n$(simulation.settings)  ")
end


"""
`struct  Cascade.LineIndex{T}`  ... defines a line index with regard to the various lineLists of data::Cascade.LineIndex.

    + lines        ::Array{T,1}              ... refers to the line list for which this index is defined.
    + process      ::Basics.AbstractProcess  ... refers to the particular lineList of cascade (data).
    + index        ::Int64                   ... index of the corresponding line.
"""
struct  LineIndex{T}
    lines          ::Array{T,1} 
    process        ::Basics.AbstractProcess
    index          ::Int64 
end 


# `Base.show(io::IO, index::Cascade.LineIndex)`  ... prepares a proper printout of the variable index::Cascade.LineIndex.
function Base.show(io::IO, index::Cascade.LineIndex) 
    println(io, "lines (typeof):        $(typeof(lines))  ")
    println(io, "process:               $(index.process)  ")
    println(io, "index:                 $(index.index)  ")
end


"""
`mutable struct  Cascade.Level`  ... defines a level specification for dealing with cascade transitions.

    + energy       ::Float64                     ... energy of the level.
    + J            ::AngularJ64                  ... total angular momentum of the level
    + parity       ::Basics.Parity               ... total parity of the level
    + NoElectrons  ::Int64                       ... total number of electrons of the ion to which this level belongs.
    + majorConfig  ::Configuration               ... major (dominant) configuration of this level.
    + relativeOcc  ::Float64                     ... relative occupation  
    + parents      ::Array{Cascade.LineIndex,1}  ... list of parent lines that (may) populate the level.     
    + daugthers    ::Array{Cascade.LineIndex,1}  ... list of daugther lines that (may) de-populate the level.     
"""
mutable struct  Level
    energy         ::Float64 
    J              ::AngularJ64 
    parity         ::Basics.Parity 
    NoElectrons    ::Int64 
    majorConfig    ::Configuration 
    relativeOcc    ::Float64 
    parents        ::Array{Cascade.LineIndex,1} 
    daugthers      ::Array{Cascade.LineIndex,1} 
end 


# `Base.show(io::IO, level::Cascade.Level)`  ... prepares a proper printout of the variable level::Cascade.Level.
function Base.show(io::IO, level::Cascade.Level) 
    println(io, "energy:        $(level.energy)  ")
    println(io, "J:             $(level.J)  ")
    println(io, "parity:        $(level.parity)  ")
    println(io, "NoElectrons:   $(level.NoElectrons)  ")
    println(io, "majorConfig:   $(level.majorConfig)  ")
    println(io, "relativeOcc:   $(level.relativeOcc)  ")
    println(io, "parents:       $(level.parents)  ")
    println(io, "daugthers:     $(level.daugthers)  ")
end


"""
`Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)`  ... returns true if both levels are equal and false otherwise.
"""
function  Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)
    if  leva.energy == levb.energy   &&   leva.J == levb.J  && leva.parity == levb.parity  && leva.NoElectrons == levb.NoElectrons 
            return( true )
    else    return( false )
    end
end

"""
`struct  Cascade.AbsorptionCrossSection`  
    ... defines the absorption cross section for a particular (incident) photon energy in terms of its discrete and 
        (direct photoionization) contributions. Of course, this absorption cross section depends on the relative population
        of the initial levels.

    + photonEnergy ::Float64              ... incident photon energy/photon-energy dependence of the absorption spectrum.
    + excitationCS ::Basics.EmProperty    ... contribution due to discrete excitation processes.
    + ionizationCS ::Basics.EmProperty    ... contribution due to contineous ionization processes.
"""
struct  AbsorptionCrossSection
    photonEnergy   ::Float64
    excitationCS   ::Basics.EmProperty
    ionizationCS   ::Basics.EmProperty
end 


# `Base.show(io::IO, cs::Cascade.AbsorptionCrossSection)`  ... prepares a proper printout of the variable cs::Cascade.AbsorptionCrossSection.
function Base.show(io::IO, cs::Cascade.AbsorptionCrossSection) 
    println(io, "photonEnergy:      $(cs.photonEnergy)  ")
    println(io, "excitationCS:      $(cs.excitationCS)  ")
    println(io, "ionizationCS:      $(cs.ionizationCS)  ")
end
