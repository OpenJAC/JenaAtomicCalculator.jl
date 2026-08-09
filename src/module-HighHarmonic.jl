
"""
`module  JAC.HighHarmonic
... a submodel of JAC that contains all methods to set-up and process high-harmonic computations. 
    It defines in particular a number of structs to deal with the observables, approaches, laser pulses and
    time-dependent fields as they frequently occur in the computation of high-harmonic spectra.
"""
module HighHarmonic


using Printf, ..Basics, ..Defaults, ..Radial, ..ManyElectron, ..Nuclear, ..TableStrings


"""
`abstract type HighHarmonic.AbstractTargetModel` 
    ... defines an abstract and a number of other types to represent an atomic target (cloud) for HHG:

    + struct LocalizedTargetModel     ... to represent a point-like localized target       
    + struct GaussianTargetModel      ... to represent a localized, Gaussian-type target (not yet)      
    + struct PlanarTargetModel        ... to represent a planar target of finite extent (not yet)      
    + struct UserTargetModel          ... to represent a user-defined 3-dimensional target (not yet)
"""
abstract type  AbstractTargetModel  end


"""
`struct  HighHarmonic.LocalizedTargetModel  <:  AbstractTargetModel`  
    ... defines a struct for point-like localized target cloud.

    + center     ::CartesianVector{Float64}   ... Center of the (point-like) target.
"""
struct  LocalizedTargetModel  <:  AbstractTargetModel
    center       ::CartesianVector{Float64}
end


"""
`struct  HighHarmonic.GaussianTargetModel  <:  AbstractTargetModel`  
    ... defines a struct for localized, Gaussian target cloud.

    + center     ::CartesianVector{Float64}   ... Center of the (point-like) target.
    + NoPoints   ::Int64                      ... Number of points to be used for the target.
    + sigma      ::Float64                    ... Gaussian half-widths of the cloud.
"""
struct  GaussianTargetModel  <:  AbstractTargetModel
    center       ::CartesianVector{Float64}
    NoPoints     ::Int64  
    sigma        ::Float64 
end


"""
`struct  HighHarmonic.TargetCloud`  
    ... defines a struct to represent a -- point-like or extended -- target cloud for the computation of high-harmonic spectra

    + model      ::AbstractTargetModel                    ... Model of the target cloud, from which the target is generated.
    + points     ::Array{WeightedCartesian{Float64},1}    ... Cartesian points and weights that define the atomic target cloud.
"""
struct  TargetCloud
    model        ::AbstractTargetModel
    points       ::Array{WeightedCartesian{Float64},1}
end



"""
`abstract type HighHarmonic.AbstractDeterctorModel` 
    ... defines an abstract and a number of other types to represent an dectector screen for high-harmonic spectra:

    + struct LocalizedDetectorModel      ... to represent a point-like localized detector       
    + struct SphericalDetectorModel      ... to represent a partial-sphere detector screen (not yet)      
    + struct PlanarDetectorModel         ... to represent a planar detector screen of finite extent (not yet)      
"""
abstract type  AbstractDetectorModel  end


"""
`struct  HighHarmonic.LocalizedDetectorModel  <:  AbstractDetectorModel`  
    ... defines a struct for point-like localized detector.

    + center     ::CartesianVector{Float64}   ... Center of the (point-like) detector.
"""
struct  LocalizedDetectorModel  <:  AbstractDetectorModel
    center       ::CartesianVector{Float64}
end


"""
`struct  HighHarmonic.SphericalDetectorModel  <:  AbstractDetectorModel`  
    ... defines a struct for a partly-spherical detector screen.

    + center     ::CartesianVector{Float64}   ... Center of the detector.
    + NoPoints   ::Int64                      ... Number of points to be used for the detector screen.
    + theta      ::Float64                    ... polar (arc) angle of the screen.
"""
struct  SphericalDetectorModel  <:  AbstractDetectorModel
    center       ::CartesianVector{Float64}
    NoPoints     ::Int64  
    theta        ::Float64 
end


"""
`struct  HighHarmonic.DetectorScreen`  
    ... defines a struct to represent a -- point-like or extended -- detector screen for the computation of high-harmonic spectra.

    + model      ::AbstractDetectorModel                    ... Model of the detector, from which the detection screen is generated.
    + points     ::Array{WeightedCartesian{Float64},1}      ... Cartesian points and weights that define the detection screen.
"""
struct  DetectorScreen
    model        ::AbstractDetectorModel
    points       ::Array{WeightedCartesian{Float64},1}
end



"""
`abstract type HighHarmonic.AbstractHhgObservable` 
    ... defines an abstract and a number of singleton types for the observables that can be computed for HHG processes.

    + struct HhgSpectrum     
    ... to compute a standard (single-color) high-harmonic spectrum (for a point-like target and point-like detector).       
        
    + struct HhgPolarizedSpectrum
    ... to compute a bi-color high-harmonic spectrum (for a point-like target and point-like detector).
        
    + struct AttoSecondPulse     (not yet)
    + struct IntensityProfile    (not yet)
    + struct PhaseProfile        (not yet)
        
    + struct HhgSpectrum
"""
abstract type  AbstractHhgObservable                              end
struct         HhgSpectrum           <:  AbstractHhgObservable   end
struct         HhgPolarizedSpectrum  <:  AbstractHhgObservable   end



"""
`abstract type HighHarmonic.AbstractHhgApproach` 
    ... defines an abstract and a number of singleton types for the computational approach/model that is applied in the
        computation of HHG spectra.

    + struct HhgHydrogenicSaddlePoint         
    ... to apply hydrogenic orbitals for the computation of the time-dependent dipole amplitude and the saddle-point
        approximation for the dipole moments.
        
    + struct HhgHydrogenicFullVolkov
    ... to apply hydrogenic orbitals for the computation of the time-dependent dipole amplitude but the full
        Volkov phase for the dipole moments.
"""
abstract type  AbstractHhgApproach                                 end
struct         HhgHydrogenicSaddlePoint  <:  AbstractHhgApproach   end
struct         HhgHydrogenicFullVolkov   <:  AbstractHhgApproach   end


"""
`abstract type HighHarmonic.AbstractHhgPulse` 
    ... defines an abstract and a number of other types to represent an incident laser pulse/field for HHG:

    + struct HhgLinearPlaneWaveLaser  ... to represent a long, linearly-polarized plane-wave laser field,
                                        propagating into z-direction.
    + struct HhgSineSquaredPulse      ... to represent a sin^2(...) pulse  (not yet)      
"""
abstract type  AbstractHhgPulse  end


"""
`struct  HighHarmonic.HhgLinearPlaneWaveLaser  <:  AbstractHhgPulse`  
    ... defines a struct for a long, plane-wave laser field of given frequency and intensity.

    + omega          ::Float64   ... Frequency of the incident laser field in [a.u.].
    + E0             ::Float64   ... Electric-field amplitude E0 in [a.u.].
    ...
    + intensity      ::Float64   ... Intensity [W/cm^2] of the incident laser field.
    + polarization   ::...       ... Polarization of the incident laser field.
    + beamProfile    ::...
    + envelope       ::...
    + length         ::...
    + centralLength  ::...
"""
struct  HhgLinearPlaneWaveLaser  <:  AbstractHhgPulse
    omega            ::Float64
    E0               ::Float64
end


"""
`HighHarmonic.HhgLinearPlaneWaveLaser()`  
    ... constructor for an (empty) instance of a linearly-polarized plane-wave pulse.
"""
function HhgLinearPlaneWaveLaser()
    HhgLinearPlaneWaveLaser(0., 0.)
end


"""
`HighHarmonic.HhgLinearPlaneWaveLaser(omega; intensity::Float64=1.0e14)`  
    ... constructor for an linearly-polarized plane-wave pulse of given frequency omega and intensity [W/cm^2].
"""
function HhgLinearPlaneWaveLaser(omega; intensity::Float64=1.0e14)
    omega_au = Defaults.convertUnits("energy: to atomic", omega);   E0_au = intensity / 3.51e16
    HhgLinearPlaneWaveLaser(omega_au, E0_au)
end


# `Base.show(io::IO, pulse::HhgLinearPlaneWaveLaser)`  ... prepares a proper printout of the variable pulse::HhgLinearPlaneWaveLaser.
function Base.show(io::IO, pulse::HhgLinearPlaneWaveLaser) 
    println(io, "omega:              $(pulse.omega)  ")
    println(io, "E0:                 $(pulse.E0)  ")
end



"""
`struct  HighHarmonic.Computation`  
    ... defines a type for the computation of high-harmonic spectra, atto-second pulses, intensity and phase profiles,
        and several others.

    + observable         ::AbstractHhgObservable    ... Specifies one of the observables to be calculated.
    + approach           ::AbstractHhgApproach      ... Specifies the approach/frame to be applied in the computations.
    + pulse              ::AbstractHhgPulse         ... Specifies the incident laser field.   
    + timeMesh           ::Array{Float64,1}         ... time mesh to be used for representing the incident and generated fields.
    + nuclearModel       ::Nuclear.Model            ... Model, charge and parameters of the nucleus.
    + initialOrbital     ::Radial.Orbital           ... (Single-electron) orbital to specify the initial state of atom.
    + target             ::TargetCloud              ... Specifies the target for generating high harmonics.
    + detector           ::DetectorScreen           ... Specifies the detector for observing high harmonics.
"""
struct  Computation
    observable           ::AbstractHhgObservable
    approach             ::AbstractHhgApproach
    pulse                ::AbstractHhgPulse
    timeMesh             ::Array{Float64,1}
    nuclearModel         ::Nuclear.Model
    initialOrbital       ::Radial.Orbital
    target               ::TargetCloud 
    detector             ::DetectorScreen 
end 


# `Base.show(io::IO, computation::HighHarmonic.Computation)`  ... prepares a proper printout of the variable computation::HighHarmonic.Computation.
function Base.show(io::IO, computation::HighHarmonic.Computation) 
    println(io, "observable:           $(computation.observable)  ")
    println(io, "approach:             $(computation.approach)  ")
    println(io, "pulse:                $(computation.pulse)  ")
    println(io, "timeMesh:             $(computation.timeMesh)  ")
    println(io, "nuclearModel:         $(computation.nuclearModel)  ")
    println(io, "initialOrbital:       $(computation.initialOrbital)  ")
    println(io, "target:               $(computation.target)  ")
    println(io, "detector:             $(computation.detector)  ")
end


"""
`HighHarmonic.computeTimeMesh(NoCycles::Int64, NoPointsPerCycle::Int64)` 
    ... to compute a time mesh for a photon field with given number of cycles (NoCycles) as well as time 
        points per cycle; a timeMesh::Array{Float64,1} is returned.
"""
function computeTimeMesh(NoCycles::Int64, NoPointsPerCycle::Int64)
    dt = 2pi / NoPointsPerCycle;    timeMesh = Float64[]
    for  k = 1:NoCycles*NoPointsPerCycle   push!( timeMesh, k*dt)   end
    
    println("Time mesh generated in the (time) interval [0., $(timeMesh[end])] with $(length(timeMesh)) points.")
    return( timeMesh )
end


"""
`HighHarmonic.computeTimeDipoleMoment(approach::HhgHydrogenicSaddlePoint, pulse::HhgLinearPlaneWaveLaser, 
                                        timeMesh::Array{Float64,1}, orbital::Radial.Orbital)` 
    ... to compute the time-dependent dipole moment D(t) for the given pulse, approximation and the orbital wave function.
        This dipole moment is calculated for all times as defined by timeMesh. A dipoleMoment::Array{Float64,1} is returned.
"""
function computeTimeDipoleMoment(approach::HhgHydrogenicSaddlePoint, pulse::HhgLinearPlaneWaveLaser, 
                                    timeMesh::Array{Float64,1}, orbital::Radial.Orbital)
    dipoleMoment = CartesianVector{ComplexF64}[]
    for  t in timeMesh    
        push!( dipoleMoment, CartesianVector{ComplexF64}(0., 0., 0.) )    
    end
    return( dipoleMoment )
end


"""
`HighHarmonic.computeElectricField(pulse::HhgLinearPlaneWaveLaser, time::Float64)` 
    ... to compute the electric field for the given pulse at time; an eField::CartesianVector{Float64} is returned.
"""
function computeElectricField(pulse::HhgLinearPlaneWaveLaser, time::Float64)
    return( CartesianVector{Float64}(0., 0., pulse.E0 * cos(pulse.omega*time) ) )
end


"""
`HighHarmonic.computeVectorPotential(pulse::HhgLinearPlaneWaveLaser, time::Float64)` 
    ... to compute the electric field for the given pulse at time; a aField::CartesianVector{Float64} is returned.
"""
function computeVectorPotential(pulse::HhgLinearPlaneWaveLaser, time::Float64)
    aField = CartesianVector{Float64}(0., 0., 0.)
    return( aField )
end


"""
`HighHarmonic.computeVolkovPhase(approach::HhgHydrogenicSaddlePoint, time::Float64, ip::Float64, aField::CartesianVector{Float64})` 
    ... to compute the Volkov phase at time for the given approximation, ionization potential and vector potential (aField);
        a Volkov phase::Float64 is returned.
"""
function computeVolkovPhase(approach::HhgHydrogenicSaddlePoint, time::Float64, ip::Float64, aField::CartesianVector{Float64})
    phase = 0.
    return( phase )
end


"""
`HighHarmonic.computeTimeDipoleAmplitude(approach::HhgHydrogenicSaddlePoint, time::Float64, aField::Array{Float64,1},
                                            orbital::Radial.Orbital)` 
    ... to compute the dipole amplitude d(t) at time for the given approximation, vector potential (aField) and the 
        orbital wave function; a dipole amplitude::CartesianVector{ComplexF64} is returned.
"""
function computeTimeDipoleAmplitude(approach::HhgHydrogenicSaddlePoint, time::Float64, aField::Array{Float64,1},
                                    orbital::Radial.Orbital)
    dipoleAmplitude = CartesianVector{ComplexF64}(0., 0., 0.)
    return( dipoleAmplitude )
end


"""
`HighHarmonic.computeFrequencyDipoleMoment(timeMesh::Array{Float64,1}, timeDipole::Array{CartesianVector{ComplexF64},1})` 
    ... to compute the frequency-dependent dipole moment D^~ (omega); 
        an omegaMoment::Array{CartesianVector{ComplexF64},1} is returned.
"""
function computeFrequencyDipoleMoment(timeMesh::Array{Float64,1}, timeDipole::Array{CartesianVector{ComplexF64},1})
    omegaMoment = CartesianVector{ComplexF64}[]
    for  t in timeMesh
        push!( omegaMoment, CartesianVector{ComplexF64}(0., 0., 0.) )
    end
    return( omegaMoment )
end


"""
`HighHarmonic.perform(comp::HighHarmonic.Computation; output::Bool=false)` 
    ... to perform a computation of (one) selected observable in HHG, such a spectra, attosecond pulses of profiles.
        Different approaches, pulse structures of the incident light as well as different target clouds and detector
        screens can also be specified together with all physical parameters: frequency, intensity and polarization of
        the incident light; size andshape of target clouds; size and shape of detector screens; ....
"""
function perform(comp::HighHarmonic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = comp.nuclearModel
    
    if       comp.observable == HighHarmonic.HhgSpectrum()
        timeDipole  = HighHarmonic.computeTimeDipoleMoment(comp.approach, comp.pulse, comp.timeMesh, comp.initialOrbital)
        omegaDipole = HighHarmonic.computeFrequencyDipoleMoment(comp.timeMesh, timeDipole)
    elseif   comp.observable == HighHarmonic.HhgPolarizedSpectrum()
    else     error("Undefined observable for HHG computations.")
    end
    
    return( 0. )
end

end # module


