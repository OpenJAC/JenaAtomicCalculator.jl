
"""
`module  JAC.Beam`  
	... a submodel of JAC that contains all methods for dealing with different beam types and polarizations.
"""
module Beam


using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..Pulse


"""
`abstract type Beam.AbstractBeamType` 
    ... defines an abstract type to distinguish different beams; see also:
    
    + struct Beam.PlaneWave       ... to model a plane-wave beam.
    + struct Beam.BesselBeam      ... to model a Bessel beam.
    + struct Beam.LaguerreGauss   ... to model a Laguerre-Gaussian beam.
    + struct Beam.SuperposedBeam  ... to model a superposition of beams.
"""
abstract type  AbstractBeamType       end


"""
`struct  Beam.PlaneWave  <:  Beam.AbstractBeamType`   
    ... to model a plane-wave beam in terms of its k-vector and polarization.

    + kx                ::Float64      ... x-component of the wave vector.          
    + ky                ::Float64      ... y-component                 
    + kz                ::Float64      ... z-component                
"""
struct   PlaneWave  <:  Beam.AbstractBeamType
    kx                  ::Float64          
    ky                  ::Float64                
    kz                  ::Float64                
end


"""
`Beam.PlaneWave()`  ... constructor for an 'default' Beam.PlaneWave.
"""
function PlaneWave()
    PlaneWave(0., 0., 1.)
end


# `Base.string(beam::PlaneWave)`  ... provides a String notation for the variable beam::PlaneWave.
function Base.string(beam::PlaneWave)
    sa = "Plane-wave beam with k-vec = ($(beam.kx), $(beam.ky), $(beam.kz))."
    return( sa )
end


# `Base.show(io::IO, beam::PlaneWave)`  ... prepares a proper printout of the beam::PlaneWave.
function Base.show(io::IO, beam::PlaneWave)
    sa = Base.string(beam);                print(io, sa)
end


"""
`struct  Beam.BesselBeam  <:  Beam.AbstractBeamType`   
    ... to model a standard Bessel beam in terms of its opening angle and polarization.

    + mOAM              ::Int64        ... m OAM-value (or sometimes l, winding number).         
    + openingAngle      ::Float64      ... opening angle                 
    + kz                ::Float64      ... z-component                
"""
struct   BesselBeam  <:  Beam.AbstractBeamType
    mOAM                ::Int64      
    openingAngle        ::Float64               
    kz                  ::Float64               
end


"""
`Beam.BesselBeam()`  ... constructor for an 'default' Beam.BesselBeam.
"""
function BesselBeam()
    BesselBeam(0, 0., 1.)
end


# `Base.string(beam::BesselBeam)`  ... provides a String notation for the variable beam::BesselBeam.
function Base.string(beam::BesselBeam)
    sa = "Bessel beam along z-axis with OAM m = $(beam.mOAM), opening angle = $(beam.openingAngle) and " * 
            "kz = $(beam.kz)."
    return( sa )
end


# `Base.show(io::IO, beam::BesselBeam)`  ... prepares a proper printout of the beam::BesselBeam.
function Base.show(io::IO, beam::BesselBeam)
    sa = Base.string(beam);                print(io, sa)
end


"""
`struct  Beam.LaguerreGauss  <:  Beam.AbstractBeamType`   
    ... to model a Laguerre-Gaussian beam LG^l_p in terms of its two mode indices, its beam waist and the
        longitudinal wave number.

    + lOAM              ::Int64        ... l OAM-value (winding number).         
    + pRadial           ::Int64        ... radial quantum number               
    + kz                ::Float64      ... z-component                
    + waist             ::Float64      ... beam waist w_0 [a_o], the radius at which the fundamental mode has
                                           fallen to 1/e of its amplitude at the focus.

        **THE WAIST IS NOT OPTIONAL:  WITHOUT IT THE MODE IS NOT DEFINED.**  Every transverse property of a
        Laguerre-Gaussian beam is carried by w_0 -- the ring radius w_0 sqrt(|l|/2) of the intensity maximum,
        the Rayleigh range z_R = k w_0^2/2 over which the beam stays collimated, the divergence 2/(k w_0), the
        Gouy phase (2p + |l| + 1) arctan(z/z_R), and the width ~ 2/w_0 of the TRANSVERSE MOMENTUM spectrum.
        That last one is what distinguishes an LG beam from a Bessel beam:  a Bessel beam has a single
        transverse momentum kappa and therefore a sharp cone, while an LG beam has a spread of order 1/w_0
        about its mean.  A struct without w_0 can name a mode but cannot evaluate it, which is why `waist` was
        added on 11-Oct-2026;  the field is appended, so a three-argument call now raises a MethodError naming
        the type rather than silently constructing a beam with no transverse scale.
"""
struct   LaguerreGauss  <:  Beam.AbstractBeamType
    lOAM                ::Int64      
    pRadial             ::Int64      
    kz                  ::Float64               
    waist               ::Float64
end


"""
`Beam.LaguerreGauss()`  ... constructor for an 'default' Beam.LaguerreGauss.
"""
function LaguerreGauss()
    LaguerreGauss(0, 0, 1., 1.)
end


"""
`Beam.rayleighRange(beam::Beam.LaguerreGauss)`
    ... returns the Rayleigh range `z_R = k w_0^2 / 2` [a_o] of the given beam, i.e. the distance from the focus
        over which the beam cross section doubles;  a `value::Float64` is returned.
"""
function rayleighRange(beam::Beam.LaguerreGauss)
    return( beam.kz * beam.waist^2 / 2 )
end


"""
`Beam.waistAt(beam::Beam.LaguerreGauss, z::Float64)`
    ... returns the beam waist `w(z) = w_0 sqrt(1 + z^2/z_R^2)` [a_o] at a distance z from the focus;  a
        `value::Float64` is returned.  This is the whole of the z-dependence of the transverse SCALE, and for a
        target displaced along the beam it is the dominant z-effect:  a displacement b_perp is felt through the
        ratio b_perp/w(z) rather than through b_perp itself.
"""
function waistAt(beam::Beam.LaguerreGauss, z::Float64)
    return( beam.waist * sqrt(1. + (z/Beam.rayleighRange(beam))^2) )
end


"""
`Beam.gouyPhase(beam::Beam.LaguerreGauss, z::Float64)`
    ... returns the Gouy phase `(2p + |l| + 1) arctan(z/z_R)` [rad] of the given mode at a distance z from the
        focus;  a `value::Float64` is returned.  It is the one place where the RADIAL index p changes the
        propagation, and it is what distinguishes LG modes of equal |l| as they propagate.
"""
function gouyPhase(beam::Beam.LaguerreGauss, z::Float64)
    return( (2*beam.pRadial + abs(beam.lOAM) + 1) * atan( z / Beam.rayleighRange(beam) ) )
end


# `Base.string(beam::LaguerreGauss)`  ... provides a String notation for the variable beam::LaguerreGauss.
function Base.string(beam::LaguerreGauss)
    sa = "Laguerre-Gaussian beam along z-axis with OAM l = $(beam.lOAM), radial index p = $(beam.pRadial), " *
            "kz = $(beam.kz) and waist w_0 = $(beam.waist) a_o."
    return( sa )
end


# `Base.show(io::IO, beam::LaguerreGauss)`  ... prepares a proper printout of the beam::LaguerreGauss.
function Base.show(io::IO, beam::LaguerreGauss)
    sa = Base.string(beam);                print(io, sa)
end


"""
`struct  Beam.Component  <:  Beam.AbstractBeamType`   
    ... to model a single component with well-defined quantum numbers in the superposition of beams; 
        such components help define vector beams and other superpositions. No test is made with regard to the complex
        amplitude/weight of the given beam component. Various procedures will be used to specify 
        well-known superpositions as list of such beam components.

    + cw            ::Complex{Float64}      ... complex weight of the beam component.                
    + beamType      ::Beam.AbstractBeamType ... to specifies the beam component in terms of its quantum numbers.       
"""
struct   Component  <:  Beam.AbstractBeamType
    cw              ::Complex{Float64}              
    beamType        ::Beam.AbstractBeamType       
end


"""
`Beam.Component()`  ... constructor for an 'default' Beam.Component.
"""
function Component()
    Component( 1.0, LaguerreGauss() )
end


# `Base.string(beam::Component)`  ... provides a String notation for the variable beam::Component.
function Base.string(beam::Component)
    sa = "Single component with cw = $(beam.cw)  of  a $(beam.beamType)"
    return( sa )
end


# `Base.show(io::IO, beam::Component)`  ... prepares a proper printout of the beam::Component.
function Base.show(io::IO, beam::Component)
    sa = Base.string(beam);                print(io, sa)
end


"""
`struct  Beam.PhotonBeam`   
    ... to model a full photon beams withwell-defined beamType, polarization, time envelope, 
        central frequency (omega), intensity and carrier-envelope phase. 
        Such photon beams can be utilized for strong-field phenomena, atomic-compass simulations 
        and other atom-beam interactions.
        
        !! This data type is still under development ... and improvements are welcome.

    + beamType      ::Beam.AbstractBeamType         ... to specifies the beam component in terms of its quantum numbers.
    + polarization  ::Basics.AbstractPolarization   ... to specify the polarization of the beam.
    + envelope      ::Pulse.AbstractEnvelope        ... to characterize the envelope-function
    + omega         ::Float64                       ... frequency (a.u.)
    + intensity     ::Float64                       ... intensity of the beam (which units ??)
    + cep           ::Float64                       ... cep-phase
"""
struct  PhotonBeam
    beamType      ::Beam.AbstractBeamType 
    polarization  ::Basics.AbstractPolarization
    envelope      ::Pulse.AbstractEnvelope
    omega         ::Float64 
    intensity     ::Float64 
    cep           ::Float64 
end


"""
`Beam.PhotonBeam()`  ... constructor for an 'default' Beam.PhotonBeam.
"""
function PhotonBeam()
    PhotonBeam( LaguerreGauss(), Basics.LinearPolarization() )
end


# `Base.show(io::IO, beam::Beam.PhotonBeam)`  ... prepares a proper printout of the beam::Beam.PhotonBeam.
function Base.show(io::IO, beam::Beam.PhotonBeam)
    println(io, "beamType:           $(beam.beamType)  ")
    println(io, "polarization:       $(beam.polarization)  ")
    println(io, "envelope:           $(beam.envelope)  ")
    println(io, "omega:              $(beam.omega)  ")
    println(io, "intensity:          $(beam.intensity)  ")
    println(io, "cep:                $(beam.cep)  ")
end


"""
`abstract type Beam.AbstractObservable` 
    ... defines an abstract type for an observable in a (twisted) beam; see also:
    
    + struct Beam.NoObservable           ... to represent no observable, usually used just for initialization.
    + struct Beam.DominantMultipoles     ... to represent regions where selected multipoles are dominant.
    + struct Beam.IntensityPattern       ... to represent the intensity pattern of a given excitation/ionization.
    + struct Beam.AnisotropyParameter    ... to represent the anisotropy A_kq for atoms at given impact vector.
"""
abstract type  AbstractObservable                               end
struct         NoObservable       <:  Beam.AbstractObservable   end


"""
`struct  Beam.DominantMultipoles  <:  Beam.AbstractObservable`   
    ... to represent regions where selected multipoles are dominant.

    + mesh              ::Basics.AbstractMesh   ... to represent a mesh at which this dominance has to be determined.         
"""
struct   DominantMultipoles  <:  Beam.AbstractObservable
    mesh                ::Basics.AbstractMesh         
end


# `Base.show(io::IO, obs::DominantMultipoles)`  ... prepares a proper printout of the obs::DominantMultipoles.
function Base.show(io::IO, obs::Beam.DominantMultipoles)
    sa = "Dominant multipole are to be determined at mesh: $(obs.mesh).";         print(io, sa)
end


"""
`struct  Beam.IntensityPattern  <:  Beam.AbstractObservable`   
    ... to represent the intensity pattern of a given excitation/ionization.

    + mesh              ::Basics.AbstractMesh   ... to represent a mesh at which the intensity pattern has to be determined.         
"""
struct   IntensityPattern  <:  Beam.AbstractObservable
    mesh                ::Basics.AbstractMesh         
end


# `Base.show(io::IO, obs::IntensityPattern)`  ... prepares a proper printout of the obs::IntensityPattern.
function Base.show(io::IO, obs::Beam.IntensityPattern)
    sa = "Intensity pattern is to be determined at mesh: $(obs.mesh).";         print(io, sa)
end


"""
`function  Beam.constructSuperposition(key::String)`   
    ... to construct various (vector) beams with well-defined polarization properties; allowed keys are:

    + "radially-polarized LG vector beam"      ... to define a radially-polarized LG vector beam.         
    + "azimuthally-polarized LG vector beam"   ... to define an azimuthally-polarized LG vector beam.   
    + "..."
    
    A superposition::Array{Beam.Component,1} of beam components is returned.
"""
function  constructSuperposition(key::String)
    wa = Beam.Component[];   cwsq2 = 1.0 / sqrt(2.0)
    #
    if      key == "radially-polarized LG vector beam"
        # Obviously, this is not correct
        push!(wa, Beam.Component(cwsq2, Beam.LaguerreGauss()) )
        push!(wa, Beam.Component(cwsq2, Beam.BesselBeam()) )
    elseif  key == "azimuthally-polarized LG vector beam"
    else    error("Unknown beam type")
    end
    
    return(wa) 
end


"""
`function  Beam.redefineEnergy(energy::Float64, beam::AbstractBeamType)`   
    ... to re-define the energy (i.e. the components of the k-vector) in Hartree for a given beam to a well-defined value.
        A newBeam of the same (concrete) type as beam is returned.
"""
function  redefineEnergy(energy::Float64, beam::AbstractBeamType)
    k = sqrt(2*energy)
    #
    if      typeof(beam) == Beam.PlaneWave
        if  beam.kx != 0  ||   beam.ky != 0   error("stop a")   end
        newBeam = Beam.PlaneWave(0., 0., k)
    elseif  typeof(beam) == Beam.BesselBeam
        kz = k * cos(beam.openingAngle)
        newBeam = Beam.BesselBeam( beam.mOAM, beam.openingAngle, kz)
    else    error("stop b")   
    end
    
    return(newBeam)
end


end # module
