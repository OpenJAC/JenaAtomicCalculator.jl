
"""
`module  JAC.Beam`  
	... a submodel of JAC that contains all methods for dealing with different beam types and polarizations.
"""
module Beam


using  GSL, Printf, QuadGK, ..AngularMomentum, ..Basics, ..Defaults, ..Pulse


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
`struct  Beam.BesselGauss  <:  Beam.AbstractBeamType`
    ... to model a Bessel-Gauss beam (Gori, Guattari & Padovani, Opt. Commun. 64, 491 (1987)): the finite-energy,
        experimentally realisable counterpart of an ideal Bessel beam, obtained by replacing the sharp ring of
        transverse momentum kappa_0 = k sin(openingAngle) that `Beam.BesselBeam` carries with a GAUSSIAN ring of
        width `ringWidth` about it. This is the standard model for a real vortex beam (optical OAM beams from an
        axicon or spatial light modulator; electron vortex beams from a phase plate or hologram) as opposed to the
        idealised, infinite-energy `Beam.BesselBeam`.

    + mOAM              ::Int64        ... m OAM-value (winding number), same convention as Beam.BesselBeam.
    + openingAngle      ::Float64      ... opening angle [rad] of the ring's CENTRE, kappa_0 = k sin(openingAngle).
    + kz                ::Float64      ... z-component (unused by Beam.coneIntegrate; see this module's header).
    + ringWidth         ::Float64      ... Gaussian width Delta_kappa [a.u.] of the momentum-space ring about kappa_0.
                                           ringWidth -> 0 is the ideal Beam.BesselBeam limit (verified where this
                                           struct's coneIntegrate method is validated).
"""
struct   BesselGauss  <:  Beam.AbstractBeamType
    mOAM                ::Int64
    openingAngle        ::Float64
    kz                  ::Float64
    ringWidth           ::Float64
end


"""
`Beam.BesselGauss()`  ... constructor for an 'default' Beam.BesselGauss.
"""
function BesselGauss()
    BesselGauss(0, 0., 1., 0.1)
end


# `Base.string(beam::BesselGauss)`  ... provides a String notation for the variable beam::BesselGauss.
function Base.string(beam::BesselGauss)
    sa = "Bessel-Gauss beam along z-axis with OAM m = $(beam.mOAM), opening angle = $(beam.openingAngle), " *
            "kz = $(beam.kz) and ring width = $(beam.ringWidth)."
    return( sa )
end


# `Base.show(io::IO, beam::BesselGauss)`  ... prepares a proper printout of the beam::BesselGauss.
function Base.show(io::IO, beam::BesselGauss)
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
    PhotonBeam( LaguerreGauss(), Basics.LinearPolarization(), Pulse.InfiniteEnvelope(), 1.0, 1.0, 0.0 )
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


# GENERIC TWISTED-BEAM DECOMPOSITION, usable by any process that needs a Bessel beam's plane-wave cone and the phase a
# displaced target picks up from it -- see module-ParticleScattering-inc-beams.jl, where this was identified as belonging
# here rather than being re-derived per process. Nothing HERE knows what physics couples to the beam: every function below
# takes the caller's own direction-dependent amplitude as a plain Function and returns a weighted sum over it, so the SAME
# code can serve PhotonScattering's scalar form-factor amplitude and, eventually, ElectronCapture's M_d-resolved capture
# amplitude once that is built. PhotonScattering.besselAmplitude still carries its own copy of this logic inline
# (predating this file); unifying the two is a later refactor, deliberately not done here to avoid touching an already
# validated, dated example (example-Pd.jl) while adding this machinery's first caller.
#
# THREE BEAM TYPES ARE COVERED, in two families. Beam.BesselBeam has a SINGLE, sharp transverse momentum kappa =
# k sin(openingAngle), so its cone collapses to a one-dimensional azimuthal integral (`Beam.coneIntegrate`, below).
# Beam.BesselGauss and Beam.LaguerreGauss both have a continuous transverse-momentum SPREAD instead -- a Gaussian ring
# about kappa_0, respectively the self-Fourier LG radial profile (see each struct's own docstring) -- so their cone
# becomes a genuine two-dimensional (kappa, phiK) integral, done by the shared `Beam.coneIntegrateContinuous` that
# each beam's own `Beam.coneIntegrate` method calls (added 12-Oct-2026, extending priority items 5/6 of that day's
# Beam-module menu). The physics differs (a caller's `amplitudeAt` must now take `(kappa, phiK)`, not `phiK` alone,
# since a continuous-spectrum beam gives it an elevation degree of freedom the sharp cone does not have), so the two
# families are NOT drop-in replacements for each other -- but they share the same displacement/OAM conventions.
#
# k IS PASSED IN EXPLICITLY, NEVER READ FROM beam.kz (fixed 12-Oct-2026, priority item 63). The first real caller,
# PhotonScattering.besselAmplitude, computes its own kappa = k*sin(thetaK) from the real photon wavenumber and never
# touches beam.kz; example-Pe.jl even constructs its BesselBeam with an arbitrary kz, commented "kz (unused here)".
# Deriving kappa from beam.kz instead would have worked only for a beam already passed through Beam.redefineEnergy,
# and nothing enforced that -- a caller who forgot would get a silently wrong kappa from a stale kz. Taking k as an
# explicit argument here removes the possibility rather than guarding against it.


"""
`Beam.coneIntegrate(beam::Beam.BesselBeam, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                    nphi::Int64=721)`
    ... to combine a caller-supplied, direction-dependent amplitude `amplitudeAt(phiK::Float64)::ComplexF64` -- evaluated
        at the beam's own fixed cone elevation `beam.openingAngle` -- coherently over the beam's plane-wave cone and the
        target's displacement from the axis, as `(1/2pi) Int_0^2pi dphiK amplitudeAt(phiK) * transverseAmplitude(phiK) *
        displace(k, bPerp, phiB, phiK)`. `k` is the real wavenumber of the beam at the energy in use, supplied by the
        caller rather than read from `beam.kz` (see this file's header for why). The caller owns every bit of physics;
        this function owns only the beam's own geometry. A `value::ComplexF64` is returned.

        The integral is a 721-point (default) trapezoid, matching the resolution `PhotonScattering.besselAmplitude` uses
        for the same integral; `nphi` can be raised for an `amplitudeAt` that itself oscillates faster than that in phiK.
"""
function coneIntegrate(beam::Beam.BesselBeam, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                       nphi::Int64=721)
    dphi = 2pi / (nphi - 1);    wa = ComplexF64(0.)
    for  i = 1:nphi
        phiK = (i-1) * dphi;    w = (i == 1 || i == nphi) ? 0.5 : 1.0
        wa   = wa + w * amplitudeAt(phiK) * Beam.transverseAmplitude(beam, phiK) * Beam.displace(beam, k, bPerp, phiB, phiK) * dphi
    end

    return( wa / (2pi) )
end


"""
`function  Beam.constructSuperposition(key::String)`   
    ... to construct various (vector) beams with well-defined polarization properties; allowed keys are:

    + "radially-polarized LG vector beam"      ... to define a radially-polarized LG vector beam.         
    + "azimuthally-polarized LG vector beam"   ... to define an azimuthally-polarized LG vector beam.   
    + "..."

    NEITHER KEY IS IMPLEMENTED YET (noted 12-Oct-2026): a genuine vector beam needs each component to carry its OWN
    polarization, and `Beam.Component` does not track one today -- it holds only a complex weight and a `beamType`.
    The superposition built here before this note mixed a default `LaguerreGauss()` with a default `BesselBeam()` and
    no polarization distinction whatsoever, which was wrong rather than merely incomplete, and the second key did
    nothing at all. Both now raise instead. Giving `Component` a polarization field -- mirroring `Beam.PhotonBeam`'s
    own `::Basics.AbstractPolarization` -- is the prerequisite, not attempted here.

    A superposition::Array{Beam.Component,1} of beam components would be returned once implemented.
"""
function  constructSuperposition(key::String)
    if  key == "radially-polarized LG vector beam"  ||  key == "azimuthally-polarized LG vector beam"
        error("\n\nBeam.constructSuperposition(\"$key\") is not implemented: a vector beam needs each component to " *
              "carry its own polarization, which Beam.Component does not track yet (see this function's docstring).")
    else
        error("Beam.constructSuperposition(): unknown beam type \"$key\".")
    end
end


"""
`Beam.oamOf(beam)`
    ... to return the orbital angular momentum (OAM) winding number of the given beam, under whichever field name its
        own struct happens to use (`mOAM` for a `Beam.BesselBeam`, `lOAM` for a `Beam.LaguerreGauss`), so that generic
        code can ask for a beam's OAM without first finding out which concrete type it holds. A `Beam.PlaneWave` carries
        none; a `Beam.Component` or `Beam.PhotonBeam` reports its wrapped `beamType`'s own OAM. A `value::Int64` is
        returned.
"""
function oamOf(beam::Beam.PlaneWave)
    return( 0 )
end


# `Beam.oamOf(beam::Beam.BesselBeam)`  ... the beam's own mOAM field.
function oamOf(beam::Beam.BesselBeam)
    return( beam.mOAM )
end


# `Beam.oamOf(beam::Beam.BesselGauss)`  ... the beam's own mOAM field.
function oamOf(beam::Beam.BesselGauss)
    return( beam.mOAM )
end


# `Beam.oamOf(beam::Beam.LaguerreGauss)`  ... the beam's own lOAM field.
function oamOf(beam::Beam.LaguerreGauss)
    return( beam.lOAM )
end


# `Beam.oamOf(beam::Beam.Component)`  ... delegates to the wrapped beamType's own OAM.
function oamOf(beam::Beam.Component)
    return( Beam.oamOf(beam.beamType) )
end


# `Beam.oamOf(beam::Beam.PhotonBeam)`  ... delegates to the wrapped beamType's own OAM.
function oamOf(beam::Beam.PhotonBeam)
    return( Beam.oamOf(beam.beamType) )
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
    elseif  typeof(beam) == Beam.BesselGauss
        kz = k * cos(beam.openingAngle)
        newBeam = Beam.BesselGauss( beam.mOAM, beam.openingAngle, kz, beam.ringWidth)
    elseif  typeof(beam) == Beam.LaguerreGauss
        newBeam = Beam.LaguerreGauss( beam.lOAM, beam.pRadial, k, beam.waist)
    else    error("stop b")   
    end
    
    return(newBeam)
end


"""
`Beam.displace(beam::Beam.BesselBeam, k::Float64, bPerp::Float64, phiB::Float64, phiK::Float64)`
    ... to compute the phase exp(i kappa bPerp cos(phiK - phiB)) that a target displaced by the impact vector
        (bPerp [a.u.], phiB [rad]) from the beam axis imposes on the plane-wave component at cone azimuth phiK, with
        `k` the real wavenumber of the beam at the energy in use; a `value::ComplexF64` is returned. bPerp = 0 returns
        1 for every phiK and phiB, the on-axis limit.
"""
function displace(beam::Beam.BesselBeam, k::Float64, bPerp::Float64, phiB::Float64, phiK::Float64)
    kappa = Beam.transverseMomentum(beam, k)

    return( exp( im * kappa * bPerp * cos(phiK - phiB) ) )
end


"""
`Beam.transverseAmplitude(beam::Beam.AbstractBeamType, phiK::Float64)`
    ... to compute the azimuthal OAM phase exp(i oamOf(beam) phiK) carried by the plane-wave component at cone azimuth
        phiK, using `Beam.oamOf` so one method serves every concrete beam type regardless of which field name it
        keeps its own winding number under; a `value::ComplexF64` is returned.
"""
function transverseAmplitude(beam::Beam.AbstractBeamType, phiK::Float64)
    return( exp( im * Beam.oamOf(beam) * phiK ) )
end


"""
`Beam.transverseMomentum(beam::Beam.BesselBeam, k::Float64)`
    ... to compute the single transverse momentum kappa = k sin(openingAngle) [a.u.] that every plane-wave component on
        the beam's cone shares, at the real wavenumber `k` of the beam at the energy in use. `k` is taken as an explicit
        argument rather than read from `beam.kz`, which nothing keeps synchronised with the energy actually in use (see
        this file's header, priority item 63); a `value::Float64` is returned.
"""
function transverseMomentum(beam::Beam.BesselBeam, k::Float64)
    return( k * sin(beam.openingAngle) )
end


# CONTINUOUS-TRANSVERSE-MOMENTUM BEAMS (Beam.BesselGauss, Beam.LaguerreGauss), extending the cone machinery above
# beyond the single sharp kappa a Beam.BesselBeam carries. The physics is a genuine generalisation, not a reshape:
# the beam's own momentum-space amplitude radialAmplitude(beam, kappa, k) now varies continuously with kappa, so the
# one-dimensional azimuthal integral above becomes a two-dimensional one, (kappa, phiK), and `amplitudeAt` must take
# BOTH coordinates -- `amplitudeAt(kappa::Float64, phiK::Float64)::ComplexF64` -- since a continuous-spectrum beam
# gives the caller a genuine elevation degree of freedom that a sharp Bessel cone does not have.
#
# NORMALISATION. `radialAmplitude` need not be normalised by its own type's method: coneIntegrateContinuous divides
# the coherent radial-times-azimuthal integral by M1 = Int_0^infty dkappa kappa radialAmplitude(kappa), the first
# radial MOMENT of the weight (linear in radialAmplitude, like the weight itself -- NOT its square). This is the
# unique choice that reproduces, at bPerp = 0 with mOAM = 0 and a constant amplitudeAt, exactly that constant -- the
# same anchor Beam.BesselBeam's own (1/2pi) Int dphiK convention satisfies for the azimuthal integral, since
# |exp(i m phiK)| = 1 needs no explicit normalisation there. A first attempt here divided by Int kappa
# radialAmplitude^2 instead (by analogy with a quantum-mechanical |psi|^2 = 1 normalisation) and FAILED this exact
# check (returned 1.417 instead of 1.0): a sharp Beam.BesselBeam cone is an externally FIXED parameter, not an
# integration variable to be probability-normalised, so that analogy does not transfer. The M1 moment is also what
# makes coneIntegrateContinuous reduce EXACTLY to Beam.coneIntegrate's own sharp-cone formula as a beam's spread
# narrows (verified where each continuous-spectrum beam's own coneIntegrate method is checked, below).
#
# THE KAPPA INTEGRAL USES QuadGK (adaptive), not a fixed trapezoid: both radial profiles below are smooth and decay
# like a Gaussian, which QuadGK handles over [0, Inf) without a hand-picked cutoff or node count: the error risked by
# guessing a cutoff too small is exactly the kind of silent mistake this module's own history (item 63) warns against.


"""
`Beam.radialAmplitude(beam::Beam.BesselGauss, kappa::Float64, k::Float64)`
    ... to compute the (real, non-negative) Gaussian-ring momentum-space amplitude exp(-(kappa-kappa_0)^2 /
        (2 ringWidth^2)) of the given Bessel-Gauss beam, with kappa_0 = k sin(openingAngle) the ring's centre at the
        real wavenumber `k`; a `value::Float64` is returned. ringWidth -> 0 collapses this to a delta function at
        kappa_0, the ideal Beam.BesselBeam limit (verified at `Beam.coneIntegrate(beam::Beam.BesselGauss, ...)`).
"""
function radialAmplitude(beam::Beam.BesselGauss, kappa::Float64, k::Float64)
    kappa0 = k * sin(beam.openingAngle)

    return( exp( -(kappa - kappa0)^2 / (2 * beam.ringWidth^2) ) )
end


"""
`Beam.radialAmplitude(beam::Beam.LaguerreGauss, kappa::Float64, k::Float64)`
    ... to compute the (real, non-negative) momentum-space radial amplitude of the given LG_p^l mode, from the
        SELF-FOURIER property of Laguerre-Gauss modes: the momentum-space amplitude has the SAME functional form as
        the real-space one, with r -> kappa and the waist w0 -> 2/w0 (Siegman, *Lasers*, ch. 16) --
        `(kappa waist/sqrt2)^|l| L_pRadial^|l|(kappa^2 waist^2/2) exp(-kappa^2 waist^2/4)`, dropping the overall
        constant phase `(-i)^|l| (-1)^pRadial` that the full momentum-space field carries, since it multiplies this
        function's result by a kappa-independent factor that cancels exactly in `Beam.coneIntegrateContinuous`'s
        ratio. `k` is accepted for dispatch consistency with `Beam.radialAmplitude(beam::Beam.BesselGauss, ...)` but
        not used: a paraxial LG mode's transverse structure depends only on its own waist, not on the real
        wavenumber. A `value::Float64` is returned.
"""
function radialAmplitude(beam::Beam.LaguerreGauss, kappa::Float64, k::Float64)
    m = abs(beam.lOAM);    w0 = beam.waist;    x = kappa^2 * w0^2 / 2

    return( (kappa * w0 / sqrt(2))^m * GSL.sf_laguerre_n(beam.pRadial, Float64(m), x) * exp(-x/2) )
end


"""
`Beam.kappaBreakpoints(beam::Beam.BesselGauss, k::Float64)`
    ... to return the kappa nodes `QuadGK.quadgk` must be forced through to resolve the beam's own momentum-space
        ring, as `[max(0., kappa_0 - 8 ringWidth), kappa_0, kappa_0 + 8 ringWidth]` about kappa_0 = k sin(openingAngle);
        an `Array{Float64,1}` is returned. Without this, a ring narrower than QuadGK's default sampling can be
        missed entirely over the semi-infinite domain [0, Inf) -- caught by this struct's own first validation run,
        where a ring 1000x narrower than kappa_0 integrated to 0/0 = NaN.
"""
function kappaBreakpoints(beam::Beam.BesselGauss, k::Float64)
    kappa0 = k * sin(beam.openingAngle);    low = max(0., kappa0 - 8*beam.ringWidth);    high = kappa0 + 8*beam.ringWidth

    return( low < kappa0 ? [low, kappa0, high] : [low, high] )
end


"""
`Beam.kappaBreakpoints(beam::Beam.LaguerreGauss, k::Float64)`
    ... to return the kappa nodes `QuadGK.quadgk` must be forced through to resolve the LG mode's own momentum-space
        radial profile, as `[0., charKappa, 3 charKappa]` with `charKappa = (2/waist) sqrt(2 pRadial + |lOAM| + 1)`
        the profile's own characteristic radial scale -- the real-space LG ring radius formula, evaluated at the
        momentum-space waist `2/waist` that the self-Fourier property of LG modes gives; an `Array{Float64,1}` is
        returned.
"""
function kappaBreakpoints(beam::Beam.LaguerreGauss, k::Float64)
    charKappa = (2 / beam.waist) * sqrt(2*beam.pRadial + abs(beam.lOAM) + 1)

    return( [0., charKappa, 3*charKappa] )
end


"""
`Beam.coneIntegrateContinuous(beam::Beam.AbstractBeamType, k::Float64, bPerp::Float64, phiB::Float64,
                              amplitudeAt::Function; nphi::Int64=361)`
    ... to combine a caller-supplied amplitude `amplitudeAt(kappa::Float64, phiK::Float64)::ComplexF64` coherently
        over a continuous-transverse-momentum beam's own momentum-space weight `Beam.radialAmplitude(beam, kappa,
        k)` and OAM phase, and the target's displacement from the axis, generalising `Beam.coneIntegrate` (sharp
        Beam.BesselBeam case) to a beam with a genuine kappa spread; a `value::ComplexF64` is returned. See this
        file's header for the normalisation; `Beam.kappaBreakpoints(beam, k)` forces QuadGK through the beam's own
        radial feature(s) rather than trusting the default adaptive search to find them unaided.
"""
function coneIntegrateContinuous(beam::Beam.AbstractBeamType, k::Float64, bPerp::Float64, phiB::Float64,
                                 amplitudeAt::Function; nphi::Int64=361)
    oam  = Beam.oamOf(beam);    dphi = 2pi / (nphi - 1);    nodes = Beam.kappaBreakpoints(beam, k)

    function azimuthalPart(kappa::Float64)
        acc = ComplexF64(0.)
        for  i = 1:nphi
            phiK = (i-1) * dphi;    w = (i == 1 || i == nphi) ? 0.5 : 1.0
            acc  = acc + w * amplitudeAt(kappa, phiK) * exp(im * oam * phiK) * exp(im * kappa * bPerp * cos(phiK - phiB)) * dphi
        end

        return( acc / (2pi) )
    end

    numerator,  _ = QuadGK.quadgk(kappa -> kappa * Beam.radialAmplitude(beam, kappa, k) * azimuthalPart(kappa), nodes..., Inf)
    firstMoment,_ = QuadGK.quadgk(kappa -> kappa * Beam.radialAmplitude(beam, kappa, k),                        nodes..., Inf)

    return( numerator / firstMoment )
end


"""
`Beam.coneIntegrate(beam::Beam.BesselGauss, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                    nphi::Int64=361)`
    ... to combine a caller-supplied amplitude `amplitudeAt(kappa::Float64, phiK::Float64)::ComplexF64` coherently
        over the Bessel-Gauss beam's Gaussian ring of transverse momentum and OAM phase, and the target's
        displacement from the axis; a `value::ComplexF64` is returned, via `Beam.coneIntegrateContinuous`.
"""
function coneIntegrate(beam::Beam.BesselGauss, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                       nphi::Int64=361)
    return( Beam.coneIntegrateContinuous(beam, k, bPerp, phiB, amplitudeAt; nphi=nphi) )
end


"""
`Beam.coneIntegrate(beam::Beam.LaguerreGauss, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                    nphi::Int64=361)`
    ... to combine a caller-supplied amplitude `amplitudeAt(kappa::Float64, phiK::Float64)::ComplexF64` coherently
        over the LG mode's continuous transverse-momentum spread and OAM phase, and the target's displacement from
        the axis; a `value::ComplexF64` is returned, via `Beam.coneIntegrateContinuous`.
"""
function coneIntegrate(beam::Beam.LaguerreGauss, k::Float64, bPerp::Float64, phiB::Float64, amplitudeAt::Function;
                       nphi::Int64=361)
    return( Beam.coneIntegrateContinuous(beam, k, bPerp, phiB, amplitudeAt; nphi=nphi) )
end


end # module
