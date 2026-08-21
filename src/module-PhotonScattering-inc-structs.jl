
# The data structures of PhotonScattering. They come first because Settings names them in its field types.
#
# The module is organized along THREE ORTHOGONAL AXES, deliberately the same arrangement that ParticleScattering uses for massive
# projectiles and that proved to carry a second projectile (the positron) without any new machinery:
#
#     beam            how the photon is prepared            plane wave, Bessel/twisted, and its polarization
#     process         what happens to the target            Rayleigh, Compton, Raman, resonant inelastic, pair creation
#     approximation   how the amplitude is evaluated        second-order Green sum, form factor, independent particle
#
# The third axis is what lets a cheap and an expensive calculation of the SAME process stand side by side and be compared, which is
# the check a scattering code most often needs and most rarely has.


"""
`abstract type PhotonScattering.AbstractPhotonProcess`
    ... defines an abstract type to distinguish what the scattered photon leaves behind; the process fixes which final levels may be
        reached and, with them, which amplitudes have to be built.

    + struct RayleighScattering    ... elastic scattering, the target returning to its initial level.
    + struct ComptonScattering     ... inelastic scattering, the target being left excited or ionized.
    + struct RamanScattering       ... inelastic scattering into a discrete final level of the same ion.
    + struct ResonantScattering    ... scattering through a resonant intermediate state (RIXS).
    + struct BoundFreePairCreation ... the photon creates an electron-positron pair, the ELECTRON being captured into a bound orbital
                                       of the ion and the POSITRON leaving as a continuum particle.
"""
abstract type  AbstractPhotonProcess                                end
struct     RayleighScattering     <:  AbstractPhotonProcess         end
struct     ComptonScattering      <:  AbstractPhotonProcess         end
struct     RamanScattering        <:  AbstractPhotonProcess         end
struct     ResonantScattering     <:  AbstractPhotonProcess         end
struct     BoundFreePairCreation  <:  AbstractPhotonProcess         end

@doc "... elastic (Rayleigh) scattering; the target returns to its initial level."                            RayleighScattering
@doc "... inelastic (Compton) scattering; the target is left excited or ionized."                             ComptonScattering
@doc "... Raman scattering, i.e. inelastic scattering into a discrete final level of the same ion."            RamanScattering
@doc "... resonant inelastic scattering (RIXS) through a resonant intermediate state."                        ResonantScattering
@doc "... bound-free pair creation, gamma + |i(N)> --> |f(N+1)> + e^+: the photon creates a pair, the ELECTRON " *
     "being captured into a bound orbital of the ion while the POSITRON leaves as a continuum particle. It is the " *
     "crossing partner of the one-photon annihilation in ParticleScattering, sharing its charge-conjugate positron " *
     "wave and its photon operator, and it has a LOWER threshold than free pair creation -- omega > 2 m c^2 - B -- " *
     "the binding energy B gained by the captured electron paying part of the rest mass."                     BoundFreePairCreation


"""
`abstract type PhotonScattering.AbstractScatteringApproximation`
    ... defines an abstract type for how a scattering amplitude is evaluated. This is the axis along which the same process may be
        computed at two very different costs, so that the cheap result can be checked against the expensive one.

    + struct SecondOrderGreen      ... the full second-order amplitude with a Green-function sum over intermediate states.
    + struct FormFactorApproximation ... the non-relativistic form-factor (and incoherent scattering function) approximation.
    + struct FirstOrderVertex      ... a single photon vertex, which is exact for the processes that have only one.
"""
abstract type  AbstractScatteringApproximation                        end
struct     SecondOrderGreen        <:  AbstractScatteringApproximation  end
struct     FormFactorApproximation <:  AbstractScatteringApproximation  end
struct     FirstOrderVertex        <:  AbstractScatteringApproximation  end

@doc "... the full second-order amplitude, summing over intermediate states with an AtomicState Green function."   SecondOrderGreen
@doc "... the non-relativistic form-factor approximation, i.e. the Rayleigh amplitude built from the atomic form " *
     "factor and the Compton one from the incoherent scattering function; cheap, and the standard against which a " *
     "full calculation is read."                                                                          FormFactorApproximation
@doc "... a single photon vertex. This is not an approximation for every process: bound-free pair creation HAS only " *
     "one vertex, so FirstOrderVertex is exact there and is its only sensible setting."                          FirstOrderVertex


"""
`struct  PhotonScattering.Channel`
    ... defines a type for one reduced amplitude of a photon-scattering event, i.e. for one combination of the quantum numbers over
        which the amplitude is decomposed. Not every field is meaningful for every process: a one-vertex process such as bound-free
        pair creation leaves outMultipole unused, carrying the incoming multipole alone.

    + kappa           ::Int64            ... kappa of the continuum electron or positron, where the process has one.
    + inMultipole     ::EmMultipole      ... multipole of the incoming photon.
    + outMultipole    ::EmMultipole      ... multipole of the outgoing photon; equal to inMultipole where there is none.
    + gauge           ::EmGauge          ... gauge in which this amplitude was evaluated.
    + totalSymmetry   ::LevelSymmetry    ... J^pi of the coupled system; amplitudes sharing it interfere and are summed coherently.
    + amplitude       ::ComplexF64       ... the reduced amplitude itself.
"""
struct  Channel
    kappa             ::Int64
    inMultipole       ::EmMultipole
    outMultipole      ::EmMultipole
    gauge             ::EmGauge
    totalSymmetry     ::LevelSymmetry
    amplitude         ::ComplexF64
end


"""
`PhotonScattering.Channel()`  ... constructor for a default PhotonScattering.Channel.
"""
function Channel()
    Channel(-1, E1, E1, Basics.Coulomb, LevelSymmetry(AngularJ64(0), Basics.plus), ComplexF64(0.))
end


# `Base.show(io::IO, ch::PhotonScattering.Channel)`  ... prepares a proper printout of the variable ch::PhotonScattering.Channel.
function Base.show(io::IO, ch::PhotonScattering.Channel)
    println(io, "kappa = $(ch.kappa),  in = $(ch.inMultipole),  out = $(ch.outMultipole),  gauge = $(ch.gauge),  " *
                "J^pi = $(ch.totalSymmetry),  amplitude = $(ch.amplitude)")
end


"""
`struct  PhotonScattering.Observables`
    ... defines a type for the observables of a photon-scattering event in one direction (theta, phi). The Stokes parameters describe
        the OUTGOING photon and are what distinguishes a photon-scattering measurement from a bare cross section.

    + theta           ::Float64    ... polar angle of the outgoing photon or positron [rad].
    + phi             ::Float64    ... azimuthal angle [rad].
    + dcs             ::EmProperty ... differential cross section d sigma / d Omega [a.u.], in both gauges.
    + P1              ::Float64    ... Stokes parameter P1, linear polarization along the reaction plane.
    + P2              ::Float64    ... Stokes parameter P2, linear polarization at 45 degrees to it.
    + P3              ::Float64    ... Stokes parameter P3, circular polarization.
"""
struct  Observables
    theta             ::Float64
    phi               ::Float64
    dcs               ::EmProperty
    P1                ::Float64
    P2                ::Float64
    P3                ::Float64
end


"""
`PhotonScattering.Observables()`  ... constructor for a default PhotonScattering.Observables.
"""
function Observables()
    Observables(0., 0., EmProperty(0., 0.), 0., 0., 0.)
end


# `Base.show(io::IO, obs::PhotonScattering.Observables)`  ... prepares a proper printout of obs::PhotonScattering.Observables.
function Base.show(io::IO, obs::PhotonScattering.Observables)
    println(io, "theta = $(obs.theta),  phi = $(obs.phi),  dcs = $(obs.dcs),  " *
                "P1 = $(obs.P1),  P2 = $(obs.P2),  P3 = $(obs.P3)")
end


"""
`struct  PhotonScattering.Line`
    ... defines a type for one photon-scattering event: one incoming photon energy, one initial and one final level, together with the
        amplitudes and observables computed for it.

    + initialLevel    ::Level            ... initial level of the target.
    + finalLevel      ::Level            ... final level reached.
    + inPhotonEnergy  ::Float64          ... energy of the incoming photon [a.u.].
    + outPhotonEnergy ::Float64          ... energy of the outgoing photon [a.u.]; zero for a process that emits none.
    + particleEnergy  ::Float64          ... kinetic energy of the outgoing continuum particle [a.u.]; zero where there is none.
    + crossSection    ::EmProperty       ... total cross section [a.u.], in both gauges.
    + channels        ::Array{PhotonScattering.Channel,1}      ... the reduced amplitudes; the primary result.
    + observables     ::Array{PhotonScattering.Observables,1}  ... derived, one per requested (theta, phi).
"""
struct  Line
    initialLevel      ::Level
    finalLevel        ::Level
    inPhotonEnergy    ::Float64
    outPhotonEnergy   ::Float64
    particleEnergy    ::Float64
    crossSection      ::EmProperty
    channels          ::Array{PhotonScattering.Channel,1}
    observables       ::Array{PhotonScattering.Observables,1}
end


"""
`PhotonScattering.Line()`  ... constructor for a default PhotonScattering.Line.
"""
function Line()
    Line( Level(), Level(), 0., 0., 0., EmProperty(0., 0.), PhotonScattering.Channel[], PhotonScattering.Observables[] )
end


# `Base.show(io::IO, line::PhotonScattering.Line)`  ... prepares a proper printout of the variable line::PhotonScattering.Line.
function Base.show(io::IO, line::PhotonScattering.Line)
    println(io, "initialLevel:       $(line.initialLevel.index)  ")
    println(io, "finalLevel:         $(line.finalLevel.index)  ")
    println(io, "inPhotonEnergy:     $(line.inPhotonEnergy)  ")
    println(io, "outPhotonEnergy:    $(line.outPhotonEnergy)  ")
    println(io, "particleEnergy:     $(line.particleEnergy)  ")
    println(io, "crossSection:       $(line.crossSection)  ")
    println(io, "channels:           $(line.channels)  ")
    println(io, "observables:        $(line.observables)  ")
end
