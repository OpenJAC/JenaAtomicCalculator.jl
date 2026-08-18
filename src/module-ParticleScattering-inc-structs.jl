
# The data structures of the ParticleScattering module, kept apart along three ORTHOGONAL axes -- the projectile,
# the interaction model and the beam -- so that none of them has to be widened when another one is. The observables
# are a functional of the amplitudes and are computed afterwards; nothing that is declared here is left unfilled.


"""
`abstract type ParticleScattering.AbstractProjectile`
    ... defines an abstract type to distinguish the (massive) projectiles that are scattered at an atom or ion. Photons are
        NOT among them: photon scattering is a second-order matrix element over intermediate atomic states and is computed
        by the RayleighCompton module, which shares no computational core with the potential scattering treated here.

    + struct Electron    ... an electron as projectile.
    + struct Positron    ... a positron as projectile; the static field changes sign and no exchange term applies.
    + struct Proton      ... a proton as projectile.
"""
abstract type  AbstractProjectile                                  end
struct     Electron            <:  AbstractProjectile              end
struct     Positron            <:  AbstractProjectile              end
struct     Proton              <:  AbstractProjectile              end

@doc "... an electron as projectile; the local exchange term of the interaction model applies to it."                    Electron
@doc "... a positron as projectile; the static field enters with the opposite sign and no exchange term applies."        Positron
@doc "... a proton as projectile."                                                                                          Proton


"""
`abstract type ParticleScattering.AbstractScatteringProcess`
    ... defines an abstract type to distinguish the scattering processes; the process fixes which final levels may be
        reached and, with it, how many independent amplitudes the scattering matrix requires.

    + struct ElasticScattering     ... the target is left in its initial level.
    + struct InelasticScattering   ... the target is left excited (not yet implemented).
"""
abstract type  AbstractScatteringProcess                           end
struct     ElasticScattering    <:  AbstractScatteringProcess      end
struct     InelasticScattering  <:  AbstractScatteringProcess      end

@doc "... elastic scattering, i.e. the target is left in its initial level."                                       ElasticScattering
@doc "... inelastic scattering, i.e. the target is left excited; not yet implemented."                           InelasticScattering


"""
`abstract type ParticleScattering.AbstractInteractionModel`
    ... defines an abstract type for how the target is represented to the projectile. This is the axis along which the
        ELSEPA code is organised, and it is kept separate from the projectile and the beam so that a further model is one
        new subtype plus one new method of ParticleScattering.scatteringPotential.

    + struct StaticField           ... the electrostatic interaction with the target density alone.
    + struct StaticFieldSlaterExchange  ... the static field plus JAC's Slater exchange term (not recommended).
    + struct StaticFieldFurnessMcCarthy ... the static field plus the energy-dependent Furness-McCarthy exchange.
"""
abstract type  AbstractInteractionModel                                        end
struct     StaticField                     <:  AbstractInteractionModel        end
struct     StaticFieldSlaterExchange       <:  AbstractInteractionModel        end
struct     StaticFieldFurnessMcCarthy      <:  AbstractInteractionModel        end

@doc "... the electrostatic interaction of the projectile with the nuclear and electronic charge density alone."       StaticField
@doc "... the static field plus the SLATER exchange term carried by JAC's DFS field. That term is designed for a " *
     "BOUND electron and is too strong for a continuum projectile: at 100 eV on He it raises the transport cross " *
     "section to 1.157 a0^2 against 0.826 without it, where the reference calculations give 0.834 (TFD) and 0.928 " *
     "(DHF). Kept because it is what the module used before 17-Aug-2026, not because it is recommended."     StaticFieldSlaterExchange
@doc "... the static field plus the energy-dependent local exchange potential of Furness and McCarthy, J. Phys. B 6, " *
     "2280 (1973), which is the term ELSEPA and the NIST database use. This is the recommended model for electrons " *
     "and the default."                                                                              StaticFieldFurnessMcCarthy


"""
`struct  ParticleScattering.PartialWave`
    ... defines a type for a single scattering partial wave of the projectile.

    + kappa           ::Int64      ... relativistic angular-momentum quantum number of the partial wave.
    + energy          ::Float64    ... kinetic energy of the projectile [a.u.].
    + phaseShift      ::Float64    ... phase shift delta_kappa of this partial wave.

        Only kappa is stored: the orbital and total angular momenta follow from it as l = -kappa-1, j = -kappa-1/2 for
        kappa < 0 and l = kappa, j = kappa-1/2 for kappa > 0, and are obtained with Basics.subshell_l(Subshell(101, kappa))
        and Basics.subshell_j(...). Carrying l alongside kappa would only give the two a chance to disagree.
"""
struct  PartialWave
    kappa             ::Int64
    energy            ::Float64
    phaseShift        ::Float64
end


"""
`ParticleScattering.PartialWave()`  ... constructor for a default ParticleScattering.PartialWave.
"""
function PartialWave()
    PartialWave(-1, 0., 0.)
end


# `Base.show(io::IO, pw::ParticleScattering.PartialWave)`  ... prepares a proper printout of pw::ParticleScattering.PartialWave.
function Base.show(io::IO, pw::ParticleScattering.PartialWave)
    sa = Basics.subshell_l(Subshell(101, pw.kappa));    sb = Basics.subshell_j(Subshell(101, pw.kappa))
    println(io, "kappa = $(pw.kappa)  (l = $sa, j = $sb),  energy = $(pw.energy),  phaseShift = $(pw.phaseShift)")
end


"""
`struct  ParticleScattering.ScatteringChannel`
    ... defines a type for one reduced scattering amplitude in a total-angular-momentum coupled basis,

            < (alpha_f J_f, eps kappa_f) J^pi || T || (alpha_i J_i, eps kappa_i) J^pi > .

        This is the GENERAL object of the module and the one that survives a many-electron target: as soon as the target
        carries J_i /= 0 or is left excited, the scattering matrix is a matrix in the magnetic sublevels and cannot be
        written with two scalar functions. The familiar direct and spin-flip amplitudes f and g are a projection of these
        channels, valid for elastic scattering from a spinless target, and are provided as such.

    + initialKappa    ::Int64            ... kappa of the incoming partial wave.
    + finalKappa      ::Int64            ... kappa of the outgoing partial wave.
    + totalSymmetry   ::LevelSymmetry    ... J^pi of the coupled (target + projectile) system.
    + amplitude       ::ComplexF64       ... the reduced amplitude itself.
"""
struct  ScatteringChannel
    initialKappa      ::Int64
    finalKappa        ::Int64
    totalSymmetry     ::LevelSymmetry
    amplitude         ::ComplexF64
end


"""
`ParticleScattering.ScatteringChannel()`  ... constructor for a default ParticleScattering.ScatteringChannel.
"""
function ScatteringChannel()
    ScatteringChannel(-1, -1, LevelSymmetry(AngularJ64(0), Basics.plus), ComplexF64(0.))
end


# `Base.show(io::IO, ch::ParticleScattering.ScatteringChannel)`  ... prepares a printout of ch::ParticleScattering.ScatteringChannel.
function Base.show(io::IO, ch::ParticleScattering.ScatteringChannel)
    println(io, "kappa_i = $(ch.initialKappa),  kappa_f = $(ch.finalKappa),  J^pi = $(ch.totalSymmetry),  " *
                "amplitude = $(ch.amplitude)")
end


"""
`struct  ParticleScattering.AngularObservables`
    ... defines a type for the observables of a scattering event at one scattering angle (theta, phi).

    + theta           ::Float64    ... polar scattering angle [rad].
    + phi             ::Float64    ... azimuthal scattering angle [rad].
    + dcs             ::Float64    ... differential cross section d sigma / d Omega [a.u.].
    + sherman         ::Float64    ... Sherman function S; defined for elastic scattering from a spinless target and
                                       returned as NaN wherever that projection does not apply.
"""
struct  AngularObservables
    theta             ::Float64
    phi               ::Float64
    dcs               ::Float64
    sherman           ::Float64
end


"""
`ParticleScattering.AngularObservables()`  ... constructor for a default ParticleScattering.AngularObservables.
"""
function AngularObservables()
    AngularObservables(0., 0., 0., NaN)
end


# `Base.show(io::IO, obs::ParticleScattering.AngularObservables)`  ... printout of obs::ParticleScattering.AngularObservables.
function Base.show(io::IO, obs::ParticleScattering.AngularObservables)
    println(io, "theta = $(obs.theta),  phi = $(obs.phi),  dcs = $(obs.dcs),  S = $(obs.sherman)")
end


"""
`struct  ParticleScattering.IntegratedObservables`
    ... defines a type for the angle-integrated cross sections of a scattering event. The two transport cross sections are
        the ones that enter transport and Monte-Carlo simulations, and are what ELSEPA reports alongside the elastic total.

    + sigmaElastic            ::Float64   ... total elastic cross section, integral of the DCS over all angles [a.u.].
    + sigmaMomentumTransfer   ::Float64   ... first transport cross section sigma_1, weighted by (1 - cos theta) [a.u.].
    + sigmaViscosity          ::Float64   ... second transport cross section sigma_2, weighted by sin^2 theta [a.u.].
"""
struct  IntegratedObservables
    sigmaElastic              ::Float64
    sigmaMomentumTransfer     ::Float64
    sigmaViscosity            ::Float64
end


"""
`ParticleScattering.IntegratedObservables()`  ... constructor for a default ParticleScattering.IntegratedObservables.
"""
function IntegratedObservables()
    IntegratedObservables(0., 0., 0.)
end


# `Base.show(io::IO, obs::ParticleScattering.IntegratedObservables)`  ... printout of obs::IntegratedObservables.
function Base.show(io::IO, obs::ParticleScattering.IntegratedObservables)
    println(io, "sigma_el = $(obs.sigmaElastic),  sigma_1 = $(obs.sigmaMomentumTransfer),  " *
                "sigma_2 = $(obs.sigmaViscosity)")
end


"""
`struct  ParticleScattering.Event`
    ... defines a type to collect the data and results of a single scattering event, i.e. of one projectile of one energy
        prepared in one beam and scattered from one initial into one final level.

    + projectile      ::ParticleScattering.AbstractProjectile        ... the scattered particle.
    + process         ::ParticleScattering.AbstractScatteringProcess ... elastic or inelastic.
    + interaction     ::ParticleScattering.AbstractInteractionModel  ... how the target is represented.
    + beamType        ::Beam.AbstractBeamType                        ... how the projectile is prepared.
    + initialLevel    ::Level                                        ... initial level of the target.
    + finalLevel      ::Level                                        ... final level of the target.
    + impactEnergy    ::Float64                                      ... kinetic energy of the projectile [a.u.].
    + partialWaves    ::Array{ParticleScattering.PartialWave,1}      ... the phase shifts that were computed.
    + channels        ::Array{ParticleScattering.ScatteringChannel,1} ... the reduced amplitudes; the primary result.
    + angular         ::Array{ParticleScattering.AngularObservables,1} ... derived, one per requested (theta, phi).
    + integrated      ::ParticleScattering.IntegratedObservables     ... derived, the angle-integrated cross sections.
"""
struct  Event
    projectile        ::ParticleScattering.AbstractProjectile
    process           ::ParticleScattering.AbstractScatteringProcess
    interaction       ::ParticleScattering.AbstractInteractionModel
    beamType          ::Beam.AbstractBeamType
    initialLevel      ::Level
    finalLevel        ::Level
    impactEnergy      ::Float64
    partialWaves      ::Array{ParticleScattering.PartialWave,1}
    channels          ::Array{ParticleScattering.ScatteringChannel,1}
    angular           ::Array{ParticleScattering.AngularObservables,1}
    integrated        ::ParticleScattering.IntegratedObservables
end


"""
`ParticleScattering.Event()`  ... constructor for a default ParticleScattering.Event.
"""
function Event()
    Event( Electron(), ElasticScattering(), StaticFieldFurnessMcCarthy(), Beam.PlaneWave(), Level(), Level(), 0.,
           PartialWave[], ScatteringChannel[], AngularObservables[], IntegratedObservables() )
end


# `Base.show(io::IO, event::ParticleScattering.Event)`  ... prepares a proper printout of event::ParticleScattering.Event.
function Base.show(io::IO, event::ParticleScattering.Event)
    println(io, "projectile:         $(event.projectile)  ")
    println(io, "process:            $(event.process)  ")
    println(io, "interaction:        $(event.interaction)  ")
    println(io, "beamType:           $(event.beamType)  ")
    println(io, "initialLevel:       $(event.initialLevel)  ")
    println(io, "finalLevel:         $(event.finalLevel)  ")
    println(io, "impactEnergy:       $(event.impactEnergy)  ")
    println(io, "partialWaves:       $(length(event.partialWaves)) computed  ")
    println(io, "channels:           $(length(event.channels)) computed  ")
    println(io, "angular:            $(length(event.angular)) angles  ")
    println(io, "integrated:         $(event.integrated)  ")
end
