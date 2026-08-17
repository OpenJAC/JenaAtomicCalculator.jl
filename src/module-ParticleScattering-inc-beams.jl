
# How the projectile is prepared. The partial waves are computed once and are independent of the beam; the beam enters
# only in how they are superposed into an observable, which is why this is a dispatch layer over the same data rather
# than a second calculation.
#
# The twisted beams are NOT implemented today, deliberately. The expressions this module carried for them were
# nonrelativistic -- one phase shift per l, sin(delta_l) exp(i delta_l) -- and have to be re-derived over both signs of
# kappa before they can be used with the Dirac partial waves. Until then they raise, rather than returning a number from
# a different theory than the one the rest of the module now uses.


"""
`ParticleScattering.beamObservables(beamType::Beam.PlaneWave, pws::Array{ParticleScattering.PartialWave,1},
                                    energy::Float64, theta::Float64, phi::Float64)`
    ... to compute the observables for a PLANE-WAVE projectile, for which the partial-wave superposition is the textbook
        one and the amplitudes are those of ParticleScattering.directAmplitude and .spinFlipAmplitude directly. An
        obs::ParticleScattering.AngularObservables is returned.
"""
function beamObservables(beamType::Beam.PlaneWave, pws::Array{ParticleScattering.PartialWave,1},
                         energy::Float64, theta::Float64, phi::Float64)
    return( ParticleScattering.angularObservables(pws, energy, theta, phi) )
end


"""
`ParticleScattering.beamObservables(beamType::Beam.BesselBeam, pws::Array{ParticleScattering.PartialWave,1},
                                    energy::Float64, theta::Float64, phi::Float64)`
    ... to compute the observables for a Bessel (vortex) projectile beam. NOT IMPLEMENTED: the partial-wave superposition
        for a twisted beam has to be re-derived over both signs of kappa before it can be combined with Dirac phase
        shifts. Nothing is returned; this method always raises.
"""
function beamObservables(beamType::Beam.BesselBeam, pws::Array{ParticleScattering.PartialWave,1},
                         energy::Float64, theta::Float64, phi::Float64)
    error("\n\nScattering of a Bessel (vortex) beam is not implemented in the Dirac treatment yet.\n" *
          "The expressions previously carried here were NONRELATIVISTIC -- a single phase shift per l -- and cannot be " *
          "combined with the two phase shifts delta_{-l-1}, delta_{+l} that the module now computes. They have to be " *
          "re-derived over both signs of kappa first. Use Beam.PlaneWave() meanwhile.")
end


"""
`ParticleScattering.beamObservables(beamType::Beam.AbstractBeamType, pws::Array{ParticleScattering.PartialWave,1},
                                    energy::Float64, theta::Float64, phi::Float64)`
    ... a fallback for every beam type without its own method, which says which one was asked for instead of silently
        falling back on the plane-wave superposition. Nothing is returned; this method always raises.
"""
function beamObservables(beamType::Beam.AbstractBeamType, pws::Array{ParticleScattering.PartialWave,1},
                         energy::Float64, theta::Float64, phi::Float64)
    error("\n\nNo partial-wave superposition is implemented for beamType = $(beamType).\n" *
          "Implemented today: Beam.PlaneWave(). The twisted beams await the Dirac re-derivation described in " *
          "module-ParticleScattering-inc-beams.jl.")
end
