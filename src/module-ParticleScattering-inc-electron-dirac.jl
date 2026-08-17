
# The Dirac partial-wave core: phase shifts for BOTH signs of kappa, the reduced scattering channels built from them,
# and the direct and spin-flip amplitudes f and g as a projection of those channels.
#
# The projection is what ELSEPA computes, and it is valid only for elastic scattering from a spinless target: the
# decomposition M = f + g (sigma . n) follows from rotational invariance and parity when the projectile is spin-1/2 and
# the target neither carries angular momentum nor changes state. The channels above it stay general, so that a target
# with J_i /= 0, or an inelastic transition, can be added without rebuilding the module around a different object.


"""
`ParticleScattering.projectileMomentum(energy::Float64)`
    ... to compute the relativistic momentum of a projectile of rest mass one and kinetic energy `energy` [a.u.], i.e.
        p = sqrt(E (E + 2c^2)) / c = sqrt(2E + E^2/c^2). A momentum::Float64 [a.u.] is returned.
"""
function projectileMomentum(energy::Float64)
    wc = Defaults.getDefaults("speed of light: c")

    return( sqrt( 2*energy + energy*energy/(wc*wc) ) )
end


"""
`ParticleScattering.computePartialWaves(projectile::ParticleScattering.AbstractProjectile,
                                        interaction::ParticleScattering.AbstractInteractionModel, level::Level,
                                        impactEnergy::Float64, nm::Nuclear.Model, grid::Radial.Grid,
                                        contSettings::Continuum.Settings, settings::ParticleScattering.Settings;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)`
    ... to compute the phase shifts delta_kappa of the projectile in the field of the given target level, for BOTH spin-orbit
        partners of every orbital angular momentum: kappa = -l-1 (j = l+1/2) and, for l >= 1, kappa = +l (j = l-1/2). The
        series is extended until the contribution of one l to the elastic cross section falls below settings.epsPartialWave
        twice in succession, or until settings.maxL is reached; the achieved l is reported so that a truncation is visible
        rather than silent. An Array{ParticleScattering.PartialWave,1} is returned.
"""
function computePartialWaves(projectile::ParticleScattering.AbstractProjectile,
                             interaction::ParticleScattering.AbstractInteractionModel, level::Level,
                             impactEnergy::Float64, nm::Nuclear.Model, grid::Radial.Grid,
                             contSettings::Continuum.Settings, settings::ParticleScattering.Settings;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)
    pot = ParticleScattering.scatteringPotential(projectile, interaction, level, nm, grid; nuclearPot=nuclearPot)
    pws = ParticleScattering.PartialWave[];    total = 0.;    quiet = 0;    lReached = 0

    for  l = 0:settings.maxL
        kappas = (l == 0) ? [-1] : [-l - 1, l]
        contribution = 0.
        for  kappa in kappas
            _, phase, _  = Continuum.generateOrbitalLocalPotential(impactEnergy, Subshell(101, kappa), pot, contSettings)
            push!( pws, ParticleScattering.PartialWave(kappa, impactEnergy, phase) )
            # The weight of a partial wave in sigma_el is (2j+1)/2 = l+1 for kappa = -l-1 and l for kappa = +l
            weight       = (kappa < 0) ? l + 1 : l
            contribution = contribution + weight * sin(phase)^2
        end
        total = total + contribution;    lReached = l
        # Two quiet l in succession, so that an accidental zero of a single phase shift does not end the series early
        if  total > 0.  &&  contribution / total < settings.epsPartialWave    quiet = quiet + 1
        else                                                                  quiet = 0
        end
        if  quiet >= 2  &&  l > 2    break    end
    end

    if  printout
        sa = (lReached == settings.maxL) ? "  <-- maxL reached, the series may NOT be converged" : ""
        println(">> ParticleScattering: $(length(pws)) partial waves up to l = $lReached" *
                " (epsPartialWave = $(settings.epsPartialWave))$sa")
    end

    return( pws )
end


"""
`ParticleScattering.phaseShift(pws::Array{ParticleScattering.PartialWave,1}, kappa::Int64)`
    ... to pick the phase shift of the partial wave with the given kappa out of the list; zero is returned if that partial
        wave was not computed, which is the correct limit for an l beyond the converged series. A phase::Float64 is returned.
"""
function phaseShift(pws::Array{ParticleScattering.PartialWave,1}, kappa::Int64)
    for  pw in pws
        if  pw.kappa == kappa    return( pw.phaseShift )    end
    end

    return( 0. )
end


"""
`ParticleScattering.scatteringChannels(pws::Array{ParticleScattering.PartialWave,1}, level::Level)`
    ... to build the reduced scattering channels from the computed partial waves. For elastic scattering the projectile
        keeps its kappa, so one channel arises per partial wave, with the reduced amplitude exp(2 i delta_kappa) - 1 and
        the total symmetry obtained by coupling the partial wave to the target level. An
        Array{ParticleScattering.ScatteringChannel,1} is returned.
"""
function scatteringChannels(pws::Array{ParticleScattering.PartialWave,1}, level::Level)
    channels = ParticleScattering.ScatteringChannel[]
    for  pw in pws
        sh   = Subshell(101, pw.kappa)
        j    = Basics.subshell_j(sh);    l = Basics.subshell_l(sh)
        # Elastic scattering from the level: the coupled symmetry runs over J_i (x) j, and the parity picks up (-1)^l
        parity = iseven(l) ? level.parity : ( level.parity == Basics.plus ? Basics.minus : Basics.plus )
        for  J in oplus(level.J, j)
            amp = exp( 2im * pw.phaseShift ) - 1.0
            push!( channels, ParticleScattering.ScatteringChannel(pw.kappa, pw.kappa, LevelSymmetry(J, parity), amp) )
        end
    end

    return( channels )
end


"""
`ParticleScattering.directAmplitude(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64)`
    ... to compute the direct scattering amplitude

            f(theta) = 1/(2ik) SUM_l { (l+1) [exp(2i delta_{-l-1}) - 1] + l [exp(2i delta_{l}) - 1] } P_l(cos theta)

        from the computed phase shifts. This is the spin-non-flip amplitude of the elastic scattering of a spin-1/2
        projectile from a SPINLESS target; see ParticleScattering.assertSpinlessTarget for where the projection applies.
        An f::ComplexF64 [a.u.] is returned.
"""
function directAmplitude(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64)
    k = ParticleScattering.projectileMomentum(energy);    x = cos(theta);    f = ComplexF64(0.)
    lMax = maximum( Basics.subshell_l(Subshell(101, pw.kappa)) for pw in pws )

    for  l = 0:lMax
        dMinus = ParticleScattering.phaseShift(pws, -l - 1)          # kappa = -l-1,  j = l+1/2
        dPlus  = (l == 0) ? 0. : ParticleScattering.phaseShift(pws, l)   # kappa = +l,   j = l-1/2
        wa     = (l + 1) * ( exp(2im*dMinus) - 1.0 ) + l * ( exp(2im*dPlus) - 1.0 )
        f      = f + wa * GSL.sf_legendre_Pl_e(l, x).val
    end

    return( f / (2im * k) )
end


"""
`ParticleScattering.spinFlipAmplitude(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64)`
    ... to compute the spin-flip scattering amplitude

            g(theta) = 1/(2ik) SUM_l { exp(2i delta_{l}) - exp(2i delta_{-l-1}) } P_l^1(cos theta)

        from the computed phase shifts. It vanishes identically when the two spin-orbit partners share a phase shift, which
        is why a treatment using only kappa = -l-1 can carry no spin dependence at all. A g::ComplexF64 [a.u.] is returned.
"""
function spinFlipAmplitude(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64)
    k = ParticleScattering.projectileMomentum(energy);    x = cos(theta);    g = ComplexF64(0.)
    lMax = maximum( Basics.subshell_l(Subshell(101, pw.kappa)) for pw in pws )

    for  l = 1:lMax
        dMinus = ParticleScattering.phaseShift(pws, -l - 1)
        dPlus  = ParticleScattering.phaseShift(pws, l)
        g      = g + ( exp(2im*dPlus) - exp(2im*dMinus) ) * GSL.sf_legendre_Plm(l, 1, x)
    end

    return( g / (2im * k) )
end


"""
`ParticleScattering.assertSpinlessTarget(event::ParticleScattering.Event)`
    ... to test whether the f/g projection applies to the given event, i.e. whether the process is elastic and the target
        level carries J_i = 0. The decomposition M = f + g (sigma . n) rests on both, and for a target with J_i /= 0 the
        scattering matrix is a matrix in the magnetic sublevels that two scalar functions cannot represent. Nothing is
        returned; the function raises where the projection would silently give a wrong number.
"""
function assertSpinlessTarget(event::ParticleScattering.Event)
    if  !(typeof(event.process) == ParticleScattering.ElasticScattering)
        error("\n\nThe direct/spin-flip amplitudes f and g are defined for ELASTIC scattering only; the present event " *
              "carries process = $(event.process). Use the reduced channels in event.channels instead.")
    end
    if  event.initialLevel.J != AngularJ64(0)
        error("\n\nThe direct/spin-flip amplitudes f and g rest on a SPINLESS target, whereas the present initial level " *
              "carries J = $(event.initialLevel.J). For J /= 0 the scattering matrix is a matrix in the magnetic " *
              "sublevels and cannot be written with two scalar functions; use the reduced channels in event.channels.")
    end

    return( nothing )
end
