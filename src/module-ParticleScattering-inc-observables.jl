
# The observables of a scattering event. They are a FUNCTIONAL of the amplitudes and are computed after them, never
# inside the partial-wave loop: that separation is what allows a further observable to be added without touching the
# machinery that produced the phase shifts.


"""
`ParticleScattering.angularObservables(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64,
                                       phi::Float64)`
    ... to compute the differential cross section and the Sherman function at one scattering angle, from the direct and
        spin-flip amplitudes,

            dsigma/dOmega = |f|^2 + |g|^2 ,        S = i (f g* - f* g) / (|f|^2 + |g|^2) .

        An obs::ParticleScattering.AngularObservables is returned.
"""
function angularObservables(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64, theta::Float64, phi::Float64)
    f = ParticleScattering.directAmplitude(pws, energy, theta)
    g = ParticleScattering.spinFlipAmplitude(pws, energy, theta)
    dcs = abs2(f) + abs2(g)
    sherman = (dcs > 0.) ? real( im * (f * conj(g) - conj(f) * g) ) / dcs : 0.

    return( ParticleScattering.AngularObservables(theta, phi, dcs, sherman) )
end


"""
`ParticleScattering.elasticCrossSectionFromPhases(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64)`
    ... to compute the total elastic cross section directly from the phase shifts,

            sigma_el = 4 pi / k^2 SUM_l [ (l+1) sin^2 delta_{-l-1} + l sin^2 delta_{l} ] ,

        i.e. without ever forming the amplitudes. Comparing it with the angular integral of the differential cross section
        is an internal consistency check that costs nothing. A sigma::Float64 [a.u.] is returned.
"""
function elasticCrossSectionFromPhases(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64)
    k = ParticleScattering.projectileMomentum(energy);    wa = 0.
    for  pw in pws
        l      = Basics.subshell_l(Subshell(101, pw.kappa))
        weight = (pw.kappa < 0) ? l + 1 : l
        wa     = wa + weight * sin(pw.phaseShift)^2
    end

    return( 4pi * wa / (k*k) )
end


"""
`ParticleScattering.integratedObservables(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64;
                                          noAngles::Int64=64)`
    ... to compute the angle-integrated cross sections by Gauss-Legendre quadrature over x = cos(theta),

            sigma_el = 2 pi INT dcs dx ,   sigma_1 = 2 pi INT (1-x) dcs dx ,   sigma_2 = 3 pi INT (1-x^2) dcs dx ,

        i.e. the elastic total together with the first (momentum-transfer) and second (viscosity) transport cross sections,
        which are the quantities transport and Monte-Carlo codes consume. An obs::ParticleScattering.IntegratedObservables
        is returned.
"""
function integratedObservables(pws::Array{ParticleScattering.PartialWave,1}, energy::Float64; noAngles::Int64=64)
    gl = Radial.GridGL(Radial.GridGaussLegendreFinite(), -1.0, 1.0, noAngles)
    sigmaEl = 0.;    sigma1 = 0.;    sigma2 = 0.

    for  i = 1:length(gl.t)
        x     = gl.t[i];    w = gl.wt[i];    theta = acos( max(-1.0, min(1.0, x)) )
        f     = ParticleScattering.directAmplitude(pws, energy, theta)
        g     = ParticleScattering.spinFlipAmplitude(pws, energy, theta)
        dcs   = abs2(f) + abs2(g)
        sigmaEl = sigmaEl + w * dcs
        sigma1  = sigma1  + w * dcs * (1.0 - x)
        sigma2  = sigma2  + w * dcs * (1.0 - x*x)
    end

    return( ParticleScattering.IntegratedObservables( 2pi * sigmaEl, 2pi * sigma1, 3pi * sigma2 ) )
end


"""
`ParticleScattering.extractCrossSections(events::Array{ParticleScattering.Event,1}, energy::Float64, phi::Float64)`
    ... to collect, for the given impact energy and azimuthal angle, the scattering angles together with the differential
        cross sections and Sherman functions of all events, in a form convenient for plotting. A named tuple
        (thetas, dcs, sherman)::NamedTuple is returned.
"""
function extractCrossSections(events::Array{ParticleScattering.Event,1}, energy::Float64, phi::Float64)
    thetas = Float64[];    dcs = Float64[];    sherman = Float64[]

    for  event in events
        if  abs(event.impactEnergy - energy) > 1.0e-10    continue    end
        for  obs in event.angular
            if  abs(obs.phi - phi) > 1.0e-10              continue    end
            push!(thetas, obs.theta);    push!(dcs, obs.dcs);    push!(sherman, obs.sherman)
        end
    end

    return( (thetas=thetas, dcs=dcs, sherman=sherman) )
end
