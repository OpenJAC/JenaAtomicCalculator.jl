
println("Oe) POSITRON elastic scattering at helium, beside the electron case: opposite-sign phase shifts and a " *
        "weaker cross section, most of all backwards.")

setDefaults("print summary: open", "zzz-ParticleScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:  17-Aug-2026
    # Last successful:  17-Aug-2026
    # FIRST STEPS WITH A SECOND PROJECTILE. The projectile is one of the module's three orthogonal axes, and a positron
    # needed no new machinery at all: the Dirac radial equation, the partial-wave series, the amplitudes f and g and
    # every observable are the same, and only the POTENTIAL differs. The electrostatic interaction changes sign
    # throughout -- the nucleus repels the positron and the electrons attract it -- and no exchange term applies at all,
    # a positron being distinguishable from the target electrons. That is one method of ParticleScattering.
    # scatteringPotential, and it is the whole implementation.
    #
    #   REPORT (He, 100 eV, rbox = 20 a0):
    #
    #                                  sigma_el   sigma_1    DCS(180 deg)   delta(kappa=-1)   partial waves
    #     electron, Furness-McCarthy   2.0852     0.94324    0.032923       +1.06447           25
    #     electron, static only        1.6311     0.82609    0.032269       +0.96543           25
    #     positron, static only        0.86029    0.32569    0.0079071      -0.60684           25
    #
    #   TWO THINGS TO READ OFF, and both are textbook.
    #
    #   (1) THE PHASE SHIFT CHANGES SIGN. The electron s-wave phase shift is +0.965 in the static field and the
    #       positron's is -0.607: an attractive potential pulls the wave in and advances its phase, a repulsive core
    #       pushes it out and retards it. Nothing enforces this in the code -- it follows from the sign flip in the
    #       potential, and it is the cleanest single signature that the positron path is doing something different.
    #
    #   (2) THE POSITRON SCATTERS LESS, and most of all backwards. In the same static field sigma_el falls from 1.63
    #       to 0.86 a0^2, a factor 1.9, while the BACKWARD cross section falls from 0.0323 to 0.0079, a factor 4.1.
    #       Large-angle scattering is what requires a close approach to the nucleus, and that is exactly what the
    #       positron is prevented from making.
    #
    #   A CLAIM WITHDRAWN, recorded because it was nearly published here. Being repelled, a positron might be expected
    #   to need FEWER partial waves; a first measurement seemed to show 25 against 71. Re-run properly, all three rows
    #   above need 25 at this energy and convergence criterion. The 71 belonged to the Slater-exchange electron, whose
    #   much stronger potential is what lengthens the series -- an interaction-model effect, not a projectile one.
    #
    #   WHAT IS NOT YET HERE. No correlation-polarization potential, which matters MORE for positrons than for
    #   electrons: there is no exchange term to compete with it, so the induced-dipole attraction is the leading
    #   correction to the static field and its omission is felt sooner. Positronium formation, which opens as a channel
    #   below the ionization threshold, is absent altogether and is not a small effect where it is open. These numbers
    #   are therefore a working static-field result, not a prediction to compare with a positron beam experiment, and
    #   no such comparison is claimed.
    #
    #   AND WHY NOT PROTONS. The same axis exists for them and the potential would be as easy -- charge +1, no exchange,
    #   i.e. the positron's potential. What blocks it is that Continuum.generateOrbitalLocalPotential solves the Dirac
    #   equation for a particle of the ELECTRON rest mass; a proton needs its own mass in both the radial equation and
    #   the momentum relation K = sqrt(E(E+2mc^2))/c. That is a change inside the Continuum module rather than here.
    #   A proton also carries ~43 times the momentum of an electron of the same energy, so the partial-wave series
    #   would need of order 43 times as many terms, and a partial-wave treatment stops being the sensible method.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 20.0)
    cases = [ ("electron, Furness-McCarthy", ParticleScattering.Electron(), ParticleScattering.StaticFieldFurnessMcCarthy()),
              ("electron, static only",      ParticleScattering.Electron(), ParticleScattering.StaticField()),
              ("positron, static only",      ParticleScattering.Positron(), ParticleScattering.StaticField()) ]
    println("\n                                sigma_el     sigma_1      DCS(180)     delta(k=-1)   partial waves")
    for  (label, projectile, model) in cases
        psSettings = ParticleScattering.Settings(ParticleScattering.Settings(),
                                                 projectile = projectile, interaction = model,
                                                 impactEnergies = [100.0], polarThetas = [Float64(pi)],
                                                 polarPhis = [0.0], printBefore = false,
                                                 epsPartialWave = 1.0e-7, maxL = 120)
        wc = Atomic.Computation(Atomic.Computation(), name="positron vs electron", grid=grid,
                                nuclearModel    = Nuclear.Model(2.0),
                                initialConfigs  = [Configuration("1s^2")],
                                finalConfigs    = [Configuration("1s^2")],
                                processSettings = psSettings )
        event = perform(wc; output=true)["particle-scattering events:"][1]
        println("  ", rpad(label, 30),
                rpad(round(event.integrated.sigmaElastic, sigdigits=5), 13),
                rpad(round(event.integrated.sigmaMomentumTransfer, sigdigits=5), 13),
                rpad(round(event.angular[1].dcs, sigdigits=5), 13),
                rpad(round(ParticleScattering.phaseShift(event.partialWaves, -1), digits=5), 14),
                length(event.partialWaves))
    end
    #
elseif  false
    # Last successful:  unknown ...
    #
end
#
setDefaults("print summary: close", "")
