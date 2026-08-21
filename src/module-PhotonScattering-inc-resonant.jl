
# Resonant inelastic scattering (RIXS):   gamma(omega_in) + |i>  -->  |f> + gamma(omega_out) ,  through a RESONANT intermediate.
#
# The same second-order amplitude as Rayleigh/Raman, and deliberately the same channel construction -- see
# PhotonScattering.determineSecondOrderChannels, which serves both -- but evaluated where Rayleigh refuses to go.
#
# THE ONE PHYSICAL DIFFERENCE, and it inverts a decision made deliberately in the Rayleigh file.
# PhotonScattering.rayleighAmplitude SKIPS an intermediate level whose energy denominator falls below settings.selfTolerance, on
# the grounds that the perturbative expression is not defined there. RIXS is exactly the case where the incoming photon is TUNED
# ONTO such a level, so skipping it would discard the entire signal. Instead the level is kept and the denominator regularised by
# the resonance width,
#
#     AbsorbThenEmit :   E_i + omega_in - E_nu + i*Gamma/2        <- resonant; the width matters and is why RIXS is finite
#     EmitThenAbsorb :   E_i - omega_out - E_nu                   <- real; this ordering never goes on shell
#
# The width sits in the resonant term ALONE. Emitting first moves AWAY from the intermediate manifold, so that denominator cannot
# vanish and attributing a decay width to it would be attributing a lifetime to a state the process never occupies.
#
# WHAT MAKES THIS A LINE AND NOT A PATHWAY. JAC's older ResonantInelastic module is pathway-shaped -- initial -> intermediate ->
# final, with an explicit intermediateMultiplet supplied by Basics.perform from intermediateConfigs. This implementation is
# LINE-shaped by the maintainer's decision of 21-Aug-2026: the intermediate levels come from settings.gMultiplet like every other
# second-order process here, so there is ONE result type and ONE entry point in the module, and PhotonScattering.Settings is
# deliberately NOT added to the intermediateMultiplet list in module-BasicsAZ-inc-perform.jl.
#
# THE OBSERVABLE IS A MAP, not a number. RIXS is measured by scanning omega_in across the resonance and dispersing omega_out, so
# the natural output is a two-dimensional intensity in (omega_in, omega_out). Here omega_out is FIXED by energy conservation once
# the final level is chosen, so the second axis is the FINAL LEVEL: one Line per (final level, omega_in), and the map is read off
# the list of lines. That is the same arrangement as Rayleigh/Raman and needs no extra machinery.


"""
`PhotonScattering.computeResonantAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                      settings::PhotonScattering.Settings; printout::Bool=true)`
    ... to compute all channel amplitudes of one resonant-scattering line and, from them, its cross section. The channels are those
        of PhotonScattering.determineSecondOrderChannels, shared with Rayleigh/Raman, and each is evaluated by
        PhotonScattering.resonantAmplitude. A newLine::PhotonScattering.Line is returned with its channels and cross section
        evaluated.
"""
function computeResonantAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                             settings::PhotonScattering.Settings; printout::Bool=true)
    newChannels = PhotonScattering.Channel[]
    for  channel in PhotonScattering.determineSecondOrderChannels(line.finalLevel, line.initialLevel, settings)
        amplitude = PhotonScattering.resonantAmplitude(channel, line.finalLevel, line.initialLevel, line.inPhotonEnergy,
                                                        line.outPhotonEnergy, grid, settings; printout=printout)
        push!( newChannels, PhotonScattering.Channel(channel.kappa, channel.inMultipole, channel.outMultipole, channel.gauge,
                                                      channel.timeOrdering, channel.totalSymmetry, amplitude) )
    end

    crossSection = PhotonScattering.rayleighCrossSection(newChannels, line.inPhotonEnergy, line.outPhotonEnergy, line.initialLevel)
    newLine      = PhotonScattering.Line( line.initialLevel, line.finalLevel, line.inPhotonEnergy, line.outPhotonEnergy,
                                          0., crossSection, newChannels, PhotonScattering.Observables[] )

    return( newLine )
end


"""
`PhotonScattering.computeResonantLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                       grid::Radial.Grid, settings::PhotonScattering.Settings; output=true)`
    ... to compute all selected resonant-scattering lines, together with their channel amplitudes and cross sections, and to display
        them. The intermediate basis is checked first, and a zero resonance width is refused rather than worked around. An
        Array{PhotonScattering.Line,1} is returned if output=true, and nothing otherwise.
"""
function computeResonantLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                              settings::PhotonScattering.Settings; output=true)
    println("")
    printstyled("PhotonScattering.computeResonantLines(): The computation of resonant (RIXS) amplitudes starts now ... \n",
                color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    if  settings.width <= 0.
        error("PhotonScattering: ResonantScattering() requires a POSITIVE settings.width. The resonant denominator " *
              "E_i + omega_in - E_nu + i*Gamma/2 is what keeps a RIXS amplitude finite on resonance; with Gamma = 0 the " *
              "result would be governed by however close the scan happens to land on an intermediate level, which is an " *
              "artefact of the energy grid rather than physics. Set the width, or use RayleighScattering() if the " *
              "non-resonant amplitude is what is wanted.")
    end
    lines = PhotonScattering.determineResonantLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotonScattering.displayResonantLines(stdout, lines)    end
    PhotonScattering.checkIntermediateBasis(lines, settings)
    newLines = PhotonScattering.Line[]
    for  line in lines
        push!( newLines, PhotonScattering.computeResonantAmplitudesProperties(line, nm, grid, settings; printout=false) )
    end
    PhotonScattering.displayResonantResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotonScattering.displayResonantResults(iostream, newLines, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotonScattering.determineResonantLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                         settings::PhotonScattering.Settings)`
    ... to determine the list of resonant-scattering lines, i.e. one line per selected pair of levels and per incoming photon energy.
        The outgoing energy follows from energy conservation, omega_out = omega_in - (E_f - E_i), and a line whose outgoing energy
        is non-positive is dropped. Unlike Rayleigh the final level is normally DIFFERENT from the initial one -- an elastic final
        level is not excluded, that being the resonant-elastic channel. An Array{PhotonScattering.Line,1} is returned with the
        amplitudes not yet evaluated.
"""
function determineResonantLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotonScattering.Settings)
    lines = PhotonScattering.Line[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            for  omega in settings.photonEnergies
                inOmega  = Defaults.convertUnits("energy: to atomic", omega)
                outOmega = inOmega - (fLevel.energy - iLevel.energy)
                if  outOmega <= 0.    continue    end
                push!( lines, PhotonScattering.Line(iLevel, fLevel, inOmega, outOmega, 0., EmProperty(0., 0.),
                                                     PhotonScattering.Channel[], PhotonScattering.Observables[]) )
            end
        end
    end

    return( lines )
end


"""
`PhotonScattering.displayResonantLines(stream::IO, lines::Array{PhotonScattering.Line,1})`
    ... to list the selected resonant-scattering lines before their amplitudes are computed. A neat table is printed but nothing is
        returned otherwise.
"""
function displayResonantLines(stream::IO, lines::Array{PhotonScattering.Line,1})
    nx = 104
    println(stream, " ")
    println(stream, "  Selected resonant (RIXS) scattering lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(16, "Energy in"; na=4) * TableStrings.center(16, "Energy out"; na=4)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy))  * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.outPhotonEnergy))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotonScattering.displayResonantResults(stream::IO, lines::Array{PhotonScattering.Line,1},
                                         settings::PhotonScattering.Settings)`
    ... to print the computed resonant cross sections. The table also carries the DETUNING of each line, i.e. how far the incoming
        photon sits from the nearest intermediate level of its own symmetry, since a RIXS number is meaningless without it: the same
        cross section means quite different things on and off resonance. A neat table is printed but nothing is returned otherwise.
"""
function displayResonantResults(stream::IO, lines::Array{PhotonScattering.Line,1}, settings::PhotonScattering.Settings)
    nx = 148
    println(stream, " ")
    println(stream, "  Resonant (RIXS) scattering cross sections   [width Gamma = " *
                    @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", settings.width)) * " " *
                    TableStrings.inUnits("energy") * "]:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(16, "Energy in"; na=4) * TableStrings.center(16, "Energy out"; na=4) *
                TableStrings.center(16, "Detuning"; na=4) * TableStrings.center(30, "Cross section  (Coulomb, Babushkin)"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        detuning = PhotonScattering.resonantDetuning(line, settings)
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy))  * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.outPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%+.6e", Defaults.convertUnits("energy: from atomic", detuning))            * "    "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "  " * @sprintf("%.6e", line.crossSection.Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotonScattering.resonantAmplitude(channel::PhotonScattering.Channel, finalLevel::Level, initialLevel::Level,
                                    inOmega::Float64, outOmega::Float64, grid::Radial.Grid,
                                    settings::PhotonScattering.Settings; printout::Bool=true)`
    ... to compute the second-order amplitude of one resonant channel. It differs from PhotonScattering.rayleighAmplitude in exactly
        one respect: the resonant (AbsorbThenEmit) denominator carries the width, E_i + omega_in - E_nu + i*Gamma/2, and NO
        intermediate level is skipped. An amplitude::ComplexF64 is returned.

        The off-shell correction of PhotonScattering.offShellFactor applies here unchanged and for the same reason: the vertices are
        evaluated at a free photon energy, where the length-velocity identity does not hold.
"""
function resonantAmplitude(channel::PhotonScattering.Channel, finalLevel::Level, initialLevel::Level,
                           inOmega::Float64, outOmega::Float64, grid::Radial.Grid,
                           settings::PhotonScattering.Settings; printout::Bool=true)
    U        = ComplexF64(0.)
    nuLevels = PhotonScattering.intermediateLevels(settings.gMultiplet, channel.totalSymmetry)

    for  nuLevel in nuLevels
        if      channel.timeOrdering == PhotonScattering.AbsorbThenEmit()
            # THE RESONANT TERM. No level is skipped; the width keeps it finite where Rayleigh would have to skip.
            denom = ComplexF64( initialLevel.energy + inOmega - nuLevel.energy, 0.5 * settings.width )
            wa = PhotoEmission.amplitude(Basics.Emission(),   channel.outMultipole, channel.gauge, outOmega,
                                         finalLevel, nuLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, outOmega, nuLevel, finalLevel)
            wb = PhotoEmission.amplitude(Basics.Absorption(), channel.inMultipole,  channel.gauge, inOmega,
                                         nuLevel, initialLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, inOmega, nuLevel, initialLevel)
        else
            # The NON-resonant ordering: emitting first moves away from the intermediate manifold, so this denominator cannot
            # vanish and carries no width. It is still guarded, since a pathological level set could bring it near zero.
            denomReal = initialLevel.energy - outOmega - nuLevel.energy
            if  abs(denomReal) < settings.selfTolerance    continue    end
            denom = ComplexF64(denomReal, 0.)
            wa = PhotoEmission.amplitude(Basics.Absorption(), channel.inMultipole,  channel.gauge, inOmega,
                                         finalLevel, nuLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, inOmega, finalLevel, nuLevel)
            wb = PhotoEmission.amplitude(Basics.Emission(),   channel.outMultipole, channel.gauge, outOmega,
                                         nuLevel, initialLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, outOmega, initialLevel, nuLevel)
        end
        U = U + wa * wb / denom
    end

    return( U )
end


"""
`PhotonScattering.resonantDetuning(line::PhotonScattering.Line, settings::PhotonScattering.Settings)`
    ... to compute how far the incoming photon of the given line sits from the NEAREST intermediate level that its channels can
        reach, i.e. min over nu of (E_i + omega_in - E_nu). A detuning::Float64 [a.u.] is returned, zero if no intermediate level is
        available at all.

        This is printed beside every resonant cross section deliberately. A RIXS number without its detuning is close to
        meaningless -- the same cross section means one thing on resonance and quite another two widths away -- and the detuning is
        the quantity a reader needs in order to know which of the two they are looking at.
"""
function resonantDetuning(line::PhotonScattering.Line, settings::PhotonScattering.Settings)
    best = 0.;   found = false
    for  channel in PhotonScattering.determineSecondOrderChannels(line.finalLevel, line.initialLevel, settings)
        for  nuLevel in PhotonScattering.intermediateLevels(settings.gMultiplet, channel.totalSymmetry)
            d = line.initialLevel.energy + line.inPhotonEnergy - nuLevel.energy
            if  !found  ||  abs(d) < abs(best)    best = d;   found = true    end
        end
    end

    return( best )
end
