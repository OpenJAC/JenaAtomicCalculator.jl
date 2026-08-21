
# Rayleigh and Raman scattering:   gamma(omega_in) + |i>  -->  |f> + gamma(omega_out) ,   omega_out fixed by energy conservation.
#
# A FRESH IMPLEMENTATION, not the older JAC.RayleighCompton moved here. That module was found on 21-Aug-2026 not to compute: its
# sum over intermediate states runs for the first level only -- `if ig != 1 continue end` at module-RayleighCompton.jl:338, under
# a comment reading "Sum over all terms of the Green function channel" -- so a second-order amplitude is reduced to its first term,
# and its pole branch throws a MethodError besides. examples/example-Dg.jl branch a carries the full diagnosis. The maintainer's
# decision was that getting this right is a new implementation rather than a repair.
#
# THE PATTERN FOLLOWED IS MultiPhotonTransition's two-photon amplitude, which does the same kind of sum, is validated across eight
# dated example branches, and differs from RayleighCompton in the one respect that matters: a RESONANT intermediate level, where
# the energy denominator vanishes, is SKIPPED by tolerance rather than integrated over as a pole. That is the honest boundary
# between non-resonant and resonant scattering -- the perturbative expression is simply not defined there -- and it removes the
# machinery that broke.
#
# THE TWO TIME-ORDERINGS are kept apart, because they carry different denominators and cannot be folded into one term:
#
#     AbsorbThenEmit :   |i> --absorb omega_in-->  |nu> --emit omega_out--> |f>     denominator  E_i + omega_in  - E_nu
#     EmitThenAbsorb :   |i> --emit omega_out-->   |nu> --absorb omega_in--> |f>    denominator  E_i - omega_out - E_nu
#
# The sign is worth care rather than confidence: MultiPhotonTransition carries a long note on a sign bug of exactly this kind,
# found 07-Aug-2026, where an absorption denominator had been used in the emission file. It produced SPURIOUS POLES, a spike where
# the true spectrum has a smooth maximum, and a gauge ratio that was Z-independent -- the signature of a structural error rather
# than a relativistic one. Both orderings appear here, so both signs must be right, and the symptom to watch for is the same.
#
# ENERGY CONSERVATION fixes the outgoing photon, so nothing is free:  omega_out = omega_in - (E_f - E_i).  For RAYLEIGH the final
# level is the initial one and omega_out = omega_in exactly; for RAMAN the final level is another discrete level of the same ion.
# The CONTINUUM Compton profile is a different calculation and is not implemented; see the module docstring.


"""
`PhotonScattering.computeRayleighAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                      settings::PhotonScattering.Settings; printout::Bool=true)`
    ... to compute all channel amplitudes of one Rayleigh/Raman line and, from them, its total cross section. Each channel is one
        combination of incoming multipole, outgoing multipole, intermediate symmetry, gauge and time ordering, and is evaluated by
        PhotonScattering.rayleighAmplitude. A newLine::PhotonScattering.Line is returned with its channels and cross section
        evaluated.
"""
function computeRayleighAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                             settings::PhotonScattering.Settings; printout::Bool=true)
    newChannels = PhotonScattering.Channel[]
    for  channel in PhotonScattering.determineRayleighChannels(line.finalLevel, line.initialLevel, settings)
        amplitude = PhotonScattering.rayleighAmplitude(channel, line.finalLevel, line.initialLevel, line.inPhotonEnergy,
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
`PhotonScattering.computeRayleighLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                       grid::Radial.Grid, settings::PhotonScattering.Settings; output=true)`
    ... to compute all selected Rayleigh/Raman lines, together with their channel amplitudes and cross sections, and to display
        them. Before any amplitude is evaluated the intermediate basis is checked against the symmetries the channels will demand,
        so that a basis which does not span them is reported at once rather than discovered deep inside a sum. An
        Array{PhotonScattering.Line,1} is returned if output=true, and nothing otherwise.
"""
function computeRayleighLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                              settings::PhotonScattering.Settings; output=true)
    println("")
    printstyled("PhotonScattering.computeRayleighLines(): The computation of Rayleigh/Raman amplitudes starts now ... \n",
                color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    lines = PhotonScattering.determineRayleighLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotonScattering.displayRayleighLines(stdout, lines)    end
    PhotonScattering.checkIntermediateBasis(lines, settings)
    newLines = PhotonScattering.Line[]
    for  line in lines
        push!( newLines, PhotonScattering.computeRayleighAmplitudesProperties(line, nm, grid, settings; printout=false) )
    end
    PhotonScattering.displayRayleighResults(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotonScattering.displayRayleighResults(iostream, newLines)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotonScattering.checkIntermediateBasis(lines::Array{PhotonScattering.Line,1}, settings::PhotonScattering.Settings)`
    ... to check, BEFORE any amplitude is evaluated, that settings.gMultiplet spans every intermediate symmetry the channels of the
        given lines will ask for, and to warn once per missing symmetry if it does not. Nothing is returned.

        This exists because of how the old module failed: RayleighCompton raised `error("stop a: Green channel not found for
        symmetry ...")` from inside computeChannelAmplitude, after all the setup and the entire channel table had been printed, so
        the user learned of an incomplete intermediate basis at the latest possible moment and in the least informative place. A
        missing symmetry is a property of the INPUT and is knowable before a single radial integral is done.
"""
function checkIntermediateBasis(lines::Array{PhotonScattering.Line,1}, settings::PhotonScattering.Settings)
    wanted = LevelSymmetry[]
    for  line in lines
        for  channel in PhotonScattering.determineRayleighChannels(line.finalLevel, line.initialLevel, settings)
            if  !(channel.totalSymmetry in wanted)    push!(wanted, channel.totalSymmetry)    end
        end
    end
    missingSyms = LevelSymmetry[]
    for  symt in wanted
        if  length( PhotonScattering.intermediateLevels(settings.gMultiplet, symt) ) == 0    push!(missingSyms, symt)    end
    end
    if  length(missingSyms) > 0
        @warn "PhotonScattering: settings.gMultiplet spans no level of the intermediate symmetries $missingSyms, which the " *
              "channels of this computation require. Those channels contribute exactly zero, which is an artefact of the " *
              "intermediate basis and NOT a selection rule. Extend gMultiplet, or read the result as a partial sum."
    else
        println(">> PhotonScattering: the intermediate basis spans all $(length(wanted)) required symmetries.")
    end

    return( nothing )
end


"""
`PhotonScattering.determineRayleighChannels(finalLevel::Level, initialLevel::Level, settings::PhotonScattering.Settings)`
    ... to determine the allowed (incoming multipole, outgoing multipole, intermediate symmetry, gauge, time ordering) channels of
        one Rayleigh/Raman line. The intermediate symmetries are those reachable from the initial symmetry by the incoming
        multipole AND from which the final symmetry is reachable by the outgoing one, which is what
        AngularMomentum.allowedTotalSymmetries computes. Both time orderings are generated for every such combination. An
        Array{PhotonScattering.Channel,1} is returned with the amplitudes not yet evaluated.
"""
function determineRayleighChannels(finalLevel::Level, initialLevel::Level, settings::PhotonScattering.Settings)
    channels = PhotonScattering.Channel[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    orderings = PhotonScattering.AbstractTimeOrdering[ PhotonScattering.AbsorbThenEmit(), PhotonScattering.EmitThenAbsorb() ]
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            for  symt in AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                for  gauge in settings.gauges
                    if      mp1.electric  &&  mp2.electric  &&  gauge == Basics.UseCoulomb      gaugex = Basics.Coulomb
                    elseif  mp1.electric  &&  mp2.electric  &&  gauge == Basics.UseBabushkin    gaugex = Basics.Babushkin
                    elseif  !mp1.electric && !mp2.electric  &&  gauge == settings.gauges[1]     gaugex = Basics.Magnetic
                    else    continue
                    end
                    for  ordering in orderings
                        push!( channels, PhotonScattering.Channel(0, mp1, mp2, gaugex, ordering, symt, ComplexF64(0.)) )
                    end
                end
            end
        end
    end

    return( channels )
end


"""
`PhotonScattering.determineRayleighLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                         settings::PhotonScattering.Settings)`
    ... to determine the list of Rayleigh/Raman lines to be computed, i.e. one line per selected pair of levels and per incoming
        photon energy of settings.photonEnergies. The OUTGOING energy is not free: it follows from energy conservation as
        omega_out = omega_in - (E_f - E_i), and a line whose outgoing energy comes out non-positive is dropped, that channel being
        closed. An Array{PhotonScattering.Line,1} is returned with the amplitudes not yet evaluated.
"""
function determineRayleighLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotonScattering.Settings)
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
`PhotonScattering.displayRayleighLines(stream::IO, lines::Array{PhotonScattering.Line,1})`
    ... to list the selected Rayleigh/Raman lines before their amplitudes are computed, so that a long run can be checked against
        what was intended. A neat table is printed but nothing is returned otherwise.
"""
function displayRayleighLines(stream::IO, lines::Array{PhotonScattering.Line,1})
    nx = 104
    println(stream, " ")
    println(stream, "  Selected Rayleigh/Raman scattering lines:")
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
`PhotonScattering.displayRayleighResults(stream::IO, lines::Array{PhotonScattering.Line,1})`
    ... to print the computed Rayleigh/Raman cross sections, one row per line and with both gauges beside each other so that their
        agreement can be read off directly. A neat table is printed but nothing is returned otherwise.
"""
function displayRayleighResults(stream::IO, lines::Array{PhotonScattering.Line,1})
    nx = 130
    println(stream, " ")
    println(stream, "  Rayleigh/Raman scattering cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(16, "Energy in"; na=4) * TableStrings.center(16, "Energy out"; na=4) *
                TableStrings.center(30, "Cross section  (Coulomb, Babushkin)"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy))  * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.outPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "  " * @sprintf("%.6e", line.crossSection.Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotonScattering.rayleighAmplitude(channel::PhotonScattering.Channel, finalLevel::Level, initialLevel::Level,
                                    inOmega::Float64, outOmega::Float64, grid::Radial.Grid,
                                    settings::PhotonScattering.Settings; printout::Bool=true)`
    ... to compute the second-order amplitude of one Rayleigh/Raman channel, i.e. the sum over those intermediate levels of
        settings.gMultiplet that carry the channel's symmetry,

            SUM_nu   <f || O(mp_out, omega_out) || nu> <nu || O(mp_in, omega_in) || i>  /  denominator(nu) ,

        the two vertices and the denominator both depending on the channel's time ordering. A resonant level, where the denominator
        falls below settings.selfTolerance, is SKIPPED -- the perturbative expression is not defined there. An
        amplitude::ComplexF64 is returned; it is exactly zero if the intermediate basis carries no level of this symmetry, which
        PhotonScattering.checkIntermediateBasis warns about before any of this is reached.
"""
function rayleighAmplitude(channel::PhotonScattering.Channel, finalLevel::Level, initialLevel::Level,
                           inOmega::Float64, outOmega::Float64, grid::Radial.Grid,
                           settings::PhotonScattering.Settings; printout::Bool=true)
    U        = ComplexF64(0.)
    nuLevels = PhotonScattering.intermediateLevels(settings.gMultiplet, channel.totalSymmetry)

    for  nuLevel in nuLevels
        if      channel.timeOrdering == PhotonScattering.AbsorbThenEmit()
            # |i> absorbs omega_in to reach |nu>, which then emits omega_out to reach |f>
            denom = initialLevel.energy + inOmega - nuLevel.energy
            if  abs(denom) < settings.selfTolerance    continue    end
            wa = PhotoEmission.amplitude(Basics.Emission(),   channel.outMultipole, channel.gauge, outOmega,
                                         finalLevel, nuLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, outOmega, nuLevel, finalLevel)
            wb = PhotoEmission.amplitude(Basics.Absorption(), channel.inMultipole,  channel.gauge, inOmega,
                                         nuLevel, initialLevel, grid; display=false, printout=printout) *
                 PhotonScattering.offShellFactor(channel.gauge, inOmega, nuLevel, initialLevel)
        else
            # |i> emits omega_out to reach |nu>, which then absorbs omega_in to reach |f>
            denom = initialLevel.energy - outOmega - nuLevel.energy
            if  abs(denom) < settings.selfTolerance    continue    end
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
`PhotonScattering.offShellFactor(gauge::EmGauge, omega::Float64, levelA::Level, levelB::Level)`
    ... to correct a photon vertex that is evaluated OFF SHELL, i.e. at a photon energy other than the transition energy of the two
        levels it connects. A factor::Float64 is returned, to multiply the amplitude by; it is exactly 1 for the Coulomb and
        Magnetic gauges and for an on-shell vertex.

        WHY THIS IS NEEDED AT ALL. Length-velocity equivalence is an ON-SHELL identity: the forms are related through
        [H, r] = -i p / m, and turning <f|p|i> into omega <f|r|i> uses (E_f - E_i) = omega. In a second-order sum the photon energy
        is a FREE parameter and the intermediate level is virtual, so both vertices are evaluated off shell, where that step is
        false. Measurement makes the consequence exact: over omega = 0.5 ... 4 eV the Coulomb amplitude is flat to four digits
        while the BABUSHKIN one is strictly proportional to omega, the ratio falling by 0.5000 per doubling. The structural source
        is InteractionStrength.multipoleTransition, whose length branch carries j_L(qr) against the velocity branch's j_L(qr)/(qr)
        and j'_L(qr): for small argument j_1(x) ~ x/3 against j_1(x)/x ~ 1/3.

        WHAT IS DONE. The length vertex is proportional to the photon energy, so evaluating it at omega rather than at its own
        transition energy dE scales it by omega/dE. Multiplying by dE/omega therefore returns it to the point where the identity
        holds, while the energy DENOMINATOR of the sum keeps the true photon energy, which is where the physics of the virtual
        intermediate state lives.

        THE APPROXIMATION IS THE LONG-WAVELENGTH LIMIT: it assumes retardation across the orbital is negligible, so that the vertex
        depends on omega only through that overall factor and not through the shape of j_L(qr). For E1 at a few eV that is amply
        true; it would NOT be for the ~1 MeV photons of the pair-creation and annihilation modules, where q a_0 ~ 274 and the
        Bessel functions oscillate across the orbital. This correction is therefore local to Rayleigh/Raman and must not be copied
        to those files without rederivation.

        A CONSEQUENCE WORTH KNOWING: after this correction the two gauges no longer differ by a power of omega, so the
        Coulomb/Babushkin comparison becomes what a gauge check is supposed to be -- a measure of WAVE-FUNCTION QUALITY. For the
        test system it should reproduce the on-shell disagreement, which was measured as 0.191 and 0.738 and is correlation-limited
        rather than a defect.
"""
function offShellFactor(gauge::EmGauge, omega::Float64, levelA::Level, levelB::Level)
    if  gauge != Basics.Babushkin    return( 1.0 )    end
    dE = abs(levelA.energy - levelB.energy)
    if  omega < 1.0e-12  ||  dE < 1.0e-12    return( 1.0 )    end

    return( dE / omega )
end


"""
`PhotonScattering.rayleighCrossSection(channels::Array{PhotonScattering.Channel,1}, inOmega::Float64, outOmega::Float64,
                                       initialLevel::Level)`
    ... to form the total Rayleigh/Raman cross section from the evaluated channel amplitudes. Amplitudes feeding the SAME
        intermediate symmetry interfere and are summed coherently first -- across both multipoles and both time orderings, which is
        the point of keeping the orderings as separate channels rather than separate results; the squared moduli are then summed
        over the symmetries, the gauges kept apart, and divided by the statistical weight of the initial level,

            sigma  =  (8 pi / 3) (omega_in omega_out^3 / c^4) / (2 J_i + 1)  *  SUM_symt | SUM_{mp_in, mp_out, ordering} A |^2 .

        A crossSection::EmProperty [a.u.] is returned.

        THE PREFACTOR IS NOT INDEPENDENTLY VERIFIED, in common with the other cross sections of this module and its crossing
        partner in ParticleScattering. What CAN be checked without it is the shape: the factor omega_in omega_out^3 is the Rayleigh
        law, so for ELASTIC scattering (omega_out = omega_in) far below any resonance the cross section must rise as omega^4, and a
        RATIO of cross sections at two energies is completely insensitive to the constant in front.

        THAT TEST HAS BEEN RUN, on beryllium-like neon at 1, 2 and 4 eV (example-Pb.jl branch a), and it SPLITS THE GAUGES. The
        Coulomb gauge reproduces the law with exponents 4.014 and 4.056 -- slightly above 4 and rising toward the 2s-2p resonance,
        which is what it should do -- so the amplitude, the sum over intermediate levels, the two time orderings and their
        denominators are all doing their job. The Babushkin gauge comes out exactly four powers steeper, 8.014 and 8.056, and the
        Coulomb/Babushkin ratio falls by 16.0 per doubling at BOTH steps, i.e. by 2^4 to three significant figures.

        An EXACT power localizes the fault, and it is not in this prefactor -- a constant cannot change a power law. It was then
        MEASURED at a single vertex, which is where it lives. For the same Be-like neon levels, |<nu||O(E1,omega)||i>| over
        omega = 0.5, 1, 2, 4 eV:

            Coulomb     7.513941e-07  7.514237e-07  7.515424e-07  7.520171e-07     flat to four digits
            Babushkin   1.825198e-07  3.650593e-07  7.301977e-07  1.460712e-06     exactly proportional to omega

        the ratio falling by 0.5000, 0.5000, 0.5002 per doubling. Two such vertices give omega^2 in the amplitude and omega^4 in
        the cross section, which is the whole of the discrepancy. The structural reason is visible in
        InteractionStrength.multipoleTransition: the length form carries j_L(qr) while the velocity form carries j_L(qr)/(qr) and
        j'_L(qr), and for small argument j_1(x) ~ x/3 against j_1(x)/x ~ 1/3.

        THE QUESTION THIS RAISED IS ANSWERED, and the answer puts the fault here rather than in PhotoEmission. Length-velocity
        equivalence is an ON-SHELL IDENTITY: the two forms are related through the commutator [H, r] = -i p / m, and turning
        <f|p|i> into omega <f|r|i> uses (E_f - E_i) = omega. OFF SHELL that step is false, so the two amplitudes have no reason to
        agree, and the Bessel asymmetry above is exactly how the discrepancy appears -- the length form picks up the power of omega
        that the commutator identity would have supplied. Nothing is wrong in InteractionStrength or PhotoEmission; this file uses
        a real-photon amplitude at a FREE omega and must supply the weighting itself. Two on-shell facts settle the direction:
        helium 1s^2 -> 1s2p at its true transition energy gives f_Coulomb = 0.378 against f_Babushkin = 0.421, a ratio of 0.897
        and nowhere near omega or 1/omega; and the GeneralizedOscillatorStrength module, which computes the LENGTH form by
        construction, reproduces PhotoExcitation's BABUSHKIN oscillator strength to 1.000035 in the K -> 0 limit.

        SO THERE ARE TWO INDEPENDENT FAULTS, and they were conflated until measured apart:

          (1) THE OFF-SHELL WEIGHTING, above: structural, understood, and belonging to this file.
          (2) BASIS QUALITY, exposed by measuring the gauge ratio AT the on-shell points of the test system rather than assuming
              it would be unity there:

                  1^- #1   dE = 10.816 eV   |A_Cou| = 7.560e-07   |A_Bab| = 3.953e-06   ratio 0.191
                  1^- #2   dE = 25.003 eV   |A_Cou| = 1.059e-03   |A_Bab| = 1.434e-03   ratio 0.738

              Neither is near unity. The initial and intermediate multiplets come from two SEPARATELY CONVERGED SCF runs, so the
              orbitals are non-orthogonal and no biorthogonal transformation is applied -- which is what PhotoEmission's
              calcBiorthogonal exists for. The two levels also differ by three orders of magnitude in amplitude, and it is the WEAK
              one whose gauges disagree worst: #1 is an intercombination line whose E1 amplitude survives only through
              singlet-triplet mixing, i.e. a small difference of large numbers. Strong cancellation is where gauge agreement fails
              first, which is why JAC computes a cancellation factor at all.

        BOTH ARE NOW RESOLVED, and neither the way it first looked.

        Fault (1) is FIXED by PhotonScattering.offShellFactor, and confirmed four ways: both gauges now give the omega^4 exponents
        4.014 and 4.056 identically; the Coulomb numbers are BIT-IDENTICAL to the pre-fix run, so only Babushkin moved; the gauge
        ratio is constant to five digits at 3.3634 where it previously varied by 256 across the scan; and that residual falls
        BETWEEN the two independently measured on-shell ratios, 5.228 and 1.354, which is where a weighted combination of the two
        intermediate levels belongs. The fourth was not designed in and came from a calculation made for another purpose.

        Fault (2) is EXONERATED. A three-point comparison -- separate DFS 0.191/0.738, common mean field 0.168/0.726, common
        nuclear field 0.0155/0.0293 -- shows that sharing the one-body Hamiltonian changes essentially nothing while destroying
        state quality changes everything. Both DFS runs used scField = DFSField(), checked rather than assumed. The on-shell
        disagreement is therefore CORRELATION, not the non-orthogonality it was blamed on, and no biorthogonal transformation is
        called for. Note also that Basics.NuclearField() is a confounded diagnostic: it shares the Hamiltonian but removes all
        screening, varying two things at once, which is why it made the ratio ten times worse instead of better.

        WHAT REMAINS. The SHAPE is verified in both gauges; the ABSOLUTE MAGNITUDE is not, this function's prefactor never having
        been derived. And nothing tested here constrains a PHASE -- a sign error in one time ordering would survive every check
        above. The Raman Hermiticity test is the missing instrument.
"""
function rayleighCrossSection(channels::Array{PhotonScattering.Channel,1}, inOmega::Float64, outOmega::Float64,
                              initialLevel::Level)
    wc = Defaults.getDefaults("speed of light: c")
    wCo = 0.;   wBa = 0.
    for  gauge in [Basics.Coulomb, Basics.Babushkin, Basics.Magnetic]
        symmetries = LevelSymmetry[]
        for  ch in channels    if  ch.gauge == gauge  &&  !(ch.totalSymmetry in symmetries)   push!(symmetries, ch.totalSymmetry)   end   end
        for  symt in symmetries
            wa = ComplexF64(0.)
            for  ch in channels
                if  ch.gauge == gauge  &&  ch.totalSymmetry == symt    wa = wa + ch.amplitude    end
            end
            # A magnetic multipole has no gauge freedom and therefore contributes to both gauge sums alike
            if      gauge == Basics.Magnetic     wCo = wCo + abs(wa)^2;    wBa = wBa + abs(wa)^2
            elseif  gauge == Basics.Coulomb      wCo = wCo + abs(wa)^2
            else                                 wBa = wBa + abs(wa)^2
            end
        end
    end
    factor = 8pi / 3 * inOmega * outOmega^3 / wc^4 / (Basics.twice(initialLevel.J) + 1)

    return( EmProperty(factor * wCo, factor * wBa) )
end
