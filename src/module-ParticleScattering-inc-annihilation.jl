
# One-photon annihilation of the projectile POSITRON with a BOUND electron:   e^+ + |i(N)>  -->  |f(N-1)> + photon.
#
# WHY THE PROCESS EXISTS AT ALL. A free electron and a positron cannot annihilate into a single photon: in their centre-of-mass
# frame the pair is at rest and the photon would have to carry zero momentum and 2 m c^2 of energy at once. A BOUND electron
# lifts that, because the nucleus takes up the recoil, and the amplitude is therefore controlled by how much of the electron's
# momentum distribution reaches p ~ m c -- which lives close to the nucleus. That is the whole physics of the process: it is a
# K-shell, high-Z process, its cross section rises roughly as Z^5 for K-shell annihilation, and it vanishes in the free limit.
#
# THE TRICK THAT MAKES IT COMPUTABLE HERE. By charge conjugation, annihilating a bound electron against an incoming positron is
# the same matrix element as the bound electron making a transition into a NEGATIVE-energy continuum state, the charge-conjugate
# of the positron wave. The operator is then literally the photon-emission operator, so this file adds no new vertex: it builds
# the charge-conjugate orbital and hands it to PhotoEmission.amplitude, exactly as InternalConversion hands its ejected-electron
# orbital to a rank-L one-particle operator. Structurally this is AutoIonization's pipeline with a photon on the outgoing side.
#
# ENERGY CONSERVATION fixes the photon energy without any freedom. Writing rest masses explicitly, the N-electron atom plus a
# positron of kinetic energy T_+ goes to an (N-1)-electron ion plus one photon, so
#
#         omega  =  2 m c^2  +  T_+  -  I ,          I = E_f - E_i  the binding energy of the annihilated electron.
#
# THE COST, and it is the reason this file carries a warning rather than a recommendation. omega ~ 2 m c^2 = 1.022 MeV gives a
# photon momentum q a_0 = omega / 3729 eV ~ 274 in atomic units, so the retardation factor exp(i q.r) oscillates ~ q <r> times
# across the orbital. For a K-shell electron <r> ~ a_0/Z, hence q <r> ~ 274/Z: about 3 for uranium, but about 137 for helium.
# The multipole series is therefore SHORT only for heavy elements and hopeless for light ones -- the opposite of the intuition
# carried over from optical transitions, and the reason a light-element test case is the wrong first target here.
#
#
# TWO-PHOTON ANNIHILATION IS NOT HERE, AND DELIBERATELY SO, decided 21-Aug-2026 by the maintainer. Do not add it to this file
# without asking him; it is a postponement, not an oversight.
#
# e^+ + |i(N)> --> |f(N-1)> + photon + photon is the dominant channel in nature -- it needs no nucleus, so it does not carry the
# Z^5 suppression of the one-photon case -- but it is SECOND order, two vertices with a sum over intermediate states, and JAC has
# no machinery for that today. The point that decided the placement is that this is the SAME missing machinery as two-photon
# (2E1) DECAY, whose absence is already a recorded defect: example-Je.jl branch a shows 1s2s ^1S_0 of helium-like carbon getting
# a 513 s lifetime from its M1 channel alone, against a true ~2.7e-5 s, i.e. wrong by ~2e7, and notes that at low density this
# would dominate. Whoever builds the second-order two-photon machinery therefore enables BOTH, and putting the annihilation half
# of it inside a particle-scattering file would hide that connection from the person who comes to fix the decay half.
#
# The relevant existing ingredient is AtomicState.GreenExpansion, the sum-over-states machinery that MultipolePolarizibility
# uses; note that that module is itself paused on a B-spline question, so the dependency is not free ground.
#
# What was retired to get here: src/module-PairAnnihilation1Photon.jl and src/module-PairAnnihilation2Photon.jl, two stubs that
# had never been compiled -- they sat inside the "#= Further processes, not yet included =#" block of JenaAtomicCalculator.jl --
# and that contained no physics at all (amplitude = 1.0, an undefined `phase`, the module used as a type). Commit 6ba9985.


"""
`struct  ParticleScattering.AnnihilationChannel`
    ... defines a type for one reduced amplitude of the one-photon annihilation of the projectile positron with a bound electron,
        i.e. for one partial wave of the incoming positron, one photon multipole and one gauge.

    + kappa           ::Int64            ... kappa of the incoming positron partial wave.
    + multipole       ::EmMultipole      ... multipole of the emitted photon.
    + gauge           ::EmGauge          ... gauge in which this amplitude was evaluated.
    + totalSymmetry   ::LevelSymmetry    ... J^pi of the coupled (ion + positron) system; amplitudes sharing it interfere.
    + amplitude       ::ComplexF64       ... the reduced amplitude itself.
"""
struct  AnnihilationChannel
    kappa             ::Int64
    multipole         ::EmMultipole
    gauge             ::EmGauge
    totalSymmetry     ::LevelSymmetry
    amplitude         ::ComplexF64
end


"""
`ParticleScattering.AnnihilationChannel()`  ... constructor for a default ParticleScattering.AnnihilationChannel.
"""
function AnnihilationChannel()
    AnnihilationChannel(-1, E1, Basics.Coulomb, LevelSymmetry(AngularJ64(0), Basics.plus), ComplexF64(0.))
end


# `Base.show(io::IO, ch::ParticleScattering.AnnihilationChannel)`  ... prepares a printout of ch::AnnihilationChannel.
function Base.show(io::IO, ch::ParticleScattering.AnnihilationChannel)
    println(io, "kappa = $(ch.kappa),  multipole = $(ch.multipole),  gauge = $(ch.gauge),  J^pi = $(ch.totalSymmetry),  " *
                "amplitude = $(ch.amplitude)")
end


"""
`struct  ParticleScattering.AnnihilationLine`
    ... defines a type for the one-photon annihilation of a positron of one impact energy with the bound electrons of one initial
        level, leaving one final level of the ion with one electron less.

    + initialLevel    ::Level            ... initial level of the N-electron target.
    + finalLevel      ::Level            ... final level of the (N-1)-electron ion.
    + positronEnergy  ::Float64          ... kinetic energy of the incoming positron [a.u.].
    + photonEnergy    ::Float64          ... energy of the emitted photon [a.u.], fixed by energy conservation, see the file header.
    + crossSection    ::EmProperty       ... total annihilation cross section, in both gauges.
    + channels        ::Array{ParticleScattering.AnnihilationChannel,1}   ... the reduced amplitudes; the primary result.
"""
struct  AnnihilationLine
    initialLevel      ::Level
    finalLevel        ::Level
    positronEnergy    ::Float64
    photonEnergy      ::Float64
    crossSection      ::EmProperty
    channels          ::Array{ParticleScattering.AnnihilationChannel,1}
end


"""
`ParticleScattering.AnnihilationLine()`  ... constructor for a default ParticleScattering.AnnihilationLine.
"""
function AnnihilationLine()
    AnnihilationLine( Level(), Level(), 0., 0., EmProperty(0., 0.), AnnihilationChannel[] )
end


# `Base.show(io::IO, line::ParticleScattering.AnnihilationLine)`  ... prepares a printout of line::AnnihilationLine.
function Base.show(io::IO, line::ParticleScattering.AnnihilationLine)
    println(io, "initialLevel:       $(line.initialLevel.index)  ")
    println(io, "finalLevel:         $(line.finalLevel.index)  ")
    println(io, "positronEnergy:     $(line.positronEnergy)  ")
    println(io, "photonEnergy:       $(line.photonEnergy)  ")
    println(io, "crossSection:       $(line.crossSection)  ")
    println(io, "channels:           $(line.channels)  ")
end


"""
`ParticleScattering.annihilationPhotonEnergy(finalLevel::Level, initialLevel::Level, positronEnergy::Float64)`
    ... to compute the energy of the single annihilation photon from energy conservation, omega = 2 m c^2 + T_+ - (E_f - E_i), where
        the last bracket is the binding energy of the annihilated electron. Nothing here is adjustable; the expression is the whole
        content of energy conservation for this process. An omega::Float64 [a.u.] is returned, which is negative -- and the channel
        therefore closed -- only for a binding energy in excess of 2 m c^2 plus the impact energy.
"""
function annihilationPhotonEnergy(finalLevel::Level, initialLevel::Level, positronEnergy::Float64)
    wc = Defaults.getDefaults("speed of light: c")

    return( 2 * wc * wc + positronEnergy - (finalLevel.energy - initialLevel.energy) )
end


"""
`ParticleScattering.chargeConjugateOrbital(orbital::Radial.Orbital)`
    ... to form the charge-conjugate of a positron continuum orbital, i.e. the negative-energy electron solution that represents the
        same physical positron. Charge conjugation C psi = i gamma^2 psi^* reverses the sign of kappa and exchanges the large and the
        small radial component,

            kappa -> -kappa ,      P(r) -> -Q(r) ,      Q(r) ->  P(r) ,

        which is what lets the annihilation matrix element be evaluated with the ordinary photon-emission operator. A
        cOrbital::Radial.Orbital is returned, carrying the conjugated subshell.

        THE OVERALL SIGN IS A CONVENTION and cannot be settled from within this file: it cancels in |amplitude|^2 for a single
        channel, but it does NOT cancel between two multipoles feeding the same total symmetry, where the amplitudes interfere.
        Any disagreement with measured multipole ratios should be checked here first.
"""
function chargeConjugateOrbital(orbital::Radial.Orbital)
    sh = Subshell(orbital.subshell.n, -orbital.subshell.kappa)

    return( Radial.Orbital( sh, orbital.isBound, orbital.useStandardGrid, orbital.energy, -orbital.Q, orbital.P,
                            orbital.Qprime, orbital.Pprime, orbital.grid ) )
end


"""
`ParticleScattering.computeAnnihilationAmplitudesProperties(line::ParticleScattering.AnnihilationLine, nm::Nuclear.Model,
                                                            grid::Radial.Grid, nrContinuum::Int64,
                                                            settings::ParticleScattering.Settings;
                                                            nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)`
    ... to compute all channel amplitudes of one annihilation line and, from them, its total cross section. For every channel the
        positron partial wave is generated in the POSITRON potential -- the electrostatic interaction with reversed sign and no
        exchange, cf. ParticleScattering.scatteringPotential -- charge-conjugated, attached to the final level as an extra electron,
        and handed to PhotoEmission.amplitude, the operator being the ordinary photon-emission operator. A
        newLine::ParticleScattering.AnnihilationLine is returned with its channels and cross section evaluated.
"""
function computeAnnihilationAmplitudesProperties(line::ParticleScattering.AnnihilationLine, nm::Nuclear.Model,
                                                 grid::Radial.Grid, nrContinuum::Int64, settings::ParticleScattering.Settings;
                                                 nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)
    contSettings = Continuum.Settings(false, nrContinuum)
    nucPot       = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    pot          = ParticleScattering.scatteringPotential(ParticleScattering.Positron(), settings.interaction, line.finalLevel,
                                                          line.positronEnergy, nm, grid; nuclearPot=nucPot)
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    redILevel    = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    newfLevel    = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
    newChannels  = ParticleScattering.AnnihilationChannel[]

    for  channel in ParticleScattering.determineAnnihilationChannels(line.finalLevel, line.initialLevel, settings)
        # channel.kappa is the kappa of the CONJUGATED orbital, i.e. of the one that enters the CSF and the angular coupling; the
        # POSITRON partial wave that has to be generated is therefore the one with the opposite sign, since charge conjugation
        # reverses kappa. Getting this the wrong way round costs no error message: the conjugated orbital then cannot couple to
        # channel.totalSymmetry, no CSF survives, and every amplitude comes out exactly zero.
        pKappa             = -channel.kappa
        pOrbital, phase, _ = Continuum.generateOrbitalLocalPotential(line.positronEnergy, Subshell(101, pKappa), pot,
                                                                     contSettings)
        cOrbital  = ParticleScattering.chargeConjugateOrbital(pOrbital)
        newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, channel.totalSymmetry, newfLevel)
        amplitude = PhotoEmission.amplitude(Basics.Emission(), channel.multipole, channel.gauge, line.photonEnergy,
                                            newcLevel, redILevel, grid; printout=printout)
        # The incoming-wave boundary condition of the positron; i^l exp(i delta) is the convention of PhotoRecombination, NOT the
        # i^l exp(-i delta) of AutoIonization, and the two differ by (-1)^l -- see AutoIonization.captureAmplitude. The l here is
        # the POSITRON's, not the conjugated orbital's.
        l         = Basics.subshell_l( Subshell(101, pKappa) )
        amplitude = amplitude * (1im)^l * exp(1im * phase)
        push!( newChannels, ParticleScattering.AnnihilationChannel(channel.kappa, channel.multipole, channel.gauge,
                                                                   channel.totalSymmetry, amplitude) )
    end

    crossSection = ParticleScattering.annihilationCrossSection(newChannels, line.positronEnergy, line.photonEnergy,
                                                               line.initialLevel)
    newLine      = ParticleScattering.AnnihilationLine( line.initialLevel, line.finalLevel, line.positronEnergy,
                                                        line.photonEnergy, crossSection, newChannels )

    return( newLine )
end


"""
`ParticleScattering.annihilationCrossSection(channels::Array{ParticleScattering.AnnihilationChannel,1}, positronEnergy::Float64,
                                             photonEnergy::Float64, initialLevel::Level)`
    ... to form the total one-photon annihilation cross section from the evaluated channel amplitudes. Amplitudes that feed the SAME
        total symmetry J^pi interfere and are therefore summed coherently first; the squared moduli are then summed over J^pi, over
        the gauges being kept apart, and divided by the statistical weight of the initial level and by the positron flux,

            sigma  =  8 pi^3 alpha omega / (k_+^2 (2 J_i + 1))  *  SUM_{J^pi} | SUM_{kappa, multipole} A |^2 ,

        with k_+ the relativistic positron momentum. A crossSection::EmProperty [a.u.] is returned.

        THE PREFACTOR IS THE ONE PIECE OF THIS FILE NOT INDEPENDENTLY CHECKED. The channel structure, the coherent sum within a total
        symmetry and the k_+^-2 flux factor are standard, but the numerical constant multiplying them has been written down rather
        than derived against a known case. Any RATIO of cross sections -- the Z-dependence above all -- is insensitive to it, which is
        why the isoelectronic scan is the meaningful first test and an absolute comparison is not.
"""
function annihilationCrossSection(channels::Array{ParticleScattering.AnnihilationChannel,1}, positronEnergy::Float64,
                                  photonEnergy::Float64, initialLevel::Level)
    wc    = Defaults.getDefaults("speed of light: c");    alpha = 1.0 / wc
    kPlus = ParticleScattering.projectileMomentum(positronEnergy)
    if  kPlus < 1.0e-12    return( EmProperty(0., 0.) )    end
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
    factor = 8 * pi^3 * alpha * photonEnergy / ( kPlus * kPlus * (Basics.twice(initialLevel.J) + 1) )

    return( EmProperty(factor * wCo, factor * wBa) )
end


"""
`ParticleScattering.computeAnnihilationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                             grid::Radial.Grid, settings::ParticleScattering.Settings; output=true)`
    ... to compute all selected one-photon annihilation lines, together with their channel amplitudes and cross sections, and to
        display them. An Array{ParticleScattering.AnnihilationLine,1} is returned if output=true, and nothing otherwise.
"""
function computeAnnihilationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                  settings::ParticleScattering.Settings; output=true)
    println("")
    printstyled("ParticleScattering.computeAnnihilationLines(): The computation of annihilation amplitudes starts now ... \n",
                color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    lines = ParticleScattering.determineAnnihilationLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    ParticleScattering.displayAnnihilationLines(stdout, lines)    end
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.positronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    nucPot      = Nuclear.nuclearPotential(nm, grid)
    newLines    = ParticleScattering.AnnihilationLine[]
    for  line in lines
        push!( newLines, ParticleScattering.computeAnnihilationAmplitudesProperties(line, nm, grid, nrContinuum, settings;
                                                                                    nuclearPot=nucPot, printout=false) )
    end
    ParticleScattering.displayAnnihilationResults(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ParticleScattering.displayAnnihilationResults(iostream, newLines)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`ParticleScattering.determineAnnihilationChannels(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)`
    ... to determine the allowed (kappa, multipole, gauge) channels of one annihilation line. The photon multipole couples the
        initial symmetry to an intermediate one, and the incoming positron's kappa then couples that to the final symmetry -- the
        same two-step construction as InternalConversion.determineChannels, a rank-L one-particle operator requiring it. Channels
        with |kappa| beyond settings.maxL + 1 are dropped. An Array{ParticleScattering.AnnihilationChannel,1} is returned with the
        amplitudes not yet evaluated.
"""
function determineAnnihilationChannels(finalLevel::Level, initialLevel::Level, settings::ParticleScattering.Settings)
    channels = ParticleScattering.AnnihilationChannel[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    for  mp in settings.multipoles
        for  symt in AngularMomentum.allowedMultipoleSymmetries(symi, mp)
            for  kappa in AngularMomentum.allowedKappaSymmetries(symt, symf)
                if  abs(kappa) > settings.maxL + 1    continue    end
                for  gauge in settings.gauges
                    if      mp.electric  &&  gauge == Basics.UseCoulomb
                        push!( channels, ParticleScattering.AnnihilationChannel(kappa, mp, Basics.Coulomb,   symt, ComplexF64(0.)) )
                    elseif  mp.electric  &&  gauge == Basics.UseBabushkin
                        push!( channels, ParticleScattering.AnnihilationChannel(kappa, mp, Basics.Babushkin, symt, ComplexF64(0.)) )
                    elseif  !mp.electric &&  gauge == settings.gauges[1]
                        push!( channels, ParticleScattering.AnnihilationChannel(kappa, mp, Basics.Magnetic,  symt, ComplexF64(0.)) )
                    end
                end
            end
        end
    end

    return( channels )
end


"""
`ParticleScattering.determineAnnihilationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                               settings::ParticleScattering.Settings)`
    ... to determine the list of annihilation lines to be computed, i.e. one line per selected pair of levels and per impact energy
        of settings.impactEnergies. A line whose photon energy comes out non-positive is dropped, that channel being closed. An
        Array{ParticleScattering.AnnihilationLine,1} is returned with the amplitudes not yet evaluated.
"""
function determineAnnihilationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                    settings::ParticleScattering.Settings)
    lines = ParticleScattering.AnnihilationLine[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            for  energy in settings.impactEnergies
                omega = ParticleScattering.annihilationPhotonEnergy(fLevel, iLevel, energy)
                if  omega <= 0.    continue    end
                push!( lines, ParticleScattering.AnnihilationLine(iLevel, fLevel, energy, omega, EmProperty(0., 0.),
                                                                   ParticleScattering.AnnihilationChannel[]) )
            end
        end
    end

    return( lines )
end


"""
`ParticleScattering.displayAnnihilationLines(stream::IO, lines::Array{ParticleScattering.AnnihilationLine,1})`
    ... to list the selected annihilation lines before their amplitudes are computed, so that a long run can be checked against what
        was intended. A neat table is printed but nothing is returned otherwise.
"""
function displayAnnihilationLines(stream::IO, lines::Array{ParticleScattering.AnnihilationLine,1})
    nx = 108
    println(stream, " ")
    println(stream, "  Selected one-photon annihilation lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(14, "Energy e+"; na=4) * TableStrings.center(14, "Energy photon"; na=4) *
                TableStrings.center(20, "No. channels"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.positronEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa  = sa * string( length(line.channels) )
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`ParticleScattering.displayAnnihilationResults(stream::IO, lines::Array{ParticleScattering.AnnihilationLine,1})`
    ... to print the computed one-photon annihilation cross sections, one row per line and with both gauges beside each other so that
        their agreement can be read off directly. A neat table is printed but nothing is returned otherwise.
"""
function displayAnnihilationResults(stream::IO, lines::Array{ParticleScattering.AnnihilationLine,1})
    nx = 122
    println(stream, " ")
    println(stream, "  One-photon annihilation cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(14, "Energy e+"; na=4) * TableStrings.center(14, "Energy photon"; na=4) *
                TableStrings.center(30, "Cross section  (Coulomb, Babushkin)"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.positronEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "  " * @sprintf("%.6e", line.crossSection.Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end
