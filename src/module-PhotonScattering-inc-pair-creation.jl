
# Bound-free pair creation:   gamma + |i(N)>  -->  |f(N+1)> + e^+ .
#
# A photon is absorbed and an electron-positron pair appears, the ELECTRON being captured directly into a bound orbital of the ion --
# so the ion gains an electron -- while the POSITRON leaves as a continuum particle. Free pair creation in vacuum is impossible for a
# single photon, for the same reason one-photon annihilation is: momentum cannot balance. Here the nucleus takes the recoil.
#
# THE THRESHOLD IS LOWERED BY THE BINDING ENERGY, and that is the cleanest qualitative signature of the process:
#
#         T_+ = omega - 2 m c^2 + B ,        so the channel opens at    omega > 2 m c^2 - B ,
#
# with B the binding energy gained by the captured electron. For a K-shell capture in a heavy ion B is of order 100 keV, so bound-free
# pair creation sets in appreciably below the 1.022 MeV of the free process. It is not an exotic correction: at the LHC it is the
# dominant beam-loss mechanism for fully stripped heavy ions, the created electron being captured by the very ion that made it.
#
# THIS IS THE CROSSING PARTNER of the one-photon annihilation in module-ParticleScattering-inc-annihilation.jl, and the two share
# their machinery deliberately: the same charge-conjugate positron wave (ParticleScattering.chargeConjugateOrbital), the same photon
# operator, run in the other direction -- absorption here, emission there. That is also the validation route this file was chosen for.
# Detailed balance relates the two cross sections at the same total energy, so the annihilation module's unverified prefactor and this
# one's can be checked against EACH OTHER without any literature at all; two prefactors that fail to satisfy detailed balance locate
# the error between them, which no single absolute number can do.
#
# The retardation caveat of the annihilation file applies here unchanged: omega > 2 m c^2 means q a_0 > 274, so the multipole series
# is short only where the captured electron sits close to the nucleus, i.e. for K-shell capture in a heavy ion.


"""
`PhotonScattering.computePairCreationAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                          nrContinuum::Int64, settings::PhotonScattering.Settings;
                                                          nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)`
    ... to compute all channel amplitudes of one bound-free pair-creation line and, from them, its total cross section. For every
        channel the positron partial wave is generated in the positron potential of the FINAL ion, charge-conjugated, attached to the
        INITIAL level as an extra electron so that both sides carry N+1 electrons, and handed to PhotoEmission.amplitude in its
        absorption form. A newLine::PhotonScattering.Line is returned with its channels and cross section evaluated.
"""
function computePairCreationAmplitudesProperties(line::PhotonScattering.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                 nrContinuum::Int64, settings::PhotonScattering.Settings;
                                                 nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)
    contSettings = Continuum.Settings(false, nrContinuum)
    nucPot       = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    pot          = ParticleScattering.scatteringPotential(ParticleScattering.Positron(), ParticleScattering.StaticField(),
                                                          line.finalLevel, line.particleEnergy, nm, grid; nuclearPot=nucPot)
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    redILevel    = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    redFLevel    = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
    newChannels  = PhotonScattering.Channel[]

    for  channel in PhotonScattering.determinePairCreationChannels(line.finalLevel, line.initialLevel, settings)
        # channel.kappa is the kappa of the CONJUGATED orbital, the one entering the CSF and the angular coupling; the POSITRON
        # partial wave carries the opposite sign, charge conjugation reversing kappa. Getting this backwards is silent -- no CSF
        # then survives the coupling and every amplitude is exactly zero. It cost a debugging round in the annihilation module.
        pKappa             = -channel.kappa
        pOrbital, phase, _ = Continuum.generateOrbitalLocalPotential(line.particleEnergy, Subshell(101, pKappa), pot, contSettings)
        cOrbital  = ParticleScattering.chargeConjugateOrbital(pOrbital)
        newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, channel.totalSymmetry, redILevel)
        amplitude = PhotoEmission.amplitude(Basics.Absorption(), channel.inMultipole, channel.gauge, line.inPhotonEnergy,
                                            redFLevel, newcLevel, grid; printout=printout)
        l         = Basics.subshell_l( Subshell(101, pKappa) )
        amplitude = amplitude * (1im)^l * exp(1im * phase)
        push!( newChannels, PhotonScattering.Channel(channel.kappa, channel.inMultipole, channel.outMultipole, channel.gauge,
                                                      channel.totalSymmetry, amplitude) )
    end

    crossSection = PhotonScattering.pairCreationCrossSection(newChannels, line.inPhotonEnergy, line.particleEnergy, line.initialLevel)
    newLine      = PhotonScattering.Line( line.initialLevel, line.finalLevel, line.inPhotonEnergy, line.outPhotonEnergy,
                                          line.particleEnergy, crossSection, newChannels, PhotonScattering.Observables[] )

    return( newLine )
end


"""
`PhotonScattering.computePairCreationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                           grid::Radial.Grid, settings::PhotonScattering.Settings; output=true)`
    ... to compute all selected bound-free pair-creation lines, together with their channel amplitudes and cross sections, and to
        display them. An Array{PhotonScattering.Line,1} is returned if output=true, and nothing otherwise.
"""
function computePairCreationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                  settings::PhotonScattering.Settings; output=true)
    println("")
    printstyled("PhotonScattering.computePairCreationLines(): The computation of pair-creation amplitudes starts now ... \n",
                color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    lines = PhotonScattering.determinePairCreationLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotonScattering.displayPairCreationLines(stdout, lines)    end
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.particleEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    nucPot      = Nuclear.nuclearPotential(nm, grid)
    newLines    = PhotonScattering.Line[]
    for  line in lines
        push!( newLines, PhotonScattering.computePairCreationAmplitudesProperties(line, nm, grid, nrContinuum, settings;
                                                                                   nuclearPot=nucPot, printout=false) )
    end
    PhotonScattering.displayPairCreationResults(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotonScattering.displayPairCreationResults(iostream, newLines)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotonScattering.determinePairCreationChannels(finalLevel::Level, initialLevel::Level, settings::PhotonScattering.Settings)`
    ... to determine the allowed (kappa, multipole, gauge) channels of one pair-creation line. The two-step construction runs the
        OPPOSITE way round from the annihilation case, and that is not a detail: the charge-conjugate positron orbital is attached to
        the INITIAL level here, so the intermediate symmetry symt is the one carried by (initial level + positron), the absorbed
        photon connects symt to the FINAL symmetry, and the positron's kappa must couple symi to symt. Taking the annihilation
        ordering instead -- symt from symi, kappa to symf -- costs no error message and yields amplitudes that are all exactly zero.
        An Array{PhotonScattering.Channel,1} is returned with the amplitudes not yet evaluated.
"""
function determinePairCreationChannels(finalLevel::Level, initialLevel::Level, settings::PhotonScattering.Settings)
    channels = PhotonScattering.Channel[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    for  mp in settings.multipoles
        for  symt in AngularMomentum.allowedMultipoleSymmetries(symf, mp)
            for  kappa in AngularMomentum.allowedKappaSymmetries(symt, symi)
                if  abs(kappa) > settings.maxKappa    continue    end
                for  gauge in settings.gauges
                    if      mp.electric  &&  gauge == Basics.UseCoulomb
                        push!( channels, PhotonScattering.Channel(kappa, mp, mp, Basics.Coulomb,   symt, ComplexF64(0.)) )
                    elseif  mp.electric  &&  gauge == Basics.UseBabushkin
                        push!( channels, PhotonScattering.Channel(kappa, mp, mp, Basics.Babushkin, symt, ComplexF64(0.)) )
                    elseif  !mp.electric &&  gauge == settings.gauges[1]
                        push!( channels, PhotonScattering.Channel(kappa, mp, mp, Basics.Magnetic,  symt, ComplexF64(0.)) )
                    end
                end
            end
        end
    end

    return( channels )
end


"""
`PhotonScattering.determinePairCreationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                             settings::PhotonScattering.Settings)`
    ... to determine the list of pair-creation lines to be computed, i.e. one line per selected pair of levels and per photon energy
        of settings.photonEnergies. A line whose positron energy comes out non-positive is dropped, the channel being below threshold
        -- which for this process means below 2 m c^2 minus the binding energy gained, not below 2 m c^2. An
        Array{PhotonScattering.Line,1} is returned with the amplitudes not yet evaluated.
"""
function determinePairCreationLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotonScattering.Settings)
    lines = PhotonScattering.Line[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            for  omega in settings.photonEnergies
                energy = PhotonScattering.positronEnergy(fLevel, iLevel, omega)
                if  energy <= 0.    continue    end
                push!( lines, PhotonScattering.Line(iLevel, fLevel, omega, 0., energy, EmProperty(0., 0.),
                                                     PhotonScattering.Channel[], PhotonScattering.Observables[]) )
            end
        end
    end

    return( lines )
end


"""
`PhotonScattering.displayPairCreationLines(stream::IO, lines::Array{PhotonScattering.Line,1})`
    ... to list the selected pair-creation lines before their amplitudes are computed, so that a long run can be checked against what
        was intended. A neat table is printed but nothing is returned otherwise.
"""
function displayPairCreationLines(stream::IO, lines::Array{PhotonScattering.Line,1})
    nx = 108
    println(stream, " ")
    println(stream, "  Selected bound-free pair-creation lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(14, "Energy photon"; na=4) * TableStrings.center(14, "Energy e+"; na=4) *
                TableStrings.center(20, "No. channels"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.particleEnergy)) * "    "
        sa  = sa * string( length(line.channels) )
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotonScattering.displayPairCreationResults(stream::IO, lines::Array{PhotonScattering.Line,1})`
    ... to print the computed bound-free pair-creation cross sections, one row per line and with both gauges beside each other so that
        their agreement can be read off directly. A neat table is printed but nothing is returned otherwise.
"""
function displayPairCreationResults(stream::IO, lines::Array{PhotonScattering.Line,1})
    nx = 122
    println(stream, " ")
    println(stream, "  Bound-free pair-creation cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(18, "i-level-f"; na=2) * TableStrings.center(18, "i--J^P--f"; na=4) *
                TableStrings.center(14, "Energy photon"; na=4) * TableStrings.center(14, "Energy e+"; na=4) *
                TableStrings.center(30, "Cross section  (Coulomb, Babushkin)"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  " * TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index) * "    "
        sa  = sa * TableStrings.symmetries_if(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity),
                                              LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.particleEnergy)) * "    "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "  " * @sprintf("%.6e", line.crossSection.Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotonScattering.pairCreationCrossSection(channels::Array{PhotonScattering.Channel,1}, inPhotonEnergy::Float64,
                                           particleEnergy::Float64, initialLevel::Level)`
    ... to form the total bound-free pair-creation cross section from the evaluated channel amplitudes. Amplitudes feeding the SAME
        total symmetry interfere and are summed coherently first; the squared moduli are then summed over the symmetries, the gauges
        being kept apart, and divided by the statistical weight of the initial level and by the incoming photon flux,

            sigma  =  8 pi^3 alpha p_+ / (omega (2 J_i + 1))  *  SUM_{J^pi} | SUM_{kappa, multipole} A |^2 ,

        with p_+ the relativistic positron momentum. A crossSection::EmProperty [a.u.] is returned.

        THE PREFACTOR IS NOT INDEPENDENTLY VERIFIED, exactly as in ParticleScattering.annihilationCrossSection, and for the same
        reason: it has been written down rather than derived against a known case. The two have been checked TOGETHER, which is what
        this process was implemented first for, and THEY ARE NOT MUTUALLY CONSISTENT.

        DETAILED BALANCE AGAINST THE CROSSING PARTNER FAILS, and the failure is not yet reduced to one cause. Running the
        time-reverse -- helium-like against hydrogen-like uranium at omega = 894.556 keV -- detailed balance requires
        sigma_ann / sigma_pc = 573.97 and JAC gives 3.1327e+05 (Coulomb) and 3.2065e+05 (Babushkin), i.e. violation factors of
        545.8 and 558.6.

        WHAT THAT DOES AND DOES NOT SHOW. It was first read here as a constant of about 4c = 548.14; that reading is WITHDRAWN,
        having rested on a single energy point. Comparing the two prefactors directly refutes it: were the amplitudes equal in
        both directions, the coded prefactors alone would give a violation of c^2/p_+ = 1326 at this point -- energy-DEPENDENT,
        going as 1/p_+, and not 548. The observed 546-559 is 0.41 of that, so the amplitudes are not equal in the two directions
        either, and BOTH the prefactor and the amplitude are implicated.

        THE SECOND-ENERGY TEST HAS NOW BEEN RUN, and it rules out both simple explanations. Doubling p_+ (T_+ = 100 -> 400 a.u.,
        p_+ = 14.1610 -> 28.4345) moves the violation from 545.8/558.6 to 448.3/462.0, i.e. by 0.821 and 0.827. A pure prefactor
        error would have given 0.498, a constant one 1.000; it is neither, the implied exponent being about p_+^(-0.28).

        The useful part is that the energy dependence is the SAME in both gauges to within half a percent, hence GAUGE-INDEPENDENT,
        while the gauge disagreement itself is nearly energy-independent (17.5 -> 17.2 across the same factor of four). So there are
        TWO separable faults: a gauge-independent, energy-dependent inconsistency between the two directions -- for which a continuum
        NORMALIZATION convention, per unit energy against per unit momentum, is the natural suspect, sitting in how each module uses
        the orbital from Continuum.generateOrbitalLocalPotential rather than in the prefactor expressions -- and a gauge-dependent,
        energy-independent amplitude defect common to both modules, which example-Of.jl already showed is not a truncation. Suspect
        the first one first. See example-Pa.jl branch b, which carries both energy points.
"""
function pairCreationCrossSection(channels::Array{PhotonScattering.Channel,1}, inPhotonEnergy::Float64, particleEnergy::Float64,
                                  initialLevel::Level)
    wc    = Defaults.getDefaults("speed of light: c");    alpha = 1.0 / wc
    pPlus = ParticleScattering.projectileMomentum(particleEnergy)
    if  inPhotonEnergy < 1.0e-12    return( EmProperty(0., 0.) )    end
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
    factor = 8 * pi^3 * alpha * pPlus / ( inPhotonEnergy * (Basics.twice(initialLevel.J) + 1) )

    return( EmProperty(factor * wCo, factor * wBa) )
end


"""
`PhotonScattering.positronEnergy(finalLevel::Level, initialLevel::Level, inPhotonEnergy::Float64)`
    ... to compute the kinetic energy of the emitted positron from energy conservation, T_+ = omega - 2 m c^2 + B, where B is the
        binding energy gained by the captured electron and is read off the two level energies. Nothing here is adjustable. A
        T::Float64 [a.u.] is returned, non-positive exactly when the channel is below its threshold -- which lies at 2 m c^2 - B and
        so BELOW the 2 m c^2 of free pair creation.
"""
function positronEnergy(finalLevel::Level, initialLevel::Level, inPhotonEnergy::Float64)
    wc = Defaults.getDefaults("speed of light: c")

    return( inPhotonEnergy + initialLevel.energy - finalLevel.energy - 2 * wc * wc )
end
