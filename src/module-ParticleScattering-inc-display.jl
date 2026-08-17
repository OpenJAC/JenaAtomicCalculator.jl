
# The display methods of the ParticleScattering module. They print, and nothing else; every number they show is taken
# from an Event that already carries it.


"""
`ParticleScattering.displayEvents(stream::IO, events::Array{ParticleScattering.Event,1})`
    ... to list the scattering events that have been selected, before the computation starts, so that a mis-specified
        selection is seen at once rather than after the partial waves have been generated. A neat table is printed but
        nothing is returned otherwise.
"""
function displayEvents(stream::IO, events::Array{ParticleScattering.Event,1})
    nx = 116
    println(stream, " ")
    println(stream, "  Selected particle-scattering events:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);          sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=4);          sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(16, "Projectile"; na=2);         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(20, "Interaction"; na=2);        sb = sb * TableStrings.hBlank(22)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  event in events
        sa  = "  "
        sa  = sa * TableStrings.center(16, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=2)
        sa  = sa * TableStrings.center(16, TableStrings.symmetries_if(LevelSymmetry(event.initialLevel.J, event.initialLevel.parity),
                                                                     LevelSymmetry(event.finalLevel.J, event.finalLevel.parity)); na=4)
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", event.impactEnergy)) * "    "
        sa  = sa * TableStrings.center(16, string(typeof(event.projectile).name.name); na=2)
        sa  = sa * TableStrings.center(20, string(typeof(event.interaction).name.name); na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`ParticleScattering.displayPhaseShifts(stream::IO, events::Array{ParticleScattering.Event,1})`
    ... to display the Dirac phase shifts of every event, for both spin-orbit partners of each orbital angular momentum.
        Showing them side by side makes the spin-orbit splitting of the phase shifts visible, which is the quantity the
        spin-flip amplitude and hence the Sherman function are built from. A neat table is printed but nothing is
        returned otherwise.
"""
function displayPhaseShifts(stream::IO, events::Array{ParticleScattering.Event,1})
    nx = 88
    println(stream, " ")
    println(stream, "  Dirac phase shifts of the scattered projectile:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "Impact energy"; na=4)
    sa = sa * TableStrings.center( 8, "l";  na=4)
    sa = sa * TableStrings.center(20, "delta (j = l+1/2)"; na=4)
    sa = sa * TableStrings.center(20, "delta (j = l-1/2)"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  event in events
        lMax = isempty(event.partialWaves) ? -1 :
               maximum( Basics.subshell_l(Subshell(101, pw.kappa)) for pw in event.partialWaves )
        for  l = 0:lMax
            dMinus = ParticleScattering.phaseShift(event.partialWaves, -l - 1)
            dPlus  = (l == 0) ? 0. : ParticleScattering.phaseShift(event.partialWaves, l)
            sa     = "  " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", event.impactEnergy)) * "     "
            sa     = sa * TableStrings.center(8, string(l); na=4)
            sa     = sa * @sprintf("%.8e", dMinus) * "        "
            sa     = sa * ( (l == 0) ? TableStrings.center(20, "--"; na=2) : @sprintf("%.8e", dPlus) )
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`ParticleScattering.displayCrossSections(stream::IO, events::Array{ParticleScattering.Event,1})`
    ... to display the differential cross sections and Sherman functions of every event, at each requested scattering
        angle. A neat table is printed but nothing is returned otherwise.
"""
function displayCrossSections(stream::IO, events::Array{ParticleScattering.Event,1})
    nx = 118
    println(stream, " ")
    println(stream, "  Differential cross sections and spin polarization:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "theta"; na=2);           sb = sb * TableStrings.center(12, "[rad]"; na=2)
    sa = sa * TableStrings.center(12, "phi";   na=2);           sb = sb * TableStrings.center(12, "[rad]"; na=2)
    sa = sa * TableStrings.center(18, "d sigma / d Omega"; na=2)
    sb = sb * TableStrings.center(18, "[a.u.]"; na=2)
    sa = sa * TableStrings.center(16, "Sherman S"; na=2);       sb = sb * TableStrings.hBlank(18)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  event in events
        for  obs in event.angular
            sa = "  "
            sa = sa * TableStrings.center(16, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=2)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", event.impactEnergy)) * "    "
            sa = sa * @sprintf("%.4e", obs.theta) * "    "
            sa = sa * @sprintf("%.4e", obs.phi)   * "      "
            sa = sa * @sprintf("%.6e", obs.dcs)   * "      "
            sa = sa * @sprintf("%.6e", obs.sherman)
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`ParticleScattering.displayIntegratedCrossSections(stream::IO, events::Array{ParticleScattering.Event,1})`
    ... to display the angle-integrated cross sections of every event: the elastic total together with the first
        (momentum-transfer) and second (viscosity) transport cross sections. A neat table is printed but nothing is
        returned otherwise.
"""
function displayIntegratedCrossSections(stream::IO, events::Array{ParticleScattering.Event,1})
    nx = 96
    println(stream, " ")
    println(stream, "  Integrated and transport cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "Impact energy"; na=4)
    sa = sa * TableStrings.center(20, "sigma_elastic"; na=2)
    sa = sa * TableStrings.center(20, "sigma_1 (moment.)"; na=2)
    sa = sa * TableStrings.center(20, "sigma_2 (viscos.)"; na=2)
    println(stream, sa)
    println(stream, "  " * TableStrings.hBlank(16) * TableStrings.center(64, "all in [a.u.]"; na=2))
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  event in events
        sa = "  " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", event.impactEnergy)) * "     "
        sa = sa * @sprintf("%.8e", event.integrated.sigmaElastic)          * "        "
        sa = sa * @sprintf("%.8e", event.integrated.sigmaMomentumTransfer) * "        "
        sa = sa * @sprintf("%.8e", event.integrated.sigmaViscosity)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end
