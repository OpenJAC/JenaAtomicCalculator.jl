#
# Two-photon emission of initially-unpolarized atoms.
#
"""
`struct  MultiPhotonDeExcitation.ReducedChannel_2pEmission`  
    ... defines a type for a two-photon emission channel for the emission of photons with well-defined 
        multipolarities.

    + K              ::AngularJ64             ... Rank K of the channel.
    + omega1         ::Float64                ... omega1.
    + omega2         ::Float64                ... omega2.
    + multipole1     ::EmMultipole            ... Multipole M1.
    + multipole2     ::EmMultipole            ... Multipole M2.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + Jsym           ::LevelSymmetry          ... Symmetry of the Green function channel used in the summation.
    + amplitude      ::Complex{Float64}       ... reduced two-photon emission amplitude U^(K, 2gamma, emission) (..)
                                                    associated with the given channel.
"""
struct  ReducedChannel_2pEmission
    K                ::AngularJ64 
    omega1           ::Float64
    omega2           ::Float64 
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge
    Jsym             ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  MultiPhotonDeExcitation.Sharing_2pEmission`  
    ... defines a type for a 2pEmission sharing to help characterize energy sharing between the two emitted photons.

    + omega1           ::Float64         ... Energy of the emitted photon 1.
    + omega2           ::Float64         ... Energy of the emitted photon 2.
    + weight           ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
    + differentialRate ::EmProperty      ... differential rate of this energy sharing.
    + channels         ::Array{MultiPhotonDeExcitation.ReducedChannel_2pEmission,1}  
                                            ... List of 2pEmission channels of this sharing.
"""
struct  Sharing_2pEmission
    omega1             ::Float64
    omega2             ::Float64
    weight             ::Float64
    differentialRate   ::EmProperty
    channels           ::Array{MultiPhotonDeExcitation.ReducedChannel_2pEmission,1}
end


"""
`struct  MultiPhotonDeExcitation.Line_2pEmission`  
    ... defines a type for a two-photon absorption line by monochromatic light that may include the definition of channels.

    + initialLevel     ::Level                ... initial-(state) level
    + finalLevel       ::Level                ... final-(state) level
    + totalRate        ::EmProperty           ... Total rate for the two-photon transition. 
    + sharings         ::Array{MultiPhotonDeExcitation.Sharing_2pEmission,1}  
                                                ... List of MultiPhotonDeExcitation.Sharing_2pEmission's of this line.
"""
struct  Line_2pEmission
    initialLevel       ::Level
    finalLevel         ::Level
    totalRate          ::EmProperty
    sharings           ::Array{MultiPhotonDeExcitation.Sharing_2pEmission,1} 
end


# `Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pEmission)`  
#   ... prepares a proper printout of the variable line::MultiPhotonDeExcitation.Line_2pEmission.
function Base.show(io::IO, line::MultiPhotonDeExcitation.Line_2pEmission) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "totalRate:         $(line.totalRate)  ")
    println(io, "sharings:          $(line.sharings)  ")
end


"""
`MultiPhotonDeExcitation.computeAmplitudes_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                                                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)` 
    ... to compute all amplitudes of the given line; a line::MultiPhotonDeExcitation.Line_2pEmission 
        is returned for which the amplitudes are now evaluated.
"""
function  computeAmplitudes_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)
    newSharings = MultiPhotonDeExcitation.Sharing_2pEmission[]
    for sharing in line.sharings
        @show sharing.omega1, sharing.omega2
        newChannels = MultiPhotonDeExcitation.ReducedChannel_2pEmission[]
        for ch in sharing.channels
            amplitude = MultiPhotonDeExcitation.computeReducedAmplitudeEmission(ch.K, line.finalLevel, ch.omega2, ch.multipole2, 
                                                            ch.Jsym, ch.omega1, ch.multipole1, line.initialLevel, ch.gauge, grid, settings.green)
            push!( newChannels, MultiPhotonDeExcitation.ReducedChannel_2pEmission(ch.K, ch.omega1, ch.omega2, ch.multipole1, ch.multipole2, 
                                                                                    ch.gauge, ch.Jsym, amplitude) )
            @show ch, amplitude
        end
        push!( newSharings, MultiPhotonDeExcitation.Sharing_2pEmission( sharing.omega1, sharing.omega2, sharing.weight, EmProperty(0.), newChannels) )
    end
    newLine = MultiPhotonDeExcitation.Line_2pEmission( line.initialLevel, line.finalLevel, EmProperty(0.), newSharings)
    return( newLine )
end


"""
`MultiPhotonDeExcitation.computeEnergyDiffCs(sharing::MultiPhotonDeExcitation.Sharing_2pEmission, line::MultiPhotonDeExcitation.Line_2pEmission, 
                                                gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)`  
    ... to compute the energy differential cross section for the given line and the energy sharing omega1. A dcs::Float64 is returned.
"""
function computeEnergyDiffCs(sharing::MultiPhotonDeExcitation.Sharing_2pEmission, line::MultiPhotonDeExcitation.Line_2pEmission, 
                                gauge::EmGauge, settings::MultiPhotonDeExcitation.Settings)
    omega1 = sharing.omega1;    omega2 = sharing.omega2;    dcs = 0.
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            symi        = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
            Klist       = oplus(line.finalLevel.J, line.initialLevel.J)
            for Jsym in symmetries
                wk = ComplexF64(0.)
                for  K in Klist
                    wk =  wk + MultiPhotonDeExcitation.getReducedAmplitudeEmission(K, line.finalLevel, omega2, mp2, Jsym, omega1, mp1, 
                                                                                        line.initialLevel, gauge, line.sharings) +
                            AngularMomentum.phaseFactor([K, +1, line.finalLevel.J, +1, line.initialLevel.J]) * 
                                MultiPhotonDeExcitation.getReducedAmplitudeEmission(K, line.finalLevel, omega1, mp1, Jsym, omega2, mp2, 
                                                                                        line.initialLevel, gauge, line.sharings)
                end
                dcs = dcs + abs( wk )^2
            end
        end
    end
    dcs = dcs * 2pi * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) * omega1 * omega2
    
    return( dcs )
end



"""
`MultiPhotonDeExcitation.computeLines(process::TwoPhotonEmission, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; output=true)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1} is returned.
"""
function  computeLines(process::TwoPhotonEmission, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings; output=true)
    println("")
    printstyled("MultiPhotonDeExcitation.computeLines(::TwoPhotonEmission): The computation of amplitudes starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    lines = MultiPhotonDeExcitation.determineLines_2pEmission(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    MultiPhotonDeExcitation.displayLines_2pEmission(lines)    end
    # Calculate all amplitudes
    newLines = MultiPhotonDeExcitation.Line_2pEmission[]
    for  line in lines
        newLine = MultiPhotonDeExcitation.computeAmplitudes_2pEmission(line, grid, settings) 
        newLine = MultiPhotonDeExcitation.computeProperties_2pEmission(newLine, grid, settings)
        push!( newLines, newLine)
    end
    # Print all results to screen
    MultiPhotonDeExcitation.displayTotalRates_2pEmission(stdout, newLines, settings)
    MultiPhotonDeExcitation.displayDifferentialRates_2pEmission(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   MultiPhotonDeExcitation.displayTotalRates_2pEmission(iostream, newLines, settings)
                        MultiPhotonDeExcitation.displayDifferentialRates_2pEmission(iostream, newLines, settings)     end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`MultiPhotonDeExcitation.computeProperties_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                                                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)` 
    ... to compute all properties of the given line; a line::MultiPhotonDeExcitation.Line_2pEmission 
        is returned for which the properties are now evaluated.
"""
function  computeProperties_2pEmission(line::MultiPhotonDeExcitation.Line_2pEmission, 
                                        grid::Radial.Grid, settings::MultiPhotonDeExcitation.Settings)
    newSharings = MultiPhotonDeExcitation.Sharing_2pEmission[];     totalRate = EmProperty(0.)
    for sharing in line.sharings
        @show sharing.omega1, sharing.omega2
        diffRate = EmProperty(0.)
        if  Basics.UseCoulomb  in  settings.gauges
                dr_Cou = MultiPhotonDeExcitation.computeEnergyDiffCs(sharing, line, EmGauge("Coulomb"), settings)
        else    dr_Cou = 0.
        end
        if  Basics.UseBabushkin  in  settings.gauges
                dr_Bab = MultiPhotonDeExcitation.computeEnergyDiffCs(sharing, line, EmGauge("Babushkin"), settings)
        else    dr_Bab = 0.
        end
        diffRate = diffRate + EmProperty(dr_Cou, dr_Bab)
        push!( newSharings, MultiPhotonDeExcitation.Sharing_2pEmission( sharing.omega1, sharing.omega2, sharing.weight, diffRate, sharing.channels) )
        totalRate = totalRate + sharing.weight * diffRate
    end
    # Calculate the totalRate 
    newLine = MultiPhotonDeExcitation.Line_2pEmission( line.initialLevel, line.finalLevel, totalRate, newSharings)
    return( newLine )
end


"""
`MultiPhotonDeExcitation.computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, 
                                                            omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                            omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                            gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})`  
    ... to compute the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) by means of the
        given Green function channels. An amplitude::Complex{Float64} is returned.
"""
function computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                            omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                        gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})
    U = Complex(0.);    found = false
    for channel in greenChannels
        if  Jsym == channel.symmetry
            found = true
            for (i, nuLevel) in enumerate(channel.gMultiplet.levels)
                U = U + PhotoEmission.amplitude(Emission(), multipole2, gauge, omega2, finalLevel, nuLevel, grid) *
                        PhotoEmission.amplitude(Emission(), multipole1, gauge, omega1, nuLevel, initialLevel, grid) / 
                        (initialLevel.energy + omega1 - nuLevel.energy)
            end
        end
    end 
    
    if    found
        wx = AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipole2.L), AngularJ64(multipole1.L), Jsym.J)
        U  = U * wx
    else  println("No Green function cannel found for amplitude U^{K, 2gamma emission} " *
                    "(K, Jf, omega2, multipole2, Jsym = $Jsym, omega1, multipole1, Ji) ")
    end 
    
    return( U )
end


"""
`MultiPhotonDeExcitation.getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, 
                                                        omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                        omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                        gauge::EmGauge, sharings:Array{MultiPhotonDeExcitation.Sharing_2pEmission,1})`  
    ... to get/return the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) 
        from the calculated list of channels. An amplitude::Complex{Float64} is returned.
"""
function getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                        omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                                        gauge::EmGauge, sharings::Array{MultiPhotonDeExcitation.Sharing_2pEmission,1})
    U = Complex(0.);    found = false
    for sharing in sharings
        for channel in sharing.channels
            ## omega1 == channel.omega1  &&  omega2 == channel.omega2  &&  
            if  K == channel.K         &&  abs(omega1 - channel.omega1) < 1.e-10  &&   abs(omega2 - channel.omega2) < 1.e-10 && 
                Jsym == channel.Jsym   &&  multipole1 == channel.multipole1  &&  multipole2 == channel.multipole2            && 
                gauge == channel.gauge
                ## (gauge == channel.gauge  ||  EmGauge("Magnetic")  == channel.gauge)
                U = channel.amplitude;    found = true;     continue
            end
        end
    end 
    
    if    found                                
            println("U^{K, 2gamma emission} ($gauge..) = $U found. ")
    else  println("No U^{K, 2gamma emission} ($gauge, $K, Jf, $omega2, $multipole2, Jsym = $Jsym, $omega1, $multipole1, Ji) amplitude found.")
    end 
    
    return( U )
end


"""
`MultiPhotonDeExcitation.determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                    settings::MultiPhotonDeExcitation.Settings)`
    ... to determine a list of MultiPhotonDeExcitation.Line_2pEmission's for transitions between the levels from the given 
        initial- and final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{MultiPhotonDeExcitation.Line_2pEmission,1} is returned. Apart from the level specification, all physical 
        properties are set to zero during this initialization process.  
"""
function  determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::MultiPhotonDeExcitation.Settings)
    lines = MultiPhotonDeExcitation.Line_2pEmission[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy   = iLevel.energy - fLevel.energy
                sharings = MultiPhotonDeExcitation.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                push!( lines, MultiPhotonDeExcitation.Line_2pEmission(iLevel, fLevel, EmProperty(0., 0.,), sharings) )
            end
        end
    end
    return( lines )
end


"""
`MultiPhotonDeExcitation.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::MultiPhotonDeExcitation.Settings)`  
    ... to determine a list of MultiPhotonDeExcitation 2pEmission Sharing's and Channel's for a transitions from the initial to 
        final level and by taking into account the particular settings of for this computation; 
        an Array{MultiPhotonDeExcitation.Sharing_2pEmission,1} is returned.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::MultiPhotonDeExcitation.Settings)
    sharings = MultiPhotonDeExcitation.Sharing_2pEmission[];    eSharings = Basics.determineEnergySharings(energy, settings.NoEnergySharings) 
    for  es in eSharings
        omega1    = es[1];    omega2 = es[2];    weight = es[3] 
        channels  = MultiPhotonDeExcitation.ReducedChannel_2pEmission[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                Klist       = oplus(symf.J, symi.J)
                for  symx in symmetries
                    for  gauge in settings.gauges
                        # Include further restrictions if appropriate
                        if     (string(mp1)[1] == 'E' || string(mp2)[1] == 'E')  &&   gauge == Basics.UseCoulomb
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Coulomb, symx, 0.) )    end 
                        elseif (string(mp1)[1] == 'E' || string(mp2)[1] == 'E')  &&   gauge == Basics.UseBabushkin    
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Babushkin, symx, 0.) )  end
                        elseif (string(mp1)[1] == 'M' && string(mp2)[1] == 'M')
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Magnetic, symx, 0.) )   end
                        end
                    end
                end
            end
        end
        push!(sharings, MultiPhotonDeExcitation.Sharing_2pEmission(omega1, omega2, weight, EmProperty(0., 0.), channels) )
    end
    return( sharings )  
end


"""
`MultiPhotonDeExcitation.displayDifferentialRates_2pEmission(stream::IO, lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1}, 
                                                                settings::MultiPhotonDeExcitation.Settings)`  
    ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayDifferentialRates_2pEmission(stream::IO, lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1}, 
                                                settings::MultiPhotonDeExcitation.Settings)
    #
    # First, print lines and sharings
    nx = 130
    println(stream, " ")
    println(stream, "  Energy-differential rates of selected two-photon emission lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=4);              
    sb = sb * TableStrings.flushleft(38, "  i -- f          omega1        omega2"; na=4)
    sa = sa * TableStrings.center(14, "Weight"; na=0);                            sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(34, "Cou -- diff. rate -- Bab"; na=3)      
    sb = sb * TableStrings.center(34, TableStrings.inUnits("rate") * "        " * 
                                        TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        #
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega1))                   * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega2))                   * "    "
            sb = sb * @sprintf("%.4e",                                              sharing.weight)                    * "       "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Coulomb))   * "   "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Babushkin)) * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`MultiPhotonDeExcitation.displayLines_2pEmission(lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1})`  
    ... to display a list of lines, sharings & (reduced) channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines_2pEmission(lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1})
    nx = 157
    println(" ")
    println("  Selected two-photon emission lines with given photon splitting:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(34, "Energies (all in " * TableStrings.inUnits("energy") *")"; na=5);              
    sb = sb * TableStrings.flushleft(34, "  i -- f      omega1      omega2  "; na=5)
    sa = sa * TableStrings.flushleft(77, "List of multipoles, gauges & intermediate level symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(77, "(K-rank, multipole_1, Jsym, multipole_2, gauge)           "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy)) * "  "
        #
        for  sharing  in  line.sharings
            sb =      @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega1))   * "  "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega2))   * "     "
            mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
            for  channel in  sharing.channels
                push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
            end
            wa = TableStrings.twoPhotonGaugeTupels(75, mpGaugeList)
            sc = sa * sb * wa[1];    println( sc )  
            for  i = 2:length(wa)
                sc = TableStrings.hBlank( length(sa*sb) ) * wa[i];    println( sc )
            end
        end
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`MultiPhotonDeExcitation.displayTotalRates_2pEmission(stream::IO, lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1}, 
                                                        settings::MultiPhotonDeExcitation.Settings)`  
    ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayTotalRates_2pEmission(stream::IO, lines::Array{MultiPhotonDeExcitation.Line_2pEmission,1}, settings::MultiPhotonDeExcitation.Settings)
    #
    # First, print lines and sharings
    nx = 88
    println(stream, " ")
    println(stream, "  Total rates of selected two-photon emission lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(21)
    sa = sa * TableStrings.center(12, "i--Energy--f"; na=3)               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(30, "Cou --   total rate   -- Bab"; na=4)      
    sb = sb * TableStrings.center(30, TableStrings.inUnits("rate") * "          " * 
                                        TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))                 * "      "
        #
        sb = sa * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Coulomb))   * "     "
        sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Babushkin)) * "   "
        println(stream,  sb )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end
