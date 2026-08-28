
"""
`module  JAC.RadiativeAuger`  
... a submodel of JAC that contains all methods for computing radiative Auger and autoionization amplitudes and rates.
"""
module RadiativeAuger

using Printf, ..AngularMomentum, ..AtomicState, ..AutoIonization, ..Basics, ..Defaults, ..ManyElectron, 
                ..Radial, ..PhotoEmission, ..TableStrings

"""
`struct  RadiativeAuger.Settings`  ... defines a type for the settings in estimating radiative-Auger and autoionization rates.

    + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
    + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
    + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
    + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + operator                ::AbstractEeInteraction   ... Auger operator that is to be used for evaluating the Auger amplitudes: 
                                                            allowed values are: CoulombInteraction(), BreitInteraction(), ...
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
"""
struct Settings
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    green                     ::Array{GreenChannel,1}
    NoEnergySharings          ::Int64 
    maxKappa                  ::Int64  
    printBefore               ::Bool   
    operator                  ::AbstractEeInteraction 
    lineSelection             ::LineSelection 
end 


"""
`RadiativeAuger.Settings()`  ... constructor for the default values of RadiativeAuger line computations
"""
function Settings()
    Settings(EmMultipole[], UseGauge[], GreenChannel[], 0, 3, false, CoulombInteraction(), LineSelection())
end


"""
`RadiativeAuger.Settings(set::RadiativeAuger.Settings;`

    multipoles..., gauges..., green..., NoEnergySharings..., maxKappa..., printBefore..., operator..., lineSelection...)

    ... constructor for modifying the given RadiativeAuger.Settings by 'overwriting' the previously selected parameters.
        A settings::RadiativeAuger.Settings is returned.
"""
function Settings(set::RadiativeAuger.Settings;
                        multipoles::Union{Nothing,Array{EmMultipole,1}} = nothing,
                        gauges::Union{Nothing,Array{UseGauge,1}} = nothing,
                        green::Union{Nothing,Array{GreenChannel,1}} = nothing,
                        NoEnergySharings::Union{Nothing,Int64} = nothing,
                        maxKappa::Union{Nothing,Int64} = nothing,
                        printBefore::Union{Nothing,Bool} = nothing,
                        operator::Union{Nothing,AbstractEeInteraction} = nothing,
                        lineSelection::Union{Nothing,LineSelection} = nothing )

    if  isnothing(multipoles)         multipolesx       = set.multipoles         else   multipolesx       = multipoles         end
    if  isnothing(gauges)             gaugesx           = set.gauges             else   gaugesx           = gauges             end
    if  isnothing(green)              greenx            = set.green              else   greenx            = green              end
    if  isnothing(NoEnergySharings)   NoEnergySharingsx = set.NoEnergySharings   else   NoEnergySharingsx = NoEnergySharings   end
    if  isnothing(maxKappa)           maxKappax         = set.maxKappa           else   maxKappax         = maxKappa           end
    if  isnothing(printBefore)        printBeforex      = set.printBefore        else   printBeforex      = printBefore        end
    if  isnothing(operator)           operatorx         = set.operator           else   operatorx         = operator           end
    if  isnothing(lineSelection)      lineSelectionx    = set.lineSelection      else   lineSelectionx    = lineSelection      end

    Settings( multipolesx, gaugesx, greenx, NoEnergySharingsx, maxKappax, printBeforex, operatorx,
              lineSelectionx )
end


# `Base.show(io::IO, settings::RadiativeAuger.Settings)`  ... prepares a proper printout of the variable settings::RadiativeAuger.Settings.
function Base.show(io::IO, settings::RadiativeAuger.Settings) 
    println(io, "multipoles:                   $(settings.multipoles)  ")
    println(io, "gauges:                       $(settings.gauges)  ")
    println(io, "green:                         (settings.green)  ")
    println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
    println(io, "maxKappa:                     $(settings.maxKappa)  ")
    println(io, "printBefore:                  $(settings.printBefore)  ")
    println(io, "operator :                    $(settings.operator )  ")
    println(io, "lineSelection:                $(settings.lineSelection)  ")
end


"""
`struct  RadiativeAuger.ReducedChannel`  
    ... defines a type for a RadiativeAuger (reduced) channel to help characterize a scattering (continuum) state of many 
        electron-states with a single free electron.

    + multipole      ::EmMultipole          ... Multipole of the photon absorption.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + omega          ::Float64              ... photon energy
    + epsilon        ::Float64              ... (free) electron energy
    + kappa          ::Int64                ... partial-wave of the free electron
    + phase          ::Float64              ... phase of the partial wave
    + tSymmetry      ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
"""
struct  ReducedChannel
    multipole        ::EmMultipole  
    gauge            ::EmGauge
    omega            ::Float64 
    epsilon          ::Float64   
    kappa            ::Int64 
    phase            ::Float64   
    tSymmetry        ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  Sharing`  
    ... defines a type for a RadiativeAuger sharing to help characterize energy sharing between the emitted photon and
        the scattering (continuum) state with a single free electron.

    + omega            ::Float64         ... Energy of the emitted photon.
    + epsilon          ::Float64         ... Energy of the (outgoing free) electron.
    + weight           ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
    + differentialRate ::EmProperty      ... differential rate of this energy sharing.
    + channels         ::Array{RadiativeAuger.ReducedChannel,1}  ... List of RadiativeAuger channels of this line.
"""
struct  Sharing
    omega              ::Float64
    epsilon            ::Float64
    weight             ::Float64
    differentialRate   ::EmProperty
    channels           ::Array{RadiativeAuger.ReducedChannel,1}
end




# RadiativeAuger line between an initial and final (bound-state) level
"""
`struct  Line`  ... defines a type for a radiative Auger line that may include the definition of sublines and their corresponding amplitudes.

    + initialLevel   ::Level          ... initial-(state) level
    + finalLevel     ::Level          ... final-(state) level
    + totalRate      ::EmProperty     ... Total rate of this line.
    + sharings       ::Array{RadiativeAuger.Sharing,1}  ... List of RadiativeAuger sharings of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    totalRate        ::EmProperty
    sharings         ::Array{RadiativeAuger.Sharing,1}
end 


# `Base.show(io::IO, line::RadiativeAuger.Line)`  ... prepares a proper printout of the variable line::RadiativeAuger.Line.
    function Base.show(io::IO, line::RadiativeAuger.Line) 
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "totalRate:              $(line.totalRate)  ")
    println(io, "sharings:               $(line.sharings)  ")
end


"""
`RadiativeAuger.computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)`  
    ... to compute all amplitudes and properties of the given line; a line::RadiativeAuger.Line is returned for which the amplitudes 
        and properties have now been evaluated.
"""
function  computeAmplitudesProperties(line::RadiativeAuger.Line, grid::Radial.Grid, settings::RadiativeAuger.Settings)
    newSharings = RadiativeAuger.Sharing[]
    for sharing in line.sharings
        newChannels = RadiativeAuger.ReducedChannel[]
        for ch in sharing.channels
            # Generate a continuum orbital
            phase     = 1.3
            amplitude = 1.0im 
            push!( newChannels, ReducedChannel(ch.multipole, ch.gauge, ch.omega, ch.epsilon, ch.kappa, phase, ch.tSymmetry, amplitude) )
        end
        # Calculate the differential rate 
        diffRate = EmProperty(-1., -1.)
        push!( newSharings, RadiativeAuger.Sharing( sharing.omega, sharing.epsilon, sharing.weight, diffRate, newChannels) )
    end
    # Calculate the totalRate 
    totalRate = EmProperty(-1., -1.)
    line = RadiativeAuger.Line( line.initialLevel, line.finalLevel, totalRate, newSharings)
    return( line )
end


"""
`RadiativeAuger.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                settings::RadiativeAuger.Settings; output=true)`  
    ... to compute the radiative Auger transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{RadiativeAuger.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                        settings::RadiativeAuger.Settings; output=true)
    println("")
    printstyled("RadiativeAuger.computeLines(): The computation of radiative Auger rates starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    lines = RadiativeAuger.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    RadiativeAuger.displayLines(lines)    end
    # Calculate all amplitudes and requested properties
    newLines = RadiativeAuger.Line[]
    for  line in lines
        newLine = RadiativeAuger.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    RadiativeAuger.displayTotalRates(stdout, lines, settings)
    RadiativeAuger.displayDifferentialRates(stdout, lines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   RadiativeAuger.displayTotalRates(iostream, lines, settings)
                        RadiativeAuger.displayDifferentialRates(iostream, lines, settings)     end
    #
    if    output    return( lines )
    else            return( nothing )
    end
end


"""
`RadiativeAuger.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)` 
    ... to determine a list of RadiativeAuger.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{RadiativeAuger.Line,1} is 
        returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RadiativeAuger.Settings)
    lines = RadiativeAuger.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy   = iLevel.energy - fLevel.energy
                channels = RadiativeAuger.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                push!( lines, RadiativeAuger.Line(iLevel, fLevel, EmProperty(0., 0.,), channels) )
            end
        end
    end
    return( lines )
end


"""
`RadiativeAuger.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)`  
    ... to determine a list of RadiativeAuger Sharing's and Channel's for a transitions from the initial to final level and by taking into 
        account the particular settings of for this computation; an Array{RadiativeAuger.Sharing,1} is returned.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::RadiativeAuger.Settings)
    sharings = RadiativeAuger.Sharing[];    eSharings = Basics.determineEnergySharings(energy, settings.NoEnergySharings) 
    for  es in eSharings
        omega     = es[1];    epsilon = es[2];    weight = es[3] 
        channels  = RadiativeAuger.ReducedChannel[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            for  gauge in settings.gauges
                symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                for  symt in symList
                    kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                    for  kappa in kappaList
                        # Include further restrictions if appropriate
                        if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                            push!(channels, ReducedChannel(mp, Basics.Coulomb,   omega, epsilon, kappa, 0., symt, Complex(0.)) )
                        elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                            push!(channels, ReducedChannel(mp, Basics.Babushkin, omega, epsilon, kappa, 0., symt, Complex(0.)) )  
                        elseif string(mp)[1] == 'M'                                
                            push!(channels, ReducedChannel(mp, Basics.Magnetic,  omega, epsilon, kappa, 0., symt, Complex(0.)) ) 
                        end 
                    end
                end
            end
        end
        push!(sharings, RadiativeAuger.Sharing(omega, epsilon, weight, EmProperty(0., 0.), channels) )
    end
    return( sharings )  
end


"""
`RadiativeAuger.displayDifferentialRates(stream::IO, lines::Array{RadiativeAuger.Line,1}, settings::RadiativeAuger.Settings)`  
    ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayDifferentialRates(stream::IO, lines::Array{RadiativeAuger.Line,1}, settings::RadiativeAuger.Settings)
    #
    # First, print lines and sharings
    nx = 130
    println(stream, " ")
    println(stream, "  Energy-differential rates of selected radiative-Auger lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=4);              
    sb = sb * TableStrings.flushleft(38, "  i -- f          omega        epsilon"; na=4)
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
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))                    * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon))                  * "    "
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
`RadiativeAuger.displayLines(lines::Array{RadiativeAuger.Line,1})`  
    ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(lines::Array{RadiativeAuger.Line,1})
    nx = 170
    println(" ")
    println("  Selected radiative-Auger lines:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(34, "Energy  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(34, "i -- f        omega      e_Auger "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "  "
        #
        for  sharing  in  line.sharings
            sb =      @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))   * "  "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon)) * "    "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(sharing.channels)
                push!( kappaMultipoleSymmetryList, (sharing.channels[i].kappa, sharing.channels[i].multipole, sharing.channels[i].gauge, 
                                                    sharing.channels[i].tSymmetry) )
            end
            wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
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
`RadiativeAuger.displayTotalRates(stream::IO, lines::Array{RadiativeAuger.Line,1}, settings::RadiativeAuger.Settings)`  
    ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayTotalRates(stream::IO, lines::Array{RadiativeAuger.Line,1}, settings::RadiativeAuger.Settings)
    #
    # First, print lines and sharings
    nx = 88
    println(stream, " ")
    println(stream, "  Total rates of selected radiative-Auger lines:")
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

end # module
