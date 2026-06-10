
"""
`module  JAC.Einstein`  
... a submodel of JAC that contains all methods for computing Einstein A and B coefficients;
    these computations are treated as level properties with regard to a single multiplet.
"""
module Einstein


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..PhotoEmission, ..Radial, ..TableStrings


"""
`struct  Einstein.Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing Einstein lines and coefficients

    + multipoles              ::Array{EmMultipoles}          ... Specifies the multipoles of the radiation field that are to be included.
    + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
    + photonEnergyShift       ::Float64                      ... An overall energy shift for all photon energies.
    + mimimumPhotonEnergy     ::Float64                      ... minimum transition energy for which (photon) transitions are included into the
                                                                    computation.
    + maximumPhotonEnergy     ::Float64                      ... maximum transition energy for which (photon) transitions are included.
"""
struct Settings  <:  AbstractPropertySettings 
    multipoles                ::Array{EmMultipole,1}
    printBefore               ::Bool 
    lineSelection             ::LineSelection
    photonEnergyShift         ::Float64
    mimimumPhotonEnergy       ::Float64
    maximumPhotonEnergy       ::Float64 
end 


"""
`Einstein.Settings(settings::Einstein.Settings;`
    
            multipoles=..,            printBefore=..,      lineSelection=..,     
            photonEnergyShift=..,     mimimumPhotonEnergy=..,         maximumPhotonEnergy=..)
            
    ... constructor for re-defining the settings::Einstein.Settings.
"""
function Settings(settings::Einstein.Settings;                            multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    
    printBefore::Union{Nothing,Bool}=nothing,                             lineSelection::Union{Nothing,LineSelection}=nothing, 
    photonEnergyShift::Union{Nothing,Float64}=nothing,  
    mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,                  maximumPhotonEnergy::Union{Nothing,Float64}=nothing)

    if  multipoles         == nothing   multipolesx          = settings.multipoles            else  multipolesx = multipoles                   end 
    if  printBefore        == nothing   printBeforex         = settings.printBefore           else  printBeforex = printBefore                 end 
    if  lineSelection      == nothing   lineSelectionx       = set.lineSelection              else  lineSelectionx = lineSelection             end 
    if  photonEnergyShift  == nothing   photonEnergyShiftx   = settings.photonEnergyShift     else  photonEnergyShiftx = photonEnergyShift     end 
    if  mimimumPhotonEnergy== nothing   mimimumPhotonEnergyx = settings.mimimumPhotonEnergy   else  mimimumPhotonEnergyx = mimimumPhotonEnergy end 
    if  maximumPhotonEnergy== nothing   maximumPhotonEnergyx = settings.maximumPhotonEnergy   else  maximumPhotonEnergyx = maximumPhotonEnergy end 
    
    Settings(multipolesx, printBeforex, selectLinesx, selectedLinesx, photonEnergyShiftx, mimimumPhotonEnergyx, maximumPhotonEnergyx)
end


"""
`Einstein.Settings()`  ... constructor for the default values of Einstein line computations.
"""
function Settings()
    Settings(EmMultipole[], false, LineSelection(), 0., 0., 0.)
end


#== """
`Einstein.Settings(multipoles::Array{EmMultipole,1};` printBefore::Bool=false, 
                    selectLines::Bool=false, selectedLines::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                    photonEnergyShift::Float64=0., mimimumPhotonEnergy::Float64=0., maximumPhotonEnergy::Float64=0.) 
    ... keyword constructor to overwrite selected value of Einstein line computations; DON'T use it anymore (11/2019).
"""
function Settings(multipoles::Array{EmMultipole,1}; printBefore::Bool=false, 
                    selectLines::Bool=false, selectedLines::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                    photonEnergyShift::Float64=0., mimimumPhotonEnergy::Float64=0., maximumPhotonEnergy::Float64=0.)
    println("!!! UPDATE: Remove this call and use Einstein.Settings(settings::Einstein.Settings) instead.")
    Settings(multipoles, printBefore, selectLines, selectedLines, photonEnergyShift, mimimumPhotonEnergy, maximumPhotonEnergy)
end  ==#


# `Base.show(io::IO, settings::Einstein.Settings)`  ... prepares a proper printout of the variable settings::Einstein.Settings.
function Base.show(io::IO, settings::Einstein.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:        $(settings.photonEnergyShift)  ")
    println(io, "mimimumPhotonEnergy:      $(settings.mimimumPhotonEnergy)  ")
    println(io, "maximumPhotonEnergy:      $(settings.maximumPhotonEnergy)  ")
end


"""
`struct  Einstein.Channel`  ... defines a type for a single Einstein emission/absorption channel that specifies the multipole, gauge and amplitude.

    + multipole         ::EmMultipole        ... Multipole of the photon emission/absorption.
    + gauge             ::EmGauge            ... Gauge for dealing with the (coupled) radiation field.
    + amplitude         ::Complex{Float64}   ... Amplitude of this multiple channel.
"""
struct Channel 
    multipole           ::EmMultipole
    gauge               ::EmGauge
    amplitude           ::Complex{Float64}
end 


# `Base.show(io::IO, channel::Einstein.Channel)`  ... prepares a proper printout of the variable channel::EinsteinChannel.
function Base.show(io::IO, channel::Einstein.Channel) 
    print(io, "$(channel.multipole) [$(channel.gauge); $(channel.amplitude)") 
end


"""
`struct  Einstein.Line`  
    ... defines a type for a Einstein line that may include the definition of sublines and their corresponding amplitudes.

    + initialLevel   ::Level            ... initial-(state) level
    + finalLevel     ::Level            ... final-(state) level
    + omega          ::Float64          ... photon energy = transition energy
    + photonRate     ::EmProperty       ... Total rate of this line.
    + angularBeta    ::EmProperty       ... Angular beta_2 coefficient.
    + channels       ::Array{Einstein.Channel,1}   ... List of radiative (photon) channels.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    omega            ::Float64
    photonRate       ::EmProperty
    angularBeta      ::EmProperty
    channels         ::Array{Einstein.Channel,1}
end 


"""
`Einstein.Line(initialLevel::Level, finalLevel::Level, photonRate::Float64)`  
    ... constructor an Einstein line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, photonRate::EmProperty)
    Line(initialLevel, finalLevel, 0, 0, photonRate, EmProperty(0., 0.), Einstein.Channel[])
end


# `Base.show(io::IO, line::Einstein.Line)`  ... prepares a proper printout of the variable line::Einstein.Line.
function Base.show(io::IO, line::Einstein.Line) 
    println(io, "initialLevel:         $(line.initialLevel)  ")
    println(io, "finalLevel:           $(line.finalLevel)  ")
    println(io, "omega:                $(line.omega)  ")
    println(io, "photonRate:           $(line.photonRate)  ")
    println(io, "angularBeta:          $(line.angularBeta)  ")
    println(io, "channels:             $(line.channels)  ")
end


"""
`Einstein.computeLines(multiplet::Multiplet, grid::Radial.Grid, settings::Einstein.Settings; output::Bool=true)`  
    ... to compute the Einstein transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{Einstein.Lines} is returned.
"""
function  computeLines(multiplet::Multiplet, grid::Radial.Grid, settings::Einstein.Settings; output::Bool=true)
    println("")
    printstyled("Einstein.computeLines(): The computation of the transition amplitudes and properties starts now ... \n", color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = Einstein.determineLines(multiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    Einstein.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = Einstein.Line[]
    for  line in lines
        newLine = Einstein.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    Einstein.displayRates(stdout, newLines)
    Einstein.displayLifetimes(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Einstein.displayRates(iostream, newLines);    Einstein.displayLifetimes(iostream, newLines)   end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`Einstein.computeAmplitudesProperties(line::Einstein.Line, grid::Radial.Grid, settings::Einstein.Settings)`  
    ... to compute all amplitudes and properties of the given line; a line::Einstein.Line is returned for which the 
        amplitudes and properties are now evaluated.
"""
function  computeAmplitudesProperties(line::Einstein.Line, grid::Radial.Grid, settings::Einstein.Settings)
    global JAC_counter
    newChannels = Einstein.Channel[];    rateC = rateB = 0.
    for channel in line.channels
        #
        amplitude = PhotoEmission.amplitude("emission", channel.multipole, channel.gauge, line.omega, line.finalLevel, line.initialLevel, grid)
        #
        push!( newChannels, Einstein.Channel( channel.multipole, channel.gauge, amplitude) )
        if       channel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    # Calculate the photonrate and angular beta if requested
    wa = 8.0pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1) * 
                                                    (Basics.twice(line.finalLevel.J) + 1)
    photonrate  = EmProperty(wa * rateC, wa * rateB)    
    angularBeta = EmProperty(-9., -9.)
    line = Einstein.Line( line.initialLevel, line.finalLevel, line.omega, photonrate, angularBeta, newChannels)
    return( line )
end



"""
`Einstein.determineChannels(finalLevel::Level, initialLevel::Level, settings::Einstein.Settings)`  
    ... to determine a list of Einstein.Channel for a transitions from the initial to final level and by taking into 
        account the particular settings of for this computation; an Array{Einstein.Channel,1} is returned.
"""
function  determineChannels(finalLevel::Level, initialLevel::Level, settings::Einstein.Settings)
    channels = Einstein.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
            if     string(mp)[1] == 'E'   push!(channels, Einstein.Channel(mp, Basics.Coulomb,   0.) )
                                            push!(channels, Einstein.Channel(mp, Basics.Babushkin, 0.) )  
            elseif string(mp)[1] == 'M'   push!(channels, Einstein.Channel(mp, Basics.Magnetic,  0.) ) 
            else   error("stop a")
            end 
        end
    end
    return( channels )  
end


"""
`Einstein.determineLines(multiplet::Multiplet, settings::Einstein.Settings)`  
    ... to determine a list of EinsteinLine's for transitions between the lines from the given multiplet and by taking 
        into account the particular selections and settings of for this computation; an Array{Einstein.Line,1} is 
        returned. Apart from the level specification, all physical properties are set to zero during the initialization 
        process.
"""
function  determineLines(multiplet::Multiplet, settings::Einstein.Settings)
    lines = Einstein.Line[]
    for  iLevel  in  multiplet.levels
        for  fLevel  in  multiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                omega    = iLevel.energy - fLevel.energy + settings.photonEnergyShift
                if  omega <= settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end  
                channels = determineChannels(fLevel, iLevel, settings) 
                push!( lines, Einstein.Line(iLevel, fLevel, omega, EmProperty(0., 0.), EmProperty(0., 0.), channels) )
            end
        end
    end
    return( lines )
end


#== """
`Einstein.determineLines(multiplet::Multiplet, settings::Einstein.Settings)`  
    ... to determine a list of EinsteinLine's for transitions between the lines from the given multiplet and by taking 
        into account the particular selections and settings of for this computation; an Array{Einstein.Line,1} is 
        returned. Apart from the level specification, all physical properties are set to zero during the initialization 
        process.
"""
function  determineLines(multiplet::Multiplet, settings::Einstein.Settings)
    if    settings.selectLines    selectLines   = true
        selectedLines = Basics.determineSelectedLines(settings.selectedLines, multiplet, multiplet)
    else                          selectLines   = false
    end

    lines = Einstein.Line[]
    for  i = 1:length(multiplet.levels)
        for  f = 1:length(multiplet.levels)
            if  selectLines  &&  !((i,f) in selectedLines )    continue   end
            omega    = multiplet.levels[i].energy - multiplet.levels[f].energy + settings.photonEnergyShift
            if  omega <= settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end  

            channels = determineChannels(multiplet.levels[f], multiplet.levels[i], settings) 
            push!( lines, Einstein.Line(multiplet.levels[i], multiplet.levels[f], omega, EmProperty(0., 0.), EmProperty(0., 0.),
                                        true, channels) )
        end
    end
    return( lines )
end ==#


"""
`Einstein.displayLifetimes(stream::IO, lines::Array{Einstein.Line,1})`  
    ... to list all lifetimes of the selected lines as can be derived from the calculated data. A neat table 
        is printed but nothing is returned otherwise.
"""
function  displayLifetimes(stream::IO, lines::Array{Einstein.Line,1})
    # Determine all initial levels (and their level information) for printing the lifetimes
    ilevels = Int64[];   istr = String[]
    for  i = 1:length(lines)
        ii = lines[i].initialLevel.index;   
        if  !(ii in ilevels)   
            sa  = "  ";    sym = LevelSymmetry( lines[i].initialLevel.J, lines[i].initialLevel.parity )
            sa = sa * TableStrings.center(10, TableStrings.level(lines[i].initialLevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", lines[i].initialLevel.energy)) * "    "
            push!( ilevels, ii);    push!( istr, sa )   
        end
    end
    # Determine the lifetime (in a.u.) of the selected initial levels
    irates = Basics.EmProperty[]
    for  ii in  ilevels
        waCoulomb = waBabushkin = 0.
        for  i = 1:length(lines)
            if   lines[i].initialLevel.index == ii    
                waCoulomb   = waCoulomb   + lines[i].photonRate.Coulomb
                waBabushkin = waBabushkin + lines[i].photonRate.Babushkin
            end
        end
        push!(irates, EmProperty( waCoulomb, waBabushkin) )
    end
    
    nx = 108
    println(stream, " ")
    println(stream, "  Radiative lifetimes (as derived from these computations):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Used Gauge"    ; na=6);                     sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(26, "Lifetime"; na=2);       
    sb = sb * TableStrings.center(26, "[a.u.]"*"          "*TableStrings.inUnits("time"); na=5)
    sa = sa * TableStrings.center(12, "Decay widths"; na=4);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #    
    for  ii = 1:length(ilevels)
        sa = istr[ii]
        sa = sa * "Coulomb          " * @sprintf("%.6e",              1.0/irates[ii].Coulomb)               * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Coulomb) )   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Coulomb) )
        println(stream, sa)
        sa = repeat(" ", length(istr[ii]) )
        sa = sa * "Babushkin        " * @sprintf("%.6e",              1.0/irates[ii].Babushkin)             * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Babushkin) ) * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Babushkin) )
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Einstein.displayLines(stream::IO, lines::Array{Einstein.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of 
        all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{Einstein.Line,1})
    nx = 115
    println(stream, " ")
    println(stream, "  Selected Einstein lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        mpGaugeList = Tuple{Basics.EmMultipole,Basics.EmGauge}[]
        for  i in 1:length(line.channels)
            push!( mpGaugeList, (line.channels[i].multipole, line.channels[i].gauge) )
        end
        sa = sa * TableStrings.multipoleGaugeTupels(50, mpGaugeList)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, " ")
    #
    return( nothing )
end


"""
`Einstein.displayRates(stream::IO, lines::Array{Einstein.Line,1})`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayRates(stream::IO, lines::Array{Einstein.Line,1})
    nx = 145
    println(stream, " ")
    println(stream, "  Einstein coefficients, transition rates and oscillator strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center( 9, "Multipole"; na=1);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=2);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"           "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(11, "Osc. strength"    ; na=3);                 sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(12, "Decay widths"; na=4);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        for  ch in line.channels
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
            sa = sa * TableStrings.center(9,  string(ch.multipole); na=4)
            sa = sa * TableStrings.flushleft(11, string(ch.gauge);  na=2)
            chRate =  8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1) * (abs(ch.amplitude)^2) * 
                                                                (Basics.twice(line.finalLevel.J) + 1)
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein A",  line, chRate)) * "  "
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein B",  line, chRate)) * "    "
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to f",           line, chRate)) * "    "
            sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to decay width", line, chRate)) * "    "
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
