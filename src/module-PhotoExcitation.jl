
"""
`module  JAC.PhotoExcitation`  
... a submodel of JAC that contains all methods for computing photo-excitation properties between some initial 
    and final-state multiplets.
"""
module PhotoExcitation

using Printf, ..AngularMomentum, ..Basics,  ..Basics,  ..Defaults, ..ManyElectron, ..Radial, ..PhotoEmission, ..TableStrings

"""
`struct  PhotoExcitation.Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing photo-excitation  lines.

    + multipoles              ::Array{EmMultipole,1}    ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge,1}       ... Specifies the gauges to be included into the computations.
    + calcForStokes           ::Bool                    ... True, if the excitation cross sections are to be calculated (and false otherwise)
                                                            for given Stokes parameter of the incident plane-wave photons.
    + calcPhotonDm            ::Bool                    ... True, if the photon density matrix of a subsequently emitted fluorescence photon 
                                                            is to be calculated and false otherwise. 
    + calcTensors             ::Bool                    ... True, if statistical tensors of the excited atom are to be calculated, false otherwise. 
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + photonEnergyShift       ::Float64                 ... An overall energy shift for all photon energies.
    + mimimumPhotonEnergy     ::Float64                 ... minimum transition energy for which (photon) transitions are included into the
                                                            computation.
    + maximumPhotonEnergy     ::Float64                 ... maximum transition energy for which (photon) transitions are included.
    + stokes                  ::ExpStokes               ... Stokes parameters of the incident radiation.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge,1}
    calcForStokes             ::Bool
    calcPhotonDm              ::Bool  
    calcTensors               ::Bool  
    printBefore               ::Bool 
    lineSelection             ::LineSelection
    photonEnergyShift         ::Float64
    mimimumPhotonEnergy       ::Float64   
    maximumPhotonEnergy       ::Float64     
    stokes                    ::ExpStokes 
end 


"""
`PhotoExcitation.Settings()`  ... 'empty' constructor for the default values of photo-excitation line computations
"""
function Settings()
    Settings(EmMultipole[], UseGauge[], false, false, false, false, LineSelection(), 0., 0., 1.0e6, Basics.ExpStokes())
end


"""
`PhotoExcitation.Settings(set::PhotoExcitation.Settings;`

        multipoles=..,          gauges=..,                  calcForStokes=..,           calcPhotonDm=..,    
        calcTensors=..,         printBefore=..,             lineSelection=..,    
        photonEnergyShift=..,   mimimumPhotonEnergy=..,     maximumPhotonEnergy=..,     stokes=..)
                    
    ... constructor for modifying the given PhotoExcitation.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoExcitation.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,        gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    calcForStokes::Union{Nothing,Bool}=nothing,                     calcPhotonDm::Union{Nothing,Bool}=nothing,    
    calcTensors::Union{Nothing,Bool}=nothing,                       printBefore::Union{Nothing,Bool}=nothing,  
    lineSelection::Union{Nothing,LineSelection}=nothing, 
    photonEnergyShift::Union{Nothing,Float64}=nothing,              mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,     
    maximumPhotonEnergy::Union{Nothing,Float64}=nothing,            stokes::Union{Nothing,ExpStokes}=nothing)  
    
    if  isnothing(multipoles)   multipolesx          = set.multipoles              else  multipolesx          = multipoles            end 
    if  isnothing(gauges)   gaugesx              = set.gauges                  else  gaugesx              = gauges                end 
    if  isnothing(calcForStokes)   calcForStokesx       = set.calcForStokes           else  calcForStokesx       = calcForStokes         end 
    if  isnothing(calcPhotonDm)   calcPhotonDmx        = set.calcPhotonDm            else  calcPhotonDmx        = calcPhotonDm          end 
    if  isnothing(calcTensors)   calcTensorsx         = set.calcTensors             else  calcTensorsx         = calcTensors           end 
    if  isnothing(printBefore)   printBeforex         = set.printBefore             else  printBeforex         = printBefore           end 
    if  isnothing(lineSelection)   lineSelectionx       = set.lineSelection           else  lineSelectionx       = lineSelection         end 
    if  isnothing(photonEnergyShift)   photonEnergyShiftx   = set.photonEnergyShift       else  photonEnergyShiftx   = photonEnergyShift     end 
    if  isnothing(mimimumPhotonEnergy)   mimimumPhotonEnergyx = set.mimimumPhotonEnergy     else  mimimumPhotonEnergyx = mimimumPhotonEnergy   end 
    if  isnothing(maximumPhotonEnergy)   maximumPhotonEnergyx = set.maximumPhotonEnergy     else  maximumPhotonEnergyx = maximumPhotonEnergy   end 
    if  isnothing(stokes)   stokesx              = set.stokes                  else  stokesx              = stokes                end 
    
    Settings( multipolesx, gaugesx, calcForStokesx, calcPhotonDmx, calcTensorsx, printBeforex, lineSelectionx,
                photonEnergyShiftx, mimimumPhotonEnergyx, maximumPhotonEnergyx, stokesx)
end


# `Base.show(io::IO, settings::PhotoExcitation.Settings)`  ... prepares a proper printout of the variable settings::PhotoExcitation.Settings.
function Base.show(io::IO, settings::PhotoExcitation.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "use-gauges:               $(settings.gauges)  ")
    println(io, "calcForStokes:            $(settings.calcForStokes)  ")
    println(io, "calcPhotonDm:             $(settings.calcPhotonDm)  ")
    println(io, "calcTensors:              $(settings.calcTensors)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:        $(settings.photonEnergyShift)  ")
    println(io, "mimimumPhotonEnergy:      $(settings.mimimumPhotonEnergy)  ")
    println(io, "maximumPhotonEnergy:      $(settings.maximumPhotonEnergy)  ")
    println(io, "stokes:                   $(settings.stokes)  ")
end


"""
`struct  PhotoExcitation.Line`  
    ... defines a type for a photo-excitation line that may include the definition of sublines and their 
        corresponding amplitudes.

    + initialLevel   ::Level                       ... initial-(state) level
    + finalLevel     ::Level                       ... final-(state) level
    + omega          ::Float64                     ... Transition frequency of this line; can be shifted w.r.t. the level energies.
    + oscStrength    ::EmProperty                  ... Absorption oscillator strength
    + crossSection   ::EmProperty                  ... Total cross section of this line.
    + staTensor      ::Array{TensorComp,1}         ... Array of statistical tensor components rho_kq
    + hasSublines    ::Bool                        ... Determines whether the individual sublines are defined in terms of their 
                                                        multipolarity, amplitude, or not; cf. PhotoEmission.Channel
    + channels       ::Array{PhotoEmission.Channel,1}  ... List of radiative (photon) channels
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    omega            ::Float64
    oscStrength      ::EmProperty
    crossSection     ::EmProperty
    staTensor        ::Array{TensorComp,1}
    hasSublines      ::Bool
    channels         ::Array{PhotoEmission.Channel,1}
end 


"""
`PhotoExcitation.Line(initialLevel::Level, finalLevel::Level, omega::Float64, crossSection::EmProperty)`  
    ... constructor an photo-excitation line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, omega::Float64, crossSection::EmProperty)
    Line(initialLevel, finalLevel, omega, EmProperty(0., 0.), crossSection, false, PhotoEmission.Channel[])
end


# `Base.show(io::IO, line::PhotoExcitation.Line)`  ... prepares a proper printout of the variable line::PhotoExcitation.Line.
function Base.show(io::IO, line::PhotoExcitation.Line) 
    println(io, "initialLevel:         $(line.initialLevel)  ")
    println(io, "finalLevel:           $(line.finalLevel)  ")
    println(io, "omega:                $(line.omega)  ")
    println(io, "oscStrength:          $(line.oscStrength)  ")
    println(io, "crossSection:         $(line.crossSection)  ")
    println(io, "staTensor:            $(line.staTensor)  ")
    println(io, "hasSublines:          $(line.hasSublines)  ")
    println(io, "channels:             $(line.channels)  ")
end


"""
`PhotoExcitation.computeAmplitudesProperties(line::PhotoExcitation.Line, grid::Radial.Grid, settings::PhotoExcitation.Settings; printout::Bool=true)`  
    ... to compute all amplitudes and properties of the given line; a line::PhotoExcitation.Line is returned for which 
        the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(line::PhotoExcitation.Line, grid::Radial.Grid, settings::PhotoExcitation.Settings; printout::Bool=true)
    newChannels = PhotoEmission.Channel[];    Ji2 = Basics.twice(line.initialLevel.J);    Jf2 = Basics.twice(line.finalLevel.J)
    for  channel  in  line.channels
        amplitude = PhotoEmission.amplitude("absorption", channel.multipole, channel.gauge, line.omega, 
                                            line.finalLevel, line.initialLevel, grid, printout=printout)
        push!( newChannels, PhotoEmission.Channel(channel.multipole, channel.gauge, amplitude) )
    end
    # Calculate the absorption oscillator strength
    oscCoulomb = oscBabushkin = 0.;   omega = line.omega;   alpha = Defaults.getDefaults("alpha")
    for  channel  in  newChannels
        wa = 8pi * alpha * line.omega / (Ji2 + 1) * (Jf2 + 1) / 2.
        #x @show  channel
        if      channel.gauge == Basics.Coulomb     oscCoulomb   = oscCoulomb    +  channel.amplitude * conj(channel.amplitude) * wa
        elseif  channel.gauge == Basics.Babushkin   oscBabushkin = oscBabushkin  +  channel.amplitude * conj(channel.amplitude) * wa
        elseif  channel.gauge == Basics.Magnetic    oscBabushkin = oscBabushkin  +  channel.amplitude * conj(channel.amplitude) * wa
                                                    oscCoulomb   = oscCoulomb    +  channel.amplitude * conj(channel.amplitude) * wa
        else    error("stop a")
        end
    end
    # Calculate the crossSection and stastistical tensors if requested
    csCoulomb = csBabushkin = 0.;   
    for  channel  in  newChannels
        if      channel.gauge == Basics.Coulomb     csCoulomb   = csCoulomb    +  channel.amplitude * conj(channel.amplitude)
        elseif  channel.gauge == Basics.Babushkin   csBabushkin = csBabushkin  +  channel.amplitude * conj(channel.amplitude)
        elseif  channel.gauge == Basics.Magnetic    csBabushkin = csBabushkin  +  channel.amplitude * conj(channel.amplitude)
                                                    csCoulomb   = csCoulomb    +  channel.amplitude * conj(channel.amplitude)
        else    error("stop b")
        end
    end
    oscFactor    = (Ji2 + 1) / (Jf2 + 1) / 2. * Defaults.getDefaults("speed of light: c")^3 / line.omega^2
    ## oscFactor    = oscFactor / (Jf2 + 1)                                                            ## Schippers (Oct. 2023)
    csFactor     = 2 * pi^3 * Defaults.getDefaults("alpha") / (line.omega * (Ji2 + 1)) ## / (Jf2 + 1)  ## Schippers (Oct. 2023)
    oscStrength  = EmProperty(oscFactor*oscCoulomb, oscFactor*oscBabushkin)
    crossSection = EmProperty(csFactor*csCoulomb, csFactor*csBabushkin)
    staTensor    = TensorComp[]
    if  settings.calcTensors    push!(staTensor, TensorComp(0, 0, 1.))      end
    
    newLine = PhotoExcitation.Line( line.initialLevel, line.finalLevel, line.omega, oscStrength, crossSection, staTensor, true, newChannels)
    return( newLine )
end


"""
`PhotoExcitation.computeCrossSection(line::PhotoExcitation.Line, stokes::ExpStokes)`  
    ... to compute the excitation cross section for the excitation of unpolarized atoms by plane-wave photons, whose polarization 
        is described by the given (experimental) Stokes parameters. A cs::EmProperty is returned.
"""
function  computeCrossSection(line::PhotoExcitation.Line, stokes::ExpStokes)
    csCoulomb = -2.;    csBabushkin = -2.
    cs = EmProperty( csCoulomb, csBabushkin )
    return( cs )
end


"""
`PhotoExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                    settings::PhotoExcitation.Settings; output=true)`  
    ... to compute the photo-excitation amplitudes and all properties as requested by the given settings. A list 
        of lines::Array{PhotoExcitation.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                        settings::PhotoExcitation.Settings; output=true)
    println("")
    printstyled("PhotoExcitation.computeLines(): The computation of the excitation cross sections, etc. starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = PhotoExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoExcitation.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = PhotoExcitation.Line[]
    for  line in lines
        newLine = PhotoExcitation.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoExcitation.displayCrossSections(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoExcitation.displayCrossSections(iostream, newLines, settings)   end
    #
    if    output    return( lines )
    else            return( nothing )
    end
end


"""
`PhotoExcitation.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                     settings::PhotoExcitation.Settings, initialLevelSelection::LevelSelection; 
                                     output::Bool=true, printout::Bool=true)`  
    ... to compute the excitation (absorption) transition amplitudes and all properties as requested by the given settings. 
        The computations and printout is adapted for larger cascade computations by including only lines with at least one channel 
        and by sending all printout to a summary file only. A list of lines::Array{PhotoExcitation.Lines} is returned.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                              settings::PhotoExcitation.Settings, initialLevelSelection::LevelSelection; output::Bool=true, printout::Bool=true) 
    # Define a common subshell list for both multiplets
    subshellList = Basics.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    lines = PhotoExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    # if  settings.printBefore    PhotoExcitation.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = PhotoExcitation.Line[]
    for  (i,line)  in  enumerate(lines)
        if  rem(i,200) == 0    println("> Excitation line $i:")   end
        if  !Basics.selectLevel(line.initialLevel, initialLevelSelection)     continue    
        ## else @show "jump to be calculated";    continue 
        end
        newLine = PhotoExcitation.computeAmplitudesProperties(line, grid, settings, printout=printout) 
        #
        # Don't add this line if it does not contribute to the decay
        wa = 0.
        for  ch in newLine.channels   wa = wa + abs(ch.amplitude)^2    end
        if   wa == 0.    continue    end
        push!( newLines, newLine)
    end
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoExcitation.displayCrossSections(iostream, newLines, settings)    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoExcitation.computeStatisticalTensor(k::Int64, q::Int64, line::PhotoExcitation.Line, stokes::ExpStokes)`  
    ... to compute the statistical tensor (component) rho_{k,q} of the final level for the excitation of unpolarized atoms by 
        plane-wave photons, whose polarization is described by the given (experimental) Stokes parameters. 
        A rho_kq::EmProperty is returned.
"""
function  computeStatisticalTensor(k::Int64, q::Int64, line::PhotoExcitation.Line, stokes::ExpStokes)
    if  k == q  ||  k == -q   rhoCoulomb = -3.;    rhoBabushkin = -3.    else   rhoCoulomb = 0.;    rhoBabushkin = 0.    end
    rho = EmProperty( rhoCoulomb, rhoBabushkin )
    return( rho )
end


"""
`PhotoExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoExcitation.Settings)`  
    ... to determine a list of PhotoExcitation.Channel for a transitions from the initial to final level and by taking 
        into account the particular settings of for this computation; an Array{PhotoExcitation.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoExcitation.Settings)
    channels = PhotoEmission.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)               push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end
    return( channels )  
end


"""
`PhotoExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoExcitation.Settings)`  
    ... to determine a list of photo-excitation Line's for transitions between the levels from the given initial- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{PhotoExcitation.Line,1} is returned. Apart from the level specification, all physical properties are set to 
        zero during the initialization process.  
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoExcitation.Settings)
    lines = PhotoExcitation.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                omega    = fLevel.energy - iLevel.energy + settings.photonEnergyShift
                if   omega == 0.  ||  omega < settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end  
                channels = PhotoExcitation.determineChannels(fLevel, iLevel, settings)
                if   length(channels) == 0   continue   end
                push!( lines, PhotoExcitation.Line(iLevel, fLevel, omega, EmProperty(0., 0.), EmProperty(0., 0.), TensorComp[], true, channels) )
            end
        end
    end
    return( lines )
end


"""
`PhotoExcitation.displayLineData(stream::IO, lines::Array{PhotoExcitation.Line,1})`  
    ... to display the calculated data, ordered by the initial levels and the photon energies involved.
        Neat tables of all initial levels and photon energies as well as all associated total cross sections are printed
        but nothing is returned otherwise.
"""
function  displayLineData(stream::IO, lines::Array{PhotoExcitation.Line,1})
    # Extract and display all initial levels by their total energy and symmetry
    energies = Float64[]
    for  line  in  lines    push!(energies, line.initialLevel.energy)  end     
    energies = unique(energies);    energies = sort(energies)
    println(stream, "\n  Initial levels, available in the given photoexcitation line data:")
    println(stream, "\n  ", TableStrings.hLine(42))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(42))
    for  energy in energies
        isNew = true
        for  line in lines
            if  line.initialLevel.energy == energy  &&  isNew
                sa  = "  ";    sym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity )
                sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * "    "
                println(stream, sa);   isNew = false
            end
        end
    end
    #
    # Extract and display photon energies for which line data are available
    omegas = Float64[]
    for  line  in  lines    push!(omegas, line.omega)   end;    
    omegas = unique(omegas);    omegas = sort(omegas)
    sa  = "    ";   
    for   omega in omegas  sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "   "      end
    #
    println(stream, "\n  Photon energies $(TableStrings.inUnits("energy")) for which given photoexcitation line data are given:")
    println(stream, "\n" * sa)
    #
    # Extract and display total cross sections, ordered by the initial levels and photon energies, for which line data are available
    println(stream, "\n  Initial level, photon energies and total cross sections for given photoexcitation line data: \n")
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Omega"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(26, "Cou--Total CS--Bab"   ; na=5);               
    sb = sb * TableStrings.center(26,TableStrings.inUnits("cross section"); na=5)
    sa = sa * TableStrings.center(20, "Shell excitations"; na=3);               
    sb = sb * TableStrings.center(20, "from  ==>  to";  na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(75))
    #
    if  length(lines) == 0    println("  >>> No photoexcitation data in lines here ... return ");    return(nothing)    end 
    #
    wLine = lines[1]
    for  energy in energies
        for  omega in omegas
            cs = Basics.EmProperty(0.)
            for  line in lines
                if  line.initialLevel.energy == energy  &&   line.omega == omega     cs = cs + line.crossSection
                    wLine = line   end
            end
            sa  = "  ";    sym = LevelSymmetry( wLine.initialLevel.J, wLine.initialLevel.parity )
            sa = sa * TableStrings.center(10, TableStrings.level(wLine.initialLevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Coulomb))   * "  " *
                        @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Babushkin)) * "     "
            println(stream, sa)
        end
    end
    
    
    return( nothing )
end


"""
`PhotoExcitation.displayLines(stream::IO, lines::Array{PhotoExcitation.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
        selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{PhotoExcitation.Line,1})
    nx = 120
    println(stream, " ")
    println(stream, "  Selected photo-excitation lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        mpGaugeList = Tuple{Basics.EmMultipole,Basics.EmGauge}[]
        for  i in 1:length(line.channels)
            push!( mpGaugeList, (line.channels[i].multipole, line.channels[i].gauge) )
        end
        sa = sa * TableStrings.multipoleGaugeTupels(50, mpGaugeList)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`PhotoExcitation.displayCrossSections(stream::IO, lines::Array{PhotoExcitation.Line,1}, settings::PhotoExcitation.Settings)`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayCrossSections(stream::IO, lines::Array{PhotoExcitation.Line,1}, settings::PhotoExcitation.Settings)
    #
    # Photoexcitation resonance strength as derived from oscillator strength
    nx = 122
    println(stream, " ")
    println(stream, "  Photoexcitation resonance strength as derived from oscillator strength:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"   ; na=3);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(24, "Cou -- f_ik -- Bab"; na=1);                        
    sb = sb * TableStrings.center(24, "     absorption   "; na=2)
    sa = sa * TableStrings.center(24, "Cou -- S -- Bab"; na=3);                        
    sb = sb * TableStrings.center(24, "   [Mb eV]     "; na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        multipoles = EmMultipole[]
        for  ch in line.channels
            multipoles = push!( multipoles, ch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
        sa = sa * @sprintf("%.4e", line.oscStrength.Coulomb)   * "   "
        sa = sa * @sprintf("%.4e", line.oscStrength.Babushkin) * "    "
        sa = sa * @sprintf("%.4e", line.oscStrength.Coulomb   * 109.761)   * "   "
        sa = sa * @sprintf("%.4e", line.oscStrength.Babushkin * 109.761)   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Photoexcitation cross sections from line strength for (completely) linearly-polarized plane-wave photons
    nx = 103
    println(stream, " ")
    println(stream, "  Photoexcitation cross sections for (completely) linearly-polarized plane-wave photons:  ... not yet correct !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"   ; na=3);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "Multipoles"; na=2);                        sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(28, "Cou -- cross section -- Bab"; na=0);       
    sb = sb * TableStrings.center(28, TableStrings.inUnits("cross section")*"       "*
                                        TableStrings.inUnits("cross section"); na=1)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        multipoles = EmMultipole[]
        for  ch in line.channels
            multipoles = push!( multipoles, ch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=5)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb))     * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin))   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    if  settings.calcForStokes
        stokes = settings.stokes
        println(stream, " ")
        println(stream, "  Photoexcitation cross sections for incident plane-wave photons with given $stokes:  ... not yet correct !!")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Energy"   ; na=3);               
        sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=2);       
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section")*"          "*
                                                TableStrings.inUnits("cross section"); na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
            crossSection = PhotoExcitation.computeCrossSection(line, stokes)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", crossSection.Babushkin))   * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTensors
        nx = 145
        stokes = settings.stokes
        println(stream, " ")
        println(stream, "  Statistical tensors rho_kq  and alignment parameters A_kq for the excitation by incident plane-wave photons")
        println(stream, "  with given $stokes:  ... not yet correct !!")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Energy"   ; na=3);               
        sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                        sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(12, "     k    q    "; na=0);                   sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(26, "Cou -- rho_kq -- Bab"; na=5);              sb = sb * TableStrings.hBlank(14)   
        sa = sa * TableStrings.center(26, "Cou --  A_kq  -- Bab"; na=2);              sb = sb * TableStrings.hBlank(14)   
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
            println(stream, sa);   sc = TableStrings.hBlank( length(sa)-1 )
            tc00   = PhotoExcitation.computeStatisticalTensor(0, 0, line, stokes)
            for  k = 2:2
                for  q = -k:k
                    tc = PhotoExcitation.computeStatisticalTensor(k, q, line, stokes)
                    if  abs(tc.Coulomb) == abs(tc.Babushkin) == 0.    continue    end
                    #
                    aCoulomb = tc.Coulomb / tc00.Coulomb;    aBabushkin = tc.Babushkin / tc00.Babushkin
                    sd = sc * TableStrings.level(k) * TableStrings.level(q) * "   "
                    sd = sd * @sprintf("%.6e", tc.Coulomb)     * "  "
                    sd = sd * @sprintf("%.6e", tc.Babushkin)   * "    "
                    sd = sd * @sprintf("%.6e", aCoulomb)       * "  "
                    sd = sd * @sprintf("%.6e", aBabushkin)     * "    "
                    println(stream, sd)
                end
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTensors   &&   settings.calcPhotonDm
        stokes = settings.stokes
        println(stream, " ")
        println(stream, "  Reduced density matrix of the fluorescence photons at selected solid angles for the prior excitation by incident plane-wave photons")
        println(stream, "  with given $stokes:")
        println(stream, " ")
        println(stream, "  Not yet implemented !!")
    end
    #
    return( nothing )
end


"""
`PhotoExcitation.estimateCrossSection(lines::Array{PhotoExcitation.Line,1}, omega::Float64, gamma::Float64, initialLevel)`  
    ... to estimate from lines the total PE cross section for any omega and for initial level. 
        The procedure assumes a Gaussian line shape for each PhotoExcitation.Line with widths gamma and distributes
        the (total cross section according to the line shape within [E_r - 10*gamma, E_r + 10*gamma].
        Other line shapes can be readily implemented; a cross section cs::EmProperty is returned.
"""
function  estimateCrossSection(lines::Array{PhotoExcitation.Line,1}, omega::Float64, gamma::Float64, initialLevel)
    function gaussianDistribution(omega::Float64, Er::Float64, gamma::Float64, gf::Basics.EmProperty)
        # Return zero if outside of interval
        if abs(Er - omega) > 10 * gamma    return(Basics.EmProperty(0.))   end
        # Account for the factor to bring 109.7617 Mb / eV first back to atomic units
        factor = 109.7617 * Defaults.convertUnits("cross section: from barn to atomic unit", 1.0e6) / 
                            Defaults.convertUnits("energy: from eV to atomic", 1.0)
        wa = (Er - omega)^2  +  gamma^2/4
        wa = factor * gamma / (2pi) / wa
        return( wa * gf )
    end
        
    cs = Basics.EmProperty(0.)
    for  line  in  lines
        if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  
            cs = cs + gaussianDistribution(omega, line.omega, gamma, line.oscStrength)   
        end  
    end
    
    
    return( cs )
end


end # module

