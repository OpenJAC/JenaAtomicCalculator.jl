
"""
`module  JAC.TwoElectronOnePhoton`  
... a submodel of JAC that contains all methods for computing two-electron-one-photon (TEOP) transition rates between 
    some initial and final-state multiplets. 
"""
module TwoElectronOnePhoton


using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, 
                ..Radial, ..TableStrings


"""
`struct  TwoElectronOnePhoton.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing radiative lines.

    + multipoles              ::Array{EmMultipoles}     ... Specifies the (radiat. field) multipoles to be included.
    + gauges                  ::Array{UseGauge}         ... Gauges to be included into the computations.
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before comput.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + photonEnergyShift       ::Float64                 ... An overall energy shift for all photon energies.
    + eeInteraction           ::AbstractEeInteraction   ... Type of the electron-electron interaction in the second-order treatment.
    + gMultiplet              ::Multiplet               ... Mean-field multiplet of intermediate levels in the computations.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge}
    printBefore               ::Bool 
    lineSelection             ::LineSelection 
    photonEnergyShift         ::Float64
    eeInteraction             ::AbstractEeInteraction
    gMultiplet                ::Multiplet
end 


"""
`TwoElectronOnePhoton.Settings()`  ... constructor for the default values of TEOP line computations
"""
function Settings()
    Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], false, LineSelection(), 0., CoulombInteraction(), Multiplet() )
end


"""
`TwoElectronOnePhoton.Settings(set::TwoElectronOnePhoton.Settings;`

        multipoles::=..,        gauges=..,                printBefore=..,         lineSelection=..,         
        photonEnergyShift=..,   eeInteraction=..,         gMultiplet=..) 
                    
    ... constructor for modifying the given TwoElectronOnePhoton.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::TwoElectronOnePhoton.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    gauges::Union{Nothing,Array{UseGauge}}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                   lineSelection::Union{Nothing,LineSelection}=nothing, 
    photonEnergyShift::Union{Nothing,Float64}=nothing,          eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,
    gMultiplet::Union{Nothing,Multiplet}=nothing)
    
    if  isnothing(multipoles)            multipolesx          = set.multipoles          else  multipolesx          = multipoles         end 
    if  isnothing(gauges)                gaugesx              = set.gauges              else  gaugesx              = gauges             end 
    if  isnothing(printBefore)           printBeforex         = set.printBefore         else  printBeforex         = printBefore        end 
    if  isnothing(lineSelection)         lineSelectionx       = set.lineSelection       else  lineSelectionx       = lineSelection      end 
    if  isnothing(photonEnergyShift)     photonEnergyShiftx   = set.photonEnergyShift   else  photonEnergyShiftx   = photonEnergyShift  end 
    if  isnothing(eeInteraction)         eeInteractionx       = set.eeInteraction       else  eeInteractionx       = eeInteraction      end 
    if  isnothing(gMultiplet)            gMultipletx          = set.gMultiplet          else  gMultipletx          = gMultiplet         end 
        
    Settings( multipolesx, gaugesx, printBeforex, lineSelectionx, photonEnergyShiftx, eeInteractionx, gMultipletx)
end


# `Base.show(io::IO, settings::TwoElectronOnePhoton.Settings)`  
#    ... prepares a proper printout of the variable settings::TwoElectronOnePhoton.Settings.
function Base.show(io::IO, settings::TwoElectronOnePhoton.Settings) 
    println(io, "multipoles:             $(settings.multipoles)  ")
    println(io, "gauges:                 $(settings.gauges)  ")
    println(io, "printBefore:            $(settings.printBefore)  ")
    println(io, "lineSelection:          $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:      $(settings.photonEnergyShift)  ")
    println(io, "gMultiplet:             $(settings.gMultiplet)  ")
end


"""
`struct  TwoElectronOnePhoton.Channel`  
    ... defines a type for a TEOP channel that specifies the multipole, gauge and amplitude. The total angular momentum and 
        parity of the intermediate Green function levels coincides with either the intial or final-level symmetry and, therefore,
        need not to be marked separately. However, these Green function levels formally occur in the computation of the amplitude.

    + multipole         ::EmMultipole        ... Multipole of the photon emission/absorption.
    + gauge             ::EmGauge            ... Gauge for dealing with the (coupled) radiation field.
    + amplitude         ::Complex{Float64}   ... Amplitude of this multiple channel.
"""
struct Channel 
    multipole           ::EmMultipole
    gauge               ::EmGauge
    amplitude           ::Complex{Float64}
end 


# `Base.show(io::IO, channel::TwoElectronOnePhoton.Channel)`  
#    ... prepares a proper printout of the variable channel::TwoElectronOnePhoton.Channel.
function Base.show(io::IO, channel::TwoElectronOnePhoton.Channel) 
    print(io, "TwoElectronOnePhoton.Channel($(channel.multipole), $(channel.gauge), amp = $(channel.amplitude)) ") 
end


"""
`struct  TwoElectronOnePhoton.Line`  
    ... defines a type for a TEOP line that may include the definition of sublines and their corresponding amplitudes.

    + initialLevel   ::Level               ... initial-(state) level
    + finalLevel     ::Level               ... final-(state) level
    + omega          ::Float64             ... Transition frequency of this line; can be shifted w.r.t. the level energies.
    + teopRate       ::EmProperty          ... Total TEOP rate of this line.
    + channels       ::Array{TwoElectronOnePhoton.Channel,1}  ... List of TEOP (photon) channels
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    omega            ::Float64
    teopRate         ::EmProperty
    channels         ::Array{TwoElectronOnePhoton.Channel,1}
end 


"""
`TwoElectronOnePhoton.Line()`  
    ... constructor an empty TEOP line.
"""
function Line()
    Line(Level(), Level(), 0., EmProperty(0., 0.), TwoElectronOnePhoton.Channel[])
end


# `Base.show(io::IO, line::TwoElectronOnePhoton.Line)`  ... prepares a proper printout of the variable line::TwoElectronOnePhoton.Line.
function Base.show(io::IO, line::TwoElectronOnePhoton.Line) 
    println(io, "initialLevel:         $(line.initialLevel)  ")
    println(io, "finalLevel:           $(line.finalLevel)  ")
    println(io, "omega:                $(line.omega)  ")
    println(io, "teopRate:             $(line.teopRate)  ")
    println(io, "channels:             $(line.channels)  ")
end




"""
`TwoElectronOnePhoton.amplitude(::Emission, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                                gMultiplet::Multiplet, grid::Radial.Grid; display::Bool=false, printout::Bool=false)`
    ... to compute the TEOP emission amplitude

                <alpha_f J_f || O^(Mp, emission) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
            +   <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, emission) || alpha_i J_i>

        emission amplitude for the interaction with the photon field of multipolarity Mp and for the given transition energy
        and gauge. A value::ComplexF64 is returned. The amplitude value is printed to screen if display=true.
"""
function amplitude(::Emission, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                    gMultiplet::Multiplet, grid::Radial.Grid; display::Bool=false, printout::Bool=false)
    #
    # Always ensure the same subshell list for all initial, intermediate and final levels. The gMultiplet must take
    # part in this merge: a physically-motivated Green-function space generally introduces orbitals that occur in
    # neither the initial nor the final configurations (e.g. the 2p of the near-degenerate 1s2p4s partner for the
    # Li-like O5+ 1s2s5p --> 1s^2 4s transition), and merging only initial and final makes those cases unusable.
    subshells  = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
    subshells  = Basics.merge(subshells, gMultiplet.levels[1].basis.subshells)
    iLevel     = Level(initialLevel, subshells)
    fLevel     = Level(finalLevel, subshells)
    nMultiplet = Multiplet(gMultiplet, subshells)
    
    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity);    eni = iLevel.energy
    nn = length(nMultiplet.levels[1].basis.csfs)
    efinal = eni - omega
    #
    # Second-order perturbation theory requires the intermediate (Green-function) level n to be genuinely
    # DIFFERENT from the initial or final level itself; a formally coincident n -- which can occur if the
    # gMultiplet's reference configurations overlap with the initial- or final-state configurations -- would
    # give a spurious, ill-defined 1/(eni-enn) or 1/(efinal-enn) singularity and must be excluded from the sum.
    selfTolerance = 1.0e-8

    if  printout   printstyled("Compute TEOP $(Mp) amplitude for the transition [$(iLevel.index)-$(fLevel.index)] ... ",
                                color=:light_green)    end
    amplitude = ComplexF64(0.)
    #
    for  r = 1:nf
        symr = LevelSymmetry(fLevel.basis.csfs[r].J, fLevel.basis.csfs[r].parity);      if  symr != symf    continue    end
        for  s = 1:ni
            syms = LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity);  if  syms != symi    continue    end
            for  nLevel in nMultiplet.levels
                symn = LevelSymmetry(nLevel.J, nLevel.parity);    enn = nLevel.energy
                #
                #   Compute <alpha_f J_f || O^(Mp, kind) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
                #   The intermediate Green-function level is, in general, multi-configurational; the two matrix
                #   elements therefore run over two INDEPENDENT CSF indices t, tp of that same level.
                if  symn == symi  &&  abs(enn - eni) > selfTolerance
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        OMp = ManyElectron.matrixElement_Mab(Mp, gauge, omega, fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            Vee       = ManyElectron.matrixElement_Vee(CoulombInteraction(), nLevel.basis, tp, iLevel.basis, s, grid)
                            amplitude = amplitude + fLevel.mc[r] * OMp * nLevel.mc[t] * nLevel.mc[tp] * Vee * iLevel.mc[s] / (eni - enn)
                        end
                    end
                end
                #
                #   Compute <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, kind) || alpha_i J_i>
                if  symn == symf  &&  abs(enn - efinal) > selfTolerance
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        Vee = ManyElectron.matrixElement_Vee(CoulombInteraction(), fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            OMp       = ManyElectron.matrixElement_Mab(Mp, gauge, omega, nLevel.basis, tp, iLevel.basis, s, grid)
                            amplitude = amplitude + fLevel.mc[r] * Vee * nLevel.mc[t] * nLevel.mc[tp] * OMp * iLevel.mc[s] / (eni - omega - enn)
                        end
                    end
                end
            end
        end
    end
    if  printout   printstyled("done. \n", color=:light_green)    end
    
    if  display  
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " TEOP^($Mp, emission) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end
    
    return( amplitude )
end
        

"""
`TwoElectronOnePhoton.computeAmplitudesProperties(line::TwoElectronOnePhoton.Line, grid::Radial.Grid, 
                                                    settings::TwoElectronOnePhoton.Settings; printout::Bool=true)`  
    ... to compute all amplitudes and properties of the given line; a line::TwoElectronOnePhoton.Line is returned for which the 
        amplitudes and properties are now evaluated.
"""
function  computeAmplitudesProperties(line::TwoElectronOnePhoton.Line, grid::Radial.Grid, settings::TwoElectronOnePhoton.Settings; 
                                        printout::Bool=true)
    newChannels = TwoElectronOnePhoton.Channel[];    rateC = 0.;    rateB = 0.
    for channel in line.channels
        #
        amplitude = TwoElectronOnePhoton.amplitude(Emission(), channel.multipole, channel.gauge, line.omega,
                                                    line.finalLevel, line.initialLevel, settings.gMultiplet, grid, printout=printout)
        #
        push!( newChannels, TwoElectronOnePhoton.Channel( channel.multipole, channel.gauge, amplitude) )
        #
        if       channel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    #     
    # Calculate the TEOP rate if requested 
    wa = 8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1) 
    photonrate  = EmProperty(wa * rateC, wa * rateB)    
    line = TwoElectronOnePhoton.Line( line.initialLevel, line.finalLevel, line.omega, photonrate, newChannels)
    return( line )
end


"""
`TwoElectronOnePhoton.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::TwoElectronOnePhoton.Settings; output=true)`  
    ... to compute the TEOP transition amplitudes and rates as requested by the given settings. A list of 
        lines::Array{TwoElectronOnePhoton.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::TwoElectronOnePhoton.Settings; output=true) 
    # Define a common subshell list for both multiplets
    subshellList = Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    #
    println("")
    printstyled("TwoElectronOnePhoton.computeLines(): The computation of the TEOP amplitudes and rates starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = TwoElectronOnePhoton.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    TwoElectronOnePhoton.displayLines(lines)    end
    # Calculate all amplitudes and requested properties
    newLines = TwoElectronOnePhoton.Line[]
    for  line in lines
        newLine = TwoElectronOnePhoton.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    TwoElectronOnePhoton.displayRates(stdout, newLines, settings)
    TwoElectronOnePhoton.displayLifetimes(stdout, newLines, settings)
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   TwoElectronOnePhoton.displayRates(iostream, newLines, settings)
                        TwoElectronOnePhoton.displayLifetimes(iostream, newLines, settings)
    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`TwoElectronOnePhoton.determineChannels(finalLevel::Level, initialLevel::Level, settings::TwoElectronOnePhoton.Settings)`  
    ... to determine a list of TwoElectronOnePhoton.Channel for a transitions from the initial to final level and by taking into 
        account the particular settings of for this computation. The definition of the channels remains the same as in single-photon
        emission, whereas the amplitudes requires special care. Since the electron-electron interaction is of zero rank, either
        the initial or final symmetry gives rises to an additional Green-function levels, together with the associated
        V^(e-e) interaction. An Array{TwoElectronOnePhoton.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::TwoElectronOnePhoton.Settings)
    channels = TwoElectronOnePhoton.Channel[];   c0   = Complex{Float64}(0.)
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb    push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Coulomb,   c0) )
                elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin  push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Babushkin, c0) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)                push!(channels, TwoElectronOnePhoton.Channel(mp, Basics.Magnetic,  c0) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end
    return( channels )  
end


"""
`TwoElectronOnePhoton.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::TwoElectronOnePhoton.Settings)`  
    ... to determine a list of TwoElectronOnePhoton Line's for transitions between the levels from the given initial- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{TwoElectronOnePhoton.Line,1} is returned. Apart from the level specification, all physical properties are set 
        to zero during the initialization process.  
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::TwoElectronOnePhoton.Settings)
    lines = TwoElectronOnePhoton.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                omega = iLevel.energy - fLevel.energy   + settings.photonEnergyShift
                if  omega <= 0.   continue   end    # no emission line for a final level above the initial one
                #
                channels = TwoElectronOnePhoton.determineChannels(fLevel, iLevel, settings)
                if   length(channels) == 0   continue   end
                push!( lines, TwoElectronOnePhoton.Line(iLevel, fLevel, omega, EmProperty(0., 0.), channels) )
            end
        end
    end
    return( lines )
end


"""
`TwoElectronOnePhoton.displayLines(lines::Array{TwoElectronOnePhoton.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
        selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(lines::Array{TwoElectronOnePhoton.Line,1})
    nx = 95
    println(" ")
    println("  Selected radiative lines:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * TableStrings.hBlank(34)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        mpGaugeList = Tuple{EmMultipole,EmGauge}[]
        for  i in 1:length(line.channels)
            push!( mpGaugeList, (line.channels[i].multipole, line.channels[i].gauge) )
        end
        sa = sa * TableStrings.multipoleGaugeTupels(50, mpGaugeList)
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    println(" ")
    #
    return( nothing )
end


"""
`TwoElectronOnePhoton.displayLifetimes(stream::IO, lines::Array{TwoElectronOnePhoton.Line,1}, settings::TwoElectronOnePhoton.Settings)`  
    ... to list all lifetimes as derived from the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayLifetimes(stream::IO, lines::Array{TwoElectronOnePhoton.Line,1}, settings::TwoElectronOnePhoton.Settings)
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
        waCoulomb = 0.;    waBabushkin = 0.
        for  i = 1:length(lines)
            if   lines[i].initialLevel.index == ii    
                waCoulomb   = waCoulomb   + lines[i].teopRate.Coulomb
                waBabushkin = waBabushkin + lines[i].teopRate.Babushkin
            end
        end
        push!(irates, EmProperty( waCoulomb, waBabushkin) )
    end
    
    nx = 105
    println(stream, " ")
    println(stream, "  Two-electron-one-photon lifetimes (as derived from these computations):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Used Gauge"    ; na=6);                     sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(26, "Lifetime"; na=3);       
    sb = sb * TableStrings.center(26, "[a.u.]"*"          "*TableStrings.inUnits("time"); na=5)
    sa = sa * TableStrings.center(12, "Decay widths"; na=4);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #    
    for  ii = 1:length(ilevels)
        sa = istr[ii]
        sa = sa * "Coulomb          " * @sprintf("%.6e",              1.0/irates[ii].Coulomb)     * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Coulomb) )   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Coulomb) )
        println(stream, sa)
        sa = repeat(" ", length(istr[ii]) )
        sa = sa * "Babushkin        " * @sprintf("%.6e",              1.0/irates[ii].Babushkin)   * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Babushkin) ) * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Babushkin) )
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`TwoElectronOnePhoton.displayRates(stream::IO, lines::Array{TwoElectronOnePhoton.Line,1}, settings::TwoElectronOnePhoton.Settings)`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayRates(stream::IO, lines::Array{TwoElectronOnePhoton.Line,1}, settings::TwoElectronOnePhoton.Settings)
    nx = 161
    println(stream, " ")
    println(stream, "  TEOP (Einstein) coefficients, transition rates and oscillator strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center( 9, "Multipole"; na=0);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=3);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"          "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(11, "Osc. strength"    ; na=3);                 sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(12, "Decay widths"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(13, "Line strength"; na=4);       
    sb = sb * TableStrings.center(12, "[a.u.]"       ; na=4)
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
            chRate =  8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1) * (abs(ch.amplitude)^2) 
            sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToEinsteinA(),    line, chRate)) * "  "
            sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToEinsteinB(),    line, chRate)) * "    "
            sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToOscillatorGf(),           line, chRate)) * "    "
            sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToDecayWidth(),   line, chRate)) * "    "
            if  ch.multipole == E1
                    sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToLineStrengthS(),     line, chRate)) * "    "
            else    sa = sa * "  --  " 
            end
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module

