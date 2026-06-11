
"""
`module  JAC.InternalRecombination`  
... a submodel of JAC that contains all methods for computing internal-recombination rates between some initial 
    (Rydberg-type) and final-state multiplets.
"""
module InternalRecombination


using  Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear,
                ..Radial, ..SpinAngular, ..TableStrings

"""
`struct  InternalRecombination.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing internal-recombination lines.

    + rydbergShells         ::Array{Shell,1}     ... List of (Rydberg) shells which are to be coupled to the initial levels.
    + printBefore           ::Bool               ... True, if all energies and lines are printed before their evaluation.
    + resonanceEnergyShift  ::Float64            ... An overall energy shift for resonance energies [Hartree].
    + gamma                 ::Float64            ... An overall widths of resonances [Hartree].
    + lineSelection         ::LineSelection      ... Specifies the selected levels, if any.
    + operator              ::AbstractEeInteraction   
        ... Electron-interaction operator that is to be used for evaluating the internal-recombination amplitudes; 
            allowed values are: CoulombInteraction(), BreitInteraction(), ...
    
    The internal recombination (stabilization) is described by the following internal transitions:

        A^q+: [Ne] 3s (10s + 10p + 10d + ...)  --> [Ne] 3d (4s + 4p + 4d + ...) 
        
    due to some (plasma-type) broadening of "widths" to make the resonance condition of the dielectronic interaction
    possible. In this process, the rydbergShells describe the (10s + 10p + 10d + ...) hollow-atom Rydberg occupation which
    is coupled "configuration-wise" to the initial N-electron level, from which it "recombines" to a doubly-exited 
    (N-1) electron level of type 3d (4s + 4p + 4d + ...), or similar. These doubly-excited levels quickly "stabilizes"
    radiatively. The process can be characterized either by a "capture rate" (nominal Auger rate) or 
    "dielectronic strength", analogue to the dielectronic recombination.
"""
struct Settings  <:  AbstractProcessSettings
    rydbergShells           ::Array{Shell,1}
    printBefore             ::Bool 
    resonanceEnergyShift    ::Float64
    gamma                   ::Float64
    lineSelection           ::LineSelection 
    operator                ::AbstractEeInteraction 
end 


"""
`InternalRecombination.Settings()`  ... constructor for the default values of InternalRecombination line computations
"""
function Settings()
    Settings(Shell[], false, 0., 0.1, LineSelection(), CoulombInteraction())
end


"""
`InternalRecombination.Settings(set::InternalRecombination.Settings;`

        rydbergShells=..,           printBefore=..,             resonanceEnergyShift=..,    
        gamma=..,                   lineSelection=..,           operator=..)
                    
    ... constructor for modifying the given InternalRecombination.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::InternalRecombination.Settings;    
    rydbergShells::Union{Nothing,Array{Shell,1}}=nothing,    printBefore::Union{Nothing,Bool}=nothing, 
    resonanceEnergyShift::Union{Nothing,Float64}=nothing,    gamma::Union{Nothing,Float64}=nothing,  
    lineSelection::Union{Nothing,LineSelection}=nothing,     operator::Union{Nothing,String}=nothing)  
    
    if  isnothing(rydbergShells)         rydbergShellsx        = set.rydbergShells        else  rydbergShellsx        = rydbergShells        end 
    if  isnothing(printBefore)           printBeforex          = set.printBefore          else  printBeforex          = printBefore          end 
    if  isnothing(resonanceEnergyShift)  resonanceEnergyShiftx = set.resonanceEnergyShift else  resonanceEnergyShiftx = resonanceEnergyShift end 
    if  isnothing(gamma)                 gammax                = set.gamma                else  gammax                = gamma                end 
    if  isnothing(lineSelection)         lineSelectionx        = set.lineSelection        else  lineSelectionx        = lineSelection        end 
    if  isnothing(operator)              operatorx             = set.operator             else  operatorx             = operator             end 

    Settings( rydbergShellsx, printBeforex, resonanceEnergyShiftx, gammax, lineSelectionx, operatorx)
end


# `Base.show(io::IO, settings::InternalRecombination.Settings)`  
#   ... prepares a proper printout of the variable settings::InternalRecombination.Settings.
function Base.show(io::IO, settings::InternalRecombination.Settings) 
    println(io, "rydbergShells:            $(settings.rydbergShells)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "resonanceEnergyShift:     $(settings.resonanceEnergyShift)  ")
    println(io, "gamma:                    $(settings.gamma)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "operator:                 $(settings.operator)  ")
end


"""
`struct  InternalRecombination.Channel`   
    ... defines a type for a InternalRecombination channel to help characterize am initial Rydberg state with 
        one electron at high n.

    + subshell       ::Subshell             ... subshell of the Rydberg electron.
    + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the overall Rydberg state.
    + amplitude      ::Complex{Float64}     ... Internal-recombination amplitude associated with the given channel.
"""
struct  Channel
    subshell         ::Subshell 
    symmetry         ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  InternalRecombination.Line`  
    ... defines a type for a InternalRecombination line that may include the definition of sublines and their 
        corresponding amplitudes.

    + initialLevel   ::Level           ... initial-(state) level
    + finalLevel     ::Level           ... final-(state) level
    + deltaEnergy    ::Float64         ... Energy difference of the initial and final levels.
    + rateZ          ::Float64         ... Rate of this line as obtained from the level energies.
    + rate           ::Float64         ... Rate of this line as obtained with shifted resonance energy.
    + channels       ::Array{InternalRecombination.Channel,1}  ... List of internal-recombination channels of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    deltaEnergy      ::Float64
    rateZ            ::Float64 
    rate             ::Float64 
    channels         ::Array{InternalRecombination.Channel,1}
end 


"""
`InternalRecombination.Line()`  
    ... constructor for an empty InternalRecombination line.
"""
function Line()
    Line(Level(), Level(), 0., 0., 0., InternalRecombination.Channel[])
end


# `Base.show(io::IO, line::InternalRecombination.Line)`  ... prepares a proper printout of the variable line::InternalRecombination.Line.
function Base.show(io::IO, line::InternalRecombination.Line) 
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "deltaEnergy:            $(line.deltaEnergy)  ")
    println(io, "rateZ:                  $(line.rateZ)  ")
    println(io, "rate:                   $(line.rate)  ")
    println(io, "channels:               $(line.channels)  ")
end


"""
`InternalRecombination.amplitude(operator::AbstractEeInteraction, finalLevel::Level, rydbergLevel::Level, grid::Radial.Grid; 
                                    display::Bool=true, printout::Bool=false)
    ... to compute the internal-recombination amplitude  
    
                <alpha_f J_f || V^(e-e) || (alpha_i J_i, n kappa) J_t>
            
        A value::ComplexF64 is returned and the amplitude value printed to screen, if display=true.
"""
function amplitude(operator::AbstractEeInteraction, finalLevel::Level, rydbergLevel::Level, grid::Radial.Grid; 
                    display::Bool=true, printout::Bool=false)
    #
    if  finalLevel.basis.subshells != rydbergLevel.basis.subshells  error("stop a")     end
    if  finalLevel.J               != rydbergLevel.J                error("stop b")     end
    fLevel = finalLevel;    iLevel = rydbergLevel
    
    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity)
    
    if  printout   printstyled("Compute internal-stabilization amplitude for the transition [$(rydbergLevel.index)-$(finalLevel.index)] ... ", 
                                color=:light_green)    end
    amplitude = ComplexF64(0.)
    #
    for  r = 1:nf
        for  s = 1:ni
            #
            #    <alpha_f J_f || V^(e-e) || (alpha_i J_i, n kappa) J_t>
            Vee       = ManyElectron.matrixElement_Vee(operator, fLevel.basis, r, iLevel.basis, s, grid)
            amplitude = amplitude + fLevel.mc[r] * Vee * iLevel.mc[s]
        end
    end
    if  printout   printstyled("done. \n", color=:light_green)    end
    
    if  display  
        println("    < level=$(finalLevel.index) [J=$symf)] ||" * " IR^($operator) ||" * " $(rydbergLevel.index) [$symi)] >  = $amplitude  ")
    end
    
    return( amplitude )
end


"""
`InternalRecombination.computeAmplitudesProperties(line::InternalRecombination.Line, nm::Nuclear.Model, 
                                                    rydbergOrbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, 
                                                    settings::InternalRecombination.Settings; printout::Bool=true)` 
    ... to compute all amplitudes and properties of the given line; a line::InternalRecombination.Line is returned 
        for which the amplitudes and properties are now evaluated.
"""
function computeAmplitudesProperties(line::InternalRecombination.Line, nm::Nuclear.Model, rydbergOrbitals::Dict{Subshell, Orbital}, 
                                        grid::Radial.Grid, settings::InternalRecombination.Settings; printout::Bool=true) 
    newChannels = InternalRecombination.Channel[];   rateZ = 0.;   rate = 0.;   gHalf = settings.gamma / 2.
    # Define a common subshell list for both multiplets
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    
    for channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
        newfLevel = Basics.generateLevelWithExtraSubshell(channel.subshell, newfLevel)
        # Generate Rydberg orbital in the field of the initial core
        newrLevel = Basics.generateLevelWithExtraElectron(rydbergOrbitals[channel.subshell], channel.symmetry, newiLevel)
        amplitude = InternalRecombination.amplitude(settings.operator, newfLevel, newrLevel, grid, printout=printout)
        # Apply energy dominator + gamma
        factorZ   = gHalf / ( (line.initialLevel.energy - line.finalLevel.energy)^2 + gHalf^2 )
        factor    = gHalf / ( line.deltaEnergy^2 + gHalf^2 )
        rateZ     = rateZ + conj(amplitude) * amplitude * factorZ
        rate      = rate  + conj(amplitude) * amplitude * factor
        push!( newChannels, InternalRecombination.Channel(channel.subshell, channel.symmetry, amplitude) )
    end
    rateZ = 2pi* rateZ;   rate = 2pi* rate
    newLine   = InternalRecombination.Line(line.initialLevel, line.finalLevel, line.deltaEnergy, rateZ, rate, newChannels)
    
    return( newLine )
end


"""
`InternalRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::InternalRecombination.Settings; output=true, printout::Bool=true)`  
    ... to compute the internal-recombination transition amplitudes and rates as requested by the given settings. A list of 
        lines::Array{InternalRecombination.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::InternalRecombination.Settings; output=true, printout::Bool=true)
    println("")
    printstyled("InternalRecombination.computeLines(): The computation of internal-recombination rates starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    # Generate orbitals for all rydberg-subshells
    rydbergSubshells = Basics.generateSubshellList(settings.rydbergShells)
    meanPot          = Basics.computePotential(Basics.DFSField(1.0), grid, initialMultiplet.levels[1].basis)
    primitives       = Bsplines.generatePrimitives(grid)
    rydbergOrbitals  = Bsplines.generateOrbitals(rydbergSubshells, meanPot, nm, primitives; printout=true)
    #
    lines = InternalRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    InternalRecombination.displayLines(stdout, lines)    end  
    # Calculate all amplitudes and requested properties
    newLines = InternalRecombination.Line[]
    for  line in lines
        newLine = InternalRecombination.computeAmplitudesProperties(line, nm, rydbergOrbitals, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    InternalRecombination.displayRates(stdout, newLines, settings)
    InternalRecombination.displayTotalRates(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   InternalRecombination.displayRates(iostream, newLines, settings)    
                        InternalRecombination.displayTotalRates(iostream, newLines, settings)        end
    #
    if    output    return( lines )
    else            return( nothing )
    end
end


"""
`InternalRecombination.determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalRecombination.Settings)`  
    ... to determine a list of internal-recombination Channel for a transitions from the initial to final level and by taking 
        into account the particular settings of for this computation; an Array{InternalRecombination.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalRecombination.Settings)
    channels  = InternalRecombination.Channel[];   
    symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    subshells = Basics.generateSubshellList(settings.rydbergShells)
    for  subsh in subshells
        tSymmetries = AngularMomentum.allowedTotalSymmetries(symi, subsh.kappa)
        for  symt in tSymmetries
            if  symt != symf      continue    end
            push!(channels, InternalRecombination.Channel(subsh, symt, Complex(0.)) )
        end
    end
    return( channels )  
end


"""
`InternalRecombination.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalRecombination.Settings)`  
    ... to determine a list of InternalRecombination.Line's for transitions between levels from the initial- and final-state 
        multiplets, and by taking into account the particular selections and settings for this computation; 
        an Array{InternalRecombination.Line,1} is returned. Apart from the level specification, all physical properties are 
        set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalRecombination.Settings)
    lines = InternalRecombination.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                dEnergy = iLevel.energy - fLevel.energy   + settings.resonanceEnergyShift
                channels = InternalRecombination.determineChannels(fLevel, iLevel, settings) 
                push!( lines, InternalRecombination.Line(iLevel, fLevel, dEnergy, 0., 0., channels) )
            end
        end
    end
    return( lines )
end


"""
`InternalRecombination.displayLines(stream::IO, lines::Array{InternalRecombination.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{InternalRecombination.Line,1})
    nx = 150
    println(stream, " ")
    println(stream, "  Selected internal-recombination lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Delta Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(37, "List of subshells and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(37, "subshell (total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.7e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * " "
        sb = "      " * @sprintf("%.7e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy))
        sa = sa * sb[end-16:end] * "    "
        subshellSymmetryList = Tuple{Subshell,LevelSymmetry}[]
        for  channel in line.channels
            push!( subshellSymmetryList, (channel.subshell, channel.symmetry) )
        end
        sa = sa * TableStrings.subshellSymmetryTupels(80, subshellSymmetryList)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`InternalRecombination.displayRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned 
        otherwise.
"""
function  displayRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)
    nx = 104
    println(stream, " ")
    println(stream, "  Internal recombination rates without and with the resonance shift: \n")
    println(stream, "  Resonance shift [Hartree]:    $(settings.resonanceEnergyShift) ")
    println(stream, "  Gamma [Hartree]:              $(settings.gamma)  \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy  " ; na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Delta energy"   ; na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(14, "Rate (dE=0.)"   ; na=0);       
    sb = sb * TableStrings.center(14, TableStrings.inUnits("rate");  na=1)
    sa = sa * TableStrings.center(15, "Rate"           ; na=2);                           
    sb = sb * TableStrings.center(14, TableStrings.inUnits("rate");  na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=2)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "  "
        sx = "     " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.deltaEnergy))
        sa = sa * sx[end-13:end]                                                                            * "   "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.rateZ))                  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.rate))                   * "  "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`InternalRecombination.displayTotalRates(stream::IO, lines::Array{InternalRecombination.Line,1}, 
                                            settings::InternalRecombination.Settings)`  
    ... to list all total rates for the initial levels. A neat table is printed but nothing is returned 
        otherwise.
"""
function  displayTotalRates(stream::IO, lines::Array{InternalRecombination.Line,1}, settings::InternalRecombination.Settings)
    nx = 106
    println(stream, " ")
    println(stream, "  Total internal-recombination rates for initial levels without and with the resonance shift: \n")
    println(stream, "  Resonance shift [Hartree]:    $(settings.resonanceEnergyShift) ")
    println(stream, "  Gamma [Hartree]:              $(settings.gamma)  \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center( 8, "i-level"; na=2);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center( 6, "J^P";     na=2);                         sb = sb * TableStrings.hBlank( 8)
    sa = sa * TableStrings.center(12, "Energy";  na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(22, "Total IDE rate (dE=0.)"; na=3);       
    sb = sb * TableStrings.center(22, TableStrings.inUnits("rate");  na=3)
    sa = sa * TableStrings.center(16, "Total IDE rate"; na=2);                           
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate");  na=2)
    sa = sa * TableStrings.center(22, "A_IDE / A_radiative"; na=3);       
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    used = falses(length(lines))
    for  (i, line) in enumerate(lines)
        if  used[i]     continue    end
        totalZ = 0.;    total = 0.
        for  j = 1:length(lines)
            if  line.initialLevel.index == lines[j].initialLevel.index      used[j] = true
                totalZ = totalZ + lines[j].rateZ
                total  = total  + lines[j].rate 
            end
        end
        #
        totalZ = 2.2e3 * totalZ;    total = 2.2e3 * total
        isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
        sa = "      " * string(line.initialLevel.index);   sb = "            " * string(isym)
        sa = sa[end-6:end] * sb[end-10:end] * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "       "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", totalZ))                      * "          "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", total))                       * "        "
        sa = sa * @sprintf("%.6e", rand(1)[1] )                       * "      "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module

