
"""
`module  JAC.ElectronCapture`  
... a submodel of JAC that contains all methods for computing electron-capture properties between some initial and final-state 
    multiplets.
"""
module ElectronCapture


using  Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, 
                ..Radial, ..TableStrings

"""
`struct  ElectronCapture.Settings`  ... defines a type for the details and parameters of computing electron-capture lines.

    + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
    + minCaptureEnergy        ::Float64                      ... Minimum energy of free (Auger) electrons to be captured.
    + maxCaptureEnergy        ::Float64                      ... Maximum energy of free (Auger) electrons to be captured.
    + maxKappa                ::Int64                        ... Maximum kappa value of partial waves to be included.
    + operator                ::AbstractEeInteraction        ... Auger/capture operator that is to be used for evaluating the capture amplitudes: 
                                                                    allowed values are: CoulombInteraction(), BreitInteraction(), ...
"""
struct Settings
    printBefore               ::Bool 
    lineSelection             ::LineSelection 
    minCaptureEnergy          ::Float64
    maxCaptureEnergy          ::Float64
    maxKappa                  ::Int64
    operator                  ::AbstractEeInteraction 
end 


"""
`ElectronCapture.Settings()`  ... constructor for the default values of ElectronCapture line computations
"""
function Settings()
    Settings(false, LineSelection(), 0., 1.0e5, 2, CoulombInteraction())
end


"""
`ElectronCapture.Settings(set::ElectronCapture.Settings;`

        printBefore=..,       minCaptureEnergy=..,       maxCaptureEnergy=..,       maxKappa=..,       operator=..)
                    
    ... constructor for modifying the given ElectronCapture.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ElectronCapture.Settings;    
    printBefore::Union{Nothing,Bool}=nothing,               lineSelection::Union{Nothing,LineSelection}=nothing, 
    minCaptureEnergy::Union{Nothing,Float64}=nothing,       maxCaptureEnergy::Union{Nothing,Float64}=nothing,
    maxKappa::Union{Nothing,Int64}=nothing,                 operator::Union{Nothing,String}=nothing)  
    
    if  printBefore      == nothing   printBeforex      = set.printBefore       else  printBeforex      = printBefore        end 
    if  lineSelection    == nothing   lineSelectionx    = set.lineSelection     else  lineSelectionx    = lineSelection      end 
    if  minCaptureEnergy == nothing   minCaptureEnergyx = set.minCaptureEnergy  else  minCaptureEnergyx = minCaptureEnergy   end 
    if  maxCaptureEnergy == nothing   maxCaptureEnergyx = set.maxCaptureEnergy  else  maxCaptureEnergyx = maxCaptureEnergy   end 
    if  maxKappa         == nothing   maxKappax         = set.maxKappa          else  maxKappax         = maxKappa           end 
    if  operator         == nothing   operatorx         = set.operator          else  operatorx         = operator           end 

    Settings( printBeforex, lineSelectionx, minCaptureEnergyx, maxCaptureEnergyx, maxKappax, operatorx)
end


# `Base.show(io::IO, settings::ElectronCapture.Settings)`  ... prepares a proper printout of the variable settings::ElectronCapture.Settings.
function Base.show(io::IO, settings::ElectronCapture.Settings) 
    println(io, "printBefore:                   $(settings.printBefore)  ")
    println(io, "lineSelection:                 $(settings.lineSelection)  ")
    println(io, "minCaptureEnergy:              $(settings.minCaptureEnergy)  ")
    println(io, "maxCaptureEnergy:              $(settings.maxCaptureEnergy)  ")
    println(io, "maxKappa:                      $(settings.maxKappa)  ")
    println(io, "operator:                      $(settings.operator)  ")
end


"""
`struct  Channel`   
    ... defines a type for a ElectronCapture channel to help characterize a scattering (continuum) state of many 
        electron-states with a single free electron.

    + kappa          ::Int64                ... partial-wave of the free electron
    + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + phase          ::Float64              ... phase of the partial wave
    + amplitude      ::Complex{Float64}     ... electron-capture amplitude associated with the given channel.
"""
struct  Channel
    kappa            ::Int64
    symmetry         ::LevelSymmetry
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  Line`  
    ... defines a type for a ElectronCapture line that may include the definition of sublines and their 
        corresponding amplitudes.

    + initialLevel   ::Level           ... initial-(state) level
    + finalLevel     ::Level           ... final-(state) level
    + electronEnergy ::Float64         ... Energy of the (incoming free) electron.
    + totalRate      ::Float64         ... Total rate of this line.
    + channels       ::Array{ElectronCapture.Channel,1}  ... List of ElectronCapture channels of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    totalRate        ::Float64
    channels         ::Array{ElectronCapture.Channel,1}
end 


"""
`ElectronCapture.Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)`  
    ... constructor for an ElectronCapture line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)
    Line(initialLevel, finalLevel, 0., totalRate, false, ElectronCapture.Channel[])
end


# `Base.show(io::IO, line::ElectronCapture.Line)`  ... prepares a proper printout of the variable line::ElectronCapture.Line.
function Base.show(io::IO, line::ElectronCapture.Line) 
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "electronEnergy:         $(line.electronEnergy)  ")
    println(io, "totalRate:              $(line.totalRate)  ")
    println(io, "channels:               $(line.channels)  ")
end


"""
`ElectronCapture.amplitude(kind::String, channel::ElectronCapture.Channel, finalLevel::Level, continuumLevel::Level,
                            grid::Radial.Grid; printout::Bool=true)`  
    ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit(), CoulombGaunt()   ElectronCapture amplitude 
        <alpha_f J_f || O^(capture, kind) || (alpha_i J_i, kappa) J_f> = <(alpha_i J_i, kappa) J_f || O^(Auger, kind) || alpha_f J_f>^*
        due to the interelectronic interaction for the given final and initial level. A value::ComplexF64 is returned.
"""
function amplitude(kind::String, channel::ElectronCapture.Channel, finalLevel::Level, continuumLevel::Level, grid::Radial.Grid; printout::Bool=true)
    aChannel  = AutoIonization.Channel(channel.kappa, channel.symmetry, channel.phase, 0.)
    amplitude = conj(AutoIonization.amplitude(kind, aChannel, continuumLevel::Level, finalLevel::Level, grid; printout=printout))
    
    return( amplitude )
end


"""
`ElectronCapture.computeAmplitudesProperties(line::ElectronCapture.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                settings::ElectronCapture.Settings; printout::Bool=true)` 
    ... to compute all amplitudes and properties of the given line; a line::ElectronCapture.Line is returned for which the amplitudes 
        and properties are now evaluated.
"""
function computeAmplitudesProperties(line::ElectronCapture.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                        settings::ElectronCapture.Settings; printout::Bool=true) 
    newChannels = ElectronCapture.Channel[];   contSettings = Continuum.Settings(false, nrContinuum);   rate = 0.
    for channel in line.channels
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, line.finalLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, newfLevel.basis.subshells)
        newfLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newfLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newiLevel, nm, grid, contSettings)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newiLevel)
        newChannel = ElectronCapture.Channel(channel.kappa, channel.symmetry, phase, 0.)
        amplitude  = ElectronCapture.amplitude(settings.operator, newChannel, newcLevel, newfLevel, grid, printout=printout)
        rate       = rate + conj(amplitude) * amplitude
        push!( newChannels, ElectronCapture.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude) )
    end
    totalRate = 2pi* rate
    newLine   = ElectronCapture.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, true, newChannels)
    
    return( newLine )
end


"""
`ElectronCapture.displayRates(stream::IO, lines::Array{ElectronCapture.Line,1}, settings::ElectronCapture.Settings)`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned 
        otherwise.
"""
function  displayRates(stream::IO, lines::Array{ElectronCapture.Line,1}, settings::ElectronCapture.Settings)
    nx = 96
    println(stream, " ")
    println(stream, "  Auger rates: \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Electron energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(14, "Capture rate"; na=2);       
    sb = sb * TableStrings.center(14, TableStrings.inUnits("rate"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.totalRate))              * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
