
"""
`module  JAC.MultipolePolarizibility`  
... a submodel of JAC that contains all methods for computing dynamic multipole polarizibilities, including 
    (dynanmic) electric- & magnetic-dipole polarizibilities and others.
"""
module MultipolePolarizibility


using Printf, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

"""
`struct  MultipolePolarizibility.Amplitude`  
    ... defines a type for an individual multipole term as defined by the multipole and frequency.

    + multipole  ::EmMultipole         ... Multipole of the dynamic polarizibility contribution.
    + omega      ::Float64             ... (dynamic) frequency of the polarizibility.
    + value      ::Complex{Float64}    ... Value of this multipole amplitude.

"""
struct  Amplitude
    multipole    ::EmMultipole
    omega        ::Float64
    value        ::Complex{Float64}
end 


"""
`MultipolePolarizibility.Amplitude()`  ... constructor for an `empty` instance of MultipolePolarizibility.Amplitude.
"""
function Amplitude()
    Amplitude(E1, 0., 0.)
end


# `Base.show(io::IO, amplitude::MultipolePolarizibility.Amplitude)`  ... prepares a proper printout of the variable amplitude.
function Base.show(io::IO, amplitude::MultipolePolarizibility.Amplitude) 
    println(io, "multipole :       $(amplitude.multipole)  ")
    println(io, "omega:            $(amplitude.omega)  ")
    println(io, "value:            $(amplitude.value)  ")
end


"""
`struct  MultipolePolarizibility.Outcome`  ... defines a type to keep the outcome of a dynamic-polarizibility computation.

    + level        ::Level                                       ... Atomic level to which the outcome refers to.
    + amplitudes   ::Array{MultipolePolarizibility.Amplitude,1}  ... multipole (and frequency-dependent) contribution to the polarizibility.
"""
struct Outcome 
    level          ::Level  
    amplitudes     ::Array{MultipolePolarizibility.Amplitude,1}
end 


"""
`MultipolePolarizibility.Outcome()`  
    ... constructor for an `empty` instance of MultipolePolarizibility.Outcome for the computation of dynamic 
        polarizibilities.
"""
function Outcome()
    Outcome(Level(), MultipolePolarizibility.Amplitude[])
end


# `Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)` 
#		 ... prepares a proper printout of the variable outcome::MultipolePolarizibility.Outcome.
function Base.show(io::IO, outcome::MultipolePolarizibility.Outcome) 
    println(io, "level:                     $(outcome.level)  ")
    println(io, "amplitudes:                $(outcome.amplitudes)  ")
end


"""
`struct  MultipolePolarizibility.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing multipolar polarizibilities.

    + multipoles                ::Array{EmMultipole,1}  ... Multipoles to be considered for the polarity.
    + nLower                    ::Int64                 ... Lower and upper indices in the configurations for the summation 
    + nUpper                    ::Int64                     over the intermediate levels.
    + omegas                    ::Array{Float64,1}      ... List of omegas (energies) of the dynamic polarizibility.
    + printBefore               ::Bool                  ... True if a list of selected levels is printed before the 
                                                            actual computations start. 
    + levelSelection            ::LevelSelection        ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    multipoles                  ::Array{EmMultipole,1}
    nLower                      ::Int64  
    nUpper                      ::Int64
    omegas                      ::Array{Float64,1}
    printBefore                 ::Bool 
    levelSelection              ::LevelSelection
end 


"""
`MultipolePolarizibility.Settings()`  
    ... constructor for an `empty` instance of MultipolePolarizibility.Settings for the computation dynamic polarizibilities.
"""
function Settings()
    Settings(EmMultipole[], 0, 0, Float64[], false, LevelSelection() )
end


"""
`MultipolePolarizibility.Settings(set::MultipolePolarizibility.Settings;`

        multipoles=.., nLower=.., nUpper=.., omegas=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::MultipolePolarizibility.Settings.
"""
function Settings(set::MultipolePolarizibility.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
        nLower::Union{Nothing,Int64}=nothing,              nUpper::Union{Nothing,Int64}=nothing,
        omegas::Union{Nothing,Array{Float64,1}}=nothing,   printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(multipoles)       multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(nLower)           nLowerx         = set.nLower         else   nLowerx         = nLower         end
    if  isnothing(nUpper)           nUpperx         = set.nUpper         else   nUpperx         = nUpper         end
    if  isnothing(omegas)           omegasx         = set.omegas         else   omegasx         = omegas         end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( multipolesx, nLowerx, nUpperx, omegasx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings:MultipolePolarizibility.Settings)`  
#		... prepares a proper printout of the variable settings::MultipolePolarizibility.Settings.
function Base.show(io::IO, settings::MultipolePolarizibility.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "nLower:                   $(settings.nLower)  ")
    println(io, "nUpper:                   $(settings.nUpper)  ")
    println(io, "omegas:                   $(settings.omegas)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`MultipolePolarizibility.computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, grid::Radial.Grid, 
                                                        settings::MultipolePolarizibility.Settings) 
    ... to compute all amplitudes and properties for a given level; an outcome::MultipolePolarizibility.Outcome is 
        returned for which the amplitudes and properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
    ## global JAC_counter
    ## return( newOutcome )
end


"""
`MultipolePolarizibility.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                                settings::MultipolePolarizibility.Settings; output=true)` 
    ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by 
        the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
            settings::MultipolePolarizibility.Settings; output=true)
    println("")
    printstyled("MultipolePolarizibility.computeOutcomes(): The computation of multipole polarizibilities starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    outcomes = MultipolePolarizibility.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    MultipolePolarizibility.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = MultipolePolarizibility.Outcome[]
    for  outcome in outcomes
        ## newOutcome = MultipolePolarizibility.computeAmplitudesProperties(outcome, nm, grid, settings) 
        newOutcome = MultipolePolarizibility.Outcome(outcome.level, [ MultipolePolarizibility.Amplitude(E1, 1., 2.)] )
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    MultipolePolarizibility.displayResults(stdout, newOutcomes)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    MultipolePolarizibility.displayResults(iostream, newOutcomes)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`MultipolePolarizibility.determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings)`  
    ... to determine a list of Outcomes's for the computation of the multipole polarizibilities for the given multiplet. It takes 
        into account the particular selections and settings. An Array{MultipolePolarizibility.Outcome,1} is returned. Apart from the 
        level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::MultipolePolarizibility.Settings) 
    outcomes = MultipolePolarizibility.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, MultipolePolarizibility.Outcome(level, MultipolePolarizibility.Amplitude[]) )
        end
    end
    return( outcomes )
end


"""
`MultipolePolarizibility.displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{MultipolePolarizibility.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected MultipolePolarizibility levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`MultipolePolarizibility.displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})`  
    ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})
    nx = 64
    println(stream, " ")
    println(stream, "  Multipole polarizibilities and amplitudes:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "xxx";     na=4)              
    sb = sb * TableStrings.center(14, "    " ; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = 1.0
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module


