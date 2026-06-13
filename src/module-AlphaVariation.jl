
"""
`module  JAC.AlphaVariation`  
	... a submodel of JAC that contains all methods for computing alpha-variation parameters for some levels.
"""
module AlphaVariation


using  Printf, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings


"""
`struct  AlphaVariation.Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing alpha-variation parameters.

    + calcK                    ::Bool             ... True if the enhancement factor need to be calculated, and false otherwise.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    calcK                      ::Bool
    printBefore                ::Bool 
    levelSelection             ::LevelSelection
end 


"""
`AlphaVariation.Settings()`  ... constructor for an `empty` instance of AlphaVariation.Settings for the computation of alpha variation parameters.
"""
function Settings()
    Settings(false, false, LevelSelection() )
end


"""
`AlphaVariation.Settings(set::AlphaVariation.Settings;`

        calcK=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::AlphaVariation.Settings.
"""
function Settings(set::AlphaVariation.Settings;
        calcK::Union{Nothing,Bool}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,          levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  calcK          == nothing   calcKx          = set.calcK          else   calcKx          = calcK          end
    if  printBefore    == nothing   printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  levelSelection == nothing   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( calcKx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::AlphaVariation.Settings)`  ... prepares a proper printout of the variable settings::AlphaVariation.Settings.
function Base.show(io::IO, settings::AlphaVariation.Settings) 
    println(io, "calcK:                    $(settings.calcK)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end



"""
`struct  AlphaVariation.Outcome`  
    ... defines a type to keep the outcome of a alpha-variation computation, such as the K enhancement factor as well other results.

    + level                     ::Level              ... Atomic level to which the outcome refers to.
    + K                         ::Float64            ... K enhancement parameter
"""
struct Outcome 
    level                       ::Level 
    K                           ::Float64
end 


"""
`AlphaVariation.Outcome()`  ... constructor for an `empty` instance of AlphaVariation..Outcome for the computation of alpha-variation properties.
"""
function Outcome()
    Outcome(Level(), 0.)
end


# `Base.show(io::IO, outcome::AlphaVariation.Outcome)`  ... prepares a proper printout of the variable outcome::AlphaVariation.Outcome.
function Base.show(io::IO, outcome::AlphaVariation.Outcome) 
    println(io, "level:                   $(outcome.level)  ")
    println(io, "K:                       $(outcome.K)  ")
end


"""
`AlphaVariation.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Nuclear.Model, grid::Radial.Grid, settings::AlphaVariation.Settings; output=true)` 
    ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by the given settings. 
        The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Nuclear.Model, grid::Radial.Grid, settings::AlphaVariation.Settings; output=true)
    println("")
    printstyled("AlphaVariation.computeOutcomes(): The computation of the alpha-variation parameters starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    outcomes = AlphaVariation.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    AlphaVariation.displayOutcomes(stdout,outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = AlphaVariation.Outcome[]
    for  outcome in outcomes
        ## newOutcome = AlphaVariation.computeAmplitudesProperties(outcome, nm, grid, settings) 
        newOutcome = AlphaVariation.Outcome(outcome.level, 2.0)
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    AlphaVariation.displayResults(stdout, newOutcomes)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    AlphaVariation.displayResults(iostream, newOutcomes)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`AlphaVariation.determineOutcomes(multiplet::Multiplet, settings::AlphaVariation.Settings)`  
    ... to determine a list of Outcomes's for the computation of the alpha-variation parameters for the given multiplet. 
        It takes into account the particular selections and settings. An Array{AlphaVariation.Outcome,1} is returned. 
        Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::AlphaVariation.Settings) 
    outcomes = AlphaVariation.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, AlphaVariation.Outcome(level, 0.) )
        end
    end
    return( outcomes )
end


"""
`AlphaVariation.displayOutcomes(stream::IO, outcomes::Array{AlphaVariation.Outcome,1})`  ... to display a list of levels that have been selected 
        for the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(stream::IO, outcomes::Array{AlphaVariation.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected AlphaVariation levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`AlphaVariation.displayResults(stream::IO, outcomes::Array{AlphaVariation.Outcome,1})`  
    ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{AlphaVariation.Outcome,1})
    nx = 64
    println(stream, " ")
    println(stream, "  Alpha variation parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "K";     na=4)              
    sb = sb * TableStrings.center(14, "    " ; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = 1.0
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "    "
        sa = sa * @sprintf("%.8e", outcome.K)                                               * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module
