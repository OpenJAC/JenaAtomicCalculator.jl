
"""
`module  JAC.AlphaVariation`
	... a submodel of JAC that contains all methods for computing alpha-variation parameters for some levels.
"""
module AlphaVariation


using  Printf, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..SelfConsistent, ..TableStrings


"""
`struct  AlphaVariation.Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing alpha-variation parameters.

    + calcQ                    ::Bool             ... True if the q (and Q) sensitivity coefficients need to be calculated, and false otherwise.
    + variationX               ::Float64          ... The finite-difference step x = (alpha/alpha_0)^2 - 1 used to shift alpha in the two
                                                        auxiliary SCF/CI computations; the physical alpha corresponds to x = 0.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start.
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcQ                      ::Bool
    variationX                 ::Float64
    printBefore                ::Bool
    levelSelection             ::LevelSelection
end


"""
`AlphaVariation.Settings()`  ... constructor for an `empty` instance of AlphaVariation.Settings for the computation of alpha variation parameters.
"""
function Settings()
    Settings(false, 0.125, false, LevelSelection() )
end


"""
`AlphaVariation.Settings(set::AlphaVariation.Settings;`

        calcQ=.., variationX=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::AlphaVariation.Settings.
"""
function Settings(set::AlphaVariation.Settings;
        calcQ::Union{Nothing,Bool}=nothing,                variationX::Union{Nothing,Float64}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,          levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcQ)            calcQx          = set.calcQ          else   calcQx          = calcQ          end
    if  isnothing(variationX)       variationXx     = set.variationX     else   variationXx     = variationX     end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( calcQx, variationXx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::AlphaVariation.Settings)`  ... prepares a proper printout of the variable settings::AlphaVariation.Settings.
function Base.show(io::IO, settings::AlphaVariation.Settings)
    println(io, "calcQ:                    $(settings.calcQ)  ")
    println(io, "variationX:               $(settings.variationX)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end



"""
`struct  AlphaVariation.Outcome`
    ... defines a type to keep the outcome of a alpha-variation computation, such as the q (and Q) sensitivity
        coefficients as well as other results.

    + level                     ::Level              ... Atomic level to which the outcome refers to.
    + omega                     ::Float64            ... Level energy at the physical value of alpha.
    + q                         ::Float64            ... Sensitivity coefficient, q = d omega / dx  at x = 0, in the same
                                                            (atomic-unit) energy scale as omega.
    + Q                         ::Float64            ... Dimensionless sensitivity coefficient, Q = q / (omega - omega_ground),
                                                            referred to the lowest level of the SAME multiplet; 0. for the
                                                            ground level itself, where this ratio is undefined.
"""
struct Outcome
    level                       ::Level
    omega                       ::Float64
    q                           ::Float64
    Q                           ::Float64
end


"""
`AlphaVariation.Outcome()`  ... constructor for an `empty` instance of AlphaVariation.Outcome for the computation of alpha-variation properties.
"""
function Outcome()
    Outcome(Level(), 0., 0., 0.)
end


# `Base.show(io::IO, outcome::AlphaVariation.Outcome)`  ... prepares a proper printout of the variable outcome::AlphaVariation.Outcome.
function Base.show(io::IO, outcome::AlphaVariation.Outcome)
    println(io, "level:                   $(outcome.level)  ")
    println(io, "omega:                   $(outcome.omega)  ")
    println(io, "q:                       $(outcome.q)  ")
    println(io, "Q:                       $(outcome.Q)  ")
end


"""
`AlphaVariation.matchLevelByOverlap(target::Level, candidates::Array{Level,1})`
    ... identifies, among candidates, the level whose mixing-coefficient vector has the LARGEST overlap with
        target's -- the state that IS the same physical state, rather than the one that merely sits at the
        same POSITION in an energy-sorted list. This is what makes alpha-variation safe across a level
        crossing: for these ions q can reach 1e5-1e6 cm^-1 against level spacings of only a few thousand, so a
        step as small as a few percent can re-order the levels, and a positional match would then silently
        difference two DIFFERENT physical states (challenge S6, closed by this function). target and every
        level in candidates must share one CSF basis, in the same order -- true whenever both were built from
        the same configs and asfSettings, since CSF generation is purely combinatorial and does not depend on
        alpha. RAISES if the best overlap is not a clear majority (< 0.5): either this level has no clean
        counterpart at the shifted alpha (the step is too large for the level spacing here, and variationX
        needs to shrink), or candidates was not built from the same basis as target.
        A tuple (level::Level, overlap::Float64) is returned.
"""
function matchLevelByOverlap(target::Level, candidates::Array{Level,1})
    pool = [ lev  for lev in candidates  if  lev.J == target.J  &&  lev.parity == target.parity ]
    if  isempty(pool)
        error("AlphaVariation.matchLevelByOverlap(): no candidate shares target's symmetry $(target.J)^" *
              "$(target.parity) -- candidates was not built from the same configs as target.")
    end
    if  length(pool[1].mc) != length(target.mc)
        error("AlphaVariation.matchLevelByOverlap(): mixing-coefficient vectors differ in length " *
              "($(length(pool[1].mc)) against $(length(target.mc))) -- candidates and target were not built " *
              "from the same CSF basis.")
    end
    overlaps = [ sum(target.mc[r] * lev.mc[r]  for r in eachindex(target.mc))^2   for lev in pool ]
    iBest    = argmax(overlaps)
    if  overlaps[iBest] < 0.5
        error("AlphaVariation.matchLevelByOverlap(): no level of symmetry $(target.J)^$(target.parity) carries " *
              "a majority overlap with the target level (the best of $(length(pool)) candidates is " *
              "$(overlaps[iBest])) -- this level has no clean counterpart at the shifted alpha, most likely " *
              "because variationX is too large relative to the level spacing here; a smaller step is needed.")
    end

    return( pool[iBest], overlaps[iBest] )
end


"""
`AlphaVariation.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, configs::Array{Configuration,1},
                                asfSettings::AsfSettings, settings::AlphaVariation.Settings; output=true)`
    ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified by
        the given settings. The multiplet is taken as the x = 0 (physical-alpha) reference; if settings.calcQ, two further
        SCF/CI computations are performed for the same configs & asfSettings at alpha shifted by x = +settings.variationX
        and x = -settings.variationX (with x = (alpha/alpha_0)^2 - 1), and q = d omega/dx is obtained from the simple
        symmetric finite-difference formula q = (omega(+x) - omega(-x)) / (2x). Levels are matched across the three
        multiplets by OVERLAP of their mixing coefficients (AlphaVariation.matchLevelByOverlap), not by position,
        so a level crossing under the shift does not silently difference two different physical states; this
        requires configs & asfSettings to be exactly the ones that produced multiplet -- a minimal consistency
        check on the electron number is done for this reason, in addition to matchLevelByOverlap's own checks.
        The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, configs::Array{Configuration,1},
                          asfSettings::AsfSettings, settings::AlphaVariation.Settings; output=true)
    println("")
    printstyled("AlphaVariation.computeOutcomes(): The computation of the alpha-variation parameters starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    # A minimal (elementary) consistency check between multiplet and configs -- just the electron number; a full CSF-by-CSF
    # comparison is not attempted here, since multiplet is trusted to be the direct performSCF(configs, ...) result.
    if  !isempty(multiplet.levels)  &&  configs[1].NoElectrons != multiplet.levels[1].basis.NoElectrons
        error("AlphaVariation.computeOutcomes(): multiplet and configs are inconsistent -- different electron numbers.")
    end
    #
    outcomes = AlphaVariation.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    AlphaVariation.displayOutcomes(stdout, outcomes)    end
    #
    if  settings.calcQ
        x      = settings.variationX
        alpha0 = Defaults.getDefaults("alpha")
        #
        printstyled(">> Compute the two auxiliary multiplets at alpha = alpha_0 * sqrt(1 +- x), x = $x ... \n", color=:light_green)
        Defaults.setDefaults("fine structure constant: alpha", alpha0 * sqrt(1.0 + x))
        multipletPlus  = SelfConsistent.performSCF(configs, nm, grid, asfSettings; printout=false)
        Defaults.setDefaults("fine structure constant: alpha", alpha0 * sqrt(1.0 - x))
        multipletMinus = SelfConsistent.performSCF(configs, nm, grid, asfSettings; printout=false)
        Defaults.setDefaults("fine structure constant: alpha", alpha0)
        #
        if  length(multipletPlus.levels) != length(multiplet.levels)  ||  length(multipletMinus.levels) != length(multiplet.levels)
            error("AlphaVariation.computeOutcomes(): the auxiliary multiplets do not have the same number of levels as multiplet.")
        end
        #
        # The lowest-energy level of the FULL multiplet is taken as the (Kozlov et al.) reference level; its own
        # q is needed below even if it is not itself among the selected outcomes. EACH LEVEL IS MATCHED BY
        # OVERLAP, NOT BY INDEX -- see matchLevelByOverlap's own docstring for why a positional match is unsafe
        # on exactly the ions this module exists to serve.
        groundLevel  = multiplet.levels[ argmin( [level.energy for level in multiplet.levels] ) ]
        groundPlus,  _ = AlphaVariation.matchLevelByOverlap(groundLevel, multipletPlus.levels)
        groundMinus, _ = AlphaVariation.matchLevelByOverlap(groundLevel, multipletMinus.levels)
        qGround      = (groundPlus.energy - groundMinus.energy) / (2x)
        #
        newOutcomes = AlphaVariation.Outcome[]
        for  outcome in outcomes
            levelPlus,  _  = AlphaVariation.matchLevelByOverlap(outcome.level, multipletPlus.levels)
            levelMinus, _  = AlphaVariation.matchLevelByOverlap(outcome.level, multipletMinus.levels)
            q           = (levelPlus.energy - levelMinus.energy) / (2x)
            excitation  = outcome.level.energy - groundLevel.energy
            Q           = excitation == 0.   ?   0.   :   (q - qGround) / excitation
            newOutcome  = AlphaVariation.Outcome(outcome.level, outcome.level.energy, q, Q)
            push!( newOutcomes, newOutcome)
        end
    else
        newOutcomes = AlphaVariation.Outcome[]
        for  outcome in outcomes
            newOutcome = AlphaVariation.Outcome(outcome.level, outcome.level.energy, 0., 0.)
            push!( newOutcomes, newOutcome)
        end
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
            push!( outcomes, AlphaVariation.Outcome(level, level.energy, 0., 0.) )
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
    ... to display the energies as well as the q- and Q-sensitivity coefficients for the selected levels. A neat table
        is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{AlphaVariation.Outcome,1})
    nx = 88
    println(stream, " ")
    println(stream, "  Alpha variation parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "omega"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "q";     na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Q";     na=4)
    sb = sb * TableStrings.center(14, "    " ; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.omega))            * "    "
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.q))                * "    "
        sa = sa * @sprintf("%.8e", outcome.Q)                                                               * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module
