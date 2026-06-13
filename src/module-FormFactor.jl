
"""
`module  JAC.FormFactor`  
... a submodel of JAC that contains all methods for computing atomic form factors for some level(s).
"""
module FormFactor


using  Printf, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals, ..TableStrings


"""
`struct  FormFactor.Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing alpha-variation parameters.

    + qList                    ::Array{Float64,1} ... List of q-values in [a.u.]
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the  
                                                        actual computations start. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    qList                      ::Array{Float64,1} 
    printBefore                ::Bool  
    levelSelection             ::LevelSelection
end 


"""
`FormFactor.Settings()`  ... constructor for an `empty` instance of FormFactor.Settings for the computation of atomic form factors.
"""
function Settings()
    Settings(Float64[0., 0.1, 1.0, 10.], false, LevelSelection() )
end


"""
`FormFactor.Settings(set::FormFactor.Settings;`

        qList=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::FormFactor.Settings.
"""
function Settings(set::FormFactor.Settings;
        qList::Union{Nothing,Array{Float64,1}}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,          levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  qList          == nothing   qListx          = set.qList          else   qListx          = qList          end
    if  printBefore    == nothing   printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  levelSelection == nothing   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( qListx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::FormFactor.Settings)`  ... prepares a proper printout of the variable settings::FormFactor.Settings.
function Base.show(io::IO, settings::FormFactor.Settings) 
    println(io, "qList:                    $(settings.qList)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end



"""
`struct  FormFactor.Outcome`  ... defines a type to keep the outcome of a form-factor computation, such as the standard form factor
                                    as well other results.

    + level                     ::Level              ... Atomic level to which the outcome refers to.
    + qValues                   ::Array{Float64,1}   ... momentum transfer |q|.
    + standardFs                ::Array{Float64,1}   ... standard atomic form factor F(q) for an (assumed) spherical charge distr.
    + modifiedFs                ::Array{Float64,1}   ... modified atomic form factor F(q) for such a charge distribution.
"""
struct Outcome 
    level                       ::Level 
    qValues                     ::Array{Float64,1}
    standardFs                  ::Array{Float64,1}
    modifiedFs                  ::Array{Float64,1} 
end 


"""
`FormFactor.Outcome()`  
    ... constructor for an `empty` instance of FormFactor.Outcome for the computation of atomic form factors.
"""
function Outcome()
    Outcome(Level(), Float64[], Float64[], Float64[])
end


# `Base.show(io::IO, outcome::FormFactor.Outcome)`  ... prepares a proper printout of the variable outcome::FormFactor.Outcome.
function Base.show(io::IO, outcome::FormFactor.Outcome) 
    println(io, "level:                   $(outcome.level)  ")
    println(io, "qValues:                 $(outcome.qValues)  ")
    println(io, "standardFs:              $(outcome.standardFs)  ")
    println(io, "modifiedFs:              $(outcome.modifiedFs)  ")
end


"""
`FormFactor.amplitude(q::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)`  
    ... to compute the momentum transfer amplitude  <(alpha_f J_f, kappa) J_i || T^(1) || alpha_i J_i>  for the given 
        final and initial level. A value::ComplexF64 is returned.
"""
function amplitude(q::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)
    nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
    printstyled("Compute (single-electron) momentum transfer matrix of dimension $nf x $ni in the final- and initial-state bases " *
                "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
    matrix = zeros(ComplexF64, nf, ni)
    #
    printstyled("done. \n", color=:light_green)
    printstyled("\nWarning -- FormFactor.amplitude():: Not yet implemented.", color=:cyan)
    #
    amplitude = 4.0 + 5.0*im
    #
    if  display
        println("   (Single-electron) Momentum transfer amplitude:  "                                                 *
                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || T^(1) ($q a.u.) ||"   *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end
    
    return( amplitude )
end


"""
`FormFactor.standardF(q::Float64, level::Level, grid::Radial.Grid)`  
        ... to compute the (real) Fourier transform of the (spherically-symmetric) charge distribution. A value::Float64 is returned.
"""
function standardF(q::Float64, level::Level, grid::Radial.Grid)
    wa = zeros( grid.NoPoints )
    # Compute the total density
    for  sh in level.basis.subshells
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])
        orb  = level.basis.orbitals[sh]
        nrho = length(orb.P)
        for    i = 1:nrho   wa[i] = wa[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Compute the full integrant; the factor 4pi * r^2 is already in the density
    if     q == 0.   
    else             for    i = 2:grid.NoPoints   wa[i] = wa[i] / grid.r[i] * sin(q*grid.r[i]) / q    end
    end
    sF = RadialIntegrals.V0(wa, grid.NoPoints, grid::Radial.Grid)
    
    return( sF )
end


"""
`FormFactor.modifiedF(q::Float64, level::Level, grid::Radial.Grid)`  
    ... to compute the (real) Fourier transform of the (spherically-symmetric) charge distribution but including a correction
        due to the local (DFS) potential, relative to the rest mass of the electron. A value::Float64 is returned.
"""
function modifiedF(q::Float64, level::Level, grid::Radial.Grid)
    wa  = zeros( grid.NoPoints )
    wc  = Defaults.getDefaults("speed of light: c")
    pot = Basics.computePotential(Basics.DFSField(), grid, level)
    # Compute the total density
    for  sh in level.basis.subshells
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])
        orb  = level.basis.orbitals[sh]
        ei   = abs(orb.energy)
        nrho = length(orb.P)
        #
        for    i = 2:nrho 
            wx = wc^2 / (wc^2 + ei + pot.Zr[i]/grid.r[i])
            if     q == 0.   wa[i] = wa[i] + occ * (orb.P[i]^2 + orb.Q[i]^2) * wx
            else             wa[i] = wa[i] + occ * (orb.P[i]^2 + orb.Q[i]^2) / grid.r[i] * sin(q*grid.r[i]) / q * wx
            end
        end
    end
    mF = RadialIntegrals.V0(wa, grid.NoPoints, grid::Radial.Grid)
    
    return( mF )
end


"""
`FormFactor.computeAmplitudesProperties(outcome::FormFactor.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::FormFactor.Settings)`  
    ... to compute all form factors, etc. for the given outcome; a newOutcome::FormFactor.Outcome is returned.
"""
function computeAmplitudesProperties(outcome::FormFactor.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::FormFactor.Settings)
    standardFs = Float64[];    modifiedFs = Float64[]
    for  i = 1:length(outcome.qValues)
        sF = FormFactor.standardF(outcome.qValues[i]::Float64, outcome.level::Level, grid::Radial.Grid)
        mF = FormFactor.modifiedF(outcome.qValues[i]::Float64, outcome.level::Level, grid::Radial.Grid)
        push!( standardFs, sF);    push!( modifiedFs, mF)
    end
    newOutcome = FormFactor.Outcome(outcome.level, outcome.qValues, standardFs, modifiedFs)
    return( newOutcome )
end


"""
`FormFactor.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::FormFactor.Settings; output=true)` 
    ... to compute (as selected) the alpha-variation parameters for the levels of the given multiplet and as specified 
        by the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::FormFactor.Settings; output=true)
    println("")
    printstyled("FormFactor.computeOutcomes(): The computation of the atomic form factors starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    outcomes = FormFactor.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    FormFactor.displayOutcomes(stdout, outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = FormFactor.Outcome[]
    for  outcome in outcomes
        newOutcome = FormFactor.computeAmplitudesProperties(outcome, nm, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    FormFactor.displayResults(stdout, newOutcomes)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    FormFactor.displayResults(iostream, newOutcomes)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`FormFactor.determineOutcomes(multiplet::Multiplet, settings::FormFactor.Settings)`  
    ... to determine a list of Outcomes's for the computation of the alpha-variation parameters for the given 
        multiplet. It takes into account the particular selections and settings. An Array{FormFactor.Outcome,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the 
        initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::FormFactor.Settings) 
    outcomes = FormFactor.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, FormFactor.Outcome(level, settings.qList, Float64[], Float64[] ) )
        end
    end
    return( outcomes )
end


"""
`FormFactor.displayOutcomes(stream::IO, outcomes::Array{FormFactor.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(stream::IO, outcomes::Array{FormFactor.Outcome,1})
    nx = 86
    println(stream, " ")
    println(stream, "  Selected FormFactor levels:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * "  q-values  ...";                                                      sb = sb * "   [a.u.]"
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        if  length(outcome.qValues) >= 1   sa = sa * @sprintf("%.2e", outcome.qValues[1] ) * ",  "   end
        if  length(outcome.qValues) >= 2   sa = sa * @sprintf("%.2e", outcome.qValues[2] ) * ",  "   end
        if  length(outcome.qValues) >= 3   sa = sa * @sprintf("%.2e", outcome.qValues[3] ) * ",  "   end
        if  length(outcome.qValues) >  3   sa = sa * "...    "   else   sa = sa * "    "             end
        println(stream,  sa ) 
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`FormFactor.displayResults(stream::IO, outcomes::Array{FormFactor.Outcome,1})`  
    ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but 
        nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{FormFactor.Outcome,1})
    nx = 88
    println(stream, " ")
    println(stream, "  Standard and modified atomic form factors [note F^(standard) (0) = N_e]:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=3)              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(14, "q [a.u.]";       na=1) 
    sa = sa * TableStrings.center(14, "standard F";     na=1) 
    sa = sa * TableStrings.center(14, "modified F";     na=4) 
    sb = sb * TableStrings.center(14, "    " ; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        sa = sa * @sprintf("%.5e", outcome.qValues[1])                                       * "    "
        sa = sa * @sprintf("%.5e", outcome.standardFs[1])                                    * "    "
        sa = sa * @sprintf("%.5e", outcome.modifiedFs[1])                                    * "    "
        println(stream, sa )
        for  i = 2:length(outcome.qValues)
            sa = repeat(" ", 47)
            sa = sa * @sprintf("%.5e", outcome.qValues[i])                                   * "    "
            sa = sa * @sprintf("%.5e", outcome.standardFs[i])                                * "    "
            sa = sa * @sprintf("%.5e", outcome.modifiedFs[i])                                * "    "
            println(stream, sa )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module
