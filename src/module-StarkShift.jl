
"""
`module  JAC.StarkShift`
... a submodel of JAC that contains all methods for computing the static electric-dipole/-quadrupole
    moments (Theta^(1), Theta^(2)) and the static Stark shift of fine-structure levels in a uniform
    electric field.

    The Stark shift of a level |J,M> in a static field E is the direct physical use of the SAME
    static E1 scalar/tensor polarizability alpha_0(J)/alpha_2(J) that module-MultipolePolarizibility.jl
    computes: DeltaE(J,M) = -(1/2) alpha_0(J) E^2 - (1/2) alpha_2(J) E^2 [3M^2-J(J+1)] / [J(2J-1)]
    (standard second-order-perturbation-theory form; the Angel & Sandars alpha_2 normalization is
    defined precisely so this simple form applies -- see van Leeuwen & Hogervorst, Z. Phys. A 316, 149
    (1984), Eq. (5), the same paper module-MultipolePolarizibility.jl's docstring already references).
    Rather than re-deriving/re-validating this sum-over-states here, StarkShift.computeStarkShifts
    calls MultipolePolarizibility.computeScalarTensorPolarizability directly, on a caller-supplied,
    independently-generated StarkShift.Settings.gMultiplet (mirrors
    MultipolePolarizibility.Settings.gMultiplet/MultiPhotonDeExcitation.Settings.gMultiplet exactly --
    see examples/example-Cl.jl), and only adds the field-strength/M-dependent energy-shift conversion
    on top. Both Coulomb and Babushkin gauge are reported throughout (EmProperty), not calibrated onto
    a common scale, following the same convention already established for MultipolePolarizibility.
"""
module StarkShift


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..MultipoleMoment,
              ..MultipolePolarizibility, ..Nuclear, ..Radial, ..TableStrings


"""
`struct  StarkShift.SublevelJ`  ... defines a type to specify a Stark-shifted sublevel with well-defined J.

    + M                      ::AngularM64        ... M_J-value
    + energy                 ::EmProperty        ... Stark energy shift of this sublevel, in both the
                                                       Coulomb and Babushkin gauge (see module docstring).
"""
struct SublevelJ
    M                        ::AngularM64
    energy                   ::EmProperty
end


# `Base.show(io::IO, Jsublevel::StarkShift.SublevelJ)`  ... prepares a proper printout of the variable Jsublevel::StarkShift.SublevelJ.
function Base.show(io::IO, Jsublevel::StarkShift.SublevelJ) 
    println(io, "Sublevel [M=$(Jsublevel.M); energy = $(Jsublevel.energy)]")
end


"""
`struct  StarkShift.Outcome`
    ... defines a type to keep the Stark shift and electric-quadrupole-moment parameters of a fine-structure level.

    + Jlevel                 ::Level             ... Fine-structure levels to which the results refer to.
    + Theta1                 ::Float64           ... Electric-dipole moment of the atom in this level.
    + Theta2                 ::Float64           ... Electric-quadrupole moment of the atom in this level.
    + alpha0                 ::EmProperty        ... Static scalar (rank 0) E1 polarizability alpha_0(J), in
                                                       a.u., that the Stark shift below was computed from
                                                       (see module docstring); EmProperty(0.) if calcStarkshifts
                                                       was not requested.
    + alpha2                 ::EmProperty        ... Static tensor (rank 2) E1 polarizability alpha_2(J), in
                                                       a.u.; identically 0. for J<1 (see module docstring).
    + Jsublevels              ::Array{StarkShift.SublevelJ,1}
        ... List of the Stark-shifted fine-structure sublevels and data with well-defined J-value
"""
struct Outcome
    Jlevel                   ::Level
    Theta1                   ::Float64
    Theta2                   ::Float64
    alpha0                   ::EmProperty
    alpha2                   ::EmProperty
    Jsublevels                ::Array{StarkShift.SublevelJ,1}
end


"""
`StarkShift.Outcome()`  ... constructor for an `empty` instance of StarkShift.Outcome.
"""
function Outcome()
    Outcome(Level(), 0., 0., EmProperty(0.), EmProperty(0.), StarkShift.SublevelJ[])
end


# `Base.show(io::IO, outcome::StarkShift.Outcome)`  ... prepares a proper printout of the variable StarkShift.Outcome.
function Base.show(io::IO, outcome::StarkShift.Outcome)
    println(io, "Jlevel:           $(outcome.Jlevel)  ")
    println(io, "Theta1:           $(outcome.Theta1)  ")
    println(io, "Theta2:           $(outcome.Theta2)  ")
    println(io, "alpha0:           $(outcome.alpha0)  ")
    println(io, "alpha2:           $(outcome.alpha2)  ")
    println(io, "Jsublevels:       $(outcome.Jsublevels)  ")
end


"""
`struct  StarkShift.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of computing the Stark shift and electric-quadrupole-moments
        of fine-structure levels.

    + calcEDM                ::Bool              ... True if the EDM of selected levels is to be calculated.
    + calcEQM                ::Bool              ... True if the EQM of selected levels is to be calculated.
    + calcStarkshifts        ::Bool              ... True if the Stark energy shifts are to be computed.
    + gMultiplet              ::Multiplet         ... A local, approximate Green multiplet of intermediate
                                                       (opposite-parity) ASF that the underlying alpha_0/alpha_2
                                                       sum-over-states runs over -- generated INDEPENDENTLY,
                                                       ahead of time, by the caller (see the module docstring
                                                       and examples/example-Cl.jl); only used if calcStarkshifts.
    + printBefore             ::Bool              ... True if a list of selected levels is printed before the actual computations start.
    + EField                 ::Float64           ... Strength of the electric field in [V/cm]
    + levelSelection         ::LevelSelection    ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcEDM                  ::Bool
    calcEQM                  ::Bool
    calcStarkshifts          ::Bool
    gMultiplet                ::Multiplet
    printBefore               ::Bool
    EField                   ::Float64
    levelSelection           ::LevelSelection
end


"""
`StarkShift.Settings()`
    ... constructor for an `empty` instance of StarkShift.Settings.
"""
function Settings()
        Settings(false, false, false, Multiplet(), false, 0., LevelSelection() )
end


"""
`StarkShift.Settings(set::StarkShift.Settings;`

        calcEDM=.., calcEQM=.., calcStarkshifts=.., gMultiplet=.., printBefore=.., EField=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::StarkShift.Settings.
"""
function Settings(set::StarkShift.Settings;
        calcEDM::Union{Nothing,Bool}=nothing,                   calcEQM::Union{Nothing,Bool}=nothing,
        calcStarkshifts::Union{Nothing,Bool}=nothing,           gMultiplet::Union{Nothing,Multiplet}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,               EField::Union{Nothing,Float64}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcEDM)           calcEDMx          = set.calcEDM          else   calcEDMx          = calcEDM          end
    if  isnothing(calcEQM)           calcEQMx          = set.calcEQM          else   calcEQMx          = calcEQM          end
    if  isnothing(calcStarkshifts)   calcStarkshiftsx  = set.calcStarkshifts  else   calcStarkshiftsx  = calcStarkshifts  end
    if  isnothing(gMultiplet)        gMultipletx       = set.gMultiplet       else   gMultipletx       = gMultiplet       end
    if  isnothing(printBefore)       printBeforex      = set.printBefore      else   printBeforex      = printBefore      end
    if  isnothing(EField)            EFieldx           = set.EField           else   EFieldx           = EField           end
    if  isnothing(levelSelection)    levelSelectionx   = set.levelSelection   else   levelSelectionx   = levelSelection   end

    Settings( calcEDMx, calcEQMx, calcStarkshiftsx, gMultipletx, printBeforex, EFieldx, levelSelectionx )
end


# `Base.show(io::IO, settings::StarkShift.Settings)`  ... prepares a proper printout of the variable settings::StarkShift.Settings.
function Base.show(io::IO, settings::StarkShift.Settings)
    println(io, "calcEDM:               $(settings.calcEDM)  ")
    println(io, "calcEQM:               $(settings.calcEQM)  ")
    println(io, "calcStarkshifts:       $(settings.calcStarkshifts)  ")
    println(io, "gMultiplet:            $(settings.gMultiplet)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "EField:                $(settings.EField)  ")
    println(io, "levelSelection:        $(settings.levelSelection)  ")
end


"""
`StarkShift.computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)`
    ... to compute all amplitudes and properties of for a given level. An outcome::StarkShift.Outcome is returned for
        which the amplitudes and all requested properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::StarkShift.Outcome, grid::Radial.Grid, settings::StarkShift.Settings)
    edm = eqm = 0.
    if  settings.calcEDM    edm = StarkShift.computeTheta1(outcome.Jlevel, grid)    end
    if  settings.calcEQM    eqm = StarkShift.computeTheta2(outcome.Jlevel, grid)    end

    alpha0 = EmProperty(0.);   alpha2 = EmProperty(0.);   Jsublevels = outcome.Jsublevels
    if  settings.calcStarkshifts
        alpha0, alpha2, Jsublevels = StarkShift.computeStarkShifts(outcome.Jlevel, settings.gMultiplet, grid, settings)
    end

    newOutcome = StarkShift.Outcome( outcome.Jlevel, edm, eqm, alpha0, alpha2, Jsublevels)
    return( newOutcome )
end


"""
`StarkShift.AU_EFIELD_IN_VCM`
    ... the CODATA atomic unit of electric field, in V/cm (1 a.u. = 5.142206751e9 V/cm; verified
        against NIST/CODATA 31-Jul-2026, not from memory alone), used to convert
        StarkShift.Settings.EField (documented in V/cm) to atomic units.
"""
const AU_EFIELD_IN_VCM = 5.142206751e9


"""
`StarkShift.computeStarkShifts(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, settings::StarkShift.Settings)`
    ... computes the static Stark shift of every M-sublevel of `level` in a field of strength
        settings.EField [V/cm], from the level's own static E1 scalar/tensor polarizability
        (MultipolePolarizibility.computeScalarTensorPolarizability, reused directly rather than
        re-derived here -- see module docstring): DeltaE(J,M) = -(1/2) alpha_0(J) E^2 -
        (1/2) alpha_2(J) E^2 [3M^2-J(J+1)] / [J(2J-1)], undefined/trivially 0. for J<1 (matching
        MultipolePolarizibility's own alpha_2 guard). A tuple (alpha0::EmProperty, alpha2::EmProperty,
        Jsublevels::Array{StarkShift.SublevelJ,1}) is returned.
"""
function computeStarkShifts(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, settings::StarkShift.Settings)
    alpha0, alpha2, _ = MultipolePolarizibility.computeScalarTensorPolarizability(level, gMultiplet, grid)
    Efield  = settings.EField / StarkShift.AU_EFIELD_IN_VCM
    Efield2 = Efield^2
    Jd      = Basics.twice(level.J) / 2.0

    Jsublevels = StarkShift.SublevelJ[]
    for  M in AngularMomentum.m_values(level.J)
        Md = Basics.twice(M) / 2.0
        scalarShift = (-0.5) * alpha0 * Efield2
        if  Jd*(2*Jd-1) == 0.0
            tensorShift = EmProperty(0.)
        else
            tensorFactor = (3*Md^2 - Jd*(Jd+1)) / (Jd*(2*Jd-1))
            tensorShift  = (-0.5*tensorFactor) * alpha2 * Efield2
        end
        push!(Jsublevels, StarkShift.SublevelJ(M, scalarShift + tensorShift))
    end
    return  (alpha0, alpha2, Jsublevels)
end


"""
`StarkShift.computeTheta1(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-dipole-moment [Theta^(1)] for the given level; an edm::Float64 is returned.
"""
function  computeTheta1(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(1) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta1 = MultipoleMoment.emmStaticAmplitude(1, level, level, grid; display=false)
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(1), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta1 = Theta1 * (-1)^(2*J - 1 + 1) * w3j * sqrt(2*J+1)
         
    Theta1JG = MultipoleMoment.amplitude(1, level, grid, display=false)
    println("   Theta^(1) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta1)  $(Theta1JG)  e a_o")
    return( Theta1 )
end



"""
`StarkShift.computeTheta2(level::Level, grid::Radial.Grid)`  
    ... to compute the electric-quadrupole-moment [Theta^(2)] for the given level; an eqm::Float64 is returned.
        Until July 2024, it remains unclear why the Clebsch-Gordan coefficient does not give a proper result
"""
function  computeTheta2(level::Level, grid::Radial.Grid)
    print(">>> Calculate the Theta^(2) (shift) coefficient for level (index=$(level.index), J=$(level.J)) ...")
    M = AngularM64(level.J.num, level.J.den)
    Theta2 = MultipoleMoment.emmStaticAmplitude(2, level, level, grid; display=false) 
             ## AngularMomentum.ClebschGordan(level.J, M, AngularJ64(2), AngularM64(0), level.J, M) / sqrt(Basics.twice(level.J) + 1)
    # Use the conversion factor by Jan Gilles (2024)
    J      = level.J.num/level.J.den
    w3j    = AngularMomentum.Wigner_3j(J, 2, J, J, 0, -J)
    Theta2 = Theta2 * (-1)^(2*J - 2 + 1) * w3j * sqrt(2*J+1)

    Theta2JG = MultipoleMoment.amplitude(2, level, grid, display=false)
    println("   Theta^(2) (level=$(level.index) [J=$(level.J)$(string(level.parity))]) = $(Theta2)  $(Theta2JG)  e a_o^2")
    return( Theta2 )
end


"""
`StarkShift.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)`  
    ... to compute (as selected) the EDM and EQM factors for the levels of the given multiplet and as specified by the given 
        settings. The results are returned in table (if required) or nothing is returned otherwise.
        The nuclear model is not used at present.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkShift.Settings; output=true)
    println("")
    printstyled("StarkShift.computeOutcomes(): The computation of the electric multipole moments starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    outcomes = StarkShift.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    StarkShift.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    newOutcomes = StarkShift.Outcome[]
    for  outcome in outcomes
        newOutcome = StarkShift.computeAmplitudesProperties(outcome, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    StarkShift.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    StarkShift.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`StarkShift.determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings)`  
    ... to determine a list of Outcomes's for the computation of electric-multipole-moments and energy shifts for 
        levels from the given multiplet. It takes into account the particular selections and settings. 
        An Array{StarkShift.Outcome,1} is returned. Apart from the level specification, all physical properties are 
        still set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::StarkShift.Settings) 
    # Define values that depend on the requested computations
    outcomes   = StarkShift.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            Jsublevels = StarkShift.SublevelJ[];   Mvalues = AngularMomentum.m_values(level.J)
            for  M in Mvalues   push!(Jsublevels, StarkShift.SublevelJ(M, EmProperty(0.)) )    end
            push!( outcomes, StarkShift.Outcome(level, 0., 0., EmProperty(0.), EmProperty(0.), Jsublevels) )
        end
    end
    return( outcomes )
end


"""
`StarkShift.displayOutcomes(outcomes::Array{StarkShift.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
        levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{StarkShift.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected levels for EMM and Stark-shift computations:")
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
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.Jlevel.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`StarkShift.displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)`  
    ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{StarkShift.Outcome,1}, settings::StarkShift.Settings)
    #
    if  settings.calcEDM  ||  settings.calcEQM
        nx = 82
        println(stream, " ")
        println(stream, "  Atomic electric-multipole-moment (EMM) shift coefficients:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=5)              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(14, "Theta^(1)/EDM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o]"; na=5)        
        sa = sa * TableStrings.center(14, "Theta^(2)/EQM"; na=5)                          
        sb = sb * TableStrings.center(14, "[e a_o^2]"; na=5)        
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            energy = outcome.Jlevel.energy
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta1) )              * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.Theta2) )              * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    if  settings.calcStarkshifts
        nx1 = 106
        println(stream, " ")
        println(stream, "  Static scalar/tensor polarizabilities entering the Stark shift (Coulomb and Babushkin gauge):")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx1))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(15, "alpha_0 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_2 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        sa = sa * TableStrings.center(15, "alpha_2 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx1))
        #
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Coulomb) )      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin) )    * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Coulomb) )      * "    "
            sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Babushkin) )    * "    "
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx1))
        #
        nx2 = 92
        println(stream, " ")
        println(stream, "  Stark-shifts (energy-shifts) of the atomic sublevels at E = $(settings.EField) V/cm:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx2))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(8,  "M";     na=4);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(17, "Delta-E [Coul]"; na=3);                    sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        sa = sa * TableStrings.center(17, "Delta-E [Bab]" ; na=3);                    sb = sb * TableStrings.center(17, "[Hz]"; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx2))
        #
        for  outcome in outcomes
            sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            for  Jsublevel in outcome.Jsublevels
                sa  = "  "
                sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * TableStrings.center(8,  string(Jsublevel.M); na=4)
                sa = sa * TableStrings.flushright(17, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Hz", Jsublevel.energy.Coulomb)) )   * "    "
                sa = sa * TableStrings.flushright(17, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Hz", Jsublevel.energy.Babushkin)) ) * "    "
                println(stream, sa )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx2))
    end
    #
    return( nothing )
end

end # module
