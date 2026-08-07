
"""
`module  JAC.MultipolePolarizibility`
... a submodel of JAC that contains all methods for computing static multipole polarizibilities.

    Stage 1 (implemented here): the static electric-dipole (E1) scalar and tensor polarizabilities
    alpha_0(J), alpha_2(J), following the classic Angel & Sandars (1968) rank-0/rank-2 decomposition
    (formula and prefactors confirmed against van Leeuwen & Hogervorst, Z. Phys. A 316, 149 (1984),
    Eq. (5)):

        alpha_0(J) = -2/(3(2J+1)) * sum_J' |<J||P||J'>|^2 / (W_J - W_J')

        alpha_2(J) = -2*sqrt(10/3)*sqrt(J(2J-1)/((2J+3)(J+1)(2J+1))) * (-1)^(2J)
                     * sum_J' (-1)^(J-J') {J J' 1; 1 2 J} |<J||P||J'>|^2 / (W_J - W_J')

    (alpha_2 is defined to be exactly 0 for J<1, where the tensor polarizability does not exist).

    The sum over J' runs, in principle, over the complete (opposite-parity) many-electron spectrum,
    including the continuum -- literally enumerating this is out of reach. Following the same
    approach already used by module-MultiPhotonTransition.jl for its own (structurally identical)
    sum-over-virtual-intermediate-states problem, the sum is instead carried out over a finite,
    user-supplied set of intermediate ASF -- a *local, approximate Green multiplet*
    (MultipolePolarizibility.Settings.gMultiplet) that brings in the physically relevant perturber
    levels (2p_1/2, 2p_3/2, ... for an ns ground level, say). This multiplet is generated
    INDEPENDENTLY, ahead of time, by the caller -- simply as a plain Atomic.Computation of energies
    & ASF (no properties) for the relevant opposite-parity configurations, e.g. 2p, 3p, 4p, ... (see
    examples/example-Cj.jl) -- and is simply consumed here; this module performs no SCF/CI generation
    of its own (mirroring MultiPhotonTransition.Settings.gMultiplet exactly). A high-n tail of such
    perturber configurations (e.g. up to 40p, 60p, ...) may also be included deliberately: these
    orbitals are no longer physically accurate bound Rydberg states at the default grid/basis size
    (see project notes), but remain properly normalized, mutually orthogonal ASF and so provide a
    simple, explicit (if crude) discrete approximation ("pseudo-continuum") to the true continuum's
    contribution -- their individual contributions can be tabulated and inspected exactly like any
    other perturber, without any change to the formula or code below.

    Every contribution (and the total alpha_0/alpha_2) is computed and reported in BOTH the Babushkin
    (length; MultipoleMoment.emmStaticAmplitude, the standard static E1 multipole-moment operator) and
    Coulomb (velocity; InteractionStrength.MabEmissionJohnsony's Johnson-formalism Coulomb-gauge
    reduced matrix element, evaluated at omega=0 where it is finite and well-defined even though the
    Babushkin-gauge form of that SAME finite-frequency machinery vanishes identically at omega=0) forms
    -- as an EmProperty, following the Coulomb/Babushkin bookkeeping already used throughout JAC (e.g.
    Einstein.jl). The two are NOT calibrated onto a common scale against one another (deliberately, for
    now -- Johnson's Coulomb-gauge reduced matrix element carries its own internal normalization
    convention, not yet independently understood/verified against the plain static multipole-moment
    operator); they are reported side by side purely for direct, qualitative gauge-agreement
    monitoring, exactly the same diagnostic role gauge comparison already plays for ordinary E1
    transition rates elsewhere in JAC.

    Stage 2 (explicitly NOT implemented here, deferred): dynamic (nonzero omega) polarizabilities,
    general multipoles beyond E1 (needs the general rank-K decomposition, not just rank-0/rank-2),
    magnetic multipoles. Stark shifts and dispersion coefficients are separate, later consumers of
    this module and are also not implemented here.
"""
module MultipolePolarizibility


using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..ManyElectron, ..MultipoleMoment, ..Nuclear, ..PhotoEmission, ..Radial, ..TableStrings


"""
`struct  MultipolePolarizibility.Contribution`
    ... defines a type to keep the individual contribution of ONE intermediate (perturber) ASF from
        gMultiplet to the scalar/tensor static E1 polarizability of some other (fixed) level.

    + nuLevel        ::Level        ... The intermediate (perturber) level nuLevel = J' from gMultiplet.
    + deltaE         ::Float64      ... Energy denominator W_J - W_J' (atomic units) entering this term.
    + contribAlpha0  ::EmProperty   ... This perturber's own additive contribution to alpha_0(J), in a.u.,
                                        in both the Coulomb and Babushkin gauge (see module docstring).
    + contribAlpha2  ::EmProperty   ... This perturber's own additive contribution to alpha_2(J), in a.u.,
                                        in both gauges; identically 0. for J<1, where the tensor
                                        polarizability does not exist.
"""
struct Contribution
    nuLevel          ::Level
    deltaE           ::Float64
    contribAlpha0    ::EmProperty
    contribAlpha2    ::EmProperty
end


"""
`struct  MultipolePolarizibility.Outcome`  ... defines a type to keep the outcome of a static-polarizability computation.

    + level          ::Level                                 ... Atomic level to which the outcome refers to.
    + alpha0         ::EmProperty                             ... static scalar (rank 0) electric-dipole
                                                                    polarizability alpha_0(J), in atomic units,
                                                                    in both the Coulomb and Babushkin gauge.
    + alpha2         ::EmProperty                             ... static tensor (rank 2) electric-dipole
                                                                    polarizability alpha_2(J), in atomic units,
                                                                    in both gauges; identically 0. for J<1,
                                                                    where the tensor polarizability does not exist.
    + contributions  ::Array{MultipolePolarizibility.Contribution,1}
                                                                ... individual, additive contribution of every
                                                                    intermediate ASF in gMultiplet that couples
                                                                    to level via E1 (see module docstring).
"""
struct Outcome
    level          ::Level
    alpha0         ::EmProperty
    alpha2         ::EmProperty
    contributions  ::Array{MultipolePolarizibility.Contribution,1}
end


"""
`MultipolePolarizibility.Outcome()`
    ... constructor for an `empty` instance of MultipolePolarizibility.Outcome for the computation of static
        polarizibilities.
"""
function Outcome()
    Outcome(Level(), EmProperty(0.), EmProperty(0.), MultipolePolarizibility.Contribution[])
end


# `Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)`
#		 ... prepares a proper printout of the variable outcome::MultipolePolarizibility.Outcome.
function Base.show(io::IO, outcome::MultipolePolarizibility.Outcome)
    println(io, "level:                     $(outcome.level)  ")
    println(io, "alpha0:                    $(outcome.alpha0)  ")
    println(io, "alpha2:                    $(outcome.alpha2)  ")
    println(io, "contributions:             $(length(outcome.contributions)) perturber level(s)  ")
end


"""
`struct  MultipolePolarizibility.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of computing multipolar polarizibilities.

    + multipoles                ::Array{EmMultipole,1}  ... Multipoles to be considered for the polarity. Stage 1
                                                              only supports E1 (the classic scalar/tensor
                                                              polarizability); any other multipole present in the
                                                              list is skipped with a warning.
    + gMultiplet                ::Multiplet             ... A local, approximate Green multiplet of intermediate
                                                              (opposite-parity) ASF that the sum-over-states runs
                                                              over -- generated INDEPENDENTLY, ahead of time, by the
                                                              caller (e.g. via AtomicState.GreenExpansion/generate(),
                                                              see examples/example-Aj.jl) and simply consumed here;
                                                              mirrors MultiPhotonTransition.Settings.gMultiplet.
    + omegas                    ::Array{Float64,1}      ... List of omegas (energies) of the dynamic polarizibility.
                                                              Stage 1 only supports the static case (omega=0);
                                                              this list is currently unused.
    + printBefore               ::Bool                  ... True if a list of selected levels is printed before the
                                                            actual computations start.
    + levelSelection            ::LevelSelection        ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    multipoles                  ::Array{EmMultipole,1}
    gMultiplet                  ::Multiplet
    omegas                      ::Array{Float64,1}
    printBefore                 ::Bool
    levelSelection              ::LevelSelection
end


"""
`MultipolePolarizibility.Settings()`
    ... constructor for an `empty` instance of MultipolePolarizibility.Settings for the computation dynamic polarizibilities.
"""
function Settings()
    Settings(EmMultipole[], Multiplet(), Float64[], false, LevelSelection() )
end


"""
`MultipolePolarizibility.Settings(set::MultipolePolarizibility.Settings;`

        multipoles=.., gMultiplet=.., omegas=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::MultipolePolarizibility.Settings.
"""
function Settings(set::MultipolePolarizibility.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
        gMultiplet::Union{Nothing,Multiplet}=nothing,
        omegas::Union{Nothing,Array{Float64,1}}=nothing,   printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(multipoles)       multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(gMultiplet)       gMultipletx     = set.gMultiplet     else   gMultipletx     = gMultiplet     end
    if  isnothing(omegas)           omegasx         = set.omegas         else   omegasx         = omegas         end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( multipolesx, gMultipletx, omegasx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings:MultipolePolarizibility.Settings)`
#		... prepares a proper printout of the variable settings::MultipolePolarizibility.Settings.
function Base.show(io::IO, settings::MultipolePolarizibility.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
    println(io, "omegas:                   $(settings.omegas)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`MultipolePolarizibility.computeScalarTensorPolarizability(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)`
    ... computes the static E1 scalar (alpha_0) and tensor (alpha_2) polarizabilities of `level` (Angel & Sandars 1968;
        prefactors as given in the module docstring) by summing the reduced E1 matrix elements over every level
        `nuLevel` in `gMultiplet`, in BOTH the Babushkin gauge (`MultipoleMoment.emmStaticAmplitude(1, nuLevel,
        level, grid)`, the standard static multipole-moment operator) and the Coulomb gauge
        (`PhotoEmission.amplitude` with `Basics.Coulomb` at omega=0., i.e. the zero-frequency limit of Johnson's
        finite-frequency radiative reduced matrix element, which -- unlike the Babushkin-gauge form of that SAME
        machinery -- remains finite and well-defined at omega=0; see module docstring). A tuple
        (alpha0::EmProperty, alpha2::EmProperty, contributions::Array{MultipolePolarizibility.Contribution,1})
        is returned, the latter giving each perturber level's own additive share of alpha0/alpha2 in both gauges.
"""
function computeScalarTensorPolarizability(level::Level, gMultiplet::Multiplet, grid::Radial.Grid)
    Jd = Basics.twice(level.J) / 2.0
    rawScalar = EmProperty(0.);   rawTensor = EmProperty(0.)
    prefactor0 = -2.0 / (3.0*(2*Jd+1))
    if  Jd < 1.0    prefactor2 = 0.0
    else            prefactor2 = -2.0*sqrt(10.0/3.0) * sqrt( Jd*(2*Jd-1) / ((2*Jd+3)*(Jd+1)*(2*Jd+1)) ) * (-1.0)^round(Int64, 2*Jd)
    end
    contributions = MultipolePolarizibility.Contribution[]

    for  nuLevel in gMultiplet.levels
        deltaE = level.energy - nuLevel.energy
        if  deltaE == 0.0   continue    end

        # Both emmStaticAmplitude's and PhotoEmission.amplitude's SpinAngular/CSF machinery need
        # BOTH levels' CSFs indexed against the SAME subshell list; nuLevel (from gMultiplet's
        # typically much larger pseudo-basis) and level (the level of interest's own, typically
        # much smaller basis) do not share one -- reduce both onto their union first (same idiom
        # as AutoIonization/InternalConversion's amplitude construction).
        subshellList = Basics.generate(OrderedSubshellList(), nuLevel.basis, level.basis)
        newNuLevel   = Basics.generateLevelWithSymmetryReducedBasis(nuLevel, subshellList)
        newLevel     = Basics.generateLevelWithSymmetryReducedBasis(level,   subshellList)

        dBabushkin = MultipoleMoment.emmStaticAmplitude(1, newNuLevel, newLevel, grid)
        dCoulomb   = PhotoEmission.amplitude(Basics.Emission(), E1, Basics.Coulomb, 0., newNuLevel, newLevel, grid)
        S = EmProperty( abs2(dCoulomb), abs2(dBabushkin) )
        if  S.Coulomb == 0.0  &&  S.Babushkin == 0.0   continue    end

        Sd = S * (1.0/deltaE)
        contribScalar = prefactor0 * Sd
        rawScalar     = rawScalar + Sd

        Jpd    = Basics.twice(nuLevel.J) / 2.0
        phase  = (-1.0)^round(Int64, Jd - Jpd)
        w6j    = AngularMomentum.Wigner_6j(level.J, nuLevel.J, AngularJ64(1), AngularJ64(1), AngularJ64(2), level.J)
        contribTensor = (prefactor2 * phase * w6j) * Sd
        rawTensor     = rawTensor + (phase * w6j) * Sd

        push!(contributions, MultipolePolarizibility.Contribution(nuLevel, deltaE, contribScalar, contribTensor))
    end

    alpha0 = prefactor0 * rawScalar
    alpha2 = prefactor2 * rawTensor

    return  (alpha0, alpha2, contributions)
end


"""
`MultipolePolarizibility.computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, gMultiplet::Multiplet,
                                                        grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
    ... to compute the static scalar/tensor E1 polarizabilities for a given level; an outcome::MultipolePolarizibility.Outcome is
        returned for which alpha0/alpha2/contributions are now evaluated explicitly. If E1 is not among
        settings.multipoles (or gMultiplet is empty), alpha0=alpha2=0. is returned without further computation.
"""
function  computeAmplitudesProperties(outcome::MultipolePolarizibility.Outcome, gMultiplet::Multiplet,
                                       grid::Radial.Grid, settings::MultipolePolarizibility.Settings)
    if  !(E1 in settings.multipoles)  ||  isempty(gMultiplet.levels)   return  outcome   end

    alpha0, alpha2, contributions = MultipolePolarizibility.computeScalarTensorPolarizability(outcome.level, gMultiplet, grid)
    return  MultipolePolarizibility.Outcome(outcome.level, alpha0, alpha2, contributions)
end


"""
`MultipolePolarizibility.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                                settings::MultipolePolarizibility.Settings; output=true)`
    ... to compute (as selected) the static E1 scalar/tensor polarizabilities for the levels of the given multiplet and
        as specified by the given settings. The intermediate (perturber) ASF are taken directly from
        settings.gMultiplet, which the caller must generate independently, ahead of time (see the module docstring
        and examples/example-Aj.jl/example-Cj.jl). The results are printed in neat tables to screen but nothing is
        returned otherwise.
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
    #
    # Calculate all amplitudes and requested properties
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    newOutcomes = MultipolePolarizibility.Outcome[]
    for  outcome in outcomes
        newOutcome = MultipolePolarizibility.computeAmplitudesProperties(outcome, settings.gMultiplet, grid, settings)
        push!( newOutcomes, newOutcome)
        MultipolePolarizibility.displayContributions(stdout, newOutcome)
        if  printSummary    MultipolePolarizibility.displayContributions(iostream, newOutcome)   end
    end
    # Print all results to screen
    MultipolePolarizibility.displayResults(stdout, newOutcomes)
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
            push!( outcomes, MultipolePolarizibility.Outcome(level, EmProperty(0.), EmProperty(0.), MultipolePolarizibility.Contribution[]) )
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
    ... to display the energies and static scalar/tensor E1 polarizabilities for the selected levels. A neat table is
        printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{MultipolePolarizibility.Outcome,1})
    nx = 122
    println(stream, " ")
    println(stream, "  Static electric-dipole scalar/tensor polarizibilities (Coulomb and Babushkin gauge, not calibrated " *
                     "onto a common scale -- see module docstring):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(15, "alpha_0 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_2 [Coul]"; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "alpha_2 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy))    * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Coulomb) )                  * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin) )                * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Coulomb) )                  * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha2.Babushkin) )                * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end


"""
`MultipolePolarizibility.displayContributions(stream::IO, outcome::MultipolePolarizibility.Outcome)`
    ... to display, for the single level outcome.level, the individual (additive) contribution of every
        intermediate ASF in gMultiplet that couples to it via E1, sorted by descending |contribAlpha0| so that
        the dominant perturber levels (e.g. 2p_1/2, 2p_3/2 for an ns ground level) appear first. A neat table is
        printed but nothing is returned otherwise. Nothing is printed if outcome.contributions is empty.
"""
function  displayContributions(stream::IO, outcome::MultipolePolarizibility.Outcome)
    if  isempty(outcome.contributions)   return( nothing )   end
    nx = 138
    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
    println(stream, " ")
    println(stream, "  Individual perturber contributions to alpha_0/alpha_2 of level $(TableStrings.level(outcome.level.index)) " *
                     "($(string(sym))), Coulomb and Babushkin gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Perturber"; na=2);                          sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J'^P'";     na=4);                          sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy";    na=4);
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Delta E";   na=4)
    sb = sb * TableStrings.center(14, "[a.u.]"; na=4)
    sa = sa * TableStrings.center(15, "d-alpha_0 [Coul]"; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_0 [Bab]" ; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_2 [Coul]"; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    sa = sa * TableStrings.center(15, "d-alpha_2 [Bab]" ; na=3);                   sb = sb * TableStrings.center(15, "[a.u.]" ; na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    ordered = sort(outcome.contributions, by = c -> -abs(c.contribAlpha0.Babushkin))
    for  c in ordered
        sa  = "  ";    nuSym = LevelSymmetry( c.nuLevel.J, c.nuLevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(c.nuLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(nuSym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", c.nuLevel.energy)) * "    "
        sa = sa * TableStrings.flushright(14, @sprintf("%.6e", c.deltaE) )                          * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha0.Coulomb) )           * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha0.Babushkin) )         * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha2.Coulomb) )           * "    "
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", c.contribAlpha2.Babushkin) )         * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end

end # module
