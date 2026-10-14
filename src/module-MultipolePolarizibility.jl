
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

    Every contribution (and the total alpha_0/alpha_2) is computed from the frequency-INDEPENDENT static
    multipole-moment operator (MultipoleMoment.emmStaticAmplitude, i.e. Johnson's r^k C_k), which is the
    operator a static polarizability -- and a C_6/C_8/C_10 dispersion coefficient -- actually requires. There
    is no gauge freedom in a static quantity, so the two components of the EmProperty are identical; they are
    retained only because Outcome carries EmProperty fields throughout JAC.

    Until 10-Aug-2026 the Coulomb component was instead taken from PhotoEmission.amplitude at omega = 0,
    i.e. from InteractionStrength.MabEmission -- the frequency-DEPENDENT transition operator -- and the two
    columns were presented side by side as a gauge-agreement diagnostic. That was a category error rather
    than a calibration one. The velocity and length forms are connected by <f|p|i> = i (E_f - E_i) <f|r|i>,
    so evaluating the velocity form at omega = 0 drops exactly the 1/omega that relates them; the result was
    alpha_0(Coulomb) ~ 1.5e-6 against 1.164 in Babushkin, six orders apart and completely insensitive to the
    radial box -- which is what identified it as the wrong operator rather than a numerical problem. There is
    no gauge diagnostic to be had from a static quantity, and none is offered here.

    Stage 2 (explicitly NOT implemented here, deferred): dynamic (nonzero omega) polarizabilities,
    general multipoles beyond E1 (needs the general rank-K decomposition, not just rank-0/rank-2),
    magnetic multipoles. Stark shifts and dispersion coefficients are separate, later consumers of
    this module and are also not implemented here.

    UNPARKED 13-Oct-2026 ON THE MAINTAINER'S INSTRUCTION, AND THE MODULE REPRODUCES EXACT HYDROGEN.
    It had been parked since 10-Aug-2026, and the whole of the reason was a stale number. The story is kept
    here because it cost two months and the mechanism is general.

      * **alpha_0(H 1s) = 4.499752 a.u. against the exact 4.500000** -- a relative error of 6e-05 -- from
        np perturbers n = 2...16 with the grid guard switched off (see below). The approach to it is what
        shows the pseudo-continuum is real rather than fortuitous: 3.4918 at 2p..4p, 3.6258 at 2p..8p --
        saturating towards the BOUND-series limit of about 3.663 -- then 4.2486 at 2p..12p as the high-n
        states cross threshold and begin to act as a continuum, 4.4997 at 2p..16p, and 4.4998 at 2p..20p,
        i.e. converged. **So the question the old block left open -- whether such a tail converges to the
        right number -- is answered, and it does.**

      * THE BOUND PART IS EXACT TERM BY TERM. A 2p+3p+4p sum gives 3.491829 against the closed-form 3.4921,
        and perturber by perturber 2.9595 / 0.4004 / 0.1320 against hydrogen's 2.9596 / 0.4005 / 0.1320, with
        each spin-orbit pair splitting 1 : 2 by its statistical weight. Both gauge columns are identical by
        construction, as a static quantity requires.

      * **THE 1.164 THAT GOT IT PARKED WAS THE RIGHT ANSWER DIVIDED BY THREE, AND THE DIVISION HAD ALREADY
        BEEN REPAIRED WHEN THE NOTE WAS WRITTEN.** Dividing today's alpha_0 by 3 reproduces 1.1639429, i.e.
        the recorded value to seven figures. The repair was `0bdff9f6` (the Mab/Mba consolidation) of
        09-Aug-2026 -- the day BEFORE the parking note -- and nobody re-ran the module afterwards. A working
        module therefore sat parked for two months behind a number that was obsolete the day it was recorded.
        `testModule_MultipolePolarizibility` now asserts the closed-form hydrogen values and has been shown to
        fail on exactly this error; the old assertion, 0 < alpha_0 < 4.5, passed 1.164 as happily as 3.4918.

    WHAT IT IS GOOD FOR NOW: static scalar and tensor polarizabilities as absolute values, not merely as lower
    bounds, PROVIDED the perturber set includes the high-n tail. What it still does not do is choose that tail
    for the caller -- see the grid note below.

    THE GRID, AND WHY THE GUARD MUST BE SWITCHED OFF FOR THE TAIL. Beyond about 5p on a 129 a.u. box,
    `Bsplines.checkGridRepresentation` refuses the perturbers, and it is RIGHT to: they are not accurate
    Rydberg states there. But that is exactly what a pseudo-continuum is -- box-confined, properly normalised,
    mutually orthogonal ASF standing in for the true continuum. So a caller who wants the converged value
    passes `AsfSettings(...; gridStopper = false)` when GENERATING THE PERTURBER MULTIPLET, deliberately and
    with this paragraph as the reason. The test of whether the basis is adequate is then not the guard but the
    CONVERGENCE OF THE SUM, which for hydrogen is checked against an exact answer above.

    THE UNITS, because comparing against a measured Stark shift without them has already produced one false alarm.
    `alpha_0` is returned and printed in ATOMIC UNITS. The Stark-shift literature quotes MHz/(kV/cm)^2, and the two
    differ by 1 a.u. = 1.244160e-04 MHz/(kV/cm)^2. Dividing the printed a.u. value by a published MHz/(kV/cm)^2 one
    yields a meaningless ratio near 2000-3000, which was recorded for eight days as "the module is too large by a
    factor of ~2600". CONVERTED, the same two Ba 6s8s cases give 1.63 against 5.68 and 1.88 against 4.68 -- a factor
    ~3 LOW, which is the lower bound this block describes, and consistent with H 1s at 0.46 of the exact value.
    So: convert before forming any ratio.

    **AND THE "NEVER ABOVE AN EXPERIMENT" THAT USED TO END THIS PARAGRAPH IS WITHDRAWN -- 14-Oct-2026, with a
    measured counterexample.** The lower-bound property is real but much narrower than that sentence claimed,
    and the two limits matter in practice:

      * it holds only for the LOWEST level of its parity, where every energy denominator carries the same sign
        so that each omitted perturber can only ADD. For an excited level such as Ba 6s8s the perturbers lying
        BELOW it contribute with the opposite sign, and the sum is then bounded in neither direction;
      * and even for a ground state it bounds only what THE GIVEN MEAN FIELD would yield with a complete
        perturber set -- it says nothing whatever about agreement with EXPERIMENT, because the field itself can
        overshoot. Measured on Li 2s against Puchalski et al.'s 164.0740 a.u., same perturber set throughout:
        DFS gives 183.7863 (12 % ABOVE), HS gives 227.9872 (39 % ABOVE), AL gives 167.6745 (2.2 % above).
        DFS underbinds the valence 2s to -0.079 Ha against the true -0.198, shrinking the 2s-2p denominator.

    **THAT IS ALSO THE USEFUL DIAGNOSTIC, so it is worth stating the other way round:** a ground-state value
    coming out ABOVE a trusted reference CANNOT be blamed on an incomplete perturber set, and points at the
    mean field instead. That is how the Li case was solved -- see `examples/example-Cj.jl` branch (C), which
    now runs on `ALField` for exactly this reason.
"""
module MultipolePolarizibility


using Printf, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..ManyElectron, ..MultipoleMoment, ..Nuclear, ..PhotoEmission, ..Radial, ..TableStrings


"""
`MultipolePolarizibility.AU_ALPHA_IN_MHZ_PER_KVCM2`
    ... the factor that converts a static E1 polarizability from atomic units into MHz/(kV/cm)^2, the unit in which
        the Stark-shift literature quotes alpha_0 and alpha_2; a `factor::Float64` is returned. It is built as
        (1/2) * (E_h/h in Hz) * (1 kV/cm / E_au)^2 * 1e-6, with E_au = 5.142206751e9 V/cm the atomic unit of electric
        field strength, and the leading 1/2 is the customary half of DeltaE = -(1/2) alpha E^2 that experimenters
        absorb into the quoted number. Its value is 1.244160e-04, so a polarizability of 13098 a.u. is 1.63 in the
        literature's unit. Dividing a printed a.u. value by a published MHz/(kV/cm)^2 one -- rather than converting
        first -- yields a spurious ratio of some thousands, which is exactly how this module was once believed to be
        too large by a factor of 2600 when it is in fact a factor ~3 LOW; see the module docstring.
"""
const AU_ALPHA_IN_MHZ_PER_KVCM2 = 1.0e-6 * Defaults.convertUnits("energy: from atomic to Hz", 1.0) *
                                  0.5 * (1000.0 / 5.142206751e9)^2


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
        `nuLevel` in `gMultiplet`, using `MultipoleMoment.emmStaticAmplitude(1, nuLevel, level, grid)` -- the
        frequency-INDEPENDENT static multipole-moment operator. A static polarizability carries no gauge freedom,
        so both components of the returned EmProperty hold the same number; see the note in the body. A tuple
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

        ## A STATIC POLARIZABILITY HAS NO GAUGE FREEDOM (corrected 10-Aug-2026), so both components of the
        ## EmProperty are filled from the same, frequency-INDEPENDENT static multipole moment. Until now the
        ## Coulomb component came from PhotoEmission.amplitude(.., Basics.Coulomb, 0., ..), i.e. MabEmission --
        ## the frequency-DEPENDENT transition operator -- evaluated at omega = 0. That is a category error, not
        ## a gauge: the velocity form is related to the length form by <f|p|i> = i (E_f - E_i) <f|r|i>, so
        ## taking it at omega = 0 drops exactly the 1/omega that connects the two. It produced alpha_0(Coulomb)
        ## ~ 1.5e-6 against 1.164 in Babushkin -- six orders, and completely insensitive to the radial box,
        ## which is what finally identified it as a wrong operator rather than a numerical problem.
        ## The operator used here, InteractionStrength.eMultipole via MultipoleMoment.emmStaticAmplitude, is
        ## Johnson's static r^k C_k: no omega, no gauge. It is also the operator that belongs in the C_6/C_8/C_10
        ## dispersion coefficients. The two columns are therefore deliberately identical; they are kept only
        ## because Outcome carries EmProperty fields throughout JAC.
        dStatic = MultipoleMoment.emmStaticAmplitude(1, newNuLevel, newLevel, grid)
        S = EmProperty( abs2(dStatic), abs2(dStatic) )
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
    nx = 141
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
    sa = sa * TableStrings.center(15, "alpha_0 [Bab]" ; na=3);                    sb = sb * TableStrings.center(15, "[MHz/(kV/cm)^2]"; na=3)
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
        sa = sa * TableStrings.flushright(15, @sprintf("%.6e", outcome.alpha0.Babushkin * MultipolePolarizibility.AU_ALPHA_IN_MHZ_PER_KVCM2) ) * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    # THE DOMINANT ENERGY DENOMINATOR, PRINTED BECAUSE IT IS THE ONE NUMBER THAT DIAGNOSES THE WHOLE FAMILY OF
    # FAILURES HERE (item 70, 14-Oct-2026).  alpha_0 ~ |<a||P||b>|^2 / (W_b - W_a), so it is INVERSELY
    # proportional to this gap, and a mean field that misplaces the lowest perturber misplaces alpha_0 by the
    # same factor -- while every internal check in this module still passes, because the formula is right and
    # only its input is wrong.  A reader who can compare ONE number against a measured level energy can catch
    # that in seconds;  a reader given only alpha_0 cannot catch it at all.
    #
    # MEASURED, Li [He]2s against Puchalski et al.'s 164.0740 a.u. with ONE perturber set throughout:
    #     DFSField  183.7863 (+12 %, gap 0.056168 Ha)    HSField  227.9872 (+39 %)    ALField  167.6745 (+2.2 %,
    #     gap 0.067688 Ha against the experimental 0.067700).  DFS puts the Li 2s at -0.079 Ha where the true
    #     value is -0.198, i.e. underbound by a factor 2.5, and the gap and alpha_0 move together exactly as
    #     the formula says they must.  `AsfSettings()` DEFAULTS to DFSField, so this is the common case, not
    #     an exotic one.
    for  outcome in outcomes
        isempty(outcome.contributions)   &&   continue
        dominant = outcome.contributions[1]
        for  c in outcome.contributions
            abs(c.contribAlpha0.Babushkin) > abs(dominant.contribAlpha0.Babushkin)   &&   (dominant = c)
        end
        share = abs(outcome.alpha0.Babushkin) > 0.0 ?
                abs(dominant.contribAlpha0.Babushkin / outcome.alpha0.Babushkin) : 0.0
        println(stream, "  Level $(outcome.level.index):  the DOMINANT perturber is " *
                        "$(LevelSymmetry(dominant.nuLevel.J, dominant.nuLevel.parity)) at " *
                        @sprintf("Delta-E = %.6f", abs(dominant.deltaE)) * " a.u. " *
                        "(" * @sprintf("%.6f", Defaults.convertUnits("energy: from atomic",
                                                                      abs(dominant.deltaE))) * " " *
                        TableStrings.inUnits("energy") * "), carrying " *
                        @sprintf("%.1f", 100*share) * " % of alpha_0.")
    end
    println(stream, "  ")
    println(stream, "  >>> alpha_0 is INVERSELY proportional to that gap, so CHECK IT against a measured level")
    println(stream, "  >>> energy before trusting alpha_0: a mean field that underbinds the valence electron")
    println(stream, "  >>> shrinks the gap and inflates alpha_0 by the same factor, with every internal check")
    println(stream, "  >>> in this module still passing.  Measured on Li [He]2s against a literature 164.0740:")
    println(stream, "  >>>     DFSField  +12 %      HSField  +39 %      ALField  +2.2 %")
    println(stream, "  >>> and `AsfSettings()` DEFAULTS to DFSField.  Prefer ALField or EOLField for a")
    println(stream, "  >>> polarizability;  see examples/example-Cj.jl branch (C).")
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
