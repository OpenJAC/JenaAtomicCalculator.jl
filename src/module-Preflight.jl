
"""
`module  Preflight`
    ... the checks that can be made on a computation BEFORE it is started, and the table of what each observable
        requires of a wave function.

        THE CHARTER, and it is the whole reason this module has a boundary rather than a topic. A check belongs here
        if it is DECIDABLE FROM THE INPUT ALONE AND COSTS NO SCF. Everything that needs a run to decide belongs in a
        verdict instead -- `SelfConsistent.ScfVerdict`, `AtomicState.RasLayerVerdict`, the CF and gauge-ratio quality
        flags -- and all of those already exist. Without that line a module of this kind receives every piece of
        advice anybody ever wants to give and becomes unsearchable.

        WHAT IT MAY SAY, AND WHAT IT MAY NOT. It reports WHAT IS MISSING FROM A SPACE and never HOW GOOD AN ANSWER
        WILL BE. The first is decidable; the second is not, and the distinction is not pedantry -- a module that
        ranked mean fields by their agreement with experiment would have recommended DFS over EOL for the Ca+ contact
        density, where DFS's closeness is a cancellation rather than accuracy. There is no INTERNAL quantity that
        ranks two mean fields against each other; only an external reference does that. For the same reason a layer
        increment is not an error bar and is never reported as one.

        WHY THE PHYSICS IS ONE TABLE. `OBSERVABLE_REQUIREMENTS` maps an observable to the radial region it weights,
        whether it needs the small component, and which classes of CSF must be present for it to mean anything. Every
        check in this module reads that table; nothing else carries physics. The table is also the reason the module
        is worth having at all: a variational calculation optimizes ONE functional, the energy, and the expectation
        value of any operator that does not commute with the Hamiltonian is wrong to FIRST order in the orbital
        displacement where the energy is wrong only to second. So a well-converged energy is no evidence whatever
        about a near-nucleus observable, and which observables a given space does constrain is knowable in advance.

        IT IS ADVISORY AND NEVER RAISES. `analyze` returns a `Preflight.Report`; a caller who wants a hard stop writes
        their own test on the findings. Refusing to start a computation is the user's decision and not this module's.
"""
module Preflight

using  Printf, ..Basics, ..ManyElectron, ..Nuclear, ..Radial, ..AtomicState, ..Atomic


"""
`struct  Preflight.ObservableRequirement`
    ... what one observable requires of a wave function; the physics content of the whole module.

    + observable       ::String             ... the observable in words, as it would be named in a paper.
    + radialWeight     ::String             ... the radial power the operator integrates, as "r^-2" or "r^+3".
    + needsSmallComponent ::Bool            ... true where the operator couples P to Q and therefore vanishes in the
                                                non-relativistic limit, as the magnetic hyperfine constants do.
    + essentialClasses ::Array{Symbol,1}    ... classes without which the observable is not constrained AT ALL.
    + refiningClasses  ::Array{Symbol,1}    ... classes that improve it; their absence degrades, it does not void.
    + reason           ::String             ... one sentence saying why, short enough to quote in a paper.

    THE TWO CLASS LISTS ARE SEPARATE BECAUSE A FLAT LIST OVER-REPORTS, and it did so the first time this module was
    run: a Lande g factor was declared UNCONSTRAINED for want of valence correlation, when a single configuration
    gives a perfectly meaningful Lande value that correlation merely refines. A severity that does not distinguish
    "meaningless without" from "better with" is not a severity.
"""
struct  ObservableRequirement
    observable            ::String
    radialWeight          ::String
    needsSmallComponent   ::Bool
    essentialClasses      ::Array{Symbol,1}
    refiningClasses       ::Array{Symbol,1}
    reason                ::String
end


"""
`Preflight.OBSERVABLE_REQUIREMENTS`
    ... the table itself, keyed by a short symbol. The radial powers are taken from the integrands in
        `InteractionStrength` rather than from the literature, so that the table cannot drift away from the code:
        `hfs_tM1` integrates r^-2 off-diagonally, `hfs_tE2` r^-3 diagonally and `hfs_tM3` r^-4 off-diagonally.
        An off-diagonal integrand is PQ+QP and so needs the small component; a diagonal one is PP+QQ and does not.
"""
const OBSERVABLE_REQUIREMENTS = Dict{Symbol,ObservableRequirement}(
    :hyperfineA     => ObservableRequirement("A, the magnetic-dipole hyperfine constant", "r^-2", true,
                            [:coreValenceS], [:valenceCorrelation],
                            "it weights the region closest to the nucleus"),
    :hyperfineB     => ObservableRequirement("B, the electric-quadrupole hyperfine constant", "r^-3", false,
                            Symbol[], [:valenceCorrelation],
                            "it weights the open shell itself rather than the nuclear region"),
    :hyperfineC     => ObservableRequirement("C, the magnetic-octupole hyperfine constant", "r^-4", true,
                            [:coreValenceS], [:valenceCorrelation],
                            "it has the steepest radial weight of any observable here and is the most sensitive " *
                            "to the inner tail of the open-shell orbital"),
    :contactDensity => ObservableRequirement("the electron density at the nucleus", "the r -> 0 limit", false,
                            [:coreValenceS], Symbol[],
                            "an energy objective does not constrain the density at the origin at all, however " *
                            "well converged it is"),
    :landeG         => ObservableRequirement("the Lande g factor", "r^0", false,
                            Symbol[], [:valenceCorrelation],
                            "it is unweighted in r and is governed by the angular composition of the level"),
    :multipoleRate  => ObservableRequirement("a bound-bound multipole transition rate", "r^+L", false,
                            [:bothStatesOneOrbitalSet], [:valenceCorrelation],
                            "a transition moment between separately optimized states is not a matrix element " *
                            "between orthonormal sets until the two bases are transformed"),
    :isotopeShift   => ObservableRequirement("a field or mass isotope shift", "the r -> 0 region", false,
                            [:coreValenceS], [:valenceCorrelation],
                            "the field shift is a contact observable and inherits the same blindness as the " *
                            "contact density") )


"""
`struct  Preflight.Finding`
    ... one thing worth saying about a computation before it is run.

    + severity    ::Symbol    ... :unconstrained, :degraded or :note; see below.
    + subject     ::Symbol    ... :grid, :csfSpace, :route, :freezing, :cost or :threading.
    + observable  ::String    ... the observable at issue, or "" where the finding is general.
    + message     ::String    ... what the input provides and what the observable needs.
    + remedy      ::String    ... the concrete change that would answer it, or "" where none is known.

    THE SEVERITIES ARE NAMED AFTER WHAT THEY SAY ABOUT THE OBSERVABLE, NOT AFTER HOW ALARMED THE MODULE IS, and that
    is deliberate. `:unconstrained` means a required class of CSF is ABSENT, so the input cannot constrain this
    observable at all -- decidable, and the strongest statement this module is entitled to make. `:degraded` means
    the observable is constrained but a known class is missing. `:note` is informational. A scale reading
    "ok / poor / bad" would be an accuracy claim in disguise, which the module docstring forbids.
"""
struct  Finding
    severity      ::Symbol
    subject       ::Symbol
    observable    ::String
    message       ::String
    remedy        ::String

    function Finding(severity::Symbol, subject::Symbol, observable::String, message::String, remedy::String)
        if  !(severity in (:unconstrained, :degraded, :note))
            error("Preflight.Finding: severity = :$severity is not recognized.  Use :unconstrained (a required " *
                  "class of CSF is absent), :degraded (constrained, but a known class is missing) or :note.")
        end
        if  !(subject in (:grid, :csfSpace, :route, :freezing, :cost, :threading))
            error("Preflight.Finding: subject = :$subject is not recognized.  Use :grid, :csfSpace, :route, " *
                  ":freezing, :cost or :threading.")
        end
        new(severity, subject, observable, message, remedy)
    end
end


"""
`struct  Preflight.Report`
    ... everything `analyze` has to say about one computation.

    + label      ::String                   ... the computation, named as Rule 21 asks: the path and the date.
    + findings   ::Array{Finding,1}         ... in the order the checks were made, not sorted by severity.
    + csfCounts  ::Dict{String,Int64}       ... per symmetry, where a count is available without building the space.
    + costNote   ::String                   ... what is known about the cost, never a promise.
"""
struct  Report
    label         ::String
    findings      ::Array{Finding,1}
    csfCounts     ::Dict{String,Int64}
    costNote      ::String
end


"""
`Base.show(io::IO, finding::Preflight.Finding)`
    ... prints one finding on two or three lines, with the severity first so that a long report can be skimmed.
"""
function Base.show(io::IO, finding::Preflight.Finding)
    tag = finding.severity == :unconstrained ? "UNCONSTRAINED" :
          finding.severity == :degraded      ? "degraded     " : "note         "
    obs = isempty(finding.observable) ? "" : "  [" * finding.observable * "]"
    println(io, "  " * tag * "  (" * string(finding.subject) * ")" * obs)
    println(io, "      " * finding.message)
    if  !isempty(finding.remedy)   println(io, "      remedy: " * finding.remedy)   end
end


"""
`Base.show(io::IO, report::Preflight.Report)`
    ... prints the whole report, with a one-line census first and the findings beneath it.
"""
function Base.show(io::IO, report::Preflight.Report)
    nUnc = count(f -> f.severity == :unconstrained, report.findings)
    nDeg = count(f -> f.severity == :degraded,      report.findings)
    nNot = count(f -> f.severity == :note,          report.findings)
    println(io, "\n  Preflight report for: " * report.label)
    println(io, "  " * "-"^104)
    @printf(io, "  %d unconstrained, %d degraded, %d notes.  This report is ADVISORY: nothing here stops a run.\n",
            nUnc, nDeg, nNot)
    if  !isempty(report.csfCounts)
        println(io, "  CSF counts: " * join([k * " = " * string(v) for (k,v) in report.csfCounts], ",  "))
    end
    if  !isempty(report.costNote)   println(io, "  cost: " * report.costNote)   end
    println(io, "  " * "-"^104)
    for f in report.findings   show(io, f)   end
end


"""
`Preflight.analyze(rep::AtomicState.Representation; observables::Array{Symbol,1}=Symbol[])`
    ... makes every check of this module that applies to a restricted-active-space representation and returns the
        `report::Preflight.Report`. Nothing is raised, whatever is found.

        `observables` MUST BE GIVEN HERE AND CANNOT BE DERIVED, which is a design point rather than an omission: a
        `Representation` says what wave function is to be built and nothing about what will be asked of it. The same
        orbitals might be wanted for a transition rate, a hyperfine constant or an energy alone, and those require
        different things. An empty list therefore means "check only what is independent of the observable".
"""
function analyze(rep::AtomicState.Representation; observables::Array{Symbol,1}=Symbol[])
    findings = Finding[]
    label    = rep.name * ", " * string(length(rep.refConfigs)) * " reference configuration(s)"

    if  rep.repType isa AtomicState.RasExpansion
        ras = rep.repType
        append!(findings, checkGrid(rep.grid, rep.refConfigs, rep.nuclearModel))
        append!(findings, checkRoute(ras))
        for obs in observables   append!(findings, checkCsfSpace(ras, rep.refConfigs, obs))   end
    else
        push!(findings, Finding(:note, :csfSpace, "",
              "this representation is not a RAS expansion, so only the observable-independent checks were made.",
              ""))
    end
    append!(findings, checkThreading())

    return( Report(label, findings, Dict{String,Int64}(), costNote()) )
end


"""
`Preflight.analyze(comp::Atomic.Computation; observables::Array{Symbol,1}=Symbol[])`
    ... makes every check that applies to an atomic computation and returns the `report::Preflight.Report`.

        Here the observables CAN be derived, from `comp.propertySettings`, and the derived list is used wherever the
        caller gives none. A caller who intends to take an expectation value that the settings do not mention must
        still name it, since nothing in the computation announces that intention.
"""
function analyze(comp::Atomic.Computation; observables::Array{Symbol,1}=Symbol[])
    obs      = isempty(observables) ? observablesOf(comp.propertySettings) : observables
    findings = Finding[]
    append!(findings, checkGrid(comp.grid, comp.configs, comp.nuclearModel))

    push!(findings, Finding(:note, :csfSpace, "",
          "a configuration list was given directly rather than a RAS ladder, so the CSF-class checks compare the " *
          "requested observables against the configurations themselves.", ""))
    for o in obs
        haskey(OBSERVABLE_REQUIREMENTS, o)  ||  continue
        req = OBSERVABLE_REQUIREMENTS[o]
        push!(findings, Finding(:note, :csfSpace, req.observable,
              "weights " * req.radialWeight * (req.needsSmallComponent ?
              " and couples P to Q, so it vanishes in the non-relativistic limit." : "."), ""))
    end
    append!(findings, checkThreading())

    return( Report(comp.name, findings, Dict{String,Int64}(), costNote()) )
end


"""
`Preflight.checkCsfSpace(ras::AtomicState.RasExpansion, refConfigs::Array{Configuration,1}, observable::Symbol)`
    ... compares what the RAS ladder contains against what one observable requires OF THIS STATE, and returns the
        `findings::Array{Preflight.Finding,1}` -- empty where the space carries every class the observable needs.

        THE REQUIREMENT DEPENDS ON THE STATE AND NOT ON THE OBSERVABLE ALONE, which is the one piece of real physics
        in this module. The magnetic hyperfine constants of a level with an UNPAIRED s electron are contact
        observables and cannot be had without core-valence s excitations; the same constants of a level whose open
        shell is an f hole are a property of that f shell, and the core-valence s class is then a refinement rather
        than a precondition. Deciding this from the reference configuration is exactly what
        `Hfs.warnOnContactObservable` does after the fact, and it is decidable before the run.
"""
function checkCsfSpace(ras::AtomicState.RasExpansion, refConfigs::Array{Configuration,1}, observable::Symbol)
    findings = Finding[]
    if  !haskey(OBSERVABLE_REQUIREMENTS, observable)
        push!(findings, Finding(:note, :csfSpace, string(observable),
              "this observable is not in Preflight.OBSERVABLE_REQUIREMENTS, so nothing is known about what it " *
              "needs.  Add it to the table rather than guessing here.", ""))
        return( findings )
    end
    req       = OBSERVABLE_REQUIREMENTS[observable]
    essential = essentialClassesFor(observable, refConfigs)
    refining  = setdiff(union(req.essentialClasses, req.refiningClasses), essential)

    for class in essential
        hasCsfClass(ras, refConfigs, class)  &&  continue
        push!(findings, Finding(:unconstrained, :csfSpace, req.observable,
              "the space carries no " * describeCsfClass(class) * ", and " * whyCsfClass(class) * ".  " *
              uppercasefirst(req.reason) * ", so this is a precondition and not a refinement.",
              remedyForCsfClass(class)))
    end
    for class in refining
        hasCsfClass(ras, refConfigs, class)  &&  continue
        push!(findings, Finding(:degraded, :csfSpace, req.observable,
              "the space carries no " * describeCsfClass(class) * ".  The observable is still constrained -- " *
              req.reason * " -- but this class would refine it.", remedyForCsfClass(class)))
    end

    return( findings )
end


"""
`Preflight.checkGrid(grid::Radial.Grid, refConfigs::Array{Configuration,1}, nm::Nuclear.Model)`
    ... compares the radial box actually given against the one `Basics.recommendedGrid` would build for the same
        reference configurations, and returns the `findings::Array{Preflight.Finding,1}`.

        THE YARDSTICK IS THE RECOMMENDED GRID AND NOT A HYDROGENIC TURNING POINT, and that is a correction rather
        than a preference. The first version of this check used Rule 12's turning point with Z_eff = Z - N + 1, and
        on Yb+ it put the 6s turning point at 36 a.u. and declared a perfectly good 52.7 a.u. box too small --
        while the measured mean radius of that 6s orbital is 3.5 a.u. A hydrogenic estimate is simply not valid for
        a PENETRATING valence orbital of a heavy ion, which sees far more than the asymptotic charge over most of
        its extent. The same over-warning is recorded in the settled register for every detector built on the ratio
        of an orbital's extent to the box: such a ratio measures where a tail EXISTS and not where its WEIGHT is.

        So this check now compares two boxes built by the same calibrated routine, which is a comparison it can
        actually make, and it reports only a box that differs from the recommendation by more than a factor of two.
"""
function checkGrid(grid::Radial.Grid, refConfigs::Array{Configuration,1}, nm::Nuclear.Model)
    findings = Finding[]
    rMax     = grid.r[end]
    wanted   = 0.0
    try
        wanted = Basics.recommendedGrid(refConfigs, nm; printout=false).r[end]
    catch
        return( findings )
    end
    sBox = @sprintf("%.1f", rMax);    sWant = @sprintf("%.1f", wanted);    sRatio = @sprintf("%.2f", rMax/wanted)

    if      rMax < 0.5 * wanted
        push!(findings, Finding(:degraded, :grid, "",
              "the box reaches " * sBox * " a.u. where the recommended grid for these configurations reaches " *
              sWant * " a.u., a factor of " * sRatio * ".  A box that is too small corrupts PROGRESSIVELY, so the " *
              "artifact grows as the correlation space grows and carries the same signature as a physical trend.",
              "use Basics.recommendedGrid for these reference configurations, or say why a smaller box is wanted"))
    elseif  rMax > 2.0 * wanted
        push!(findings, Finding(:degraded, :grid, "",
              "the box reaches " * sBox * " a.u. where the recommended grid reaches " * sWant * " a.u., a factor " *
              "of " * sRatio * ".  The spline count is fixed, so a box this wide spends splines on empty space " *
              "and misrepresents the high-n members of each symmetry -- often by returning a DIFFERENT state " *
              "rather than an inaccurate one.",
              "use Basics.recommendedGrid unless the extra range is needed for a continuum orbital"))
    end

    return( findings )
end


"""
`Preflight.checkRoute(ras::AtomicState.RasExpansion)`
    ... checks the self-consistent-field route against the ladder it will be applied to, and returns the
        `findings::Array{Preflight.Finding,1}`. Only settings that cannot do what they appear to do are reported.
"""
function checkRoute(ras::AtomicState.RasExpansion)
    findings = Finding[]
    route    = ras.settings.scfRoute
    if  route isa Basics.RotationRoute  &&  Basics.incrementTolerance(route) > 0.  &&  route.maxIterations < 6
        sTol = @sprintf("%.3f", Basics.incrementTolerance(route))
        push!(findings, Finding(:degraded, :route, "",
              "incrementTolerance = " * sTol * " was requested with maxIterations = " * string(route.maxIterations) *
              ", and the doubling test cannot fire before iteration 6 -- at a budget of 3 a reference step was " *
              "measured to return an energy BELOW its own converged value, so a test allowed to fire that early " *
              "stops on noise.  The tolerance is therefore a silent no-op and the run will stop on the ceiling.",
              "raise maxIterations to at least 12, so that the test has a baseline to compare against"))
    end
    if  ras.settings.freezing == :allFree  &&  length(ras.steps) > 2
        push!(findings, Finding(:note, :freezing, "",
              "freezing = :allFree releases every orbital in every layer, so each layer re-optimizes what the " *
              "previous ones settled.  That is a legitimate choice and it makes the layer increments no longer " *
              "additive, so an increment sequence from this ladder cannot be compared with one from a frozen " *
              "ladder.", ""))
    end

    return( findings )
end


"""
`Preflight.checkThreading()`
    ... reports the one threading hazard that is known and decidable without a run, and returns the
        `findings::Array{Preflight.Finding,1}` -- empty on a single-threaded session.
"""
function checkThreading()
    findings = Finding[]
    if  Threads.nthreads() > 1
        push!(findings, Finding(:note, :threading, "",
              "this session runs with " * string(Threads.nthreads()) * " threads, and the memoisation caches in " *
              "SpinAngular are not thread-safe -- a real eight-thread run has died with a concurrent-write " *
              "assertion.  This is a crash and not a wrong number, so a run that COMPLETES is unaffected by it.",
              "run single-threaded if the computation must not be interrupted"))
    end

    return( findings )
end


"""
`Preflight.costNote()`
    ... says what is known about the cost of a computation without estimating one, and returns a `note::String`.
        The estimate itself is made by the optimized-level solver at the start of each layer, where the CSF count
        is known; repeating it here from the input alone would be a second, worse estimate of the same quantity.
"""
function costNote()
    return( "the optimized-level solver prints a measured cost and time estimate at the start of every layer, " *
            "once the CSF count is known;  it is anchored on a 69-electron ion and is good to a factor of three " *
            "either way across elements." )
end


"""
`Preflight.describeCsfClass(class::Symbol)`
    ... names a class of configuration state function in words a reader of a paper would use, and returns the
        `description::String`.
"""
function describeCsfClass(class::Symbol)
    class == :coreValenceS            &&  return( "core-valence single excitations out of an s subshell" )
    class == :valenceCorrelation      &&  return( "double excitations among the valence shells" )
    class == :fSymmetryVirtuals       &&  return( "virtual orbitals of f or g symmetry" )
    class == :bothStatesOneOrbitalSet &&  return( "statement that both states share one orthonormal orbital set" )

    return( "class :" * string(class) )
end


"""
`Preflight.essentialClassesFor(observable::Symbol, refConfigs::Array{Configuration,1})`
    ... returns the `classes::Array{Symbol,1}` without which THIS state's value of the observable is not constrained
        at all. It is the table's `essentialClasses` for every observable but the two magnetic hyperfine constants,
        whose contact character depends on whether the level has an unpaired s electron.
"""
function essentialClassesFor(observable::Symbol, refConfigs::Array{Configuration,1})
    haskey(OBSERVABLE_REQUIREMENTS, observable)  ||  return( Symbol[] )
    classes = copy(OBSERVABLE_REQUIREMENTS[observable].essentialClasses)
    if  observable in (:hyperfineA, :hyperfineC)  &&  !isContactLike(refConfigs)
        classes = filter(c -> c != :coreValenceS, classes)
    end

    return( classes )
end


"""
`Preflight.hasCsfClass(ras::AtomicState.RasExpansion, refConfigs::Array{Configuration,1}, class::Symbol)`
    ... returns `true::Bool` where the ladder's steps contain the named class of configuration state function, and
        `false` where they do not. A class this function does not recognize returns `false`, so that an unknown
        requirement is reported as missing rather than silently satisfied.

        THE CORE-VALENCE s TEST MUST EXCLUDE THE VALENCE s ITSELF, and the first version of this module did not --
        which made it answer the one question the whole module was built for wrongly. Asked whether a space
        correlating 4f and 6s carried core-valence s excitations, it said yes, because `6s` has l = 0 and appeared
        in the excitation list. An s shell qualifies only if it is CLOSED in the reference configuration and lies
        BELOW the outermost occupied s shell, i.e. only if it is genuinely a core shell that has been released.
"""
function hasCsfClass(ras::AtomicState.RasExpansion, refConfigs::Array{Configuration,1}, class::Symbol)
    if      class == :coreValenceS
        nValence = outermostSShell(refConfigs)
        for step in ras.steps
            isempty(step.seTo)  &&  isempty(step.deTo)  &&  continue
            for sh in union(step.seFrom, step.deFrom)
                sh.l == 0  &&  sh.n < nValence  &&  isClosedIn(sh, refConfigs)  &&  return( true )
            end
        end
    elseif  class == :valenceCorrelation
        for step in ras.steps
            !isempty(step.deTo)  &&  return( true )
        end
    elseif  class == :fSymmetryVirtuals
        for step in ras.steps
            any(sh -> sh.l >= 3, step.seTo)  &&  return( true )
            any(sh -> sh.l >= 3, step.deTo)  &&  return( true )
        end
    elseif  class == :bothStatesOneOrbitalSet
        return( false )
    end

    return( false )
end


"""
`Preflight.isClosedIn(shell::Shell, refConfigs::Array{Configuration,1})`
    ... returns `true::Bool` where the given shell is fully occupied in every reference configuration, which is the
        test of whether releasing it amounts to a core excitation.
"""
function isClosedIn(shell::Shell, refConfigs::Array{Configuration,1})
    full = 2*(2*shell.l + 1)
    for conf in refConfigs
        get(conf.shells, shell, 0) == full  ||  return( false )
    end

    return( true )
end


"""
`Preflight.isContactLike(refConfigs::Array{Configuration,1})`
    ... returns `true::Bool` where some reference configuration carries an s subshell of ODD occupation, so that the
        level has an unpaired s electron and its magnetic hyperfine constants are contact observables.
"""
function isContactLike(refConfigs::Array{Configuration,1})
    for conf in refConfigs
        for (sh, occ) in conf.shells
            sh.l == 0  &&  isodd(occ)  &&  return( true )
        end
    end

    return( false )
end


"""
`Preflight.observablesOf(settings::Array{Basics.AbstractPropertySettings,1})`
    ... maps the property settings of an atomic computation onto the keys of `OBSERVABLE_REQUIREMENTS`, and returns
        the `observables::Array{Symbol,1}`. A settings type with no entry in the table contributes nothing, which is
        why a caller who wants a quantity the table does not know must name it.
"""
function observablesOf(settings::Array{Basics.AbstractPropertySettings,1})
    observables = Symbol[]
    for set in settings
        name = string(typeof(set))
        occursin("Hfs",         name)  &&  append!(observables, [:hyperfineA, :hyperfineB, :hyperfineC])
        occursin("LandeZeeman", name)  &&  push!(observables, :landeG)
        occursin("IsotopeShift", name) &&  push!(observables, :isotopeShift)
    end

    return( unique(observables) )
end


"""
`Preflight.outermostSShell(refConfigs::Array{Configuration,1})`
    ... returns the principal quantum number `n::Int64` of the outermost occupied s shell over all reference
        configurations, or zero where none is occupied. It is the boundary between a core s shell and the valence one.
"""
function outermostSShell(refConfigs::Array{Configuration,1})
    n = 0
    for conf in refConfigs
        for (sh, occ) in conf.shells
            sh.l == 0  &&  occ > 0  &&  sh.n > n  &&  (n = sh.n)
        end
    end

    return( n )
end


"""
`Preflight.remedyForCsfClass(class::Symbol)`
    ... returns the `remedy::String`, a concrete change to the input that would supply the named class, or an empty
        string where the class cannot be supplied by a change to the ladder alone.
"""
function remedyForCsfClass(class::Symbol)
    class == :coreValenceS            &&  return( "release the outermost core s subshell into `fromShells` and add " *
                                                  "a singles-only layer of one or two further s shells" )
    class == :valenceCorrelation      &&  return( "add a layer with de = true over the open shells" )
    class == :fSymmetryVirtuals       &&  return( "add virtual shells of f and g symmetry to a layer" )
    class == :bothStatesOneOrbitalSet &&  return( "compute both states in one orbital set, or apply the " *
                                                  "biorthogonal transformation before the matrix element" )

    return( "" )
end


"""
`Preflight.whyCsfClass(class::Symbol)`
    ... says what a class of configuration state function does for an observable, and returns the `why::String`.
        The reason belongs to the CLASS and not to the observable: the first version of this module pasted the
        observable's reason onto whichever class happened to be missing, which produced sentences that did not
        follow from one another.
"""
function whyCsfClass(class::Symbol)
    class == :coreValenceS            &&  return( "only an excitation out of a closed core s shell can relax the " *
                                                  "density in the region closest to the nucleus" )
    class == :valenceCorrelation      &&  return( "the open shells are otherwise treated as a single " *
                                                  "configuration" )
    class == :fSymmetryVirtuals       &&  return( "two f electrons couple to L = 0 ... 6, so d, f, g and h " *
                                                  "virtuals are each a separate correlation channel" )
    class == :bothStatesOneOrbitalSet &&  return( "a matrix element between two independently optimized orbital " *
                                                  "sets is not a matrix element between orthonormal states" )

    return( "the reason for class :" * string(class) * " is not recorded" )
end


end # module
