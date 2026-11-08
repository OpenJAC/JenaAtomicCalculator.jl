
"""
`module  JAC.AtomicState`  
	... a submodel of JAC that contains all methods to set-up and process atomic representations.
"""
module AtomicState


# using Interact
using  Printf, ..Basics, ..Hamiltonian, ..ManyElectron, ..Nuclear, ..Radial, ..SelfConsistent, ..SpinAngular, ..ResultKeys

export  MeanFieldSettings, MeanFieldBasis, MeanFieldMultiplet, OneElectronSettings, OneElectronSpectrum, CiSettings, CiExpansion,
        RasSettings, RasStep, RasLayer, RasExpansion, GreenSettings, GreenChannel, GreenExpansion, Representation,
        RasAutoSettings, FixedLayerCount, IncrementBelow, IncrementRatioBelow, CostCeiling, CombinedPolicy,
        RasCompleteShellPolicy, RasLimitedAngularMomentumPolicy, RasSpectrumOrderedShellPolicy


"""
`abstract type AtomicState.AbstractRepresentationType` 
    ... defines an abstract type and a number of data types to work with and distinguish different atomic representations; see also:
    
    + struct MeanFieldBasis       ... to represent (and generate) a mean-field basis and, especially, a set of 
                                        (mean-field) orbitals.
    + struct MeanFieldMultiplet   ... to represent (and generate) a mean-field multiplet, based on (mean-field) orbitals.
    + struct OneElectronSpectrum  ... to represent (and generate) a one-electron spectrum for the mean-field 
                                        potential of refConfigs.
    + struct CiExpansion          ... to represent (and generate) a configuration-interaction representation.
    + struct RasExpansion         ... to represent (and generate) a restricted active-space representation.
    + struct GreenExpansion       ... to represent (and generate) an approximate (many-electron) Green functions.
"""
abstract type  AbstractRepresentationType                           end


### Mean-field basis ############################################################################################
"""
`struct  AtomicState.MeanFieldSettings`  
    ... a struct for defining the settings for a mean-field basis (orbital) representation.

    + scField           ::AbstractScField   
        ... Specify the (mean) self-consistent field as DFSField() or HSField(); note that not all AbstractScField's 
            allowed here.
    + gridStopper       ::Bool
        ... False to let `Bsplines.checkGridRepresentation` WARN instead of raising, so that deliberately
            box-confined high-n orbitals can be generated.  **Added 15-Oct-2026, and it was a real gap rather
            than a convenience:** this struct carried `scField` alone, so the whole
            `Representation`/`MeanFieldMultiplet` route could not express a PSEUDO-CONTINUUM intermediate set
            at all -- and that is what kept `example-Dh.jl` branch (b) on a bound-only sum for months, giving
            10.104 /s where the exact two-photon rate is 8.22906.  The recipe that converges it (n_max ~ 100,
            the grid guard relaxed, the box held fixed) had to be driven through an `Atomic.Computation`
            instead.  The default is TRUE, i.e. the guard raises as before, so no existing caller changes;
            set it false only where the high-n states are MEANT to be box-confined stand-ins for the
            continuum, and say so at the call site.
"""
struct  MeanFieldSettings
    scField              ::AbstractScField
    gridStopper           ::Bool
end


"""
`AtomicState.MeanFieldSettings()`  ... constructor for setting the default values.
        A `meanFieldSettings::AtomicState.MeanFieldSettings` is returned.
"""
function MeanFieldSettings()
    MeanFieldSettings(DFSField(), true)
end


"""
`AtomicState.MeanFieldSettings(scField::AbstractScField)`
    ... constructor from a field alone, keeping `gridStopper = true`.  It exists so that the one-argument form
        every existing caller uses goes on meaning what it meant;  a `settings::MeanFieldSettings` is returned.
"""
function MeanFieldSettings(scField::AbstractScField)
    MeanFieldSettings(scField, true)
end


"""
`AtomicState.MeanFieldSettings(settings::AtomicState.MeanFieldSettings;
                               scField::Union{Nothing,AbstractScField}=nothing,              gridStopper::Union{Nothing,Bool}=nothing)`
    ... the standard JAC keyword copy-constructor; a `settings::MeanFieldSettings` is returned.
"""
function MeanFieldSettings(settings::AtomicState.MeanFieldSettings;
    scField::Union{Nothing,AbstractScField}=nothing,              gridStopper::Union{Nothing,Bool}=nothing)

    if  isnothing(scField)          scFieldx        = settings.scField          else   scFieldx        = scField          end
    if  isnothing(gridStopper)      gridStopperx    = settings.gridStopper      else   gridStopperx    = gridStopper      end

    MeanFieldSettings( scFieldx, gridStopperx )
end


function Base.show(io::IO, settings::MeanFieldSettings)
        println(io, "scField:                  $(settings.scField)  ")
        println(io, "gridStopper:              $(settings.gridStopper)  ")
end


"""
`struct  AtomicState.MeanFieldBasis  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a mean-field orbital basis.

    + settings         ::AtomicState.MeanFieldSettings      ... Settings for the given mean-field orbital basis
"""
struct   MeanFieldBasis  <:  AbstractRepresentationType
    settings           ::AtomicState.MeanFieldSettings
end


"""
`AtomicState.MeanFieldBasis()`  
    ... constructor for an 'empty' instance of a variable::AtomicState.MeanFieldBasis, i.e. a mean-field basis built with
        the default `MeanFieldSettings()`.  It completes the set of zero-argument forms of the representation types, so that
        each of them can be written down before its settings are chosen.  A basis::AtomicState.MeanFieldBasis is returned.
"""
function MeanFieldBasis()
    MeanFieldBasis( AtomicState.MeanFieldSettings() )
end


function Base.string(basis::MeanFieldBasis)
    if !(basis.settings.scField  in [DFSField(), HSField()])
        error("A MeanFieldBasis presently supports only a DFSField() or HSField() but received:  $(basis.settings.scField)")
    end
    sa = "Mean-field orbital basis for a $(basis.settings.scField) SCF field:"
    return( sa )
end


function Base.show(io::IO, basis::MeanFieldBasis)
    sa = Base.string(basis);       print(io, sa, "\n")
end


"""
`struct  AtomicState.MeanFieldMultiplet  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a mean-field orbital basis and multiplet.

    + settings         ::AtomicState.MeanFieldSettings      ... Settings for the given mean-field orbital basis and multiplet.
"""
struct   MeanFieldMultiplet  <:  AbstractRepresentationType
    settings           ::AtomicState.MeanFieldSettings
end


"""
`AtomicState.MeanFieldMultiplet()`  
    ... constructor for an 'empty' instance of a variable::AtomicState.MeanFieldMultiplet, i.e. a mean-field multiplet
        built with the default `MeanFieldSettings()`.  It completes the set of zero-argument forms of the representation
        types, so that each of them can be written down before its settings are chosen.  A
        multiplet::AtomicState.MeanFieldMultiplet is returned.
"""
function MeanFieldMultiplet()
    MeanFieldMultiplet( AtomicState.MeanFieldSettings() )
end


function Base.string(basis::MeanFieldMultiplet)
    sa = "Mean-field multiplet for a $(basis.settings.scField) SCF field:"
    return( sa )
end


function Base.show(io::IO, basis::MeanFieldMultiplet)
    sa = Base.string(basis);       print(io, sa, "\n")
end



### Spectrum: One-electron spectrum for given (reference) configurations ####################################################
"""
`struct  AtomicState.OneElectronSettings`  
    ... a struct for defining the settings for a one-electron spectrum.

    + nMax                 ::Int64            ... maximum principal quantum number.
    + lValues              ::Array{Int64,1}   ... l-values (partial waves) for which orbitals are to be generated.
    + levelSelectionMean   ::LevelSelection   ... Level(s) of the mean-field to generate the underlying atomic potential
"""
struct  OneElectronSettings
    nMax                   ::Int64
    lValues                ::Array{Int64,1}
    levelSelectionMean     ::LevelSelection
end


"""
`AtomicState.OneElectronSettings()`  ... constructor for setting the default values.
        A `oneElectronSettings::AtomicState.OneElectronSettings` is returned.
"""
function OneElectronSettings()
    OneElectronSettings(10, [0, 1], LevelSelection(true, indices=[1]))
end


"""
`AtomicState.OneElectronSettings(settings::AtomicState.OneElectronSettings;
    nMax::Union{Nothing,Int64}=nothing,                             lValues::Union{Nothing,Array{Int64,1}}=nothing,
    levelSelectionMean::Union{Nothing,LevelSelection}=nothing)`

    ... constructor for modifying the given OneElectronSettings by 'overwriting' the explicitly selected parameters.
        A `oneElectronSettings::AtomicState.OneElectronSettings` is returned.
"""
function OneElectronSettings(settings::AtomicState.OneElectronSettings;
    nMax::Union{Nothing,Int64}=nothing,                             lValues::Union{Nothing,Array{Int64,1}}=nothing,
    levelSelectionMean::Union{Nothing,LevelSelection}=nothing)
    
    if  isnothing(nMax)                 nMaxx               = settings.nMax                else   nMaxx               = nMax               end 
    if  isnothing(lValues)              lValuesx            = settings.lValues             else   lValuesx            = lValues            end 
    if  isnothing(levelSelectionMean)   levelSelectionMeanx = settings.levelSelectionMean  else   levelSelectionMeanx = levelSelectionMean end 
    
    OneElectronSettings( nMaxx, lValuesx, levelSelectionMeanx)
end


function Base.show(io::IO, settings::OneElectronSettings)
        println(io, "nMax:                 $(settings.nMax)  ")
        println(io, "lValues:              $(settings.lValues)  ")
        println(io, "levelSelectionMean:   $(settings.levelSelectionMean)  ")
end


"""
`struct  AtomicState.OneElectronSpectrum  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a one-electron spectrum for the mean-field potential of given levels 
        from a set of reference configurations.

    + settings         ::AtomicState.OneElectronSettings         ... Settings for the given OneElectronSpectrum
"""
struct  OneElectronSpectrum     <:  AbstractRepresentationType
    settings           ::AtomicState.OneElectronSettings
end


"""
`AtomicState.OneElectronSpectrum()`  
    ... constructor for an 'empty' instance of a variable::AtomicState.OneElectronSpectrum, i.e. a one-electron spectrum
        built with the default `OneElectronSettings()`.  It completes the set of zero-argument forms of the representation
        types, so that each of them can be written down before its settings are chosen.  A
        spectrum::AtomicState.OneElectronSpectrum is returned.
"""
function OneElectronSpectrum()
    OneElectronSpectrum( AtomicState.OneElectronSettings() )
end


function Base.string(spectrum::OneElectronSpectrum)
    sa = "One-electron spectrum for given atomic potential:"
    return( sa )
end


function Base.show(io::IO, spectrum::OneElectronSpectrum)
    sa = Base.string(spectrum);       print(io, sa, "\n")
    println(io, "... and the current settings:")
    println(io, "$(spectrum.settings)  ")
end



### RAS: Restricted-Active-Space expansions ##################################################################################
"""
`struct  AtomicState.RasSettings`  
    ... a struct for defining the settings for a restricted active-space computations.

    + levelsScf            ::Array{Int64,1}         ... Levels on which the optimization is carried out, BY INDEX.
        Read only when `levelSelectionCI` is INACTIVE; when both are given, `levelSelectionCI` wins and the RAS
        driver says so. Until 01-Sep-2026 this field was declared, documented and printed but NEVER READ (priority
        item 22), so a user who wrote `RasSettings([1], ...)` to optimize a layer on the ground level alone
        silently got an EOL functional over every level `levelSelectionCI` had selected -- thirteen of them in one
        measured case -- while the printout confirmed the wrong belief.
        PREFER `levelSelectionCI` with `configurations=`: selecting an EOL target by INDEX is unstable, because
        the indices refer to the energy-sorted multiplet and a correlation configuration sinking below the
        reference silently changes which levels are optimized.
    + scfRoute             ::Basics.AbstractScfRoute   ... the SCF route by which each RAS step is solved, carrying its own
        iteration budget and whatever else that solver needs. **This replaced `maxIterationsScf::Int64` on 04-Oct-2026 and the
        change is not cosmetic.** The old field was converted, one line deep in the driver, into
        `Basics.RotationRoute(maxIterationsScf)` -- so a RAS ladder could name ONLY the rotation route and ONLY its iteration
        count, while `RotationRoute` also carries `nVirtual` and `stepping`, `FockRoute` carries `unscaledOffDiagonal`, and the
        route hierarchy exists precisely so that a computation may be repeated on a more expensive route when a cheaper one does
        not converge. Four of the six things a route can say were unreachable from a RAS computation. Naming the route here makes
        all of them reachable and gives the quantity ONE spelling, the same as `AsfSettings.scfRoute`.
    + accuracyScf          ::Float64                ... convergence criterion for the SCF field.  This stays a settings field
        rather than moving into the route, exactly as in `AsfSettings`, and the split is deliberate: an iteration BUDGET is a
        property of the solver -- a mean field converges in a handful of iterations where the rotation route needs thousands --
        whereas a convergence TOLERANCE is a property of the answer the caller wants, and is the same question whichever solver
        is asked.
    + eeInteractionCI      ::AbstractEeInteraction  ... logical flag to include Breit interactions.
    + levelSelectionCI     ::LevelSelection         ... Specifies the selected levels, if any; also the EOL target
        of every step when active, and then takes precedence over `levelsScf`.
    + freezing             ::Symbol                 ... WHICH ORBITALS A LAYER MAY STILL MOVE, one of
        `:previousLayersFrozen` (the default), `:coreFrozen` or `:allFree`; an unrecognized symbol RAISES.  These are
        Rule 21's three freezing policies, and naming the policy is what makes a layer increment comparable between
        runs -- so it belongs in the settings that are printed with the computation rather than in a script.
        * `:previousLayersFrozen` ... every shell an earlier layer optimized is held, and only the shells this layer
          introduces are varied.  GRASP's default, and the safe one: the layer's space contains the previous layer's,
          and with the earlier orbitals held it cannot help but descend.
        * `:coreFrozen` ... only the `coreShells` of the expansion are held; everything else is released, so the
          valence shells may relax against the correlation this layer adds.
        * `:allFree` ... nothing is held.  The richest policy and the most expensive.
        **THE RELEASED POLICIES RUN IN TWO PHASES, and that is what makes them trustworthy.**  A layer first WARMS its
        new shells with the earlier ones frozen, and only then RELEASES them, starting the optimized-level solver from
        the warm orbitals through `ManyElectron.StartFromPrevious`.  Because the warming phase leaves a set that COVERS
        the basis, the average-level pre-pass is then skipped, and so does not re-derive the inherited orbitals for the
        configuration AVERAGE rather than for the level energy.  (That skip follows from the COVERAGE and not from
        `StartFromPrevious` alone;  a handed-in set that leaves subshells uncovered would leave them hydrogenic, so
        `SelfConsistent.averageLevelPrePass` pre-passes those alone with the handed-in ones frozen.)
        Without that two-step entry a released layer begins far ABOVE the layer it sits on: measured on 43Ca+ [Ar] 4s
        with 3s -> 5s,6s singles, 36.1 mHa above, which is more than the whole correlation effect being measured, and
        the resulting layer increment had a BUDGET-DEPENDENT SIGN (+11.2, +5.9, +1.0 mHa at budgets 20, 60, 120).
        Repaired, the same layer gives -0.013, -0.053, -0.063 mHa -- monotone, correctly signed, and saturating.
        **WHAT THE RELEASE IS WORTH, measured on three ions with the same two-layer structure and `freezing` as the
        only variable** -- a reference layer, then ONE correlation layer with singles and doubles, every box checked
        adequate and every correlation orbital contracted into the valence region, so no grid confound:

            ion      reference        CSFs   frozen incr.   allFree incr.   deeper   doubling test   cost
            Ca+      [Ar] 4s             5   -0.0234 mHa    -0.0604 mHa     +158 %   29.5 % / 7.3 %  ~2.0x
            C II     1s^2 2s^2 2p       41  -77.1587 mHa   -80.4217 mHa      +4.2 %  0.0001%/ 0.31%   1.65x
            Ti III   [Ar] 3d^2        2264 -150.8840 mHa  -161.8682 mHa      +7.3 %  0.025 % / 0.25%  1.83x

        The release is deeper on every case above -- but **READ THE NEXT PARAGRAPH BEFORE CHOOSING IT, because that
        advantage is almost entirely the release buying back something the REFERENCE layer was forbidden to do, and
        it is far better not to forbid it.**
        **WHERE THE CORE IS ALLOWED TO RELAX MATTERS MORE THAN THIS POLICY DOES, and earlier is better.** Isolated
        07-Oct-2026 on Ti III [Ar] 3d^2 + 4s4p4d(SD), one code path, adaptive stop, `coreShells` the only variable:

            layer 1            step 1         :previousLayersFrozen    :allFree        policy gap
            core FROZEN      -852.160048      -852.310901 (305 s)    -852.321827 (526 s)   10.93 mHa
            core FREE        -852.176877      -852.327971 (321 s)    -852.328019 (544 s)    0.048 mHa

        Letting the core relax at the REFERENCE layer is worth **16.8 mHa** on step 1, and the best total improves by
        **6.2 mHa** -- so no later freezing policy recovers what a frozen core costs: the best frozen-core run still
        ends 6.2 mHa ABOVE the worst core-free one. And once the core has relaxed there, **this policy is worth
        0.048 mHa, a factor 228 less than it appeared**, at which point `:previousLayersFrozen` is the right choice
        because it is 1.7x faster.
        **SO `:allFree` IS A REMEDY FOR A FROZEN-CORE REFERENCE, NOT A GENERAL IMPROVEMENT.** Declare `coreShells`
        only when cost forces it; if it is forced, `:allFree` wins back about 11 mHa of the ~17 mHa given away.
        The same ordering holds on C II, less sharply: the gap falls from 3.26 mHa to 0.61 mHa when its 1s is
        released, a 1s having little to relax.  It costs 1.65x to 1.83x rather than
        the 2x two phases suggest, because the warm-up leaves the release a good start.  Cross-checked under both
        `stepping = :curvature` and `:newton`, so it is a property of the handover and not of the solver.
        **AND A CLAIM MADE ON THE Ca+ CASE ALONE IS WITHDRAWN HERE:** that the released policy SATURATES in the
        iteration budget while the frozen one does not.  That is true only of Ca+.  On both open-shell ions BOTH
        policies converge and the FROZEN one converges better.  **The argument for the release is its DEPTH, not its
        convergence.**
        **WHAT MAKES A LAYER HARD TO CONVERGE IS THE LAYER'S OWN SHELL STRUCTURE, not the ion's.**  Ca+'s layer adds
        5s and 6s -- high-n orbitals of the SAME symmetry as the occupied 4s, diffuse and nearly degenerate with it;
        C II adds 3s,3p,3d and Ti III adds 4s,4p,4d, compact orbitals mostly in l channels the reference does not
        occupy.  Measured, the drift of the increment when the budget doubles: 0.0069 mHa for Ca+ against 0.0001 mHa
        for C II, i.e. 69x more drift on a 3000x smaller increment.  A layer of high-n orbitals in an
        already-occupied symmetry is the hard case AND earns little;  a layer opening new l channels is the easy case
        AND earns the real correlation energy.  Budget accordingly, and test it by DOUBLING the budget and comparing
        the layer INCREMENT -- under about 10 % the layer is converged.  That test replaces the "about twelve
        iterations" rule of thumb, which is strongly Z-dependent (45:1 at Z = 92 against 7:1 at Z = 8) and so cannot
        be carried between elements.
"""
struct  RasSettings
    levelsScf              ::Array{Int64,1}
    scfRoute               ::Basics.AbstractScfRoute
    accuracyScf            ::Float64 
    eeInteractionCI        ::AbstractEeInteraction 
    levelSelectionCI       ::LevelSelection
    freezing               ::Symbol

    # The policy is validated HERE and not in the driver, so that a mistyped name is named where the caller wrote it
    # rather than silently selecting a different physics -- the same reasoning as for RotationRoute's stepping.
    function RasSettings(levelsScf::Array{Int64,1}, scfRoute::Basics.AbstractScfRoute, accuracyScf::Float64,
                         eeInteractionCI::AbstractEeInteraction, levelSelectionCI::LevelSelection, freezing::Symbol)
        if  !(freezing in (:previousLayersFrozen, :coreFrozen, :allFree))
            error("AtomicState.RasSettings: freezing = :$freezing is not recognized.  Use :previousLayersFrozen " *
                  "(the default: only this layer's new shells are varied), :coreFrozen (the coreShells are held, " *
                  "everything else is released) or :allFree (nothing is held).  See the docstring for what each is " *
                  "worth and what it costs.")
        end
        new(levelsScf, scfRoute, accuracyScf, eeInteractionCI, levelSelectionCI, freezing)
    end
end


"""
`AtomicState.RasSettings(levelsScf::Array{Int64,1}, scfRoute::Basics.AbstractScfRoute, accuracyScf::Float64,
                         eeInteractionCI::AbstractEeInteraction, levelSelectionCI::LevelSelection)`  
    ... the five-argument form every caller written before 06-Oct-2026 uses;  it keeps the freezing policy that was
        the only one available then, so no existing computation changes.
"""
function RasSettings(levelsScf::Array{Int64,1}, scfRoute::Basics.AbstractScfRoute, accuracyScf::Float64,
                     eeInteractionCI::AbstractEeInteraction, levelSelectionCI::LevelSelection)
    RasSettings(levelsScf, scfRoute, accuracyScf, eeInteractionCI, levelSelectionCI, :previousLayersFrozen)
end


"""
`AtomicState.RasSettings()`  ... constructor for setting the default values.
        A `rasSettings::AtomicState.RasSettings` is returned.
"""
function RasSettings()
    # RotationRoute(24) is what the driver built from the old maxIterationsScf = 24, so the default is unchanged in behaviour.
    RasSettings(Int64[1], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(), LevelSelection(), :previousLayersFrozen )
end


"""
`AtomicState.RasSettings(set::AtomicState.RasSettings;
    levelsScf::Union{Nothing,Array{Int64,1}}=nothing,                      scfRoute::Union{Nothing,Basics.AbstractScfRoute}=nothing,
    accuracyScf::Union{Nothing,Float64}=nothing,                           eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,
    levelSelectionCI::Union{Nothing,LevelSelection}=nothing,               freezing::Union{Nothing,Symbol}=nothing,
    maxIterationsScf::Union{Nothing,Int64}=nothing)`

    ... constructor for modifying the given AtomicState.RasSettings by 'overwriting' the previously selected parameters; a
        settings::AtomicState.RasSettings is returned.

        IT ALSO REFUSES `maxIterationsScf` BY NAME, and that is the point of having it. `AsfSettings` retired the same field
        into the routes and refuses the keyword with a message naming the replacement, so a caller is told what to write;
        `RasSettings` carried it on as its SECOND POSITIONAL argument, so the one quantity had two spellings and only one of
        them announced itself. A user who learned the new convention from the `AsfSettings` message and wrote it into a
        `RasSettings` call got a `MethodError` about positional arguments instead -- which is what this refusal replaces.
"""
function RasSettings(set::AtomicState.RasSettings;
    levelsScf::Union{Nothing,Array{Int64,1}}=nothing,                      scfRoute::Union{Nothing,Basics.AbstractScfRoute}=nothing,
    accuracyScf::Union{Nothing,Float64}=nothing,                           eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,
    levelSelectionCI::Union{Nothing,LevelSelection}=nothing,               freezing::Union{Nothing,Symbol}=nothing,
    maxIterationsScf::Union{Nothing,Int64}=nothing)

    if  !isnothing(maxIterationsScf)
        error("RasSettings no longer carries maxIterationsScf; the iteration budget belongs to the scf route, as it does in " *
              "AsfSettings.  Write, for instance, scfRoute = Basics.RotationRoute($maxIterationsScf) for an optimized level " *
              "by orbital rotation, Basics.FockRoute($maxIterationsScf) for the Fock route, or " *
              "Basics.MeanFieldRoute($maxIterationsScf) for a local screened potential.  Naming the route also reaches what " *
              "the old field could not: nVirtual and stepping for the rotation route, unscaledOffDiagonal for the Fock route.")
    end
    if  isnothing(levelsScf)          levelsScfx          = set.levelsScf          else  levelsScfx          = levelsScf          end 
    if  isnothing(scfRoute)           scfRoutex           = set.scfRoute           else  scfRoutex           = scfRoute           end 
    if  isnothing(accuracyScf)        accuracyScfx        = set.accuracyScf        else  accuracyScfx        = accuracyScf        end 
    if  isnothing(eeInteractionCI)    eeInteractionCIx    = set.eeInteractionCI    else  eeInteractionCIx    = eeInteractionCI    end 
    if  isnothing(levelSelectionCI)   levelSelectionCIx   = set.levelSelectionCI   else  levelSelectionCIx   = levelSelectionCI   end 
    if  isnothing(freezing)           freezingx           = set.freezing           else  freezingx           = freezing           end 

    RasSettings( levelsScfx, scfRoutex, accuracyScfx, eeInteractionCIx, levelSelectionCIx, freezingx )
end


function Base.show(io::IO, settings::RasSettings)
        println(io, "levelsScf:            $(settings.levelsScf)  " *
                    (settings.levelSelectionCI.active ? "  (NOT in force: levelSelectionCI is active and wins)" :
                                                        "  (in force as the EOL target of each step)"))
        println(io, "scfRoute:             $(settings.scfRoute)   (budget $(Basics.maxIterations(settings.scfRoute)))  ")
        println(io, "accuracyScf:          $(settings.accuracyScf)  ")
        println(io, "eeInteractionCI:      $(settings.eeInteractionCI)  ")
        println(io, "levelSelectionCI:     $(settings.levelSelectionCI)  ")
        println(io, "freezing:             $(settings.freezing)  " *
                    (settings.freezing == :previousLayersFrozen ? "  (only each layer's new shells are varied)" :
                     settings.freezing == :coreFrozen           ? "  (coreShells held, the rest released; two SCF phases per layer)" :
                                                                  "  (nothing held; two SCF phases per layer)"))
end


"""
`struct  AtomicState.RasStep`  
    ... specifies an individual step of a (relativistic) restricted active space computation for a set of levels. This struct 
        comprises all information to generate the orbital basis and to perform the associated SCF and multiplet computations for a 
        selected number of levels.

    + seFrom            ::Array{Shell,1}        ... Single-excitations from shells   [sh_1, sh_2, ...]
    + seTo              ::Array{Shell,1}        ... Single-excitations to shells  [sh_1, sh_2, ...]
    + deFrom            ::Array{Shell,1}        ... Double-excitations from shells   [sh_1, sh_2, ...]
    + deTo              ::Array{Shell,1}        ... Double-excitations to shells  [sh_1, sh_2, ...]
    + teFrom            ::Array{Shell,1}        ... Triple-excitations from shells   [sh_1, sh_2, ...]
    + teTo              ::Array{Shell,1}        ... Triple-excitations to shells  [sh_1, sh_2, ...]
    + qeFrom            ::Array{Shell,1}        ... Quadrupole-excitations from shells   [sh_1, sh_2, ...]
    + qeTo              ::Array{Shell,1}        ... Quadrupole-excitations to shells  [sh_1, sh_2, ...]
    + shapedTo          ::Array{Shell,1}        ... Which of this step's OWN *To shells get a separate,
        independently-contracted `Basics.CoreValenceCorrelationShape` orbital rather than the ordinary
        `Basics.SCFShape` -- see that type's own docstring. Empty (the default) leaves every shell this step
        generates as `SCFShape`, exactly as before this field existed. When set: any OTHER shell this step's
        own *To lists name (a core-valence DOUBLE needs the full valence+new-shell list to generate at all,
        exactly as an ordinary second-order step already requires) must be in `frozenShells`, since this
        step's own SCF re-optimizes only what is unfrozen and nothing may vary except the shaped shell(s);
        and its `treatment` must not be `Basics.Variational()`: the refined, shaped orbital this step produces
        is kept as a separate P/Q junction (`AtomicState.RasLayerVerdict.shapedBasis`) rather than merged into
        the ladder's own running multiplet, since ordinary Slater-Condon CI assumes every subshell of a kappa
        is mutually orthonormal, which a shaped realization of an existing label deliberately is not.
    + frozenShells      ::Array{Shell,1}        ... List of shells that are kept 'frozen' in this step.
    + constraints       ::Array{String,1}       ... List of Strings to define 'constraints/restrictions' 
                                                    to the generated CSF basis.
    + scfRoute          ::Union{Nothing,Basics.AbstractScfRoute}   ... an SCF route for THIS step alone, or `nothing` to inherit
        the ladder's `RasSettings.scfRoute`.  `nothing` is the default and is what nearly every ladder wants.  The override exists
        because a ladder's steps are not alike: measured, a REFERENCE layer is converged by about twelve iterations while a
        CORRELATION layer never converges at all -- it descends indefinitely with increments that grow -- so one budget for the
        whole ladder cannot express both, and spending the correlation layer's budget on the reference one buys nothing.  This
        mirrors `treatment`, which is likewise a per-step policy rather than a ladder-wide one.
    + treatment         ::Basics.AbstractQTreatment  ... how the configurations this step generates are to be treated:
                                                    Basics.Variational() puts every one of them into the CI, which is
                                                    what a step has always done and remains the default, while
                                                    Basics.SecondOrder(promoteAbove, discardBelow) ranks each one by
                                                    the fraction of the wave function it carries and promotes, folds
                                                    in perturbatively or discards it accordingly.
"""
struct  RasStep
    seFrom              ::Array{Shell,1}
    seTo                ::Array{Shell,1}
    deFrom              ::Array{Shell,1}
    deTo                ::Array{Shell,1}
    teFrom              ::Array{Shell,1}
    teTo                ::Array{Shell,1}
    qeFrom              ::Array{Shell,1}
    qeTo                ::Array{Shell,1}
    shapedTo            ::Array{Shell,1}
    frozenShells        ::Array{Shell,1}
    constraints         ::Array{String,1}
    treatment           ::Basics.AbstractQTreatment
    scfRoute            ::Union{Nothing,Basics.AbstractScfRoute}
end


"""
`AtomicState.RasStep()`  ... constructor for an 'empty' instance of a variable::AtomicState.RasStep
        A `rasStep::AtomicState.RasStep` is returned.
"""
function RasStep()
    RasStep(Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[],    Shell[], Shell[], String[],
            Basics.Variational(), nothing)
end


"""
`AtomicState.RasStep(rasStep::AtomicState.RasStep;
                     seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[],
                     deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[],
                     teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[],
                     qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[],
                     shapedTo::Array{Shell,1}=Shell[],
                     frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[],
                     treatment::Union{Nothing,Basics.AbstractQTreatment}=nothing,
                     scfRoute::Union{Nothing,Basics.AbstractScfRoute}=nothing)`

                    seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[],
                    deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[],
                    teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[],
                    qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[],
                    shapedTo::Array{Shell,1}=Shell[],
                    frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[],
                    treatment::Union{Nothing,Basics.AbstractQTreatment}=nothing,
                    scfRoute::Union{Nothing,Basics.AbstractScfRoute}=nothing

    ... constructor for modifying the given rasStep by specifying all excitations, frozen shells and
        constraints optionally.

        NOTE THAT `scfRoute` CANNOT BE CLEARED BACK TO `nothing` THROUGH THIS CONSTRUCTOR, since `nothing` is also how the
        keyword says "leave this field alone" -- which is the price of the keyword-copy pattern and is the same for every other
        field here.  Build a fresh `RasStep()` to drop an override.
        A `rasStep::AtomicState.RasStep` is returned.
"""
function RasStep(rasStep::AtomicState.RasStep;
                    seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[],
                    deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[],
                    teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[],
                    qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[],
                    shapedTo::Array{Shell,1}=Shell[],
                    frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[],
                    treatment::Union{Nothing,Basics.AbstractQTreatment}=nothing,
                    scfRoute::Union{Nothing,Basics.AbstractScfRoute}=nothing)
    if  seFrom == Shell[]   sxFrom = rasStep.seFrom   else      sxFrom = seFrom      end
    if  seTo   == Shell[]   sxTo   = rasStep.seTo     else      sxTo   = seTo        end
    if  deFrom == Shell[]   dxFrom = rasStep.deFrom   else      dxFrom = deFrom      end
    if  deTo   == Shell[]   dxTo   = rasStep.deTo     else      dxTo   = deTo        end
    if  teFrom == Shell[]   txFrom = rasStep.teFrom   else      txFrom = teFrom      end
    if  teTo   == Shell[]   txTo   = rasStep.teTo     else      txTo   = teTo        end
    if  qeFrom == Shell[]   qxFrom = rasStep.qeFrom   else      qxFrom = qeFrom      end
    if  qeTo   == Shell[]   qxTo   = rasStep.qeTo     else      qxTo   = qeTo        end
    if  shapedTo == Shell[]      shapedx= rasStep.shapedTo      else      shapedx= shapedTo           end
    if  frozen == Shell[]        frozx  = rasStep.frozenShells   else      frozx  = frozen             end
    if  constraints ==String[]   consx  = rasStep.constraints    else      consx  = constraints        end
    if  isnothing(treatment)     treatx = rasStep.treatment      else      treatx = treatment          end
    if  isnothing(scfRoute)      routex = rasStep.scfRoute       else      routex = scfRoute           end

    RasStep( sxFrom, sxTo, dxFrom, dxTo, txFrom, txTo, qxFrom, qxTo, shapedx, frozx, consx, treatx, routex)
end


function Base.string(step::AtomicState.RasStep)
    sa = "\nCI or RAS step with $(length(step.frozenShells)) (explicitly) frozen shell(s): $(step.frozenShells)  ... and virtual excitations"
    return( sa )
end


function Base.show(io::IO, step::AtomicState.RasStep)
    sa = Base.string(step);                 print(io, sa, "\n")
    if  length(step.seFrom) > 0
        sa = "   Singles from:          { ";      for  sh in step.seFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.seTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.deFrom) > 0
        sa = "   Doubles from:          { ";      for  sh in step.deFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.deTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.teFrom) > 0
        sa = "   Triples from:          { ";      for  sh in step.teFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.teTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.qeFrom) > 0
        sa = "   Quadruples from:       { ";      for  sh in step.qeFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.qeTo     sa = sa * string(sh) * ", "  end;
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end
    if  length(step.shapedTo) > 0
        sa = "   Shaped to CoreValenceCorrelationShape: { ";      for  sh in step.shapedTo   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }";             print(io, sa, "\n")
    end
    # A variational step is the default and says nothing; a perturbative one must announce itself, since the
    # CSF count it reports is then no longer the size of the space it accounts for.
    if  typeof(step.treatment) != Basics.Variational
        sa = "   Q-space treatment:     " * string(step.treatment);    print(io, sa, "\n")
    end
end


"""
`struct  AtomicState.RasLayer`
    ... a struct to specify a single layer of a multi-layer restricted active-space computation in a simplified,
        non-cumulative form: only the shells that are NEW at this layer need to be given; the frozen and active
        shell sets of the equivalent AtomicState.RasStep are derived automatically layer by layer, see
        AtomicState.RasExpansion(..., layers::Array{RasLayer,1}, ...).

    + newShells        ::Array{Shell,1}        ... Shells introduced (as virtual/correlation shells) at this layer.
    + se               ::Bool                  ... Include single excitations from the reference into all active
                                                    shells (up to and including this layer) if true.
    + de               ::Bool                  ... Include double excitations from the reference into all active
                                                    shells (up to and including this layer) if true.
    + treatment        ::Basics.AbstractQTreatment  ... how this layer's configurations are treated; see
                                                    AtomicState.RasStep. Variational() by default, so a layered
                                                    expansion written before 11-Sep-2026 is unaffected.
"""
struct  RasLayer
    newShells          ::Array{Shell,1}
    se                 ::Bool
    de                 ::Bool
    treatment          ::Basics.AbstractQTreatment
end


"""
`AtomicState.RasLayer(newShells::Array{Shell,1}; se::Bool=true, de::Bool=true,`
                     `treatment::Basics.AbstractQTreatment=Basics.Variational())`
    ... constructor for a variable::AtomicState.RasLayer with the given new shells and, by default, both single
        and double excitations enabled.
        A `rasLayer::AtomicState.RasLayer` is returned.
"""
function RasLayer(newShells::Array{Shell,1}; se::Bool=true, de::Bool=true,
                  treatment::Basics.AbstractQTreatment=Basics.Variational())
    RasLayer(newShells, se, de, treatment)
end


function Base.string(layer::AtomicState.RasLayer)
    sa = "RAS layer, adding shells $(layer.newShells)  (single excitations = $(layer.se), double excitations = $(layer.de))"
    return( sa )
end


function Base.show(io::IO, layer::AtomicState.RasLayer)
    sa = Base.string(layer);       print(io, sa, "\n")
end


"""
`struct  AtomicState.RasExpansion    <:  AbstractRepresentationType`
    ... a struct to represent (and generate) a restricted active-space representation.

    + symmetries       ::Array{LevelSymmetry,1}         ... Symmetries of the levels/CSF in the many-electron basis.
    + NoElectrons      ::Int64                          ... Number of electrons.
    + steps            ::Array{AtomicState.RasStep,1}   ... List of SCF steps that are to be done in this model 
                                                            computation.
    + settings         ::AtomicState.RasSettings        ... Settings for the given RAS computation
"""
struct         RasExpansion    <:  AbstractRepresentationType
    symmetries         ::Array{LevelSymmetry,1}
    NoElectrons        ::Int64
    steps              ::Array{AtomicState.RasStep,1}
    settings           ::AtomicState.RasSettings 
    end


"""
`AtomicState.RasExpansion()`  ... constructor for an 'empty' instance of the a variable::AtomicState.RasExpansion
        A `rasExpansion::AtomicState.RasExpansion` is returned.
"""
function RasExpansion()
    RasExpansion([Basics.LevelSymmetry(0, Basics.plus)], 0, AtomicState.RasStep[], AtomicState.RasSettings())
end


"""
`AtomicState.RasExpansion(symmetries::Array{LevelSymmetry,1}, NoElectrons::Int64, coreShells::Array{Shell,1},
                          fromShells::Array{Shell,1}, layers::Array{AtomicState.RasLayer,1}, settings::AtomicState.RasSettings)`
    ... constructor for a variable::AtomicState.RasExpansion from a simplified, non-cumulative list of
        AtomicState.RasLayer's; internally translated into the equivalent Array{AtomicState.RasStep,1}. For each
        layer, only the newly-introduced shells need to be given: `coreShells` (always frozen, never excited from)
        and `fromShells` (the reference valence shells that single/double excitations originate from) stay the
        same for every layer; `seTo`/`deTo` accumulate `fromShells` plus every newShells seen so far (including
        the current layer). `fromShells` are optimized together with the FIRST layer's new shells (their one and
        only chance to vary, exactly as for a plain reference SCF), then frozen from the second layer onward,
        alongside every earlier layer's own newShells (never the current layer's).
        A `rasExpansion::AtomicState.RasExpansion` is returned.
"""
function RasExpansion(symmetries::Array{LevelSymmetry,1}, NoElectrons::Int64, coreShells::Array{Shell,1},
                       fromShells::Array{Shell,1}, layers::Array{AtomicState.RasLayer,1}, settings::AtomicState.RasSettings)
    steps  = AtomicState.RasStep[]
    frozen = deepcopy(coreShells)
    to     = deepcopy(fromShells)
    prior  = AtomicState.RasStep()
    for  (i, layer)  in  enumerate(layers)
        append!(to, layer.newShells)
        sxFrom = layer.se ? deepcopy(fromShells) : Shell[];    sxTo = layer.se ? deepcopy(to) : Shell[]
        dxFrom = layer.de ? deepcopy(fromShells) : Shell[];    dxTo = layer.de ? deepcopy(to) : Shell[]
        step   = AtomicState.RasStep(prior; seFrom=sxFrom, seTo=sxTo, deFrom=dxFrom, deTo=dxTo, frozen=deepcopy(frozen),
                                             treatment=layer.treatment)
        push!(steps, step)
        if  i == 1   append!(frozen, fromShells)   end
        append!(frozen, layer.newShells)
        prior  = step
    end

    RasExpansion(symmetries, NoElectrons, steps, settings)
end


function Base.string(expansion::RasExpansion)
    sa = "RAS expansion for symmétry $(expansion.symmetries) and with $(length(expansion.steps)) steps:"
    return( sa )
end


function Base.show(io::IO, expansion::RasExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "$(expansion.steps)  ")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end



### RAS automation: layer verdicts, and the policies that decide an automatic ladder ##########################################
"""
`struct  AtomicState.RasLayerVerdict`
    ... what one completed RAS layer is worth, in numbers a stopping or shell policy can read, modelled on
        SelfConsistent.ScfVerdict: advice about what has already happened, never a gate on it. It is built by
        the driver AFTER a layer has run and is never consulted while that layer is still being computed.

    + istep         ::Int64          ... which layer this is; 1 is the reference layer, which carries no increment.
    + newShells     ::Array{Shell,1} ... the shells this layer introduced, for the printed record.
    + csfCount      ::Int64          ... size of this layer's CI space -- the cheap, reliable cost signal
        (reference_ras_cost_scaling found no clean per-CSF law, so this stays a raw count rather than a derived cost).
    + wallTime      ::Float64        ... seconds this layer took to run.
    + energy        ::Float64        ... the lowest target-level energy at this layer.
    + increment     ::Float64        ... `energy` minus the PREVIOUS layer's energy; `NaN` for the reference layer,
        which has no previous layer to compare against. This is Rule 21's own quantity -- the one a policy must
        look at, never the total.
    + scfVerdict    ::Union{Nothing,SelfConsistent.ScfVerdict}   ... this layer's own SCF verdict, so a policy that
        wants the doubling-test drift, not only the increment, can read it; `nothing` for a layer treated to second
        order, since that layer's orbitals are never optimized and no SCF verdict exists for it.
    + shapedBasis   ::Union{Nothing,Basis}   ... the refined, separate Q-space `Basis` a `RasStep.shapedTo`
        layer built -- its own CSF space (the reference plus that layer's shaped excitation alone), with every
        OTHER subshell's orbital the literal SAME object `priorMultiplet` already used, and the new shell(s)
        tagged `Basics.CoreValenceCorrelationShape` and SCF-refined on this small space. `nothing` for every
        ordinary layer. `energy` and `increment` above are NOT about this space -- they stay the ladder's own
        P-progression, unchanged by a shaped layer, which is why this field exists separately rather than
        folding Q's energy into them.
    + shapedPartition ::Union{Nothing,Hamiltonian.ShapeCrossingPartition}   ... the second-order P/Q coupling
        this shape-crossing layer is actually FOR -- `Hamiltonian.secondOrderShapeCrossingCorrection` applied
        between `priorMultiplet` (P) and `shapedBasis`'s own multiplet (Q), biorthogonalized first. `nothing`
        for every ordinary layer, and whenever `shapedBasis` is `nothing`.
"""
struct  RasLayerVerdict
    istep           ::Int64
    newShells       ::Array{Shell,1}
    csfCount        ::Int64
    wallTime        ::Float64
    energy          ::Float64
    increment       ::Float64
    scfVerdict      ::Union{Nothing,SelfConsistent.ScfVerdict}
    shapedBasis     ::Union{Nothing,Basis}
    shapedPartition ::Union{Nothing,Hamiltonian.ShapeCrossingPartition}
end


# `Base.show(io::IO, verdict::RasLayerVerdict)`  ... prints one line summarizing the layer verdict.
function Base.show(io::IO, verdict::RasLayerVerdict)
    sa = "RAS layer $(verdict.istep), adding $(verdict.newShells): $(verdict.csfCount) CSF, " *
         @sprintf("%.1f s", verdict.wallTime) * ", energy = $(verdict.energy) Ha, increment = " *
         (isnan(verdict.increment) ? "n/a (reference layer)" : @sprintf("%+.4f mHa", 1000*verdict.increment)) *
         (isnothing(verdict.shapedBasis) ? "" : "  (shape-crossing leaf: the ladder's own P-progression is untouched)") *
         (isnothing(verdict.shapedPartition) ? "" :
          @sprintf("; P/Q second order %+.4f mHa over %d target level(s)", 1000*sum(verdict.shapedPartition.deltaE),
                   length(verdict.shapedPartition.targetLevels)))
    println(io, sa)
end


"""
`struct  AtomicState.RasStoppingVerdict`
    ... what a stopping policy concluded about whether a RAS ladder should propose another layer, in the same
        advisory spirit as SelfConsistent.ScfVerdict: it is read by the driver and PRINTED, never used to
        silently decide anything a human cannot see.

    + stop      ::Bool      ... whether the policy says to stop proposing further layers.
    + reason    ::String    ... one sentence naming WHY, printed verbatim; empty when `stop` is false.
"""
struct  RasStoppingVerdict
    stop        ::Bool
    reason      ::String
end


# `Base.show(io::IO, verdict::RasStoppingVerdict)`  ... prints the policy's verdict.
function Base.show(io::IO, verdict::RasStoppingVerdict)
    println(io, verdict.stop ? "STOP: $(verdict.reason)" : "continue")
end


"""
`AtomicState.collectLadderShells(coreShells::Array{Shell,1}, fromShells::Array{Shell,1}, layersSoFar::Array{RasLayer,1})`  
    ... gathers every shell that the RAS ladder has already used -- the core, the reference valence, and the new
        shells of each layer built so far -- so that a shell policy can exclude them from what it proposes next.  A
        shells::Array{Shell,1} is returned, with duplicates removed and the order of first appearance kept.  This is
        plumbing shared by every concrete shell policy rather than something a caller builds or inspects, and it is
        therefore not exported;  a shell present anywhere in the ladder is never proposed a second time by any policy.
"""
function collectLadderShells(coreShells::Array{Shell,1}, fromShells::Array{Shell,1}, layersSoFar::Array{RasLayer,1})
    used = vcat(deepcopy(coreShells), deepcopy(fromShells))
    for layer in layersSoFar    append!(used, layer.newShells)    end

    return( unique(used) )
end


"""
`abstract type AtomicState.AbstractRasStoppingPolicy`
    ... defines an abstract and a number of concrete types to decide whether another RAS layer should be
        PROPOSED, given the record of every layer run so far. A policy is read-only advice: it answers one
        question (stop, or go on) and says WHY, and it cannot abort a layer that is already running. The
        hierarchy is intentionally small and is meant to grow -- a new stopping rule is a new struct plus one
        `AtomicState.evaluateStoppingPolicy` method, never a change to the driver that calls it, and never a
        change to any existing policy.
"""
abstract type  AbstractRasStoppingPolicy    end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.AbstractRasStoppingPolicy, history::Array{AtomicState.RasLayerVerdict,1})`
    ... asks `policy` whether the ladder recorded in `history` (oldest first, one entry per layer already run)
        should stop here. A `verdict::AtomicState.RasStoppingVerdict` is returned; see the concrete policies
        for what each one actually tests.
"""
function evaluateStoppingPolicy end


"""
`struct  AtomicState.FixedLayerCount  <:  AtomicState.AbstractRasStoppingPolicy`
    ... stops once exactly `n` layers have run. This is the escape hatch that makes the automatic driver a
        strict superset of today's hand-built one: wrapping a fully-specified, hand-written layer list in
        `FixedLayerCount(length(layers))` reproduces exactly the existing behaviour.

    + n     ::Int64   ... the number of layers to run, reference layer included.
"""
struct  FixedLayerCount  <:  AtomicState.AbstractRasStoppingPolicy
    n       ::Int64
end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.FixedLayerCount, history::Array{AtomicState.RasLayerVerdict,1})`
    ... stops once `length(history) >= policy.n`. A `verdict::AtomicState.RasStoppingVerdict` is returned.
"""
function evaluateStoppingPolicy(policy::FixedLayerCount, history::Array{RasLayerVerdict,1})
    if  length(history) >= policy.n
        return( RasStoppingVerdict(true, "FixedLayerCount($(policy.n)): $(length(history)) layer(s) have run.") )
    end

    return( RasStoppingVerdict(false, "") )
end


"""
`struct  AtomicState.IncrementBelow  <:  AtomicState.AbstractRasStoppingPolicy`
    ... stops once `consecutive` layer increments IN A ROW fall below `tolerance`. A single small increment
        proves nothing on its own -- `reference_ras_practical_recipe` measured a case where ONE frozen-looking
        increment understated the true remaining error by a factor of 100 -- so this policy refuses to trust
        one data point; `consecutive` is 2 by default and should not be lowered to 1 without reading why.

    + tolerance     ::Float64   ... the absolute increment (in Hartree) a layer must fall below.
    + consecutive   ::Int64     ... how many layers in a row must fall below `tolerance` before stopping.
"""
struct  IncrementBelow  <:  AtomicState.AbstractRasStoppingPolicy
    tolerance       ::Float64
    consecutive     ::Int64
end


"""
`AtomicState.IncrementBelow(tolerance::Float64; consecutive::Int64=2)`
    ... constructor taking the tolerance positionally and `consecutive` as a keyword, defaulted to the safe
        value of 2. An `policy::AtomicState.IncrementBelow` is returned.
"""
function IncrementBelow(tolerance::Float64; consecutive::Int64=2)
    IncrementBelow(tolerance, consecutive)
end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.IncrementBelow, history::Array{AtomicState.RasLayerVerdict,1})`
    ... stops once the last `policy.consecutive` layers all have `abs(increment) < policy.tolerance`. A
        `verdict::AtomicState.RasStoppingVerdict` is returned.
"""
function evaluateStoppingPolicy(policy::IncrementBelow, history::Array{RasLayerVerdict,1})
    if  length(history) < policy.consecutive + 1    # +1: the reference layer itself carries no increment
        return( RasStoppingVerdict(false, "") )
    end
    tail = history[end-policy.consecutive+1:end]
    if  all( v -> !isnan(v.increment)  &&  abs(v.increment) < policy.tolerance,  tail )
        return( RasStoppingVerdict(true, "IncrementBelow($(policy.tolerance)): the last $(policy.consecutive) " *
                "layer increment(s) were all below tolerance -- a HINT, not proof the ladder has converged; " *
                "see reference_ras_practical_recipe for a case where this understated the remaining error by " *
                "up to 100x.") )
    end

    return( RasStoppingVerdict(false, "") )
end


"""
`struct  AtomicState.IncrementRatioBelow  <:  AtomicState.AbstractRasStoppingPolicy`
    ... stops once the newest layer increment's magnitude falls below `fraction` times the PREVIOUS layer's,
        i.e. a decaying-ladder test. Useful across different elements/ions without re-tuning an absolute
        tolerance per case, since Rule 21's own warning is that the increment is not an error bar in absolute
        terms -- but its rate of decay is at least comparable between elements without rescaling.

    + fraction  ::Float64   ... the ratio |increment_n| / |increment_(n-1)| must fall below.
"""
struct  IncrementRatioBelow  <:  AtomicState.AbstractRasStoppingPolicy
    fraction    ::Float64
end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.IncrementRatioBelow, history::Array{AtomicState.RasLayerVerdict,1})`
    ... stops once `abs(history[end].increment) / abs(history[end-1].increment) < policy.fraction`. A
        `verdict::AtomicState.RasStoppingVerdict` is returned.
"""
function evaluateStoppingPolicy(policy::IncrementRatioBelow, history::Array{RasLayerVerdict,1})
    if  length(history) < 3    return( RasStoppingVerdict(false, "") )    end   # need two non-NaN increments
    last = history[end];    prev = history[end-1]
    if  isnan(last.increment)  ||  isnan(prev.increment)  ||  prev.increment == 0.
        return( RasStoppingVerdict(false, "") )
    end
    ratio = abs(last.increment) / abs(prev.increment)
    if  ratio < policy.fraction
        return( RasStoppingVerdict(true, "IncrementRatioBelow($(policy.fraction)): the latest increment is " *
                @sprintf("%.1f %%", 100*ratio) * " of the previous one.") )
    end

    return( RasStoppingVerdict(false, "") )
end


"""
`struct  AtomicState.CostCeiling  <:  AtomicState.AbstractRasStoppingPolicy`
    ... refuses to propose a further layer once the LAST layer that actually ran already met or exceeded the
        given cost, in CSF count or wall time. This is a safety valve, not a convergence judgement, and it is
        deliberately phrased against the last MEASURED layer rather than a predicted next one, since the
        existing `tryRun` cost estimator is known to be off by 3-4x (`reference_ras_practical_recipe`) and no
        clean per-CSF cost law fits any case measured so far (`reference_ras_cost_scaling`).

    + maxCsf        ::Int64     ... the CSF count a layer may not meet or exceed.
    + maxWallTime   ::Float64   ... the wall time (seconds) a layer may not meet or exceed.
"""
struct  CostCeiling  <:  AtomicState.AbstractRasStoppingPolicy
    maxCsf          ::Int64
    maxWallTime     ::Float64
end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.CostCeiling, history::Array{AtomicState.RasLayerVerdict,1})`
    ... stops once the last layer's own `csfCount` or `wallTime` meets or exceeds the ceiling. A
        `verdict::AtomicState.RasStoppingVerdict` is returned.
"""
function evaluateStoppingPolicy(policy::CostCeiling, history::Array{RasLayerVerdict,1})
    if  isempty(history)    return( RasStoppingVerdict(false, "") )    end
    last = history[end]
    if  last.csfCount >= policy.maxCsf  ||  last.wallTime >= policy.maxWallTime
        return( RasStoppingVerdict(true, "CostCeiling: the last layer already used $(last.csfCount) CSF in " *
                @sprintf("%.1f s", last.wallTime) * ", at or beyond the ceiling ($(policy.maxCsf) CSF, " *
                @sprintf("%.1f s", policy.maxWallTime) * ").") )
    end

    return( RasStoppingVerdict(false, "") )
end


"""
`struct  AtomicState.CombinedPolicy  <:  AtomicState.AbstractRasStoppingPolicy`
    ... stops as soon as ANY of `policies` says to stop (a logical OR), and is the actual seam that keeps the
        hierarchy expandable: a brand-new policy type needs no change here, it simply becomes one more entry
        in somebody's list. A convergence rule alone, trusted without a cost ceiling and a hard layer count
        behind it, is exactly the silent over-trust this item must not build -- so a `CombinedPolicy` such as
        `CombinedPolicy([IncrementBelow(1.0e-4), CostCeiling(5000, 1800.), FixedLayerCount(6)])` is the
        recommended way to use this hierarchy, rather than any one policy alone.

    + policies  ::Array{AtomicState.AbstractRasStoppingPolicy,1}   ... the sub-policies, tested in order.
"""
struct  CombinedPolicy  <:  AtomicState.AbstractRasStoppingPolicy
    policies    ::Array{AtomicState.AbstractRasStoppingPolicy,1}
end


"""
`AtomicState.evaluateStoppingPolicy(policy::AtomicState.CombinedPolicy, history::Array{AtomicState.RasLayerVerdict,1})`
    ... tests every sub-policy in order and stops on the first one that does, carrying its reason forward. A
        `verdict::AtomicState.RasStoppingVerdict` is returned.
"""
function evaluateStoppingPolicy(policy::CombinedPolicy, history::Array{RasLayerVerdict,1})
    for  sub in policy.policies
        verdict = AtomicState.evaluateStoppingPolicy(sub, history)
        if  verdict.stop    return( verdict )    end
    end

    return( RasStoppingVerdict(false, "") )
end


"""
`abstract type AtomicState.AbstractRasShellPolicy`
    ... defines an abstract and a number of concrete types to decide WHICH new shells the next RAS layer
        introduces, given what the expansion already has (coreShells, fromShells, every earlier layer's own
        newShells). AtomicState.AbstractRasStoppingPolicy decides WHETHER another layer runs; this family
        decides WHAT it would contain if it did -- a separate question, asked first. Like the stopping
        hierarchy, this one is deliberately small and meant to grow by adding a struct and one
        `AtomicState.proposeNextShells` method, never by changing the driver or an existing policy.
"""
abstract type  AbstractRasShellPolicy    end


"""
`AtomicState.proposeNextShells(policy::AtomicState.AbstractRasShellPolicy, coreShells::Array{Shell,1},`
                               `fromShells::Array{Shell,1}, layersSoFar::Array{AtomicState.RasLayer,1},`
                               `candidateOrbitals::Dict{Subshell,Orbital})`
    ... asks `policy` which shells the NEXT layer should introduce; `candidateOrbitals` is the empty dict for a
        policy that does not need one (see AtomicState.needsEnergySpectrum). An `Array{Shell,1}` is returned,
        which may be empty if the policy has nothing left to propose -- an empty proposal is read by the
        driver as "stop", independent of whatever the stopping policy says, since there is nothing left to run.
"""
function proposeNextShells end


"""
`AtomicState.needsEnergySpectrum(policy::AtomicState.AbstractRasShellPolicy)`
    ... whether `policy` needs a precomputed, named one-electron spectrum to propose shells. Defaults to
        `false`, so a new, purely combinatorial policy needs no change here; a new policy that DOES read
        energies overrides this method, the same one-method-to-add discipline as the rest of this hierarchy.
        A `value::Bool` is returned.
"""
needsEnergySpectrum(policy::AbstractRasShellPolicy) = false


"""
`struct  AtomicState.RasCompleteShellPolicy  <:  AtomicState.AbstractRasShellPolicy`
    ... proposes every `l = 0, ..., n-1` of the next unused principal quantum number `n` -- GRASP's own
        default layering convention, and the one every RAS example in this codebase has used by hand so far.
"""
struct  RasCompleteShellPolicy  <:  AtomicState.AbstractRasShellPolicy    end


"""
`AtomicState.proposeNextShells(policy::AtomicState.RasCompleteShellPolicy, coreShells::Array{Shell,1},`
                               `fromShells::Array{Shell,1}, layersSoFar::Array{AtomicState.RasLayer,1},`
                               `candidateOrbitals::Dict{Subshell,Orbital})`
    ... proposes the complete next principal shell, `n = 1 + maximum(n of every shell already in the ladder)`.
        `candidateOrbitals` is ignored. An `Array{Shell,1}` is returned.
"""
function proposeNextShells(policy::RasCompleteShellPolicy, coreShells::Array{Shell,1}, fromShells::Array{Shell,1},
                            layersSoFar::Array{RasLayer,1}, candidateOrbitals::Dict{Subshell,Orbital})
    used  = AtomicState.collectLadderShells(coreShells, fromShells, layersSoFar)
    nextN = isempty(used) ? 1 : maximum(sh.n for sh in used) + 1

    return( Shell[ Shell(nextN, l)  for l in 0:nextN-1 ] )
end


"""
`struct  AtomicState.RasLimitedAngularMomentumPolicy  <:  AtomicState.AbstractRasShellPolicy`
    ... caps how far the active space reaches in `l`, INDEPENDENTLY of how far `n` climbs -- the remedy for
        AtomicState.RasCompleteShellPolicy being "rather extreme" once `n` exceeds about 5 (at `n=6` it would
        propose `6s,6p,6d,6f,6g,6h` in one layer, mostly high-`l` partial waves far from anything the
        reference occupies, and `reference_ras_practical_recipe` measured that a layer of high-`n`,
        already-occupied-symmetry orbitals is both the hardest to converge and the least rewarding).

        For the next unused `n`, this policy proposes only `l = 0, ..., min(n-1, lMaxRef + extraL)`, where
        `lMaxRef` is the largest `l` occupied by `coreShells` or `fromShells`. So with `extraL = 2` on a
        `3d^2` reference (`lMaxRef = 2`), layer `n=6` proposes only `6s,6p,6d,6f` -- `6g` and `6h` are never
        opened, however high `n` goes, because they sit beyond `lMaxRef+extraL = 4`. The ladder keeps
        climbing in `n` without ever again reaching for a new, likely-unrewarding high-`l` partial wave.

    + extraL    ::Int64   ... how many `l`-values beyond the reference's own maximum may ever be opened.
"""
struct  RasLimitedAngularMomentumPolicy  <:  AtomicState.AbstractRasShellPolicy
    extraL      ::Int64
end


"""
`AtomicState.proposeNextShells(policy::AtomicState.RasLimitedAngularMomentumPolicy, coreShells::Array{Shell,1},`
                               `fromShells::Array{Shell,1}, layersSoFar::Array{AtomicState.RasLayer,1},`
                               `candidateOrbitals::Dict{Subshell,Orbital})`
    ... proposes the next unused `n`'s shells with `l` capped at `lMaxRef + policy.extraL`, as described in
        the struct's own docstring. `candidateOrbitals` is ignored. An `Array{Shell,1}` is returned.
"""
function proposeNextShells(policy::RasLimitedAngularMomentumPolicy, coreShells::Array{Shell,1}, fromShells::Array{Shell,1},
                            layersSoFar::Array{RasLayer,1}, candidateOrbitals::Dict{Subshell,Orbital})
    used    = AtomicState.collectLadderShells(coreShells, fromShells, layersSoFar)
    nextN   = isempty(used) ? 1 : maximum(sh.n for sh in used) + 1
    lMaxRef = maximum( sh.l  for sh in vcat(coreShells, fromShells) )
    lCap    = min(nextN - 1, lMaxRef + policy.extraL)

    return( Shell[ Shell(nextN, l)  for l in 0:lCap ] )
end


"""
`struct  AtomicState.RasSpectrumOrderedShellPolicy  <:  AtomicState.AbstractRasShellPolicy`
    ... for medium and heavy elements, where neither AtomicState.RasCompleteShellPolicy nor
        AtomicState.RasLimitedAngularMomentumPolicy is reliable, because both assume shells group cleanly by
        `n` -- which stops holding once a real atomic potential, rather than a bare hydrogenic one, is
        involved. For a `3d^3 4s^2`-type reference, the next worthwhile correlation orbitals are not "finish
        n=4, then start n=5" but something that may read `5s, 4p, 4d, 4f, ...`: shells from DIFFERENT
        principal quantum numbers interleaved by where they actually sit, because screening and `3d`/`4s`
        mixing move `4p`, `4d`, `4f` and `5s` out of their light-atom, hydrogenic order.

        This policy answers the question numerically rather than combinatorially, reusing
        AtomicState.OneElectronSpectrum -- generated once by the driver into `candidateOrbitals`, see
        AtomicState.needsEnergySpectrum -- rather than inventing new orbital generation: for every partial
        wave `l` present in `candidateOrbitals`, it finds the lowest-`n` shell of that `l` not already in the
        ladder, ranks those one-candidate-per-`l` shells by their ACTUAL one-electron energy, and proposes
        the `count` cheapest.

        THIS IS THE ONE PLACE THIS HIERARCHY DOES SOMETHING GRASP'S OWN WORKFLOW NEVER HAS TO, AND IT IS
        WORTH KNOWING WHY. GRASP's correlation virtuals are generated per `kappa` and referred to only by
        their sequential index within that `kappa` -- "the 3rd d_5/2 virtual" -- never by a principal
        quantum number, because once an orbital is correlated and polarized, the node-counting convention
        that even DEFINES `n` for a hydrogenic function stops applying cleanly, so GRASP never has to decide
        whether that 3rd virtual "is" 4d or 5d. JAC's RAS API cannot take that way out: `RasLayer.newShells`,
        every `se/deFrom/To` list, `frozenShells` and this very policy's own candidates are all referred to
        by name, because that is how a user reads and writes a RAS ladder. So where GRASP can simply ask its
        spectrum for the next cheapest orbital of a given `kappa` and never decide what to call it, this
        policy keeps names attached to a measured energy ordering that may cross `n`.

        A diffuse, box-truncated state can rank artificially low in energy (Rule 12's territory), so this
        policy inherits a risk the other two never face, since they never consult a computed energy at all --
        its first real use on any one case is a judgement call for a human, not an automatic pass/fail.

    + count             ::Int64                            ... how many new shells to propose per layer.
    + spectrumSettings  ::AtomicState.OneElectronSettings   ... the candidate pool (nMax, lValues) the
        driver generates, once, via the existing `Basics.generate(::AtomicState.OneElectronSpectrum, rep)`.
"""
struct  RasSpectrumOrderedShellPolicy  <:  AtomicState.AbstractRasShellPolicy
    count               ::Int64
    spectrumSettings    ::AtomicState.OneElectronSettings
end


"""
`AtomicState.needsEnergySpectrum(policy::AtomicState.RasSpectrumOrderedShellPolicy)`
    ... always `true` -- this policy cannot propose anything without the generated candidate spectrum. A
        `value::Bool` is returned.
"""
needsEnergySpectrum(policy::RasSpectrumOrderedShellPolicy) = true


"""
`AtomicState.proposeNextShells(policy::AtomicState.RasSpectrumOrderedShellPolicy, coreShells::Array{Shell,1},`
                               `fromShells::Array{Shell,1}, layersSoFar::Array{AtomicState.RasLayer,1},`
                               `candidateOrbitals::Dict{Subshell,Orbital})`
    ... ranks the lowest-`n` unused shell of every partial wave present in `candidateOrbitals` by its actual
        one-electron energy there, and proposes the `policy.count` cheapest, as described in the struct's own
        docstring. An `Array{Shell,1}` is returned, shorter than `policy.count` once the candidate pool itself
        is exhausted, and empty once nothing more can be proposed from it.
"""
function proposeNextShells(policy::RasSpectrumOrderedShellPolicy, coreShells::Array{Shell,1}, fromShells::Array{Shell,1},
                            layersSoFar::Array{RasLayer,1}, candidateOrbitals::Dict{Subshell,Orbital})
    used      = AtomicState.collectLadderShells(coreShells, fromShells, layersSoFar)
    lValues   = unique( Basics.subshell_l(ss)  for ss in keys(candidateOrbitals) )
    candidates = Tuple{Shell,Float64}[]
    for  l in lValues
        nsOfL = sort(unique( ss.n  for ss in keys(candidateOrbitals)  if Basics.subshell_l(ss) == l ))
        for  n in nsOfL
            sh = Shell(n, l)
            if  !(sh in used)
                energy = minimum( orb.energy  for (ss, orb) in candidateOrbitals
                                              if  ss.n == n  &&  Basics.subshell_l(ss) == l )
                push!(candidates, (sh, energy))
                break    # the lowest unused n of this l only; a higher n of the SAME l is not offered yet
            end
        end
    end
    sort!(candidates, by = c -> c[2])

    return( Shell[ c[1]  for c in candidates[1:min(policy.count, length(candidates))] ] )
end


"""
`struct  AtomicState.RasAutoSettings`
    ... bundles the choices an automatic RAS ladder needs: WHICH shells to propose next, WHETHER to propose
        another layer at all, and a hard ceiling that applies no matter what either policy says.

    + stoppingPolicy    ::AtomicState.AbstractRasStoppingPolicy   ... decides whether another layer runs.
    + shellPolicy       ::AtomicState.AbstractRasShellPolicy      ... decides what it would contain if it did.
    + layerTemplate     ::AtomicState.RasLayer                    ... the `se`/`de`/`treatment` to reuse for
        every auto-proposed layer; its own `newShells` field is ignored, since the shell policy fills that in.
    + maxLayers         ::Int64                                   ... a hard ceiling that ALWAYS applies,
        independent of `stoppingPolicy`, so an automatic ladder can never run unbounded even if a
        AtomicState.CombinedPolicy forgets a count-based term.
"""
struct  RasAutoSettings
    stoppingPolicy      ::AtomicState.AbstractRasStoppingPolicy
    shellPolicy         ::AtomicState.AbstractRasShellPolicy
    layerTemplate       ::AtomicState.RasLayer
    maxLayers           ::Int64
end


"""
`AtomicState.RasAutoSettings(stoppingPolicy::AtomicState.AbstractRasStoppingPolicy,`
                             `shellPolicy::AtomicState.AbstractRasShellPolicy;`
                             `layerTemplate::AtomicState.RasLayer=AtomicState.RasLayer(Shell[]), maxLayers::Int64=10)`
    ... constructor taking the two policies positionally and the template/ceiling as keywords with safe
        defaults. A `settings::AtomicState.RasAutoSettings` is returned.
"""
function RasAutoSettings(stoppingPolicy::AbstractRasStoppingPolicy, shellPolicy::AbstractRasShellPolicy;
                          layerTemplate::RasLayer=RasLayer(Shell[]), maxLayers::Int64=10)
    RasAutoSettings(stoppingPolicy, shellPolicy, layerTemplate, maxLayers)
end



### CI: Configuration Interaction expansions ##################################################################################
"""
`struct  AtomicState.CiSettings`  
    ... a struct for defining the settings for a configuration-interaction (CI) expansion.

    + eeInteractionCI      ::AbstractEeInteraction   ... Specifies the treatment of the e-e interaction.
    + levelSelectionCI     ::LevelSelection          ... Specifies the selected levels, if any.
"""
struct  CiSettings
    eeInteractionCI        ::AbstractEeInteraction 
    levelSelectionCI       ::LevelSelection
end


"""
`AtomicState.CiSettings()`  ... constructor for setting the default values.
        A `ciSettings::AtomicState.CiSettings` is returned.
"""
function CiSettings()
    CiSettings(CoulombInteraction(), LevelSelection() )
end


"""
`AtomicState.CiSettings(settings::AtomicState.CiSettings;
    eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,        levelSelectionCI::Union{Nothing,LevelSelection}=nothing)`

    ... constructor for modifying the given CiSettings by 'overwriting' the explicitly selected parameters.
        A `ciSettings::AtomicState.CiSettings` is returned.
"""
function CiSettings(settings::AtomicState.CiSettings;
    eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,        levelSelectionCI::Union{Nothing,LevelSelection}=nothing)
    
    if  isnothing(eeInteractionCI)       eeInteractionCIx      = settings.eeInteractionCI       else   eeInteractionCIx      = eeInteractionCI      end 
    if  isnothing(levelSelectionCI)      levelSelectionCIx     = settings.levelSelectionCI      else   levelSelectionCIx     = levelSelectionCI     end 
    
    CiSettings( eeInteractionCIx, levelSelectionCIx)
end


function Base.show(io::IO, settings::CiSettings)
        println(io, "eeInteractionCI:          $(settings.eeInteractionCI)  ")
        println(io, "levelSelectionCI:         $(settings.levelSelectionCI)  ")
end


"""
`struct  AtomicState.CiExpansion  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a configuration-interaction representation.

    + applyOrbitals    ::Dict{Subshell, Orbital}
    + excitations      ::AtomicState.RasStep            ... Excitations beyond refConfigs.
    + settings         ::AtomicState.CiSettings         ... Settings for the given CI expansion
"""
struct         CiExpansion     <:  AbstractRepresentationType
    applyOrbitals      ::Dict{Subshell, Orbital}
    excitations        ::AtomicState.RasStep
    settings           ::AtomicState.CiSettings
end


"""
`AtomicState.CiExpansion()`  
    ... constructor for an 'empty' instance of a variable::AtomicState.CiExpansion, i.e. one with no orbitals, no excitations beyond
        the reference configurations and the default settings.  It is also what `AtomicState.Representation()` is built from, so the
        zero-argument form of a representation depends on this one existing.  A ciExpansion::AtomicState.CiExpansion is returned.
"""
function CiExpansion()
    CiExpansion( Dict{Subshell, Orbital}(), AtomicState.RasStep(), AtomicState.CiSettings() )
end


function Base.string(expansion::CiExpansion)
    sa = "CI expansion with (additional) excitations:"
    return( sa )
end


function Base.show(io::IO, expansion::CiExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "$(expansion.excitations)  ")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end


    
### Green function representation ##################################################################################

"""
`abstract type AtomicState.AbstractGreenApproach` 
    ... defines an abstract and a number of singleton types for approximating a many-electron Green
        function expansion for calculating second-order processes.

    + struct SingleCSFwithoutCI        
    ... to approximate the many-electron multiplets (gMultiplet) for every chosen level symmetry by dealing with each CSF 
        independently, and without any configuration interaction. This is a fast but also very rough approximation.
        
    + struct CoreSpaceCI                
    ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction between the 
        bound-state orbitals into account.
        
    + struct DampedSpaceCI                
    ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction for all, the 
        bound and free-electron, orbitals into account but by including a damping factor e^{tau*r}, tau > 0 into the 
        electron densities rho_ab (r) --> rho_ab (r) * e^{tau*r}
"""
abstract type  AbstractGreenApproach                                  end
struct         SingleCSFwithoutCI  <:  AtomicState.AbstractGreenApproach   end
struct         CoreSpaceCI         <:  AtomicState.AbstractGreenApproach   end
struct         DampedSpaceCI       <:  AtomicState.AbstractGreenApproach   end


"""
`struct  AtomicState.GreenSettings`  
    ... defines a type for defining the details and parameters of the approximate Green (function) expansion.

    + nMax                     ::Int64            ... maximum principal quantum numbers of (single-electron) 
                                                        excitations that are to be included into the representation.
    + lValues                  ::Array{Int64,1}   ... List of (non-relativistic) orbital angular momenta for which
                                                        (single-electron) excitations are to be included.
    + dampingTau               ::Float64          ... factor tau (> 0.) that is used to 'damp' the one- and two-electron
                                                        interactions strength: exp( - tau * r)
    + printBefore              ::Bool             ... True if a short overview is to be printed before. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
    + scField                  ::AbstractScField
        ... the self-consistent field in which the Green space is generated.

        WHY THIS FIELD EXISTS (added 07-Aug-2026). A Green expansion is used as the INTERMEDIATE spectrum of a
        second-order calculation, and gauge invariance there requires the intermediate states to be eigenstates
        of the SAME one-body Hamiltonian as the initial and final states. Until now `generate(GreenExpansion,...)`
        hardcoded `AsfSettings()` for the reference SCF and `Basics.DFSField(1.0)` for the mean potential, so a
        user had NO way to match it to their own computation -- and a mismatch breaks gauge invariance
        systematically without preventing either gauge from converging on its own, i.e. invisibly. Measured on
        H 2s -> 1s two-photon decay, moving from a mismatched to a matched potential took the length gauge from
        35 % error to 5.5 %.
"""
struct GreenSettings 
    nMax                       ::Int64
    lValues                    ::Array{Int64,1}
    dampingTau                 ::Float64
    printBefore                ::Bool 
    levelSelection             ::LevelSelection
    scField                    ::AbstractScField
end 


"""
`AtomicState.GreenSettings()`  ... constructor for an `empty` instance of AtomicState.GreenSettings.
        A `greenSettings::AtomicState.GreenSettings` is returned.
"""
function GreenSettings()
    # `GreenSettings`, not `Settings` (fixed 07-Aug-2026): this constructor called `Settings(...)`, so the
    # zero-argument form either raised or returned the wrong type. Every caller had to use the positional form.
    GreenSettings( 0, Int64[], 0., false, LevelSelection(), Basics.DFSField() )
end


"""
`AtomicState.GreenSettings(nMax::Int64, lValues::Array{Int64,1}, dampingTau::Float64, printBefore::Bool,
                        levelSelection::LevelSelection)`
    ... backward-compatible five-argument constructor; the self-consistent field defaults to Basics.DFSField(),
        which is what `generate(GreenExpansion, ...)` used unconditionally before `scField` existed. Supply the
        sixth argument (or use the keyword constructor) whenever the Green space must match the potential of a
        surrounding computation -- see the note in the struct docstring.
        A `greenSettings::AtomicState.GreenSettings` is returned.
"""
function GreenSettings(nMax::Int64, lValues::Array{Int64,1}, dampingTau::Float64, printBefore::Bool,
                       levelSelection::LevelSelection)
    GreenSettings( nMax, lValues, dampingTau, printBefore, levelSelection, Basics.DFSField() )
end


"""
`AtomicState.GreenSettings(settings::AtomicState.GreenSettings;
    nMax::Union{Nothing,Int64}=nothing,                     lValues::Union{Nothing,Array{Int64,1}}=nothing,
    dampingTau::Union{Nothing,Float64}=nothing,             printBefore::Union{Nothing,Bool}=nothing,
    levelSelection::Union{Nothing,LevelSelection}=nothing,  scField::Union{Nothing,AbstractScField}=nothing)`

    ... the standard JAC keyword copy-constructor.
        A `greenSettings::AtomicState.GreenSettings` is returned.
"""
function GreenSettings(settings::AtomicState.GreenSettings;
    nMax::Union{Nothing,Int64}=nothing,                     lValues::Union{Nothing,Array{Int64,1}}=nothing,
    dampingTau::Union{Nothing,Float64}=nothing,             printBefore::Union{Nothing,Bool}=nothing,
    levelSelection::Union{Nothing,LevelSelection}=nothing,  scField::Union{Nothing,AbstractScField}=nothing)

    if  isnothing(nMax)            nMaxx           = settings.nMax            else  nMaxx           = nMax            end
    if  isnothing(lValues)         lValuesx        = settings.lValues         else  lValuesx        = lValues         end
    if  isnothing(dampingTau)      dampingTaux     = settings.dampingTau      else  dampingTaux     = dampingTau      end
    if  isnothing(printBefore)     printBeforex    = settings.printBefore     else  printBeforex    = printBefore     end
    if  isnothing(levelSelection)  levelSelectionx = settings.levelSelection  else  levelSelectionx = levelSelection  end
    if  isnothing(scField)         scFieldx        = settings.scField         else  scFieldx        = scField         end

    GreenSettings( nMaxx, lValuesx, dampingTaux, printBeforex, levelSelectionx, scFieldx )
end


function Base.show(io::IO, settings::AtomicState.GreenSettings) 
    println(io, "nMax:                     $(settings.nMax)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
    println(io, "dampingTau:               $(settings.dampingTau)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  AtomicState.GreenChannel`  
    ... defines a type for a single symmetry channel of an (approximate) Green (function) expansion.

    + symmetry          ::LevelSymmetry    ... Level symmetry of this part of the representation.
    + gMultiplet        ::Multiplet        ... Multiplet of (scattering) levels of this symmetry.
"""
struct GreenChannel 
    symmetry            ::LevelSymmetry
    gMultiplet          ::Multiplet
end   


"""
`AtomicState.GreenChannel()`  ... constructor for an `empty` instance of AtomicState.GreenChannel.
        A `greenChannel::AtomicState.GreenChannel` is returned.
"""
function GreenChannel()
    GreenChannel( LevelSymmetry(0, Basics.plus), ManyElectron.Multiplet() )
end


function Base.show(io::IO, channel::AtomicState.GreenChannel) 
    println(io, "symmetry:                $(channel.symmetry)  ")
    println(io, "gMultiplet:              $(channel.gMultiplet)  ")
end


"""
`struct  AtomicState.GreenExpansion  <:  AbstractRepresentationType`  
    ... defines a type to keep an (approximate) Green (function) expansion that is associated with a given set of reference
        configurations.

    + approach          ::AtomicState.AbstractGreenApproach  ... Approach used to approximate the representation.
    + excitationScheme  ::Basics.AbstractExcitationScheme    ... Applied excitation scheme w.r.t. refConfigs. 
    + levelSymmetries   ::Array{LevelSymmetry,1}             ... Total symmetries J^P to be included into Green expansion.
    + NoElectrons       ::Int64                              ... Number of electrons.
    + settings          ::AtomicState.GreenSettings          ... settings for the Green (function) expansion.
"""
struct GreenExpansion  <:  AbstractRepresentationType
    approach            ::AtomicState.AbstractGreenApproach
    excitationScheme    ::Basics.AbstractExcitationScheme 
    levelSymmetries     ::Array{LevelSymmetry,1}
    NoElectrons         ::Int64 
    settings            ::AtomicState.GreenSettings
end   


"""
`AtomicState.GreenExpansion()`  ... constructor for an `empty` instance of AtomicState.GreenExpansion.
        A `greenExpansion::AtomicState.GreenExpansion` is returned.
"""
function GreenExpansion()
    GreenExpansion( AtomicState.SingleCSFwithoutCI(), Basics.NoExcitationScheme(), LevelSymmetry[], 0, AtomicState.GreenSettings())
end


function Base.string(expansion::GreenExpansion)
    sa = "Green (function) expansion in $(expansion.approach) approach and for excitation scheme  $(expansion.excitationScheme)," *
            "\nincluding (Green function channels with) symmetries $(expansion.levelSymmetries):"
    return( sa )
end


function Base.show(io::IO, expansion::GreenExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end



### Representation ##################################################################################
"""
`struct  AtomicState.Representation`  
    ... a struct for defining an atomic state representation. Such representations often refer to approximate wave function approximations
    of
        one or several levels but may concern also a mean-field basis (for some multiplet of some given configurations) or Green functions,
        etc.

    + name             ::String                      ... to assign a name to the given model.
    + nuclearModel     ::Nuclear.Model               ... Model, charge and parameters of the nucleus.
    + grid             ::Radial.Grid                 ... The radial grid to be used for the computation.
    + refConfigs       ::Array{Configuration,1}      ... List of references configurations, at least 1.
    + repType          ::AbstractRepresentationType  ... Specifies the particular representation.
    + orbitalSource    ::Basics.AbstractOrbitalSource
        ... says WHICH one-electron orbital set the representation is generated on.  `SeparateOrbitals()` is the default and is what JAC has
            always done:  the generator runs its own SCF on `refConfigs`.  `GivenOrbitals(orbitals)` instead hands the generator a set the
            CALLER already holds, and every multiplet of the representation is then a CI on it, built by
            `Hamiltonian.performCIwithFrozenOrbitals` with no further orbital optimization.

            WHY IT IS HERE AS WELL AS ON `Atomic.Computation` (added 05-Nov-2026).  An explicit SECOND-ORDER amplitude is valid only
            on ONE orthonormal set, because `ManyElectron.matrixElement_HminusE` drops the overlap term `-E <r|s>`, which vanishes
            identically only then.  `Atomic.Computation` has carried `orbitalSource` since 02-Nov-2026, so the PROCESS half of such a
            computation could be put on one set -- but its Green space is a `Representation`, which could not be, and a caller who
            wanted both halves consistent had to leave the `generate` route and call `Hamiltonian.performCIwithFrozenOrbitals` by
            hand.  That is what every second-order harness in `tools/` and both second-order branches of `examples/example-Dj.jl`
            did.  Both halves are now reachable through the ordinary route.

            ONLY `SeparateOrbitals()` AND `GivenOrbitals(...)` MEAN ANYTHING HERE.  A representation has no initial/final pair, so
            `InitialStateOrbitals()` and `FinalStateOrbitals()` have nothing to refer to, and `Basics.generate` raises on them
            rather than quietly picking one of the two.  Not every representation can honour a given set either -- one whose whole
            purpose is to OPTIMIZE orbitals cannot -- and those generators raise as well;  see `Basics.generate` for which do what.
"""
struct  Representation
    name               ::String  
    nuclearModel       ::Nuclear.Model
    grid               ::Radial.Grid
    refConfigs         ::Array{Configuration,1} 
    repType            ::AbstractRepresentationType
    orbitalSource      ::Basics.AbstractOrbitalSource
end


"""
`AtomicState.Representation()`  ... constructor for an 'empty' instance of the a variable::AtomicState.Representation
        A `representation::AtomicState.Representation` is returned.
"""
function Representation()
    Representation("", Nuclear.Model(1.0), Radial.Grid(), ManyElectron.Configuration[], CiExpansion(), Basics.SeparateOrbitals())
end


"""
`AtomicState.Representation(name::String, nuclearModel::Nuclear.Model, grid::Radial.Grid, refConfigs::Array{Configuration,1},
                            repType::AbstractRepresentationType; orbitalSource::Basics.AbstractOrbitalSource=Basics.SeparateOrbitals())`
    ... the five-argument constructor, which is how a representation has always been written down and goes on meaning exactly what it
        meant: the generator runs its own SCF on `refConfigs`.  The sixth field is reached by the keyword, so that a caller who holds an
        orbital set already -- the one case where a second-order amplitude is valid at all -- can hand it over without the positional
        form.  A `rep::AtomicState.Representation` is returned.

            orbitals = res[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
            rep      = Representation(name, nm, grid, [Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                      MeanFieldMultiplet(MeanFieldSettings()); orbitalSource = Basics.GivenOrbitals(orbitals) )
"""
function Representation(name::String, nuclearModel::Nuclear.Model, grid::Radial.Grid, refConfigs::Array{Configuration,1},
                        repType::AbstractRepresentationType; orbitalSource::Basics.AbstractOrbitalSource=Basics.SeparateOrbitals())

    return( AtomicState.Representation(name, nuclearModel, grid, refConfigs, repType, orbitalSource) )
end


"""
`AtomicState.Representation(wa::Bool)`  
    ... to return a default representation::AtomicState.Representation, which the keyword copy-constructor above then
        modifies;  the argument wa is not read.  The examples below show the typical use.

  Example for the generation of a mean-field basis:

        name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
        grid        = Radial.Grid(true)
        nuclearM    = Nuclear.Model(8.)
        refConfigs  = [Configuration("[He] 2s^2 2p^4")]
        mfSettings  = MeanFieldSettings()
        Representation(name, nuclearM, grid, refConfigs, MeanFieldBasis(mfSettings) )
    
`AtomicState.Representation( ... example for the computation of a configuration-interaction (CI) expansion)`  

        name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
        grid        = Radial.Grid(true)
        nuclearM    = Nuclear.Model(8.)
        refConfigs  = [Configuration("[He] 2s^2 2p^4")]
        orbitals    = wb[ResultKeys.MeanFieldBasis].orbitals #   get a proper set of orbitals
        ciSettings  = CiSettings(true, false, Int64[], false, LevelSymmetry[] )
        from        = [Shell("2s")]
        to          = [Shell("2s"), Shell("2p")]
        excitations = RasStep(RasStep(), seFrom=from, seTo=to, deFrom=from, deTo=to, frozen=[Shell("1s")])
        Representation(name, nuclearM, grid, refConfigs, CiExpansion(orbitals, excitations, ciSettings) )
    
`AtomicState.Representation( ... example for the computation of a restricted-active-space (RAS) expansion)`  

        name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
        refConfigs  = [Configuration("[He] 2s^2")]
        rasSettings = RasSettings([1], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(), LevelSelection() )
        from        = [Shell("2s")]
        
        frozen      = [Shell("1s")]
        to          = [Shell("2s"), Shell("2p")]
        step1       = RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), 
                                frozen=deepcopy(frozen))

        append!(frozen, [Shell("2s"), Shell("2p")])
        append!(to,     [Shell("3s"), Shell("3p"), Shell("3d")])
        step2       = RasStep(step1; seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))

        append!(frozen, [Shell("3s"), Shell("3p"), Shell("3d")])
        append!(to,     [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
        step3       = RasStep(step2, seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))

        Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs, 
                        RasExpansion(LevelSymmetry(0, Basics.plus), 4, [step1, step2, step3], rasSettings) )
    
`AtomicState.Representation( ... example for the computation of Green(function) expansion)`  

        name            = "Lithium 1s^2 2s ground configuration"
        refConfigs      = [Configuration("[He] 2s")]
        levelSymmetries = [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)]
        greenSettings   = GreenSettings(5, [0, 1, 2], 0.01, true, false, Int64[])
        Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                        GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                        levelSymmetries, 3, greenSettings) ) 
                        
    Either example is adapted by overwriting the parameters that matter for the case at hand.
"""
function Representation(wa::Bool)    
    AtomicState.Representation()    
end


"""
`AtomicState.Representation(name::String, nuclearModel::Nuclear.Model, refConfigs::Array{Configuration,1},
                            repType::AbstractRepresentationType; printout::Bool=true,
                            orbitalSource::Basics.AbstractOrbitalSource=Basics.SeparateOrbitals())`
    ... constructor for which NO grid is given, so that the radial box is derived from the REFERENCE
        configurations by Basics.recommendedGrid; a rep::AtomicState.Representation is returned.

        THE REFERENCE CONFIGURATIONS ARE THE RIGHT THING TO SIZE THE BOX FROM, and the alternative -- handing
        the estimate every configuration the expansion will contain -- is actively wrong.  The estimate reads a
        shell by its (n,l) and cannot see the ROLE it plays: a spectroscopic 3s, occupied in a real 1s^2 2s 3s
        state, is diffuse and is bound by its own ionization potential, whereas a CORRELATION 3s added to
        1s^2 2s^2 is never occupied, describes the short-range correlation hole, and contracts onto the valence
        region.  Same quantum numbers, extent differing by an order of magnitude.  Fed the correlation layers of
        a beryllium RAS expansion, the estimate returns 66 a.u. where the reference alone asks 20.5, and a box
        that large starves the fixed B-spline basis precisely on the high-n orbitals the layers introduce.

        MEASURED on a three-layer beryllium expansion (reference, +2p, +3s3p3d) with a full EOL optimisation per
        layer, one grid per run:

            box sized from   rbox     final energy      dE(layer 3 - layer 2)
            refConfigs       20.9     -14.61925410      -0.000135
            all shells       66.7     -14.61911933      -0.001730
            Radial.Grid(true) 614     -14.61805822      +0.001205

        The reference-sized box gives the LOWEST energy, which at fixed layer structure is the variational
        criterion that decides this.  The uncorrelated first layer is box-independent to eight digits across
        that whole range, so nothing is being traded away.  And the layer-to-layer differences are not
        comparable across grids at all -- the last column changes sign, and a correlation layer that RAISES the
        energy is not physics -- so a layer difference should never be quoted without stating the grid.

        The correlation orbitals are still checked, not merely assumed to fit: Bsplines.checkOrbitalBox runs on
        the converged orbitals of every layer and reports if any of them reaches the wall.

    + name          ::String                             ... name of the representation.
    + nuclearModel  ::Nuclear.Model                      ... nuclear model.
    + refConfigs    ::Array{Configuration,1}             ... reference configurations, which set the box.
    + repType       ::AbstractRepresentationType         ... the expansion to be generated.
"""
function Representation(name::String, nuclearModel::Nuclear.Model, refConfigs::Array{Configuration,1},
                        repType::AbstractRepresentationType; printout::Bool=true,
                        orbitalSource::Basics.AbstractOrbitalSource=Basics.SeparateOrbitals())
    grid = Basics.recommendedGrid(refConfigs, nuclearModel, printout=printout)

    return( AtomicState.Representation(name, nuclearModel, grid, refConfigs, repType, orbitalSource) )
end



function Base.string(rep::Representation)
    sa = "Atomic representation:   $(rep.name) for Z = $(rep.nuclearModel.Z) and with reference configurations: \n   "
    for  refConfig  in  rep.refConfigs     sa = sa * string(refConfig) * ",  "     end
    return( sa )
end


function Base.show(io::IO, rep::Representation)
    sa = Base.string(rep);            print(io, sa, "\n")
    println(io, "representation type:   $(rep.repType)  ")
    println(io, "nuclearModel:          $(rep.nuclearModel)  ")
    println(io, "grid:                  $(rep.grid)  ")
    if  rep.orbitalSource != Basics.SeparateOrbitals()
        println(io, "orbitalSource:         $(rep.orbitalSource)  ")
    end
end


# The measured anchors used by `sizeRasStep` to BRACKET a layer: (nCsf, peak RSS in GB, wall clock in s; 0 = not
# measured).  Ti III with 3p opened, one process each under /usr/bin/time -v, 14-Sep-2026.  THE 11 298 ENTRY IS
# PROVISIONAL -- that run did not complete, and it is the single point that makes the series look unfittable, so it
# should be re-measured before any law is rejected or built on it.
const ANCHORS = [ (887, 1.6, 98.0), (3167, 2.1, 777.0), (11298, 12.8, 0.0), (25085, 29.7, 36000.0) ]

fmtTime(t::Float64) = t < 120. ? @sprintf("%.0f s", t) : (t < 7200. ? @sprintf("%.0f min", t/60) : @sprintf("%.1f h", t/3600))


"""
`AtomicState.sizeRasStep(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1},
                         step::AtomicState.RasStep; printout::Bool=true)`
    ... counts what a RAS step would cost BEFORE any of it is run, so that a layer can be judged against the machine
        it is meant to run on. Only the CSF expansion is built -- no orbitals, no SCF, no angular coefficients -- so
        this is seconds even for a layer that would take hours. A NamedTuple
        `(nCsf, nSub, lMax, blocks, sumN2, maxN2, blockShare)` is returned.

        WHAT IT REPORTS, AND WHY NOT MORE. Everything here is either an exact COUNT or a RATIO of two counts, and
        that restraint is deliberate. Everything an EOL/RAS layer builds per CSF PAIR lives inside one symmetry
        block, so the pair work is `sum_b n_b^2` over the blocks and NOT `nCsf^2`; with the stores of all blocks
        resident, as they are today, that sum is what is paid, while a route holding ONE block at a time would pay
        `max_b n_b^2` instead. Both are counted here, and `blockShare = max_b n_b^2 / sum_b n_b^2` is what such a
        route would save -- 31 % on a Ti III layer of 25 085 CSFs, measured 14-Sep-2026.

        **A RATIO IS QUOTED WHERE AN ABSOLUTE FIGURE IS NOT, BECAUSE THE UNKNOWN CONSTANT CANCELS IN IT.** Four peak-RSS
        anchors on Ti III layers with 3p opened -- 887 CSFs at 1.647 GB, 3 167 at 2.149, 11 298 at 12.784, 25 085 at
        29.73 -- are fitted by NO law in the CSF count alone: the space-dependent cost per CSF runs 1.7e-4, 2.1e-4,
        1.0e-3, 1.1e-3 GB, a fivefold jump between the second and third and then flat, because those four are not four
        sizes of one thing but layers that add 4p, then 4d, then 4f. Two successive fitted laws disagreed by 24x in
        slope for exactly that reason. So the anchors are printed as MEASUREMENTS to interpolate between, and this
        routine does not pretend to a law it does not have.

        THE TWO DESIGN LEVERS IT IS MEANT TO INFORM, both measured rather than argued:
        - the angular coefficient count per CSF is set mainly by l, at roughly 2.7x per step -- about 55 entries per
          CSF for p, 130-172 for d, 340-389 for f -- and grows with the subshell count as well, so `lMax` and `nSub`
          are reported beside the CSF count;
        - CONCENTRATING the open electrons beats splitting them: 3d^5 gives 37 CSFs and 6 363 entries against
          3p^3 3d^2 with 141 and 33 545, a factor 5.3 for the same five open-shell electrons.

        AND THE J DISTRIBUTION IS NOT A LEVER, though it looks like one in the table below: it follows from the
        configurations and the angular coupling and is not chosen. The only control there is which symmetries are
        asked for at all, which drops whole blocks.
"""
function sizeRasStep(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1},
                     step::AtomicState.RasStep; printout::Bool=true)
    basis  = Basics.generateBasis(refConfigs, symmetries, step)
    nCsf   = length(basis.csfs);    nSub = length(basis.subshells)
    lMax   = isempty(basis.subshells) ? 0 : maximum( Basics.subshell_l(sh)  for sh in basis.subshells )

    blocks = Dict{LevelSymmetry,Int64}()
    for  csf  in  basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity);    blocks[sym] = get(blocks, sym, 0) + 1
    end
    sumN2 = sum( Float64(n)^2  for (_, n) in blocks; init=0.0 )
    maxN2 = isempty(blocks) ? 0.0 : Float64(maximum(values(blocks)))^2
    share = sumN2 > 0. ? maxN2 / sumN2 : 0.

    if  printout
        println("\n>> RAS step sizing:  $nCsf CSFs over $nSub subshells (l_max = $lMax), $(length(blocks)) symmetry blocks.")
        println(">>   J^P            CSFs           n_b^2      share of the pair work")
        for  (sym, n)  in  sort(collect(blocks), by = x -> -x[2])
            println(">>   " * rpad(string(sym), 12) * lpad(string(n), 8) * lpad(@sprintf("%14.4g", Float64(n)^2), 16) *
                    lpad(@sprintf("%8.1f %%", 100*Float64(n)^2/max(sumN2,1.0)), 16))
        end
        println(">>   sum n_b^2 = " * @sprintf("%.4g", sumN2) * ",  largest block = " * @sprintf("%.4g", maxN2) *
                @sprintf(" (%.1f %% of the sum)", 100*share))
        println(">>")
        println(">>   ROUTE COMPARISON, as a ratio, which is the part that does not depend on an unmeasured constant:")
        println(">>     all symmetry blocks resident (today)        1.00 x the pair work")
        println(">>     one block at a time                         " * @sprintf("%.2f", share) *
                " x the pair work, for 26-32 % more run time")
        println(">>")
        println(">>   MEASURED ANCHORS -- peak RSS and wall clock of Ti III layers with 3p opened, one process each")
        println(">>   (14-Sep-2026).  These are MEASUREMENTS, not a fitted curve:")
        println(">>        nCsf      peak RSS     wall clock")
        println(">>         887       1.6 GB          98 s")
        println(">>       3 167       2.1 GB         777 s")
        println(">>      11 298      12.8 GB          --        (PROVISIONAL: that run did not complete)")
        println(">>      25 085      29.7 GB       > 10 h       (stopped by its own time limit, not by memory)")
        println(">>")
        (loA, hiA) = (0, 0)
        for (n, g, t) in ANCHORS
            n <= nCsf  &&  (loA = findfirst(x -> x[1] == n, ANCHORS))
            hiA == 0  &&  n >= nCsf  &&  (hiA = findfirst(x -> x[1] == n, ANCHORS))
        end
        if      loA != 0  &&  loA == hiA
            # the step lands ON an anchor, so quote what was actually measured rather than a bracket around it
            (n1, g1, t1) = ANCHORS[loA]
            println(">>   THIS STEP IS THE ANCHOR AT $n1 CSFs, which was measured at " * @sprintf("%.1f GB", g1) *
                    (t1 > 0. ? " in " * fmtTime(t1) : "") * ".")
        elseif  loA != 0  &&  hiA != 0  &&  loA != hiA
            (n1, g1, t1) = ANCHORS[loA];    (n2, g2, t2) = ANCHORS[hiA]
            println(">>   THIS STEP SITS BETWEEN TWO OF THEM, so expect roughly")
            println(">>     memory   between " * @sprintf("%.1f", g1) * " and " * @sprintf("%.1f GB", g2) *
                    "   (bracketed by the anchors at $n1 and $n2 CSFs)")
            t1 > 0. && t2 > 0. &&
                println(">>     time     between " * fmtTime(t1) * " and " * fmtTime(t2))
        elseif  hiA == 0
            (n1, g1, t1) = ANCHORS[end]
            println(">>   THIS STEP IS LARGER THAN EVERY ANCHOR (the largest is $n1 CSFs at " *
                    @sprintf("%.1f GB", g1) * "), so the figures above are a LOWER bound and nothing here")
            println(">>   should be read as a prediction:  extrapolating these anchors is exactly what produced two")
            println(">>   successive cost laws that disagreed by 24x in slope.")
        else
            println(">>   THIS STEP IS SMALLER THAN EVERY ANCHOR, so expect the ~1.5 GB baseline to dominate:  a JAC")
            println(">>   process costs about that before any layer is built.")
        end
        println(">>")
        println(">>   AND THE RATIO ABOVE IS THE FIRM PART.  An absolute figure needs a constant that has NOT been")
        println(">>   pinned down -- the cost per CSF runs 1.7e-4, 2.1e-4, 1.1e-3 GB across these anchors -- whereas")
        println(">>   the one-block-at-a-time ratio carries that same constant in numerator and denominator, where it")
        println(">>   cancels.  Trust the " * @sprintf("%.2f", share) * " x;  treat the GB as a bracket.")
    end

    return( (nCsf=nCsf, nSub=nSub, lMax=lMax, blocks=blocks, sumN2=sumN2, maxN2=maxN2, blockShare=share) )
end


"""
`AtomicState.tryRun(rep::AtomicState.Representation; nSample::Int64=2000, printout::Bool=true)`
    ... a TRY RUN of a RAS ladder: for every step it reports what the production job would cost in memory and in
        time, WITHOUT running any of it, so that a ladder can be discussed and decided before a machine is
        committed to it. An `Array{NamedTuple,1}`, one entry per step, is returned.

        WHY IT TAKES THE REPRESENTATION AND NOT A SETTINGS FLAG. It is handed THE VERY OBJECT that `generate`
        would be handed, so the estimate cannot describe a different calculation from the production run -- which
        is the one advantage a `tryRun` field in `AsfSettings` would have had, obtained here without any of its
        costs. A flag there would have to make `performSCF` return something other than a `Multiplet` (a fake one,
        or `nothing`, or a union every caller must then handle); `AsfSettings` is threaded through every module,
        so a stray `true` left in a copy-constructor would silently produce no physics everywhere and look like a
        convergence failure; and "estimate instead of compute" is a verb, not a property of the Hamiltonian. With
        a separate entry point there is no flag to default, and the production path is untouched.

        HOW THE NUMBERS ARE OBTAINED -- COUNTED WHERE POSSIBLE, SAMPLED WHERE NOT. The CSF expansion of each step
        is generated (cheap: no orbitals, no SCF) and its symmetry blocks are counted exactly. Everything an
        EOL/RAS step builds per CSF PAIR lives inside one block, so the pair work is `sum_b n_b^2` and NOT
        `nCsf^2`; with the stores of all blocks resident, as they are today, that sum is what is paid, while a
        route holding one block at a time would pay `max_b n_b^2`. Those are exact. The only unmeasured quantity
        is the cost of ONE pair, and that is SAMPLED: `nSample` pairs are drawn at RANDOM from each block and
        their angular coefficients actually computed, on this system's own shells and occupations.

        RANDOM PAIRS, NOT A CONTIGUOUS SUB-BLOCK, AND THE DIFFERENCE MATTERS. Cost per pair FALLS as a block grows
        -- measured 14-Sep-2026, 61.5 coefficients per pair at 3d^2 against 29.1 at 3d^5 -- because most pairs of
        a large block are far apart in excitation and carry nothing at all (61 % were empty in one measured
        block). A contiguous subset over-samples the near-diagonal pairs and would over-predict; a random draw
        samples the true distribution, so the mean is unbiased. The scatter over the sample is reported beside the
        estimate, since 60 % of draws returning zero makes a small sample noisy.

        WHAT IS SOLID AND WHAT IS NOT. The per-ITERATION figures rest on exact counts times a measured per-pair
        cost, and are the trustworthy part. The TOTAL needs the iteration count, which is genuinely
        unpredictable -- so `maxIterationsScf` from the ladder's own settings is printed as the multiplier it is,
        and the total is given as "up to", not as a forecast.

        AND THE ITERATION BUDGET IS PRINTED AS A DIAL, because for a large open-shell system it may be the only
        one left to turn. The cost is linear in it, and the rotation route descends MONOTONICALLY, so a run
        stopped early returns a variational UPPER BOUND on the energy: the orbitals are less converged, not
        wrong. Trading convergence for feasibility is therefore a legitimate choice, and the row exists so that
        the choice can be made on numbers -- "60 iterations is 22 h but 15 is 5 h" -- rather than discovered
        after a job has been queued.
"""
function tryRun(rep::AtomicState.Representation; nSample::Int64=2000, printout::Bool=true)
    if  !(rep.repType isa AtomicState.RasExpansion)
        error("AtomicState.tryRun is implemented for a RasExpansion; the representation given carries " *
              "$(typeof(rep.repType)).  Other representation types follow by the same pattern.")
    end
    # Through the accessor, so that a new route type needs no change here; every route carries its own natural budget.
    repType = rep.repType;    maxIter = Basics.maxIterations(repType.settings.scfRoute)
    results = NamedTuple[]

    if  printout
        println("\n", "="^124)
        println("TRY RUN -- what this RAS ladder would cost, WITHOUT running it.   $(length(repType.steps)) step(s), " *
                "up to $maxIter SCF iterations each.")
        println("="^124)
    end

    for  (istep, step)  in  enumerate(repType.steps)
        tGen  = @elapsed (basis = Basics.generateBasis(rep.refConfigs, repType.symmetries, step))
        nCsf  = length(basis.csfs);    nSub = length(basis.subshells)
        blocks = Dict{LevelSymmetry,Int64}()
        for  csf  in  basis.csfs
            sym = LevelSymmetry(csf.J, csf.parity);    blocks[sym] = get(blocks, sym, 0) + 1
        end
        sumN2 = sum( Float64(n)^2  for (_, n) in blocks; init=0.0 )
        maxN2 = isempty(blocks) ? 0.0 : Float64(maximum(values(blocks)))^2

        # SAMPLE the per-pair cost on this system's own CSFs.  THE DRAW MUST BE WITHIN ONE BLOCK, and getting
        # that wrong is the obvious trap: `sum_b n_b^2` counts only pairs that SHARE a symmetry, so sampling (r,s)
        # from the whole CSF list scores zero on every cross-block draw and dilutes the mean by whatever fraction
        # those are.  Measured before the fix on a 25 085-CSF step: 0.1 coefficients per pair, against about 3
        # when the draw is done properly -- a thirtyfold under-estimate of the store.  A block is therefore chosen
        # first, with probability proportional to n_b^2 (its share of the work), and r and s are drawn inside it.
        idxOfBlock = Dict{LevelSymmetry,Array{Int64,1}}()
        for  (i, csf)  in  enumerate(basis.csfs)
            sym = LevelSymmetry(csf.J, csf.parity);    push!(get!(idxOfBlock, sym, Int64[]), i)
        end
        blockList = collect(keys(blocks));    wBlock = [ Float64(blocks[b])^2 / max(sumN2,1.0)  for b in blockList ]
        cumW = cumsum(wBlock)
        drawPair() = begin
            u = rand();    ib = findfirst(x -> x >= u, cumW);    ib === nothing && (ib = length(cumW))
            idx = idxOfBlock[blockList[ib]]
            (basis.csfs[idx[rand(1:length(idx))]], basis.csfs[idx[rand(1:length(idx))]])
        end
        countPair(csfR, csfS) = begin
            m = 0
            for  cf  in  SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus),
                                                               csfR, csfS, basis.subshells)      m += 1   end
            for  cf  in  SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus),
                                                         csfR, csfS, basis.subshells)            m += 1   end
            m
        end
        # WARM UP FIRST, AND TIME THE LOOP AS A WHOLE.  Wrapping `@elapsed` around each individual pair charges
        # Julia's JIT compilation of the angular routines to whichever pair happened to be first -- seconds of it
        # -- and that one pair then sets the per-pair cost.  Measured 14-Sep-2026: two runs of the SAME step gave
        # 22 min and 7.5 h per iteration, a twentyfold swing, from exactly this.  One warm-up call outside the
        # timing and one `@elapsed` around the whole loop removes it.
        entries = Int64[];    tSample = 0.0
        # AND THE LOOP IS TIMED TWICE, THE FIRST PASS DISCARDED.  A handful of warm-up draws is NOT enough here:
        # `SpinAngular` specialises on the occupation pattern of the pair it is given, so fresh specialisations go
        # on being compiled well into the sample and land in the timing.  Measured 14-Sep-2026, the same 887-CSF
        # step timed 2 s and then 27 s of angular build in two runs -- a thirteenfold swing -- while the STORE,
        # which does not depend on timing, reproduced at 8 MB both times.  Running the whole loop twice and
        # keeping the second pass removes it;  the counts come from the second pass too, so both are consistent.
        for  pass = 1:2
            passEntries = Int64[]
            tPass = @elapsed begin
                for  k = 1:nSample
                    (csfR, csfS) = drawPair();    push!(passEntries, countPair(csfR, csfS))
                end
            end
            entries = passEntries;    tSample = tPass
        end
        nDrawn   = max(length(entries), 1)
        perPair  = sum(entries) / nDrawn
        sdPair   = nDrawn > 1 ? sqrt(sum((e - perPair)^2 for e in entries) / (nDrawn - 1)) : 0.0
        # THE STANDARD ERROR OF THE MEAN IS WHAT THE ESTIMATE INHERITS, not the spread of the sample.  With most
        # pairs empty and a heavy tail on the rest, the two differ by sqrt(nSample) and only the first says how
        # well the total is known;  it is reported so that a noisy estimate announces itself.
        sePair   = sdPair / sqrt(nDrawn)
        relErr   = perPair > 0. ? sePair / perPair : 0.
        sPerPair = tSample / nDrawn
        fillPct  = 100 * count(!iszero, entries) / nDrawn

        # 16 bytes per stored entry is the flat-CSR cost of the present store: an Int32 label index, an Int32
        # value index, and the amortised share of the interning tables and the pair pointer.
        bytesPerEntry = 16.0
        memAll = sumN2 * perPair * bytesPerEntry / 1024^3
        memOne = maxN2 * perPair * bytesPerEntry / 1024^3
        fmtMem(g) = g >= 1.0 ? @sprintf("%.2f GB", g) : (g >= 1.0e-3 ? @sprintf("%.0f MB", g*1024) : @sprintf("%.1f kB", g*1024^2))
        # THE ANGULAR BUILD IS THE SAME WORK ON EITHER ROUTE -- every block's coefficients are needed each
        # iteration whichever way they are stored.  WHAT THE ROUTE CHANGES IS HOW OFTEN IT IS PAID: holding all
        # blocks pays it ONCE and reuses them, holding one at a time pays it EVERY iteration.  So that route buys
        # memory WITH time and can NEVER be faster;  an earlier version of this routine applied the memory
        # formula (max_b n_b^2) to the time column and so showed the cheaper-memory route as also the quicker,
        # which is the trade backwards.
        tBuild = sumN2 * sPerPair
        # A FULL SCF ITERATION IS LARGER THAN ITS ANGULAR PART, by a measured factor: 15.26 s of angular build
        # against 47.7 s per iteration at 3167 CSFs, and 4.78 s against 8.4 s at 887 -- so the build is roughly a
        # third of an iteration, and the penalty for holding one block at a time is about +30 % wall clock rather
        # than a doubling.  This routine times only the angular work and cannot see the radial integrals, the
        # orbital update or the diagonalisation, so that factor is carried explicitly.
        iterOverAngular = 3.0
        tIter  = tBuild * iterOverAngular
        tAllT  = tBuild + tIter * maxIter                      # build once, then iterate
        tOneT  = (tBuild + tIter) * maxIter                    # rebuild every iteration

        if  printout
            println("\n>> step $istep:  $nCsf CSFs, $nSub subshells, $(length(blocks)) blocks, largest " *
                    @sprintf("%.0f %%", 100*maxN2/max(sumN2,1.0)) * " of the pair work;  CSF list built in " *
                    @sprintf("%.1f s", tGen) * ";  sampled $nDrawn random pairs, " *
                    @sprintf("%.2f +- %.2f", perPair, sePair) * " coefficients each (s.e. of the mean, " *
                    @sprintf("%.0f %%", 100*relErr) * "), " * @sprintf("%.0f %%", fillPct) * " non-empty." *
                    (relErr > 0.20 ? "  ** NOISY: raise nSample **" : ""))
            println("   " * rpad("route", 38) * rpad("store", 12) * rpad("angular build", 14) *
                    rpad("paid", 18) * "total, up to $maxIter iterations")
            println("   " * rpad("all blocks resident (today)", 38) * rpad(fmtMem(memAll), 12) *
                    rpad(fmtTime(tBuild), 14) * rpad("once", 18) * fmtTime(tAllT))
            println("   " * rpad("one block at a time (not yet built)", 38) * rpad(fmtMem(memOne), 12) *
                    rpad(fmtTime(tBuild), 14) * rpad("every iteration", 18) * fmtTime(tOneT))
            println("   " * rpad("", 38) *
                    @sprintf("-> one block at a time would cost %.0f %% more time and keep %.0f %% of the store",
                             100*(tOneT/max(tAllT,1e-30) - 1), 100*memOne/max(memAll,1e-30)))
            # THE ITERATION BUDGET IS A DIAL, AND FOR A LARGE SYSTEM IT MAY BE THE ONLY ONE LEFT.  The rotation
            # route descends monotonically, so a run stopped early returns a VARIATIONAL UPPER BOUND on the
            # energy -- the orbitals are less good, not wrong -- and trading convergence for feasibility is a
            # legitimate choice rather than a fudge.  The cost is linear in the budget, so it is shown as a dial.
            local budgets = sort(unique(Int64[ b  for b in [maxIter, maxIter ÷ 2, maxIter ÷ 4, 10, 5]  if 1 <= b <= maxIter ]), rev=true)
            println("   " * rpad("iteration budget (today's route)", 38) *
                    join([ @sprintf("%d -> %s", b, fmtTime(tBuild + tIter*b))  for b in budgets ], "   |   "))
        end

        push!(results, (step=istep, nCsf=nCsf, nSub=nSub, blocks=length(blocks), sumN2=sumN2, maxN2=maxN2,
                        perPair=perPair, sdPerPair=sdPair, fillPercent=fillPct, memAllGB=memAll, memOneGB=memOne,
                        secAngularBuild=tBuild, secPerIteration=tIter, secTotalAll=tAllT,
                        secTotalOne=tOneT, maxIterations=maxIter))
    end

    if  printout
        println("\n>> THE PER-ITERATION FIGURES ARE THE SOLID ONES -- exact pair counts times a per-pair cost measured")
        println(">> on this very system.  The totals assume the full iteration budget and are an UPPER bound, not a")
        println(">> forecast:  a step that converges in five iterations costs a quarter of what is shown.  The store")
        println(">> figure covers the angular coefficient store only;  a running job also carries the ~1.5 GB JAC")
        println(">> baseline, the CI matrix of one block at a time, and the radial caches.")
        println(">> THE TIME COLUMN IS BUILT FROM THE ANGULAR WORK, sampled here, times a measured factor of 3 for the")
        println(">> rest of an SCF iteration.  Against the one step run to completion it comes out LOW by about 1.75x")
        println(">> -- helmtop did 16 iterations of the 25 085-CSF step in under ten hours -- so read the totals as a")
        println(">> FLOOR.  Item 30 can never be faster:  it does the same angular work and merely pays it per")
        println(">> iteration instead of once, which is how it buys the store back.  That route is NOT IMPLEMENTED yet;")
        println(">> its row says what it WOULD cost, so that the trade can be judged before anyone builds it.")
        println("")
        println(">> " * "-"^116)
        println(">> HOW TO MAKE A LARGE RUN FIT, in the order worth trying.  These are Julia and JAC facts that are easy")
        println(">> to forget between calculations, so they are printed where a large job is being planned:")
        println(">>")
        println(">>  1. CUT THE PROBLEM, NOT THE HEAP.  A JAC EOL solve ALLOCATES far more than it holds -- measured")
        println(">>     14-Sep-2026, 133.5 GB allocated over a 130 s solve whose live data is about 2 GB, with 12 % of")
        println(">>     the wall clock already in garbage collection -- so what a job occupies is mostly heap grown to")
        println(">>     absorb churn, which Julia does not hand back.  THAT MAKES `--heap-size-hint` LOOK LIKE THE")
        println(">>     ANSWER AND IT IS NOT:  measured 14/15-Sep on the 25 085-CSF Ti III row, `--heap-size-hint=12G`")
        println(">>     still peaked at 26.6 GB (against 29.7 GB uncapped, a 10 % saving) and the job was then KILLED")
        println(">>     at two hours.  The flag is a HINT -- the collector tries to stay under it and grows anyway when")
        println(">>     the allocation rate demands.  Try it if you like, but do NOT size a machine on it.")
        println(">>")
        println(">>  2. CUT THE ITERATION BUDGET, as the dial above shows.  The rotation route descends monotonically,")
        println(">>     so a run stopped early returns a variational UPPER BOUND -- less converged orbitals, not wrong")
        println(">>     ones -- and the cost is linear in the budget.  For a large open-shell system this is often the")
        println(">>     only practical route, and it is an honest one.")
        println(">>")
        println(">>  3. CONCENTRATE THE OPEN ELECTRONS rather than splitting them over shells.  Measured: 3d^5 gives 37")
        println(">>     CSFs and 6 363 angular coefficients against 3p^3 3d^2 with 141 and 33 545 -- a factor 5.3 for the")
        println(">>     same five electrons, because two open shells in jj-coupling generate far more distinct subshell")
        println(">>     quadruples than their lower ranks save.")
        println(">>")
        println(">>  4. ASK FOR FEWER SYMMETRIES if the physics allows:  each J^P is an independent block, and dropping")
        println(">>     one removes its whole n_b^2 share.  The J DISTRIBUTION itself is not a lever -- it follows from")
        println(">>     the configurations -- but which symmetries are requested is.")
        println(">>")
        println(">>  5. RUN THIS TRY RUN FIRST, on any ladder that looks expensive.  It costs a minute against hours.")
        println(">> " * "-"^116)
    end

    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end

end # module