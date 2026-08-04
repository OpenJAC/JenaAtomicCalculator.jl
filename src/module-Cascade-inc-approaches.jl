
## Cascade approaches: the CONDITIONS that define them.
##
## A cascade approach fixes exactly three things -- how the atomic levels are represented, which bound orbitals
## they are built from, and how the continuum orbitals of the free electron are generated. It fixes NOTHING
## else. Which processes are included, which multipoles are applied, whether shake configurations are added,
## and whether Breit/QED enter the CI matrix all remain properties of the cascade SCHEME or of AsfSettings,
## where they already live. Keeping those axes apart is deliberate: the historical ShakeMCA approach duplicated
## the shakeFromShells/shakeToShells fields that Cascade.StepwiseDecayScheme already carried, and "only E1
## transitions" was written into an approach while being a property of the requested processes -- in both cases
## two places could disagree about the same question, with nothing to arbitrate.
##
## The three axes are declared below as separate type hierarchies, so that an approach is a NAME for one
## particular combination rather than a thing in its own right. Every dispatch site asks for the conditions,
## never for the approach name; a combination that no approach happens to name is therefore still expressible
## and still runs, which is what makes controlled one-axis-at-a-time experiments possible.


"""
`abstract type Cascade.AbstractLevelRepresentation`
    ... defines an abstract type to specify how the atomic levels of a cascade block are represented.

    + struct SingleCSF        ... every level is a single CSF; no configuration mixing of any kind.
    + struct BlockCI          ... configuration mixing WITHIN each block (intermediate coupling), but no
                                  mixing between different blocks.
    + struct GroupedBlockCI   ... configuration mixing within blocks that may themselves comprise SEVERAL
                                  configurations, as grouped by the user; still no mixing between blocks.
"""
abstract type  AbstractLevelRepresentation                                    end
struct         SingleCSF        <:  Cascade.AbstractLevelRepresentation       end
struct         BlockCI          <:  Cascade.AbstractLevelRepresentation       end
struct         GroupedBlockCI   <:  Cascade.AbstractLevelRepresentation       end


"""
`abstract type Cascade.AbstractBoundOrbitals`
    ... defines an abstract type to specify which set of bound one-electron orbitals the cascade blocks are
        built from.

    + struct GlobalOrbitals       ... ONE common set for the whole cascade, generated from the mean field of
                                      the initial atom or ion and re-used for every charge state.
    + struct ChargeStateOrbitals  ... one common set per charge state, generated from the mean field of that
                                      charge state and shared by all of its blocks.
"""
abstract type  AbstractBoundOrbitals                                          end
struct         GlobalOrbitals       <:  Cascade.AbstractBoundOrbitals         end
struct         ChargeStateOrbitals  <:  Cascade.AbstractBoundOrbitals         end


"""
`abstract type Cascade.AbstractContinuumOrbitals`
    ... defines an abstract type to specify how the continuum orbital of the emitted (Auger) electron is
        generated.

    + struct StepAveragedContinuum       ... ONE continuum orbital set per cascade step, generated at the
                                             step's mean transition energy and shared by all fine-structure
                                             transitions of that step.
    + struct TransitionResolvedContinuum ... one continuum orbital per fine-structure transition, at that
                                             transition's own energy. This matters beyond accuracy: near
                                             threshold, a step-averaged energy can leave a channel open that
                                             is in fact closed (or vice versa), which adds or removes a whole
                                             branch of the decay tree rather than just mis-scaling a rate.
"""
abstract type  AbstractContinuumOrbitals                                              end
struct         StepAveragedContinuum        <:  Cascade.AbstractContinuumOrbitals     end
struct         TransitionResolvedContinuum  <:  Cascade.AbstractContinuumOrbitals     end


"""
`struct  Cascade.RefinedSCA  <:  Basics.AbstractCascadeApproach`
    ... a cascade approach that keeps SCA's level representation and orbitals but resolves the continuum
        orbitals per fine-structure transition; see Cascade.conditions for the full set of conditions.

        Note this concrete type is declared HERE and not in Basics: only the abstract supertype
        Basics.AbstractCascadeApproach has to live in Basics, because DecayYield.Settings needs it as a
        compile-time field type and module-DecayYield.jl is included before module-Cascade.jl. Any subtype may
        be declared wherever it belongs, which keeps Basics from accumulating cascade vocabulary.
"""
struct   RefinedSCA  <:  Basics.AbstractCascadeApproach   end


"""
`struct  Cascade.CascadeApproachConditions`
    ... the three conditions that together define a cascade approach.

    + levels     ::Cascade.AbstractLevelRepresentation  ... representation of the levels of a block.
    + bound      ::Cascade.AbstractBoundOrbitals        ... origin of the bound one-electron orbitals.
    + continuum  ::Cascade.AbstractContinuumOrbitals    ... generation of the free-electron orbitals.
"""
struct  CascadeApproachConditions
    levels          ::Cascade.AbstractLevelRepresentation
    bound           ::Cascade.AbstractBoundOrbitals
    continuum       ::Cascade.AbstractContinuumOrbitals
end


# `Base.show(io::IO, conditions::Cascade.CascadeApproachConditions)`  ... prepares a proper printout.
function Base.show(io::IO, conditions::Cascade.CascadeApproachConditions)
    println(io, "levels:                     $(conditions.levels)  ")
    println(io, "bound:                      $(conditions.bound)  ")
    println(io, "continuum:                  $(conditions.continuum)  ")
end


"""
`Cascade.conditions(approach::Basics.AbstractCascadeApproach)`
    ... returns the Cascade.CascadeApproachConditions that define the given cascade approach. Each approach
        escalates exactly one axis with respect to the previous one, so that the cost and the accuracy of a
        cascade can be raised one property at a time:

        + AverageSCA  ... SingleCSF      / GlobalOrbitals      / StepAveragedContinuum
        + SCA         ... BlockCI        / ChargeStateOrbitals / StepAveragedContinuum
        + RefinedSCA  ... BlockCI        / ChargeStateOrbitals / TransitionResolvedContinuum
        + UserMCA     ... GroupedBlockCI / ChargeStateOrbitals / TransitionResolvedContinuum

        The SCA -> RefinedSCA boundary is the one real cost cliff of this hierarchy. Switching SingleCSF to
        BlockCI was measured at only ~1.1x for a 34-CSF configuration (Fe 3d^6, 04-Aug-2026), since the CSFs
        spread over many small J^P blocks, the upper-triangle shortcut halves the work and the radial-integral
        cache absorbs much of the rest. Resolving the continuum per transition instead of per step, by
        contrast, replaces one continuum-orbital generation per step by one per fine-structure transition and
        is expected to dominate any Auger-rich cascade. Separating the two lets the cheap accuracy be taken
        without paying for the expensive one -- they were bundled into a single step of the published
        hierarchy, which is why its "each approach costs almost nothing more than the previous" claim holds
        for one half and not for the other.

        A fifth tier, RelaxedMCA, is deliberately NOT defined: it would relax the bound orbitals
        independently for every block and connect initial and final states through the (validated)
        BiOrthogonal module, together with a transition-resolved continuum. It is the natural end point of
        this hierarchy and the enabling piece exists, but the cost of a biorthogonal transformation at cascade
        scale has never been measured, so no name is offered that cannot yet be selected.

        An  conditions::Cascade.CascadeApproachConditions  is returned.
"""
function conditions(approach::Basics.AbstractCascadeApproach)
    if      approach == Basics.AverageSCA()
        return( CascadeApproachConditions( SingleCSF(),      GlobalOrbitals(),      StepAveragedContinuum()       ) )
    elseif  approach == Basics.SCA()
        return( CascadeApproachConditions( BlockCI(),        ChargeStateOrbitals(), StepAveragedContinuum()       ) )
    elseif  approach == Cascade.RefinedSCA()
        return( CascadeApproachConditions( BlockCI(),        ChargeStateOrbitals(), TransitionResolvedContinuum() ) )
    elseif  approach == Basics.UserMCA()
        return( CascadeApproachConditions( GroupedBlockCI(), ChargeStateOrbitals(), TransitionResolvedContinuum() ) )
    else
        error("Unsupported cascade approach: $approach.")
    end
end


"""
`Cascade.asfSettingsForConditions(conditions::Cascade.CascadeApproachConditions, settings::AsfSettings)`
    ... returns the AsfSettings with which a cascade block must be represented under the given conditions.
        The level representation is a property of the approach and not a free setting, so the returned
        settings override settings.eeInteractionCI accordingly: DiagonalCoulomb() makes Hamiltonian.setupMatrix
        skip every off-diagonal element -- before any spin-angular or radial-integral work, so it is also the
        cheap path that SingleCSF is meant to be -- while CoulombInteraction() gives the intermediate coupling
        of BlockCI/GroupedBlockCI. An  AsfSettings  is returned.
"""
function asfSettingsForConditions(conditions::Cascade.CascadeApproachConditions, settings::AsfSettings)
    if      conditions.levels == SingleCSF()    return( AsfSettings(settings; eeInteractionCI=DiagonalCoulomb())    )
    else                                        return( AsfSettings(settings; eeInteractionCI=CoulombInteraction()) )
    end
end


"""
`Cascade.validateConditions(conditions::Cascade.CascadeApproachConditions)`
    ... checks, at cascade set-up time, that the requested conditions can actually be realized, and raises an
        error otherwise. This is deliberately done up front: a cascade that is going to fail should fail
        before hours of computation, not in the middle of some step. Nothing is returned.

        Exchange between the free and the bound electrons is NOT checked here, and deliberately so: every
        Continuum.Settings built anywhere in this module hardwires includeExchange = false, and a user who
        edits that is refused by Continuum.jl itself ("Continuum orbital for local potential does not allow
        'exchange'"). Re-checking it here would only pretend to a plumbing that does not exist. Lifting that
        restriction needs a non-local (Hartree-Fock) continuum solver and is out of scope for every approach
        above.
"""
function validateConditions(conditions::Cascade.CascadeApproachConditions)
    ## GroupedBlockCI needs cascade blocks that may hold several configurations; the block generation still
    ## puts exactly one configuration into every block, so the grouping has nothing to act upon yet.
    if  conditions.levels == GroupedBlockCI()
        error("Cascade.UserMCA() is not yet available: it needs cascade blocks that comprise several " *
              "configurations, while the block generation currently places exactly one configuration into " *
              "each block. Use Cascade.RefinedSCA() or Basics.SCA() instead.")
    end

    return( nothing )
end
