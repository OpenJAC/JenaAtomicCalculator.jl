
## Cascade approaches: the CONDITIONS that define them.
##
## A cascade approach fixes FOUR conditions and nothing else. Which processes are included, which multipoles
## are applied, whether shake configurations are added, whether the Breit interaction or QED enter the
## Hamiltonian -- all of that stays with the cascade SCHEME and with AsfSettings, where it already lives.
## Keeping those apart is deliberate: the historical ShakeMCA approach duplicated the shakeFromShells/
## shakeToShells fields that Cascade.StepwiseDecayScheme already carried, and "only E1 transitions" was
## written into an approach while being a property of the requested processes; in both cases two places could
## disagree about the same question with nothing to arbitrate.
##
## Each approach differs from the one below it in exactly TWO of the conditions, and every approach is
## reported to the user as the same six lines in the same order, stated absolutely (never as "as SCA, but
## ..."), so that reading one block tells you what you got without holding another approach in your head:
##
##      (1) level representation      <- from .levels
##      (2) configurations per block  <- from .levels
##      (3) bound orbitals            <- from .bound
##      (4) continuum orbitals        <- from .continuum
##      (5) continuum potential       <- from .continuumPotential
##      (6) e-e interaction           <- from .levels plus the caller's AsfSettings
##
## Cascade.displayApproach generates all six lines FROM the conditions. That is the point: the previous
## free-text blocks in each -inc- file kept announcing "single-CSF levels and without any configuration
## mixing" for months while the code underneath was doing full configuration interaction.


"""
`abstract type Cascade.AbstractLevelRepresentation`
    ... specifies how the atomic levels of a cascade block are represented.

    + struct SingleCSF        ... every level is a single CSF; no configuration mixing of any kind.
    + struct BlockCI          ... configuration mixing WITHIN each block (intermediate coupling); one
                                  configuration per block, so no mixing between configurations.
    + struct GroupedBlockCI   ... as BlockCI, but a block may comprise SEVERAL configurations as grouped by
                                  the user, which admits correlation between those configurations.
"""
abstract type  AbstractLevelRepresentation                                    end
struct         SingleCSF        <:  Cascade.AbstractLevelRepresentation       end
struct         BlockCI          <:  Cascade.AbstractLevelRepresentation       end
struct         GroupedBlockCI   <:  Cascade.AbstractLevelRepresentation       end


"""
`abstract type Cascade.AbstractBoundOrbitals`
    ... specifies which set of bound one-electron orbitals the cascade blocks are built from.

    + struct GlobalOrbitals       ... ONE set for the whole cascade, from the mean field of the initial ion.
                                      The decayed states are then described with orbitals that never relax
                                      after the hole has been filled.
    + struct ChargeStateOrbitals  ... one self-consistent field per charge state, shared by all blocks with
                                      the same number of electrons.
    + struct MultipletOrbitals    ... one self-consistent field per multiplet, i.e. per block; every block's
                                      levels relax on their own orbitals.
"""
abstract type  AbstractBoundOrbitals                                          end
struct         GlobalOrbitals       <:  Cascade.AbstractBoundOrbitals         end
struct         ChargeStateOrbitals  <:  Cascade.AbstractBoundOrbitals         end
struct         MultipletOrbitals    <:  Cascade.AbstractBoundOrbitals         end


"""
`abstract type Cascade.AbstractContinuumOrbitals`
    ... specifies at which energy the orbital of the emitted (Auger) electron is generated.

    + struct StepAveragedContinuum       ... ONE set per cascade step, at that step's mean transition energy,
                                             shared by all fine-structure transitions of the step.
    + struct TransitionResolvedContinuum ... one orbital per fine-structure transition, at that transition's
                                             own energy.
"""
abstract type  AbstractContinuumOrbitals                                              end
struct         StepAveragedContinuum        <:  Cascade.AbstractContinuumOrbitals     end
struct         TransitionResolvedContinuum  <:  Cascade.AbstractContinuumOrbitals     end


"""
`abstract type Cascade.AbstractContinuumPotential`
    ... specifies the potential in which the continuum orbital is generated.

    + struct LocalPotential   ... a local (Dirac-Fock-Slater) potential; no exchange with the bound electrons.
    + struct WithExchange     ... exchange between the free and the bound electrons included. NOT available:
                                  JAC generates continuum orbitals in a local potential only, and lifting that
                                  needs a non-local (Hartree-Fock) continuum solver. Selecting it is reported
                                  as a warning rather than an error, so that an approach which merely ASKS for
                                  exchange stays usable for everything else it offers.
"""
abstract type  AbstractContinuumPotential                                     end
struct         LocalPotential  <:  Cascade.AbstractContinuumPotential         end
struct         WithExchange    <:  Cascade.AbstractContinuumPotential         end


"""
`struct  Cascade.RefinedSCA  <:  Basics.AbstractCascadeApproach`
    ... a cascade approach that refines SCA by relaxing the bound orbitals per multiplet and by resolving the
        continuum orbitals per fine-structure transition; see Cascade.conditions.

        Declared here rather than in Basics: only the abstract supertype Basics.AbstractCascadeApproach has to
        live in Basics, because DecayYield.Settings needs it as a compile-time field type and
        module-DecayYield.jl is included before module-Cascade.jl. Subtypes may be declared where they belong.
"""
struct   RefinedSCA  <:  Basics.AbstractCascadeApproach   end


"""
`struct  Cascade.CascadeApproachConditions`
    ... the four conditions that together define a cascade approach.

    + levels             ::Cascade.AbstractLevelRepresentation  ... representation of the levels of a block.
    + bound              ::Cascade.AbstractBoundOrbitals        ... origin of the bound one-electron orbitals.
    + continuum          ::Cascade.AbstractContinuumOrbitals    ... energy at which free-electron orbitals are
                                                                    generated.
    + continuumPotential ::Cascade.AbstractContinuumPotential   ... potential in which they are generated.
"""
struct  CascadeApproachConditions
    levels              ::Cascade.AbstractLevelRepresentation
    bound               ::Cascade.AbstractBoundOrbitals
    continuum           ::Cascade.AbstractContinuumOrbitals
    continuumPotential  ::Cascade.AbstractContinuumPotential
end


"""
`Cascade.conditions(approach::Basics.AbstractCascadeApproach)`
    ... returns the Cascade.CascadeApproachConditions that define the given cascade approach. Each approach
        changes exactly TWO conditions with respect to the one below it:

        AverageSCA   SingleCSF      / GlobalOrbitals      / StepAveraged       / LocalPotential
        SCA          BlockCI        / ChargeStateOrbitals / StepAveraged       / LocalPotential
        RefinedSCA   BlockCI        / MultipletOrbitals   / TransitionResolved / LocalPotential
        UserMCA      GroupedBlockCI / MultipletOrbitals   / TransitionResolved / WithExchange

        An  conditions::Cascade.CascadeApproachConditions  is returned.
"""
function conditions(approach::Basics.AbstractCascadeApproach)
    if      approach == Basics.AverageSCA()
        return( CascadeApproachConditions( SingleCSF(),      GlobalOrbitals(),      StepAveragedContinuum(),       LocalPotential() ) )
    elseif  approach == Basics.SCA()
        return( CascadeApproachConditions( BlockCI(),        ChargeStateOrbitals(), StepAveragedContinuum(),       LocalPotential() ) )
    elseif  approach == Cascade.RefinedSCA()
        return( CascadeApproachConditions( BlockCI(),        MultipletOrbitals(),   TransitionResolvedContinuum(), LocalPotential() ) )
    elseif  approach == Basics.UserMCA()
        return( CascadeApproachConditions( GroupedBlockCI(), MultipletOrbitals(),   TransitionResolvedContinuum(), WithExchange()   ) )
    else
        error("Unsupported cascade approach: $approach.")
    end
end


"""
`Cascade.displayApproach(stream::IO, approach::Basics.AbstractCascadeApproach, settings::AsfSettings)`
    ... writes the six lines that state, absolutely and always in the same order, what the given approach
        actually does. Nothing is returned.
"""
function displayApproach(stream::IO, approach::Basics.AbstractCascadeApproach, settings::AsfSettings)
    cd = Cascade.conditions(approach)

    if      cd.levels == SingleCSF()        sLevels = "single CSF; no configuration mixing"
    else                                    sLevels = "CI within each block (intermediate coupling)"       end
    if      cd.levels == GroupedBlockCI()   sGroup  = "several, as grouped by the user"
    else                                    sGroup  = "one"                                                end
    if      cd.bound == GlobalOrbitals()    sBound  = "one set for the whole cascade, from the initial ion"
    elseif  cd.bound == ChargeStateOrbitals()  sBound = "one self-consistent field per charge state"
    else                                    sBound  = "one self-consistent field per multiplet"            end
    if      cd.continuum == StepAveragedContinuum()
                                            sCont   = "one set per cascade step, at the step's mean energy"
    else                                    sCont   = "one per fine-structure transition, at its own energy" end
    if      cd.continuumPotential == LocalPotential()
                                            sPot    = "local (DFS); no exchange with the bound electrons"
    else                                    sPot    = "exchange requested -- NOT available, see the warning"  end
    if      cd.levels == SingleCSF()        sEe     = "Coulomb only (DiagonalCoulomb); Breit cannot be selected here"
    else                                    sEe     = "$(typeof(settings.eeInteractionCI)), as given by AsfSettings"  end

    println(stream, "\n* Cascade approach:  $approach")
    println(stream, "    (1) level representation ..... $sLevels")
    println(stream, "    (2) configurations per block . $sGroup")
    println(stream, "    (3) bound orbitals ........... $sBound")
    println(stream, "    (4) continuum orbitals ....... $sCont")
    println(stream, "    (5) continuum potential ...... $sPot")
    println(stream, "    (6) e-e interaction in H ..... $sEe\n")

    return( nothing )
end


"""
`Cascade.asfSettingsForApproach(approach::Basics.AbstractCascadeApproach, settings::AsfSettings)`
    ... returns the AsfSettings with which the cascade blocks must be represented for the given approach. The
        approach decides ONLY whether the Hamiltonian matrix is taken diagonal or is fully diagonalized; the
        e-e interaction itself (Coulomb, Gaunt, Breit) stays the caller's choice through AsfSettings.

        JAC packs both questions into the single field eeInteractionCI, whose values are DiagonalCoulomb,
        CoulombInteraction, CoulombGaunt, BreitInteraction and CoulombBreit -- there is no "DiagonalBreit".
        Consequently SingleCSF can only be realized as DiagonalCoulomb, and a Breit request cannot be honoured
        together with AverageSCA; that is a limitation of the type, not a choice made here. From SCA upwards
        the caller's setting is passed through untouched, and is only promoted when it says DiagonalCoulomb,
        which would contradict the configuration mixing the approach asks for. An  AsfSettings  is returned.
"""
function asfSettingsForApproach(approach::Basics.AbstractCascadeApproach, settings::AsfSettings)
    cd = Cascade.conditions(approach)
    if      cd.levels == SingleCSF()
        return( AsfSettings(settings; eeInteractionCI=DiagonalCoulomb()) )
    elseif  settings.eeInteractionCI == DiagonalCoulomb()
        return( AsfSettings(settings; eeInteractionCI=CoulombInteraction()) )
    else
        return( settings )
    end
end


"""
`Cascade.generateBoundOrbitals(approach::Basics.AbstractCascadeApproach, comp::Cascade.Computation,
                               confs::Array{Configuration,1}; printout::Bool=true)`
    ... generates the bound one-electron orbitals prescribed by the approach, as a dictionary that maps the
        number of electrons onto the orbital set to be used for every configuration with that many electrons.
        An empty dictionary is returned for MultipletOrbitals, where the orbitals cannot be prepared in
        advance because each block generates its own; the caller then runs the self-consistent field itself.
        A  Dict{Int64, Dict{Subshell, Orbital}}  is returned.

        For GlobalOrbitals the single set is generated from a mean-field basis spanning ALL of the cascade's
        configurations at once. Seeding it from the initial configurations alone -- the literal wording of the
        published AverageSCA -- is not sufficient in general: an excitation or ionization cascade reaches
        configurations with subshells the initial one never had, and the frozen-orbital CI would then fail on a
        missing subshell. Spanning all configurations keeps one common set while guaranteeing that every
        subshell the cascade can reach is present.
"""
function generateBoundOrbitals(approach::Basics.AbstractCascadeApproach, comp::Cascade.Computation,
                               confs::Array{Configuration,1}; printout::Bool=true)
    cd         = Cascade.conditions(approach)
    orbitalSets = Dict{Int64, Dict{Subshell, Orbital}}()
    if      cd.bound == MultipletOrbitals()
        return( orbitalSets )
    elseif  cd.bound == GlobalOrbitals()
        allConfigs = unique( append!(Configuration[], confs, comp.initialConfigs) )
        orbitals   = Cascade.generateMeanFieldOrbitals(comp, allConfigs)
        for  conf in allConfigs    orbitalSets[conf.NoElectrons] = orbitals    end
    else
        allConfigs = unique( append!(Configuration[], confs, comp.initialConfigs) )
        for  NoElectrons  in  unique( [conf.NoElectrons for conf in allConfigs] )
            sameConfigs                = filter(c -> c.NoElectrons == NoElectrons, allConfigs)
            orbitalSets[NoElectrons]   = Cascade.generateMeanFieldOrbitals(comp, sameConfigs)
        end
    end

    return( orbitalSets )
end


"""
`Cascade.generateMeanFieldOrbitals(comp::Cascade.Computation, confs::Array{Configuration,1})`
    ... generates one common set of bound orbitals from a mean-field basis spanning the given configurations.
        A  Dict{Subshell, Orbital}  is returned.
"""
function generateMeanFieldOrbitals(comp::Cascade.Computation, confs::Array{Configuration,1})
    repBasis  = AtomicState.Representation("Cascade mean-field basis", comp.nuclearModel, comp.grid, confs,
                                           AtomicState.MeanFieldBasis( AtomicState.MeanFieldSettings() ))
    repOutput = Basics.generate(repBasis, output=true)

    return( repOutput["mean-field basis"].orbitals )
end


"""
`Cascade.validateConditions(approach::Basics.AbstractCascadeApproach)`
    ... checks, at cascade set-up time, that the requested approach can actually be realized. Conditions that
        are not implemented raise an error here, before any orbital or amplitude work is done, so that a
        cascade which cannot succeed fails in the first second rather than in the middle of some step. The one
        exception is the continuum exchange, which only warns: an approach that merely ASKS for exchange stays
        usable for everything else it offers. Nothing is returned.
"""
function validateConditions(approach::Basics.AbstractCascadeApproach)
    cd = Cascade.conditions(approach)

    if  cd.levels == GroupedBlockCI()
        error("$approach is not yet available: it needs cascade blocks that comprise several configurations, " *
              "while the block generation currently places exactly one configuration into each block.")
    end
    ## No check is needed for TransitionResolvedContinuum: it IS implemented, in the Auger branch of
    ## Cascade.computeSteps that calls AutoIonization.computeLinesCascade and generates one continuum orbital
    ## per fine-structure transition. It had simply never been reachable on purpose -- the dispatch keyed on
    ## the approach NAME, so every approach other than AverageSCA fell into it by accident rather than by
    ## request. Both values of condition (4) are now selected by the condition itself.
    if  cd.continuumPotential == WithExchange()
        @warn("$approach asks for exchange between the free and the bound electrons, which JAC cannot yet " *
              "provide: continuum orbitals are generated in a local potential only, and lifting that needs a " *
              "non-local (Hartree-Fock) continuum solver. The computation proceeds WITHOUT exchange; every " *
              "other condition of this approach is honoured.")
    end

    return( nothing )
end
