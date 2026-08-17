
"""
`module  JAC.ParticleScattering`
... a submodel of JAC that contains all methods for computing scattering amplitudes and cross sections for the scattering of
    MASSIVE particles -- electrons, positrons, protons, ... -- at atoms and ions. The module is organised along three
    orthogonal axes, so that none of them has to be widened when another one is:

    + the PROJECTILE       ... Electron, Positron, Proton;
    + the INTERACTION MODEL ... how the target is represented, StaticField, StaticFieldExchange, ...;
    + the BEAM             ... how the projectile is prepared, Beam.PlaneWave, Beam.BesselBeam, ...

    The scattering amplitudes are computed as Dirac partial waves for BOTH spin-orbit partners of each orbital angular
    momentum, kappa = -l-1 and kappa = +l, and are kept as reduced amplitudes in a total-angular-momentum coupled basis.
    The familiar direct and spin-flip amplitudes f and g are a projection of these, valid for elastic scattering from a
    spinless target, and are provided as such.

    Photon scattering is NOT part of this module: it is a second-order matrix element over intermediate atomic states
    (Kramers-Heisenberg) and shares no computational core with the potential scattering treated here. It is computed by
    JAC.RayleighCompton. What the two do share -- the preparation of a twisted beam -- lives in JAC.Beam.
"""
module ParticleScattering


using  Printf, GSL, QuadGK, SpecialFunctions,
        ..AngularMomentum, ..Basics, ..Beam, ..Bsplines, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron,
        ..Nuclear, ..Radial, ..RadialIntegrals, ..SpinAngular, ..TableStrings


# The data structures come first: Settings below refers to them in its field types, and a struct cannot name a type
# that does not exist yet.
include("module-ParticleScattering-inc-structs.jl")


"""
`struct  ParticleScattering.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of computing scattering amplitudes and cross sections.

    + projectile          ::ParticleScattering.AbstractProjectile         ... the scattered particle.
    + process             ::ParticleScattering.AbstractScatteringProcess  ... elastic or inelastic scattering.
    + interaction         ::ParticleScattering.AbstractInteractionModel   ... how the target is represented.
    + beamType            ::Beam.AbstractBeamType                         ... how the projectile is prepared.
    + polarization        ::Basics.AbstractPolarization                   ... polarization of the beam.
    + impactEnergies      ::Array{Float64,1}                              ... impact energies of the projectile.
    + polarThetas         ::Array{Float64,1}                              ... polar scattering angles [rad].
    + polarPhis           ::Array{Float64,1}                              ... azimuthal scattering angles [rad].
    + printBefore         ::Bool               ... True, if all events are printed before their evaluation.
    + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
    + epsPartialWave      ::Float64
        ... convergence criterion for the partial-wave series; it is ended once the contribution of one l to the elastic
            cross section falls below this fraction of the running total, twice in succession.
    + maxL                ::Int64
        ... hard upper bound on the orbital angular momentum, a backstop rather than the normal way of ending the series;
            a computation that reaches it says so.
"""
struct Settings  <:  AbstractProcessSettings
    projectile            ::ParticleScattering.AbstractProjectile
    process               ::ParticleScattering.AbstractScatteringProcess
    interaction           ::ParticleScattering.AbstractInteractionModel
    beamType              ::Beam.AbstractBeamType
    polarization          ::Basics.AbstractPolarization
    impactEnergies        ::Array{Float64,1}
    polarThetas           ::Array{Float64,1}
    polarPhis             ::Array{Float64,1}
    printBefore           ::Bool
    lineSelection         ::LineSelection
    epsPartialWave        ::Float64
    maxL                  ::Int64
end


"""
`ParticleScattering.Settings()`  ... constructor for the default ParticleScattering.Settings.
"""
function Settings()
    Settings( ParticleScattering.Electron(), ParticleScattering.ElasticScattering(),
              ParticleScattering.StaticFieldExchange(), Beam.PlaneWave(), Basics.LinearPolarization(),
              Float64[], Float64[], Float64[], false, LineSelection(), 1.0e-6, 200 )
end


"""
`ParticleScattering.Settings(set::ParticleScattering.Settings;`

        projectile=..,          process=..,                 interaction=..,         beamType=..,
        polarization=..,        impactEnergies=..,          polarThetas=..,         polarPhis=..,
        printBefore=..,         lineSelection=..,           epsPartialWave=..,      maxL=.. )

    ... constructor for modifying the given ParticleScattering.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ParticleScattering.Settings;
    projectile::Union{Nothing,ParticleScattering.AbstractProjectile}=nothing,
    process::Union{Nothing,ParticleScattering.AbstractScatteringProcess}=nothing,
    interaction::Union{Nothing,ParticleScattering.AbstractInteractionModel}=nothing,
    beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,
    polarization::Union{Nothing,Basics.AbstractPolarization}=nothing,
    impactEnergies::Union{Nothing,Array{Float64,1}}=nothing,     polarThetas::Union{Nothing,Array{Float64,1}}=nothing,
    polarPhis::Union{Nothing,Array{Float64,1}}=nothing,          printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,         epsPartialWave::Union{Nothing,Float64}=nothing,
    maxL::Union{Nothing,Int64}=nothing)

    if  isnothing(projectile)       projectilex     = set.projectile      else  projectilex     = projectile      end
    if  isnothing(process)          processx        = set.process         else  processx        = process         end
    if  isnothing(interaction)      interactionx    = set.interaction     else  interactionx    = interaction     end
    if  isnothing(beamType)         beamTypex       = set.beamType        else  beamTypex       = beamType        end
    if  isnothing(polarization)     polarizationx   = set.polarization    else  polarizationx   = polarization    end
    if  isnothing(impactEnergies)   impactEnergiesx = set.impactEnergies  else  impactEnergiesx = impactEnergies  end
    if  isnothing(polarThetas)      polarThetasx    = set.polarThetas     else  polarThetasx    = polarThetas     end
    if  isnothing(polarPhis)        polarPhisx      = set.polarPhis       else  polarPhisx      = polarPhis       end
    if  isnothing(printBefore)      printBeforex    = set.printBefore     else  printBeforex    = printBefore     end
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection   else  lineSelectionx  = lineSelection   end
    if  isnothing(epsPartialWave)   epsPartialWavex = set.epsPartialWave  else  epsPartialWavex = epsPartialWave  end
    if  isnothing(maxL)             maxLx           = set.maxL            else  maxLx           = maxL            end

    Settings( projectilex, processx, interactionx, beamTypex, polarizationx, impactEnergiesx, polarThetasx, polarPhisx,
              printBeforex, lineSelectionx, epsPartialWavex, maxLx )
end


# `Base.show(io::IO, settings::ParticleScattering.Settings)`
#       ... prepares a proper printout of the variable settings::ParticleScattering.Settings.
function Base.show(io::IO, settings::ParticleScattering.Settings)
    println(io, "projectile:            $(settings.projectile)  ")
    println(io, "process:               $(settings.process)  ")
    println(io, "interaction:           $(settings.interaction)  ")
    println(io, "beamType:              $(settings.beamType)  ")
    println(io, "polarization:          $(settings.polarization)  ")
    println(io, "impactEnergies:        $(settings.impactEnergies)  ")
    println(io, "polarThetas:           $(settings.polarThetas)  ")
    println(io, "polarPhis:             $(settings.polarPhis)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "lineSelection:         $(settings.lineSelection)  ")
    println(io, "epsPartialWave:        $(settings.epsPartialWave)  ")
    println(io, "maxL:                  $(settings.maxL)  ")
end


"""
`ParticleScattering.determineEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                    settings::ParticleScattering.Settings)`
    ... to determine the list of scattering events that are to be computed, i.e. one event per selected pair of levels and
        per impact energy, with the amplitudes and observables still empty. For elastic scattering the initial and final
        level must be the same one. An Array{ParticleScattering.Event,1} is returned.
"""
function determineEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)
    events = ParticleScattering.Event[]

    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            if  typeof(settings.process) == ParticleScattering.ElasticScattering  &&  iLevel.index != fLevel.index
                continue
            end
            for  en  in  settings.impactEnergies
                en < 0.  &&  error("Negative impact energy $en in ParticleScattering.Settings.")
                # The impact energies are given in the CURRENTLY selected unit, as everywhere else in JAC, and are
                # converted here once; everything downstream works in atomic units.
                energy = Defaults.convertUnits("energy: to atomic", en)
                push!( events, ParticleScattering.Event(settings.projectile, settings.process, settings.interaction,
                                                        settings.beamType, iLevel, fLevel, energy,
                                                        ParticleScattering.PartialWave[], ParticleScattering.ScatteringChannel[],
                                                        ParticleScattering.AngularObservables[],
                                                        ParticleScattering.IntegratedObservables()) )
            end
        end
    end

    return( events )
end


"""
`ParticleScattering.computeAmplitudesProperties(event::ParticleScattering.Event, nm::Nuclear.Model, grid::Radial.Grid,
                                                nrContinuum::Int64, settings::ParticleScattering.Settings;
                                                nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)`
    ... to compute the partial waves, the reduced scattering channels and all observables of the given event. A new
        event::ParticleScattering.Event is returned, for which these are now evaluated.
"""
function computeAmplitudesProperties(event::ParticleScattering.Event, nm::Nuclear.Model, grid::Radial.Grid,
                                     nrContinuum::Int64, settings::ParticleScattering.Settings;
                                     nuclearPot::Union{Nothing,Radial.Potential}=nothing, printout::Bool=true)
    contSettings = Continuum.Settings(false, nrContinuum)
    pws = ParticleScattering.computePartialWaves(event.projectile, event.interaction, event.finalLevel, event.impactEnergy,
                                                 nm, grid, contSettings, settings; nuclearPot=nuclearPot, printout=printout)
    channels = ParticleScattering.scatteringChannels(pws, event.finalLevel)
    #
    angular = ParticleScattering.AngularObservables[]
    for  theta  in  settings.polarThetas,  phi  in  settings.polarPhis
        push!( angular, ParticleScattering.beamObservables(event.beamType, pws, event.impactEnergy, theta, phi) )
    end
    integrated = ParticleScattering.integratedObservables(pws, event.impactEnergy)

    return( ParticleScattering.Event(event.projectile, event.process, event.interaction, event.beamType,
                                     event.initialLevel, event.finalLevel, event.impactEnergy, pws, channels,
                                     angular, integrated) )
end


"""
`ParticleScattering.computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                  grid::Radial.Grid, settings::ParticleScattering.Settings; output=true)`
    ... to compute all selected scattering events, together with their partial waves, reduced channels and observables, and
        to display them. An Array{ParticleScattering.Event,1} is returned if output=true, and nothing otherwise.
"""
function  computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::ParticleScattering.Settings; output=true)
    println("")
    printstyled("ParticleScattering.computeEvents(): The computation of scattering amplitudes starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    events = ParticleScattering.determineEvents(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    ParticleScattering.displayEvents(stdout, events)    end
    # Determine the maximum energy and check the consistency of the grid
    maxEnergy = 0.;   for  event in events   maxEnergy = max(maxEnergy, event.impactEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # The nuclear potential depends only on (nm, grid) and is built once for all events
    nucPot    = Nuclear.nuclearPotential(nm, grid)
    newEvents = ParticleScattering.Event[]
    for  event in events
        push!( newEvents, ParticleScattering.computeAmplitudesProperties(event, nm, grid, nrContinuum, settings;
                                                                        nuclearPot=nucPot) )
    end
    # Print all results to screen
    ParticleScattering.displayPhaseShifts(stdout, newEvents)
    ParticleScattering.displayCrossSections(stdout, newEvents)
    ParticleScattering.displayIntegratedCrossSections(stdout, newEvents)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ParticleScattering.displayPhaseShifts(iostream, newEvents)
                       ParticleScattering.displayCrossSections(iostream, newEvents)
                       ParticleScattering.displayIntegratedCrossSections(iostream, newEvents)     end

    if    output    return( newEvents )
    else            return( nothing )
    end
end


include("module-ParticleScattering-inc-interaction.jl")
include("module-ParticleScattering-inc-electron-dirac.jl")
include("module-ParticleScattering-inc-beams.jl")
include("module-ParticleScattering-inc-observables.jl")
include("module-ParticleScattering-inc-display.jl")

end # module
