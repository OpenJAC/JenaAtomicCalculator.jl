
"""
`module  JAC.PhotonScattering`
... a submodel of JAC for the SCATTERING OF PHOTONS at atoms and ions, i.e. for those processes that have a photon on BOTH sides of
    the reaction, or that convert an incoming photon into particles -- as distinct from the absorption and emission processes of
    PhotoIonization, PhotoExcitation, PhotoEmission and PhotoRecombination, which have a photon on one side only.

    The module is organized along three orthogonal axes, `beamType` x `process` x `approximation`, the same arrangement that
    ParticleScattering uses for massive projectiles; see module-PhotonScattering-inc-structs.jl for what each axis carries and why
    the third one is worth having.

    IMPLEMENTED TODAY: bound-free pair creation, gamma + |i(N)> --> |f(N+1)> + e^+, in
    module-PhotonScattering-inc-pair-creation.jl; and Rayleigh/Raman scattering in
    module-PhotonScattering-inc-rayleigh.jl.

    ON RAYLEIGH: it is a FRESH implementation, not the older module JAC.RayleighCompton moved here. That module was found on
    21-Aug-2026 not to compute -- its sum over intermediate states runs for the first level only (`if ig != 1 continue end` at
    module-RayleighCompton.jl:338), so a second-order amplitude is reduced to its first term, and the pole branch throws a
    MethodError besides. The maintainer's decision was that getting it right is a new implementation rather than a repair; see
    examples/example-Dg.jl branch a for the full diagnosis. This file follows MultiPhotonTransition's two-photon amplitudes
    instead, which perform the same kind of sum, are validated across eight dated example branches, and skip resonant
    denominators by tolerance rather than integrating over a pole.

    ON COMPTON: `ComptonScattering()` presently means Raman-type inelastic scattering into DISCRETE final levels, which is what
    a final-state Multiplet can express. The CONTINUUM Compton profile -- an ejected electron and the doubly differential
    d^2 sigma / dOmega dOmega_out -- is NOT implemented in JAC under any name, and no current `Line` can express it, since a
    line carries one outgoing energy fixed by energy conservation. The token is kept so that the eventual implementation has
    its place named; what is not kept is the silent over-promise.

    PLANNED: resonant inelastic scattering (RIXS), which exists today as the separate module ResonantInelastic and whose working
    state has not been established, it carrying no example file at all; and scattering of TWISTED (Bessel) photon beams, for
    which Beam.BesselBeam already provides the beam side. Neither is done.
"""
module PhotonScattering

using  Printf,
        ..AngularMomentum, ..AtomicState, ..Basics, ..Beam, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear,
        ..ParticleScattering, ..PhotoEmission, ..Radial, ..TableStrings


include("module-PhotonScattering-inc-structs.jl")


"""
`struct  PhotonScattering.Settings  <:  AbstractProcessSettings`
    ... defines a type for the settings of a photon-scattering computation.

    + process             ::PhotonScattering.AbstractPhotonProcess          ... what happens to the target.
    + approximation       ::PhotonScattering.AbstractScatteringApproximation ... how the amplitude is evaluated.
    + beamType            ::Beam.AbstractBeamType                           ... how the incoming photon is prepared.
    + photonEnergies      ::Array{Float64,1}    ... energies of the incoming photon [a.u.].
    + multipoles          ::Array{EmMultipole,1}
        ... multipoles of the radiation field to be included. For pair creation the incoming photon carries omega > 2 m c^2, i.e.
            q a_0 > 274, so the series is short only where the captured electron sits close to the nucleus; cf. the note on q <r> in
            module-ParticleScattering-inc-annihilation.jl, which applies here unchanged.
    + gauges              ::Array{UseGauge,1}   ... gauges in which the amplitudes are evaluated; two that disagree are the standard
                                                    internal check on a photon amplitude, so both are the default.
    + polarThetas         ::Array{Float64,1}    ... polar angles at which the differential observables are wanted [rad].
    + polarPhis           ::Array{Float64,1}    ... azimuthal angles [rad].
    + incidentStokes      ::ExpStokes           ... Stokes parameters of the incident radiation.
    + maxKappa            ::Int64               ... largest |kappa| of the continuum particle to be included.
    + gMultiplet          ::Union{Multiplet, Array{AtomicState.GreenChannel,1}}
        ... the intermediate levels over which a SECOND-order amplitude is summed, needed by every SecondOrderGreen process and
            ignored by a first-order one. Either a plain Multiplet -- typically a short, explicitly chosen list of 2-4 levels,
            which is what makes the sum inspectable -- or a Green expansion, exactly as MultiPhotonTransition accepts both.
    + selfTolerance       ::Float64
        ... a resonant intermediate level, where the energy denominator vanishes, is SKIPPED once |denominator| falls below this
            value. That is the boundary between non-resonant and resonant scattering, and it is guarded rather than hidden: the
            perturbative expression is simply not defined there.
    + printBefore         ::Bool                ... True, if all lines are printed before their evaluation.
    + lineSelection       ::LineSelection       ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractProcessSettings
    process               ::PhotonScattering.AbstractPhotonProcess
    approximation         ::PhotonScattering.AbstractScatteringApproximation
    beamType              ::Beam.AbstractBeamType
    photonEnergies        ::Array{Float64,1}
    multipoles            ::Array{EmMultipole,1}
    gauges                ::Array{UseGauge,1}
    polarThetas           ::Array{Float64,1}
    polarPhis             ::Array{Float64,1}
    incidentStokes        ::ExpStokes
    maxKappa              ::Int64
    gMultiplet            ::Union{Multiplet, Array{AtomicState.GreenChannel,1}}
    selfTolerance         ::Float64
    printBefore           ::Bool
    lineSelection         ::LineSelection
end


"""
`PhotonScattering.Settings()`  ... constructor for the default PhotonScattering.Settings.
"""
function Settings()
    Settings( PhotonScattering.BoundFreePairCreation(), PhotonScattering.FirstOrderVertex(), Beam.PlaneWave(),
              Float64[], EmMultipole[E1], UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Float64[],
              Basics.ExpStokes(), 4, Multiplet(), 1.0e-8, false, LineSelection() )
end


"""
`PhotonScattering.Settings(set::PhotonScattering.Settings;`

        process=..,             approximation=..,           beamType=..,            photonEnergies=..,
        multipoles=..,          gauges=..,                  polarThetas=..,         polarPhis=..,
        incidentStokes=..,      maxKappa=..,                gMultiplet=..,          selfTolerance=..,
        printBefore=..,         lineSelection=.. )

    ... constructor for modifying the given PhotonScattering.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotonScattering.Settings;
    process::Union{Nothing,PhotonScattering.AbstractPhotonProcess}=nothing,
    approximation::Union{Nothing,PhotonScattering.AbstractScatteringApproximation}=nothing,
    beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,      photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,     gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    polarThetas::Union{Nothing,Array{Float64,1}}=nothing,        polarPhis::Union{Nothing,Array{Float64,1}}=nothing,
    incidentStokes::Union{Nothing,ExpStokes}=nothing,            maxKappa::Union{Nothing,Int64}=nothing,
    gMultiplet::Union{Nothing,Multiplet,Array{AtomicState.GreenChannel,1}}=nothing,
    selfTolerance::Union{Nothing,Float64}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                    lineSelection::Union{Nothing,LineSelection}=nothing)

    if  isnothing(process)          processx        = set.process         else  processx        = process         end
    if  isnothing(approximation)    approximationx  = set.approximation   else  approximationx  = approximation   end
    if  isnothing(beamType)         beamTypex       = set.beamType        else  beamTypex       = beamType        end
    if  isnothing(photonEnergies)   photonEnergiesx = set.photonEnergies  else  photonEnergiesx = photonEnergies  end
    if  isnothing(multipoles)       multipolesx     = set.multipoles      else  multipolesx     = multipoles      end
    if  isnothing(gauges)           gaugesx         = set.gauges          else  gaugesx         = gauges          end
    if  isnothing(polarThetas)      polarThetasx    = set.polarThetas     else  polarThetasx    = polarThetas     end
    if  isnothing(polarPhis)        polarPhisx      = set.polarPhis       else  polarPhisx      = polarPhis       end
    if  isnothing(incidentStokes)   incidentStokesx = set.incidentStokes  else  incidentStokesx = incidentStokes  end
    if  isnothing(maxKappa)         maxKappax       = set.maxKappa        else  maxKappax       = maxKappa        end
    if  isnothing(gMultiplet)       gMultipletx     = set.gMultiplet      else  gMultipletx     = gMultiplet      end
    if  isnothing(selfTolerance)    selfTolerancex  = set.selfTolerance   else  selfTolerancex  = selfTolerance   end
    if  isnothing(printBefore)      printBeforex    = set.printBefore     else  printBeforex    = printBefore     end
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection   else  lineSelectionx  = lineSelection   end

    Settings( processx, approximationx, beamTypex, photonEnergiesx, multipolesx, gaugesx, polarThetasx, polarPhisx,
              incidentStokesx, maxKappax, gMultipletx, selfTolerancex, printBeforex, lineSelectionx )
end


# `Base.show(io::IO, settings::PhotonScattering.Settings)`  ... prepares a proper printout of settings::PhotonScattering.Settings.
function Base.show(io::IO, settings::PhotonScattering.Settings)
    println(io, "process:               $(settings.process)  ")
    println(io, "approximation:         $(settings.approximation)  ")
    println(io, "beamType:              $(settings.beamType)  ")
    println(io, "photonEnergies:        $(settings.photonEnergies)  ")
    println(io, "multipoles:            $(settings.multipoles)  ")
    println(io, "gauges:                $(settings.gauges)  ")
    println(io, "polarThetas:           $(settings.polarThetas)  ")
    println(io, "polarPhis:             $(settings.polarPhis)  ")
    println(io, "incidentStokes:        $(settings.incidentStokes)  ")
    println(io, "maxKappa:              $(settings.maxKappa)  ")
    println(io, "gMultiplet:            $(typeof(settings.gMultiplet))  ")
    println(io, "selfTolerance:         $(settings.selfTolerance)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "lineSelection:         $(settings.lineSelection)  ")
end


"""
`PhotonScattering.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                               settings::PhotonScattering.Settings; output=true)`
    ... the single entry point of the module, reached from Basics.perform(::Atomic.Computation). It dispatches on settings.process,
        every process having its own pipeline in its own -inc- file, and a process that is named but not yet implemented saying so
        rather than failing obscurely. An Array{PhotonScattering.Line,1} is returned if output=true, and nothing otherwise.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                      settings::PhotonScattering.Settings; output=true)
    if      settings.process == PhotonScattering.BoundFreePairCreation()
        return( PhotonScattering.computePairCreationLines(finalMultiplet, initialMultiplet, nm, grid, settings; output=output) )
    elseif  settings.process in [PhotonScattering.RayleighScattering(), PhotonScattering.ComptonScattering(),
                                 PhotonScattering.RamanScattering()]
        return( PhotonScattering.computeRayleighLines(finalMultiplet, initialMultiplet, nm, grid, settings; output=output) )
    elseif  settings.process == PhotonScattering.ResonantScattering()
        error("PhotonScattering: resonant inelastic scattering is not yet part of this module -- it is still the separate " *
              "module ResonantInelastic, whose working state has not been established (it carries no example file).")
    else
        error("PhotonScattering: no pipeline for process $(settings.process).")
    end
end


"""
`PhotonScattering.intermediateLevels(gMultiplet::Multiplet, Jsym::LevelSymmetry)`
    ... to select those intermediate levels of the given multiplet that carry the symmetry Jsym, i.e. the levels over which a
        second-order amplitude of that intermediate symmetry is summed. An Array{Level,1} is returned, empty if the multiplet
        does not span Jsym at all -- which the caller must treat as a reason to warn rather than as a zero amplitude.

        This mirrors MultiPhotonTransition.intermediateLevels, deliberately and by the same name, but is implemented locally:
        MultiPhotonTransition is included AFTER this module in JenaAtomicCalculator.jl and is therefore not visible here, and a
        module-level dependency of photon scattering on multiphoton transitions would be a heavy price for a filter.
"""
function intermediateLevels(gMultiplet::Multiplet, Jsym::LevelSymmetry)
    return( filter(lev -> LevelSymmetry(lev.J, lev.parity) == Jsym, gMultiplet.levels) )
end


"""
`PhotonScattering.intermediateLevels(gChannels::Array{AtomicState.GreenChannel,1}, Jsym::LevelSymmetry)`
    ... to select those intermediate levels of the given Green-function channels that carry the symmetry Jsym; the levels of every
        channel whose own symmetry matches are collected. An Array{Level,1} is returned, empty if no channel carries Jsym.
"""
function intermediateLevels(gChannels::Array{AtomicState.GreenChannel,1}, Jsym::LevelSymmetry)
    levels = Level[]
    for  ch in gChannels
        if  Jsym == ch.symmetry    append!(levels, ch.gMultiplet.levels)    end
    end

    return( levels )
end


include("module-PhotonScattering-inc-pair-creation.jl")
include("module-PhotonScattering-inc-rayleigh.jl")

end # module
