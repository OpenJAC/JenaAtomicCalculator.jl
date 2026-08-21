
"""
`module  JAC.PhotonScattering`
... a submodel of JAC for the SCATTERING OF PHOTONS at atoms and ions, i.e. for those processes that have a photon on BOTH sides of
    the reaction, or that convert an incoming photon into particles -- as distinct from the absorption and emission processes of
    PhotoIonization, PhotoExcitation, PhotoEmission and PhotoRecombination, which have a photon on one side only.

    The module is organized along three orthogonal axes, `beamType` x `process` x `approximation`, the same arrangement that
    ParticleScattering uses for massive projectiles; see module-PhotonScattering-inc-structs.jl for what each axis carries and why
    the third one is worth having.

    IMPLEMENTED TODAY: bound-free pair creation, gamma + |i(N)> --> |f(N+1)> + e^+, in
    module-PhotonScattering-inc-pair-creation.jl.

    PLANNED, and named here so that the intended scope is visible rather than inferred: Rayleigh, Compton and Raman scattering,
    which exist today as the separate module RayleighCompton and are to be absorbed as -inc-rayleigh-compton.jl; resonant inelastic
    scattering (RIXS), which exists as the separate module ResonantInelastic and is to be absorbed as -inc-resonant.jl once its
    working state has been established, it carrying no example file at all today; and scattering of TWISTED (Bessel) photon beams,
    for which Beam.BesselBeam already provides the beam side. None of that is done, and no part of this module should be read as
    claiming it is.
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
    printBefore           ::Bool
    lineSelection         ::LineSelection
end


"""
`PhotonScattering.Settings()`  ... constructor for the default PhotonScattering.Settings.
"""
function Settings()
    Settings( PhotonScattering.BoundFreePairCreation(), PhotonScattering.FirstOrderVertex(), Beam.PlaneWave(),
              Float64[], EmMultipole[E1], UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Float64[],
              Basics.ExpStokes(), 4, false, LineSelection() )
end


"""
`PhotonScattering.Settings(set::PhotonScattering.Settings;`

        process=..,             approximation=..,           beamType=..,            photonEnergies=..,
        multipoles=..,          gauges=..,                  polarThetas=..,         polarPhis=..,
        incidentStokes=..,      maxKappa=..,                printBefore=..,         lineSelection=.. )

    ... constructor for modifying the given PhotonScattering.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotonScattering.Settings;
    process::Union{Nothing,PhotonScattering.AbstractPhotonProcess}=nothing,
    approximation::Union{Nothing,PhotonScattering.AbstractScatteringApproximation}=nothing,
    beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,      photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,     gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    polarThetas::Union{Nothing,Array{Float64,1}}=nothing,        polarPhis::Union{Nothing,Array{Float64,1}}=nothing,
    incidentStokes::Union{Nothing,ExpStokes}=nothing,            maxKappa::Union{Nothing,Int64}=nothing,
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
    if  isnothing(printBefore)      printBeforex    = set.printBefore     else  printBeforex    = printBefore     end
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection   else  lineSelectionx  = lineSelection   end

    Settings( processx, approximationx, beamTypex, photonEnergiesx, multipolesx, gaugesx, polarThetasx, polarPhisx,
              incidentStokesx, maxKappax, printBeforex, lineSelectionx )
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
        error("PhotonScattering: Rayleigh, Compton and Raman scattering are not yet part of this module -- they are still " *
              "the separate module RayleighCompton, which is to be absorbed here as -inc-rayleigh-compton.jl. Use " *
              "RayleighCompton.Settings for now.")
    elseif  settings.process == PhotonScattering.ResonantScattering()
        error("PhotonScattering: resonant inelastic scattering is not yet part of this module -- it is still the separate " *
              "module ResonantInelastic, whose working state has not been established (it carries no example file).")
    else
        error("PhotonScattering: no pipeline for process $(settings.process).")
    end
end


include("module-PhotonScattering-inc-pair-creation.jl")

end # module
