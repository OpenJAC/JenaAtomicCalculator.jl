
"""
`module  JAC.PhotonScattering`
... a submodel of JAC for the SCATTERING OF PHOTONS at atoms and ions, i.e. for those processes that have a photon on BOTH sides of
    the reaction, or that convert an incoming photon into particles -- as distinct from the absorption and emission processes of
    PhotoIonization, PhotoExcitation, PhotoEmission and PhotoRecombination, which have a photon on one side only.

    The module is organized along three orthogonal axes, `beamType` x `process` x `approximation`, the same arrangement that
    ParticleScattering uses for massive projectiles; see module-PhotonScattering-inc-structs.jl for what each axis carries and why
    the third one is worth having.

    IMPLEMENTED TODAY: bound-free pair creation, gamma + |i(N)> --> |f(N+1)> + e^+, in
    module-PhotonScattering-inc-pair-creation.jl; Rayleigh/Raman scattering in
    module-PhotonScattering-inc-rayleigh.jl; and resonant inelastic scattering (RIXS) in
    module-PhotonScattering-inc-resonant.jl. The last two share their channel construction
    (determineSecondOrderChannels) and differ in one physical respect: Rayleigh SKIPS a resonant intermediate level, RIXS keeps it
    and regularises the denominator with settings.width.

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

    ON RIXS: it is likewise a FRESH implementation and is LINE-shaped, not pathway-shaped. The older module ResonantInelastic
    carries a Pathway type and takes its intermediate levels from the intermediateMultiplet that Basics.perform builds from
    intermediateConfigs; here they come from settings.gMultiplet like every other second-order process, so the module keeps ONE
    result type and ONE entry point. PhotonScattering.Settings is therefore deliberately NOT added to the intermediateMultiplet
    list in module-BasicsAZ-inc-perform.jl.

    PLANNED: scattering of TWISTED (Bessel) photon beams, for which Beam.BesselBeam already provides the beam side. Not done.
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
    + width               ::Float64
        ... resonance width Gamma of the intermediate levels [a.u.], read ONLY by ResonantScattering(). It regularises the
            resonant denominator as E_i + omega_in - E_nu + i*Gamma/2, which is what lets RIXS be evaluated ON resonance where
            Rayleigh must skip. It is zero by default, and the resonant pipeline refuses to run with a zero width rather than
            dividing by something arbitrarily small.
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
    width                 ::Float64
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
              Basics.ExpStokes(), 4, Multiplet(), 0., 1.0e-8, false, LineSelection() )
end


"""
`PhotonScattering.Settings(set::PhotonScattering.Settings;`

        process=..,             approximation=..,           beamType=..,            photonEnergies=..,
        multipoles=..,          gauges=..,                  polarThetas=..,         polarPhis=..,
        incidentStokes=..,      maxKappa=..,                gMultiplet=..,          width=..,
        selfTolerance=..,
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
    width::Union{Nothing,Float64}=nothing,                       selfTolerance::Union{Nothing,Float64}=nothing,
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
    if  isnothing(width)            widthx          = set.width           else  widthx          = width           end
    if  isnothing(selfTolerance)    selfTolerancex  = set.selfTolerance   else  selfTolerancex  = selfTolerance   end
    if  isnothing(printBefore)      printBeforex    = set.printBefore     else  printBeforex    = printBefore     end
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection   else  lineSelectionx  = lineSelection   end

    Settings( processx, approximationx, beamTypex, photonEnergiesx, multipolesx, gaugesx, polarThetasx, polarPhisx,
              incidentStokesx, maxKappax, gMultipletx, widthx, selfTolerancex, printBeforex, lineSelectionx )
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
    println(io, "width:                 $(settings.width)  ")
    println(io, "selfTolerance:         $(settings.selfTolerance)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "lineSelection:         $(settings.lineSelection)  ")
end


"""
`PhotonScattering.checkCombination(settings::PhotonScattering.Settings)`
    ... to check, BEFORE anything is computed, that the requested (beamType, process, approximation) triple is one the module can
        actually evaluate. Nothing is returned if it can; otherwise this raises. The three axes are orthogonal by design, which
        means most of the grid is EMPTY, and a Settings object can be constructed for any point on it -- so without this check a
        request for a Bessel beam or a form-factor amplitude would be accepted and silently answered with plane-wave second-order
        numbers. That is the failure mode this exists to prevent, and it follows the fallback pattern of
        ParticleScattering.scatteringPotential.

        TWO KINDS OF REFUSAL, deliberately worded differently, because they tell the user to do different things. A combination
        that is NOT YET IMPLEMENTED is a feature to wait for. A combination that is MEANINGLESS is a request to rethink: bound-free
        pair creation has ONE vertex, so SecondOrderGreen is not a more expensive version of it but a category error, and a
        second-order process cannot be evaluated with FirstOrderVertex for the same reason in reverse.
"""
function checkCombination(settings::PhotonScattering.Settings)
    secondOrder = [PhotonScattering.RayleighScattering(), PhotonScattering.ComptonScattering(),
                   PhotonScattering.RamanScattering(),    PhotonScattering.ResonantScattering()]
    firstOrder  = [PhotonScattering.BoundFreePairCreation()]

    # (a) meaningless combinations -- the order of the process and the order of the approximation disagree
    if      settings.process in firstOrder   &&  !(settings.approximation == PhotonScattering.FirstOrderVertex())
        error("\n\nThe combination\n    process       = $(settings.process)\n" *
              "    approximation = $(settings.approximation)\n\n" *
              "does not describe anything. Bound-free pair creation has a SINGLE photon vertex, so FirstOrderVertex() is not an " *
              "approximation to it -- it is exact, and is the only setting that means anything. SecondOrderGreen() and " *
              "FormFactorApproximation() both presuppose two vertices. Use FirstOrderVertex().")
    elseif  settings.process in secondOrder  &&  settings.approximation == PhotonScattering.FirstOrderVertex()
        error("\n\nThe combination\n    process       = $(settings.process)\n" *
              "    approximation = $(settings.approximation)\n\n" *
              "does not describe anything. A scattering process has TWO photon vertices -- one for the incoming photon and one " *
              "for the outgoing -- so it cannot be evaluated with a single-vertex amplitude. Use SecondOrderGreen().")
    end

    # (b) not yet implemented
    if      !(settings.beamType isa Beam.PlaneWave)
        error("\n\nNo photon-scattering process is implemented yet for\n    beamType = $(settings.beamType)\n\n" *
              "Only Beam.PlaneWave() is available today. Scattering of TWISTED beams -- Beam.BesselBeam and " *
              "Beam.LaguerreGauss -- is a planned stage of this module and the beam side already exists in the Beam module, but " *
              "no amplitude here uses it: settings.beamType is not yet read by any pipeline. Ask for a plane wave, or wait.")
    elseif  settings.approximation == PhotonScattering.FormFactorApproximation()
        error("\n\nThe form-factor approximation is not implemented yet for\n    process = $(settings.process)\n\n" *
              "Only SecondOrderGreen() is available for a scattering process today. The form-factor route -- the Rayleigh " *
              "amplitude from the atomic form factor and the Compton one from the incoherent scattering function -- is planned " *
              "as the cheap counterpart to the full sum, and JAC's FormFactor module already supplies the ingredient; nothing " *
              "here calls it. Use SecondOrderGreen().")
    end

    return( nothing )
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
    PhotonScattering.checkCombination(settings)
    if      settings.process == PhotonScattering.BoundFreePairCreation()
        return( PhotonScattering.computePairCreationLines(finalMultiplet, initialMultiplet, nm, grid, settings; output=output) )
    elseif  settings.process in [PhotonScattering.RayleighScattering(), PhotonScattering.ComptonScattering(),
                                 PhotonScattering.RamanScattering()]
        return( PhotonScattering.computeRayleighLines(finalMultiplet, initialMultiplet, nm, grid, settings; output=output) )
    elseif  settings.process == PhotonScattering.ResonantScattering()
        return( PhotonScattering.computeResonantLines(finalMultiplet, initialMultiplet, nm, grid, settings; output=output) )
    else
        error("PhotonScattering: no pipeline for process $(settings.process).")
    end
end


"""
`PhotonScattering.displaySupportedCombinations(stream::IO)`
    ... to print which (beamType, process, approximation) combinations this module can evaluate and which it cannot, so that a user
        can ask BEFORE running rather than discover it by being refused. A neat table is printed but nothing is returned otherwise.

        The three axes are orthogonal by design and most of the grid is therefore empty. That is a feature -- the empty places are
        named so the intended scope is visible -- but it does mean the supported set is much smaller than the type system allows,
        and it is worth being able to see it at a glance.
"""
function displaySupportedCombinations(stream::IO)
    nx = 118
    println(stream, " ")
    println(stream, "  PhotonScattering: which combinations can be computed today")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    beamType        process                  approximation             state")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    PlaneWave       BoundFreePairCreation    FirstOrderVertex          YES")
    println(stream, "    PlaneWave       RayleighScattering       SecondOrderGreen          YES")
    println(stream, "    PlaneWave       RamanScattering          SecondOrderGreen          YES")
    println(stream, "    PlaneWave       ResonantScattering       SecondOrderGreen          YES  (needs settings.width > 0)")
    println(stream, "    PlaneWave       ComptonScattering        SecondOrderGreen          PARTLY -- see the note below")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    PlaneWave       any scattering process   FormFactorApproximation   not implemented")
    println(stream, "    BesselBeam      any                      any                       not implemented")
    println(stream, "    LaguerreGauss   any                      any                       not implemented")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    any             BoundFreePairCreation    SecondOrderGreen          MEANINGLESS -- one vertex only")
    println(stream, "    any             any scattering process   FirstOrderVertex          MEANINGLESS -- two vertices")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, " ")
    println(stream, "  ComptonScattering presently means RAMAN-type inelastic scattering into DISCRETE final levels, which is what")
    println(stream, "  a final-state Multiplet can express. The CONTINUUM Compton profile -- an ejected electron and the doubly")
    println(stream, "  differential d^2 sigma / dOmega dOmega_out -- is not implemented in JAC under any name, and no Line can carry")
    println(stream, "  it, a line holding one outgoing energy fixed by energy conservation.")
    println(stream, " ")

    return( nothing )
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
include("module-PhotonScattering-inc-resonant.jl")

end # module
