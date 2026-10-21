
"""
`module  JAC.PhotoIonizationInterference`
... a submodel of JAC that combines the DIRECT photoionization (PI) amplitude and the RESONANT amplitude -- excitation into an
    autoionizing level followed by its autoionization decay -- COHERENTLY, and that analyses the interference between them in the
    cross section (the Fano line shape) and, where requested, in the angular distribution of the photoelectron.

    Both routes lead to the SAME final state -- a one-electron-less ion plus a free electron of a given partial wave -- so their
    amplitudes have to be added before, and not after, they are squared. Per channel, i.e. per partial wave kappa of the outgoing
    electron and per total symmetry J_t of the (ion + electron) scattering state, the module forms

        A_total = A_PI  +  Sum_d  A_excite(i-->d) * A_decay(d-->f; kappa) / (E_i + omega - E_d + i Gamma_d/2) ,

    where the sum runs over those intermediate (autoionizing, above-threshold) levels d whose symmetry matches the channel. The
    direct term A_PI is the one of `PhotoIonization`, the excitation amplitude is a plain bound-bound E1 (or higher-multipole)
    amplitude from `PhotoEmission`'s own `Absorption()` kind -- the SAME amplitude `PhotoExcitation` itself calls -- and the decay
    amplitude is `AutoIonization`'s own Auger amplitude, converted to the continuum-phase convention `PhotoIonization` uses.

    THIS IS THE TIME-REVERSED COUNTERPART OF `PhotoRecombinationInterference`, not a copy of it: there RR/DR share a continuum
    electron on the INITIAL side of both terms; here PI and the resonant route share NEITHER state with a continuum electron
    attached on the same side as the OTHER term's non-continuum state -- the excitation step is bound-bound throughout (i and d are
    both ordinary bound levels, exactly as in a standalone `PhotoExcitation` computation) and the decay step ends in the SAME
    (f, kappa) continuum channel the direct term reaches. This structural symmetry is why NO extra degeneracy-ratio normalization
    constant is needed here, unlike the sqrt((2Jf+1)/(2Ji+1)) factor `PhotoRecombinationInterference.dielectronicNormalization`
    derives -- see `PhotoIonizationInterference.photoionizationNormalization`'s own docstring for the argument and for how it is
    checked rather than merely asserted.

    THE CONTINUUM-PHASE-CONVENTION MISMATCH, FOUND BY READING THE CODE, NOT BY ASSUMING THE TWO MODULES AGREE.
    `PhotoIonization.amplitude` multiplies its real reduced matrix element by i^(-l) exp(-i phase); `AutoIonization.amplitude`
    multiplies its own by i^(+l) exp(-i phase) -- the SAME phase, but the OPPOSITE l-prefactor, differing by i^(2l) = (-1)^l. That
    factor is harmless in anything that squares the amplitude (which is why a plain Auger rate never reveals it), but not in a
    coherent sum, where it reweights every partial wave against the others. `AutoIonization.amplitudeForPhotoionizationConvention`
    supplies the Auger decay amplitude already converted to `PhotoIonization`'s own convention, so that no (-1)^l is left implicit
    at any call site here.

    SCOPE, DELIBERATELY NARROW, the same choice `PhotoRecombinationInterference` made for the analogous reason: a resonant-decay
    width built only from the configurations the computation itself names (a LOWER bound, reported as such); E1 for the
    angular-distribution formula (`calcAnisotropy`); and NO linear-polarization observable, because -- unlike recombination, where
    the resonant step ends in a radiated photon whose polarization is a real observable -- the resonant step here ends in
    AUTOIONIZATION, not photon emission, so there is no analogous polarization to compute.
"""
module PhotoIonizationInterference


using  Printf, ..AutoIonization, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..PhotoEmission,
       ..PhotoIonization, ..Radial, ..TableStrings


"""
`struct  PhotoIonizationInterference.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of computing the interference of direct photoionization with resonant
        (autoionization-mediated) photoionization.

    + multipoles               ::Array{EmMultipole,1}   ... Multipoles of the absorbed photon that are taken into account.
    + gauges                   ::Array{UseGauge,1}      ... Gauges to be included into the computations.
    + photonEnergies           ::Array{Float64,1}       ... Photon energies, i.e. the scan ACROSS the resonance.
    + lValues                  ::Array{Int64,1}         ... Orbital angular momenta of the outgoing electron that are considered.
    + augerOperator            ::AbstractEeInteraction  ... Operator for the autoionization-decay amplitude.
    + includePI                ::Bool                   ... True, if the direct (photoionization) term is included.
    + includeResonance         ::Bool                   ... True, if the resonant (excitation + autoionization) term is included.
    + intermediateWidths       ::Dict{Int64,Float64}    ... Optional per-level override of the total width Gamma_d [a.u.], keyed on
                                                            the index of the intermediate level; where an entry exists it wins over
                                                            the width that is computed from the configurations named in the
                                                            computation.
    + calcAnisotropy           ::Bool                   ... True, if the anisotropy parameters beta_nu of the photoelectron are
                                                            computed.
    + printBefore              ::Bool                   ... True, if all selected pathways are printed before their evaluation.
    + pathwaySelection         ::PathwaySelection       ... Specifies the selected pathways, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                 ::Array{EmMultipole,1}
    gauges                     ::Array{UseGauge,1}
    photonEnergies             ::Array{Float64,1}
    lValues                    ::Array{Int64,1}
    augerOperator              ::AbstractEeInteraction
    includePI                  ::Bool
    includeResonance           ::Bool
    intermediateWidths         ::Dict{Int64,Float64}
    calcAnisotropy             ::Bool
    printBefore                ::Bool
    pathwaySelection           ::PathwaySelection
end


"""
`PhotoIonizationInterference.Settings()`
    ... constructor for an `empty` instance of PhotoIonizationInterference.Settings; a settings::Settings is returned.
"""
function Settings()
    Settings( EmMultipole[E1], UseGauge[UseCoulomb], Float64[], [0,1,2,3,4,5], CoulombInteraction(), true, true,
              Dict{Int64,Float64}(), false, false, PathwaySelection() )
end


"""
`PhotoIonizationInterference.Settings(set::PhotoIonizationInterference.Settings;`

        multipoles=..,          gauges=..,              photonEnergies=..,     lValues=..,
        augerOperator=..,       includePI=..,           includeResonance=..,   intermediateWidths=..,
        calcAnisotropy=..,      printBefore=..,          pathwaySelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::PhotoIonizationInterference.Settings; a
        settings::Settings is returned.
"""
function Settings(set::PhotoIonizationInterference.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,         gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,         lValues::Union{Nothing,Array{Int64,1}}=nothing,
    augerOperator::Union{Nothing,AbstractEeInteraction}=nothing,     includePI::Union{Nothing,Bool}=nothing,
    includeResonance::Union{Nothing,Bool}=nothing,                   intermediateWidths::Union{Nothing,Dict{Int64,Float64}}=nothing,
    calcAnisotropy::Union{Nothing,Bool}=nothing,                     printBefore::Union{Nothing,Bool}=nothing,
    pathwaySelection::Union{Nothing,PathwaySelection}=nothing)

    if  isnothing(multipoles)          multipolesx         = set.multipoles         else   multipolesx         = multipoles         end
    if  isnothing(gauges)              gaugesx             = set.gauges             else   gaugesx             = gauges             end
    if  isnothing(photonEnergies)      photonEnergiesx     = set.photonEnergies     else   photonEnergiesx     = photonEnergies     end
    if  isnothing(lValues)             lValuesx            = set.lValues            else   lValuesx            = lValues            end
    if  isnothing(augerOperator)       augerOperatorx      = set.augerOperator      else   augerOperatorx      = augerOperator      end
    if  isnothing(includePI)           includePIx          = set.includePI          else   includePIx          = includePI          end
    if  isnothing(includeResonance)    includeResonancex   = set.includeResonance   else   includeResonancex   = includeResonance   end
    if  isnothing(intermediateWidths)  intermediateWidthsx = set.intermediateWidths else   intermediateWidthsx = intermediateWidths end
    if  isnothing(calcAnisotropy)      calcAnisotropyx     = set.calcAnisotropy     else   calcAnisotropyx     = calcAnisotropy     end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end
    if  isnothing(pathwaySelection)    pathwaySelectionx   = set.pathwaySelection   else   pathwaySelectionx   = pathwaySelection   end

    Settings( multipolesx, gaugesx, photonEnergiesx, lValuesx, augerOperatorx, includePIx, includeResonancex, intermediateWidthsx,
              calcAnisotropyx, printBeforex, pathwaySelectionx )
end


# `Base.show(io::IO, settings::PhotoIonizationInterference.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::PhotoIonizationInterference.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
    println(io, "augerOperator:            $(settings.augerOperator)  ")
    println(io, "includePI:                $(settings.includePI)  ")
    println(io, "includeResonance:         $(settings.includeResonance)  ")
    println(io, "intermediateWidths:       $(settings.intermediateWidths)  ")
    println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "pathwaySelection:         $(settings.pathwaySelection)  ")
end


"""
`struct  PhotoIonizationInterference.Channel`
    ... defines a type for ONE total symmetry J_t of the (ion + photoelectron) scattering state, and carries the direct, the
        resonant and the summed amplitude separately, so that the three cross sections and hence the interference term can be
        formed afterwards.

    + symmetry                 ::LevelSymmetry                  ... Total J^parity of the scattering state.
    + amplitudesPI              ::Array{MultipoleAmplitude,1}    ... Direct (photoionization) amplitudes, one entry per multipole,
                                                                      both gauges.
    + amplitudesResonant        ::Array{MultipoleAmplitude,1}    ... Resonant amplitudes, one entry per multipole, both gauges.
    + amplitudes                ::Array{MultipoleAmplitude,1}    ... The COHERENT SUM, one entry per multipole, both gauges.
"""
struct  Channel
    symmetry                   ::LevelSymmetry
    amplitudesPI                ::Array{MultipoleAmplitude,1}
    amplitudesResonant           ::Array{MultipoleAmplitude,1}
    amplitudes                  ::Array{MultipoleAmplitude,1}
end


"""
`PhotoIonizationInterference.Channel()`
    ... constructor for an `empty` instance of PhotoIonizationInterference.Channel; a channel::Channel is returned.
"""
function Channel()
    Channel( LevelSymmetry(AngularJ64(0), Basics.plus), MultipoleAmplitude[], MultipoleAmplitude[], MultipoleAmplitude[] )
end


# `Base.show(io::IO, channel::PhotoIonizationInterference.Channel)`  ... prepares a proper printout of channel.
function Base.show(io::IO, channel::PhotoIonizationInterference.Channel)
    println(io, "symmetry:                 $(channel.symmetry)  ")
    println(io, "amplitudesPI:             $(channel.amplitudesPI)  ")
    println(io, "amplitudesResonant:       $(channel.amplitudesResonant)  ")
    println(io, "amplitudes:               $(channel.amplitudes)  ")
end


"""
`struct  PhotoIonizationInterference.PartialWave`
    ... defines a type for ONE partial wave of the outgoing photoelectron together with the total symmetries it serves; it mirrors
        `PhotoIonization.PartialWave`, since the radial orbital and the scattering phase belong to (energy, kappa) and not to any
        one total symmetry.

    + kappa                    ::Int64                  ... Partial wave of the outgoing photoelectron.
    + energy                   ::Float64                ... Energy of the photoelectron [a.u.].
    + phase                    ::Float64                ... Scattering phase of this partial wave.
    + channels                 ::Array{Channel,1}       ... The total symmetries that this partial wave serves.
"""
struct  PartialWave
    kappa                      ::Int64
    energy                     ::Float64
    phase                      ::Float64
    channels                   ::Array{Channel,1}
end


"""
`PhotoIonizationInterference.PartialWave()`
    ... constructor for an `empty` instance of PhotoIonizationInterference.PartialWave; a pw::PartialWave is returned.
"""
function PartialWave()
    PartialWave( 0, 0., 0., PhotoIonizationInterference.Channel[] )
end


# `Base.show(io::IO, pw::PhotoIonizationInterference.PartialWave)`  ... prepares a proper printout of pw.
function Base.show(io::IO, pw::PhotoIonizationInterference.PartialWave)
    println(io, "kappa:                    $(pw.kappa)  ")
    println(io, "energy:                   $(pw.energy)  ")
    println(io, "phase:                    $(pw.phase)  ")
    println(io, "channels:                 $(pw.channels)  ")
end


"""
`struct  PhotoIonizationInterference.Pathway`
    ... defines a type for the coherent photoionization of an initial ion level into a final level plus photoelectron, for ONE
        photon energy and through a set of intermediate (autoionizing) resonances.

    + initialLevel             ::Level                  ... Initial-(state) level of the atom or ion.
    + finalLevel               ::Level                  ... Final-(state) level of the ion.
    + intermediateLevels       ::Array{Level,1}         ... The resonances that were actually included in the resonant term.
    + photonEnergy             ::Float64                ... Energy of the absorbed photon [a.u.]; ONE value for both terms.
    + electronEnergy           ::Float64                ... Energy of the outgoing photoelectron [a.u.].
    + crossSection             ::EmProperty             ... Cross section from the SUMMED amplitude.
    + crossSectionPI           ::EmProperty             ... Cross section from the direct term alone.
    + crossSectionResonant     ::EmProperty             ... Cross section from the resonant term alone.
    + interference             ::EmProperty             ... crossSection - crossSectionPI - crossSectionResonant.
    + betaNu                   ::Array{EmProperty,1}    ... Anisotropy parameters beta_nu, nu = 2 and 4, from the summed
                                                            amplitude (E1 only).
    + partialWaves             ::Array{PartialWave,1}   ... Partial waves, each with the channels it serves.
"""
struct  Pathway
    initialLevel               ::Level
    finalLevel                 ::Level
    intermediateLevels         ::Array{Level,1}
    photonEnergy                ::Float64
    electronEnergy              ::Float64
    crossSection                ::EmProperty
    crossSectionPI              ::EmProperty
    crossSectionResonant         ::EmProperty
    interference                ::EmProperty
    betaNu                      ::Array{EmProperty,1}
    partialWaves                 ::Array{PartialWave,1}
end


"""
`PhotoIonizationInterference.Pathway()`
    ... constructor for an `empty` instance of PhotoIonizationInterference.Pathway; a pathway::Pathway is returned.
"""
function Pathway()
    Pathway( Level(), Level(), Level[], 0., 0., EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.),
             EmProperty[], PhotoIonizationInterference.PartialWave[] )
end


# `Base.show(io::IO, pathway::PhotoIonizationInterference.Pathway)`  ... prepares a proper printout of pathway.
function Base.show(io::IO, pathway::PhotoIonizationInterference.Pathway)
    println(io, "initialLevel:             $(pathway.initialLevel)  ")
    println(io, "finalLevel:               $(pathway.finalLevel)  ")
    println(io, "intermediateLevels:       $(pathway.intermediateLevels)  ")
    println(io, "photonEnergy:             $(pathway.photonEnergy)  ")
    println(io, "electronEnergy:           $(pathway.electronEnergy)  ")
    println(io, "crossSection:             $(pathway.crossSection)  ")
    println(io, "crossSectionPI:           $(pathway.crossSectionPI)  ")
    println(io, "crossSectionResonant:     $(pathway.crossSectionResonant)  ")
    println(io, "interference:             $(pathway.interference)  ")
    println(io, "betaNu:                   $(pathway.betaNu)  ")
    println(io, "partialWaves:             $(pathway.partialWaves)  ")
end


"""
`PhotoIonizationInterference.amplitudeResonant(mp::EmMultipole, kappa::Int64, phase::Float64, dLevel::Level, newcLevel::Level,
        initialLevel::Level, totalEnergy::Float64, gammaD::Float64, grid::Radial.Grid, settings::PhotoIonizationInterference.Settings,
        exciteCache::Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC} = Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC}())`
    ... to compute the resonant contribution of ONE intermediate level d to the photoionization amplitude of one channel, i.e.

            A_excite(i-->d) * A_decay(d-->f,kappa) / (E_i + omega - E_d + i Gamma_d/2) ,

        with NO extra normalization constant -- see `PhotoIonizationInterference.photoionizationNormalization`'s own docstring for
        why none is needed here, unlike the sqrt((2Jf+1)/(2Ji+1)) factor `PhotoRecombinationInterference.dielectronicNormalization`
        derives for the structurally different RR/DR case. A value::EmPropertyC is returned.

        The excitation amplitude is a plain bound-bound `PhotoEmission.amplitude(Absorption(),...)` call -- exactly the amplitude
        `PhotoExcitation` itself uses -- carrying no continuum-phase ambiguity at all, since both i and d are ordinary bound
        levels. The decay amplitude comes from `AutoIonization.amplitudeForPhotoionizationConvention`, which is
        `AutoIonization.amplitude` converted to `PhotoIonization`'s own i^(-l) exp(-i phase) convention; using the plain
        `AutoIonization.amplitude` here instead would reweight the partial waves against one another by (-1)^l, exactly the
        mistake `AutoIonization.captureAmplitude`'s own docstring warns against for the recombination case.

        THE EXCITATION AMPLITUDE DOES NOT DEPEND ON kappa, phase OR newcLevel -- only on (mp, dLevel, initialLevel) -- so the SAME
        bound-bound radial integral is otherwise recomputed once per partial wave and per channel that shares this d-level, for
        every photon energy of the scan. `exciteCache` memoizes it; the caller supplies ONE dict shared across the whole
        `computePathways` run (dLevel and initialLevel repeat across the entire photon-energy scan, so this is the dominant
        saving), and an isolated call may simply omit it.
"""
function amplitudeResonant(mp::EmMultipole, kappa::Int64, phase::Float64, dLevel::Level, newcLevel::Level, initialLevel::Level,
                           totalEnergy::Float64, gammaD::Float64, grid::Radial.Grid, settings::PhotoIonizationInterference.Settings,
                           exciteCache::Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC} = Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC}())
    denominator = Complex(totalEnergy - dLevel.energy, 0.5 * gammaD)
    if  abs(denominator) < 1.0e-30    return( EmPropertyC(0.0im) )    end
    #
    cacheKey = (mp, dLevel.index, initialLevel.index)
    if  haskey(exciteCache, cacheKey)
        ampExcite = exciteCache[cacheKey]
    else
        omega_id = dLevel.energy - initialLevel.energy
        if  string(mp)[1] == 'E'
            aC = PhotoEmission.amplitude(Absorption(), mp, Basics.Coulomb,   omega_id, dLevel, initialLevel, grid, display=false, printout=false)
            aB = PhotoEmission.amplitude(Absorption(), mp, Basics.Babushkin, omega_id, dLevel, initialLevel, grid, display=false, printout=false)
            ampExcite = EmPropertyC(aC, aB)
        else
            aM = PhotoEmission.amplitude(Absorption(), mp, Basics.Magnetic,  omega_id, dLevel, initialLevel, grid, display=false, printout=false)
            ampExcite = EmPropertyC(aM)
        end
        exciteCache[cacheKey] = ampExcite
    end
    if  ampExcite == EmPropertyC(0.0im)    return( EmPropertyC(0.0im) )    end
    #
    ampDecay = AutoIonization.amplitudeForPhotoionizationConvention(settings.augerOperator, kappa, phase, newcLevel, dLevel, grid;
                                                                    printout=false)

    return( (ampDecay / denominator) * ampExcite )
end


"""
`PhotoIonizationInterference.photoionizationNormalization(initialLevel::Level, dLevel::Level)`
    ... supplies the constant N that would put the resonant amplitude on the same footing as the direct one, analogous to
        `PhotoRecombinationInterference.dielectronicNormalization`. STRUCTURAL ARGUMENT FOR N = 1, not assumed: the direct
        amplitude's bra is the BOUND initial level i and its ket is the (f, kappa) J_t CONTINUUM state; the resonant amplitude's
        two halves share exactly those same two ends -- `PhotoEmission.amplitude(Absorption(),...)` has bra i and ket d (both
        bound, d carrying the SAME J_t the channel itself carries), and `AutoIonization.amplitudeForPhotoionizationConvention` has
        bra d and ket (f, kappa) J_t -- so NEITHER half introduces a continuum electron on the side the direct term does not
        already carry one. This is the opposite of the RR/DR case, where the capture amplitude's bra (i, kappa) J_t already
        carries a continuum electron that the direct term's OWN bra does too, while the stabilization amplitude's bra/ket (d, f)
        are both bound -- an asymmetric structure that is exactly what forces `dielectronicNormalization`'s derived factor.

        THIS VALUE IS NOT TRUSTED ON THE STRENGTH OF THE ARGUMENT ALONE: example-Dy.jl's own branch measures the ratio of the
        integrated resonant-alone cross section (interference and the direct term both switched off) to the independently
        computed sum rule `2 pi^2 alpha f_exc BR_auger` -- `f_exc` from `PhotoExcitation.computeAmplitudesProperties` and
        `BR_auger` from `AutoIonization.computeLines`/`displayLifetimes`, both already-validated existing functions -- and reports
        the measured ratio, exactly the discipline `dielectronicNormalization`'s own docstring describes for its case. A
        value::Float64 is returned.
"""
function photoionizationNormalization(initialLevel::Level, dLevel::Level)

    return( 1.0 )
end


"""
`PhotoIonizationInterference.computeAmplitudesProperties(pathway::PhotoIonizationInterference.Pathway,
        intermediateMultiplet::Multiplet, widths::Dict{Int64,Float64}, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
        settings::PhotoIonizationInterference.Settings,
        exciteCache::Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC} = Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC}())`
    ... to compute, for the given pathway, the direct and the resonant amplitude of every channel and to add them coherently; the
        three cross sections, the interference term and the anisotropy parameters follow from them. A newPathway::Pathway is
        returned.

        ONE continuum orbital is generated per partial wave, exactly as in `PhotoIonization.computeAmplitudesProperties`, and the
        excitation amplitude of the resonant term is evaluated at the resonance's OWN energy E_d - E_i, not at the photon energy
        of the overall process -- the two coincide only exactly at resonance, same distinction
        `PhotoRecombinationInterference.amplitudeDR`'s own docstring draws for its stabilizing amplitude.

        `exciteCache` is passed straight through to `PhotoIonizationInterference.amplitudeResonant`, which explains what it saves
        and why; `computePathways` is the caller that shares ONE dict across every pathway of the scan.
"""
function computeAmplitudesProperties(pathway::PhotoIonizationInterference.Pathway, intermediateMultiplet::Multiplet,
                                     widths::Dict{Int64,Float64}, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                     settings::PhotoIonizationInterference.Settings,
                                     exciteCache::Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC} = Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC}())
    contSettings = Continuum.Settings(false, nrContinuum)
    redILevel    = Basics.generateLevelWithSymmetryReducedBasis(pathway.initialLevel, pathway.initialLevel.basis.subshells)
    newfLevel    = Basics.generateLevelWithSymmetryReducedBasis(pathway.finalLevel, redILevel.basis.subshells)
    totalEnergy  = pathway.initialLevel.energy + pathway.photonEnergy
    newPartialWaves = PhotoIonizationInterference.PartialWave[]

    for  pw  in  pathway.partialWaves
        cSubshell = Subshell(101, pw.kappa)
        newiLevel = Basics.generateLevelWithExtraSubshell(cSubshell, redILevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(pathway.electronEnergy, cSubshell, newfLevel, nm, grid, contSettings)
        newChannels = PhotoIonizationInterference.Channel[]
        #
        for  ch  in  pw.channels
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newfLevel)
            ampsPI = MultipoleAmplitude[];    ampsRes = MultipoleAmplitude[];    ampsTot = MultipoleAmplitude[]
            #
            for  mp  in  settings.multipoles
                # The direct, photoionization term
                if  settings.includePI
                    if  string(mp)[1] == 'E'
                        aC = PhotoIonization.amplitude("photoionization", mp, Basics.Coulomb,   pw.kappa, phase,
                                                       pathway.photonEnergy, newcLevel, newiLevel, grid)
                        aB = PhotoIonization.amplitude("photoionization", mp, Basics.Babushkin, pw.kappa, phase,
                                                       pathway.photonEnergy, newcLevel, newiLevel, grid)
                        ampPI = EmPropertyC(aC, aB)
                    else
                        aM = PhotoIonization.amplitude("photoionization", mp, Basics.Magnetic,  pw.kappa, phase,
                                                       pathway.photonEnergy, newcLevel, newiLevel, grid)
                        ampPI = EmPropertyC(aM)
                    end
                else
                    ampPI = EmPropertyC(0.0im)
                end
                # The resonant, excitation-then-autoionization term
                ampRes = EmPropertyC(0.0im)
                if  settings.includeResonance
                    for  dLevel  in  pathway.intermediateLevels
                        if  LevelSymmetry(dLevel.J, dLevel.parity) != ch.symmetry    continue    end
                        nFactor = PhotoIonizationInterference.photoionizationNormalization(pathway.initialLevel, dLevel)
                        ampRes  = ampRes + nFactor * PhotoIonizationInterference.amplitudeResonant(mp, pw.kappa, phase, dLevel,
                                            newcLevel, pathway.initialLevel, totalEnergy, get(widths, dLevel.index, 0.), grid,
                                            settings, exciteCache)
                    end
                end
                push!( ampsPI,  MultipoleAmplitude(mp, ampPI) )
                push!( ampsRes, MultipoleAmplitude(mp, ampRes) )
                push!( ampsTot, MultipoleAmplitude(mp, ampPI + ampRes) )
            end
            push!( newChannels, PhotoIonizationInterference.Channel(ch.symmetry, ampsPI, ampsRes, ampsTot) )
        end
        push!( newPartialWaves, PhotoIonizationInterference.PartialWave(pw.kappa, pathway.electronEnergy, phase, newChannels) )
    end
    #
    # The three cross sections, each through the SAME function of PhotoIonization, so that they are commensurable by construction
    csTot = PhotoIonization.computeCrossSection( PhotoIonizationInterference.toPhotoIonizationWaves(newPartialWaves, :total),
                                                 pathway.photonEnergy, pathway.initialLevel )
    csPI  = PhotoIonization.computeCrossSection( PhotoIonizationInterference.toPhotoIonizationWaves(newPartialWaves, :pi),
                                                 pathway.photonEnergy, pathway.initialLevel )
    csRes = PhotoIonization.computeCrossSection( PhotoIonizationInterference.toPhotoIonizationWaves(newPartialWaves, :resonant),
                                                 pathway.photonEnergy, pathway.initialLevel )
    interference = csTot - csPI - csRes
    #
    betaNu = EmProperty[]
    if  settings.calcAnisotropy
        totWaves = PhotoIonizationInterference.toPhotoIonizationWaves(newPartialWaves, :total)
        for  nu in [2, 4]
            push!( betaNu, PhotoIonization.computeAngularBeta(pathway.initialLevel, pathway.finalLevel, totWaves) )
        end
    end

    newPathway = PhotoIonizationInterference.Pathway( pathway.initialLevel, pathway.finalLevel, pathway.intermediateLevels,
                     pathway.photonEnergy, pathway.electronEnergy, csTot, csPI, csRes, interference, betaNu, newPartialWaves )

    return( newPathway )
end


"""
`PhotoIonizationInterference.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
        initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationInterference.Settings;
        output::Bool=true)`
    ... to compute the coherent photoionization amplitudes and all requested properties for the selected pathways and for every
        photon energy of the settings; this is the standard entry point of this module. A list of
        pathways::Array{PhotoIonizationInterference.Pathway,1} is returned if output = true, and nothing otherwise.
"""
function computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                         grid::Radial.Grid, settings::PhotoIonizationInterference.Settings; output::Bool=true)
    println("")
    printstyled("PhotoIonizationInterference.computePathways(): The computation of the coherent PI + resonant amplitudes starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------------------ \n",
                color=:light_green)
    println("")
    widths   = PhotoIonizationInterference.computeTotalWidths(intermediateMultiplet, initialMultiplet, finalMultiplet, nm, grid, settings)
    pathways = PhotoIonizationInterference.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoIonizationInterference.displayPathways(stdout, pathways, widths)    end
    maxEnergy = 0.;    for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    # dLevel and initialLevel repeat across the WHOLE photon-energy scan, so one cache shared by every pathway removes almost
    # all of the redundancy; see PhotoIonizationInterference.amplitudeResonant's own docstring for why this exists.
    exciteCache = Dict{Tuple{EmMultipole,Int64,Int64},EmPropertyC}()
    newPathways = PhotoIonizationInterference.Pathway[]
    for  pathway in pathways
        push!( newPathways, PhotoIonizationInterference.computeAmplitudesProperties(pathway, intermediateMultiplet, widths, nm, grid,
                                                                                    nrContinuum, settings, exciteCache) )
    end
    PhotoIonizationInterference.displayResults(stdout, newPathways, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoIonizationInterference.displayResults(iostream, newPathways, settings)    end

    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`PhotoIonizationInterference.computeTotalWidths(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
        finalMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationInterference.Settings)`
    ... to determine the total width Gamma_d of every intermediate level, as the sum of its Auger rates into the levels of the
        final multiplet and of its radiative rates into the levels of the initial multiplet. A widths::Dict{Int64,Float64} is
        returned, keyed on the index of the intermediate level.

        This width is built ONLY from the configurations that the computation names, and is therefore a LOWER bound: any decay
        channel outside them is missing, which makes the resonance too narrow and its peak too high. The breakdown is printed for
        exactly that reason, and `settings.intermediateWidths` overrides the computed value wherever the user supplies one.

        THE TWO DECAY DIRECTIONS, named explicitly because they are the reverse of `PhotoRecombinationInterference`'s own
        `computeTotalWidths`: there, Gamma_a decays INTO the INITIAL (bare ion) multiplet and Gamma_r INTO the FINAL (recombined)
        multiplet, because recombination's bare ion sits on the initial side. Here Gamma_a decays INTO the FINAL (one-electron-
        less ion) multiplet and Gamma_r INTO the INITIAL (bound, lower-lying) multiplet, because photoionization's bare ion sits
        on the final side.
"""
function computeTotalWidths(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, finalMultiplet::Multiplet, nm::Nuclear.Model,
                            grid::Radial.Grid, settings::PhotoIonizationInterference.Settings)
    widths = Dict{Int64,Float64}()
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(); printBefore=false, operator=settings.augerOperator)
    photonSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=settings.multipoles,
                                            gauges=[UseCoulomb, UseBabushkin], printBefore=false)
    augerLines  = AutoIonization.computeLines(finalMultiplet, intermediateMultiplet, nm, grid, augerSettings; output=true)
    photonLines = PhotoEmission.computeLines(initialMultiplet, intermediateMultiplet, grid, photonSettings; output=true)
    #
    println("\n  Total widths Gamma_d of the intermediate levels, from the configurations named in this computation only:")
    println("  ", TableStrings.hLine(96))
    println("     level     J^P            Gamma_a [a.u.]      Gamma_r [a.u.]      Gamma_d [a.u.]     source")
    println("  ", TableStrings.hLine(96))
    for  dLevel in intermediateMultiplet.levels
        gammaA = 0.;    gammaR = 0.
        for  ln in augerLines    if  ln.initialLevel.index == dLevel.index   gammaA = gammaA + ln.totalRate            end   end
        for  ln in photonLines   if  ln.initialLevel.index == dLevel.index   gammaR = gammaR + ln.photonRate.Coulomb   end   end
        gammaD = gammaA + gammaR;    source = "computed"
        if  haskey(settings.intermediateWidths, dLevel.index)
            gammaD = settings.intermediateWidths[dLevel.index];    source = "user override"
        end
        widths[dLevel.index] = gammaD
        sym = LevelSymmetry(dLevel.J, dLevel.parity)
        println("  " * TableStrings.center(8, string(dLevel.index); na=2) * TableStrings.center(10, string(sym); na=2) *
                @sprintf("%.6e", gammaA) * "        " * @sprintf("%.6e", gammaR) * "        " * @sprintf("%.6e", gammaD) *
                "     " * source)
    end
    println("  ", TableStrings.hLine(96))
    println("  Decay channels outside the named configurations are NOT included; the widths above are lower bounds.\n")

    return( widths )
end


"""
`PhotoIonizationInterference.determineChannels(finalLevel::Level, initialLevel::Level,
        settings::PhotoIonizationInterference.Settings)`
    ... to determine the partial waves of the outgoing photoelectron and, for each of them, the total symmetries of the
        (ion + electron) scattering state that contribute; the list is that of `PhotoIonization.determineChannels` and is obtained
        from it, so that the direct term is bookkept exactly as in that module. An Array{PhotoIonizationInterference.PartialWave,1}
        is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonizationInterference.Settings)
    piSettings = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=settings.multipoles, gauges=settings.gauges,
                                          lValues=settings.lValues)
    piWaves    = PhotoIonization.determineChannels(finalLevel, initialLevel, piSettings)
    partialWaves = PhotoIonizationInterference.PartialWave[]
    for  pw in piWaves
        channels = PhotoIonizationInterference.Channel[]
        for  ch in pw.channels
            push!( channels, PhotoIonizationInterference.Channel(ch.symmetry, MultipoleAmplitude[], MultipoleAmplitude[],
                                                                  MultipoleAmplitude[]) )
        end
        push!( partialWaves, PhotoIonizationInterference.PartialWave(pw.kappa, pw.energy, pw.phase, channels) )
    end

    return( partialWaves )
end


"""
`PhotoIonizationInterference.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
        initialMultiplet::Multiplet, settings::PhotoIonizationInterference.Settings)`
    ... to determine the list of pathways for which the coherent photoionization amplitude is to be computed, one per (initial
        level, final level, photon energy). Each pathway carries the intermediate levels that SURVIVE the pathway selection for
        its own (i, f) pair, so that a single resonance -- or a chosen subset of them -- can be studied in isolation. Which of the
        surviving levels then actually contributes is decided channel by channel through the symmetry condition J_d = J_t. An
        Array{PhotoIonizationInterference.Pathway,1} is returned.
"""
function determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                           settings::PhotoIonizationInterference.Settings)
    pathways = PhotoIonizationInterference.Pathway[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            dLevels = Level[]
            for  dLevel in intermediateMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, dLevel, fLevel, settings.pathwaySelection)   push!( dLevels, dLevel )   end
            end
            if  length(dLevels) == 0    continue    end
            for  omega in settings.photonEnergies
                omegaAu = Defaults.convertUnits("energy: to atomic", omega)
                electronEnergy = omegaAu - (fLevel.energy - iLevel.energy)
                if  electronEnergy < 0.    continue    end
                pWaves = PhotoIonizationInterference.determineChannels(fLevel, iLevel, settings)
                push!( pathways, PhotoIonizationInterference.Pathway(iLevel, fLevel, dLevels, omegaAu, electronEnergy,
                           EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty[], pWaves) )
            end
        end
    end

    return( pathways )
end


"""
`PhotoIonizationInterference.displayPathways(stream::IO, pathways::Array{PhotoIonizationInterference.Pathway,1},
        widths::Dict{Int64,Float64})`
    ... to display a list of the pathways that have been selected for the computations, together with the photon and electron
        energies and the number of partial waves. A neat table is printed but nothing is returned otherwise.
"""
function displayPathways(stream::IO, pathways::Array{PhotoIonizationInterference.Pathway,1}, widths::Dict{Int64,Float64})
    nx = 118
    println(stream, " ")
    println(stream, "  Selected photoionization-interference pathways:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2)
    sa = sa * TableStrings.center(16, "E_photon"; na=2)
    sa = sa * TableStrings.center(16, "E_electron"; na=2)
    sa = sa * TableStrings.center(14, "no. of kappa"; na=2)
    sa = sa * TableStrings.center(16, "no. of resonances"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  pathway in pathways
        isym = LevelSymmetry(pathway.initialLevel.J, pathway.initialLevel.parity)
        fsym = LevelSymmetry(pathway.finalLevel.J,   pathway.finalLevel.parity)
        sa   = "  " * TableStrings.center(18, TableStrings.levels_if(pathway.initialLevel.index, pathway.finalLevel.index); na=2)
        sa   = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=2)
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "  "
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "  "
        sa   = sa * TableStrings.center(14, string(length(pathway.partialWaves)); na=2)
        sa   = sa * TableStrings.center(16, string(length(pathway.intermediateLevels)); na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotoIonizationInterference.displayResults(stream::IO, pathways::Array{PhotoIonizationInterference.Pathway,1},
        settings::PhotoIonizationInterference.Settings)`
    ... to display the three cross sections, the interference term and -- where requested -- the anisotropy parameter beta_2, for
        every pathway and photon energy. A neat table is printed but nothing is returned otherwise.
"""
function displayResults(stream::IO, pathways::Array{PhotoIonizationInterference.Pathway,1},
                        settings::PhotoIonizationInterference.Settings)
    nx = 123
    println(stream, " ")
    println(stream, "  Coherent PI + resonant photoionization:  cross sections [a.u.] and the interference term, Coulomb gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=1)
    sa = sa * TableStrings.center(15, "E_photon"; na=1)
    sa = sa * TableStrings.center(15, "sigma total"; na=1)
    sa = sa * TableStrings.center(15, "sigma PI"; na=1)
    sa = sa * TableStrings.center(15, "sigma resonant"; na=1)
    sa = sa * TableStrings.center(15, "interference"; na=1)
    sa = sa * TableStrings.center(13, "beta_2"; na=1)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  pathway in pathways
        sa = "  " * TableStrings.center(16, TableStrings.levels_if(pathway.initialLevel.index, pathway.finalLevel.index); na=1)
        sa = sa * @sprintf("%14.6e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy)) * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSection.Coulomb)         * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSectionPI.Coulomb)       * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSectionResonant.Coulomb) * " "
        sa = sa * @sprintf("%14.6e", pathway.interference.Coulomb)        * " "
        if  length(pathway.betaNu) > 0   sa = sa * @sprintf("%12.5e", pathway.betaNu[1].Coulomb) * " "
        else                             sa = sa * TableStrings.center(13, "--"; na=1)
        end
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end


"""
`PhotoIonizationInterference.toPhotoIonizationWaves(partialWaves::Array{PhotoIonizationInterference.PartialWave,1}, which::Symbol)`
    ... to package the amplitudes of the given partial waves into `Array{PhotoIonization.PartialWave,1}`, so that the cross section
        and the angular-distribution formula can be taken from that module unchanged. `which` selects which set of amplitudes
        travels: `:total` the coherent sum, `:pi` the direct term alone and `:resonant` the resonant term alone. An
        Array{PhotoIonization.PartialWave,1} is returned.

        This is what makes the limit includeResonance = false reproduce `PhotoIonization` by CONSTRUCTION: in that limit the very
        same numbers are handed to the very same functions.
"""
function toPhotoIonizationWaves(partialWaves::Array{PhotoIonizationInterference.PartialWave,1}, which::Symbol)
    piWaves = PhotoIonization.PartialWave[]
    for  pw in partialWaves
        piChannels = PhotoIonization.Channel[]
        for  ch in pw.channels
            if      which == :total      amps = ch.amplitudes
            elseif  which == :pi         amps = ch.amplitudesPI
            elseif  which == :resonant   amps = ch.amplitudesResonant
            else    error("stop a: unknown selector $which")
            end
            push!( piChannels, PhotoIonization.Channel(ch.symmetry, amps) )
        end
        push!( piWaves, PhotoIonization.PartialWave(pw.kappa, pw.energy, pw.phase, piChannels) )
    end

    return( piWaves )
end

end # module
