
"""
`module  JAC.PhotoRecombinationInterference`
... a submodel of JAC that combines the RADIATIVE (RR) and the DIELECTRONIC (DR) recombination amplitudes COHERENTLY, and that analyses
    the interference between them in the cross section as well as in the angular distribution and the linear polarization of the emitted
    photon.

    Both processes lead to the same final state -- a recombined ion plus one photon -- so their amplitudes have to be added before, and not
    after, they are squared.  Per channel, i.e. per partial wave kappa of the incoming electron, per total symmetry J_t of the (ion +
    electron) scattering state and per multipole of the emitted photon, the module forms

        A_total = A_RR  +  N * Sum_d  A_capture(d; kappa, J_t) * A_radiative(f <-- d; multipole) / (E_i + E - E_d + i Gamma_d/2)

    where the sum runs over those intermediate (doubly-excited, autoionizing) levels d whose symmetry matches the channel; the electron-
    electron interaction is a scalar, so only J_d = J_t contributes.  The direct term A_RR is the one of `PhotoRecombination`, the capture
    amplitude is built from the Auger amplitude of `AutoIonization`, and the stabilizing amplitude is the one of `PhotoEmission`.

    THE CAPTURE AMPLITUDE IS NOT conj(Auger), even though `ElectronCapture.amplitude` forms it that way.  The two source modules attach
    DIFFERENT phase conventions to the continuum state -- `PhotoRecombination` uses i^l exp(+i phi) and `AutoIonization` i^l exp(-i phi) --
    so the amplitude that belongs to PhotoRecombination's state is `AutoIonization.amplitude * exp(2 i phi)`, which differs from the
    conjugate by (-1)^l.  `ElectronCapture` is not wrong for its own purpose, since it only ever squares the result; but a coherent sum
    does not square first, and (-1)^l depends on the partial wave, so the difference is not even an overall sign.  See `amplitudeDR`.

    Both terms are evaluated at ONE AND THE SAME photon energy, omega = E + E_i - E_f, which is what energy conservation for the overall
    process requires and what makes the two amplitudes coherent at all.  It is NOT E_d - E_f; the two coincide only at exact resonance.

    Because the summed amplitudes are written into the very slots of a `PhotoRecombination.Line`, the cross section and every anisotropy
    parameter follow from `PhotoRecombination.computeCrossSectionForMultipoles` and `PhotoRecombination.computeAnisotropyParameter` without
    any new angular-momentum algebra -- and the limit `includeDR = false` reproduces that module by construction rather than by accident.
    Only the linear polarization is new here; see `computeLinearPolarization`.

    Reference for the physics: X.-M. Tong, Phys. Rev. A 107 (2023) 052801, for Be-like ions produced from Li-like ones through the KLL
    resonances, where the polarization is found to change sharply as the electron energy crosses the resonance; and A. J. Gonzalez Martinez
    et al. / Knapp et al., Phys. Rev. Lett. 74 (1995) 54, for the original observation of the Fano asymmetry in uranium.
"""
module PhotoRecombinationInterference


using  Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear,
       ..PhotoEmission, ..PhotoRecombination, ..Radial, ..TableStrings


"""
`struct  PhotoRecombinationInterference.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of computing the interference of radiative and dielectronic recombination.

    + multipoles               ::Array{EmMultipole,1}   ... Multipoles of the emitted photon that are taken into account.
    + gauges                   ::Array{UseGauge,1}      ... Gauges to be included into the computations.
    + electronEnergies         ::Array{Float64,1}       ... Energies of the incoming electron, i.e. the scan ACROSS the resonance.
    + maxKappa                 ::Int64                  ... Maximum kappa of the partial waves of the incoming electron.
    + augerOperator            ::AbstractEeInteraction  ... Operator for the dielectronic-capture amplitudes.
    + includeRR                ::Bool                   ... True, if the direct (radiative-recombination) term is included.
    + includeDR                ::Bool                   ... True, if the resonant (dielectronic-recombination) term is included.
    + intermediateWidths       ::Dict{Int64,Float64}    ... Optional per-level override of the total width Gamma_d [a.u.], keyed on the
                                                            index of the intermediate level; where an entry exists it wins over the width
                                                            that is computed from the configurations named in the computation.
    + calcAnisotropy           ::Bool                   ... True, if the anisotropy parameters beta_nu of the photon are computed.
    + calcPolarization         ::Bool                   ... True, if the linear polarization of the photon is computed.
    + printBefore              ::Bool                   ... True, if all selected pathways are printed before their evaluation.
    + pathwaySelection         ::PathwaySelection       ... Specifies the selected pathways, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                 ::Array{EmMultipole,1}
    gauges                     ::Array{UseGauge,1}
    electronEnergies           ::Array{Float64,1}
    maxKappa                   ::Int64
    augerOperator              ::AbstractEeInteraction
    includeRR                  ::Bool
    includeDR                  ::Bool
    intermediateWidths         ::Dict{Int64,Float64}
    calcAnisotropy             ::Bool
    calcPolarization           ::Bool
    printBefore                ::Bool
    pathwaySelection           ::PathwaySelection
end


"""
`PhotoRecombinationInterference.Settings()`
    ... constructor for an `empty` instance of PhotoRecombinationInterference.Settings; a settings::Settings is returned.
"""
function Settings()
    Settings( EmMultipole[E1], UseGauge[UseCoulomb], Float64[], 2, CoulombInteraction(), true, true,
              Dict{Int64,Float64}(), true, true, false, PathwaySelection() )
end


"""
`PhotoRecombinationInterference.Settings(set::PhotoRecombinationInterference.Settings;`

        multipoles=..,          gauges=..,              electronEnergies=..,    maxKappa=..,
        augerOperator=..,       includeRR=..,           includeDR=..,           intermediateWidths=..,
        calcAnisotropy=..,      calcPolarization=..,    printBefore=..,         pathwaySelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::PhotoRecombinationInterference.Settings; a
        settings::Settings is returned.
"""
function Settings(set::PhotoRecombinationInterference.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,         gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    electronEnergies::Union{Nothing,Array{Float64,1}}=nothing,       maxKappa::Union{Nothing,Int64}=nothing,
    augerOperator::Union{Nothing,AbstractEeInteraction}=nothing,     includeRR::Union{Nothing,Bool}=nothing,
    includeDR::Union{Nothing,Bool}=nothing,                          intermediateWidths::Union{Nothing,Dict{Int64,Float64}}=nothing,
    calcAnisotropy::Union{Nothing,Bool}=nothing,                     calcPolarization::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                        pathwaySelection::Union{Nothing,PathwaySelection}=nothing)

    if  isnothing(multipoles)          multipolesx         = set.multipoles         else   multipolesx         = multipoles         end
    if  isnothing(gauges)              gaugesx             = set.gauges             else   gaugesx             = gauges             end
    if  isnothing(electronEnergies)    electronEnergiesx   = set.electronEnergies   else   electronEnergiesx   = electronEnergies   end
    if  isnothing(maxKappa)            maxKappax           = set.maxKappa           else   maxKappax           = maxKappa           end
    if  isnothing(augerOperator)       augerOperatorx      = set.augerOperator      else   augerOperatorx      = augerOperator      end
    if  isnothing(includeRR)           includeRRx          = set.includeRR          else   includeRRx          = includeRR          end
    if  isnothing(includeDR)           includeDRx          = set.includeDR          else   includeDRx          = includeDR          end
    if  isnothing(intermediateWidths)  intermediateWidthsx = set.intermediateWidths else   intermediateWidthsx = intermediateWidths end
    if  isnothing(calcAnisotropy)      calcAnisotropyx     = set.calcAnisotropy     else   calcAnisotropyx     = calcAnisotropy     end
    if  isnothing(calcPolarization)    calcPolarizationx   = set.calcPolarization   else   calcPolarizationx   = calcPolarization   end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end
    if  isnothing(pathwaySelection)    pathwaySelectionx   = set.pathwaySelection   else   pathwaySelectionx   = pathwaySelection   end

    Settings( multipolesx, gaugesx, electronEnergiesx, maxKappax, augerOperatorx, includeRRx, includeDRx, intermediateWidthsx,
              calcAnisotropyx, calcPolarizationx, printBeforex, pathwaySelectionx )
end


# `Base.show(io::IO, settings::PhotoRecombinationInterference.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::PhotoRecombinationInterference.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "electronEnergies:         $(settings.electronEnergies)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "augerOperator:            $(settings.augerOperator)  ")
    println(io, "includeRR:                $(settings.includeRR)  ")
    println(io, "includeDR:                $(settings.includeDR)  ")
    println(io, "intermediateWidths:       $(settings.intermediateWidths)  ")
    println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
    println(io, "calcPolarization:         $(settings.calcPolarization)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "pathwaySelection:         $(settings.pathwaySelection)  ")
end


"""
`struct  PhotoRecombinationInterference.Channel`
    ... defines a type for ONE total symmetry J_t of the (ion + free electron) scattering state, and carries the direct, the resonant and
        the summed amplitude separately, so that the three cross sections and hence the interference term can be formed afterwards.

    + symmetry                 ::LevelSymmetry                  ... Total J^parity of the scattering state.
    + amplitudesRR             ::Array{MultipoleAmplitude,1}    ... Direct (RR) amplitudes, one entry per multipole, both gauges.
    + amplitudesDR             ::Array{MultipoleAmplitude,1}    ... Resonant (DR) amplitudes, one entry per multipole, both gauges.
    + amplitudes               ::Array{MultipoleAmplitude,1}    ... The COHERENT SUM, one entry per multipole, both gauges.
"""
struct  Channel
    symmetry                   ::LevelSymmetry
    amplitudesRR               ::Array{MultipoleAmplitude,1}
    amplitudesDR               ::Array{MultipoleAmplitude,1}
    amplitudes                 ::Array{MultipoleAmplitude,1}
end


"""
`PhotoRecombinationInterference.Channel()`
    ... constructor for an `empty` instance of PhotoRecombinationInterference.Channel; a channel::Channel is returned.
"""
function Channel()
    Channel( LevelSymmetry(AngularJ64(0), Basics.plus), MultipoleAmplitude[], MultipoleAmplitude[], MultipoleAmplitude[] )
end


# `Base.show(io::IO, channel::PhotoRecombinationInterference.Channel)`  ... prepares a proper printout of channel.
function Base.show(io::IO, channel::PhotoRecombinationInterference.Channel)
    println(io, "symmetry:                 $(channel.symmetry)  ")
    println(io, "amplitudesRR:             $(channel.amplitudesRR)  ")
    println(io, "amplitudesDR:             $(channel.amplitudesDR)  ")
    println(io, "amplitudes:               $(channel.amplitudes)  ")
end


"""
`struct  PhotoRecombinationInterference.PartialWave`
    ... defines a type for ONE partial wave of the incoming free electron together with the total symmetries it serves; it mirrors
        `PhotoRecombination.PartialWave`, since the radial orbital and the scattering phase belong to (energy, kappa) and not to any one
        total symmetry.

    + kappa                    ::Int64                  ... Partial wave of the incoming free electron.
    + energy                   ::Float64                ... Energy of the free electron [a.u.].
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
`PhotoRecombinationInterference.PartialWave()`
    ... constructor for an `empty` instance of PhotoRecombinationInterference.PartialWave; a pw::PartialWave is returned.
"""
function PartialWave()
    PartialWave( 0, 0., 0., PhotoRecombinationInterference.Channel[] )
end


# `Base.show(io::IO, pw::PhotoRecombinationInterference.PartialWave)`  ... prepares a proper printout of pw.
function Base.show(io::IO, pw::PhotoRecombinationInterference.PartialWave)
    println(io, "kappa:                    $(pw.kappa)  ")
    println(io, "energy:                   $(pw.energy)  ")
    println(io, "phase:                    $(pw.phase)  ")
    println(io, "channels:                 $(pw.channels)  ")
end


"""
`struct  PhotoRecombinationInterference.Pathway`
    ... defines a type for the coherent recombination of an electron with an initial ion level into a final level, for ONE energy of the
        incoming electron and through a set of intermediate resonances.

    + initialLevel             ::Level                  ... Initial-(state) level of the ion.
    + finalLevel               ::Level                  ... Final-(state) level of the recombined ion.
    + intermediateLevels       ::Array{Level,1}         ... The resonances that were actually included in the resonant term.
    + electronEnergy           ::Float64                ... Energy of the incoming electron [a.u.].
    + photonEnergy             ::Float64                ... Energy of the emitted photon [a.u.]; ONE value for both terms.
    + betaGamma2               ::Float64                ... beta^2 * gamma^2 of the incoming electron.
    + crossSection             ::EmProperty             ... Cross section from the SUMMED amplitude.
    + crossSectionRR           ::EmProperty             ... Cross section from the direct term alone.
    + crossSectionDR           ::EmProperty             ... Cross section from the resonant term alone.
    + interference             ::EmProperty             ... crossSection - crossSectionRR - crossSectionDR.
    + betaNu                   ::Array{EmPropertyC,1}   ... Anisotropy parameters beta_nu, nu = 2 and 4, from the summed amplitude.
    + linearPolarization       ::EmPropertyC            ... Linear polarization of the photon emitted at 90 degrees.
    + partialWaves             ::Array{PartialWave,1}   ... Partial waves, each with the channels it serves.
"""
struct  Pathway
    initialLevel               ::Level
    finalLevel                 ::Level
    intermediateLevels         ::Array{Level,1}
    electronEnergy             ::Float64
    photonEnergy               ::Float64
    betaGamma2                 ::Float64
    crossSection               ::EmProperty
    crossSectionRR             ::EmProperty
    crossSectionDR             ::EmProperty
    interference               ::EmProperty
    betaNu                     ::Array{EmPropertyC,1}
    linearPolarization         ::EmPropertyC
    partialWaves               ::Array{PartialWave,1}
end


"""
`PhotoRecombinationInterference.Pathway()`
    ... constructor for an `empty` instance of PhotoRecombinationInterference.Pathway; a pathway::Pathway is returned.
"""
function Pathway()
    Pathway( Level(), Level(), Level[], 0., 0., 1., EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.),
             EmPropertyC[], EmPropertyC(0.0im), PhotoRecombinationInterference.PartialWave[] )
end


# `Base.show(io::IO, pathway::PhotoRecombinationInterference.Pathway)`  ... prepares a proper printout of pathway.
function Base.show(io::IO, pathway::PhotoRecombinationInterference.Pathway)
    println(io, "initialLevel:             $(pathway.initialLevel)  ")
    println(io, "finalLevel:               $(pathway.finalLevel)  ")
    println(io, "intermediateLevels:       $(pathway.intermediateLevels)  ")
    println(io, "electronEnergy:           $(pathway.electronEnergy)  ")
    println(io, "photonEnergy:             $(pathway.photonEnergy)  ")
    println(io, "betaGamma2:               $(pathway.betaGamma2)  ")
    println(io, "crossSection:             $(pathway.crossSection)  ")
    println(io, "crossSectionRR:           $(pathway.crossSectionRR)  ")
    println(io, "crossSectionDR:           $(pathway.crossSectionDR)  ")
    println(io, "interference:             $(pathway.interference)  ")
    println(io, "betaNu:                   $(pathway.betaNu)  ")
    println(io, "linearPolarization:       $(pathway.linearPolarization)  ")
    println(io, "partialWaves:             $(pathway.partialWaves)  ")
end


"""
`PhotoRecombinationInterference.computeAmplitudesProperties(pathway::PhotoRecombinationInterference.Pathway,`
        `intermediateMultiplet::Multiplet, widths::Dict{Int64,Float64}, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,`
        `settings::PhotoRecombinationInterference.Settings)`
    ... to compute, for the given pathway, the direct and the resonant amplitude of every channel and to add them coherently; the three
        cross sections, the interference term, the anisotropy parameters and the linear polarization follow from them.  A
        newPathway::Pathway is returned.

        ONE continuum orbital is generated per partial wave, exactly as in `PhotoRecombination.computeAmplitudesProperties`, and the
        radiative amplitude of the resonant term is evaluated at the SAME photon energy as the direct term.
"""
function computeAmplitudesProperties(pathway::PhotoRecombinationInterference.Pathway, intermediateMultiplet::Multiplet,
                                     widths::Dict{Int64,Float64}, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                     settings::PhotoRecombinationInterference.Settings)
    contSettings = Continuum.Settings(false, nrContinuum)
    redFLevel    = Basics.generateLevelWithSymmetryReducedBasis(pathway.finalLevel, pathway.finalLevel.basis.subshells)
    newiLevel    = Basics.generateLevelWithSymmetryReducedBasis(pathway.initialLevel, redFLevel.basis.subshells)
    nFactor      = PhotoRecombinationInterference.dielectronicNormalization(pathway.initialLevel, pathway.finalLevel)
    totalEnergy  = pathway.initialLevel.energy + pathway.electronEnergy
    newPartialWaves = PhotoRecombinationInterference.PartialWave[]

    for  pw  in  pathway.partialWaves
        cSubshell = Subshell(101, pw.kappa)
        newfLevel = Basics.generateLevelWithExtraSubshell(cSubshell, redFLevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(pathway.electronEnergy, cSubshell, newiLevel, nm, grid, contSettings)
        newChannels = PhotoRecombinationInterference.Channel[]
        #
        for  ch  in  pw.channels
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newiLevel)
            ampsRR    = MultipoleAmplitude[];    ampsDR = MultipoleAmplitude[];    ampsTot = MultipoleAmplitude[]
            #
            for  mp  in  settings.multipoles
                # The direct, radiative-recombination term
                if  settings.includeRR
                    if  string(mp)[1] == 'E'
                        aC = PhotoRecombination.amplitude("photorecombination", mp, Basics.Coulomb,   pw.kappa, phase,
                                                          pathway.photonEnergy, newfLevel, newcLevel, grid)
                        aB = PhotoRecombination.amplitude("photorecombination", mp, Basics.Babushkin, pw.kappa, phase,
                                                          pathway.photonEnergy, newfLevel, newcLevel, grid)
                        ampRR = EmPropertyC(aC, aB)
                    else
                        aM = PhotoRecombination.amplitude("photorecombination", mp, Basics.Magnetic,  pw.kappa, phase,
                                                          pathway.photonEnergy, newfLevel, newcLevel, grid)
                        ampRR = EmPropertyC(aM)
                    end
                else
                    ampRR = EmPropertyC(0.0im)
                end
                # The resonant, dielectronic-recombination term
                ampDR = EmPropertyC(0.0im)
                if  settings.includeDR
                    for  dLevel  in  pathway.intermediateLevels
                        if  LevelSymmetry(dLevel.J, dLevel.parity) != ch.symmetry    continue    end
                        ampDR = ampDR + nFactor * PhotoRecombinationInterference.amplitudeDR(mp, pw.kappa, phase, dLevel, newfLevel,
                                            newcLevel, redFLevel, totalEnergy, pathway.photonEnergy,
                                            get(widths, dLevel.index, 0.), grid, settings)
                    end
                end
                push!( ampsRR,  MultipoleAmplitude(mp, ampRR) )
                push!( ampsDR,  MultipoleAmplitude(mp, ampDR) )
                push!( ampsTot, MultipoleAmplitude(mp, ampRR + ampDR) )
            end
            push!( newChannels, PhotoRecombinationInterference.Channel(ch.symmetry, ampsRR, ampsDR, ampsTot) )
        end
        push!( newPartialWaves, PhotoRecombinationInterference.PartialWave(pw.kappa, pathway.electronEnergy, phase, newChannels) )
    end
    #
    # The three cross sections, each through the SAME function of PhotoRecombination, so that they are commensurable by construction
    csTot = PhotoRecombination.computeCrossSectionForMultipoles( settings.multipoles,
                PhotoRecombinationInterference.toPhotoRecombinationLine(pathway, newPartialWaves, :total) )
    csRR  = PhotoRecombination.computeCrossSectionForMultipoles( settings.multipoles,
                PhotoRecombinationInterference.toPhotoRecombinationLine(pathway, newPartialWaves, :rr) )
    csDR  = PhotoRecombination.computeCrossSectionForMultipoles( settings.multipoles,
                PhotoRecombinationInterference.toPhotoRecombinationLine(pathway, newPartialWaves, :dr) )
    interference = csTot - csRR - csDR
    #
    betaNu = EmPropertyC[];    linPol = EmPropertyC(0.0im)
    totLine = PhotoRecombinationInterference.toPhotoRecombinationLine(pathway, newPartialWaves, :total)
    if  settings.calcAnisotropy
        for  nu in [2, 4]    push!( betaNu, PhotoRecombination.computeAnisotropyParameter(nu, totLine) )    end
    end
    if  settings.calcPolarization   linPol = PhotoRecombinationInterference.computeLinearPolarization(totLine, settings.multipoles)   end

    newPathway = PhotoRecombinationInterference.Pathway( pathway.initialLevel, pathway.finalLevel, pathway.intermediateLevels,
                     pathway.electronEnergy, pathway.photonEnergy, pathway.betaGamma2, csTot, csRR, csDR, interference,
                     betaNu, linPol, newPartialWaves )

    return( newPathway )
end


"""
`PhotoRecombinationInterference.amplitudeDR(mp::EmMultipole, kappa::Int64, phase::Float64, dLevel::Level, newfLevel::Level,`
        `newcLevel::Level, redFLevel::Level, totalEnergy::Float64, omega::Float64, gammaD::Float64, grid::Radial.Grid,`
        `settings::PhotoRecombinationInterference.Settings)`
    ... to compute the resonant contribution of ONE intermediate level d to the recombination amplitude of one channel, i.e.

            A_capture(d) * A_radiative(f <-- d) / (E_total - E_d + i Gamma_d / 2) ,

        without the overall normalization, which `dielectronicNormalization` supplies once per pathway.  The stabilizing amplitude is
        evaluated at the photon energy of the overall process rather than at E_d - E_f.  A value::EmPropertyC is returned.

        THE CAPTURE PHASE IS NOT conj(Auger), and this is the one place where the two source modules cannot simply be combined.  Both carry
        a factor for the continuum state, but not the same one: `PhotoRecombination.amplitude` multiplies by i^l exp(+i phi), so the
        state it uses is |c> = i^l exp(+i phi) |real orbital>, whereas `AutoIonization.amplitude` multiplies by i^l exp(-i phi) and returns
        i^l exp(-i phi) M_A with M_A REAL (its matrix is filled from real interaction strengths and real spin-angular coefficients).  The
        capture amplitude that belongs to PhotoRecombination's state is therefore

            <d|V|c> = i^l exp(+i phi) M_A = AutoIonization.amplitude * exp(2 i phi) ,

        while conj(AutoIonization.amplitude) = (-i)^l exp(+i phi) M_A differs from it by (-1)^l.  `ElectronCapture.amplitude` does take the
        conjugate, and is right to for its own purpose -- it only ever feeds |amplitude|^2, in which any such factor cancels -- but that
        makes it no guide here.  Since (-1)^l depends on the PARTIAL WAVE, using it would not even be an overall sign: it would reweight the
        partial waves against each other inside a coherent sum, and so change the interference itself.
"""
function amplitudeDR(mp::EmMultipole, kappa::Int64, phase::Float64, dLevel::Level, newfLevel::Level, newcLevel::Level,
                     redFLevel::Level, totalEnergy::Float64, omega::Float64, gammaD::Float64, grid::Radial.Grid,
                     settings::PhotoRecombinationInterference.Settings)
    denominator = Complex(totalEnergy - dLevel.energy, 0.5 * gammaD)
    if  abs(denominator) < 1.0e-30    return( EmPropertyC(0.0im) )    end
    #
    ampCapture = AutoIonization.amplitude(settings.augerOperator, kappa, phase, newcLevel, dLevel, grid; printout=false) *
                 exp( 2im * phase )
    if  ampCapture == ComplexF64(0.)    return( EmPropertyC(0.0im) )    end
    #
    if  string(mp)[1] == 'E'
        aC = PhotoEmission.amplitude(Emission(), mp, Basics.Coulomb,   omega, redFLevel, dLevel, grid, display=false, printout=false)
        aB = PhotoEmission.amplitude(Emission(), mp, Basics.Babushkin, omega, redFLevel, dLevel, grid, display=false, printout=false)
        ampRad = EmPropertyC(aC, aB)
    else
        aM = PhotoEmission.amplitude(Emission(), mp, Basics.Magnetic,  omega, redFLevel, dLevel, grid, display=false, printout=false)
        ampRad = EmPropertyC(aM)
    end

    return( (ampCapture / denominator) * ampRad )
end


"""
`PhotoRecombinationInterference.computeLinearPolarization(line::PhotoRecombination.Line, multipoles::Array{EmMultipole,1})`
    ... to compute the degree of linear polarization P = (I_parallel - I_perpendicular)/(I_parallel + I_perpendicular) of the photon
        emitted perpendicular to the electron beam, for an initially unpolarized ion; a value::EmPropertyC is returned that holds both
        gauges.

        For PURE ELECTRIC-DIPOLE radiation the polarization is not an independent quantity: the same alignment of the emitting system fixes
        both the angular distribution and the polarization, and at 90 degrees they are related rigorously by

            P = -3 beta_2 / (2 - beta_2) ,

        which is what is evaluated here from the beta_2 of `PhotoRecombination.computeAnisotropyParameter`, i.e. from the COHERENTLY SUMMED
        amplitudes.  Nothing is assumed about the relative size of the direct and the resonant term, so the interference enters the
        polarization exactly as it enters beta_2, which is the effect reported by Tong (2023).

        RESTRICTION, stated rather than hidden: this relation holds for E1 only.  Once a magnetic or a higher electric multipole
        contributes, the nu,2 component of the photon density matrix stops being proportional to the nu,0 component and a full density-
        matrix treatment is needed.  Rather than return a number that would silently be wrong, this function returns NaN in that case, and
        the display method prints it as such.
"""
function computeLinearPolarization(line::PhotoRecombination.Line, multipoles::Array{EmMultipole,1})
    if  multipoles != [E1]    return( EmPropertyC(Complex(NaN, 0.), Complex(NaN, 0.)) )    end
    beta2 = PhotoRecombination.computeAnisotropyParameter(2, line)
    polC  = -3.0 * beta2.Coulomb   / (2.0 - beta2.Coulomb)
    polB  = -3.0 * beta2.Babushkin / (2.0 - beta2.Babushkin)

    return( EmPropertyC(polC, polB) )
end


"""
`PhotoRecombinationInterference.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,`
        `initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,`
        `settings::PhotoRecombinationInterference.Settings; output::Bool=true)`
    ... to compute the coherent recombination amplitudes and all requested properties for the selected pathways and for every electron
        energy of the settings; this is the standard entry point of this module.  A list of
        pathways::Array{PhotoRecombinationInterference.Pathway,1} is returned if output = true, and nothing otherwise.
"""
function computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                         grid::Radial.Grid, settings::PhotoRecombinationInterference.Settings; output::Bool=true)
    println("")
    printstyled("PhotoRecombinationInterference.computePathways(): The computation of the coherent RR + DR amplitudes starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------------------ \n",
                color=:light_green)
    println("")
    widths   = PhotoRecombinationInterference.computeTotalWidths(intermediateMultiplet, initialMultiplet, finalMultiplet, nm, grid,
                                                                 settings)
    pathways = PhotoRecombinationInterference.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoRecombinationInterference.displayPathways(stdout, pathways, widths)    end
    maxEnergy = 0.;    for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    newPathways = PhotoRecombinationInterference.Pathway[]
    for  pathway in pathways
        push!( newPathways, PhotoRecombinationInterference.computeAmplitudesProperties(pathway, intermediateMultiplet, widths, nm, grid,
                                                                                      nrContinuum, settings) )
    end
    PhotoRecombinationInterference.displayResults(stdout, newPathways, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoRecombinationInterference.displayResults(iostream, newPathways, settings)    end

    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`PhotoRecombinationInterference.computeTotalWidths(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,`
        `finalMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoRecombinationInterference.Settings)`
    ... to determine the total width Gamma_d of every intermediate level, as the sum of its Auger rates into the levels of the initial
        multiplet and of its radiative rates into the levels of the final multiplet.  A widths::Dict{Int64,Float64} is returned, keyed on
        the index of the intermediate level.

        This width is built ONLY from the configurations that the computation names, and is therefore a LOWER bound: any decay channel
        outside them is missing, which makes the resonance too narrow and its peak too high.  The breakdown is printed for exactly that
        reason, and `settings.intermediateWidths` overrides the computed value wherever the user supplies one.
"""
function computeTotalWidths(intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, finalMultiplet::Multiplet, nm::Nuclear.Model,
                            grid::Radial.Grid, settings::PhotoRecombinationInterference.Settings)
    widths = Dict{Int64,Float64}()
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(); printBefore=false, maxKappa=settings.maxKappa,
                                            operator=settings.augerOperator)
    photonSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=settings.multipoles,
                                            gauges=[UseCoulomb, UseBabushkin], printBefore=false)
    augerLines  = AutoIonization.computeLines(initialMultiplet, intermediateMultiplet, nm, grid, augerSettings; output=true)
    photonLines = PhotoEmission.computeLines(finalMultiplet, intermediateMultiplet, grid, photonSettings; output=true)
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
`PhotoRecombinationInterference.determineChannels(finalLevel::Level, initialLevel::Level,`
        `settings::PhotoRecombinationInterference.Settings)`
    ... to determine the partial waves of the incoming electron and, for each of them, the total symmetries of the (ion + electron)
        scattering state that contribute; the list is that of `PhotoRecombination.determineChannels` and is obtained from it, so that the
        direct term is bookkept exactly as in that module.  An Array{PhotoRecombinationInterference.PartialWave,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombinationInterference.Settings)
    prSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(); multipoles=settings.multipoles, gauges=settings.gauges,
                                             electronEnergies=settings.electronEnergies, maxKappa=settings.maxKappa)
    prWaves    = PhotoRecombination.determineChannels(finalLevel, initialLevel, prSettings)
    partialWaves = PhotoRecombinationInterference.PartialWave[]
    for  pw in prWaves
        channels = PhotoRecombinationInterference.Channel[]
        for  ch in pw.channels
            push!( channels, PhotoRecombinationInterference.Channel(ch.symmetry, MultipoleAmplitude[], MultipoleAmplitude[],
                                                                    MultipoleAmplitude[]) )
        end
        push!( partialWaves, PhotoRecombinationInterference.PartialWave(pw.kappa, pw.energy, pw.phase, channels) )
    end

    return( partialWaves )
end


"""
`PhotoRecombinationInterference.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,`
        `initialMultiplet::Multiplet, settings::PhotoRecombinationInterference.Settings)`
    ... to determine the list of pathways for which the coherent recombination amplitude is to be computed, one per (initial level, final
        level, electron energy).  Each pathway carries the intermediate levels that SURVIVE the pathway selection for its own (i, f) pair,
        so that a single resonance -- or a chosen subset of them -- can be studied in isolation, which is what an interference study
        normally wants.  Which of the surviving levels then actually contributes is decided channel by channel through the symmetry
        condition J_d = J_t.  An Array{PhotoRecombinationInterference.Pathway,1} is returned.
"""
function determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                           settings::PhotoRecombinationInterference.Settings)
    pathways = PhotoRecombinationInterference.Pathway[]
    wc       = Defaults.getDefaults("speed of light: c")
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            dLevels = Level[]
            for  dLevel in intermediateMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, dLevel, fLevel, settings.pathwaySelection)   push!( dLevels, dLevel )   end
            end
            if  length(dLevels) == 0    continue    end
            for  en in settings.electronEnergies
                enAu = Defaults.convertUnits("energy: to atomic", en)
                if  enAu < 0.    continue    end
                omega  = enAu + iLevel.energy - fLevel.energy
                if  omega <= 0.    continue    end
                gamma  = 1.0 + enAu / wc^2
                beta   = sqrt( 1.0 - 1.0/gamma^2 )
                pWaves = PhotoRecombinationInterference.determineChannels(fLevel, iLevel, settings)
                push!( pathways, PhotoRecombinationInterference.Pathway(iLevel, fLevel, dLevels, enAu, omega,
                           beta^2*gamma^2, EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.),
                           EmPropertyC[], EmPropertyC(0.0im), pWaves) )
            end
        end
    end

    return( pathways )
end


"""
`PhotoRecombinationInterference.dielectronicNormalization(initialLevel::Level, finalLevel::Level)`
    ... to supply the constant N that puts the resonant amplitude on the same footing as the direct one, so that the two may be added.

        N is fixed by requiring that the resonant term ALONE, integrated over one isolated resonance, reproduce the standard dielectronic
        resonance strength S = (pi^2/E_res) (2J_d+1)/(2(2J_i+1)) Gamma_a A_r / Gamma_d.  Carrying JAC's own three normalizations through
        that condition -- Gamma_a = 2 pi Sum |A_Auger|^2, A_r = 8 pi alpha omega/(2J_d+1) Sum |A_rad|^2, and
        sigma_RR = (1/beta gamma^2) 8 pi^3 alpha^3 omega/(2J_f+1) Sum |A_RR|^2 -- everything cancels except

            N = sqrt( (2J_f+1) / (2J_i+1) ) .

        This value is DERIVED, and the derivation is not trusted on its own: branch b of example-Dvnew.jl measures the ratio of the
        integrated resonant cross section to the resonance strength of the DielectronicRecombination module and reports it.  A
        value::Float64 is returned.
"""
function dielectronicNormalization(initialLevel::Level, finalLevel::Level)
    wa = sqrt( (Basics.twice(finalLevel.J) + 1.0) / (Basics.twice(initialLevel.J) + 1.0) )

    return( wa )
end


"""
`PhotoRecombinationInterference.displayPathways(stream::IO, pathways::Array{PhotoRecombinationInterference.Pathway,1},`
        `widths::Dict{Int64,Float64})`
    ... to display a list of the pathways that have been selected for the computations, together with the electron and photon energies and
        the number of partial waves.  A neat table is printed but nothing is returned otherwise.
"""
function displayPathways(stream::IO, pathways::Array{PhotoRecombinationInterference.Pathway,1}, widths::Dict{Int64,Float64})
    nx = 118
    println(stream, " ")
    println(stream, "  Selected photorecombination-interference pathways:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2)
    sa = sa * TableStrings.center(16, "E_electron"; na=2)
    sa = sa * TableStrings.center(16, "E_photon"; na=2)
    sa = sa * TableStrings.center(14, "no. of kappa"; na=2)
    sa = sa * TableStrings.center(16, "no. of resonances"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  pathway in pathways
        isym = LevelSymmetry(pathway.initialLevel.J, pathway.initialLevel.parity)
        fsym = LevelSymmetry(pathway.finalLevel.J,   pathway.finalLevel.parity)
        sa   = "  " * TableStrings.center(18, TableStrings.levels_if(pathway.initialLevel.index, pathway.finalLevel.index); na=2)
        sa   = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=2)
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "  "
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.photonEnergy))   * "  "
        sa   = sa * TableStrings.center(14, string(length(pathway.partialWaves)); na=2)
        sa   = sa * TableStrings.center(16, string(length(pathway.intermediateLevels)); na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`PhotoRecombinationInterference.displayResults(stream::IO, pathways::Array{PhotoRecombinationInterference.Pathway,1},`
        `settings::PhotoRecombinationInterference.Settings)`
    ... to display the three cross sections, the interference term and -- where they were requested -- the anisotropy parameter beta_2 and
        the linear polarization, for every pathway and electron energy.  A neat table is printed but nothing is returned otherwise.
"""
function displayResults(stream::IO, pathways::Array{PhotoRecombinationInterference.Pathway,1},
                        settings::PhotoRecombinationInterference.Settings)
    nx = 136
    println(stream, " ")
    println(stream, "  Coherent RR + DR recombination:  cross sections [a.u.] and the interference term, Coulomb gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "i-level-f"; na=1)
    sa = sa * TableStrings.center(15, "E_electron"; na=1)
    sa = sa * TableStrings.center(15, "sigma total"; na=1)
    sa = sa * TableStrings.center(15, "sigma RR"; na=1)
    sa = sa * TableStrings.center(15, "sigma DR"; na=1)
    sa = sa * TableStrings.center(15, "interference"; na=1)
    sa = sa * TableStrings.center(13, "beta_2"; na=1)
    sa = sa * TableStrings.center(13, "lin. polar."; na=1)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  pathway in pathways
        sa = "  " * TableStrings.center(16, TableStrings.levels_if(pathway.initialLevel.index, pathway.finalLevel.index); na=1)
        sa = sa * @sprintf("%14.6e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSection.Coulomb)   * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSectionRR.Coulomb) * " "
        sa = sa * @sprintf("%14.6e", pathway.crossSectionDR.Coulomb) * " "
        sa = sa * @sprintf("%14.6e", pathway.interference.Coulomb)   * " "
        if  length(pathway.betaNu) > 0   sa = sa * @sprintf("%12.5e", real(pathway.betaNu[1].Coulomb)) * " "
        else                             sa = sa * TableStrings.center(13, "--"; na=1)
        end
        if  settings.calcPolarization    sa = sa * @sprintf("%12.5e", real(pathway.linearPolarization.Coulomb)) * " "
        else                             sa = sa * TableStrings.center(13, "--"; na=1)
        end
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end


"""
`PhotoRecombinationInterference.toPhotoRecombinationLine(pathway::PhotoRecombinationInterference.Pathway,`
        `partialWaves::Array{PhotoRecombinationInterference.PartialWave,1}, which::Symbol)`
    ... to package the amplitudes of the given partial waves into a `PhotoRecombination.Line`, so that the cross section and the anisotropy
        parameters can be taken from that module unchanged.  `which` selects which set of amplitudes travels: `:total` the coherent sum,
        `:rr` the direct term alone and `:dr` the resonant term alone.  A line::PhotoRecombination.Line is returned.

        This is what makes the limit includeDR = false reproduce `PhotoRecombination` by CONSTRUCTION: in that limit the very same numbers
        are handed to the very same functions.
"""
function toPhotoRecombinationLine(pathway::PhotoRecombinationInterference.Pathway,
                                  partialWaves::Array{PhotoRecombinationInterference.PartialWave,1}, which::Symbol)
    prWaves = PhotoRecombination.PartialWave[]
    for  pw in partialWaves
        prChannels = PhotoRecombination.Channel[]
        for  ch in pw.channels
            if      which == :total   amps = ch.amplitudes
            elseif  which == :rr      amps = ch.amplitudesRR
            elseif  which == :dr      amps = ch.amplitudesDR
            else    error("stop a: unknown selector $which")
            end
            push!( prChannels, PhotoRecombination.Channel(ch.symmetry, amps) )
        end
        push!( prWaves, PhotoRecombination.PartialWave(pw.kappa, pw.energy, pw.phase, prChannels) )
    end

    return( PhotoRecombination.Line(pathway.initialLevel, pathway.finalLevel, pathway.electronEnergy, pathway.photonEnergy,
                                    pathway.betaGamma2, 0., EmProperty(0., 0.), prWaves) )
end

end # module
