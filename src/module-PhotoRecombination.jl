
"""
`module  JAC.PhotoRecombination`  
... a submodel of JAC that contains all methods for computing photo-recomnbination, i.e. radiative recombination 
    and radiative electron capture properties between some initial and final-state multiplets.
"""
module PhotoRecombination


using Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..HydrogenicIon, ..InteractionStrength,
                ..ManyElectron, ..PhotoEmission, ..Radial, ..Nuclear, ..TableStrings

"""
`struct  PhotoRecombination.Settings  <:  AbstractProcessSettings` ... defines a type for the details and parameters of computing photo recombination lines.

    + multipoles          ::Array{EmMultipole}  ... Multipoles of the radiation field that are to be included.
    + gauges              ::Array{UseGauge}     ... Gauges to be included into the computations.
    + electronEnergies    ::Array{Float64,1}    ... List of electron energies [in default units].
    + ionEnergies         ::Array{Float64,1}    ... List of ion energies [in MeV/u].
    + useIonEnergies      ::Bool                ... Make use of ion energies in [MeV/u] to obtain the electron energies.
    + calcTotalCs         ::Bool                ... True, if the total cross sections is to be calculated/displayed for all initial levels.
    + calcAnisotropy      ::Bool                ... True, if the overall anisotropy is to be calculated.
    + calcTensors         ::Bool                ... True, if the statistical tensors are to be calculated and 
                                                    false otherwise.
    + printBefore         ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + maxKappa            ::Int64               ... Maximum kappa value of partial waves to be included.
    + lineSelection       ::LineSelection       ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles            ::Array{EmMultipole}
    gauges                ::Array{UseGauge}
    electronEnergies      ::Array{Float64,1} 
    ionEnergies           ::Array{Float64,1}
    useIonEnergies        ::Bool
    calcTotalCs           ::Bool
    calcAnisotropy        ::Bool
    calcTensors           ::Bool 
    printBefore           ::Bool 
    maxKappa              ::Int64 
    lineSelection         ::LineSelection 
end 


"""
`PhotoRecombination.Settings(set::PhotoRecombination..Settings;`

        multipoles=..,          gauges=..,              electronEnergies=..,          ionEnergies=..,     
        useIonEnergies=..,      calcTotalCs..,          calcAnisotropy=..,            calcTensors=..,             
        printBefore=..,         maxKappa=..,            lineSelection=..)
                    
    ... constructor for modifying the given PhotoRecombination..Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoRecombination.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    electronEnergies::Union{Nothing,Array{Float64,1}}=nothing,              ionEnergies::Union{Nothing,Array{Float64,1}}=nothing,       
    useIonEnergies::Union{Nothing,Bool}=nothing,                            calcTotalCs::Union{Nothing,Bool}=nothing, 
    calcAnisotropy::Union{Nothing,Bool}=nothing,                            calcTensors::Union{Nothing,Bool}=nothing,   
    printBefore::Union{Nothing,Bool}=nothing,                               maxKappa::Union{Nothing,Int64}=nothing, 
    lineSelection::Union{Nothing,LineSelection}=nothing)  
    
    if  isnothing(multipoles)          multipolesx        = set.multipoles        else  multipolesx        = multipoles         end 
    if  isnothing(gauges)              gaugesx            = set.gauges            else  gaugesx            = gauges             end 
    if  isnothing(electronEnergies)    electronEnergiesx  = set.electronEnergies  else  electronEnergiesx  = electronEnergies   end 
    if  isnothing(ionEnergies)         ionEnergiesx       = set.ionEnergies       else  ionEnergiesx       = ionEnergies        end 
    if  isnothing(useIonEnergies)      useIonEnergiesx    = set.useIonEnergies    else  useIonEnergiesx    = useIonEnergies     end 
    if  isnothing(calcTotalCs)         calcTotalCsx       = set.calcTotalCs       else  calcTotalCsx       = calcTotalCs        end 
    if  isnothing(calcAnisotropy)      calcAnisotropyx    = set.calcAnisotropy    else  calcAnisotropyx    = calcAnisotropy     end 
    if  isnothing(calcTensors)         calcTensorsx       = set.calcTensors       else  calcTensorsx       = calcTensors        end 
    if  isnothing(printBefore)         printBeforex       = set.printBefore       else  printBeforex       = printBefore        end 
    if  isnothing(maxKappa)            maxKappax          = set.maxKappa          else  maxKappax          = maxKappa           end 
    if  isnothing(lineSelection)       lineSelectionx     = set.lineSelection     else  lineSelectionx     = lineSelection      end 

    Settings( multipolesx, gaugesx, electronEnergiesx, ionEnergiesx, useIonEnergiesx, calcTotalCsx, calcAnisotropyx, 
              calcTensorsx, printBeforex, maxKappax, lineSelectionx)
end


"""
`PhotoRecombination.Settings()`  ... constructor for the default values of photo recombination line computations
"""
function Settings()
    Settings(EmMultipole[], UseGauge[], Float64[], Float64[], false, false, false, false, false, 0, LineSelection() )
end


# `Base.show(io::IO, settings::PhotoRecombination.Settings)`  
#		... prepares a proper printout of the variable settings::PhotoRecombination.Settings.
function Base.show(io::IO, settings::PhotoRecombination.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "electronEnergies:         $(settings.electronEnergies)  ")
    println(io, "ionEnergies:              $(settings.ionEnergies)  ")
    println(io, "useIonEnergies:           $(settings.useIonEnergies)  ")
    println(io, "calcTotalCs:              $(settings.calcTotalCs)  ")
    println(io, "calcAnisotropy:           $(settings.calcAnisotropy)  ")
    println(io, "calcTensors:              $(settings.calcTensors)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end








#####################################################################################################################
## THE PHYSICAL CHANNEL, built BESIDE the flat one above -- the second module to carry it, after PhotoIonization.
##
## A PhotoRecombination.Channel labels (multipole, gauge, kappa, symmetry).  Of those four, only kappa and the
## total symmetry describe the scattering STATE: a multipole is a term in the expansion of the interaction
## OPERATOR and a gauge is a representation of that same operator.  Here the state is described once and the
## operator hangs off it -- a partial wave carries kappa, the electron energy and the scattering phase, which
## depend on nothing else; a channel carries one total J^parity; and an amplitude carries BOTH gauges in a
## Basics.EmPropertyC, so that the gauge is never a label anywhere.
##
## What this module tests that PhotoIonization could not: its consumers are different.  The anisotropy
## parameter needs PAIRS of amplitudes from DIFFERENT partial waves, and needed the same gauge guard written
## three times -- which is where all three of them had the wrong operand until 91a91c0.
#####################################################################################################################




"""
`struct  PhotoRecombination.Channel`
    ... ONE asymptotic scattering state: the initial ion plus a free electron, coupled to a total symmetry.

    + symmetry       ::LevelSymmetry                            ... total J^parity of the scattering state.
    + amplitudes     ::Array{MultipoleAmplitude,1}        ... one entry per contributing multipole.
"""
struct  Channel
    symmetry         ::LevelSymmetry
    amplitudes       ::Array{MultipoleAmplitude,1}
end


"""
`struct  PhotoRecombination.PartialWave`
    ... ONE partial wave of the incoming free electron, and the channels it serves.

    + kappa          ::Int64                            ... partial wave of the free electron.
    + energy         ::Float64                          ... energy of the free electron.
    + phase          ::Float64                          ... scattering phase; a property of (energy, kappa).
    + channels       ::Array{Channel,1}           ... the total symmetries this partial wave serves.

        The radial orbital and the phase belong HERE, not to a channel: they do not depend on the total
        symmetry, which is precisely why one kappa can serve several of them. computeAmplitudesProperties
        needed a per-kappa dictionary to exploit that; here it follows from the shape.
"""
struct  PartialWave
    kappa            ::Int64
    energy           ::Float64
    phase            ::Float64
    channels         ::Array{Channel,1}
end


"""
`struct  PhotoRecombination.Line`
    ... as PhotoRecombination.Line, but carrying partial waves instead of flat channels; every other field is
        the same, so that retiring the flat form later is a deletion and a rename rather than a rewrite.

    + initialLevel   ::Level                            ... initial-(state) level
    + finalLevel     ::Level                            ... final-(state) level
    + electronEnergy ::Float64                          ... Energy of the (incoming free) electron.
    + photonEnergy   ::Float64                          ... Energy of the emitted photon.
    + betaGamma2     ::Float64                          ... beta^2 * gamma^2.
    + weight         ::Float64                          ... weight of line in the integration over electron energies.
    + crossSection   ::EmProperty                       ... Cross section for this electron capture.
    + partialWaves   ::Array{PartialWave,1}       ... partial waves, each with the channels it serves.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    photonEnergy     ::Float64
    betaGamma2       ::Float64
    weight           ::Float64
    crossSection     ::EmProperty
    partialWaves     ::Array{PartialWave,1}
end


# `Base.show(io::IO, line::PhotoRecombination.Line)`  ... prepares a proper printout of the variable line::PhotoRecombination.Line.
function Base.show(io::IO, line::PhotoRecombination.Line)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "electronEnergy:    $(line.electronEnergy)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "betaGamma2:        $(line.betaGamma2)  ")
    println(io, "weight:            $(line.weight)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "partialWaves:      $(line.partialWaves)  ")
end




"""
`PhotoRecombination.amplitude(kind::String, mp::EmMultipole, gauge::EmGauge, kappa::Int64, phase::Float64,
                                energy::Float64, finalLevel::Level, continuumLevel::Level, grid::Radial.Grid)`

    The gauge is passed EXPLICITLY rather than read off a channel: it is a representation of the interaction
    operator, and since the retirement of the flat channel form it is no longer a label carried by any state.
    kappa and the phase likewise belong to the partial wave, and are handed over as the numbers they are.
    ... to compute the kind = (photorecombination) amplitude  
        < alpha_f J_f || O^(photorecombination) || (alpha_i J_i, epsilon kappa) J_t>  due to the electron-photon 
        interaction for the given final and continuum level, the partial wave of the outgoing electron as well as 
        the given multipole and gauge. A value::ComplexF64 is returned.
"""
function amplitude(kind::String, mp::EmMultipole, gauge::EmGauge, kappa::Int64, phase::Float64, energy::Float64,
                    finalLevel::Level, continuumLevel::Level, grid::Radial.Grid)
    if      kind in [ "photorecombination"]
    #-----------------------------------
        amplitude = PhotoEmission.amplitude(Emission(), mp, gauge, energy, finalLevel, continuumLevel, grid)
        amplitude = im^Basics.subshell_l(Subshell(101, kappa)) * exp( im*phase ) * amplitude
    else    error("stop b")
    end
    
    return( amplitude )
end


"""
`PhotoRecombination.checkConsistentMultiplets(finalMultiplet::Multiplet, initialMultiplet::Multiplet)`  
    ... to check that the given initial- and final-state levels and multiplets are consistent to each other and
        to avoid later problems with the computations. An error message is issued if an inconsistency occurs,
        and nothing is returned otherwise.
"""
function  checkConsistentMultiplets(finalMultiplet::Multiplet, initialMultiplet::Multiplet)
    initialSubshells      = initialMultiplet.levels[1].basis.subshells;             ni = length(initialSubshells)
    finalSubshells        = finalMultiplet.levels[1].basis.subshells
    
    if initialSubshells[1:end] == finalSubshells[1:ni]
    else
        error("\nThe order of subshells must be equal for the initial- and final states. \n" *
                "However, the initial states can have less subshells; this limitation arises from the angular coefficients.")
    end
        
    return( nothing )
end













"""
`PhotoRecombination.crossSectionStobbe(energy::Float64, Z::Float64)`  
    ... to evaluate the (non-relativistic) Stobbe cross section for the RR of a free electron with energy into the 1s state of 
        initially bare ions with nuclear charge Z; an cs::Float64 is returned.
"""
function crossSectionStobbe(energy::Float64, Z::Float64)
    eta = sqrt( Z^2 / 2. / energy)
    cs  = 2^8 * pi^2 * Defaults.getDefaults("alpha")^3 / 3.
    cs  = cs * eta^6 * exp(-4 * eta * atan(1/eta)) / (1 - exp(-2 * pi * eta)) / (eta^2 + 1)^2
    return( cs )
end



"""
`PhotoRecombination.crossSectionKramers(energy::Float64, Z::Float64, nLowUp::Tuple{Int64,Int64})`  
    ... to evaluate the (non-relativistic) Kramers cross section for the RR of a free electron with energy into all shells with
        n = n_Low ... n_Up of initially bare ions with nuclear charge Z; an cs::Float64 is returned.
"""
function crossSectionKramers(energy::Float64, Z::Float64, nLowUp::Tuple{Int64,Int64})
    eta = sqrt( Z^2 / 2. / energy);     wa = 32 * pi * Defaults.getDefaults("alpha")^3 / 3. / sqrt(3.)
    cs  = 0.;       for  n=nLowUp[1]:nLowUp[2]  cs = cs + wa * eta^4 / n / (eta^2 + n^2)    end
    return( cs )
end


"""
`PhotoRecombination.crossSectionKramersTotal(energy::Float64, Z::Float64)`  
    ... to evaluate the (non-relativistic) Kramers cross section for the RR of a free electron with energy into any shell
        of initially bare ions with nuclear charge Z; an cs::Float64 is returned.
"""
function crossSectionKramersTotal(energy::Float64, Z::Float64)
    eta = sqrt( Z^2 / 2. / energy);     wa = 16 * pi * Defaults.getDefaults("alpha")^3 / 3. / sqrt(3.)
    cs  = wa * eta^2 * log(1 + eta^2)
    return( cs )
end


"""
`PhotoRecombination.crossSectionBellTotal(energy::Float64, Z::Float64)`  
    ... to evaluate the (non-relativistic) total Bell & Bell cross section for the RR of a free electron with energy 
        into any shell of initially bare ions with nuclear charge Z; an cs::Float64 is returned.
"""
function crossSectionBellTotal(energy::Float64, Z::Float64)
    eta = sqrt( Z^2 / 2. / energy);     wa = 32 * pi * Defaults.getDefaults("alpha")^3 / 3. / sqrt(3.)
    cs = log(eta) + 0.1492 + 0.5250 / eta^(2/3)
    cs  = wa * eta^2 * cs
    return( cs )
end











"""
`PhotoRecombination.displayRateCoefficients(stream::IO, isym::LevelSymmetry, temperatures::Array{Float64,1}, alphaRR::Array{EmProperty,1})`  
    ... to print all rate coefficients for the selected temperatures in neat tables, 
        though nothing is returned otherwise.
"""
function  displayRateCoefficients(stream::IO, isym::LevelSymmetry, temperatures::Array{Float64,1}, alphaRR::Array{EmProperty,1})
    #
    ntemps = length(temperatures)
    nx = 54 + 17 * min(ntemps, 7)
    println(stream, " ")
    println(stream, "  Rate coefficients [cm^3/s] for initial level with symmetry J^P = $isym:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = sb = "                                   "
    for  nt = 1:min(ntemps, 7)
        sa = sa * TableStrings.center(14, "T = " * @sprintf("%.2e", temperatures[nt]); na=3);       
        sb = sb * TableStrings.center(14, "[K]"; na=3)
    end
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    println(stream, "  ")
    sa = "  alpha^RR (T, i; Coulomb gauge):    " 
    sb = "  alpha^RR (T, i; Babushkin gauge):  " 
    for  nt = 1:min(ntemps, 7) 
        sa = sa * @sprintf("%.4e", alphaRR[nt].Coulomb)    * "       "
        sb = sb * @sprintf("%.4e", alphaRR[nt].Babushkin)  * "       "
    end
    println(stream, sa);    println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`PhotoRecombination.plasmaRateKotelnikov(Te::Float64, Z::Float64)`  
    ... to evaluate the RR plasma rate coefficient at temperature Te and for a Maxwellian distribution of electron
        energies and an initially bare ions with nuclear charge Z; an alpha::Float64 is returned.
"""
function plasmaRateKotelnikov(Te::Float64, Z::Float64)
    wx    = 2*Te / Z^2
    alpha = 8.414 * Defaults.getDefaults("alpha")^3 * Z
    alpha = alpha * (log(1.0 + 1.0/wx) + 3.499)
    alpha = alpha / (sqrt(wx)  +  0.6517 * wx  +  0.2138 * wx^1.5)
    return( alpha )
end


"""
`PhotoRecombination.plasmaRateKotelnikov_1s(Te::Float64, Z::Float64)`  
    ... to evaluate the RR plasma rate coefficient at temperature Te and for a Maxwellian distribution of electron
        energies and an initially bare ions with nuclear charge Z; here only the capture into 1s is taken into account.
        An alpha::Float64 is returned.
"""
function plasmaRateKotelnikov_1s(Te::Float64, Z::Float64)
    wx    = 2*Te / Z^2
    alpha = 8.414 * Defaults.getDefaults("alpha")^3 * Z
    alpha = alpha / (sqrt(wx)  +  0.3593 * wx^(7/6)  +  0.1471 * wx^1.5)
    return( alpha )
end


"""
`PhotoRecombination.plasmaRatePartialSeaton(Te::Float64, Z::Float64, n::Int64)`  
    ... to evaluate the RR plasma rate coefficient at temperature Te and for a Maxwellian distribution of electron
        energies and an initially bare ions with nuclear charge Z.
        An alpha::Float64 is returned.
"""
function plasmaRatePartialSeaton(Te::Float64, Z::Float64, n::Int64)
    wx    = Z^2 / (2*n^2 * Te)
    alpha = 64 * Defaults.getDefaults("alpha")^3 / 3. * sqrt(pi/3.) * wx^(3/2) * exp(wx)
    alpha = 11 * alpha  # fudge factor
    # Formal integral is still missing.
    return( alpha )
end


"""
`PhotoRecombination.plasmaRateSeaton(Te::Float64, Z::Float64)`  
    ... to evaluate the RR plasma rate coefficient at temperature Te and for a Maxwellian distribution of electron
        energies and an initially bare ions with nuclear charge Z; this formula is valid for 2*Te / Z^2 << 1.
        An alpha::Float64 is returned.
"""
function plasmaRateSeaton(Te::Float64, Z::Float64)
    wx    = 2*Te / Z^2
    alpha = 32 * sqrt(pi) * Defaults.getDefaults("alpha")^3 * Z / (3. * sqrt(3.))
    alpha = alpha * sqrt(1/wx) * (log(1/wx) + 0.8576 + 0.9380 * (1/wx)^(-1/3))
    return( alpha )
end


#####################################################################################################################
## STAGE A of the physical channel: the core.  Everything below is ADDITIVE -- nothing above is altered, and the
## flat path is bitwise unchanged, which work/diag-photorec-gate.jl checks.
#####################################################################################################################


"""
`PhotoRecombination.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)`
    ... as PhotoRecombination.determineChannels, but returning the partial waves of the physical form; an
        Array{PhotoRecombination.PartialWave,1} is returned with all amplitudes still zero.

        The SELECTION RULES are exactly those of determineChannels above -- the same
        AngularMomentum.allowedMultipoleSymmetries(symf, mp) reaching the total symmetry from the FINAL
        (recombined) level, then allowedKappaSymmetries(symi, symt) coupling the free electron to the INITIAL
        ion -- and the same `abs(kappa) > settings.maxKappa` restriction. Only the nesting is inverted: there,
        multipole and gauge are outermost and the physical pair (kappa, symt) is regenerated for each of them;
        here the distinct pairs are collected first and every multipole reaching a pair is attached to it.
        The gauge is not iterated over at all.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoRecombination.Settings)
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    ## (1) Collect, for every physical pair (kappa, symt), the multipoles that can reach it.
    mpsFor = Dict{Tuple{Int64,LevelSymmetry}, Array{EmMultipole,1}}()
    order  = Tuple{Int64,LevelSymmetry}[]                    ## to keep a reproducible sequence
    for  mp in settings.multipoles
        for  symt in AngularMomentum.allowedMultipoleSymmetries(symf, mp)
            for  kappa in AngularMomentum.allowedKappaSymmetries(symi, symt)
                if  abs(kappa) > settings.maxKappa      continue    end
                key = (kappa, symt)
                if  haskey(mpsFor, key)   push!(mpsFor[key], mp)
                else                      mpsFor[key] = EmMultipole[mp];   push!(order, key)
                end
            end
        end
    end
    ## (2) Group the pairs by kappa: one partial wave per kappa, carrying the symmetries it serves.
    kappas = Int64[];   for (kappa, symt) in order    if !(kappa in kappas)   push!(kappas, kappa)   end    end
    partialWaves = PartialWave[]
    for  kappa in kappas
        channels = Channel[]
        for  (ka, symt) in order
            ka == kappa   ||   continue
            amps = MultipoleAmplitude[]
            for  mp in mpsFor[(ka, symt)]
                push!(amps, MultipoleAmplitude(mp, EmPropertyC(Complex(0.), Complex(0.))))
            end
            push!(channels, Channel(symt, amps))
        end
        push!(partialWaves, PartialWave(kappa, 0., 0., channels))
    end
    return( partialWaves )
end


"""
`PhotoRecombination.computeAmplitudesProperties(line::PhotoRecombination.Line, nm::Nuclear.Model,
        grid::Radial.Grid, nrContinuum::Int64, settings::PhotoRecombination.Settings;
        nuclearPot::Union{Nothing,Radial.Potential}=nothing, primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... as PhotoRecombination.computeAmplitudesProperties, but on partial waves; a Line with all
        amplitudes and the cross section evaluated is returned. The two keywords carry quantities that are
        CONSTANT for a whole computation, exactly as in the flat version.

        ONE continuum orbital and ONE extended final level per PARTIAL WAVE -- structurally, not by a
        dictionary. The flat function needs `cOrbitals`, `cPhases` and `fLevels` keyed on the continuum
        subshell precisely because its channels repeat kappa once per (multipole, gauge); here the loop nesting
        says it, and there is nothing to key on.

        An electric multipole is evaluated twice, once per gauge, into one EmPropertyC; a magnetic multipole
        once, into an EmPropertyC with equal components.
"""
function computeAmplitudesProperties(line::PhotoRecombination.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                           nrContinuum::Int64, settings::PhotoRecombination.Settings;
                                           nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                           primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    contSettings = Continuum.Settings(false, nrContinuum)
    redFLevel    = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, line.finalLevel.basis.subshells)
    newiLevel    = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, redFLevel.basis.subshells)
    newPartialWaves = PhotoRecombination.PartialWave[]

    for  pw in line.partialWaves
        cSubshell = Subshell(101, pw.kappa)
        newfLevel = Basics.generateLevelWithExtraSubshell(cSubshell, redFLevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, cSubshell, newiLevel, nm, grid,
                                                            contSettings; nuclearPot=nuclearPot, primitives=primitives)
        newChannels = PhotoRecombination.Channel[]
        for  ch in pw.channels
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newiLevel)
            newAmps   = MultipoleAmplitude[]
            for  ma in ch.amplitudes
                mp = ma.multipole
                if  string(mp)[1] == 'E'
                    ampC = PhotoRecombination.amplitude("photorecombination", mp, Basics.Coulomb, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                    ampB = PhotoRecombination.amplitude("photorecombination", mp, Basics.Babushkin, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                    push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampC, ampB)))
                else
                    ampM = PhotoRecombination.amplitude("photorecombination", mp, Basics.Magnetic, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                    push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampM)))
                end
            end
            push!(newChannels, PhotoRecombination.Channel(ch.symmetry, newAmps))
        end
        push!(newPartialWaves, PhotoRecombination.PartialWave(pw.kappa, line.electronEnergy, phase, newChannels))
    end
    newLine      = PhotoRecombination.Line(line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                                 line.betaGamma2, line.weight, line.crossSection, newPartialWaves)
    crossSection = PhotoRecombination.computeCrossSectionForMultipoles(settings.multipoles, newLine)
    newLine      = PhotoRecombination.Line(line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                                 line.betaGamma2, line.weight, crossSection, newPartialWaves)
    return( newLine )
end


"""
`PhotoRecombination.computeCrossSectionForMultipoles(multipoles::Array{EmMultipole,1},
                                                           line::PhotoRecombination.Line)`
    ... as PhotoRecombination.computeCrossSectionForMultipoles: the cross section from the amplitudes of the
        given line, restricted to the given multipoles; a cs::EmProperty is returned.

        The three-way `if` on the gauge label is GONE. abs2 of an EmPropertyC is an EmProperty, so the sum is
        gauge-paired by construction, and a magnetic amplitude -- equal components -- enters both gauge sums
        by itself, which is exactly what the flat form's third branch had to say by hand.

        The summation is INCOHERENT over multipoles as well as over channels, exactly as in the flat version.
        The new form would make a coherent multipole sum easy to write; that would be a physics change and is
        deliberately not made here.
"""
function computeCrossSectionForMultipoles(multipoles::Array{EmMultipole,1}, line::PhotoRecombination.Line)
    cs = EmProperty(0., 0.)
    for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
        if  !(ma.multipole  in  multipoles)     continue    end
        cs = cs + abs2(ma.amplitude)
    end
    csFactor = 8 * pi^3 * Defaults.getDefaults("alpha")^3 * line.photonEnergy / (Basics.twice(line.finalLevel.J)+1)
    return( 1.0 / line.betaGamma2 * (csFactor * cs) )
end


"""
`PhotoRecombination.computeAnisotropyParameter(nu::Int64, line::PhotoRecombination.Line)`
    ... as PhotoRecombination.computeAnisotropyParameter: the anisotropy parameter beta_nu of the emitted
        photons for an initially unpolarized ion. An EmPropertyC is returned, holding BOTH gauges.

        NOTE ON THE SIGNATURE: no `gauge` argument. In the flat form a gauge is a channel label that has to be
        selected for, so the display calls the function eight times -- twice per nu. Here one call returns both.

        THIS IS THE FUNCTION THE NEW FORM HAS TO JUSTIFY IN THIS MODULE, and the case for it is not
        hypothetical. The flat version needs the same gauge guard written THREE times -- once for the
        normalization wn, once for the outer channel, once for the inner -- and until 91a91c0 all three tested
        the requested gauge instead of the channel's own, so every magnetic multipole was silently dropped and
        the odd-nu parameters came out identically zero. Here the product conj(amp) * ampp is componentwise,
        so Coulomb meets Coulomb and Babushkin meets Babushkin, and a magnetic amplitude has equal components
        and so enters both. There is no guard to write and none to get wrong.

        Pairs are formed ACROSS partial waves: kappa and kappa' are independent, and that interference is the
        physical content of the anisotropy. The parity rule and every angular factor are reproduced verbatim.
"""
function computeAnisotropyParameter(nu::Int64, line::PhotoRecombination.Line)
    Ji = line.initialLevel.J;    Jf = line.finalLevel.J
    ## Flatten once: one entry per (partial wave, total symmetry, multipole), each carrying both gauges.
    entries = Tuple{Int64,LevelSymmetry,EmMultipole,EmPropertyC}[]
    for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
        push!(entries, (pw.kappa, ch.symmetry, ma.multipole, ma.amplitude))
    end
    wn = EmProperty(0., 0.)
    for  (kapa, syma, mpa, ampa)  in entries      wn = wn + abs2(ampa)      end
    #
    wa = EmPropertyC(0.0im)
    for  (kapa, syma, mpa, ampa)  in entries
        J = syma.J;    L = mpa.L;    if  mpa.electric   p = 1   else    p = 0   end
        j = AngularMomentum.kappa_j(kapa);    l = AngularMomentum.kappa_l(kapa)
        #
        for  (kapb, symb, mpb, ampb)  in entries
            Jp = symb.J;   Lp = mpb.L;   if  mpb.electric   pp = 1   else    pp = 0   end
            jp = AngularMomentum.kappa_j(kapb);     lp = AngularMomentum.kappa_l(kapb)
            #
            if  1 + (-1)^(L + p + Lp + pp - nu) == 0    continue    end
            wa = wa + (1.0im)^(L + p - Lp - pp) * AngularMomentum.phaseFactor([Ji, -1, AngularJ64(1//2), -1, Jf]) *
                                sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), l, lp, j, jp, J, Jp]) ) *
                                AngularMomentum.ClebschGordan( l, AngularM64(0), lp, AngularM64(0),  AngularJ64(nu),  AngularM64(0)) *
                                AngularMomentum.ClebschGordan( AngularJ64(L), AngularM64(1), AngularJ64(Lp), AngularM64(-1),
                                                            AngularJ64(nu), AngularM64(0)) *
                                AngularMomentum.Wigner_6j(J, Jp, AngularJ64(nu), AngularJ64(Lp), AngularJ64(L), Jf) *
                                AngularMomentum.Wigner_6j(J, Jp, AngularJ64(nu), jp, j, Ji) *
                                AngularMomentum.Wigner_6j(j, jp, AngularJ64(nu), lp, l, AngularJ64(1//2)) *
                                (conj(ampa) * ampb)
        end
    end
    return( EmPropertyC(- 0.5 * wa.Coulomb / wn.Coulomb, - 0.5 * wa.Babushkin / wn.Babushkin) )
end


#####################################################################################################################
## STAGE B: the bridge back to the flat form, and the displays.
#####################################################################################################################






"""
`PhotoRecombination.displayLines(stream::IO, lines::Array{PhotoRecombination.Line,1})`
`PhotoRecombination.displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings)`
    ... the display layer.  Both are the tables that were printed before the flat channel form was retired,
        reading partial waves instead of channels; the printed output is unchanged, and
        work/diag-retire-baseline.jl holds the pre-retirement text it is diffed against.

        TWO PLACES NEEDED MORE THAN A RENAME.  displayLines prints one row per (multipole, GAUGE), which is no
        longer how an amplitude is stored, so that expansion is written out where the table is built -- the
        gauge is a property of the PRESENTATION here, not of the channel.  And displayResults called
        computeAnisotropyParameter eight times per line, twice per rank; one call per rank now returns both
        gauges, so there are four.
"""
function  displayLines(stream::IO, lines::Array{PhotoRecombination.Line,1})
    nx = 181
    println(stream, " ")
    println(stream, "  Selected photorecombination lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(10, "Energy_if"; na=2);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "Energy e_r"; na=1);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(10, "omega"; na=5);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center( 7, "beta^2*"; na=2)
    sb = sb * TableStrings.center( 7, "gamma^2"; na=2)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa = sa * @sprintf("%.2e", line.betaGamma2)                                                   * "  "
        ## The table shows one row per (multipole, GAUGE), which the amplitudes no longer carry as a label:
        ## an electric multipole is one amplitude holding two gauges, a magnetic one has no gauge freedom.
        ## Expanding it here keeps this table exactly as it always was; the gauge is a property of the
        ## PRESENTATION, not of the channel.
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            if  string(ma.multipole)[1] == 'E'
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Coulomb,   ch.symmetry) )
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Babushkin, ch.symmetry) )
            else
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Magnetic,  ch.symmetry) )
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        if  length(wa) > 0  sb = sa * wa[1]   else    sb = sa    end;    println(stream,  sb )  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


function  displayResults(stream::IO, lines::Array{PhotoRecombination.Line,1}, settings::PhotoRecombination.Settings)
    if  settings.useIonEnergies   nx = 170    else  nx = 146   end
    println(stream, " ")
    println(stream, "  Photorecombination cross sections:")
    println(stream, "  ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(16, "i--J^P--f"   ; na=3);                       sb = sb * TableStrings.hBlank(19)
    sa = sa * TableStrings.center(12, "i--Energy--f"; na=4)               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "omega"     ; na=4)             
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "Energy e_r"; na=3)             
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    if  settings.useIonEnergies
        sa = sa * TableStrings.center(12, "Ion energies"; na=3)             
        sb = sb * TableStrings.center(12, "[MeV/u]";      na=3)
    end
    sa = sa * TableStrings.center( 7, "beta^2*"; na=4)
    sb = sb * TableStrings.center( 7, "gamma^2"; na=4)
    sa = sa * TableStrings.flushleft(16, "Multipoles"; na=1);                      sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                            TableStrings.inUnits("cross section"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(17, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, fsym); na=3)
        en = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        if  settings.useIonEnergies 
            # E^(electron) [keV] = E^(projectile) [MeV/u] / 1.8228885
            enIon = 1.8228885 * Defaults.convertUnits("energy: from atomic to eV", line.electronEnergy) / 1000.
            sa = sa * @sprintf("%.6e", enIon) * "    "
        end
        sa = sa * @sprintf("%.2e", line.betaGamma2)                                                   * "   "
        multipoles = EmMultipole[]
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            multipoles = push!( multipoles, ma.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "                              "
        sa = sa * TableStrings.flushleft(16, mpString[1:16];  na=2)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb))     * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin))   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    if  settings.calcTotalCs 
        nx = 143
        println(stream, " ")
        println(stream, "  Total photorecombination cross sections for the intial levels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level"   ; na=0);                       sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(12, "i--J^P "   ; na=3);                       sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(12, "i--Energy "; na=3)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=5)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_r"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=17);                      sb = sb * TableStrings.hBlank(15)
        sa = sa * TableStrings.flushleft(57, "Cou -- Total cross section -- Bab"; na=4)  
        sb = sb * TableStrings.center(57, TableStrings.inUnits("cross section") * "          " * 
                                            TableStrings.inUnits("cross section");  na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #
        index = Int64[];    eEnergy = Float64[]
        for  line in lines
            if  line.initialLevel.index  in  index  &&  line.electronEnergy  in  eEnergy   continue  end
            push!(index, line.initialLevel.index);    push!(eEnergy, line.electronEnergy) 
            sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
            sa = sa * TableStrings.center(17, TableStrings.level(line.initialLevel.index))
            sa = sa * TableStrings.center(17, string(isym))
            en = line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
                multipoles = push!( multipoles, ma.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "                              "
            sa = sa * TableStrings.flushleft(26, mpString[1:26];  na=3)
            tcs = Basics.EmProperty(0.)
            for  lineb  in  lines
                if  lineb.initialLevel.index == line.initialLevel.index  &&
                    lineb.electronEnergy     == line.electronEnergy      tcs = tcs + lineb.crossSection   end
            end
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", tcs.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", tcs.Babushkin))   * "    "
            println(stream, sa);   sa = TableStrings.hBlank(102)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcAnisotropy  
        nx = 133
        println(stream, " ")
        println(stream, "  Anisotropy angular parameters beta_nu of the emitted photons for initially unpolarized ions:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(16, "i--J^P--f"   ; na=3);                       sb = sb * TableStrings.hBlank(19)
        sa = sa * TableStrings.center(12, "i--Energy--f"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_r"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=5);                         sb = sb * TableStrings.hBlank(15)
        sa = sa * TableStrings.flushleft(57, "beta's"; na=4)  
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = " ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(17, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.initialLevel.energy - line.finalLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
                multipoles = push!( multipoles, ma.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
            be = PhotoRecombination.computeAnisotropyParameter(1, line)
            sa = sa * @sprintf("%.4e", real(be.Coulomb))   * " (1, C);  "
            sa = sa * @sprintf("%.4e", real(be.Babushkin)) * " (1, B);  "
            println(stream, sa);   sa = TableStrings.hBlank(102)
            be = PhotoRecombination.computeAnisotropyParameter(2, line)
            sa = sa * @sprintf("%.4e", real(be.Coulomb))   * " (2, C);  "
            sa = sa * @sprintf("%.4e", real(be.Babushkin)) * " (2, B);  "
            println(stream, sa);   sa = TableStrings.hBlank(102)
            be = PhotoRecombination.computeAnisotropyParameter(3, line)
            sa = sa * @sprintf("%.4e", real(be.Coulomb))   * " (3, C);  "
            sa = sa * @sprintf("%.4e", real(be.Babushkin)) * " (3, B);  "
            println(stream, sa);   sa = TableStrings.hBlank(102)
            be = PhotoRecombination.computeAnisotropyParameter(4, line)
            sa = sa * @sprintf("%.4e", real(be.Coulomb))   * " (4, C);  "
            sa = sa * @sprintf("%.4e", real(be.Babushkin)) * " (4, B);  "
            println(stream, sa);   sa = TableStrings.hBlank(102)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTensors   
        println(stream, " ")
        println(stream, "  Reduced statistical tensors of the recombined ion ... not yet implemented !!")
        println(stream, " ")
    end
    #
    return( nothing )
end


#####################################################################################################################
## STAGE C: the two drivers, so that a whole computation can be run either way and compared end to end.
#####################################################################################################################


"""
`PhotoRecombination.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                         settings::PhotoRecombination.Settings)`
    ... as PhotoRecombination.determineLines, but producing Line with partial waves; an
        Array{PhotoRecombination.Line,1} is returned.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoRecombination.Settings)
    lines = PhotoRecombination.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                #
                electronEnergies = Float64[]
                if  settings.useIonEnergies
                    for en in settings.ionEnergies
                        # E^(electron) [keV] = E^(projectile) [MeV/u] / 1.8228885
                        en_au = 1000. / 1.8228885 * Defaults.convertUnits("energy: from eV to atomic", en);      push!(electronEnergies, en_au)
                    end
                else
                    for en in settings.electronEnergies
                        en_au = Defaults.convertUnits("energy: to atomic", en);      push!(electronEnergies, en_au)
                    end
                end
                #
                for  en in electronEnergies
                    betaGamma2 = 1.0
                    wc         = Defaults.getDefaults("speed of light: c")
                    Gamma      = 1.0 + en / wc^2
                    beta       = sqrt( 1.0 - 1.0/Gamma^2)
                    betaGamma2 = beta^2 * Gamma^2
                    #
                    if  en < 0    continue   end
                    omega        = en + iLevel.energy - fLevel.energy
                    partialWaves = PhotoRecombination.determineChannels(fLevel, iLevel, settings)
                    push!( lines, PhotoRecombination.Line(iLevel, fLevel, en, omega, betaGamma2, 0.,
                                                                EmProperty(0., 0.), partialWaves) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                       grid::Radial.Grid, settings::PhotoRecombination.Settings; output::Bool=true)`
    ... as PhotoRecombination.computeLines, but on the physical form; a list of
        lines::Array{PhotoRecombination.Line,1} is returned.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::PhotoRecombination.Settings; output::Bool=true)
    println("")
    printstyled("PhotoRecombination.computeLines(): The computation of photo-recombination properties starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    PhotoRecombination.checkConsistentMultiplets(finalMultiplet, initialMultiplet)
    #
    lines = PhotoRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoRecombination.displayLines(stdout, lines)    end
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    ## Both are constants of the whole computation and are READ-ONLY below, so they are safely shared.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    newLines = Vector{PhotoRecombination.Line}(undef, length(lines))
    Threads.@threads for  i  in  eachindex(lines)
        newLines[i] = PhotoRecombination.computeAmplitudesProperties(lines[i], nm, grid, nrContinuum, settings;
                                                                           nuclearPot=nuclearPot, primitives=primitives)
    end
    PhotoRecombination.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoRecombination.displayResults(iostream, newLines, settings)    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoRecombination.computeLinesWithContinuumOrbital(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
        nm::Nuclear.Model, grid::Radial.Grid, cOrbitals::Dict{Subshell, Orbital}, energyGrid::Radial.GridGL,
        settings::PhotoRecombination.Settings; output::Bool=true)`
    ... as PhotoRecombination.computeLinesWithContinuumOrbital, but on the physical form; a list of
        lines::Array{PhotoRecombination.Line,1} is returned. This is the driver that
        Cascade.computeSteps(::Cascade.RadiativeRecombinationScheme, ...) uses.

        MIRRORED AS IT STANDS, including two things that differ from computeLines and are NOT changed here:
        the phase is set to 0. rather than obtained from the orbital, and the cross section is summed inline
        over ALL multipoles instead of going through computeCrossSectionForMultipoles, so `settings.multipoles`
        does not filter it. Both are properties of the flat function; this is a translation.

        The continuum orbital is looked up ONCE per partial wave rather than once per channel, which is the
        same saving the physical form gives everywhere else -- here it also means the energy-grid index search
        runs once per partial wave instead of once per (multipole, gauge, kappa).
"""
function computeLinesWithContinuumOrbital(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                                grid::Radial.Grid, cOrbitals::Dict{Subshell, Orbital}, energyGrid::Radial.GridGL,
                                                settings::PhotoRecombination.Settings; output::Bool=true)
    PhotoRecombination.checkConsistentMultiplets(finalMultiplet, initialMultiplet);     ie = 0
    #
    lines = PhotoRecombination.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoRecombination.displayLines(stdout, lines)    end
    newLines = PhotoRecombination.Line[]
    for  (i,line)  in  enumerate(lines)
        if  rem(i,500) == 0    println("> PhotoRecombination line $i:  ... calculated ")    end
        redFLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, line.finalLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, redFLevel.basis.subshells)
        en        = line.electronEnergy
        # Find index of en in energyGrid.t with given accuracy; terminate if nothing is found
        ie = 0
        for  it = 1:length(energyGrid.t)   if   abs( (energyGrid.t[it]-en)/en ) < 0.0001   ie = it;   break   end   end
        if  ie == 0   stop("a")     end
        #
        newPartialWaves = PhotoRecombination.PartialWave[];    cs = EmProperty(0., 0.)
        for  pw in line.partialWaves
            cSubsh    = Subshell(100+ie, pw.kappa)
            newfLevel = Basics.generateLevelWithExtraSubshell(cSubsh, redFLevel)
            cOrbital  = cOrbitals[cSubsh];      phase = 0.
            newChannels = PhotoRecombination.Channel[]
            for  ch in pw.channels
                newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newiLevel)
                newAmps   = MultipoleAmplitude[]
                for  ma in ch.amplitudes
                    mp = ma.multipole
                    if  string(mp)[1] == 'E'
                        ampC = PhotoRecombination.amplitude("photorecombination", mp, Basics.Coulomb, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                        ampB = PhotoRecombination.amplitude("photorecombination", mp, Basics.Babushkin, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                        amp  = EmPropertyC(ampC, ampB)
                    else
                        ampM = PhotoRecombination.amplitude("photorecombination", mp, Basics.Magnetic, pw.kappa, phase, line.photonEnergy, newfLevel, newcLevel, grid)
                        amp  = EmPropertyC(ampM)
                    end
                    cs = cs + abs2(amp)                     ## over ALL multipoles, as the flat version does
                    push!(newAmps, MultipoleAmplitude(mp, amp))
                end
                push!(newChannels, PhotoRecombination.Channel(ch.symmetry, newAmps))
            end
            push!(newPartialWaves, PhotoRecombination.PartialWave(pw.kappa, en, phase, newChannels))
        end
        csFactor     = 8 * pi^3 * Defaults.getDefaults("alpha")^3 * line.photonEnergy / (Basics.twice(line.finalLevel.J)+1)
        crossSection = 1.0 / line.betaGamma2 * (csFactor * cs)
        push!( newLines, PhotoRecombination.Line(line.initialLevel, line.finalLevel, line.electronEnergy,
                                                       line.photonEnergy, line.betaGamma2, energyGrid.wt[ie],
                                                       crossSection, newPartialWaves) )
    end
    PhotoRecombination.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoRecombination.displayResults(iostream, newLines, settings)    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end

end # module
