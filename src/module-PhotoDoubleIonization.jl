
"""
`module  JAC.PhotoDoubleIonization`  
... a submodel of JAC that computes photo-double ionization, i.e. the absorption of ONE photon with TWO electrons emitted into the
    continuum. The final state is built as a pair of outgoing partial waves: the first electron couples to the final ion to give an
    intermediate symmetry, the second couples to that to give the total symmetry of the scattering state, and a multipole must reach that
    total symmetry from the initial level. The amplitude is second order, with both time orderings of the knock-out (TS1) mechanism, and
    the intermediate states are levels of the (N-1)-electron ion with ONE continuum partial wave added.

    STATUS, 15-Aug-2026: POSTPONED BY THE MAINTAINER. Do NOT pick this module up as a task without saying so
    first. It was rebuilt on this date from the parked quasi-shell version, and the rebuild is UNFINISHED:
    the structure works and passes several internal tests, but the ABSOLUTE SCALE IS WRONG BY THREE TO FOUR
    ORDERS OF MAGNITUDE and no number it produces means anything yet.

    DO NOT MISREAD THE STATE. That the module now runs end to end, returns finite numbers and passes its
    symmetry tests is NOT evidence that it is nearly right. The three self-tests below are necessary and
    cheap, and they were all passing while the result was still wrong by 10^3-10^4. In particular, a
    beautifully symmetric energy-sharing distribution says only that the pair set is exchange-closed.

    WHAT WORKS, and was verified on He (Z=2, 1s^2 -> He^2+ + 2e) at 200/400/800 eV:
      * The five-layer structure Line -> Sharing -> PartialWavePair -> Channel -> MultipoleAmplitude.
        For He at E1 every pair correctly serves the single total symmetry 1-, and the eight pairs at
        maxKappa=2 are the relativistic decomposition of the textbook eps-s eps-p and eps-p eps-d channels.
      * The double-ionization threshold comes out at 77.77 eV against 79.0 eV measured.
      * The energy-sharing distribution is mirror-symmetric to 1e-10 .. 1e-14, and EACH TIME ORDERING is
        separately symmetric, which is the stronger test: a defect in one can otherwise be masked by the other.
      * The distribution is U-shaped and deepens with photon energy, as it must far above threshold.
      * The total cross section reproduces the weighted sum of the differential ones to 0.00%.

    WHAT IS WRONG. sigma(2+)/sigma(+) comes out 0.001% at 200 eV rising to 0.006% at 800 eV, against a
    measured 3-4% near 200 eV falling to the asymptotic 1.68%. So the result is ~3600x too small AND trends
    the wrong way. The single-ionization cross section computed from the SAME orbitals, normalization and
    prefactor lands within 15% of the literature, so the deficit is specific to the second-order
    double-ionization amplitude, not to the machinery underneath it.

    THE REFERENCE THAT SETTLES THE NORMALIZATION is Kornberg & Miraglia, Phys. Rev. A 48, 3714 (1993),
    "Double photoionization of helium: use of a correlated two-electron continuum wave function"
    (examples/papers/1993.pra-kornberg-pdi-helium.pdf). Their Eq. (2) reads

        d^5 sigma(2+) / (d eps1 d Omega1 d Omega2)  =  4 pi^2 alpha a0^2 k1 k2 C^(G) |e.T_fi^(G)(k1,k2)|^2

    with C^(L) = E_gamma, C^(V) = 1/E_gamma, C^(A) = E_gamma^(-3). THREE FACTORS ARE MISSING HERE:
      1. the momentum product k1 k2, which this module does not apply at all. It is O(1) -- 1.90 at the end
         sharing of the 200 eV case and 4.49 at the midpoint -- so it CANNOT close the gap, but it flattens
         the U by a factor 2.4 and is therefore a shape prediction to test.
      2. the gauge-dependent prefactor C^(G): length and velocity differ by E_gamma^2, whereas this module
         applies ONE gauge-independent prefactor to both Coulomb and Babushkin.
      3. a factor 1/2: "the correctly normalized total cross section is one half of the integral in the
         shake-off region". This is NOT the same double counting that Basics.determineEnergySharings handles.
    Their final state is normalized to a delta function in MOMENTUM space while JAC normalizes its continuum
    orbitals per unit ENERGY; since dE = k dk that conversion generates k factors of its own and must be done
    deliberately, not assumed. Doing it may absorb part or all of item 1.

    THE OTHER SUSPECT is the intermediate space, which is severely truncated: the numbers above used four
    Green levels, maxKappa=1 and a 5-point Gauss-Legendre integral over the WHOLE excess-energy range, for an
    integrand peaked near threshold. The rise of the ratio with photon energy is the signature of exactly
    that. A convergence scan was started and did not finish; it is the first thing to run on resuming.

    WHAT IS NOT A PROBLEM. The Babushkin/Coulomb ratio of ~2.7 is NOT evidence of a defect: Kornberg &
    Miraglia report a large length-velocity discrepancy at all energies, with the two forms separated by
    about an order of magnitude over much of the range and velocity the trustworthy one at high energy.
    Do not spend time hunting it.

    QUANTITATIVE TARGETS for the next attempt, both from that paper: Fig. 2 gives d sigma / d eps1 in Mb/eV
    at 120, 150 and 200 eV and at 1, 2 and 3 keV; Fig. 3 gives sigma(2+)/sigma(+) against photon energy with
    the 1.68% asymptote marked. Note their total cross section cost about 20 h on a 10-Mflop machine, so this
    is an expensive calculation even when done correctly.

    The error() in computeLines stays until a number is earned, and is to be removed in the same commit that
    earns it -- not before.
"""
module PhotoDoubleIonization


using Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron, ..PhotoEmission,
                ..TableStrings

"""
`struct  PhotoDoubleIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing photo-double ionization lines.

    + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
    + photonEnergies          ::Array{Float64,1}        ... List of photon energies [in user-selected units].  
    + electronEnergyShift     ::Float64                 ... An overall energy shift of the two emitted electrons [in user-selected units];
        it corrects the total electron energy for a shift of the i-f transition energy and thereby separates the energies at which the
        continuum orbitals are generated from the sharing coordinates.
    + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
    + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
    + calcDifferentialCs      ::Bool                    ... True, if the energy-differential cs are to be calculated and false otherwise.  
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + eeInteraction           ::AbstractEeInteraction   ... Type of the electron-electron interaction in the second-order treatment.
    + gMultiplet              ::Multiplet               ... Mean-field multiplet of intermediate levels in the computations.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    photonEnergies            ::Array{Float64,1} 
    electronEnergyShift       ::Float64 
    NoEnergySharings          ::Int64         
    maxKappa                  ::Int64 
    calcDifferentialCs        ::Bool 
    printBefore               ::Bool
    lineSelection             ::LineSelection
    eeInteraction             ::AbstractEeInteraction
    gMultiplet                ::Multiplet
end 


"""
`PhotoDoubleIonization.Settings()`  ... constructor for the default values of PhotoDoubleIonization line computations
"""
function Settings()
    Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], 0., 0, 0, false, false, 
                LineSelection(), CoulombInteraction(), Multiplet())
end


"""
`PhotoDoubleIonization.Settings(set::PhotoDoubleIonization.Settings;`

        multipoles=..,          gauges=..,                  photonEnergies=..,          
        electronEnergyShift=.., NoEnergySharings=..,     
        maxKappa=..,            calcDifferentialCs..,       printBefore=..,             lineSelection=..,       
        eeInteraction=..,       gMultiplet=..)
                    
    ... constructor for modifying the given PhotoDoubleIonization.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoDoubleIonization.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergyShift::Union{Nothing,Float64}=nothing,
    NoEnergySharings::Union{Nothing,Int64}=nothing,       
    maxKappa::Union{Nothing,Int64}=nothing,                                 calcDifferentialCs::Union{Nothing,Bool}=nothing,      
    printBefore::Union{Nothing,Bool}=nothing,                               lineSelection::Union{Nothing,LineSelection}=nothing, 
    eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,            gMultiplet::Union{Nothing,Multiplet}=nothing)  
    
    if  isnothing(multipoles)           multipolesx         = set.multipoles         else  multipolesx         = multipoles          end 
    if  isnothing(gauges)               gaugesx             = set.gauges             else  gaugesx             = gauges              end 
    if  isnothing(photonEnergies)       photonEnergiesx     = set.photonEnergies     else  photonEnergiesx     = photonEnergies      end 
    if  isnothing(electronEnergyShift)  electronEnergyShiftx= set.electronEnergyShift else electronEnergyShiftx= electronEnergyShift end 
    if  isnothing(NoEnergySharings)     NoEnergySharingsx   = set.NoEnergySharings   else  NoEnergySharingsx   = NoEnergySharings    end 
    if  isnothing(maxKappa)             maxKappax           = set.maxKappa           else  maxKappax           = maxKappa            end 
    if  isnothing(calcDifferentialCs)   calcDifferentialCsx = set.calcDifferentialCs else  calcDifferentialCsx = calcDifferentialCs  end 
    if  isnothing(printBefore)          printBeforex        = set.printBefore        else  printBeforex        = printBefore         end 
    if  isnothing(lineSelection)        lineSelectionx      = set.lineSelection      else  lineSelectionx      = lineSelection       end 
    if  isnothing(eeInteraction)        eeInteractionx      = set.eeInteraction      else  eeInteractionx      = eeInteraction       end 
    if  isnothing(gMultiplet)           gMultipletx         = set.gMultiplet         else  gMultipletx         = gMultiplet          end 

    Settings( multipolesx, gaugesx, photonEnergiesx, electronEnergyShiftx, NoEnergySharingsx, maxKappax, calcDifferentialCsx, 
                printBeforex, lineSelectionx, eeInteractionx, gMultipletx)
end


# `Base.show(io::IO, settings::PhotoDoubleIonization.Settings)`  
#   ... prepares a proper printout of the variable settings::PhotoDoubleIonization.Settings.
function Base.show(io::IO, settings::PhotoDoubleIonization.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "electronEnergyShift:      $(settings.electronEnergyShift)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "NoEnergySharings:         $(settings.NoEnergySharings)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "calcDifferentialCs:       $(settings.calcDifferentialCs)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "eeInteraction:            $(settings.eeInteraction)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
end


"""
`struct  PhotoDoubleIonization.Channel`
    ... defines a type for one total symmetry of the complete final state, i.e. of the final ion together with both emitted electrons, and
        collects the multipole amplitudes that reach this symmetry from the initial level.

    + symmetry       ::LevelSymmetry                    ... total angular momentum and parity of the complete scattering state.
    + amplitudes     ::Array{MultipoleAmplitude,1}      ... one entry per contributing multipole, each holding BOTH gauges.
"""
struct  Channel
    symmetry         ::LevelSymmetry
    amplitudes       ::Array{MultipoleAmplitude,1}
end


# `Base.show(io::IO, channel::PhotoDoubleIonization.Channel)`  ... prepares a proper printout of channel::PhotoDoubleIonization.Channel.
function Base.show(io::IO, channel::PhotoDoubleIonization.Channel)
    println(io, "symmetry:               $(channel.symmetry)  ")
    println(io, "amplitudes:             $(channel.amplitudes)  ")
end


"""
`struct  PhotoDoubleIonization.PartialWavePair`
    ... defines a type for one pair of outgoing partial waves, i.e. for the two electrons that leave the ion together. The two electrons are
        coupled in sequence: the first couples to the final ion to give xSymmetry, and the second couples to that intermediate symmetry to
        give the total symmetry of a channel. One pair generally serves several total symmetries, which is why the channels are held here
        rather than the pair being repeated for each of them; the two continuum orbitals are therefore generated once per pair.

    + kappa1         ::Int64                ... partial wave of the first emitted electron.
    + energy1        ::Float64              ... energy at which the orbital of the first electron is generated.
    + phase1         ::Float64              ... scattering phase of the first partial wave.
    + xSymmetry      ::LevelSymmetry        ... intermediate symmetry of (final ion + electron 1).
    + kappa2         ::Int64                ... partial wave of the second emitted electron.
    + energy2        ::Float64              ... energy at which the orbital of the second electron is generated.
    + phase2         ::Float64              ... scattering phase of the second partial wave.
    + channels       ::Array{PhotoDoubleIonization.Channel,1}   ... total symmetries that this pair of partial waves serves.
"""
struct  PartialWavePair
    kappa1           ::Int64
    energy1          ::Float64
    phase1           ::Float64
    xSymmetry        ::LevelSymmetry
    kappa2           ::Int64
    energy2          ::Float64
    phase2           ::Float64
    channels         ::Array{PhotoDoubleIonization.Channel,1}
end


# `Base.show(io::IO, pair::PhotoDoubleIonization.PartialWavePair)`  ... prepares a proper printout of pair::PhotoDoubleIonization.PartialWavePair.
function Base.show(io::IO, pair::PhotoDoubleIonization.PartialWavePair)
    println(io, "kappa1:                 $(pair.kappa1)  ")
    println(io, "energy1:                $(pair.energy1)  ")
    println(io, "phase1:                 $(pair.phase1)  ")
    println(io, "xSymmetry:              $(pair.xSymmetry)  ")
    println(io, "kappa2:                 $(pair.kappa2)  ")
    println(io, "energy2:                $(pair.energy2)  ")
    println(io, "phase2:                 $(pair.phase2)  ")
    println(io, "channels:               $(pair.channels)  ")
end


"""
`struct  PhotoDoubleIonization.Sharing`
    ... defines a type for one division of the excess energy between the two emitted electrons. The list of sharings of a line IS the
        energy-sharing distribution, which is the primary observable of this process, and each sharing owns two quantities that belong
        nowhere else: the quadrature weight of this sharing point, and the differential cross section summed over all pairs and channels
        at this point.

        The sharing coordinates epsilon1, epsilon2 are what Basics.determineEnergySharings produces and define the quadrature point. They
        are NOT the same as the energies at which the orbitals are generated, which additionally carry settings.electronEnergyShift and
        are held by the partial-wave pair.

    + omega          ::Float64         ... energy of the incident photon.
    + epsilon1       ::Float64         ... sharing coordinate of (free) electron 1.
    + epsilon2       ::Float64         ... sharing coordinate of (free) electron 2.
    + weight         ::Float64         ... Gauss-Legendre weight of this sharing for energy-integrated quantities.
    + differentialCs ::EmProperty      ... energy-differential cross section dsigma/depsilon1 at this sharing.
    + partialWavePairs ::Array{PhotoDoubleIonization.PartialWavePair,1}   ... pairs of partial waves at this sharing.
"""
struct  Sharing
    omega            ::Float64
    epsilon1         ::Float64
    epsilon2         ::Float64
    weight           ::Float64
    differentialCs   ::EmProperty
    partialWavePairs ::Array{PhotoDoubleIonization.PartialWavePair,1}
end


# `Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of sharing::PhotoDoubleIonization.Sharing.
function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)
    println(io, "omega:                  $(sharing.omega)  ")
    println(io, "epsilon1:               $(sharing.epsilon1)  ")
    println(io, "epsilon2:               $(sharing.epsilon2)  ")
    println(io, "weight:                 $(sharing.weight)  ")
    println(io, "differentialCs:         $(sharing.differentialCs)  ")
    println(io, "partialWavePairs:       $(sharing.partialWavePairs)  ")
end


"""
`struct  PhotoDoubleIonization.Line`
    ... defines a type for a photo-double ionization line between an initial and a final level, for one photon energy, in which a single
        photon is absorbed and TWO electrons are emitted into the continuum. The cross section is the sharing-integrated total.

    + initialLevel   ::Level                  ... initial-(state) level.
    + finalLevel     ::Level                  ... final-(state) level of the doubly-ionized ion.
    + photonEnergy   ::Float64                ... energy of the absorbed photon.
    + crossSection   ::EmProperty             ... total cross section, integrated over all energy sharings.
    + sharings       ::Array{PhotoDoubleIonization.Sharing,1}  ... the energy-sharing distribution of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    photonEnergy     ::Float64
    crossSection     ::EmProperty
    sharings         ::Array{PhotoDoubleIonization.Sharing,1}
end


"""
`PhotoDoubleIonization.Line()`  ... constructor for an empty PhotoDoubleIonization line; a line::PhotoDoubleIonization.Line is returned.
"""
function Line()
    Line(Level(), Level(), 0., EmProperty(0., 0.), PhotoDoubleIonization.Sharing[])
end


# `Base.show(io::IO, line::PhotoDoubleIonization.Line)`  ... prepares a proper printout of the variable line::PhotoDoubleIonization.Line.
function Base.show(io::IO, line::PhotoDoubleIonization.Line)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "sharings:          $(line.sharings)  ")
end


"""
`PhotoDoubleIonization.amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
        nLevels::Array{Level,1}, nWeights::Array{Float64,1}, grid::Radial.Grid; display::Bool=false, printout::Bool=false)`
    ... computes the second-order photo-double ionization amplitude

                <alpha_f J_f || O^(Mp, absorption) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>  / (E_i - E_n)
            +   <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, absorption) || alpha_i J_i>  / (E_i - omega - E_n)

        for the interaction with a photon of frequency omega, multipolarity Mp and the given gauge; a value::ComplexF64 is returned.
        The first term is the ordering in which the two electrons interact first and the photon is absorbed afterwards, the second the
        ordering in which the photon is absorbed first. Together they are the knock-out (TS1) mechanism.

        The intermediate levels are handed in rather than taken from a multiplet, because they are not bound states: each of them is a
        level of the (N-1)-electron ion with ONE continuum partial wave added, and such a state has one continuum electron exactly as the
        final state has two. A Green multiplet of bound levels cannot contribute at all, since it differs from the final state in two
        orbitals and the one-body photon operator connects only one.

        Because those partial waves are normalized PER ENERGY INTERVAL, so are the intermediate levels, and the sum over them is really an
        integral over the intermediate electron energy. nWeights carries the quadrature weight of each intermediate level and must have the
        same length as nLevels; passing ones would silently turn the integral back into a plain sum.

        The intermediate level is in general multi-configurational, so the two matrix elements of each ordering run over two INDEPENDENT
        CSF indices t and tp; a single index with mc[t]^2 would keep only the diagonal terms. Each ordering has its own symmetry condition
        and the two are mutually exclusive whenever symi != symf, which is why neither may skip the loop over intermediate levels.
"""
function amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                   nLevels::Array{Level,1}, nWeights::Array{Float64,1}, grid::Radial.Grid; display::Bool=false, printout::Bool=false)
    length(nLevels) == length(nWeights)  ||  error("nLevels and nWeights must have the same length.")
    if  length(nLevels) == 0    return( ComplexF64(0.) )    end

    # Always ensure the same subshell list for all initial, intermediate and final levels.  Each intermediate
    # level carries a continuum subshell of its own kappa, so every one of them has to enter the merge.
    subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
    for  nLevel in nLevels   subshells = Basics.merge(subshells, nLevel.basis.subshells)   end
    iLevel    = Level(initialLevel, subshells)
    fLevel    = Level(finalLevel,   subshells)

    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity);    eni = iLevel.energy

    if  printout   printstyled("Compute photo-double $(Mp) ionization amplitude for the transition " *
                               "[$(iLevel.index)-$(fLevel.index)] ... ", color=:light_green)    end
    amplitude = ComplexF64(0.)

    for  r = 1:nf
        symr = LevelSymmetry(fLevel.basis.csfs[r].J, fLevel.basis.csfs[r].parity);      if  symr != symf    continue    end
        for  s = 1:ni
            syms = LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity);  if  syms != symi    continue    end
            for  (k, nLevel) in enumerate(nLevels)
                nLevel = Level(nLevel, subshells)
                symn   = LevelSymmetry(nLevel.J, nLevel.parity);    enn = nLevel.energy;    wn = nWeights[k]
                nn     = length(nLevel.basis.csfs)

                #   Compute <alpha_f J_f || O^(Mp, kind) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
                if  symn == symi
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        OMp = ManyElectron.matrixElement_Mab(Mp, gauge, omega, fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            Vee       = ManyElectron.matrixElement_Vee(CoulombInteraction(), nLevel.basis, tp, iLevel.basis, s, grid)
                            amplitude = amplitude + wn * fLevel.mc[r] * OMp * nLevel.mc[t] * nLevel.mc[tp] * Vee *
                                                    iLevel.mc[s] / (eni - enn)
                        end
                    end
                end

                #   Compute <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, kind) || alpha_i J_i>
                if  symn == symf
                    for  t = 1:nn
                        if  nLevel.mc[t] == 0.  continue    end
                        Vee = ManyElectron.matrixElement_Vee(CoulombInteraction(), fLevel.basis, r, nLevel.basis, t, grid)
                        for  tp = 1:nn
                            if  nLevel.mc[tp] == 0.  continue    end
                            OMp       = ManyElectron.matrixElement_Mab(Mp, gauge, omega, nLevel.basis, tp, iLevel.basis, s, grid)
                            amplitude = amplitude + wn * fLevel.mc[r] * Vee * nLevel.mc[t] * nLevel.mc[tp] * OMp *
                                                    iLevel.mc[s] / (eni - omega - enn)
                        end
                    end
                end
            end
        end
    end
    if  printout   printstyled("done. \n", color=:light_green)    end

    if  display
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " PhotoDouble^($Mp, absorption) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end

    return( amplitude )
end


"""
`PhotoDoubleIonization.generateIntermediateLevels(symn::LevelSymmetry, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
        energyGrid::Radial.GridGL, contSettings::Continuum.Settings, maxKappa::Int64;
        nuclearPot::Union{Nothing,Radial.Potential}=nothing, primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... builds the intermediate levels of the second-order sum for ONE required total symmetry symn, by adding a continuum partial wave to
        every level of the (N-1)-electron Green multiplet. Every kappa with abs(kappa) <= maxKappa that couples a Green level to symn is
        taken, at every energy of energyGrid.

        These states are not bound: they carry one continuum electron, which is what lets the one-body photon operator connect them to a
        final state carrying two. Because the partial waves are normalized per energy interval, the corresponding quadrature weight of
        energyGrid is returned with each level, and the caller must use it -- the sum over intermediate states is an integral over the
        intermediate electron energy.

        A tuple (nLevels, nWeights)::Tuple{Array{Level,1}, Array{Float64,1}} is returned.
"""
function generateIntermediateLevels(symn::LevelSymmetry, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                          energyGrid::Radial.GridGL, contSettings::Continuum.Settings, maxKappa::Int64;
                                          nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                          primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                          spectators::Union{Nothing,Array{Tuple{Subshell,Float64},1}}=nothing)
    nLevels = Level[];    nWeights = Float64[]

    # A ONE-BODY operator connects this intermediate level to the final state, so the continuum electron that
    # the photon does NOT act upon must be the very same orbital in both -- same kappa, same energy AND the same
    # subshell label, since two continuum subshells that differ only in their label count as different orbitals
    # and make the one-body matrix element vanish.  The caller therefore passes the pair's own subshells.
    if  !isnothing(spectators)
        for  gLevel in gMultiplet.levels
            gLevel = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
            symg   = LevelSymmetry(gLevel.J, gLevel.parity)
            for  (subsh, en) in spectators
                if  !(subsh.kappa in AngularMomentum.allowedKappaSymmetries(symg, symn))    continue    end
                nOrbital, nPhase = Continuum.generateOrbitalForLevel(en, subsh, gLevel, nm, grid, contSettings;
                                                                    nuclearPot=nuclearPot, primitives=primitives)
                push!(nLevels,  Basics.generateLevelWithExtraElectron(nOrbital, symn, gLevel))
                push!(nWeights, 1.0)
            end
        end
        return( nLevels, nWeights )
    end

    for  gLevel in gMultiplet.levels
        # The Green multiplet's basis holds CSFs of every symmetry it was built from; only the CSFs of this
        # level's own symmetry may be coupled to a partial wave, so the basis is reduced first.
        gLevel = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
        symg   = LevelSymmetry(gLevel.J, gLevel.parity)
        for  kappa in AngularMomentum.allowedKappaSymmetries(symg, symn)
            if  abs(kappa) > maxKappa    continue    end
            shn = Subshell(103, kappa)
            for  (ie, en)  in  enumerate(energyGrid.t)
                nOrbital, nPhase = Continuum.generateOrbitalForLevel(en, shn, gLevel, nm, grid, contSettings;
                                                                    nuclearPot=nuclearPot, primitives=primitives)
                push!(nLevels,  Basics.generateLevelWithExtraElectron(nOrbital, symn, gLevel))
                push!(nWeights, energyGrid.wt[ie])
            end
        end
    end

    return( nLevels, nWeights )
end


"""
`PhotoDoubleIonization.computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
        nrContinuum::Int64, settings::PhotoDoubleIonization.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing,
        primitives::Union{Nothing,Bsplines.Primitives}=nothing, printout::Bool=false)`
    ... computes the amplitudes of every partial-wave pair of the given line, at every energy sharing, and from them the energy-differential
        cross section of each sharing. TWO continuum orbitals are generated per pair -- one for each emitted electron, at the energies the
        pair carries -- and they are shared by all channels of that pair, which is what the pair layer is for.

        The two electrons are coupled in sequence: the first to the final ion, giving the intermediate xSymmetry, and the second to that,
        giving the total symmetry of the channel. An electric multipole is evaluated once per gauge into a single EmPropertyC, a magnetic
        one once into an EmPropertyC with equal components, so that no gauge is iterated over here.

        The energy-differential cross section dsigma/depsilon1 of each sharing is the INCOHERENT sum of abs2 over pairs, channels and
        multipoles, times the same prefactor 4 pi^2 alpha omega / (2 (2J_i+1)) that photoionization uses -- the operator and the per-energy
        normalization of the continuum orbitals are the same, and here BOTH outgoing electrons carry it. The total cross section is the
        sharing-integrated differential one; the exchange double counting is already handled by Basics.determineEnergySharings and is not
        applied again. A newLine::PhotoDoubleIonization.Line with all amplitudes and both cross sections evaluated is returned.
"""
function  computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                            nrContinuum::Int64, settings::PhotoDoubleIonization.Settings;
                                            nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                            primitives::Union{Nothing,Bsplines.Primitives}=nothing, printout::Bool=false)
    newSharings  = PhotoDoubleIonization.Sharing[]
    contSettings = Continuum.Settings(false, nrContinuum)
    symi         = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    # The sum over intermediate states is an integral over the intermediate electron energy, because those states
    # are normalized per energy interval.  The same number of points is used as for the energy sharings.
    Ji2          = Basics.twice(line.initialLevel.J)
    csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
    maxIntEnergy = maximum([sh.epsilon1 + sh.epsilon2  for sh in line.sharings])
    intermediateGrid = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.01, maxIntEnergy, settings.NoEnergySharings; printout=false)

    for  sharing in line.sharings
        newPairs = PhotoDoubleIonization.PartialWavePair[]
        for  pw in sharing.partialWavePairs
            sh1 = Subshell(101, pw.kappa1);    sh2 = Subshell(102, pw.kappa2)
            # A fresh subshell list: appending to the level's own list would mutate the level that was passed in.
            subshells = copy(line.initialLevel.basis.subshells);    push!(subshells, sh1);    push!(subshells, sh2)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshells)
            redFLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshells[1:end-2])

            # Electron 1 couples to the final ion and gives the intermediate symmetry of the pair.
            cOrbital1, phase1 = Continuum.generateOrbitalForLevel(pw.energy1, sh1, line.finalLevel, nm, grid, contSettings;
                                                                 nuclearPot=nuclearPot, primitives=primitives)
            xLevel            = Basics.generateLevelWithExtraElectron(cOrbital1, pw.xSymmetry, redFLevel)
            # Electron 2 couples to that intermediate symmetry and gives the total symmetry of each channel.
            cOrbital2, phase2 = Continuum.generateOrbitalForLevel(pw.energy2, sh2, line.finalLevel, nm, grid, contSettings;
                                                                 nuclearPot=nuclearPot, primitives=primitives)

            newChannels = PhotoDoubleIonization.Channel[]
            for  ch in pw.channels
                cLevel  = Basics.generateLevelWithExtraElectron(cOrbital2, ch.symmetry, xLevel)
                # The intermediate states of the second-order sum carry ONE continuum electron and are built here,
                # for the two symmetries the two time orderings require: symi for the ordering in which the photon
                # is absorbed last, and the total symmetry of this channel for the ordering in which it is absorbed
                # first.  They come with the quadrature weights of their per-energy normalization.
                nLevelsI, nWeightsI = PhotoDoubleIonization.generateIntermediateLevels(symi, settings.gMultiplet, nm, grid,
                                            intermediateGrid, contSettings, settings.maxKappa;
                                            nuclearPot=nuclearPot, primitives=primitives,
                                            spectators=[(sh1, pw.energy1), (sh2, pw.energy2)])
                nLevelsF, nWeightsF = PhotoDoubleIonization.generateIntermediateLevels(ch.symmetry, settings.gMultiplet, nm, grid,
                                            intermediateGrid, contSettings, settings.maxKappa;
                                            nuclearPot=nuclearPot, primitives=primitives)
                nLevels  = vcat(nLevelsI,  nLevelsF)
                nWeights = vcat(nWeightsI, nWeightsF)
                newAmps  = MultipoleAmplitude[]
                for  ma in ch.amplitudes
                    mp = ma.multipole
                    if  string(mp)[1] == 'E'
                        ampC = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Coulomb,   line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout)
                        ampB = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Babushkin, line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout)
                        push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampC, ampB)))
                    else
                        ampM = PhotoDoubleIonization.amplitude(Absorption(), mp, Basics.Magnetic,  line.photonEnergy, cLevel,
                                                               newiLevel, nLevels, nWeights, grid; printout=printout)
                        push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampM, ampM)))
                    end
                end
                push!(newChannels, PhotoDoubleIonization.Channel(ch.symmetry, newAmps))
            end
            push!(newPairs, PhotoDoubleIonization.PartialWavePair(pw.kappa1, pw.energy1, phase1, pw.xSymmetry,
                                                                 pw.kappa2, pw.energy2, phase2, newChannels))
        end

        # Incoherent over pairs, channels and multipoles; abs2 of an EmPropertyC keeps the two gauges apart by itself.
        dcs = EmProperty(0., 0.)
        for  pw in newPairs,  ch in pw.channels,  ma in ch.amplitudes    dcs = dcs + abs2(ma.amplitude)    end
        # The same prefactor as for photoionization: the electron-photon operator and the per-energy
        # normalization of the continuum orbitals are the same, and BOTH outgoing electrons already carry it.
        dcs = csFactor * dcs
        push!(newSharings, PhotoDoubleIonization.Sharing(sharing.omega, sharing.epsilon1, sharing.epsilon2,
                                                         sharing.weight, dcs, newPairs))
    end

    # The total is the sharing-integrated differential cross section.  The exchange double counting of
    # eps1 <-> eps2 is already handled inside Basics.determineEnergySharings and must NOT be applied again.
    crossSection = EmProperty(0., 0.)
    for  sh in newSharings   crossSection = crossSection + sh.weight * sh.differentialCs   end

    newLine = PhotoDoubleIonization.Line(line.initialLevel, line.finalLevel, line.photonEnergy, crossSection, newSharings)

    return( newLine )
end


"""
`PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::PhotoDoubleIonization.Settings; output::Bool=true)`  
    ... to compute the photo-double ionization transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{PhotoDoubleIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::PhotoDoubleIonization.Settings; output::Bool=true)
    println("")
    printstyled("PhotoDoubleIonization.computeLines(): The computation of photo-double ionization properties starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    error("\n\nPhotoDoubleIonization is POSTPONED (15-Aug-2026) and does not produce meaningful numbers.\n" *
          ">>> The module was rebuilt on that date onto partial-wave PAIRS and it now runs end to end, but\n"   *
          ">>> its ABSOLUTE SCALE IS WRONG BY 3-4 ORDERS OF MAGNITUDE: sigma(2+)/sigma(+) comes out 0.001%\n"   *
          ">>> at 200 eV against 3-4% measured, and it trends the wrong way with photon energy.\n"              *
          ">>> Do NOT read the internal tests it passes (mirror symmetry to 1e-10, threshold 77.77 eV) as\n"    *
          ">>> evidence that it is nearly right; they were all passing while the result was this wrong.\n"      *
          ">>> The prefactor of Kornberg & Miraglia, Phys. Rev. A 48, 3714 (1993), Eq. (2) is not yet applied\n"*
          ">>> -- see the STATUS block in the module docstring for the three missing factors and the\n"         *
          ">>> truncated intermediate space. Remove this error only together with that work.\n")

    # The nuclear potential and the B-spline primitives are constant for the whole computation and are built once.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    #
    lines = PhotoDoubleIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoDoubleIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.
    for  line in lines,  sh in line.sharings   maxEnergy = max(maxEnergy, sh.epsilon1, sh.epsilon2)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoDoubleIonization.Line[]
    for  line in lines
        println("\n>> Calculate photo-double ionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
        newLine = PhotoDoubleIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings;
                                                                    nuclearPot=nuclearPot, primitives=primitives) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoDoubleIonization.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoDoubleIonization.displayResults(iostream, newLines, settings)     end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoDoubleIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)`  
    ... to determine a list of PhotoDoubleIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{PhotoDoubleIonization.Line,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)
    lines = PhotoDoubleIonization.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                # Add lines for all photon energies
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    excess   = omega_au - (fLevel.energy - iLevel.energy)
                    if  excess < 0.    continue   end  
                    sharings = PhotoDoubleIonization.determineSharingsAndChannels(fLevel, iLevel, omega_au, settings) 
                    push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, omega_au, EmProperty(0., 0.), sharings) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoDoubleIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64,
                                                          settings::PhotoDoubleIonization.Settings)`
    ... determines the energy sharings of the two emitted electrons and, for each of them, the pairs of outgoing partial waves together
        with the total symmetries they serve. The coupling is read off the types: the first electron couples to the final ion to give
        xSymmetry, the second couples to that to give the total symmetry of a channel, and that total symmetry is what a multipole must
        reach from the initial level. Every kappa with abs(kappa) > settings.maxKappa is discarded.

        The distinct pairs (kappa1, xSymmetry, kappa2) are collected FIRST and every multipole that reaches one of their total symmetries
        is then attached, so that a pair which serves several symmetries appears once and its two continuum orbitals are generated once.
        The gauge is not iterated over at all: an electric multipole yields one EmPropertyC holding both gauges, a magnetic one an
        EmPropertyC with equal components.

        An Array{PhotoDoubleIonization.Sharing,1} is returned, with all amplitudes and cross sections still zero.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64,
                                            settings::PhotoDoubleIonization.Settings)
    symi     = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf     = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    shift_au = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    sharings = PhotoDoubleIonization.Sharing[]

    # The excess energy is what remains of the photon energy after the two electrons have been removed.
    eSharings = Basics.determineEnergySharings(omega - (finalLevel.energy - initialLevel.energy), settings.NoEnergySharings)

    for  es in eSharings
        epsilon1 = es[1];    epsilon2 = es[2];    weight = es[3]

        # Collect, for every distinct pair of partial waves, which total symmetry is reached by which multipole.
        pairs = Dict{Tuple{Int64,LevelSymmetry,Int64}, Dict{LevelSymmetry,Array{EmMultipole,1}}}()
        for  mp in settings.multipoles
            for  symt in AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                for  kappa1 = -settings.maxKappa-1:settings.maxKappa
                    if  kappa1 == 0  ||  abs(kappa1) > settings.maxKappa    continue    end
                    for  symx in AngularMomentum.allowedTotalSymmetries(symf, kappa1)
                        for  kappa2 in AngularMomentum.allowedKappaSymmetries(symt, symx)
                            if  abs(kappa2) > settings.maxKappa    continue    end
                            key = (kappa1, symx, kappa2)
                            haskey(pairs, key)         ||  (pairs[key] = Dict{LevelSymmetry,Array{EmMultipole,1}}())
                            haskey(pairs[key], symt)   ||  (pairs[key][symt] = EmMultipole[])
                            mp in pairs[key][symt]     ||  push!(pairs[key][symt], mp)
                        end
                    end
                end
            end
        end

        newPairs = PhotoDoubleIonization.PartialWavePair[]
        for  (key, symtDict) in pairs
            kappa1, symx, kappa2 = key
            newChannels = PhotoDoubleIonization.Channel[]
            for  (symt, mpList) in symtDict
                amplitudes = MultipoleAmplitude[]
                for  mp in mpList   push!(amplitudes, MultipoleAmplitude(mp, EmPropertyC(0.)))    end
                push!(newChannels, PhotoDoubleIonization.Channel(symt, amplitudes))
            end
            push!(newPairs, PhotoDoubleIonization.PartialWavePair(kappa1, epsilon1 + shift_au, 0., symx,
                                                                 kappa2, epsilon2 + shift_au, 0., newChannels))
        end
        push!(sharings, PhotoDoubleIonization.Sharing(omega, epsilon1, epsilon2, weight, EmProperty(0., 0.), newPairs))
    end

    return( sharings )
end


"""
`PhotoDoubleIonization.displayLines(stream::IO, lines::Array{PhotoDoubleIonization.Line,1})`
    ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{PhotoDoubleIonization.Line,1})
    #
    # First, print lines and sharings
    nx = 94
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & sharings:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
    sb = sb * TableStrings.flushleft(54, "  i -- f         omega        epsilon_1    epsilon_2"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        #
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))    * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "   "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Second, print lines and channles
    nx = 130
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & channels:   ... channels are shown only for the first sharing")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(100, "Channels (all energies in " * TableStrings.inUnits("energy") * ")" ; na=5);              
    sb = sb * TableStrings.flushleft(100, "Multipole  Gauge      quasi-Subshell   J^P_x   kappa   -->    J^P_t"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        for (ic, ch) in enumerate(line.sharings[1].channels)
            if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * string(ch.multipole) * "         " * string(ch.gauge)[1:3] * "            " 
            sb = sb * string(ch.quasiSubshell) * "       " * string(ch.xSymmetry) * "   "
            sb = sb * string(Subshell(1,ch.kappa))[end-4:end] * "   -->    "
            sb = sb * string(ch.tSymmetry)
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
