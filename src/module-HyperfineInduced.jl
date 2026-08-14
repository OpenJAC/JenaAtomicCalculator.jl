
"""
`module  JAC.HyperfineInduced`
... a submodel of JAC for computing HYPERFINE-INDUCED transitions, i.e. radiative transitions that are forbidden
    for the isolated nucleus or the isolated electron cloud, and become possible only because the hyperfine
    interaction mixes the IJF-coupled levels. Three kinds are treated by one and the same amplitude:

    * ELECTRONIC   -- the nucleus stays put, the electronic level changes, and the decay proceeds through the
                      hyperfine admixture of a nearby electronic level. This is the classic hyperfine quenching
                      of the 3P_0,2 -- 1S_0 lines in beryllium-like ions, and of the 27Al+ clock transition.
    * NUCLEAR      -- the electronic level stays put and the NUCLEAR level changes; the transition borrows its
                      strength from the mixing of two nuclear states of equal F. This is nuclear hyperfine mixing
                      (NHM), the mechanism that shortens the 229Th isomer lifetime from hours to tens of ms and
                      the 205Pb one from ~15 min to ~32 ms.
    * MIXED        -- both change, i.e. the hyperfine electronic bridge.

    REWRITTEN 06-Aug-2026, and the module's own note of September 2025 asked for exactly this: the IJF_Vector /
    IJF_Level / IJF_Multiplet defined here duplicated the Hfs family, and "it is neither recommended nor useful to
    define analogue data structures in different modules". They are retired; everything now runs on
    Hfs.HfBasisVector / Hfs.HfLevel / Hfs.HfMultiplet. The previous version could not be used in any case: FIVE of
    its function bodies were dead through identical-signature duplicates that Julia resolves silently by load
    order (computeLines, determineIJFlevels, displayHfMultiplet, displayRates, amplitude), only one of its three
    `amplitude` methods was reachable, and the two unreachable ones called `amplitude_Wu` with a
    String argument for which no method exists.

    The structure now mirrors PhotoEmission one-for-one -- Settings, Channel, Line, amplitude,
    computeAmplitudesProperties, computeLines, determineChannels, determineLines, displayLines, displayRates,
    displayLifetimes -- with Hfs.HfLevel in place of Level, so that a reader who knows the one knows the other.
"""
module HyperfineInduced


using Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..Hfs, ..ManyElectron, ..Radial, ..Nuclear,
              ..TableStrings, ..PhotoEmission


"""
`struct  HyperfineInduced.Settings  <:  AbstractProcessSettings`
    ... defines the details and parameters for computing hyperfine-induced transitions.

    + multipoles         ::Array{EmMultipole,1}
        ... multipoles of the RADIATION FIELD to be included, i.e. of the emitted photon.
    + hfMultipoles       ::Array{EmMultipole,1}
        ... multipoles of the HYPERFINE INTERACTION that mixes the IJF-coupled levels [M1, E2, E3, ...]. These two
            lists are independent: the 235U case, for instance, emits an E1 photon while the mixing is carried by
            the nuclear E3 moment.
    + gauges             ::Array{UseGauge,1}      ... gauges to be included into the computations.
    + isomers            ::Array{Nuclear.Isomer,1}
        ... the nuclear states to be included in the hyperfine basis. ONE isomer gives ordinary hyperfine
            quenching; TWO or more admit nuclear hyperfine mixing, since two nuclear levels of equal F then
            appear as separate basis vectors and are mixed by the diagonalization.
    + mixingLevels       ::Array{Int64,1}
        ... electronic levels admitted into the hyperfine basis, so that they can LEND their strength through the
            hyperfine admixture. Empty -- the default -- means the WHOLE multiplet mixes, which is both the
            physically correct choice and a cheap one, since a hyperfine matrix has the dimension of a handful of
            levels times a handful of F values. Restricting it is an optimization, never physics: it changes the
            AMPLITUDE, not merely what is printed.
    + initialFvalues     ::Array{AngularJ64,1}    ... F values retained in the initial basis; all, if empty.
    + finalFvalues       ::Array{AngularJ64,1}    ... as above, for the final hyperfine levels.
    + calcOverview       ::Bool
        ... run in OVERVIEW mode: build the multiplets and the full hyperfine representation, print the electronic
            levels with their roles, the hyperfine levels with their stable indices, and a RANKED list of the
            admixture channels that will carry the transition -- then stop, without computing any amplitude. This
            is the cheap first pass that tells a user which levels matter, so that lineSelection can be written
            from evidence instead of guessed at and corrected.
    + printBefore        ::Bool                   ... print the selected lines before computing them.
    + calcLifetimes      ::Bool                   ... compute and tabulate level lifetimes from the summed rates.
    + lineSelection      ::LineSelection          ... selected lines, if any.
    + photonEnergyShift  ::Float64                ... an overall shift of all photon energies.

    HOW TRANSITIONS ARE SELECTED -- read this before setting anything (redesigned 06-Aug-2026).

    THERE IS EXACTLY ONE SELECTOR, `lineSelection`, and it works on HYPERFINE levels:

        lineSelection = LineSelection(true, indexPairs = [(3,1), (3,2)])
        lineSelection = LineSelection(true, symmetryPairs = [(LevelSymmetry(AngularJ64(1), Basics.minus),
                                                             LevelSymmetry(AngularJ64(2), Basics.minus))])

    * `indexPairs` are (initial, final) indices of HYPERFINE levels -- the stable, energy-ordered index that
      Hfs.computeHyperfineRepresentation assigns and that every table in this module prints. They are NOT
      electronic level indices.
    * `symmetryPairs` are (initial, final) LevelSymmetry built from F and the parity -- the HYPERFINE F, not the
      electronic J.
    * Empty / inactive means every energetically allowed line is computed.

    GET THE INDICES FROM `calcOverview = true`, which builds everything, prints each hyperfine level with its
    index, its F^P, its dominant electronic parent and its dominant nuclear parent, ranks the admixture channels
    that will carry the transition, and stops without computing an amplitude. That is the intended first step on
    any new system; writing indices without it is guesswork.

    WHY THERE IS ONLY ONE. Earlier versions also selected by electronic parent (`iLevelIndex`/`iAddIndices`,
    later `initialLevels`/`finalLevels`) and briefly by nuclear parent. Both were wrong in practice:

    * selecting by electronic parent CANNOT express the question an isomer lifetime asks. In B-like 205Pb77+ the
      isomer-dominated hyperfine levels sit on BOTH 2p_1/2 and 2p_3/2, so any electronic selection that keeps
      them also keeps the ordinary 2p_3/2 -> 2p_1/2 M1 line at 1e11 /s and buries the hyperfine-induced effect
      by ten orders of magnitude.
    * `iAddIndices` additionally widened the BASIS and the set of REPORTED levels at the same time, so naming a
      perturber silently reported its own allowed decay as though it were the line of interest.

    A hyperfine level is identified by neither parent alone; it is identified by its index. Hence one selector,
    on that index, and `mixingLevels` -- which is about the basis and therefore about the amplitude -- kept
    strictly separate from it.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                  ::Array{EmMultipole,1}
    hfMultipoles                ::Array{EmMultipole,1}
    gauges                      ::Array{UseGauge,1}
    isomers                     ::Array{Nuclear.Isomer,1}
    mixingLevels                ::Array{Int64,1}
    initialFvalues              ::Array{AngularJ64,1}
    finalFvalues                ::Array{AngularJ64,1}
    calcOverview                ::Bool
    printBefore                 ::Bool
    calcLifetimes               ::Bool
    lineSelection               ::LineSelection
    photonEnergyShift           ::Float64
end


"""
`HyperfineInduced.Settings()`  ... constructor for the default values.
"""
function Settings()
    Settings(EmMultipole[E1], EmMultipole[M1], UseGauge[Basics.UseCoulomb], Nuclear.Isomer[],
             Int64[], AngularJ64[], AngularJ64[], false, false, true, LineSelection(), 0.)
end


"""
`HyperfineInduced.Settings(set::HyperfineInduced.Settings;`

        multipoles=..,      hfMultipoles=..,    gauges=..,          isomers=..,
        mixingLevels=..,    initialFvalues=..,
        finalFvalues=..,    calcOverview=..,    printBefore=..,     calcLifetimes=..,
        lineSelection=..,   photonEnergyShift=..)

    ... the standard JAC keyword copy-constructor, which this Settings previously lacked.
"""
function Settings(set::HyperfineInduced.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    hfMultipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
    gauges::Union{Nothing,Array{UseGauge,1}}=nothing,           isomers::Union{Nothing,Array{Nuclear.Isomer,1}}=nothing,
    mixingLevels::Union{Nothing,Array{Int64,1}}=nothing,        initialFvalues::Union{Nothing,Array{AngularJ64,1}}=nothing,
    finalFvalues::Union{Nothing,Array{AngularJ64,1}}=nothing,   calcOverview::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,
    calcLifetimes::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,        photonEnergyShift::Union{Nothing,Float64}=nothing)

    if  isnothing(multipoles)         multipolesx        = set.multipoles         else  multipolesx        = multipoles         end
    if  isnothing(hfMultipoles)       hfMultipolesx      = set.hfMultipoles       else  hfMultipolesx      = hfMultipoles       end
    if  isnothing(gauges)             gaugesx            = set.gauges             else  gaugesx            = gauges             end
    if  isnothing(isomers)            isomersx           = set.isomers            else  isomersx           = isomers            end
    if  isnothing(mixingLevels)       mixingLevelsx      = set.mixingLevels       else  mixingLevelsx      = mixingLevels       end
    if  isnothing(initialFvalues)     initialFvaluesx    = set.initialFvalues     else  initialFvaluesx    = initialFvalues     end
    if  isnothing(finalFvalues)       finalFvaluesx      = set.finalFvalues       else  finalFvaluesx      = finalFvalues       end
    if  isnothing(calcOverview)       calcOverviewx      = set.calcOverview       else  calcOverviewx      = calcOverview       end
    if  isnothing(printBefore)        printBeforex       = set.printBefore        else  printBeforex       = printBefore        end
    if  isnothing(calcLifetimes)      calcLifetimesx     = set.calcLifetimes      else  calcLifetimesx     = calcLifetimes      end
    if  isnothing(lineSelection)      lineSelectionx     = set.lineSelection      else  lineSelectionx     = lineSelection      end
    if  isnothing(photonEnergyShift)  photonEnergyShiftx = set.photonEnergyShift  else  photonEnergyShiftx = photonEnergyShift  end

    Settings(multipolesx, hfMultipolesx, gaugesx, isomersx, mixingLevelsx, initialFvaluesx, finalFvaluesx,
             calcOverviewx, printBeforex, calcLifetimesx, lineSelectionx, photonEnergyShiftx)
end


# `Base.show(io::IO, settings::HyperfineInduced.Settings)`  ... prepares a proper printout.
function Base.show(io::IO, settings::HyperfineInduced.Settings)
    println(io, "multipoles (radiation):    $(settings.multipoles)  ")
    println(io, "hfMultipoles (interaction):$(settings.hfMultipoles)  ")
    println(io, "gauges:                    $(settings.gauges)  ")
    println(io, "isomers:                   $(length(settings.isomers)) nuclear state(s)  ")
    println(io, "mixingLevels (admix only): $(isempty(settings.mixingLevels)  ? "all of the multiplet" :
                                              settings.mixingLevels)  ")
    println(io, "initialFvalues:            $(isempty(settings.initialFvalues) ? "all" : settings.initialFvalues)  ")
    println(io, "finalFvalues:              $(isempty(settings.finalFvalues)   ? "all" : settings.finalFvalues)  ")
    println(io, "calcOverview:              $(settings.calcOverview)  ")
    println(io, "printBefore:               $(settings.printBefore)  ")
    println(io, "calcLifetimes:             $(settings.calcLifetimes)  ")
    println(io, "lineSelection:             $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:         $(settings.photonEnergyShift)  ")
end


"""
`struct  HyperfineInduced.Channel`
    ... one channel of a hyperfine-induced transition, i.e. one (multipole, gauge) of the radiation field.

    + multipole  ::EmMultipole        ... multipole of the emitted photon.
    + gauge      ::EmGauge            ... gauge.
    + amplitude  ::Complex{Float64}   ... transition amplitude.
    + rate       ::Float64            ... rate of this channel.
"""
struct Channel
    multipole           ::EmMultipole
    gauge               ::EmGauge
    amplitude           ::Complex{Float64}
    rate                ::Float64
end


# `Base.show(io::IO, channel::HyperfineInduced.Channel)`  ... prepares a proper printout.
function Base.show(io::IO, channel::HyperfineInduced.Channel)
    println(io, "multipole: $(channel.multipole), gauge: $(channel.gauge), amplitude: $(channel.amplitude), " *
                "rate: $(channel.rate)")
end


"""
`struct  HyperfineInduced.Line`
    ... one hyperfine-induced transition between two IJF-coupled levels.

    + initialLevel ::Hfs.HfLevel                        ... initial hyperfine level.
    + finalLevel   ::Hfs.HfLevel                        ... final hyperfine level.
    + omega        ::Float64                            ... photon energy.
    + photonRate   ::EmProperty                         ... total rate, summed over the channels.
    + channels     ::Array{HyperfineInduced.Channel,1}  ... the individual (multipole, gauge) channels.
"""
struct  Line
    initialLevel     ::Hfs.HfLevel
    finalLevel       ::Hfs.HfLevel
    omega            ::Float64
    photonRate       ::EmProperty
    channels         ::Array{HyperfineInduced.Channel,1}
end


"""
`HyperfineInduced.Line(initialLevel::Hfs.HfLevel, finalLevel::Hfs.HfLevel, omega::Float64)`
    ... constructor for a line without channels or rate yet.
"""
function Line(initialLevel::Hfs.HfLevel, finalLevel::Hfs.HfLevel, omega::Float64)
    Line(initialLevel, finalLevel, omega, EmProperty(0., 0.), HyperfineInduced.Channel[])
end


# `Base.show(io::IO, line::HyperfineInduced.Line)`  ... prepares a proper printout.
function Base.show(io::IO, line::HyperfineInduced.Line)
    println(io, "initialLevel (F):   $(line.initialLevel.F)  ")
    println(io, "finalLevel (F):     $(line.finalLevel.F)  ")
    println(io, "omega:              $(line.omega)  ")
    println(io, "photonRate:         $(line.photonRate)  ")
    println(io, "channels:           $(length(line.channels))  ")
end


"""
`HyperfineInduced.doubleFactorial(n::Int64)`
    ... to return the double factorial n!! for n >= -1, as a Float64.

        Written out because the previous implementation was a three-branch lookup that returned values for
        n = 3, 5, 7 and called error("stop a") for anything else -- which silently confined the whole module to
        multipoles L = 1, 2, 3 and would have aborted on the E3 hyperfine interaction of 235U at L = 3 combined
        with any higher radiation multipole.
"""
function  doubleFactorial(n::Int64)
    if  n <= 0    return( 1.0 )    end
    wa = 1.0;   k = n
    while  k > 1    wa = wa * k;   k = k - 2    end
    return( wa )
end


"""
`HyperfineInduced.amplitude(::Emission, mp::EmMultipole, gauge::EmGauge, omega::Float64,
                            finalLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel, grid::Radial.Grid; printout::Bool=true)`
    ... to compute the hyperfine-induced transition amplitude between two IJF-coupled levels; an
        amplitude::ComplexF64 is returned.

        THE STRUCTURE. Each hyperfine level is a superposition sum_k mc[k] |(I_k J_k) F>, so the amplitude is a
        double sum over the components of the two levels, weighted by their mixing coefficients. Within one pair
        of components exactly one of two things can happen, since the radiation operator acts either on the
        nucleus or on the electrons but never on both:

        * the ELECTRONIC level changes while the isomer does not -- the ordinary multipole amplitude, recoupled
          with Hfs.recouplingElectronicOperator (nucleus as spectator);
        * the NUCLEAR level changes while the electronic level does not -- the nuclear multipole matrix element
          <I_p||W^L||I_q>, recoupled with Hfs.recouplingNuclearOperator (electrons as spectator), and multiplied
          by the multipole-field factor below so that it carries the same normalization as the electronic one.

        Both are summed here, which is what lets ONE routine cover purely electronic, purely nuclear and mixed
        hyperfine-induced transitions. The mixing coefficients themselves already contain all the hyperfine
        amplitudes and energy denominators, having come out of the diagonalization in
        Hfs.computeHyperfineRepresentation.

        REPLACES an amplitude that omitted the sqrt((2F+1)(2F'+1)) of the recoupling, gated the nuclear term on a
        magic `abs(mc*mc) > 0.8`, and used a double factorial defined only for L = 3, 5, 7. The `0.8` gate in
        particular silenced the nuclear term for exactly the strongly mixed cases the module exists to describe.
"""
function  amplitude(::Emission, mp::EmMultipole, gauge::EmGauge, omega::Float64,
                    finalLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel, grid::Radial.Grid; printout::Bool=true)
    amp   = ComplexF64(0.)
    alpha = Defaults.getDefaults("alpha")
    L     = mp.L
    ## THE MULTIPOLE FIELD FACTOR, derived 06-Aug-2026 rather than assumed, and verified against the bare
    ## nucleus. It puts a nuclear reduced matrix element on the same footing as PhotoEmission's electronic
    ## amplitude, i.e. it is fixed by demanding that the rate below reproduce the standard result
    ##
    ##     A(sigma L)  =  8 pi (L+1) / (L [(2L+1)!!]^2) * (omega/c)^(2L+1) * B(sigma L)
    ##
    ## Inserting A = 8 pi alpha omega/(2 F_i + 1) |elementM * fld|^2 with elementM = sqrt((2 I_i + 1) B), the
    ## degeneracies cancel identically and everything collapses to
    ##
    ##     fld  =  sqrt( (L+1)/L ) * (alpha omega)^L / (2L+1)!!
    ##
    ## with NO 4 pi, no (2L+1), and -- the point -- no dependence on the nuclear spin. An earlier version had
    ## sqrt((L+1)(2L+1)/(4 pi L)) and was too fast by ~9300; the spurious spin dependence needed to repair it
    ## was what revealed that the manuscript's elementM convention, not this factor, was the real error.
    ##
    ## THE EXTRA alpha FOR MAGNETIC MULTIPOLES is the 1/c of the magnetic multipole radiation operator. JAC
    ## stores nuclear moments with mu_N = 5.446170e-4/2 (the mu_B = 1/2 convention), which is correct for the
    ## hyperfine INTERACTION -- V11 and V22 match Shabaev to 5 % -- whereas the radiation formula is written in
    ## the Gaussian convention mu_B = alpha/2. The factor belongs to the photon operator, so it is applied here
    ## and not to the moment, which must stay consistent with Isomer.mu and with Hfs.
    fld   = sqrt( (L+1)/L ) * (alpha * omega)^L / HyperfineInduced.doubleFactorial(2L+1)
    if  !mp.electric    fld = fld * alpha    end
    #
    for  (q, iVec)  in  enumerate(initialLevel.hfBasisVectors)
        if  abs(initialLevel.mc[q]) < 1.0e-12    continue    end
        for  (p, fVec)  in  enumerate(finalLevel.hfBasisVectors)
            if  abs(finalLevel.mc[p]) < 1.0e-12    continue    end
            weight = finalLevel.mc[p] * initialLevel.mc[q]
            #
            ## (1) NUCLEAR radiation: the electronic level is untouched, the isomer changes.
            if  fVec.levelJ.index == iVec.levelJ.index
                wn = Hfs.computeInteractionAmplitudeM(mp, fVec.isomer, iVec.isomer)
                if  fVec.isomer != iVec.isomer  &&  wn != 0.
                    wa = Hfs.recouplingNuclearOperator(iVec.levelJ.J, iVec.isomer.spinI, iVec.F,
                                                       fVec.isomer.spinI, fVec.F, L)
                    amp = amp + weight * wa * wn * fld
                end
            end
            #
            ## (2) ELECTRONIC radiation: the isomer is untouched, the electrons radiate.
            ##
            ## NO CONDITION on the electronic level, and this matters (corrected 06-Aug-2026): requiring the
            ## electronic level to CHANGE looks natural, but it silently discards the dominant channel of nuclear
            ## hyperfine mixing. In H-like 229Th89+ every basis vector carries the same 1s level, so that
            ## condition switched the electronic operator off altogether and left only the slow nuclear M1 -- an
            ## enhancement of 4 orders instead of the ~6 that Shabaev reports. The fast path is the ordinary
            ## HYPERFINE M1 transition F -> F' within ONE nuclear state, driven by the electron's magnetic moment;
            ## the isomer-dominated level borrows it through its ground-state admixture. Same electronic level,
            ## same isomer, different F: non-zero precisely because the recoupling factor below carries the F
            ## dependence, and PhotoEmission.amplitude then supplies the level's own magnetic moment.
            if  fVec.isomer == iVec.isomer
                wa = Hfs.recouplingElectronicOperator(iVec.isomer.spinI, iVec.levelJ.J, iVec.F,
                                                      fVec.levelJ.J, fVec.F, L)
                if  wa != 0.
                    we = PhotoEmission.amplitude(Emission(), mp, gauge, omega, fVec.levelJ, iVec.levelJ, grid,
                                                 display=false, printout=false)
                    amp = amp + weight * wa * we
                end
            end
        end
    end
    return( amp )
end


"""
`HyperfineInduced.computeAmplitudesProperties(line::HyperfineInduced.Line, grid::Radial.Grid,
                            settings::HyperfineInduced.Settings; printout::Bool=true)`
    ... to compute all amplitudes and the rate of the given line; a new HyperfineInduced.Line is returned.

        THE STATISTICAL FACTOR, corrected 06-Aug-2026. The rate is

            A = 8 pi alpha omega / (2 F_initial + 1) * |amplitude|^2

        exactly as PhotoEmission.computeAmplitudesProperties forms it for an electronic line, with the hyperfine
        weight (2F+1) of the EMITTING level in place of (2J+1). The previous version instead multiplied by
        `Basics.twice(line.finalLevel.F)`, i.e. by 2F of the FINAL level -- the wrong level, and 2F rather than
        2F+1, differing by about a factor twenty for F = 2. Using PhotoEmission's own prefactor also makes the
        I -> 0 limit exact rather than approximate: with no nuclear spin every hyperfine level coincides with its
        electronic parent, the recoupling factors reduce to 1, and this rate must reproduce PhotoEmission's for
        the same transition. That is the check to run first if the absolute scale is ever in doubt.
"""
function  computeAmplitudesProperties(line::HyperfineInduced.Line, grid::Radial.Grid,
                                      settings::HyperfineInduced.Settings; printout::Bool=true)
    newChannels = HyperfineInduced.Channel[];    rateC = 0.;    rateB = 0.
    wa = 8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.F) + 1)
    for  channel in line.channels
        amplitude = HyperfineInduced.amplitude(Emission(), channel.multipole, channel.gauge, line.omega,
                                               line.finalLevel, line.initialLevel, grid, printout=printout)
        rate      = wa * abs(amplitude)^2
        push!( newChannels, HyperfineInduced.Channel(channel.multipole, channel.gauge, amplitude, rate) )
        if       channel.gauge == Basics.Coulomb     rateC = rateC + rate
        elseif   channel.gauge == Basics.Babushkin   rateB = rateB + rate
        elseif   channel.gauge == Basics.Magnetic    rateB = rateB + rate;   rateC = rateC + rate
        end
    end
    return( HyperfineInduced.Line(line.initialLevel, line.finalLevel, line.omega,
                                  EmProperty(rateC, rateB), newChannels) )
end


"""
`HyperfineInduced.determineChannels(finalLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel,
                            settings::HyperfineInduced.Settings)`
    ... to determine the (multipole, gauge) channels of one hyperfine-induced transition; an
        Array{HyperfineInduced.Channel,1} is returned. The selection rules are those of the TOTAL angular momenta
        F and the TOTAL parity, nuclear x electronic -- which is the point of working in the IJF basis: a
        transition forbidden for the electrons alone, or for the nucleus alone, may still be allowed for F.
"""
function determineChannels(finalLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel, settings::HyperfineInduced.Settings)
    channels = HyperfineInduced.Channel[]
    symi = LevelSymmetry(initialLevel.F, initialLevel.parity)
    symf = LevelSymmetry(finalLevel.F,   finalLevel.parity)
    for  mp in settings.multipoles
        if  AngularMomentum.isAllowedMultipole(symi, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                if      gauge == Basics.UseCoulomb    &&  mp.electric
                    push!(channels, HyperfineInduced.Channel(mp, Basics.Coulomb,   ComplexF64(0.), 0.) )
                elseif  gauge == Basics.UseBabushkin  &&  mp.electric
                    push!(channels, HyperfineInduced.Channel(mp, Basics.Babushkin, ComplexF64(0.), 0.) )
                elseif  !(mp.electric)  &&  !hasMagnetic
                    push!(channels, HyperfineInduced.Channel(mp, Basics.Magnetic,  ComplexF64(0.), 0.) )
                    hasMagnetic = true
                end
            end
        end
    end
    return( channels )
end


"""
`HyperfineInduced.determineLines(finalMultiplet::Hfs.HfMultiplet, initialMultiplet::Hfs.HfMultiplet,
                            settings::HyperfineInduced.Settings)`
    ... to determine all hyperfine-induced lines with a positive photon energy and at least one open channel;
        an Array{HyperfineInduced.Line,1} is returned.
"""
function  determineLines(finalMultiplet::Hfs.HfMultiplet, initialMultiplet::Hfs.HfMultiplet,
                         settings::HyperfineInduced.Settings)
    lines  = HyperfineInduced.Line[]
    eShift = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    ## The selection is done on indices and F^P here, NOT through Basics.selectLevelPair: that method dispatches
    ## on ManyElectron.Level, and an Hfs.HfLevel is not one.
    ##
    ## THE INDEX USED IS `lev.index`, the stable energy-ordered label that Hfs.computeHyperfineRepresentation
    ## assigns to the FULL diagonalized multiplet -- not the position in this loop. The two differ whenever
    ## restricting mixingLevels or the F values changes the size of the multiplet, and it is the stable index
    ## that a user reads off the printed tables and puts into a LineSelection. Selecting by loop position would
    ## silently mean a different pair as soon as the basis changed.
    for  iLevel  in  initialMultiplet.hfLevels
        for  fLevel  in  finalMultiplet.hfLevels
            if  settings.lineSelection.active
                isym = LevelSymmetry(iLevel.F, iLevel.parity);   fsym = LevelSymmetry(fLevel.F, fLevel.parity)
                if  !( (iLevel.index, fLevel.index) in settings.lineSelection.indexPairs  ||
                       (isym,fsym) in settings.lineSelection.symmetryPairs )    continue    end
            end
            omega = iLevel.energy - fLevel.energy + eShift
            if  omega <= 0.    continue    end
            channels = HyperfineInduced.determineChannels(fLevel, iLevel, settings)
            if  length(channels) == 0    continue    end
            push!( lines, HyperfineInduced.Line(iLevel, fLevel, omega, EmProperty(0., 0.), channels) )
        end
    end
    return( lines )
end


"""
`HyperfineInduced.selectHyperfineBasisVectors(basis::Array{Hfs.HfBasisVector,1}, Fvalues::Array{AngularJ64,1})`
    ... to keep only those basis vectors whose total angular momentum F appears in Fvalues; all vectors are kept
        if Fvalues is empty. An Array{Hfs.HfBasisVector,1} is returned.

        Renamed from `generateBasis` (06-Aug-2026): it never generated anything, it selected from a basis that
        Hfs.defineHyperfineBasis had already built.
"""
function  selectHyperfineBasisVectors(basis::Array{Hfs.HfBasisVector,1}, Fvalues::Array{AngularJ64,1})
    if  length(Fvalues) == 0    return( basis )    end
    return( filter(v -> v.F in Fvalues, basis) )
end


"""
`HyperfineInduced.determineHyperfineMultiplet(multiplet::Multiplet, isomers::Array{Nuclear.Isomer,1},
                            levelIndex::Int64, addIndices::Array{Int64,1}, Fvalues::Array{AngularJ64,1},
                            grid::Radial.Grid, settings::HyperfineInduced.Settings; printout::Bool=false)`
    ... to build and diagonalize the IJF-coupled hyperfine representation for the selected electronic levels and
        nuclear states; an Hfs.HfMultiplet is returned.

        Renamed from `determineIJFlevels` (06-Aug-2026), together with the retirement of the IJF_* types: it
        returns an Hfs.HfMultiplet, and the name should say so. The electronic levels admitted are levelIndex
        together with addIndices -- the level whose hyperfine expansion is wanted, plus the neighbours whose
        admixture makes the transition possible in the first place.
"""
function  determineHyperfineMultiplet(multiplet::Multiplet, isomers::Array{Nuclear.Isomer,1},
                                      mixingLevels::Array{Int64,1}, Fvalues::Array{AngularJ64,1},
                                      grid::Radial.Grid, settings::HyperfineInduced.Settings; printout::Bool=false)
    nl = length(multiplet.levels)
    for  k in mixingLevels
        if  k < 1  ||  k > nl
            error("\n\nHyperfineInduced.determineHyperfineMultiplet():  STOP -- level index $k was requested in " *
                  "mixingLevels, \nbut the electronic multiplet has only $nl levels.\n")
        end
    end
    ## An empty mixingLevels means the WHOLE multiplet mixes -- the physical default, and a cheap one. Note that
    ## this list governs the BASIS and therefore the amplitude itself; which transitions are reported is a
    ## separate question, settled solely by settings.lineSelection on the hyperfine indices below.
    if  length(mixingLevels) == 0    basisIndices = collect(1:nl)
    else                             basisIndices = sort(unique(mixingLevels))
    end
    subMultiplet = Multiplet("selected levels", [multiplet.levels[k] for k in basisIndices])
    basis        = Hfs.defineHyperfineBasis(subMultiplet, isomers; printout=printout)
    basis        = HyperfineInduced.selectHyperfineBasisVectors(basis, Fvalues)
    if  length(basis) == 0
        error("\n\nHyperfineInduced.determineHyperfineMultiplet():  STOP -- no hyperfine basis vector survived " *
              "the selection.\n>>> The requested F values $Fvalues are not compatible with I + J for any of the " *
              "selected levels.\n")
    end
    return( Hfs.computeHyperfineRepresentation(basis, grid; hfMultipoles=settings.hfMultipoles,
                                               printout=printout) )
end


"""
`HyperfineInduced.computeHfsCoefficientsAB(multiplet::Multiplet, isomer::Nuclear.Isomer, index::Int64,
                            grid::Radial.Grid)`
    ... to compute the ordinary hyperfine A and B coefficients of one electronic level, as a diagnostic on the
        nuclear input BEFORE any transition amplitude is evaluated. A tuple (A, B)::Tuple{Float64,Float64} in
        atomic units is returned.

        Renamed from `splitting_AB` (06-Aug-2026). It is worth running: A and B depend on the same nuclear moments
        and the same electronic matrix elements as the hyperfine mixing does, but are directly comparable with
        tabulated values, so a wrong mu, Q or level index shows up here rather than in a final rate.
"""
function  computeHfsCoefficientsAB(multiplet::Multiplet, isomer::Nuclear.Isomer, index::Int64, grid::Radial.Grid)
    if  index < 1  ||  index > length(multiplet.levels)    return( (0., 0.) )    end
    level = multiplet.levels[index]
    J     = AngularMomentum.oneJ(level.J);    Ix = AngularMomentum.oneJ(isomer.spinI)
    if  J == 0.  ||  Ix == 0.    return( (0., 0.) )    end
    ampM1 = real( Hfs.amplitude(Basics.M1, level, level, grid; printout=false) )
    wx    = Defaults.convertUnits("moment: from nuclear magneton to atomic", 1.0)
    A     = ampM1 / sqrt(J*(J+1)) * wx * isomer.mu / Ix
    B     = 0.
    if  !(isomer.spinI in [AngularJ64(0), AngularJ64(1//2)])  &&  level.J != AngularJ64(1//2)
        ampE2 = real( Hfs.amplitude(Basics.E2, level, level, grid; printout=false) )
        wy    = Defaults.convertUnits("cross section: from barn to atomic unit", 1.0)
        B     = 2 * ampE2 * sqrt( (2J-1) / ((J+1)*(2J+3)) ) * wy * isomer.Q
    end
    return( (A, B) )
end


"""
`HyperfineInduced.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                            grid::Radial.Grid, settings::HyperfineInduced.Settings; output=true)`
    ... to compute all hyperfine-induced transition amplitudes and properties as requested by the settings; a list
        of lines::Array{HyperfineInduced.Line,1} is returned if output = true.

        The two multiplets always refer to the ELECTRONIC system; every piece of nuclear information comes from
        settings.isomers, which replaces the isotopic data of an ordinary structure computation. Passing two or
        more isomers is what admits nuclear hyperfine mixing.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                       settings::HyperfineInduced.Settings; output=true)
    subshellList = Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    println("")
    printstyled("HyperfineInduced.computeLines(): The computation of hyperfine-induced transitions starts now ... \n",
                color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    if  length(settings.isomers) == 0
        error("\n\nHyperfineInduced.computeLines():  STOP -- settings.isomers is empty, so there is no nuclear "     *
              "state to couple to \nand nothing hyperfine-induced can happen.\n>>> Supply at least one "             *
              "Nuclear.Isomer; supply two (ground + isomeric) for nuclear hyperfine mixing.\n")
    end
    println(">>> $(length(settings.isomers)) nuclear state(s); radiation multipoles $(settings.multipoles), " *
            "hyperfine-interaction multipoles $(settings.hfMultipoles)")
    for  isomer in settings.isomers
        println("    I = $(isomer.spinI)$(isomer.parity),  excitation energy $(isomer.energy),  mu = $(isomer.mu)," *
                "  Q = $(isomer.Q),  nuclear transitions $(isomer.multipoleM)")
    end
    #
    ## The A and B coefficients first: they rest on the same nuclear moments and electronic matrix elements as the
    ## mixing does, but are directly comparable with tabulated values, so bad nuclear input shows up here.
    A, B = HyperfineInduced.computeHfsCoefficientsAB(initialMultiplet, settings.isomers[1], 1, grid)
    println(">>> initial level 1:  A = $(A) a.u.,  B = $(B) a.u.")
    A, B = HyperfineInduced.computeHfsCoefficientsAB(finalMultiplet, settings.isomers[1], 1, grid)
    println(">>> final level 1:    A = $(A) a.u.,  B = $(B) a.u.")
    #
    ## Show which electronic levels admix. With mixingLevels defaulting to the whole multiplet, this is what
    ## tells a user what the default actually did -- it governs the AMPLITUDE, not merely the printout.
    HyperfineInduced.displayElectronicBasis(stdout, initialMultiplet, settings.mixingLevels, "initial")
    HyperfineInduced.displayElectronicBasis(stdout, finalMultiplet,   settings.mixingLevels, "final")
    #
    ## OVERVIEW MODE: build everything, show what is there and what will carry the transition, then stop. No
    ## amplitude is evaluated, so this is the cheap pass that lets lineSelection be written from evidence.
    if  settings.calcOverview
        for  (mp, Fv, sRole) in [ (initialMultiplet, settings.initialFvalues, "initial"),
                                  (finalMultiplet,   settings.finalFvalues,   "final") ]
            hfAll = HyperfineInduced.determineHyperfineMultiplet(mp, settings.isomers, settings.mixingLevels,
                                                                 Fv, grid, settings)
            HyperfineInduced.displayHyperfineComposition(stdout, hfAll, sRole)
            HyperfineInduced.displayAdmixtureRanking(stdout, hfAll, sRole)
        end
        println("\n>>> Overview only (calcOverview = true); no amplitude was computed.")
        println(">>> Put the HYPERFINE level indices above into lineSelection, e.g.")
        println(">>>     lineSelection = LineSelection(true, indexPairs = [(3,1), (3,2)])\n")
        if  output    return( HyperfineInduced.Line[] )
        else          return( nothing )
        end
    end
    #
    initialHfMultiplet = HyperfineInduced.determineHyperfineMultiplet(initialMultiplet, settings.isomers,
                                settings.mixingLevels, settings.initialFvalues, grid, settings)
    finalHfMultiplet   = HyperfineInduced.determineHyperfineMultiplet(finalMultiplet, settings.isomers,
                                settings.mixingLevels, settings.finalFvalues, grid, settings)
    println("\n>>> hyperfine levels: $(length(initialHfMultiplet.hfLevels)) initial, " *
            "$(length(finalHfMultiplet.hfLevels)) final")
    HyperfineInduced.displayHyperfineComposition(stdout, initialHfMultiplet, "initial")
    #
    lines = HyperfineInduced.determineLines(finalHfMultiplet, initialHfMultiplet, settings)
    if  settings.printBefore    HyperfineInduced.displayLines(stdout, lines)    end
    #
    newLines = HyperfineInduced.Line[]
    for  line in lines
        push!( newLines, HyperfineInduced.computeAmplitudesProperties(line, grid, settings; printout=false) )
    end
    #
    HyperfineInduced.displayRates(stdout, newLines, settings)
    if  settings.calcLifetimes    HyperfineInduced.displayLifetimes(stdout, newLines, settings)    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        HyperfineInduced.displayRates(iostream, newLines, settings)
        if  settings.calcLifetimes    HyperfineInduced.displayLifetimes(iostream, newLines, settings)    end
    end
    #
    if  output    return( newLines )
    else          return( nothing )
    end
end


"""
`HyperfineInduced.dominantElectronicLevel(level::Hfs.HfLevel)`
    ... to return the electronic level that carries the largest weight in the given hyperfine level; a Level is
        returned, and its `.index` is the index in the original electronic multiplet.

        This is what decides whether a hyperfine level is REPORTED as decaying or is present only to lend its
        strength: a level is "the 3P_0 level" if 3P_0 dominates it, whatever small admixtures it also carries.
"""
function  dominantElectronicLevel(level::Hfs.HfLevel)
    p = 1;   wmax = -1.
    for  (i, mc)  in  enumerate(level.mc)
        if  abs(mc) > wmax   wmax = abs(mc);   p = i    end
    end
    return( level.hfBasisVectors[p].levelJ )
end


"""
`HyperfineInduced.displayElectronicBasis(stream::IO, multiplet::Multiplet, mixingLevels::Array{Int64,1},
                            role::String)`
    ... to list the electronic levels that enter the hyperfine basis, with index, symmetry and energy; nothing is
        returned.

        Every level listed here contributes to the AMPLITUDE through its hyperfine admixture, whether or not any
        transition involving it is finally reported -- reporting is settled separately, by lineSelection on the
        hyperfine indices. The table exists so that a user can see what `mixingLevels` did, since its default is
        to admit the whole multiplet.
"""
function  displayElectronicBasis(stream::IO, multiplet::Multiplet, mixingLevels::Array{Int64,1}, role::String)
    nl = length(multiplet.levels)
    if  length(mixingLevels) == 0    basisIndices = collect(1:nl)
    else                             basisIndices = sort(unique(mixingLevels))
    end
    nx = 76
    println(stream, " ")
    println(stream, "  Electronic levels in the $role hyperfine basis  (all of them admix; " *
            (length(mixingLevels) == 0 ? "mixingLevels = all" : "restricted by mixingLevels") * "):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "     level      J^P         energy " * TableStrings.inUnits("energy"))
    println(stream, "  ", TableStrings.hLine(nx))
    e0 = multiplet.levels[1].energy
    for  k in basisIndices
        lev = multiplet.levels[k]
        println(stream, "      " * TableStrings.center(6, string(k); na=3) *
                TableStrings.center(10, string(LevelSymmetry(lev.J, lev.parity)); na=3) *
                @sprintf("%14.6e", Defaults.convertUnits("energy: from atomic", lev.energy - e0)))
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`HyperfineInduced.displayHyperfineComposition(stream::IO, hfMultiplet::Hfs.HfMultiplet, role::String;
                            threshold::Float64=1.0e-3)`
    ... to list, for every reported hyperfine level, the electronic levels and nuclear states that compose it,
        with their mixing coefficients; nothing is returned.

        THIS IS THE TABLE THAT ANSWERS "which additional levels are of interest". The mixing coefficients are the
        physics: a hyperfine-induced transition happens because the decaying level is not pure, and the size of
        each admixture says directly how much each perturber contributes. Reading it also shows whether the basis
        was large enough -- if the largest admixture sits on the level highest in the basis, the basis is probably
        too small.
"""
function  displayHyperfineComposition(stream::IO, hfMultiplet::Hfs.HfMultiplet, role::String;
                                      threshold::Float64=1.0e-7, nAlways::Int64=3)
    ## THE THRESHOLD MUST BE SMALL, and the first version's 1.0e-3 was not. The two mechanisms this module treats
    ## differ by two orders of magnitude in their admixtures: nuclear hyperfine mixing in 229Th89+ gives 1.2e-2,
    ## whereas the electronic quenching of the Be-like 3P_0 level rests entirely on 1.8e-4 (33S12+, from the
    ## neighbouring 3P_1 at 0.59 eV). A threshold tuned to the first printed "-1.000000  <-- parent" for the
    ## second and hid the very admixture that makes the transition possible. Hence 1e-7 -- and, in addition, the
    ## nAlways largest contributions are printed whatever their size, because a composition table that shows only
    ## the parent has failed at its one job.
    ##
    ## What the Be-like table then shows is worth the space: 3P_1 contributes 1.8e-4, 1P_1 only 3e-6 (sixty
    ## times less, being ~25 eV away instead of 0.59 eV), and 3P_2 appears at 1e-6 solely because hfMultipoles
    ## includes E2, which connects J = 0 to J = 2 where M1 cannot. None of that is guessable from the input.
    nx = 104
    println(stream, " ")
    println(stream, "  Composition of the reported $role hyperfine levels  (|c| > $threshold, or among the " *
            "$nAlways largest):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    for  hfLevel in hfMultiplet.hfLevels
        sym = LevelSymmetry(hfLevel.F, hfLevel.parity)
        println(stream, "     hyperfine level " * string(hfLevel.index) * ":   F^P = " * string(sym) * ",  E = " *
                @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", hfLevel.energy)) * " " *
                TableStrings.inUnits("energy"))   ## absolute; differences are what matter
        perm = sortperm(abs.(hfLevel.mc), rev=true)
        for  (rank, p)  in  enumerate(perm)
            if  abs(hfLevel.mc[p]) <= threshold  &&  rank > nAlways    continue    end
            vec  = hfLevel.hfBasisVectors[p]
            jsym = LevelSymmetry(vec.levelJ.J, vec.levelJ.parity)
            isym = LevelSymmetry(vec.isomer.spinI, vec.isomer.parity)
            sa   = p == perm[1]  ?  "  <-- parent"  :  "  <-- admixture"
            println(stream, "        " * @sprintf("%+9.6f", hfLevel.mc[p]) *
                    "   electronic level " * TableStrings.center(4, string(vec.levelJ.index); na=1) *
                    " (J^P = " * TableStrings.center(7, string(jsym); na=1) * "),  nuclear I^P = " *
                    TableStrings.center(7, string(isym); na=1) *
                    (vec.isomer.energy != 0. ? "*" : " ") * sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`HyperfineInduced.displayAdmixtureRanking(stream::IO, hfMultiplet::Hfs.HfMultiplet, role::String;
                            nShow::Int64=20)`
    ... to list, ranked by size, every admixture channel that can carry a hyperfine-induced transition out of the
        given hyperfine multiplet, labelled as ELECTRONIC or NUCLEAR; nothing is returned.

        THIS IS THE TABLE THE OVERVIEW MODE EXISTS FOR. A hyperfine-induced transition happens only because the
        decaying level is not pure, so the admixture coefficients ARE the physics: their squares set the rate, and
        their ranking says which perturbers matter and which may be ignored.

        ELECTRONIC and NUCLEAR channels are ranked TOGETHER, deliberately. Both are the same thing -- a coupling
        matrix element over an energy denominator -- and treating them separately would make 229Th89+ look like a
        special case. It is not: there the electronic list is simply EMPTY, because the ion has a single 1s level
        and nothing electronic to mix with, and the whole effect is carried by one nuclear channel between the two
        isomers. A user reading that learns the true situation rather than meeting a fallback rule.

        The coefficients are taken from the DIAGONALIZED representation, so they are exact rather than the
        first-order V/dE estimate; dE is printed alongside, since a small denominator is what a user scanning for
        candidate perturbers is really looking for.
"""
function  displayAdmixtureRanking(stream::IO, hfMultiplet::Hfs.HfMultiplet, role::String; nShow::Int64=20)
    ## collect every non-parent contribution as one channel
    chans = NamedTuple{(:c,:lev,:kind,:from,:to,:dE),
                       Tuple{Float64,Int64,String,String,String,Float64}}[]
    for  hfLevel in hfMultiplet.hfLevels
        perm   = sortperm(abs.(hfLevel.mc), rev=true)
        parent = hfLevel.hfBasisVectors[perm[1]]
        for  p in perm[2:end]
            vec = hfLevel.hfBasisVectors[p]
            if  abs(hfLevel.mc[p]) < 1.0e-12    continue    end
            ## ELECTRONIC: the nuclear state is the same and another electronic level admixes.
            ## NUCLEAR:    the electronic level is the same and another nuclear state admixes.
            if      vec.isomer == parent.isomer  &&  vec.levelJ.index != parent.levelJ.index   kind = "electronic"
            elseif  vec.isomer != parent.isomer  &&  vec.levelJ.index == parent.levelJ.index   kind = "nuclear"
            else                                                                               kind = "mixed"
            end
            sFrom = "el " * string(parent.levelJ.index) * " (" *
                    string(LevelSymmetry(parent.levelJ.J, parent.levelJ.parity)) * "), nuc " *
                    string(LevelSymmetry(parent.isomer.spinI, parent.isomer.parity)) *
                    (parent.isomer.energy != 0. ? "*" : "")
            sTo   = "el " * string(vec.levelJ.index) * " (" *
                    string(LevelSymmetry(vec.levelJ.J, vec.levelJ.parity)) * "), nuc " *
                    string(LevelSymmetry(vec.isomer.spinI, vec.isomer.parity)) *
                    (vec.isomer.energy != 0. ? "*" : "")
            dE    = (vec.levelJ.energy + Defaults.convertUnits("energy: to atomic", vec.isomer.energy)) -
                    (parent.levelJ.energy + Defaults.convertUnits("energy: to atomic", parent.isomer.energy))
            ## DEGENERATE PAIRS ARE NOT CHANNELS, and must be dropped rather than ranked. If two basis vectors
            ## are exactly degenerate -- which happens whenever the multipoles in hfMultipoles produce no
            ## DIAGONAL hyperfine splitting, e.g. an E3-only calculation, or a nuclear state with mu = 0 -- then
            ## any superposition of them is an equally valid eigenvector and the eigensolver returns an
            ## arbitrary one. The tell-tale signature is dE = 0 with a coefficient near 1/sqrt(2) = 0.707. Such a
            ## pair carries no transition in any case, since a channel with dE = 0 emits no photon. The cure for
            ## the underlying degeneracy is to include M1 (and E2) in hfMultipoles alongside the mixing
            ## multipole, so that the hyperfine levels are split as they physically are.
            if  abs(dE) < 1.0e-12    continue    end
            push!(chans, (c = hfLevel.mc[p], lev = hfLevel.index, kind = kind, from = sFrom, to = sTo,
                          dE = Defaults.convertUnits("energy: from atomic", dE)))
        end
    end
    sort!(chans, by = ch -> -abs(ch.c))
    nx = 132
    println(stream, " ")
    println(stream, "  Admixture channels of the $role hyperfine levels, ranked -- these carry the transition:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "     coefficient    kind         hf level    parent                        admixed with" *
            "                  dE " * TableStrings.inUnits("energy"))
    println(stream, "  ", TableStrings.hLine(nx))
    if  length(chans) == 0
        println(stream, "     -- none: every hyperfine level is pure, so no hyperfine-induced transition is possible.")
    end
    for  (n, ch) in enumerate(chans)
        if  n > nShow
            println(stream, "     ... and $(length(chans)-nShow) further channels, all smaller.")
            break
        end
        println(stream, "     " * @sprintf("%+11.3e", ch.c) * "   " * TableStrings.flushleft(12, ch.kind; na=1) *
                TableStrings.center(8, string(ch.lev); na=3) * TableStrings.flushleft(30, ch.from; na=1) *
                TableStrings.flushleft(30, ch.to; na=1) * @sprintf("%+11.4e", ch.dE))
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> The rate scales as the SQUARE of these coefficients; a channel ten times smaller " *
            "contributes a hundred times less.")
    return( nothing )
end


"""
`HyperfineInduced.dominantNuclearState(level::Hfs.HfLevel)`
    ... to return the isomer that carries the largest weight in the given hyperfine level; a Nuclear.Isomer is
        returned.

        Every level of a hyperfine-mixed computation is a superposition of nuclear states, so strictly speaking
        no level "is" the isomer. It is nevertheless always dominated by one of them -- the admixture that drives
        the whole effect is typically of order a per cent -- and naming that one is what makes an output table
        readable.
"""
function  dominantNuclearState(level::Hfs.HfLevel)
    p = 1;   wmax = -1.
    for  (i, mc)  in  enumerate(level.mc)
        if  abs(mc) > wmax   wmax = abs(mc);   p = i    end
    end
    return( level.hfBasisVectors[p].isomer )
end


"""
`HyperfineInduced.nuclearStateLabel(level::Hfs.HfLevel)`
    ... to return a short label "I^P" for the nuclear state dominating the given hyperfine level, with an
        asterisk marking an excited isomer; a String is returned.

        WHY THIS IS NEEDED (06-Aug-2026): F alone does not identify a level once more than one nuclear state is
        in the basis. In H-like 229Th89+ the ground state (I = 5/2) gives F = 2, 3 and the isomer (I = 3/2) gives
        F = 1, 2, so F = 2 occurs TWICE -- and two physically quite different transitions, one that de-excites
        the nucleus and one that merely relaxes its hyperfine state, both printed as "1+ --> 2+". Tables that
        show only F^P are genuinely ambiguous and were misread as soon as they were produced.
"""
function  nuclearStateLabel(level::Hfs.HfLevel)
    iso = HyperfineInduced.dominantNuclearState(level)
    return( string(LevelSymmetry(iso.spinI, iso.parity)) * (iso.energy != 0. ? "*" : " ") )
end


"""
`HyperfineInduced.displayLines(stream::IO, lines::Array{HyperfineInduced.Line,1})`
    ... to list the selected hyperfine-induced lines and their channels before the computation; nothing is returned.
"""
function  displayLines(stream::IO, lines::Array{HyperfineInduced.Line,1})
    nx = 96
    println(stream, " ")
    println(stream, "  Selected hyperfine-induced lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(20, "F^P (i --> f)"; na=4);   sb = sb * TableStrings.center(20, "            "; na=4)
    sa = sa * TableStrings.center(18, "omega  " * TableStrings.inUnits("energy"); na=4)
    sb = sb * TableStrings.center(18, "                  "; na=4)
    sa = sa * TableStrings.flushleft(40, "multipoles and gauges"; na=2)
    sb = sb * TableStrings.flushleft(40, "                     "; na=2)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        isym = LevelSymmetry(line.initialLevel.F, line.initialLevel.parity)
        fsym = LevelSymmetry(line.finalLevel.F,   line.finalLevel.parity)
        sa   = "  " * TableStrings.center(20, string(isym) * "  -->  " * string(fsym); na=4)
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega)) * "      "
        mpg  = Tuple{EmMultipole,EmGauge}[]
        for  ch in line.channels    push!(mpg, (ch.multipole, ch.gauge))    end
        sa   = sa * TableStrings.multipoleGaugeTupels(50, mpg)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> $(length(lines)) hyperfine-induced lines will be computed.")
    return( nothing )
end


"""
`HyperfineInduced.displayRates(stream::IO, lines::Array{HyperfineInduced.Line,1},
                            settings::HyperfineInduced.Settings)`
    ... to list the transition rates of all computed hyperfine-induced lines; nothing is returned.
"""
function  displayRates(stream::IO, lines::Array{HyperfineInduced.Line,1}, settings::HyperfineInduced.Settings)
    nx = 140
    println(stream, " ")
    println(stream, "  Hyperfine-induced transition rates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(12, "level i-f"; na=2);       sb = sb * TableStrings.center(12, "         "; na=2)
    sa = sa * TableStrings.center(20, "F^P (i --> f)"; na=3);   sb = sb * TableStrings.center(20, "            "; na=3)
    sa = sa * TableStrings.center(22, "nuclear I^P (i --> f)"; na=3)
    sb = sb * TableStrings.center(22, "  (* = excited isomer)"; na=3)
    sa = sa * TableStrings.center(18, "omega  " * TableStrings.inUnits("energy"); na=3)
    sb = sb * TableStrings.center(18, "                  "; na=3)
    sa = sa * TableStrings.center(12, "multipole"; na=3);       sb = sb * TableStrings.center(12, "         "; na=3)
    sa = sa * TableStrings.center(30, "Cou -- rate -- Bab"; na=3)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("rate"); na=3)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        isym = LevelSymmetry(line.initialLevel.F, line.initialLevel.parity)
        fsym = LevelSymmetry(line.finalLevel.F,   line.finalLevel.parity)
        mps  = unique([ch.multipole for ch in line.channels])
        inuc = HyperfineInduced.nuclearStateLabel(line.initialLevel)
        fnuc = HyperfineInduced.nuclearStateLabel(line.finalLevel)
        sa   = "  " * TableStrings.center(12, string(line.initialLevel.index) * "-" *
                                              string(line.finalLevel.index); na=2)
        sa   = sa * TableStrings.center(20, string(isym) * "  -->  " * string(fsym); na=3)
        sa   = sa * TableStrings.center(22, inuc * " --> " * fnuc; na=3)
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega)) * "     "
        sa   = sa * TableStrings.center(12, join([string(mp) for mp in mps], ","); na=3)
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.photonRate.Coulomb))   * "   "
        sa   = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.photonRate.Babushkin))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`HyperfineInduced.displayLifetimes(stream::IO, lines::Array{HyperfineInduced.Line,1},
                            settings::HyperfineInduced.Settings)`
    ... to list, for every initial hyperfine level, the summed hyperfine-induced decay rate and the resulting
        lifetime; nothing is returned.

        This is the quantity the whole module exists to produce: for 229Th89+ it is the isomeric lifetime that
        nuclear hyperfine mixing shortens from hours to tens of milliseconds, and for 205Pb77+ the one that drops
        from about fifteen minutes to some tens of milliseconds. Note that only the hyperfine-induced channels
        computed here enter -- any competing decay path (internal conversion, ordinary electronic transitions) is
        NOT included, so a lifetime printed here is an upper bound on the true one.
"""
function  displayLifetimes(stream::IO, lines::Array{HyperfineInduced.Line,1}, settings::HyperfineInduced.Settings)
    ## Group by the initial hyperfine level, keyed on its INDEX. Until 06-Aug-2026 Hfs.HfLevel carried no index
    ## and this had to key on the surrogate pair (2F, energy) -- which is fragile precisely where it matters, since
    ## two levels of equal F differing only in their nuclear parent are exactly what this module computes.
    total = Dict{Int64,EmProperty}();   symOf = Dict{Int64,LevelSymmetry}();   nucOf = Dict{Int64,String}()
    for  line in lines
        k = line.initialLevel.index
        total[k] = get(total, k, EmProperty(0., 0.)) + line.photonRate
        symOf[k] = LevelSymmetry(line.initialLevel.F, line.initialLevel.parity)
        nucOf[k] = HyperfineInduced.nuclearStateLabel(line.initialLevel)
    end
    nx = 124
    println(stream, " ")
    println(stream, "  Hyperfine-induced decay rates and lifetimes of the initial levels:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(8, "level"; na=2);           sb = sb * TableStrings.center(8, "     "; na=2)
    sa = sa * TableStrings.center(14, "F^P"; na=4);            sb = sb * TableStrings.center(14, "   "; na=4)
    sa = sa * TableStrings.center(18, "nuclear I^P"; na=4)
    sb = sb * TableStrings.center(18, "(* = isomer)"; na=4)
    sa = sa * TableStrings.center(30, "Cou -- total rate -- Bab"; na=4)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("rate"); na=4)
    sa = sa * TableStrings.center(30, "Cou -- lifetime -- Bab"; na=2)
    sb = sb * TableStrings.center(30, "[s]"; na=2)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    for  k in sort(collect(keys(total)))
        rC = Defaults.convertUnits("rate: from atomic", total[k].Coulomb)
        rB = Defaults.convertUnits("rate: from atomic", total[k].Babushkin)
        sa = "  " * TableStrings.center(8, string(k); na=2)
        sa = sa * TableStrings.center(14, string(symOf[k]); na=4)
        sa = sa * TableStrings.center(18, nucOf[k]; na=4)
        sa = sa * @sprintf("%.6e", rC) * "   " * @sprintf("%.6e", rB) * "      "
        sa = sa * (rC > 0. ? @sprintf("%.6e", 1/rC) : "     ---    ") * "   "
        sa = sa * (rB > 0. ? @sprintf("%.6e", 1/rB) : "     ---    ")
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> Only the hyperfine-induced channels computed here are summed; competing decay paths are not.")
    return( nothing )
end

end # module
