
"""
`module  JAC.DielectronicRecombination`  
... a submodel of JAC that contains all methods for computing dielectronic recombination strength & rate coefficients 
    for some given initial, intermediate and final-state multiplets. Special code has been implemented to account for
    the capture into high-n Rydberg shells by incorporating, in addition to the explicitly calculated Auger and radiative
    rates, also several hydrogenic and empirical corrections to the DR resonance strength. Further corrections can be
    readily added if the needs arise.

    DR is represented by the two physical steps it consists of: a CaptureLine holds the capture i + e- --> m together
    with the two total widths Gamma_a(m) and Gamma_r(m) that decide its fate, and a PhotonLine holds one radiative
    stabilization m --> f + hv. Combining them gives both the resonance strengths and the DR satellite spectrum.

    The shells that cannot be treated explicitly are covered by corrections, see
    DielectronicRecombination.AbstractCorrections: the radiative decay into final shells above n^(final) by
    HydrogenicCorrections, and the capture into Rydberg shells above n^(lowest-captured) by RydbergTailCorrection.

    Before any dielectronic-recombination calculations are done, it is generally suggested to carefully check the given
    lists of initial, intermediate and final configuration since missing levels in these lists are "hard" to detect 
    automatically.
"""
module DielectronicRecombination


using Base.Threads, Distributed, Printf, ProgressMeter, SpecialFunctions,
        ..AngularMomentum, ..AutoIonization, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..Hfs, ..ManyElectron,
        ..Nuclear, ..PhotoEmission, ..Radial, ..TableStrings


"""
`abstract type DielectronicRecombination.AbstractCorrections` 
    ... defines an abstract type to distinguish different types of corrections to the decay rates and strength.
        These corrections are based on the classification of shells:

        n^(core)  <   n^(final)  <  n^(hydrogenic)  <  n^(lowest-captured)  <=  n^(max)      ... where

        n^(core)            ... refers to the (maximum) principal quantum number to which initial core electrons are excited;
        n^(final)           ... the maximum number for which shells are treated explicitly in the representation of the final
                                levels f, EXCLUDING the shell the electron was captured into (which is always present in the
                                final basis and would otherwise swallow this number);
        n^(hydrogenic)      ... to the maximum n-shell, to which the radiative decay is modeled by scaled-hydrogenic rates,
                            ... and which can be omitted also from the list.
        n^(lowest-captured) ... is the lowest, high-n shell, into which the additional electron is captured and which must
                                (of course) occur explicitly in the basis of the intermediate and final levels.
        n^(max)             ... the highest shell for which the DR resonances are still estimated, by extrapolating the
                                explicitly computed Rydberg series; all shells with n > n^(max) are neglected completely.

    + struct DielectronicRecombination.RydbergTailCorrection
        ... to extrapolate the explicitly computed Rydberg series to the shells n^(lowest-captured) < n <= n^(max) that are
            too numerous to be treated explicitly, and -- optionally -- to the orbital angular momenta l that exist at those
            shells but were never computed. The capture and Auger rates are scaled as n^(-3), the core radiative rate is
            held n-independent, and the radiative decay of the Rydberg spectator itself is computed (not scaled) by the
            hydrogenic machinery. Replaces the former EmpiricalCorrections, whose scalings could not be reconstructed.
    + struct DielectronicRecombination.HydrogenicCorrections
        ... to add for missing final decay levels to the (total) photon decay rates by scaling the corresponding rates
            of non-relativistic hydrogenic ions with a suitable effective charge (Zeff); these hydrogenic corrections improve
            goth, the total photon rate as well as the resonance strength.
    + struct DielectronicRecombination.MaximumlCorrection
        ... to exclude all subshells with l > l_max in the hydrogenic corrections; this restriction does not apply to the
            given resonance levels, which can be controlled (and are specified) by the list of intermediate configurations.
    + struct DielectronicRecombination.ResonanceWindowCorrection
        ... to keep only those resonances whose capture energy falls inside a given window [E_min, E_max], for instance the
            range actually scanned in an EBIT or storage-ring measurement.
"""
abstract type  AbstractCorrections       end


"""
`struct  DielectronicRecombination.RydbergTailCorrection  <:  DielectronicRecombination.AbstractCorrections`
    ... to estimate the DR resonances of the shells  n^(lowest-captured) < n <= nMax  by extrapolating the Rydberg
        series that WAS computed explicitly, and -- optionally -- the orbital angular momenta that exist at those
        shells but were never computed. This matters because a Rydberg series does not stop where one runs out of
        computer time: while the resonances remain Auger dominated their strength per (n,l) is nearly independent
        of n, so the shells just above the explicit ones can carry a large part of alpha(T).

        The extrapolation rests on ONE scaling assumption, A_a ~ n^(-3), which follows from the normalization of a
        Rydberg orbital near the core and which this module VERIFIES against the user's own data whenever two
        Rydberg shells were computed explicitly. Everything else is either exact (the Rydberg energy shift) or
        computed rather than scaled (the radiative decay of the spectator, via HydrogenicCorrections).

    + nMax          ::Int64
        ... highest principal quantum number included in the tail sum; shells above it are neglected completely.
            Physically the series is cut by field ionization or by collisions in the plasma, so this is a property
            of the experiment and not of the atom.
    + effectiveZ    ::Union{Float64,Missing}
        ... effective charge seen by the Rydberg electron; derived as Z - n^(core) if missing.
    + lMax          ::Union{Int64,Missing}
        ... highest l to extrapolate to at each shell. The l values present at n^(lowest-captured) are extrapolated
            in n in any case; this field additionally extrapolates in l, up to lMax or until the contribution has
            fallen below 1e-4 of the l = 0 term. Set to a value <= the largest computed l, or leave missing, to
            suppress the l extrapolation entirely.
    + nExponent     ::Union{Float64,Missing}
        ... exponent of the n-scaling of the capture and Auger rates; 3.0 if missing.

        SETTLED 02-Sep-2026, AND THE DEFAULT 3.0 IS RIGHT. Measured on H-like C5+ (the system of
        examples/example-Df.jl branches g and h) with each Rydberg shell computed in a basis OF ITS OWN, so that
        no two shells can mix:

            l = 1, the only l whose level content is the same at every shell (10 levels at n = 8, 9 and 10):
                 p(8->9) = 3.066     p(9->10) = 3.064     p(8->10) = 3.065

        Three independent pairs agreeing to +-0.001. The residual 2 % above 3 is small, systematic and not fully
        explained -- a quantum defect of the right sign accounts for at most 3.01 -- and it costs the tail sum
        (n = 10..30 scaled from n0 = 9) only 2.7 %, which is far inside any other uncertainty a DR rate
        coefficient carries. So `nExponent = 3.0` stands.

        WHAT THE EARLIER LOW VALUES WERE. Branch h reported p = 1.27 for l = 1, and that number is an ARTIFACT OF
        PUTTING n = 8 AND n = 9 IN ONE CI BASIS, not a property of the law. Measured both ways on the same system:
        W(8, l=1) falls 8.137e-4 -> 7.296e-4 (-10.3 %) and W(9, l=1) rises 5.671e-4 -> 6.279e-4 (+10.7 %) when the
        two shells share a basis, while their SUM is conserved to 1.7 % and 72 % of what n = 8 loses reappears at
        n = 9. That is redistribution between the shells, and the ratio of two shells is exactly what the exponent
        is formed from -- so a shared basis flattens it towards zero. See the warning in `measureRydbergExponent`.
"""
struct   RydbergTailCorrection               <:  DielectronicRecombination.AbstractCorrections
    nMax              ::Int64
    effectiveZ        ::Union{Float64,Missing}
    lMax              ::Union{Int64,Missing}
    nExponent         ::Union{Float64,Missing}
end


"""
`Base.show(io::IO, corr::RydbergTailCorrection)`
    ... prepares a proper printout of the corr::RydbergTailCorrection; nothing is returned.
"""
function Base.show(io::IO, corr::RydbergTailCorrection)
    println(io, "RydbergTailCorrection() with: ")
    println(io, "nMax:            $(corr.nMax)  ")
    println(io, "effectiveZ:      $(corr.effectiveZ)  ")
    println(io, "lMax:            $(corr.lMax)  ")
    println(io, "nExponent:       $(corr.nExponent)  ")
end


"""
`struct  DielectronicRecombination.HydrogenicCorrections  <:  DielectronicRecombination.AbstractCorrections`  
    ... to add for missing final decay levels the photon decay rates for non-relativistic hydrogenic ions;
        this improves the total photon rate as well as the resonance strength. These corrections are taken into
        account for all shells with n^{final}+1 <= n <= nHydrogenic

    + nHydrogenic       ::Union{Int64,Missing}   
        ... upper principal quantum number nHydrogenic for which hydrogenic corrections to the radiative photon rates are 
            calculated explicitly; the photon rates are further scaled if some proper effectiveZ and/or rateScaling
            is provided.
    + effectiveZ      ::Union{Float64,Missing}   ... effective charge Z_eff for the hydrogenic correction.
    + rateScaling     ::Union{Float64,Missing}   ... scaling factor to scale the photon rates
"""
struct   HydrogenicCorrections               <:  DielectronicRecombination.AbstractCorrections
    nHydrogenic       ::Union{Int64,Missing}  
    effectiveZ        ::Union{Float64,Missing}
    rateScaling       ::Union{Float64,Missing}
end


"""
`Base.show(io::IO, corr::HydrogenicCorrections)`
    ... prepares a proper printout of the corr::HydrogenicCorrections; nothing is returned.
"""
function Base.show(io::IO, corr::HydrogenicCorrections)
    println(io, "HydrogenicCorrections() with: ")
    println(io, "nHydrogenic:     $(corr.nHydrogenic)  ")
    println(io, "effectiveZ:      $(corr.effectiveZ)  ")
    println(io, "rateScaling:     $(corr.rateScaling)  ")
end     


"""
`struct  DielectronicRecombination.MaximumlCorrection  <:  DielectronicRecombination.AbstractCorrections`  
    ... to exclude all subshells with l > l_max, both in the treatment of the corrections shells.

    + maximum_l    ::Union{Int64,Missing}   
        ... maximum orbital angular momentum quantum number for which contributions to the DR strengths are 
            taken into account. This number applies for all subshells for which other corrections are 
            requested, whereas the "physical subshells" are defined by the configuration lists.
"""
struct   MaximumlCorrection                  <:  DielectronicRecombination.AbstractCorrections
    maximum_l      ::Union{Int64,Missing} 
end


"""
`Base.show(io::IO, corr::MaximumlCorrection)`
    ... prepares a proper printout of the corr::MaximumlCorrection; nothing is returned.
"""
function Base.show(io::IO, corr::MaximumlCorrection)
    println(io, "MaximumlCorrection(lmax = $(corr.maximum_l)): ")
end     


"""
`struct  DielectronicRecombination.ResonanceWindowCorrection  <:  DielectronicRecombination.AbstractCorrections`  
    ... to exclude all DR resonances outside of a given "window [E_min, E_max]" of resonance energies with
        regard to the initial level.

    + energyMin  ::Float64   ... minimum energy [Hartree] of the resonances to be considered.  
    + energyMax  ::Float64   ... maximum energy [Hartree] of the resonances to be considered.   
"""
struct   ResonanceWindowCorrection           <:  DielectronicRecombination.AbstractCorrections
    energyMin    ::Float64  
    energyMax    ::Float64   
end


"""
`Base.show(io::IO, corr::ResonanceWindowCorrection)`
    ... prepares a proper printout of the corr::ResonanceWindowCorrection; nothing is returned.
"""
function Base.show(io::IO, corr::ResonanceWindowCorrection)
    println(io, "ResonanceWindowCorrection() with: ")
    println(io, "energyMin:  $(corr.energyMin)  ")
    println(io, "energyMax:  $(corr.energyMax)  ")
end     



"""
`struct  DielectronicRecombination.EmpiricalTreatment`  
    ... defines an (internal) type to communicate and distribute the physical (and technical) parameters
        that are utilized to make the requested empirical corrections or just nothing. This data type should
        not be applied by the user but is initialized by the given (set of) corrections.
        Otherwise, it is treated like any other type in JAC. All parameters are made physically "explicit",
        even if they were "missing" originally, and can be directly applied in the empirical treatment of
        the DR process. The following hierarchy of shells is used:

        n^(core)  <   n^(final)  <  n^(hydrogenic)  <  n^(lowest-captured)  <=  n^(max)

    + doRydbergTailCorrection     ::Bool    ... True, if the Rydberg series is to be extrapolated, false o/w.
    + doHydrogenicCorrections     ::Bool    ... True, if hydrogenic corrections are needed, false o/w.
    + doMaximumlCorrection        ::Bool    ... True, if a maximum l values is used, false o/w.
    + doResonanceWindowCorrection ::Bool    ... True, if a window of resonances is specified, false o/w.
    + nCore                       ::Int64
        ... (maximum) principal quantum number to which initial core electrons are excited;
    + nFinal                      ::Int64
        ... the maximum number for which shells are treated explicitly in the representation of the final levels f,
            EXCLUDING the shell the electron was captured into.
    + nHydrogenic                 ::Int64
        ... maximum n-shell, to which the radiative decay is modeled by scaled-hydrogenic rates.
    + nLowestCaptured             ::Int64
        ... lowest, high-n shell, into which the additional electron is captured and which must (of course) occur
            explicitly in the basis of the intermediate and final levels.
    + nMax                        ::Int64
        ... highest n-shell for which the DR resonances are still estimated by extrapolating the explicitly computed
            Rydberg series; all shells with n > nMax are neglected completely.
    + maximum_l                   ::Int64    ... maximum l value; is set to a large value if not specified by the user.
    + lMaxTail                    ::Int64
        ... highest l to extrapolate to in the tail; equal to the largest computed l when no l extrapolation is wanted.
    + nExponent                   ::Float64  ... exponent of the n-scaling of the capture and Auger rates.
    + hydrogenicEffectiveZ        ::Float64  ... effective charge Z_eff for the hydrogenic correction.
    + hydrogenicRateScaling       ::Float64  ... scaling factor to modify the estimated hydrogenic rates.
    + tailEffectiveZ              ::Float64  ... effective Z seen by the Rydberg electron in the tail extrapolation.
    + resonanceEnergyMin:         ::Float64  ... minimum energy [Hartree] of the resonances to be considered.
    + resonanceEnergyMax:         ::Float64  ... maximum energy [Hartree] of the resonances to be considered.
"""
struct   EmpiricalTreatment
    doRydbergTailCorrection     ::Bool
    doHydrogenicCorrections     ::Bool
    doMaximumlCorrection        ::Bool
    doResonanceWindowCorrection ::Bool
    nCore                       ::Int64
    nFinal                      ::Int64
    nHydrogenic                 ::Int64
    nLowestCaptured             ::Int64
    nMax                        ::Int64
    maximum_l                   ::Int64
    lMaxTail                    ::Int64
    nExponent                   ::Float64
    hydrogenicEffectiveZ        ::Float64
    hydrogenicRateScaling       ::Float64
    tailEffectiveZ              ::Float64
    resonanceEnergyMin          ::Float64
    resonanceEnergyMax          ::Float64
end


"""
`Base.show(io::IO, tr::EmpiricalTreatment)`
    ... prepares a proper printout of the tr::EmpiricalTreatment; nothing is returned.
"""
function Base.show(io::IO, tr::EmpiricalTreatment)
    println(io, "doRydbergTailCorrection:  $(tr.doRydbergTailCorrection)  ")
    println(io, "doHydrogenicCorrections:  $(tr.doHydrogenicCorrections)  ")
    println(io, "doMaximumlCorrection:     $(tr.doMaximumlCorrection)  ")
    println(io, "doResonanceWindowCorrect: $(tr.doResonanceWindowCorrection)  ")
    println(io, "nCore:                    $(tr.nCore)  ")
    println(io, "nFinal:                   $(tr.nFinal)  ")
    println(io, "nHydrogenic:              $(tr.nHydrogenic)  ")
    println(io, "nLowestCaptured:          $(tr.nLowestCaptured)  ")
    println(io, "nMax:                     $(tr.nMax)  ")
    println(io, "maximum_l:                $(tr.maximum_l)  ")
    println(io, "lMaxTail:                 $(tr.lMaxTail)  ")
    println(io, "nExponent:                $(tr.nExponent)  ")
    println(io, "hydrogenicEffectiveZ:     $(tr.hydrogenicEffectiveZ)  ")
    println(io, "hydrogenicRateScaling:    $(tr.hydrogenicRateScaling)  ")
    println(io, "tailEffectiveZ:           $(tr.tailEffectiveZ)  ")
    println(io, "resonanceEnergyMin:       $(tr.resonanceEnergyMin)  ")
    println(io, "resonanceEnergyMax:       $(tr.resonanceEnergyMax)  ")
end


"""
`struct  DielectronicRecombination.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing dielectronic recombination pathways.

    + multipoles            ::Array{EmMultipoles}  ... Multipoles of the radiation field that are to be included.
    + gauges                ::Array{UseGauge}      ... Specifies the gauges to be included into the computations.
    + calcOnlyPassages      ::Bool                 
        ... Only compute resonance strength but without making all the pathways explicit. This option is useful
            for the capture into high-n shells or if the photons are not considered explicit. It also treats the 
            shells differently due to the given core shells < final-state shells < hydrogenically-scaled shells <
            capture-shells < asymptotic-shells. Various correction and multi-threading techniques can be applied
            to deal with or omit different classes of these shells.
    + calcRateAlpha         ::Bool                 
        ... True, if the DR rate coefficients are to be calculated, and false o/w.
    + calcHyperfineResolved ::Bool                 
        ... True, if the DR resonance strength are calculated for hyperfine-resolved levels, and false o/w.
            If true, it need to come together with calcOnlyPassages = true, and no fine-structure resolved rates
            and strength are computed in this case.
    + calcPhotonSpectrum    ::Bool
        ... True, if the individual PhotonLine's (m --> f) are to be RETAINED and displayed, i.e. if the DR satellite
            spectrum is wanted, and false o/w. Note that this flag controls RETENTION, not computation: the total
            radiative width Gamma_r(m) = sum_f A_r(m,f) is needed for every resonance strength, so each m --> f rate
            is evaluated in either case. Setting it false accumulates those rates into Gamma_r and discards them,
            which saves memory (the point for large Rydberg manifolds) but essentially no computing time.
            Only used by the fine-structure resolved route.
    + printBefore           ::Bool
        ... True, if all energies and pathways are printed before their evaluation.
    + pathwaySelection      ::PathwaySelection     ... Specifies the selected levels/pathways, if any.
    + electronEnergyShift   ::Float64              
        ... An overall energy shift for all electron energies (i.e. from the initial to the resonance levels [Hartree].
    + photonEnergyShift     ::Float64              
        ... An overall energy shift for all photon energies (i.e. from the resonance to the final levels.
    + mimimumPhotonEnergy   ::Float64              
        ... minimum transition energy for which photon transitions are  included into the evaluation.
    + temperatures          ::Array{Float64,1}     
        ... list of temperatures for which plasma rate coefficients are displayed; however, these rate coefficients
            only include the contributions from those pathways that are calculated here explicitly.
    + corrections           ::Array{DielectronicRecombination.AbstractCorrections,1}
        ... Specify, if appropriate, the inclusion of additional corrections to the rates and DR strengths.
    + augerOperator         ::AbstractEeInteraction 
        ... Auger operator that is to be used for evaluating the Auger amplitude's; the allowed values are: 
            CoulombInteraction(), BreitInteration(), CoulombBreit(), CoulombGaunt().
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles              ::Array{EmMultipole,1}
    gauges                  ::Array{UseGauge}
    calcOnlyPassages        ::Bool
    calcRateAlpha           ::Bool
    calcHyperfineResolved   ::Bool
    calcPhotonSpectrum      ::Bool
    printBefore             ::Bool
    pathwaySelection        ::PathwaySelection
    electronEnergyShift     ::Float64
    photonEnergyShift       ::Float64
    mimimumPhotonEnergy     ::Float64
    temperatures            ::Array{Float64,1}
    corrections             ::Array{DielectronicRecombination.AbstractCorrections,1}
    augerOperator           ::AbstractEeInteraction
end 


"""
`DielectronicRecombination.Settings()`  
    ... constructor for the default values of dielectronic recombination pathway computations.
"""
function Settings()
    Settings([E1], UseGauge[], false, false, false, false, false, PathwaySelection(), 0., 0., 0., Float64[],
             DielectronicRecombination.AbstractCorrections[], CoulombInteraction())
end


"""
` (set::DielectronicRecombination.Settings;`

        multipoles=..,             gauges=..,                  
        calcOnlyPassages=..,       calcRateAlpha=..,         calcHyperfineResolved=..,
        calcPhotonSpectrum=..,
        printBefore=..,            pathwaySelection=..,      electronEnergyShift=..,   photonEnergyShift=..,       
        mimimumPhotonEnergy=..,    temperatures=..,          corrections=..,           augerOperator=..)
                    
    ... constructor for modifying the given DielectronicRecombination.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::DielectronicRecombination.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,               gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    calcOnlyPassages::Union{Nothing,Bool}=nothing,                         calcRateAlpha::Union{Nothing,Bool}=nothing,  
    calcHyperfineResolved::Union{Nothing,Bool}=nothing,
    calcPhotonSpectrum::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                              pathwaySelection::Union{Nothing,PathwaySelection}=nothing,
    electronEnergyShift::Union{Nothing,Float64}=nothing,                   photonEnergyShift::Union{Nothing,Float64}=nothing, 
    mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,                   temperatures::Union{Nothing,Array{Float64,1}}=nothing,    
    corrections::Union{Nothing,Array{AbstractCorrections,1}}=nothing,      augerOperator::Union{Nothing,AbstractEeInteraction}=nothing)
    
    if  isnothing(multipoles)            multipolesx           = set.multipoles            else  multipolesx           = multipoles            end 
    if  isnothing(gauges)                gaugesx               = set.gauges                else  gaugesx               = gauges                end 
    if  isnothing(calcOnlyPassages)      calcOnlyPassagesx     = set.calcOnlyPassages      else  calcOnlyPassagesx     = calcOnlyPassages      end 
    if  isnothing(calcRateAlpha)         calcRateAlphax        = set.calcRateAlpha         else  calcRateAlphax        = calcRateAlpha         end 
    if  isnothing(calcHyperfineResolved) calcHyperfineResolvedx= set.calcHyperfineResolved else  calcHyperfineResolvedx= calcHyperfineResolved end
    if  isnothing(calcPhotonSpectrum)    calcPhotonSpectrumx   = set.calcPhotonSpectrum    else  calcPhotonSpectrumx   = calcPhotonSpectrum    end 
    if  isnothing(printBefore)           printBeforex          = set.printBefore           else  printBeforex          = printBefore           end 
    if  isnothing(pathwaySelection)      pathwaySelectionx     = set.pathwaySelection      else  pathwaySelectionx     = pathwaySelection      end 
    if  isnothing(electronEnergyShift)   electronEnergyShiftx  = set.electronEnergyShift   else  electronEnergyShiftx  = electronEnergyShift   end 
    if  isnothing(photonEnergyShift)     photonEnergyShiftx    = set.photonEnergyShift     else  photonEnergyShiftx    = photonEnergyShift     end 
    if  isnothing(mimimumPhotonEnergy)   mimimumPhotonEnergyx  = set.mimimumPhotonEnergy   else  mimimumPhotonEnergyx  = mimimumPhotonEnergy   end 
    if  isnothing(temperatures)          temperaturesx         = set.temperatures          else  temperaturesx         = temperatures          end 
    if  isnothing(corrections)           correctionsx          = set.corrections           else  correctionsx          = corrections           end 
    if  isnothing(augerOperator)         augerOperatorx        = set.augerOperator         else  augerOperatorx        = augerOperator         end 

    Settings( multipolesx, gaugesx, calcOnlyPassagesx, calcRateAlphax, calcHyperfineResolvedx, calcPhotonSpectrumx, printBeforex,
              pathwaySelectionx, electronEnergyShiftx, photonEnergyShiftx, mimimumPhotonEnergyx, temperaturesx,
              correctionsx, augerOperatorx )
end


"""
`Base.show(io::IO, settings::DielectronicRecombination.Settings)`
    ... prepares a proper printout of the variable settings::DielectronicRecombination.Settings; nothing is returned.
"""
function Base.show(io::IO, settings::DielectronicRecombination.Settings) 
    println(io, "multipoles:                 $(settings.multipoles)  ")
    println(io, "use-gauges:                 $(settings.gauges)  ")
    println(io, "calcOnlyPassages:           $(settings.calcOnlyPassages)  ")
    println(io, "calcRateAlpha:              $(settings.calcRateAlpha)  ")
    println(io, "calcHyperfineResolved:      $(settings.calcHyperfineResolved)  ")
    println(io, "calcPhotonSpectrum:         $(settings.calcPhotonSpectrum)  ")
    println(io, "printBefore:                $(settings.printBefore)  ")
    println(io, "pathwaySelection:           $(settings.pathwaySelection)  ")
    println(io, "electronEnergyShift:        $(settings.electronEnergyShift)  ")
    println(io, "photonEnergyShift:          $(settings.photonEnergyShift)  ")
    println(io, "mimimumPhotonEnergy:        $(settings.mimimumPhotonEnergy)  ")
    println(io, "temperatures:               $(settings.temperatures)  ")
    println(io, "corrections:                $(settings.corrections)  ")
    println(io, "augerOperator:              $(settings.augerOperator)  ")
end



"""
`struct  DielectronicRecombination.ResonanceSelection`  
    ... defines a type for selecting classes of resonances in terms of leading configurations.

    + active          ::Bool              ... initial-(state) level
    + fromShells      ::Array{Shell,1}    ... List of shells from which excitations are to be considered.
    + toShells        ::Array{Shell,1}    ... List of shells to which (core-shell) excitations are to be considered.
    + intoShells      ::Array{Shell,1}    ... List of shells into which electrons are initially placed (captured).
"""
struct  ResonanceSelection
    active            ::Bool  
    fromShells        ::Array{Shell,1} 
    toShells          ::Array{Shell,1} 
    intoShells        ::Array{Shell,1} 
end 


"""
`DielectronicRecombination.ResonanceSelection()`  
    ... constructor for an 'empty' instance of a ResonanceSelection()
"""
function ResonanceSelection()
    ResonanceSelection(false, Shell[], Shell[], Shell[] )
end


"""
`Base.show(io::IO, rSelection::DielectronicRecombination.ResonanceSelection)`
    ... prepares a proper printout of resonance::DielectronicRecombination.ResonanceSelection; nothing is returned.
"""
function Base.show(io::IO, rSelection::DielectronicRecombination.ResonanceSelection) 
    println(io, "active:           $(rSelection.active)  ")
    println(io, "fromShells:       $(rSelection.fromShells)  ")
    println(io, "toShells:         $(rSelection.toShells)  ")
    println(io, "intoShells:       $(rSelection.intoShells)  ")
end


include("module-DielectronicRecombination-inc-FS-resolved.jl")
include("module-DielectronicRecombination-inc-HF-resolved.jl")

end # module
