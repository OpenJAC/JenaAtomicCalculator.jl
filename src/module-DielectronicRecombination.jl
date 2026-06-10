
"""
`module  JAC.DielectronicRecombination`  
... a submodel of JAC that contains all methods for computing dielectronic recombination strength & rate coefficients 
    for some given initial, intermediate and final-state multiplets. Special code has been implemented to account for 
    the capture into high-n Rydberg shells by incorporating, in addition to the explicitly calculated Auger and radiative
    rates, also several hyrogenic and empirical corrections to the DR resonance strength. Further corrections can be 
    readily added if the needs arise.
        
    The program distinguishes between the prior set-up of all pathways (i -- m -- f) as well as of passages (i --m).
    For passages, the radiative stabilization into any lower-lying level is calculated "on fly" by dividing the 
    possible final levels into two groups: (i) The explicitly calculated final levels with a principal quantum number
    n <= nFinal and (ii) a set of final levels with principal quantum numbers nFinal < n <= nHydrogenic < nLowestCaptured.
    See DielectronicRecombination.AbstractCorrections for possible hydrogenic, empirical and other corrections that 
    can be taken into account.
    
    Before any dielectronic-recombination calculations are done, it is generally suggested to carefully check the given
    lists of initial, intermediate and final configuration since missing levels in these lists are "hard" to detect 
    automatically.
"""
module DielectronicRecombination


using Base.Threads, Distributed, Printf, ProgressMeter, SpecialFunctions,
        ..AngularMomentum, ..AutoIonization, ..Basics, ..Continuum, ..Defaults, ..Hfs, ..ManyElectron, ..Nuclear, 
        ..PhotoEmission, ..Radial, ..TableStrings


"""
`abstract type DielectronicRecombination.AbstractCorrections` 
    ... defines an abstract type to distinguish different types of corrections to the decay rates and strength.
        These corrections are based on the classification of shell:
        
        n^(core)  <   n^(final)  <  n^(hydrogenic)  <  n^(lowest-captured)  <  n^(lower-empirical)    
                  <=  n^(upper-empirical)              ... where
                  
        n^(core)            ... refers to the (maximum) principal quantum number to which initial core electrons are excited;
        n^(final)           ... the maximum number for which shells are treated explicitly in the representation of the final levels f;
        n^(hydrogenic)      ... to the maximum n-shell, to which the radiative decay is modeled by scaled-hydrogenic rates, 
                            ... and which can be omitted also from the list. 
        n^(lowest-captured) ... is the lowest, high-n shell, into which the additional electron is captured and which must
                                (of course) occur explicitly in the basis of the intermediate and final levels. 
        [n^(lower-empirical)  <=  n^(upper-empirical)]  
                            ... designates additional (empirical) high-n shells for which the contributions to the DR resonances 
                                are still estimated empirically by using arguments from quantum-defect theory. All shells with 
                                n > n^(upper-empirical) are neglected completely for their contributions to the DR spectra; see also:
    
    + struct DielectronicRecombination.EmpiricalCorrections  
        ... to estimate empirically the contributions of additional resonances for the capture of an electron into
            shells with [n^(lower-empirical)  <=  n^(upper-empirical)]. A simple scaling of the rates, calculated initially
            for n^(lowest-captured), ... only, is utilized for estimating the associated strength for these additional
            resonances.
    + struct DielectronicRecombination.HydrogenicCorrections  
        ... to add for missing final decay levels to the (total) photon decay rates by scaling the corresponding rates
            of non-relativistic hydrogenic ions with a suitable effective charge (Zeff); these hydrogenic corrections improve
            goth, the total photon rate as well as the resonance strength.
    + struct DielectronicRecombination.MaximumlCorrection  
        ... to exclude all subshells with l > l_max in the hydrogenic corrections; this restriction does not apply to the 
            given resonance levels, which can be controlled (and are specified) by the list of intermediate configurations.
"""
abstract type  AbstractCorrections       end


"""
`struct  DielectronicRecombination.EmpiricalCorrections  <:  DielectronicRecombination.AbstractCorrections`  
    ... to include empirical corrections for the shells with [n^(lower-empirical)  <=  n^(upper-empirical)].
        A rather rude model is used so far.

    + nUpperEmpirical ::Union{Int64,Missing}   
        ... The upper-empirical shell for which rate contributions are estimated; the lower-empirical shell = n^(captured-max + 1) 
            is derived from the given configuration lists. No corrections are made for nUpperEmpirical <= n^(captured-max + 1).
    + effectiveZ      ::Union{Float64,Missing}  ... effective charge Z_eff for the hydrogenic correction (inactive).
    + rateScaling     ::Union{Float64,Missing}  ... scaling factor to modify the estimated rates.
"""
struct   EmpiricalCorrections                <:  DielectronicRecombination.AbstractCorrections
    nUpperEmpirical   ::Union{Int64,Missing}
    effectiveZ        ::Union{Float64,Missing}
    rateScaling       ::Union{Float64,Missing} 
end


# `Base.show(io::IO, corr::EmpiricalCorrections)`  ... prepares a proper printout of the corr::EmpiricalCorrections.
function Base.show(io::IO, corr::EmpiricalCorrections)
    println(io, "EmpiricalCorrections() with: ")
    println(io, "nUpperEmpirical: $(corr.nUpperEmpirical)  ")
    println(io, "effectiveZ:      $(corr.effectiveZ)  ")
    println(io, "rateScaling:     $(corr.rateScaling)  ")
end     


"""
`struct  DielectronicRecombination.HydrogenicCorrections  <:  DielectronicRecombination.AbstractCorrections`  
    ... to add for missing final decay levels the photon decay rates for non-relativistic hydrogenic ions;
        this improves the total photon rate as well as the resonance strength. These corrections are taken into
        account for all shells with n^{final}+1 <= n <= nHydrogenic

    + nHydrogenic       ::Union{Int64,Missing}   
        ... upper principal quantum number nHydrogenic for which hydrogenic correctios to the radiative photon rates are 
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


# `Base.show(io::IO, corr::HydrogenicCorrections)`  ... prepares a proper printout of the corr::HydrogenicCorrections.
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


# `Base.show(io::IO, corr::MaximumlCorrection)`  ... prepares a proper printout of the corr::MaximumlCorrection.
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


# `Base.show(io::IO, corr::ResonanceWindowCorrection)`  ... prepares a proper printout of the corr::ResonanceWindowCorrection.
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
        
        n^(core)  <   n^(final)  <  n^(hydrogenic)  <  n^(lowest-captured)  <  n^(lower-empirical)    
                  <=  n^(upper-empirical) 
        
    + doEmpiricalCorrections      ::Bool    ... True, if empirical corrections are needed, false o/w.
    + doHydrogenicCorrections     ::Bool    ... True, if hydrogenic corrections are needed, false o/w.
    + doMaximumlCorrection        ::Bool    ... True, if a maximum l values is used, false o/w.
    + doResonanceWindowCorrection ::Bool    ... True, if a window of resonances is specified, false o/w.
    + nCore                       ::Int64   
        ... (maximum) principal quantum number to which initial core electrons are excited;
    + nFinal                      ::Int64   
        ... the maximum number for which shells are treated explicitly in the representation of the final levels f;
    + nHydrogenic                 ::Int64   
        ... maximum n-shell, to which the radiative decay is modeled by scaled-hydrogenic rates. 
    + nLowestCaptured             ::Int64   
        ... lowest, high-n shell, into which the additional electron is captured and which must (of course) occur 
            explicitly in the basis of the intermediate and final levels.
    + nLowerEmpirical             ::Int64   
        ... maximum n-shell, to which the radiative decay is modeled by scaled-hydrogenic rates. 
    + nUpperEmpirical             ::Int64   
        ... additional (empirical) high-n shells for which the contributions to the DR resonances are still 
            estimated empirically by using arguments from quantum-defect theory.
    + maximum_l                   ::Int64    ... maximum l value; is set to a large value if not specified by the user.
    + hydrogenicEffectiveZ        ::Float64  ... effective charge Z_eff for the hydrogenic correction (inactive).
    + hydrogenicRateScaling       ::Float64  ... scaling factor to modify the estimated hydrogenic rates.
    + empiricalEffectiveZ         ::Float64  ... effective Z for empirical estimates
    + empiricalRateScaling        ::Float64  ... scaling factor to modify the empirical rates.
    + resonanceEnergyMin:         ::Float64  ... minimum energy [Hartree] of the resonances to be considered.
    + resonanceEnergyMax:         ::Float64  ... maximum energy [Hartree] of the resonances to be considered.
"""
struct   EmpiricalTreatment
    doEmpiricalCorrections      ::Bool 
    doHydrogenicCorrections     ::Bool
    doMaximumlCorrection        ::Bool  
    doResonanceWindowCorrection ::Bool  
    nCore                       ::Int64   
    nFinal                      ::Int64   
    nHydrogenic                 ::Int64   
    nLowestCaptured             ::Int64   
    nLowerEmpirical             ::Int64   
    nUpperEmpirical             ::Int64   
    maximum_l                   ::Int64 
    hydrogenicEffectiveZ        ::Float64 
    hydrogenicRateScaling       ::Float64
    empiricalEffectiveZ         ::Float64
    empiricalRateScaling        ::Float64
    resonanceEnergyMin          ::Float64
    resonanceEnergyMax          ::Float64
end


# `Base.show(io::IO, tr::EmpiricalTreatment)`  ... prepares a proper printout of the tr::EmpiricalTreatment.
function Base.show(io::IO, tr::EmpiricalTreatment)
    println(io, "doEmpiricalCorrections:   $(tr.doEmpiricalCorrections)  ")
    println(io, "doHydrogenicCorrections:  $(tr.doHydrogenicCorrections)  ")
    println(io, "doMaximumlCorrection:     $(tr.doMaximumlCorrection)  ")
    println(io, "nCore:                    $(tr.nCore)  ")
    println(io, "nFinal:                   $(tr.nFinal)  ")
    println(io, "nHydrogenic:              $(tr.nHydrogenic)  ")
    println(io, "nLowestCaptured:          $(tr.nLowestCaptured)  ")
    println(io, "nUpperEmpirical:          $(tr.nUpperEmpirical)  ")
    println(io, "maximum_l:                $(tr.maximum_l)  ")
    println(io, "hydrogenicEffectiveZ:     $(tr.hydrogenicEffectiveZ)  ")
    println(io, "hydrogenicRateScaling:    $(tr.hydrogenicRateScaling)  ")
    println(io, "empiricalEffectiveZ:      $(tr.empiricalEffectiveZ)  ")
    println(io, "empiricalRateScaling:     $(tr.empiricalRateScaling)  ")
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
            capture-shells < asymptotic-shells. Various correction and multi-threading techiques can be applied
            to deal with or omit different classes of these shells.
    + calcRateAlpha         ::Bool                 
        ... True, if the DR rate coefficients are to be calculated, and false o/w.
    + calcHyperfineResolved ::Bool                 
        ... True, if the DR resonance strength are calculated for hyperfine-resolved levels, and false o/w.
            If true, it need to come together with calcOnlyPassages = true, and no fine-structure resolved rates 
            and strength are computed in this case.
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
            only include the contributions from those pathsways that are calculated here explicitly.
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
    Settings([E1], UseGauge[], false, false, false, false, PathwaySelection(), 0., 0., 0., Float64[],  
             DielectronicRecombination.AbstractCorrections[], CoulombInteraction())
end


"""
` (set::DielectronicRecombination.Settings;`

        multipoles=..,             gauges=..,                  
        calcOnlyPassages=..,       calcRateAlpha=..,         calcHyperfineResolved=..,         
        printBefore=..,            pathwaySelection=..,      electronEnergyShift=..,   photonEnergyShift=..,       
        mimimumPhotonEnergy=..,    temperatures=..,          corrections=..,           augerOperator=..)
                    
    ... constructor for modifying the given DielectronicRecombination.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::DielectronicRecombination.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,               gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    calcOnlyPassages::Union{Nothing,Bool}=nothing,                         calcRateAlpha::Union{Nothing,Bool}=nothing,  
    calcHyperfineResolved::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                              pathwaySelection::Union{Nothing,PathwaySelection}=nothing,
    electronEnergyShift::Union{Nothing,Float64}=nothing,                   photonEnergyShift::Union{Nothing,Float64}=nothing, 
    mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,                   temperatures::Union{Nothing,Array{Float64,1}}=nothing,    
    corrections::Union{Nothing,Array{AbstractCorrections,1}}=nothing,      augerOperator::Union{Nothing,AbstractEeInteraction}=nothing)
    
    if  multipoles           == nothing  multipolesx           = set.multipoles            else  multipolesx           = multipoles            end 
    if  gauges               == nothing  gaugesx               = set.gauges                else  gaugesx               = gauges                end 
    if  calcOnlyPassages     == nothing  calcOnlyPassagesx     = set.calcOnlyPassages      else  calcOnlyPassagesx     = calcOnlyPassages      end 
    if  calcRateAlpha        == nothing  calcRateAlphax        = set.calcRateAlpha         else  calcRateAlphax        = calcRateAlpha         end 
    if  calcHyperfineResolved== nothing  calcHyperfineResolvedx= set.calcHyperfineResolved else  calcHyperfineResolvedx= calcHyperfineResolved end 
    if  printBefore          == nothing  printBeforex          = set.printBefore           else  printBeforex          = printBefore           end 
    if  pathwaySelection     == nothing  pathwaySelectionx     = set.pathwaySelection      else  pathwaySelectionx     = pathwaySelection      end 
    if  electronEnergyShift  == nothing  electronEnergyShiftx  = set.electronEnergyShift   else  electronEnergyShiftx  = electronEnergyShift   end 
    if  photonEnergyShift    == nothing  photonEnergyShiftx    = set.photonEnergyShift     else  photonEnergyShiftx    = photonEnergyShift     end 
    if  mimimumPhotonEnergy  == nothing  mimimumPhotonEnergyx  = set.mimimumPhotonEnergy   else  mimimumPhotonEnergyx  = mimimumPhotonEnergy   end 
    if  temperatures         == nothing  temperaturesx         = set.temperatures          else  temperaturesx         = temperatures          end 
    if  corrections          == nothing  correctionsx          = set.corrections           else  correctionsx          = corrections           end 
    if  augerOperator        == nothing  augerOperatorx        = set.augerOperator         else  augerOperatorx        = augerOperator         end 

    Settings( multipolesx, gaugesx, calcOnlyPassagesx, calcRateAlphax, calcHyperfineResolvedx, printBeforex, 
              pathwaySelectionx, electronEnergyShiftx, photonEnergyShiftx, mimimumPhotonEnergyx, temperaturesx,
              correctionsx, augerOperatorx )
end


# `Base.show(io::IO, settings::DielectronicRecombination.Settings)`  ... prepares a proper printout of the variable settings::DielectronicRecombination.Settings.
function Base.show(io::IO, settings::DielectronicRecombination.Settings) 
    println(io, "multipoles:                 $(settings.multipoles)  ")
    println(io, "use-gauges:                 $(settings.gauges)  ")
    println(io, "calcOnlyPassages:           $(settings.calcOnlyPassages)  ")
    println(io, "calcRateAlpha:              $(settings.calcRateAlpha)  ")
    println(io, "calcHyperfineResolved:      $(settings.calcHyperfineResolved)  ")
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
`struct  DielectronicRecombination.Pathway`  
    ... defines a type for a dielectronic recombination pathways that may include the definition of channels and 
        their corresponding amplitudes.

    + initialLevel      ::Level                   ... initial-(state) level
    + intermediateLevel ::Level                   ... intermediate-(state) level
    + finalLevel        ::Level                   ... final-(state) level
    + electronEnergy    ::Float64                 ... energy of the (incoming, captured) electron
    + photonEnergy      ::Float64                 ... energy of the (emitted) photon
    + captureRate       ::Float64                 ... rate for the electron capture (Auger rate)
    + photonRate        ::EmProperty              ... rate for the photon emission
    + angularBeta       ::EmProperty              ... beta parameter of the photon emission
    + reducedStrength   ::EmProperty              ... reduced resonance strength S(i -> d -> f) * Gamma_d of this pathway;
                                                        this reduced strength does not require the knowledge of Gamma_d for each pathway.
    + captureChannels   ::Array{AutoIonization.Channel,1}   ... List of |i> -->  |n>   dielectronic (Auger) capture channels.
    + photonChannels    ::Array{PhotoEmission.Channel,1}    ... List of |n> -->  |f>   radiative stabilization channels.
"""
struct  Pathway
    initialLevel        ::Level
    intermediateLevel   ::Level
    finalLevel          ::Level
    electronEnergy      ::Float64
    photonEnergy        ::Float64 
    captureRate         ::Float64
    photonRate          ::EmProperty
    angularBeta         ::EmProperty
    reducedStrength     ::EmProperty
    captureChannels     ::Array{AutoIonization.Channel,1} 
    photonChannels      ::Array{PhotoEmission.Channel,1} 
end 


"""
`DielectronicRecombination.Pathway()`  
    ... constructor for an 'empty' instance of a dielectronic recombination pathway between a specified 
        initial, intermediate and final level.
"""
function Pathway()
    em = EmProperty(0., 0.)
    Pathway(initialLevel, intermediateLevel, finalLevel, 0., 0., 0., em, em, em, AutoIonization.Channel[], PhotoEmission.Channel[])
end


# `Base.show(io::IO, pathway::DielectronicRecombination.Pathway)`  ... prepares a proper printout of the variable pathway::DielectronicRecombination.Pathway.
function Base.show(io::IO, pathway::DielectronicRecombination.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "electronEnergy:             $(pathway.electronEnergy)  ")
    println(io, "photonEnergy:               $(pathway.photonEnergy)  ")
    println(io, "captureRate:                $(pathway.captureRate)  ")
    println(io, "photonRate:                 $(pathway.photonRate)  ")
    println(io, "angularBeta:                $(pathway.angularBeta)  ")
    println(io, "reducedStrength:            $(pathway.reducedStrength)  ")
    println(io, "captureChannels:            $(pathway.captureChannels)  ")
    println(io, "photonChannels:             $(pathway.photonChannels)  ")
end


"""
`struct  DielectronicRecombination.Passage`  
    ... defines a type for a dielectronic recombination passage, i.e. a (reduced) pathways, that include the 
        definition of channels and their corresponding amplitudes for the individual i --> m resonances, whereas
        the subsequent radiative stabilization is considered only later.

    + initialLevel      ::Level                   ... initial-(state) level
    + intermediateLevel ::Level                   ... intermediate-(state) level
    + electronEnergy    ::Float64                 ... energy of the (incoming, captured) electron
    + captureRate       ::Float64                 ... rate for the electron capture (Auger rate)
    + photonRate        ::EmProperty              ... rate for the photon emission
    + reducedStrength   ::EmProperty              
        ... reduced resonance strength Sum_f S(i -> d -> f) * Gamma_d of this passage; this reduced strength does 
            not require the knowledge of Gamma_d for the individual passage.
    + captureChannels   ::Array{AutoIonization.Channel,1}   ... List of |i> -->  |n>   dielectronic (Auger) capture channels.
"""
struct  Passage
    initialLevel        ::Level
    intermediateLevel   ::Level
    electronEnergy      ::Float64
    captureRate         ::Float64
    photonRate          ::EmProperty
    reducedStrength     ::EmProperty
    captureChannels     ::Array{AutoIonization.Channel,1} 
end 


"""
`DielectronicRecombination.Passage()`  
    ... constructor for an 'empty' instance of a dielectronic recombination passage between a specified 
        initial and intermediate level.
"""
function Passage()
    em = EmProperty(0., 0.)
    Passage(Level(), Level(), 0., 0., em, em, AutoIonization.Channel[])
end


# `Base.show(io::IO, pathway::DielectronicRecombination.Passage)`  
#   ... prepares a proper printout of the variable pathway::DielectronicRecombination.Passage.
function Base.show(io::IO, pathway::DielectronicRecombination.Passage) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "electronEnergy:             $(pathway.electronEnergy)  ")
    println(io, "captureRate:                $(pathway.captureRate)  ")
    println(io, "photonRate:                 $(pathway.photonRate)  ")
    println(io, "reducedStrength:            $(pathway.reducedStrength)  ")
    println(io, "captureChannels:            $(pathway.captureChannels)  ")
end


"""
`struct  DielectronicRecombination.Resonance`  
    ... defines a type for a dielectronic resonance as defined by a given initial and resonance level but by summing over all final levels

    + initialLevel      ::Level             ... initial-(state) level
    + intermediateLevel ::Level             ... intermediate-(state) level
    + resonanceEnergy   ::Float64           ... energy of the resonance w.r.t. the inital-state
    + resonanceStrength ::EmProperty        ... strength of this resonance due to the stabilization into any of the allowed final levels.
    + captureRate       ::Float64           ... capture (Auger) rate to form the intermediate resonance, starting from the initial level.
    + augerRate         ::Float64           ... total (Auger) rate for an electron emission of the intermediate resonance
    + photonRate        ::EmProperty        ... total photon rate for a photon emission, i.e. for stabilization.
"""
struct  Resonance
    initialLevel        ::Level
    intermediateLevel   ::Level
    resonanceEnergy     ::Float64 
    resonanceStrength   ::EmProperty
    captureRate         ::Float64
    augerRate           ::Float64
    photonRate          ::EmProperty
end 


"""
`DielectronicRecombination.Resonance()`  
    ... constructor for an 'empty' instance of a dielectronic resonance as defined by a given initial and resonance 
        level but by summing over all final levels.
"""
function Resonance()
    em = EmProperty(0., 0.)
    Resonance(initialLevel, intermediateLevel, 0., em, 0., 0., em)
end


# `Base.show(io::IO, resonance::DielectronicRecombination.Resonance)`  ... prepares a proper printout of the variable resonance::DielectronicRecombination.Resonance.
function Base.show(io::IO, resonance::DielectronicRecombination.Resonance) 
    println(io, "initialLevel:               $(resonance.initialLevel)  ")
    println(io, "intermediateLevel:          $(resonance.intermediateLevel)  ")
    println(io, "resonanceEnergy:            $(resonance.resonanceEnergy)  ")
    println(io, "resonanceStrength:          $(resonance.resonanceStrength)  ")
    println(io, "captureRate:                $(resonance.captureRate)  ")
    println(io, "augerRate:                  $(resonance.augerRate)  ")
    println(io, "photonRate:                 $(resonance.photonRate)  ")
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


# `Base.show(io::IO, resonance::DielectronicRecombination.ResonanceSelection)`  ... prepares a proper printout of resonance::DielectronicRecombination.ResonanceSelection.
function Base.show(io::IO, rSelection::DielectronicRecombination.ResonanceSelection) 
    println(io, "active:           $(rSelection.active)  ")
    println(io, "fromShells:       $(rSelection.fromShells)  ")
    println(io, "toShells:         $(rSelection.toShells)  ")
    println(io, "intoShells:       $(rSelection.intoShells)  ")
end


include("module-DielectronicRecombination-inc-FS-resolved.jl")
include("module-DielectronicRecombination-inc-HF-resolved.jl")

end # module
