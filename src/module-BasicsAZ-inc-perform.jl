
export  perform

# 7Apr25

"""
`BasicsAZ.resultKeyStrings(T::Type)`
    ... gives the string key(s) under which `Basics.perform` stores the result belonging to `T`, where `T` is
        either a settings type or one of the `ResultKeys`; an `Array{String,1}` is returned, EMPTY where the type
        is not one that produces a result.

        IT RETURNS A LIST, AND THERE ARE TWO REASONS.  The first is that JAC stores some quantities under MORE THAN
        ONE STRING: the Cascade modules write `"autoionization lines:"`, `"photoemission lines:"` and
        `"photoexcitation lines:"` where `Basics.perform` writes `"AutoIonization lines:"`, `"radiative lines:"`
        and `"photo-excitation lines:"` -- the same physical quantities under keys differing in case, in spacing
        or in the word chosen, and `"name"`/`"name:"` differ by a colon.  **One typed key lists every variant, so
        the caller never meets the second vocabulary at all.**  That is the whole of the maintainer's instruction
        of 05-Oct-2026 that there should not be two vocabularies, and it needed no renaming to carry out.

        The second reason: `DielectronicRecombination.Settings` produces
        either `"dielectronic recombination lines:"` or `"hyperfine-resolved dielectronic recombination lines:"`
        depending on a field of those very settings.  Asking by TYPE makes that ambiguity disappear for the
        caller -- they ask for the DR result and are handed whichever was computed -- where a string key forced
        them to know which branch had run.  It is the one place where the typed key is not merely safer but says
        something the string could not.

        THIS TABLE IS THE MIGRATION, and it is derived from `Basics.perform` itself: forty settings types, one key
        each but for the DR pair, and six `ResultKeys` for the results that have no settings type.  While the
        strings remain accepted, this is the only place the correspondence is written down.
"""
function  resultKeyStrings(T::Type)
    # --- the results that have no settings type of their own
    T === ResultKeys.Multiplet              &&  return( ["multiplet:"] )
    T === ResultKeys.Grid                   &&  return( ["grid:"] )
    T === ResultKeys.InitialMultiplet       &&  return( ["initialMultiplet"] )
    T === ResultKeys.FinalMultiplet         &&  return( ["finalMultiplet"] )
    T === ResultKeys.IntermediateMultiplet  &&  return( ["intermediateMultiplet"] )
    T === ResultKeys.IjfMultiplet           &&  return( ["IJF multiplet:"] )
    # --- level properties, several of which one computation may request
    T === Einstein.Settings                 &&  return( ["Einstein lines:"] )
    T === Hfs.Settings                      &&  return( ["HFS outcomes:"] )
    T === LandeZeeman.Settings              &&  return( ["Zeeman parameter outcomes:"] )
    T === StarkShift.Settings               &&  return( ["Stark-shift outcomes:"] )
    T === StarkZeeman.Settings              &&  return( ["Stark-Zeeman outcomes:"] )
    T === IsotopeShift.Settings             &&  return( ["Isotope parameter outcomes:"] )
    T === AlphaVariation.Settings           &&  return( ["alpha variation parameter outcomes:"] )
    T === FormFactor.Settings               &&  return( ["Form factor outcomes:"] )
    T === DecayYield.Settings               &&  return( ["Fluorescence and AutoIonization yield outcomes:"] )
    T === MultipolePolarizibility.Settings  &&  return( ["Polarizibility outcomes:"] )
    T === ReducedDensityMatrix.Settings     &&  return( ["RDM outcomes:"] )
    T === WeakInteractionEnhancement.Settings  &&  return( ["Weak-interaction enhancement outcomes:"] )
    # --- the one process of the computation
    T === AutoIonization.Settings           &&  return( ["AutoIonization lines:", "autoionization lines:"] )
    T === RayleighCompton.Settings          &&  return( ["Rayleigh-Compton lines:"] )
    T === ElectronCapture.Settings          &&  return( ["electron-capture lines:"] )
    T === DoubleAutoIonization.Settings     &&  return( ["Double-Auger lines:"] )
    T === DielectronicRecombination.Settings   &&
        return( ["dielectronic recombination lines:", "hyperfine-resolved dielectronic recombination lines:"] )
    T === MultiPhotonTransition.Settings    &&  return( ["multi-photon transition lines:"] )
    T === PhotoIonization.Settings          &&  return( ["photoionization lines:"] )
    T === PhotoDoubleIonization.Settings    &&  return( ["Single-photon double-ionization lines:"] )
    T === PhotoExcitation.Settings          &&  return( ["photo-excitation lines:", "photoexcitation lines:"] )
    T === PhotoExcitationAutoion.Settings   &&  return( ["photo-excitation-autoionization pathways:"] )
    T === PhotoExcitationFluores.Settings   &&  return( ["photo-excitation-fluorescence pathways:"] )
    T === PhotoEmission.Settings            &&  return( ["radiative lines:", "photoemission lines:"] )
    T === CoulombExcitation.Settings        &&  return( ["Coulomb excitation lines:"] )
    T === CoulombIonization.Settings        &&  return( ["Coulomb ionization lines:"] )
    T === RadiativeAuger.Settings           &&  return( ["radiative Auger sharings:"] )
    T === PhotoRecombination.Settings       &&  return( ["photo recombination lines:"] )
    T === ImpactExcitation.Settings         &&  return( ["impact-excitation lines:"] )
    T === InternalRecombination.Settings    &&  return( ["internal-recombination lines:"] )
    T === InternalConversion.Settings       &&  return( ["internal conversion lines:"] )
    T === TwoElectronOnePhoton.Settings     &&  return( ["two-electron-one-photon lines:"] )
    T === ParticleScattering.Settings       &&  return( ["particle-scattering events:"] )
    T === PhotonScattering.Settings         &&  return( ["photon-scattering lines:"] )
    T === BeamPhotoExcitation.Settings      &&  return( ["beam-assisted photo-excitation:"] )
    T === HyperfineInduced.Settings         &&  return( ["hyperfine-induced transitions:"] )
    T === MultiPhotonIonization.Settings    &&  return( ["multi-photon single ionization:"] )
    T === CrystalFieldEmission.Settings     &&  return( ["crystal-field-resolved emission lines:"] )
    T === PhotoRecombinationInterference.Settings  &&  return( ["photorecombination-interference pathways:"] )
    T === GeneralizedOscillatorStrength.Settings   &&  return( ["generalized oscillator strengths:"] )

    # --- the cascade and simulation results, which have a scheme rather than a settings type
    T === ResultKeys.CascadeData            &&  return( ["cascade data:"] )
    T === ResultKeys.CascadeScheme          &&  return( ["cascade scheme"] )
    T === ResultKeys.Name                   &&  return( ["name", "name:"] )
    T === ResultKeys.DataFormat             &&  return( ["data format:"] )
    T === ResultKeys.InitialMultiplets      &&  return( ["initial multiplets:"] )
    T === ResultKeys.GeneratedMultiplets    &&  return( ["generated multiplets:"] )
    T === ResultKeys.DielectronicMultiplets &&  return( ["dielectronic multiplets:"] )
    T === ResultKeys.ImpactExcitedMultiplets   &&  return( ["impact-excited multiplets:"] )
    T === ResultKeys.PhotoExcitedMultiplets &&  return( ["photoexcited multiplets:"] )
    T === ResultKeys.PhotoIonizedMultiplets &&  return( ["photoionized multiplets:"] )
    T === ResultKeys.DielectronicCaptureLines  &&  return( ["dielectronic-capture lines:"] )
    T === ResultKeys.HollowIonLineData      &&  return( ["hollow-ion line data:"] )
    T === ResultKeys.PhotoRecombinationLineData &&  return( ["photo-recombination line data:"] )
    T === ResultKeys.PhotoExcitationLineData   &&  return( ["photoexcitation line data:"] )
    T === ResultKeys.SimulationData         &&  return( ["data:"] )
    T === ResultKeys.SimulationProperty     &&  return( ["property:"] )
    # --- the generate(::Representation) results
    T === ResultKeys.MeanFieldBasis         &&  return( ["mean-field basis"] )
    T === ResultKeys.MeanFieldMultiplet     &&  return( ["mean-field multiplet"] )
    T === ResultKeys.MeanPotential          &&  return( ["mean potential"] )
    T === ResultKeys.Orbitals               &&  return( ["orbitals"] )
    T === ResultKeys.ReferenceMultiplet     &&  return( ["reference multiplet"] )
    T === ResultKeys.CiMultiplet            &&  return( ["CI multiplet"] )
    T === ResultKeys.GreenChannels          &&  return( ["Green channels"] )
    T === ResultKeys.EiiCrossSections       &&  return( ["EII cross sections:"] )

    return( String[] )
end


"""
`Basics.typedKeyFor(key::AbstractString)`
    ... names the typed key that replaces a deprecated string key; a `String` is returned, empty where the string
        is not one JAC ever used.  It is the migration aid for `apps/` scripts: the deprecation warning points at
        it, and it answers in one call what the table below would otherwise have to be read for.
"""
function  Basics.typedKeyFor(key::AbstractString)
    for  T  in  BasicsAZ.allResultKeyTypes()
        key in BasicsAZ.resultKeyStrings(T)   &&   return( string(T) )
    end
    return( "" )
end


"""
`BasicsAZ.allResultKeyTypes()`
    ... lists every type that labels a result -- the forty settings types and the `ResultKeys` -- so that a string
        can be mapped back to the typed key that replaces it; an `Array{Type,1}` is returned.
"""
function  allResultKeyTypes()
    return( Type[ ResultKeys.Multiplet, ResultKeys.Grid, ResultKeys.InitialMultiplet, ResultKeys.FinalMultiplet,
            ResultKeys.IntermediateMultiplet, ResultKeys.IjfMultiplet, ResultKeys.CascadeData,
            ResultKeys.CascadeScheme, ResultKeys.Name, ResultKeys.DataFormat, ResultKeys.InitialMultiplets,
            ResultKeys.GeneratedMultiplets, ResultKeys.DielectronicMultiplets, ResultKeys.ImpactExcitedMultiplets,
            ResultKeys.PhotoExcitedMultiplets, ResultKeys.PhotoIonizedMultiplets,
            ResultKeys.DielectronicCaptureLines, ResultKeys.HollowIonLineData,
            ResultKeys.PhotoRecombinationLineData, ResultKeys.PhotoExcitationLineData, ResultKeys.SimulationData,
            ResultKeys.SimulationProperty, ResultKeys.MeanFieldBasis, ResultKeys.MeanFieldMultiplet,
            ResultKeys.MeanPotential, ResultKeys.Orbitals, ResultKeys.ReferenceMultiplet,
            ResultKeys.CiMultiplet, ResultKeys.GreenChannels, ResultKeys.EiiCrossSections,
            Einstein.Settings, Hfs.Settings, LandeZeeman.Settings, StarkShift.Settings, StarkZeeman.Settings,
            IsotopeShift.Settings, AlphaVariation.Settings, FormFactor.Settings, DecayYield.Settings,
            MultipolePolarizibility.Settings, ReducedDensityMatrix.Settings, WeakInteractionEnhancement.Settings,
            AutoIonization.Settings, RayleighCompton.Settings, ElectronCapture.Settings,
            DoubleAutoIonization.Settings, DielectronicRecombination.Settings, MultiPhotonTransition.Settings,
            PhotoIonization.Settings, PhotoDoubleIonization.Settings, PhotoExcitation.Settings,
            PhotoExcitationAutoion.Settings, PhotoExcitationFluores.Settings, PhotoEmission.Settings,
            CoulombExcitation.Settings, CoulombIonization.Settings, RadiativeAuger.Settings,
            PhotoRecombination.Settings, ImpactExcitation.Settings, InternalRecombination.Settings,
            InternalConversion.Settings, TwoElectronOnePhoton.Settings, ParticleScattering.Settings,
            PhotonScattering.Settings, BeamPhotoExcitation.Settings, HyperfineInduced.Settings,
            MultiPhotonIonization.Settings, CrystalFieldEmission.Settings,
            PhotoRecombinationInterference.Settings, GeneralizedOscillatorStrength.Settings ] )
end


"""
`Base.haskey(r::Basics.PerformResults, T::Type)`
    ... true where this computation produced a result for the settings type or `ResultKeys` key `T`; a
        `value::Bool` is returned.  It completes the typed interface: a caller that must branch on whether a
        result exists -- `Cascade.Simulation` does, since a cascade may or may not carry photoexcitation lines --
        can now ask that question in the same vocabulary it uses to fetch the answer.
"""
function  Base.haskey(r::Basics.PerformResults, T::Type)
    for  k  in  BasicsAZ.resultKeyStrings(T)     haskey(r.dict, k)  &&  return( true )     end
    return( false )
end


"""
`Base.getindex(r::Basics.PerformResults, T::Type)`
    ... returns the result belonging to the settings type or `ResultKeys` key `T` -- `res[PhotoEmission.Settings]`,
        `res[Hfs.Settings]`, `res[ResultKeys.FinalMultiplet]`.  A wrong name is an `UndefVarError` where it is
        WRITTEN, naming the name, rather than a `KeyError` where the result is used.

        Where `T` is a type this computation did not produce a result for, the error says so and lists the settings
        types and keys that ARE present, translated back from the strings they are still stored under.
"""
function  Base.getindex(r::Basics.PerformResults, T::Type)
    cands = BasicsAZ.resultKeyStrings(T)
    if  isempty(cands)
        error("Basics.perform(): $T is not a type that labels a result.  Ask by the SETTINGS TYPE that produced " *
              "the result -- res[PhotoEmission.Settings], res[Hfs.Settings] -- or, for a multiplet or the grid, " *
              "by one of the ResultKeys:  Multiplet, Grid, InitialMultiplet, FinalMultiplet, " *
              "IntermediateMultiplet, IjfMultiplet.")
    end
    for  k  in  cands     haskey(r.dict, k)  &&  return( r.dict[k] )     end
    sa = "Basics.perform(): this computation produced no result for $T.\n"
    sa = sa * "   It produced " * string(length(r.dict)) * " result(s), under these keys:\n"
    for  k  in  sort(collect(keys(r.dict)))     sa = sa * "       \"" * k * "\"\n"     end
    error(sa)
end


"""
`Base.getindex(r::Basics.PerformResults, key::AbstractString)`
    ... returns the result stored under `key`; where the key is absent the error NAMES the keys that are present,
        and names the intended one where the request differs from it only in spacing, case or punctuation.
"""
function  Base.getindex(r::PerformResults, key::AbstractString)
    if  haskey(r.dict, key)
        # AND THERE IS A BOUNDARY THE TYPED KEYS DO NOT CROSS, found by the suite on 05-Oct-2026 and worth stating
        # where someone will meet it.  A typed key indexes what `perform` RETURNS.  It does NOT index cascade data
        # that has been SERIALISED and read back: `Cascade.Simulation` consumes
        # `simulation.computationData[i]["results"]`, a plain `Dict{String,Any}` restored from a JLD2 file, and
        # indexing that with a type raises `KeyError: key ResultKeys.InitialMultiplets not found`.  Making it work
        # would mean defining `getindex(::Dict{String,Any}, ::Type)`, which is type piracy on Base and would reach
        # every Dict in every package.  So the stored path keeps its strings deliberately, and the warnings it
        # raises here are correct rather than a conversion anybody forgot.
        #
        # STEP 4 OF THE MIGRATION, 06-Oct-2026: a string key that HAS a typed replacement now RAISES.
        #   It warned between 05 and 06-Oct, and the deprecation period ended when the evidence it was waiting on
        # arrived: two `/testExamples` sweeps, 254 branch runs, with ZERO KeyErrors and ZERO UndefVarErrors naming
        # a ResultKeys type (priority item 49).  The warning also earned its keep on the way out -- it found nine
        # sites in `examples/` that a grep had missed, because nobody knew that key family was unconverted.
        #   THE MESSAGE IS THE SAME ONE, raised rather than warned, so an `apps/` script that was never swept is
        # told exactly what to write at the moment it stops working.
        replacement = Basics.typedKeyFor(key)
        # A STRING WITH NO TYPED REPLACEMENT IS NOT DEPRECATED -- it is the only way to ask, and it keeps working
        # silently.  The RAS ladder's "step1", "step2", ... are the case: one per layer, their number unknown
        # until the ladder is built, so no fixed set of types can name them.  Raising on those would be telling
        # the user to do something that cannot be done.
        if  replacement != ""
            error("Basics.perform(): the string key \"" * key * "\" has been REPLACED BY A TYPED KEY and is no " *
                  "longer accepted.  Write  res[" * replacement * "]  instead.  A typed key is checked where it " *
                  "is WRITTEN, so a misspelling is an UndefVarError naming the name rather than a KeyError three " *
                  "call levels below it;  see priority item 39.")
        end
        return( r.dict[key] )
    end
    near = [ k  for k in keys(r.dict)  if  Basics.normalizeResultKey(k) == Basics.normalizeResultKey(key) ]
    sa   = "Basics.perform(): there is no result under the key \"" * key * "\".\n"
    if  !isempty(near)
        sa = sa * "   DID YOU MEAN  \"" * first(near) * "\" ?  It differs from what you asked for only in " *
                  "spacing, case or punctuation.\n"
    end
    sa = sa * "   This computation produced " * string(length(r.dict)) * " result(s), under these keys:\n"
    for  k  in  sort(collect(keys(r.dict)))     sa = sa * "       \"" * k * "\"\n"     end
    sa = sa * "   A key is a literal string and is not checked until it is used, which is why this message lists " *
              "them;  see priority item 39."
    error(sa)
end


"""
`Basics.perform(computation::Atomic.Computation)`  
    ... to perform the computation as prescribed by comp. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.perform(computation::Atomic.Computation; output::Bool=false)`  
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(computation::Atomic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = computation.nuclearModel

    Atomic.warnAboutUnusedSettings(computation)

    # Distinguish between the computation of level energies and properties and the simulation of atomic processes
    if   length(computation.configs) != 0
        multiplet = SelfConsistent.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings)
        LSjj.expandLevelsIntoLS(multiplet, computation.asfSettings.jjLS)
        #
        if output    results = Base.merge( results, Dict("multiplet:" => multiplet) ) 
                        results = Base.merge( results, Dict("grid:"      => computation.grid) )  end
        
        # Now compute all requested properties
        for settings  in computation.propertySettings
            if      typeof(settings) == Einstein.Settings
                outcome = Einstein.computeLines(multiplet,        computation.grid, settings)    
                if output    results = Base.merge( results, Dict("Einstein lines:" => outcome) )                  end
                #
            ## elseif  typeof(settings) == Hfs.Settings    && settings.calcIJFexpansion  
            ##     outcome = Hfs.computeHyperfineMultiplet(multiplet, nModel, computation.grid, settings)         
            ##     if output    results = Base.merge( results, Dict("IJF multiplet:" => outcome) )                   end
            ##     #
            elseif  typeof(settings) == Hfs.Settings
                outcome = Hfs.computeOutcomes(multiplet, nModel,  computation.grid, settings)         
                if output    results = Base.merge( results, Dict("HFS outcomes:" => outcome) )                    end
                #
            elseif  typeof(settings) == LandeZeeman.Settings 
                outcome = LandeZeeman.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Zeeman parameter outcomes:" => outcome) )       end
                #
            elseif  typeof(settings) == StarkShift.Settings
                outcome = StarkShift.computeOutcomes(multiplet, nModel,  computation.grid, settings)
                if output    results = Base.merge( results, Dict("Stark-shift outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == StarkZeeman.Settings
                outcome = StarkZeeman.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Stark-Zeeman outcomes:" => outcome) )           end
                #
            elseif  typeof(settings) == IsotopeShift.Settings 
                outcome = IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == AlphaVariation.Settings
                outcome = AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, computation.configs,
                                                          computation.asfSettings, settings)
                if output    results = Base.merge( results, Dict("alpha variation parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == FormFactor.Settings 
                outcome = FormFactor.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Form factor outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == DecayYield.Settings 
                outcome = DecayYield.computeOutcomes(computation.configs, computation.asfSettings, 
                                                        multiplet, nModel, computation.grid, settings)     
                if output    results = Base.merge( results, Dict("Fluorescence and AutoIonization yield outcomes:" => outcome) )   end
                #
            elseif  typeof(settings) == MultipolePolarizibility.Settings
                outcome = MultipolePolarizibility.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Polarizibility outcomes:" => outcome) )         end
                #
            elseif  typeof(settings) == ReducedDensityMatrix.Settings
                outcome = ReducedDensityMatrix.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("RDM outcomes:" => outcome) )         end
                #
            elseif  typeof(settings) == WeakInteractionEnhancement.Settings
                outcome = WeakInteractionEnhancement.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Weak-interaction enhancement outcomes:" => outcome) )   end
                #
            end
        end
        
    else
        initialMultiplet = SelfConsistent.performSCF(computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
        finalMultiplet   = SelfConsistent.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        if  output   results["initialMultiplet"] = initialMultiplet;   results["finalMultiplet"] = finalMultiplet    end 
        #
        if typeof(computation.processSettings) in [PhotoExcitationFluores.Settings, PhotoExcitationAutoion.Settings,
                                                   DielectronicRecombination.Settings, ResonantInelastic.Settings,
                                                   PhotoRecombinationInterference.Settings]
            intermediateMultiplet = SelfConsistent.performSCF(computation.intermediateConfigs, nModel, computation.grid, computation.intermediateAsfSettings)
            if  output   results["intermediateMultiplet"] = intermediateMultiplet    end 
        end
        #
        if      typeof(computation.processSettings) == AutoIonization.Settings 
            outcome = AutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == RayleighCompton.Settings 
            outcome = RayleighCompton.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Rayleigh-Compton lines:" => outcome) )                end
        elseif  typeof(computation.processSettings) == ElectronCapture.Settings
            ## The capture is the time reverse of an Auger transition, so `initialConfigs` are the ion BEFORE the
            ## capture and `finalConfigs` the doubly excited resonance it forms.
            outcome = ElectronCapture.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                   computation.processSettings)
            if output    results = Base.merge( results, Dict("electron-capture lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == DoubleAutoIonization.Settings   
            outcome = DoubleAutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Double-Auger lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == DielectronicRecombination.Settings
            ## settings.calcHyperfineResolved selects the hyperfine-resolved route, which builds the hyperfine
            ## multiplets itself from these same three ELECTRONIC multiplets and recouples the electronic
            ## amplitudes; it is a strict analogue of the fine-structure one and needs no separate entry.
            if  computation.processSettings.calcHyperfineResolved
                outcome = DielectronicRecombination.computeHfCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet,
                                                                          nModel, computation.grid, computation.processSettings)
                if output    results = Base.merge( results, Dict("hyperfine-resolved dielectronic recombination lines:" => outcome) )   end
            else
                outcome = DielectronicRecombination.computeCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                                                                        computation.grid, computation.processSettings)
                if output    results = Base.merge( results, Dict("dielectronic recombination lines:" => outcome) )   end
            end
        elseif  typeof(computation.processSettings) == MultiPhotonTransition.Settings
            outcome = MultiPhotonTransition.computeLines(computation.processSettings.scheme,
                                                            finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon transition lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == PhotoIonization.Settings   
            outcome = PhotoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photoionization lines:" => outcome) )                 end
        elseif  typeof(computation.processSettings) == PhotoDoubleIonization.Settings   
            outcome = PhotoDoubleIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Single-photon double-ionization lines:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoExcitation.Settings
            outcome = PhotoExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation lines:" => outcome) )                     end
        elseif  typeof(computation.processSettings) == PhotoExcitationAutoion.Settings  
            outcome = PhotoExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-autoionization pathways:" => outcome) )   end
        elseif  typeof(computation.processSettings) == PhotoExcitationFluores.Settings
            outcome = PhotoExcitationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-fluorescence pathways:" => outcome) )     end
        elseif  typeof(computation.processSettings) == PhotoEmission.Settings
            outcome = PhotoEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative lines:" => outcome) )                            end
        elseif  typeof(computation.processSettings) == ResonantInelastic.Settings 
            outcome = ResonantInelastic.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                        computation.grid, computation.processSettings) 
        elseif  typeof(computation.processSettings) == CoulombExcitation.Settings
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )                   end
        elseif  typeof(computation.processSettings) == RadiativeAuger.Settings
            outcome = RadiativeAuger.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative Auger sharings:" => outcome) )                   end
        elseif  typeof(computation.processSettings) == PhotoRecombination.Settings 
            outcome = PhotoRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo recombination lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == ImpactExcitation.Settings 
            outcome = ImpactExcitation.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == InternalRecombination.Settings 
            outcome = InternalRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal-recombination lines:" => outcome) )          end
        elseif  typeof(computation.processSettings) == TwoElectronOnePhoton.Settings
            ## The Green-function (intermediate) multiplet is NOT generated here; it is prepared by the user in a
            ## separate, standard Atomic.Computation -- whose configurations are chosen to include just those levels
            ## that are expected to contribute strongly to the second-order amplitude -- and handed over explicitly
            ## as settings.gMultiplet.
            outcome = TwoElectronOnePhoton.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("two-electron-one-photon lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == ParticleScattering.Settings 
            outcome = ParticleScattering.computeEvents(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("particle-scattering events:" => outcome) )         end
        elseif  typeof(computation.processSettings) == PhotonScattering.Settings
            outcome = PhotonScattering.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("photon-scattering lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == BeamPhotoExcitation.Settings 
            outcome = BeamPhotoExcitation.computeOutcomes(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("beam-assisted photo-excitation:" => outcome) )         end
        elseif  typeof(computation.processSettings) == HyperfineInduced.Settings 
            outcome = HyperfineInduced.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("hyperfine-induced transitions:" => outcome) )         end
            #
            #
        ## THE SEVEN TESTS BELOW USED TO COMPARE A TYPE WITH AN INSTANCE -- `typeof(x) == Coulex()` and its six
        ## siblings -- which is false for every x, so all seven processes fell through to error("stop b") and were
        ## unreachable.  They also tested the wrong THING: computation.processSettings is declared
        ## ::Basics.AbstractProcessSettings, so it is a Settings object and never one of the AbstractProcess
        ## singletons.  Each now tests the Settings type of the module whose function the branch calls, which is
        ## the form the rest of this chain already uses.  Fixed 24-Aug-2026.
        elseif  typeof(computation.processSettings) == CoulombExcitation.Settings
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == CoulombIonization.Settings
            ## Reserved and empty until 25-Aug-2026, when the module was written; this branch used to raise and say
            ## that CoulombIonization held no computation to reach.  It now does.
            outcome = CoulombIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                     computation.processSettings)
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )           end
        elseif  typeof(computation.processSettings) == MultiPhotonIonization.Settings
            outcome = MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                         computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == InternalConversion.Settings
            outcome = InternalConversion.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        elseif  typeof(computation.processSettings) == CrystalFieldEmission.Settings
            outcome = CrystalFieldEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("crystal-field-resolved emission lines:" => outcome) )    end
        elseif  typeof(computation.processSettings) == PhotoRecombinationInterference.Settings
            outcome = PhotoRecombinationInterference.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                                                                     computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("photorecombination-interference pathways:" => outcome) )    end
        elseif  typeof(computation.processSettings) == GeneralizedOscillatorStrength.Settings
            outcome = GeneralizedOscillatorStrength.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("generalized oscillator strengths:" => outcome) )    end
        else
            error("stop b")
        end
    end
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    # WRAPPED AT THE RETURN, not at creation: every `results = Base.merge(results, Dict(...))` above keeps working
    # on a plain Dict, so none of the forty-odd accumulation sites in this function had to be touched.
    return( isnothing(results) ? results : Basics.PerformResults(results) )
end
