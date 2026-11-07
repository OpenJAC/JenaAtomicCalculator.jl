
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
    T === HfsPaschenBack.Settings           &&  return( ["Hfs Paschen-Back outcomes:"] )
    T === IsotopeShift.Settings             &&  return( ["Isotope parameter outcomes:"] )
    T === AlphaVariation.Settings           &&  return( ["alpha variation parameter outcomes:"] )
    T === FormFactor.Settings               &&  return( ["Form factor outcomes:"] )
    T === DecayYield.Settings               &&  return( ["Fluorescence and AutoIonization yield outcomes:"] )
    T === MultipolePolarizibility.Settings  &&  return( ["Polarizibility outcomes:"] )
    T === DispersionCoefficient.Settings    &&  return( ["Dispersion coefficient outcomes:"] )
    T === ReducedDensityMatrix.Settings     &&  return( ["RDM outcomes:"] )
    T === WeakInteractionEnhancement.Settings  &&  return( ["Weak-interaction enhancement outcomes:"] )
    T === CrystalField.Settings             &&  return( ["Crystal-field splitting outcomes:"] )
    T === NuclearShielding.Settings         &&  return( ["Nuclear-shielding outcomes:"] )
    # --- the one process of the computation
    T === AutoIonization.Settings           &&  return( ["AutoIonization lines:", "autoionization lines:"] )
    T === RayleighCompton.Settings          &&  return( ["Rayleigh-Compton lines:"] )
    T === ElectronCapture.Settings          &&  return( ["electron-capture lines:"] )
    T === ThreeBodyRecombination.Settings   &&  return( ["three-body recombination lines:"] )
    T === DoubleAutoIonization.Settings     &&  return( ["Double-Auger lines:"] )
    T === DielectronicRecombination.Settings   &&
        return( ["dielectronic recombination lines:", "hyperfine-resolved dielectronic recombination lines:"] )
    T === MultiPhotonTransition.Settings    &&  return( ["multi-photon transition lines:"] )
    T === PhotoIonization.Settings          &&  return( ["photoionization lines:"] )
    T === PhotoDoubleIonization.Settings    &&  return( ["Single-photon double-ionization lines:"] )
    T === PhotoExcitation.Settings          &&  return( ["photo-excitation lines:", "photoexcitation lines:"] )
    T === StarkBroadening.Settings          &&  return( ["Stark-broadening lines:"] )
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
    T === PhotoIonizationInterference.Settings     &&  return( ["photoionization-interference pathways:"] )
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
    # --- Plasma.SuperConfigurationOpacityScheme's results
    T === ResultKeys.SuperConfigurationGroups               &&  return( ["superconfiguration groups"] )
    T === ResultKeys.SuperConfigurationPopulations          &&  return( ["superconfiguration populations"] )
    T === ResultKeys.SuperConfigurationTransitionArrays     &&  return( ["superconfiguration transition arrays"] )
    T === ResultKeys.SuperConfigurationTransitionGf         &&  return( ["superconfiguration transition gf"] )
    T === ResultKeys.SuperConfigurationOpacityContributions &&
        return( ["superconfiguration opacity contributions"] )
    T === ResultKeys.SuperConfigurationOpacitySpectrum      &&  return( ["superconfiguration opacity spectrum"] )
    T === ResultKeys.SuperConfigurationIonizationGroups     &&  return( ["superconfiguration ionization groups"] )
    T === ResultKeys.SuperConfigurationIonizationCharges    &&  return( ["superconfiguration ionization charges"] )
    T === ResultKeys.SuperConfigurationIonizationStagePopulations   &&
        return( ["superconfiguration ionization stage populations"] )
    T === ResultKeys.SuperConfigurationIonizationGroupPopulations   &&
        return( ["superconfiguration ionization group populations"] )
    T === ResultKeys.SuperConfigurationIonizationElectronDensity    &&
        return( ["superconfiguration ionization electron density"] )
    T === ResultKeys.SuperConfigurationMixtureSpeciesOpacities      &&
        return( ["superconfiguration mixture species opacities"] )
    T === ResultKeys.SuperConfigurationMixtureOpacitySpectrum       &&
        return( ["superconfiguration mixture opacity spectrum"] )
    T === ResultKeys.SuperConfigurationMixturePlanckOpacity         &&
        return( ["superconfiguration mixture planck opacity"] )
    T === ResultKeys.SuperConfigurationMixtureCharges               &&
        return( ["superconfiguration mixture charges"] )
    T === ResultKeys.SuperConfigurationMixtureStagePopulations      &&
        return( ["superconfiguration mixture stage populations"] )
    T === ResultKeys.SuperConfigurationMixtureElectronDensity       &&
        return( ["superconfiguration mixture electron density"] )
    T === ResultKeys.SuperConfigurationMixtureBoundFreeOpacity      &&
        return( ["superconfiguration mixture bound-free opacity"] )
    T === ResultKeys.SuperConfigurationMixtureFreeFreeOpacity       &&
        return( ["superconfiguration mixture free-free opacity"] )
    T === ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum  &&
        return( ["superconfiguration mixture total opacity spectrum"] )
    T === ResultKeys.SuperConfigurationMixtureTotalPlanckOpacity    &&
        return( ["superconfiguration mixture total planck opacity"] )
    T === ResultKeys.SuperConfigurationMixtureRosselandOpacity      &&
        return( ["superconfiguration mixture rosseland opacity"] )
    T === ResultKeys.SuperConfigurationMixtureResonanceCorrections  &&
        return( ["superconfiguration mixture resonance corrections"] )

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
            ResultKeys.SuperConfigurationGroups, ResultKeys.SuperConfigurationPopulations,
            ResultKeys.SuperConfigurationTransitionArrays, ResultKeys.SuperConfigurationTransitionGf,
            ResultKeys.SuperConfigurationOpacityContributions, ResultKeys.SuperConfigurationOpacitySpectrum,
            ResultKeys.SuperConfigurationIonizationGroups, ResultKeys.SuperConfigurationIonizationCharges,
            ResultKeys.SuperConfigurationIonizationStagePopulations,
            ResultKeys.SuperConfigurationIonizationGroupPopulations,
            ResultKeys.SuperConfigurationIonizationElectronDensity,
            ResultKeys.SuperConfigurationMixtureSpeciesOpacities, ResultKeys.SuperConfigurationMixtureOpacitySpectrum,
            ResultKeys.SuperConfigurationMixturePlanckOpacity, ResultKeys.SuperConfigurationMixtureCharges,
            ResultKeys.SuperConfigurationMixtureStagePopulations, ResultKeys.SuperConfigurationMixtureElectronDensity,
            ResultKeys.SuperConfigurationMixtureBoundFreeOpacity, ResultKeys.SuperConfigurationMixtureFreeFreeOpacity,
            ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum,
            ResultKeys.SuperConfigurationMixtureTotalPlanckOpacity, ResultKeys.SuperConfigurationMixtureRosselandOpacity,
            ResultKeys.SuperConfigurationMixtureResonanceCorrections,
            Einstein.Settings, Hfs.Settings, LandeZeeman.Settings, StarkShift.Settings, StarkZeeman.Settings,
            HfsPaschenBack.Settings,
            IsotopeShift.Settings, AlphaVariation.Settings, FormFactor.Settings, DecayYield.Settings,
            MultipolePolarizibility.Settings, DispersionCoefficient.Settings, ReducedDensityMatrix.Settings,
            WeakInteractionEnhancement.Settings, CrystalField.Settings, NuclearShielding.Settings,
            AutoIonization.Settings, RayleighCompton.Settings, ElectronCapture.Settings,
            DoubleAutoIonization.Settings, DielectronicRecombination.Settings, ThreeBodyRecombination.Settings,
            MultiPhotonTransition.Settings,
            PhotoIonization.Settings, PhotoDoubleIonization.Settings, PhotoExcitation.Settings,
            StarkBroadening.Settings,
            PhotoExcitationAutoion.Settings, PhotoExcitationFluores.Settings, PhotoEmission.Settings,
            CoulombExcitation.Settings, CoulombIonization.Settings, RadiativeAuger.Settings,
            PhotoRecombination.Settings, ImpactExcitation.Settings, InternalRecombination.Settings,
            InternalConversion.Settings, TwoElectronOnePhoton.Settings, ParticleScattering.Settings,
            PhotonScattering.Settings, BeamPhotoExcitation.Settings, HyperfineInduced.Settings,
            MultiPhotonIonization.Settings, CrystalFieldEmission.Settings,
            PhotoRecombinationInterference.Settings, PhotoIonizationInterference.Settings,
            GeneralizedOscillatorStrength.Settings ] )
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
`BasicsAZ.checkFrozenOrbitalsSpan(frozenOrbitals::Dict{Subshell,Radial.Orbital}, configs::Array{Configuration,1},
                                    configName::String, sourceName::String)`
    ... tests that `frozenOrbitals` holds every subshell the given configurations occupy and RAISES an error naming the absent ones
        if it does not;  `nothing` is returned.  It exists because the failure it replaces is otherwise a bare
        `KeyError: key 2s_1/2 not found` thrown by `Dict` inside the CI matrix setup and wrapped in a `TaskFailedException` --
        loud, but telling a user nothing about what to change.  MEASURED 02-Nov-2026: an SCF on `[1s^2]` yields the single
        orbital `1s_1/2`, so a Green space of `[1s, 2s, 2p]` frozen onto it lacks `2s_1/2, 2p_1/2, 2p_3/2`, while the same SCF
        over the union `[1s^2, 1s 2s, 1s 2p]` lacks none.  That is the whole content of the check:  ONE orbital set is usable
        only where it SPANS every subshell the multiplets built on it occupy.
"""
function checkFrozenOrbitalsSpan(frozenOrbitals::Dict{Subshell,Radial.Orbital}, configs::Array{Configuration,1},
                                    configName::String, sourceName::String)
    relconfs = ConfigurationR[]
    for  conf in configs
        append!( relconfs, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf) )
    end
    absent = setdiff( Basics.generateSubshellList(relconfs), collect(keys(frozenOrbitals)) )
    if  !isempty(absent)
        error("BasicsAZ.checkFrozenOrbitalsSpan(): the orbital set taken from $sourceName does not span $configName -- "    *
              "it has no $(join(string.(absent), ", ")).  ONE shared orbital set is usable only where it SPANS every "       *
              "subshell the other multiplets occupy, so widen the configuration list that PRODUCED $sourceName until it "   *
              "generates those subshells (for a 1s^2 initial state and a 2s/2p Green space that is "                        *
              "[1s^2, 1s 2s, 1s 2p] rather than [1s^2]), or use Basics.SeparateOrbitals() and accept that a "               *
              "second-order amplitude through these multiplets carries the dropped -E <r|s> term.  "                        *
              "See Basics.AbstractOrbitalSource.")
    end

    return( nothing )
end


"""
`Basics.checkGridConvergence(computation::Atomic.Computation, extract::Function; factor::Float64=2.0, 
                                tolerance::Float64=1.0e-2, printout::Bool=true)`  
    ... re-performs `computation` on a twin grid whose outer step `hp` is finer by `factor`, with the radial BOX held, and
        compares whatever number(s) `extract` pulls out of the two result dictionaries.  It is the check that both measured
        tables in `Basics.recommendedGrid` end on -- halve `hp` until the answer stops moving -- and it is here because
        nothing else in JAC performs it:  every other grid and continuum guard is A PRIORI, read off the grid before the
        answer exists, while the failure this one catches is POSTERIOR.  `extract` is handed the `Basics.PerformResults` of
        `Basics.perform(...; output=true)` -- an `AbstractDict{String,Any}` read by TYPED keys -- and may return a number or a
        collection of numbers.  A value::NamedTuple with the fields `converged`, `relative`, `index`, `coarse`, `fine`,
        `coarseGrid` and `fineGrid` is returned.

            Basics.checkGridConvergence(comp, r -> [lev.energy for lev in r[ResultKeys.Multiplet].levels])
            Basics.checkGridConvergence(comp, r -> [ln.crossSection.Coulomb for ln in r[PhotoIonization.Settings]])

        WHY THE MEAN FIELD IS RE-DERIVED RATHER THAN THE ORBITALS CARRIED OVER.  A halved `hp` shares NO mesh with the
        original -- measured 06-Nov-2026, 2 of 1330 points survive and both lie at r ~ 1e-8 -- because the harmonic blend of
        `h*r` and `hp` carries the step change inward to every radius, not only to the outer linear part.  So nothing can be
        restricted onto the finer mesh, and JAC has no interpolation either (`Basics.interpolateOnGridGrasp92` is a stub that
        raises).  Re-performing the whole computation is therefore the only honest route, and it is what both of those tables
        actually did.  The cost is that the computation runs TWICE.

        THE DEFAULT `tolerance` IS MEASURED, and the calibration is the second table in that same docstring:  RR into 1s 4s of
        H-like Bi with the box fixed at 142 a.u. gives sigma_C = 3.8338e-10, 5.0133e-09, 9.7472e-03, 1.1324e-02, 1.1324e-02 b
        at hp = 0.472, 0.200, 0.100, 0.050, 0.025.  Read as this check reads it -- each grid against its own halving -- those
        are relative changes of 1.00, 1.00, 0.139 and 0.000, so EVERY threshold in (0, 0.139) separates the three unconverged
        grids from the converged one and none outside it does.  1.0e-2 sits inside that interval.  Note what the same table
        says about the alternatives:  the two gauges agree to better than 2 % at every row, so a gauge comparison is BLIND
        here while the answer is wrong by seven orders.

        AND THE BOX IS HELD ONLY TO WITHIN ONE OUTER STEP, since `Radial.Grid` lands on the first mesh point beyond the
        requested `rbox`:  20.1267 against 20.0046 a.u. for hp = 2e-2 against 1e-2, i.e. 0.6 %.  That is harmless under
        `CoulombFunctionNorm`, measured insensitive to the box, and is NOT harmless under a sine convention, where the number
        of oscillations in the box is what the normalization fits -- so under those conventions a failed verdict may mean
        "the box moved" rather than "the answer has not converged".

    + computation   ::Atomic.Computation   ... the computation to be repeated; its own grid is the coarse one of the pair.
    + extract       ::Function             ... pulls the number(s) to be compared out of the result dictionary of `Basics.perform`.
    + factor        ::Float64              ... factor by which the outer step `hp` is refined; it must exceed 1.
    + tolerance     ::Float64              ... largest relative change that still counts as converged.
    + printout      ::Bool                 ... whether the two grids and the comparison are printed.
"""
function Basics.checkGridConvergence(computation::Atomic.Computation, extract::Function; factor::Float64=2.0, 
                                        tolerance::Float64=1.0e-2, printout::Bool=true)
    if  factor <= 1.0      error("Basics.checkGridConvergence(): factor = $factor must exceed 1;  it is the factor by "  *
                                 "which the outer step hp is REFINED.")                                                   end
    if  tolerance <= 0.0   error("Basics.checkGridConvergence(): tolerance = $tolerance must be positive.")               end

    coarseGrid = computation.grid
    rbox       = coarseGrid.r[coarseGrid.NoPoints]
    fineGrid   = Radial.Grid(coarseGrid; hp = coarseGrid.hp/factor, rbox = rbox)

    if  printout
        println("\n>> Basics.checkGridConvergence():  hp = "       * @sprintf("%.4e", coarseGrid.hp)                      *
                " -> "                                             * @sprintf("%.4e", fineGrid.hp)                        *
                " a.u., NoPoints = "                               * string(coarseGrid.NoPoints)                          *
                " -> "                                             * string(fineGrid.NoPoints)                            *
                ", nsL = "                                         * string(coarseGrid.nsL)                                *
                " -> "                                             * string(fineGrid.nsL)                                 *
                ", box "                                           * @sprintf("%.4f", rbox)                               *
                " -> "                                             * @sprintf("%.4f", fineGrid.r[fineGrid.NoPoints])      *
                " a.u.\n>> The computation now runs TWICE.")
    end

    coarseRes = Basics.perform(computation; output=true)
    fineRes   = Basics.perform( Atomic.Computation(computation; grid=fineGrid); output=true )

    wa      = extract(coarseRes);                                 wb    = extract(fineRes)
    vCoarse = wa isa Number ? Float64[wa] : collect(Float64, wa); vFine = wb isa Number ? Float64[wb] : collect(Float64, wb)
    if  length(vCoarse) == 0
        error("Basics.checkGridConvergence(): `extract` returned nothing on the given grid, so there is nothing to "      *
              "compare;  check the key it reads against the result dictionary of Basics.perform.")
    end
    if  length(vCoarse) != length(vFine)
        error("Basics.checkGridConvergence(): `extract` returned $(length(vCoarse)) number(s) on the given grid and "     *
              "$(length(vFine)) on the refined one, so the two cannot be compared term by term.  A finer grid may "       *
              "legitimately yield a different NUMBER of lines or levels, in which case `extract` must select the ones "   *
              "to be compared rather than return all of them.")
    end

    relative = zeros( length(vCoarse) )
    for  k = 1:length(vCoarse)
        den         = max( abs(vFine[k]), abs(vCoarse[k]) )
        relative[k] = den > 0.0  ?  abs(vFine[k] - vCoarse[k]) / den  :  0.0
    end
    worst = maximum(relative);    index = argmax(relative);    converged = (worst <= tolerance)

    if  printout
        for  k = 1:length(vCoarse)
            if  length(vCoarse) <= 10  ||  k == index
                println("   [" * string(k) * "]  " * @sprintf("%16.8e", vCoarse[k]) * "  ->  "                            *
                        @sprintf("%16.8e", vFine[k]) * "    relative change " * @sprintf("%.3e", relative[k])             *
                        (k == index  ?  "   <-- worst"  :  ""))
            end
        end
        println(">> Worst relative change " * @sprintf("%.3e", worst) * " against a tolerance of "                        *
                @sprintf("%.3e", tolerance) * ":  "                                                                       *
                (converged  ?  "CONVERGED on the given grid."                                                             :
                               "NOT CONVERGED -- the answer is still moving with hp, so the given grid is too coarse "    *
                               "for the quantity compared.  Refine it and repeat."))
    end
    if  !converged
        Defaults.warn(AddWarning(), "Basics.checkGridConvergence(): the quantity compared moved by "                      *
                      @sprintf("%.3e", worst) * " when hp was refined by " * @sprintf("%.2f", factor)                     *
                      ", above the tolerance of " * @sprintf("%.3e", tolerance) * ";  the given grid is too coarse "      *
                      "for it.")
    end

    return( (converged=converged, relative=worst, index=index, coarse=vCoarse, fine=vFine,
             coarseGrid=coarseGrid, fineGrid=fineGrid) )
end


"""
`Basics.perform(computation::Atomic.Computation; output::Bool=false)`
    ... to perform the computation as prescribed by computation, i.e. the SCF and CI steps followed by whichever properties and
        processes its settings request; all relevant intermediate and final results are printed to screen (stdout). With
        output = false, nothing is returned. With output = true, a results::Basics.PerformResults is returned, an
        AbstractDict{String,Any} whose keys depend on the type and the specifications of the computation.
"""
function Basics.perform(computation::Atomic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = computation.nuclearModel

    Atomic.warnAboutUnusedSettings(computation)

    # Distinguish between the computation of level energies and properties and the simulation of atomic processes
    if   length(computation.configs) != 0
        multiplet = SelfConsistent.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings)
        LSjj.expandLevelsIntoLS(multiplet, computation.asfSettings.jjLS)

        if output    results = Base.merge( results, Dict("multiplet:" => multiplet) ) 
                        results = Base.merge( results, Dict("grid:"      => computation.grid) )  end
        
        # Now compute all requested properties
        for settings  in computation.propertySettings
            if      typeof(settings) == Einstein.Settings
                outcome = Einstein.computeLines(multiplet,        computation.grid, settings)    
                if output    results = Base.merge( results, Dict("Einstein lines:" => outcome) )                  end

            # elseif  typeof(settings) == Hfs.Settings    && settings.calcIJFexpansion  
            #     outcome = Hfs.computeHyperfineMultiplet(multiplet, nModel, computation.grid, settings)         
            #     if output    results = Base.merge( results, Dict("IJF multiplet:" => outcome) )                   end
            #     #
            elseif  typeof(settings) == Hfs.Settings
                outcome = Hfs.computeOutcomes(multiplet, nModel,  computation.grid, settings)         
                if output    results = Base.merge( results, Dict("HFS outcomes:" => outcome) )                    end

            elseif  typeof(settings) == LandeZeeman.Settings 
                outcome = LandeZeeman.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Zeeman parameter outcomes:" => outcome) )       end

            elseif  typeof(settings) == StarkShift.Settings
                outcome = StarkShift.computeOutcomes(multiplet, nModel,  computation.grid, settings)
                if output    results = Base.merge( results, Dict("Stark-shift outcomes:" => outcome) )            end

            elseif  typeof(settings) == StarkZeeman.Settings
                outcome = StarkZeeman.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Stark-Zeeman outcomes:" => outcome) )           end

            elseif  typeof(settings) == HfsPaschenBack.Settings
                outcome = HfsPaschenBack.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Hfs Paschen-Back outcomes:" => outcome) )       end

            elseif  typeof(settings) == IsotopeShift.Settings 
                outcome = IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end

            elseif  typeof(settings) == AlphaVariation.Settings
                outcome = AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, computation.configs,
                                                          computation.asfSettings, settings)
                if output    results = Base.merge( results, Dict("alpha variation parameter outcomes:" => outcome) )      end

            elseif  typeof(settings) == FormFactor.Settings 
                outcome = FormFactor.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Form factor outcomes:" => outcome) )            end

            elseif  typeof(settings) == DecayYield.Settings 
                outcome = DecayYield.computeOutcomes(computation.configs, computation.asfSettings, 
                                                        multiplet, nModel, computation.grid, settings)     
                if output    results = Base.merge( results, Dict("Fluorescence and AutoIonization yield outcomes:" => outcome) )   end

            elseif  typeof(settings) == MultipolePolarizibility.Settings
                outcome = MultipolePolarizibility.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Polarizibility outcomes:" => outcome) )         end

            elseif  typeof(settings) == DispersionCoefficient.Settings
                outcome = DispersionCoefficient.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Dispersion coefficient outcomes:" => outcome) )  end

            elseif  typeof(settings) == ReducedDensityMatrix.Settings
                outcome = ReducedDensityMatrix.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("RDM outcomes:" => outcome) )         end

            elseif  typeof(settings) == WeakInteractionEnhancement.Settings
                outcome = WeakInteractionEnhancement.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Weak-interaction enhancement outcomes:" => outcome) )   end

            elseif  typeof(settings) == CrystalField.Settings
                outcome = CrystalField.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Crystal-field splitting outcomes:" => outcome) )        end

            elseif  typeof(settings) == NuclearShielding.Settings
                outcome = NuclearShielding.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Nuclear-shielding outcomes:" => outcome) )              end

            end
        end
        
    else
        # WHICH ORBITAL SET THE TWO MULTIPLETS ARE BUILT ON -- `computation.orbitalSource`, see Basics.AbstractOrbitalSource.
        # SeparateOrbitals() is the default and its branch below is character for character what this code always did,
        # two independent SCF runs;  the other two run ONE SCF and build the remaining multiplet as a CI on those frozen
        # orbitals.  That is what an explicit SECOND-ORDER amplitude needs: ManyElectron.matrixElement_HminusE drops the
        # overlap term -E <r|s>, which vanishes identically only on one orthonormal set.  The statement ORDER inside the
        # default branch is preserved deliberately -- performSCF resets the global relativistic subshell list, so moving
        # an expandLevelsIntoLS across an SCF call is not obviously inert and is not worth risking for tidiness.
        if      computation.orbitalSource isa Basics.SeparateOrbitals
            initialMultiplet = SelfConsistent.performSCF(computation.initialConfigs, nModel, computation.grid,
                                                         computation.initialAsfSettings)
            LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
            finalMultiplet   = SelfConsistent.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
            LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        elseif  computation.orbitalSource isa Basics.InitialStateOrbitals
            initialMultiplet = SelfConsistent.performSCF(computation.initialConfigs, nModel, computation.grid,
                                                         computation.initialAsfSettings)
            LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
            sharedOrbitals   = initialMultiplet.levels[1].basis.orbitals
            BasicsAZ.checkFrozenOrbitalsSpan(sharedOrbitals, computation.finalConfigs, "finalConfigs", "initialConfigs")
            finalMultiplet   = Hamiltonian.performCIwithFrozenOrbitals(computation.finalConfigs, sharedOrbitals, nModel,
                                                                       computation.grid, computation.finalAsfSettings)
            LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        elseif  computation.orbitalSource isa Basics.FinalStateOrbitals
            finalMultiplet   = SelfConsistent.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
            LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
            sharedOrbitals   = finalMultiplet.levels[1].basis.orbitals
            BasicsAZ.checkFrozenOrbitalsSpan(sharedOrbitals, computation.initialConfigs, "initialConfigs", "finalConfigs")
            initialMultiplet = Hamiltonian.performCIwithFrozenOrbitals(computation.initialConfigs, sharedOrbitals, nModel,
                                                                       computation.grid, computation.initialAsfSettings)
            LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
        elseif  computation.orbitalSource isa Basics.GivenOrbitals
            # THE ONE THAT COMPLETES A SECOND-ORDER COMPUTATION.  The Green space is built by the CALLER before `perform`
            # runs -- it sits in the process settings -- so the caller already holds the orbitals, and the SCF behind them
            # may have run over a WIDER configuration list than any multiplet's CI.  That is what a 1s^2 initial state
            # needs when its Green space occupies 2s and 2p:  measured 02-Nov-2026 on He-like Be, one SCF over
            # [1s^2, 1s 2s, 1s 2p] carries a 1s^2 initial state and a [1s, 2s, 2p] Green space with condition (b) met, at
            # a price of 3.3e-02 Ha (0.24 %) on the initial-state energy against its own dedicated SCF.
            sharedOrbitals   = computation.orbitalSource.orbitals
            BasicsAZ.checkFrozenOrbitalsSpan(sharedOrbitals, computation.initialConfigs, "initialConfigs", "orbitalSource")
            BasicsAZ.checkFrozenOrbitalsSpan(sharedOrbitals, computation.finalConfigs,   "finalConfigs",   "orbitalSource")
            initialMultiplet = Hamiltonian.performCIwithFrozenOrbitals(computation.initialConfigs, sharedOrbitals, nModel,
                                                                       computation.grid, computation.initialAsfSettings)
            LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
            finalMultiplet   = Hamiltonian.performCIwithFrozenOrbitals(computation.finalConfigs, sharedOrbitals, nModel,
                                                                       computation.grid, computation.finalAsfSettings)
            LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        else
            error("Basics.perform(): computation.orbitalSource is $(computation.orbitalSource), which this method does not "  *
                  "know how to honour.  The supported values are Basics.SeparateOrbitals(), Basics.InitialStateOrbitals(), "  *
                  "Basics.FinalStateOrbitals() and Basics.GivenOrbitals(orbitals);  see Basics.AbstractOrbitalSource.")
        end
        if  output   results["initialMultiplet"] = initialMultiplet;   results["finalMultiplet"] = finalMultiplet    end 
        # HOW A SECOND-ORDER CALLER REACHES THE SHARED SET, since its Green space must sit on these same orbitals and
        # `settings.gMultiplet` is built by the caller before `perform` ever runs:  it is already returned, as
        # `res[ResultKeys.InitialMultiplet].levels[1].basis.orbitals` for InitialStateOrbitals() and the FinalMultiplet
        # key for FinalStateOrbitals().  Hand that dictionary to Hamiltonian.performCIwithFrozenOrbitals to build the
        # Green space on it.  NO separate result entry is added on purpose:  it would need its own ResultKeys type to be
        # readable at all, and that is a new promise in the public surface for something already reachable.

        if typeof(computation.processSettings) in [PhotoExcitationFluores.Settings, PhotoExcitationAutoion.Settings,
                                                   DielectronicRecombination.Settings, ResonantInelastic.Settings,
                                                   PhotoRecombinationInterference.Settings, PhotoIonizationInterference.Settings]
            intermediateMultiplet = SelfConsistent.performSCF(computation.intermediateConfigs, nModel, computation.grid,
                                                              computation.intermediateAsfSettings)
            if  output   results["intermediateMultiplet"] = intermediateMultiplet    end 
        end

        if      typeof(computation.processSettings) == AutoIonization.Settings
            # continuumExchange is resolved HERE, where computation.finalAsfSettings.scField is in scope, and travels on as an
            # ordinary (now concrete) field of the settings object handed to computeLines -- see Continuum.resolveExchangeDefault,
            # item 99.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = AutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           AutoIonization.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("AutoIonization lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == RayleighCompton.Settings 
            outcome = RayleighCompton.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Rayleigh-Compton lines:" => outcome) )                end
        elseif  typeof(computation.processSettings) == ElectronCapture.Settings
            # The capture is the time reverse of an Auger transition, so `initialConfigs` are the ion BEFORE the
            # capture and `finalConfigs` the doubly excited resonance it forms. The free electron approaches that
            # INITIAL ion (see ElectronCapture.computeAmplitudesProperties), so it is initialAsfSettings.scField
            # that decides the automatic default here, unlike every OUTGOING-electron process above.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.initialAsfSettings.scField)
            outcome = ElectronCapture.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           ElectronCapture.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("electron-capture lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == DoubleAutoIonization.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = DoubleAutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           DoubleAutoIonization.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("Double-Auger lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == ThreeBodyRecombination.Settings
            # initialConfigs are the ion BEFORE capture (one electron fewer); finalConfigs the bound level the
            # two free electrons recombine into -- the same direction as ElectronCapture, above.
            outcome = ThreeBodyRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                           computation.processSettings)
            if output    results = Base.merge( results, Dict("three-body recombination lines:" => outcome) )        end
        elseif  typeof(computation.processSettings) == DielectronicRecombination.Settings
            # settings.calcHyperfineResolved selects the hyperfine-resolved route, which builds the hyperfine
            # multiplets itself from these same three ELECTRONIC multiplets and recouples the electronic
            # amplitudes; it is a strict analogue of the fine-structure one and needs no separate entry.
            # The captured electron approaches the INITIAL ion in both routes (as in ElectronCapture and
            # PhotoRecombination), so initialAsfSettings.scField decides the automatic default.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.initialAsfSettings.scField)
            dcrSettings = DielectronicRecombination.Settings(computation.processSettings; continuumExchange=exchangex)
            if  computation.processSettings.calcHyperfineResolved
                outcome = DielectronicRecombination.computeHfCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet,
                                                                          nModel, computation.grid, dcrSettings)
                if output    results = Base.merge( results, Dict("hyperfine-resolved dielectronic recombination lines:" => outcome) )   end
            else
                outcome = DielectronicRecombination.computeCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                                                                        computation.grid, dcrSettings)
                if output    results = Base.merge( results, Dict("dielectronic recombination lines:" => outcome) )   end
            end
        elseif  typeof(computation.processSettings) == MultiPhotonTransition.Settings
            outcome = MultiPhotonTransition.computeLines(computation.processSettings.scheme,
                                                            finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon transition lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == PhotoIonization.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = PhotoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           PhotoIonization.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("photoionization lines:" => outcome) )                 end
        elseif  typeof(computation.processSettings) == PhotoDoubleIonization.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = PhotoDoubleIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           PhotoDoubleIonization.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("Single-photon double-ionization lines:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoExcitation.Settings
            outcome = PhotoExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("photo-excitation lines:" => outcome) )                     end
        elseif  typeof(computation.processSettings) == StarkBroadening.Settings
            outcome = StarkBroadening.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("Stark-broadening lines:" => outcome) )                     end
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
            # The captured electron approaches the INITIAL ion (see PhotoRecombination.computeAmplitudesProperties), so
            # initialAsfSettings.scField decides the automatic default here, not finalAsfSettings.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.initialAsfSettings.scField)
            outcome = PhotoRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           PhotoRecombination.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("photo recombination lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == ImpactExcitation.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = ImpactExcitation.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           ImpactExcitation.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("impact-excitation lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == InternalRecombination.Settings 
            outcome = InternalRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal-recombination lines:" => outcome) )          end
        elseif  typeof(computation.processSettings) == TwoElectronOnePhoton.Settings
            # The Green-function (intermediate) multiplet is NOT generated here; it is prepared by the user in a
            # separate, standard Atomic.Computation -- whose configurations are chosen to include just those levels
            # that are expected to contribute strongly to the second-order amplitude -- and handed over explicitly
            # as settings.gMultiplet.
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


        # THE SEVEN TESTS BELOW USED TO COMPARE A TYPE WITH AN INSTANCE -- `typeof(x) == Coulex()` and its six
        # siblings -- which is false for every x, so all seven processes fell through to error("stop b") and were
        # unreachable.  They also tested the wrong THING: computation.processSettings is declared
        # ::Basics.AbstractProcessSettings, so it is a Settings object and never one of the AbstractProcess
        # singletons.  Each now tests the Settings type of the module whose function the branch calls, which is
        # the form the rest of this chain already uses.  Fixed 24-Aug-2026.
        elseif  typeof(computation.processSettings) == CoulombExcitation.Settings
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == CoulombIonization.Settings
            # Reserved and empty until 25-Aug-2026, when the module was written; this branch used to raise and say
            # that CoulombIonization held no computation to reach.  It now does.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = CoulombIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           CoulombIonization.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )           end
        elseif  typeof(computation.processSettings) == MultiPhotonIonization.Settings
            outcome = MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                         computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == InternalConversion.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = InternalConversion.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                           InternalConversion.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        elseif  typeof(computation.processSettings) == CrystalFieldEmission.Settings
            outcome = CrystalFieldEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("crystal-field-resolved emission lines:" => outcome) )    end
        elseif  typeof(computation.processSettings) == PhotoRecombinationInterference.Settings
            # As in PhotoRecombination: the captured electron approaches the INITIAL ion, so initialAsfSettings.scField
            # decides the automatic default here.
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.initialAsfSettings.scField)
            outcome = PhotoRecombinationInterference.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                           computation.grid, PhotoRecombinationInterference.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("photorecombination-interference pathways:" => outcome) )    end
        elseif  typeof(computation.processSettings) == PhotoIonizationInterference.Settings
            exchangex = Continuum.resolveExchangeDefault(computation.processSettings.continuumExchange, computation.finalAsfSettings.scField)
            outcome = PhotoIonizationInterference.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                           computation.grid, PhotoIonizationInterference.Settings(computation.processSettings; continuumExchange=exchangex))
            if output    results = Base.merge( results, Dict("photoionization-interference pathways:" => outcome) )    end
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
