
"""
`module  ResultKeys`
    ... a small namespace holding the keys for those results of `Basics.perform` that have NO settings type of
        their own: the multiplets and the radial grid.  Every other result is asked for by the settings type that
        produced it -- `res[PhotoEmission.Settings]`, `res[Hfs.Settings]` -- which needs no new name at all, since
        the caller has just written that type three lines above as `processSettings = PhotoEmission.Settings(...)`.

        WHY A MODULE RATHER THAN SIX NAMES IN `Basics`, which is the question this file exists to answer.  The two
        natural names are already taken by the DATA: `ManyElectron.Multiplet` and `Radial.Grid`.  A key called
        `Basics.Multiplet` would therefore sit one qualifier away from the datum of the same name -- two different
        things, one word apart -- which is exactly the near-miss priority item 39 was filed to remove.  The only
        alternative inside `Basics` is to coin names nobody will guess (`LevelMultiplet`, `ComputationGrid`), and
        a name that must be looked up defeats the purpose.  A namespace disambiguates without inventing anything:
        `ResultKeys.Multiplet` is a key, `Multiplet` is the datum.

        The keys are TYPES, not instances, so that the whole interface reads by one rule -- a key is a bare type
        qualified by the module that owns it.
"""
module ResultKeys

using  ..Basics

"""
`struct  ResultKeys.Multiplet`             ... the multiplet of a structure computation;  was `"multiplet:"`.
"""
struct  Multiplet              <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.Grid`                  ... the radial grid used;  was `"grid:"`.
"""
struct  Grid                   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.InitialMultiplet`      ... the initial-state multiplet of a process computation;  was
                                               `"initialMultiplet"` -- note that this one carried NO colon.
"""
struct  InitialMultiplet       <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.FinalMultiplet`        ... the final-state multiplet of a process computation;  was
                                               `"finalMultiplet"`, likewise without a colon.
"""
struct  FinalMultiplet         <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.IntermediateMultiplet` ... the intermediate-state multiplet, where the process has one;  was
                                               `"intermediateMultiplet"`, likewise without a colon.
"""
struct  IntermediateMultiplet  <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.IjfMultiplet`          ... the IJF-coupled multiplet of a hyperfine computation;  was
                                               `"IJF multiplet:"`.  That branch is presently disabled in
                                               `Basics.perform`, and the key is carried so that re-enabling it
                                               needs no new name.
"""
struct  IjfMultiplet           <: Basics.AbstractResultKey   end


# ---------------------------------------------------------------------------------------------------------------
# THE CASCADE AND SIMULATION RESULTS.  A Cascade computation is driven by a SCHEME rather than by a settings type,
# so its results have nothing to be keyed by and need names here -- unlike a process result, which is asked for by
# the settings that produced it.  Where a Cascade result IS the same quantity as a process result -- its
# photoemission, autoionization and photoexcitation lines -- it keeps the process key and is NOT given one here;
# see `BasicsAZ.resultKeyStrings`, which lists every string variant a quantity is stored under.
# ---------------------------------------------------------------------------------------------------------------

"""
`struct  ResultKeys.CascadeData`            ... the cascade's blocks and lines;  was `"cascade data:"`.
"""
struct  CascadeData            <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.CascadeScheme`          ... the scheme the cascade was run under;  was `"cascade scheme"`.
"""
struct  CascadeScheme          <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.Name`                   ... the computation's name;  was `"name"` OR `"name:"` -- the same
                                                quantity under two strings, which is why it is one key here.
"""
struct  Name                   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.DataFormat`             ... the format the cascade data was written in;  was `"data format:"`.
"""
struct  DataFormat             <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.InitialMultiplets`      ... the multiplets a cascade starts from;  was `"initial multiplets:"`.
                                                PLURAL, and distinct from `InitialMultiplet`, which is the single
                                                initial-state multiplet of a process computation.
"""
struct  InitialMultiplets      <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.GeneratedMultiplets`    ... the multiplets a cascade generated;  was
                                                `"generated multiplets:"`.
"""
struct  GeneratedMultiplets    <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.DielectronicMultiplets` ... was `"dielectronic multiplets:"`.
"""
struct  DielectronicMultiplets <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.ImpactExcitedMultiplets` ... was `"impact-excited multiplets:"`.
"""
struct  ImpactExcitedMultiplets <: Basics.AbstractResultKey  end

"""
`struct  ResultKeys.PhotoExcitedMultiplets` ... was `"photoexcited multiplets:"`.
"""
struct  PhotoExcitedMultiplets <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.PhotoIonizedMultiplets` ... was `"photoionized multiplets:"`.
"""
struct  PhotoIonizedMultiplets <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.DielectronicCaptureLines` ... was `"dielectronic-capture lines:"`.  NOT the same quantity as
                                                `DielectronicRecombination.Settings`, which is the full
                                                recombination;  this is the capture step alone.
"""
struct  DielectronicCaptureLines <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.HollowIonLineData`      ... was `"hollow-ion line data:"`.
"""
struct  HollowIonLineData      <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.PhotoRecombinationLineData` ... was `"photo-recombination line data:"`.  LINE DATA, i.e. the
                                                cascade's own record, not the `PhotoRecombination.Settings` lines.
"""
struct  PhotoRecombinationLineData <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.PhotoExcitationLineData` ... was `"photoexcitation line data:"`.
"""
struct  PhotoExcitationLineData <: Basics.AbstractResultKey  end

"""
`struct  ResultKeys.SimulationData`         ... the data a cascade simulation consumed;  was `"data:"`.
"""
struct  SimulationData         <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SimulationProperty`     ... the property a cascade simulation computed;  was `"property:"`.
"""
struct  SimulationProperty     <: Basics.AbstractResultKey   end


# ---------------------------------------------------------------------------------------------------------------
# THE `generate(::Representation)` RESULTS.  A representation is not driven by settings either, so these need
# names here, like the cascade's.
# ---------------------------------------------------------------------------------------------------------------

"""
`struct  ResultKeys.MeanFieldBasis`         ... was `"mean-field basis"`.
"""
struct  MeanFieldBasis         <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.MeanFieldMultiplet`     ... was `"mean-field multiplet"`.
"""
struct  MeanFieldMultiplet     <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.MeanPotential`          ... was `"mean potential"`.
"""
struct  MeanPotential          <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.Orbitals`               ... was `"orbitals"`.
"""
struct  Orbitals               <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.ReferenceMultiplet`     ... was `"reference multiplet"`.
"""
struct  ReferenceMultiplet     <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.CiMultiplet`            ... was `"CI multiplet"`.
"""
struct  CiMultiplet            <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.GreenChannels`          ... was `"Green channels"`.
"""
struct  GreenChannels          <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.EiiCrossSections`       ... was `"EII cross sections:"`.
"""
struct  EiiCrossSections       <: Basics.AbstractResultKey   end


# ---------------------------------------------------------------------------------------------------------------
# PLASMA.SUPERCONFIGURATIONOPACITYSCHEME'S RESULTS.  Driven by a scheme, not a settings type, exactly the
# Cascade/generate reason above -- so these need names here too, rather than a sixth vocabulary of raw strings.
# ---------------------------------------------------------------------------------------------------------------

"""
`struct  ResultKeys.SuperConfigurationGroups`  ... the `SuperConfigurations.Group` per reference configuration;
                                                     was `"superconfiguration groups"`.
"""
struct  SuperConfigurationGroups             <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationPopulations`  ... the Boltzmann population fraction per group, matching
                                                          `SuperConfigurationGroups` in order; was
                                                          `"superconfiguration populations"`.
"""
struct  SuperConfigurationPopulations        <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationTransitionArrays`  ... one `SuperConfigurations.TransitionArray` per
                                                               E1-connected pair of reference configurations;
                                                               was `"superconfiguration transition arrays"`.
"""
struct  SuperConfigurationTransitionArrays   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationTransitionGf`  ... each array's own oscillator-strength sum
                                                           `gf = (2/3)*meanEnergy*totalStrength`, matching
                                                           `SuperConfigurationTransitionArrays` in order; was
                                                           `"superconfiguration transition gf"`.
"""
struct  SuperConfigurationTransitionGf       <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationOpacityContributions`  ... that `gf` weighted by the lower configuration's
                                                                   population; was
                                                                   `"superconfiguration opacity contributions"`.
"""
struct  SuperConfigurationOpacityContributions <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.SuperConfigurationOpacitySpectrum`  ... the Stage-3 summed per-absorber cross section
                                                              kappa(omega), present only when
                                                              `scheme.opacityOmegas` was non-empty; was
                                                              `"superconfiguration opacity spectrum"`.
"""
struct  SuperConfigurationOpacitySpectrum    <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationIonizationGroups`  ... one `SuperConfigurations.Group` array per
                                                               ionization stage; was
                                                               `"superconfiguration ionization groups"`.
"""
struct  SuperConfigurationIonizationGroups   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationIonizationCharges`  ... each stage's charge, matching
                                                                `SuperConfigurationIonizationGroups` in
                                                                order; was
                                                                `"superconfiguration ionization charges"`.
"""
struct  SuperConfigurationIonizationCharges  <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationIonizationStagePopulations`  ... the equilibrium population fraction
                                                                         per stage, summing to 1; was
                                                                         `"superconfiguration ionization
                                                                         stage populations"`.
"""
struct  SuperConfigurationIonizationStagePopulations <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.SuperConfigurationIonizationGroupPopulations`  ... the same equilibrium population,
                                                                         resolved per GROUP within each
                                                                         stage, summing to 1 over ALL stages
                                                                         together; was
                                                                         `"superconfiguration ionization
                                                                         group populations"`.
"""
struct  SuperConfigurationIonizationGroupPopulations <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.SuperConfigurationIonizationElectronDensity`  ... the converged self-consistent
                                                                        electron number density [a.u.]; was
                                                                        `"superconfiguration ionization
                                                                        electron density"`.
"""
struct  SuperConfigurationIonizationElectronDensity <: Basics.AbstractResultKey end

"""
`struct  ResultKeys.SuperConfigurationMixtureSpeciesOpacities`  ... one Stage-3 `kappa(omega)` array per
                                                                        species, in `scheme`'s own order; was
                                                                        `"superconfiguration mixture species
                                                                        opacities"`.
"""
struct  SuperConfigurationMixtureSpeciesOpacities   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureOpacitySpectrum`  ... the number-fraction-weighted SUM of every
                                                                       species' own spectrum; was
                                                                       `"superconfiguration mixture opacity
                                                                       spectrum"`.
"""
struct  SuperConfigurationMixtureOpacitySpectrum    <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixturePlanckOpacity`  ... the Planck mean of the BOUND-BOUND-ONLY
                                                                     mixture spectrum (unaffected by
                                                                     `boundFreeApprox`, so its meaning never
                                                                     changes); was `"superconfiguration
                                                                     mixture planck opacity"`. See
                                                                     `ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum`
                                                                     for the Planck mean of the FULL
                                                                     (bound-bound+bound-free) spectrum once
                                                                     bound-free is turned on.
"""
struct  SuperConfigurationMixturePlanckOpacity      <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureBoundFreeOpacity`  ... added 18-Oct-2026: the mixture's own
                                                                        bound-free (photoionization) kappa(omega),
                                                                        on the SAME `opacityOmegas` grid, summed
                                                                        over every adjacent-stage pair of every
                                                                        species and weighted the same way the
                                                                        bound-bound spectrum is. Present ONLY
                                                                        when `scheme.boundFreeApprox !== nothing`;
                                                                        was `"superconfiguration mixture
                                                                        bound-free opacity"`. See
                                                                        `ResultKeys.SuperConfigurationMixtureFreeFreeOpacity`
                                                                        for the SEPARATE free-free contribution.
"""
struct  SuperConfigurationMixtureBoundFreeOpacity   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureFreeFreeOpacity`  ... added 19-Oct-2026: the mixture's own
                                                                        free-free (inverse-bremsstrahlung)
                                                                        kappa(omega), on the SAME `opacityOmegas`
                                                                        grid, summed over EVERY stage (no
                                                                        adjacency needed -- a stage's own ion
                                                                        charge alone sets its free-free cross
                                                                        section) of every species and weighted
                                                                        the same way the bound-bound spectrum
                                                                        is. Present ONLY when
                                                                        `scheme.freeFreeIncluded`; was
                                                                        `"superconfiguration mixture free-free
                                                                        opacity"`.
"""
struct  SuperConfigurationMixtureFreeFreeOpacity    <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum`  ... added 18-Oct-2026, GENERALIZED
                                                                            19-Oct-2026: bound-bound PLUS
                                                                            whichever of bound-free/free-free
                                                                            are switched on, i.e.
                                                                            `SuperConfigurationMixtureOpacitySpectrum`
                                                                            plus `...BoundFreeOpacity` (if
                                                                            `scheme.boundFreeApprox !== nothing`)
                                                                            plus `...FreeFreeOpacity` (if
                                                                            `scheme.freeFreeIncluded`), on the
                                                                            same `opacityOmegas` grid. Present
                                                                            whenever EITHER continuum switch is
                                                                            on; this is the spectrum both
                                                                            `Plasma.planckMeanOpacity` and
                                                                            `Plasma.rosselandMeanOpacity` are
                                                                            evaluated on for
                                                                            `...TotalPlanckOpacity` and
                                                                            `...RosselandOpacity` below.
"""
struct  SuperConfigurationMixtureTotalOpacitySpectrum <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureTotalPlanckOpacity`  ... added 18-Oct-2026, GENERALIZED
                                                                           19-Oct-2026: the Planck mean of the
                                                                           TOTAL spectrum (bound-bound plus
                                                                           whichever of bound-free/free-free
                                                                           are on) -- distinct from
                                                                           `SuperConfigurationMixturePlanckOpacity`,
                                                                           which stays bound-bound-only so its
                                                                           own meaning never changes. Present
                                                                           whenever EITHER continuum switch is on.
"""
struct  SuperConfigurationMixtureTotalPlanckOpacity <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureRosselandOpacity`  ... added 18-Oct-2026, GENERALIZED
                                                                         19-Oct-2026: the Rosseland mean of the
                                                                         TOTAL spectrum (bound-bound plus
                                                                         whichever of bound-free/free-free are
                                                                         on) -- MUCH LESS floor-sensitive than
                                                                         bound-bound alone once a continuum
                                                                         contribution exists (measured stable
                                                                         across 9 orders of magnitude of the
                                                                         internal floor choice on a real case,
                                                                         vs. 12 orders of movement for
                                                                         bound-bound alone), but NOT
                                                                         unconditionally floor-free -- see
                                                                         `Plasma.displayMixtureSpectrumBreakdown`'s
                                                                         own docstring. Present whenever EITHER
                                                                         continuum switch is on; remains a
                                                                         LOWER BOUND whenever the OTHER one is
                                                                         still off.
"""
struct  SuperConfigurationMixtureRosselandOpacity   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureCharges`  ... one charge array per species (matching
                                                               `refConfigsPerSpeciesStages`'s own stage
                                                               order), added 18-Oct-2026 when the mixture
                                                               scheme gained its own stage axis; was
                                                               `"superconfiguration mixture charges"`.
"""
struct  SuperConfigurationMixtureCharges            <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureStagePopulations`  ... one stage-population array per species,
                                                                        summing to 1 WITHIN that species
                                                                        (not across the whole mixture); was
                                                                        `"superconfiguration mixture stage
                                                                        populations"`.
"""
struct  SuperConfigurationMixtureStagePopulations   <: Basics.AbstractResultKey   end

"""
`struct  ResultKeys.SuperConfigurationMixtureElectronDensity`  ... the ONE converged electron number density
                                                                       [a.u.] shared by the WHOLE mixture
                                                                       (every species, every stage); was
                                                                       `"superconfiguration mixture electron
                                                                       density"`.
"""
struct  SuperConfigurationMixtureElectronDensity    <: Basics.AbstractResultKey   end


# ---------------------------------------------------------------------------------------------------------------
# WHAT IS DELIBERATELY NOT HERE: the RAS ladder's per-step results, stored as "step1", "step2", ... one per layer.
# Their number is not known until the ladder is built, so no fixed set of types can name them and they KEEP THEIR
# STRINGS.  That is why the deprecation warning fires only where a typed key actually replaces the string: a
# string with no replacement is not deprecated, it is the only way to ask.
# ---------------------------------------------------------------------------------------------------------------

# NOTHING IS EXPORTED, AND THAT IS THE POINT OF THE MODULE.  Exporting `Multiplet` and `Grid` would let a
# `using ResultKeys` shadow `ManyElectron.Multiplet` and `Radial.Grid` in the caller's own namespace -- which is
# precisely the collision this module exists to prevent.  The keys are meant to be written qualified, and the
# qualifier is what makes them unambiguous:  res[ResultKeys.Multiplet]  against  Multiplet the datum.

end # module
