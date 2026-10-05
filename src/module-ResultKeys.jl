
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


# NOTHING IS EXPORTED, AND THAT IS THE POINT OF THE MODULE.  Exporting `Multiplet` and `Grid` would let a
# `using ResultKeys` shadow `ManyElectron.Multiplet` and `Radial.Grid` in the caller's own namespace -- which is
# precisely the collision this module exists to prevent.  The keys are meant to be written qualified, and the
# qualifier is what makes them unambiguous:  res[ResultKeys.Multiplet]  against  Multiplet the datum.

end # module
