
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

# NOTHING IS EXPORTED, AND THAT IS THE POINT OF THE MODULE.  Exporting `Multiplet` and `Grid` would let a
# `using ResultKeys` shadow `ManyElectron.Multiplet` and `Radial.Grid` in the caller's own namespace -- which is
# precisely the collision this module exists to prevent.  The keys are meant to be written qualified, and the
# qualifier is what makes them unambiguous:  res[ResultKeys.Multiplet]  against  Multiplet the datum.

end # module
