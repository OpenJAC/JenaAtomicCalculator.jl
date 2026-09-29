
println("Cp) Apply & test the NuclearShielding module: how the electrons of a closed-shell ion modify, at the")
println("    nucleus, an electric field gradient applied from outside.")

using Printf

# WRITTEN 29-Sep-2026, together with module-NuclearShielding.jl.
#
# WHAT THE MODULE ANSWERS, in one sentence: what does the nucleus actually FEEL, given what is applied from outside?
# It is the mirror image of `Hfs`.  There a NUCLEAR moment couples to the electrons, giving the hyperfine constants
# A, B and C.  Here an EXTERNAL field of the same multipolarity pushes on the electrons, the cloud is distorted, and
# the distorted cloud makes a field of its own AT THE NUCLEUS.  The ratio is the shielding factor, which is why the
# two modules carry the same multipole selectors.
#
# WHY IT MATTERS.  For a heavy ion the electrons do not screen the applied gradient, they AMPLIFY AND INVERT it -- by
# a factor near 185 for Th(4+).  So a nuclear quadrupole splitting measured in a crystal is about 99 % an atomic
# quantity and 1 % crystallography, and the factor is the whole bridge between a lattice geometry and an observable.
# The same bridge, at rank 1 and magnetic, converts an NMR frequency into a nuclear magnetic moment -- which is how a
# wrong shielding constant once produced an apparent 7-sigma failure of strong-field QED in bismuth.
#
# WHAT IS AND IS NOT COMPUTED.  A shielding factor is a LINEAR response of the FREE ION to a field assumed UNIFORM
# over it.  This module returns exactly that, relativistically, with no crystal or molecular content: no charge
# transfer, no bonding, and (so far) no self-consistency -- the induced field is not yet allowed to act back.
# `displayResults` states all of this on every run, beside the number, rather than leaving it to a docstring.


if  true
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- Th(4+) [Rn] gives gamma_inf = -184.95, hence 1 - gamma_inf = 185.95.
    #   VERIFIED AGAINST AN INDEPENDENT RELATIVISTIC CALCULATION: Feiock & Johnson, Phys. Rev. 187, 39 (1969),
    #   quoted by Peik & Okhapkin, C. R. Physique 16, 516 (2015), give about -178 for Rn-like Th(4+) -- so 4 %,
    #   from a different method, different code and a different decade.  The breakdown printed below also shows
    #   the outermost p shell dominating, which is the expected structure for a closed-shell ion.
    println("\nCp-a)  Th(4+): the quadrupole antishielding factor, and where it comes from.")
    setDefaults("print summary: open", "zzz-NuclearShielding.sum")
    nm   = Nuclear.Model(90.0)
    conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false)
    setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    wa   = Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                              configs=[conf], asfSettings=asf)
    mp   = perform(wa; output=true)["multiplet:"]
    NuclearShielding.computeOutcomes(mp, nm, grid, NuclearShielding.Settings())
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- the three CLOSED-SHELL ions of Sternheimer, Phys. Rev. 159, 266 (1967),
    #   Table I.  This work against his RAW values (his published ones carry a hand-applied factor 1/1.1 for his
    #   use of neutral-atom orbitals, which we do not need since we use ion orbitals):
    #       Y(3+)   -36.2  against  -34.1    +6 %
    #       In(3+)  -27.4  against  -25.3    +8 %
    #       Bi(3+)  -53.4  against  -42.9   +24 %
    #   THE TREND IS THE RESULT, not the individual numbers: the excess grows with Z, which is the shape of a
    #   relativistic correction -- his calculation is non-relativistic and ours is not.  Note also that his
    #   l -> l+-2 channels were ESTIMATED from a Thomas-Fermi model and are 5-10 % of his total, while we compute
    #   them properly, so agreement much better than ~15 % would not be meaningful in either direction.
    println("\nCp-b)  The validation set: three closed-shell ions with published values.")
    for (name, Z, cstr) in (("Y^3+", 39.0, "[Kr]"), ("In^3+", 49.0, "[Kr] 4d^10"),
                            ("Bi^3+", 83.0, "[Xe] 4f^14 5d^10 6s^2"))
        conf = Configuration(cstr);   nm = Nuclear.Model(Z)
        grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
        asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
        mp   = perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                       configs=[conf], asfSettings=asf); output=true)["multiplet:"]
        out  = NuclearShielding.computeOutcomes(mp, nm, grid, NuclearShielding.Settings(); output=true)
        @printf("\n  >>> %-8s gamma_inf = %9.2f\n\n", name, out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- A NEGATIVE RESULT, AND THAT IS WHY THE BRANCH EXISTS.  The far-field
    #   coefficient assumes the external charges are infinitely distant.  In CaF2 the nearest neighbours sit
    #   2.4-3.9 A away, comparable with the ion, so the uniform-field assumption was expected to fail and to
    #   explain a factor 1.55-1.69 between the computed 186 and the 110 extracted from measurements.
    #   IT DOES NOT.  Replacing r^2 by the exact q r_<^2 / r_>^3 of a charge at the real distance changes the
    #   answer by 0.0-0.5 %, because Th(4+) is small and tightly bound and almost none of its density lies
    #   beyond 2.4 A.  So the factorisation into lattice x atomic factor is SOUND for this system, and the
    #   110-vs-186 discrepancy must be sought elsewhere -- in the assigned defect charges or in the screening
    #   this module does not yet include.  The sanity check (a charge at 1000 A) reproduces the far field exactly.
    println("\nCp-c)  Does the near neighbour break the uniform-field assumption?  (It does not.)")
    nm   = Nuclear.Model(90.0);   conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    mp   = perform(Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                   configs=[conf], asfSettings=asf); output=true)["multiplet:"]
    ANG  = 1.0/0.529177210903
    for model in (NuclearShielding.UniformField(),
                  NuclearShielding.NeighbourField([(1.0, 1000.0*ANG)]),
                  NuclearShielding.NeighbourField([( 2.0, 3.8627*ANG)]),
                  NuclearShielding.NeighbourField([( 1.0, 2.3654*ANG)]))
        out = NuclearShielding.computeOutcomes(mp, nm, grid,
                  NuclearShielding.Settings(NuclearShielding.Settings(); fieldModel=model); output=true)
        @printf("\n  >>> %-46s gamma = %9.2f\n\n", string(model), out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- the sum over intermediate states converges FAST and on BOUND virtuals:
    #   keeping only the largest 1, 2, 5, 10, 20 terms per channel reaches the full answer to better than 1 %
    #   by about 20 states out of ~180 available.  That is the practical argument for the sum-over-states
    #   formulation: the answer can be SEEN to come from a handful of low-lying orbitals rather than from the
    #   continuum, so the quality of the pseudo-continuum is not the limiting factor here.
    println("\nCp-d)  How many intermediate states does the response actually need?")
    nm   = Nuclear.Model(90.0);   conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    mp   = perform(Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                   configs=[conf], asfSettings=asf); output=true)["multiplet:"]
    for nv in (1, 2, 5, 10, 20, 0)
        out = NuclearShielding.computeOutcomes(mp, nm, grid,
                  NuclearShielding.Settings(NuclearShielding.Settings(); nVirtualMax=nv); output=true)
        @printf("\n  >>> nVirtualMax = %3d   gamma_inf = %9.2f\n\n", nv, out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- A GUARD TEST: both calls below MUST raise, and raising is the pass.
    #   (i)  an OPEN-SHELL ion.  A shielding factor is the SCALAR response of a spherical ion; for an open shell
    #        the response is not one number and the factorisation does not apply.  Returning something plausible
    #        would be worse than refusing, so the module refuses and says why.
    #   (ii) calcM1.  The magnetic-dipole shielding sigma -- the NMR analogue, and the quantity behind the
    #        bismuth hyperfine puzzle -- is NOT implemented.  Both rank-1 operators already exist in JAC
    #        (InteractionStrength.zeeman_n1 and hfs_tM1) but the response is not wired up or validated, so the
    #        settings refuse rather than return a number that would be wrong rather than approximate.
    println("\nCp-e)  Two things the module refuses to do, and says why.")
    nm   = Nuclear.Model(90.0)
    for (what, cstr, sets) in (("an OPEN-SHELL ion", "[Rn] 5f^2", NuclearShielding.Settings()),
                               ("calcM1 (not implemented)", "[Rn]",
                                NuclearShielding.Settings(NuclearShielding.Settings(); calcM1=true)))
        conf = Configuration(cstr)
        grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
        asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
        try
            mp = perform(Atomic.Computation(Atomic.Computation(); name="x", grid=grid, nuclearModel=nm,
                         configs=[conf], asfSettings=asf); output=true)["multiplet:"]
            NuclearShielding.computeOutcomes(mp, nm, grid, sets)
            @printf("\n  >>> %-28s DID NOT RAISE -- this is a FAILURE of the guard.\n\n", what)
        catch err
            @printf("\n  >>> %-28s raised as designed:\n       %s\n\n", what,
                    first(split(sprint(showerror, err), "\n"), 4) |> x -> join(x, "\n       "))
        end
    end
    #
end
