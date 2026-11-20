#
# crossCheckBenchmark.jl -- TIER 0 of the standing cross-code benchmark (priority item 23).
#
#   Run from the JAC root:   julia --project=. tools/crossCheckBenchmark.jl
#
# WHAT THIS IS FOR.  Every check JAC applies to itself is self-referential: the suite compares against its own
# approved files, the cancellation factor and the gauge ratio are computed from the amplitudes they judge, and a
# `Last successful` date records one person's reading of JAC's own output.  The remedy is a standing comparison
# against ANOTHER code, and GRASP2018 is the first partner -- same formalism, so it can be given the SAME CSF
# list, which makes a disagreement a BUG rather than a modelling difference.
#
# BUT A CODE-AGAINST-CODE COMPARISON NEEDS A TOLERANCE, AND THE TOLERANCE MUST NOT BE GUESSED.  That is what this
# script measures, and it is why it comes FIRST and runs JAC alone.  For a ONE-ELECTRON ion in a POINT-nucleus
# field the Dirac energy is known in closed form, so each code can be measured against THEORY instead of against
# the other.  The worst deviation found here is the floor below which an agreement between two codes proves
# nothing:  if JAC is 1e-6 from the exact value, then matching GRASP to 1e-6 is matching two discretizations and
# not two physics implementations.
#
# WHY A ONE-ELECTRON CASE AND NOT A GRID CHECK.  `Bsplines.checkGridRepresentation` already compares JAC's
# B-spline solution with `Basics.computeDiracEnergy` and is the Rule 12 guard;  it returns a worst deviation but
# tests hydrogenic proxies INSIDE an SCF, not the end-to-end answer a user gets.  This script goes through
# `Atomic.Computation` + `perform`, which is the route every example and every application uses, so what it
# measures is what a user would publish.
#
# THE FIELD MUST BE THE BARE NUCLEUS, AND THIS IS THE ONE TRAP IN THE SETUP.  A DFS or HS mean field builds its
# exchange from the density, so on a SINGLE electron it is a self-interaction and the answer is NOT the exact
# Dirac energy -- it would look like a discretization error of the wrong size.  `Basics.NuclearField()` is the
# only correct choice here, with `NoneQed()` and a point nucleus, because those are the terms the closed form
# contains and nothing else.
#
# WHAT THIS WOULD MISS.  It measures the radial discretization of a ONE-ELECTRON problem and nothing else.  It is
# blind to every error in the electron-electron interaction, in the angular algebra, and in the CSF machinery --
# which is most of what a many-electron code can get wrong.  So it sets a FLOOR for Tier 1 and is not itself a
# validation of JAC:  a clean table here is necessary and nowhere near sufficient.
#
# TWO CALIBRATIONS, MEASURED 20-Nov-2026, AND BOTH CHANGE HOW THE TABLE MUST BE READ.
#
#  1. THE REFERENCE VALUE IS THE WEAKER SIDE AT LOW Z.  `Basics.computeDiracEnergy` ends on
#     `m c^2 (1/wa - 1)`, and at Z = 1 that subtracts 0.9999734 from 1, discarding five significant digits:
#     against a 200-bit evaluation of the SAME formula the Float64 one carries 2.62e-12 at Z = 1, 1.02e-14 at
#     Z = 18, and 1.85e-16 at Z = 92 -- it IMPROVES with Z as the cancellation eases, while the code's own error
#     grows.  The two cross between Z = 1 and Z = 18.  So this script compares against a BigFloat evaluation and
#     prints the Float64 formula's own error beside it;  reading the raw Z = 1 agreement as JAC's accuracy would
#     credit the code with the formula's round-off, which is precisely the comparison-tuned-to-agree failure.
#
#  2. WHAT IS LEFT IS GOVERNED BY `rnt`, AND NO GRID GUARD LOOKS AT IT.  At Z = 92 the 1s came out 1.275e-4 Ha
#     (2.62e-8) above the exact value, and that offset is INVARIANT under a factor 8 in `hp` and under a tenfold
#     box with nsL going 82 -> 476 -- it is converged, not unresolved.  It is controlled by the INNERMOST grid
#     point: rnt/10 cuts it by a factor 30, to 4.3e-6 Ha.  The large component goes as r^gamma with
#     gamma = sqrt(1-(Z alpha)^2) = 0.7411 at Z = 92, a NON-INTEGER power that a polynomial B-spline basis cannot
#     follow at the origin, so the error lives at r -> 0 where neither `hp` nor `rbox` refines.
#     **CONSEQUENCE FOR `Basics.checkGridConvergence`: it refines `hp` with the box held, so it is BLIND to this
#     class and will report `converged` on a value carrying a systematic offset.**  Not a defect of it -- it was
#     built for the outer region -- but it is what this instrument would otherwise have been trusted to catch.
#     NOT yet established: where the optimum `rnt` lies.  Shrinking `rnt` at fixed `NoPoints` also shrinks the
#     box (0.196 -> 0.113 a.u. over the decades tried), so those rows are confounded exactly as the register
#     records for `hp`;  settling it needs `NoPoints` raised to hold the box.
#
# THE TIER-1 TOLERANCE IS A SEPARATE MEASUREMENT AND THE INSTRUMENT ALREADY EXISTS:
# `Basics.checkGridConvergence(comp, extract)` re-performs a computation with `hp` refined and THE BOX HELD.  Do
# not measure convergence by refining `hp` alone (`Radial.Grid(grid; hp=)` keeps `NoPoints` fixed, so the box
# SHRINKS), and do not use points-per-oscillation as a criterion -- it is measured non-monotone.
#
using JenaAtomicCalculator, Printf
setprecision(BigFloat, 200)

"""
`diracEnergyBig(sh::Subshell, Z::Float64)`
    ... evaluates the SAME closed-form point-nucleus Dirac expression as `Basics.computeDiracEnergy`, in 200-bit
        arithmetic, and returns the energy::BigFloat in atomic units without the rest energy.  It exists because
        the Float64 form loses five digits to the `1/wa - 1` cancellation at low Z (2.62e-12 at Z = 1), which is
        larger than the deviation being measured there;  the reference must not be the weaker side of its own
        comparison.
"""
function diracEnergyBig(sh::Subshell, Z::Float64)
    A = BigFloat(Defaults.getDefaults("alpha"));   C = BigFloat(Defaults.getDefaults("speed of light: c"))
    jPlusHalf = BigFloat(Basics.subshell_2j(sh) + 1) / 2;    nr = BigFloat(sh.n) - jPlusHalf;   ZZ = BigFloat(Z)
    g  = sqrt(jPlusHalf^2 - ZZ^2 * A^2)
    wa = sqrt(1 + ZZ^2 * A^2 / (nr + g)^2)

    return( C^2 * (1/wa - 1) )
end

# Each case is one configuration of a ONE-ELECTRON ion, so every level in its multiplet corresponds to exactly
# one subshell and the pairing is by the angular momentum j, not by energy order.
cases = [ ("H      Z =  1",   1.0, "1s", ["1s_1/2"]),
          ("H      Z =  1",   1.0, "2p", ["2p_1/2", "2p_3/2"]),
          ("Ar17+  Z = 18",  18.0, "1s", ["1s_1/2"]),
          ("Ar17+  Z = 18",  18.0, "2s", ["2s_1/2"]),
          ("Ar17+  Z = 18",  18.0, "2p", ["2p_1/2", "2p_3/2"]),
          ("Xe53+  Z = 54",  54.0, "1s", ["1s_1/2"]),
          ("Xe53+  Z = 54",  54.0, "2p", ["2p_1/2", "2p_3/2"]),
          ("U91+   Z = 92",  92.0, "1s", ["1s_1/2"]),
          ("U91+   Z = 92",  92.0, "2p", ["2p_1/2", "2p_3/2"]) ]

settings = AsfSettings(AsfSettings(); generateScf=true, scField=Basics.NuclearField(),
                       eeInteraction=CoulombInteraction(), eeInteractionCI=CoulombInteraction(),
                       qedModel=NoneQed())

function runTier0(cases, settings)
    println("\n" * "="^104)
    println("  TIER 0 of the cross-code benchmark:  JAC through perform() against the EXACT point-nucleus Dirac energy")
    println("  one electron, bare nuclear field, no Breit, no QED -- so the closed form is the whole answer")
    println("="^104)
    @printf("  %-14s %-9s %18s %18s %11s %11s   %s\n",
            "ion", "subshell", "E(JAC) [a.u.]", "E(Dirac) [a.u.]", "JAC dev.", "formula", "rbox, rnt, NoPoints")
    println("  " * "-"^102)

    worst = 0.0;   worstWhat = "";   rows = 0

    for (label, Z, confString, subshellStrings) in cases
        configs = [Configuration(confString)]
        nm      = Nuclear.Model(Z, Nuclear.PointNucleus())
        grid    = Basics.recommendedGrid(configs, nm)
        comp    = Atomic.Computation(Atomic.Computation(); name="tier0", nuclearModel=nm, grid=grid,
                                     configs=configs, asfSettings=settings)
        multiplet = nothing
        try
            results   = redirect_stdout(devnull) do;  Basics.perform(comp; output=true)   end
            multiplet = results[ResultKeys.Multiplet]
        catch err
            @printf("  %-14s %-9s   RAISED: %s\n", label, confString, first(sprint(showerror, err), 60))
            continue
        end

        for shString in subshellStrings
            sh    = Subshell(shString)
            jSh   = Basics.subshell_j(sh)
            exact   = diracEnergyBig(sh, Z)
            f64     = Basics.computeDiracEnergy(sh, Z)
            devF64  = Float64(abs(BigFloat(f64) - exact) / abs(exact))
            # Pair the level to the subshell by j.  A one-electron configuration has at most one level per j, so
            # a missing match is a real failure and is reported rather than silently skipped.
            idx   = findall(lev -> lev.J == jSh, multiplet.levels)
            if  length(idx) != 1
                @printf("  %-14s %-9s   NO UNIQUE LEVEL with J = %s (%d found)\n",
                        label, shString, string(jSh), length(idx))
                continue
            end
            en  = multiplet.levels[idx[1]].energy
            dev = Float64(abs(BigFloat(en) - exact) / abs(exact))
            rows += 1
            if  dev > worst   worst = dev;   worstWhat = "$label  $shString"   end
            @printf("  %-14s %-9s %18.10f %18.10f %11.2e %11.2e   %.2f, %.0e, %d\n",
                    label, shString, en, Float64(exact), dev, devF64, grid.r[end], grid.rnt, grid.NoPoints)
        end
    end

    println("  " * "-"^102)
    @printf("  %d rows.  WORST relative deviation %.2e, at %s.\n", rows, worst, worstWhat)
    println("  The `formula` column is the Float64 reference's OWN error against 200-bit arithmetic;  where it is")
    println("  comparable to the JAC column, that row measures the REFERENCE and not the code -- see the header.")
    println("  THAT IS THE TIER-1 FLOOR: an agreement between two codes closer than this compares discretizations,")
    println("  not implementations.  It bounds a ONE-ELECTRON radial problem only -- see the header.")
    println("="^104 * "\n")

    return( (rows = rows, worst = worst, worstWhat = worstWhat) )
end

runTier0(cases, settings)
