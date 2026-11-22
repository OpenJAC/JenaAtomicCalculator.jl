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
#
# TIER 1 -- ADDED 20-Nov-2026:  the same MANY-ELECTRON Dirac-Hartree-Fock level against GRASP2018.
#
# This is the first test that exercises the ELECTRON-ELECTRON INTERACTION and the ANGULAR ALGEBRA,
# to which Tier 0 is completely blind.  Everything that could differ between the two codes is
# pinned:  a POINT nucleus (two finite-nucleus models are not identical, and the difference would
# read as a code difference), Coulomb only with NO Breit and NO QED, and a SINGLE reference
# configuration so that the CSF list is trivially identical in both codes.  JAC uses
# `Basics.ALField()`, the average-level Dirac-Fock field -- NOT `DFSField()`, which is the local
# Dirac-Fock-Slater approximation and a different operator -- and for one CSF that is the same
# functional GRASP's `rmcdhf` minimises.
#
# THE GRIDS ARE DELIBERATELY *NOT* MADE IDENTICAL, AND THAT IS A CORRECTION TO THE FIRST PLAN.
# GRASP's default point-nucleus grid is RNT = exp(-65/16)/Z, H = 0.0625, HP = 0 and N = 220, which
# on a purely exponential mesh puts its outer boundary near 175 a.u. at Z = 92 -- against the 0.196
# a.u. that `Basics.recommendedGrid` chooses there.  Copying GRASP's four numbers into JAC would
# STARVE its fixed B-spline basis, which is Rule 12's own warning about a box far larger than the
# orbitals.  The four parameters simply do not mean the same thing to a finite-difference
# integrator, which carries a point at every radius, as they do to a basis-set solver.  So each
# code is run on a grid adequate FOR IT and the comparison rests on the Hamiltonian and the CSF
# list being identical, which is what makes a disagreement a bug.
#
# THE GRASP NUMBERS ARE NOT TYPED IN.  `tools/graspLevels.sh` is called for every row, so the
# partner's value is re-derived on every run together with the oracle's commit -- a hand-recorded
# number that nothing re-measures is this project's most repeated failure.
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

#####  THE ROW SCHEMA  #####

"""
`struct BenchmarkRow`
    ... ONE comparison, in the form the whole benchmark is keyed by.  The field names are part of the
        deliverable and not an implementation detail:  the same row is printed for a reader, written to
        `docs/src/validation-rows.csv` for any other harness to read, and rendered into the published page, so
        they are expected to stay stable.

        THE ORDER OF THE KEYS IS THE ARRANGEMENT, and it is deliberate.  `tier` comes first because it says HOW
        a row must be read -- the same 1e-7 difference means "a bug" at tier 1 and "expected, and explained" at
        tier 2 -- and a table that mixes the two says nothing.  `quantity` comes second because that is what a
        reader actually asks about ("can I trust JAC's E1 rates?"), not which partner was used.  And the PARTNER
        IS A FIELD RATHER THAN A SECTION, because the most valuable rows are those where two partners cover one
        quantity, and grouping by code would separate exactly those.

    + tier           ::Int64     ... 0 against a closed form, 1 same formalism (a difference is a BUG), 2 a
                                     different method (a difference is informative, never decisive).
    + quantity       ::String    ... what is compared, e.g. "level-energy".
    + system         ::String    ... the ion, e.g. "U91+".
    + Z              ::Float64
    + configuration  ::String
    + level          ::String    ... the level or subshell the row is about.
    + model          ::String    ... nucleus, mean field, ee-interaction and QED, in one string.
    + grid           ::String    ... rnt, h, hp, NoPoints, rbox -- the discretization the JAC value was got on.
    + jacValue       ::Float64
    + partner        ::String    ... "exact-dirac-formula" at tier 0, otherwise the code's name.
    + partnerVersion ::String    ... the pinned commit, or "closed form".  A number whose code state is unknown
                                     is not evidence.
    + partnerValue   ::Float64
    + relDiff        ::Float64
    + tolerance      ::Float64
    + verdict        ::String    ... AGREES, DIFFERS, REFERENCE-LIMITED, PARTNER-LIMITED, or RAISED.
                                     The two "-LIMITED" words both say that a deviation is understood and that the
                                     OTHER side carries it:  REFERENCE-LIMITED at tier 0, where the Float64 closed
                                     form is as inaccurate as the deviation, and PARTNER-LIMITED at tier 1, where
                                     the partner code is.  Neither is ever DERIVED by this harness -- see
                                     `PARTNER_LIMITED` below for why attribution is not a tool's to make.
    + note           ::String    ... one sentence where the verdict needs one;  empty otherwise.
"""
struct BenchmarkRow
    tier           ::Int64
    quantity       ::String
    system         ::String
    Z              ::Float64
    configuration  ::String
    level          ::String
    model          ::String
    grid           ::String
    jacValue       ::Float64
    partner        ::String
    partnerVersion ::String
    partnerValue   ::Float64
    relDiff        ::Float64
    tolerance      ::Float64
    verdict        ::String
    note           ::String
end

const CSV_FIELDS = ["tier", "quantity", "system", "Z", "configuration", "level", "model", "grid",
                    "jac_value", "partner", "partner_version", "partner_value", "rel_diff", "tolerance",
                    "verdict", "note"]

"""
`csvCell(x)`
    ... renders one field for the machine-readable file, quoting a string that contains a comma or a quote and
        giving a float enough digits to be re-checked rather than merely read;  a String is returned.
"""
csvCell(x::AbstractString) = occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : String(x)
csvCell(x::Integer)        = string(x)
csvCell(x::AbstractFloat)  = @sprintf("%.12g", x)


function runTier0(cases, settings)
    println("\n" * "="^104)
    println("  TIER 0 of the cross-code benchmark:  JAC through perform() against the EXACT point-nucleus Dirac energy")
    println("  one electron, bare nuclear field, no Breit, no QED -- so the closed form is the whole answer")
    println("="^104)
    @printf("  %-14s %-9s %18s %18s %11s %11s   %s\n",
            "ion", "subshell", "E(JAC) [a.u.]", "E(Dirac) [a.u.]", "JAC dev.", "formula", "rbox, rnt, NoPoints")
    println("  " * "-"^102)

    worst = 0.0;   worstWhat = "";   rows = 0;   out = BenchmarkRow[];   tol = 1.0e-7

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
            # TIER 0 MEASURES AND DOES NOT GRADE, and the first version of this got it wrong.  Giving tier 0 a
            # pass/fail threshold is circular: the tier-1 tolerance is READ OFF this tier's worst row, so
            # judging tier 0 by that same number stamps the one informative row `AGREES` and hides it.  The
            # U91+ 1s offset of 1.275e-4 Ha was reported as agreement exactly so, on 20-Nov-2026.  A tolerance
            # of 0 in a tier-0 row therefore means "no threshold: this row is a measurement".
            # REFERENCE-LIMITED is the separate statement that the Float64 reference carries as much error as
            # the deviation, which is true below about Z = 18 and is NOT a defect of the code.
            verdict = devF64 >= dev ? "REFERENCE-LIMITED" : "MEASURED"
            push!(out, BenchmarkRow(0, "level-energy", replace(strip(label), r"\s+" => " "), Z, confString, shString,
                                    "point nucleus; bare nuclear field; Coulomb; no QED",
                                    @sprintf("rnt=%.1e h=%.4g hp=%.4g N=%d rbox=%.3f",
                                             grid.rnt, grid.h, grid.hp, grid.NoPoints, grid.r[end]),
                                    en, "exact-dirac-formula", "closed form", Float64(exact), dev, 0.0, verdict,
                                    verdict == "REFERENCE-LIMITED" ?
                                        "the Float64 reference formula's own error is " *
                                        @sprintf("%.2e", devF64) * ", so this row bounds the reference" : ""))
        end
    end

    println("  " * "-"^102)
    # The floor is the worst row the REFERENCE can still resolve;  a reference-limited row bounds the formula.
    floorRows = filter(r -> r.verdict != "REFERENCE-LIMITED", out)
    floor     = isempty(floorRows) ? 0.0 : maximum(r -> r.relDiff, floorRows)
    floorAt   = isempty(floorRows) ? "" :
                floorRows[argmax([r.relDiff for r in floorRows])].system * " " *
                floorRows[argmax([r.relDiff for r in floorRows])].level
    @printf("  %d rows.  WORST relative deviation %.2e, at %s.\n", rows, worst, worstWhat)
    @printf("  MEASURED FLOOR, excluding the %d reference-limited rows: %.2e, at %s.\n",
            length(out) - length(floorRows), floor, floorAt)
    println("  Tier 0 MEASURES and does not grade:  there is no pass/fail here, because tier 1's tolerance is")
    println("  read off this floor and judging tier 0 by it would be circular.")
    println("  The `formula` column is the Float64 reference's OWN error against 200-bit arithmetic;  where it is")
    println("  comparable to the JAC column, that row measures the REFERENCE and not the code -- see the header.")
    println("  THAT IS THE TIER-1 FLOOR: an agreement between two codes closer than this compares discretizations,")
    println("  not implementations.  It bounds a ONE-ELECTRON radial problem only -- see the header.")
    println("="^104 * "\n")

    return( (rows = rows, worst = worst, worstWhat = worstWhat, out = out) )
end


#####  TIER 1:  many-electron Dirac-Hartree-Fock levels against GRASP2018  #####

# Each case is ONE reference configuration, so the CSF list is a single CSF in both codes.  The
# GRASP columns name the configuration in GRASP's own nonrelativistic notation and its active set.
tier1cases = [ ("He-like  Z = 18", 18.0, "1s^2",            "1s(2,*)",                  "1s",       "0,0", 0, Basics.plus),
               ("He-like  Z = 54", 54.0, "1s^2",            "1s(2,*)",                  "1s",       "0,0", 0, Basics.plus),
               ("Be-like  Z = 18", 18.0, "1s^2 2s^2",       "1s(2,*)2s(2,*)",           "1s,2s",    "0,0", 0, Basics.plus),
               ("Be-like  Z = 54", 54.0, "1s^2 2s^2",       "1s(2,*)2s(2,*)",           "1s,2s",    "0,0", 0, Basics.plus),
               ("Ne-like  Z = 26", 26.0, "1s^2 2s^2 2p^6",  "1s(2,*)2s(2,*)2p(6,*)",    "1s,2s,2p", "0,0", 0, Basics.plus),
               ("Ne-like  Z = 54", 54.0, "1s^2 2s^2 2p^6",  "1s(2,*)2s(2,*)2p(6,*)",    "1s,2s,2p", "0,0", 0, Basics.plus),
               # THE OPEN-SHELL ROWS, added 21-Nov-2026, and they are the point rather than more of the same.
               # Every case above is a CLOSED shell with J = 0, so until now the oracle had never compared an
               # open-shell level at all -- and an open shell is where the AVERAGE-LEVEL functional, and which
               # levels it averages over, can differ between two codes that agree perfectly on a closed one.
               # A 2J = 1 ground term of 1s^2 2s is still a SINGLE CSF, so tier 1's "the same CSF list" stays
               # automatic and the mean field is the only thing left that can differ.  This is the cheap test
               # of whether the velocity-gauge E1 disagreement at high Z is a setup difference.
               ("Li-like  Z = 26", 26.0, "1s^2 2s",         "1s(2,*)2s(1,*)",           "1s,2s",    "1,1", 1, Basics.plus),
               ("Li-like  Z = 54", 54.0, "1s^2 2s",         "1s(2,*)2s(1,*)",           "1s,2s",    "1,1", 1, Basics.plus) ]

"""
`graspOracleCommit()`
    ... reads the pinned commit of the GRASP oracle clone, so that every partner value carries the code state
        it came from;  a short hash::String is returned, or "unknown" where the clone cannot be read.
"""
function graspOracleCommit()
    dir = get(ENV, "GRASP_BINDIR", joinpath(homedir(), "grasp-oracle", "bin"))
    return( try
        strip(read(`git -C $(dirname(dir)) rev-parse --short HEAD`, String))
    catch
        "unknown"
    end )
end


"""
`graspLevelEnergy(Z::Float64, conf::String, active::String, jrange::String)`
    ... runs `tools/graspLevels.sh` for one case and returns the first level energy::Float64 in
        Hartree, or `nothing` where GRASP raised or printed nothing parsable.  The partner's value
        is re-derived on every run rather than stored, so it cannot go stale against the oracle.
"""
function graspLevelEnergy(Z::Float64, conf::String, active::String, jrange::String)
    script = joinpath(@__DIR__, "graspLevels.sh")
    isfile(script)  ||  return( nothing )
    out = try
        read(`bash $script $(Int(Z)) $conf $active $jrange`, String)
    catch
        return( nothing )
    end
    for line in split(out, '\n')
        startswith(strip(line), "#")  &&  continue
        m = match(r"(-?\d+\.\d+)[DdEe]([+-]\d+)", line)
        isnothing(m)  &&  continue

        return( parse(Float64, m.captures[1]) * 10.0^parse(Int, m.captures[2]) )
    end

    return( nothing )
end


function runTier1(cases, settings)
    println("\n" * "="^104)
    println("  TIER 1:  many-electron Dirac-Hartree-Fock levels, JAC against GRASP2018")
    println("  point nucleus, Coulomb only, no Breit, no QED, ONE reference configuration -- so the CSF list")
    println("  is a single CSF in both codes and only the implementations differ")
    println("="^104)
    @printf("  %-16s %-16s %19s %19s %12s %10s\n",
            "ion", "configuration", "E(JAC) [a.u.]", "E(GRASP) [a.u.]", "rel. diff", "verdict")
    println("  " * "-"^102)

    worst = 0.0;   worstWhat = "";   rows = 0;   tol = 1.0e-7;   out = BenchmarkRow[]
    oracle = graspOracleCommit()

    for (label, Z, jacConf, graspConf, active, jrange, twoJ, parity) in cases
        configs = [Configuration(jacConf)]
        nm      = Nuclear.Model(Z, Nuclear.PointNucleus())
        grid    = Basics.recommendedGrid(configs, nm)
        comp    = Atomic.Computation(Atomic.Computation(); name="tier1", nuclearModel=nm, grid=grid,
                                     configs=configs, asfSettings=settings)
        enJac = try
            results = redirect_stdout(devnull) do;  Basics.perform(comp; output=true)   end
            results[ResultKeys.Multiplet].levels[1].energy
        catch err
            @printf("  %-16s %-16s   JAC RAISED: %s\n", label, jacConf, first(sprint(showerror, err), 50))
            continue
        end
        enGrasp = graspLevelEnergy(Z, graspConf, active, jrange)
        if  isnothing(enGrasp)
            @printf("  %-16s %-16s %19.10f   GRASP gave no parsable energy\n", label, jacConf, enJac)
            continue
        end
        dev  = abs(enJac - enGrasp) / abs(enGrasp)
        rows += 1
        if  dev > worst   worst = dev;   worstWhat = "$label  $jacConf"   end
        verdict = dev <= tol ? "AGREES" : "DIFFERS"
        @printf("  %-16s %-16s %19.10f %19.10f %12.2e %10s\n",
                label, jacConf, enJac, enGrasp, dev, verdict)
        push!(out, BenchmarkRow(1, "level-energy", replace(strip(label), r"\s+" => " "), Z, jacConf, "2J = " * jrange,
                                "point nucleus; average-level Dirac-Fock; Coulomb; no Breit; no QED; 1 CSF",
                                @sprintf("JAC rnt=%.1e h=%.4g hp=%.4g N=%d rbox=%.3f; GRASP on its own default",
                                         grid.rnt, grid.h, grid.hp, grid.NoPoints, grid.r[end]),
                                enJac, "GRASP2018", oracle, enGrasp, dev, tol, verdict,
                                verdict == "DIFFERS" ?
                                    "the Hamiltonian and the CSF list are identical by construction, so this " *
                                    "is a defect in one of the two codes" : ""))
    end

    println("  " * "-"^102)
    @printf("  %d rows at a tolerance of %.0e, read off the Tier-0 floor.  WORST %.2e, at %s.\n",
            rows, tol, worst, worstWhat)
    println("  A row marked DIFFERS is a BUG in one code or the other, NOT a modelling difference --")
    println("  the Hamiltonian and the CSF list are identical by construction.")
    println("="^104 * "\n")

    return( (rows = rows, worst = worst, worstWhat = worstWhat, out = out) )
end

# accuracyScf IS SET EXPLICITLY, AND LEAVING IT AT THE DEFAULT WAS A REAL DEFECT OF THIS HARNESS
# (found 21-Nov-2026).  `SelfConsistent.solveAverageLevelField`'s own docstring says the default 1e-6
# is an overlap defect, QUADRATIC in the orbital change, so it "accepts orbitals that are still moving
# by 1.4e-3" and leaves the argon energy 6.2e-05 Ha short.  A cross-code comparison must not be limited
# by one side's stopping criterion:  until this was set, the oracle was measuring JAC's SCF truncation
# and calling it a disagreement with GRASP.  Measured, going from 1e-6 to 1e-10:
#     He-like Z = 18   transition energy  3.34e-06 -> 1.32e-08,  A_Coulomb 2.2e-04 -> 2.5e-05,
#                      A_Babushkin 6.6e-05 -> 1.000000 exactly
#     Li-like Z = 26   transition energy  3.59e-04 -> 5.1e-05,   A_Coulomb 6.6e-02 -> 1.2e-02,
#                      A_Babushkin 1.6e-03 -> 2.5e-04
# THIS IS NOT TUNING JAC TOWARD GRASP, which Rule 17 step 3a forbids:  the energy and BOTH gauges
# improve together, which a fitted knob does not do, and converging an SCF is the right choice here
# whatever the partner says.  1e-12 buys nothing over 1e-10 (measured, both cases).
tier1settings = AsfSettings(AsfSettings(); generateScf=true, scField=Basics.ALField(),
                            eeInteraction=CoulombInteraction(), eeInteractionCI=CoulombInteraction(),
                            qedModel=NoneQed(), accuracyScf=1.0e-10,
                            scfRoute=Basics.AverageLevelRoute(2000))


#####  TIER 1:  E1 transition rates against GRASP2018  #####

# WHY THE FIGURE OF MERIT IS THE CODE-AGAINST-CODE RATIO PER GAUGE, AND NOT THE GAUGE RATIO.
# The register's section G settles this and it is worth not re-deriving: the Babushkin/Coulomb
# ratio "warns on a neutral system and is BLIND on a highly charged one", it is undefined for a
# transition INSIDE one configuration (the velocity form of an electric multipole is structurally
# zero between two orbitals of the same subshell), and a detector watching it move is reading
# noise.  So each gauge is compared SEPARATELY between the two codes, and the gauge ratio is
# carried only as a secondary note.
#
# AND THE TRANSITIONS MUST CONNECT DIFFERENT CONFIGURATIONS, for the same reason.
#
# TWO TRAPS IN ONE FIELD, BOTH MEASURED HERE ON 21-Nov-2026.
#
#  (1) `Einstein.Line.photonRate` is in ATOMIC UNITS, not inverse seconds -- read as 1/s it is low
#      by 1/AU_TIME = 4.13e16.
#
#  (2) AND IT ONCE WAS NOT THE EINSTEIN A COEFFICIENT -- A DEFECT THIS HARNESS FOUND IN `src/`.
#      Until 21-Nov-2026 `module-Einstein.jl` formed
#          wa = 8 pi alpha omega / (2J_initial+1) * (2J_final+1)
#      while `module-PhotoEmission.jl` formed the same quantity WITHOUT the last factor.  The
#      second is the standard Einstein A, so one field name carried two quantities differing by
#      (2J_final+1), and photonRate/AU_TIME equalled A only when the LOWER level had J = 0.  That
#      is why this harness looked correct on the He-like cases (J_final = 0) and reported a
#      spurious factor 2.13 on the Li-like ones (J_final = 1/2).
#      FIXED in commit 24469284:  the two modules now agree and photonRate IS the A coefficient in
#      atomic units.  So do NOT divide by (2J_final+1) here again -- doing so would restore the
#      very error this paragraph records, and on the He-like rows it would again be invisible.
const AU_TIME = 2.4188843265e-17        # s;  one atomic unit of time

"""
`einsteinA(line)`
    ... returns the Einstein A coefficients of `line` in 1/s as the NamedTuple `(coulomb, babushkin)`.
        `Einstein.Line.photonRate` is the A coefficient in ATOMIC UNITS, so only the unit conversion
        is applied;  see trap (2) in the note above for why nothing else may be.
"""
function einsteinA(line)
    return( (coulomb   = line.photonRate.Coulomb  /AU_TIME,
             babushkin = line.photonRate.Babushkin/AU_TIME) )
end

#  THE PARTNER-LIMITED ATTRIBUTIONS -- AN EXPLICIT LIST, AND DELIBERATELY NOT A CRITERION.
#
#  A row is PARTNER-LIMITED when the deviation has been traced to the PARTNER code.  That is a physics
#  judgement about somebody else's code, and this harness must not make it:  the same reasoning as Rule 7's,
#  where `Last successful` is never machine-written because reproduction proves stability and not correctness.
#  A mechanical rule -- "length agrees and velocity does not" -- would be exactly the self-certifying detector
#  that cost a day on 21-Nov-2026:  that pattern LOCALISES a difference to the velocity amplitude, it does NOT
#  say whose.  So each entry below is added by hand, with the evidence that earned it, and nothing is inferred.
#
#  WHAT THIS WOULD MISS, stated because a quiet detector gets believed:  it says nothing at all about rows not
#  listed, and it CANNOT notice that an attribution has gone stale -- if the partner is ever upgraded past the
#  pinned commit, every entry here must be re-derived rather than trusted.  The oracle commit is recorded on
#  every row for exactly that reason.
#
#  Keyed by (quantity, system);  both fine-structure lines of a case are covered by one entry.
const PARTNER_LIMITED = [
    (quantity = "E1-rate-coulomb", system = "He-like Z = 54", oracle = "9006157",
     why = "traced to GRASP2018, 21-Nov-2026, by four measurements: (1) the LENGTH form agrees to 5e-05 and is "   *
           "FLAT in Z (1.000047/1.000055 at Z = 54, 1.000049/1.000054 at Z = 92) -- and it uses the SAME orbitals " *
           "and the SAME CI mixing as the velocity form, so the orbitals agree and the many-electron assembly is "  *
           "validated; (2) JAC's one-body velocity form satisfies the EXACT gauge identity to 4.3e-08 at Z = 92; "  *
           "(3) the two codes agree in the velocity gauge to 2e-05 at Z = 18, 26 and 36, the divergence setting "   *
           "in progressively from Z ~ 47, which is numerics and not algebra; (4) (A_C/A_B - 1)*Z must be smooth "   *
           "and vanish as Z -> infinity, the orbitals becoming exactly hydrogenic -- JAC gives 0.861, 0.835, "      *
           "0.823, 0.824, 0.829, 0.847, 0.878, 0.924 over Z = 18..92 while GRASP gives 0.862, 0.835, 0.823, "       *
           "0.772, 0.631, -0.646, -5.509, -13.331, reversing sign above Z ~ 60. GRASP's SCF is converged (11 "      *
           "digits). The mechanism -- a finite-difference solver losing the SMALL components as 1s becomes "        *
           "compact, which only the velocity form reads at full weight -- is INFERRED FROM BEHAVIOUR; GRASP's "     *
           "source was not read."),
    (quantity = "E1-rate-coulomb", system = "He-like Z = 66", oracle = "9006157",
     why = "the same measurement set as the Z = 54 entry above, which swept Z = 18..92 in both codes: the "      *
           "length form agrees to 5e-05 and is flat in Z while the velocity form departs progressively from "    *
           "Z ~ 47, and (A_C/A_B - 1)*Z -- which must vanish as Z -> infinity, the orbitals becoming exactly "   *
           "hydrogenic -- reverses SIGN in the partner above Z ~ 60 while JAC stays smooth at ~0.85/Z."),
    (quantity = "E1-rate-coulomb", system = "He-like Z = 92", oracle = "9006157",
     why = "the same measurement set as the Z = 54 entry above.  At this Z the partner's velocity rate differs "  *
           "from JAC's by 5.5 % and 18 % on the two lines while the LENGTH rates agree to 4.9e-05 and 5.4e-05 "   *
           "-- and the length form uses the SAME orbitals and the SAME CI mixing, so the orbitals agree and only " *
           "the velocity amplitude does not.  JAC's own one-body gauge identity holds to 4.3e-08 here."),
                   ]

#  A CAUSE WE HAVE TRACED TO OUR OWN SIDE.  The verdict stays DIFFERS, because the DEFAULT behaviour really
#  does differ and a validation page must describe the default;  what changes is the NOTE, which says which code
#  carries it and what to switch on.  Leaving the generic "a defect in one of the two codes" on a row whose cause
#  is known would be true and useless.  Hand-maintained for the same reason as `PARTNER_LIMITED`:  attributing a
#  deviation is a judgement and not something a harness may infer.
const KNOWN_CAUSE = [
    (quantity = "E1-rate-coulomb", system = "Li-like Z = 26",
     why = "traced to JAC, 21/22-Nov-2026, and remediable:  the average-level field does not optimise the "      *
           "ROTATION between two subshells of the same kappa carrying unequal generalised occupations (here 1s " *
           "at 2.0 and 2s at 1/3), so its orbitals come out orthonormal but not stationary in that angle.  With " *
           "`SelfConsistent.setAlSameKappaRotation(true)` this row reads 0.999998 instead, the functional "      *
           "falling 1.75e-05 Ha to within 4.3e-07 Ha of the partner.  The default is OFF because switching it "  *
           "on moves the eighth or ninth digit of many approved references at once;  see the note at "           *
           "`SelfConsistent.GBL_AL_SAME_KAPPA_ROTATION`."),
                    ]

"""
`knownCause(quantity::AbstractString, system::AbstractString)`
    ... returns the recorded cause of a DIFFERS row where it has been traced BY HAND, or `nothing` where it has
        not.  The verdict is deliberately left as DIFFERS:  the default behaviour does differ, and a validation
        page describes the default rather than the best setting reachable.
"""
function knownCause(quantity::AbstractString, system::AbstractString)
    for  e  in  KNOWN_CAUSE
        if  e.quantity == quantity  &&  e.system == system    return( e.why )    end
    end

    return( nothing )
end


"""
`partnerLimited(quantity::AbstractString, system::AbstractString, oracle::AbstractString)`
    ... says whether a DIFFERS row has been attributed BY HAND to the partner code, returning the justification
        string, or `nothing` where it has not.  An attribution is only honoured for the oracle commit it was
        made against, so upgrading the partner retires every attribution instead of silently carrying it over.
"""
#  AbstractString throughout:  `oracle` arrives as a SubString from `strip(read(...))`.
function partnerLimited(quantity::AbstractString, system::AbstractString, oracle::AbstractString)
    for  e  in  PARTNER_LIMITED
        if  e.quantity == quantity  &&  e.system == system  &&  e.oracle == oracle    return( e.why )    end
    end

    return( nothing )
end


e1cases = [("He-like  Z = 18", 18.0, ["1s^2", "1s 2p"], "1s(2,*);1s(1,*)2p(1,*)", "1s,2p",    "0,4"),
           ("He-like  Z = 26", 26.0, ["1s^2", "1s 2p"], "1s(2,*);1s(1,*)2p(1,*)", "1s,2p",    "0,4"),
           ("He-like  Z = 54", 54.0, ["1s^2", "1s 2p"], "1s(2,*);1s(1,*)2p(1,*)", "1s,2p",    "0,4"),
           # Z = 66 AND Z = 92 ARE HERE TO SHOW A TREND, not to lengthen the table.  The partner's
           # velocity-gauge E1 rate degrades PROGRESSIVELY above Z ~ 40, and a single high-Z row would read as
           # an isolated oddity;  the onset and the endpoint together show it is a trend with a direction.  The
           # length-gauge rows beside them agree to 5e-05 at every Z and are what prove the orbitals agree, so
           # each pair of rows carries its own control.
           ("He-like  Z = 66", 66.0, ["1s^2", "1s 2p"], "1s(2,*);1s(1,*)2p(1,*)", "1s,2p",    "0,4"),
           ("He-like  Z = 92", 92.0, ["1s^2", "1s 2p"], "1s(2,*);1s(1,*)2p(1,*)", "1s,2p",    "0,4"),
           ("Li-like  Z = 26", 26.0, ["1s^2 2s", "1s^2 2p"], "1s(2,*)2s(1,*);1s(2,*)2p(1,*)",
                                                             "1s,2s,2p", "1,3")]

"""
`graspE1Lines(Z::Float64, confs::String, active::String, jrange::String)`
    ... runs `tools/graspTransitions.sh` and returns the E1 lines it reports as a vector of
        NamedTuples `(energy, aCoulomb, aBabushkin)` with the energy in cm^-1 and the rates in
        1/s;  an empty vector where GRASP produced nothing parsable.
"""
function graspE1Lines(Z::Float64, confs::String, active::String, jrange::String)
    script = joinpath(@__DIR__, "graspTransitions.sh")
    isfile(script)  ||  return( NamedTuple[] )
    out = try
        read(`bash $script $(Int(Z)) $confs $active $jrange`, String)
    catch
        return( NamedTuple[] )
    end
    rows = NamedTuple[]
    for line in split(out, '\n')
        startswith(strip(line), "#")  &&  continue
        f = split(strip(line))
        length(f) == 7  ||  continue
        push!(rows, (energy = parse(Float64, f[5]), aCoulomb = parse(Float64, f[6]),
                     aBabushkin = parse(Float64, f[7])))
    end

    return( rows )
end


"""
`runTier1E1(cases, settings)`
    ... computes the E1 lines of each case in JAC, pairs them with GRASP's by TRANSITION ENERGY,
        and compares the two rates in each gauge separately;  returns the rows for the schema.

        The pairing is by energy because the two codes order and index their levels differently.  The
        window is wide (2e-3) rather than tight, and that is deliberate: measured 21-Nov-2026, the two
        codes' transition energies agree to 9e-10 for the He-like cases but only to 3.5e-04 for the
        Li-like ones, so a window set from the closed-shell agreement would drop the open-shell lines
        entirely.  The window is still far tighter than the spacing between any two lines of these
        cases, so the match stays unambiguous -- and the energy itself is now COMPARED as its own
        `transition-energy` row rather than merely used for pairing and discarded, which is what let a
        3.5e-04 disagreement sit unseen underneath the rate rows it partly explains.
"""
function runTier1E1(cases, settings)
    println("\n" * "="^104)
    println("  TIER 1:  E1 transition rates, JAC against GRASP2018")
    println("  one orbital set in both codes, so no biorthogonal transformation enters either side;")
    println("  each GAUGE compared separately -- the Babushkin/Coulomb ratio is NOT the figure of merit")
    println("="^104)
    @printf("  %-16s %11s %11s %9s %12s %12s %9s %12s %12s %9s\n", "ion",
            "E JAC", "E GRASP", "dE", "A_C JAC", "A_C GRASP", "ratio", "A_B JAC", "A_B GRASP", "ratio")
    println("  " * "-"^102)

    out = BenchmarkRow[];   tol = 1.0e-3;   worst = 0.0;   worstWhat = "";   oracle = graspOracleCommit()
    # THE TRANSITION ENERGY IS COMPARED TOO, AND AT THE LEVEL TIER'S OWN 1e-7 -- not at the rate tolerance.
    # It is reachable: the tier-1 levels agree to ~2e-10 relative, i.e. ~6e-08 Ha absolute, so a transition
    # energy of 1.8 Ha inherits ~6e-08 relative and one of 115 Ha far less.  Where this row says DIFFERS the
    # two codes' ORBITALS differ, and that must be settled before the rate rows beside it can be read at all,
    # since A_B ~ omega^3 while A_C ~ omega^1 -- the same energy error moves the two gauges by different
    # amounts and so shows up as a false gauge-selective disagreement.  Do NOT widen it to match the rates.
    tolE = 1.0e-7

    for (label, Z, jacConfs, graspConfs, active, jrange) in cases
        configs = [Configuration(c) for c in jacConfs]
        nm      = Nuclear.Model(Z, Nuclear.PointNucleus())
        grid    = Basics.recommendedGrid(configs, nm)
        comp    = Atomic.Computation(Atomic.Computation(); name="tier1e1", nuclearModel=nm, grid=grid,
                                     configs=configs, asfSettings=settings,
                                     propertySettings=[Einstein.Settings([E1], false, LineSelection(),
                                                                         0., 0., 1.0e6)])
        jacLines = try
            res = redirect_stdout(devnull) do;  Basics.perform(comp; output=true)   end
            res[Einstein.Settings]
        catch err
            @printf("  %-16s   JAC RAISED: %s\n", label, first(sprint(showerror, err), 50))
            continue
        end
        gLines = graspE1Lines(Z, graspConfs, active, jrange)
        if  isempty(gLines)
            @printf("  %-16s   GRASP produced no parsable E1 line\n", label);   continue
        end

        for g in gLines
            # pair by transition energy;  1e-4 relative is far tighter than any line spacing here
            best = nothing;  bestd = Inf
            for l in jacLines
                e = l.omega*219474.6313702
                d = abs(e - g.energy)/g.energy
                if  d < bestd   bestd = d;  best = l   end
            end
            # A PAIRING FAILURE IS REPORTED, NOT SKIPPED.  The first version of this dropped any
            # GRASP line it could not match and printed nothing, so the Li-like doublet vanished
            # from the table and the run looked clean -- a detector that goes quiet when it fails
            # is the worst kind.  The window is 2e-3 rather than 1e-4 because a 2s-2p transition
            # energy is a SMALL DIFFERENCE of two large totals:  the levels agree to 1e-6, the
            # splitting between them to only a few 1e-4.  It is still far tighter than the
            # spacing between any two lines of these cases, so the match stays unambiguous.
            if  best === nothing  ||  bestd > 2.0e-3
                @printf("  %-16s %10.0f   NO JAC LINE within 2e-3 of this GRASP line (closest %.1e)\n",
                        label, g.energy, bestd)
                continue
            end
            aj = einsteinA(best);   jc = aj.coulomb;   jb = aj.babushkin
            (jc <= 0.0 || g.aCoulomb <= 0.0)  &&  continue
            rc = jc/g.aCoulomb;   rb = jb/g.aBabushkin
            jacE = best.omega * 219474.6313702
            devE = abs(jacE/g.energy - 1.0)
            @printf("  %-16s %11.0f %11.0f %9.1e %12.5e %12.5e %9.6f %12.5e %12.5e %9.6f\n",
                    label, jacE, g.energy, devE, jc, g.aCoulomb, rc, jb, g.aBabushkin, rb)
            if  devE > worst   worst = devE;   worstWhat = "$label  transition energy"   end
            push!(out, BenchmarkRow(1, "transition-energy", replace(strip(label), r"\s+" => " "), Z,
                                    join(jacConfs, " + "), @sprintf("%.0f cm^-1 (JAC)", jacE),
                                    "point nucleus; average-level Dirac-Fock; Coulomb; no Breit; no QED; one orbital set",
                                    @sprintf("JAC rnt=%.1e h=%.4g hp=%.4g N=%d rbox=%.3f; GRASP on its own default",
                                             grid.rnt, grid.h, grid.hp, grid.NoPoints, grid.r[end]),
                                    jacE, "GRASP2018", oracle, g.energy, devE, tolE,
                                    devE <= tolE ? "AGREES" : "DIFFERS",
                                    devE > tolE ? "a transition energy is a SMALL DIFFERENCE of two large " *
                                        "totals, so it is the sharpest one-CSF test of the two codes' orbitals; " *
                                        "settle this row before reading the rate rows beside it" : ""))
            for (gauge, jv, gv, rr) in (("coulomb", jc, g.aCoulomb, rc), ("babushkin", jb, g.aBabushkin, rb))
                dev = abs(rr - 1.0)
                if  dev > worst   worst = dev;   worstWhat = "$label  $gauge"   end
                quantity = "E1-rate-" * gauge
                system   = replace(strip(label), r"\s+" => " ")
                #  A DIFFERS row may have been attributed BY HAND to the partner;  nothing here infers it.
                why      = dev > tol ? partnerLimited(quantity, system, oracle) : nothing
                verdict  = dev <= tol ? "AGREES" : (isnothing(why) ? "DIFFERS" : "PARTNER-LIMITED")
                mine     = dev > tol ? knownCause(quantity, system) : nothing
                note     = dev <= tol ? "" :
                           !isnothing(why)  ? why  :
                           !isnothing(mine) ? mine :
                           "the CSF list and the Hamiltonian are identical by construction, " *
                           "so this is a defect in one of the two codes"
                push!(out, BenchmarkRow(1, quantity, system, Z,
                                        join(jacConfs, " + "), @sprintf("%.0f cm^-1", g.energy),
                                        "point nucleus; average-level Dirac-Fock; Coulomb; no Breit; no QED; one orbital set",
                                        @sprintf("JAC rnt=%.1e h=%.4g hp=%.4g N=%d rbox=%.3f; GRASP on its own default",
                                                 grid.rnt, grid.h, grid.hp, grid.NoPoints, grid.r[end]),
                                        jv, "GRASP2018", oracle, gv, dev, tol, verdict, note))
            end
        end
    end

    println("  " * "-"^102)
    @printf("  %d rows;  rates at a tolerance of %.0e, transition energies at %.0e.  WORST deviation = %.2e, at %s.\n",
            length(out), tol, tolE, worst, worstWhat)
    println("  The JAC column is the Einstein A from `einsteinA()`:  photonRate is the A coefficient")
    println("  in atomic units, so the unit conversion is the only thing applied to it.")

    return( (out = out, worst = worst) )
end


#####  OUTPUTS:  the machine-readable file, and the published page  #####

"""
`writeBenchmarkOutputs(rows::Vector{BenchmarkRow}; write::Bool=false)`
    ... with `write = true`, writes `docs/src/validation-rows.csv` (one line per comparison, with the field
        names of `CSV_FIELDS`) and rewrites the region between the markers in `docs/src/validation.md`;  with
        `write = false` it forms the rows and writes NEITHER, so a diagnostic run leaves `docs/` untouched.
        Nothing outside the markers is touched;  that text is prose and stays under human control.  The path of
        the file written is returned.

        WHY A MACHINE-READABLE FILE AND NOT ONLY A TABLE.  A page can be read by a person and by nobody else.
        The point of a standing cross-code comparison is that ANOTHER project's harness can pick ours up, add
        its own column and compare without asking us anything, so the row carries its own model, its own grid,
        its own tolerance and the partner's pinned commit, under names that are expected not to change.
"""
function writeBenchmarkOutputs(rows::Vector{BenchmarkRow}; write::Bool=false)
    docs = normpath(joinpath(@__DIR__, "..", "docs", "src"))
    csv  = joinpath(docs, "validation-rows.csv")
    # THE CSV IS GATED BEHIND `--write` TOO, AND IT WAS NOT UNTIL 21-Nov-2026.  It is a TRACKED file, so an
    # ungated write meant that every diagnostic run of this harness -- and there were a dozen in one afternoon
    # -- left `docs/` dirty with numbers nobody had decided to publish, inviting exactly the accidental commit
    # that "docs only at release" exists to prevent.  The page was already gated;  the two now behave alike,
    # which is also what keeps the CSV and the page from disagreeing with each other.
    if  write
        open(csv, "w") do io
            println(io, join(CSV_FIELDS, ","))
            for r in rows
                println(io, join(csvCell.([r.tier, r.quantity, r.system, r.Z, r.configuration, r.level, r.model,
                                           r.grid, r.jacValue, r.partner, r.partnerVersion, r.partnerValue,
                                           r.relDiff, r.tolerance, r.verdict, r.note]), ","))
            end
        end
        println("\n>> ", csv, "  --  ", length(rows), " rows")
    else
        println("\n>> ", length(rows), " rows formed;  ", csv, " NOT written (pass --write).")
    end

    commit = try   strip(read(`git rev-parse --short HEAD`, String))   catch;  "unknown"   end
    md = IOBuffer()
    println(md, "<!-- BEGIN generated: cross-code validation -->")
    for tier in sort(unique(r.tier for r in rows))
        these = filter(r -> r.tier == tier, rows)
        println(md, "\n#### Tier ", tier, " — ",
                tier == 0 ? "against a closed form (this tier measures the floor; it does not grade)" :
                tier == 1 ? "same formalism, same CSF list: a difference is a defect" :
                            "a different method: a difference is informative, not decisive", "\n")
        println(md, "| system | level | quantity | JAC | partner | partner value | rel. diff | verdict |")
        println(md, "|---|---|---|---|---|---|---|---|")
        for r in these
            println(md, @sprintf("| %s | %s | %s | %.10g | %s `%s` | %.10g | %.2e | %s |",
                                 r.system, r.level, r.quantity, r.jacValue, r.partner, r.partnerVersion,
                                 r.partnerValue, r.relDiff, r.verdict))
        end
    end
    #  PARTNER-LIMITED IS COUNTED SEPARATELY FROM DIFFERS, because the two say opposite things to a reader:
    #  DIFFERS is an open question, PARTNER-LIMITED is a closed one whose answer was that the partner carries
    #  the deviation.  Folding them together would hide the open rows among the settled ones, which is the
    #  failure the tiers exist to prevent.
    nDiff = count(r -> r.verdict == "DIFFERS", rows)
    nPart = count(r -> r.verdict == "PARTNER-LIMITED", rows)
    print(md, "\n", length(rows), " comparisons, ", nDiff, " marked `DIFFERS`")
    nPart > 0  &&  print(md, " and ", nPart, " `PARTNER-LIMITED` (understood, and traced to the partner -- ",
                         "see the note on each such row)")
    println(md, ".  Measured ", Libc.strftime("%d %B %Y", time()), " at commit `", commit,
            "`; regenerate with `julia --project=. tools/crossCheckBenchmark.jl --write`.")
    println(md, "<!-- END generated: cross-code validation -->")
    block = String(take!(md))

    if  write
        page = joinpath(docs, "validation.md")
        isfile(page)  ||  error("writeBenchmarkOutputs(): $page does not exist;  create the page with its two " *
                                "markers before running with --write.")
        text = read(page, String)
        a = findfirst("<!-- BEGIN generated: cross-code validation -->", text)
        b = findfirst("<!-- END generated: cross-code validation -->", text)
        if  a === nothing  ||  b === nothing
            error("writeBenchmarkOutputs(): the markers are missing from $page;  add them around the table " *
                  "before running with --write.")
        end
        Base.write(page, text[1:first(a)-1] * block * text[last(b)+1:end])
        println(">> ", page, "  updated between the markers.")
    else
        println("\n", block)
        println("\n>> run with --write to replace the region between the markers in docs/src/validation.md.")
    end

    return( csv )
end


t0 = runTier0(cases, settings)
t1 = runTier1(tier1cases, tier1settings)
t2 = runTier1E1(e1cases, tier1settings)
writeBenchmarkOutputs(vcat(t0.out, t1.out, t2.out); write = ("--write" in ARGS))
