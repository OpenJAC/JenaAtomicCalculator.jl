#
# tools/probeTeaRegularIntegral.jl
#
# STEP 2 of building R(eps), the regular part of the continuum branch of the two-electron Auger amplitude
# (priority items 93 and 98).  Written 26-Oct-2026.  A diagnostic, not a test.
#
# Step 1 (tools/probeTeaRegularPart.jl) established that the two-sided continuum matrix element
# <f_core, eps_f kappa_f | V | g_core, eps kappa> can be formed, is hermitian to 3.8e-10, and genuinely depends
# on eps.  THIS probe asks the next question, and only that one:
#
#     WHAT DOES THE PRINCIPAL-VALUE INTEGRAL OVER R(eps) DO TO THE RATE, and is it converged in the mesh?
#
# THREE RUNS, ONE VARIABLE EACH.
#
#     A   NoIntermediateEnergies = 0, gCoreMultiplet set   -> MUST reproduce the committed delta-piece result.
#                                                             This is the Rule 19 baseline for the two NEW
#                                                             Settings fields: with R(eps) switched off they
#                                                             must change nothing at all.
#     B   NoIntermediateEnergies = n, for a few n          -> R(eps) added.  The SEQUENCE in n is the result;
#                                                             a single value of n is not a number worth quoting.
#
# WHY THE SEQUENCE AND NOT A VALUE.  A plain Gauss-Legendre sum of 1/x through a pole HAS NO LIMIT -- refining
# the mesh merely puts nodes nearer the pole -- and PhotoDoubleIonization measured exactly that divergence on
# helium before the principal value was subtracted: 4.35 -> 4.18 -> 4.44 -> 5.15 -> 6.35 -> 7.99 -> 9.64e-05 for
# 3...61 points, with the increments NOT shrinking (priority item 33, closed 05-Oct-2026).  The subtraction is
# implemented here from the start, so the sequence is expected to SETTLE; if it marches instead, the cause is
# the assembly and not the physics, and this probe is what distinguishes the two.
#
# AND THE RANGE IS A SECOND CONVERGENCE PARAMETER, not a technicality.  An intermediate state is VIRTUAL, so its
# energy is not restricted by energy conservation; the mesh must extend BEYOND the energetically allowed maximum
# and the pole must lie strictly INSIDE it, which is what intermediateEnergyFactor > 1 secures.  PDI's own
# default of 3 is 1.7 % below its factor-8 value, so it is scanned here too rather than trusted.
#
# COST, MEASURED FIRST AND DELIBERATELY.  A naive loop is 41 cores x ~4 kappas x (mesh + 1 pole) intermediate
# levels, each needing a continuum orbital and up to two matrix elements, PER observed partial wave and PER
# final level.  The probe therefore runs ONE final level while the cost is unknown, and reports the measured
# seconds so that the full five-level run can be sized rather than guessed.
#
# Run:  julia --project=<JAC-root> tools/probeTeaRegularIntegral.jl
#
# WHAT IT ESTABLISHED, 26-Oct-2026, final level 1 (4S), factor 3, measured at -t 1:
#
#     nInt        rate [1e-7]     vs. baseline    increment
#     0 (off)        130.3710          1.0000         --
#     3              188.2904          1.4443         --
#     5              193.9500          1.4877       +5.66
#     9              191.7247          1.4706       -2.23
#
# (A) THE BASELINE HOLDS, AND IT WAS PREDICTED BEFORE THE RUN: 130.3710 against the committed five-level
#     run-B value of 130.371 for this term.  The two new Settings fields are inert when switched off.
#
# (B) THE SEQUENCE SETTLES, WHICH IS THE WHOLE POINT OF SCANNING IT.  The increments shrink by a factor 2.5
#     AND CHANGE SIGN -- oscillatory quadrature convergence.  The failure mode this scan exists to exclude
#     looks quite different: PhotoDoubleIonization, with the pole unregularized, gave
#     4.35 -> 4.18 -> 4.44 -> 5.15 -> 6.35 -> 7.99 -> 9.64e-05 for 3...61 points, increments NOT shrinking
#     (item 33).  So the principal-value subtraction is doing its job here.
#
# (C) AND THE PHYSICS IS THE SURPRISE: R(eps) PUSHES THE RATE BACK UP BY ~47 %, i.e. in the OPPOSITE
#     direction to the delta piece of the same branch.  For this term: bound-only 180.711, bound + delta
#     130.371, bound + delta + R(eps) ~191.7.  THE TWO HALVES OF THE CONTINUUM BRANCH LARGELY CANCEL EACH
#     OTHER, and the net effect of the whole branch on 4S is a ~6 % INCREASE over bound-only rather than the
#     28 % decrease the delta piece alone suggested.
#
#     DO NOT EXTRAPOLATE THAT TO A TOTAL FROM ONE LEVEL.  The delta piece was strongly term-dependent
#     (0.62 to 0.72 across 4S/2D/2P), so R(eps) may be too; the five-level number is a separate measurement.
#
# COST, measured: 58 s with R(eps) off, 395 s at 9 points for one final level, i.e. ~37 s per mesh point.
# The full five-level run at 9 points is therefore about half an hour -- affordable, and worth doing rather
# than extrapolating.
#
# THE FIVE-LEVEL TOTAL, 26-Oct-2026, nInt 9 / factor 3 / maxKappa 4:
#
#     term     bound    +delta   +delta+R(eps)    R/delta   full/bound
#     4S     180.711   130.371         191.725     1.4706       1.0609
#     2D     266.586   176.857         298.969     1.6905       1.1215
#     2P     138.089    86.004         159.074     1.8496       1.1520
#     TOTAL  585.385   393.232         649.767     1.6524       1.1100
#
# i.e. 2.850x Zitnik's 228e-7 and 2.589x the experiment-anchored 251(46)e-7.  **THE CONTINUUM INTERMEDIATES DO
# NOT EXPLAIN THE REMAINING DISCREPANCY.**  They were the candidate held by elimination after five others were
# closed on measurement, and the branch nearly cancels itself -- ending 11 % ABOVE bound-only, so the rate sits
# where the H - E repair alone had left it.  That is the result, and it closes off the intermediate-state space.
# R(eps) is also term-dependent in the OPPOSITE order to the delta piece (1.47/1.69/1.85 against 0.72/0.66/0.62),
# so the two halves compensate term by term as well as in total and neither alone predicts the other.
#
# TWO CONVERGENCE AXES, BOTH MEASURED RATHER THAN ASSUMED.
#
#   maxKappa: FLAT -- 191.7301 / 191.7247 / 191.7247 for 2 / 4 / 6.  This matters because
#     generateIntermediateLevels truncates the INTERMEDIATE electron's partial waves at the same maxKappa as the
#     observed electron's, and that sum did not exist when the observed one was closed.  The confound is
#     deliberate: maxKappa moves BOTH sums, so a flat result proves both converged in one scan.
#
#   intermediateEnergyFactor: PARTLY RESOLUTION, PARTLY REAL, and only a density-matched test separates them.
#     At a fixed 9 points the factor scan marches -4 % per step, but widening the range at fixed point count
#     also coarsens the mesh the pole sits in.  At MATCHED density of 3.00 points per unit range: factor 3 ->
#     191.7247, factor 8 -> 178.8695, so the range costs a real 6.7 %.  A 7 % correction moves the five-level
#     total to ~600e-7, i.e. ~2.39x -- so THE VERDICT ABOVE SURVIVES THE WHOLE REMAINING UNCERTAINTY, and the
#     range scan is about quoting a number rather than about the physics.
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

# [AR-TEA-RPART-STEP2]  ONE question: what does the PV integral over R(eps) do, and does it converge in the mesh?
# The setup is the branch-f setup verbatim, including both continuum defaults -- see the note in
# tools/probeTeaContinuumIntermediates.jl for what omitting them costs (1.1 %, uniformly, and it looks physical).

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
asfSettings = AsfSettings(AsfSettings(), scField=Basics.HSField())
nm          = Nuclear.Model(18.)

gMultiplet = perform(Atomic.Computation(Atomic.Computation(), name="bound Green multiplet", grid=grid, nuclearModel=nm,
                        configs=[Configuration("1s^2 2s^2 2p^5 3s^2 3p^5"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^4")],
                        asfSettings=asfSettings); output=true)[ResultKeys.Multiplet]
cMultiplet = perform(Atomic.Computation(Atomic.Computation(), name="continuum-intermediate ionic cores", grid=grid,
                        nuclearModel=nm, configs=[Configuration("1s^2 2s^2 2p^5 3s^2 3p^4"),
                        Configuration("1s^2 2s^2 2p^5 3s 3p^5"), Configuration("1s^2 2s^2 2p^5 3p^6")],
                        asfSettings=asfSettings); output=true)[ResultKeys.Multiplet]
@printf("\n>>> core multiplet: %d levels\n", length(cMultiplet.levels))

# ONE final level while the cost is being measured; indexPairs=[(1,0)] would take all five.
const FINAL_LEVEL = 1
const TERM_OF_1   = "4S"

function runOne(label::String, nInt::Int64, factor::Float64)
    st = AutoIonization.Settings(AutoIonization.Settings(), maxAugerEnergy=1.0e6, maxKappa=4,
                                 operator=CoulombInteraction(), calcTeAuger=true,
                                 lineSelection=LineSelection(true, indexPairs=[(1,FINAL_LEVEL)]),
                                 gMultiplet=gMultiplet, gCoreMultiplet=cMultiplet,
                                 NoIntermediateEnergies=nInt, intermediateEnergyFactor=factor)
    wa = Atomic.Computation(Atomic.Computation(), name=label, grid=grid, nuclearModel=nm,
                            initialConfigs=[Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")], initialAsfSettings=asfSettings,
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")], finalAsfSettings=asfSettings,
                            processSettings=st)
    t0    = time()
    lines = perform(wa; output=true)[AutoIonization.Settings]
    dt    = time() - t0
    rate  = sum(l.totalRate for l in lines) * 1.0e7

    return( dt, rate )
end

results = []
println("\n" * "="^104)
@printf("  AR TWO-ELECTRON AUGER -- R(eps) STEP 2:  the principal-value integral, final level %d (%s)\n",
        FINAL_LEVEL, TERM_OF_1)
println("="^104)
@printf("  %-34s %6s %8s %12s %14s %10s\n", "run", "nInt", "factor", "time [s]", "rate [1e-7]", "vs. A")
println("  " * "-"^100)

dtA, rA = runOne("A  R(eps) OFF -- the baseline", 0, 3.)
push!(results, (0, 3., dtA, rA))
@printf("  %-34s %6d %8.1f %12.1f %14.4f %10s\n", "A  R(eps) OFF (Rule 19 baseline)", 0, 3., dtA, rA, "1.0000")

for nInt in [3, 5, 9]
    dt, r = runOne("B  R(eps) ON, $nInt points", nInt, 3.)
    push!(results, (nInt, 3., dt, r))
    @printf("  %-34s %6d %8.1f %12.1f %14.4f %10.4f\n", "B  R(eps) ON", nInt, 3., dt, r, r/rA)
end

println("\n  " * "="^100)
println("  WHAT TO READ FROM THIS")
println("  " * "="^100)
println("  A must reproduce the committed delta-piece value for this final level (4S = 130.371e-7 of the")
println("  five-level run of 26-Oct-2026).  If it does not, the two new Settings fields are not inert when")
println("  switched off, and nothing below means anything.")
println()
println("  The B sequence in nInt is the result.  SETTLING means the principal-value subtraction is doing its")
println("  job; a MONOTONE MARCH with non-shrinking increments is the signature of an unregularized pole, which")
println("  is what PhotoDoubleIonization measured on helium before item 33 was closed -- and it would mean the")
println("  assembly here is wrong rather than the physics being interesting.")
@printf("\n  cost: %.1f s with R(eps) off, %.1f s at %d points;  per-point cost ~ %.1f s\n",
        dtA, results[end][3], results[end][1],
        length(results) > 1 ? (results[end][3] - dtA) / max(results[end][1], 1) : 0.)
@printf("  naive size of the FULL five-level run at 9 points: %.0f s = %.1f h\n",
        5 * (dtA + (results[end][3] - dtA) * 9 / max(results[end][1], 1)),
        5 * (dtA + (results[end][3] - dtA) * 9 / max(results[end][1], 1)) / 3600)

println("\n=== TEA R(eps) STEP-2 PROBE FINISHED ===")
