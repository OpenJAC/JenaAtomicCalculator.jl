#
# tools/probeTeaContinuumIntermediates.jl
#
# Runs the argon two-electron Auger case of examples/example-De.jl branch f TWICE in one process -- once with
# the BOUND intermediate sum alone, exactly as the branch stands, and once with the delta piece of the
# CONTINUUM intermediate sum added coherently -- and reports both.  Written 26-Oct-2026 for priority items 93
# and 98.  It is a diagnostic, not a test.
#
# WHY TWO RUNS IN ONE PROCESS.  The first run is the BASELINE required before any convention-touching change
# (Rule 19): AutoIonization.Settings gained a gCoreMultiplet field, and with that field left empty the branch
# must reproduce its recorded 585.39e-7 a.u. EXACTLY.  A second number is only worth reading once the first
# one is confirmed, and running both in one process puts them on one compilation and one set of orbitals.
#
# WHAT THE SECOND RUN ADDS.  The sum over intermediate states of a second-order amplitude runs over a COMPLETE
# set, which partitions by bound-electron count: 16 bound electrons (what gMultiplet carries) and 15 bound + 1
# free (what gCoreMultiplet carries, with a continuum spectator attached).  The two are disjoint and add
# COHERENTLY.  Zitnik et al., PRA 93, 021401(R) (2016), compute ONLY the continuum part and omit the bound part
# as negligible; JAC until today computed ONLY the bound part.  Neither is the whole amplitude.
#
# THE PREDICTION, WRITTEN BEFORE THE RUN.  tools/probeTeaDenominators.jl measured the two propagators and found
# them to be the same quantity with OPPOSITE SIGN -- the dominant bound family at +252.7 to +256.7 eV against
# the whole continuum family at -234.0 to -246.8 eV, agreeing to about 6 %.  So the two parts should interfere
# DESTRUCTIVELY and run B should come out BELOW run A.  The numerators are a different matter: they carry the
# two operators in opposite order, V^(e-e) then H - E on the bound path and the reverse on the continuum one,
# so the cancellation is not guaranteed to be as complete as the denominators alone suggest.  That is exactly
# what this run measures, and the amount by which it falls short is the useful number.
#
# WHAT IT ESTABLISHED, 26-Oct-2026.  THE INTERFERENCE IS DESTRUCTIVE, and worth about a third of the rate:
#
#     term            A (bound)   B (bound+cont)    B / A
#     4S                180.711          130.371   0.7214
#     2D                266.586          176.857   0.6634
#     2P                138.089           86.004   0.6228
#     TOTAL             585.385          393.232   0.6717     [1e-7 a.u.]
#
#     against Zitnik Table I theory 228e-7:      2.57x  ->  1.72x
#     against the experiment-anchored 251(46):   2.33x  ->  1.57x
#
# Run A reproduced the recorded 585.39e-7 to 1e-5 and term by term, so the gCoreMultiplet field changes nothing
# when it is empty.  In AMPLITUDE the continuum delta piece is sqrt(0.6717) = 0.82, i.e. about 18 % of the bound
# amplitude with opposite sign -- a real contribution and NOT the near-cancellation the propagators alone
# suggested.  The reason is visible in the per-intermediate output of run A: the deep 3p^-2 bound family
# dominates (4.7e-4 against ~6e-5 for 2p^-1 3p^-1 on the p_1/2 wave) despite its LARGER denominator, and that
# family is exactly the one with no delta-piece partner.
#
# AND THE CANCELLATION IS TERM-DEPENDENT, 0.62 to 0.72, which is the new information.  The bound sum alone was
# non-uniform across the terms after the H - E repair (3.06/2.67/2.00 against Zitnik); the continuum piece pulls
# 2P down hardest and 4S least, which flattens that pattern in the right direction without closing it.
#
# ALSO WORTH KNOWING: the junctions come out at 0.005 to 0.027 with MIXED SIGNS from one core to the next, so
# the continuum contribution partly cancels WITHIN ITSELF before it ever meets the bound sum.  Its net size is
# not the sum of its parts, and a truncated core space can therefore err in either direction.
#
# WHAT IS STILL MISSING, and it is the larger half: R(eps), the regular part of the continuum branch, carried as
# a principal-value integral with an on-shell pole at eps = E_i - E_g.  That is what Zitnik et al. actually
# compute, so the 1.57x remaining is NOT a verdict on their number -- the two calculations still do not contain
# the same terms.
#
# CERTIFIED: src/ stood still for the whole of the A/B run (the peer session's last edit to
# module-SelfConsistent-inc-optimizedlevel.jl was 09:10:47, the run began after 09:13).
#
# Run:  julia --project=<JAC-root> tools/probeTeaContinuumIntermediates.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

# [AR-TEA-CONTINUUM]  ONE question: what does the delta piece of the continuum intermediate sum do to the total,
# in sign and in size?  ONE variable between the two runs: gCoreMultiplet, empty against the three one-2p-hole
# configurations of Ar^3+.  Everything else is the branch-f setup verbatim.

# THE CONTINUUM DEFAULTS ARE PART OF THE SETUP, NOT DECORATION.  examples/example-De.jl sets both at the top of
# the file, and a first attempt at this probe omitted them: the run then used Coulomb-function normalization
# instead of pure sine and every rate came out 1.1 % LOW -- uniformly, because all five final levels share the
# one continuum orbital machinery.  That uniformity is the signature of a normalization difference, and it is
# what the baseline check exists to catch.
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
asfSettings = AsfSettings(AsfSettings(), scField=Basics.HSField())
nm          = Nuclear.Model(18.)

# The BOUND intermediate space, as branch f supplies it: 16 electrons.
gComp      = Atomic.Computation(Atomic.Computation(), name="bound Green multiplet", grid=grid, nuclearModel=nm,
                                configs = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^5"),
                                           Configuration("1s^2 2s^2 2p^6 3s^2 3p^4")], asfSettings=asfSettings)
gMultiplet = perform(gComp; output=true)[ResultKeys.Multiplet]

# The IONIC CORES of the continuum intermediate space: 15 electrons, one 2p hole.  These are exactly the final
# configurations of the first-order L23^2 width anchor in example-De.jl branch (g), which is why the first step
# of this path is machinery already measured against experiment (361 meV against 323(15), i.e. 1.12x).
cComp      = Atomic.Computation(Atomic.Computation(), name="continuum-intermediate ionic cores", grid=grid, nuclearModel=nm,
                                configs = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^4"),
                                           Configuration("1s^2 2s^2 2p^5 3s 3p^5"),
                                           Configuration("1s^2 2s^2 2p^5 3p^6")], asfSettings=asfSettings)
cMultiplet = perform(cComp; output=true)[ResultKeys.Multiplet]

# The five 3p^-3 final levels group into the three LS terms by index, as branch f records them.
termOf(idx) = idx == 1 ? "4S" : (idx in [2,3] ? "2D" : "2P")

function runOne(label::String, coreMultiplet::Multiplet)
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), printBefore=true, maxAugerEnergy=1.0e6,
                                            maxKappa=4, operator=CoulombInteraction(),
                                            lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                                            calcTeAuger=true, gMultiplet=gMultiplet, gCoreMultiplet=coreMultiplet)
    wa = Atomic.Computation(Atomic.Computation(), name=label, grid=grid, nuclearModel=nm,
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")], initialAsfSettings=asfSettings,
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")], finalAsfSettings=asfSettings,
                            processSettings = augerSettings)
    return( perform(wa; output=true)[AutoIonization.Settings] )
end

linesA = runOne("RUN A: bound intermediates only (the baseline)",         Multiplet())
linesB = runOne("RUN B: bound + continuum delta piece, added coherently", cMultiplet)

function report(label, lines)
    println("\n  " * label)
    println("  " * "-"^104)
    @printf("  %-6s %-9s %14s %18s %12s\n", "f", "J^P", "eps_f [eV]", "rate [1e-7 a.u.]", "term")
    terms = Dict("4S" => 0., "2D" => 0., "2P" => 0.);   total = 0.
    for  line in lines
        idx = line.finalLevel.index;    r = line.totalRate * 1.0e7
        @printf("  %-6d %-9s %14.3f %18.3f %12s\n", idx, string(LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)),
                Defaults.convertUnits("energy: from atomic to eV", line.electronEnergy), r, termOf(idx))
        terms[termOf(idx)] += r;    total += r
    end
    @printf("\n  term sums [1e-7 a.u.]:   4S %10.3f    2D %10.3f    2P %10.3f\n", terms["4S"], terms["2D"], terms["2P"])
    @printf("  TOTAL [1e-7 a.u.]:    %12.3f\n", total)

    return( total, terms )
end

println("\n" * "="^112)
totA, trmA = report("RUN A -- bound intermediates only;  MUST reproduce the recorded 585.39e-7 (4S 180.7, 2D 266.6, 2P 138.1)", linesA)
totB, trmB = report("RUN B -- bound + continuum delta piece", linesB)

println("\n" * "="^112)
@printf("  BASELINE CHECK:   run A total %12.3fe-7 against the recorded 585.390e-7   ->  %s\n", totA,
        abs(totA - 585.39)/585.39 < 1.0e-4 ? "REPRODUCES" : "DOES NOT REPRODUCE -- stop and read before anything else")
println("  " * "-"^108)
@printf("  %-30s %14s %14s %12s\n", "quantity", "A (bound)", "B (bound+cont)", "B / A")
for  t in ["4S", "2D", "2P"]
    @printf("  %-30s %14.3f %14.3f %12.4f\n", "term " * t, trmA[t], trmB[t], trmA[t] == 0. ? NaN : trmB[t]/trmA[t])
end
@printf("  %-30s %14.3f %14.3f %12.4f\n", "TOTAL", totA, totB, totA == 0. ? NaN : totB/totA)
println("  " * "-"^108)
@printf("  against Zitnik Table I theory 228e-7:      A %8.2fx    B %8.2fx\n", totA/228., totB/228.)
@printf("  against the experiment-anchored 251(46):   A %8.2fx    B %8.2fx\n", totA/251., totB/251.)
println("\n  B < A would confirm the destructive interference the denominators predict; B > A would refute it.")
println("\n=== TEA CONTINUUM-INTERMEDIATE PROBE FINISHED ===")
