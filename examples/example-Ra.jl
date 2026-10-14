
println("Ra) Apply & test the ExportDataSet/ExportAdas modules: write a BASIC adf04 file from JAC's own output.")
println("    Demonstrates the pipeline, not a converged physics result -- see this file's own comments.")

setDefaults("print summary: open", "zzz-ExportAdas.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

if  true
    # Last successful:  14-Oct-2026
    #
    # Branch a: He-like C V (Z=6), 1s^2 1S_0 (ground) -> 1s 2p 1P_1 -- THE SAME transition example-Dl.jl's
    #   branch c already exercised for ImpactExcitation.computeEffStrengths. Reused here rather than
    #   re-run: that branch's own 14-min REDUCED calculation gave Upsilon(T) machinery-verified to 0.3%
    #   against the standard Maxwellian formula (NOT accuracy-dated -- see its own "Last visit", not "Last
    #   successful"), and re-running it would cost the same 14 minutes for a plumbing demonstration that
    #   does not need a fresh number. Only the A-value (cheap, bound-state only, no continuum) is computed
    #   fresh here.
    #
    # WHAT THIS BRANCH IS FOR: showing that JAC's own computed levels, A-value and effective collision
    # strengths assemble into a real external file with NO new physics -- not a new, dated physics result.
    # The ionisation potential (489.99 eV, NIST ASD) is quoted, not computed, since it plays no role in
    # whether the level/transition data itself is right.
    ni = Nuclear.Model(6.)
    wa = Atomic.Computation(Atomic.Computation(), name="Ra-a: He-like C V, 1s^2 -> 1s2p A-value", grid=Radial.Grid(true),
                            nuclearModel=ni,
                            configs=[Configuration("1s^2"), Configuration("1s 2p")],
                            propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.) ] )
    wb    = perform(wa; output=true)
    lines = wb[Einstein.Settings]
    ## Two J=1(odd) -> J=0(even, ground) lines exist here: the spin-forbidden intercombination ^3P_1 and
    ## the E1-allowed, strongly dominant ^1P_1 "resonance line" -- picked out by A-value, not list order,
    ## which is the physically meaningful criterion and not sensitive to how Einstein happens to order lines.
    candidates = [ln for ln in lines if ln.initialLevel.J == AngularJ64(1) && ln.finalLevel.J == AngularJ64(0)]
    line       = candidates[ argmax([ln.photonRate.Babushkin for ln in candidates]) ]

    ## Einstein.Line.photonRate is in ATOMIC UNITS regardless of `setDefaults("unit: rate", ...)`, which
    ## only affects later DISPLAY formatting, not the stored field -- the same kind of gotcha the IsotopeShift
    ## King-plot work ran into overnight. Convert explicitly.
    aValue  = Defaults.convertUnits("rate: from atomic to 1/s", line.photonRate.Babushkin)
    deltaE  = Defaults.convertUnits("energy: from atomic to Kayser", line.initialLevel.energy - line.finalLevel.energy)
    println("\n>>> Computed A(1s2p 1P_1 -> 1s^2 1S_0) = $aValue /s,  deltaE = $deltaE cm^-1\n")

    levels = [ ExportDataSet.LevelRecord(1, "1s2", 1, 0, 0.0, 0.0),
               ExportDataSet.LevelRecord(2, "1s2p", 1, 1, 1.0, deltaE) ]

    ## Upsilon(T), T = [1e5, 5e5, 1e6, 2e6, 5e6] K, quoted from example-Dl.jl branch c's own report
    ## (04-Aug-2026): machinery-verified to 0.3% against the standard Maxwellian formula, not accuracy-dated.
    temperatures = [1.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6]
    upsilons     = [0.042018, 0.049975, 0.059768, 0.077571, 0.11769]
    transitions  = [ ExportDataSet.TransitionRecord(2, 1, aValue, upsilons) ]

    ds = ExportDataSet.DataSet("C", 4, 6, 489.99*8065.544, levels, temperatures, transitions)

    println(">>> adf04 output:\n")
    ExportAdas.writeAdf04(stdout, ds)
    #
end
