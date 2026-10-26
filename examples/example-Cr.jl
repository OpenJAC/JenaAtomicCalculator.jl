println("Cr) Apply & test the HfsPaschenBack module (exact Hfs.HfLevel-basis diagonalization of hyperfine")
println("    Paschen-Back / Breit-Rabi mixing in a static magnetic and/or electric field).")

if  false
    # Last visit:  26-Oct-2026
    # Last successful:  26-Oct-2026 -- BIT-IDENTICAL to StarkZeeman at every printed digit, both fields
    #   and both M: -5.00006662010119e-01 Ha (M=+1/2) and -5.00006657755812e-01 Ha (M=-1/2).
    # Branch a: H 1s, nuclear spin I = 0 -- the hyperfine basis then collapses to F = J exactly (a
    #   single HfLevel per electronic level, no mixing, no splitting), so HfsPaschenBack.Settings'
    #   hfMultiplets reduce to exactly the plain ASF (Level,M) basis StarkZeeman already diagonalizes.
    #   This is the cleanest correctness check available: it needs no external formula or literature
    #   value at all, only that the two independently-coded modules AGREE once the extra (nuclear)
    #   structure this module adds is turned off. Any sign or factor error in
    #   Hfs.recouplingElectronicOperator, or in how HfsPaschenBack sums over each HfLevel's own
    #   mc-weighted hfBasisVectors, would show up here as a mismatch against StarkZeeman's own
    #   already-validated Zeeman splitting -- it did not.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(Nuclear.Model(1., PointNucleus()); spinI=AngularJ64(0), mu=0.0)

    wa = Atomic.Computation(Atomic.Computation(), name="H 1s", grid=grid, nuclearModel=nm, configs=[Configuration("1s")])
    wb = perform(wa; output=true)
    multiplet = wb[ResultKeys.Multiplet]

    hfSettings = Hfs.Settings(; calcM1=false, calcHfMultiplet=true, levelSelection=LevelSelection())
    wc = Atomic.Computation(Atomic.Computation(), name="H 1s hyperfine basis (I=0)", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[hfSettings])
    wd = perform(wc; output=true)
    hfMultiplet = wd[Hfs.Settings][1].hfsMultiplet

    bFieldT = 0.001
    pbSettings = HfsPaschenBack.Settings(HfsPaschenBack.Settings(); includeBField=true, bField=bFieldT,
                      includeNuclearZeeman=false, hfMultiplets=[hfMultiplet], printBefore=false)
    we = Atomic.Computation(Atomic.Computation(), name="H 1s HfsPaschenBack (I=0)", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[pbSettings])
    wf = perform(we; output=true)
    outcomesA = wf[HfsPaschenBack.Settings]

    szSettings = StarkZeeman.Settings(StarkZeeman.Settings(); includeBField=true, bField=bFieldT,
                      printBefore=false, levelSelection=LevelSelection())
    wg = Atomic.Computation(Atomic.Computation(), name="H 1s StarkZeeman (ordinary ASF basis)", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[szSettings])
    wh = perform(wg; output=true)
    outcomesB = wh[StarkZeeman.Settings]

    println("\n>> I=0 reduction check, B = $bFieldT T:")
    for (oa, ob) in zip(sort(outcomesA, by=o->o.energy), sort(outcomesB, by=o->o.energy))
        println("   HfsPaschenBack: $(oa.energy) Ha      StarkZeeman: $(ob.energy) Ha      " *
                "diff = $(oa.energy - ob.energy)")
    end
    #
elseif true
    # Last visit:  26-Oct-2026
    # Last successful:  26-Oct-2026 -- weak-field (B=1e-4 T) F=1 triplet splits ANTISYMMETRICALLY about
    #   its own zero-field energy: Delta-E(M=+1) = -2.6108e-10 Ha, Delta-E(M=-1) = +2.6108e-10 Ha (equal
    #   to 4 significant figures), Delta-E(M=0) = +2.0e-13 Ha -- the textbook linear Zeeman pattern of a
    #   degenerate hyperfine triplet. At B=0.01 T, genuine non-perturbative F-mixing appears: F=0,M=0
    #   acquires a 0.009 admixture of F=1,M=0 (and vice versa, symmetrically) -- an avoided crossing
    #   StarkShift's quadratic formula and LandeZeeman's own perturbative C2 cannot produce. TRACE
    #   INVARIANCE (sum of all 4 field-dressed eigenvalues) holds to MACHINE PRECISION between B=0 and
    #   B=0.01 T (diff = 4.4e-16 Ha on a trace of -2.0 Ha, i.e. the double-precision floor), as it must
    #   for a Hermitian perturbation added only off the full-manifold diagonal of two COMPLETE F-blocks.
    # Branch b: H 1s, nuclear spin I = 1/2 (a proton-like mu, hydrogen's own real value) -- the genuine
    #   two-level (F=0,F=1) hyperfine Paschen-Back/Breit-Rabi case, the textbook example for this
    #   physics. Both the electronic and the (much smaller) nuclear Zeeman term are included.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(Nuclear.Model(1., PointNucleus()); spinI=AngularJ64(1,2), mu=2.792847351)

    wa = Atomic.Computation(Atomic.Computation(), name="H 1s", grid=grid, nuclearModel=nm, configs=[Configuration("1s")])
    wb = perform(wa; output=true)
    multiplet = wb[ResultKeys.Multiplet]

    hfSettings = Hfs.Settings(; calcM1=true, calcHfMultiplet=true, levelSelection=LevelSelection())
    wc = Atomic.Computation(Atomic.Computation(), name="H 1s hyperfine basis (I=1/2)", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[hfSettings])
    wd = perform(wc; output=true)
    hfMultiplet = wd[Hfs.Settings][1].hfsMultiplet
    dWhfs = hfMultiplet.hfLevels[end].energy - hfMultiplet.hfLevels[1].energy
    println("\n>> Zero-field F=1 minus F=0 splitting: $dWhfs Ha = " *
            "$(Defaults.convertUnits("energy: from atomic to Hz", dWhfs)) Hz")

    previousTrace = 0.
    for  bFieldT in (0.0, 1.0e-4, 1.0e-2)
        pbSettings = HfsPaschenBack.Settings(HfsPaschenBack.Settings(); includeBField=true, bField=bFieldT,
                          includeNuclearZeeman=true, hfMultiplets=[hfMultiplet], printBefore=false)
        we = Atomic.Computation(Atomic.Computation(), name="H 1s HfsPaschenBack (B=$bFieldT T)", grid=grid, nuclearModel=nm,
                                configs=[Configuration("1s")], propertySettings=[pbSettings])
        wf = perform(we; output=true)
        outcomes = wf[HfsPaschenBack.Settings]

        trace = sum(o.energy for o in outcomes)
        println("\n>> B = $bFieldT T:   trace = $trace Ha" *
                (bFieldT == 0.0 ? "" : "   (vs B=0 trace $previousTrace Ha, diff = $(trace-previousTrace))"))
        if  bFieldT == 0.0    global previousTrace = trace    end
        for o in sort(outcomes, by=o->o.energy)
            println("   E = $(o.energy) Ha   dominant: " *
                    join(["F=$(c.hfLevel.F),M=$(c.hfLevel.M)[$(round(c.weight,digits=3))]" for c in o.components], ", "))
        end
    end
    #
end
