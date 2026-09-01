#
println("Ai) Apply & test for restricted-active-space (RAS) expansions.")

if  true
    # Last visit:  01-Sep-2026
    # Last successful:  01-Sep-2026 -- step1 -14.57132546, step2 -14.61403156 Ha, w_ref 0.90451, 2p at
    #                   2.5276 / 2.5289 a.u., 22.0 s.  See the 01-Sep survey below.
    #
    # ===== SURVEY, 01-Sep-2026 (priority items 7 and 10): ALL THREE BRANCHES RE-RUN AND DATED. =====
    #   branch 1 (2-layer EOL) : -14.57132546 / -14.61403156 Ha            22.0 s
    #   branch 2 (3-layer EOL) : -14.57132546 / -14.61403156 / -14.61475586 Ha   16.7 s
    #   branch 3 (3-layer DFS) : -14.570964036 / -14.605737784 / -14.616949323 Ha  0.6 s
    #
    # TWO THINGS CHANGED SINCE THE 27-Aug SURVEY, and they pull in opposite directions.
    #  (i)  THE ENERGY ROSE, and that is correct.  27-Aug recorded step1 -14.575891726 and step2
    #       -14.618826816; step2 is now 4.8 mHa HIGHER.  Its own note says why: "until then no RAS layer
    #       froze anything ... BOTH layers here re-optimized every orbital".  Freezing is now real in the AL
    #       pass as well as the rotation (01-Sep), so this is a CONSTRAINED minimum and must lie above the
    #       unconstrained one.  It sits within 3e-05 Ha of the converged -14.614058864 that
    #       TestFrames.testRepresentation_RasExpansion asserts, the residual being this branch's 614 a.u. box
    #       and 24-iteration budget.
    #  (ii) BRANCH 2 NOW CONVERGES IN SECONDS.  Its own comment records that the third layer "did NOT converge
    #       within 40 minutes ... energy still dropping steadily after 19 iterations, no plateau".  It takes
    #       16.7 s.  That is the largest single before/after in this survey, and it is the four EOL
    #       convergence fixes of 29-Aug to 01-Sep rather than anything about this case.
    #
    # WHAT THE THIRD LAYER IS WORTH, AND WHY -- read the radii, which each step now prints (item 25).
    # The n = 3 layer earns only -0.00072 Ha, and its orbitals say why: 3s 6.9036, 3p 8.1255 / 8.4324,
    # 3d 7.0131 / 7.1395 a.u., against a 2s at 2.6398 and a 2p at 2.5276.  The layer sits THREE TIMES further
    # out than the electrons it is meant to correlate.  Compare the high-Z branch (d) below, where the layer
    # that contracts ONTO the valence region earns sixty times more than the one that does not.
    #
    # AND A DISCREPANCY WORTH KEEPING: branch 3, the hand-written DFS loop, reaches -14.616949 -- 2.2 mHa
    # BELOW the EOL branch 2 on the same three layers.  Part of it is that branch 3 is a COPY of the driver
    # loop and still seeds its start orbitals from the bare `nuclearPot`, where the driver itself now uses the
    # screened mean field (priority item 23, fixed 01-Sep); the rest is not yet explained.  A field that
    # optimizes a configuration average should not beat one that optimizes the very level being reported.
    # 27-Jul-2026 (the last verified run): -14.57047 (ref) -> -14.57133 (step1, EOL reproduces reference) ->
    #                   -14.58990 (step2, 2p correlation) Hartree; monotonic lowering, ~0.017 Ha 2p-correlation
    #                   contribution called the right order of magnitude for Be's dominant 2s^2<->2p^2 channel.
    # SURVEY, 27-Aug-2026 (priority item 50).  Re-run on current code this branch gives step1 -14.575891726 and
    # step2 -14.618826816, i.e. 4.6 mHa and 28.9 mHa BELOW the recorded values, and the 2p-correlation
    # contribution is now 0.043 Ha rather than 0.017.  The direction is right -- every change since July was a
    # convergence fix, and a variational energy that falls has found more correlation, not less -- but the
    # branch is NOT re-dated, for two reasons.  Step 2 still ends on "STOPPED after 24 iterations" with
    # |grad| = 1.7e-03, so it is an upper bound and not a converged number; and it runs on Radial.Grid(true),
    # 614 a.u. for a four-electron atom, the very grid on which testRepresentation_RasExpansion records that a
    # correlation layer can RAISE the energy.  Re-dating wants a matched box and a larger iteration budget
    # first.  The same physical case on the reference-sized box, at 24/60/120 iterations, is converged to
    # 6.7e-09 Ha -- see the value recorded in TestFrames.testRepresentation_RasExpansion.
    #   Part of the move is the frozen-orbital fix of the same day: until then no RAS layer froze anything,
    # despite the description below, so BOTH layers here re-optimized every orbital.
    # Scenario A: 2-layer Be RAS. Reference 1s^2 2s^2, core 1s always frozen. Layer 1 = reference SCF only
    # (no correlation shells, no excitations -- just re-optimizes 1s/2s on the single reference CSF). Layer 2
    # adds 2p as a single+double-excitation correlating shell; 1s and 2s stay frozen (already optimized in
    # layer 1), only 2p is variationally optimized (EOL, target = lowest ^1S_0 level).
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state -- 2-layer RAS (reference, then 2p correlation)"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1]) )
    coreShells  = [Shell("1s")]
    fromShells  = [Shell("2s")]
    layers      = [ RasLayer(Shell[]; se=false, de=false),
                     RasLayer([Shell("2p")]) ]

    wa          = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs,
                                 RasExpansion([LevelSymmetry(0, Basics.plus)], 4, coreShells, fromShells, layers, rasSettings) )
    println("wa = $wa")

    wb = generate(wa, output=true)
    #
elseif  false
    # Last visit:  01-Sep-2026
    # Last successful:  01-Sep-2026 -- -14.57132546 / -14.61403156 / -14.61475586 Ha in 16.7 s, w_ref
    #                   0.90211 at step 3, 3l orbitals at 6.90-8.43 a.u.  THE OPEN ISSUE BELOW IS CLOSED:
    #                   what took more than 40 minutes without converging now takes seventeen seconds.
    # OPEN ISSUE (28-Jul-2026, root cause CONFIRMED, fix deferred; RESOLVED 01-Sep-2026 by the four EOL
    # convergence fixes): this 3rd layer (3s,3p,3d, 11 CSF) did
    # NOT converge within 40 minutes even at a single EOL target level -- energy still dropping steadily
    # after 19 iterations, no plateau. This is NOT primarily a performance problem: it is the same
    # winner-take-all instability confirmed for Scenario A's own step 2 (see the KNOWN LIMITATION note on
    # SelfConsistent.solveOptimizedLevelField and project_eol_implementation.md) -- competing near-degenerate
    # correlating CSFs (here across 3s/3p/3d as well as 2p) fight for the same correlation channel, and the
    # (1.0/occ) Fock-matrix scaling diverges as a losing CSF's mixing coefficient shrinks, rather than
    # settling. More competing shells (this layer) means more such instabilities at once, hence even slower/
    # less likely to plateau than Scenario A's single 2p_1/2-vs-2p_3/2 case. Needs the real
    # DA/inhomogeneous-term mechanism (deferred) before this scenario can be trusted or dated. Steps 1-2
    # below are shared with (and reproduce exactly) the working, dated Scenario A above.
    # Scenario B: 3-layer Be RAS, extending Scenario A. Layer 1 = reference SCF only; layer 2 adds 2p;
    # layer 3 adds 3s, 3p, 3d as a further correlation layer (1s, 2s, 2p frozen throughout layer 3). This is
    # the "2-3 layers" scope the RAS scheme is meant for -- a small, transparent convergence-vs-layer-count
    # test case (each layer's lowest-level energy should decrease monotonically as correlation is added).
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state -- 3-layer RAS (reference, 2p, then 3s3p3d)"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1]) )
    coreShells  = [Shell("1s")]
    fromShells  = [Shell("2s")]
    layers      = [ RasLayer(Shell[]; se=false, de=false),
                     RasLayer([Shell("2p")]),
                     RasLayer([Shell("3s"), Shell("3p"), Shell("3d")]) ]

    wa          = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs,
                                 RasExpansion([LevelSymmetry(0, Basics.plus)], 4, coreShells, fromShells, layers, rasSettings) )
    println("wa = $wa")

    wb = generate(wa, output=true)
    #
elseif  false
    # Last visit:  01-Sep-2026
    # Last successful:  01-Sep-2026 -- -14.570964036 / -14.605737784 / -14.616949323 Ha, total 0.6 s.
    #                   FASTER THAN THE EOL ROUTE BY A FACTOR OF THIRTY, as this branch expects -- and, less
    #                   expectedly, 2.2 mHa LOWER than it on the same three layers.  A field that optimizes a
    #                   configuration AVERAGE should not beat one that optimizes the very level being
    #                   reported, so this is a real discrepancy and not a win.  Part of it: this branch is a
    #                   COPY of the driver loop and still seeds its start orbitals from the bare `nuclearPot`,
    #                   whereas the driver itself now uses the screened mean field (item 23, fixed 01-Sep).
    #                   Bring the two into line before drawing any conclusion from the comparison.
    # Scenario C: DFS-Field analog of Scenario B (SAME 3-layer structure: reference, then 2p, then
    # 3s3p3d), for a rough timing/accuracy comparison against the EOL mode above. DFS-Field is a much
    # cheaper local (Slater-exchange) mean field -- a single potential build per layer, no per-CSF-pair
    # Fock matrix, no outer CI+SCF loop -- so it should be far faster, at the cost of not optimizing
    # orbitals against the CI-mixed target level(s) the way EOL does. Hand-rolled here (does not call
    # Basics.generate(::RasExpansion,...), which hardcodes scField=EOLField()) since this is a one-off
    # diagnostic comparison, not a permanent second mode.
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state -- 3-layer RAS, DFS-Field (timing/accuracy reference)"
    refConfigs  = [Configuration("[He] 2s^2")]
    nModel      = Nuclear.Model(4.)
    grid        = Radial.Grid(true)
    coreShells  = [Shell("1s")]
    fromShells  = [Shell("2s")]
    layers      = [ RasLayer(Shell[]; se=false, de=false),
                     RasLayer([Shell("2p")]),
                     RasLayer([Shell("3s"), Shell("3p"), Shell("3d")]) ]
    rasExpansion = RasExpansion([LevelSymmetry(0, Basics.plus)], 4, coreShells, fromShells, layers, RasSettings())
    wa          = Representation(name, nModel, grid, refConfigs, rasExpansion)
    println("wa = $wa")

    t0             = time()
    priorMultiplet = SelfConsistent.performSCF(refConfigs, nModel, grid, AsfSettings(); printout=true)
    nuclearPot     = Nuclear.nuclearPotential(nModel, grid)
    subshellList   = Basics.extractRelativisticSubshellList(wa)
    primitives     = Bsplines.generatePrimitives(grid)
    startOrbitals  = Bsplines.generateOrbitals(subshellList, nuclearPot, nModel, primitives, printout=true)
    println(">> [DFS] reference multiplet took $(round(time()-t0, digits=1)) s,  lowest level energy = $(priorMultiplet.levels[1].energy)")

    for  (istep, step)  in  enumerate(rasExpansion.steps)
        global priorMultiplet
        tstep      = time()
        println("")
        printstyled("++ [DFS] Compute the orbitals and multiplet for step $istep ... \n", color=:light_yellow)
        basis      = Basics.generateBasis(wa.refConfigs, rasExpansion.symmetries, step)
        orbitals   = Basics.generateOrbitalsForBasis(basis, step.frozenShells, priorMultiplet.levels[1].basis, startOrbitals)
        basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )

        frozenSubshellsThisStep = [ sh  for  shell in step.frozenShells  for sh in basis.subshells
                                        if  sh.n == shell.n  &&  Basics.subshell_l(sh) == shell.l ]

        stepSettings = AsfSettings( AsfSettings();  scField=Basics.DFSField(),  frozenSubshells=frozenSubshellsThisStep )

        multiplet  = SelfConsistent.performSCF(basis, nModel, grid, stepSettings; printout=true)
        println(">> [DFS] step $istep took $(round(time()-tstep, digits=1)) s,  lowest level energy = $(multiplet.levels[1].energy)")
        priorMultiplet = multiplet
    end
    println(">> [DFS] total wall time = $(round(time()-t0, digits=1)) s")
    #
elseif  false
    # Last visit:  01-Sep-2026
    # Last successful:  01-Sep-2026 -- numbers below, produced by the run this branch performs.
    #
    # d) A HIGH-Z FOUR-LAYER REFERENCE: Be-like at Z = 92, the 2s^2 (0+) -> 2s2p (1-) resonance line.
    #    Added as a reference case because the low-Z branches above cannot show what a RAS expansion does
    #    when relativity dominates, and because four layers is where this driver currently stops being
    #    trustworthy (see the NOTE ON DEPTH below).
    #
    #    THE MODEL.  Reference = the full n = 2 COMPLEX {2s^2, 2s2p, 2p^2}, not 2s^2 alone: the ground state's
    #    dominant correlation is 2s^2 <-> 2p^2, so a single-configuration reference would leave it out while
    #    the 1- state carries 2p as a spectroscopic shell -- two DIFFERENT correlation models, whose energies
    #    cannot be compared.  A multireference also avoids a zero-occupation crash: a single-configuration 0+
    #    basis still carries 2p in its subshell list with mean occupation 0, and the AL pass rightly refuses it.
    #    ONE EXPANSION PER SYMMETRY, the transition formed afterwards from two absolute totals, each
    #    variational for its own level.
    #
    #    MEASURED, 01-Sep-2026, four layers (reference / +{3s,3p,3d} / +{4s..4f} / +{5s..5g}), EOL target by
    #    reference weight, Coulomb and Coulomb+Breit:
    #
    #      layer     Coulomb [cm^-1]   +Breit [cm^-1]   Breit shift    (2s^2 -> 2s2p)
    #        1          2 523 378.55     2 748 515.60    225 137.06
    #        2          2 499 736.60     2 728 708.61    228 972.01
    #        3          2 499 876.70     2 728 932.92    229 056.22
    #        4          2 500 032.86     2 729 178.36    229 145.50
    #
    #      0+ totals [Ha]  C : -12040.8608630542 / -12040.8626981076 / -12040.8637114179 / -12040.8647392645
    #                      CB: -12025.9649174697 / -12025.9687581045 / -12025.9705520430 / -12025.9723015723
    #      1- totals [Ha]  C : -12029.3635044239 / -12029.4730601271 / -12029.4734351020 / -12029.4737514117
    #                      CB: -12013.4417589374 / -12013.5358468782 / -12013.5366187808 / -12013.5372500104
    #      wall time       0+ 0.11 h / 0.14 h      1- 0.39 h / 0.53 h      (Coulomb / +Breit)
    #
    #    WHAT THE NUMBERS SAY, and why this case is worth keeping.
    #    (i)   CORRELATION IS WORTH -23 346 cm^-1 (0.9 %) AND 99 % OF IT ARRIVES IN THE FIRST CORRELATION
    #          LAYER.  Layer 2 moves the transition by -23 642; layers 3 and 4 move it BACK by +140 and +156.
    #          Every total energy falls monotonically, but the DIFFERENCE overshoots and returns.
    #    (ii)  BREIT AND CORRELATION ARE SEPARABLE -- BUT NOT AT LAYER 1.  The Breit shift changes by 3 835
    #          cm^-1 (1.7 %) from the bare reference to one correlation layer, then by 84 and 89 cm^-1
    #          (0.04 %).  A Breit correction taken from the UNCORRELATED reference is wrong by ~2 %; taken
    #          with a single layer it is good to 0.04 %.
    #    (iii) THE NON-ADDITIVITY IS IN THE CORRELATION ENERGY ITSELF: the same layer buys -0.00184 Ha under
    #          Coulomb and -0.00384 under Breit for 0+ (2.1x MORE), and -0.10956 against -0.09409 for 1-
    #          (14 % LESS).  Breit does not scale correlation uniformly.
    #    (iv)  BREIT COSTS ~30 % MORE WALL-CLOCK AND NOTHING IN THE SCF: the orbital radii come out identical
    #          to four digits with and without it, because the EOL path applies it once at the final CI and
    #          leaves the variational loop pure Coulomb.
    #    (v)   WHICH LAYER PAYS IS VISIBLE IN THE RADII, which each step now prints.  The n = 3 layer contracts
    #          onto the valence region for 1- (3d at 0.0613 a.u., against 2p at 0.0439/0.0532 and 2s at 0.0545)
    #          and earns -0.10956 Ha; for 0+ it sits at twice that radius (3d 0.1105/0.1149) and earns
    #          -0.00184.  SIXTY TIMES the correlation energy, and the radius predicts it before the energy is
    #          looked at.  Successive layers then expand and saturate: 0.05 -> 0.12 -> 0.22 -> 0.27 a.u.,
    #          +80 % from n=3 to n=4 but only +20 % from n=4 to n=5, with every l of the fifth layer between
    #          0.260 and 0.274 -- one shell of a single characteristic size.
    #    (vi)  AT HIGH Z THE REFERENCE WEIGHT IS A BLUNT INSTRUMENT: it stays at 1.00000 through layer 2 here,
    #          where Be at Z = 4 drops to 0.905 on its first layer.  Correlation energy is roughly constant in
    #          Hartree while the total scales like Z^2, so the RELATIVE admixture collapses with Z.  Use the
    #          radii instead.
    #
    #    NOTE ON DEPTH -- WHY FOUR AND NOT FIVE.  A fifth layer runs inside the same time budget and returns a
    #    WRONG number: the energy RISES by 568 Ha (0+) and 535 Ha (1-), which a larger basis cannot do to a
    #    variational minimum.  The fault is the DEPTH and not the layer: the same {6s..6h} layer is clean when
    #    it is step 4 instead of step 5, and the failure survives a change of reference model and does not need
    #    an h shell.  Priority item 32.  Until that is closed, four layers is the honest limit of this driver.
    #
    Z      = 92.0
    refs   = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s^1 2p^1"), Configuration("1s^2 2p^2")]
    layers = [ RasLayer(Shell[]; se=false, de=false),
               RasLayer([Shell("3s"), Shell("3p"), Shell("3d")]),
               RasLayer([Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]),
               RasLayer([Shell("5s"), Shell("5p"), Shell("5d"), Shell("5f"), Shell("5g")]) ]
    # The 0+ symmetry alone, so that the branch stays around ten minutes; swap the symmetry below for the 1-
    # partner and difference the two lowest totals to reproduce the table above.
    rasSettings = RasSettings([1], 60, 1.0e-6, Basics.CoulombBreit(0.0), LevelSelection(true, configurations=refs))
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rnt = 2.0e-7)
    wa     = Representation("Be-like Z=92 -- 4-layer RAS, Coulomb+Breit", Nuclear.Model(Z), grid, refs,
                            RasExpansion([LevelSymmetry(0, Basics.plus)], 4, [Shell("1s")],
                                         [Shell("2s"), Shell("2p")], layers, rasSettings) )
    println("wa = $wa")
    wb = generate(wa, output=true)
    for i = 1:length(layers)
        k = "step" * string(i)
        haskey(wb, k)  &&  println(">> step $i : lowest level = $(sort(wb[k].levels, by = l -> l.energy)[1].energy) Ha")
    end
    #
end
