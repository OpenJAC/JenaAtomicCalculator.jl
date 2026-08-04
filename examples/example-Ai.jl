#
println("Ai) Apply & test for restricted-active-space (RAS) expansions.")

if  true
    # Last successful:  27-Jul-2026: -14.57047 (ref) -> -14.57133 (step1, EOL reproduces reference) ->
    #                   -14.58990 (step2, 2p correlation) Hartree; monotonic lowering, ~0.017 Ha 2p-correlation
    #                   contribution is the right order of magnitude for Be's dominant 2s^2<->2p^2 channel.
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
    # Last successful:  unknown ...
    # OPEN ISSUE (28-Jul-2026, root cause CONFIRMED, fix deferred): this 3rd layer (3s,3p,3d, 11 CSF) did
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
    # Last successful:  unknown ...
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
end
