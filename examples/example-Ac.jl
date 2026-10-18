#
println("Ac) Apply & test the CI part for open-shell multiplets of increasing complexity, without Breit interaction.")

if  false
    # Last visit:  18-Oct-2026
    # Last reproduced:  18-Oct-2026 (1 of 2 to 1e-5, 1 more to 1e-3)
    # Last successful:  31-Aug-2026: E1 = -54.414213603268 Ha, nlev=5.
    #   THE NUMBER MOVED BY 1.43 mHa FROM THE 28-Jul-2026 VALUE OF -54.412785002103, AND THAT IS NOT A
    #   REGRESSION. Bisected 31-Aug over the 469 commits since: the shift is `4cc94eb` (15-Aug-2026), "The
    #   kink-aware Slater integral becomes the standard" -- XL_Coulomb integrating the Slater kernel across the
    #   r_< / r_> cusp instead of through it. That commit validated the new rule against the ANALYTIC
    #   F^0(1s,1s) = 5Z/8 (7.85e-5 against the old 2.59e-4 on this very grid) and re-approved twelve reference
    #   files. So the old number was the less accurate one; this branch's comment simply was never updated.
    #   A SECOND, MUCH SMALLER STEP of 9.7e-8 Ha sits at `3876ab0` (31-Jul), the kappa-sign B-spline
    #   boundary-condition fix -- also deliberate, and far too small to be what was noticed.
    # Earlier timing note, kept: CI-only (isolated) 0.194 s BEFORE the Hamiltonian.jl caching fix (see
    #   module-Hamiltonian.jl/module-InteractionStrength.jl, 28-Jul-2026) -> 0.142 s AFTER (1.37x, bit-identical
    #   energy at the time) -- the smallest branch (3 blocks), so the smallest speedup of the four.
    # Branch 1 (smallest): N I, single reference 1s^2 2s^2 2p^3 -- the richest SINGLE p-shell case (2p^2 and
    # 2p^4 are identical in CSF count/structure by particle-hole symmetry, both only 5 CSF/2 blocks; 2p^3
    # alone already reaches 5 CSF across 3 symmetry blocks). Pure Coulomb (AsfSettings() defaults already
    # give eeInteractionCI=CoulombInteraction(), qedModel=NoneQed()), standard DFSField (not Claude2/EOL).
    name        = "N I 1s^2 2s^2 2p^3 -- Branch 1 (single reference)"
    refConfigs  = [Configuration("[He] 2s^2 2p^3")]
    nModel      = Nuclear.Model(7.)
    grid        = Radial.Grid(true)
    settings    = AsfSettings()

    t0 = time()
    mp = SelfConsistent.performSCF(refConfigs, nModel, grid, settings; printout=true)
    println(">> [Branch 1] SCF+CI (bundled) took $(round(time()-t0, digits=3)) s,  lowest level energy = $(mp.levels[1].energy),  nlev = $(length(mp.levels))")

    t1  = time()
    mp2 = Hamiltonian.performCI(mp.levels[1].basis, nModel, grid, settings; printout=false)
    println(">> [Branch 1] CI-only (isolated re-diagonalization) took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mp2.levels[1].energy)")

elseif  false
    # Last visit:  18-Oct-2026
    # Last successful:  unknown ... THE 28-Jul-2026 CLAIM IS WITHDRAWN, 18-Oct-2026. Its E1 = -54.406266419642 Ha
    #   was TOO BOUND, and not by a grid or a convergence margin: `SpinAngular.twoParticleMoveOne` attached the
    #   direct and crossed quadruple labels the wrong way round whenever the spectator subshell lies BETWEEN the
    #   donor and the acceptor, so weight landed on ranks the conversion to R^k annihilates and the rank that
    #   survives was left at zero. Fixed in `1d05ab3b` (13-Sep-2026), validated against GRASP2018 to nine
    #   figures. This branch's 2p^2 3s / 2p^2 3p references ARE that geometry, which is why it moved while
    #   Branch 1, a single 2p^3 reference, is bit-identical to 1e-10.
    #   MEASURED ON A CLEAN HEAD, 18-Oct-2026: E1 = -54.399823189844 Ha, nlev = 34 (unchanged). The whole chain,
    #   each point re-run from an archived code state: 28-Jul -54.406266419640 (reproduced EXACTLY), -1.51 mHa at
    #   `4cc94eb` (15-Aug, the kink-aware Slater integral), then +7.956 mHa at the label fix. Nothing since has
    #   moved it. RE-DATING IS A HUMAN'S PHYSICAL JUDGEMENT (Rule 7); the above is a measurement, not that.
    # Earlier timing note, kept: CI-only (isolated) 0.685 s BEFORE the caching fix -> 0.401 s AFTER (1.71x
    #   speedup, bit-identical energy at the time).
    # Branch 2: N I + one correlating layer -- reference 2p^3 plus single/double excitations represented
    # here as explicit correlating reference configurations (2p^2 3s, 2p^2 3p). Genuine "increasing
    # complexity" needs growing the CI/correlation space, not just relabeling a single p-shell (2p^n alone
    # tops out at 5 CSF regardless of n) -- this is what actually stresses the O(n^2) CSF-pair loop in
    # Hamiltonian.setupMatrix and multi-block cache reuse.
    name        = "N I -- Branch 2 (+ one correlating layer: 2p^2 3s, 2p^2 3p)"
    refConfigs  = [Configuration("[He] 2s^2 2p^3"), Configuration("[He] 2s^2 2p^2 3s"), Configuration("[He] 2s^2 2p^2 3p")]
    nModel      = Nuclear.Model(7.)
    grid        = Radial.Grid(true)
    settings    = AsfSettings()

    t0 = time()
    mp = SelfConsistent.performSCF(refConfigs, nModel, grid, settings; printout=true)
    println(">> [Branch 2] SCF+CI (bundled) took $(round(time()-t0, digits=3)) s,  lowest level energy = $(mp.levels[1].energy),  nlev = $(length(mp.levels))")

    t1  = time()
    mp2 = Hamiltonian.performCI(mp.levels[1].basis, nModel, grid, settings; printout=false)
    println(">> [Branch 2] CI-only (isolated re-diagonalization) took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mp2.levels[1].energy)")

elseif  false
    # Last visit:  18-Oct-2026
    # Last successful:  unknown ... THE 28-Jul-2026 CLAIM IS WITHDRAWN, 18-Oct-2026, for the same reason as
    #   Branch 2 and from the same commit: its E1 = -54.400809719941 Ha was too bound by the interchanged
    #   two-particle quadruple labels of `SpinAngular.twoParticleMoveOne`, fixed in `1d05ab3b` (13-Sep-2026).
    #   MEASURED ON A CLEAN HEAD, 18-Oct-2026: E1 = -54.394094072958 Ha, nlev = 62 (unchanged). Chain: 28-Jul
    #   -54.400809721682 (reproduced to 1.7e-09), -1.52 mHa at `4cc94eb`, then +8.233 mHa at the label fix --
    #   slightly MORE than Branch 2's +7.956, as expected, since 2p^2 3d adds further spectator-between pairs.
    #   RE-DATING IS A HUMAN'S PHYSICAL JUDGEMENT (Rule 7).
    # Earlier timing note, kept: CI-only (isolated) 1.128 s BEFORE the caching fix -> 0.586 s AFTER (1.93x
    #   speedup, energy matching to 9 significant figures -- the ~2e-9 Ha residual is floating-point
    #   summation-order noise from the cache reuse, not a physics change).
    # Branch 3: N I + two correlating layers -- Branch 2 plus 2p^2 3d, growing the CI space further.
    # Branches 1-4 are independent, self-contained computations, NOT a nested convergence series -- energies
    # are NOT expected to lower monotonically branch to branch, since DFS jointly re-optimizes orbitals
    # against each branch's own full config list rather than freezing+growing like EOL/RAS.
    name        = "N I -- Branch 3 (+ two correlating layers: adds 2p^2 3d)"
    refConfigs  = [Configuration("[He] 2s^2 2p^3"), Configuration("[He] 2s^2 2p^2 3s"), Configuration("[He] 2s^2 2p^2 3p"),
                   Configuration("[He] 2s^2 2p^2 3d")]
    nModel      = Nuclear.Model(7.)
    grid        = Radial.Grid(true)
    settings    = AsfSettings()

    t0 = time()
    mp = SelfConsistent.performSCF(refConfigs, nModel, grid, settings; printout=true)
    println(">> [Branch 3] SCF+CI (bundled) took $(round(time()-t0, digits=3)) s,  lowest level energy = $(mp.levels[1].energy),  nlev = $(length(mp.levels))")

    t1  = time()
    mp2 = Hamiltonian.performCI(mp.levels[1].basis, nModel, grid, settings; printout=false)
    println(">> [Branch 3] CI-only (isolated re-diagonalization) took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mp2.levels[1].energy)")

elseif  true
    # Last visit:  18-Oct-2026
    # Last successful:  unknown ... THE 28-Jul-2026 CLAIM IS WITHDRAWN, 18-Oct-2026, and THIS IS THE BRANCH THAT
    #   SHOWS THE SIZE OF THE DEFECT. Its E1 = -74.960229746756 Ha was too bound by 142.8 mHa -- about
    #   31000 cm^-1 -- through the interchanged two-particle quadruple labels of
    #   `SpinAngular.twoParticleMoveOne`, fixed in `1d05ab3b` (13-Sep-2026) and validated against GRASP2018.
    #   MEASURED ON A CLEAN HEAD, 18-Oct-2026: E1 = -74.820550955084 Ha, nlev = 139 (unchanged). Chain: 28-Jul
    #   -74.960229746759 (reproduced EXACTLY), -3.04 mHa at `4cc94eb`, then +142.836 mHa at the label fix, which
    #   alone accounts for the whole move -- the fix commit and HEAD agree to every digit.
    #   WHY IT IS THE LARGEST OF THE FOUR: six reference configurations, 2p^3(3s,3p,3d) and 2s 2p^4(3s,3p), so
    #   far more CSF pairs sit in the spectator-between geometry than in Branches 2-3. The defect was INVISIBLE
    #   to any single-reference test and left the CSF count untouched, which is why nlev = 139 both sides.
    #   RE-DATING IS A HUMAN'S PHYSICAL JUDGEMENT (Rule 7).
    # Earlier timing note, kept: CI-only (isolated) 2.017 s BEFORE the caching fix -> 0.939 s AFTER (2.15x
    #   speedup -- the LARGEST of the four branches, confirming the predicted pattern that speedup grows with
    #   block count: 1.37x (3 blocks) -> 1.71x (7 blocks) -> 1.93x (9 blocks) -> 2.15x (11 blocks)).
    # Branch 4 (largest): O I, a richer multi-reference case -- reference 2p^4 plus correlating 2p^3(3s,3p,3d)
    # and 2s2p^4(3s,3p). "Open-shell SYSTEMS" (plural) supports varying the physical system, not only
    # resizing one fixed system's CI space -- this also exercises a genuinely multi-reference CI setup
    # (6 reference configurations at once) reaching well beyond Branches 1-3's single-open-shell scope.
    name        = "O I -- Branch 4 (richer multi-reference: 2p^4, 2p^3(3s,3p,3d), 2s2p^4(3s,3p))"
    refConfigs  = [Configuration("[He] 2s^2 2p^4"), Configuration("[He] 2s^2 2p^3 3s"), Configuration("[He] 2s^2 2p^3 3p"),
                   Configuration("[He] 2s^2 2p^3 3d"), Configuration("[He] 2s 2p^4 3s"), Configuration("[He] 2s 2p^4 3p")]
    nModel      = Nuclear.Model(8.)
    grid        = Radial.Grid(true)
    settings    = AsfSettings()

    t0 = time()
    mp = SelfConsistent.performSCF(refConfigs, nModel, grid, settings; printout=true)
    println(">> [Branch 4] SCF+CI (bundled) took $(round(time()-t0, digits=3)) s,  lowest level energy = $(mp.levels[1].energy),  nlev = $(length(mp.levels))")

    t1  = time()
    mp2 = Hamiltonian.performCI(mp.levels[1].basis, nModel, grid, settings; printout=false)
    println(">> [Branch 4] CI-only (isolated re-diagonalization) took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mp2.levels[1].energy)")

elseif  false
    # Last successful:  29-Jul-2026, RE-VERIFIED after the ALField promotion/rename.
    # NEW (promoted, bVector-native) implementation: DFS-Field converges fully in 18.414 s,
    # E1 = -54.412785002103 Ha (matches Branch 1 exactly, as it must -- same system). AL-Field, capped at
    # the SAME maxIterationsScf=5 as before, now takes only 41.218 s (~11.7x FASTER than the old
    # implementation's 481.766 s for the identical 5-iteration cap) and reaches E1 = -54.425890 Ha -- only
    # 0.013 Ha (0.02%) away from the converged DFS answer after just 5 iterations, vs. the OLD implementation
    # still being 3.69 Ha away at the same iteration count. Both a large speedup AND a qualitatively better
    # per-iteration convergence rate -- direct, concrete confirmation that the ALField promotion fixed both
    # the correctness bug AND the extreme per-iteration cost documented below (old numbers, for the historical
    # record): DFS 16.639 s / E1=-54.412785 Ha; OLD buggy AL-Field (same 5-iteration cap) 481.766 s
    # (~96 s/iteration) / E1=-50.723 Ha, still far from converged.
    # SEEN AGAIN 12-Oct-2026 (/testExamples): DFS-Field now gives E1 = -54.414213603268 Ha, 1.43 mHa below
    #   the number quoted above -- NOT A REGRESSION, the same `4cc94eb` (15-Aug) kink-aware-Slater-integral
    #   shift that Branch 1's own comment already explains; this branch's text was simply never updated to
    #   match. AL-Field (still capped at 5 iterations, still correctly reported NOT CONVERGED) now lands at
    #   E1 = -54.430615937669245 Ha, moved by the same cause plus its own iteration-count sensitivity.
    # Branch 5: DFS-Field vs AL-Field comparison, reusing Branch 1's small (5 CSF) N I system.
    name        = "N I -- Branch 5 (DFS-Field vs AL-Field, Branch 1's system)"
    refConfigs  = [Configuration("[He] 2s^2 2p^3")]
    nModel      = Nuclear.Model(7.)
    grid        = Radial.Grid(true)

    t0 = time()
    settingsDFS = AsfSettings()
    mpDFS       = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsDFS; printout=true)
    println(">> [Branch 5, DFS-Field] took $(round(time()-t0, digits=3)) s,  lowest level energy = $(mpDFS.levels[1].energy)")

    t1 = time()
    settingsAL = AsfSettings(AsfSettings(); scField=Basics.ALField(), scfRoute=Basics.AverageLevelRoute(5))
    mpAL       = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsAL; printout=true)
    println(">> [Branch 5, AL-Field, maxIterationsScf=5] took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mpAL.levels[1].energy)")

end
