#
println("Ac) Apply & test the CI part for open-shell multiplets of increasing complexity, without Breit interaction.")

if  false
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
    # Last successful:  28-Jul-2026: E1 = -54.406266419642 Ha, nlev=34. CI-only (isolated) timing:
    # 0.685 s BEFORE the caching fix -> 0.401 s AFTER (1.71x speedup, bit-identical energy).
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
    # Last successful:  28-Jul-2026: E1 = -54.400809719941 Ha, nlev=62. CI-only (isolated) timing:
    # 1.128 s BEFORE the caching fix -> 0.586 s AFTER (1.93x speedup, energy matches to 9 significant
    # figures -- the ~2e-9 Ha residual is floating-point summation-order noise from the cache reuse, not a
    # physics change).
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
    # Last successful:  28-Jul-2026: E1 = -74.960229746756 Ha, nlev=139. CI-only (isolated) timing:
    # 2.017 s BEFORE the caching fix -> 0.939 s AFTER (2.15x speedup -- the LARGEST of the four branches,
    # confirming the predicted pattern that speedup grows with block count: 1.37x (3 blocks) -> 1.71x
    # (7 blocks) -> 1.93x (9 blocks) -> 2.15x (11 blocks). Energy matches to 9 significant figures (same
    # floating-point summation-order noise as Branch 3).
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
    settingsAL = AsfSettings(AsfSettings(); scField=Basics.ALField(), maxIterationsScf=5)
    mpAL       = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsAL; printout=true)
    println(">> [Branch 5, AL-Field, maxIterationsScf=5] took $(round(time()-t1, digits=3)) s,  lowest level energy = $(mpAL.levels[1].energy)")

end
