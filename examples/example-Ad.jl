#
println("Ad) Apply & test for the frequency-independent Breit interaction for an internally generated neon multiplet.")

if  true
    # Last successful:  29-Jul-2026: E1(Coulomb) = -6331.332822 Ha; Delta E(Breit) = 5.245391 Ha,
    # Delta E(Gaunt) = 4.880632 Ha, Gaunt/Breit ratio = 0.93 (literature: ~90% for inner shells -- close
    # match). Timing (CI-only): Breit/Coulomb = 1.38x, Gaunt/Coulomb = 1.07x, Gaunt/Breit = 0.78x (Gaunt
    # genuinely cheaper than full Breit, as expected from skipping the retardation term).
    # Branch 1 (small): Cl-like Xe^35+ (1s^2 2s^2 2p^6 3s^2 3p^5, single reference) -- Z=54, matching the
    # actual configuration this file originally used (its own comment mistakenly said "Cl-like Fe^10+";
    # the configs are Ne-core 3p^5, i.e. Cl-like Xe, not Fe -- corrected here). Compares pure Coulomb,
    # full zero-frequency Breit (CoulombBreit(0.)), and Gaunt-only (magnetic part, CoulombGaunt()) at one
    # small system, reporting energies, timings, and the empirical Gaunt/full-Breit ratio -- literature
    # (arxiv.org/abs/2512.03179, GRASP2K studies) reports ~90% for inner/core shells.
    name       = "Cl-like Xe -- Branch 1 (single reference)"
    refConfigs = [Configuration("[Ne] 3s^2 3p^5")]
    nModel     = Nuclear.Model(54.)
    grid       = Radial.Grid(true)

    settingsCoulomb = AsfSettings()
    settingsBreit   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.))
    settingsGaunt   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombGaunt())

    t0 = time()
    mpCoulomb = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
    tCoulomb  = time() - t0
    basis     = mpCoulomb.levels[1].basis
    t1 = time()
    mpCoulombCI = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb; printout=false)
    tCoulombCI  = time() - t1
    println(">> [Branch 1] Coulomb:        SCF+CI = $(round(tCoulomb,digits=3)) s,  CI-only = $(round(tCoulombCI,digits=3)) s,  E1 = $(mpCoulombCI.levels[1].energy)")

    t2 = time()
    mpBreitCI = Hamiltonian.performCI(basis, nModel, grid, settingsBreit; printout=false)
    tBreitCI  = time() - t2
    println(">> [Branch 1] CoulombBreit:   CI-only = $(round(tBreitCI,digits=3)) s,  E1 = $(mpBreitCI.levels[1].energy)")

    t3 = time()
    mpGauntCI = Hamiltonian.performCI(basis, nModel, grid, settingsGaunt; printout=false)
    tGauntCI  = time() - t3
    println(">> [Branch 1] CoulombGaunt:   CI-only = $(round(tGauntCI,digits=3)) s,  E1 = $(mpGauntCI.levels[1].energy)")

    dBreit = mpBreitCI.levels[1].energy - mpCoulombCI.levels[1].energy
    dGaunt = mpGauntCI.levels[1].energy - mpCoulombCI.levels[1].energy
    println(">> [Branch 1] Delta E (Breit) = $dBreit Ha,  Delta E (Gaunt) = $dGaunt Ha,  Gaunt/Breit ratio = $(round(dGaunt/dBreit,digits=3))")
    println(">> [Branch 1] Timing ratio CI-only:  Breit/Coulomb = $(round(tBreitCI/tCoulombCI,digits=2))x,  Gaunt/Coulomb = $(round(tGauntCI/tCoulombCI,digits=2))x,  Gaunt/Breit = $(round(tGauntCI/tBreitCI,digits=2))x")

elseif  false
    # Last successful:  29-Jul-2026: E1(Coulomb) = -6331.333055 Ha, nlev=31; Delta E(Breit) = 5.246673 Ha,
    # Delta E(Gaunt) = 4.881829 Ha, Gaunt/Breit ratio = 0.93 (matches Branch 1 -- expected, correlating
    # shells stay within n=3). Timing (CI-only): Breit/Coulomb = 1.22x, Gaunt/Coulomb = 1.14x, Gaunt/Breit
    # = 0.93x (ratios shrink vs. Branch 1 as the shared angular-coefficient overhead grows with CSF count,
    # diluting Breit/Gaunt's relative extra cost).
    # Branch 2 (medium): same Cl-like Xe ion family, with a correlating layer added (reference + 3s 3p^6 +
    # 3s^2 3p^4 3d, matching Ad.jl's own original 3-config list) -- a larger CSF space for the timing-factor
    # question at bigger scale, same Coulomb/Breit/Gaunt comparison and Gaunt-ratio check as Branch 1.
    name       = "Cl-like Xe -- Branch 2 (+ correlating layer: 3s 3p^6, 3s^2 3p^4 3d)"
    refConfigs = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")]
    nModel     = Nuclear.Model(54.)
    grid       = Radial.Grid(true)

    settingsCoulomb = AsfSettings()
    settingsBreit   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.))
    settingsGaunt   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombGaunt())

    t0 = time()
    mpCoulomb = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
    tCoulomb  = time() - t0
    basis     = mpCoulomb.levels[1].basis
    t1 = time()
    mpCoulombCI = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb; printout=false)
    tCoulombCI  = time() - t1
    println(">> [Branch 2] Coulomb:        SCF+CI = $(round(tCoulomb,digits=3)) s,  CI-only = $(round(tCoulombCI,digits=3)) s,  E1 = $(mpCoulombCI.levels[1].energy),  nlev = $(length(mpCoulombCI.levels))")

    t2 = time()
    mpBreitCI = Hamiltonian.performCI(basis, nModel, grid, settingsBreit; printout=false)
    tBreitCI  = time() - t2
    println(">> [Branch 2] CoulombBreit:   CI-only = $(round(tBreitCI,digits=3)) s,  E1 = $(mpBreitCI.levels[1].energy)")

    t3 = time()
    mpGauntCI = Hamiltonian.performCI(basis, nModel, grid, settingsGaunt; printout=false)
    tGauntCI  = time() - t3
    println(">> [Branch 2] CoulombGaunt:   CI-only = $(round(tGauntCI,digits=3)) s,  E1 = $(mpGauntCI.levels[1].energy)")

    dBreit = mpBreitCI.levels[1].energy - mpCoulombCI.levels[1].energy
    dGaunt = mpGauntCI.levels[1].energy - mpCoulombCI.levels[1].energy
    println(">> [Branch 2] Delta E (Breit) = $dBreit Ha,  Delta E (Gaunt) = $dGaunt Ha,  Gaunt/Breit ratio = $(round(dGaunt/dBreit,digits=3))")
    println(">> [Branch 2] Timing ratio CI-only:  Breit/Coulomb = $(round(tBreitCI/tCoulombCI,digits=2))x,  Gaunt/Coulomb = $(round(tGauntCI/tCoulombCI,digits=2))x,  Gaunt/Breit = $(round(tGauntCI/tBreitCI,digits=2))x")

elseif  false
    # Last successful:  29-Jul-2026: Delta E(Breit) grows monotonically with Z as expected (5.245 -> 14.872
    # -> 31.358 Ha for Z=54,74,92 -- roughly Z^3.3 empirically between these points, consistent with an
    # inner-shell-dominated relativistic correction). Gaunt/Breit ratio drifts DOWN slightly with Z (0.930
    # -> 0.925 -> 0.919) -- matches the literature nuance that the retardation term matters proportionally
    # more for heavier atoms, even though Gaunt still dominates throughout this range.
    # Branch 3 (Z-scaling sanity check): same single-reference Cl-like configuration across three Z values
    # along the isoelectronic sequence (Xe Z=54, W Z=74, U Z=92) -- confirms the Breit/Gaunt correction
    # grows with Z (a physical sanity check needing no external literature reference, just internal
    # monotonicity), and tracks how the Gaunt/Breit ratio itself evolves with Z.
    refConfigs = [Configuration("[Ne] 3s^2 3p^5")]
    grid       = Radial.Grid(true)
    settingsCoulomb = AsfSettings()
    settingsBreit   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.))
    settingsGaunt   = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombGaunt())

    for  Z  in  [54., 74., 92.]
        nModel = Nuclear.Model(Z)
        mpCoulomb   = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
        basis       = mpCoulomb.levels[1].basis
        mpCoulombCI = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb; printout=false)
        mpBreitCI   = Hamiltonian.performCI(basis, nModel, grid, settingsBreit;   printout=false)
        mpGauntCI   = Hamiltonian.performCI(basis, nModel, grid, settingsGaunt;   printout=false)
        dBreit = mpBreitCI.levels[1].energy - mpCoulombCI.levels[1].energy
        dGaunt = mpGauntCI.levels[1].energy - mpCoulombCI.levels[1].energy
        println(">> [Branch 3] Z=$Z:  E(Coulomb) = $(mpCoulombCI.levels[1].energy) Ha,  Delta E(Breit) = $dBreit Ha,  " *
                "Delta E(Gaunt) = $dGaunt Ha,  Gaunt/Breit ratio = $(round(dGaunt/dBreit,digits=3))")
    end

end
