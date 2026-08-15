#
println("Ad) Apply & test for the Breit interaction, frequency-independent and frequency-dependent.")

if  true
    # Last successful:  14-Aug-2026: E1(Coulomb) = -6331.33287512 Ha; Delta E(Breit) = 5.24539117 Ha,
    # Delta E(Gaunt) = 4.88063173 Ha, Gaunt/Breit ratio = 0.930 (literature: ~90% for inner shells).
    #   THIS BRANCH CAUGHT A REAL BUG, and is the reason to keep running it. Between 10-Aug and 14-Aug
    # it read Delta E(Breit) = 1.82701667 Ha and a ratio of 0.80. Bisection put the change at 8f0930b,
    # which gave AngularMomentum.CL_reduced_me the parity rule it had genuinely been missing -- correct
    # in itself, but it silently zeroed the nu = L block of XL_Breit_coefficients, whose guard demands
    # l_a+l_c+L ODD while the shared prefactor now required EVEN. That block carries the dominant
    # magnetic term, so Gaunt came out a factor 3.3 too small; the retardation part, which wants even
    # parity, was untouched. Fixed by giving the nu = L block its own -kappa prefactor. The whole JAC
    # test suite stayed at 45/45 across both the breakage and the repair, so this file, not the suite,
    # is what covers the Breit interaction.
    #   The 29-Jul-2026 reading is reproduced to 8 significant figures (5.245391 / 4.880632 / 0.93); the
    # residual 3e-8 tracks the 2e-8 that E1 itself has moved since, from unrelated SCF changes. Timing
    # then (CI-only): Breit/Coulomb = 1.38x, Gaunt/Coulomb = 1.07x, Gaunt/Breit = 0.78x (Gaunt genuinely
    # cheaper than full Breit, as expected from skipping the retardation term).
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
    # Last successful:  14-Aug-2026: E1(Coulomb) = -6331.333108 Ha, nlev=31; Delta E(Breit) = 5.246673 Ha,
    # Delta E(Gaunt) = 4.881829 Ha, Gaunt/Breit ratio = 0.930 (matches Branch 1 -- expected, correlating
    # shells stay within n=3). Both Breit numbers reproduce the 29-Jul-2026 reading to every digit recorded
    # then, which is how the CL_reduced_me parity repair of branch 1 was confirmed to be complete rather
    # than merely to have restored the single case it was diagnosed on. Timing (CI-only) on 29-Jul:
    # Breit/Coulomb = 1.22x, Gaunt/Coulomb = 1.14x, Gaunt/Breit = 0.93x (ratios shrink vs. Branch 1 as the
    # shared angular-coefficient overhead grows with CSF count, diluting Breit/Gaunt's relative extra cost).
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
    # Last successful:  14-Aug-2026: Delta E(Breit) grows monotonically with Z as expected (5.245391 ->
    # 14.872162 -> 31.357648 Ha for Z=54,74,92 -- roughly Z^3.3 empirically between these points, consistent
    # with an inner-shell-dominated relativistic correction). Gaunt/Breit ratio drifts DOWN slightly with Z
    # (0.930 -> 0.925 -> 0.919) -- matches the literature nuance that the retardation term matters
    # proportionally more for heavier atoms, even though Gaunt still dominates throughout this range.
    # All six numbers reproduce the 29-Jul-2026 reading exactly, across three widely separated Z; together
    # with branch 2 this is what establishes that the parity repair of branch 1 restored the nu = L block
    # generally, and not just for the one case it was diagnosed on.
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

elseif  false
    # Last successful:  14-Aug-2026: the correction is quadratic in `factor` to 0.04%. A two-term fit
    # A*f^2 + B*f^4 through the five points gives A = -6.386609e-02, B = +2.8208e-03 and reproduces every
    # one of them to ~2e-5 of the leading term, so the ratio drifts only from -0.063822 (f=0.125) to
    # -0.061045 (f=1) over an eightfold range of f. That is the signature the implementation must have and
    # is the sharpest check in this file. Delta E(Breit, omega=0) = 5.24539117 Ha.
    #   Note that these constants are NOT the ones first recorded on 14-Aug-2026 (A = -7.222400e-02):
    # those were measured while the nu = L Gaunt block was zeroed by the CL_reduced_me parity rule, see
    # branch 1. The quadratic LAW held in both cases -- as it must, since the retardation part was never
    # affected -- which is worth knowing about this test: it checks the frequency machinery, not the
    # magnitude of the Breit interaction it multiplies.
    # Branch 4 (frequency dependence -- INTERNAL consistency, no external reference needed): retardation
    # enters at O(omega^2) and CoulombBreit(factor) scales omega, so [dE(factor) - dE(0)] / factor^2 has to
    # be constant up to an O(factor^2) remainder. This tests the W kernel of Grant & Pyper eq. (6), the
    # collapse of the 'S' coefficients onto the multipole L, and the radial wiring all at once -- a wrong
    # kernel, a wrong power of omega or a mis-paired coefficient all break the quadratic law.
    refConfigs = [Configuration("[Ne] 3s^2 3p^5")]
    nModel     = Nuclear.Model(54.)
    grid       = Radial.Grid(true)
    settingsCoulomb = AsfSettings()

    mpCoulomb = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
    basis     = mpCoulomb.levels[1].basis
    eCoulomb  = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb; printout=false).levels[1].energy
    d0        = Hamiltonian.performCI(basis, nModel, grid,
                    AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)); printout=false).levels[1].energy - eCoulomb
    println(">> [Branch 4] E(Coulomb) = $eCoulomb Ha,  Delta E(Breit, omega=0) = $d0 Ha")

    for  f  in  [0.125, 0.25, 0.5, 0.75, 1.0]
        dF = Hamiltonian.performCI(basis, nModel, grid,
                 AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(f)); printout=false).levels[1].energy - eCoulomb
        println(">> [Branch 4] factor = $f:  Delta E(Breit) = $dF Ha,  frequency correction = $(dF-d0) Ha,  " *
                "correction/factor^2 = $((dF-d0)/f^2)")
    end

elseif  false
    # Last successful:  15-Aug-2026: the relative frequency correction is negative throughout and grows
    # monotonically in magnitude, -0.25 % at Z=26, -1.16 % at Z=54 and -2.30 % at Z=79, i.e. roughly as
    # Z^2, which is what omega ~ Delta E / c demands.
    #   THIS BRANCH IS A RULE 12 CASE STUDY, and is worth reading before believing any Breit number on a
    # default grid. Run on the standard grid the Z=79 point came out at +12.71 %, i.e. with the WRONG SIGN,
    # and the SCF reported sign changes for 2p_3/2 and 3p_3/2. Nothing was wrong with the angular or
    # frequency machinery: a Cl-like Z=79 ion sees Z_eff = Z - 16 = 63, so its 3p outer turning point is at
    # r_+ = 0.27 a.u. and the standard grid is some three orders of magnitude wider than the orbitals need.
    # Boxes of 2.0 and 5.0 a.u. agree with the matched box below to every digit, so the plateau is broad --
    # it is only the enormous default that fails. Z=26 and Z=54 are unaffected, being far less overshot.
    #   The earlier readings of -0.94 %, -3.82 % and +8.22 % on this branch predate the CL_reduced_me parity
    # repair of branch 1, which roughly trebled the Breit energy in the denominator.
    # Branch 5 (frequency dependence -- Z scaling): omega = |E_a - E_c| / c grows with Z, so the frequency
    # correction to the Breit interaction must grow along an isoelectronic sequence. Note that it is fed
    # ONLY by exchange-type contributions: for a direct matrix element a = c and b = d, hence omega = 0
    # identically, which is why the correction is much smaller than the Breit energy itself.
    refConfigs = [Configuration("[Ne] 3s^2 3p^5")]
    settingsCoulomb = AsfSettings()
    settingsBreit0  = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.))
    settingsBreit1  = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(1.))

    for  Z  in  [26., 54., 79.]
        ## Rule 12: match the box to the 3p valence orbital, r_+ = (n^2/Z_eff)(1 + sqrt(1 - l(l+1)/n^2))
        ## with n = 3, l = 1 and Z_eff = Z - 16 for a Cl-like ion; take rbox = 2.5 r_+.
        rPlus     = (9.0/(Z - 16.0)) * (1.0 + sqrt(1.0 - 2.0/9.0))
        grid      = Radial.Grid(Radial.Grid(false); rnt = 1.0e-7, h = 5.0e-2, hp = 1.0e-2, rbox = 2.5*rPlus)
        nModel    = Nuclear.Model(Z)
        mpCoulomb = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
        basis     = mpCoulomb.levels[1].basis
        eCoulomb  = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb; printout=false).levels[1].energy
        d0 = Hamiltonian.performCI(basis, nModel, grid, settingsBreit0; printout=false).levels[1].energy - eCoulomb
        d1 = Hamiltonian.performCI(basis, nModel, grid, settingsBreit1; printout=false).levels[1].energy - eCoulomb
        println(">> [Branch 5] Z=$Z:  r_+ = $(round(rPlus,digits=3)) a.u., rbox = $(round(2.5*rPlus,digits=3)) a.u.;  " *
                "Delta E(Breit, omega=0) = $d0 Ha,  Delta E(Breit, omega/=0) = $d1 Ha,  " *
                "frequency correction = $(d1-d0) Ha  ($(round(100*(d1-d0)/abs(d0),digits=2)) %)")
    end

end
