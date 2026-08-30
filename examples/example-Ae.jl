#
println("Ae) Apply & test the QED model corrections to the level structure of atoms and ions.")

if  true
    # Last successful:  29-Jul-2026: E1(Coulomb) = -1519.352175 Ha. Delta E(Petersburg) = +1.609769 Ha
    # (~43.8 eV -- a physically plausible 1s Lamb-shift magnitude for Z=54, correct sign: self-energy
    # dominates and RAISES the energy). Delta E(Sydney) = +1.241483 Ha -- SAME SIGN as Petersburg, ratio
    # 0.771 (previously -0.359361 Ha opposite sign; then +0.212221 Ha, ratio 0.132, after the first Sydney
    # completion). QedSydney was completed by implementing the 3 previously-missing terms from Flambaum &
    # Ginges, PRA 72, 052115 (2005): qedElectricFormFactor (eq. 10), qedMagneticFormFactor (eq. 7),
    # qedWichmannKrollSimple (eq. 12) -- all in module-RadialIntegrals.jl, wired into
    # InteractionStrengthQED.qedLocal. Also fixed a real transcription bug in qedLowFrequency's B(Z)
    # coefficient (0.035 -> 0.35). A SECOND round of fixes then resolved a Z-dependent sign instability
    # (see Branch 3): qedElectricFormFactor's small-distance damping used FG's stated alpha^2 exponent,
    # which over-suppresses the term at Z gtrsim 50 where it matters most; switched to Rci-Q's stated
    # alpha^3, and replaced FG's generic (n=5-calibrated, per Rci-Q's own text) A(Z,r) coefficient with
    # Rci-Q's Table I per-n fit (arXiv:2512.01515) for s-orbitals. Sign is now correct and STABLE across
    # Z=30-100 (see Branch 3), though the magnitude ratio still declines from 0.89 (Z=30) to 0.40 (Z=100) --
    # a smooth, understood residual gap between two independent approximate models, not an instability.
    # Branch 1 (small, technical-smoothness check): H-like Xe (Z=54, 1s), matching example-Ad.jl's Z choice
    # for a coherent cross-reference. Compares NoneQed() (pure Coulomb baseline), QedPetersburg()
    # (Shabaev/Volotka model self-energy + Uehling VP), and QedSydney() (now-complete Flambaum-Ginges
    # radiative potential: Uehling + Wichmann-Kroll + electric form factor + magnetic form factor +
    # low-frequency self-energy).
    name       = "H-like Xe -- Branch 1 (1s)"
    refConfigs = [Configuration("1s")]
    nModel     = Nuclear.Model(54.)
    grid       = Radial.Grid(true)

    settingsNone       = AsfSettings()
    settingsPetersburg = AsfSettings(AsfSettings(); qedModel=QedPetersburg())
    settingsSydney     = AsfSettings(AsfSettings(); qedModel=QedSydney())

    mpNone = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsNone; printout=false)
    basis  = mpNone.levels[1].basis
    mpNoneCI       = Hamiltonian.performCI(basis, nModel, grid, settingsNone;       printout=false)
    mpPetersburgCI = Hamiltonian.performCI(basis, nModel, grid, settingsPetersburg; printout=false)
    mpSydneyCI     = Hamiltonian.performCI(basis, nModel, grid, settingsSydney;     printout=false)

    dPetersburg = mpPetersburgCI.levels[1].energy - mpNoneCI.levels[1].energy
    dSydney     = mpSydneyCI.levels[1].energy     - mpNoneCI.levels[1].energy
    println(">> [Branch 1] E1(Coulomb) = $(mpNoneCI.levels[1].energy) Ha")
    println(">> [Branch 1] Delta E(Petersburg) = $dPetersburg Ha,  Delta E(Sydney) = $dSydney Ha,  Sydney/Petersburg ratio = $(round(dSydney/dPetersburg,digits=3))")
    println(">> [Branch 1] QED(Petersburg)/E(Coulomb) = $(round(dPetersburg/mpNoneCI.levels[1].energy,digits=6))")

elseif  false
    # Last successful:  29-Jul-2026: E1(Coulomb) = -3363.898765 Ha. Delta E(Breit) = 2.365959 Ha,
    # Delta E(QED, Petersburg) = 1.763950 Ha, Delta E(Breit+QED combined) = 4.129909 Ha -- EXACTLY additive
    # (residual = 0.0 Ha). Caveat: this [He] 2s system is a single-CSF case (no CI mixing), so exact
    # additivity here is largely structural -- both corrections are simply added into the same trivial 1x1
    # Hamiltonian "matrix" before diagonalization. A genuinely multi-CSF system could show a nonzero
    # residual from CI-mixing cross-terms; not tested here.
    # Branch 2 (medium, Breit+QED combined): Li-like Xe (Z=54, [He] 2s valence -- n=2, kappa=-1, well
    # within QedPetersburg's documented n<=4 validity range). Four-way comparison: Coulomb-only,
    # Coulomb+Breit, Coulomb+QED(Petersburg), Coulomb+Breit+QED(Petersburg) together -- answers whether
    # Breit+QED combine additively (Delta E(Breit)+Delta E(QED) ~= Delta E(combined)) or show real
    # cross-terms, directly addressing "Breit+QED are often discussed together."
    name       = "Li-like Xe -- Branch 2 ([He] 2s, Breit+QED combined)"
    refConfigs = [Configuration("[He] 2s")]
    nModel     = Nuclear.Model(54.)
    grid       = Radial.Grid(true)

    settingsCoulomb      = AsfSettings()
    settingsBreit        = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.))
    settingsQed          = AsfSettings(AsfSettings(); qedModel=QedPetersburg())
    settingsBreitAndQed  = AsfSettings(AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.), qedModel=QedPetersburg())

    mpCoulomb = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsCoulomb; printout=false)
    basis     = mpCoulomb.levels[1].basis
    mpCoulombCI     = Hamiltonian.performCI(basis, nModel, grid, settingsCoulomb;     printout=false)
    mpBreitCI       = Hamiltonian.performCI(basis, nModel, grid, settingsBreit;       printout=false)
    mpQedCI         = Hamiltonian.performCI(basis, nModel, grid, settingsQed;         printout=false)
    mpBreitAndQedCI = Hamiltonian.performCI(basis, nModel, grid, settingsBreitAndQed; printout=false)

    dBreit    = mpBreitCI.levels[1].energy       - mpCoulombCI.levels[1].energy
    dQed      = mpQedCI.levels[1].energy         - mpCoulombCI.levels[1].energy
    dCombined = mpBreitAndQedCI.levels[1].energy - mpCoulombCI.levels[1].energy
    println(">> [Branch 2] E1(Coulomb) = $(mpCoulombCI.levels[1].energy) Ha")
    println(">> [Branch 2] Delta E(Breit) = $dBreit Ha,  Delta E(QED) = $dQed Ha,  Delta E(Breit+QED combined) = $dCombined Ha")
    println(">> [Branch 2] Sum of separate shifts = $(dBreit+dQed) Ha,  vs. combined = $dCombined Ha,  " *
            "non-additive residual = $(dCombined-(dBreit+dQed)) Ha")

elseif  false
    # Last successful:  29-Jul-2026: widened Z-scan (1s), AFTER the second round of QedSydney fixes (see
    # Branch 1 comment -- alpha^3 damping + Rci-Q Table I coefficient in qedElectricFormFactor):
    # Delta E(Petersburg)/Delta E(Sydney)/ratio = Z=30: +0.230363/+0.205793/+0.893;
    # Z=54: +1.609769/+1.241483/+0.771;  Z=74: +4.572163/+2.897395/+0.634;
    # Z=92: +9.640799/+4.613263/+0.479;  Z=100: +12.914564/+5.183285/+0.401.
    # RESOLVED: the sign is now correct and STABLE across this entire Z=30-100 range -- no more flips (the
    # first-fix version flipped negative at Z=74/92, see git history / memory for that intermediate result).
    # The ratio still declines smoothly with Z (0.89 -> 0.40), i.e. Sydney increasingly UNDER-predicts
    # Petersburg's self-energy at higher Z -- but this is now a well-behaved, monotonic, understood residual
    # gap between two independent approximate models (each individually only claimed accurate to a few
    # percent against rigorous QED, not against each other), not the earlier sign-flipping instability.
    # Further improving the high-Z agreement would need Rci-Q's remaining, un-implemented pieces (magnetic/
    # low-frequency per-shell re-fits for l>0, the full Wichmann-Kroll correction, finite-nuclear-size
    # correction) -- a bigger transcription project, deferred (see memory), not needed to resolve the sign
    # instability that was the original concern. Validity-limit probe (Z=54, 5s, n=5, outside Petersburg's
    # n<=4 self-energy range): Delta E(Petersburg) = -0.002334 Ha -- small, negative, physically sensible
    # (self-energy forced to exactly 0 by the code's own guard; only a small VP/Uehling contribution
    # survives). Graceful degradation, not an error -- confirms Petersburg is technically safe to use
    # outside its stated validity range, it just silently drops to VP-only there.
    # Branch 3 (Z-scan + validity-limit probe): (a) H-like ions across Z=30-100 (widened from the original
    # 54/74/92 grid to give genuine confidence after the second fix), tracking how the Petersburg/Sydney
    # agreement evolves with Z; (b) one explicit case OUTSIDE QedPetersburg's documented validity range
    # (H-like Xe in a 5s configuration, n=5 > the n<=4 limit) to concretely show what happens there --
    # InteractionStrengthQED.selfEnergyVolotka returns EXACTLY 0. for n>4 by explicit code guard (not an
    # error, not a silently wrong nonzero value), so only the Uehling vacuum-polarization piece survives for
    # Petersburg at n=5; this directly answers "which features can be used technically smoothly."
    # OWN GRID for this branch: Bsplines.checkGridRepresentation refused the one used here and named
    # about 2.4 a.u. as the box these subshells want. A local grid is used rather than
    # changing a shared one, so the branches that are already matched are not disturbed -- a box too
    # LARGE starves the fixed B-spline basis just as badly as one too small (Rule 12). The
    # non-exponential family is used because the exponential one quantises r_max coarsely and cannot
    # be tuned to a target; hp = rbox/300 is Basics.recommendedGrid's own recipe.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-07, h = 5.0e-2, hp = 0.0080, rbox = 2.4)
    settingsNone       = AsfSettings()
    settingsPetersburg = AsfSettings(AsfSettings(); qedModel=QedPetersburg())
    settingsSydney     = AsfSettings(AsfSettings(); qedModel=QedSydney())

    for  Z  in  [30., 54., 74., 92., 100.]
        nModel = Nuclear.Model(Z)
        mpNone = SelfConsistent.performSCF([Configuration("1s")], nModel, grid, settingsNone; printout=false)
        basis  = mpNone.levels[1].basis
        mpNoneCI       = Hamiltonian.performCI(basis, nModel, grid, settingsNone;       printout=false)
        mpPetersburgCI = Hamiltonian.performCI(basis, nModel, grid, settingsPetersburg; printout=false)
        mpSydneyCI     = Hamiltonian.performCI(basis, nModel, grid, settingsSydney;     printout=false)
        dPetersburg = mpPetersburgCI.levels[1].energy - mpNoneCI.levels[1].energy
        dSydney     = mpSydneyCI.levels[1].energy     - mpNoneCI.levels[1].energy
        println(">> [Branch 3] Z=$Z (1s):  Delta E(Petersburg) = $dPetersburg Ha,  Delta E(Sydney) = $dSydney Ha,  " *
                "Sydney/Petersburg ratio = $(round(dSydney/dPetersburg,digits=3))")
    end

    # Validity-limit probe: H-like Xe in a 5s (n=5) configuration -- outside selfEnergyVolotka's n<=4 range
    nModel   = Nuclear.Model(54.)
    mpNone5s = SelfConsistent.performSCF([Configuration("5s")], nModel, grid, settingsNone; printout=false)
    basis5s  = mpNone5s.levels[1].basis
    mpNoneCI5s       = Hamiltonian.performCI(basis5s, nModel, grid, settingsNone;       printout=false)
    mpPetersburgCI5s = Hamiltonian.performCI(basis5s, nModel, grid, settingsPetersburg; printout=false)
    dPetersburg5s = mpPetersburgCI5s.levels[1].energy - mpNoneCI5s.levels[1].energy
    println(">> [Branch 3] Z=54 (5s, OUTSIDE n<=4 validity range):  Delta E(Petersburg) = $dPetersburg5s Ha " *
            "(self-energy term forced to exactly 0 by code guard; only Uehling VP survives)")

elseif  false
    # Last successful:  29-Jul-2026: Be-like Fe (Z=26), 6-CSF basis spanning J=0(x3)/1/2(x2), even parity.
    # QED (Petersburg) shifts EVERY level's energy by ~0.14-0.16 Ha (self-energy-dominated, matches the 1s
    # weight in each CSF), AND shifts the CSF-mixing coefficients in every block that has real mixing:
    # Level 1 (J=0): mc shift up to 4.03e-4;  Level 2 (J=0): up to 2.89e-4;  Level 4/5 (J=2): up to 2.25e-4;
    # Level 6 (J=0): up to 4.06e-4. Level 3 (J=1) shows EXACTLY ZERO mc shift -- this is a clean built-in
    # control, not a null result: that CSF block has only ONE contributing CSF (mc=[0,0,0,1,0,0]), so there
    # is nothing to mix, and QED (a diagonal, per-orbital correction) cannot create mixing where there is
    # no partner CSF to mix with. CONFIRMS: QED enters the CI matrix as a genuine additive diagonal term
    # (same architecture as Breit), which measurably reshapes the eigenvectors (CSF mixing / "representation")
    # in any block with real near-degenerate mixing, not merely a rigid post-hoc shift of level energies --
    # directly answering whether QED "affects level energies + representation in a good sense" (it does,
    # and the zero-shift control case confirms the effect isn't a numerical artifact).
    # Branch 4 (multi-CSF representation check): every branch above uses a single-orbital/single-CSF system,
    # so QED's effect on CSF MIXING (not just a scalar energy shift) had never actually been demonstrated --
    # Branch 2's own comment already flagged this gap. Be-like Fe (Z=26), 2s^2/2p^2 CSFs mixing within the
    # same J=0 even-parity block: QedPetersburg's self-energy is diagonal per orbital (selfEnergyVolotka
    # only contributes for a.subshell==b.subshell) and differs substantially between s- and p-type orbitals
    # (per the fez_qed table), so it shifts the 2s^2 and 2p^2 diagonal CI-matrix entries by DIFFERENT
    # amounts -- this should measurably change the CSF-mixing coefficients of the resulting level(s), not
    # just their energies, directly confirming QED enters the "representation" (eigenvectors), not just a
    # rigid post-hoc energy shift.
    name       = "Be-like Fe -- Branch 4 (2s^2/2p^2 mixing, QED effect on representation)"
    refConfigs = [Configuration("[He] 2s^2"), Configuration("[He] 2p^2")]
    nModel     = Nuclear.Model(26.)
    grid       = Radial.Grid(true)

    settingsNone       = AsfSettings()
    settingsPetersburg = AsfSettings(AsfSettings(); qedModel=QedPetersburg())

    mpNone = SelfConsistent.performSCF(refConfigs, nModel, grid, settingsNone; printout=false)
    basis  = mpNone.levels[1].basis
    mpNoneCI       = Hamiltonian.performCI(basis, nModel, grid, settingsNone;       printout=false)
    mpPetersburgCI = Hamiltonian.performCI(basis, nModel, grid, settingsPetersburg; printout=false)

    for  i = 1:length(mpNoneCI.levels)
        lvNone = mpNoneCI.levels[i];   lvQed = mpPetersburgCI.levels[i]
        println(">> [Branch 4] Level $i, J=$(lvNone.J), parity=$(lvNone.parity):")
        println(">>   E(no QED) = $(lvNone.energy) Ha,  mc(no QED) = $(lvNone.mc)")
        println(">>   E(QED)    = $(lvQed.energy) Ha,   mc(QED)    = $(lvQed.mc)")
        println(">>   mc shift  = $(lvQed.mc - lvNone.mc)")
    end

end
