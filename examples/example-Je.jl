#
println("Je) Plasma.CollisionalRadiativeScheme: kinetic (non-LTE) level populations for a small level set.")


# ORIGINAL DESIGN SKETCH (18-Jul-2026), condensed: this branch started as a pure "ReducedCRScheme" design
# note -- never implemented, never runnable, deliberately just a considered-but-undecided idea for a SMALL,
# closed-form (7x7-style) collisional-radiative population balance among a hand-picked handful of levels
# within one ion, as the honest kinetic replacement for an LTE-Boltzmann placeholder used elsewhere. It was
# judged "completely unrealistic" at the time ("EI rates don't exist at the needed level of detail").
# 02-Aug-2026: re-evaluated after Jd.jl's SatelliteDiagnosticScheme work gave direct, hands-on evidence that
# ImpactExcitation genuinely works in moderate-threshold/moderate-Z regimes (it has a real, but localized,
# bug at high-keV/high-Z thresholds, deferred to its own session) -- the judgment no longer holds. Redesigned
# and implemented as Plasma.CollisionalRadiativeScheme (module-Plasma.jl / module-Plasma-inc-collisional-
# radiative.jl); full design discussion and scope decisions live in that session's plan, not repeated here.
# Key differences from the original sketch: no manual levelConfigs list (generated automatically from
# refConfigs via NoExcitations/upperShellNo, mirroring SahaBoltzmannScheme's own generation policy); no
# automatic Saha-coupling for absolute normalization and no multi-charge-state/ionization coupling (both
# explicitly out of scope -- this scheme answers "how is this ion's population distributed among its
# levels", not "how much of this ion is present"); no isotopicMixture/qRange batching (each ion is solved
# independently, so batching would add convenience only, not capability) -- refConfigs (against
# nuclearModel.Z) alone fixes the element and charge state, exactly as for every other scheme.


if true
    #
    # Last visit:      04-Aug-2026
    # Last successful: unknown ... BLOCKED, see the 04-Aug-2026 report below (no-CI level set)
    #
    # Branch a: 7-level He-like carbon (1s^2, 1s2s, 1s2p), the scheme's own first validation target --
    #   the same small w/x/y/z-adjacent level set the original design sketch above was built around.
    #   NoExcitations=1, upperShellNo=2 regenerates exactly the 3 configurations (1s^2 -> 1S0; 1s2s -> 1S0,
    #   3S1; 1s2p -> 1P1, 3P0,1,2 -- 7 levels total). No level here lies anywhere near the next charge
    #   state's ionization threshold, so the scheme's known AutoIonization-competition limitation (see the
    #   @warn in the driver) is not exercised by this branch.
    #
    # REPORT (03-Aug-2026, re-run from cache after the balance-solve reformulation): unchanged and healthy.
    # Correct coronal-limit behaviour -- the ground state holds essentially the entire population (1.0000 up
    # to Te=3e5 K, 9.9752e-1 at Te=1e6 K) and every n=2 level rises monotonically with Te. Gauge agreement is
    # 4-5 significant figures across the whole table (level 2 at Te=1e4 K: 4.1994e-150 Coulomb vs 4.1959e-150
    # Babushkin; level 5 at Te=1e6 K: 4.4188e-4 vs 4.4185e-4). No negative populations. This branch is also
    # what pins down the cleanup policy in Plasma.solveCollisionalRadiativeBalance: its populations span ~150
    # orders of magnitude (4.2e-150 at Te=1e4 K to 1.8e-3 at Te=1e6 K), all of it genuine exp(-dE/kT)
    # Boltzmann-driven physics, so the solver must never floor on |population| -- only on negatives.
    #
    # REPORT (04-Aug-2026, comparison against known C V values): the CR machinery itself checks out, but the
    # LEVEL SET it is built on does not, so this branch CANNOT yet be dated successful.
    #
    # What passes:
    #   * Oscillator-strength sum rule. f(6->1) + f(7->1) = 0.2339 + 0.4654 = 0.6993 against the known
    #     f(w) = 0.6967 for C V, i.e. +0.4%. The radial integrals and the E1 operator are therefore sound.
    #   * The collisional-radiative balance itself. Level 7 is firmly coronal (A_7 = 6.274e11 /s vs a
    #     collisional depopulation ~1e7 /s), so n7/n1 = C_17*Ne/A_7 inverts cleanly: at Te=1e6 K it gives
    #     C_17 = 1.389e-11 cm^3/s, i.e. an effective collision strength Upsilon_7 = 0.056 (Upsilon_6 = 0.045,
    #     total 0.101). The van Regemorter estimate for this transition is Omega = (8pi/sqrt 3)*f/dE_Ry*gbar
    #     = 0.452*gbar = 0.090 for the usual gbar ~ 0.2. Agreement to ~12% -- well inside van Regemorter's own
    #     factor-of-2. This validates ImpactExcitation's rates, the thermal average, the detailed-balance
    #     de-excitation and the balance solve as one chain.
    #
    # What fails -- the level set carries NO configuration interaction:
    #   * Every level of the cached repMultiplet is a single CSF (mc = 1.0 exactly). Verified by loading the
    #     .jld cache directly and by re-running the generation path. The cause is NOT the per-configuration
    #     loop in Plasma.generateCollisionalRadiativeLevels: Hamiltonian.performCIwithFrozenOrbitals uses, by
    #     its own docstring, "just the diagonal part of the Hamiltonian matrix" and prints "no CI, CSF
    #     diagonal". Its name notwithstanding, it performs no CI. Feeding it ALL configurations in ONE call
    #     changes nothing.
    #   * Consequence for the two J=1 odd levels. SelfConsistent.performSCF on the same configurations gives
    #     mixing coefficients (-0.5741, +0.8188) and (-0.8188, -0.5741), essentially the exact LS vectors
    #     (-sqrt(1/3), +sqrt(2/3)) -- i.e. near-pure 1P1 and 3P1. The no-CI set instead returns the bare jj
    #     CSFs (1s2p_1/2)_1 and (1s2p_3/2)_1. The trace is preserved (-1150.8935 eV either way) but the
    #     splitting collapses from 3.902 eV to 1.330 eV, which puts the missing off-diagonal element at
    #     V = sqrt((3.902^2 - 1.330^2))/2 = 1.83 eV. The J=0 even block loses V = 8.85 eV the same way,
    #     leaving the ground state 0.26 eV too high (-880.05 vs -880.31 eV).
    #   * Observable damage. The w-line strength is split 1:2.007 across the two J=1 levels -- the statistical
    #     jj ratio -- instead of being concentrated in 1P1, so the intercombination line y is not reproduced
    #     at all and any w/x/y/z ratio from this branch is meaningless. Full CI also recovers the known
    #     near-degeneracy of x and y (302.938 vs 302.915 eV, 0.023 eV apart, against 0.02 eV in nature);
    #     the no-CI set puts them 1.33 eV apart.
    #   * Excitation energies improve correspondingly. No-CI: z -0.59%, x -0.57%, w -0.85% (scattered).
    #     Full CI: z -0.50%, x -0.49%, y -0.49%, w -0.35% -- uniform, as a DF calculation missing only
    #     ground-state correlation should be.
    #
    # Second, independent gap: level 5 (1s2s 1S0) can only decay by two-photon 2E1 emission, which JAC does
    # not compute. Its only listed channel is the M1 5->2 at 1.95e-3 /s, giving a 513 s lifetime against a
    # true ~2.7e-5 s (A(2E1) ~ 3.7e4 /s, scaling A ~ Z^6 from the He value 51.3 /s). At Ne=1e15 cm^-3 the
    # collisional coupling out of level 5 is faster than either, so branch a's populations are probably not
    # much affected -- but at low density this would dominate, and the printed lifetime is wrong by ~2e7.
    #
    # CACHE HAZARD found while doing this: collisionalRadiativeLevelFingerprint covers the grid, nuclear
    # model, asfSettings, NoExcitations and upperShellNo -- but NOT the level-generation ALGORITHM. Switching
    # generateCollisionalRadiativeLevels to a genuine CI would leave the fingerprint unchanged, so these
    # .jld caches would be silently reused and the fix would appear to do nothing. RESOLVED 04-Aug-2026 by
    # adding a "generationMethod" version tag to the fingerprint; bump it whenever the generation changes.
    #
    # REPORT (04-Aug-2026, after the frozen-orbital CI fix, commit 7cc164b): the no-CI defect described above
    # is FIXED and this branch now reproduces the C V K-alpha structure properly.
    #   * The resonance line w is back where it belongs. Level 7 (1P1) carries A = 9.45e11 /s (Coulomb) /
    #     8.22e11 /s (Babushkin) with f = 0.6954, against the known f(w) = 0.6967 for C V -- 0.19%. The
    #     literature A(w) ~ 8.9e11 /s sits BETWEEN the two gauges. Previously that strength was split 1:2
    #     across two levels.
    #   * The intercombination line y is now genuinely weak: level 4 (3P1) gives A = 1.59e7 /s, i.e. w/y
    #     ~ 6e4, as it must be for a low-Z He-like ion. Previously y and w differed by a factor 2.
    #   * The x-y near-degeneracy is reproduced: levels 4 and 5 lie 0.023 eV apart (0.02 eV in nature).
    #   * Excitation energies are now uniformly ~0.35-0.50% low (w 306.82 vs 307.90 eV; z 297.47 vs 298.96;
    #     x 302.94 vs 304.42; y 302.91 vs 304.40), the expected signature of a DF calculation missing only
    #     ground-state correlation. Before the fix they were scattered over 0.57-0.85%.
    #   * The populations change physically: with the triplets no longer drained by a spurious E1, the
    #     metastable 1s2s 3S1 rises from 1.78e-3 to 1.09e-2 at Te=1e6 K.
    #   * Coronal cross-check still holds. Inverting n7/n1 = C_17*Ne/A_7 at Te=1e6 K gives Upsilon = 0.068
    #     against a van Regemorter estimate of 0.448*gbar = 0.090 (gbar ~ 0.2) -- 24%, inside van
    #     Regemorter's own factor-of-2.
    #
    # STILL MISSING -- deliberately, this being a first version of CR modelling (4-Aug-2026 decision):
    #   (a) INTER-CONFIGURATION CI. Plasma.generateCollisionalRadiativeLevels diagonalizes ONE configuration
    #       at a time, which is the intended simplification here: it keeps the scheme element-agnostic, keeps
    #       the cost linear in the number of configurations, and lets every level be identified with the
    #       configuration it came from. The price for He-like carbon is the 1s^2 - 1s2s J=0 mixing (an 8.85 eV
    #       off-diagonal element), which leaves the ground state 0.26 eV too high (-880.05 vs -880.31 eV) and
    #       is the main reason the excitation energies are uniformly ~0.4% low. Note this does NOT affect the
    #       w/y structure above: that mixing is INTRA-configuration (both J=1 CSFs live in 1s2p) and is now
    #       correctly included. Options for lifting the restriction are discussed but not implemented; a
    #       user-defined grouping of configurations into CI blocks is the most likely route, and the CR
    #       scheme is not the right place to prototype it.
    #   (b) TWO-PHOTON (2E1) DECAY. Level 6 (1s2s 1S0) cannot decay by any single-photon channel to the
    #       ground state (0 -> 0 is strictly forbidden), so JAC gives it only the M1 6->2 route and a 513 s
    #       lifetime, against a true ~2.7e-5 s (A(2E1) ~ 3.7e4 /s, scaling A ~ Z^6 from the He value 51.3 /s)
    #       -- wrong by ~2e7. At Ne=1e15 cm^-3 collisional coupling out of that level is faster than either
    #       rate, so branch a's populations are probably little affected; at low density this would dominate
    #       and the level would be badly over-populated. This belongs to a multi-photon excitation/decay
    #       module, not here.
    # Both are why this branch stays "Last visit" rather than "Last successful": the agreement above is
    # genuinely good, but two known channels are absent by construction.
    #
    # hp is set by the FASTEST continuum electron this branch asks for, not by the bound orbitals.
    # Continuum.gridConsistency refused hp = 2.0e-2 at 327.33 Hartree: the de Broglie wavelength is
    # 0.245564 a.u., and the guard requires 15 points per oscillation, i.e. hp <= 1.637e-2.
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-5, h=5.0e-2, hp=1.5e-2, rbox=20.0)
    nm   = Nuclear.Model(6.)

    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings(); maxKappa=5, numElectronEnergies=6,
                                           maxEnergyMultiplier=10.0, temperatures=[1e4, 3e4, 1e5, 3e5, 1e6])
    peSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=[E1, M1])
    aiSettings = AutoIonization.Settings()

    # cacheDirectory keeps the auto-written level-/rates-cache files out of examples/ itself; the two
    # filenames below are the caches generated by this very branch on 03-Aug-2026 and are reused on every
    # later run. Both are validated by an EXACT fingerprint match, so any change to the grid, nuclear model,
    # asfSettings, NoExcitations/upperShellNo (level cache) or to the ImpactExcitation settings (rates cache)
    # makes them be skipped with a printed note and silently regenerated -- they can never go stale unnoticed.
    # The rates cache stores the RAW ImpactExcitation lines, i.e. before the thermal average, so changing
    # ieSettings.temperatures alone still reuses it. Note the paths are relative to the working directory the
    # example is run from (examples/), exactly as cacheDirectory itself is.
    scheme = Plasma.CollisionalRadiativeScheme(Plasma.CollisionalRadiativeScheme(); NoExcitations=1, upperShellNo=2,
                                               ieSettings=ieSettings, peSettings=peSettings, aiSettings=aiSettings,
                                               levelsFilenames=["example-Je.dat/newCRLevelsZ6A12-2026-08-04T09.jld"],
                                               ratesFilenames =["example-Je.dat/newCRRatesZ6A12-2026-08-04T09.jld"],
                                               cacheDirectory="example-Je.dat")

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                     refConfigs=[Configuration("1s^2")], asfSettings=AsfSettings(),
                                     settings=Plasma.Settings(0., 1.0e15, true))
    wb = perform(computation, output=true)
    #
elseif false
    #
    # Last visit:      03-Aug-2026
    # Last successful: unknown ... not yet run/verified
    #
    # Branch b: boron-like Ne5+ (1s^2 2s^2 2p^1), NoExcitations=1, upperShellNo=2 -- a genuinely open-shell
    #   (2p^1, richer term structure than branch a's closed/near-closed configs) level set whose fromShells
    #   set is {1s,2s,2p} (all three shells of refConfigs), unlike branch a's single 1s shell. Several of
    #   the generated 1s/2s->nl single excitations sit at inner-shell threshold energies far above what the
    #   chosen temperatures populate, so a number of the resulting IE rate coefficients/effective collision
    #   strengths are expected to underflow toward ~0 (not necessarily structurally/exactly zero) rather
    #   than be forbidden outright -- branch a's fully-connected, closed-shell 7-level system never
    #   exercised this. Chosen (over the same idea at upperShellNo=3, or the same upperShellNo=2 idea for
    #   carbon-like Ne 1s^2 2s^2 2p^2) via Plasma.estimateCollisionalRadiativeCost: N=18 levels/153 pairs/
    #   ~1.4h extrapolated IE cost here, vs. N=119/7021 pairs/~56.5h for upperShellNo=3 (not worth it -- the
    #   near-zero-inner-shell-line behavior is already fully exercised at upperShellNo=2; going to n=3 just
    #   multiplies the pair count) and N=25/300 pairs/~2.3h for carbon-like Ne at upperShellNo=2 (same idea,
    #   higher cost, no added value over boron-like).
    #
    # REPORT (03-Aug-2026, re-run from cache after the balance-solve reformulation): the tiny NEGATIVE
    # populations previously seen for the inner-shell block 11-18 at low Te (e.g. -2.89e-27) are GONE -- all
    # 90 rows (18 levels x 5 temperatures) are now non-negative. Root cause was NOT the rates but the balance
    # solve itself: it sacrificed the LAST level's equation to the sum(n)=1 constraint, leaving the
    # least-populated level determined only indirectly, through a near-singular combination that amplified
    # roundoff into signed noise. Plasma.solveCollisionalRadiativeBalance now sacrifices the GROUND STATE's
    # equation instead (see the long comment there); every weakly-populated level thereby keeps its own
    # well-conditioned equation. Four checks on the new run:
    #   (i)   the 2p ground term relaxes to its statistical ratio -- n(3/2-)/n(1/2-) = 0.657/0.333 = 1.97 at
    #         Te=1e6 K against the exact g-ratio 4/2 = 2, and is already 0.619/0.381 at Te=1e4 K;
    #   (ii)  the block 11-18 is exactly 0.0000e+00 at Te=1e4 K ONLY -- there every feeding rate has genuinely
    #         underflowed to exact 0.0, so an exactly-zero population is the correct answer, not a floored one;
    #   (iii) from Te=3e4 K upward that same block carries real, strictly positive, monotonically rising
    #         population (~1e-155 -> 1e-13), i.e. the fix removed noise WITHOUT erasing physics. This is why
    #         the cleanup in solveCollisionalRadiativeBalance floors negatives only and never |population|:
    #         branch a shows genuine coronal populations as small as 4.2e-150, which any |.|-floor would kill;
    #   (iv)  both gauges agree to 3-4 significant figures throughout the block (level 11 at Te=1e6 K:
    #         8.8301e-13 Coulomb vs 8.8727e-13 Babushkin).
    # Left as "Last visit": these are internal-consistency checks, not a comparison against an independent
    # CR code or literature value for this ion, which is what Rule 7 asks for before dating a branch.
    #
    # SUPERSEDED (04-Aug-2026): every number in the REPORT above was produced BEFORE the frozen-orbital CI
    # fix (commit 7cc164b) and no longer describes what this branch computes. The two .jld files named below
    # predate that fix too; they are now correctly REJECTED by the fingerprint's "generationMethod" tag, so
    # the next run of this branch regenerates them at the full ~1.4 h ImpactExcitation cost and should have
    # its filenames updated here afterwards. The negative-population conclusion itself is unaffected -- that
    # was verified independently on branch a and by the synthetic decoupled-block test -- but the population
    # values, the level ordering and the Te at which the inner-shell block lights up all have to be re-read.
    #
    # KNOWN LIMITATION, unchanged by any of the above and the reason this branch cannot become "Last
    # successful" as it stands: levels 11-18 sit ~880 eV above the ground state, whereas the ionization
    # potential of Ne5+ is only ~158 eV. That whole block is therefore ~720 eV ABOVE the ionization
    # threshold -- these are autoionizing states whose Auger rates (~1e14 /s) would dwarf the radiative decay
    # that the balance currently gives them as their only loss channel. The driver's own @warn says exactly
    # this. Their populations are consequently wrong by construction, not merely unvalidated. Fixing it needs
    # either a level set kept below the ionization threshold, or genuine AutoIonization competition wired
    # into the scheme -- note that scheme.aiSettings is currently ACCEPTED BUT NEVER READ by the driver,
    # so setting it below buys nothing today.
    #
    # hp is set by the FASTEST continuum electron this branch asks for, not by the bound orbitals.
    # Continuum.gridConsistency refused hp = 2.0e-2 at 327.33 Hartree: the de Broglie wavelength is
    # 0.245564 a.u., and the guard requires 15 points per oscillation, i.e. hp <= 1.637e-2.
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-5, h=5.0e-2, hp=1.5e-2, rbox=20.0)
    nm   = Nuclear.Model(10.)

    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings(); maxKappa=5, numElectronEnergies=6,
                                           maxEnergyMultiplier=10.0, temperatures=[1e4, 3e4, 1e5, 3e5, 1e6])
    peSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=[E1, M1])
    aiSettings = AutoIonization.Settings()

    # Cache files generated by this branch on 03-Aug-2026; see branch a for the exact-fingerprint reuse rules.
    # Reusing the rates cache here is what turns the ~1.4 h ImpactExcitation run into a few seconds.
    scheme = Plasma.CollisionalRadiativeScheme(Plasma.CollisionalRadiativeScheme(); NoExcitations=1, upperShellNo=2,
                                               ieSettings=ieSettings, peSettings=peSettings, aiSettings=aiSettings,
                                               levelsFilenames=["example-Je.dat/newCRLevelsZ10A20-2026-08-03T11.jld"],
                                               ratesFilenames =["example-Je.dat/newCRRatesZ10A20-2026-08-03T12.jld"],
                                               cacheDirectory="example-Je.dat")

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                     refConfigs=[Configuration("1s^2 2s^2 2p^1")], asfSettings=AsfSettings(),
                                     settings=Plasma.Settings(0., 1.0e15, true))
    wb = perform(computation, output=true)
end
