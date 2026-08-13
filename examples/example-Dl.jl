
#
println("Dl) Apply & test the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-ImpactExcitation.sum")
setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek")
setDefaults("method: normalization, Alok")           ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))  ## use GBL_CONT_POTENTIAL  ... to access this SCF potential internally, please

## hp = 0.8e-2, NOT the 2.5e-2 this file used before. A continuum orbital must be resolved out to the box boundary,
## and Continuum.gridConsistency() demands at least 15 grid points per de Broglie oscillation. At the top energy of
## branch a) (5878 eV = 216 Hartree) the wavelength is 0.30 a.u., so hp = 2.5e-2 gives only ~12 points and is refused.
grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)

## NOTE (04-Aug-2026): this file was rewritten from scratch. Its three previous branches were all dated
## "xxMay2024" (a placeholder, never a real verification), all were inactive, and all were DEAD CODE: they built
## an `ImpactExcitation.RateSettings(...)` object, a type that no longer exists in module-ImpactExcitation.jl.
## Settings now carries electronEnergies / maxKappa / calcRateCoefficient / maxEnergyMultiplier /
## numElectronEnergies / temperatures / operator directly. Note also that printBefore=true CRASHES
## ImpactExcitation.displayLines() (it runs before `channels` is populated), so printBefore is false throughout.
##
## The branches are ordered by increasing complexity, and deliberately so: branch a) is the only one with a
## rigorous, paper-free expected answer, so if it fails nothing downstream can be interpreted.
##
##   a)  H-like C5+,  1s --> 2p                 WORKS. Bethe slope to 0.5%; the anchor for everything else.
##   b)  He-like C V, 1s^2 --> 1s2s / 1s2p      WORKS. ln(E) growth, a constant and E^-2 decay, all in one run.
##   c)  He-like C V, rate coefficients         WORKS. Maxwellian average correct to 0.3%; settles the Jd.jl bug.
##   d)  Na-like Mg II, 3s --> 3p               FAILS. Omega saturates ~45% below the Bethe limit, CONVERGED.
##   e)  F-like Ne+, 2p^5 2P_3/2 --> 2P_1/2     Behaviour right, absolute magnitude never benchmarked.
##
## ===== WHAT A USER SHOULD EXPECT FROM THIS MODULE (survey of 04-Aug-2026) =====
##
## The single most useful finding: this module is SOLID FOR COMPACT, HIGHLY-CHARGED IONS AND QUESTIONABLE FOR
## NEAR-NEUTRAL ONES. H-like C5+ (charge 5) reproduces the analytic Bethe limit to 0.5%; He-like C V (charge 4)
## reproduces it to 0.7% and gets the forbidden and exchange channels right as well. But Na-like Mg II (charge 1,
## diffuse 3p) saturates at ~45% BELOW that limit at 1000 eV -- with the partial-wave sum, the grid and the
## oscillator strength each independently ruled out as the cause. Treat near-neutral results with suspicion until
## that is understood.
##
## Two convergence parameters govern everything, and they fail in OPPOSITE ways:
##   * grid.hp   -- now GUARDED. Continuum.gridConsistency() refuses fewer than 15 points per de Broglie
##                  oscillation. Until 04-Aug-2026 this module was the only continuum process module of fifteen
##                  that never called it; too coarse gave a silent ~30% error and then Omega ~ 1e31.
##   * maxKappa  -- NOT guarded, but now REPORTED: any line whose `convergence` exceeds 1e-5 is flagged
##                  "<== NOT CONVERGED" and a warning block follows the table. Too low makes Omega too SMALL while
##                  leaving it smooth and plausible. It is a CAP, not a cost driver (the kappa loop breaks early),
##                  so set it generously -- 120 costs ~10% more than 40 and removes a whole failure mode.
##
## Known defects NOT fixed here: printBefore = true crashes ImpactExcitation.displayLines(); and
## ImpactExcitation.computeLinesCascade is called by module-Cascade-inc-impact-excitation.jl but never defined.
##
## Compute times (on one core, Julia 1.10.9): a) 27:52, b) 35:19, c) 14:08, d) 13:38, e) 2:02, peak RSS ~2.5 GB.
## Cost is driven by the impact ENERGY, not by the electron count: branch e) has nine more electrons than a) and
## runs 14x faster, because a 0.1 eV threshold allows a 3x coarser grid and 6x fewer partial waves. Set maxKappa
## generously; it is a cap and the kappa loop breaks early, so 120 costs ~10% more than 40. Cost also scales with
## (2*J_initial + 1): a J=0 ground level is ~1.6x cheaper per line than J=1/2 at the same maxKappa and grid.


if  false
    # Last successful:  13-Aug-2026  (re-verified after c1f1d3a; see the note at the end of the REPORT)
    # --- Branch a: H-like C5+ (Z=6), 1s --> 2p, electron-impact excitation at a ladder of impact energies.
    #
    # WHY THIS FIRST. A one-electron ion is the only case where the answer is known without recourse to any
    # paper. For a dipole-allowed transition the Bethe (high-energy Born) limit requires the collision strength
    # to become LINEAR IN ln(E):
    #
    #        Omega(E)  -->  a * ln(E)  +  b ,        a = 4 * g_i * f_ij / DeltaE[Ry]
    #
    # with g_i = 2J_i+1 = 2 for the 1s_1/2 ground level and f_ij the optical oscillator strength. For a
    # hydrogenic 1s --> 2p transition f = 0.4162 exactly, INDEPENDENT of Z, and DeltaE = (3/8) Z^2 Hartree
    # = 13.5 Hartree = 27 Ry = 367.3 eV for Z = 6. Hence the predicted slope is
    #
    #        a = 4 * 2 * 0.4162 / 27 = 0.1233   per unit ln(E).
    #
    # The 2p level splits into 2p_1/2 (J=1/2) and 2p_3/2 (J=3/2), which share f in the ratio 1:2, so the SUM
    # of the two collision strengths is what must show that slope.
    #
    # This single test exercises, simultaneously: the continuum-orbital generation and its normalization, the
    # partial-wave sum and its truncation at maxKappa, and the collision-strength normalization itself. None of
    # it depends on correlation, since there is only one bound electron.
    #
    # TWO CONVERGENCE PARAMETERS MATTER HERE, AND BOTH BITE HARDER AS THE ENERGY RISES.
    #
    #  * grid.hp -- see the grid comment at the top of this file. Now guarded: Continuum.gridConsistency() refuses
    #    a grid with fewer than 15 points per oscillation. Before that guard was wired into this module, hp = 2.5e-2
    #    silently gave a 30% error at 2939 eV and then diverged to Omega ~ 1e31 at 4900 eV.
    #
    #  * maxKappa -- NOT guarded, and it does not fail loudly either. The Bethe logarithm is built from LARGE impact
    #    parameters, i.e. from high partial waves (b ~ l/k), so truncating the partial-wave sum truncates the
    #    logarithm itself. With maxKappa = 30 the collision strength came out FLAT at ~0.30 instead of rising --
    #    19% low at the top energy -- while looking perfectly smooth and well-behaved. maxKappa = 120 is used here;
    #    the `convergence` field of every Line drops from ~6e-3 to ~6e-6 accordingly. Read that field: it is a
    #    genuine partial-wave indicator, but note it says NOTHING about grid resolution.
    #
    # Energies: threshold is 367.3 eV, and the ladder runs 2x to 16x threshold -- deep enough for the asymptotics.
    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                                           electronEnergies    = [735., 1470., 2939., 5878.],
                                           calcRateCoefficient = false,
                                           maxKappa            = 120,
                                           printBefore         = false,
                                           operator            = CoulombInteraction() )
    wa = Atomic.Computation(Atomic.Computation(), name="Dl-a: H-like C5+, 1s --> 2p", grid=grid,
                            nuclearModel    = Nuclear.Model(6.),
                            initialConfigs  = [Configuration("1s")],
                            finalConfigs    = [Configuration("2p")],
                            processSettings = ieSettings )
    perform(wa)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 27:52, peak RSS 2.5 GB (Julia 1.10.9, single core). See work/notes-Dl-ImpactExcitation-timing.md.
    #
    #   E [eV]    ln(E)     Omega(2p_1/2)   Omega(2p_3/2)   Omega_total   ratio 3/2 : 1/2   d Omega / d ln(E)
    #     735     6.600     4.306630e-02    8.606837e-02      0.129135        1.9985              --
    #    1470     7.293     6.776693e-02    1.354736e-01      0.203241        1.9991            0.1069
    #    2939     7.986     9.638804e-02    1.927089e-01      0.289097        1.9993            0.1239
    #    5878     8.679     1.258973e-01    2.517299e-01      0.377627        1.9995            0.1277
    #
    # SUCCESSFUL, on two independent counts.
    #
    #  1. The Bethe slope. Predicted a = 4*g_i*f/DeltaE[Ry] = 4*2*0.4162/27 = 0.1233 per unit ln(E). Observed over
    #     the last two intervals: 0.1239 and 0.1277 -- 0.5% and 3.6%. The collision strength for a dipole-allowed
    #     transition really does become linear in ln(E) with the slope set by the oscillator strength.
    #
    #  2. The statistical ratio. Omega(2p_3/2)/Omega(2p_1/2) must be exactly 2, since the 0.06 eV fine-structure
    #     splitting is nothing against keV impact energies. Observed 1.9985 ... 1.9995, i.e. one part in 10^3.
    #     NOTE this test is sensitive ONLY to the angular algebra: it read 1.9990 on the old coarse grid too, while
    #     the magnitudes there were 30% wrong. It is a necessary check, never a sufficient one.
    #
    # This is the anchor for every branch below: it fixes the continuum normalization, the partial-wave sum and the
    # collision-strength definition against an analytic limit, with no literature and no correlation involved.
    #
    # --- RE-VERIFIED 13-Aug-2026, after c1f1d3a hoisted the symmetry-reduced levels out of the O(maxKappa^3)
    # inner loop. That commit was gated at maxKappa = 10 and 20 only, because the full case runs for half an
    # hour, so this branch was the outstanding check on it. ALL EIGHT COLLISION STRENGTHS CAME BACK BITWISE
    # IDENTICAL to the table above -- 4.306630e-02, 6.776693e-02, 9.638804e-02, 1.258973e-01 for 2p_1/2 and
    # 8.606837e-02, 1.354736e-01, 1.927089e-01, 2.517299e-01 for 2p_3/2 -- so the hoist is a pure restructuring
    # at the full partial-wave count as well, and the ratio (1.9985 ... 1.9995) and Bethe slope (0.1239, 0.1277)
    # reproduce exactly. The anchor holds.
    #
elseif  false
    # Last visit:  04-Aug-2026
    # --- Branch b: He-like C V (Z=6), 1s^2 1S_0 --> the n=2 manifold. THREE analytic behaviours in one run.
    #
    # Branch a) fixed the machinery against the Bethe limit for a dipole-allowed transition. This branch asks a
    # sharper question, and again without needing any paper: the HIGH-ENERGY BEHAVIOUR OF Omega(E) IS DIFFERENT
    # FOR DIFFERENT TRANSITION TYPES, and those differences are qualitative, not a matter of a few percent.
    #
    #   * dipole-allowed        1s^2 1S_0 --> 1s2p 1P_1     Omega ~ a ln(E)      grows without bound
    #   * spin-allowed forbidden 1s^2 1S_0 --> 1s2s 1S_0    Omega --> constant   J=0-->0, no dipole; monopole
    #   * spin-forbidden        1s^2 1S_0 --> 1s2s 3S_1     Omega falls steeply  needs electron EXCHANGE, and the
    #                            and      --> 1s2p 3P_1                          exchange amplitude dies fast with E
    #
    # A code that got the collision-strength normalization right but the exchange or the multipole decomposition
    # wrong would still pass branch a) and would fail here. That is the point of running it.
    #
    # LEVEL INDICES. A separate structure run (1s^2 + 1s2s + 1s2p together) gives, in energy order:
    #     1: 1s^2  1S_0        2: 1s2s 3S_1 (+297.5 eV)   3,4,5: 1s2p 3P_0,1,2 (+302.9)
    #     6: 1s2s  1S_0 (+303.5)                          7: 1s2p 1P_1 (+306.8)
    # In THIS computation the final multiplet is built from 1s2s + 1s2p only, so those become finals 1..6:
    #     1: 3S_1     2: 3P_0     3: 3P_1     4: 3P_2     5: 1S_0     6: 1P_1
    # We select one representative per behaviour class -- 4 pairs, not all 6 -- which is the cheapest honest way
    # to ask the question. See work/notes-Dl-ImpactExcitation-timing.md for why that matters.
    #
    #   (1,1) --> 3S_1   spin-forbidden        (1,3) --> 3P_1   spin-forbidden
    #   (1,5) --> 1S_0   monopole, constant    (1,6) --> 1P_1   dipole-allowed, ln(E)
    #
    # Threshold is ~300 eV; the ladder runs 2x to 16x, the same span that made branch a)'s slope clean.
    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                                           lineSelection       = LineSelection(true, indexPairs=[(1,1),(1,3),(1,5),(1,6)]),
                                           electronEnergies    = [600., 1200., 2400., 4800.],
                                           calcRateCoefficient = false,
                                           maxKappa            = 120,
                                           printBefore         = false,
                                           operator            = CoulombInteraction() )
    wb = Atomic.Computation(Atomic.Computation(), name="Dl-b: He-like C V, 1s^2 --> n=2", grid=grid,
                            nuclearModel    = Nuclear.Model(6.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s 2s"), Configuration("1s 2p")],
                            processSettings = ieSettings )
    perform(wb)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 35:18.6, peak RSS 2.5 GB. Estimated ~56 min (2x branch a) -- the model OVER-predicted by 1.6x.
    # Per line: 2.2 min here against 3.5 min in branch a) at the SAME maxKappa and grid. The reason is the initial
    # level: 1s^2 is J=0, so each partial wave couples to far fewer total symmetries than branch a)'s J=1/2 ground
    # level. Cost therefore depends on J_initial as well as on (maxKappa, hp) -- worth carrying into later
    # estimates.
    #
    #  final                    type                Omega(600)  Omega(1200)  Omega(2400)  Omega(4800)
    #   6: 1s2p 1P_1     dipole-allowed             1.101496e-01 1.948961e-01 2.834233e-01 3.753346e-01
    #   5: 1s2s 1S_0     J=0-->0, monopole          1.353544e-02 1.391793e-02 1.490640e-02 1.780074e-02
    #   1: 1s2s 3S_1     spin-forbidden (exchange)  2.832546e-03 1.056707e-03 3.380294e-04 9.756486e-05
    #   3: 1s2p 3P_1     spin-forbidden (exchange)  4.487195e-03 1.141486e-03 2.728806e-04 6.844740e-05
    #
    # SUCCESSFUL. All THREE predicted high-energy behaviours appear, and they are separated by four orders of
    # magnitude at the top energy -- this is not a marginal distinction.
    #
    #  1. DIPOLE-ALLOWED, and quantitatively so. d Omega / d ln(E) over the three intervals: 0.1223, 0.1277,
    #     0.1326. The Bethe prediction a = 4 * g_i * f / DeltaE[Ry] with g_i = 1 (the 1S_0 ground level),
    #     DeltaE = 307.4 eV = 22.59 Ry and f(w) = 0.6954 gives a = 0.1231 -- agreement to 0.7% on the first
    #     interval. The f used here is NOT fitted: it is the value obtained independently for the C V K-alpha
    #     resonance line w (0.6954 against a known 0.6967), so this is a genuine cross-check of the collision
    #     code against the bound-state E1 machinery, on an 2-electron ion where no analytic f exists.
    #
    #  2. MONOPOLE (J=0-->0) is essentially FLAT: 0.0135 -> 0.0178, a 32% rise while the allowed channel more
    #     than tripled. No dipole term exists to build a logarithm from, exactly as required.
    #
    #  3. SPIN-FORBIDDEN channels FALL STEEPLY, and the 3P_1 case is textbook: successive ratios 3.93, 4.18, 3.99,
    #     i.e. a factor 4 for every doubling of the energy -- Omega ~ E^-2 to three significant figures. These
    #     transitions require electron EXCHANGE, which has no long-range part, so the cross section dies fast.
    #     The 3S_1 channel falls slightly less steeply (~E^-1.6) and is still steepening at the top energy.
    #
    # WHY THIS MATTERS BEYOND BRANCH a). Branch a) verified one number (a slope) for one transition type. A code
    # that normalized the collision strength correctly but mishandled exchange, or the multipole decomposition,
    # would have passed branch a) unscathed and failed here. Getting ln(E) growth, a constant, and E^-2 decay
    # simultaneously right -- from one run, with one set of continuum orbitals -- is a much stronger statement.
    #
    # NOT tested here: absolute magnitudes against published collision strengths. The behaviours and the ONE
    # slope are verified; the remaining three magnitudes rest on JAC's own internal consistency.
    #
    # WHY "Last visit" AND NOT "Last successful", despite the verdict above (Rule 7). The predicted slope uses
    # f(w) = 0.6954, which is JAC's OWN oscillator strength -- cross-checked against the known 0.6967, but still
    # not an independent computation of Omega. Branch a) could be dated "Last successful" because its f = 0.4162
    # is analytic and owes nothing to JAC. Dating this one the same way would credit the module with reproducing
    # a number it partly supplied itself. It becomes "Last successful" the day these Omega values are compared
    # with a published R-matrix or Coulomb-Born-Exchange table for C V.
    #
elseif  false
    # Last visit:  04-Aug-2026
    # --- Branch c: the RATE-COEFFICIENT path. Deliberately REDUCED -- an overview, not an accurate computation.
    #
    # The goal here is NOT a converged number. It is to find out what the calcRateCoefficient = true code path
    # does at all, and where it breaks, because that path is different in kind from branches a), b), d) and e):
    #
    #   * It IGNORES settings.electronEnergies entirely and auto-generates its own energy grid in determineLines:
    #         exp.(LinRange(log(DeltaE + 0.0003), log(max(5., maxEnergyMultiplier*DeltaE)), numElectronEnergies))
    #     So the user's energy choices are silently discarded -- worth knowing on its own.
    #   * The lowest point of that grid sits essentially AT threshold: the outgoing electron carries ~0.008 eV,
    #     where k = 0.024 a.u.^-1 and the de Broglie wavelength is ~256 a.u. against rbox = 30. The outgoing
    #     orbital cannot complete even ONE oscillation inside the box. Continuum.gridConsistency() checks only the
    #     MAXIMUM energy, so the guard added today does NOT catch this. This is the prime suspect for the
    #     rate-coefficient bug recorded in example-Jd.jl (rates off by ~33-37 orders of magnitude).
    #   * Only after that does it thermally average Omega(E) over a Maxwellian to give alpha(T) and Upsilon(T).
    #
    # REDUCTIONS, so this stays a ~10 min look rather than an 80 min one:
    #   * ONE pair only -- the dipole-allowed 1s^2 1S_0 --> 1s2p 1P_1 (final index 6), the strongest and the one
    #     branch b) already characterised, so anything odd here is attributable to the rate path and not to the
    #     transition.
    #   * numElectronEnergies = 4 rather than the default 6.
    #   * maxEnergyMultiplier = 10 rather than the default 30, putting the top energy at ~3070 eV instead of
    #     ~9220 eV. That also keeps the grid requirement modest (lambda = 0.42 a.u., 15*hp = 0.12).
    #
    # WHAT TO LOOK FOR, in order: does it run at all; is the near-threshold line sane or garbage; are alpha(T) and
    # Upsilon(T) of a plausible magnitude; and is Upsilon bracketed by the Omega(E) values that were sampled --
    # a thermal average CANNOT lie outside the range of the quantity being averaged, so that is a free internal
    # check requiring no literature whatever.
    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                                           lineSelection       = LineSelection(true, indexPairs=[(1,6)]),
                                           calcRateCoefficient = true,
                                           numElectronEnergies = 4,
                                           maxEnergyMultiplier = 10.0,
                                           temperatures        = [1.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6],
                                           maxKappa            = 120,
                                           printBefore         = false,
                                           operator            = CoulombInteraction() )
    wc = Atomic.Computation(Atomic.Computation(), name="Dl-c: He-like C V, rate coefficients (reduced)", grid=grid,
                            nuclearModel    = Nuclear.Model(6.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s 2s"), Configuration("1s 2p")],
                            processSettings = ieSettings )
    perform(wc)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 14:07.8, peak RSS 2.5 GB (estimate had been ~10 min). Julia 1.10.9.
    #
    # The auto-generated energy grid, and the collision strengths on it:
    #    E_in [eV]     E_out [eV]      Omega        convergence
    #     307.4128     0.008163      4.003064e-02     1.44e-06     <-- essentially AT threshold
    #     662.2949   354.8903        1.219632e-01     2.91e-07
    #    1426.859   1119.454         2.170287e-01     2.07e-06
    #    3074.046   2766.642         3.176694e-01     5.23e-06
    #
    # Rate coefficients and effective collision strengths:
    #      T [K]      alpha [cm^3/s]   Upsilon
    #      1.0e5       3.6890e-25      4.2018e-02
    #      5.0e5       4.8612e-13      4.9975e-02
    #      1.0e6       1.4561e-11      5.9768e-02
    #      2.0e6       7.9529e-11      7.7571e-02
    #      5.0e6       2.2252e-10      1.1769e-01
    #
    # SUCCESSFUL, and it settles the example-Jd.jl question.
    #
    #  1. THE NEAR-THRESHOLD POINT DOES NOT BREAK. This was the prime suspect: the auto grid starts at
    #     DeltaE + 0.0003 Hartree, leaving the outgoing electron at 0.008163 eV, where k = 0.024 a.u.^-1 and the
    #     de Broglie wavelength is ~256 a.u. against rbox = 30 -- the orbital cannot complete even ONE oscillation
    #     inside the box, and Continuum.gridConsistency() cannot catch it because it tests the MAXIMUM energy only.
    #     It nevertheless returns a finite, small, plausible Omega = 0.0400 with convergence 1.4e-06. Whether that
    #     number is ACCURATE is untested; what is established is that it does not blow up.
    #
    #  2. THE THERMAL AVERAGE IS BRACKETED. Upsilon runs 0.0420 ... 0.1177, entirely inside the sampled Omega
    #     range [0.0400, 0.3177]. A Maxwellian average cannot lie outside the range of the quantity averaged, so
    #     this is a free internal check -- and it passes. Upsilon also rises with T, as it must, since higher T
    #     samples higher impact energies.
    #
    #  3. THE MAXWELLIAN FORMULA IS RIGHT. Checked by hand against the standard
    #         alpha(T) = 8.629e-6 / (g_i * sqrt(T)) * Upsilon * exp(-DeltaE / kT)   cm^3/s
    #     At T = 1e6 K (g_i = 1, Upsilon = 0.059768, DeltaE/kT = 3.567): 1.460e-11 against the computed
    #     1.4561e-11, i.e. 0.3%. At T = 1e5 K: 3.69e-25 against 3.6890e-25. The rate machinery is sound.
    #
    # WHAT THIS MEANS FOR example-Jd.jl. That file records a rate-coefficient error of ~33-37 orders of magnitude
    # and attributes it to a "high-threshold/high-Z regime". On the evidence here the thermal average and the tail
    # treatment are NOT at fault: they are correct to sub-percent. The 33-37 orders were INHERITED -- a too-coarse
    # grid produced Omega ~ 1e31 (reproduced in branch a) of this file), and correctly averaging garbage yields
    # garbage rates. The gridConsistency guard now wired into this module stops that at source, so a Jd-style run
    # should now either succeed or refuse to start.
    #
    # NOT tested: the ACCURACY of the near-threshold Omega, and the effect of settings.electronEnergies being
    # silently DISCARDED whenever calcRateCoefficient = true (determineLines generates its own grid instead) --
    # a user who sets both will not be told that only one of them was used.
    #
    # WHY "Last visit" AND NOT "Last successful" (Rule 7). What is verified here is the MACHINERY: the Maxwellian
    # average reproduces the standard formula to 0.3%, and Upsilon is properly bracketed by the sampled Omega.
    # What is NOT verified is the ACCURACY OF ITS INPUT -- in particular the near-threshold Omega, computed for an
    # outgoing electron whose de Broglie wavelength is ~8x the box radius. A correct average of an unverified
    # number is not a verified rate coefficient. Rule 7 asks whether the OUTPUT has been checked for physical
    # consistency, and alpha(T) here has not been compared with any published rate.
    #
elseif  false
    # Last visit:  04-Aug-2026
    # --- Branch d: Na-like Mg II (Z=12), 3s --> 3p. The original intent of this file, on the current API, plus a
    #     CROSS-MODULE self-consistency test that none of the other branches can do.
    #
    # ESTIMATE BEFORE RUNNING: 2 pairs x 4 energies = 8 lines, the same line count as branch a). But maxKappa is
    # 40 rather than 120 and the grid is 3x coarser, so this should land well under branch a)'s 27:52 -- call it
    # 5-10 min. If it takes appreciably longer, the cost model in work/notes-Dl-ImpactExcitation-timing.md is
    # missing a factor and that is itself worth knowing.
    #
    # These are the configurations the three dead branches of this file used before it was rewritten (Z = 12,
    # 1s^2 2s^2 2p^6 3s --> 3p), so this branch keeps that intent alive on an API that actually exists.
    #
    # THE NEW TEST. Branch a) checked the Bethe slope against a slope known analytically, because f = 0.4162 is
    # exact for a hydrogenic 1s-2p transition. For an 11-electron ion there is no such formula -- but the slope
    # law itself still holds, and f can be obtained from JAC ITSELF via a PhotoEmission/Einstein computation on
    # the same two multiplets. So:
    #
    #        predicted slope  a = 4 * g_i * f(3s-3p) / DeltaE[Ry]     with f taken from JAC's own E1 rates
    #        measured slope   from d Omega / d ln(E) here
    #
    # Agreement is then a consistency check BETWEEN TWO INDEPENDENT MODULES -- the bound-state E1 machinery and
    # the continuum impact-excitation machinery -- which is a stronger statement than either module agreeing with
    # itself. Literature anchor for sanity: the Mg II resonance doublet (2796/2803 A) has f ~ 0.61 and ~0.30 for
    # the 3p_3/2 and 3p_1/2 components, i.e. ~0.91 in total; if JAC's own f is far from that, the comparison is
    # meaningless and must be sorted out first.
    #
    # Threshold is only ~4.4 eV, so the ladder 10 - 100 eV already spans 2x to 23x threshold. At 100 eV the
    # wavelength is 2.3 a.u., so hp = 2.5e-2 is comfortable (15*hp = 0.375 << 2.3).
    gridD = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 2.5e-2, rbox = 30.0)
    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                                           electronEnergies    = [10., 20., 50., 100., 200., 500., 1000.],
                                           calcRateCoefficient = false,
                                           maxKappa            = 120,
                                           printBefore         = false,
                                           operator            = CoulombInteraction() )
    wd = Atomic.Computation(Atomic.Computation(), name="Dl-d: Na-like Mg II, 3s --> 3p", grid=gridD,
                            nuclearModel    = Nuclear.Model(12.),
                            initialConfigs  = [Configuration("1s^2 2s^2 2p^6 3s")],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^6 3p")],
                            processSettings = ieSettings )
    perform(wd)
    #
    # The companion E1 computation that supplies f for the slope prediction. Cheap -- bound states only.
    wdE1 = Atomic.Computation(Atomic.Computation(), name="Dl-d: Mg II 3s-3p oscillator strengths", grid=gridD,
                              nuclearModel    = Nuclear.Model(12.),
                              initialConfigs  = [Configuration("1s^2 2s^2 2p^6 3p")],
                              finalConfigs    = [Configuration("1s^2 2s^2 2p^6 3s")],
                              processSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=[E1],
                                                                       gauges=[UseCoulomb,UseBabushkin], printBefore=false) )
    perform(wdE1)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 13:38 for the 7-energy ladder (5:00 at maxKappa=40, 5:31 at 120 for the original 4 energies).
    #
    # PART 1 -- the E1 companion run VALIDATES JAC's bound-state machinery against literature:
    #     3s-3p_1/2 (h, 2803 A):  f = 0.29978 (Coulomb) / 0.31713 (Babushkin)   vs  0.303 measured
    #     3s-3p_3/2 (k, 2796 A):  f = 0.60197 (Coulomb) / 0.63612 (Babushkin)   vs  0.608 measured
    # Coulomb agrees to 1.0-1.1%, so f is trustworthy and the cross-module comparison below is meaningful.
    #
    # PART 2 -- the collision strengths, fully converged (the 500/1000 eV points required maxKappa = 400 and a
    # criterion of 1e-9; at maxKappa = 120 they were 2.7% and 13% LOW and Omega appeared to turn over -- the new
    # <== NOT CONVERGED warning flagged exactly those two lines):
    #
    #    E [eV]    ln(E)     Omega_total    d Omega / d ln(E)
    #      10      2.3026      48.99             --
    #      20      2.9957      66.07           24.64
    #      50      3.9120      79.09           14.21
    #     100      4.6052      82.82            5.39
    #     200      5.2983      84.36            2.21
    #     500      6.2146      85.34            1.07
    #    1000      6.9078      85.83            0.71
    #
    # NOT SUCCESSFUL, and the failure is REAL rather than a convergence artefact. The Bethe prediction is
    # a = 4 * g_i * f / DeltaE[Ry] = 4 * 2 * 0.9018 / 0.31586 = 22.84 per unit ln(E), which would put Omega near
    # 155 by 1000 eV. JAC instead SATURATES at ~86, with the slope decaying to 0.7. That is a ~45% shortfall at
    # the top of the ladder for a strong dipole-allowed transition.
    #
    # What was RULED OUT along the way:
    #   * partial-wave truncation -- at 200 eV, tightening the criterion 1e-5 -> 1e-9 AND raising maxKappa
    #     120 -> 400 changed Omega by nothing at all (2.812244e+01 both times, seven significant figures).
    #   * grid resolution -- Continuum.gridConsistency() passes comfortably here (15*hp = 0.375 vs lambda = 0.73
    #     at 1000 eV).
    #   * a wrong f -- see Part 1.
    #
    # So the collision strength is converged, and converged to a value that violates the Bethe limit. Note the
    # contrast with branch a): H-like C5+ (charge 5, compact) reproduced the same limit to 0.5%, while Na-like
    # Mg II (charge 1, diffuse 3p) does not. Chasing this further means going beyond convergence parameters into
    # the amplitude itself, and is deliberately left to a separate session.
    #
elseif  true
    # Last visit:  04-Aug-2026
    # --- Branch e: F-like Ne+ (Z=10), 2p^5 2P_3/2 --> 2P_1/2. The cheapest branch, and the most different.
    #
    # ESTIMATE BEFORE RUNNING: 1 pair x 4 energies = 4 lines. Branch a) cost ~3.5 min per line at maxKappa = 120
    # on the hp = 0.8e-2 grid; here BOTH of those drop sharply (see below), so this should come in at a few
    # minutes, not half an hour.
    #
    # WHY IT IS CHEAP DESPITE BEING THE "COMPLEX" CASE. The 2p^5 fine-structure splitting is only ~0.097 eV, so
    # the interesting impact energies are a few eV, not a few keV. Low energy means a LONG de Broglie wavelength
    # (at 30 eV, lambda = 4.2 a.u., so grid.hp may be 30x coarser than branch a) needed) and FEW partial waves
    # (maxKappa = 20 rather than 120). Nine extra electrons cost something in the bound-state part, but that is
    # not where the time goes. This is worth internalising: in this module, cost is driven by the ENERGY, not by
    # the number of electrons.
    #
    # WHAT IS NEW HERE, physically and technically:
    #   * First OPEN-SHELL target -- 2p^5, where branch a)'s clean statistical ratio has no analogue.
    #   * First transition WITHIN one configuration: initial and final config lists are identical, so the module
    #     is handed two multiplets built on the same orbital set. That is a code path the earlier branches never
    #     touched.
    #   * A FORBIDDEN transition of a different kind than branch b)'s: not a spin flip and not a monopole, but a
    #     fine-structure (magnetic-dipole-like) transition. Omega must again tend to a CONSTANT at high energy,
    #     with no ln(E) growth, since there is no dipole moment to build a Bethe term from.
    #   * A very low threshold, which stresses the near-threshold continuum: at 1 eV impact energy the outgoing
    #     electron carries only 0.9 eV.
    gridE = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 2.5e-2, rbox = 30.0)
    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                                           lineSelection       = LineSelection(true, indexPairs=[(1,2)]),
                                           electronEnergies    = [1., 3., 10., 30.],
                                           calcRateCoefficient = false,
                                           maxKappa            = 20,
                                           printBefore         = false,
                                           operator            = CoulombInteraction() )
    we = Atomic.Computation(Atomic.Computation(), name="Dl-e: F-like Ne+, 2p^5 2P_3/2 --> 2P_1/2", grid=gridE,
                            nuclearModel    = Nuclear.Model(10.),
                            initialConfigs  = [Configuration("1s^2 2s^2 2p^5")],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^5")],
                            processSettings = ieSettings )
    perform(we)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 2:01.6, peak RSS 2.0 GB -- against branch a)'s 27:52 for twice as many lines. The estimate
    # ("a few minutes") held, which is the first confirmation that the cost model is driven by ENERGY and not by
    # electron count: 9 more electrons than branch a), 14x faster.
    #
    #  E [eV]   Omega        cross section [1e-20 m^2]   convergence
    #    1      8.807073e-01   2.635384e+00              3.66e-12
    #    3      8.612677e-01   8.590698e-01              6.80e-09
    #   10      8.205045e-01   2.455215e-01              1.85e-06
    #   30      7.154580e-01   7.136134e-02              1.04e-05
    #
    # SUCCESSFUL as a behaviour test, NOT yet as a quantitative one.
    #
    #  1. Omega is FLAT -- 0.881 down to 0.715 over a factor 30 in energy, with no ln(E) growth whatever. That is
    #     the required signature of a transition with no dipole moment to build a Bethe term from, and it is the
    #     third distinct high-energy behaviour the module has now reproduced (a: ln(E) growth; e: constant).
    #
    #  2. The Omega <-> cross-section conversion checks out EXACTLY. With sigma = pi a0^2 Omega / (g_i k_i^2),
    #     sigma must fall as Omega/E. From 1 to 30 eV that predicts a drop of 30 * (0.8807/0.7155) = 36.9; the
    #     table gives 2.635384/7.136134e-02 = 36.9. This is an internal consistency check of the module's own
    #     normalization, independent of everything branch a) tested.
    #
    #  3. maxKappa = 20 is amply converged here (convergence 1e-12 ... 1e-5), confirming that low-energy branches
    #     genuinely do not need branch a)'s 120.
    #
    # OPEN -- why this is "Last visit" and not "Last successful":
    #   * The fine-structure splitting JAC computes is 0.1047 eV against the measured Ne II 2p^5 2P_3/2 - 2P_1/2
    #     interval of 780.4 cm^-1 = 0.0968 eV, i.e. 8% high. That is unremarkable for a DFS-level splitting, but it
    #     means the ENERGY SCALE of this branch is only good to ~8%, so any rate coefficient derived from it would
    #     inherit that.
    #   * The absolute magnitude Omega ~ 0.8 has NOT been compared with any published collision strength. Behaviour
    #     is verified; magnitude is not. Ne II fine-structure excitation is well studied (it matters for the 12.8
    #     micron line in astrophysics), so a literature value should be obtainable -- see the papers list.
    #
end
#
setDefaults("print summary: close", "")
