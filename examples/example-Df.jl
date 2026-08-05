
println("Df) Apply & test the DielectronicRecombination module with ASF from internally generated initial-, intermediate and final-state multiplets.")

setDefaults("print summary: open", "zzz-DielectronicRecombination.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")
setDefaults("unit: strength", "cm^2 eV")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, asymptotic Coulomb")
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")

## NOTE (04-Aug-2026): this file was rewritten from scratch. Its seven previous branches were ALL dead code:
## every one of them called `Dielectronic.Settings(...)`, but there is no module or alias `Dielectronic`
## anywhere in src/ -- the module is `DielectronicRecombination`. They also passed 8 or 9 positional arguments
## (two branches used mutually inconsistent signatures) to a Settings struct that has THIRTEEN fields, used the
## removed `process = Dierec()` keyword, and relied on a `JAC.` binding that is not exported. Two branches were
## simultaneously `true`, so the second could never run. Five carried "Last successful: unknown", two claimed
## "13May2024". The physics INTENT is preserved below -- the Tu et al. and Xu et al. comparisons are kept as
## targets -- but none of the old code was salvageable.
##
##   a)  He-like C4+ + e- --> Li-like C3+, K-LL          -- the small anchor; every check is internal
##   b)  the same system, plus the rate coefficient alpha(T)
##   c)  branch a)'s system through the RE-STRUCTURED route -- equivalence + the DR satellite spectrum
##   d)  branch b)'s Xe53+ case through the re-structured route -- the Harman regression
##   e)  He-like W72+ K-LL -- the Tu et al. (2016) benchmark, wide J_f spread
##   f)  Be-like Au75+ with an n=19 Rydberg spectator -- the HydrogenicCorrections, ported and repaired
##   g)  H-like C5+ with TWO explicit Rydberg shells -- validates the RydbergTailCorrection's n-scaling
##   h)  K-LL DR of Li-like C3+ with a HYPERFINE-RESOLVED initial ion -- the hyperfine route
##
## THE TWO STRATEGIES OF THIS MODULE, for orientation:
##   * FINE-STRUCTURE resolved (module-DielectronicRecombination-inc-FS-resolved.jl) -- the working path. Two
##     sub-modes: full pathway computation (default), and `calcOnlyPassages`, which collects each resonance's
##     contribution directly and is meant for high-n Rydberg intermediates and semi-empirical corrections.
##   * HYPERFINE resolved (module-DielectronicRecombination-inc-HF-resolved.jl) -- REWRITTEN 05-Aug-2026 and now
##     a peer of the fine-structure path, selected simply by `calcHyperfineResolved = true`. It gives the INITIAL
##     (recombining) ion the full nuclear moments and sets mu = Q = 0 for the intermediate and final levels, whose
##     hyperfine structure is negligible; their F is retained because it carries the angular factors and summed
##     over on display. No electronic amplitude is recomputed -- the route runs the fine-structure machinery and
##     RECOUPLES its amplitudes, so the two cannot drift apart. Every run prints the F-sum rule, which must hold
##     exactly. (Before the rewrite this path could not run at all: three of its functions were called but never
##     defined, and it passed a Multiplet where an Hfs.HfMultiplet was required.)
##
## Note that Settings() defaults to gauges = UseGauge[] -- EMPTY. Set the gauges explicitly or no radiative
## channel is computed at all.

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


if  false
    # Last visit:  05-Aug-2026
    # --- Branch a: He-like C4+ (Z=6) + e- --> Li-like C3+, K-LL dielectronic recombination.
    #
    # WHY THIS FIRST. It is the smallest system in which DR is a genuine three-step process, and small enough
    # that every quantity can be checked by hand from the module's own printout, with no literature at all:
    #
    #   1. ENERGY CONSISTENCY. electronEnergy must equal E(intermediate) - E(initial), and photonEnergy must
    #      equal E(intermediate) - E(final). Both are readable off the level tables that precede the results.
    #
    #   2. THE RESONANCE-STRENGTH FORMULA. The Resonance table prints resonanceEnergy, captureRate, augerRate
    #      and photonRate SEPARATELY from resonanceStrength, so the strength can be recomputed by hand from
    #
    #           S = (h^2 / (2 m_e E_res)) * (g_d / (2 g_i)) * A_a * A_r / (A_a + A_r)
    #
    #      and compared. This is the same style of check that confirmed the Maxwellian average in example-Dl.jl.
    #
    #   3. THE SUM RULE. Summing the per-pathway reducedStrength over all final levels of one intermediate level
    #      must reproduce that intermediate level's resonanceStrength.
    #
    #   4. GAUGE AGREEMENT on everything radiative (photonRate, resonanceStrength). The capture side is
    #      gauge-independent by construction, so a gauge spread here localises to the radiative step.
    #
    # The K-LL resonances are the 1s2l2l' doubly-excited levels of Li-like carbon: a free electron is captured
    # into a 2l orbital while a 1s electron is simultaneously promoted to 2l', and the resulting autoionizing
    # level then stabilises by emitting a photon down to 1s^2 2l. Literature target, not yet compared here:
    # Xu et al., PRA 2016.
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles  = [E1],
                                                    gauges      = [UseCoulomb, UseBabushkin],
                                                    printBefore = false )
    wa = Atomic.Computation(Atomic.Computation(), name="Df-a: K-LL DR into He-like C4+", grid=grid,
                            nuclearModel        = Nuclear.Model(6.),
                            initialConfigs      = [Configuration("1s^2")],
                            intermediateConfigs = [Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            processSettings     = drSettings )
    perform(wa)
    #
elseif  false
    # Last visit:  04-Aug-2026
    # --- Branch b: K-LL DR of H-like Xe53+ (Z=54). THE LITERATURE ANCHOR of this file.
    #
    # This is the one case for which independent reference data sit inside this repository:
    #   apps/b24-apps-schippers-dr-breit-xenon/b24.01-harmann-Xe53DRtheo4.dat  -- Harman's theoretical resonance
    #   energies, strengths and widths, and apps-dr-xenon.report-01..03, which hold earlier JAC output for the
    #   same system. (The app's own script b24.01-Xe53KLLDR.jl is itself dead under the current API: it calls
    #   `Dielectronic.Settings` with 10 positional arguments.)
    #
    # It is also STRUCTURALLY SIMPLER than branch a): the recombining ion is H-like, so ONE electron is present
    # before capture and the resonance is a doubly-excited He-like 2l2l' state with an EMPTY K shell. The K-LL
    # resonances therefore sit near 20.6-21.6 keV, which is what makes the branch expensive rather than difficult.
    #
    # DELIBERATE DIFFERENCE FROM THE APP SCRIPT. Its finalConfigs were {1s^2, 1s2s, 2s2p} -- which omits 1s2p and
    # instead lists 2s2p, an autoionizing configuration identical to one of the intermediates. The stabilizing
    # transition is 2p --> 1s, so 2p^2 --> 1s2p is the main radiative channel of the 2p^2 resonances and must be
    # present; 2s^2 has no E1 channel at all (2s --> 1s is forbidden), and 1s^2 cannot be reached from 2l2l' by a
    # single E1 at all. The final set used here is therefore {1s^2, 1s2s, 1s2p}: 1s^2 is kept only so that its
    # absence of channels is visible rather than assumed.
    #
    # WHAT THIS BRANCH IS FOR. Branch a) verified ratios, a sum rule and the width relation -- all of which are
    # blind to an overall scale factor. The module's reducedStrength carries (2J_m+1)/(2J_i+1) with a factor 2
    # explicitly commented out and marked "factor 2 is not really clear" in
    # module-DielectronicRecombination-inc-FS-resolved.jl. Only an ABSOLUTE comparison can decide that, and this
    # is it.
    grid54     = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    # QED matters at Z = 54: it shifts the K-shell binding by tens of eV and hence the resonance energies. The
    # standard self-consistent field is used; no special SCF field is needed or wanted here.
    asfSet54   = AsfSettings(AsfSettings(), qedModel=QedPetersburg(), eeInteraction=CoulombBreit(1.0))
    #
    # THE BREIT INTERACTION IS NOT OPTIONAL HERE. Run first with the default augerOperator = CoulombInteraction()
    # and the RADIATIVE rates come out within 0.1-1.6% of the 2024 reference, while the AUGER rates fall short by
    # 19%, 39%, 73% and 93% for m = 1, 5, 2 and 6 respectively -- the shortfall growing as the Auger rate itself
    # gets smaller. That is the signature of Breit opening capture channels that the Coulomb interaction alone
    # barely supports, and at Z = 54 it dominates the weak resonances. The application directory this benchmark
    # comes from is named "schippers-dr-breit-xenon" precisely because that is the effect it was built to study.
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles    = [E1],
                                                    gauges        = [UseCoulomb, UseBabushkin],
                                                    augerOperator = CoulombBreit(1.0),
                                                    printBefore   = false )
    wb = Atomic.Computation(Atomic.Computation(), name="Df-b: K-LL DR of H-like Xe53+", grid=grid54,
                            nuclearModel        = Nuclear.Model(54., "Fermi"),
                            # NOTE the explicit ZERO occupations. DielectronicRecombination.checkConsistentMultiplets
                            # requires the initial, intermediate and final states to be built over the SAME ordered
                            # subshell list (the initial one may be shorter), a restriction coming from the angular
                            # coefficients. Writing Configuration("2s^2") would drop 1s from the intermediate subshell
                            # list and abort the run; "1s^0 2s^2 2p^0" keeps it. This is why the app script uses this
                            # notation throughout -- it is load-bearing, not decoration.
                            initialConfigs      = [Configuration("1s^1 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^0 2s^2 2p^0"), Configuration("1s^0 2s^1 2p^1"),
                                                   Configuration("1s^0 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^0 2p^0"), Configuration("1s^1 2s^1 2p^0"),
                                                   Configuration("1s^1 2s^0 2p^1")],
                            initialAsfSettings  = asfSet54, intermediateAsfSettings = asfSet54, finalAsfSettings = asfSet54,
                            processSettings     = drSettings )
    perform(wb)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 1:50, peak RSS 2.5 GB. Ten K-LL resonances, 20.6 - 21.6 keV.
    #
    # COMPARISON WITH HARMAN (apps/b24-apps-schippers-dr-breit-xenon/b24.01-harmann-Xe53DRtheo4.dat):
    #
    #   Harman E [eV]   Harman S      this run (m)   S [cm^2 eV], Coulomb   ratio
    #     20632.4       2.28e-20        m1, 20625      2.2381e-20           0.982
    #     20631.7       7.6e-21         m2, 20628      7.7012e-21           1.013
    #     20655.6       1.75e-20        m3, 20653      1.7102e-20           0.977
    #     21051         7.36e-21        m5, 21048      7.3938e-21           1.005
    #     21095.8       1.05e-21        m6, 21093      1.0826e-21           1.031
    #     21112.3       6.44e-20        m7, 21113      6.3310e-20           0.983
    #     20726.8       1.16e-21        m4, 20721      8.3063e-22           0.716   <-- outlier
    #
    # SUCCESSFUL. Six of seven resonances agree with independent theory to 0.5 - 3.1%.
    #
    #  1. THE FACTOR-2 QUESTION IS SETTLED. reducedStrength carries (2J_m+1)/(2J_i+1) with a factor 2 commented
    #     out and marked "factor 2 is not really clear" in module-DielectronicRecombination-inc-FS-resolved.jl.
    #     Agreement at the 1-3% level across six independent resonances cannot coexist with a missing factor of
    #     two. The factor must stay commented out. Do not reopen this without an equally strong comparison.
    #
    #  2. BREIT IS ESSENTIAL AT Z = 54, and the evidence is clean because this branch was first run WITHOUT it.
    #     With augerOperator = CoulombInteraction() the radiative rates were already within 0.07-1.6% of the
    #     reference while the Auger rates fell short by 19% (m1), 39% (m5), 73% (m2) and 93% (m6) -- the deficit
    #     growing as the Auger rate itself gets smaller. Breit opens capture channels that Coulomb alone barely
    #     supports. Switching to CoulombBreit(1.0) moved the strengths from 18-93% off to 0.5-3.1%.
    #
    #  3. THE OUTLIER, m4 at 20721 eV, 28% low. m4 is the 2p^2 3P_0 level: J = 0 and even parity, so
    #     autoionization into 1s (J = 1/2) admits a SINGLE channel, eps s_1/2, and as a pure triplet its Auger
    #     amplitude is small (1.42e13 against 4.49e14 for m1). Its width is radiative-dominated, so S is very
    #     nearly proportional to A_a and a 28% deficit in the strength is a 28% deficit in that one
    #     near-cancelling single-channel amplitude. Adding 1s2p to the final configurations -- absent from the
    #     2024 app run -- improved it from 32% to 28%, so the missing configuration was NOT the cause. This is a
    #     plausible mechanism consistent with everything else seen here, NOT a demonstration; it remains open.
    #
    # WHY "Last visit" AND NOT "Last successful" (Rule 7). The agreement is real, but the reference is a single
    # theoretical dataset, not experiment, and the one outlier is unexplained. A "Last successful" date should
    # wait until the m4 amplitude is understood.
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch c: branch a)'s system through the RE-STRUCTURED route (CaptureLine + PhotonLine).
    #
    # module-DielectronicRecombination-inc-FS-resolved.jl is exercised instead of the Pathway/Passage/Resonance
    # code. The point is EQUIVALENCE: the "Total Auger rates, radiative rates and resonance strengths" table must
    # come out bit-identical to branch a). Any difference is a defect in the new route, since the physics, the
    # grid and the settings are otherwise the same.
    #
    # calcPhotonSpectrum = true additionally prints the DR satellite spectrum, which the old route cannot produce:
    # it computes the per-final-level photon energies and then discards them. Two identities are checkable there
    # without any reference data, and they are the unit test for setTotalRates():
    #     sum_f A_r(m,f)      must equal the totalPhotonRate column of every capture line sharing that m
    #     sum_f intensity     must equal that capture line's resonance strength
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    printBefore        = false,
                                                    calcPhotonSpectrum = true )
    wc = Atomic.Computation(Atomic.Computation(), name="Df-c: K-LL DR into He-like C4+, CaptureLine route", grid=grid,
                            nuclearModel        = Nuclear.Model(6.),
                            initialConfigs      = [Configuration("1s^2 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^1 2s^2 2p^0"), Configuration("1s^1 2s^1 2p^1"),
                                                   Configuration("1s^1 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^1 2p^0"), Configuration("1s^2 2s^0 2p^1")],
                            processSettings     = drSettings )
    perform(wc)
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch d: branch b)'s Xe53+ literature case through the RE-STRUCTURED route.
    #
    # The regression that matters most: the six resonances matched against Harman in branch b) must stay within
    # 0.5-3.1%, and the m4 outlier must stay at 28% -- neither improved nor worsened by the rewrite. A rewrite
    # that changed the agreement would be changing physics, which it must not.
    grid54     = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    asfSet54   = AsfSettings(AsfSettings(), qedModel=QedPetersburg(), eeInteraction=CoulombBreit(1.0))
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    augerOperator      = CoulombBreit(1.0),
                                                    printBefore        = false,
                                                    calcPhotonSpectrum = true )
    wd = Atomic.Computation(Atomic.Computation(), name="Df-d: K-LL DR of H-like Xe53+, CaptureLine route", grid=grid54,
                            nuclearModel        = Nuclear.Model(54., "Fermi"),
                            initialConfigs      = [Configuration("1s^1 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^0 2s^2 2p^0"), Configuration("1s^0 2s^1 2p^1"),
                                                   Configuration("1s^0 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^0 2p^0"), Configuration("1s^1 2s^1 2p^0"),
                                                   Configuration("1s^1 2s^0 2p^1")],
                            initialAsfSettings  = asfSet54, intermediateAsfSettings = asfSet54, finalAsfSettings = asfSet54,
                            processSettings     = drSettings )
    perform(wd)
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch e: K-LL DR of He-like W72+ (Z=74) -- the Tu et al. benchmark, and the first case with a WIDE
    #     spread of J_f. Structurally identical to branch a): initial 1s^2, intermediate 1s2l2l', final 1s^2 2l,
    #     but at Z = 74 instead of Z = 6.
    #
    # REFERENCE: B. Tu, J. Xiao, Y. Shen, Y. Yang, D. Lu, T. H. Xu, W. X. Li, C. Y. Chen, Y. Fu, B. Wei,
    #   C. Zheng, L. Y. Huang, R. Hutton, X. Wang, K. Yao, Y. Zou, B. H. Zhang, Y. J. Tang,
    #   "KLL dielectronic recombination resonant strengths of He-like up to O-like tungsten ions",
    #   Phys. Plasmas 23, 053301 (2016).      examples/papers/2016.pp-tu-dr-kll.pdf
    #
    # Their Table II lists, for He-like W, the resonance energy, the natural width and the resonance strength of
    # nine strong KLL resonances, computed with a fully relativistic configuration-interaction method in FAC and
    # including both Coulomb and Breit effects. Experimental strengths were measured at the Shanghai EBIT with
    # uncertainties below 11%.
    #
    # THE FORMULA AGREES WITH JAC'S. Their Eq. (1) is
    #     S = (g_d / 2g_i) * (pi^2 hbar^3 / (m_e E_res)) * A_r A_a / (sum A_r + sum A_a)
    # while JAC forms  pi^2/k^2 * A_a * g_d/g_i  with k^2 = 2E in atomic units. Since
    # pi^2/(2E) * g_d/g_i == pi^2/E * g_d/(2g_i), the two are identical -- the factor 2 that is commented out in
    # the source ("factor 2 is not really clear") is absorbed in k^2 = 2E and must stay commented out. This is
    # now settled from the published formula, not only from numerical agreement.
    #
    # WHY THIS BRANCH MATTERS FOR THE RADIATIVE-PREFACTOR FIX. Carbon (branch a) is Auger dominated and xenon
    # (branch b) radiatively dominated, and both were used to establish the correction. Tungsten K-LL spans a
    # much wider range of J_f -- the final levels run over 1s^2 2s_1/2, 1s^2 2p_1/2 and 1s^2 2p_3/2, and the
    # middle states from J = 1/2 to 5/2 -- so the (2J_f+1)/pi error of the old prefactor would show up here as a
    # RESONANCE-DEPENDENT scatter rather than an overall shift. A corrected formula should agree across the whole
    # table; a fitted one should not.
    #
    # FINAL STATES UP TO n = 4. Tu's CI carries 1s(2l)^(p+1)nl' intermediates and 1s^2(2l)^p nl' finals with n up
    # to 7; a first run with n = 2 finals only gave three strong resonances within 1-8% but a total strength 14%
    # high, and energies 0.6% high. The FINAL set is extended here rather than the intermediate one, for a specific
    # reason: Gamma_r = sum_f A_r(m,f) is what the missing final states suppress, and it enters every resonance
    # strength through the branching ratio, whereas the capture side is untouched -- there are still 16 capture
    # lines, so the expensive continuum-orbital work does not grow at all. Extending the intermediate set instead
    # would add whole new resonance series (KLM, KLN) rather than improve the K-LL ones, at large cost.
    # This is therefore a targeted improvement, NOT Tu's full CI, and the comparison remains approximate.
    #
    # ============================== REPORT (05-Aug-2026) ==============================
    #
    # Wall time 1:01, unchanged from the n=2 run -- the capture side dominates the cost and it did not grow.
    #
    # THE n=4 FINAL STATES CHANGE ESSENTIALLY NOTHING, and the reason is instructive. Gamma_r(m=1) moved from
    # 8.3297e+14 to 8.3271e+14 (0.03%), the total strength stayed at 24.6e-20 cm^2 eV, and the three resonances
    # that match Tu kept their ratios to four figures (0.952, 1.083, 0.989). The extra final levels DO appear in
    # the satellite table, but only through CI admixture with 1s^2 2s / 1s^2 2p: a direct 1s2l2l' --> 1s^2 nl
    # transition with n >= 3 is a TWO-electron jump (2l --> 1s together with 2l' --> nl) and a one-body E1
    # operator cannot drive it. Extending the final set could therefore never have repaired Gamma_r.
    #
    # This is a useful NEGATIVE result: it localises the disagreement with Tu to the INTERMEDIATE-state CI. Tu's
    # 1s(2l)^(p+1)nl' configurations mix into the 1s2l2l' levels themselves, changing both the K-LL energies and
    # their wavefunctions -- consistent with the energies here sitting 0.6% high, and with the strengths being
    # sensitive to those energies through the branching ratio. Closing the gap means adding 1s2lnl' to the
    # INTERMEDIATE configurations, which in JAC also generates the KLM/KLN resonance series and multiplies the
    # expensive continuum work; that is a separate, much larger computation.
    #
    # STATUS vs Tu Table II (He-like W, FAC + Breit): three strong resonances matched by J^P,
    #     He_1 (1/2+) 3.099 vs m1  2.950   ratio 0.952
    #     He_3 (1/2-) 4.504 vs m4  4.879   ratio 1.083
    #     He_4 (3/2-) 3.499 vs m9  3.461   ratio 0.989
    # total strength 24.6 against Tu's 21.5 (14% high), energies systematically +0.6%. Dated "Last visit": the
    # agreement is encouraging but the CI is demonstrably smaller than the reference's, so this is not a
    # validation.
    grid74     = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    asfSet74   = AsfSettings(AsfSettings(), qedModel=QedPetersburg(), eeInteraction=CoulombBreit(1.0))
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    augerOperator      = CoulombBreit(1.0),
                                                    printBefore        = false,
                                                    calcPhotonSpectrum = true )
    we = Atomic.Computation(Atomic.Computation(), name="Df-e: K-LL DR into He-like W72+ (Tu et al. 2016)", grid=grid74,
                            nuclearModel        = Nuclear.Model(74., "Fermi"),
                            # All configurations are written over the SAME ordered subshell list, 1s..4f, with explicit
                            # zero occupations (checkConsistentMultiplets requires it).
                            initialConfigs      = [Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0")],
                            intermediateConfigs = [Configuration("1s^1 2s^2 2p^0 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0"), Configuration("1s^1 2s^1 2p^1 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^1 2s^0 2p^2 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0")],
                            finalConfigs        = [Configuration("1s^2 2s^1 2p^0 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0"), Configuration("1s^2 2s^0 2p^1 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^1 3p^0 3d^0 4s^0 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^1 3d^0 4s^0 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^1 4s^0 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^0 4s^1 4p^0 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^0 4s^0 4p^1 4d^0 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^0 4s^0 4p^0 4d^1 4f^0"),
                                                   Configuration("1s^2 2s^0 2p^0 3s^0 3p^0 3d^0 4s^0 4p^0 4d^0 4f^1")],
                            initialAsfSettings  = asfSet74, intermediateAsfSettings = asfSet74, finalAsfSettings = asfSet74,
                            processSettings     = drSettings )
    perform(we)
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch f: Delta-n = 0 DR of Be-like Au75+ (Z=79) with a high-n Rydberg spectator -- the first branch
    #     that exercises the HydrogenicCorrections, and the reason they were ported to the re-structured route.
    #
    # THE SYSTEM, following apps/apps-schippers-dr-belike/job-belike-a.jl (S. Schippers, 11/2023):
    #
    #     initial        1s^2 2s^2                  (Be-like Au75+)
    #     capture        1s^2 2s 2p 19l             (2s -> 2p core excitation + capture into n = 19)
    #     stabilize (a)  --> 1s^2 2s^2 19l          core decay 2p -> 2s;  EXPLICIT below
    #     stabilize (b)  --> 1s^2 2s 2p n'l'        decay of the Rydberg SPECTATOR;  NOT explicit
    #
    # Channel (b) is the point. The captured electron sits at n = 19 and can itself radiate down to any
    # n' < 19; each such step lowers the ion below the autoionization threshold and stabilizes the resonance
    # exactly as the core decay does. Representing all of 1s^2 2s 2p n'l' for n' = 4...18 explicitly is out of
    # reach, so those channels are estimated with hydrogenic rates. Without them Gamma_r(m) is too small, and
    # since the strength carries the branching ratio Gamma_r/(Gamma_a+Gamma_r), every resonance strength is too
    # small with it. Note the regime: the ratio (Rydberg binding)/(2s-2p excitation energy) grows like Z, so it
    # is at HIGH Z that most of the n' range is genuinely bound and the correction is physically meaningful.
    # Z = 79 is not decoration here.
    #
    # WHY n = 19 AND NOT SOMETHING CHEAPER. The intermediate level must lie ABOVE 1s^2 2s^2 + a free electron,
    # i.e. the Rydberg binding must be smaller than the 2s -> 2p excitation energy. For Be-like gold the latter
    # is a few hundred eV while a captured electron seeing charge 75 is bound by 75^2*13.6/n^2 eV, which does not
    # drop below the excitation energy until n ~ 18. The high shell is forced by the physics, not chosen.
    # The Rydberg l values are cut at 19d (the app used up to 19g) purely to keep the run affordable.
    #
    # WHAT IS BEING TESTED, in order of strength:
    #   1. computeHydrogenicRate itself -- ALREADY VALIDATED independently (05-Aug-2026; the check is not yet part
    #      of the committed suite, it belongs in TestFrames as a known-answer test and needs no approved file):
    #      140 of 140 E1-allowed transitions with n <= 8 reproduce a direct numerical evaluation of the radial
    #      dipole integral to better than 0.001%, the Z^4 scaling is exact to 1e-14, and seven tabulated H rates
    #      agree to 0.076%. The Infeld-Hull recursion is correct; nothing here needs to re-establish that.
    #   2. the WIRING -- that the correction is evaluated once per intermediate level m (not once per (i,m) as in
    #      the old route), lands in Gamma_r(m) where the strength is formed, and is reported rather than absorbed.
    #   3. the SUMMATION RANGE, which was wrong before -- see below.
    #
    # THE BUG THIS BRANCH EXISTS TO CATCH. In the old route the correction was identically zero for this very
    # system, by two independent defects: the live loop ran nf = nFinal+1 : nHydrogenic, i.e. UPWARD, while a
    # hydrogenic rate vanishes unless nf < n_Rydberg; and nFinal was read off the whole final basis, which always
    # contains the captured shell itself, so nFinal = 19 rather than the intended 3. Both are corrected in
    # module-DielectronicRecombination-inc-FS-resolved.jl. The correction range here should therefore print as
    # 4 <= n_f <= 18, and the added rates must be NON-ZERO. If the table shows zeros, the fix has regressed.
    #
    # WHAT TO LOOK AT IN THE OUTPUT. The "Hydrogenic corrections to the total photon rate" table prints, per
    # intermediate level, the explicitly computed Gamma_r, the hydrogenic amount added, and their ratio. The last
    # column is the honest one: effectiveZ = 75 and rateScaling = 1.0 are sensitivity knobs, not physics, and a
    # result carried by them is a result that needs a larger explicit final set instead.
    #
    # ============================== REPORT (05-Aug-2026) ==============================
    #
    # RUNS CLEAN, and the correction is non-zero for the first time: the table prints 4 <= n_f <= 18 (the app's
    # nHydrogenic = 22 capped at n_R - 1) and adds 1.3695e12, 1.7905e12 and 2.4452e12 1/s for Rydberg l_R = 0, 1
    # and 2. Those three numbers reproduce an INDEPENDENT summation over the same range, term by term, exactly --
    # so the wiring, the cap and the l-resolution are all doing what they claim. Before the repair every one of
    # them was 0.0.
    #
    # THE PHYSICS COMES OUT RIGHT WHERE IT CAN BE CHECKED. Two resonance groups appear, at electron energies
    # 12.6-16.4 eV and ~2000 eV. Adding back the n=19 binding of an electron seeing charge 75 (75^2*13.6/361 =
    # 212 eV) gives 2s -> 2p excitation energies of ~225 eV and ~2212 eV, i.e. the 2p_1/2 and 2p_3/2 branches
    # separated by the ~2 keV fine-structure splitting expected for Be-like gold. The explicit Gamma_r is
    # 2.82e9 1/s and is constant to four digits across all thirteen levels of the low group regardless of J or
    # l_R -- exactly the signature of a pure SPECTATOR decay, where the 2p -> 2s core rate cannot depend on how
    # the Rydberg electron is coupled.
    #
    # THE CORRECTION SUPPLIES 99.8-100% OF Gamma_r, and the built-in warning fires. This is not a malfunction, it
    # is the model reporting on itself: the entire Rydberg-decay channel is unrepresented here (the finals are
    # 1s^2 2s^2 nl, whereas the spectator decay ends in 1s^2 2s 2p n'l'), so there is nothing for the estimate to
    # supplement -- it IS the answer. Reducing that fraction would mean carrying 1s^2 2s 2p n'l' finals for
    # n' = 4...18 explicitly, which is precisely what the correction exists to avoid. Read the strengths of this
    # branch as model-dominated.
    #
    # THE RE-AUTOIONIZATION CUTOFF, added 05-Aug-2026 and exercised here for the first time.
    #
    # Nothing prevents the Rydberg electron from falling to a lower n. It does, the photon is emitted, and the
    # rate is real. The issue is a different one: the strength formula S = C * Gamma_r/(Gamma_a + Gamma_r) treats
    # EVERY radiative decay as TERMINATING the process, and that holds only if the state one lands in lies below
    # the autoionization threshold. The core stays excited during a spectator decay, so after 19l -> n'l' the ion
    # is 1s^2 2s 2p n'l', sitting at  E_electron + B(19) - B(n')  above 1s^2 2s^2 + e-. Where that is still
    # positive the ion can Auger straight back: the photon was emitted and the electron then left again, so
    # counting the rate in Gamma_r overestimates the DR yield. The sum is therefore cut at
    #
    #     n' <= n^(bound) = floor( 1 / sqrt( 2*E_electron/Zeff^2 + 1/n_R^2 ) ),      B(n) = Zeff^2/(2n^2).
    #
    # THE CUTOFF IS SELF-CONSISTENT WITH THE CAPTURE STEP, and this run demonstrates it. For the 2p_1/2 group
    # n^(bound) comes out as 18 = n_R - 1 and EXACTLY ZERO rate is removed -- every level prints dropped = 0.0.
    # That is not luck: binding(19) = 212 eV < 225 eV is precisely why n = 19 is the lowest shell this resonance
    # can be captured into, so the same inequality that opens the capture channel closes the bound-state question.
    # For the 2p_3/2 group, built on a 2211 eV core excitation and sitting 2000 eV above threshold, n^(bound) = 5:
    # falling 19 -> 10 sheds only 765 eV of the 2211 eV needed, so those steps do not stabilize.
    #
    # WHAT IT CHANGES, measured:
    #     hydrogenic rate removed, summed over all levels      35.9 %
    #     2p_1/2 group  (13 resonances), summed strength       +0.00 %   -- nothing cut, as predicted
    #     2p_3/2 group  (30 resonances), summed strength       -2.26 %
    #     total DR strength                                    -0.03 %
    #     largest single change: level 20                      x0.593
    #
    # So the cut is large in Gamma_r (a factor ~2 for the 2p_3/2 group) and small in the strengths, because
    # Gamma_a << Gamma_r for nearly all of those levels and the branching ratio is already saturated at ~1. It
    # bites only where the two widths are comparable -- level 20 has Gamma_a = 1.06e12 against Gamma_r = 1.37e12.
    # That is specific to this system: in an AUGER-DOMINATED regime S is proportional to Gamma_r and the full
    # factor would propagate straight into the strengths.
    #
    # STILL A LOWER BOUND, and the table says so. A state landing above threshold can radiate a SECOND time --
    # the core 2p -> 2s, or the spectator falling further -- and stabilize then. Resolving that needs a radiative
    # cascade with branching at every step. The dropped rate is printed per level, so the uncut upper bound is
    # always recoverable by adding it back; the truth lies between the two.
    #
    # (Related, and harmless only by accident: the model also ignores the Pauli principle. A hydrogenic 19p -> 1s
    # term would be included if the range reached n' = 1, although 1s^2 2s 2p 1s does not exist. Here n^(final) = 3
    # keeps the sum clear of it.)
    #
    # DATED "Last visit", NOT "Last successful". The machinery is verified -- summation range, cap, per-l
    # resolution, once-per-m evaluation, reporting, and the underlying rate to 0.001% against an independent
    # implementation. The PHYSICS is not: there is no literature comparison for this system in the repository,
    # the strengths are dominated by an adjustable model, and the bound-state limitation above is known and
    # uncorrected. A "Last successful" date here would claim more than has been shown.
    grid79     = Radial.Grid(Radial.Grid(false), rnt = 6.0e-6, h = 5.0e-2, hp = 6.0e-3, rbox = 10.0)
    asfSet79   = AsfSettings(AsfSettings(), eeInteraction=CoulombBreit(1.0))
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    augerOperator      = CoulombBreit(1.0),
                                                    printBefore        = false,
                                                    calcPhotonSpectrum = true,
                                                    # nHydrogenic = 22 is the app's own value; it is capped
                                                    # internally at n_Rydberg - 1 = 18. effectiveZ = 75 is the
                                                    # charge a fifth electron sees outside the Be-like core.
                                                    corrections        = DielectronicRecombination.AbstractCorrections[
                                                                            DielectronicRecombination.HydrogenicCorrections(22, 75., 1.0)] )
    wf = Atomic.Computation(Atomic.Computation(), name="Df-f: Delta-n=0 DR of Be-like Au75+ with n=19 spectator", grid=grid79,
                            nuclearModel        = Nuclear.Model(79., "Fermi"),
                            # Same ordered subshell list everywhere, with explicit zero occupations.
                            initialConfigs      = [Configuration("1s^2 2s^2 2p^0 3s^0 3p^0 3d^0 19s^0 19p^0 19d^0")],
                            intermediateConfigs = [Configuration("1s^2 2s^1 2p^1 3s^0 3p^0 3d^0 19s^1 19p^0 19d^0"),
                                                   Configuration("1s^2 2s^1 2p^1 3s^0 3p^0 3d^0 19s^0 19p^1 19d^0"),
                                                   Configuration("1s^2 2s^1 2p^1 3s^0 3p^0 3d^0 19s^0 19p^0 19d^1")],
                            # The 1s^2 2s^2 19l finals carry the core-decay channel (a). The 1s^2 2s^2 3l finals are
                            # NOT reachable from 1s^2 2s 2p 19l by a one-body E1 -- that would be a two-electron jump --
                            # but they are what fixes n^(final) = 3, i.e. they declare which low-n shells are treated
                            # explicitly and hence where the hydrogenic estimate is allowed to start.
                            finalConfigs        = [Configuration("1s^2 2s^2 2p^0 3s^1 3p^0 3d^0 19s^0 19p^0 19d^0"),
                                                   Configuration("1s^2 2s^2 2p^0 3s^0 3p^1 3d^0 19s^0 19p^0 19d^0"),
                                                   Configuration("1s^2 2s^2 2p^0 3s^0 3p^0 3d^1 19s^0 19p^0 19d^0"),
                                                   Configuration("1s^2 2s^2 2p^0 3s^0 3p^0 3d^0 19s^1 19p^0 19d^0"),
                                                   Configuration("1s^2 2s^2 2p^0 3s^0 3p^0 3d^0 19s^0 19p^1 19d^0"),
                                                   Configuration("1s^2 2s^2 2p^0 3s^0 3p^0 3d^0 19s^0 19p^0 19d^1")],
                            initialAsfSettings  = asfSet79, intermediateAsfSettings = asfSet79, finalAsfSettings = asfSet79,
                            processSettings     = drSettings )
    perform(wf)
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch g: DR of H-like C5+ (Z=6) into He-like C4+ with the core excited to 2p and the electron captured
    #     into a RYDBERG shell -- the validation vehicle for the RydbergTailCorrection's n-scaling.
    #
    # WHAT THIS BRANCH IS FOR. The whole Rydberg-tail extrapolation rests on ONE scaling assumption,
    #
    #     A_a(n,l)  ~  n^(-3),
    #
    # which follows from the normalization of a Rydberg orbital at the core -- the capture amplitude samples the
    # wavefunction where the core sits, and |psi_nl(r->0)|^2 goes as n^(-3). That is textbook, but it is an
    # assumption about THIS code until checked against it, and the check is free in any calculation carrying TWO
    # Rydberg shells: the ratio of the weighted capture strengths W(n,l) = sum_m A_a(m)(2J_m+1) of two adjacent
    # shells at the same l measures the exponent directly. So n = 5 AND n = 6 are computed explicitly here.
    #
    # WHY AN H-LIKE CORE, after a He-like one failed. The first version of this branch used He-like C4+ with the
    # core excited to 1s2p. It could not do the job, for a reason worth recording: the He-like core's 3P - 1P
    # splitting is about 4 eV, while the 5 -> 6 Rydberg spacing for an electron seeing charge 4 is only
    # 16*13.6*(1/25 - 1/36) = 2.66 eV. The two Rydberg shells therefore INTERLEAVE -- the run showed four clusters
    # of capture energies, 3P.5l at 294.4, 3P.6l at 297.1, 1P.5l at 298.5 and 1P.6l at 301.1 eV -- so a W(n,l)
    # summed over levels of a given (n,l) mixes two different core terms whose relative weights differ between the
    # shells. The measured "exponent" then reports that difference and not the n-scaling; it came out 1.92 with a
    # per-l scatter from 0.20 to 3.08, which is a property of the chosen system, not of the code.
    #
    # An H-like core has no such structure: the 2p_1/2 - 2p_3/2 splitting of C is alpha^2 Z^4/32 = 0.06 eV, forty
    # times SMALLER than the Rydberg spacing, so the two shells are cleanly separated and each is a single series.
    # THE GENERAL RULE, which the module now enforces by counting levels: an n-scaling can only be measured within
    # a series whose core is fixed, so the core structure must be small against the Rydberg spacing.
    #
    # THE SYSTEM.
    #     initial        1s                         (H-like C5+)
    #     capture        2p nl,  n = 5, 6           (1s -> 2p core excitation + capture into nl)
    #     stabilize      --> 1s nl                  (core decay 2p -> 1s)
    #
    # The 1s -> 2p excitation costs about 367 eV (Lyman-alpha of C VI) while an electron seeing charge 5 is bound
    # by 25*13.6/n^2 eV (13.6 eV at n = 5), so the resonances sit near 353 eV and are comfortably autoionizing.
    # Two electrons only, so this is also cheaper than the He-like version it replaces.
    #
    # NOTE ON effectiveZ. The internal default is Z - n^(core), which for this configuration list resolves to 4;
    # the Rydberg electron outside a single 2p core actually sees 5. It is therefore given explicitly -- a good
    # illustration that the derived default is a convenience and not always right.
    #
    # WHAT THE OUTPUT SHOULD SHOW, in order of importance:
    #   1. THE MEASURED EXPONENT, per l, with the number of levels entering each W. Equal level counts are the
    #      precondition for the comparison to mean anything; unequal ones are now excluded and flagged.
    #   2. THE l-PROFILE FIT over l = 0..3, its slope, largest residual, and the share of the weighted capture
    #      strength carried by the extrapolated l = 4..8.
    #   3. CONVERGENCE: each shell's share of the tail, and Gamma_a/Gamma_r. While that ratio is >> 1 the
    #      resonances are Auger dominated and the strength per shell falls slowly; below 1 the n^(-3) descent has
    #      begun and the remaining tail is small.
    #   4. THE alpha(T) SPLIT, explicit versus extrapolated -- whether the tail mattered at all.
    #
    # ============================== REPORT (05-Aug-2026) ==============================
    #
    # THE MACHINERY IS VERIFIED. THE n^(-3) LAW IS NOT. Both halves matter; do not read the first as the second.
    #
    # WHAT PASSED.
    #   * The H-like core did what it was chosen for: the level counts now match at every l (4 vs 4, 10 vs 10,
    #     12 vs 12, 12 vs 12), so the two shells really are one series sampled twice. In the earlier He-like
    #     version they did not, and the module now says so instead of averaging over the mismatch.
    #   * Convergence is real: 26 shells, 1070 extrapolated lines, monotonic decline, last shell 0.57 % of the
    #     tail. Gamma_a/Gamma_r falls from 16.4 at n = 6 through 1 near n = 16 to 0.13 at n = 30, i.e. the series
    #     is followed from the Auger-dominated regime into the radiative one where the n^(-3) descent sets in.
    #   * The l-profile fit is well behaved: ln W(l) = -6.640 - 0.318 l, decaying, and the extrapolated l = 4..8
    #     carry only 6.6 % of the weighted capture strength at n0 -- small enough that the fit is a correction
    #     and not the answer.
    #   * alpha(T) splits cleanly, and reports something uncomfortable but true: the tail carries 51-70 % of
    #     alpha^DR between 1e5 and 3e6 K. For THIS system the plasma rate coefficient is dominated by
    #     extrapolation, and the printed warning says exactly that.
    #
    # WHAT FAILED, and it is the central assumption. The measured exponents are
    #       l = 0 :  p =  2.798        l = 1 :  p =  1.979
    #       l = 2 :  p =  1.720        l = 3 :  p = -1.598
    # The last is unphysical -- the weighted capture strength GROWS from n = 5 to n = 6 -- and the sequence
    # degrades monotonically with l. The leading suspicion is the accuracy of the high-l Rydberg orbitals rather
    # than the law: a Rydberg orbital's amplitude near the core goes as r^l, the capture rate samples precisely
    # that region, and by l = 3 the rate is already 15x below l = 1, so it is a small difference of small
    # quantities. Consistent with this, l = 0 -- the largest core amplitude and the most reliable -- lands within
    # 7 % of 3. But this is a SUSPICION and is not established here; n = 5, 6 is also simply not asymptotic.
    #
    # THE CODE DOES NOT ADOPT THE BAD NUMBER. It warns and builds the tail on the assumed p = 3.0, so the
    # extrapolation rests on the textbook law rather than on a measurement that cannot be defended. That is the
    # right default, and it is why this branch is dated "Last visit" and not "Last successful".
    #
    # OPEN, TO BE SETTLED ON REAL SYSTEMS (user's instruction, 05-Aug-2026). Whether p = 3 holds in JAC should be
    # decided where DR is actually being computed for a purpose -- a real high-n case with converged orbitals,
    # higher n0 (where the law is asymptotic), and l restricted to the range whose capture rates are numerically
    # solid. Until then, treat nExponent = 3.0 as an assumption carried on textbook authority, and read the
    # measured value the module prints as a diagnostic of the calculation rather than of the law.
    grid6      = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 45.0)
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles    = [E1],
                                                    gauges        = [UseCoulomb, UseBabushkin],
                                                    printBefore   = false,
                                                    calcRateAlpha = true,
                                                    temperatures  = [1.0e5, 3.0e5, 1.0e6, 3.0e6],
                                                    corrections   = DielectronicRecombination.AbstractCorrections[
                                                        # spectator decay into the un-represented n = 3, 4
                                                        DielectronicRecombination.HydrogenicCorrections(4, 5.0, 1.0),
                                                        # extrapolate n = 6..30 and l = 4..8; nExponent missing -> 3.0
                                                        # is ASSUMED but also MEASURED, since two shells are present
                                                        DielectronicRecombination.RydbergTailCorrection(30, 5.0, 8, missing)] )
    # One ordered subshell list for all three multiplets, with explicit zero occupations.
    wg = Atomic.Computation(Atomic.Computation(), name="Df-g: Rydberg-tail validation, H-like C5+", grid=grid6,
                            nuclearModel        = Nuclear.Model(6.),
                            initialConfigs      = [Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0")],
                            intermediateConfigs = [Configuration("1s^0 2s^0 2p^1 5s^1 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^1 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^1 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^0 5f^1 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^0 5f^0 6s^1 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^0 5f^0 6s^0 6p^1 6d^0 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^1 6f^0"),
                                                   Configuration("1s^0 2s^0 2p^1 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^1")],
                            finalConfigs        = [Configuration("1s^1 2s^0 2p^0 5s^1 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^1 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^1 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^1 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^0 6s^1 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^0 6s^0 6p^1 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^1 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^0 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^1"),
                                                   # the low-n finals fix n^(final) = 2 and hence where the
                                                   # hydrogenic spectator estimate is allowed to start
                                                   Configuration("1s^1 2s^1 2p^0 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0"),
                                                   Configuration("1s^1 2s^0 2p^1 5s^0 5p^0 5d^0 5f^0 6s^0 6p^0 6d^0 6f^0")],
                            processSettings     = drSettings )
    perform(wg)
    #
elseif  true
    # Last visit:  05-Aug-2026
    # --- Branch h: K-LL DR of Li-like C3+ (Z=6) with the initial ion HYPERFINE-RESOLVED -- the first branch that
    #     exercises module-DielectronicRecombination-inc-HF-resolved.jl on the dimension where the physics lives.
    #
    # WHY A Li-LIKE INITIAL ION. The hyperfine route splits the capture resonances by F_i, and F_i is only
    # interesting when the initial level actually has more than one. Branch a's He-like 1s^2 ground state has
    # J_i = 0, so F_i = I is unique and its statistical weight is exactly 1: the F-sum rule still tests the whole
    # F_m/F_f recoupling machinery (50 intermediate hyperfine levels there), but the initial-state dimension is
    # trivial. Li-like 1s^2 2s has J_i = 1/2, so with I = 3/2 there are TWO initial hyperfine levels, F_i = 1 and
    # F_i = 2, and every K-LL resonance appears twice.
    #
    #     initial        1s^2 2s                     (Li-like C3+, J_i = 1/2)
    #     capture        1s 2s^2 2p, 1s 2s 2p^2, 1s 2p^3     (K-LL, one 1s promoted + one captured)
    #     stabilize      --> 1s^2 2s^2, 1s^2 2s 2p
    #
    # A DELIBERATELY UNPHYSICAL mu. The 2s hyperfine splitting of Li-like carbon is ~1.8 GHz per unit nuclear
    # magneton, i.e. some 7 micro-eV -- far below anything visible beside a 200 eV resonance energy, and far below
    # what any DR measurement resolves. mu is therefore set to 200 nuclear magnetons, which is nonsense as nuclear
    # physics but turns the splitting into a few meV so that it can be READ OFF the table and checked. Nothing
    # else in the calculation depends on mu, so this rescales the effect without distorting anything around it.
    #
    # WHAT TO CHECK, in order:
    #   1. THE F-SUM RULE, printed by the route itself: the statistically weighted hyperfine strengths must
    #      reproduce the fine-structure total EXACTLY (identity among recoupling coefficients, not an
    #      approximation). It held to 4.4e-16 for the He-like case; here the F_i weights are 3/8 and 5/8 rather
    #      than 1, so this run tests them as well.
    #   2. THE SPLITTING of the capture energies over F_i, printed as "shift from lowest F_i" in meV. It must
    #      equal the hyperfine splitting of the INITIAL ion with reversed sign, and for J_i = 1/2, I = 3/2 the
    #      interval rule gives E(F=2) - E(F=1) = 2A, so the two resonance groups must be separated by exactly 2A.
    #   3. TWO ROWS, F_i = 1 and F_i = 2, with strengths in the ratio of their statistical weights when the
    #      splitting is small compared with everything else.
    #
    # ============================== REPORT (05-Aug-2026) ==============================
    #
    # ALL THREE CHECKS PASS, and the third is confirmed by three independent routes agreeing to six digits.
    #
    # 1. TWO initial hyperfine levels, F_i = 1 and F_i = 2, as J_i = 1/2 with I = 3/2 requires; 93 intermediate
    #    and 12 final hyperfine levels, recoupled from 16 electronic capture and 22 electronic photon lines.
    #
    # 2. THE F-SUM RULE holds to 6.4e-11:
    #        sum_F (2F_i+1)/((2I+1)(2J_i+1)) * S(hyperfine)   Coulomb 5.29406373e-05
    #        sum   S(fine structure)                          Coulomb 5.29406373e-05
    #    Note this is a stronger test than the He-like case (branch a's system with I = 3/2), where J_i = 0 gives
    #    a single F_i of weight exactly 1; here the two weights are 3/8 and 5/8 and both must be right. The
    #    tolerance is looser than the 4.4e-16 seen there for a PHYSICAL reason, not a numerical one: with a real
    #    initial splitting the hyperfine and fine-structure resonances sit at slightly different energies, so
    #    their strengths are not exactly proportional. 6e-11 is the size of that effect, not of an error.
    #
    # 3. THE SPLITTING of the capture energies over F_i comes out as 4.212978 meV, and
    #        A of the initial 1s^2 2s level, from Hfs.amplitude   ->  2A = 4.212978 meV
    #        the repaired hyperfine matrix, F=2 minus F=1         ->       4.212978 meV
    #        this branch's capture-energy spread                  ->       4.212978 meV
    #    i.e. ratio 1.000000. The interval rule E(F) = A/2 [F(F+1) - I(I+1) - J(J+1)] gives exactly 2A for
    #    J = 1/2, I = 3/2, so the resonance splitting IS the initial hyperfine splitting, with reversed sign.
    #    That chain ties the DR route, the hyperfine representation and the A constant together in one number.
    #
    # 4. The two strengths are equal to four digits (4.0340e-20 and 4.0341e-20 in Coulomb gauge), as they must
    #    be when the splitting is tiny beside every other scale: the F_i dependence of the strength enters only
    #    through the (2F_m+1)/(2F_i+1) weight and the tiny energy shift.
    #
    # DATED "Last visit", not "Last successful": mu = 200 is deliberately unphysical, so nothing here can be
    # compared with a measurement. What is established is that the machinery is internally exact.
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    nmHf   = Nuclear.Model(Nuclear.Model(6.); spinI = AngularJ64(3//2), mu = 200.0, Q = 0.0)
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles            = [E1],
                                                    gauges                = [UseCoulomb, UseBabushkin],
                                                    printBefore           = false,
                                                    calcHyperfineResolved = true )
    wh = Atomic.Computation(Atomic.Computation(), name="Df-h: hyperfine-resolved K-LL DR of Li-like C3+", grid=grid,
                            nuclearModel        = nmHf,
                            initialConfigs      = [Configuration("1s^2 2s")],
                            intermediateConfigs = [Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2"),
                                                   Configuration("1s 2p^3")],
                            finalConfigs        = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p")],
                            processSettings     = drSettings )
    perform(wh)
    #
end
#
setDefaults("print summary: close", "")
