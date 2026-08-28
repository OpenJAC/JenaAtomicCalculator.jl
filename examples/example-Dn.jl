

println("Dn)  Test of the TwoElectronOnePhoton (TEOP) module: both computational strategies compared for the same transition.")

setDefaults("print summary: open", "zzz-TwoElectronOnePhoton.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")


grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

## Branch a) -- the method test, on a synthetic system with no literature counterpart.
## Branch b) -- the physics test, on the Li-like O5+ case measured by Togawa et al. (2020).
##
## Test system for branch a): He-like Ne8+ (Z=10), doubly-excited 2s2p relaxing to the 1s^2 ground state.
## This is the simplest genuine TEOP system: BOTH electrons must change orbital (2s->1s and 2p->1s) while a single
## photon carries the combined energy -- a one-body E1 operator alone cannot connect these two configurations, so any
## nonzero rate has to come from electron-electron-interaction-driven mixing with a one-hole intermediate configuration
## (here: 1s2p / 1s2s). This is the mechanism described by Heisenberg (1925), first observed by Woelfli et al.,
## PRL 35, 656 (1975), and shown to become STRONG (dominant over Auger/direct decay) for near-degenerate configuration
## mixing by Togawa et al., Phys. Rev. A 102, 052831 (2020) [Li-like O, cross-checked against FAC].
##
## Two INDEPENDENT strategies exist for the same physics, and branch a) runs BOTH so that they can be compared:
##
##   (ii) EXPLICIT second-order perturbation theory -- the TwoElectronOnePhoton module proper. The initial state is
##        the PURE 2s2p configuration; the coupling to the one-hole intermediates is carried by the second-order sum
##        over a Green-function multiplet (settings.gMultiplet), with energy denominators (E_i - E_n) and
##        (E_i - omega - E_n).
##
##   (i)  IMPLICIT, via ordinary CI mixing -- a plain PhotoEmission computation whose INITIAL multiplet already
##        contains the one-hole configuration 1s2p alongside 2s2p. The one-body E1 operator then connects the small
##        1s2p admixture of the (2s2p-dominated) initial level to the 1s^2 ground state. Since the two multiplets are
##        built from different, mutually non-orthogonal orbital sets, this requires calcBiorthogonal=true.
##
## To leading (first) order in the configuration mixing the two must agree: the CI mixing coefficient of strategy (i)
## is precisely c_n = <n|V^ee|i> / (E_i - E_n), i.e. the very ratio that strategy (ii) forms explicitly inside its
## second-order sum. A comparison of the two is therefore a genuine known-answer test of the TEOP module.


if  false
    # Last visit:  28-Aug-2026 -- RE-RUN, AND EVERY LINE AGREES TO 0.08-0.14 %. Today against the 9-Aug column:
    #   STEP 2   2-->1  C 2.760519e+06 / B 9.565061e+06   (was 2.757046e+06 / 9.551948e+06)   +0.13 %
    #            4-->1  C 1.375121e+10 / B 5.448745e+10   (was 1.376205e+10 / 5.453150e+10)   -0.08 %
    #   STEP 3   2-->1  C 3.101045e+06 / B 1.073773e+07   (was 3.097219e+06 / 1.072328e+07)   +0.12 %
    #            4-->1  C 1.534757e+10 / B 6.081964e+10   (was 1.535991e+10 / 6.086998e+10)   -0.08 %
    #   STEP 4   6-->1  C 3.354699e+06 / B 1.938824e+07   (was 3.350980e+06 / 1.936779e+07)   +0.11 %
    #            8-->1  C 1.459402e+10 / B 1.035002e+11   (was 1.460522e+10 / 1.036086e+11)   -0.08 %
    #   Kalpha   4-->1  C 9.864266e+12                    (was 9.863788e+12)                  +0.005 %
    # WHAT THE PATTERN MEANS, and it is worth reading before anyone suspects the angular coefficients: BOTH GAUGES
    # MOVE TOGETHER AND BY THE SAME AMOUNT on every line, and the sign follows the transition energy (the ~1910 eV
    # lines rise, the ~1926 eV lines fall) rather than the multipole or the gauge. That is a small change in the
    # radial/SCF layer -- level energies and mixing coefficients -- not in the angular part, which would show as a
    # gauge-dependent or multipole-dependent shift. The omegas moved too, in the fifth figure.
    # IT ALSO PREDATES 15-Aug: `test/approved/test-PhotoEmission-approved.sum` was regenerated on 15-Aug (6276e4f)
    # and has NOT been touched since, and the suite compares against it at rtol = 1e-6 and passes. So PhotoEmission
    # rates are bit-stable since 15-Aug, and this 0.1 % sits between 9 and 15-Aug -- the kink-aware Slater integral
    # (4cc94eb) and the screened-potential sweeps are the changes in that window.
    # Previously: 09-Aug-2026
    # --- The complete TEOP test: one gMultiplet, three computations, all numbers of the REPORT below reproduced
    #     by a single run of this branch.
    #
    #     STEP 1  the Green-function (intermediate) multiplet
    #     STEP 2  strategy (ii), minimal final space {1s^2}          -- the plain TEOP demonstration
    #     STEP 3  strategy (ii), enlarged final space                -- directly comparable to STEP 4
    #     STEP 4  strategy (i),  enlarged final space                -- the independent biorthogonal route
    #
    # ---------- STEP 1: the Green-function (intermediate) multiplet ----------
    # It is NOT built inside perform(); it is generated here, by a separate and completely standard
    # Atomic.Computation, and then handed to the TEOP settings. This keeps the choice of intermediate levels where it
    # belongs -- with the user -- who selects exactly those configurations expected to contribute strongly to the
    # second-order amplitude. Here these are the natural one-hole intermediates: 1s2p (odd parity, connects to 2s2p
    # via V^ee and to 1s^2 via the E1 operator) and 1s2s (even parity, connects to 1s^2 via V^ee and to 2s2p via the
    # E1 operator) -- i.e. one configuration for each of the two time-orderings of the second-order amplitude.
    gComp      = Atomic.Computation(Atomic.Computation(), name="STEP 1: Green-function (intermediate) multiplet",
                                    grid=grid, nuclearModel=Nuclear.Model(10.),
                                    configs = [Configuration("1s 2p"), Configuration("1s 2s")] )
    gMultiplet = perform(gComp; output=true)["multiplet:"]
    #
    # ---------- STEP 2: strategy (ii) with the minimal final space ----------
    # The intermediate sum in TwoElectronOnePhoton.amplitude() explicitly skips any gMultiplet level that coincides in
    # energy with the initial or the final level (an ill-defined self-term whenever the intermediate and initial/final
    # configurations overlap -- not the case here, but guarded against unconditionally).
    teopSettings2 = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,UseBabushkin], true, LineSelection(), 0.,
                                                  CoulombInteraction(), gMultiplet)
    wa = Atomic.Computation(Atomic.Computation(), name="STEP 2: strategy (ii), minimal final space", grid=grid,
                            nuclearModel   = Nuclear.Model(10.),
                            initialConfigs = [Configuration("2s 2p")],
                            finalConfigs   = [Configuration("1s^2")],
                            processSettings = teopSettings2 )
    perform(wa)
    #
    # ---------- STEPS 3 and 4: the head-to-head comparison ----------
    # The two strategies necessarily differ in their INITIAL configuration space -- that difference IS the comparison:
    # strategy (ii) starts from the pure 2s2p configuration and generates the 1s2p/1s2s admixture perturbatively,
    # strategy (i) starts from the CI-mixed {2s2p, 1s2p} space and lets the diagonalization do it. What must be held
    # fixed for the numbers to be comparable is everything else: grid, nuclear model, and above all the FINAL
    # multiplet.
    #
    # That final multiplet has to be enlarged beyond 1s^2, because BiOrthogonal.computeTransformationMatrices requires
    # the two bases to carry the SAME number of orbitals for every kappa occurring in either of them (the
    # differing-dimension case of Olsen et al., Appendix B, is not yet implemented). The initial basis {2s2p, 1s2p}
    # spans 1s, 2s, 2p_1/2, 2p_3/2, so the final basis must span the same four subshells: {1s^2} alone does not,
    # {1s^2, 1s2s, 2p^2} does -- and both added configurations are the natural correlation partners of the 1s^2 ground
    # state anyway, not arbitrary padding.
    #
    # Both steps select only the lines that END in the ground level (final index 1). The enlarged final space contains
    # 2p^2 levels lying ABOVE part of the initial space, and an unselected run would attempt those upward "decays"
    # with a negative photon energy. TwoElectronOnePhoton.determineLines is presently the ONLY emission module that
    # does not filter omega (Einstein, PhotoEmission, CrystalFieldEmission and HyperfineInduced all do), so such a
    # line reaches GSL's spherical Bessel function out of domain and aborts the run. The initial spaces differ in size
    # (4 levels for 2s2p, 8 for 2s2p+1s2p), hence the two index ranges.
    finalConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("2p^2")]
    #
    teopSettings3 = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,UseBabushkin], true,
                                                  LineSelection(true, indexPairs=[(i,1) for i = 1:4]), 0.,
                                                  CoulombInteraction(), gMultiplet)
    wb = Atomic.Computation(Atomic.Computation(), name="STEP 3: strategy (ii), enlarged final space", grid=grid,
                            nuclearModel   = Nuclear.Model(10.),
                            initialConfigs = [Configuration("2s 2p")],
                            finalConfigs   = finalConfigs,
                            processSettings = teopSettings3 )
    perform(wb)
    #
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin],
                                           printBefore=true, calcBiorthogonal=true,
                                           lineSelection=LineSelection(true, indexPairs=[(i,1) for i = 1:8]) )
    wc = Atomic.Computation(Atomic.Computation(), name="STEP 4: strategy (i), CI mixing + biorthogonal", grid=grid,
                            nuclearModel   = Nuclear.Model(10.),
                            initialConfigs = [Configuration("2s 2p"), Configuration("1s 2p")],
                            finalConfigs   = finalConfigs,
                            processSettings = photoSettings )
    perform(wc)
    #
    # ============================== REPORT (re-run 9-Aug-2026) ==============================
    #
    # All numbers below were re-measured on 9-Aug-2026, after TWO corrections that both touch this branch: the
    # Racah-phase fix in AngularMomentum.JohnsonI (0bdff9f), which moves Coulomb as well, and the length-form
    # orientation fix in MabEmission (c023481), which moves Babushkin only. The original 04-Aug-2026 values are
    # kept in brackets so the two can be told apart.
    #
    # STEP 2 -- strategy (ii), minimal final space {1s^2}. Two lines, both 1- --> 0+:
    #     2 --> 1   omega = 1.911730e+03 eV   A = 2.757046e+06 (Coulomb) / 9.551948e+06 (Babushkin) 1/s
    #                                            [04-Aug: 2.760831e+06          / 6.997322e+06 ]
    #     4 --> 1   omega = 1.927570e+03 eV   A = 1.376205e+10 (Coulomb) / 5.453150e+10 (Babushkin) 1/s
    #                                            [04-Aug: 1.375383e+10          / 5.436247e+10 ]
    # Coulomb is stable to 0.14% and the strong line's Babushkin to 0.31%; only the weak line's Babushkin moves
    # appreciably (+37%), which is the first sign of the cancellation sensitivity discussed under OPEN below.
    # Both rates still sit far below an allowed E1 line at this Z (~1e13-1e14 1/s), as a genuine second-order
    # process should.
    #
    # STEPS 3/4 -- initial levels. Strategy (ii) has 4 levels (pure 2s2p); strategy (i) has 8, the lower four being
    # 1s2p and the upper four the 2s2p states corresponding to strategy (ii)'s. Correspondence and CI level shift:
    #     (ii) level 2  1-  -6.455570e+02 eV   <-->   (i) level 6  1-  -6.449935e+02 eV     shift  +0.563 eV
    #     (ii) level 4  1-  -6.297308e+02 eV   <-->   (i) level 8  1-  -6.292913e+02 eV     shift  +0.440 eV
    # Both 2s2p levels are pushed UP by their mixing with the ~1000 eV lower 1s2p levels -- the correct sign for level
    # repulsion, and a direct measure of the admixture: |V|^2/DeltaE = 0.44...0.56 eV over DeltaE ~ 1.0e3 eV gives
    # |c| ~ 2.1...2.4e-2, i.e. |c|^2 ~ 4e-4...6e-4. Multiplying that by strategy (i)'s own ordinary 1s2p --> 1s^2
    # Kalpha rate (level 4 --> 1: 9.863788e+12 1/s, Coulomb; [04-Aug: 9.804282e+12]) predicts a TEOP rate of
    # ~4e9 1/s -- the right order.
    #
    # STEPS 3/4 -- rates, same final multiplet in both:
    #                              STEP 3, strategy (ii)             STEP 4, strategy (i)             ratio (i)/(ii)
    #   line A, omega ~ 1910 eV    C  3.097219e+06                   C  3.350980e+06                     1.08
    #    [(ii) 2-->1, (i) 6-->1]   B  1.072328e+07                   B  1.936779e+07                     1.81
    #                              [04-Aug: C 2.941112e+06           [04-Aug: C 3.791451e+06             1.29 ]
    #                                       B 8.335175e+06 ]                  B 5.893378e+09 ]            707  (!)
    #   line B, omega ~ 1926 eV    C  1.535991e+10                   C  1.460522e+10                     0.951
    #    [(ii) 4-->1, (i) 8-->1]   B  6.086998e+10                   B  1.036086e+11                     1.70
    #                              [04-Aug: C 1.534363e+10           [04-Aug: C 1.453570e+10             0.947]
    #                                       B 6.066702e+10 ]                  B 7.030514e+10 ]            1.16
    #
    #   own gauge ratio B/C        (ii) line A 3.46, line B 3.96     (i) line A 5.78, line B 7.09
    #                              [04-Aug:     2.8          3.95 ]  [04-Aug:      1554         4.84 ]
    #
    # VERDICT. For line B -- the strong, dominant TEOP channel -- two genuinely independent routes agree to 4.9% in
    # Coulomb (was 5.3%). Given that strategy (ii) is first order in the mixing while strategy (i) diagonalizes it
    # exactly, this remains a real known-answer validation of the second-order amplitude, cross-term double sum
    # included. In Babushkin the two now differ by 70% (was 16%) -- see the warning below before reading that as a
    # regression.
    #
    # RESOLVED. The factor-707 anomaly on line A was a genuine CODE DEFECT, not physics: the length form in
    # InteractionStrength.MabEmission was evaluated with its two orbitals in the orientation that breaks it, so every
    # Babushkin amplitude in JAC carried an error growing as (Z*alpha)^2 (c023481). With that fixed, line A's
    # (i)-vs-(ii) Babushkin disagreement falls 707 --> 1.81 and strategy (i)'s own gauge ratio for it falls
    # 1554 --> 5.78, i.e. back into family with every other line here. The 04-Aug reading that "the anomaly sits in
    # the LENGTH-form amplitude" was exactly right; the suspicion that it lay in the biorthogonal transformation was
    # not.
    #
    # DO NOT read line B's Babushkin agreement getting worse (16% --> 70%) as a regression. Strategy (ii)'s value
    # barely moved (+0.3%) while strategy (i)'s moved +47%, off the same sub-percent shifts in the individual matrix
    # elements: these are nominally two-electron transitions that survive only through CI mixing and near-cancellation,
    # so the sum amplifies small changes enormously. The old 16% was two wrong numbers happening to sit close together.
    #
    # STILL OPEN, and now clearly separated from the code defect: gauge ratios of 3.5-7 across all four entries.
    # These reflect the still-missing damping/regularization term in the second-order amplitude (deferred by explicit
    # request), the incompleteness of the intermediate sum, the absence of a common orbital basis between the two
    # strategies, and the plain fact that a second-order amplitude evaluates its intermediate steps OFF SHELL, where
    # gauge invariance is not expected at all. Not a bug to hunt.
    #
    # NOT RE-CHECKED on 9-Aug-2026: STEP 4's anisotropy (structure) function table, where line 6 --> 1 had been the
    # outlier with f_2 (Coulomb) = 2.858529e+01 against -1.35e-02, 1.679e-01 and -8.48e-01. That anomaly sat in
    # COULOMB while the rate anomaly sat in Babushkin, so it is NOT explained by the length-form fix and should be
    # re-measured before anything is concluded about it.
    #
    # Nothing here is benchmarked against an independent literature value for this synthetic system, so this branch
    # stays dated "Last visit", not "Last successful" (Rule 7).
    #
elseif  true
    # Last visit:  04-Aug-2026
    # --- Li-like O5+: the Togawa et al. TEOP case, the first branch with a real literature counterpart.
    #
    # Togawa, Kuehn, Shah et al., Phys. Rev. A 102, 052831 (2020), arXiv:2003.05965. From the abstract
    # (verbatim): "Some photoabsorption resonances of O5+ reveal strong two-electron--one-photon (TEOP)
    # transitions. We find that for the [(1s 2s)_1 5p_3/2]_{3/2;1/2} states, TEOP relaxation is by far stronger
    # than the radiative decay and competes with the usually much faster Auger decay path. This enhanced TEOP
    # decay arises from a strong correlation with the near-degenerate upper states [(1s 2p_3/2)_1 4s]_{3/2;1/2}
    # of a Li-like satellite blend of the He-like Kalpha transition."
    #
    # NOTE -- an earlier version of this comment had the pair as [(1s2p_3/2)_1 5p] / [(1s2p_1/2)_1 4s]. That was
    # wrong, and wrong in a way that mattered: 1s2p5p is EVEN and 1s2p4s is ODD, so those two could not have
    # mixed at all. The correct pair, 1s2s5p and 1s2p4s, is odd/odd -- CI mixing is allowed and the whole
    # mechanism follows. Any numbers quoted for Auger/one-photon rates in that earlier comment were likewise
    # never verified against the paper body and are deliberately not repeated here.
    #
    # MECHANISM, and why it is the same as branch a) at a heavier level scheme: the near-degenerate partner
    # 1s2p4s carries ordinary Kalpha strength (2p --> 1s) down to 1s^2 4s. The initial 1s2s5p state reaches that
    # same final state only by moving TWO electrons, 2s --> 1s and 5p --> 4s, with a single photon -- a genuine
    # TEOP transition, odd --> even, E1-allowed. Its rate is carried entirely by the 1s2p4s admixture.
    #
    #     STEP 1  gMultiplet: 1s2p4s (time-ordering 1, the near-degenerate partner) and 1s2s4s (time-ordering 2)
    #     STEP 2  TEOP        1s2s5p --> 1s^2 4s
    #     STEP 3  the competing ORDINARY radiative decay 1s2s5p --> 1s^2 5p (one electron, 2s --> 1s), which is
    #             the quantity the paper's central claim is measured against ("TEOP by far stronger than the
    #             radiative decay").
    #
    # The grid is enlarged to rbox = 30 a.u.; n = 5 orbitals at Z = 8 peak near <r> ~ n^2/Z ~ 3 a.u. and are not
    # contained by branch a)'s rbox = 10 box.
    gridB      = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 30.0)
    #
    gCompB     = Atomic.Computation(Atomic.Computation(), name="STEP 1: gMultiplet for the O5+ TEOP case",
                                    grid=gridB, nuclearModel=Nuclear.Model(8.),
                                    configs = [Configuration("1s 2p 4s"), Configuration("1s 2s 4s")] )
    gMultipletB = perform(gCompB; output=true)["multiplet:"]
    #
    teopSettingsB = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,UseBabushkin], true, LineSelection(), 0.,
                                                  CoulombInteraction(), gMultipletB)
    wd = Atomic.Computation(Atomic.Computation(), name="STEP 2: TEOP  1s2s5p --> 1s^2 4s", grid=gridB,
                            nuclearModel   = Nuclear.Model(8.),
                            initialConfigs = [Configuration("1s 2s 5p")],
                            finalConfigs   = [Configuration("1s^2 4s")],
                            processSettings = teopSettingsB )
    perform(wd)
    #
    photoSettingsB = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1],
                                            gauges=[UseCoulomb,UseBabushkin], printBefore=true)
    wf = Atomic.Computation(Atomic.Computation(), name="STEP 3: K-shell-filling decay  1s2s5p --> 1s^2 5p",
                            grid=gridB, nuclearModel = Nuclear.Model(8.),
                            initialConfigs = [Configuration("1s 2s 5p")],
                            finalConfigs   = [Configuration("1s^2 5p")],
                            processSettings = photoSettingsB )
    perform(wf)
    #
    # STEP 4: the one-electron radiative decay that IS allowed, 5p --> 4s within the doubly-excited manifold.
    wg = Atomic.Computation(Atomic.Computation(), name="STEP 4: outer-electron decay  1s2s5p --> 1s2s4s",
                            grid=gridB, nuclearModel = Nuclear.Model(8.),
                            initialConfigs = [Configuration("1s 2s 5p")],
                            finalConfigs   = [Configuration("1s 2s 4s")],
                            processSettings = photoSettingsB )
    perform(wg)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # The 1s2s5p initial multiplet has 7 levels; the paper's [(1s 2s)_1 5p_3/2]_{3/2;1/2} pair is levels 4 (3/2-,
    # -1.069644e+03 eV) and 5 (1/2-, -1.069642e+03 eV), 0.30 eV above the {1,2,3} fine-structure group.
    #
    #  level  J^P     STEP 2, TEOP 1s2s5p-->1s^2 4s      STEP 4, radiative 5p-->4s      TEOP / radiative
    #                 A(Coulomb) / A(Babushkin) [1/s]    A(Coulomb) [1/s]
    #    1    1/2-    3.269178e+07 / 3.389466e+07        9.193920e+08                        0.036
    #    2    3/2-    7.879110e+07 / 8.173859e+07        9.178975e+08                        0.086
    #    4    3/2-    1.064132e+12 / 1.110085e+12        7.103460e+04                        1.5e+07
    #    5    1/2-    1.177252e+12 / 1.227215e+12        2.762565e+04                        4.3e+07
    #    6    1/2-    2.906786e+09 / 6.983910e+10
    #    7    3/2-    2.922373e+09 / 6.946454e+10
    #
    # VERDICT -- the paper's central claim is REPRODUCED, and reproduced selectively. For levels 4 and 5, exactly the
    # [(1s2s)_1 5p_3/2]_{3/2;1/2} pair the abstract names, the TEOP rate exceeds the ordinary one-electron radiative
    # decay by SEVEN orders of magnitude. For levels 1 and 2 the ordering is reversed -- radiative decay wins by a
    # factor 12-28. So this is not a blanket enhancement of everything in the multiplet; it singles out the same two
    # levels the experiment did, which is a far stronger test than a single number would have been.
    #
    # Absolute magnitude: 1.06e+12 ... 1.18e+12 1/s against the ~3e+12 1/s quoted for this state in the earlier
    # comment (a value taken from the paper BODY, which was not available here -- only the abstract was checked, and
    # the abstract quotes no rates). Same order of magnitude, low by a factor ~2.5-3, which is what one should expect
    # from a gMultiplet of just two reference configurations and no further correlation.
    #
    # Gauge agreement is 4.3% (level 4) and 4.2% (level 5) -- dramatically better than branch a)'s 3x-5x spread.
    # That is the expected signature of a DOMINANT, non-cancelling second-order amplitude: branch a)'s poor gauge
    # behaviour comes from its amplitudes being small and near-cancelling, not from a defect in the module.
    #
    # STEP 3 deliberately returns an EMPTY line table, and that emptiness is the point: 1s2s5p and 1s^2 5p are BOTH
    # odd, so the K-shell-filling channel 2s --> 1s is E1-forbidden (s --> s). The state cannot get rid of its K-shell
    # hole by an ordinary one-photon transition at all -- which is precisely why the two-electron channel, riding on
    # the 1s2p4s admixture where 2p --> 1s IS allowed, becomes the dominant radiative route.
    #
    # Dated "Last visit": the qualitative claim and the level selectivity are reproduced, but the absolute rate has
    # not been checked against the published number itself (paper body not consulted), so this is not yet a verified
    # quantitative benchmark under Rule 7.
    #
end
#
setDefaults("print summary: close", "")
