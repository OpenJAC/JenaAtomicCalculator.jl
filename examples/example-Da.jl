#
println("Da) Apply & test the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("unit: energy", "eV")   ## setDefaults("unit: rate", "a.u.")
setDefaults("unit: rate", "1/s")
setDefaults("print summary: open", "zzz-radiative.sum")


if  false
    # Last successful:  28-Aug-2026 -- RE-RUN UNCHANGED. Every number recorded below reproduces to all four
    # figures it was quoted to: 4.585787e12 / 4.755319e12 (Kalpha1, Coulomb/Babushkin), 2.305319e12 / 2.390639e12
    # (Kalpha2), at 849.9454 and 849.8406 eV. The three derived checks also hold: intensity ratio 1.9892, gauge
    # gap 3.70 % on BOTH fine-structure components, radiative widths 4.536e-3 eV (Coulomb) / 4.703e-3 eV
    # (Babushkin). Run as one of three PhotoEmission checks on 28-Aug (with Dn and Ba); see Dn for the one place
    # where anything moved at all.
    # Previously: 9-Aug-2026
    # Branch 0: Ne K-alpha, 1s-hole "1s 2s^2 2p^6" -> 2p-hole "1s^2 2s^2 2p^5" (2p fills the K-hole).
    #   Set up per explicit user request/hypothesis: our SCF procedure may not produce sufficiently relaxed
    #   orbitals for shallow valence-hole cases (like the Fe X 3p-hole tested above) for a bi-orthogonal
    #   transformation to show a large effect; a K-hole (1s-hole) should perturb the WHOLE electron cloud's
    #   screening far more, giving substantially larger orbital non-orthogonality to test against.
    #   PRE-CHECK (orbital overlap matrix, diagonal <initial orbital|final orbital>): 1s: 1-<i|f>=3.16e-4,
    #   2s: 1.73e-3, 2p_1/2: 4.00e-3, 2p_3/2: 3.97e-3 -- THREE TO FOUR ORDERS OF MAGNITUDE larger than the Fe
    #   X 3p-hole case (which was ~1e-8). Full transformation matrix Cleft[1s,2s]=0.0147 (vs ~1e-4 for Fe X).
    #   Confirms the hypothesis at the orbital level; the actual PhotoEmission rate effect (with
    #   calcBiorthogonal) was ~4.56% Coulomb-gauge shift, ~38x larger than Fe X's ~0.12% -- see
    #   project_biorthogonal_transformation.md for the full BiOrthogonal validation follow-up.
    #
    #   REPORT (WITHOUT biorthogonal, this run): two lines, 1->1 (2p_3/2 hole, "Kalpha1") and 1->2
    #   (2p_1/2 hole, "Kalpha2"):
    #     Line   Gauge       A [1/s]        Energy [eV]
    #     1->1   Coulomb     4.586e12       849.94
    #     1->1   Babushkin   4.755e12
    #     1->2   Coulomb     2.305e12       849.84
    #     1->2   Babushkin   2.391e12
    #   Three independent literature/physical checks, all consistent: (1) transition energy 849.9 eV vs.
    #   tabulated Ne Kalpha ~848.6-849.8 eV (X-ray Data Booklet) -- agreement <0.2%; (2) intensity ratio
    #   A(1->1)/A(1->2) = 1.99 vs. the textbook Kalpha1:Kalpha2 statistical-weight ratio 2p_3/2:2p_1/2 =
    #   4:2 = 2:1 -- essentially exact; (3) total radiative width (sum of both lines) = 4.54e-3 eV (Coulomb)
    #   / 4.70e-3 eV (Babushkin) vs. the width implied by the well-known Ne K-shell fluorescence yield
    #   omega_K ~= 1.7-1.8% (Bambynek et al., Rev. Mod. Phys. 44, 716 (1972); Krause, J. Phys. Chem. Ref.
    #   Data 8, 307 (1979)) times the tabulated total K-width Gamma_K ~= 0.24 eV, i.e. Gamma_rad ~= 4.1e-3
    #   eV -- agreement ~10-15%, good for an uncorrelated single-configuration DHF calculation. Gauge
    #   disagreement here (3.70% on BOTH lines) is much SMALLER than Fe X branch 1's (~55%, see below) -- not simply
    #   explained by "more relaxation = more gauge disagreement": a strongly-bound, deeply-relaxed K-hole
    #   transition can still be well described by a single configuration, while the shallow Fe X valence
    #   transition is more correlation-sensitive despite its much smaller absolute non-orthogonality.
    #   RE-RUN 9-Aug-2026, after the MabEmission length-form orientation fix (c023481). Coulomb is unchanged
    #   to every digit; the two Babushkin rates moved by +0.5% = (Z*alpha)^2 at Z=10, exactly as that fix
    #   predicts. All three checks above are unaffected: ratio 1.9892 in BOTH gauges, widths 4.536e-3 eV
    #   (Coulomb) / 4.704e-3 eV (Babushkin). The 3.7% gauge gap is therefore NOT the code defect -- it is
    #   genuine correlation/relaxation, and the earlier suspicion that it was "not a module bug" was right.
    #   What the fix DID repair is an inconsistency between the two components: the gauge gap used to read
    #   3.7% on Kalpha1 but 3.2% on Kalpha2, and two fine-structure components of one multiplet must share
    #   it. They now both read 3.70%.
    #   IMPORTANT: run useBiorthogonal=false and =true in SEPARATE Julia sessions, as in the branches above.
    useBiorthogonal = false
    grid = Radial.Grid(true)
    setDefaults("standard grid", grid)
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin],
                                            printBefore=true, calcBiorthogonal=useBiorthogonal)
    comp = Atomic.Computation(Atomic.Computation(), name="Ne K-alpha 1s-hole -> 2p-hole",
              grid=grid, nuclearModel=Nuclear.Model(10.),
              initialConfigs = [Configuration("1s 2s^2 2p^6")],
              finalConfigs   = [Configuration("1s^2 2s^2 2p^5")],
              processSettings = photoSettings );
    perform(comp)
    #
elseif false
    # Last successful:  1-Aug-2026
    # Branch 1: Fe X (Ne-like) [Ne]3s 3p^6 -> [Ne]3s^2 3p^5 ONLY (the simplest possible subset of the
    #   user's 20-year-old test case below) -- WITHOUT vs. WITH the bi-orthogonal transformation
    #   (calcBiorthogonal), explored in full detail per explicit request. This restricted config was chosen
    #   deliberately: the FULL test case (branch below, with an added [Ne]3s^2 3p^4 3d initial config) hits a
    #   genuine, current BiOrthogonal limitation -- see that branch's REPORT -- so this simpler subset is
    #   used to actually exercise and compare the transformation in detail, without deleting or replacing
    #   the original branch.
    #   IMPORTANT: run useBiorthogonal=false and =true in SEPARATE Julia sessions (flip the flag below and
    #   rerun `include(...)` fresh each time) -- an internal caching path makes a same-session before/after
    #   comparison unreliable.
    #
    #   REPORT: this config's [Ne]3s 3p^6 -> [Ne]3s^2 3p^5 transition set is small -- only 2 distinct lines
    #   (1->1 and 1->2), 4 rows counting both gauges; there is no way to reach anywhere near 20 lines from
    #   this particular pair of configurations (noted honestly, not padded). Both settings ran cleanly
    #   (unlike branch 2 below, this config has matching kappa symmetries on both sides -- no BiOrthogonal
    #   orbital-mismatch issue). Einstein A-values [1/s], WITHOUT vs. WITH bi-orthogonal:
    #
    #     Line   Gauge       WITHOUT              WITH                rel. change
    #     1->1   Coulomb     2.264040e+10         2.261340e+10        -0.1193%
    #     1->1   Babushkin   3.504448e+10         3.504464e+10        +0.0005%
    #     1->2   Coulomb     9.515701e+09         9.503589e+09        -0.1273%
    #     1->2   Babushkin   1.514383e+10         1.514386e+10        +0.0002%
    #
    #   Both transitions show the SAME pattern: the Coulomb (velocity) gauge shifts by ~0.12-0.13%, while
    #   the Babushkin (length) gauge barely moves at all (~0.0002-0.0005%, two orders of magnitude smaller).
    #   This is consistent across both independent lines, not a one-off -- physically expected, since the
    #   length-gauge matrix element is dominated by the LARGE-r part of the wavefunction (far less sensitive
    #   to the small orbital non-orthogonality the bi-orthogonal transformation corrects for), while the
    #   velocity-gauge matrix element (involving a radial derivative) is comparatively more sensitive to
    #   exactly this kind of small basis imperfection. This is not, by itself, a proof that the
    #   bi-orthogonal result is "more correct" (both gauges should agree exactly only in the limit of an
    #   exact/complete basis) -- but it is a clean, physically sensible, reproducible signature of what the
    #   transformation actually does to this calculation.
    useBiorthogonal = false
    grid = Radial.Grid(true)
    setDefaults("standard grid", grid)
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin],
                                            printBefore=true, calcBiorthogonal=useBiorthogonal)
    comp = Atomic.Computation(Atomic.Computation(), name="Fe X 3s 3p^6 -> 3s^2 3p^5",
              grid=grid, nuclearModel=Nuclear.Model(26.),
              initialConfigs = [Configuration("[Ne] 3s 3p^6")],
              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")],
              processSettings = photoSettings );
    perform(comp)
    #
elseif false
    # Last successful:  1-Aug-2026
    # Branch 2: Fe X (Ne-like) [Ne]3s 3p^6 / [Ne]3s^2 3p^4 3d -> [Ne]3s^2 3p^5 / [Ne]3s^2 3p^4 3d -- the
    #   user's own long-standing test case, basis for a paper roughly 20 years ago. finalConfigs EXPANDED
    #   this session (per explicit user request) to also include "[Ne] 3s^2 3p^4 3d", closing the orbital
    #   set under single excitations so BOTH bases now contain a kappa=+2 (d_3/2) orbital -- this resolves
    #   the earlier dimensional mismatch that made calcBiorthogonal=true fail outright (see prior REPORT,
    #   now superseded). The genuine 3p^5 <- 3d decays remain the physically meaningful ones (lower final
    #   energy); the "3d <- 3d" pairs are the same nominal configuration on both sides and are expected to
    #   be irrelevant/filtered by energy ordering, included only to give BiOrthogonal a well-defined
    #   same-dimension problem to solve.
    #
    #   REPORT (WITHOUT biorthogonal): 18 distinct E1 lines (36 rows incl. both gauges), all positive
    #   rates, all transition energies in the expected 32-53 eV Fe X M-shell/EUV range. Two groups, as
    #   anticipated: direct 3p^5<-3s3p^6 decays within the low-lying manifold (levels 1,3,4,5, ~33-47 eV --
    #   note levels 3,4,5 are NOT a pure 3p^5-hole 2P term, since this is now a genuinely 2-configuration
    #   (3s3p^6 / 3s^2 3p^4 3d) initial CI mix, so some low-lying 3d character is already present here, not
    #   cleanly separated by energy as the branch-2 setup note above anticipated -- an honest correction to
    #   that expectation) and higher, more clearly 3d-dominated levels (8,9,10,11,... ~50-53 eV) decaying
    #   down into the low manifold. Gauge (Coulomb vs. Babushkin) disagreement ranges widely across lines,
    #   roughly 1-35% -- comparable to, not worse than, branch 1's own ~55% baseline disagreement for the
    #   simplest single-configuration subset, consistent with this being a more correlated (2-config)
    #   description of closely related physics. No unphysical values (no negative rates, no NaN/zero where
    #   a nonzero rate is expected). No external literature table available this session for this specific
    #   2-configuration calculation (the user's own ~20-year-old paper would be the natural reference, not
    #   accessible here) -- dated on the physical-consistency criterion (Rule 7: no zero rates, wrong
    #   units, or clearly wrong magnitudes found), not an independent numerical comparison.
    #   IMPORTANT: run useBiorthogonal=false and =true in SEPARATE Julia sessions, as in branch 1.
    useBiorthogonal = false
    grid = Radial.Grid(true)
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5, h = 1.0e-2, hp = 0., rbox = 5.0)
    setDefaults("standard grid", grid)
    defaultsSettings = PhotoEmission.Settings()
    photoSettings = PhotoEmission.Settings(defaultsSettings, multipoles=[E1], gauges=[UseCoulomb,UseBabushkin], printBefore=true,
                                            calcBiorthogonal=useBiorthogonal)

    comp = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for the spectrum Fe X",
              grid=grid, nuclearModel=Nuclear.Model(26.);
              initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s^2 3p^4 3d")],
              processSettings = photoSettings );
    @show comp
    perform(comp)
    #
elseif false
    # Last successful:  1-Aug-2026
    # Compute transition probabilities among low-lying levels of He-like Si^12+ ions
    # Test: Lifetimes of various helium- and lithium-like ions; cf. Figure 8.1 in section 8.1.a
    #
    #   REPORT: 3 lines, exactly the expected He-like w/y/z triad: 1->1 (1s2s 3S1 -> 1s^2 1S0, M1
    #   "forbidden" line z) A=2.92e7 /s; 3->1 (1s2p 3P1 -> 1s^2 1S0, E1 intercombination line y)
    #   A=1.23e11(C)/1.14e11(B) /s, gauge disagreement ~7.3%; 6->1 (1s2p 1P1 -> 1s^2 1S0, E1 resonance
    #   line w) A=3.91e13(C)/3.98e13(B) /s, gauge disagreement ~1.9%. Hierarchy z << y << w (ratios
    #   ~1:4000:300000) and resonance-line lifetime ~1.0e-14 s (~10 fs) both match well-known He-like
    #   ion systematics for this Z (Si, Z=14) -- forbidden lines strongly suppressed relative to the fully
    #   allowed resonance line, gauge agreement tightest for the strongest (best single-configuration-
    #   converged) line. No external "Figure 8.1" table available this session to check absolute numbers
    #   against; dated on this qualitative-systematics + physical-consistency basis (Rule 7).
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles = [E1,M1], gauges = [UseCoulomb,JAC.UseBabushkin], 
                                       printBefore = true )
    wa = Atomic.Computation(Atomic.Computation(), name="xx",  grid=grid, nuclearModel=Nuclear.Model(14.),
                            initialConfigs = [Configuration("1s 2s"), Configuration("1s 2p")],
                            finalConfigs   = [Configuration("1s^2")], 
                            processSettings= pSettings)
    wb = @time( perform(wa) )
    #
elseif false
    # Last successful:  1-Aug-2026
    # Compute transition probabilities among low-lying levels of Li-like Mg^9+ ions
    # Test: Lifetimes for lithium-like ions; cf. Figure 8.1 in section 8.1.a
    #
    #   REPORT: 15 levels, 39 rows (E1 + M1 lines mixed). All nonzero rates positive and sensibly scaled
    #   (E1 resonance-like lines ~1e13, weaker E1 satellites ~1e9, M1 lines ~1e5-1e7 /s). Levels 6,7 (from
    #   the "1s 2p^2" even configuration) show EXACTLY ZERO M1 rate to level 1 ("1s^2 2s") despite parity
    #   allowing it -- this is the expected LS-coupling M1 selection rule (Delta-L=0, Delta-S=0): only an
    #   L=0,S=1/2 component of 1s2p^2 can couple to the pure 2S term of 1s^2 2s via M1, so a component with
    #   L!=0 or S!=1/2 gives an exact zero, not a bug. Gauge (Coulomb/Babushkin) agreement mostly within
    #   ~5-20% across the E1 lines. No external "Figure 8.1" table available this session; dated on the
    #   physical-consistency + selection-rule basis (Rule 7).
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles = [E1,M1], gauges = [UseCoulomb,JAC.UseBabushkin], 
                                       printBefore = true )
    wa = Atomic.Computation(Atomic.Computation(), name="Lifetimes for lithium-like ions", 
                            grid=grid, nuclearModel=Nuclear.Model(12.), 
                            initialConfigs = [Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                            processSettings = pSettings );
    wb = @time( perform(wa) )
    #
elseif true
    # Last visit:  18-Aug-2026
    # Last successful:  18-Aug-2026
    # THE QUALITY INDICATORS, on a case chosen so that the two of them disagree.
    #
    # Every rate this module prints carries two cheap indicators of how far it may be trusted, and they are
    # cheap because both are by-products of numbers already computed:
    #
    #   * Cowan's CANCELLATION FACTOR, CF = |SUM_rs c_f(r) M(r,s) c_i(s)| / SUM_rs |c_f(r) M(r,s) c_i(s)| --
    #     the coherent sum that IS the amplitude, over the sum of the magnitudes of its terms.  CF = 1 means
    #     nothing cancels; CF near 0 means the amplitude is a small residue of large opposing contributions, so
    #     a small error in any mixing coefficient moves it by a large factor.  Cowan's rule of thumb: below
    #     0.10 unreliable, below 0.05 no information.  It falls out of the very sum that forms the amplitude.
    #
    #   * The BABUSHKIN/COULOMB RATE RATIO.  Length and velocity forms agree exactly only for exact wave
    #     functions, so a ratio far from 1 says the description is incomplete.
    #
    # WHY BOTH, and this branch is the demonstration.  They fail independently: CF is a property of ONE gauge's
    # amplitude and says nothing about the orbitals; the ratio is a property of the orbitals and says nothing
    # about cancellation.  Be-like Fe shows one of each in a single run.
    #
    #   REPORT (Fe XXIII, 1s^2 2s 2p -> 1s^2 2s^2 + 1s^2 2p^2, so the even J=0 levels are genuine 2s^2/2p^2
    #   mixtures and the sums really do cancel):
    #
    #     transition                    CF Coulomb   CF Babushkin   B/C ratio   verdict
    #     2 -> 1  (1- -> 0+)  E1          0.0672        0.1911        4.4895    CF low
    #     4 -> 1  (1- -> 0+)  E1          1.0000        1.0000        1.9169    gauges 50%+
    #
    #   Line 2 is the INTERCOMBINATION line, 2s2p 3P_1 -> 2s^2 1S_0: it is weak precisely because it is a small
    #   residue -- CF = 0.067 says that 93 % of the amplitude has cancelled -- and the two gauges then differ by
    #   a factor 4.5, which is what such a residue does.
    #
    #   Line 4 is the RESONANCE line, 2s2p 1P_1 -> 2s^2 1S_0: CF = 1.0000 in BOTH gauges, i.e. nothing cancels
    #   at all, and the rates still differ by 92 %.  Nothing about the cancellation factor could have predicted
    #   that; it is the orbitals, not the sum.  This is why one indicator is not enough.
    #
    #   NOT claimed: that either flag is an error bar.  A high CF and a ratio near 1 do not prove a rate right;
    #   they only remove two known ways for it to be wrong.
    pSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1, M1, E2], printBefore=false)
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    wa = Atomic.Computation(Atomic.Computation(), name="Be-like Fe: quality indicators", grid=grid,
                            nuclearModel = Nuclear.Model(26.),
                            initialConfigs = [Configuration("1s^2 2s 2p")],
                            finalConfigs   = [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")],
                            processSettings = pSettings );
    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
