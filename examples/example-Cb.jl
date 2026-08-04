
println("Cb) Apply & test the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
println("    Branches follow Andersson & Jonsson, CPC 178 (2008) 156-170 (examples/papers/")
println("    2008.cpc-andersson-jonsson-hfs-zeemann.pdf); see project memory project_zeeman_hfs_bugs.md.")

Defaults.setDefaults("unit: energy", "Hz")
setDefaults("print summary: open", "zzz-Hfs.sum")

if  false
    # Last successful:  26-Jul-2026
    # Branch a: HFS A parameter for hydrogen 1s (the 21 cm line). Proton: I=1/2, mu=2.7928 nmu, Q=0,
    #   rms radius=0.8797 fm. "uniform" nuclear model, not "Fermi": for Z=1, JAC's 2-parameter Fermi
    #   charge distribution (fixed skin thickness) cannot represent an rms radius below ~1.86 fm --
    #   far above hydrogen's physical 0.88 fm -- and raises an explicit error (module-Nuclear.jl, fixed
    #   25-Jul-2026; see memory project_isotope_shift_ris3ris4.md). NuclearField avoids DFS
    #   self-interaction error for single-electron systems.
    #   TWO real bugs found and fixed here (26-Jul-2026):
    #   (1) InteractionStrength.hfs_tM1 was missing the alpha (fine-structure constant) prefactor
    #   required by Andersson & Jonsson (2008), CPC 178, Eq. (49) -- confirmed by direct re-reading of
    #   the paper (p. 161): <n_a kappa_a||t^(1)||n_b kappa_b> = -alpha(kappa_a+kappa_b)
    #   <-kappa_a||C^(1)||kappa_b>[r^-2]. This alone brought A(1s) from ~2.76e5 MHz down to 2012.65 MHz
    #   -- much better, but still ~1.42x (~sqrt(2)) too large.
    #   (2) Hfs.amplitude and Hfs.computeInteractionMatrix used coeff.T directly from
    #   SpinAngular.computeCoefficientsNonScalar, which (for rank>0 one-particle operators) applies an
    #   internal "GRASP-like" sqrt(2*j_a+1) factor that the paper's own Eq. (48) reduced-matrix-element
    #   sum does NOT want (the "pure" coefficient). By contrast, IsotopeShift.amplitude's rank-0 case
    #   uses computeCoefficientsScalar, where the equivalent conversion step is deliberately commented
    #   out in the source, and IsotopeShift.amplitude re-applies sqrt(2j_a+1) itself externally --
    #   confirming this is a real, load-bearing distinction between the two SpinAngular code paths, not
    #   a coincidence. Fixed by dividing each coeff.T by sqrt(2*j_a+1) in both Hfs functions.
    #   Verified: A(1s) = 1423.157 MHz vs the famous 1420.405751 MHz (21 cm line) -- 0.19% agreement.
    #   Cross-checked independently with H(2p): A(2p_1/2)=43.359 MHz, A(2p_3/2)=8.671 MHz -- internally
    #   consistent (ratio exactly 0.2), though no independently pre-verified literature target was used
    #   for this second case (see branch b for that instead, a genuinely independent system).
    wa = Atomic.Computation(Atomic.Computation(), name="Cb-a-H1s", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8797, AngularJ64(1//2), 2.7928, 0.0, 0.0),
                            configs=[Configuration("1s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.NuclearField()),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  26-Jul-2026
    # Branch b: HFS A parameter for the Na (Z=11) [Ne]3s ^2S_1/2 ground state -- the paper's own second
    #   test case (Section 6: "excellent agreement with the Breit-Rabi formula", no explicit number given
    #   in the paper text). Na-23: I=3/2, mu=+2.2176 n.m., Q=+0.10 barn (Q plays no role for calcM1), rms
    #   radius estimated 2.98 fm (not independently verified). Single configuration, no core-polarization
    #   correlation.
    #   Target: the well-known Na ground-state hyperfine splitting frequency 1771.626 MHz (F=2 to F=1)
    #   corresponds to A = splitting/2 = 885.813 MHz (via Delta E(F) = A/2[F(F+1)-I(I+1)-J(J+1)], I=3/2,
    #   J=1/2).
    #   Verified: A(3s) = 871.219 MHz vs 885.8 MHz -- 1.65% low. A larger residual than H(1s)'s 0.19% is
    #   expected and not concerning: Na's 3s hyperfine constant is well known to receive a significant
    #   contribution from CORE POLARIZATION (the valence 3s electron polarizing the closed
    #   1s^2 2s^2 2p^6 core), which this single-configuration calculation does not include at all. A
    #   second, independent confirmation (different Z, different electron count, same underlying fix)
    #   landing within ~2% using only a bare single configuration is a good result.
    wa = Atomic.Computation(Atomic.Computation(), name="Cb-b-Na3s", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., "uniform", 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  29-Jul-2026
    # Branch c (AL-Field / Breit revisit): reuses branch b's Na [Ne] 3s system, now that the AL-Field SCF
    #   (Basics.ALField()) has been root-cause-fixed and promoted to the standard implementation (see
    #   project_df_al_kink_bug.md) and Breit has been checked to run cleanly through the CI step for other
    #   properties this session. Four settings combinations, same nuclear model as branch b:
    #   DFS+Coulomb (reproduces branch b's baseline, A(3s)=871.219 MHz), AL+Coulomb, DFS+Breit, AL+Breit.
    #   STRUCTURAL PREDICTION for Breit (stated before running, then checked, not asserted after the fact):
    #   Na [Ne] 3s is a SINGLE-CSF system -- its CI "matrix" is a trivial 1x1 block, so the mixing
    #   coefficient is exactly 1.0 regardless of eeInteractionCI, and neither module-Hfs.jl nor its
    #   operators reference Breit/eeInteraction at all (confirmed by grep) -- Breit should shift only the
    #   total energy, leaving A(3s) EXACTLY unchanged.
    #   RESULT (1), AL-Field: A(3s) = 640.967 MHz (AL) vs 871.219 MHz (DFS) -- a LARGE, ~26% DECREASE, moving
    #   FURTHER from the 885.8 MHz literature target, not closer. This is a real, somewhat unexpected finding
    #   (AL-Field's genuine non-local exchange was hoped to partially compensate for the missing core-
    #   polarization contribution, but instead makes agreement markedly worse for this valence-orbital-shape-
    #   sensitive property) -- not investigated further here (would need comparing AL vs DFS 3s orbital
    #   shapes directly near the nucleus, out of scope for this revisit). A concrete, reportable answer to
    #   question (1), just not the hoped-for direction.
    #   RESULT (2), Breit: A(3s) = 871.2193 MHz (DFS+Breit) and 640.9670 MHz (AL+Breit) -- EXACTLY identical
    #   (to the printed precision) to the corresponding Coulomb-only runs, precisely confirming the
    #   structural prediction above: Breit enters the CI step cleanly (no crash, sensible small energy shift,
    #   ~1.6e10 Hz shift in the level's absolute energy -- see the "Energy [Hz]" column), but has ZERO effect
    #   on this single-CSF system's A-constant, exactly as predicted from JAC's architecture, not a surprise.
    wa1 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-DFS-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., "uniform", 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-AL-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., "uniform", 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField()),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa2)

    wa3 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-DFS-Breit", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., "uniform", 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-AL-Breit", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., "uniform", 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField(),
                                                                  eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa4)
    #
elseif  true
    # Last successful:  30-Jul-2026 (J=3/2/kappa=+2 level ONLY -- see OPEN BUG below for J=5/2/kappa=-3)
    # Branch e (M3 fix verification): Sc (Z=21) [Ar] 3d 4s^2 ^2D_3/2,5/2 -- a single-open-subshell (3d^1)
    #   system, same structural class as branch b/c's Na 3s, chosen specifically to verify TWO bugs just
    #   fixed in the M3 (magnetic-octupole) HFS pathway: (1) Hfs.computeInteractionMatrix's calcM3 block
    #   called InteractionStrength.hfs_tE2 instead of hfs_tM3 (a copy-paste dispatch bug, module-Hfs.jl);
    #   (2) hfs_tM3 itself was missing the alpha (fine-structure constant) prefactor required by Andersson
    #   & Jonsson (2008), CPC 178, Eq. (49) -- the SAME bug already found+fixed in hfs_tM1 on 26-Jul-2026,
    #   recurring here since Eq. (49) is generic in the multipole rank L (module-InteractionStrength.jl).
    #   Literature target: de Groote et al., Phys. Lett. B 810 (2022) 135790 (arXiv:2005.00414), "Precision
    #   measurement of the magnetic octupole moment in 45Sc" -- reports A, B, AND C hyperfine constants for
    #   the SAME 3d4s^2 ^2D_3/2,5/2 states, at multiple theory levels including bare (uncorrelated,
    #   single-configuration) Dirac-Fock -- the honest apples-to-apples comparison for this single-CSF JAC
    #   calculation (per user's explicit request: "a comparison at DF/DFS level is enough ... just to know
    #   that this branch is running and (very) reasonable", NOT a precision nuclear-moment match).
    #   Dirac-Fock targets from the paper's Table 1 (A/mu, B/Q in MHz/mu_nuc, MHz/barn; C/Omega converted
    #   here from the paper's 10^-2 kHz/(mu_nuc x barn) into MHz/(mu_nuc x fm^2) -- see note below):
    #     D3/2:  A/mu = 49.520 MHz/mu_nuc   B/Q = 107.037 MHz/barn   C/Omega = 1.91e-7 MHz/(mu_nuc fm^2)
    #     D5/2:  A/mu = 21.066 MHz/mu_nuc   B/Q = 151.39  MHz/barn   C/Omega = 7.8e-8  MHz/(mu_nuc fm^2)
    #   Nuclear.Model's mu/Q/Omega are all left at 0. here on purpose -- Hfs.jl's printed "A/mu", "B/Q",
    #   "C/Omega" columns are computed directly from the electronic reduced matrix elements alone (see
    #   Hfs.computeAmplitudesProperties, module-Hfs.jl) and do not depend on the actual nm.mu/Q/Omega
    #   values at all, so this side-steps entirely a separate, NOT-fixed-here units inconsistency noticed
    #   while setting this branch up: Nuclear.Model.Omega's docstring says "[nuclear magnetons x barn]",
    #   but Hfs.jl's internal wx conversion factor for C/Omega (module-Defaults.jl, "moment: from nuclear
    #   magneton x fm^2 to atomic") is confirmed (by independently re-deriving the nuclear-magneton-times-
    #   area atomic-unit conversion factor by hand) to actually be in mu_nuc x fm^2, not barn -- a 100x
    #   mismatch if a caller ever sets nm.Omega following its own documented barn convention. Flagged for a
    #   future look, not fixed in this pass (module-Defaults.jl, a different module again).
    #   45Sc rms charge radius (3.6 fm) is an estimate, not independently verified (same caveat as Cd.jl's
    #   Ge radius) -- not expected to matter much for this DF-level sanity check.
    #   RESULT for J=3/2 (valence 3d_3/2, kappa=+2): A/mu=58.41 (target 49.52, +18%), B/Q=103.11 (target
    #   107.04, -3.7%), C/Omega=2.429e-7 (target 1.91e-7, +27%, SAME SIGN) -- all three constants land
    #   within the expected ~20-30% "single-CSF vs Dirac-Fock" ballpark (comparable to the ~18-27% scatter
    #   already seen between A/mu and B/Q themselves), confirming BOTH M3 fixes above are working correctly
    #   and producing physically sensible numbers, exactly what was asked for.
    #   OPEN BUG found while verifying (NOT fixed here, out of scope for this pass): for J=5/2 (valence
    #   3d_5/2, kappa=-3), A/mu=25.94 is still reasonable (target 21.07, +23%), but B/Q=2.647e+7 (target
    #   151.4 -- off by ~1.7e5x) and C/Omega=2.114e+4 (target 7.8e-8 -- off by ~2.7e11x) are catastrophically
    #   wrong, driven by amplitudeE2/amplitudeM3 themselves blowing up (1.49e5 and -9.58e7 respectively, vs
    #   0.60 and -0.00225 for the J=3/2 level) while amplitudeM1 stays normal for both levels. So M1 (rank 1)
    #   is fine for kappa=-3, but E2/M3 (rank>=2) are not -- a DIFFERENT symptom pattern from, but the same
    #   TRIGGER CONDITION (kappa<=-3) as, the still-open LandeZeeman N1 bug (see example-Cd.jl branch c and
    #   memory project_zeeman_hfs_bugs.md, which affects a rank-1 term instead) -- flagged as a related but
    #   presumably DIFFERENT bug (different code path: Hfs's E2/M3 coefficient/amplitude machinery, not
    #   LandeZeeman's zeeman_n1/CL_reduced_me_rb), not investigated further per the user's explicit request
    #   to keep this pass at "runs and is reasonable", not a deep dive. Only the J=3/2 result above is dated;
    #   J=5/2's E2/M3 numbers are NOT to be trusted (Rule 7).
    nm = Nuclear.Model(21., "Fermi", 45., 3.6, AngularJ64(7//2), 0.0, 0.0, 0.0)
    wa = Atomic.Computation(Atomic.Computation(), name="Cb-e-Sc3d4s2-DFS", grid=Radial.Grid(true),
                            nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d 4s^2")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=true, calcM3=true, printBefore=true) ] )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

