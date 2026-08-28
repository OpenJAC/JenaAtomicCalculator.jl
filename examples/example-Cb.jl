
println("Cb) Apply & test the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
println("    Branches follow Andersson & Jonsson, CPC 178 (2008) 156-170 (examples/papers/")
println("    2008.cpc-andersson-jonsson-hfs-zeemann.pdf); see project memory project_zeeman_hfs_bugs.md.")

Defaults.setDefaults("unit: energy", "Hz")
setDefaults("print summary: open", "zzz-Hfs.sum")

if  false
    # Last successful:  28-Aug-2026 -- A(1s) now reads 1.421215e+03 MHz where this branch recorded 1423.157 MHz, a
    # move of 0.14 % and TOWARD the measured value: agreement with the famous 1420.405751 MHz 21 cm line improves
    # from 0.19 % to 0.057 %. The record is from 26-Jul, so the shift is not attributable to any one change; it is
    # the same small radial/SCF drift seen in example-Dn.jl on 28-Aug.
    # Previously:  26-Jul-2026
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
                            nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8797, AngularJ64(1//2), 2.7928, 0.0, 0.0),
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
                            nuclearModel=Nuclear.Model(11., UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
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
                            nuclearModel=Nuclear.Model(11., UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-AL-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField()),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa2)

    wa3 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-DFS-Breit", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cb-c-Na3s-AL-Breit", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(11., UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0),
                            configs=[Configuration("[Ne] 3s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField(),
                                                                  eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )
    perform(wa4)
    #
elseif  false
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
    nm = Nuclear.Model(21., FermiNucleus(), 45., 3.6, AngularJ64(7//2), 0.0, 0.0, 0.0)
    wa = Atomic.Computation(Atomic.Computation(), name="Cb-e-Sc3d4s2-DFS", grid=Radial.Grid(true),
                            nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d 4s^2")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=true, calcM3=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  true
    # Last successful:  28-Aug-2026 -- REPRODUCES EXACTLY, and being dated 18-Aug it is the tighter of this file's
    # two checks: Ho I A = 8.957463e+02 MHz and B = -9.516364e+02 MHz against the recorded 895.75 and -951.6;
    # Pr I A = 9.301445e+02 MHz and B = -5.849614e+00 MHz against the recorded 930.14 and -5.85.
    # Previously:  18-Aug-2026 (Ho I; the Pr I comparison is NOT dated -- see the citation caveat)
    # Branch f (open f shell): the hyperfine constants of two NEUTRAL lanthanides, Ho I [Xe] 4f^11 6s^2
    #   ^4I_15/2 and Pr I [Xe] 4f^3 6s^2 ^4I_9/2.  These are the branches that were blocked from 11-Aug-2026
    #   until 18-Aug-2026, and it is worth recording why, because the reason was not what it looked like.
    #
    #   WHAT BLOCKED THEM.  Both atoms came out with an INVERTED 4f spin-orbit splitting -- 4f_7/2 more
    #   bound than 4f_5/2 -- which inverts the level ordering directly, so Ho I gave a J = 9/2 ground level
    #   where it must be 15/2, and Pr I gave 15/2 where it must be 9/2.  Comparing A and B with measurement
    #   is meaningless while the ground level is the wrong one, so the HFS work was postponed with it.
    #   Re-measured on 18-Aug the inversion is GONE: Pr I gives J = 9/2 and Ho I J = 15/2, both correct, and
    #   with the right multiplet DIRECTION -- regular (ascending) for the less-than-half-filled 4f^3, and
    #   inverted (descending) for the more-than-half-filled 4f^11, which is what Hund's third rule requires.
    #   The 4f_5/2 - 4f_7/2 splitting is now normal at 2714.65 cm^-1, i.e. zeta_4f = 776 cm^-1 against the
    #   ~750 cm^-1 expected.  NOTE the trap that made this look worse than it was: the splitting is
    #   zeta*(l+1/2) = 3.5 zeta for l = 3, so comparing 2715 cm^-1 against zeta itself suggests a 4x error
    #   where there is none.  It was NOT the radial box either -- rbox = 25 gives the same normal ordering
    #   today, so the fix is in the SCF; which commit fixed it has not been established.
    #
    #   NUCLEAR DATA.  165Ho: I = 7/2, mu = +4.173 mu_N, Q = +3.49 b, rms radius 5.21 fm.
    #                  141Pr: I = 5/2, mu = +4.2754 mu_N, Q = -0.0589 b, rms radius 4.99 fm.
    #   The box is taken from Basics.recommendedGrid, which sizes it from the configuration (58.95 a.u. for
    #   neutral Pr I, set by the 6s seeing Zeff = 2.85).  Neither atom is expensive: 41 CSFs over 9 J^P
    #   symmetries, largest block 7x7, with the runtime almost entirely SCF.
    #
    #   RESULT, ground level of each:
    #                        JAC            measured        deviation
    #       Ho I  A       895.75 MHz       800.583 MHz       +11.9 %
    #       Ho I  B      -951.6  MHz     -1668     MHz       43 % low in magnitude
    #       Pr I  A       930.14 MHz     (see caveat)
    #       Pr I  B        -5.85 MHz     (see caveat)
    #
    #   The Ho I values are Dankwort et al. (1974), the ABMR ground-state measurement.  THE Pr I COMPARISON
    #   IS DELIBERATELY NOT QUOTED and this branch is not dated on it: the value ~926 MHz that comes to mind
    #   for a Pr ground state belongs, as far as can be checked here, to Pr II 4f^3 6s ^5I_4 and not to
    #   neutral Pr I, and Ginibre, Physica Scripta 39, 694 / 39, 710 (1989) -- the obvious source -- is Pr II
    #   in its first part.  A 0.4 % "agreement" with a value taken from the wrong charge state would be worse
    #   than no comparison at all.  The neutral Pr I ground-state constants need the ABMR literature and a
    #   verified citation before this half means anything.
    #
    #   WHY A IS GOOD AND B IS NOT, which is the physics of this branch rather than a defect.  A and B are
    #   built from the SAME 4f orbital, so a simple error in its radial extent would move both the same way
    #   by a similar factor.  Instead A is ~12 % HIGH and B is ~43 % LOW -- opposite directions -- so this is
    #   not a bad <r^-3>.  It is a many-body contribution missing from B specifically: a single-configuration
    #   calculation has no core polarization at all, and for an f electron the polarized core ANTISHIELDS the
    #   field gradient, enhancing |B|.  Leaving it out must make |B| too small, which is the sign obtained.
    #   The magnetic constant, dominated by the 4f orbital itself, survives the same omission far better.
    #   So ~12 % on A from a bare single configuration on an open-f-shell neutral is a respectable result,
    #   and B needs correlation before it can be compared at all -- the standard situation for lanthanide
    #   quadrupole constants, not a fault of this calculation.
    #
    #   DEFECT FOUND WHILE RUNNING THIS, not fixed here (a different module):  the "g_J" column of Hfs's
    #   "HFS amplitudes and g_J factors" table is never computed.  Hfs.Outcome is constructed with its gJ
    #   field hard-wired to 1. (module-Hfs.jl:464) and nothing ever assigns it, so every HFS table JAC has
    #   printed shows g_J = 1.00000000 for every level, as if computed.  The correct route exists and is
    #   validated -- LandeZeeman, cf. example-Cd.jl -- so this is a routing fix, not new physics.
    for  (label, Z, A, rad, I, mu, Q, conf)  in
            [ ("Ho I 4f^11 6s^2", 67., 165., 5.21, 7//2, 4.173,  3.49,    Configuration("[Xe] 4f^11 6s^2")),
              ("Pr I 4f^3  6s^2", 59., 141., 4.99, 5//2, 4.2754, -0.0589, Configuration("[Xe] 4f^3 6s^2")) ]
        println("\n\n***  $label  ***\n")
        nm   = Nuclear.Model(Z, UniformNucleus(), A, rad, AngularJ64(I), mu, Q, 0.0)
        grid = Basics.recommendedGrid([conf], nm; printout=true)
        wa   = Atomic.Computation(Atomic.Computation(), name="Cb-f-$label", grid=grid, nuclearModel=nm,
                                  configs=[conf],
                                  propertySettings=[ Hfs.Settings(calcM1=true, calcE2=true, printBefore=false) ] )
        perform(wa)
    end
    #
end
#
setDefaults("print summary: close", "")

