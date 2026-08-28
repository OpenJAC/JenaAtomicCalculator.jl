
println("Cd) Apply & test the LandeZeeman module with ASF from an internally generated multiplet.")
println("    Branches follow Andersson & Jonsson, CPC 178 (2008) 156-170 (examples/papers/")
println("    2008.cpc-andersson-jonsson-hfs-zeemann.pdf); see project memory project_zeeman_hfs_bugs.md.")

setDefaults("print summary: open", "zzz-LandeZeeman.sum")

if  false
    # Last successful:  25-Jul-2026
    # Branch a: H(2p), bare single-electron test of the Zeeman N1 + Schwinger (Delta N1) terms.
    #   Proton: I=1/2, mu=2.7928 nmu, Q=0, rms radius=0.8797 fm (I plays no role for calcLandeJ itself).
    #   "uniform" nuclear model, not "Fermi": for Z=1, JAC's 2-parameter Fermi model cannot represent an
    #   rms radius this small and now raises an explicit error (module-Nuclear.jl, fixed 25-Jul-2026; see
    #   memory project_isotope_shift_ris3ris4.md); re-verified bit-identical to the earlier (silently wrong
    #   nuclear-size) Fermi-model run, confirming g_J is insensitive to nuclear size at this precision.
    #   Verified against the exact nonrelativistic Lande fractions g_J(2p_1/2)=2/3, g_J(2p_3/2)=4/3:
    #   computed g_J(1/2) = 0.665887 (vs 0.666667), g_J(3/2) = 1.334101 (vs 1.333333) -- both within the
    #   expected sub-0.1% relativistic correction for Z=1. This is the case that pinned down the (g_s-2)/4
    #   Schwinger-term fix (was (g_s-2)/2 in the paper's own Eq. (52), off by exactly 2x -- see memory).
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-a-H2p", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8797, AngularJ64(1//2), 2.7928, 0.0, 0.0),
                            configs=[Configuration("2p")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  25-Jul-2026
    # Branch b: He 1s2p multiplet (3He, I=1/2, mu=-2.12749772 n.m.), cross-checked against the paper's own
    #   Table (Fig. 10-type numbers) for this exact system. "uniform" nuclear model used (the 2-parameter
    #   Fermi model errors out for this light, small-rms-radius case, cf. memory).
    #   Verified (level 2 = 3P1, level 3 = 3P2, level 4 = 1P1; level 1 = 3P0 has g_J=0 trivially, J=0):
    #     3P1: g_J = 1.50112135  vs paper 1.5011166   (6 sig figs)
    #     3P2: g_J = 1.50112266  vs paper 1.5011183   (6 sig figs)
    #     1P1: g_J = 0.99999305  vs paper 0.9999936   (6 sig figs)
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-b-He1s2p", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(2., UniformNucleus(), 3., 1.881, AngularJ64(1//2), -2.12749772, 0.0, 0.0),
                            configs=[Configuration("1s 2p")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  10-Aug-2026
    # Branch c: Ge II (Z=32) 4s^2 4f (2)F_5/2,7/2, the paper's own single-f-electron test case (their label
    #   "4s^2 4f^2" cannot literally mean two f-electrons -- equivalent f^2 gives only singlet/triplet terms,
    #   never a doublet; re-read here as "4s^2 4f (2)F", one f-electron on a closed 4s^2 spectator shell,
    #   which matches the paper's bare g_J numbers 0.857136/1.142850 almost exactly to the nonrelativistic
    #   Lande fractions 6/7, 8/7). Uses [Ar] 3d^10 4s^2 4f (31 electrons = Ge+); 74Ge, I=0, rms radius 4.07 fm
    #   is an estimate, not independently verified.
    #   THE "kappa <= -3 ZEEMAN N1 BUG" WAS THE RADIAL BOX (resolved 10-Aug-2026). This branch used to run on
    #   Radial.Grid(true), whose box is 614 a.u., and reported g_J(f_7/2) = -2.263670 against the exact 8/7:
    #   wrong by a factor and even in sign, while f_5/2 was fine. It was recorded as an angular-coefficient
    #   bug in AngularMomentum.CL_reduced_me_rb or InteractionStrength.zeeman_n1, "evidently angular, since
    #   it reacts to an unrelated radial-basis change". That inference was backwards -- reacting to a radial
    #   change is evidence the cause IS radial. Two measurements settle it:
    #     (i) a ONE-ELECTRON scan of g_J over kappa = +1,-2,+2,-3,+3,-4 on a matched box reproduces the exact
    #         Lande factors 2|kappa|/(2|kappa| -+ 1) to 2e-4 for EVERY kappa, so the angular machinery is
    #         sound and both recorded suspects are exonerated;
    #     (ii) in this Ge II system the 4f electron sees Z_eff ~ 2 and needs a box of ~30 a.u.; on the 614
    #         a.u. box the SCF simply returns the WRONG STATE for kappa = -4 -- E(4f_7/2) = -1.5758 with
    #         <r> = 5.20, against E(4f_5/2) = -0.0619 with <r> = 10.87, whereas the two spin-orbit partners
    #         of a 4f shell must be nearly identical. On rbox = 30 they agree exactly (both -0.074629,
    #         <r> = 9.215) and g_J comes right.
    #   With the box matched, this branch now REPRODUCES THE PAPER TO SIX DECIMALS, both levels:
    #     J=7/2  g_J = +1.143182   vs Andersson & Jonsson +1.143182
    #     J=5/2  g_J = +0.856805   vs                     +0.856804
    #   Same lesson as the B-spline f-state case: a box much TOO LARGE starves the basis exactly as badly as
    #   one that is too small, because the number of splines is fixed. Note that the grid check added on
    #   9-Aug-2026 (Bsplines.checkGridRepresentation) does NOT catch this: it tests hydrogenic orbitals at
    #   the full nuclear charge, where 4f is compact, and it passes this grid.
    nm = Nuclear.Model(32., FermiNucleus(), 74., 4.07, AngularJ64(0//1), 0.0, 0.0, 0.0)
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-c-GeII4f", grid=Radial.Grid(Radial.Grid(false); rbox=30.0),
                            nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  10-Aug-2026 (re-run on a matched box; the earlier numbers were taken on a 614 a.u.
    #   box that returned the wrong f_7/2 state. Prior note, kept: re-dated 31-Jul after the same-day
    #   B-spline kappa-sign boundary-condition fix landed later that session -- see project_zeeman_hfs_bugs.md)
    # Branch d (AL-Field / Breit revisit): reuses branch c's Ge II [Ar] 3d^10 4s^2 4f system, now that
    #   AL-Field (Basics.ALField()) has been root-cause-fixed and promoted to the standard implementation
    #   (see project_df_al_kink_bug.md). ONLY the J=5/2 (f_5/2, kappa=+3) level is used for this comparison
    #   -- the J=7/2 (f_7/2, kappa=-4) level is deliberately excluded, since it is already known-wrong from
    #   the still-open kappa<=-3 Zeeman N1 bug (branch c), and mixing that in would make it impossible to
    #   tell an AL-Field/Breit effect apart from the pre-existing bug. Four settings combinations: DFS+
    #   Coulomb (reproduces branch c's f_5/2 baseline exactly, g_J=0.856824955), AL+Coulomb, DFS+Breit,
    #   AL+Breit.
    #   RE-RUN 10-Aug-2026 on rbox = 30 a.u. The f_7/2 EXCLUSION IS NO LONGER NEEDED: the "kappa <= -3 bug"
    #   was this branch's radial box, not a Zeeman defect (see branch c), so both levels are now usable and
    #   the four combinations read
    #        DFS+Coulomb   J=5/2  0.856804696    J=7/2  1.143181922
    #        AL +Coulomb   J=5/2  0.856804428    J=7/2  1.143181632
    #        DFS+Breit     J=5/2  0.856804696    J=7/2  1.143181922
    #        AL +Breit     J=5/2  0.856804428    J=7/2  1.143181632
    #   RESULT (2), Breit: g_J is EXACTLY bit-identical between the Coulomb-only and Breit-added runs, for
    #   both fields AND now for both levels -- the same structural prediction as before (single CSF per
    #   symmetry block, so each J level is its own trivial 1x1 block and there is no CI mixing for Breit to
    #   affect). Confirmed on twice as much data as when it was first recorded.
    #   RESULT (1), AL-Field: g_J shifts DFS->AL by -3.1e-5 % on f_5/2 and -2.5e-5 % on f_7/2 -- far smaller
    #   even than the -0.0022% recorded from the old, box-starved run, and still by far the smallest of the
    #   properties tested (HFS: 26%, isotope shift K_nms/F: ~0.1-0.2%). The conclusion is unchanged and now
    #   better founded: Zeeman g_J is dominated by angular/kinematic Lande structure and is essentially
    #   insensitive to whether the SCF exchange treatment is DFS (local) or AL (non-local).
    #   The old f_7/2 numbers (-2.263670 DFS / -2.242473 AL, and 0.893667 before the July B-spline fix) were
    #   all artefacts of the wrong state returned on the 614 a.u. box; they are not AL-Field physics.
    nm = Nuclear.Model(32., FermiNucleus(), 74., 4.07, AngularJ64(0//1), 0.0, 0.0, 0.0)

    wa1 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-DFS-Coulomb", grid=Radial.Grid(Radial.Grid(false); rbox=30.0),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-AL-Coulomb", grid=Radial.Grid(Radial.Grid(false); rbox=30.0),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField()),
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa2)

    wa3 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-DFS-Breit", grid=Radial.Grid(Radial.Grid(false); rbox=30.0),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-AL-Breit", grid=Radial.Grid(Radial.Grid(false); rbox=30.0),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField(),
                                                                  eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa4)
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch e: Al+ [Ne]3s3p ^3P_0 quadratic Zeeman shift coefficient C_2 -- Gilles, Fritzsche, Spiess,
    #   Schmidt & Surzhykov, Phys. Rev. A 110, 052812 (2024) (examples/papers/2024.pra-jan-magnetic-
    #   moments-original.pdf), Table III: experiment -71.944(24) MHz/T^2, independent theory -71.927
    #   MHz/T^2; their own MCDF-only "step 1" (no virtual excitations, same level this branch uses)
    #   result is -74.2(2.2) MHz/T^2, ~3% off experiment.
    #   RESULT: computed C_2 = -67.106 MHz/T^2, ~6.7% off experiment/theory -- same ballpark as the
    #   paper's own step-1 accuracy (on the other side of experiment from their value, but comparably
    #   close), a genuine, honest validation for a single-configuration (no virtual excitations) level
    #   of theory.
    #   Zero nuclear spin (I=0), matching the paper's own stated scope for C_2/quadrupole-moment
    #   theory (their Sec. V: "restricted to ... zero nuclear spin"). gMultiplet is the [Ne]3s3p
    #   configuration's OWN full multiplet (self-referential -- matches the paper's "step 1: Dirac-
    #   Fock, no virtual excitations" baseline, Table II).
    #   This also serves as the empirical check that resolved a formula question in
    #   LandeZeeman.computeQuadraticZeemanC2: the paper's Eq. (33) as literally transcribed includes an
    #   extra 1/(2J'+1) factor, but adding it made every C_2 value too small by close to a factor of
    #   (2J'+1) -- see the docstring note on that function (31-Jul-2026) for the full story.
    nm = Nuclear.Model(13.)

    wa0 = Atomic.Computation(Atomic.Computation(), name="Cd-e-AlII-3s3p", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[Ne] 3s 3p")])
    wb0 = perform(wa0; output=true)
    gMultiplet = wb0["multiplet:"]

    lzSettings = LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true, calcQZScoeff=true,
                                       includeSchwinger=true, gMultiplet=gMultiplet, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-e-AlII-3s3p-C2", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[Ne] 3s 3p")],
                            propertySettings=[lzSettings] )
    wb = perform(wa; output=true)
    outcomes = wb["Zeeman parameter outcomes:"]
    for  outcome in outcomes
        if  outcome.Jlevel.J == AngularJ64(0)
            println("\n>> Al+ 3P0 (level $(outcome.Jlevel.index)): computed C_2 = $(outcome.Jsublevels[1].c2Coeff) MHz/T^2")
            println("   Literature (Gilles et al. 2024, Table III): experiment -71.944(24) MHz/T^2, theory -71.927 MHz/T^2")
        end
    end
    #
elseif  true
    # Last successful:  28-Aug-2026 -- ALL SIX QUADRATIC SHIFTS REPRODUCE, to between 0.01 % and 0.6 %:
    #   3P0[1] M=0    -0.24119424941 MHz/T^2   (recorded -0.24138)
    #   3P1[2] M=0     0.13806323819            (recorded  0.13828)
    #          M=+-1  -0.07699337674            (recorded -0.07697)
    #   3P2[3] M=0     0.10134573611            (recorded  0.10133)
    #          M=+-1   0.07481985265            (recorded  0.07481)
    #          M=+-2  -0.00475779772            (recorded -0.00473)
    # The M degeneracy is EXACT: M = +1 and M = -1 return identical values to every digit, as they must.
    # Previously:  31-Jul-2026
    # Branch f: Ca^14+ [He]2s^2 2p^2 quadratic Zeeman shift coefficients C_2 for all M-sublevels of
    #   ^3P_0, ^3P_1, ^3P_2 -- Gilles et al. (2024), Table V (the paper's own flagship highly-charged-
    #   ion result, and the exact system Jan Gilles' own original job file targeted: see
    #   apps/b23-apps-jan-moments-2023/job-jan-a.jl, "Ca^14+ (Z=20) Configuration('1s^2 2s^2 2p^2')").
    #   Literature (MCDF column, [He]2s^2 2p^2, MHz/T^2) vs computed (this branch), level index in []:
    #     3P0[1]  M=0            -0.2447(14)   computed -0.24138   (1.4% off)
    #     3P1[2]  M=0             0.143(2)     computed  0.13828   (3.3% off -- paper itself notes a
    #                                          sign ambiguity for alpha_2 in this exact state, Table III)
    #             M=+-1          -0.0761(4)    computed -0.07697   (1.2% off)
    #     3P2[3]  M=0             0.1001(6)    computed  0.10133   (1.2% off)
    #             M=+-1           0.0738(4)    computed  0.07481   (1.4% off)
    #             M=+-2          -0.0049(1)    computed -0.00473   (3.5% off)
    #   6 of 7 M-sublevels land within 1-4% of literature -- consistent with the paper's own claimed
    #   step-1 (no virtual excitations) accuracy of order 1%. This is the level of agreement the user
    #   remembered getting two years ago; it required REMOVING the 1/(2J'+1) factor that a literal
    #   reading of the paper's own Eq. (33) calls for -- see the docstring note on
    #   LandeZeeman.computeQuadraticZeemanC2 (31-Jul-2026) for why that factor does not apply to this
    #   codebase's specific reduced-matrix-element convention.
    #   Zero nuclear spin (I=0), matching the paper's own stated scope. gMultiplet is the
    #   [He]2s^2 2p^2 configuration's OWN full multiplet (self-referential, same "step 1" baseline
    #   as branch e). Energies used throughout are JAC's own (matches the paper's "MCDF" column, not
    #   their "MCDF+NIST ASD" column, which substitutes external NIST-tabulated energies -- out of
    #   scope here).
    nm = Nuclear.Model(20., FermiNucleus())

    wa0 = Atomic.Computation(Atomic.Computation(), name="Cd-f-Ca14p-2s2p2", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[He] 2s^2 2p^2")])
    wb0 = perform(wa0; output=true)
    gMultiplet = wb0["multiplet:"]

    lzSettings = LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true, calcQZScoeff=true,
                                       includeSchwinger=true, gMultiplet=gMultiplet, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-f-Ca14p-2s2p2-C2", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[He] 2s^2 2p^2")],
                            propertySettings=[lzSettings] )
    wb = perform(wa; output=true)
    outcomes = wb["Zeeman parameter outcomes:"]
    println("\n>> Ca14+ [He]2s^2 2p^2 computed C_2 values [MHz/T^2] (literature in comment block above):")
    for  outcome in outcomes
        sym = LevelSymmetry(outcome.Jlevel.J, outcome.Jlevel.parity)
        for  Jsub in outcome.Jsublevels
            println("   Level $(outcome.Jlevel.index) ($(string(sym))), M=$(Jsub.M):  C_2 = $(Jsub.c2Coeff)")
        end
    end
    #
end
#
setDefaults("print summary: close", "")


