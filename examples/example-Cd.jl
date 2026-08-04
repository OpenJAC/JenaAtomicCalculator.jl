
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
                            nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8797, AngularJ64(1//2), 2.7928, 0.0, 0.0),
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
                            nuclearModel=Nuclear.Model(2., "uniform", 3., 1.881, AngularJ64(1//2), -2.12749772, 0.0, 0.0),
                            configs=[Configuration("1s 2p")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown -- deliberately left undated, see below (known bug).
    # Branch c: Ge II (Z=32) 4s^2 4f (2)F_5/2,7/2, the paper's own single-f-electron test case (their label
    #   "4s^2 4f^2" cannot literally mean two f-electrons -- equivalent f^2 gives only singlet/triplet terms,
    #   never a doublet; re-read here as "4s^2 4f (2)F", one f-electron on a closed 4s^2 spectator shell,
    #   which matches the paper's bare g_J numbers 0.857136/1.142850 almost exactly to the nonrelativistic
    #   Lande fractions 6/7, 8/7). Uses [Ar] 3d^10 4s^2 4f (31 electrons = Ge+); 74Ge, I=0, rms radius 4.07 fm
    #   is an estimate, not independently verified.
    #   OPEN BUG (found 25-Jul-2026, RE-VERIFIED 31-Jul-2026 after the unrelated B-spline kappa-sign fix --
    #   see memory project_zeeman_hfs_bugs.md): only the J=5/2 level (f_5/2, kappa=+3) is trustworthy here --
    #   computed g_J = 0.856825 vs paper's total (with Schwinger) 0.856804, a good match (even slightly
    #   closer than the pre-fix 0.857095). The J=7/2 level (f_7/2, kappa=-4) is WRONG (computed -2.263670 vs
    #   paper's 1.143182) -- a newly found Zeeman N1 bug affects every kappa<=-3 subshell (d5/2, f7/2, ...),
    #   not just this one -- confirmed via a kappa scan across p1/2,p3/2,d3/2,d5/2,f5/2,f7/2 in this same
    #   ion, all at the same Z. NOTE: the WRONG value itself is not stable -- before the 31-Jul B-spline fix
    #   this same J=7/2 level gave +0.893667 (smaller magnitude, same sign as the exact target); after the
    #   fix it flipped sign and grew to -2.263670. Root cause not yet isolated (AngularMomentum.CL_reduced_me_rb
    #   or its use in InteractionStrength.zeeman_n1, not yet narrowed further) -- evidently in the angular
    #   coefficient machinery, since it reacts to an unrelated radial-basis change. NOT fixed in code. Left
    #   undated on purpose (Rule 7): the J=7/2 half of this branch's output is known-wrong.
    nm = Nuclear.Model(32., "Fermi", 74., 4.07, AngularJ64(0//1), 0.0, 0.0, 0.0)
    wa = Atomic.Computation(Atomic.Computation(), name="Cd-c-GeII4f", grid=Radial.Grid(true),
                            nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  31-Jul-2026 (re-dated: the 30-Jul-2026 numbers below were stale after the same-day
    #   B-spline kappa-sign boundary-condition fix landed later that session -- see project_zeeman_hfs_bugs.md)
    # Branch d (AL-Field / Breit revisit): reuses branch c's Ge II [Ar] 3d^10 4s^2 4f system, now that
    #   AL-Field (Basics.ALField()) has been root-cause-fixed and promoted to the standard implementation
    #   (see project_df_al_kink_bug.md). ONLY the J=5/2 (f_5/2, kappa=+3) level is used for this comparison
    #   -- the J=7/2 (f_7/2, kappa=-4) level is deliberately excluded, since it is already known-wrong from
    #   the still-open kappa<=-3 Zeeman N1 bug (branch c), and mixing that in would make it impossible to
    #   tell an AL-Field/Breit effect apart from the pre-existing bug. Four settings combinations: DFS+
    #   Coulomb (reproduces branch c's f_5/2 baseline exactly, g_J=0.856824955), AL+Coulomb, DFS+Breit,
    #   AL+Breit.
    #   RESULT (2), Breit: g_J(f_5/2) is EXACTLY bit-identical between the Coulomb-only and Breit-added runs
    #   (both for DFS and for AL) -- the third and final confirmation this session of the same structural
    #   prediction (single-CSF-per-symmetry-block system here too: each J level is its own trivial 1x1
    #   block, no CI mixing for Breit to affect).
    #   RESULT (1), AL-Field: g_J(f_5/2) shifts DFS->AL by only -0.0022% (0.856824955 -> 0.856806270) -- BY
    #   FAR the smallest of the three properties tested this session (HFS: 26%, isotope shift K_nms/F:
    #   ~0.1-0.2%, Zeeman: 0.002%). This is exactly the predicted "or nothing significant is changing" result
    #   -- Zeeman g_J is dominated by angular/kinematic Lande structure, essentially insensitive to whether
    #   the SCF exchange treatment is DFS (local) or AL (non-local). (For reference/curiosity only, NOT a
    #   trustworthy data point: g_J(f_7/2) is -2.263670 (DFS) vs -2.242473 (AL) -- but f_7/2 is already wrong
    #   from the separate, pre-existing kappa<=-3 bug, so this reflects how that bug interacts with
    #   different orbital shapes, not a genuine AL-Field physics effect; not used for any conclusion here,
    #   per the exclusion stated above. Also NOTE: these f_7/2 numbers themselves shifted substantially from
    #   the 30-Jul-2026 run (then: 0.893667/1.142966) after the same-day B-spline fix -- see branch c and
    #   the KNOWN LIMITATION docstring on LandeZeeman.amplitude(::ZeemanN1,...) for the full picture.)
    nm = Nuclear.Model(32., "Fermi", 74., 4.07, AngularJ64(0//1), 0.0, 0.0, 0.0)

    wa1 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-DFS-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-AL-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField()),
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa2)

    wa3 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-DFS-Breit", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d^10 4s^2 4f")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true,
                                                includeSchwinger=true, printBefore=true) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cd-d-GeII4f-AL-Breit", grid=Radial.Grid(true),
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
    # Last successful:  31-Jul-2026
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
    nm = Nuclear.Model(20., "Fermi")

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


