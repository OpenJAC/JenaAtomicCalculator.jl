
println("Dv) Apply & test the ElectronCapture module: capture rates and the ALIGNMENT of the captured state.")
#
# WHAT THIS MODULE IS FOR, since the question is a fair one.  A capture rate cannot be measured on its own: the
# resonance always decays, and what an experiment sees is the decay.  The DR resonance STRENGTH is not the answer
# either -- it contains the branching ratio A_r/(A_a + A_r), so it is not a capture-only quantity, and it is what
# DielectronicRecombination already produces.  What DOES belong to the capture alone is the ALIGNMENT of the
# doubly excited state: the capture of an electron from a definite direction populates the magnetic sublevels of
# the resonance unequally, and in a pure two-step process every angular and polarization property of what follows
# -- the angular distribution of the Auger electron, the anisotropy and linear polarization of the DR satellite
# photon -- is fixed by that alignment times the decay amplitudes.
#
#   Balashov, Grum-Grzhimailo & Kabachnik, "Polarization and Correlation Phenomena in Atomic Collisions", Sec. 4.3.2
#   Fritzsche, Kabachnik & Surzhykov, Phys. Rev. A 78, 032703 (2008), Eqs. (3) and (4)
#
# The module implements Eq. (4) of the latter for A_k0, and Eq. (3) for the capture rate itself.

if  false
    # Last visit:  16-Oct-2026
    # Last successful:  unknown ...
    #
    # a) K-LL DIELECTRONIC CAPTURE into helium-like carbon: e- + C(5+) 1s -> C(4+) 2l2l'.
    #    The point of the branch is the ALIGNMENT column, and the check to make on it is a selection rule rather
    #    than a reference number: A_k0 is nonzero only for EVEN k <= 2J_d, so every J_d = 0 resonance must show
    #    NO alignment at all and every J_d = 1 resonance exactly one number, A_20.  An unpolarized electron beam
    #    on an unpolarized ion can ALIGN the resonance but never ORIENT it, which is why the odd ranks are absent.
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 12.0)
    ecSettings = ElectronCapture.Settings(ElectronCapture.Settings(); maxKappa = 5, calcAlignment = true,
                                           printBefore = true)
    wa = Atomic.Computation(Atomic.Computation(); name = "K-LL capture into He-like C", grid = grid,
                            nuclearModel   = Nuclear.Model(6.),
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2s^2"), Configuration("2s 2p"), Configuration("2p^2")],
                            processSettings = ecSettings )
    perform(wa)
    #
elseif  false
    # Last visit:  02-Sep-2026
    # Last successful:  unknown ...
    #
    # b) DETAILED BALANCE, which is the strongest test in this file because it is an EXACT identity between two
    #    modules rather than a comparison with a stored number.  Eq. (3) of the 2008 paper reads
    #
    #        P_cap(0 -> d)  =  (2J_d + 1) / (2 (2J_0 + 1))  *  P_A(d -> 0)
    #
    #    so the capture rate this module computes and the Auger rate AutoIonization computes for the REVERSED
    #    pair must stand in exactly that ratio, level by level.  Any deviation is a defect in one of the two.
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 12.0)
    ionConf   = [Configuration("1s")]
    resConfs  = [Configuration("2s^2"), Configuration("2s 2p"), Configuration("2p^2")]
    wcap = perform( Atomic.Computation(Atomic.Computation(); name = "capture", grid = grid,
                    nuclearModel = Nuclear.Model(6.), initialConfigs = ionConf, finalConfigs = resConfs,
                    processSettings = ElectronCapture.Settings(ElectronCapture.Settings();
                                        maxKappa = 5, calcAlignment = false) ); output = true )
    waug = perform( Atomic.Computation(Atomic.Computation(); name = "Auger (the reverse)", grid = grid,
                    nuclearModel = Nuclear.Model(6.), initialConfigs = resConfs, finalConfigs = ionConf,
                    processSettings = AutoIonization.Settings(AutoIonization.Settings(); maxKappa = 5) );
                    output = true )
    capLines = wcap[ElectronCapture.Settings];    augLines = waug[AutoIonization.Settings]
    println("\n  Detailed balance,  P_cap / P_A  against  (2J_d+1) / (2(2J_0+1)) :")
    for  cl in capLines
        for  al in augLines
            if  cl.initialLevel.index != al.finalLevel.index  ||  cl.finalLevel.index != al.initialLevel.index
                continue
            end
            al.totalRate == 0.  &&  continue
            Jd = Basics.twice(cl.finalLevel.J);    J0 = Basics.twice(cl.initialLevel.J)
            expected = (Jd + 1) / (2 * (J0 + 1))
            @printf("    d=%d -> 0=%d :  P_cap/P_A = %.6f   expected = %.6f   ratio = %.6f\n",
                    cl.finalLevel.index, cl.initialLevel.index, cl.totalRate/al.totalRate, expected,
                    (cl.totalRate/al.totalRate)/expected)
        end
    end
    #
elseif  false
    # Last visit:  02-Sep-2026
    # Last successful:  unknown ...
    #
    # c) THE SAME AT HIGH Z, where the 2008 paper puts its emphasis: K-LL capture into helium-like ions of a
    #    heavy element, for which the alignment of the L_1/2 L_3/2 resonances is what makes the subsequent
    #    K-alpha hypersatellite anisotropic.  Xenon (Z = 54) is used here rather than the paper's uranium so
    #    that the branch stays a test rather than an overnight run; the qualitative feature to look for is that
    #    |A_20| GROWS with Z, because the capture increasingly favours particular partial waves.
    #    NOTE Rule 12: the box is set by the n = 2 orbitals of a highly charged ion, so it is SMALL.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 5.0e-3, rbox = 3.0)
    ecSettings = ElectronCapture.Settings(ElectronCapture.Settings(); maxKappa = 6, calcAlignment = true)
    wa = Atomic.Computation(Atomic.Computation(); name = "K-LL capture into He-like Xe", grid = grid,
                            nuclearModel   = Nuclear.Model(54.),
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2s^2"), Configuration("2s 2p"), Configuration("2p^2")],
                            processSettings = ecSettings )
    perform(wa)
    #
    #
elseif  true
    # Last visit:  21-Oct-2026
    # Last successful:  21-Oct-2026   ## A_20 = -0.79999 (J = 2) and -0.93832 (J = 1) for 2s_1/2 2p_3/2, with the
    #                                 ## six-resonance table below;  CoulombBreit(1.), maxKappa = 6, rbox = 5 a.u.
    #
    #   THE REFERENCE FOR THIS BRANCH IS NOW JAC'S OWN RESULT, on the maintainer's instruction of 21-Oct-2026, and
    #   no longer the published pair.  The published values stay below as the physical comparison they are, but what
    #   a later run must REPRODUCE is this table, computed twice from independent scripts and agreeing to all five
    #   digits:
    #
    #        level   J^P      A_20       A_40      resonance
    #          3     1 -    +0.33084      --       2s_1/2 2p_1/2
    #          5     2 -    -0.79999   -0.16401    2s_1/2 2p_3/2
    #          6     1 +    -0.22995      --       2p^2
    #          7     2 +    -1.19175   +1.58802    2p^2
    #          8     1 -    -0.93832      --       2s_1/2 2p_3/2
    #          9     2 +    -1.10437   +1.19722    2p^2
    #
    #   Levels 1, 2, 4 and 10 have J_d < 1 and carry no alignment.  Two companion numbers are worth having beside
    #   the table, because they are what make it interpretable: the capture rates run 1.3e+13 to 3.0e+14 /s, and
    #   each of these captures proceeds through exactly TWO partial waves, so there is no partial-wave truncation
    #   error in any of these values.
    #
    #   THE 02-Sep-2026 CLAIM IS WITHDRAWN, and the reason is written here rather than deleted with it.  That
    #   claim read "both published alignment parameters REPRODUCED TO 0.9 %", on the values -0.89776 (J = 2)
    #   and -0.92599 (J = 1).  Both came from a DEFECTIVE Breit kernel: InteractionStrength.besselPhiPsi
    #   evaluated the spherical Bessel functions from a twelve-term small-argument series with no large-argument
    #   branch, and in a CAPTURE calculation the continuum electron carries full amplitude out to the box wall,
    #   where the argument omega*r reaches ~130.  The attribution is not argued but measured: an isolated clone
    #   of the repaired code with only that one route removed returns -0.89776 and -0.92599 to all five digits.
    #
    #   AND THE REFUTATION WAS PRINTED IN THE SAME TABLE ALL ALONG, which is the part worth remembering: in that
    #   clone the capture RATES come out at 1e+35 to 1e+44 /s, against 1e+13 to 1e+14 here and a physical ceiling
    #   near one per atomic time unit, ~1e+18.  The amplitudes were twelve orders of magnitude too large -- and
    #   A_20 = rho_20/rho_00 is a RATIO, which divides that common factor straight out, so it went on looking
    #   civilised and inside its physical range.  When a claim rests on a ratio, check one ABSOLUTE number from
    #   the same run against a physical bound.
    #
    #   WHAT NOW HOLDS DIFFERS BY LEVEL, and the reason is structural rather than numerical.  The selection rules
    #   from a 1s initial ion allow at most TWO partial waves into each final symmetry, so A_20 here is a
    #   function of exactly one modulus ratio x = |c2/c1|^2 and one relative phase, with no partial-wave
    #   truncation error at all, and the two single-wave values below are pure angular geometry:
    #
    #        resonance                  published    CoulombBreit(1.)   CoulombGaunt()   one wave alone
    #        2s_1/2 2p_3/2  J = 1 (1-)    -0.918         -0.93832          -0.92138      -0.70711 (p_3/2)
    #        2s_1/2 2p_3/2  J = 2 (2-)    -0.890         -0.79999          -0.80719      -0.83666 (p_3/2)
    #
    #   THE J = 1 COMPARISON WITH THE PUBLISHED VALUE HOLDS.  Its two channels, p_3/2 and p_1/2, are both Coulomb-
    #   allowed, and they interfere with the SAME sign the published value needs.  Holding the relative phase at
    #   its measured value, -0.918 corresponds to a channel ratio x = |c2/c1|^2 of 0.0162 against our 0.0198 --
    #   22 % in one small amplitude, which is what a single-configuration basis against the paper's MCDF should
    #   look like.  Here the magnetic interaction is decisive rather than a correction: Coulomb alone gives
    #   -1.27103, which is 38.5 % out.
    #
    #   THE J = 2 COMPARISON DOES NOT HOLD, AND ITS 10 % IS A SIGN RATHER THAN AN ACCURACY.  The VALUE above is this
    #   branch's reference like any other;  what is not claimed is agreement with -0.890.  Its second channel is f_5/2,
    #   which needs rank >= 3 and is therefore COULOMB-DARK -- the Coulomb-only run returns -0.83666, the pure p_3/2
    #   geometry, to five decimals -- so the whole alignment rides on one magnetically created amplitude entering as
    #   sqrt(x).  Every magnetic option available (Gaunt, static Breit, retarded Breit at four frequencies) puts the
    #   interference ABOVE that geometric value, while -0.890 lies BELOW it: the two sit on opposite sides of -0.83666.
    #   Put most sharply, and this is the strongest statement available because the relative phase is the one quantity
    #   a kernel error cannot touch: AT THE MEASURED PHASE, -0.890 IS UNREACHABLE AT ANY CHANNEL RATIO WHATEVER.  A
    #   larger basis changes x and not a phase, so this is not a correlation shortfall, and the level must not be
    #   quoted AS AGREEING with the published value in either direction until the relative sign of the p_3/2 and
    #   f_5/2 capture amplitudes is settled.
    #
    #   The other levels printed are different resonances: level 3 (1-) is 2s_1/2 2p_1/2 and level 9 (2+) is 2p^2.
    #
    # d) RECONSTRUCT THE PUBLISHED ALIGNMENTS.  Fritzsche, Kabachnik & Surzhykov, PRA 78, 032703 (2008), Sec. IV B,
    #    for the K-LL dielectronic recombination of hydrogenlike U(91+) into helium-like uranium, quote
    #
    #        A_20( 2s_1/2 2p_3/2, J = 2 ) = -0.890        A_20( 2s_1/2 2p_3/2, J = 1 ) = -0.918
    #
    #    obtained with MCDF wave functions and with BOTH the static Coulomb repulsion and the Breit interaction in
    #    the capture amplitudes -- hence CoulombBreit() here rather than the default CoulombInteraction().
    #
    #    WHAT AGREEMENT WOULD MEAN, and what it would not.  This branch uses a single-configuration basis where the
    #    paper used MCDF, so the two are not the same calculation and exact agreement is not the target; landing
    #    near -0.89 and -0.92 would show that the angular algebra of Eq. (4) and the capture amplitudes are right,
    #    since the alignment is a RATIO in which the common radial factors cancel.  A large discrepancy would point
    #    at the amplitudes, not at the formula, because the formula is already verified independently: for a single
    #    partial wave it reproduces the pure geometric values (-0.70711 for J_d=1 p_3/2, -0.83666 and -0.95618 for
    #    J_d=2), and over all two-wave mixtures A_20 stays inside its physical range [-sqrt(2), +1/sqrt(2)].
    #
    #    THE GRID IS A COMPROMISE, and the reason is worth stating: the n = 2 orbitals of uranium turn over near
    #    4/92 = 0.04 a.u. and want a TINY box, while Continuum.gridConsistency refuses anything whose normalisation
    #    point sits inside 2 a.u.  The box and the inner resolution are fixed by DIFFERENT physics, so the box is
    #    set by the continuum (5 a.u.) and the bound orbitals are carried by a very fine logarithmic inner mesh
    #    (rnt = 1e-7, h = 0.03) rather than by a small box.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-7, h = 3.0e-2, hp = 2.0e-3, rbox = 5.0)
    ecSettings = ElectronCapture.Settings(ElectronCapture.Settings(); maxKappa = 6, calcAlignment = true,
                                           operator = CoulombBreit(1.0))
    wa = Atomic.Computation(Atomic.Computation(); name = "K-LL capture into He-like U (the 2008 paper)", grid = grid,
                            nuclearModel   = Nuclear.Model(92.),
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2s^2"), Configuration("2s 2p"), Configuration("2p^2")],
                            processSettings = ecSettings )
    wb = perform(wa; output = true)
    println("\n  This branch's reference (21-Oct-2026, CoulombBreit(1.), maxKappa = 6):  " *
            "level 5 (2-) = -0.79999,  level 8 (1-) = -0.93832")
    println("  Published (PRA 78, 032703 (2008), Sec. IV B):  A_20(2s2p_3/2, J=2) = -0.890,  " *
            "A_20(2s2p_3/2, J=1) = -0.918  -- the J = 1 comparison holds, the J = 2 one does not;  see above")
    for  ln in wb[ElectronCapture.Settings]
        length(ln.alignment) == 0  &&  continue
        @printf("    level %2d   J^P = %s   A_20 = %9.5f\n", ln.finalLevel.index,
                string(LevelSymmetry(ln.finalLevel.J, ln.finalLevel.parity)), ln.alignment[1])
    end
    #
end
