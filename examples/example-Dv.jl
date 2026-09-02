
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
    # Last visit:  02-Sep-2026
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
    capLines = wcap["electron-capture lines:"];    augLines = waug["AutoIonization lines:"]
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
    # Last visit:  02-Sep-2026
    # Last successful:  02-Sep-2026 -- the two published alignment parameters are REPRODUCED TO 0.9 %:
    #
    #        resonance                    published        JAC         difference
    #        2s_1/2 2p_3/2  J = 2  (2-)     -0.890       -0.89776        0.9 %
    #        2s_1/2 2p_3/2  J = 1  (1-)     -0.918       -0.92599        0.9 %
    #
    #   with the right sign, the right magnitude AND the right ordering -- J = 1 is more strongly aligned than
    #   J = 2 in both.  About 1 % is what a single-configuration basis should give against the paper's MCDF, and
    #   the alignment is a RATIO in which the common radial factors cancel, which is why it comes out this close
    #   from a much simpler wave function.  The other levels printed are different resonances: level 3 (1-) is
    #   2s_1/2 2p_1/2 and level 9 (2+) is 2p^2.
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
    println("\n  Published (PRA 78, 032703 (2008), Sec. IV B):  A_20(2s2p_3/2, J=2) = -0.890,  " *
            "A_20(2s2p_3/2, J=1) = -0.918")
    for  ln in wb["electron-capture lines:"]
        length(ln.alignment) == 0  &&  continue
        @printf("    level %2d   J^P = %s   A_20 = %9.5f\n", ln.finalLevel.index,
                string(LevelSymmetry(ln.finalLevel.J, ln.finalLevel.parity)), ln.alignment[1])
    end
    #
end
