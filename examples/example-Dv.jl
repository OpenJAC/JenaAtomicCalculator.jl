
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

if  true
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
end
