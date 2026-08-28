#
println("Ds) Test of the CoulombExcitation module with ASF from an internally generated initial- and final-state multiplet.")


if  true
    # Last successful:  28-Aug-2026 -- RE-RUN. This branch's hard, quantitative claim is a SYMMETRY, and it still
    # holds bit-for-bit: sigma(Mi,Mf) = sigma(-Mi,-Mf) with no post-hoc averaging. At Tp = 600 MeV/u the 0+ -> 1-
    # line gives Mf = -1 and Mf = +1 both exactly 5.809794e-01 barn, against 1.431879e-01 for Mf = 0, total
    # 1.305147e+00 barn. The 0+ -> 0- line is identically zero at every energy, as it must be.
    # NO ABSOLUTE CROSS SECTION WAS EVER RECORDED HERE, so nothing else could be compared; the A2 discussion below
    # is qualitative and is unchanged.
    # Previously:  23-Jul-2026 -- He-like U90+, K-shell -> L-shell (1s^2 -> 1s2p) Coulomb excitation, scanned
    # over Tp=10..600 MeV/u, following the K->L benchmark case of Surzhykov, Jentschura, Stohlker, Gumberidze,
    # Fritzsche, Phys. Rev. A 77, 042722 (2008) (papers/b07.pra-excitation-alignment-original.pdf).
    #
    # ROOT CAUSE FOUND (previously worked around by a pragmatic Mi,Mf<->-Mi,-Mf symmetrization, now removed):
    # computeAmplitude() was missing the cmplx(zero,-one) = -i prefactor that RATIP's coulex_pure_matrix() applies
    # to every magnetic (F23-based) contribution but that Eq. (8) of the paper does not show explicitly. With this
    # factor restored, sigma(Mi,Mf) = sigma(-Mi,-Mf) now holds EXACTLY, with no post-hoc averaging, confirmed both
    # algebraically (the L=t bracket becomes A + i*beta*CG(Mf)*B, so bracket(-Mf)=conj(bracket(Mf)), giving equal
    # magnitudes) and empirically. Two further real, independent bugs were found and fixed along the way in shared
    # JAC utilities (both confirmed against textbook values, both verified not to regress test/runtests.jl, 40/40
    # pass): AngularMomentum.sphericalYlm's missing (-1)^m phase for negative m; AngularMomentum.ClebschGordan's
    # phase using Jab instead of the correct Mab.
    #
    # A2 alignment comparison against Fig. 2 (computed by hand from the sigma(Mf) table printed below; the
    # displayCrossSections table does not yet print alignment directly -- see settings.calcAlignment output on
    # the returned Line objects instead): the QUALITATIVE trend matches the paper (negative A2 at low Tp,
    # crossing to positive, increasing with Tp -- i.e. relativistic/magnetic effects favouring |Mf|=1 at high
    # energy, exactly as Sec. IV of the paper describes). Level 4 (1s2p, J=1) crosses zero between Tp=100-165
    # MeV/u (paper: 165 MeV/u for 1s2p1/2 3P1) but reaches A2=+1.00 at 600 MeV/u vs. the paper's +0.45; level 2
    # (also J=1) crosses zero between Tp=50-100 MeV/u (paper: 324 MeV/u for 1s2p3/2 1P1) and reaches +1.10 at 600
    # MeV/u vs. the paper's +0.27. So the CROSSING ORDER (one state crosses earlier than the other) and the
    # SHAPE match, but the crossing energies and asymptotic magnitudes do not yet match quantitatively (roughly
    # 2x too large at 600 MeV/u). Most likely explanation: this branch uses only Configuration("1s 2p") as the
    # final-state basis (missing 1s2s, no Breit interaction), while the paper uses a fuller MCDF treatment --
    # NOT yet confirmed; a fairer quantitative test would add 1s2s to finalConfigs and enable the Breit
    # interaction before concluding there is a remaining bug.
    setDefaults("print summary: open", "zzz-CoulombExcitation.sum")
    ceSettings  = CoulombExcitation.Settings(CoulombExcitation.Settings(),
                                             ionEnergies=[10., 50., 100., 165., 200., 250., 300., 324., 400., 500., 600.],
                                             calcAlignment=true, printBefore=false, zerosGL=10)
    wa = Atomic.Computation(Atomic.Computation(), name = "U90-KL", grid = Radial.Grid(true), nuclearModel = Nuclear.Model(92.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s 2p")],
                            processSettings = ceSettings )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end

