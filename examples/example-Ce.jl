
println("Ce) Apply & test the FormFactor module with ASF from an internally generated multiplet.")
println("    Branch a: exact sum-rule check F^(standard)(q=0) = N_e (electron number).")
println("    Branch b: exact closed-form comparison against the hydrogenic 1s form factor")
println("    F(q) = 16 Z^4 / (4 Z^2 + q^2)^2 [a.u.], derived analytically from the exact nonrelativistic")
println("    1s density -- see project memory project_formfactor_scattering.md.")

if  false
    # Last successful:  31-Jul-2026
    # Branch a: Ne-like Fe (Z=26) [Ne] 3s^2 3p^5 / [Ne] 3s 3p^6, N_e = 17 electrons for both configurations.
    #   Tests the ONE exact, model-independent identity a standard form factor must satisfy: F^(standard)(0)
    #   is simply the total electron density integrated over all space, i.e. F(0) = N_e exactly, for ANY
    #   level of ANY atom/ion -- see FormFactor.displayResults's own header note "F^(standard)(0) = N_e".
    #   This requires no external literature at all (a hard normalization sum rule, not a fitted/tabulated
    #   number) and is a good first sanity check before trusting F(q) at q != 0.
    #   VERIFIED: F^(standard)(0) = 1.70000e+01 = 17 exactly (to displayed precision) for every level. F(q)
    #   falls off monotonically and smoothly with increasing q, as expected (17.0 -> 16.99 -> 16.2 -> 4.76 ->
    #   0.111 -> 0.00314 for q=0,0.1,1,10,100,1000 a.u.). F^(modified)(q) stays close to F^(standard)(q)
    #   throughout (differing by only ~0.05-0.3%), consistent with relativistic/binding-energy corrections
    #   being a small effect for Z=26 valence electrons.
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Ce-a-FeNeLike-sumrule", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings=[ FormFactor.Settings([0., 0.1, 1.0, 10., 100., 1000.], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch b: H(1s), bare single-electron test against the EXACT closed-form hydrogenic form factor.
    #   For a nonrelativistic hydrogenic 1s density n(r) = (Z^3/pi) exp(-2Zr), the standard form factor is
    #   F(q) = integral n(r) exp(i q.r) d^3r = 16 Z^4 / (4 Z^2 + q^2)^2 [a.u.] (derived by hand: F(q) =
    #   (4Z^3/q) integral_0^infinity r exp(-2Zr) sin(qr) dr, and integral_0^infinity r exp(-ar) sin(br) dr =
    #   2ab/(a^2+b^2)^2 with a=2Z, b=q). For Z=1: F(q) = 16/(4+q^2)^2, F(0)=1=N_e (consistent with branch a's
    #   sum rule). "uniform" nuclear model and scField=Basics.NuclearField() are used, as elsewhere in this
    #   codebase for single-electron H systems, to avoid the Z=1 Fermi-model floor (module-Nuclear.jl) and
    #   the DFS self-interaction error for a one-electron system, respectively.
    #   VERIFIED, exact vs. computed: F(0)=1.0/1.00000 (exact); F(0.5)=0.885813/0.885817 (4.5e-6 rel.);
    #   F(1)=0.64/0.640009 (1.4e-5); F(2)=0.25/0.250012 (4.8e-5); F(4)=0.04/0.0400046 (1.2e-4);
    #   F(10)=0.00147920/0.00148009 (6e-4). The residual grows smoothly and monotonically with q -- exactly
    #   the expected signature of Z=1 relativistic corrections (~(Z*alpha)^2 ~ 5x10^-5) plus a small
    #   finite-nuclear-size contribution, not numerical noise. A genuine, independent, non-tabulated external
    #   check (no web/literature dependency, unlike Hubbell/Waasmaier-Kirfel which the user is separately
    #   trying to obtain).
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Ce-b-H1s-exact", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8797, AngularJ64(1//2), 0.0, 0.0, 0.0),
                            configs=[Configuration("1s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.NuclearField()),
                            propertySettings=[ FormFactor.Settings([0., 0.5, 1.0, 2.0, 4.0, 10.0], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch c: neutral Ne (Z=10), [Ne] closed shell -- comparison against Waasmaier & Kirfel (1995),
    #   Acta Cryst. A51, 416, examples/papers/1995.f0_WaasKirf.dat, a 5-Gaussian analytical fit
    #   f0(k) = c + sum_i a_i exp(-b_i k^2), k = sin(theta)/lambda [Angstrom^-1] (identical to the
    #   "x = sin(theta/2)/lambda" convention of Hubbell et al. (1975), since crystallographic theta IS the
    #   Bragg half-angle). Converted to JAC's atomic-unit momentum transfer via
    #   q [a.u.] = 4 pi a_0[Angstrom] k = 6.64983 k (standard elastic-scattering relation, independent of any
    #   possibly-mistranscribed formula in the source PDFs -- see project memory).
    #   Ne coefficients: a=[4.183749, 2.905726, 0.520513, 1.135641, 1.228065], c=0.025576,
    #   b=[8.175457, 3.252536, 0.063295, 21.813910, 0.224952]. Target f0(k=0)=9.999270 (WK's own ~7e-5
    #   relative fit residual vs the exact Z=10).
    #   VERIFIED, WK target vs. JAC-computed: q=0: 9.999270/10.00000 (7.3e-5, matches WK's own fit residual);
    #   q=0.665: 9.352189/9.31237 (0.43%); q=1.995: 6.078899/6.00742 (1.18%); q=3.990: 2.788834/2.80107
    #   (0.44%); q=6.650: 1.608401/1.61805 (0.60%). The small (<1.2%), q-dependent residual is physically
    #   explicable: WK's underlying tables ultimately derive from (largely non-relativistic) Hartree-Fock
    #   densities, whereas JAC uses a relativistic Dirac-Fock-Slater (DFS) mean field with local exchange --
    #   different many-electron treatments are expected to differ at this level, not a bug.
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Ce-c-Ne-WaasKirf", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne]")],
                            propertySettings=[ FormFactor.Settings([0., 0.665, 1.995, 3.9899, 6.6498], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch d: neutral Ar (Z=18), [Ar] closed shell -- comparison against BOTH Waasmaier & Kirfel (1995)
    #   AND Hubbell et al. (1975), examples/papers/1975.jpcrd-hubbel-form-factors.pdf, Table I (p. 498):
    #   F(x=0,Z) = 1.8000e+01 = 18 EXACTLY, an independent (different source, different method) confirmation
    #   of the same sum rule tested in branch a, now cross-checked against a second, completely separate
    #   40-year-older tabulation. Same q-conversion as branch c.
    #   Ar coefficients: a=[7.188004, 6.638454, 0.454180, 1.929593, 1.523654], c=0.265954,
    #   b=[0.956221, 15.339877, 15.339862, 39.043823, 0.062409]. Target f0(k=0)=17.999839 (WK) vs. Hubbell's
    #   exact 18.0 -- the two independent sources agree at q=0 to 5 significant figures.
    #   VERIFIED, WK target vs. JAC-computed: q=0: 17.999839/18.00000 (matches BOTH WK and Hubbell's exact
    #   18 to <1e-5); q=0.665: 16.298101/16.2623 (0.22%); q=1.995: 10.217116/10.2217 (0.045%);
    #   q=3.990: 6.878667/6.85767 (0.31%); q=6.650: 4.460082/4.43964 (0.46%). Even tighter agreement than
    #   Ne (all residuals <0.5%) -- consistent with Ar's larger, smoother closed-shell density being less
    #   sensitive to the DFS-vs-Hartree-Fock exchange-treatment difference discussed in branch c.
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Ce-d-Ar-WaasKirfHubbell", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(18.),
                            configs=[Configuration("[Ar]")],
                            propertySettings=[ FormFactor.Settings([0., 0.665, 1.995, 3.9899, 6.6498], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last successful:  31-Jul-2026
    # Branch e: Fe3+ ION (Z=26, 23 electrons), [Ar] 3d^5 -- comparison against Waasmaier & Kirfel (1995),
    #   which (unlike the older Cromer-Mann fit) explicitly covers IONS, not just neutral atoms -- the
    #   deliberate point of this branch, complementing branches c/d's neutral-atom checks.
    #   Fe3+ coefficients: a=[9.721638, 63.403847, 2.141347, 2.629274, 7.033846], c=-61.930725,
    #   b=[4.869297, 0.000293, 4.867602, 13.539076, 0.338520]. Target f0(k=0)=22.999227 (vs exact 23,
    #   WK's own ~3e-5 relative fit residual). Same q-conversion as branches c/d.
    #   VERIFIED, WK target vs. JAC-computed (level 1, 5/2+; all 5 levels of the 3d^5 configuration agree to
    #   5+ sig figs with each other, as expected -- the form factor depends on the shared radial density, not
    #   the term/J coupling): q=0: 22.999227/23.00000 (3.4e-5); q=0.665: 22.078566/22.0548 (0.11%);
    #   q=1.995: 16.725511/16.6681 (0.34%); q=3.990: 9.768990/9.78509 (0.17%); q=6.650: 6.559570/6.54864
    #   (0.17%). Excellent agreement (all <0.4%), comparable to the neutral Ar case -- confirms JAC's
    #   FormFactor machinery works well for an ION too, not just neutral atoms.
    setDefaults("print summary: open", "zzz-FormFactor.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Ce-e-Fe3plus-WaasKirf", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ar] 3d^5")],
                            propertySettings=[ FormFactor.Settings([0., 0.665, 1.995, 3.9899, 6.6498], true, LevelSelection())] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


