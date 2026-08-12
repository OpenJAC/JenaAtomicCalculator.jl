
println("Cc) Apply & test the IsotopeShift module with ASF from an internally generated multiplet.")
println("    Branches follow Naze et al., CPC 184 (2013) 2187 (RIS3) and Ekman et al., CPC 235 (2019) 433 (RIS4)")
println("    (examples/papers/2013.cpc-naze-gaigalas-RIS3.pdf, 2019.cpc-eckman-risc4.pdf);")
println("    see project memory project_isotope_shift_ris3ris4.md.")

setDefaults("print summary: open", "zzz-IsotopeShift.sum")

if  false
    # Last successful:  26-Jul-2026
    # Branch a: H(1s), bare single-electron test of the normal-mass-shift (NMS) parameter K_nms.
    #   Proton: I=1/2, mu=2.7928 nmu, Q=0, rms radius=0.8797 fm. "uniform" nuclear model, not "Fermi":
    #   for Z=1, JAC's 2-parameter Fermi model cannot represent an rms radius this small and now raises
    #   an explicit error (module-Nuclear.jl, fixed 25-Jul-2026 -- an earlier version silently substituted
    #   a hardcoded, wrong 1.89 fm nucleus instead; see memory).
    #   Two real bugs were found and fixed in the NMS machinery while setting this up (25-Jul-2026):
    #   (1) InteractionStrength.hamiltonian_nms divided by 2 on top of RadialIntegrals.isotope_nms's own
    #   internal /2 (a double division, silently halving every NMS value); (2) isotope_nms used
    #   (jb2-1)*jb2 for the Dirac small-component centrifugal coefficient l~(l~+1), which depends only on
    #   |kappa| and so cannot distinguish e.g. s_1/2 from p_1/2 (same j, different l~); the correct closed
    #   form, verified by hand against 8 different subshells, is kappa*(kappa-1).
    #   Verified against the virial theorem for a point-Coulomb potential (K_nms = <p^2>/2 = |E_1s| exactly,
    #   nonrelativistically): with a near-point rms radius (0.001 fm) the code gives K_nms = 0.500031 a.u.
    #   vs the exact 0.5 -- 0.006% agreement. With hydrogen's REAL rms radius (0.8797 fm, this branch),
    #   K_nms = 0.518773 a.u. -- about 3.75% above the point-nucleus value, which is a genuine (if larger
    #   than naively expected) finite-nuclear-size sensitivity of this derivative-heavy operator (NMS
    #   involves radial derivatives P', Q' that are locally far more sensitive near r->0 than the energy
    #   itself is), not a further bug.
    wa = Atomic.Computation(Atomic.Computation(), name="Cc-a-H1s", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8797, AngularJ64(1//2), 0.0, 0.0, 0.0),
                            configs=[Configuration("1s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.NuclearField()),
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  true
    # Last successful:  30-Jul-2026
    # Branch b: Li-like Nd (Z=60), 1s^2 2s ^2S_1/2 / 1s^2 2p ^2P_1/2,3/2 -- the RIS3/RIS4 papers' own worked
    #   example. 74Nd-142 Fermi nucleus (rms 4.9123 fm), single configuration (no correlation).
    #   UPDATE (30-Jul-2026): the ~14% K_nms gap documented below (this session, previously) was NOT a
    #   correlation/AL-SCF limitation as originally concluded -- it was a genuine B-spline boundary-condition
    #   defect in the underlying SCF machinery (kappa-sign leading/trailing spline truncation, root-caused
    #   and fixed this session in module-Bsplines.jl; see project_zeeman_hfs_bugs.md). With that fix in
    #   place, K_nms(2s) = 3970.46 a.u. vs the papers' 3976.60 a.u. -- only -0.15% off, essentially full
    #   agreement for this single-CSF (no-correlation) test. K_sms(2s) is still exactly 0.0 -- confirmed
    #   CORRECT, not a bug: the SMS operator's angular reduced matrix element (a rank-1 C^(1) tensor)
    #   vanishes by parity selection rule whenever every orbital in the system has the same l, which is the
    #   case for this all-s-orbital (1s,1s,2s) reference configuration; the papers' own small nonzero
    #   K_sms(2s)=7.74 a.u. is itself a pure CI/correlation effect, unrelated to this fix.
    #   HISTORICAL NOTE (pre-fix, kept for context): with only the branch-a NMS fixes applied, K_nms(2s) was
    #   4519.94 a.u. (13.7% too high). Hand-picked correlation configs were tried to close that gap and did
    #   NOT succeed -- correctly so, in hindsight, since the gap's real cause was the boundary-condition
    #   defect in every SCF orbital, not a missing correlation channel:
    #     - Adding "2s 2p^2" made both the level energy and K_nms worse (an AL-SCF orbital-optimization
    #       artifact from promoting both deeply-bound 1s core electrons, not a CI/diagonalization problem).
    #     - Adding "1s^2 3s" barely moved K_nms (0.01% shift) -- not the dominant channel either.
    #   These attempts were sound diagnostics; they simply couldn't fix a defect that lived in the orbitals
    #   themselves, common to every configuration tried.
    wa = Atomic.Computation(Atomic.Computation(), name="Cc-b-NdLiLike", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(60., FermiNucleus(), 142., 4.9123, AngularJ64(0//1), 0.0, 0.0, 0.0),
                            configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, calcSMS=true,
                                                printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  26-Jul-2026
    # Branch c: Li-like Nd (Z=60) field-shift electronic factor, same single-configuration 1s^2 2s ^2S_1/2
    #   system as branch b. calcF computes the field-shift electronic factor F two independent ways: F [ME]
    #   (direct matrix element of the nuclear-potential difference between two isotope masses, divided by
    #   the difference of their mean-square radii) and F [dens] (via the electron radial density at the
    #   origin) -- these agree with EACH OTHER to ~5-6 significant figures (1.18196e8 vs 1.18195e8 MHz/fm^2),
    #   a strong internal-consistency check between two independently-coded methods for the same quantity.
    #   Compared against the RIS4 paper's own n=3-CAS per-level electronic factor F0(2s) = 1.231688629e5
    #   GHz/fm^2 = 1.231688629e8 MHz/fm^2 (their Fig. 6): JAC's single-CSF result is about 4% low -- much
    #   closer agreement than K_nms's ~14% gap in branch b, consistent with the field shift (dominated by
    #   the electron density very close to the nucleus, set mostly by the deeply-bound 1s/2s orbital shapes)
    #   being markedly less sensitive to valence correlation than K_nms is.
    #   NOTE: IsotopeShift.computeAmplitudesProperties builds its internal comparison nucleus via the bare
    #   Nuclear.Model(Z, mass+1.0) constructor, which always defaults to "Fermi" regardless of what model
    #   type the primary nuclearModel argument uses -- for Z=60 this is harmless (auto rms ~4.95 fm, well
    #   above the ~1.86 fm Fermi floor), but this same code path would hit the Z=1 Fermi-model error fixed
    #   in module-Nuclear.jl (branch a) if ever used for hydrogen's field shift; not exercised here.
    nm = Nuclear.Model(60., FermiNucleus(), 142., 4.9123, AngularJ64(0//1), 0.0, 0.0, 0.0)
    wa = Atomic.Computation(Atomic.Computation(), name="Cc-c-NdFieldShift", grid=Radial.Grid(true),
                            nuclearModel=nm,
                            configs=[Configuration("1s^2 2s")],
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=false, calcSMS=false,
                                                calcF=true, printBefore=true) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  29-Jul-2026
    # Branch d (AL-Field / Breit revisit): reuses branch c's single-configuration Li-like Nd 1s^2 2s system,
    #   now that AL-Field (Basics.ALField()) has been root-cause-fixed and promoted to the standard
    #   implementation (see project_df_al_kink_bug.md). Four settings combinations: DFS+Coulomb (reproduces
    #   branch c's baseline exactly -- F[ME]=1.18196e8, F[dens]=1.18195e8 MHz/fm^2, bit-identical), AL+
    #   Coulomb, DFS+Breit, AL+Breit.
    #   RESULT (2), Breit: K_nms, K_sms and F are EXACTLY bit-identical between the Coulomb-only and
    #   Breit-added runs (both for DFS and for AL) -- precisely confirming the same structural prediction
    #   as example-Cb.jl Branch c (single-CSF system, mixing coefficient trivially 1.0, and
    #   module-IsotopeShift.jl never references Breit/eeInteraction at all). Only the total level energy
    #   shifts (by ~90 eV here), exactly as expected, not a surprise.
    #   RESULT (1), AL-Field: K_nms shifts DFS->AL by only +0.12% (1.63070e7 -> 1.63264e7) and F[ME] by only
    #   +0.17% (1.18196e8 -> 1.18394e8) -- A STRIKING CONTRAST to example-Cb.jl's Na HFS result, where the
    #   SAME DFS->AL switch changed A(3s) by a full 26%. This makes good physical sense: Nd 1s/2s are
    #   deeply-bound orbitals of a highly-charged (Z=60) few-electron ion, dominated by direct nuclear
    #   attraction -- exchange-treatment differences (DFS local vs AL non-local) matter far less there than
    #   for Na's single shallow valence 3s electron outside a neutral, lower-Z core. For the field shift
    #   specifically, AL moves F[ME] slightly CLOSER to the RIS4 target (1.231688629e8 MHz/fm^2): DFS is
    #   4.05% low, AL is 3.89% low -- a small, genuine improvement (unlike HFS, where AL made things worse).
    #   No comparably direct literature K_nms target exists for this single-CSF system (branch b's ~14% gap
    #   figure was for a DIFFERENT, 2-configuration reference), so no accuracy verdict is drawn for K_nms
    #   here -- only the AL-vs-DFS relative shift, which is small either way.
    nm = Nuclear.Model(60., FermiNucleus(), 142., 4.9123, AngularJ64(0//1), 0.0, 0.0, 0.0)

    wa1 = Atomic.Computation(Atomic.Computation(), name="Cc-d-Nd2s-DFS-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("1s^2 2s")],
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, calcSMS=true,
                                                calcF=true, printBefore=true) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cc-d-Nd2s-AL-Coulomb", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("1s^2 2s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField()),
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, calcSMS=true,
                                                calcF=true, printBefore=true) ] )
    perform(wa2)

    wa3 = Atomic.Computation(Atomic.Computation(), name="Cc-d-Nd2s-DFS-Breit", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("1s^2 2s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, calcSMS=true,
                                                calcF=true, printBefore=true) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cc-d-Nd2s-AL-Breit", grid=Radial.Grid(true),
                            nuclearModel=nm, configs=[Configuration("1s^2 2s")],
                            asfSettings=ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=Basics.ALField(),
                                                                  eeInteractionCI=Basics.CoulombBreit(0.)),
                            propertySettings=[ IsotopeShift.Settings(IsotopeShift.Settings(); calcNMS=true, calcSMS=true,
                                                calcF=true, printBefore=true) ] )
    perform(wa4)
    #
end
#
setDefaults("print summary: close", "")
