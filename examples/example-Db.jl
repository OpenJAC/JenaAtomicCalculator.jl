
println("Db) Apply & test the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoExcitation.sum")

if  true
    # Last successful:  28-Aug-2026 -- REPRODUCES EXACTLY. Energies 1.84650119 and 1.84657903 eV (record: 1.8464),
    # summed absorption oscillator strength f = 0.65048 (Coulomb) / 0.88455 (Babushkin) against the recorded
    # 0.6505 / 0.8845, still bracketing the literature f(2s-2p) = 0.7414.
    # THIS BRANCH IS UNAFFECTED BY THE 6-Aug OSCILLATOR-STRENGTH FIX and the two branches below are not; see there.
    # The reason is exact: 03cc54f corrected f by (2J_i+1)/2, and here J_i = 1/2, so the factor is 1.
    # Previously: 1-Aug-2026
    # Compute the photoexcitation cross sections for neutral lithium
    #
    #   REPORT: two fine-structure lines, 1->1 (2s J=1/2 -> 2p_1/2) and 1->2 (2s J=1/2 -> 2p_3/2), both at
    #   1.8464 eV -- matches the famous Li resonance ("D") line at 670.8 nm = 1.848 eV to <0.2%. Summed
    #   oscillator strength over both lines: f = 0.6505 (Coulomb) / 0.8845 (Babushkin), bracketing the
    #   well-known literature value f(2s-2p) = 0.7414 (e.g. Wiese, Smith & Miles 1969; NIST ASD) -- the
    #   ~30% Coulomb/Babushkin spread and the literature value falling neatly between them is the expected
    #   signature of a bare single-configuration (no core-polarization) calculation of an alkali resonance
    #   line, which is well known to be sensitive to core-polarization corrections not included here.
    #   NOTE: JAC itself flags the separate "photoexcitation cross sections" table (barns, as opposed to
    #   the oscillator-strength/resonance-strength table used above) as still under development at runtime
    #   -- a genuine, pre-existing incompleteness in PhotoExcitation (the printout wording was cleaned up
    #   1-Aug-2026, but the underlying cross-section formula itself was NOT touched/fixed); the
    #   oscillator-strength table is unaffected and was used for this report.
    pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(3.),
                            initialConfigs  = [Configuration("1s^2 2s")],
                            finalConfigs    = [Configuration("1s^2 2p")],
                            processSettings = pSettings )

    wb = perform(wa)
    #
elseif false
    # Last successful:  28-Aug-2026 -- THE RECORDED f-VALUES BELOW WERE PRE-FIX AND ARE NOW DOUBLED, CORRECTLY.
    # Today: f(1->2) = 1.8309e-02 (Coulomb) / 1.8578e-02 (Babushkin), f(1->4) = 7.6852e-01 / 7.8003e-01, at
    # 3.88624723e+03 and 3.90608068e+03 eV. The REPORT below quotes 0.0091/0.0086 and 0.3843/0.3886, i.e. almost
    # exactly half (ratio 1.99979 on the Coulomb resonance line).
    # THIS IS NOT A REGRESSION. `03cc54f` (6-Aug-2026) fixed the absorption oscillator strength, which "came out
    # too large by (2J_i+1)/2" -- so for J_i = 0 the old value was too SMALL by two, and the fix doubles it. The
    # branch is dated 1-Aug, five days BEFORE that commit, so its numbers are simply pre-fix. The fix was found by
    # checking f_ik against PhotoEmission's A_ki through the identity A_ki = 2 alpha^3 omega^2 (g_i/g_k) f_ik.
    # The new f(1->4) ~ 0.77 for the He-like resonance line is also the physically expected size.
    # Previously: 1-Aug-2026 (pre-fix numbers; the REPORT text below still quotes them)
    # Test of absorption f-values for 1s^2 --> 1s2p ^1P_1 resonance  in the helium isoelectronic sequence; cf. Table 8.1 in Section 8.1.a
    #
    #   REPORT: two E1 lines from the "1s 2p" multiplet's two J=1(-) levels: 1->2 (weak, f=0.0091(C)/
    #   0.0086(B), the relativistically-mixed 3P_1 "intercombination-like" line) and 1->4 (resonance line,
    #   the 1P_1 term) f=0.3843(C)/0.3886(B), gauge disagreement only ~1.3% -- both magnitude and gauge
    #   agreement are exactly what's expected for a strongly E1-allowed He-like resonance line at this Z
    #   (Ca, Z=20; known He-like 1s^2-1s2p ^1P_1 oscillator strengths run roughly 0.28-0.39 across the
    #   isoelectronic sequence). No external "Table 8.1" available this session to check the exact
    #   tabulated number; dated on this systematics + gauge-agreement basis (Rule 7). Same
    #   still-under-development cross-section caveat as branch 0 applies here too.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(20.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s 2p")],
                            processSettings = pSettings )
    wb = perform(wa)
    #
elseif false
    # Last successful:  28-Aug-2026 -- SAME PRE-FIX HALVING AS THE BRANCH ABOVE, and for the same reason.
    # Today: f(1->2) = 1.8213e-02 / 1.7616e-02, f(1->4) = 7.6512e-01 / 7.3954e-01 at 3.88546592e+03 and
    # 3.90529938e+03 eV, against the recorded 0.3826(C)/0.3684(B) for 1->4 -- ratio 2.0003.
    # See the branch above: `03cc54f`, 6-Aug-2026, and this branch is dated 1-Aug.
    # Previously: 1-Aug-2026 (pre-fix numbers)
    # Branch 2: Ca^18+ (He-like) 1s^2 --> 1s2p, WITHOUT vs. WITH the bi-orthogonal transformation
    #   (calcBiorthogonal), added 1-Aug-2026 together with the new calcBiorthogonal field on
    #   PhotoExcitation.Settings (analog to PhotoEmission.Settings, code in PhotoExcitation.computeLines).
    #   initialConfigs EXPANDED to also include "1s 2p" (unoccupied for the "1s^2" ground level) so the
    #   initial and final bases have matching orbital counts per kappa symmetry -- otherwise
    #   BiOrthogonal.computeTransformationMatrices errors out on the known differing-dimension limitation
    #   (0 vs 1 kappa=+1 orbitals), exactly as for the Fe X case in Da.jl branch 2. This expansion also
    #   surfaced and fixed a real, previously-undiscovered bug in
    #   BiOrthogonal.checkClosureUnderDeexcitation: it read the mutable global "standard subshell list"
    #   default instead of the basis's own (fixed) subshell list, causing a BoundsError whenever the
    #   initial and final CSF lists have different sizes (any excitation/ionization-type process is a
    #   natural trigger for this, not just a contrived edge case) -- fixed to always use basis.subshells.
    #   IMPORTANT: run WITHOUT/WITH in SEPARATE Julia sessions, as in all other biorthogonal branches.
    #
    #   REPORT: with the expanded config, the "1->2"/"1->4" oscillator strengths shift slightly from
    #   branch 1's simpler single-configuration numbers (f(1->4) = 0.3826(C)/0.3684(B) here vs.
    #   0.3843(C)/0.3886(B) there) -- expected, since the "1s^2" level's own orbital is now optimized
    #   jointly with the added "1s 2p" virtual rather than alone. The WITHOUT vs. WITH comparison itself,
    #   however, gave IDENTICAL results to the last printed digit for every line -- a genuine, well-
    #   understood null effect, not a bug: every kappa block in this minimal 2-configuration setup
    #   (kappa=-1: just 1s; kappa=+1: just 2p_1/2; kappa=-2: just 2p_3/2) is exactly 1-dimensional, so
    #   computeTransformationMatrices's LU decomposition is trivial (a 1x1 "matrix"), reducing Cleft/Cright
    #   to a pure orbital RESCALING with no other same-kappa orbital to mix with -- and a pure rescaling
    #   has, by construction, zero net effect on a properly normalized transition matrix element. This is
    #   consistent with (and a further confirmation of) the BiOrthogonal validation from the Da.jl session:
    #   the module correctly recognizes "nothing to correct" when there is no genuine multi-orbital
    #   non-orthogonality to act on. Seeing a NONTRIVIAL effect here would require a genuinely correlated
    #   setup with 2+ orbitals per kappa (e.g. an added 2s/3s core-polarization orbital) -- not attempted
    #   this session, a natural next step given that He-like resonance f-values ARE known to be
    #   correlation-sensitive (see branch 1's report).
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    useBiorthogonal = false
    pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin],
                                          printBefore=true, calcBiorthogonal=useBiorthogonal)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(20.),
                            initialConfigs  = [Configuration("1s^2"), Configuration("1s 2p")],
                            finalConfigs    = [Configuration("1s 2p")],
                            processSettings = pSettings )
    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
