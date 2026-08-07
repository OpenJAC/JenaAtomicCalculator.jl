
println("Dh) Apply & test the MultiPhotonTransition module: two-photon excitation and decay between bound levels.")

setDefaults("print summary: open", "zzz-MultiPhotonTransition.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

## REBUILT 06-Aug-2026, together with the module, which was renamed from MultiPhotonDeExcitation -- a name that
## described only half of what it does, since it has always carried absorption processes too.
##
## THE PREVIOUS VERSION OF THIS FILE COULD NOT RUN. Its branches 4 and 5 passed seven positional arguments to an
## eight-field Settings, and both used `mpGreen`, which was assigned in branch 3 -- a branch guarded by `false`.
## Since the branches form one if/elseif chain only ONE ever executes, so branch 2 likewise referred to `wb`
## from branch 1 without it existing. Every branch below is therefore SELF-CONTAINED: it builds what it needs.
##
## THE THREE PROCESSES, and what distinguishes them:
##  DIAGNOSTICS
##   a)  overview of the intermediate-state sum      -- which levels carry it, which denominators are dangerous
##  TWO-PHOTON EMISSION  (total energy fixed, shared between the photons -> continuous spectrum, energy sharings)
##   b)  H 2s -> 1s                                  -- THE ANCHOR: exact rate 8.2206 /s
##   c)  H-like Z-scan, Z = 1..54                    -- reproduces the relativistic suppression of the Z^6 law
##   d)  energy-differential spectrum                -- symmetry + sum rule; this is what found the sign bug
##   e)  He-like 1s2s 1S_0 -> 1s^2, He and Ar16+     -- J = 0 -> J = 0, so only K = 0
##   f)  H 2s -> 1s via a GREEN-FUNCTION sum         -- the alternative intermediate representation
##  TWO-PHOTON ABSORPTION  (photon energies imposed, transition driven -> cross sections, not rates)
##   g)  H 1s -> 2s, every polarization from one run -- J = 1/2 -> 1/2, so K = {0,1}: SEE THE WARNING BELOW
##   h)  Mg 3s^2 -> 3s4s via the 3s3p 1P_1 resonance -- a real atom; J = 0 -> J = 0, only K = 0
##   i)  bichromatic, two beams of different colour  -- scaffolded; fixes the single-beam convention
##  BEYOND TWO PHOTONS
##   j)  three-photon                                -- not implemented; fails with what would be needed
##
## ABSORPTION POLARIZATION OBSERVABLES -- two defects found and fixed 07-Aug-2026 (blocker A2).
##   (1) a typo: the K-weight read `2*Basics.twice(K) + 1`, and since twice(K) = 2K that is 4K+1, not 2K+1.
##       Harmless for K = 0 (both give 1), which is why Mg looked correct; for K = 1 it used 5 instead of 3.
##   (2) ODD K survived for a SINGLE BEAM, where it must vanish. Two photons from one beam are identical
##       bosons, so their polarization state must be SYMMETRIC; exchanging the multipoles in the 3-j carries
##       (-1)^(L1+L2+K) = (-1)^K for E1E1, so odd K is antisymmetric and forbidden. LINEAR light got this right
##       by accident -- its helicity sum is coherent, so the (+,-) and (-,+) terms cancel on their own -- while
##       the UNPOLARIZED sum is incoherent and accumulated |+-> and |-+> as separate states, neither of which is
##       individually symmetric. K = 1 was supplying 92 % of the unpolarized cross section for H 1s -> 2s.
##       The restriction applies to the MONOCHROMATIC scheme only: with two distinguishable beams the photons
##       are not identical and every K contributes -- one more reason to complete the bichromatic case.
##   VERIFIED: H 1s -> 2s now gives sigma_linear/sigma_unpolarized = 2.0000 in BOTH gauges, identical to Mg --
##   correct, since H 1s -> 2s is s -> s and therefore K = 0 only, exactly like a J = 0 -> 0 transition. Mg is
##   bit-identical to before (it was already K = 0 only), and its density-matrix factor 0.72 is unchanged.
##
##
## SELECTION is by `lineSelection` on level indices, and every branch carries an explicit `calcOverview` flag.
##
## UNITS, settled 06-Aug-2026 -- nothing had been defined or documented before. The generalized two-photon cross
## section is in cm^4 s, from W [1/s] = sigma^(2) * F^2 with F the photon flux density in photons cm^-2 s^-1,
## and is also reported in GM (1 GM = 1e-50 cm^4 s). For a SINGLE beam the two photons are indistinguishable and
## the combinatorial factor is a convention; it is fixed by requiring agreement with the bichromatic case,
## W = sigma^(2) * F_1 * F_2, in the limit omega_1 -> omega_2. See MultiPhotonTransition.Settings.
##
## STATUS, 07-Aug-2026. Two structural bugs were found and fixed (blockers A1, A2), both the same root cause --
## identical-boson exchange symmetry not enforced: the EMISSION exchange phase was (-1)^(K+J_f+J_i) where it
## must be (-1)^(L1+L2-K), and the single-beam ABSORPTION sum admitted ODD K, which is forbidden for two
## identical photons. Only branch d is dated "Last successful", and only for what it verifies -- the spectral
## SHAPE and the photon-exchange SYMMETRY, both independent of the open overall constant.
##
## THE ONE OPEN ISSUE THAT GATES EVERY EMISSION MAGNITUDE: a constant factor 6.679 (was 12.98 before A1),
## measured to 0.05 % over Z = 1..3 and therefore a prefactor, not a physics omission. It is NOT calibrated
## away -- see the note at the prefactor in module-MultiPhotonTransition-inc-2p-emission.jl. Until it is
## derived, emission MAGNITUDES from this file must not be quoted; SHAPES, SYMMETRIES, RATIOS and Z-TRENDS may
## be, since the constant divides out of all of them.
##
## A NOTE ON THE NUCLEAR MODEL. `Nuclear.Model(1.0)` defaults to Fermi and RAISES for hydrogen: with JAC's fixed
## skin thickness the 2-parameter Fermi distribution has a minimum rms radius near 1.86 fm, far above hydrogen's
## 0.88 fm. That is a good error -- it refuses rather than returning a silently wrong charge distribution -- so
## the hydrogen branches below use the point model, which is in any case the right choice for the benchmark:
## 8.2206 s^-1 is a point-nucleus number.
##
## REFERENCES for the verification stage:
##   Drake, PRA 34, 2871 (1986)                        -- 2-photon rates, H-like AND He-like; the standard table
##   Labzowsky et al.                                  -- A(2s->1s) = 8.2206 /s for neutral hydrogen
##   Amaro, Surzhykov, Parente, Indelicato & Santos,   -- relativistic 2s->1s, B-polynomials vs B-splines
##       J. Phys. A 44, 245302 (2011); arXiv:1104.4818
##   Fritzsche, Approximate Atomic Green Functions,    -- the intermediate-state summation itself
##       Molecules 26, 2660 (2021)

## THE RADIAL BOX MATTERS ENORMOUSLY HERE, and rbox = 10 a.u. -- the value this file used to carry -- is far
## too small. The intermediate np states are diffuse, and truncating them at 10 a.u. gave 357 /s for H 2s -> 1s
## where the same calculation at rbox = 80 gives 11.2 /s. That factor of 32 was a BOX artifact, not physics, and
## it very nearly led to "correcting" a prefactor that turned out to be right. See the convergence table in
## branch a and work/diag-2p-convergence.jl.
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 80.0)


##
## =====================================================================================================
##  DIAGNOSTICS -- run this first on any new system
## =====================================================================================================
if  false
    # Last visit:  06-Aug-2026
    #
    # --- Branch a: OVERVIEW of the intermediate-state sum for the H-like 2s -> 1s two-photon decay.
    #
    # RUN THIS FIRST on any new system. A two-photon amplitude is a SUM over intermediate levels |nu>, and that
    # sum is both the expensive part of the calculation and the delicate one. This branch builds everything,
    # ranks the intermediate levels by their actual contribution
    #
    #     | <f|O(mp2)|nu> <nu|O(mp1)|i> / (E_i + omega1 - E_nu) |
    #
    # and prints each energy denominator alongside -- then stops, without forming a single amplitude.
    #
    # IT ANSWERS THREE QUESTIONS that no configuration list can:
    #   * which intermediate levels actually carry the transition (usually a handful; the rest is noise);
    #   * whether any denominator is NEAR-RESONANT. `selfTolerance` removes an exactly singular term, but a
    #     merely small denominator is far more dangerous -- it yields a large finite number that looks like a
    #     result while the non-resonant perturbative treatment has quietly stopped applying;
    #   * whether the basis is big enough. If the largest contributions sit on the highest levels included, the
    #     sum is truncated and the answer is not converged, however smooth it looks.
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    interConfs  = [Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p")]
    interRep    = Representation("intermediate np levels", ni, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings()))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.EnergyDiffCs()], 4 ),
                        multipoles = [E1], gauges = [UseCoulomb],
                        intermediateStates = interMp,
                        calcOverview  = true,             ## <-- the whole point of this branch
                        lineSelection = LineSelection() )
    wo = Atomic.Computation(Atomic.Computation(), name="Dh-a: overview of the 2s -> 1s intermediate sum", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("2s")],
                            finalConfigs   = [Configuration("1s")],
                            processSettings= mpSettings )
    perform(wo)
    #
##
## =====================================================================================================
##  TWO-PHOTON EMISSION  (decay: the total energy E_i - E_f is FIXED and SHARED between the two photons,
##  so a line has a continuous spectrum and is resolved by energy sharings)
## =====================================================================================================
elseif  false
    # Last visit:  07-Aug-2026        ## NOT successful: the open prefactor, and this basis is not converged
    #
    # RE-RUN 07-Aug-2026 AFTER THE EXCHANGE-PHASE FIX (blocker A1). Babushkin = 1.461 /s with the 2p..8p basis
    # written below; the CONVERGED value (2p..45p, same box) is 1.231 /s, so this branch as configured is 19 %
    # off its own converged limit -- n_max = 8 is not enough, and the gauge ratio 0.153 here against 0.907 when
    # converged says so plainly. Against the exact 8.2206 /s the converged result is low by the open constant
    # 6.679 documented at the prefactor in module-MultiPhotonTransition-inc-2p-emission.jl.
    # NOT DATED SUCCESSFUL: neither the basis nor the absolute scale is settled. What this branch DID establish
    # is the exchange-phase bug -- K-resolving it showed K = 1 carrying 48 % of a transition where it must
    # vanish, which is how blocker A1 was found.
    #
    # --- Branch b: TWO-PHOTON EMISSION, 2s -> 1s in neutral hydrogen. THE ANCHOR OF THIS WHOLE FILE.
    #
    # A(2s -> 1s) = 8.2206 s^-1 is the single most-checked number in two-photon atomic physics: Breit & Teller,
    # Klarsfeld, Goldman & Drake, Labzowsky and others all agree on it. The transition is strictly forbidden for
    # ONE photon -- 2s and 1s have the SAME parity, so E1 cannot connect them -- which is precisely why the
    # two-photon channel governs the 2s lifetime.
    #
    # WHY IT IS THE RIGHT ANCHOR: one electron, so it is cheap; the answer is unambiguous to five digits; and
    # ANY error in the normalisation of the second-order amplitude shows up here immediately and undisguised.
    # No other branch of this file should be believed until this one reproduces 8.2206 s^-1.
    #
    # The intermediate states are the np levels: E1 twice, i.e. 2s -> np -> 1s. The sum formally runs over the
    # COMPLETE np spectrum INCLUDING THE CONTINUUM; truncating it at bound states up to some n_max is the
    # leading approximation here, and branch o is what shows whether that truncation has converged.
    #
    # THE INTERMEDIATE STATES MUST COME FROM THE SAME ONE-BODY HAMILTONIAN as the initial and final ones. This
    # is why `initialAsfSettings`, `finalAsfSettings` and the MeanFieldSettings below all carry
    # Basics.NuclearField(): -Z/r has no self-consistency, so it is identical for all three by construction, and
    # for a one-electron system it gives EXACT hydrogenic states. Gauge invariance in second-order perturbation
    # theory requires exactly this; a separately-converged DFS potential for the intermediate multiplet breaks
    # it systematically, and does so WITHOUT preventing either gauge from converging on its own.
    #
    # CONVERGENCE, measured 06/07-Aug-2026 (work/diag-2p-convergence.jl, diag-samepotential.jl, diag-nf-push.jl):
    #
    #     potential          n_max, rbox    Coulomb   Babushkin   Cou/Bab   Bab/8.2206
    #     DFS (mismatched)      8,  40       24.31      13.03      1.865      1.586
    #     DFS (mismatched)     12,  80       11.17      11.13      1.0035     1.354
    #     NuclearField          8,  40       22.09       9.726     2.272      1.183
    #     NuclearField         12,  80       15.09       8.693     1.736      1.057
    #     NuclearField         16, 120       14.90       8.661     1.720      1.054
    #     NuclearField         20, 160       16.09       8.774     1.834      1.067
    #
    # THREE LESSONS, and the first two each nearly cost something.
    #
    # (i) THE BOX DOMINATES. At rbox = 10 -- the value this file first carried -- the branch returned 357 /s, a
    #     factor 32 above the rbox = 80 result, purely from cutting off the diffuse np states. That artifact
    #     looked exactly like a wrong normalisation and would have been "fixed" by inserting a bogus factor into
    #     correct code. The prefactor 2pi*alpha^2*omega1*omega2/(2J_i+1) is RIGHT; do not touch it.
    #
    # (ii) THE "Cou/Bab = 1.0035" AT DFS n_max = 12 WAS A COINCIDENTAL CROSSING, NOT CONVERGENCE. Read as
    #     evidence of gauge invariance it supported a wrong story for a day. Two things exposed it: the DFS
    #     ratio is non-monotonic (2.08, 1.86, 1.0035, 1.17 -- it dips through 1 and returns), and at the SAME
    #     basis size NuclearField gives 1.736. Never trust a single point of a ratio; trust a trend.
    #
    # (iii) POTENTIAL CONSISTENCY IS WORTH A FACTOR OF SIX IN THE ERROR. Moving the intermediate states into the
    #     same Hamiltonian took the LENGTH gauge from 35 % high (DFS, 11.13) to 5.5 % high (8.66).
    #
    # WHERE IT NOW STANDS. Babushkin (length) is CONVERGED at ~8.7 /s, a stable 5-7 % above the exact 8.2206;
    # Coulomb (velocity) is stable at ~1.75x that. Both facts have ONE cause: the intermediate sum contains only
    # BOUND states. The velocity form weights the inner region and high-lying intermediates far more heavily
    # than the length form, so omitting the continuum hurts it disproportionately -- which is why the gauge split
    # persists instead of closing. The length gauge is the trustworthy number here.
    #
    # NOT DATED. 5.5 % with a 75 % gauge split is an understood discrepancy, not agreement. Closing it needs a
    # genuine continuum in the intermediate spectrum (a working Green expansion -- branch g -- or explicit
    # continuum orbitals), NOT more bound states: n_max = 12, 16, 20 already agree among themselves.
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    interConfs  = [Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p"),
                   Configuration("6p"), Configuration("7p"), Configuration("8p")]
    scf         = Basics.NuclearField()          ## the SAME one-body Hamiltonian for all three; see above
    asfA        = AsfSettings(AsfSettings(); scField = scf)
    interRep    = Representation("intermediate np levels", ni, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings(scf)))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.EnergyDiffCs()], 8 ),
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp,
                        calcOverview  = false,
                        lineSelection = LineSelection(), printBefore = true )
    wa = Atomic.Computation(Atomic.Computation(), name="Dh-b: 2s -> 1s two-photon decay of H", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("2s")],
                            finalConfigs   = [Configuration("1s")],
                            initialAsfSettings = asfA, finalAsfSettings = asfA,
                            processSettings= mpSettings )
    perform(wa)
    #
elseif  false
    # Last visit:  07-Aug-2026        ## NOT successful: n_max = 12 is far from converged
    #
    # RE-RUN 07-Aug-2026 AFTER THE EXCHANGE-PHASE FIX. Babushkin, and the deviation from the non-relativistic
    # Z^6 law normalised at Z = 1:
    #
    #     Z              1        2        6       18       36       54
    #     Bab [1/s]    1.341   85.82    6.246e4  4.482e7  2.715e9  2.826e10
    #     A/(A_1 Z^6)  1.000   0.9998   0.998    0.982    0.930    0.850     <- relativistic suppression
    #     Cou/Bab      0.302   0.302    0.301    0.293    0.275    0.245
    #
    # THE TREND IS THE POINT, and it is prefactor-independent: A/(A_1 Z^6) divides the unknown constant out, so
    # the relativistic suppression is measurable even while the absolute scale is open. Before blocker A1 the
    # same column read 1.000, 0.999, 0.992, 0.932, 0.769, 0.595 -- the old exchange phase was manufacturing a
    # spurious Z-dependence on top of the real one.
    #
    # NOT DATED: n_max = 12 gives a gauge ratio of 0.30, against 0.907 for the converged 2p..45p basis, so this
    # is nowhere near converged and the suppression curve above will move. Raise n_max (and keep rbox ~ 80/Z)
    # before quoting it against Drake, PRA 34, 2871 (1986).
    #
    # --- Branch c: the same 2s -> 1s decay along the H-LIKE ISOELECTRONIC SEQUENCE.
    #
    # A single number can be reproduced by a compensating pair of errors; a TREND is much harder to fake. The
    # non-relativistic rate scales as Z^6, so the interesting quantity is the deviation from that scaling --
    # which is where the relativistic and retardation corrections live. Drake, PRA 34, 2871 (1986) tabulates
    # exactly this for hydrogen-like ions and is the comparison.
    #
    # Watch the GAUGE AGREEMENT as Z grows: Coulomb and Babushkin must approach each other, and their spread is
    # an honest measure of how well the truncated intermediate sum is doing.
    ## THE BOX MUST SCALE AS 1/Z: a hydrogenic orbital shrinks like 1/Z, so a box fixed in absolute units would
    ## be far too small at Z = 1 and absurdly wasteful at Z = 54. rbox = 80/Z keeps the SAME physical coverage at
    ## every Z, which is what makes the trend meaningful rather than a box artifact -- and a box artifact is
    ## exactly what produced the bogus factor 32 at the start of this work.
    for  Z in [1.0, 2.0, 6.0, 18.0, 36.0, 54.0]
        println("\n", "="^110);   println("  H-like Z = $Z: 2s -> 1s two-photon decay");   println("="^110)
        ni          = Nuclear.Model(Z, Z < 2.0 ? "point" : "Fermi")
        scfZ        = Basics.NuclearField()          ## same one-body Hamiltonian for initial, final, intermediate
        asfZ        = AsfSettings(AsfSettings(); scField = scfZ)
        gridZ       = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 80.0/Z)
        interConfs  = [Configuration("$(n)p") for n = 2:12]
        interRep    = Representation("intermediate np", ni, gridZ, interConfs, MeanFieldMultiplet(MeanFieldSettings(scfZ)))
        interMp     = generate(interRep, output=true)["mean-field multiplet"]
        #
        mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                            scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                         MultiPhotonTransition.AbstractMultiPhotonProperty[
                                             MultiPhotonTransition.EnergyDiffCs()], 6 ),
                            multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                            intermediateStates = interMp, calcOverview = false,
                            lineSelection = LineSelection() )
        wb = Atomic.Computation(Atomic.Computation(), name="Dh-c: H-like Z = $Z", grid=gridZ,
                                nuclearModel   = ni,
                                initialConfigs = [Configuration("2s")],
                                finalConfigs   = [Configuration("1s")],
                                initialAsfSettings = asfZ, finalAsfSettings = asfZ,
                                processSettings= mpSettings )
        perform(wb)
    end
    #
elseif  false
    # Last successful:  07-Aug-2026
    #
    # DATED SUCCESSFUL for what it actually verifies, which is the SHAPE and the SYMMETRY -- both independent
    # of the open overall constant. Re-run 07-Aug-2026 after the exchange-phase fix (Babushkin):
    #
    #     0.201  0.864  1.573  2.084  2.367  2.481 | 2.481  2.367  2.084  1.573  0.864  0.201
    #
    # THREE CHECKS PASS:
    #   * PHOTON-EXCHANGE SYMMETRY is exact to all printed digits -- every mirror pair identical. The two
    #     emitted photons are indistinguishable, so this is forced, and it is the check that first exposed the
    #     K-sum-inside-the-modulus bug (06-Aug) when it failed by 60-90 %.
    #   * THE SHAPE is now a smooth dome PEAKED AT THE CENTRE, as the true H 2s -> 1s spectrum is. Before the
    #     exchange-phase fix it had a MINIMUM at the centre, and before the denominator-sign fix it carried a
    #     spurious spike at omega1 = 3.225 eV from poles that cannot exist.
    #   * THE INTEGRAL over the Gauss-Legendre sharings reproduces the total rate of branch b.
    #
    # WHAT IS *NOT* CLAIMED: the absolute scale, which carries the open constant 6.679. This branch is dated for
    # the spectral shape and the symmetry, not for the magnitude -- see the prefactor note in the module.
    #
    # --- Branch d: the ENERGY-DIFFERENTIAL SPECTRUM of one two-photon line.
    #
    # Unlike a one-photon line, a two-photon decay has a CONTINUOUS spectrum: only the sum omega1 + omega2 is
    # fixed, and the two photons share it in any proportion. That spectrum is the characteristic signature of
    # the process, and it is what makes two-photon decay visible as a continuum in astrophysical plasmas.
    #
    # TWO CHECKS FOR FREE, needing no literature at all, which is why this deserves its own branch:
    #   * the spectrum must be SYMMETRIC about omega1 = omega2 = (E_i - E_f)/2. The two photons are
    #     indistinguishable, so the differential rate cannot tell them apart; any asymmetry is a bug in the
    #     handling of the two multipole orderings, not physics.
    #   * its INTEGRAL over the sharings must reproduce the total rate of branch a. The sharings sit at
    #     Gauss-Legendre nodes precisely so that the weighted sum performs that integral exactly.
    # Together these test the sharing machinery and the amplitude symmetry independently of any normalisation.
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    interConfs  = [Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p"),
                   Configuration("6p"), Configuration("7p"), Configuration("8p")]
    scfC        = Basics.NuclearField()
    asfC        = AsfSettings(AsfSettings(); scField = scfC)
    interRep    = Representation("intermediate np levels", ni, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings(scfC)))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.EnergyDiffCs()], 12 ),  ## many sharings: resolve the shape
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp, calcOverview = false,
                        lineSelection = LineSelection() )
    wc = Atomic.Computation(Atomic.Computation(), name="Dh-d: differential 2-photon spectrum of H 2s -> 1s", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("2s")],
                            finalConfigs   = [Configuration("1s")],
                            initialAsfSettings = asfC, finalAsfSettings = asfC,
                            processSettings= mpSettings )
    perform(wc)
    #
elseif  false
    # Last visit:  07-Aug-2026        ## NOT successful: correlation-limited, plus the open prefactor
    #
    # RE-RUN 07-Aug-2026 AFTER THE EXCHANGE-PHASE FIX: He Bab = 2.670 /s, Ar16+ Bab = 1.543e7 /s.
    # Scaling by the open constant 6.679 gives He ~ 17.8 /s against the accepted 51.3 /s -- still a factor ~3
    # low, which is the CORRELATION deficiency this branch was built to expose, not a machinery problem: the
    # f-value check (work/diag-ca-fvalue.jl) puts Mg and Ca one-photon oscillator strengths 25-45 % from their
    # measured values at this level of CI, and a two-photon amplitude carries two such matrix elements.
    # NOT DATED: two independent open items (correlation, prefactor) sit between this and a claim.
    #
    # --- Branch e: HELIUM-LIKE two-photon emission, 1s2s 1S_0 -> 1s^2 1S_0, for He (Z=2) and Ar16+ (Z=18).
    #
    # WHY THIS IS THE RIGHT SECOND SYSTEM, and the reason is sharper than "another benchmark".
    #
    # J_i = J_f = 0, SO Klist = oplus(0,0) = {0} -- ONE value of K. That matters three times over:
    #   (1) every K-coupling complication disappears. The interference bug fixed on 06-Aug-2026 (the K-sum
    #       taken inside the modulus) cannot even ARISE here, so if a discrepancy remains it lies in the radial
    #       or normalisation part, cleanly separated from the angular part;
    #   (2) K = 0 IS the scalar part, so this branch validates precisely the quantity the ABSORPTION side reports
    #       as TotalAlpha0. Doing this first de-risks branch d rather than competing with it;
    #   (3) Drake, PRA 34, 2871 (1986) tabulates hydrogen-like AND helium-like, so one reference covers both.
    #
    # THE DISCRIMINATING VALUE, which is the real point. Branch a converged to a gauge ratio of 1.0035 but landed
    # at 11.2 /s against 8.2206, and that 1.35 excess was ATTRIBUTED to the missing continuum in the intermediate
    # sum. That is a hypothesis, not a result. This branch tests it on an independent system:
    #     * a SIMILAR ~1.35 excess here  -> the continuum diagnosis is confirmed, and building continuum
    #       machinery (a working Green route) is clearly the right investment;
    #     * a DIFFERENT factor           -> the diagnosis is wrong, and it is far better to learn that on a cheap
    #       two-electron system than after building that machinery on a false premise.
    #
    # BENCHMARKS:
    #     He   (Z=2):   51.3 /s      Drake, Victor & Dalgarno (1969); also Drake PRA 34, 2871 (1986)
    #     Ar16+ (Z=18): see Drake 1986, which tabulates to Z = 36 -- NOT quoted here from memory.
    # One expectation that needs no table: as Z grows the 2s electron sees an increasingly unscreened nucleus, so
    # the He-like rate must APPROACH the H-like rate at the same Z from below. For Z = 18 the H-like value is
    # 8.2206 * 18^6 = 2.80e8 /s, so the Ar16+ result should be somewhat under that -- a sanity bound, not a
    # benchmark.
    #
    # RESULT, 07-Aug-2026 -- THE DISCRIMINATING TEST ANSWERED, BUT NOT AS PREDICTED:
    #
    #     system        omega [eV]    Coulomb     Babushkin    Cou/Bab    reference
    #     He   (Z=2)      19.77        6.74         5.53        1.22      51.3 /s     -> 8-9x LOW
    #     Ar16+ (Z=18)   3126.4        2.72e8       1.47e8      1.85      H-like 2.80e8 -> 0.97
    #
    # THE CONTINUUM HYPOTHESIS IS NOT CONFIRMED. Branch a (H-like) came out 1.35x HIGH; this comes out 8-9x LOW.
    # Opposite directions, so a single cause cannot explain both -- which is exactly what this branch was built
    # to find out, and precisely why it was worth running BEFORE investing in continuum machinery.
    #
    # BUT THE TEST IS CONFOUNDED, and that must be said rather than over-read:
    #   * HE IS CORRELATION-LIMITED, not machinery-limited. Its 1s2s 1S_0 is genuinely correlated and a
    #     single-configuration treatment underestimates it. Read as an effective charge: 51.3 /s corresponds to
    #     Z_eff = 1.36, our 6.0 /s to Z_eff = 0.95, i.e. the screening is overestimated. The energy says the
    #     same -- 19.77 eV against the true 20.62 eV, 4 % low.
    #   * AR16+ IS THE CLEAN POINT, because at Z = 18 the system is nearly hydrogenic and correlation barely
    #     matters. Energy 3126.4 eV vs ~3123 eV is 0.1 %, and the Coulomb rate sits at 0.97 of the H-like
    #     expectation 8.2206*18^6 = 2.80e8 -- just below it, where it should be.
    #   * NEITHER IS CONVERGED. Gauge ratios 1.22 and 1.85 are both far worse than branch a's 1.0035, so
    #     1snp to n = 10 is too small a basis in both cases. Converge AR FIRST: it is the case where a residual
    #     discrepancy can be attributed to the method rather than to correlation.
    #
    # SO THE MACHINERY LOOKS RIGHT WHERE CORRELATION IS UNIMPORTANT, and the "discrepancies" have different
    # causes -- correlation at low Z, basis truncation everywhere, and something still unexplained in H-like.
    # No date is written until Ar16+ is converged in gauge and compared with Drake's actual tabulated value.
    #
    # SETUP NOTES. 1s2s produces BOTH 3S_1 and 1S_0; only the singlet is wanted, and it must not be confused with
    # 2 3S_1 -> 1 1S_0, which is a different process dominated by M1 rather than 2E1. The selection is therefore
    # made on SYMMETRY, 0+ -> 0+, which picks the singlet robustly however the levels happen to be ordered.
    # The intermediate states are 1snp 1P_1 (J = 1, odd), reached by E1 from the 1S_0.
    for  (Z, rbox, nmax) in [(2.0, 60.0, 10), (18.0, 20.0, 10)]
        println("\n", "="^110);   println("  He-like Z = $Z: 1s2s 1S_0 -> 1s^2 1S_0 two-photon decay");   println("="^110)
        ni          = Nuclear.Model(Z, Z < 2.5 ? "point" : "Fermi")
        gridH       = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = rbox)
        interConfs  = [Configuration("1s $(n)p") for n = 2:nmax]
        interRep    = Representation("1snp intermediate levels", ni, gridH, interConfs,
                                     MeanFieldMultiplet(MeanFieldSettings()))
        interMp     = generate(interRep, output=true)["mean-field multiplet"]
        #
        mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                            scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                         MultiPhotonTransition.AbstractMultiPhotonProperty[
                                             MultiPhotonTransition.EnergyDiffCs()], 8 ),
                            multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                            intermediateStates = interMp, calcOverview = false,
                            ## 0+ -> 0+ only: the 1S_0 singlet, NOT the 3S_1 triplet (a different, M1 process)
                            lineSelection = LineSelection(true,
                                symmetryPairs = [(LevelSymmetry(0, Basics.plus), LevelSymmetry(0, Basics.plus))]),
                            printBefore = true )
        wh = Atomic.Computation(Atomic.Computation(), name="Dh-e: He-like Z = $Z", grid=gridH,
                                nuclearModel   = ni,
                                initialConfigs = [Configuration("1s 2s")],
                                finalConfigs   = [Configuration("1s^2")],
                                processSettings= mpSettings )
        perform(wh)
    end
    #
elseif  false
    # Last visit:  07-Aug-2026        ## NOT successful: the Green route remains worse than the multiplet one
    #
    # RE-RUN 07-Aug-2026 AFTER THE EXCHANGE-PHASE FIX: Bab = 1.144 /s, Cou/Bab = 0.349, against the multiplet
    # route's converged 1.231 /s at Cou/Bab = 0.907. So the Green expansion is still the WORSE representation
    # here -- its gauge ratio is far from 1 -- even though giving it a matched potential (the scField added to
    # AtomicState.GreenSettings on 07-Aug-2026) improved it by a factor 15 earlier. The damped CI space simply
    # does not span the spectrum this transition needs; see the notes in the branch body.
    # NOT DATED.
    #
    # --- Branch f: the same H 2s -> 1s decay, but with a GREEN-FUNCTION intermediate representation.
    #
    # THIS IS THE HONEST NEXT STEP AFTER BRANCH a, and the reason is specific rather than general. Branch a,
    # once its radial box was large enough, converged to a gauge ratio of 1.0035 -- the amplitude is
    # gauge-invariant, so the second-order machinery and the prefactor are right -- yet it lands at 11.2 /s
    # against the accepted 8.2206 /s, a residual of 1.35.
    #
    # THAT RESIDUAL HAS ONE IDENTIFIABLE CAUSE. Branch a sums over BOUND np states only. The exact second-order
    # sum runs over the COMPLETE spectrum, and for 2s -> 1s the continuum contribution is both substantial and
    # of the sign that LOWERS the rate -- consistent with 11.2 being high rather than low. No number of extra
    # bound states can fix that: 12p, 16p and beyond are already converged among themselves (see branch a), and
    # what is missing is not a further bound state but a different part of the spectrum.
    #
    # A GREEN EXPANSION IS THE RIGHT INSTRUMENT because it represents the resolvent rather than enumerating
    # eigenstates, so the continuum enters through the pseudo-states of a damped, finite-box CI space instead of
    # being truncated away. `Settings.intermediateStates` accepts either representation and dispatches on the
    # type through MultiPhotonTransition.intermediateLevels, so NOTHING in the amplitude changes between branch a
    # and this one -- which is exactly what makes the comparison clean: any change in the answer comes from the
    # intermediate spectrum alone.
    #
    # WHAT TO WATCH, in this order:
    #   (1) does the gauge ratio STAY near 1? If a Green expansion breaks the gauge invariance that branch a
    #       achieved, the expansion itself is at fault, not the physics.
    #   (2) does the rate MOVE DOWNWARD from 11.2 toward 8.2206? Direction first, magnitude second.
    #   (3) is it stable against nMax and dampingTau? A result that drifts with the damping parameter is a
    #       numerical artifact of the expansion, not a continuum contribution.
    # Only if all three hold does 8.2206 s^-1 become a claim rather than a coincidence -- and only then does any
    # branch of this file earn a "Last successful" date.
    #
    # RESULT, 06-Aug-2026 -- AN HONEST NEGATIVE. The branch runs end to end, so the Green route is wired and
    # dispatched correctly (Settings.intermediateStates accepts the channels and MultiPhotonTransition.
    # intermediateLevels selects them by symmetry, with nothing else in the amplitude changed). But:
    #
    #     route                       Coulomb    Babushkin    Cou/Bab
    #     mean-field  (branch a)       11.17       11.13       1.0035
    #     Green       (this branch)    20.78       97.68       0.213
    #
    # CRITERION (1) FAILS. The gauge invariance that branch a had achieved is destroyed, so by the rule set out
    # above the expansion itself is at fault and the rates must NOT be interpreted -- neither the direction of
    # the move nor its size means anything while the two gauges disagree by a factor of five.
    #
    # WHAT TO TRY NEXT, in order of suspicion:
    #   * dampingTau = 0.01: the damping regularises the resolvent but also distorts it. Scan it, and treat any
    #     result that drifts with it as an artifact (criterion 3), not as a continuum contribution.
    #   * the reference space: one configuration (2p) with DeExciteSingleElectron may simply not generate a
    #     pseudo-spectrum that spans the region the 2s -> 1s sum needs.
    #   * the approach: DampedSpaceCI vs CoreSpaceCI vs SingleCSFwithoutCI. For a ONE-electron system the CI is
    #     trivial, so a difference between these would itself be diagnostic.
    #   * nMax = 12 and lValues = [1]: too small a pseudo-spectrum would leave the same gap as branch a rather
    #     than closing it, but it should not break the gauge ratio -- so if the ratio stays bad as nMax grows,
    #     the problem is in the expansion's construction, not its size.
    #
    # ROUTE 1 APPLIED, 07-Aug-2026 -- THE CAUSE IS CONFIRMED, BUT IT IS NOT THE CURE:
    #
    #     Green space potential        Coulomb   Babushkin   Cou/Bab   Bab/8.2206
    #     DFS (as this branch had it)   0.335      3.396      0.0987      0.413
    #     NuclearField (matched)        6.416      4.399      1.459       0.535
    #
    # MATCHING THE POTENTIAL MOVED THE GAUGE RATIO BY A FACTOR OF 15, with nothing else changed. Together with
    # the same fix improving branch a's LENGTH gauge from 35 % to 5.5 % error, two INDEPENDENT code paths now
    # confirm that a mismatched intermediate potential was the dominant defect. That is established, not
    # conjectured.
    #
    # BUT THE GREEN EXPANSION IS STILL NOT A GOOD CONTINUUM: Babushkin 4.40 against 8.2206, i.e. 47 % LOW, where
    # the bound-state route with the SAME matched potential gives 8.66, 5.5 % high. Note its gauge ratio is
    # actually BETTER (1.46 vs 1.72) while its absolute value is far worse -- the damped CI space is more
    # internally consistent but spans less of the spectrum that matters. With nMax = 12, lValues = [1], ONE
    # reference configuration and dampingTau = 0.01 it is simply not complete.
    #
    # SO THE REAL FIX IS COMPLETENESS, not this expansion: build the intermediate multiplet from the FULL
    # B-spline eigenspectrum of the same potential, bound AND positive-energy pseudo-states. Summing over a
    # complete B-spline eigenbasis in a box is mathematically IDENTICAL to solving the Dalgarno-Lewis
    # inhomogeneous equation (H - E_i - omega)|chi> = O|i>, so it is not an approximation to the exact method --
    # it IS the exact method written as a sum, and completeness makes gauge invariance a pass/fail test rather
    # than a judgement call. Bsplines.generateOrbitals already takes an explicit potential and
    # Bsplines.findPositiveBranchStart already guards the Dirac negative-energy sea, which is not optional here.
    ni          = Nuclear.Model(1.0, "point")
    ## The intermediate states of a 2s -> 1s two-photon decay are p levels: E1 twice, so J^P = 1/2- and 3/2-.
    levelSyms   = [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(3//2, Basics.minus)]
    ## GreenSettings(nMax, lValues, dampingTau, printBefore, levelSelection, scField). The SIXTH argument was
    ## added on 07-Aug-2026 for exactly this branch: until then generate(GreenExpansion, ...) hardcoded
    ## AsfSettings() and Basics.DFSField(1.0), so a Green space could NEVER be matched to the potential of the
    ## computation consuming it -- and for a second-order intermediate spectrum that mismatch breaks gauge
    ## invariance systematically. (The zero-argument GreenSettings() was also broken, calling Settings(...);
    ## fixed at the same time.)
    scfG        = Basics.NuclearField()
    asfG        = AsfSettings(AsfSettings(); scField = scfG)
    greenSet    = GreenSettings(12, [1], 0.01, false, LevelSelection(), scfG)
    greenRep    = Representation("Green expansion for the np intermediate spectrum", ni, grid,
                                 [Configuration("2p")],
                                 GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(),
                                                 levelSyms, 1, greenSet) )
    greenChs    = generate(greenRep, output=true)["Green channels"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.EnergyDiffCs()], 8 ),
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = greenChs,      ## <-- the Green route, not a mean-field multiplet
                        calcOverview  = false,
                        lineSelection = LineSelection(), printBefore = true )
    wg = Atomic.Computation(Atomic.Computation(), name="Dh-f: 2s -> 1s via a Green-function intermediate sum", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("2s")],
                            finalConfigs   = [Configuration("1s")],
                            initialAsfSettings = asfG, finalAsfSettings = asfG,
                            processSettings= mpSettings )
    perform(wg)
    #
##
## =====================================================================================================
##  TWO-PHOTON ABSORPTION  (excitation: the photon energies are IMPOSED, and the transition is driven
##  rather than spontaneous; the observables are cross sections, not rates)
## =====================================================================================================
elseif  true
    # Last visit:  06-Aug-2026
    #
    # --- Branch g: TWO-PHOTON ABSORPTION, monochromatic, with EVERY polarization from ONE calculation.
    #
    # The 1s -> 2s excitation is the time-reverse of branch a and equally forbidden for one photon. Both photons
    # come from the same beam and each carries omega = (E_f - E_i)/2.
    #
    # ONE CALCULATION, FIVE OBSERVABLES. The amplitude is decomposed into irreducible parts of rank K, and every
    # polarization observable is a fixed algebraic combination of the SAME |M_K|^2 -- so requesting all of them
    # costs nothing beyond requesting one:
    #     TotalAlpha0            the K = 0 part alone; the ONLY surviving term when J_i = J_f = 0
    #     TotalCsLinear          Stokes (1,0,0)
    #     TotalCsRightCircular   Stokes (0,0,1)
    #     TotalCsUnpolarized     Stokes (0,0,0)
    #     TotalCsDensityMatrix   the user's own Stokes vector, of which the three above are special values
    #
    # THEIR RATIOS ARE A CHECK IN THEMSELVES. Angular algebra alone fixes them, independently of any overall
    # normalisation -- so this branch can be tested even before an absolute two-photon absorption benchmark for
    # a general atom is available, which is the one reference still missing.
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    interConfs  = [Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p"),
                   Configuration("6p"), Configuration("7p")]
    interRep    = Representation("intermediate np levels", ni, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings()))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonAbsorptionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.TotalAlpha0(), MultiPhotonTransition.TotalCsLinear(),
                                         MultiPhotonTransition.TotalCsRightCircular(),
                                         MultiPhotonTransition.TotalCsUnpolarized(),
                                         MultiPhotonTransition.TotalCsDensityMatrix()] ),
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp, calcOverview = false,
                        stokes = ExpStokes(0.6, 0.0, 0.8),      ## partially polarized, for TotalCsDensityMatrix
                        lineSelection = LineSelection(), printBefore = true )
    wd = Atomic.Computation(Atomic.Computation(), name="Dh-g: 1s -> 2s two-photon absorption of H", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2s")],
                            processSettings= mpSettings )
    perform(wd)
    #
elseif  false
    # Last visit:  07-Aug-2026
    #
    # --- Branch h: TWO-PHOTON ABSORPTION in MAGNESIUM, 3s^2 1S_0 -> 3s4s 1S_0, via the 3s3p 1P_1 resonance.
    #
    # THE FIRST REAL-ATOM CASE, and the first ABSORPTION case ever run in this module. It was chosen because it
    # fits four constraints at once, and because hydrogen turned out to be a poor model of what this module is
    # for -- see "WHY NOT HYDROGEN" below.
    #
    #   (1) J = 0 -> J = 0, so Klist = oplus(0,0) = {0}. ONE value of K, hence no K-coupling at all, and K = 0
    #       IS the scalar part -- so TotalAlpha0 is the entire answer rather than one component of it.
    #   (2) ONE DOMINANT INTERMEDIATE, the 3s3p 1P_1 resonance: strong, and close in energy. This is the
    #       multiplet approach at its best, and calcOverview should show 3s3p standing far above 3s4p, 3s5p.
    #   (3) ABSORPTION, which is the direction of interest and which had never been executed.
    #   (4) A REAL ATOM, not a noble gas -- the gap in the literature that motivates this work.
    #
    # WHY NOT HYDROGEN, and this is the lesson of 06/07-Aug-2026. For H 2s -> 1s the CONTINUUM carries 43.5 % of
    # the dipole sum rule (the Thomas-Reiche-Kuhn bound part is 0.565, measured 0.5599 in work/diag-trk.jl):
    # there is no intermediate resonance at all, the strength is smeared over the whole spectrum, and a
    # bound-state multiplet can never represent it. That made hydrogen an EXCELLENT bug-finder -- it exposed a
    # factor-32 radial-box artifact, a K-interference that broke photon-exchange symmetry, and a potential
    # mismatch worth a factor 6 -- but a poor model of the intended use. For a real atom the intermediate states
    # are genuine, well-separated levels and a handful dominate, so no pseudo-continuum is needed.
    #
    # WHAT CAN AND CANNOT BE CHECKED HERE. Mg 3s^2 and 3s4s are correlated states, so the ABSOLUTE cross section
    # is limited by the CI, exactly as helium showed (branch h, 8-9x low from correlation alone). But the
    # POLARIZATION RATIOS -- linear vs right-circular vs unpolarized -- are fixed by angular algebra and are
    # INDEPENDENT of any overall normalisation or of the quality of the wave functions. They are therefore a
    # real check available today, without the absolute two-photon absorption benchmark that is still missing.
    # For a 0 -> 0 transition through a J = 1 intermediate the ratios are strongly constrained; in particular
    # RIGHT-CIRCULAR light cannot drive 0 -> 0 by two E1 photons at all (the two photons would have to supply
    # 2 units of angular momentum along the beam, which no J = 0 -> J = 0 transition can absorb), so
    # TotalCsRightCircular MUST come out at or near zero. That is a sharp, parameter-free test.
    #
    # ON THE POTENTIAL: the same scField is used for the initial, final AND intermediate multiplets, following
    # the 07-Aug-2026 finding that a mismatch there breaks gauge invariance systematically. Note the mismatch is
    # far milder here than in hydrogen -- 3s^2, 3s3p and 3s4s share the same [Ne] core, so their self-consistent
    # potentials are close, whereas DFS-vs-NuclearField for a single electron was a large difference.
    ni          = Nuclear.Model(12.0)
    gridM       = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 40.0)
    scfM        = Basics.DFSField()
    asfM        = AsfSettings(AsfSettings(); scField = scfM)
    ## the intermediate 3snp 1P_1 levels; 3s3p is expected to dominate
    ## MODERATELY ENLARGED, 07-Aug-2026: 3snp up to n = 7 plus 3p^2 correlation. Not a pseudo-continuum -- for a
    ## real atom the sum is carried by a few genuine resonances, which is precisely why the multiplet approach
    ## suits real atoms and did not suit hydrogen. calcOverview should show 3s3p dominating by a wide margin.
    interConfs  = [Configuration("[Ne] 3s 3p"), Configuration("[Ne] 3s 4p"), Configuration("[Ne] 3s 5p"),
                   Configuration("[Ne] 3s 6p"), Configuration("[Ne] 3s 7p"), Configuration("[Ne] 3p^2")]
    interRep    = Representation("3snp intermediate levels", ni, gridM, interConfs,
                                 MeanFieldMultiplet(MeanFieldSettings(scfM)))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonAbsorptionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.TotalAlpha0(),
                                         MultiPhotonTransition.TotalCsLinear(),
                                         MultiPhotonTransition.TotalCsRightCircular(),
                                         MultiPhotonTransition.TotalCsUnpolarized(),
                                         MultiPhotonTransition.TotalCsDensityMatrix()] ),
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp, calcOverview = false,
                        stokes = ExpStokes(0.6, 0.0, 0.8),   ## partially polarized, for TotalCsDensityMatrix
                        ## 0+ -> 0+ only: the 1S_0 of 3s4s, NOT its 3S_1
                        lineSelection = LineSelection(true,
                            symmetryPairs = [(LevelSymmetry(0, Basics.plus), LevelSymmetry(0, Basics.plus))]),
                        printBefore = true )
    wi = Atomic.Computation(Atomic.Computation(), name="Dh-h: Mg 3s^2 -> 3s4s two-photon absorption", grid=gridM,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("[Ne] 3s^2")],
                            finalConfigs   = [Configuration("[Ne] 3s 4s")],
                            initialAsfSettings = asfM, finalAsfSettings = asfM,
                            processSettings= mpSettings )
    perform(wi)
    #
elseif  false
    # Last visit:  06-Aug-2026
    #
    # --- Branch i: TWO-PHOTON ABSORPTION, BICHROMATIC -- two beams of different frequency.
    #
    # THE SCHEME EXISTED BUT HAD NO computeLines AT ALL until 06-Aug-2026: it was declared, documented in the
    # abstract type and given a default constructor, yet had no `-inc-` file, so selecting it died in a bare
    # MethodError with nothing to indicate why. It now has its own file and reports precisely what is missing.
    #
    # WHY IT MATTERS BEYOND ITS OWN PHYSICS. With two distinguishable beams the rate is W = sigma^(2) * F_1 * F_2
    # with NO combinatorial factor to argue about, whereas the single-beam case involves indistinguishable
    # photons and therefore a convention. Requiring the two to agree in the limit omega_1 -> omega_2 is what
    # FIXES that convention instead of leaving it to be guessed -- exactly the class of silent factor that has
    # cost this project real time before.
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    interConfs  = [Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p")]
    interRep    = Representation("intermediate np levels", ni, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings()))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.TotalCsUnpolarized()], 3.4 ),  ## omegaLess
                        multipoles = [E1], gauges = [UseCoulomb],
                        intermediateStates = interMp, calcOverview = false,
                        lineSelection = LineSelection() )
    we = Atomic.Computation(Atomic.Computation(), name="Dh-i: bichromatic two-photon absorption of H", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2s")],
                            processSettings= mpSettings )
    perform(we)
    #
##
## =====================================================================================================
##  BEYOND TWO PHOTONS -- not implemented; the schemes exist and say so
## =====================================================================================================
elseif  false
    # Last visit:  06-Aug-2026
    #
    # --- Branch j: THREE-PHOTON -- the scheme exists and says clearly that it is not implemented.
    #
    # This branch is deliberately kept, and deliberately expected to FAIL. A scheme that is declared but has no
    # method behind it is exactly what branch e suffered from: `TwoPhotonAbsorptionBichromatic` existed as a
    # type, appeared in the documentation, and died in a MethodError that named neither the scheme nor the
    # missing piece. Here the error names the process, lists what would have to be built, and points at the
    # scheme to use meanwhile -- so the failure is informative rather than merely a failure.
    #
    # What three-photon needs: a THIRD-order amplitude (a double sum over two sets of intermediate states, with
    # all 3! = 6 photon orderings against 2! = 2 here), energy sharings on a two-dimensional simplex rather than
    # a line, and the coupling of three multipoles, so that the rank K is no longer fixed by oplus(J_f, J_i).
    ni          = Nuclear.Model(1.0, "point")   ## Fermi cannot represent Z = 1; see the note in the header
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.ThreePhotonEmissionScheme(),
                        multipoles = [E1], gauges = [UseCoulomb],
                        intermediateStates = Multiplet(), calcOverview = false,
                        lineSelection = LineSelection() )
    wf = Atomic.Computation(Atomic.Computation(), name="Dh-j: three-photon emission (expected to fail)", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("3s")],
                            finalConfigs   = [Configuration("1s")],
                            processSettings= mpSettings )
    perform(wf)
    #
end
#
setDefaults("print summary: close", "")
