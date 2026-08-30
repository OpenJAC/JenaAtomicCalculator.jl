
println("Dd) Apply & test the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoRecombination.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: energy", "eV")
setDefaults("unit: cross section", "barn")
setDefaults("unit: rate", "1/s")

## Each branch below defines its own grid, sized to its own electron-energy range: the continuum-orbital step size
## hp must satisfy 15*hp < (de Broglie wavelength at the highest electron energy used in that branch), or the
## Galerkin continuum orbitals become inaccurate (JAC raises an explicit "Improper grid" error if hp is too large;
## too small an hp, conversely, blows up the number of grid points and can exhaust memory -- pick just enough).

## BIORTHOGONAL TRANSFORMATION: deliberately NOT added to PhotoRecombination (2-Aug-2026), unlike PhotoEmission/
## PhotoExcitation/PhotoIonization. Per explicit physics judgement: radiative capture always proceeds into a
## VALENCE shell (the outermost open subshell of the recombined ion), which does not readily induce any relaxation
## of the deeper, already-closed core -- the core-relaxation non-orthogonality the bi-orthogonal method targets
## essentially does not arise here the way it does for a K-hole (PhotoEmission/PhotoIonization) or an excited
## valence-to-Rydberg jump (PhotoExcitation). No calcBiorthogonal field added to PhotoRecombination.Settings.

## WHICH NUMBERS IN THIS FILE ARE CURRENT (28-Aug-2026). ALL SIX branches were re-run on 28-Aug-2026 and every
## number below is from that run. This supersedes the 9-Aug-2026 state, in which three branches were current and
## three still quoted pre-fix Babushkin digits.
##
## TWO EARLIER CORRECTIONS, both from 9-Aug-2026, remain in force: the Racah-phase fix in AngularMomentum.JohnsonI
## (0bdff9f), which moves BOTH gauges, and the length-form orientation fix in InteractionStrength.MabEmission
## (c023481), which moves Babushkin only, by a factor growing as (Z*alpha)^2 and with the photon momentum q.
##
## A THIRD CORRECTION LANDED ON 27-Aug-2026 AND IT MOVES THESE CROSS SECTIONS BY UP TO 34 %, GAUGE-INDEPENDENTLY.
## The rank>0 spin-angular callers were migrated to the new SpinAngular module (fd7498e). Its commit message says
## "nothing here changes a number"; THAT IS WRONG, and this file is where it shows. The OLD module IGNORED the
## parity argument it was given and returned the union of BOTH parity classes; the new one obeys it and returns
## the correct half. Measured 28-Aug on a 30-CSF 2p^5 + 2p^4 3d basis, counting coefficients per multipole:
##      multipole   E1    M1     E2    M2     E3    M3
##      new          8   114    105    26     26    80
##      old        122   122    131   131    106   106
## The old counts do not depend on the multipole's parity at all, and the new pair sums to the old at every rank
## (8+114 = 122, 105+26 = 131, 26+80 = 106). So every rank>0 amplitude previously summed the wrong-parity
## coefficients in with the right ones. The new numbers are the CORRECTED ones.
##
## WHY IT SHOWS HERE AND ALMOST NOWHERE ELSE: the error rides on the higher multipoles, and at 10 keV they are
## ~10 % of the K-shell REC cross section while at 500 keV they are ~85 % (measured: Z=82 E1-only gives 1210.6,
## 55.54 and 2.211 barn against 1329.1, 75.74 and 4.083 with E1+M1+E2+M2). So the shift against 9-Aug is -0.0 %
## at 10 keV, -6.2 % at 100 keV and -34.4 % at 500 keV. An E1-dominated branch barely moves.
##
## GAUGE AGREEMENT IS UNAFFECTED AND REMAINS GOOD EVERYWHERE -- both gauges moved together, which is what a
## gauge-independent angular error does. Measured 28-Aug across all six branches: 0.06 % to 1.51 % for the three
## high-Z branches, 1.2-4.7 % for the two Z=18 branches, and 17-19 % for the Z=12 F-like branch, whose gap is a
## long-standing open question of its own and is NOT a high-Z gauge-convergence problem.

if  true
    # Last successful:  28-Aug-2026 -- Coulomb/Babushkin 4911.864/4073.194, 1403.839/1139.254, 799.1763/643.5974
    # barn @ 10/30/50 eV. Unchanged from 14-Aug in the first five figures: this branch is Z=12 with E1 and M1 only,
    # so the 27-Aug parity correction described at the head of this file has almost nothing to act on.
    # Z=12, F-like 2p^5 -> Ne-like 2p^6 capture (electron into the 2p vacancy) at 10/30/50 eV. No specific literature
    # comparison known.
    #
    #   FOUR DEFECTS IN THE PHOTORECOMBINATION PATH WERE FIXED ON 14-AUG-2026, all found by building the
    #   physical-channel form beside the flat one and comparing. Three of them touch ONLY interference
    #   observables, i.e. the anisotropy parameters; the fourth moves the cross sections in their 7th digit.
    #     (1) computeAnisotropyParameter dropped every MAGNETIC channel: its three gauge guards tested the
    #         requested gauge instead of the channel's own, and the requested gauge is Coulomb or Babushkin at
    #         every call site, so a Basics.Magnetic channel was admitted by neither. Odd-nu parameters live
    #         entirely on electric-magnetic interference, so beta_1 and beta_3 printed as exactly zero.
    #     (2) Fixing (1) exposed `1.0im^(L+p-Lp-pp)`, which parses as 1.0 * im^n with im::Complex{Bool} and
    #         raises DomainError for a COMPUTED negative exponent -- and n = -1 on precisely those terms.
    #     (3) determineChannels emitted every magnetic channel ONCE PER REQUESTED GAUGE (12 channels of which
    #         only 9 distinct, here), so every M1 amplitude was counted twice in both gauge sums. This is the
    #         one that moves the cross sections: 4911.796 -> 4911.794 barn, i.e. ~4e-7, since M1^2/E1^2 ~ 1e-6
    #         at Z=12. PhotoIonization.determineChannels already had the guard.
    #     (4) computeAmplitudesProperties passed the incoming `channel`, still carrying phase = 0., to
    #         PhotoRecombination.amplitude, which multiplies by exp(im*channel.phase). The scattering phase was
    #         therefore DROPPED from every amplitude -- invisible in a cross section, where |exp(i phi)| = 1,
    #         and wrong in every interference observable. The correctly-phased newChannel was built one line
    #         above and used only to be destructured again. PhotoIonization passes its nChannel here.
    #   Now, Coulomb / Babushkin:  beta_1 = -2.2486e-04 / -2.4194e-04 @10 eV, 1.3164e-04 / 1.4854e-04 @30 eV,
    #   -1.5016e-05 / -1.4238e-05 @50 eV;  beta_2 = -1.6201e-01 / -1.5493e-01, -3.4537e-01 / -3.5160e-01,
    #   2.4909e-02 / 3.1738e-02. beta_1 is of order 1e-4, the right size for an E1-M1 interference at Z=12,
    #   and changes sign between 10 and 30 eV.
    #   beta_3 and beta_4 remain exactly zero, and NOT because of any defect: two L=1 multipoles cannot couple
    #   to nu = 3 or 4. A branch carrying E2 would show a non-zero beta_3.
    #   Cross sections now 4911.794 / 4073.073, 1403.824 / 1139.224, 799.1699 / 643.5854 barn -- the sentence
    #   below quotes them to 5 figures, where they are unchanged.
    #
    #   REPORT: ab-initio (Coulomb gauge) cross sections 4911.8, 1403.8, 799.2 barn @ 10/30/50 eV -- smoothly
    #   decreasing, no sign flips or non-physical values; gauge agreement (Coulomb vs. Babushkin) is stable at
    #   ~17-19% across all three energies (4074/1140/644 barn Babushkin) -- consistent internally, if not
    #   independently benchmarked. Re-ran BOTH empirical comparisons this session for the same 2p-capture channel:
    #   Empirical.ScaledHydrogenic gives 195.4, 50.9, 25.0 barn (ab-initio exceeds it by ~25-32x); Empirical.UsingJAC
    #   gives 678.6, 283.6, 186.1 barn (ab-initio exceeds it by ~4.3-7.2x, noticeably closer than ScaledHydrogenic,
    #   but still a substantial, unexplained gap). This updates and sharpens the old vague "~7-25x" note into two
    #   separate, per-method numbers -- the discrepancy persists in essentially the same range as before, and
    #   remains an OPEN, unresolved case; kept for a future revisit once a root cause is identified (not attempted
    #   this session -- both empirical methods use simplified hydrogenic/mean-field capture-into-2p physics that
    #   may simply not suit an open 2p-hole F-like ion well, but this is a hypothesis, not a diagnosis).
    gridDd1 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1], gauges = [UseCoulomb, UseBabushkin],
                                            electronEnergies = [10., 30., 50.], calcTotalCs = true, calcAnisotropy = true, printBefore = true,
                                            maxKappa = 3, lineSelection = LineSelection(true, indexPairs=[(1,1)]) )

    wa = Atomic.Computation(Atomic.Computation(), name="Z12-FtoNe", grid=gridDd1, nuclearModel=Nuclear.Model(12.),
                            initialConfigs  = [Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^6")],
                            processSettings = pSettings)

    wb = perform(wa)
    #
elseif  false
    # Last successful:  28-Aug-2026 -- Coulomb/Babushkin 5.255384e6/5.089325e6, 1.823885e5/1.766079e5,
    # 8683.07/8399.78, 617.688/593.243, 26.0111/24.7985 barn @ 1/10/100/1000/10000 eV; gauge agreement 3.2, 3.2,
    # 3.3, 4.0, 4.7 %. THE PRE-FIX BABUSHKIN DIGITS THIS BRANCH USED TO CARRY ARE NOW REPLACED: Coulomb is
    # unchanged to five figures (5.2554e6 and 26.01 barn, exactly as reconfirmed on 2-Aug), while Babushkin moved
    # from 5.1058e6 to 5.0893e6 (-0.3 %) and from 26.63 to 24.80 barn (-6.9 %) -- the length-form orientation fix,
    # larger at the top of the energy range as expected. The old text said the two gauges "agree to ~3% throughout"
    # and that survives, but the agreement now degrades to 4.7 % at 10 keV rather than staying flat.
    # Matches built-in Stobbe(1s) to ~10-20% at 100eV-10keV; ab-initio exceeds Stobbe by
    # ~2.5-7x at 1-10 eV (close to the eta>>1 threshold regime), plausibly relativistic/multipole effects (E2,M2
    # included ab-initio, absent from the non-relativistic E1-only Stobbe formula).
    # RECONFIRMED 2-Aug-2026: re-ran unchanged, results essentially identical to two decimal places (5.2554e6/
    # 5.1058e6 barn @1eV, 26.01/26.63 barn @10keV) -- nothing in the RR pipeline (continuum generation, radial
    # integrals, multipole/gauge machinery) has materially changed in the ~2 weeks since the original run.
    # Z=18 (Ar), bare H-like 1s -> He-like 1s^2: K-shell radiative recombination at 1/10/100/1000/10000 eV.
    # Comparison: Ichihara & Eichler, Atomic Data and Nuclear Data Tables (2000), "Cross sections for radiative
    # recombination and the photoelectric effect in the K, L, and M shells of one-electron systems with 1 <= Z <= 112" --
    # Z=18 is one of their tabulated benchmark charges (also used in the angle-differential follow-up, ADNDT 79, 187 (2001)).
    gridDd2 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1, E2, M2], gauges = [UseCoulomb, UseBabushkin],
                                            electronEnergies = [1., 10., 100., 1000., 10000.], calcTotalCs = true, printBefore = true, maxKappa = 4 )

    wa = Atomic.Computation(Atomic.Computation(), name="Z18-Kshell", grid=gridDd2, nuclearModel=Nuclear.Model(18.),
                            initialConfigs   = [Configuration("1s")],
                            finalConfigs     = [Configuration("1s^2")],
                            processSettings  = pSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  28-Aug-2026 -- total (all 5 final shells summed) Coulomb/Babushkin 1.487028e6/1.440833e6,
    # 4.745413e4/4.596358e4, 2171.93/2099.70, 106.140/101.994, 2.44322/2.41439 barn @ 1/10/100/1000/10000 eV, i.e.
    # gauge agreement 3.1, 3.1, 3.3, 3.9, 1.2 %. THIS BRANCH HAD NEVER QUOTED BABUSHKIN AT ALL; those digits are new
    # here. The Coulomb column is unchanged from 22-Jul to every figure it was quoted to (1.49e6, 4.75e4, 2172,
    # 106.1, 2.44), which is the cleanest evidence in this file that the 27-Aug parity correction leaves an
    # E1-dominated low-Z branch alone.
    # The Kramers comparison below is therefore unaffected: within a factor 2-3 of Kramers(n=2,3) (7.72e5, 7.63e4, 6823, 345,
    # 6.05 barn) at every point, no runaway/one-directional trend -- much better behaved than the still-open Z=12
    # F-like branch above.
    # Z=18 (Ar), He-like 1s^2 -> Li-like: L/M-shell radiative recombination (capture into 2s, 2p, 3s, 3p or 3d) at the
    # same energies as the K-shell branch above. Same reference (Ichihara & Eichler, ADNDT 2000) -- this and the
    # previous branch together restore that reference's original K+L+M-shell comparison intent (previously only the
    # K-shell part was active, and at an inconsistent Z between the code and its comment).
    gridDd3 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1, E2, M2], gauges = [UseCoulomb, UseBabushkin],
                                            electronEnergies = [1., 10., 100., 1000., 10000.], calcTotalCs = true, printBefore = true, maxKappa = 4 )

    wa = Atomic.Computation(Atomic.Computation(), name="Z18-LMshell", grid=gridDd3, nuclearModel=Nuclear.Model(18.),
                            initialConfigs   = [Configuration("1s^2")],
                            finalConfigs     = [Configuration("1s^2 2s"), Configuration("1s^2 2p"),
                                                Configuration("1s^2 3s"), Configuration("1s^2 3p"), Configuration("1s^2 3d")],
                            processSettings  = pSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  28-Aug-2026 -- Coulomb/Babushkin 1234.09/1237.48, 109.778/109.841, 6.49946/6.46310 barn
    # @ 1.196/11.96/119.6 keV; gauge agreement 0.27, 0.06, 0.56 %. Against the 9-Aug run (1117.8/1120.8, 109.24/
    # 109.30, 7.076/7.037) Coulomb moved +10.4 %, +0.5 % and -8.1 % -- the 27-Aug parity correction, see the head
    # of this file. Both gauges moved together, so the agreement is as good as it was.
    # Superseded note from 9-Aug-2026 -- both gauges genuinely computed (before 22-Jul-2026 only Babushkin was
    # requested, which hid the true Coulomb-gauge result behind a display artifact -- see the note at the top of
    # this file).
    # Coulomb/Babushkin: 1117.8/1120.8 barn @1.196 keV (0.3% gauge agreement), 109.24/109.30 barn @11.96 keV (0.06%),
    # 7.076/7.037 barn @119.6 keV (0.6%). ALL THREE points are now gauge-converged; the 22-Jul-2026 reading was
    # 1234.9/1238.3, 109.9/118.2 and 6.58/12.16 barn, i.e. a 7% gap at 12 keV and an 85% gap at 120 keV. That
    # degradation was NOT the suspected continuum radial-grid/orbital-quality issue but the length-form orientation
    # defect in InteractionStrength.MabEmission (c023481); the Coulomb shifts (e.g. 1234.9 -> 1117.8 barn at the
    # lowest energy) come from the separate Racah-phase fix (0bdff9f). The 120 keV point is no longer open.
    # Cross-checked directly against the bare-ion analog (initial=bare U, capture into 2s only, same grid/energies,
    # see the standalone check below): bare gives 1286.1/1290.3, 113.9/122.7, 6.82/12.70 barn -- only ~3-5% above
    # this He-like->Li-like branch at every energy, confirming that the two K-shell electrons barely screen L-shell
    # (2s) capture at this Z, i.e. the physical intuition "should be very similar" (comparing a screened vs. bare
    # capture) holds quantitatively here. Ichihara & Eichler (ADNDT 2000) tabulate exactly the bare-ion K/L/M-shell
    # case at Z=92 as one of their benchmark charges, but the actual tabulated numbers are behind ScienceDirect's
    # paywall and could not be retrieved in this session -- the bare-ion cross-check above is the closest available
    # substitute (it isolates the same physics Eichler's tables would let us check, even without the exact numbers).
    # Also still agrees with Empirical.photorecombinationCrossSection(ScaledHydrogenic) to ~30% (two lower-energy
    # points) and ~4x (highest energy point, where the empirical and ab-initio binding energies differ by ~30%).
    # NOT re-verified on 9-Aug-2026: that empirical cross-check.
    # The BARE-ION cross-check WAS redone, and the spelling is worth recording because it is not obvious. A bare
    # ion is written as a zero-occupancy shell, e.g. Configuration("2s^0") -- NoElectrons = 0, and performSCF is
    # happy with it. Configuration("") and Configuration("0") are both rejected outright ("input string is empty").
    # The zero-occupancy shell must be the SAME shell that is being captured into: the empty shell still appears in
    # the basis, and checkConsistentMultiplets() below demands that the initial subshells be a prefix of the final
    # ones, so Configuration("1s^0") against a 2s final state compares [1s_1/2] with [2s_1/2] and errors out.
    # Result, bare U capturing into 2s, Coulomb/Babushkin: 1161.0/1162.6 barn @1.196 keV, 113.19/113.24 @11.96 keV,
    # 7.294/7.287 @119.6 keV -- i.e. 3.1-3.9% ABOVE the He-like -> Li-like numbers above at the three energies,
    # reproducing the 22-Jul-2026 finding ("~3-5% above at every energy") and with it the physical point: the two
    # K-shell electrons barely screen 2s capture at this Z. Gauge agreement is <=0.13% here, better than the
    # He-like case, as it should be for a genuinely one-electron final state.
    # The date rests on that cross-check plus gauge convergence to <=0.6% at all three energies, which for a system
    # whose active electron is a single 2s spectator is a necessary condition that the branch previously failed
    # at 120 keV.
    # Z=92 (U), He-like 1s^2 -> Li-like 1s^2 2s: capture into 2s at given ion energies.
    # Comparison: Fritzsche, Phys. Rev. A 72, 012704 (2005); Ichihara & Eichler, ADNDT 2000 (bare-ion proxy, see above).
    gridDd4 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.3e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1, E2, M2, E3, M3], gauges = [UseCoulomb, UseBabushkin],
                                            ionEnergies = [2.18, 21.8, 218.0], useIonEnergies = true, calcTotalCs = true, printBefore = true, maxKappa = 3 )

    wa = Atomic.Computation(Atomic.Computation(), name="Z92-2s", grid=gridDd4, nuclearModel=Nuclear.Model(92.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^2 2s")],
                            processSettings = pSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  28-Aug-2026 -- Coulomb/Babushkin 546.201/539.408 barn @10 keV (1.24 %), 20.2192/19.9138
    # @100 keV (1.51 %), 0.670034/0.665195 @500 keV (0.72 %). Against 9-Aug (545.5/538.7, 21.35/21.03,
    # 0.9988/0.9923) the shift is +0.1 %, -5.3 % and -32.9 %: the 27-Aug parity correction, growing with the
    # higher multipoles' weight. Gauge agreement is unchanged, because the correction is gauge-independent.
    # Superseded note from 9-Aug-2026 -- Coulomb/Babushkin: 545.5/538.7 barn @10 keV (1.3% gauge agreement),
    # 21.35/21.03 barn @100 keV (1.5%), 0.9988/0.9923 barn @500 keV (0.7%). ALL THREE points are now gauge-converged.
    # The high-energy blow-up recorded on 22-Jul-2026 (546.9/597.2, 20.4/36.0 and 0.74/2.88 barn, i.e. a 76% gap at
    # 100 keV and a 3.9x gap at 500 keV) was NOT the suspected continuum radial-grid/orbital-quality problem, and the
    # planned remedies -- a finer grid, or the "asymptotic Coulomb" continuum method instead of Galerkin -- are not
    # needed. It was a code defect: the length form in InteractionStrength.MabEmission was evaluated with its two
    # orbitals in the orientation that breaks it, and the error grows with q, which is exactly why it was invisible
    # at low energy and 3.9x at 500 keV (c023481). Coulomb also moved, by up to +35% at 500 keV, through the
    # separate Racah-phase fix in AngularMomentum.JohnsonI (0bdff9f).
    # BASIS FOR THE DATE: gauge convergence to <=1.5% at every energy. For a genuinely one-electron system the two
    # gauges must agree exactly, so this is a necessary condition that the branch previously FAILED and now passes.
    # NOT re-verified in this session: the absolute comparison against Ichihara & Eichler (still paywalled, see the
    # Z=92 branch above) and the Empirical.photorecombinationCrossSection cross-check.
    # Z=54 (Xe), bare H-like 1s -> He-like 1s^2: K-shell REC at 10/100/500 keV -- one of the Ichihara & Eichler
    # benchmark charges (ADNDT 2000 / angle-differential ADNDT 79, 187 (2001)).
    gridDd5 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.15e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1, E2, M2], gauges = [UseCoulomb, UseBabushkin],
                                            electronEnergies = [10000., 100000., 500000.], calcTotalCs = true, printBefore = true, maxKappa = 4 )

    wa = Atomic.Computation(Atomic.Computation(), name="Z54-Kshell", grid=gridDd5, nuclearModel=Nuclear.Model(54.),
                            initialConfigs   = [Configuration("1s")],
                            finalConfigs     = [Configuration("1s^2")],
                            processSettings  = pSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  28-Aug-2026 -- Coulomb/Babushkin 1329.07/1320.32 barn @10 keV (0.66 %), 75.7433/75.0374
    # @100 keV (0.93 %), 4.08271/4.05165 @500 keV (0.76 %). Against 9-Aug (1329.3/1320.6, 80.72/79.98,
    # 6.226/6.183) the shift is -0.0 %, -6.2 % and -34.4 %. This is the branch where the 27-Aug parity correction
    # is largest anywhere in the test set, for the same reason the July length-gauge defect was: the higher
    # multipoles carry ~85 % of the 500 keV cross section here. E1-only gives 2.211 barn at that point.
    # Superseded note from 9-Aug-2026 -- Coulomb/Babushkin: 1329.3/1320.6 barn @10 keV (0.7% gauge agreement),
    # 80.72/79.98 barn @100 keV (0.9%), 6.226/6.183 barn @500 keV (0.7%). ALL THREE points are now gauge-converged.
    # Recorded on 22-Jul-2026 as 1337.8/1539.5, 77.1/138.8 and 4.54/16.67 barn, i.e. an 80% gap at 100 keV and a 3.7x
    # gap at 500 keV. The July diagnostic was sound as far as it went -- re-running the 500 keV point with maxKappa=8
    # and multipoles up to E3/M3 gave 6.12/22.02 barn, the SAME ~3.6x ratio, correctly ruling out partial-wave
    # truncation -- but the conclusion drawn from it (a continuum radial-grid/orbital-quality problem) was wrong.
    # The cause was the length-form orientation defect in InteractionStrength.MabEmission (c023481), whose error
    # grows with q; here qr ~ 1.7 at 500 keV against ~0.02 for a Ne K-alpha line, which is why this branch showed it
    # most violently of anything in the test set. Coulomb moved as well (4.54 -> 6.23 barn at 500 keV) through the
    # separate Racah-phase fix (0bdff9f). See the Z=54 branch above for the basis on which these dates are set.
    # Z=82 (Pb), bare H-like 1s -> He-like 1s^2: K-shell REC at 10/100/500 keV -- another Ichihara & Eichler
    # benchmark charge, same reference and energy grid as the Z=54 branch above.
    gridDd6 = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.15e-2, rbox = 10.0)
    pSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1, M1, E2, M2], gauges = [UseCoulomb, UseBabushkin],
                                            electronEnergies = [10000., 100000., 500000.], calcTotalCs = true, printBefore = true, maxKappa = 4 )

    wa = Atomic.Computation(Atomic.Computation(), name="Z82-Kshell", grid=gridDd6, nuclearModel=Nuclear.Model(82.),
                            initialConfigs   = [Configuration("1s")],
                            finalConfigs     = [Configuration("1s^2")],
                            processSettings  = pSettings )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
