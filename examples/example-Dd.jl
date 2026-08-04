
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

if  true
    # Last successful:  2-Aug-2026
    # Z=12, F-like 2p^5 -> Ne-like 2p^6 capture (electron into the 2p vacancy) at 10/30/50 eV. No specific literature
    # comparison known.
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
    # Last successful:  22-Jul-2026 -- Coulomb/Babushkin agree to ~3% throughout (5.26e6/5.11e6 @1eV down to
    # 26.0/26.6 barn @10keV). Matches built-in Stobbe(1s) to ~10-20% at 100eV-10keV; ab-initio exceeds Stobbe by
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
    # Last successful:  22-Jul-2026 -- total (all 5 final shells summed) Coulomb cross sections 1.49e6, 4.75e4, 2172,
    # 106.1, 2.44 barn @ 1/10/100/1000/10000 eV; within a factor 2-3 of Kramers(n=2,3) (7.72e5, 7.63e4, 6823, 345,
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
    # Last successful:  22-Jul-2026 -- both gauges now genuinely computed (previously only Babushkin was requested,
    # which hid the true Coulomb-gauge result behind a display artifact -- see the note at the top of this file).
    # Coulomb/Babushkin: 1234.9/1238.3 barn @1.196 keV (~0.3% gauge agreement), 109.9/118.2 barn @11.96 keV (~7%),
    # 6.58/12.16 barn @119.6 keV (~85% -- gauge convergence clearly degrades at the highest energy). Same pattern
    # seen in the Z=54/Z=82 K-shell branches below, and CONFIRMED there not to be a maxKappa truncation issue
    # (increasing maxKappa 4->8 plus higher multipoles left the gap essentially unchanged) -- most likely a
    # continuum radial-grid/orbital-quality issue at these relativistic energies rather than an angular-momentum
    # problem; treat the 120 keV point here as open/unconverged, not a validated result, until revisited.
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
    # Last successful:  22-Jul-2026 -- Coulomb/Babushkin: 546.9/597.2 barn @10 keV (~9% gauge agreement, fine),
    # 20.4/36.0 barn @100 keV (~76% gap), 0.74/2.88 barn @500 keV (~3.9x gap). The 10 keV point is trustworthy; the
    # 100/500 keV points are NOT well gauge-converged -- confirmed NOT a maxKappa truncation issue (Z=82 diagnostic
    # below, same symptom, maxKappa 4->8 plus E3/M3 left the gap essentially unchanged), most likely a continuum
    # radial-grid/orbital-quality issue at these relativistic energies (500 keV is comparable to the electron rest
    # mass) rather than an angular-momentum or gauge-formalism problem. Treat the two higher-energy points here as
    # open/unconverged, not as a validated comparison, until revisited (e.g. finer grid or the "asymptotic Coulomb"
    # continuum method instead of Galerkin).
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
    # Last successful:  22-Jul-2026 -- Coulomb/Babushkin: 1337.8/1539.5 barn @10 keV (~15% gauge agreement, fine),
    # 77.1/138.8 barn @100 keV (~80% gap), 4.54/16.67 barn @500 keV (~3.7x gap). Same pattern as the Z=54 branch
    # above (gap grows sharply with energy, not with Z alone) -- diagnosed here directly: re-running the 500 keV
    # point alone with maxKappa=8 and multipoles up to E3/M3 (vs. maxKappa=4, up to E2/M2 in the main branch) gave
    # 6.12/22.02 barn, essentially the SAME ~3.6x gauge ratio -- so this is not a partial-wave truncation problem.
    # Treat the 100/500 keV points as open/unconverged (see the Z=54 branch comment above for the likely cause);
    # only the 10 keV point is a trustworthy comparison for now.
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
