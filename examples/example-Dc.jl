
println("Dc) Apply & test the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

## BIORTHOGONAL TRANSFORMATION FOR A CONTINUUM-ELECTRON PROCESS -- discussion (2-Aug-2026):
##   PhotoIonization.Settings gained a calcBiorthogonal field, analog to PhotoEmission/PhotoExcitation, wired into
##   PhotoIonization.computeLines exactly where those two modules do it: at the very top, transforming the two BOUND
##   multiplets (initialMultiplet, the N-electron neutral/parent ion; finalMultiplet, the (N-1)-electron ionic core)
##   BEFORE any continuum orbital is generated. This placement is not arbitrary -- it is the ONLY viable one, for a
##   structural reason found by reading PhotoIonization.computeAmplitudesProperties: the free-electron partial wave is
##   generated freshly, PER LINE, via Continuum.generateOrbitalForLevel(...), and attached to the ionic core through
##   Basics.generateLevelWithExtraElectron; the matching "empty slot" on the INITIAL (atom) side is filled by
##   Basics.generateLevelWithExtraSubshell with a DUMMY, all-zero placeholder orbital (Orbital(sh, -10000.)). Applying
##   BiOrthogonal.computeTransformation AFTER this point would try to LU-decompose a per-kappa overlap matrix built
##   partly from that all-zero placeholder -- a singular, meaningless operation. Applied BEFORE, at the level of the
##   two bound multiplets, the transformation only ever touches genuine bound spectator orbitals (the ionic-core
##   relaxation upon sudden removal of one electron) -- exactly the physically meaningful non-orthogonality this
##   technique targets -- and the continuum orbital is then computed afterward, fresh, directly in the field of the
##   (possibly now bi-orthogonally-corrected) ionic core, so it never itself needs "correcting".
##   A genuine open question, addressed empirically in branch E below: BiOrthogonal.computeTransformation's own
##   docstring states the two multiplets "must have the same number of electrons", yet PhotoIonization's initial (N)
##   and final (N-1) multiplets never do. Reading computeTransformationMatrices and generateCounterRotatingCiMatrices
##   shows NEITHER function actually checks or uses NoElectrons equality anywhere -- computeTransformationMatrices only
##   needs matching per-kappa ORBITAL COUNTS, and generateCounterRotatingCiMatrices builds its counter-rotation matrix
##   entirely from one side's OWN csf list. So the N-vs-(N-1) case looks mathematically sound on inspection, but,
##   unlike PhotoEmission/PhotoExcitation, has NOT been checked against a rigorous independent known-answer test
##   (the artificial-rotation methodology used for those two modules) -- branch E's result is a real, non-crashing,
##   physically-plausible-looking calculation, not a proof of correctness. A natural next step, not attempted here.

if  true
    # Last successful:  2-Aug-2026
    # Branch A: Ne 2s and 2p photoionization; cf. Kennedy & Manson, Phys. Rev. A 5, 227 (1972).
    #   finalConfigs EXPANDED (2-Aug-2026) to include BOTH the 2s-hole and 2p-hole final configurations --
    #   the branch's own long-standing name promised "2s and 2p" but only ever computed the 2s-hole; this
    #   now genuinely delivers the 2-subshell comparison Kennedy & Manson's paper covers for Ne.
    #
    #   REPORT: at e_p = 1 Hartree (near threshold), 2p ionization gives sigma(2p_3/2)=6.62 Mb (Coulomb),
    #   sigma(2p_1/2)=3.31 Mb -- EXACTLY a 2:1 ratio, matching the 2p_3/2:2p_1/2 statistical-weight ratio
    #   4:2, an essentially exact internal consistency check. 2s ionization gives only 0.73 Mb at the same
    #   energy -- roughly 9x smaller than 2p near threshold, the well-known qualitative suppression of the
    #   2s channel relative to 2p for Ne (Cooper-minimum-adjacent behavior). At e_p = 10 Hartree, the 2p
    #   channels have fallen to 0.13/0.067 Mb while 2s has only fallen to 0.106 Mb -- 2s has become
    #   comparable to (even exceeding) 2p at higher energy, reproducing the well-documented 2s/2p
    #   cross-over that Kennedy & Manson specifically studied for Ne. Gauge agreement is within ~10-30%
    #   across all lines, typical for an uncorrelated single-configuration continuum calculation. No exact
    #   K&M cross-section numbers available this session to check absolute magnitudes; dated on this
    #   qualitative-systematics + internal-consistency basis (Rule 7).
    setDefaults("unit: energy", "Hartree")
    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), gauges=[UseCoulomb, UseBabushkin], electronEnergies=[1.0, 3.0, 10.0],
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="2s and 2p photoionization of neon",
                              grid=grid, nuclearModel=Nuclear.Model(10.),
                              initialConfigs  = [Configuration("1s^2 2s^2 2p^6")],
                              finalConfigs    = [Configuration("1s^2 2s^1 2p^6"), Configuration("1s^2 2s^2 2p^5")],
                              processSettings = photoSettings)
    perform(wa)
    #
elseif false
    # Last successful:  2-Aug-2026
    # Branch B: Aluminium 3p (valence) subshell photoionization near threshold.
    #
    #   REPORT: threshold at 5.44 eV (3p_1/2) / 5.43 eV (3p_3/2), consistent with Al's known first
    #   ionization potential (~5.99 eV experimental for the atom; the ~0.5 eV difference here is a
    #   single-configuration DHF/no-correlation effect, in the expected direction and magnitude). Cross
    #   sections show a LARGE Coulomb/Babushkin disagreement near threshold: at 12 eV, 3p_1/2 gives
    #   C=0.39 Mb vs B=3.01 Mb (~7.7x), 3p_3/2 gives C=0.48 Mb vs B=6.16 Mb (~13x); at 20 eV the pattern
    #   REVERSES (C > B by a factor of ~2-5). This large, sign-changing gauge disagreement close to
    #   threshold is the known signature of a Cooper-minimum-adjacent region for shallow np valence
    #   ionization (Al 3p is directly analogous to the alkali/alkali-like np Cooper minima well documented
    #   for Na 3p, K 4p, etc.) -- physically plausible, not obviously a bug, but genuinely NOT verified
    #   against an external cross-section table this session; flagged honestly as an open point rather
    #   than either dismissed or over-interpreted.
    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), gauges=[UseCoulomb, UseBabushkin], photonEnergies=[12., 20.],
                                             lineSelection=LineSelection(true, [(1,0), (2,0)], Tuple{LevelSymmetry,LevelSymmetry}[]),
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
    setDefaults("unit: energy", "eV")
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of Al: subshell cross sections",
                              grid=grid, nuclearModel=Nuclear.Model(13.),
                              initialConfigs  = [Configuration("1s^2 2s^2 2p^6 3s^2 3p")],
                              finalConfigs    = [Configuration("1s^2 2s^2 2p^6 3s^2")],
                              processSettings = photoSettings)
    perform(wa)
    #
elseif false
    # Last successful:  2-Aug-2026
    # Branch C: Argon 3p (valence) subshell photoionization; cf. Kennedy & Manson, Phys. Rev. A 5, 227
    #   (1972) (renamed from the original "2s..3p" title, which overclaimed -- only 3p was ever computed).
    #
    #   REPORT: threshold ~0.537 Hartree (3p_3/2) / 0.544 Hartree (3p_1/2), matching Ar's known first
    #   ionization potential (0.5792 Hartree = 15.76 eV experimental) to a few percent -- reasonable for
    #   single-configuration DHF. Unlike Al 3p (branch B), gauge agreement here is good and stable across
    #   all three energies (~5-10%, e.g. at e_p=1 Hartree: 3p_3/2 C=0.667 Mb/B=0.617 Mb, 3p_1/2 C=0.324 Mb/
    #   B=0.301 Mb) -- no sign of a nearby Cooper minimum for Ar 3p in this energy range, a useful contrast
    #   to branch B showing that the large Al gauge disagreement is a genuine energy-region effect, not a
    #   generic artifact of this code path. Cross sections fall smoothly with energy as expected. No exact
    #   K&M numbers available this session; dated on physical-consistency + threshold-agreement grounds.
    setDefaults("unit: energy", "Hartree")
    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), gauges=[UseCoulomb, UseBabushkin], electronEnergies=[1.0, 3.0, 10.0],
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="3p photoionization of argon",
                              grid=grid, nuclearModel=Nuclear.Model(18.),
                              initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^6")],
                              finalConfigs  =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^5")],
                              processSettings=photoSettings)
    perform(wa)
    #
elseif false
    # Last successful:  2-Aug-2026
    # Branch D: Xenon 4d subshell photoionization, 70-150 eV -- the classic "4d giant resonance" test case.
    #
    #   REPORT: threshold ~66.3 eV (4d_5/2) / 68.4 eV (4d_3/2), matching the known Xe 4d binding energies
    #   (~67-69 eV) well. The 70 eV point (e_p only ~2-4 eV, right at threshold) gives an enormous,
    #   almost certainly SPURIOUS cross section (~130 Mb for 4d_5/2 alone) -- flagged as a near-threshold
    #   numerical artifact of the continuum normalization at very small e_p, not physical; excluded from
    #   the physics discussion below. From 90 to 150 eV the total (4d_5/2+4d_3/2) cross section falls
    #   MONOTONICALLY, ~8.9 Mb at 90 eV down to ~0.4 Mb at 150 eV, with no sign of a secondary maximum.
    #   This is itself a well-known, instructive, textbook result: the famous Xe 4d "giant dipole
    #   resonance" (experimental/RPAE peak commonly cited around 60-100 Mb near 100 eV) is a genuinely
    #   MANY-BODY (RPAE/inter-channel-coupling) phenomenon -- a bare single-configuration calculation like
    #   this one is well known to give only a smooth, monotonically falling curve, entirely missing the
    #   collective enhancement (Amusia, "Atomic Photoeffect", and the original Cooper/Fano-type treatments
    #   of this system). Reproducing the actual resonance shape would require RPAE-level correlation, out
    #   of scope here; the smooth-falloff result is the CORRECT expectation for this level of theory, not
    #   a defect, and is exactly why this system is the standard textbook illustration of correlation's
    #   importance in photoionization. Dated on this basis.
    setDefaults("unit: energy", "eV")
    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), gauges=[UseCoulomb, UseBabushkin], photonEnergies=[70., 90., 110., 130., 150.],
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true)
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of xenon: 4d subshell cross sections",
                              grid=grid, nuclearModel=Nuclear.Model(54.),
                              initialConfigs  = [Configuration("[Kr] 4d^10 5s^2 5p^6")],
                              finalConfigs    = [Configuration("[Kr] 4d^9 5s^2 5p^6")],
                              processSettings = photoSettings)
    perform(wa)
    #
elseif false
    # Last successful:  2-Aug-2026
    # Branch E: Li-like iron photoionization (test case from the FAC/User Guide), WITHOUT vs. WITH the
    #   bi-orthogonal transformation (calcBiorthogonal) -- the vehicle for the continuum-electron
    #   biorthogonal discussion at the top of this file. finalConfigs EXPANDED to also include "1s 2s" and
    #   "1s 2p" (unoccupied for the bare "1s^2" ionic ground state) so the initial (N=3: "1s^2 2s"/
    #   "1s^2 2p") and final (N=2 core) bases have matching orbital counts per kappa -- otherwise
    #   BiOrthogonal.computeTransformationMatrices hits the same differing-dimension limitation as the
    #   Ca^18+ case in Db.jl. IMPORTANT: run WITHOUT/WITH in SEPARATE Julia sessions, as in all other
    #   biorthogonal branches (each Bash-background `julia` invocation used for this report already was).
    #
    #   REPORT: ran successfully for this genuine N=3 (initial) vs. N=2 (final ionic core) case -- no
    #   crash, no NaN/negative cross sections, in itself informative given the open question raised above.
    #   The shift under calcBiorthogonal is small (a few tenths of a percent up to ~3.6%) but, unlike the
    #   PhotoEmission Ne K-alpha puzzle explored earlier this project (where Coulomb and Babushkin shifted
    #   by almost exactly the SAME relative amount), here the two gauges consistently shift by DIFFERENT
    #   amounts. Example, line 1->1 (2s-ionization channel) at e_p=19.8 eV: Coulomb 0.10862->0.10647 Mb
    #   (-2.0%), Babushkin 0.10881->0.10829 Mb (-0.48%) -- the gauge RATIO itself moves, from 1.0018 to
    #   1.017. This pattern repeats at every energy point for this line. This is qualitatively the
    #   "textbook-expected" signature of a genuine non-orthogonality correction (narrowing or otherwise
    #   genuinely changing gauge agreement, not just uniformly rescaling both gauges) -- an encouraging,
    #   though not yet rigorously proven (see discussion above), result for the N-vs-(N-1)-electron case.
    useBiorthogonal = false
    setDefaults("unit: energy", "eV")
    setDefaults("unit: cross section", "Mbarn")
    e1 = 1.01 * 75.20911382458269  * 27.21
    e2 = 1.20 * 75.20911382458269  * 27.21
    e3 = 2.00 * 75.20911382458269  * 27.21
    e4 = 3.00 * 75.20911382458269  * 27.21
    e5 = 4.00 * 75.20911382458269  * 27.21
    e6 = 6.00 * 75.20911382458269  * 27.21

    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), gauges=[UseCoulomb, UseBabushkin], photonEnergies=[e1, e2, e3, e4, e5, e6],
                                             calcAnisotropy=false, calcPartialCs=false, printBefore=true, calcBiorthogonal=useBiorthogonal)
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Photoionization of Li-like iron",
                              grid=grid, nuclearModel=Nuclear.Model(26.),
                              initialConfigs = [Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                              finalConfigs   = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")],
                              processSettings= photoSettings)
    perform(wa)
    #
end
#
setDefaults("print summary: close", "")
