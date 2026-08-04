
println("Dm)  Test of the InternalRecombination module with ASF from an internally generated initial and final-state multiplet.")

setDefaults("print summary: open", "zzz-InternalRecombination.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

## ABOUT THIS PROCESS (2-Aug-2026, revisited together with S. Fritzsche): "internal recombination" (also called
## "internal excitation"/"internal stabilization"/"IDE" interchangeably in the prior correspondence with
## M. Pajek & L. Michel, see apps/apps-pajek-internal-excitation/) is NOT a standard, widely-published named
## process -- a genuine web-literature search this session found no paper using this exact term; the closest
## hit was Fritzsche, Huang & Huang, Eur. Phys. J. D 79, 22 (2025) on DR into high-n Rydberg shells, which the
## user confirmed is relevant background for DR examples but not for this module specifically. The scarce,
## non-quantitative literature situation is expected and was confirmed by the user directly.
##
## Physically: an N-electron core with an ALREADY-BOUND, very-high-n Rydberg spectator electron (n=10 or
## higher) is formally degenerate with a doubly-excited (N-electron) configuration with all-low-n occupation
## -- e.g. "3s 10d" and "3p^2" for the same ion. An exact, energy-conserving transition between these two
## discrete bound levels is formally impossible unless SOME level broadening (natural, plasma/Stark, or
## similar) gives the two Lorentzian profiles finite overlap; InternalRecombination.Settings' `gamma` field is
## exactly this broadening width, and computeAmplitudesProperties' rate formula is an explicit Lorentzian,
## `rate ~ (gamma/2) / (deltaEnergy^2 + (gamma/2)^2)`, which -> 0 as gamma->0 unless deltaEnergy is exactly
## zero. This matches the user's own framing precisely: "formally a fully resonant process which can hardly
## happen in atoms/ions... becomes possible by line-broadening ... or similar."
##
## The user's own motivating physics question (this session): can this mechanism explain why some
## experimentally-produced ions with a high-n Rydberg spectator (n=10+) are nonetheless observed to decay
## FAST via photoemission? If the "3s 10d -> 3p^2"-type internal-recombination channel is strong, the ion can
## jump directly from a high-n Rydberg configuration to an all-low-n doubly-excited one, which then decays via
## an ORDINARY, fast allowed radiative transition (e.g. 3p->3s) -- potentially much faster overall than the
## Rydberg electron's own slow, many-step radiative cascade down from n=10. This is exactly the scenario
## explored (for Xe^35+) in the 2023-2024 correspondence with M. Pajek's group, motivated by their measured
## fluorescence yields for highly-charged Xe ions with high-n Rydberg spectators.
##
## NO BIORTHOGONAL TRANSFORMATION NEEDED here (confirmed again this session): capture is always into a
## valence/Rydberg-adjacent orbital, never a core hole, so no core-relaxation non-orthogonality of the kind
## PhotoEmission/PhotoIonization/AutoIonization's calcBiorthogonal addresses arises. There is also a
## structural reason specific to this module: initialLevel and finalLevel (as passed to Atomic.Computation)
## represent the SAME electron count (the Rydberg electron is attached dynamically to BOTH sides via
## rydbergOrbitals/generateLevelWithExtraElectron resp. the dummy-placeholder generateLevelWithExtraSubshell
## trick) -- unlike PhotoIonization/AutoIonization's genuine N-vs-(N-1) mismatch.
##
## *** BUG FOUND, ROOT CAUSE CONFIRMED, FIXED -- BUT A SECOND, OPEN ISSUE REMAINS (2/3-Aug-2026) ***
## Branch 1 below is a faithful reproduction of a real, previously-validated case from
## apps/apps-pajek-internal-excitation/job-internal-a.jl and report-06 (5 May 2024): Xe^35+ K-like
## [Ar]3d(10g) --> [Ar](4s^2+4s4p+...+4f^2), same grid, same settings, same configs. The companion radiative
## calculation (job branch 2, [Ar](4s+4p+4d+4f)-->[Ar]3d) reproduces report-06's numbers EXACTLY, to every
## printed digit -- confirming the SCF/grid/orbital machinery in general has not drifted; the discrepancy is
## specific to InternalRecombination's own Rydberg-orbital generation.
##
## ROOT CAUSE #1, CONFIRMED + FIXED (3-Aug-2026): InternalRecombination.computeLines built its Rydberg-orbital
## mean-field potential as `Basics.computePotential(Basics.DFSField(1.0), grid, basis)` alone -- this returns
## ONLY the electronic screening potential, NOT the nuclear attraction (confirmed directly: potential.Zr at
## r=0.053 a.u. for Xe was -17.17, vs. the correct Nuclear.nuclearPotential(...).Zr=54.01 at the same point --
## the nuclear term was simply absent). Diagonalizing the Dirac equation with no real potential well removes
## the normal energetic gap between the unphysical negative-energy "Dirac sea" branch and genuine atomic bound
## states, so Bsplines.findPositiveBranchStart returned orbitals from right at the -1.999*c^2 Dirac-sea
## threshold (~-37530 Hartree) instead of real bound Rydberg states (which should be at ~-6 to -1500 Hartree --
## confirmed via the SCF-generated 3d orbital, -164 Hartree, sane). FIXED: meanPot now explicitly adds
## Nuclear.nuclearPotential(nm, grid), matching the pattern already used correctly in
## module-BasicsAZ-inc-generate.jl. Verified directly: re-diagonalizing with the fix immediately restores a
## physically sane bound-state spectrum (e.g. g_9/2: 1 near-threshold "bound" state then straight into
## positive/unbound energies -> 7 properly bound states at -5.14/-3.57/-2.62 Hartree, a clean decreasing
## Rydberg series). Full diagnostic + a code-level warning in Bsplines.findPositiveBranchStart's docstring:
## see project_bsplines_spurious_dirac_sea_bug.md. This is a GENUINE, INDEPENDENTLY-VERIFIED bug fix, kept in
## module-InternalRecombination.jl regardless of the item below.
##
## ROOT CAUSE #2, STILL OPEN: fixing the above does NOT restore agreement with report-06. Rates now OVERSHOOT
## by ~1e7-1e9x (previously they UNDERSHOT by ~1000-3700x) -- e.g. level 1 total rate-with-shift: 2.215e19 /s
## now vs. 1.284e7 /s before the fix vs. 4.716e10 /s in report-06. Since line.deltaEnergy (the Lorentzian
## rate-formula input) comes only from the ordinary SCF initial/final level energies -- untouched by the
## potential fix, which only changes Rydberg-orbital generation -- the entire further swing must be in the
## Vee matrix-element amplitude itself. NOT investigated further this session (deliberately, per explicit
## user decision to keep the verified fix and flag the remainder as open rather than keep digging in what was
## meant to be a wrap-up session) -- a genuinely separate, unexplained issue for a dedicated follow-up.
## Candidates for that follow-up: compare individual channel amplitudes (not just totals) across pre-fix/
## post-fix/report-06 to see whether the blowup is a uniform global scale factor or concentrated in specific
## channels; also consider that report-06's own 2024 numbers were computed with the since-deleted
## Basics.generateOrbitalsForPotential, which may not have been a clean ground truth either.
##
## Until root cause #2 is resolved, treat ALL absolute InternalRecombination rate magnitudes (both branches
## below) with real caution -- the qualitative/comparative conclusions drawn in each branch's REPORT are still
## meaningful, but the absolute numbers are not yet trustworthy.
##
## Also noted in passing: InternalRecombination.displayTotalRates' "A_IDE / A_radiative" column is computed
## via `rand(1)[1]` -- literally a random number, not a real ratio. This matches report-06's own account that
## no clean, single "mean ratio" was ever actually automated in the code; the user manually combined separate
## IR-rate and radiative-rate runs by hand. Both branches below follow that same manual-comparison approach
## (see each REPORT) rather than relying on this non-functional display column.

if  true
    # Last visit:  3-Aug-2026
    # Last successful:  unknown ... (see REPORT -- root cause #2 above still open, absolute magnitude not
    #   yet trustworthy)
    # Branch 1: reproduction of the Xe^35+ K-like case from apps/apps-pajek-internal-excitation/
    #   job-internal-a.jl and report-06 (5-May-2024, correspondence with M. Pajek & L. Michel):
    #   [Ar] 3d (10g) --> [Ar] (4s^2 + 4s4p + 4s4d + 4s4f + 4p^2 + 4p4d + 4p4f + 4d^2 + 4d4f + 4f^2)
    #
    #   REPORT: run three times this session, tracking the investigation above. (1) BEFORE the nuclear-
    #   potential fix: total IR rates ~1000-3700x SMALLER than report-06 (level 1 rate-with-shift 1.284e7 /s
    #   vs. report-06's 4.716e10 /s) -- traced to Rydberg orbitals landing in the spurious Dirac-sea branch.
    #   (2) AFTER the fix (module-InternalRecombination.jl, 3-Aug-2026): the Rydberg-orbital spectrum is now
    #   independently confirmed sane (real bound states at the right energy scale, see the file-level note),
    #   but the resulting IR rates now OVERSHOOT report-06 by ~1e7-1e9x instead (level 1 rate-with-shift
    #   2.215e19 /s) -- root cause #2, still open (see above). The companion radiative calculation (branch
    #   below this comment, run separately) reproduces report-06's numbers EXACTLY in all three runs (never
    #   affected) -- confirming the SCF/grid/orbital machinery in general is fine and the issue is confined
    #   to InternalRecombination's own Rydberg-orbital/matrix-element pathway. Companion radiative rates
    #   (Xe^35+, [Ar](4s+4p+4d+4f)-->3d, matches report-06 part H exactly): odd-parity (4p/4f-derived) levels
    #   decay FAST via E1 (up to 3.6e13 /s, sub-attosecond-to-few-attosecond lifetimes), even-parity
    #   (4s/4d-derived) levels only via forbidden M1/E2 (as slow as ~1e6-1e9 /s) -- a clean selection-rule
    #   explanation for why only a few of the many doubly-excited configurations actually contribute to the
    #   observed fast X-ray emission, reproducing report-06's own qualitative conclusion. NOT dated pending
    #   resolution of root cause #2.
    rydbergShells = Basics.generateShellList(10, 10, "g")
    irSettings    = InternalRecombination.Settings(rydbergShells, true, -0.5, 0.05, LineSelection(), CoulombInteraction())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(54.01),
                            initialConfigs = [Configuration("[Ar] 3d")],
                            finalConfigs   = [Configuration("[Ar] 4s^2"), Configuration("[Ar] 4s 4p"), Configuration("[Ar] 4s 4d"),
                                              Configuration("[Ar] 4s 4f"), Configuration("[Ar] 4p^2"), Configuration("[Ar] 4p 4d"),
                                              Configuration("[Ar] 4p 4f"), Configuration("[Ar] 4d^2"), Configuration("[Ar] 4d 4f"),
                                              Configuration("[Ar] 4f^2")],
                            processSettings= irSettings )

    wb = perform(wa)
    #
elseif  false
    # Last visit:  3-Aug-2026
    # Last successful:  unknown ... (see REPORT -- absolute magnitude not yet trustworthy, see finding above)
    # Branch 2: new, minimal illustrative case for Na-like Fe (Fe XVI, Z=26), directly testing the user's own
    #   hypothesis: [Ne] 3s (10d) --> [Ne] 3p^2 -- does this internal-recombination channel compete with (or
    #   even outrun) the subsequent ordinary 3p^2 -> 3s3p E1 decay, i.e. can it explain fast photoemission
    #   from an ion carrying a high-n Rydberg spectator?
    #
    #   REPORT: total IR rate (level 1, J=1/2+, no resonance shift applied) = 1.298e12 /s (both "dE=0" and
    #   "with shift" columns coincide exactly since resonanceEnergyShift=0 here, as expected from the rate
    #   formula) -- this number predates the nuclear-potential fix above and has NOT been rerun with it;
    #   given root cause #2's rate swing is orders of magnitude, this specific number should now be treated
    #   as stale, not just "uncertain". Companion radiative calc (run separately, unaffected by the IR fix):
    #   [Ne]3p^2 -> [Ne]3s3p E1 rate = 1.457e10 (Coulomb) / 1.960e10 (Babushkin) /s for the dominant line,
    #   i.e. lifetime ~51-68 ps -- this part remains valid. The qualitative comparison drawn before the fix
    #   (IR step ~66-89x faster than the subsequent radiative step, i.e. IR not the bottleneck, supportive of
    #   the user's hypothesis) should be re-run once root cause #2 is resolved before drawing any conclusion
    #   with confidence -- NOT re-verified this session. NOT dated pending resolution of root cause #2.
    rydbergShells = Basics.generateShellList(10, 10, "d")
    irSettings    = InternalRecombination.Settings(rydbergShells, true, 0., 0.05, LineSelection(), CoulombInteraction())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=0.6e-2, rbox=15.0),
                            nuclearModel=Nuclear.Model(26.),
                            initialConfigs = [Configuration("[Ne] 3s")],
                            finalConfigs   = [Configuration("[Ne] 3p^2")],
                            processSettings= irSettings )

    wb = perform(wa)
end

setDefaults("print summary: close", "")
