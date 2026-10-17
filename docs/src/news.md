

# Recently implemented features and up-dates



# 2026

## Version 0.7.0

!!! warning "Four breaking changes"
    * **String keys into a `perform` result now raise.** A result is asked for by the settings type that
      produced it, `res[PhotoEmission.Settings]`, or by one of the `ResultKeys` types for results that have no
      settings. The old strings were deprecated in 0.6.0 and have been removed; a retired key raises an error
      that names what to write instead.
    * **`TestFrames` is no longer exported.** From 1.0.0 an exported name is a promise, and this one is the
      test suite. Call it as `JenaAtomicCalculator.TestFrames` if you need it.
    * **The `SuperConfiguration` configuration theme has been retired**, having never had a caller. Super
      configurations are now reached through `Plasma.SuperConfigurationOpacityScheme` and the
      `SuperConfigurations` module.
    * **An unrecognized `stepping` no longer passes silently** but names the options it accepts.

* **You can now compute static polarizabilities without choosing a set of intermediate states.** The new
    [`DalgarnoLewis`](api-atomic.md) module solves for the perturbed orbital directly, which saturates the
    channel: there is no perturber sum to truncate and no box to match. Hydrogen's dipole polarizability comes
    out **4.499752 a.u. against the exact 4.5**, from two solves, and unchanged to nine digits over boxes from
    40 to 120 a.u. *(October'26)*

* **And you can now SCREEN them.** The new [`RandomPhase`](api-atomic.md) module couples the response of every
    electron to the induced field of all the others. The screened core polarizability of Ba⁺ comes out
    **9.944 a.u. against the published RPA value 10.6(5)** -- inside that value's own uncertainty, and the
    first screened polarizability in JAC that can be quoted against the literature. Exchange is not yet
    included, which is consistent with the residual. *(October'26)*

* **A caveat that matters more than a feature: do NOT compute a polarizability in the default field.**
    `AsfSettings()` defaults to `DFSField`, and measured on Li 2s against 164.07 a.u. the three fields give
    **+12 % (DFS), +39 % (Hartree-Slater) and +2.2 % (average-level)** -- same code, same basis. The response
    machinery is exact to six parts in 10⁵; the accuracy you get is set almost entirely by the mean field.
    Use `scField = Basics.ALField()` or an optimized level. *(October'26)*

* **Level energies can be fitted to measurement.** The new [`AdjustedLevels`](api-atomic.md) module scales the
    Coulomb interaction to reproduce a set of known levels and merges the result back into a multiplet, with a
    second, orthogonal (Uylings/Raassen-style) route that fits parameters independently. Its companion
    [`NistCalibration`](api-infrastructure.md) replaces computed level energies by measured ones explicitly,
    so a cascade or rate built on them says which energies it used. *(October'26)*

* **Plasma opacities from super configurations.** `Plasma.SuperConfigurationOpacityScheme` groups the
    enormous level spaces of a hot plasma into super configurations and builds the unresolved transition arrays
    between them -- a Gaussian mean and variance per array rather than a line list -- together with
    [average-atom optical properties](api-plasma.md): the Planck and Rosseland means, and a DC conductivity
    taken from the scattering phase shifts with nothing fitted. *(September/October'26)*

* **Nuclear magnetic shielding, self-consistently.** The new [`NuclearShielding`](api-properties.md) module
    computes the shielding of the nuclear magnetic moment by the electron cloud, with the induced electron
    response solved rather than estimated. *(September'26)*

* **Atomic data can be exported for the data community.** `ExportAdas` writes ADAS `adf04`, `ExportChianti`
    writes CHIANTI-style files, and `ExportVamdc` writes a VAMDC/XSAMS sketch. These are demonstrations rather
    than complete writers, and each says what it does and does not cover. *(September'26)*

* **An uncertainty you can actually quote.** `UncertaintyEstimate` runs a quantity across a few named variants
    -- a different mean field, a wider box, one more correlation layer -- and reports the spread, which is a
    defensible statement about a computed number where a single value is not. *(October'26)*

* **Isotope shifts gained a King-plot family**, so field and mass shifts can be separated from measured
    isotope shifts the way the experimental literature does it. *(September'26)*

* **A second-order treatment of a RAS Q space.** A correlation layer can now promote the strongly coupled part
    of its Q space into the CI and fold the rest in at second order, reporting its own sum of weights and
    refusing when the fold-in is not justified. It has a domain: a deep, weakly coupled core folds well, a
    same-`n` core does not, and the step says which it saw. *(September'26)*

* **Three-photon transitions are NOT in this release, and the path is still in flight.** What exists is the
    energy-sharing machinery only -- `computeLines(::ThreePhotonEmissionScheme)` says so when called -- while
    the amplitudes are not implemented. Work on that path was in progress when 0.7.0 was snapshotted and is
    deliberately not part of it. Two-photon emission and absorption are unaffected and complete. *(October'26)*

* **A caveat for two-photon rates at high nuclear charge.** The sum-over-states route of
    `MultiPhotonTransition` reproduces the exact hydrogenic Z⁶ law at Z = 1 and **fails at Z = 18** -- the two
    gauges diverge as states are added rather than converging. Use `DalgarnoLewis.twoPhotonAmplitude` instead,
    which gives 2.7678e8 /s at Z = 18 against the exact 2.7989e8. The limitation is documented at the function
    itself. *(October'26)*

## Version 0.6.0

!!! warning "One breaking change"
    `AsfSettings.maxIterationsScf` has been **removed** and replaced by `AsfSettings.scfRoute`. The iteration
    budget now belongs to the route that actually runs, because one number could not serve four solvers whose
    natural counts differ by two orders of magnitude -- which is why its default of 24 was simultaneously right
    for a mean field and far too small for an optimized level. Passing the old keyword raises an error that
    names the replacement: write, for instance, `scfRoute = Basics.RotationRoute(400)`.

* **You can now choose HOW the orbitals are found, not only what they are optimized for.** `AsfSettings` gained
    `scfRoute`, a family of routes that made explicit what was previously fixed: `MeanFieldRoute`,
    `AverageLevelRoute`, `RotationRoute` (direct minimization of the level energy) and `FockRoute` (the Fock
    equations solved to self-consistency). The default, `Basics.AutomaticRoute()`, takes the standard procedure
    for whichever field you name, so nothing changes unless you ask it to. **The new page
    [Self-consistent fields](scf-routes.md) says which to use and what each costs** -- including the measured
    accuracy of the fast route against the accurate one, from 6e-4 of the total energy near neutrality to 6e-7
    for a highly charged ion. *(September'26)*

* **The optimized-level solvers now report on themselves.** A run says whether it converged and at which
    iteration, or that it exhausted its budget while still moving; the Fock route additionally states that it is
    approximate, by how much, and checks a signature of its own failure -- the two spin-orbit partners of a
    subshell must agree in mean radius to within a bound that grows with `Z`. It also returns the best orbital
    set it saw rather than the last. Previously such a run ended in silence. *(September'26)*

* **Four process modules now print how far their two gauges agree.** `Einstein`, `PhotoExcitation`,
    `PhotoIonization` and `PhotoRecombination` add a short table giving the relative Babushkin/Coulomb deviation
    and a verdict, beside the results it qualifies, in both the log and the summary file. It costs nothing --
    the quantity was already held in both gauges -- and it is honest about its limits: a large deviation is good
    evidence that a number is wrong, a small one is NOT evidence that it is right, since both gauges can miss
    the same correlation. *(September'26)*

* **A frequency-dependent Breit request on an EOL or RAS basis now refuses**, instead of silently returning the
    static limit. Orbitals from those bases carry no orbital energy, so the frequency came out zero and
    `CoulombBreit(1.)` quietly gave the `CoulombBreit(0.)` answer. *(September'26)*

* **The Breit interaction is substantially faster.** `CoulombBreit(factor, :swept)` computes the same radial
    integrals in one pass rather than a double sum, exploiting the fact that every Breit kernel factorizes;
    measured at 116x on the kernel and 4.4x on the surrounding work. `:direct` remains the default and the
    reference form. *(September'26)*

* **Oscillator strengths and Einstein A coefficients: four defects, and why they hid each other.** The printed
    oscillator strengths and, in `Einstein.displayRates`, the A and B coefficients and decay widths are
    corrected in this release, so values from these columns change. Each defect was
    settled against hydrogen Lyman-alpha, where `A = 6.2649e8` 1/s and `f = 0.27750` (2p_3/2) are known exactly,
    rather than by comparing JAC against itself. `Basics.recast(RecastRateToOscillatorF, ...)` returned a
    quantity that is an oscillator strength in no convention -- 1.5e-07 against the tabulated 0.27750, low by
    six orders; `RecastRateToOscillatorGf` returned the UNWEIGHTED f while its name and both its callers said
    gf; `Einstein.displayRates` carried a spurious (2J_f+1), and since A, B, gf and the decay width are all
    built from that one expression, all four of its columns were affected. `Einstein`, `PhotoEmission` and
    `TwoElectronOnePhoton` had used one header for two different quantities; all three now print `g_l f` under a
    header that says `gf`, that being the tabulated convention and the one obeying the Thomas-Reiche-Kuhn sum
    rule directly.

    **Why it survived so long is worth knowing**, because it is the shape of check most people would run:
    Einstein's f and PhotoEmission's gf used to agree numerically for a shared line -- the missing `g_l` and the
    extra (2J_f+1) cancelled. Two wrong quantities produced one agreeing number, so cross-checking the two
    modules confirmed the error rather than exposing it, and repairing either defect alone looked like a
    regression. *(September'26)*

* **New modules**: `GeneralizedOscillatorStrength` (Bethe surface, bound-bound), `ResonantImpactIonization`,
    `PhotonScattering`, `PhotoRecombinationInterference` (coherent RR + DR), `WeakInteractionEnhancement` and
    `WeakInteractionMoment`. *(August/September'26)*

* **The documentation has a front door.** The home page now opens with a complete, runnable example -- the
    K-alpha lines of He-like iron -- together with the output it actually produces; and a new page,
    [How do I ...?](how-to.md), answers sixteen common questions with one worked example each, rather than
    leaving a newcomer to navigate several hundred example branches. *(September'26)*

# Earlier

!!! note "The EOL path is under active development"
    The optimized-level field is being worked on continuously, and 0.5.0 is a snapshot of it rather than a
    finished state. It converges and reports honestly whether it did -- that reporting is itself part of the
    current work -- but the optimizer, its preconditioning and its exit criteria are still changing between
    releases. **Use `Basics.EOLField()` with care for production numbers, and prefer `Basics.ALField()` where an
    individual orbital energy matters**; check the printed convergence line rather than assuming it. The
    average-level field is settled by comparison.

* **Elastic scattering of electrons and positrons, on Dirac partial waves:** The `ParticleScattering` module has
    been rebuilt around a relativistic partial-wave analysis that solves for BOTH spin-orbit partners, `kappa = -l-1`
    and `kappa = +l`. It now yields the direct and spin-flip amplitudes, the differential cross section, the Sherman
    function and the elastic, momentum-transfer and viscosity cross sections, for electrons and positrons and for
    several interaction models -- static field, Slater exchange and the energy-dependent Furness-McCarthy exchange
    used by ELSEPA. The amplitudes reproduce Eqs. (1)-(3) and (13) of the NIST Electron Elastic-Scattering
    Cross-Section Database expression for expression; the transport cross section of helium agrees with the published
    Dirac-Hartree-Fock values to 1.6 %, and the backward cross section joins the exact Rutherford limit to 0.3 % at
    2 keV. See `? ParticleScattering` and `examples/example-Ob.jl` .. `Oe.jl`. *(August'26)*

* **The radial box chooses itself:** `Basics.recommendedGrid` sizes the radial grid from the configurations under
    study rather than from a fixed default, and `Bsplines.checkOrbitalBox` judges the result from the CONVERGED
    orbitals instead of a hydrogenic stand-in. A box far too large starves the B-spline basis exactly as badly as one
    too small, and had previously been mis-diagnosed as an angular-momentum defect. *(August'26)*

* **The optimized-level (EOL) field is solved by orbital rotation:** The former solver optimized a degenerate
    stationary point. The rotation now searches conjugate directions and, since `fcbb9eb`, L-BFGS; it reports whether
    it CONVERGED or merely stopped moving, and two real defects were found on the way -- an orbital gradient that was
    not the gradient of the energy being minimized, and a line search comparing two different functions. The
    average-level field gained Anderson acceleration. *(August'26)*

* **Resonant electron-impact excitation, end to end:** `Cascade.ElectronExcitationScheme` now accounts for the
    resonant contribution -- dielectronic capture followed by radiative or Auger decay -- alongside the direct one,
    with the capture step delegated to `Cascade.DielectronicCaptureScheme`. *(August'26)*

* **Bound-free thresholds of neutral atoms:** The X-ray compilations behind `Empirical.scaledBindingEnergy` are
    core-level tables and carry no valence shell, so the threshold of every atom whose outermost shell is an s shell
    fell back on a Slater-screened estimate -- which errs in BOTH directions: Na 3s +42 %, Ca 4s +13 %, Sr 5s -22 %.
    For the ground-state valence shell of a neutral atom the binding energy IS the first ionization potential, which
    is now used and makes all three exact. *(August'26)*

* **Two-electron Auger amplitudes:** Five defects were repaired in `AutoIonization.computeTeaAmplitude`, of which the
    decisive one confined every possible intermediate space to the initial state's own orbitals -- so the feature
    threw the moment anyone tried to converge it. The rate is now computable and is being validated against the
    argon L23^2-M^3 measurement; it is not yet quantitative. *(August'26)*

* **Hyperfine multiplets are reachable again:** `Hfs` computed the hyperfine structure correctly but an `error(...)`
    stood in front of the working code, and four call sites still read field names of a type retired long ago.
    *(August'26)*

* **What the documentation publishes is now a rule, and it is checked:** The API pages selected content with
    hand-maintained allowlists and the build downgraded every consistency check to a warning, so whole modules were
    absent and nothing said so. Rule 16 of `CLAUDE.md` states what is published -- a physics module qualifies when an
    example file that uses it carries at least two verified `Last successful` dates -- and `docs/checkCoverage.jl`
    re-derives the sets rather than reading a list. Modules in the published API: 35 -> 51. *(August'26)*

* **Crystal-field splittings and crystal-field-resolved emission:** Two new modules `CrystalField` and
    `CrystalFieldEmission` support the (point-charge) splitting of atomic levels in a crystalline environment
    as well as the transitions between the split sublevels; see `? CrystalField`. *(August'26)*

* **Hyperfine-resolved dielectronic recombination:** Resonance strengths can now be resolved into hyperfine
    levels, by re-coupling the fine-structure amplitudes; this route required (and came with) a repair of the
    hyperfine representation itself. See `? DielectronicRecombination`. *(August'26)*

* **Restructured dielectronic recombination:** The DR computations were placed on explicit capture and photon
    lines, which also makes the DR satellite spectrum directly available. Note that this work uncovered and
    corrected an erroneous prefactor in the radiative rate; DR strengths obtained with earlier versions should
    be re-computed, especially in the Auger-dominated regime. *(August'26)*

* **Collisional-radiative level populations:** A new `Plasma.CollisionalRadiativeScheme` solves the kinetic
    balance of level populations of a given ion in a plasma of given temperature and density;
    see `? Plasma`. *(August'26)*

* **Corrected boundary conditions of the B-spline basis:** A defect in the (kappa-sign dependent) boundary
    conditions of the B-spline orbitals has been identified and repaired. This affects all properties that are
    sensitive to the orbitals near the nucleus; for the hyperfine, isotope-shift and Zeeman parameters, the
    agreement with published values improved considerably, and results from earlier versions should be
    re-computed. *(July'26)*

* **Bi-orthogonal transformation:** Transition amplitudes can now be computed from *independently optimized*
    initial- and final-state multiplets, by first transforming the two orbital sets to bi-orthogonality;
    this is supported for photoexcitation, photoionization and autoionization. See `? BiOrthogonal`. *(July'26)*

* **General jj-LS recoupling:** The transformation of atomic levels into a LS-coupled basis has been
    generalized to an arbitrary number of open shells; see `? LSjj`. *(July'26)*

* **AlphaVariation (q-factors):** A working implementation now computes the sensitivity of atomic levels to a
    hypothetical variation of the fine-structure constant, by re-running the SCF and CI at a slightly shifted
    alpha and taking the finite difference. See `? AlphaVariation`. *(July'26)*

* **Plasma computations made real:** The `SahaBoltzmannScheme`, `LineShiftScheme` and `AverageAtomScheme` of
    `Plasma.Computation()` were fixed and completed, replacing several placeholder or broken computation chains
    by working, typed-dispatch implementations. *(July'26)*

* **More ForPedestrians functions:** `computeChargeStateDistribution`, `computeBranchingFractions` and
    `computeLevelEnergies` (now also for isoelectronic sequences) were added to the simple-man's interface;
    see `? ForPedestrians`. *(June'26)*


# 2025

* **Simplified language to deal with electron configurations:** A number of functions have been worked out to simplify
    the generation, manipulation, extraction and display of (large lists of) electron configuration; these 
    functions also help to extract different information about leading configuration, coupling schemes, etc. 
    See, for example, `? AbstractConfigurationThemes` or `? extractFromConfiguration`. *(October'25)* 

* **ForPedestrians:** A new module `ForPedestrians` provides *simple-man's functions* to compute low-lying level
    energies of atoms and ions as well as transition rates, cross sections and selected estimates with minimum
    input but with a number of simplifying assumptions; `? ForPedestrians`. *(September'25)* 

* **Documentation:** An efficient scheme has been established to provide both, stable and development, version for the
    JenaAtomicCalculator package. *(July'25)* 

* **New B-Spline module:** The B-Spline bases of atomic orbitals has been re-implemented to make self-consistent-field
    computations more efficient. *(March'25)* 

* **Re-organized Plasma.Computation():** The computation of plasma properties have been expanded and re-organized 
    in order to support the `SahaBoltzmannScheme` and the `LineShiftScheme`. *(January'25)* 


# 2024

* **Dielectronic recombination into high-h shells:** New empirical and run-time features now support the efficient 
    computation of DR resonances strengths if the electron capture occurs into high-n (n > 15) shells.  *(October'24)* 

* **New and re-organized basic data types:** Several new abstract and concrete data types have been implemented
    (and re-organized) in the module `Basics`. *(August'24)* 

* **First design of empirical computations:** A new kind Empirical.Computation() has been established to support
    electron-impact ionization (EII) cross sections. *(March'24)* 
