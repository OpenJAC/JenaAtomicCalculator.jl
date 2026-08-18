

# Recently implemented features and up-dates



# 2026

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
