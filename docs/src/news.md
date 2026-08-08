

# Recently implemented features and up-dates



# 2026

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
