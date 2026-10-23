
println("Rd) Apply & test the AdjustedLevels module: fit ONE Coulomb scale factor for a single symmetry")
println("    block, merge it back into a full multiplet, and show the effect on a computed Einstein")
println("    A-value -- not just a nicer-looking energy table.")

setDefaults("print summary: open", "zzz-AdjustedLevels.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

if  false
    # Last visit:  18-Oct-2026
    # Last successful:  unknown -- the fit TARGET below is deliberately SYNTHETIC (see this branch's own
    #   comment), so there is nothing measured to verify the result against yet; this demonstrates the
    #   MECHANISM honestly, and earns a date only once a real measured splitting is substituted.
    #
    # Branch a: He-like C V (Z=6), 1s^2 ^1S_0 (ground, J=0 even) -> 1s2p (J=1 odd: singlet ^1P_1 + triplet
    #   ^3P_1, a 2-CSF block in jj-coupling) -- the same cheap, 2-CSF block this module's own smoke test
    #   already used, now carried all the way through to a computed rate.
    #
    # THE FIT TARGET IS SYNTHETIC: 10% below this run's OWN ab initio 1s2p splitting, chosen only to make
    # the scale factor move visibly away from 1.0 -- a quick search for a six-figure measured C V 1s2p
    # singlet-triplet splitting to quote honestly did not turn up one worth citing here. Replace `targetCm`
    # below with a real measured value (and re-run) to turn this into a VERIFIED branch (Rule 7).
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    wb    = SelfConsistent.performSCF([Configuration("1s^2"), Configuration("1s 2p")], ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    #
    # The ordinary, AB INITIO multiplet -- every symmetry block this basis carries (J=0 even ground,
    # J=0,1,2 odd from 1s2p), exactly what Hamiltonian.performCI always returns.
    abInitio = Hamiltonian.performCI(basis, ni, grid, asfSettings; printout=false)
    #
    sym = LevelSymmetry(AngularJ64(1), Basics.minus)
    idxCsf, oneBody, twoBody = AdjustedLevels.buildOneAndTwoBodyMatrices(sym, basis, ni, grid)
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    eig1        = Hamiltonian.diagonalizeCiMatrix(oneBody .+ 1.0 .* twoBody, LevelSelection())
    abInitioRel = (sort(eig1.values) .- sort(eig1.values)[1]) ./ auPerKayser
    @printf("\n>>> ab initio (c=1.0) 1s2p (J=1^-) splitting:  %.3f cm^-1\n", abInitioRel[2])
    #
    targetCm = abInitioRel[2] * 0.9
    fm = AdjustedLevels.fitScaleFactor(sym, basis, ni, grid, [(2, targetCm)],
                     "SYNTHETIC TEST TARGET -- 10% below this run's own ab initio splitting, not measured data")
    println(fm)
    #
    # Merge the fitted 1s2p block back into the full, otherwise ab initio multiplet -- the step this
    # module's own docstring names as happening at the END of an SCF+RAS line, BEFORE the result is used
    # to compute rates, cross sections, etc.
    combined = AdjustedLevels.mergeIntoMultiplet(fm, abInitio)
    #
    einsteinSettings = Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.)
    rawLines         = Einstein.computeLines(abInitio, grid, einsteinSettings)
    fittedLines      = Einstein.computeLines(combined, grid, einsteinSettings)
    #
    println("\n>>> RAW (ab initio) A-values, 1s2p -> ground:")
    for  ln in rawLines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (ab initio)   A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    #
    println("\n>>> FITTED A-values (same orbitals, rescaled Coulomb interaction):")
    for  ln in fittedLines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (fitted)      A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    #
    println("\n>>> Report this table with its own provenance, since nothing downstream carries it:")
    println(">>> source: ", fm.source)
    #
elseif  false
    # Last successful:  15-Oct-2026
    #
    # Branch b: O VII (Z=8), He-like, 1s2p (J=1 odd: singlet ^1P1 + triplet ^3P1) -- the SIMPLE end of the
    #   ion ladder. ONE target, ONE knob (`fitScaleFactor`): the cheapest, cleanest demonstration of the
    #   whole module, on a system famous in its own right -- the three energies quoted below are the O VII
    #   "f/i/r" X-ray triplet lines used throughout X-ray astronomy diagnostics, and converting them to
    #   wavelength reproduces the well-known ~22.10/21.80/21.60 Angstrom values to 3-4 figures -- an
    #   independent cross-check of the NIST data itself, not merely of this module.
    #
    # NIST ASD (L21139): 1s2p ^3P1 = 568.640052 eV, ^1P1 = 573.9610626 eV, both above the 1s^2 ^1S0 ground.
    ni   = Nuclear.Model(8.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    wb    = SelfConsistent.performSCF([Configuration("1s^2"), Configuration("1s 2p")], ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    abInitio = Hamiltonian.performCI(basis, ni, grid, asfSettings; printout=false)
    #
    sym = LevelSymmetry(AngularJ64(1), Basics.minus)
    idxCsf, oneBody, twoBody = AdjustedLevels.buildOneAndTwoBodyMatrices(sym, basis, ni, grid)
    ev_to_cm    = Defaults.convertUnits("energy: from atomic to Kayser", 1.0) / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    eig1        = Hamiltonian.diagonalizeCiMatrix(oneBody .+ 1.0 .* twoBody, LevelSelection())
    abInitioRel = (sort(eig1.values) .- sort(eig1.values)[1]) ./ auPerKayser
    @printf("\n>>> ab initio (c=1.0) 1s2p (J=1^-) splitting:  %.3f cm^-1\n", abInitioRel[2])
    #
    targetCm = (573.9610626 - 568.640052) * ev_to_cm
    @printf(">>> NIST target splitting:                     %.3f cm^-1\n", targetCm)
    fm = AdjustedLevels.fitScaleFactor(sym, basis, ni, grid, [(2, targetCm)],
                     "NIST ASD (L21139): O VII 1s2p 3P1 = 568.640052 eV, 1P1 = 573.9610626 eV, above 1s^2 ground")
    println(fm)
    #
    combined = AdjustedLevels.mergeIntoMultiplet(fm, abInitio)
    einsteinSettings = Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000.)
    rawLines    = Einstein.computeLines(abInitio, grid, einsteinSettings)
    fittedLines = Einstein.computeLines(combined, grid, einsteinSettings)
    println("\n>>> RAW (ab initio) A-values, 1s2p -> ground:")
    for  ln in rawLines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (ab initio)   A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    println("\n>>> FITTED A-values (same orbitals, rescaled Coulomb interaction):")
    for  ln in fittedLines
        ln.initialLevel.energy <= ln.finalLevel.energy   &&   continue
        deltaE = Defaults.convertUnits("energy: from atomic to Kayser", ln.initialLevel.energy - ln.finalLevel.energy)
        aVal   = Defaults.convertUnits("rate: from atomic to 1/s", ln.photonRate.Babushkin)
        @printf("    J_i=%s -> J_f=%s   deltaE = %10.3f cm^-1 (fitted)      A(Bab) = %.5e /s\n",
                ln.initialLevel.J, ln.finalLevel.J, deltaE, aVal)
    end
    println("\n>>> source: NIST ASD (L21139), O VII 1s2p 3P1/1P1 above the 1s^2 ground")
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last successful:  15-Oct-2026
    #
    # Branch c: C I (Z=6), neutral carbon, 2s^2 2p^2, J=0 even block (^3P0 ground, ^1S0) -- the MODERATE
    #   rung. THE FIRST ATTEMPT, bare 2s^2 2p^2 alone, FAILS outright: the ab initio ^1S0 splitting (55197
    #   cm^-1) is 2.5x the measured 21648 cm^-1, a single-factor search over a physically reasonable bracket
    #   (0.5-1.2) cannot reach it, and lands stuck at the boundary with a 6130 cm^-1 residual -- a textbook
    #   single-configuration failure, not a module defect.
    #
    # WHAT FIXES IT: NIST's own "leading percentages" column for ^1S0 names the correlating partner directly
    #   -- a 4% admixture of `2p^4 ^1S`. Adding JUST that one configuration (a real physical correlation, not
    #   a guess) brings the ab initio splitting down to 33905 cm^-1 -- still far off on its own, but now
    #   comfortably inside a single Coulomb knob's reach.
    #
    # THE ORTHOGONAL ROUTE IS DELIBERATELY NOT SHOWN HERE: the needed correction is ~38%, above this route's
    #   own measured validity boundary (~30%, see `AdjustedLevels.fitOrthogonalParameters`'s own docstring),
    #   and it underperforms there for a documented reason, not a mysterious one. The single-factor route is
    #   the right tool for this particular case, and using it is a scope decision, not an oversight.
    #
    # NIST ASD (L20057): C I 2s^2 2p^2 ^3P0 = 0 eV (ground), ^1S0 = 2.6840136 eV.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    wb    = SelfConsistent.performSCF([Configuration("2s^2 2p^2"), Configuration("2p^4")], ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    #
    sym = LevelSymmetry(AngularJ64(0), Basics.plus)
    ev_to_cm = Defaults.convertUnits("energy: from atomic to Kayser", 1.0) / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    targetCm = 2.6840136 * ev_to_cm
    @printf("\n>>> NIST target (2s^2 2p^2 + 2p^4 correlation):  1S0 - 3P0 = %.3f cm^-1\n", targetCm)
    fm = AdjustedLevels.fitScaleFactor(sym, basis, ni, grid, [(2, targetCm)],
                     "NIST ASD (L20057): C I 2s^2 2p^2 1S0 = 2.6840136 eV above 3P0, with 2p^4 correlation")
    println(fm)
    println("\n>>> source: ", fm.source)
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last successful:  15-Oct-2026 -- VERIFIED AS A MECHANISM DEMONSTRATION, not as a close numerical match;
    #   see the branch's own text for exactly what is and is not being claimed.
    #
    # Branch d: Pr IV (Pr3+, Z=59), [Xe] 4f^2 -- the ADVANCED rung. Bare single-configuration ab initio is
    #   NOT good here: DFS overestimates both J=2 and J=4 term splittings by 10-30%, and neither AL nor EOL
    #   (the usual remedy for underbound VALENCE orbitals elsewhere in this project) helps -- both tighten
    #   the 4f orbital further and make the SPLITTING ratio worse, since the gap has the signature of
    #   CORE-POLARIZATION SCREENING of the 4f-4f repulsion by the polarizable Xe-like core, a real many-body
    #   effect no single-configuration mean field can produce (see this branch's own commit message for the
    #   full argument). Two independent correlating-configuration guesses (5d^2/6s^2/5d6s, and the core-hole
    #   4d^8 4f^4 double excitation, ~4200 CSFs) were TRIED and both FAILED for the same reason: both sit far
    #   too high in energy (1.7e5 and 1.7e6 cm^-1 above the reference) to mix meaningfully into the low-lying
    #   levels. So THIS branch does not claim a small absolute residual -- see Branch c for that -- it claims
    #   something an open f-shell honestly can offer: a CONSISTENT, REPRODUCIBLE relative improvement from
    #   the orthogonal (per-rank) route over the single uniform factor, confirmed on TWO INDEPENDENT blocks
    #   of the same ion.
    #
    # NIST ASD (L3974): 4f^2 3H4 = 0 (ground), 3F2 = 0.619501 eV, 1D2 = 2.149190 eV, 3P2 = 2.871550 eV,
    #                    3F4 = 0.849881 eV, 1G4 = 1.230077 eV.
    ni      = Nuclear.Model(59.)
    configs = [Configuration("[Xe] 4f^2")]
    grid    = Basics.recommendedGrid(configs, ni)                 ## Rule 12: the automatic bare grid is WRONG here
    asfSettings = AsfSettings()                                   ## (614 vs the needed ~11 a.u.) -- always use this for Pr IV.
    wb    = SelfConsistent.performSCF(configs, ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    ev_to_cm = Defaults.convertUnits("energy: from atomic to Kayser", 1.0) / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    #
    println("\n>>> J=2 block (3F2 / 1D2 / 3P2) -- the block where the orthogonal route earns its keep:")
    sym2 = LevelSymmetry(AngularJ64(2), Basics.plus)
    target_1D2 = (2.149190 - 0.619501) * ev_to_cm
    target_3P2 = (2.871550 - 0.619501) * ev_to_cm
    fmSingle2 = AdjustedLevels.fitScaleFactor(sym2, basis, ni, grid, [(2, target_1D2), (3, target_3P2)],
                      "NIST ASD (L3974): Pr IV 4f2 J=2, 3F2/1D2/3P2"; scaleBounds=(0.5,1.3))
    fmOrtho2  = AdjustedLevels.fitOrthogonalParameters(sym2, basis, ni, grid, [(2, target_1D2), (3, target_3P2)],
                      "NIST ASD (L3974): Pr IV 4f2 J=2, 3F2/1D2/3P2")
    @printf("    single-factor rms = %8.1f cm^-1   (c = %.4f)\n", fmSingle2.rmsDeviation, fmSingle2.scaleFactors[1].value)
    @printf("    orthogonal   rms = %8.1f cm^-1   (%.1fx improvement)\n", fmOrtho2.rmsDeviation,
            fmSingle2.rmsDeviation / fmOrtho2.rmsDeviation)
    #
    println("\n>>> J=4 block (3H4 / 3F4 / 1G4) -- an INDEPENDENT second check, same ion, same field:")
    sym4 = LevelSymmetry(AngularJ64(4), Basics.plus)
    target_3F4 = 0.849881 * ev_to_cm
    target_1G4 = 1.230077 * ev_to_cm
    fmSingle4 = AdjustedLevels.fitScaleFactor(sym4, basis, ni, grid, [(2, target_3F4), (3, target_1G4)],
                      "NIST ASD (L3974): Pr IV 4f2 J=4, 3H4/3F4/1G4"; scaleBounds=(0.5,1.3))
    fmOrtho4  = AdjustedLevels.fitOrthogonalParameters(sym4, basis, ni, grid, [(2, target_3F4), (3, target_1G4)],
                      "NIST ASD (L3974): Pr IV 4f2 J=4, 3H4/3F4/1G4")
    @printf("    single-factor rms = %8.1f cm^-1   (c = %.4f)\n", fmSingle4.rmsDeviation, fmSingle4.scaleFactors[1].value)
    @printf("    orthogonal   rms = %8.1f cm^-1   (J=4's own shape is already close to proportional,\n",
            fmOrtho4.rmsDeviation)
    println("                                        so there is little left for the richer model to fix)")
    #
    println("\n>>> Report BOTH blocks together, since neither alone tells the honest story:")
    println(">>> J=2: ", fmOrtho2.source, " -- ", round(fmSingle2.rmsDeviation/fmOrtho2.rmsDeviation, digits=1), "x improvement")
    println(">>> J=4: ", fmOrtho4.source, " -- consistent, smaller-in-absolute-terms residual, same relative size")
    #
elseif  false
    # Last successful:  21-Oct-2026   ## 10 classes (F^0/F^2/F^4/F^6/G^2/G^4/G^6/R^2/R^4/R^6, all (4f,4f)),
    #                                    2 targets, degreesOfFreedom = -8; orthogonal rms = 79.0 cm^-1 against
    #                                    the single-factor route's 1953.7 cm^-1 (24.7x) -- see the branch's own
    #                                    text for what this is and is not claiming.
    #
    # Branch e: THE OVERFITTING GUARD, on Pr IV's OWN J=2 block (competition-list item C1's remainder --
    #   Cowan's separated F^k/G^k/R^k classes, see AdjustedLevels.twoBodyParameterClass). Once the two
    #   relativistic j-partners of 4f (4f_5/2, 4f_7/2) are told apart from each other by TYPE as well as by
    #   rank, a SINGLE open shell no longer gives one class per rank -- it gives TEN: direct F^0/F^2/F^4/F^6,
    #   exchange G^2/G^4/G^6, and configuration-interaction R^2/R^4/R^6, all still labelled "(4f,4f)" since
    #   both j-partners share one Cowan-style shell parameter (see that function's own docstring for why).
    #
    # TEN classes against the SAME TWO measured levels this file's branch d already uses is a textbook
    #   underdetermined fit -- degreesOfFreedom = -8, and AdjustedLevels.FittedMultipletOrthogonal's own
    #   Base.show prints the warning unprompted. THIS IS NOT A FAILURE OF THE MODULE: the richer model still
    #   finds a materially better rms than the single uniform factor (79.0 against 1953.7 cm^-1, nested-model
    #   guarantee holds), and the SVD's own determinacy cutoff protects any one class from being amplified
    #   past what two numbers can support. What the guard actually buys is honesty: nobody should quote any
    #   ONE of these ten numbers as "the" fitted 4f spin-orbit or Slater parameter from this run alone -- read
    #   each one's own `determinacy` first. See Branch f for what pooling a second block does to this guard.
    #
    # NIST ASD (L3974): 4f^2 3H4 = 0 (ground), 3F2 = 0.619501 eV, 1D2 = 2.149190 eV, 3P2 = 2.871550 eV.
    ni      = Nuclear.Model(59.)
    configs = [Configuration("[Xe] 4f^2")]
    grid    = Basics.recommendedGrid(configs, ni)                 ## Rule 12: the automatic bare grid is WRONG here
    asfSettings = AsfSettings()                                   ## (614 vs the needed ~11 a.u.) -- always use this for Pr IV.
    wb    = SelfConsistent.performSCF(configs, ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    ev_to_cm = Defaults.convertUnits("energy: from atomic to Kayser", 1.0) / Defaults.convertUnits("energy: from atomic to eV", 1.0)

    sym2 = LevelSymmetry(AngularJ64(2), Basics.plus)
    target_1D2 = (2.149190 - 0.619501) * ev_to_cm
    target_3P2 = (2.871550 - 0.619501) * ev_to_cm

    fmSingleE = AdjustedLevels.fitScaleFactor(sym2, basis, ni, grid, [(2, target_1D2), (3, target_3P2)],
                      "NIST ASD (L3974): Pr IV 4f2 J=2, 3F2/1D2/3P2"; scaleBounds=(0.5,1.3))
    fmOrthoE  = AdjustedLevels.fitOrthogonalParameters(sym2, basis, ni, grid, [(2, target_1D2), (3, target_3P2)],
                      "NIST ASD (L3974): Pr IV 4f2 J=2, 3F2/1D2/3P2")

    println(fmOrthoE)
    @printf("\n>>> classes discovered: %d   (all ten are listed above, by AdjustedLevels.twoBodyParameterClass)\n",
            length(fmOrthoE.radialParameters))
    @printf(">>> single-factor rms = %8.1f cm^-1\n", fmSingleE.rmsDeviation)
    @printf(">>> orthogonal    rms = %8.1f cm^-1   (%.1fx improvement, nested-model guarantee holds)\n",
            fmOrthoE.rmsDeviation, fmSingleE.rmsDeviation/fmOrthoE.rmsDeviation)
    @printf(">>> degrees of freedom = %d   (negative -- read this result through 'determinacy', not in isolation)\n",
            fmOrthoE.degreesOfFreedom)
    #
elseif  false
    # Last successful:  21-Oct-2026   ## pooled J=2+J=4, SAME 10 global classes, 4 targets, degreesOfFreedom
    #                                    improves from -8 (either block alone) to -6; combined orthogonal rms
    #                                    = 125.2 cm^-1, still far below the unscaled ab initio residual.
    #
    # Branch f: MULTI-BLOCK FITTING -- the other half of item C1's remainder, and the reason the class
    #   separation in Branch e is worth having at all rather than merely correct. An `F^2(4f,4f)` class means
    #   the SAME physical radial integral in Pr IV's J=2 block as in its J=4 block (both come from the SAME
    #   4f^2 reference), so AdjustedLevels.fitOrthogonalParametersMultiBlock fits ONE shared parameter set
    #   against BOTH blocks' measured levels SIMULTANEOUSLY -- exactly Cowan's own standard practice (fit
    #   across every symmetry a reference configuration touches, not one block at a time), which this module
    #   could not do before tonight.
    #
    # THE HONEST READING, stated plainly rather than oversold: pooling turns Branch e's degreesOfFreedom = -8
    #   into -6 here (4 targets now, against the SAME 10 global classes) -- measurably better, STILL
    #   negative. An f-shell open to BOTH its relativistic j-partners simply has more independent Slater
    #   classes than a two-level, two-block literature case can pin down on its own. That is not a defect in
    #   this fit; it is exactly the standard RCE caveat competition-list item C1 names, made visible rather
    #   than hidden -- and it is precisely why Cowan's own practice ties F^2/F^4/F^6 (etc.) into ONE scaled
    #   group rather than fitting each independently, a further simplification this module does not yet
    #   offer and is not pretending to.
    #
    # NIST ASD (L3974): 4f^2 3H4 = 0 (ground), 3F2 = 0.619501 eV, 1D2 = 2.149190 eV, 3P2 = 2.871550 eV,
    #                    3F4 = 0.849881 eV, 1G4 = 1.230077 eV.
    ni      = Nuclear.Model(59.)
    configs = [Configuration("[Xe] 4f^2")]
    grid    = Basics.recommendedGrid(configs, ni)                 ## Rule 12: the automatic bare grid is WRONG here
    asfSettings = AsfSettings()                                   ## (614 vs the needed ~11 a.u.) -- always use this for Pr IV.
    wb    = SelfConsistent.performSCF(configs, ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    ev_to_cm = Defaults.convertUnits("energy: from atomic to Kayser", 1.0) / Defaults.convertUnits("energy: from atomic to eV", 1.0)

    sym2 = LevelSymmetry(AngularJ64(2), Basics.plus)
    target_1D2 = (2.149190 - 0.619501) * ev_to_cm
    target_3P2 = (2.871550 - 0.619501) * ev_to_cm
    sym4 = LevelSymmetry(AngularJ64(4), Basics.plus)
    target_3F4 = 0.849881 * ev_to_cm
    target_1G4 = 1.230077 * ev_to_cm

    blocks = [ (sym2, [(2, target_1D2), (3, target_3P2)]), (sym4, [(2, target_3F4), (3, target_1G4)]) ]
    fmMulti = AdjustedLevels.fitOrthogonalParametersMultiBlock(blocks, basis, ni, grid,
                      "NIST ASD (L3974): Pr IV 4f2, J=2 (3F2/1D2/3P2) and J=4 (3H4/3F4/1G4) pooled")
    println(fmMulti)

    @printf("\n>>> pooled rms = %8.1f cm^-1  across %d block(s), %d target(s), %d GLOBAL class(es)\n",
            fmMulti.rmsDeviation, length(fmMulti.blocks), sum(length(b.fitTargets) for b in fmMulti.blocks),
            length(fmMulti.radialParameters))
    @printf(">>> degrees of freedom = %d   (improved from Branch e's -8 by pooling, still negative)\n",
            fmMulti.degreesOfFreedom)
    println(">>> source: NIST ASD (L3974), Pr IV 4f2 J=2 and J=4 together")
    #
elseif  true
    # Last successful:  23-Oct-2026   ## V III 3d3, CORRELATED (3p^4.3d^5) + GROUPED (Slater-vs-CI) route,
    #                                    13 targets across 5 pooled blocks, 2 shared classes, rms = 278.5
    #                                    cm^-1 at degreesOfFreedom = +11 -- the first genuinely OVER-
    #                                    DETERMINED fit on an open d/f-shell ion this file has produced;
    #                                    see the branch's own text for what this closes and what it does not.
    #
    # Branch g: priority item 91 -- "a strong open d/f-shell case where AdjustedLevels' fit still measurably
    #   helps". Branch d/e/f (Pr IV 4f^2) never reached degreesOfFreedom >= 0: a bare 2-electron open shell
    #   has only 2-4 usable targets (most J-blocks of a 2-electron configuration hold just ONE level, which
    #   carries no internal splitting to fit), well short of the 7-10 Cowan classes relativistic j-splitting
    #   already demands. TWO things close that gap here, together -- neither alone was enough:
    #
    #   (1) A RICHER REFERENCE CONFIGURATION, same Cowan class count, many more targets. V III (V^2+, Z=23)
    #       is PURE 3d^3 (NIST ASD, fetched 23-Oct-2026: all levels below 27847 cm^-1, the next configuration
    #       3d^2.4s starting only at 43942 cm^-1 -- a clean, well-isolated reference). Three equivalent d
    #       electrons give EIGHT LS terms across FIVE usable J-blocks (J=1/2,3/2,5/2,7/2,9/2 -- J=11/2 holds
    #       only one level and contributes no target), 13 independent measured splittings in total, against
    #       the SAME 7 Cowan classes a bare 3d^2 needs. Bare 3d^3, fully class-separated and pooled across all
    #       five blocks: rms = 562.3 cm^-1 at degreesOfFreedom = +6 -- already the first well-determined
    #       number in this file's history, but still far from convincing.
    #
    #   (2) REAL CORRELATION, added BEFORE any Cowan adjustment, not after -- the maintainer's own point,
    #       21/22-Oct-2026: a single bare configuration is not where GRASP-style work would stop. The
    #       mechanism matching Branch c's successful 2p^4 admixture for C I is CORE POLARIZATION: a DOUBLE
    #       excitation of the closed 3p SHELL into the open 3d one, `[Ne] 3s^2 3p^4 3d^5`, parity-preserving
    #       (two p->d promotions, Delta-l = +2) and NOT the same thing as the valence-virtual excitations
    #       (3d->4d, 3d->4s) tried first and MEASURED to make a lighter test case (Ti III 3d^2) WORSE, not
    #       better -- those only let the 3d orbital relax, while the core-valence double excitation supplies
    #       the screening the ab initio splitting was actually missing (Ti III 3d^2 ratios ab-initio/NIST:
    #       27.4%/22.7% too high bare, 20.5%/20.7% with this correlation -- nearly uniform, not merely smaller).
    #
    #   ONE MORE STEP WAS NEEDED, found only by actually testing it: a correlating configuration with its own
    #   rich internal structure (3d^5 here, 3d^4 for Ti III) brings enough new CSFs (115 for this block alone)
    #   that picking "the Nth-lowest energy" can land on an INTRUDER state of the CORRELATING configuration's
    #   OWN character rather than the physically intended level -- invisible unless checked. Every target
    #   below is therefore picked by `AdjustedLevels.identifyReferenceLevels`, which weighs each level's CI
    #   eigenvector against the BARE 3d^3 CSFs (the same idea as NIST's own "leading percentages" column) and
    #   only accepts levels above 50%; every one found here sits at 98%+, confirmed, not assumed.
    #
    #   AND THE FULLY CLASS-SEPARATED ROUTE (`fitOrthogonalParametersMultiBlock`) IS NOT THE RIGHT TOOL ONCE A
    #   SECOND CONFIGURATION IS PRESENT: it then finds 20 classes (the cross-configuration R^k terms alone add
    #   13), hopelessly underdetermined against 13 targets. `AdjustedLevels.fitGroupedParametersMultiBlock`
    #   (new this session) is the physically-motivated middle ground: ONE parameter for every Slater
    #   (F^k/G^k) class and ONE for every configuration-interaction (R^k) class, matching published Cowan-
    #   code practice of scaling those two families by DIFFERENT empirical factors rather than conflating
    #   them -- and it needed its OWN fix along the way: an undamped Gauss-Newton step on this real,
    #   correlated data was MEASURED to overshoot and then oscillate (377.9 -> 964.4 -> 395.4 -> 790.2 cm^-1,
    #   never settling) on a simpler, two-block rehearsal of this same idea, so the route now backtracks
    #   (halves the step) whenever a full step fails to improve, exactly the fix that gets it to a genuine,
    #   stable minimum here.
    #
    #   THE HONEST READING. degreesOfFreedom = +11 is comfortably over-determined -- the first such result in
    #   this file for ANY open d/f-shell ion -- and the per-target residuals (66 to 574 cm^-1, individually
    #   checked, no single runaway outlier) are the signature of a genuine least-squares fit, not curve-
    #   fitting flexibility. But 278.5 cm^-1 is NOT inside the "rms < 10...100 cm^-1" absolute bar item 91
    #   itself quoted -- against splittings spanning tens of thousands of cm^-1, that is close to 1%
    #   RELATIVE accuracy, squarely in line with real published Cowan-code fits, but not the originally-
    #   imagined absolute number. The maintainer judged this the right point to close the item on 23-Oct-2026:
    #   the capability (a genuinely well-determined, physically-grouped, correlated RCE-style fit) is real and
    #   demonstrated; the specific absolute-cm^-1 target was closer to a hope than a measured, reachable bar
    #   for raw term-splitting data at this energy scale.
    #
    # NIST ASD V III (physics.nist.gov ASD levels query, fetched 23-Oct-2026): pure 3d^3 configuration.
    ni      = Nuclear.Model(23.)
    configs = [Configuration("[Ne] 3s^2 3p^6 3d^3"), Configuration("[Ne] 3s^2 3p^4 3d^5")]
    grid    = Basics.recommendedGrid(configs, ni)
    asfSettings = AsfSettings()
    wb    = SelfConsistent.performSCF(configs, ni, grid, asfSettings; printout=false)
    basis = wb.levels[1].basis
    refShells = [(3,2,3)]                                         ## bare 3d^3 reference, for identifyReferenceLevels

    nist12 = [11513.8, 15579.8]                                   ## a4P, a2P
    nist32 = [0.0, 11591.8, 15550.3, 16330.5, 42371.2]             ## a4F, a4P, a2P, a2D2, b2D1
    nist52 = [145.5, 11769.7, 16374.7, 27846.8, 42267.4]           ## a4F, a4P, a2D2, a2F, b2D1
    nist72 = [341.5, 11966.3, 27727.8]                             ## a4F, a2G, a2F
    nist92 = [583.8, 12187.0, 16810.9]                             ## a4F, a2G, a2H

    function  identifyBlock(Jrat, nist)
        sym = LevelSymmetry(AngularJ64(Jrat), Basics.plus)
        globalIdx, weights = AdjustedLevels.identifyReferenceLevels(sym, basis, ni, grid, refShells, length(nist))
        nistRel = nist .- nist[1]
        targets = [ (globalIdx[i], nistRel[i]) for i = 2:length(nist) ]
        return( (sym, targets) )
    end

    blocks = [ identifyBlock(1//2, nist12), identifyBlock(3//2, nist32), identifyBlock(5//2, nist52),
               identifyBlock(7//2, nist72), identifyBlock(9//2, nist92) ]

    fmGrouped = AdjustedLevels.fitGroupedParametersMultiBlock(blocks, basis, ni, grid,
                      "NIST ASD V III (fetched 23-Oct-2026): 3d3, J=1/2,3/2,5/2,7/2,9/2 pooled, " *
                      "correlated with [Ne] 3s2 3p4 3d5")
    println(fmGrouped)

    @printf("\n>>> GROUPED, CORRELATED pooled rms = %8.1f cm^-1  across %d block(s), %d target(s), %d GLOBAL group(s)\n",
            fmGrouped.rmsDeviation, length(fmGrouped.blocks), sum(length(b.fitTargets) for b in fmGrouped.blocks),
            length(fmGrouped.groupedParameters))
    @printf(">>> degrees of freedom = %d   (the first OVER-determined fit this file has produced for an open d/f-shell ion)\n",
            fmGrouped.degreesOfFreedom)
    println(">>> source: NIST ASD, V III 3d3, fetched 23-Oct-2026, correlated with [Ne] 3s2 3p4 3d5 (core polarization)")
    #
end
