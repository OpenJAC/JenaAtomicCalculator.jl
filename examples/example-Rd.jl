
println("Rd) Apply & test the AdjustedLevels module: fit ONE Coulomb scale factor for a single symmetry")
println("    block, merge it back into a full multiplet, and show the effect on a computed Einstein")
println("    A-value -- not just a nicer-looking energy table.")

setDefaults("print summary: open", "zzz-AdjustedLevels.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

if  true
    # Last visit:  14-Oct-2026
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
end
