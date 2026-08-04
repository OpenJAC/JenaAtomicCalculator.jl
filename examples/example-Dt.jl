
println("Dt) Apply & test the CrystalFieldEmission module (crystal-field-resolved multipole")
println("    transitions) with initial- and final-state multiplets declared directly via")
println("    Atomic.Computation's initialConfigs/finalConfigs/processSettings, exactly as for any")
println("    other JAC process (e.g. PhotoEmission, see example-Da.jl). This is the transition/")
println("    process-side companion to the level/property-side module-CrystalField.jl (see")
println("    example-Ck.jl); the two are deliberately kept in separate modules and separate example")
println("    branches (C = properties, D = processes), matching how JAC itself always separates")
println("    AbstractPropertySettings/Outcome from AbstractProcessSettings/Line and never mixes")
println("    the two within one module or one example branch.")

if  true
    # Last successful:  29-Jul-2026
    # Branch a: empty-lattice cross-check against a textbook, implementation-independent fact --
    #   NOT just "re-derive my own formula". With an empty CrystalField.Lattice, every CfLevel
    #   collapses to a single, definite (parentLevel, M) basis vector, so the crystal-field-resolved
    #   amplitudes must reduce exactly to ordinary (unsplit) E1 amplitudes. The textbook fact checked
    #   here is the classic "sodium D-line" degeneracy ratio: for np_(3/2) -> ns_(1/2) vs.
    #   np_(1/2) -> ns_(1/2), the total M-summed line strength ratio must be exactly
    #   (2*3/2+1)/(2*1/2+1) = 2.0, independent of Z or of any detail of this implementation.
    #   System: H(2p) [J=1/2,3/2] -> H(1s) [J=1/2], E1, Coulomb gauge, empty lattice.
    setDefaults("print summary: open", "zzz-CrystalFieldEmission.sum")
    #
    emptyLattice = CrystalField.Lattice()
    ceSettings = CrystalFieldEmission.Settings(CrystalFieldEmission.Settings(); multipoles=[E1], gauges=[UseCoulomb],
                                                lattice=emptyLattice, maxRank=1, includeJmixing=false, printBefore=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dt-a-H2p-to-H1s", grid=Radial.Grid(true),
                               nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
                               initialConfigs=[Configuration("2p")], finalConfigs=[Configuration("1s")],
                               processSettings=ceSettings)
    wb = perform(comp; output=true)
    lines = wb["crystal-field-resolved emission lines:"]
    CrystalFieldEmission.displayLines(stdout, lines)
    #
    lineStrength(line) = sum( abs(ch.amplitude)^2 for ch in line.channels )
    initialJ2(line)    = Basics.twice(line.initialCfLevel.cfBasis[1].parentLevel.J)
    strength32 = sum( lineStrength(line) for line in lines if initialJ2(line) == 3; init=0.0 )
    strength12 = sum( lineStrength(line) for line in lines if initialJ2(line) == 1; init=0.0 )
    println(">> Summed |amplitude|^2, 2p_3/2 -> 1s_1/2:  $strength32")
    println(">> Summed |amplitude|^2, 2p_1/2 -> 1s_1/2:  $strength12")
    println(">> Ratio (expected exactly 2.0, the textbook np_3/2/np_1/2 degeneracy ratio):  $(strength32/strength12)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  29-Jul-2026
    # Branch b: crystal-field-modified case, reusing the octahedral (Oh) lattice already validated
    #   in example-Ck.jl branches a/c. System: H(3d) -> H(2p), E1. Since 3d has l=2 and 2p has l=1,
    #   this is a real (allowed) E1 transition, now computed between Stark-split sublevels on both
    #   sides. No literature reference exists for this synthetic system; what is checked instead is
    #   a basis-rotation invariant that must hold regardless: the total transition strength, summed
    #   over ALL final Stark sublevels and starting from one fixed initial Stark sublevel, must equal
    #   the corresponding total computed directly between the UNSPLIT parent levels (since CfLevel.mc
    #   are unitary-eigenvector coefficients of a Hermitian matrix, summing |amplitude|^2 over a
    #   complete final basis is invariant under that unitary change of basis -- Parseval's theorem).
    setDefaults("print summary: open", "zzz-CrystalFieldEmission.sum")
    #
    charge = -2.0;   rho = 200.0
    ions   = [ CrystalField.PointCharge(charge, rho, 0.0,     0.0),
               CrystalField.PointCharge(charge, rho, pi,      0.0),
               CrystalField.PointCharge(charge, rho, pi/2,    0.0),
               CrystalField.PointCharge(charge, rho, pi/2,    pi),
               CrystalField.PointCharge(charge, rho, pi/2,    pi/2),
               CrystalField.PointCharge(charge, rho, pi/2,    3pi/2) ]
    lattice = CrystalField.Lattice(ions, "Oh")
    #
    ceSettings = CrystalFieldEmission.Settings(CrystalFieldEmission.Settings(); multipoles=[E1], gauges=[UseCoulomb],
                                                lattice=lattice, maxRank=4, includeJmixing=false, printBefore=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dt-b-H3d-to-H2p", grid=Radial.Grid(true),
                               nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
                               initialConfigs=[Configuration("3d")], finalConfigs=[Configuration("2p")],
                               processSettings=ceSettings)
    wb = perform(comp; output=true)
    lines = wb["crystal-field-resolved emission lines:"]
    CrystalFieldEmission.displayLines(stdout, lines)
    #
    # Basis-rotation-invariant cross-check (Parseval-type): FIX one specific, unambiguous initial
    # state -- one particular unsplit parent (level,M) sublevel of 3d, obtained from an EMPTY-lattice
    # initial CfMultiplet (so each of its CfLevels is trivially a single one-hot (parentLevel,M)
    # state) -- and compare the total strength to ALL final sublevels computed two ways: (A) against
    # the REAL (crystal-field-split) final CfMultiplet, and (B) against the EMPTY-lattice final
    # CfMultiplet. Since CfLevel.mc are unitary-eigenvector coefficients of a Hermitian matrix, the
    # final CfLevels in case (A) and the final (parentLevel,M) states in case (B) span the SAME final
    # Hilbert space related by a unitary rotation; summing |amplitude|^2 to a complete final basis
    # for one FIXED initial state must therefore give the identical total in both cases. Only the
    # final side's splitting differs between (A) and (B) -- the initial state is the same fixed
    # vector both times, which is essential for the invariant to hold.
    initialMultiplet = wb["initialMultiplet"];   finalMultiplet = wb["finalMultiplet"];   grid = comp.grid
    emptyLat = CrystalField.Lattice()
    emptyInitialSettings = CrystalField.Settings(CrystalField.Settings(); lattice=emptyLat, maxRank=1, includeJmixing=false)
    emptyInitialOutcomes = CrystalField.computeOutcomes(initialMultiplet, emptyLat, grid, emptyInitialSettings)
    fixedInitialCfLevel  = emptyInitialOutcomes[1].cfMultiplet.cfLevels[1]     # one specific, unsplit 3d (parentLevel,M) state
    #
    # both 2p levels (J=1/2 and J=3/2) must be included on the final side for the invariant to hold
    # -- otherwise the "complete final basis" requirement of the Parseval argument is violated.
    finalLevelSelection = LevelSelection(true, indices=[lev.index for lev in finalMultiplet.levels])
    finalOutcomesSplit = CrystalField.computeOutcomes(finalMultiplet, lattice, grid, CrystalField.Settings(CrystalField.Settings();
                                                        lattice=lattice, model=CrystalField.PointChargeModel(), maxRank=4,
                                                        levelSelection=finalLevelSelection))
    finalOutcomesEmpty = CrystalField.computeOutcomes(finalMultiplet, emptyLat, grid, CrystalField.Settings(CrystalField.Settings();
                                                        lattice=emptyLat, maxRank=1, levelSelection=finalLevelSelection))
    #
    strengthA = sum( abs(CrystalFieldEmission.computeLineAmplitude(E1, Basics.Coulomb, fCfLevel, fixedInitialCfLevel, grid))^2
                      for outcome in finalOutcomesSplit  for fCfLevel in outcome.cfMultiplet.cfLevels )
    strengthB = sum( abs(CrystalFieldEmission.computeLineAmplitude(E1, Basics.Coulomb, fCfLevel, fixedInitialCfLevel, grid))^2
                      for outcome in finalOutcomesEmpty  for fCfLevel in outcome.cfMultiplet.cfLevels )
    println(">> Total strength from one fixed initial (3d) sublevel to ALL final sublevels:")
    println("     (A) real crystal-field-split final (2p) manifold:   $strengthA")
    println("     (B) unsplit (empty-lattice) final (2p) manifold:    $strengthB")
    println("   These MUST match to numerical precision (Parseval/basis-rotation invariance).")
    #
    setDefaults("print summary: close", "")
    #
end
