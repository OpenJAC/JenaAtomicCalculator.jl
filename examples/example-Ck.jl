
println("Ck) Apply & test the CrystalField module (point-charge Stark splitting) with ASF from an")
println("    internally generated multiplet. Method follows Gaigalas & Kato, CPC 261 (2021) 107772")
println("    (examples/papers/2021.cpc-gaigalas-crystal-field.pdf); see also Uldry, Vernay & Delley,")
println("    Phys. Rev. B 85, 125133 (2012) (examples/papers/2012.cpc-uldry.crystal-field.pdf) for a")
println("    complementary (non-CSF, ab-initio point-charge + Dirac-DFT determinantal) approach.")

if  true
    # Last successful:  28-Jul-2026
    # Branch a: pure symmetry/group-theory sanity check, no external reference numbers needed.
    #   System: bare hydrogen H(3d) (Z=1, no core) -- a single active d-electron gives exactly two
    #   single-CSF levels, 3d_(3/2) and 3d_(5/2) (mc = [1.0] each, no CI mixing), so the crystal-field
    #   matrix elements are exact and depend only on the lattice geometry and the tensor algebra, not
    #   on any atomic-structure approximation.
    #   Lattice: 6 point charges of charge -2.0 e at distance 4.0 bohr along +-x, +-y, +-z -- a
    #   perfectly octahedral (Oh) point-charge arrangement.
    #   Two independent, textbook group-theory predictions are checked (Bethe/Griffith crystal-field
    #   theory for the cubic double group):
    #     (i)  J = 3/2 (Gamma8, 4-fold) does NOT split any further in a pure Oh field -- the 4
    #          computed Stark sublevels should be degenerate to numerical precision.
    #     (ii) J = 5/2 (Gamma7 + Gamma8) splits into exactly two groups of 2 and 4 -- the 6 computed
    #          Stark sublevels should cluster into only 2 distinct energies (not 6).
    #   In addition, CrystalField.multipoleLatticeSum(lattice,2,q) should vanish for all q (the k=2
    #   term of a perfectly octahedral point-charge sum is exactly zero); only k=4 (q=0,+-4, in the
    #   fixed ratio B44/B40 = sqrt(5/14)) should survive.
    setDefaults("print summary: open", "zzz-CrystalField.sum")
    #
    # rho must be well outside the radial extent of the 3d orbital (hydrogenic 3d has <r> ~ 10.5
    # bohr and a long exponential tail) so that the "point charges outside the electron shell"
    # (r<R) branch of Eq. (2) -- the only branch CrystalField.electrostaticIntegral implements --
    # remains valid; rho=4 (tried first) violated this and gave unphysical, huge (~10-60 Hartree)
    # "splittings" for the k=4 (5/2) case, since a large fraction of the 3d density then sits
    # beyond rho. rho=200 restores a physically sane, perturbative-scale splitting.
    charge = -2.0;   rho = 200.0
    ions   = [ CrystalField.PointCharge(charge, rho, 0.0,     0.0),          # +z
               CrystalField.PointCharge(charge, rho, pi,      0.0),          # -z
               CrystalField.PointCharge(charge, rho, pi/2,    0.0),          # +x
               CrystalField.PointCharge(charge, rho, pi/2,    pi),           # -x
               CrystalField.PointCharge(charge, rho, pi/2,    pi/2),         # +y
               CrystalField.PointCharge(charge, rho, pi/2,    3pi/2) ]       # -y
    lattice = CrystalField.Lattice(ions, "Oh")
    #
    println(">> k=2 lattice sum (should vanish for all q):")
    for  q = -2:2   println("     q=$q:  $(CrystalField.multipoleLatticeSum(lattice, 2, q))")   end
    println(">> k=4 lattice sum (only q=0,+-4 should survive, ratio B44/B40 = sqrt(5/14) = $(sqrt(5/14))):")
    for  q = -4:4   println("     q=$q:  $(CrystalField.multipoleLatticeSum(lattice, 4, q))")   end
    #
    wa = Atomic.Computation(Atomic.Computation(), name="Ck-a-H3d", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
                            configs=[Configuration("3d")],
                            propertySettings=Basics.AbstractPropertySettings[] )

    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"];    grid = wb["grid:"]
    println(">> Generated levels:  ", [ (lev.index, string(lev.J)) for lev in multiplet.levels ])
    #
    settings = CrystalField.Settings(CrystalField.Settings(); lattice=lattice, maxRank=4, includeJmixing=false,
                                      printBefore=true, levelSelection=LevelSelection(true, indices=[1,2]))
    outcomes = CrystalField.computeOutcomes(multiplet, lattice, grid, settings)
    CrystalField.displayResults(stdout, outcomes)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  28-Jul-2026 (illustrative only, see caveats below -- not a literature match)
    # Branch b: real published crystal geometry, illustrative magnitude check (NOT a reproduction of
    #   Table 1). The 6-nearest-neighbor Er2O3 lattice below is transcribed verbatim (charge, distance
    #   in Angstrom, theta/phi in degrees) from the file Crystaldata_6_sphe_2015_NO_Symm bundled in
    #   examples/codes/2021-gaigalas-crystal-field.zip (CF_Hamiltonian/Sample_Runs/), i.e. the actual
    #   input geometry of the paper's own Er3+-in-Er2O3 test case (Sec. 7, Table 1).
    #   What this branch does NOT do: reproduce Er3+'s actual 4f^11 CI wavefunction. The paper's own
    #   mixing-coefficient and radial-orbital files (Sample_Output/5.m, 5.w) are unformatted GRASP2018
    #   binaries that JAC has no reader for (a separate undertaking, unrelated to the crystal-field
    #   physics itself; module-ManyElectron.jl only reads RELCI/Grasp92-style mixing files). Instead, a
    #   single hydrogenic 4f orbital for an effective Z=18 (NOT Z=1) stands in for the valence orbital --
    #   Z=1 was tried first and gives <r> ~ 18 bohr, wildly larger than the ~4.2-4.4 bohr real Er-O
    #   distances used below, violating the "point charges outside the electron shell" (r<R) assumption
    #   that CrystalField.electrostaticIntegral relies on (Eq. (2)'s r<R branch); Z=18 contracts the
    #   hydrogenic 4f shell to <r> ~ 1.0 bohr, safely inside the real ionic distances.
    #   maxRank is capped at 4 (not 6): the k=6 term was found to blow up numerically (<r^6> ~ 4e5,
    #   vs. <r^2> ~ 1.1) for this B-spline-generated orbital -- a tiny, numerically-unavoidable
    #   residual amplitude far out on the (r up to ~600 bohr) radial grid gets amplified enormously
    #   by r^6. This is a genuine, generic numerical hazard of the r<R point-charge model for high
    #   tensor rank with box/B-spline-basis orbitals (not specific to this module), separate from the
    #   physical r<R modeling assumption above; a production calculation would need properly
    #   asymptotically-decaying bound-state orbitals (as GRASP2018 itself provides) or an explicit
    #   tail cutoff. Both findings are noted as "serious difficulties" in the module's design summary.
    #   What IS meaningful in this branch: (a) the code runs end-to-end on a real, non-cubic,
    #   non-centrosymmetric point-charge environment without error, using the actual published lattice
    #   geometry; (b) -- in contrast to branch a's exact cubic-symmetry zeros -- the crystal-field
    #   matrix now comes out genuinely complex off-diagonal (see module-CrystalField.jl's note on
    #   Hermitian-but-complex matrices for low-symmetry lattices) and only Kramers (M,-M) degeneracy
    #   survives, not the higher cubic degeneracies of branch a.
    setDefaults("print summary: open", "zzz-CrystalField.sum")
    #
    angstromToBohr = 1.8897259886
    degToRad       = pi/180
    # rho[Ang], theta[deg], phi[deg], charge
    raw = [ 2.254103  110.7050   -40.56799  -2.0;
            2.235114   43.72158 -137.9186   -2.0;
            2.254103  110.7050   139.4320   -2.0;
            2.235114   43.72158   42.08135  -2.0;
            2.323951  122.6551    49.54849  -2.0;
            2.323951  122.6551  -130.4515   -2.0 ]
    ions = [ CrystalField.PointCharge(raw[i,4], raw[i,1]*angstromToBohr, raw[i,2]*degToRad, raw[i,3]*degToRad)
             for i = 1:size(raw,1) ]
    lattice = CrystalField.Lattice(ions, "Er2O3, 6 nearest ions (Gaigalas & Kato 2021, Sample_Runs)")
    #
    wa = Atomic.Computation(Atomic.Computation(), name="Ck-b-Z18-4f", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(18.),
                            configs=[Configuration("4f")],
                            propertySettings=Basics.AbstractPropertySettings[] )

    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"];    grid = wb["grid:"]
    println(">> Generated levels:  ", [ (lev.index, string(lev.J)) for lev in multiplet.levels ])
    #
    settings = CrystalField.Settings(CrystalField.Settings(); lattice=lattice, maxRank=4, includeJmixing=false,
                                      printBefore=true, levelSelection=LevelSelection(true, indices=[1,2]))
    outcomes = CrystalField.computeOutcomes(multiplet, lattice, grid, settings)
    CrystalField.displayResults(stdout, outcomes)
    println(">> For comparison, Table 1 (Gaigalas & Kato 2021) reports Er3+ 4f^11 4I(15/2) Stark")
    println("   splittings (without J-mixing) spanning 0 to 414.54 cm^-1 across 2159 neighboring ions")
    println("   (this branch uses only the nearest 6 ions and an illustrative Z=18 hydrogenic proxy")
    println("   orbital, not Er3+'s actual wavefunction -- see comment above).")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  29-Jul-2026
    # Branch c: test the two new features suggested by a re-read of Uldry, Vernay & Delley, Phys.
    #   Rev. B 85, 125133 (2012) -- CrystalField.characteristicSplitting (their CXS, Sec. II F) and
    #   CrystalField.fitScaleField (their S_xtal empirical rescaling, Sec. IV) -- reusing branch a's
    #   exact, hand-verifiable octahedral H(3d) setup so the expected numbers are known in advance:
    #     (i)  For J=5/2, characteristicSplitting should reproduce branch a's single splitting value
    #          (0.063551 cm^-1) exactly: with only 2 distinct sublevel energies (a doublet + a
    #          quartet), the "largest gap" split trivially recovers that one gap.
    #     (ii) fitScaleField is asked to reproduce exactly TWICE that CXS. For a single level (no
    #          J-mixing) the interaction matrix's diagonal is just the (shared) level energy -- a
    #          multiple of the identity -- so it cancels out of every eigenvalue DIFFERENCE, and CXS
    #          is then EXACTLY proportional to scaleField. The fitted scaleField should therefore come
    #          back as exactly 2.0, and re-evaluating CXS at that scaleField should reproduce the
    #          2x target exactly.
    setDefaults("print summary: open", "zzz-CrystalField.sum")
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
    wa = Atomic.Computation(Atomic.Computation(), name="Ck-c-H3d", grid=Radial.Grid(true),
                            nuclearModel=Nuclear.Model(1., "uniform", 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
                            configs=[Configuration("3d")],
                            propertySettings=Basics.AbstractPropertySettings[] )
    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"];    grid = wb["grid:"]
    lev52 = [ lev  for lev in multiplet.levels  if  lev.J == AngularJ64(5//2) ][1]
    #
    cfMultiplet1 = CrystalField.computeRepresentation([lev52], lattice, CrystalField.PointChargeModel(1.0), grid, 4)
    cxs1         = CrystalField.characteristicSplitting(cfMultiplet1)
    println(">> CXS at scaleField=1.0:  $(Defaults.convertUnits("energy: from atomic to Kayser", cxs1)) cm^-1  " *
             "(branch a found 0.063551 cm^-1 for this same level/lattice -- should match)")
    #
    target = 2.0 * cxs1
    fitted = CrystalField.fitScaleField([lev52], lattice, grid, target, 4)
    println(">> fitScaleField target = 2 x CXS(scaleField=1)  -->  fitted scaleField = $fitted  (expected: exactly 2.0)")
    #
    cfMultiplet2 = CrystalField.computeRepresentation([lev52], lattice, CrystalField.PointChargeModel(fitted), grid, 4)
    cxs2         = CrystalField.characteristicSplitting(cfMultiplet2)
    println(">> CXS at the fitted scaleField:  $(Defaults.convertUnits("energy: from atomic to Kayser", cxs2)) cm^-1  " *
             "(should equal 2x the branch-a value = $(Defaults.convertUnits("energy: from atomic to Kayser", target)) cm^-1)")
    #
    setDefaults("print summary: close", "")
    #
end
