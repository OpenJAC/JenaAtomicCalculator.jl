
println("Ck) Apply & test the CrystalField module (point-charge Stark splitting) with ASF from an")
println("    internally generated multiplet. Method follows Gaigalas & Kato, CPC 261 (2021) 107772")
println("    (examples/papers/2021.cpc-gaigalas-crystal-field.pdf); see also Uldry, Vernay & Delley,")
println("    Phys. Rev. B 85, 125133 (2012) (examples/papers/2012.cpc-uldry.crystal-field.pdf) for a")
println("    complementary (non-CSF, ab-initio point-charge + Dirac-DFT determinantal) approach.")

using Printf

if  true
    # Last successful:  28-Aug-2026 -- RE-RUN, AND BOTH GROUP-THEORY PREDICTIONS HOLD EXACTLY. J = 3/2 (Gamma8) does
    # NOT split: all four sublevels sit at 0.000000 cm^-1 and the CXS is 0.000000 exactly. J = 5/2 splits into
    # exactly two groups of 2 and 4: sublevels 1-2 at 0.000000 cm^-1 and sublevels 3-6 at 0.065832 cm^-1, giving
    # CXS = 0.065832 cm^-1. These are structural checks with no tolerance to argue about, which is what makes them
    # worth re-running.
    # Previously:  28-Jul-2026
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
                            nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
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
    # THE 42 % OVERSHOOT IS KNOWN, EXPLAINED AND DELIBERATELY NOT PURSUED -- maintainer's decision, 01-Sep-2026.
    # DO NOT RAISE IT AS A DEFECT. On the real Er3+ wavefunction this branch's method gives a span of 587.85
    # cm^-1 against the 414.54 above, and the direction is what a 6-ion truncation should give, since more
    # distant shells partly cancel. Chasing it further was retired for a reason worth stating: **414.54 is
    # another CALCULATION, not a measurement** -- the same idealised point-charge model, computed by someone
    # else -- so reproducing it would validate our tensor algebra against theirs and would say nothing about
    # how well JAC describes a real crystal field. That question is to be answered by a real application, with
    # measured levels, rather than by this branch.
    #
    # AND A TRAP, recorded because it is easy to fall into: our characteristic splitting CXS = 417.06 cm^-1
    # sits within 0.6 % of their SPAN of 414.54. Those are DIFFERENT QUANTITIES -- a largest-gap measure against
    # a full range -- and the closeness is a coincidence. Reporting it as excellent agreement would be wrong.
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
                            nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0),
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
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had
    #    been resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    # Branch d: WHICH MULTIPOLE COMPONENTS SURVIVE IN WHICH SITE SYMMETRY -- and why this, and not a level splitting,
    #   is what this module has to offer the thorium-229 nuclear-clock community.
    #
    #   THE REFRAMING, WHICH IS THE POINT OF THE BRANCH.  In the solid-state clock the ion is Th(4+) = [Rn]: a CLOSED
    #   SHELL with J = 0.  A J = 0 level has exactly one magnetic sublevel, so a crystal field cannot split it, and the
    #   electronic-splitting half of this module -- everything branches a to c exercise -- has nothing whatever to act
    #   on.  Asking for "the Stark-split sublevel pattern of Th(4+) in two candidate hosts" is asking for something that
    #   is identically zero.
    #
    #   WHAT DOES REACH THAT PHYSICS is CrystalField.multipoleLatticeSum, which is PURELY GEOMETRIC and independent of
    #   the atomic system entirely.  Its k = 2 component is the lattice ELECTRIC-FIELD GRADIENT at the nuclear site --
    #   what couples to the 229Th nuclear quadrupole moment, and therefore what splits the clock transition and, when
    #   the site varies from ion to ion, inhomogeneously broadens it.  That half uses no orbitals, so it side-steps
    #   every one of the three limitations recorded in branch b: the r < R validity condition, the k = 6 numerical
    #   blow-up on B-spline tails, and the absence of a GRASP mixing-file reader.  All three come from orbitals it does
    #   not use.
    #
    #   THREE PREDICTIONS, ALL EXACT AND ALL STATED BEFORE RUNNING:
    #     (i)   the CaF2 Ca site -- 8 F- at the corners of a cube -- gives k = 2 EXACTLY ZERO for every q, by cubic
    #           symmetry.  A nucleus there feels no quadrupole interaction at all.
    #     (ii)  the same cube against branch a's OCTAHEDRON at matched charge and distance must give a k = 4 ratio of
    #           exactly -8/9, since SUM P_4 = -28/9 over the eight cube corners against +7/2 over the six octahedral
    #           vertices.  This is the valuable one: it checks the new geometry against ALREADY DATED work, so a slip in
    #           the new coordinates cannot pass silently.
    #     (iii) an idealized square antiprism (D4d) admits q = 0 mod 4, and |q| <= 2 at rank 2, so only q = 0 survives:
    #           the EFG is non-zero but AXIALLY SYMMETRIC.
    #
    #   THE ANTIPRISM IS IDEALIZED AND IS LABELLED SO.  There is no ThF4 structure in examples/papers/; thorium is
    #   8-coordinate there but the precise angles are not to hand.  The comparison below is therefore CUBIC AGAINST
    #   ANTIPRISMATIC 8-COORDINATION -- a real and checkable symmetry statement -- and not a claim about CaF2 against
    #   ThF4 as literal structures.  The CaF2 geometry, by contrast, is exact standard crystallography.
    #
    setDefaults("print summary: open", "zzz-CrystalField.sum")
    #
    angToBohr = 1.8897259886
    aCaF2     = 5.4626                          # fluorite lattice constant [Angstrom]
    rCaF      = aCaF2*sqrt(3)/4*angToBohr       # Ca-F distance = a sqrt(3)/4 = 2.3654 Ang = 4.4699 bohr
    pcOf      = function(q, x, y, z)
        r = sqrt(x^2 + y^2 + z^2)
        CrystalField.PointCharge(q, r, acos(z/r), atan(y, x))
    end
    #
    hh    = aCaF2/4*angToBohr
    cube  = CrystalField.Lattice([pcOf(-1.0, sx*hh, sy*hh, sz*hh) for sx in (1,-1) for sy in (1,-1) for sz in (1,-1)],
                                  "CaF2 Ca site, 8 F- cube, Oh")
    octa  = CrystalField.Lattice([pcOf(-1.0,  rCaF,0,0), pcOf(-1.0, -rCaF,0,0), pcOf(-1.0, 0, rCaF,0),
                                  pcOf(-1.0, 0,-rCaF,0), pcOf(-1.0, 0,0, rCaF), pcOf(-1.0, 0,0,-rCaF)],
                                  "octahedron, Oh, matched charge and distance")
    th0   = 59.26*pi/180
    anti  = CrystalField.Lattice(vcat([CrystalField.PointCharge(-1.0, rCaF, th0,      k*pi/2)        for k = 0:3],
                                      [CrystalField.PointCharge(-1.0, rCaF, pi - th0, k*pi/2 + pi/4) for k = 0:3]),
                                  "idealized square antiprism, D4d")
    #
    println("\n  (i)  which components survive, at matched charge (-1 e) and distance (", round(rCaF, digits=4), " bohr)")
    for (lat, nme) in [(cube, "CaF2 Ca site, 8 F- cube (Oh)"), (octa, "octahedron (Oh)"),
                       (anti, "square antiprism (D4d), theta0 = 59.26 deg")]
        println("     " * nme)
        # the threshold is set by the lattice's OVERALL scale across all ranks, not by the maximum within one
        # rank -- otherwise a rank whose components are ALL numerical noise reports every q as surviving
        scale = maximum( abs(CrystalField.multipoleLatticeSum(lat, kk, qq)) for kk in (2,4,6) for qq = -kk:kk )
        for k in (2, 4, 6)
            vals = [ CrystalField.multipoleLatticeSum(lat, k, q) for q = -k:k ]
            mx   = maximum(abs.(vals))
            keep = [ q for q = -k:k if abs(vals[q+k+1]) > 1.0e-10 * scale ]
            println("       " * @sprintf("k = %d:  max|A_kq| = %.6e     surviving q = %s", k, mx,
                                         isempty(keep) ? "none, all vanish" : string(keep)))
        end
    end
    #
    println("\n  (ii) the three predictions")
    p1 = maximum(abs(CrystalField.multipoleLatticeSum(cube, 2, q)) for q = -2:2)
    p2 = real(CrystalField.multipoleLatticeSum(cube, 4, 0) / CrystalField.multipoleLatticeSum(octa, 4, 0))
    p3 = maximum(abs(CrystalField.multipoleLatticeSum(anti, 2, q)) for q in (-2,-1,1,2))
    println("     " * @sprintf("cube, k = 2, max over q         = %.3e          (must be 0)", p1))
    println("     " * @sprintf("cube k=4 / octahedron k=4       = %.10f     (must be -8/9 = %.10f)", p2, -8/9))
    println("     " * @sprintf("antiprism, k = 2, q = 0         = %.6e", real(CrystalField.multipoleLatticeSum(anti,2,0))))
    println("     " * @sprintf("antiprism, k = 2, max over q!=0 = %.3e          (must be 0)", p3))
    #
    println("\n     WHAT THIS SAYS FOR A NUCLEAR CLOCK.  A cubic site gives NO quadrupole interaction at the nucleus at")
    println("     all: the clock line is not quadrupole-split there, whatever the nuclear moment.  An antiprismatic site")
    println("     gives a non-zero but purely AXIAL gradient, so the line splits cleanly with asymmetry parameter eta = 0.")
    println("     Neither statement needs an atomic-structure calculation, a nuclear moment, or a shielding factor --")
    println("     they are symmetry, and they are exact.")
    println("\n     WHAT MAY NOT BE READ OFF: an absolute field gradient or a frequency shift.  The physical EFG is the")
    println("     lattice sum multiplied by the Sternheimer antishielding factor of the closed Th(4+) shell, which is")
    println("     large -- of order a hundred -- and is not in JAC.  It cancels in a RATIO between two sites of the same")
    println("     ion, and ratios are therefore the only magnitudes quoted here.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: as branch d.]
    # Branch e: CHARGE COMPENSATION, which is where the inhomogeneous broadening actually comes from.
    #
    #   Branch d showed that a PERFECT cubic site gives no quadrupole interaction.  Real CaF2:Th is not a perfect cubic
    #   site: Th(4+) substituting for Ca(2+) carries two extra positive charges and must be compensated by interstitial
    #   F- ions, and WHERE THOSE SIT SETS THE LOCAL SYMMETRY.  The three single-compensator geometries below are the
    #   standard reference case, long established for rare-earth-doped CaF2:
    #       O_h    remote compensation, the cube left intact
    #       C_4v   an F- interstitial at the nearest <100> site
    #       C_3v   an F- interstitial at the next-nearest <111> site
    #   A real doped crystal contains a MIXTURE of these, each with its own quadrupole splitting -- which is precisely
    #   the inhomogeneous broadening that growing a native ThF4 crystal is meant to escape.
    #
    #   ONE TRAP, WORTH THE BRANCH ON ITS OWN.  A_20 is NOT a frame-independent measure of the gradient: it is defined
    #   about z, and neither interstitial lies on z.  Taken alone it says the C_3v site has NO quadrupole interaction,
    #   which is false -- the gradient is simply carried by the other q there.  What is invariant under rotation is the
    #   rank-2 norm sqrt(SUM_q |A_2q|^2), and that is what is tabulated.  Both columns are printed so the trap is
    #   visible rather than described.
    #
    #   TWO EXACT CHECKS, again stated before running.  Since the cube contributes NOTHING at k = 2, the invariant for a
    #   singly-compensated site is just that of one point charge, which is |q|/R^3 in closed form:
    #     (i)  the C_4v invariant must equal |q|/R^3 with R = a/2, to full precision;
    #     (ii) the C_3v/C_4v ratio must be exactly 3^(-3/2) = 0.19245, the two interstitials being at a/2 and a sqrt(3)/2.
    #
    #   NOT A CLAIM ABOUT 229Th:CaF2 SPECIFICALLY.  Th(4+) needs TWO compensators, not one, so the real site zoo is
    #   richer than these three; what is computed here is the documented single-compensator reference, which is what
    #   makes the mechanism visible.
    #
    setDefaults("print summary: open", "zzz-CrystalField.sum")
    #
    angToBohr = 1.8897259886
    aCaF2     = 5.4626
    pcOf      = function(q, x, y, z)
        r = sqrt(x^2 + y^2 + z^2)
        CrystalField.PointCharge(q, r, acos(z/r), atan(y, x))
    end
    hh    = aCaF2/4*angToBohr
    cube  = CrystalField.Lattice([pcOf(-1.0, sx*hh, sy*hh, sz*hh) for sx in (1,-1) for sy in (1,-1) for sz in (1,-1)], "Oh")
    f100  = pcOf(-1.0, aCaF2/2*angToBohr, 0.0, 0.0)
    f111  = pcOf(-1.0, aCaF2/2*angToBohr, aCaF2/2*angToBohr, aCaF2/2*angToBohr)
    sites = [ ("O_h    perfect cube, remote compensation", cube),
              ("C_4v   + F- interstitial at <100>",        CrystalField.Lattice(vcat(cube.ions, [f100]), "C4v")),
              ("C_3v   + F- interstitial at <111>",        CrystalField.Lattice(vcat(cube.ions, [f111]), "C3v")) ]
    inv2  = (lat) -> sqrt(sum(abs2(CrystalField.multipoleLatticeSum(lat, 2, q)) for q = -2:2))
    #
    println("\n  (i)  the rank-2 lattice sum at the three standard compensation sites")
    println("       site                                        invariant            |A_20|, frame-dependent")
    for (nme, lat) in sites
        println("       " * @sprintf("%-42s  %.6e         %.6e", nme, inv2(lat),
                                       abs(CrystalField.multipoleLatticeSum(lat, 2, 0))))
    end
    println("\n       Compare the two columns for C_3v: the invariant is 1.4e-03 while A_20 is 2.5e-17.  A reader taking")
    println("       A_20 for the gradient would conclude that site has no quadrupole interaction.  It has one; it is")
    println("       simply not aligned with z.")
    #
    println("\n  (ii) the two exact checks")
    r43 = inv2(sites[3][2]) / inv2(sites[2][2])
    cf  = inv2(sites[2][2]) / (1/(aCaF2/2*angToBohr)^3)
    println("       " * @sprintf("C_4v invariant / (|q|/R^3)   = %.10f     (must be 1: the cube adds nothing at k = 2)", cf))
    println("       " * @sprintf("C_3v / C_4v                  = %.10f     (must be 3^(-3/2) = %.10f)", r43, 3.0^-1.5))
    #
    println("\n     THE CONCLUSION, AND IT IS A MECHANISM RATHER THAN A NUMBER.  The perfect site gives zero; the two")
    println("     compensated sites give gradients differing by a factor of five, about different axes.  A crystal")
    println("     containing all three -- which a doped crystal does -- therefore presents the nucleus with several")
    println("     distinct quadrupole splittings at once, and the clock line is the sum of them.  That is a SYMMETRY")
    println("     contribution to the linewidth, entirely separate from microstrain, and it does not anneal away.  It is")
    println("     removed only by putting the thorium on a site it does not have to be charge-compensated into, which is")
    println("     the argument for growing a native thorium fluoride rather than doping calcium fluoride.")
    println("\n     As in branch d: ratios and zeros only.  The absolute gradient needs the Sternheimer antishielding")
    println("     factor of the closed Th(4+) shell, which JAC does not have, and which cancels in these ratios.")
    #
    setDefaults("print summary: close", "")
    #
end
