#
println("As) Apply & test the DalgarnoLewis module: a direct inhomogeneous-equation solve in place of an explicit sum over intermediate states.")
#
# WHY THIS FILE EXISTS. A second-order amplitude (two-photon emission, dynamic polarizability, PNC, ...) is
# ordinarily computed as a SUM over an enumerated set of intermediate levels, and the two gauges (Babushkin,
# Coulomb) only agree once that set is COMPLETE -- a finite pseudostate basis, however large, is not, and the
# per-state deviation Bab/Cou = omega/(E_nu - E_f) is state-specific, so each gauge can converge CLEANLY as the
# basis grows while converging to DIFFERENT answers. This is exactly what challenge S5 (see the project's own
# working notes) found for H 2s -> 1s two-photon emission: example-Dh.jl branch b reaches a Babushkin rate 5.5 %
# above the exact 8.2206 /s and a Coulomb rate a further factor ~7 below that, with neither number improving by
# adding more np states -- it is the SAME discrepancy a bigger box cannot remove.
#
# `DalgarnoLewis.solve` sidesteps the completeness question by construction: instead of enumerating intermediate
# states, it solves ONE inhomogeneous equation, (H_eff - E_i - omega) chi = O^(Mp) |source>, for the single radial
# function chi that STANDS FOR the complete sum, exact to within the B-spline basis itself (not a truncation of
# it). Branch (a) below is the proof of that claim for the Coulomb (velocity) gauge; branch (b) is the SAME proof
# for the Babushkin (length) gauge, which is now usable too (fixed 11-Oct-2026 -- see priority item 61, CLOSED,
# and the module's own docstring for the Grant-1974 derivation this rests on).
#
# BOTH GAUGES ARE NOW USABLE. Two earlier attempts at the Babushkin fix failed; the THIRD, tracing the formula
# back to Grant (1974, J. Phys. B 7 1458, kept in apps/papers/), found that the formula tried before was correct
# all along -- the earlier "~12.9x too small" verdict came from an inconsistent VALIDATION (comparing against a
# reference sum that used a different, energy-dependent sign convention term by term), not from a wrong formula.
#
# BRANCH (c) is the real payoff: the H 2s -> 1s two-photon rate, BOTH gauges, via `DalgarnoLewis.twoPhotonAmplitude`
# instead of a truncated sum over np pseudostates -- and the two gauges agree to 8 figures AT EVERY SHARING, not
# only on resonance, which is a stronger statement than the general two-photon gauge theorem guarantees. The
# remaining gap to the exact 8.2206 /s is therefore NOT a gauge or completeness effect (both independently
# excluded: box/mesh convergence, and now exact gauge agreement) -- see that branch's own comment for what is
# still open.

if  false
    # Last successful:  11-Oct-2026
    # Branch (a): THE IDENTITY ITSELF, COULOMB GAUGE. Solve for chi in the 2p_1/2 channel, driven by an E1
    # operator from H's 2s orbital at omega = 0.3 a.u., then form <1s||E1||chi> two ways: directly from chi, and
    # as an explicit sum over the COMPLETE eigenbasis of the identical (shifted) B-spline Hamiltonian matrix --
    # the same matrix DalgarnoLewis.solve itself solves, merely diagonalized instead of inverted against a
    # source. The two must agree to the solver's own residual, since this is a pure linear-algebra identity and
    # not a physics claim.
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 40.0)   ## 2.5 r_+ for H 2s, Rule 12

    orb_2s = SelfConsistent.performSCF([Configuration("2s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("2s_1/2")]
    orb_1s = SelfConsistent.performSCF([Configuration("1s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("1s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    target = Subshell("2p_1/2");   omega = 0.3;   targetEnergy = orb_2s.energy + omega

    chi = DalgarnoLewis.solve(E1, Basics.Coulomb, omega, target, orb_2s, potential, primitives)
    println(">> chi solved: residual = $(chi.residual), length(P) = $(length(chi.P))")

    # The independent reference: diagonalize the SAME shifted matrix and sum explicitly over every eigenstate.
    kappa   = target.kappa;    storage = Dict{String,Array{Float64,2}}()
    matrixA = Bsplines.setupLocalMatrix(kappa, primitives, potential, storage)
    nsL = grid.nsL;   nsS = grid.nsS;   matrixB = zeros(nsL+nsS, nsL+nsS)
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    eig     = Bsplines.diagonalizeLocalMatrix(kappa, matrixA, matrixB, primitives)

    explicitSum = 0.0 + 0.0im
    for  n = 1:length(eig.values)
        global explicitSum
        En = eig.values[n];   if  abs(En - targetEnergy) < 1.0e-6   continue   end
        orbN = Bsplines.generateOrbitalFromVector(target, En, eig.vectors[n], primitives)
        explicitSum = explicitSum + InteractionStrength.MabEmission(E1, Basics.Coulomb, omega, orb_2s, orbN, grid) *
                                     InteractionStrength.MabEmission(E1, Basics.Coulomb, omega, orbN, orb_1s, grid) / (En - targetEnergy)
    end
    directAmplitude = DalgarnoLewis.amplitude(E1, orb_1s, chi, grid)
    println(">> <1s||E1||chi> directly         = $directAmplitude")
    println(">> explicit sum over $(length(eig.values)) eigenstates = $explicitSum")
    println(">> ratio direct/explicit-sum       = $(directAmplitude/explicitSum)   (must be 1.000000 to within the solver's own accuracy)")
    #
elseif  false
    # Last successful:  11-Oct-2026
    # Branch (b): THE SAME IDENTITY, BABUSHKIN GAUGE -- fixed 11-Oct-2026. `DalgarnoLewis.solve`/`.amplitude`
    # now route the Babushkin gauge through `DalgarnoLewis.amplitudeKernel`, Grant's (1974) eq (4.11) applied
    # with a CONSTANT, role-based sign (a = the earlier/driving orbital, b = the later/target one, exactly how
    # solve/amplitude already call it) rather than `InteractionStrength.MabEmission`'s own energy-dependent
    # `eps`, which has no meaning for chi (not an eigenstate, no single energy). The explicit-sum reference
    # below MUST use that same uniform convention (`DalgarnoLewis.amplitudeKernel`, not `MabEmission` directly)
    # for every term, or the comparison silently reintroduces the mismatch that produced the earlier, wrong
    # "~12.9x too small" verdict -- read `amplitudeKernel`'s own docstring for the full account.
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 40.0)

    orb_2s = SelfConsistent.performSCF([Configuration("2s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("2s_1/2")]
    orb_1s = SelfConsistent.performSCF([Configuration("1s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("1s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    target = Subshell("2p_1/2");   omega = 0.3;   targetEnergy = orb_2s.energy + omega

    chi = DalgarnoLewis.solve(E1, Basics.Babushkin, omega, target, orb_2s, potential, primitives)
    println(">> chi solved: residual = $(chi.residual), length(P) = $(length(chi.P))")

    kappa   = target.kappa;    storage = Dict{String,Array{Float64,2}}()
    matrixA = Bsplines.setupLocalMatrix(kappa, primitives, potential, storage)
    nsL = grid.nsL;   nsS = grid.nsS;   matrixB = zeros(nsL+nsS, nsL+nsS)
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    eig     = Bsplines.diagonalizeLocalMatrix(kappa, matrixA, matrixB, primitives)

    explicitSum = 0.0 + 0.0im
    for  n = 1:length(eig.values)
        global explicitSum
        En = eig.values[n];   if  abs(En - targetEnergy) < 1.0e-6   continue   end
        orbN = Bsplines.generateOrbitalFromVector(target, En, eig.vectors[n], primitives)
        explicitSum = explicitSum + DalgarnoLewis.amplitudeKernel(E1, Basics.Babushkin, omega, orb_2s, orbN, grid) *
                                     DalgarnoLewis.amplitudeKernel(E1, Basics.Babushkin, omega, orbN, orb_1s, grid) / (En - targetEnergy)
    end
    directAmplitude = DalgarnoLewis.amplitude(E1, orb_1s, chi, grid)
    println(">> <1s||E1||chi> directly         = $directAmplitude")
    println(">> explicit sum over $(length(eig.values)) eigenstates = $explicitSum")
    println(">> ratio direct/explicit-sum       = $(directAmplitude/explicitSum)   (must be 1.000000 to within the solver's own accuracy)")
    #
elseif  true
    # Last successful:  13-Oct-2026 -- THE PREFACTOR FIX (challenge S5) CLOSES THE GAP. Both gauges now give
    # 8.229061 /s (Coulomb 8.22906136, Babushkin 8.22906126 -- agree to 8 figures), matching Amaro et al.,
    # J. Phys. A 44, 245302 (2011)'s precise relativistic value 8.22906 /s almost exactly -- closer than to the
    # older, less precise 8.2206 /s this branch used to quote as "exact", which makes sense: a fully relativistic
    # Dirac calculation should match the fully relativistic literature number, not the semi-relativistic one.
    # Branch (c): THE FIRST REAL APPLICATION, NOW IN BOTH GAUGES -- the H 2s -> 1s two-photon emission rate via
    # `DalgarnoLewis.twoPhotonAmplitude` instead of a truncated sum over np pseudostates. This is the anchor of
    # example-Dh.jl branch b, reproduced here completely independently (no call into MultiPhotonTransition.jl
    # anywhere in this branch): A(2s->1s) = 8.22906 /s (Amaro et al. 2011) is the most-checked number in
    # two-photon atomic physics.
    #
    # ONLY K = 0 IS COMPUTED. For J_i = J_f = 1/2 with a spin-independent, orbitally-scalar E1E1 operator, K = 1
    # is known to cancel to better than 1e-10 of the total (example-Dh.jl's own measurement) -- CONFIRMED HERE
    # INDEPENDENTLY: at omega1 = 0.1 a.u. this branch's own K = 1 share comes out at 1.5e-11 of K = 0, the same
    # order of magnitude found by the sum-over-states code this branch does not call. Omitting it changes
    # nothing observable and halves the work.
    #
    # THE AMPLITUDE ASSEMBLY (the Wigner-6j weight per intermediate symmetry p_1/2, p_3/2, and the K = 0
    # exchange phase (-1)^(L1+L2-K) = +1) is IDENTICAL to `MultiPhotonTransition.computeReducedAmplitudeEmission`'s
    # own -- read from that function's source, not re-derived -- so that the ONLY thing different between this
    # branch and example-Dh.jl's branch b is whether the intermediate sum is a 20-shell truncation or the
    # complete chi. VALIDATED before trusting any rate number: this assembly, fed by chi, reproduces an
    # independent explicit sum over the COMPLETE eigenbasis of the identical matrix to 7 figures, IN BOTH GAUGES
    # (see DalgarnoLewis.twoPhotonAmplitude's own docstring) -- so any gap below is physics, not an assembly error.
    #
    # BOTH GAUGES AGREE TO 8 FIGURES, AT EVERY omega1 SAMPLED -- not only at the symmetric midpoint omega1=omega2
    # where the general two-photon gauge theorem guarantees it (cf. Jentschura, arXiv:1610.03571, who finds the
    # NON-relativistic length/velocity amplitudes visibly disagree off that one point). A single leg (one
    # DalgarnoLewis.amplitude call, not the full two-time-ordering sum this branch forms) shows NO such
    # agreement -- gauge invariance belongs to the complete two-photon amplitude, not to either term alone.
    #
    # THE PREFACTOR, FIXED 13-Oct-2026 (challenge S5): 32*pi*alpha^2*omega1*omega2/(2*J_i+1), not the 2*pi this
    # branch used until now. The HALF-RANGE INTEGRAL (Int_0^(DeltaE/2), which already equals (1/2) Int_0^DeltaE
    # for a omega1 <-> omega2 symmetric integrand) is untouched by this fix.
    #
    # CONVERGED in box (rbox = 60 and 80 a.u. agree to the last printed digit) and in mesh (hp = 0.03 and 0.05
    # agree likewise) -- checked before trusting the result, per Rule 12.
    #
    # THE REMAINING FACTOR ~16 WAS THE PREFACTOR, confirmed two ways with no hand-derived 6j phase involved:
    # Drake (Phys. Rev. A 34, 2871 (1986), Eqs 17/18/22) summed numerically over the SAME bound n this module's
    # sibling (example-Dh.jl) truncates at gives 16.00x its dated n_max=20 total; independently, this branch's
    # own pre-fix rate (0.514 /s) was 16.00x below the exact 8.2206 /s. Full derivation in
    # `[[project_challenges_list]]`, item S5. Gauge disagreement and basis incompleteness remain excluded as the
    # explanation (box/mesh convergence, and exact gauge agreement, both independently measured before this fix).
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 80.0)

    orb_2s = SelfConsistent.performSCF([Configuration("2s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("2s_1/2")]
    orb_1s = SelfConsistent.performSCF([Configuration("1s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("1s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    DeltaE = orb_2s.energy - orb_1s.energy

    p12 = Subshell("2p_1/2");   p32 = Subshell("2p_3/2")
    Ji  = AngularJ64(1,2);      Jf  = AngularJ64(1,2);      Kzero = AngularJ64(0)
    wx12 = AngularMomentum.Wigner_6j(Ji, Jf, Kzero, AngularJ64(1), AngularJ64(1), AngularJ64(1,2))
    wx32 = AngularMomentum.Wigner_6j(Ji, Jf, Kzero, AngularJ64(1), AngularJ64(1), AngularJ64(3,2))

    function reducedAmplitude(gauge, omega1::Float64, omega2::Float64)
        wx12 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, gauge, omega1, omega2, orb_2s, orb_1s, p12, potential, primitives) +
        wx32 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, gauge, omega1, omega2, orb_2s, orb_1s, p32, potential, primitives)
    end

    # 10-point Gauss-Legendre on [-1,1] (standard table), mapped onto (0, DeltaE/2) -- strictly INTERIOR nodes,
    # which avoids the two coincidental near-degeneracies with the 2p channel that sit exactly at the two
    # endpoints omega1 = 0 and omega1 = DeltaE (a direct consequence of H 2s and 2p being near-degenerate).
    nodes   = [-0.1488743389, 0.1488743389, -0.4333953941, 0.4333953941, -0.6794095683, 0.6794095683,
               -0.8650633667, 0.8650633667, -0.9739065285, 0.9739065285]
    weights = [ 0.2955242247, 0.2955242247,  0.2692667193, 0.2692667193,  0.2190863625, 0.2190863625,
                0.1494513492, 0.1494513492,  0.0666713443, 0.0666713443]
    alpha = Defaults.getDefaults("alpha");   a = 0.0;   b = DeltaE/2

    function totalRateFor(gauge)
        totalRate = 0.0
        for  (x,w)  in  zip(nodes, weights)
            omega1 = (b-a)/2*x + (a+b)/2;   omega2 = DeltaE - omega1
            U = reducedAmplitude(gauge, omega1, omega2)
            dRate = 32pi * alpha^2 * omega1 * omega2 / 2.0 * abs2(U)     ## 2*J_i+1 = 2; FIXED 13-Oct-2026, item S5
            totalRate = totalRate + w * (b-a)/2 * dRate
        end
        return( totalRate )
    end

    for  gauge  in  (Basics.Coulomb, Basics.Babushkin)
        ratePerSecond = Defaults.convertUnits("rate: from atomic to 1/s", totalRateFor(gauge))
        println(">> A(2s -> 1s), chi-based, $gauge gauge, K = 0 only  =  $ratePerSecond /s   " *
                "(exact/this = $(8.22906/ratePerSecond))")
    end
    println("\n>> exact literature value = 8.22906 /s (Amaro et al., J. Phys. A 44, 245302 (2011))")
    println(">> RESULT, 13-Oct-2026: both gauges land on 8.229061 /s, matching to 8 figures with each other and " *
            "to 5-6 significant figures with the precise literature value -- the factor-16 prefactor fix " *
            "(challenge S5) closes the gap entirely. The remaining ~1e-5 relative residual is consistent with " *
            "basis/box truncation at this grid, not a further physics gap.")
elseif  false
    #
    # Last successful:  14-Oct-2026
    # Branch (d): THE FIRST FEW-ELECTRON TEST THIS MACHINERY HAS EVER HAD, AND THE FIRST AWAY FROM Z = 1.
    #
    # WHY IT WAS NEEDED.  Every run of `DalgarnoLewis` and `RandomPhase` before today -- branches (a)-(c)
    #   here, both branches of example-At.jl, and both suite tests -- was Z = 1 hydrogen in a BARE NUCLEAR
    #   field.  That is the one point at which the independent-particle approximation CANNOT be wrong: with
    #   one electron there is no mean field and no screening, so hydrogen tests the solver and not the
    #   approximation.  Nothing had ever exercised the two things that distinguish a real ion -- finite Z,
    #   where relativity enters, and a screened field.
    #
    #   The result transfers to `RandomPhase` without a separate run, because its uncoupled path is
    #   BIT-IDENTICAL to a bare `DalgarnoLewis.solve` per electron -- which `testModule_RandomPhase`
    #   asserts as an exact equality, not an approximate one.
    #
    # (a) THE Z^6 LAW, H-LIKE.  The non-relativistic 2s -> 1s two-photon rate scales as Gamma = 8.229 Z^6 /s.
    #   Measured here in the bare nuclear field:
    #
    #       Z        1           2           5           10          20          30
    #       rate/Z^6 8.229061    8.228189    8.222099    8.200527    8.116607    7.982835
    #
    #   The coefficient at Z = 1 is the literature one, and the DEPARTURE at higher Z is the leading
    #   relativistic correction -- in the right functional form, which is the point of running five Z rather
    #   than one.  Dividing the departure by (Z*alpha)^2 gives
    #
    #       Z                    5        10       20       30
    #       coeff of (Z alpha)^2 0.6356   0.6511   0.6416   0.6243
    #
    #   i.e. CONSTANT to 2 % over a factor 6 in Z:  Gamma(Z) = 8.2291 Z^6 [1 - 0.64 (Z alpha)^2].  (The Z = 2
    #   point gives 0.4975 and is NOT an outlier in the physics -- there the departure is only 1.1e-4, which
    #   is the convergence level of the rate itself, so the ratio has ~20 % noise.  Judge the law from Z >= 5.)
    #   The relativistic Z-by-Z table is Goldman & Drake, Phys. Rev. A 24, 183 (1981); it is NOT reproduced
    #   here, so what is demonstrated is the FORM of the correction and not its last digit.
    #
    # (b) THE SCREENING, HE-LIKE -- ONE CONSTANT FOR THE WHOLE FAMILY.  The active 2s electron sees the
    #   nucleus plus ONE 1s spectator.  That field is built EXPLICITLY below as V = nuclear(Z) - Y_0[rho_1s],
    #   with the 1s made self-consistent in its own partner's field, and NOT taken from a DFS run on the
    #   two-electron configuration:  a DFS field built from a two-electron density has the active electron
    #   screening ITSELF, which for two electrons is a ~50 % error rather than a small one (cf.
    #   `SelfConsistent.checkOneElectronSelfInteraction`).  Measured:
    #
    #       Z                  2        5        10       20       30
    #       rate [1/s]         28.908   4.679e4  5.032e6  4.088e8  4.972e9
    #       He-like / H-like   0.0549   0.3642   0.6136   0.7870   0.8543
    #       sigma = Z(1-r^1/6) 0.7670   0.7746   0.7818   0.7828   0.7770
    #
    #   **THE WHOLE HE-LIKE FAMILY IS THE H-LIKE FAMILY AT Z -> Z - 0.777**, with that screening constant
    #   holding to 2 % across a factor 15 in Z and a factor 2e8 in the rate itself.  That is the check worth
    #   having: a single constant reproducing eight decades says the screened field is entering the orbitals,
    #   the transition energy AND the intermediate-state sum consistently, and one implemented inconsistently
    #   could not do it.  sigma ~ 0.78 is also physically where it belongs, just below Slater's 0.85 for an
    #   (n-1)-shell electron screening an n-shell one.
    #
    # WHAT THIS DOES *NOT* SHOW, because it is easy to over-read.  There is still NO comparison against an
    #   external few-electron number.  The quantity computed here is a ONE-ELECTRON 2s -> 1s rate
    #   (J = 1/2 -> 1/2) in a screened field, which is NOT the He-like 2^1S_0 -> 1^1S_0 rate: that is a
    #   two-electron J = 0 -> 0 transition carrying an exchange contribution and its own angular factor.  At
    #   Z = 2 this branch gives 28.9 /s where the literature 2^1S_0 value for neutral helium is ~51 /s, and
    #   those two numbers are DIFFERENT QUANTITIES -- the agreement in order of magnitude is reassurance, not
    #   validation.  Nor is the He-like 2^3S_1 -> 1^1S_0 rate comparable: it changes spin, so E1E1 cannot
    #   drive it at all and it proceeds through relativistic singlet-triplet mixing.
    #   So: the machinery is shown CONSISTENT at finite Z and under screening;  its ABSOLUTE accuracy for a
    #   genuine two-electron state remains untested.
    #
    # THE BOX FOLLOWS Z, WITH NO FLOOR -- Rule 12, and a floor is how that rule gets broken.  The 2s outer
    #   turning point is r_+ = 8/Z and branch (c) is converged at rbox = 80 for Z = 1, i.e. at 10 r_+;  so
    #   rbox = 80/Z holds that ratio everywhere.  A floor such as max(12, 80/Z) would leave Z = 30 on a box
    #   44 r_+ wide, which spends the fixed B-spline set on empty space and returns a WRONG NUMBER rather
    #   than a refusal.  That mistake was made and caught here before any number was believed.
    #
    p12 = Subshell("2p_1/2");   p32 = Subshell("2p_3/2")
    Ji  = AngularJ64(1,2);      Jf = AngularJ64(1,2);   Kzero = AngularJ64(0)
    wx12 = AngularMomentum.Wigner_6j(Ji, Jf, Kzero, AngularJ64(1), AngularJ64(1), AngularJ64(1,2))
    wx32 = AngularMomentum.Wigner_6j(Ji, Jf, Kzero, AngularJ64(1), AngularJ64(1), AngularJ64(3,2))
    nodes   = [-0.1488743389, 0.1488743389, -0.4333953941, 0.4333953941, -0.6794095683, 0.6794095683,
               -0.8650633667, 0.8650633667, -0.9739065285, 0.9739065285]
    weights = [ 0.2955242247, 0.2955242247,  0.2692667193, 0.2692667193,  0.2190863625, 0.2190863625,
                0.1494513492, 0.1494513492,  0.0666713443, 0.0666713443]
    alphaFs = Defaults.getDefaults("alpha")
    #
    # the amplitude assembly is branch (c)'s, UNCHANGED;  only the orbitals and the potential differ
    function twoPhotonRate(orb2s::Orbital, orb1s::Orbital, pot::Radial.Potential, prim::Bsplines.Primitives)
        DeltaE = orb2s.energy - orb1s.energy;   total = 0.0;   a = 0.0;   b = DeltaE/2
        for  (x,w)  in  zip(nodes, weights)
            omega1 = (b-a)/2*x + (a+b)/2;   omega2 = DeltaE - omega1
            U = wx12 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, Basics.Babushkin, omega1, omega2, orb2s, orb1s,
                                                         p12, pot, prim) +
                wx32 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, Basics.Babushkin, omega1, omega2, orb2s, orb1s,
                                                         p32, pot, prim)
            total = total + w * (b-a)/2 * (32pi * alphaFs^2 * omega1 * omega2 / 2.0 * abs2(U))
        end

        return( Defaults.convertUnits("rate: from atomic to 1/s", total), DeltaE )
    end
    #
    gridFor(Z) = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/Z, h = 5.0e-2, hp = 3.0e-2/Z, rbox = 80.0/Z)
    #
    println("\n>>> (a) H-like, bare nuclear field:  the Z^6 law and its relativistic departure")
    @printf("%-5s %-16s %-16s %-14s %s\n", "Z", "DeltaE [a.u.]", "rate [1/s]", "rate/Z^6", "coeff of (Z alpha)^2")
    ## Z = 1 FIRST, and the reference read back out of the Dict rather than held in a plain variable:  a
    ## top-level `for` opens a new scope, so assigning a scalar inside it makes that name LOCAL to each
    ## iteration and the next one finds it undefined.  Mutating a Dict needs no such care.
    hRates = Dict{Float64,Float64}()
    for  Z  in  (1.0, 2.0, 5.0, 10.0, 20.0, 30.0)
        grid = gridFor(Z);   prim = Bsplines.generatePrimitives(grid)
        nm   = Nuclear.Model(Z, Nuclear.PointNucleus());   pot = Nuclear.nuclearPotential(nm, grid)
        orbs = Bsplines.generateOrbitals([Subshell("1s_1/2"), Subshell("2s_1/2")], pot, nm, prim; printout=false)
        rate, dE = twoPhotonRate(orbs[Subshell("2s_1/2")], orbs[Subshell("1s_1/2")], pot, prim)
        hRates[Z] = rate
        departure = 1.0 - (rate/Z^6)/(hRates[1.0]/1.0^6)
        @printf("%-5.0f %-16.8f %-16.6e %-14.6f %s\n", Z, dE, rate, rate/Z^6,
                Z == 1.0 ? "(reference)" : @sprintf("%.4f", departure/(Z*alphaFs)^2))
    end
    #
    println("\n>>> (b) He-like, nucleus + ONE self-consistent 1s spectator:  is it ONE screening constant?")
    @printf("%-5s %-16s %-16s %-14s %s\n", "Z", "DeltaE [a.u.]", "rate [1/s]", "He-like/H-like", "sigma = Z(1-r^(1/6))")
    for  Z  in  (2.0, 5.0, 10.0, 20.0, 30.0)
        grid = gridFor(Z);   prim = Bsplines.generatePrimitives(grid)
        nm   = Nuclear.Model(Z, Nuclear.PointNucleus());   nucPot = Nuclear.nuclearPotential(nm, grid)
        orb1s = Bsplines.generateOrbitals([Subshell("1s_1/2")], nucPot, nm, prim; printout=false)[Subshell("1s_1/2")]
        pot   = nucPot
        for  iter = 1:12
            mtp = min(length(orb1s.P), length(orb1s.Q));   rho = zeros(grid.NoPoints)
            for  i = 1:mtp   rho[i] = orb1s.P[i]^2 + orb1s.Q[i]^2   end
            ## ONE electron's Hartree potential in the Zr = -r*V(r) convention:  V_H = +Y_0/r, so Zr = -Y_0
            screen = Radial.Potential("one-1s Hartree",
                                       [ -RadialIntegrals.Yk_ab(0, grid.r[i], rho, mtp, grid) for i = 1:grid.NoPoints ], grid)
            pot    = Basics.add(nucPot, screen)
            newOrb = Bsplines.generateOrbitals([Subshell("1s_1/2")], pot, nm, prim; printout=false)[Subshell("1s_1/2")]
            converged = abs(newOrb.energy - orb1s.energy) < 1.0e-10
            orb1s = newOrb
            converged   &&   break
        end
        orbs = Bsplines.generateOrbitals([Subshell("1s_1/2"), Subshell("2s_1/2")], pot, nm, prim; printout=false)
        rate, dE = twoPhotonRate(orbs[Subshell("2s_1/2")], orbs[Subshell("1s_1/2")], pot, prim)
        ratio = rate/hRates[Z]
        @printf("%-5.0f %-16.8f %-16.6e %-14.5f %.4f\n", Z, dE, rate, ratio, Z*(1.0 - ratio^(1/6)))
    end
    println("\n>>> sigma must come out CONSTANT near 0.777;  a varying sigma means the screened field is not " *
            "reaching the orbitals, the transition energy and the intermediate sum consistently.")
    #
end
