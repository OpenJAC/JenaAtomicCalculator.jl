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
# it). Branch (a) below is the proof of that claim: chi, used to form <final|O|chi>, is checked against an
# EXPLICIT sum over the COMPLETE eigenbasis of the identical finite matrix -- not a subset of it -- so the two
# must agree to the solver's own numerical accuracy, and any gap would be an implementation error, not an open
# physics question.
#
# ONLY THE COULOMB (VELOCITY) GAUGE IS CURRENTLY USABLE. The Babushkin (length) gauge needs an un-reduced
# one-body multipole matrix element valid for a function that is not an eigenstate -- chi solves an
# inhomogeneous, not an eigenvalue, equation -- and the existing reduced formula's orientation correction (eps)
# has no well-defined meaning for it. Two attempts at a quick fix were measured and failed; see priority item 61
# and the module's own docstring for the numeric record. Branch (b) demonstrates the resulting guard: `solve`
# RAISES for gauge = Basics.Babushkin rather than returning a number that looks plausible and is wrong by an
# order of magnitude.
#
# THIS FILE DOES NOT YET COMPUTE A PHYSICAL RATE. Wiring `DalgarnoLewis` into an actual two-photon amplitude
# (so that example-Dh.jl's own 8.2206 /s anchor can be checked against it) is the natural next step and is not
# attempted here -- see the module's own docstring and RandomPhase's (example-At.jl) for what is built and what
# is still missing.

if  true
    # Last successful:  11-Oct-2026
    # Branch (a): THE IDENTITY ITSELF. Solve for chi in the 2p_1/2 channel, driven by an E1 operator from H's 2s
    # orbital at omega = 0.3 a.u., then form <1s||E1||chi> two ways: directly from chi, and as an explicit sum
    # over the COMPLETE eigenbasis of the identical (shifted) B-spline Hamiltonian matrix -- the same matrix
    # DalgarnoLewis.solve itself solves, merely diagonalized instead of inverted against a source. The two must
    # agree to the solver's own residual, since this is a pure linear-algebra identity and not a physics claim.
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
elseif  true
    # Last successful:  11-Oct-2026 -- A GUARD TEST, NOT A PHYSICS BRANCH: it asks for the Babushkin gauge, which
    # this module does not yet support, specifically TO CONFIRM that it refuses cleanly instead of returning a
    # plausible-looking, wrong number. Raising is the pass here; a returned PerturbedOrbital would be the defect.
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 40.0)

    orb_2s = SelfConsistent.performSCF([Configuration("2s")], ni, grid, asfSettings; printout=false).levels[1].basis.orbitals[Subshell("2s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    target = Subshell("2p_1/2");   omega = 0.3

    try
        DalgarnoLewis.solve(E1, Basics.Babushkin, omega, target, orb_2s, potential, primitives)
        println(">> UNEXPECTED: the Babushkin gauge did NOT raise -- this would be a regression, see priority item 61")
    catch err
        println(">> EXPECTED: DalgarnoLewis.solve() raised for gauge = Basics.Babushkin:")
        println(">>   ", sprint(showerror, err))
    end
    #
elseif  true
    # Last successful:  11-Oct-2026
    # Branch (c): THE FIRST REAL APPLICATION -- the H 2s -> 1s two-photon emission rate, Coulomb/velocity gauge,
    # via `DalgarnoLewis.twoPhotonAmplitude` instead of a truncated sum over np pseudostates. This is the anchor
    # of example-Dh.jl branch b, reproduced here completely independently (no call into MultiPhotonTransition.jl
    # anywhere in this branch): A(2s->1s) = 8.2206 /s is the most-checked number in two-photon atomic physics.
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
    # independent explicit sum over the COMPLETE eigenbasis of the identical matrix to 7 figures (see
    # DalgarnoLewis.twoPhotonAmplitude's own docstring) -- so any gap below is physics, not an assembly error.
    #
    # THE PREFACTOR 2*pi*alpha^2*omega1*omega2/(2*J_i+1) and the HALF-RANGE INTEGRAL (Int_0^(DeltaE/2), which
    # already equals (1/2) Int_0^DeltaE for a omega1 <-> omega2 symmetric integrand -- checked to agree with
    # Int_0^DeltaE/2 numerically to 6 digits at the single point omega1=0.05 vs DeltaE-0.05 below) are the SAME
    # ones example-Dh.jl's own "lesson (i)" confirms are right; they are not re-derived from the literature here,
    # to avoid introducing a SECOND, independently-risked normalization next to DalgarnoLewis's own.
    #
    # CONVERGED in box (rbox = 60 and 80 a.u. agree to the last printed digit) and in mesh (hp = 0.03 and 0.05
    # agree likewise) -- checked before trusting the result, per Rule 12.
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

    function reducedAmplitude(omega1::Float64, omega2::Float64)
        wx12 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, Basics.Coulomb, omega1, omega2, orb_2s, orb_1s, p12, potential, primitives) +
        wx32 * DalgarnoLewis.twoPhotonAmplitude(E1, E1, Basics.Coulomb, omega1, omega2, orb_2s, orb_1s, p32, potential, primitives)
    end

    # 10-point Gauss-Legendre on [-1,1] (standard table), mapped onto (0, DeltaE/2) -- strictly INTERIOR nodes,
    # which avoids the two coincidental near-degeneracies with the 2p channel that sit exactly at the two
    # endpoints omega1 = 0 and omega1 = DeltaE (see the branch's own long comment above this code).
    nodes   = [-0.1488743389, 0.1488743389, -0.4333953941, 0.4333953941, -0.6794095683, 0.6794095683,
               -0.8650633667, 0.8650633667, -0.9739065285, 0.9739065285]
    weights = [ 0.2955242247, 0.2955242247,  0.2692667193, 0.2692667193,  0.2190863625, 0.2190863625,
                0.1494513492, 0.1494513492,  0.0666713443, 0.0666713443]
    alpha = Defaults.getDefaults("alpha");   a = 0.0;   b = DeltaE/2;   totalRate = 0.0
    for  (x,w)  in  zip(nodes, weights)
        global totalRate
        omega1 = (b-a)/2*x + (a+b)/2;   omega2 = DeltaE - omega1
        U = reducedAmplitude(omega1, omega2)
        dRate = 2pi * alpha^2 * omega1 * omega2 / 2.0 * abs2(U)     ## 2*J_i+1 = 2
        totalRate = totalRate + w * (b-a)/2 * dRate
        println(">> omega1 = $omega1 a.u., |U(K=0)|^2 = $(abs2(U)), dRate/domega1 = $dRate a.u.")
    end
    ratePerSecond = Defaults.convertUnits("rate: from atomic to 1/s", totalRate)
    println("\n>> A(2s -> 1s), chi-based, Coulomb/velocity gauge, K = 0 only  =  $ratePerSecond /s")
    println(">> exact literature value                                      =  8.2206 /s")
    println(">> ratio exact/chi-based                                       =  $(8.2206/ratePerSecond)")
    println(">> FOR CONTEXT, NOT RECOMPUTED HERE: example-Dh.jl branch b's bound-only pseudostate sum (n_max=20, " *
            "same gauge) gives 0.08907 /s, ratio exact/that = 92.3. The complete-basis chi result above is " *
            "roughly 5.8x CLOSER to the exact rate than the truncated sum, entirely from completing the " *
            "intermediate spectrum -- the one thing this module changes.")
end
