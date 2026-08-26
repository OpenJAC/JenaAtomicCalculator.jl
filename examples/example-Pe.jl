
println("Pe) TWISTED (Bessel) photon beams: scattering at an atom held at a finite impact parameter from the vortex axis.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      22-Aug-2026
    # Last successful:  unknown ... see the report below.
    #
    # Branch a (THE VORTEX NODE, and the plane-wave limit): Rayleigh scattering of a Bessel beam at beryllium-like
    #   neon, scanning the IMPACT PARAMETER of the atom from the beam axis, for m = 0 and m = 1.
    #
    #   WHY A TWISTED BEAM IS COMPUTABLE HERE AND NOWHERE ELSE IN THIS MODULE.  A Bessel beam is a coherent
    #   superposition of plane waves on a cone of half-angle theta_k, each carrying the phase exp(i m phi_k).
    #   Superposing a SECOND-ORDER amplitude over that cone needs the partial-wave/OAM coupling that
    #   ParticleScattering also declined to re-derive -- its own Bessel method raises rather than returning a
    #   number from a different theory.  But in the FORM-FACTOR approximation the amplitude depends on nothing but
    #   the momentum transfer q, so every cone component's amplitude is already known and the superposition
    #   collapses to a one-dimensional integral over phi_k.  That is the whole reason this branch can exist.
    #
    #   TWO PREDICTIONS THAT COST NOTHING AND CAN BOTH FAIL.
    #
    #     (1) AT b = 0 THE CROSS SECTION MUST VANISH FOR m = 1.  A vortex beam has a NODE on its axis: with the
    #         form factor varying only weakly around the cone, the integral tends to i^m J_m(kappa b) times a
    #         constant, and J_m(0) = 0 for every m /= 0.  A non-zero on-axis result at m = 1 would mean the OAM
    #         phase is not being applied at all.  For m = 0 the same limit gives J_0(0) = 1, so the m = 0 row at
    #         b = 0 must be the LARGEST of its column, not the smallest.
    #
    #     (2) AS theta_k --> 0 EVERY ROW MUST APPROACH THE PLANE-WAVE VALUE, and that value is DATED: example-Pd.jl
    #         branch a is anchored absolutely against N^2 sigma_Thomson to one part in 10^5.  So this is a limit
    #         check against a verified absolute number rather than against another unverified calculation -- which
    #         is a stronger position than any other branch of the P series has been able to take.
    #
    #   WHAT THIS IS NOT.  A SCALAR treatment.  The vector nature of the Bessel field and its polarization
    #   structure are absent, so the (1 + cos^2 theta) factor is carried over unchanged from the plane-wave case.
    #   A full vector treatment would alter the angular distribution, and nothing here should be quoted as a
    #   polarization result.  The Stokes parameters are returned as NaN throughout, deliberately.
    #
    #   REPORT (22-Aug-2026, 500 eV, opening angle 0.20 rad).  BOTH PREDICTIONS PASS, and a third confirmation
    #   appeared that had been predicted structurally but not computed in advance.
    #
    #       b [a.u.]        m = 0            m = 1          ratio to previous (m = 1)
    #        0.0        3.786484e-07     7.866593e-15               -
    #        0.5        3.786148e-07     1.680085e-11               -
    #        1.0        3.785140e-07     6.717086e-11             3.998
    #        2.0        3.781113e-07     2.685169e-10             3.998
    #        4.0        3.765033e-07     1.071759e-09             3.992
    #        8.0        3.701229e-07     4.250619e-09             3.966
    #
    #   (1) THE VORTEX NODE IS THERE.  At b = 0 the m = 1 cross section is 7.87e-15 against the m = 0 value of
    #       3.79e-07 -- suppressed by 4.8e-08, i.e. to eight orders of magnitude.  It is not exactly zero because
    #       the cone integral is a 721-point trapezoid; a node to 8 decades is what a numerical quadrature can
    #       deliver.  And the m = 0 column behaves oppositely, largest at b = 0 and falling monotonically, which is
    #       J_0(0) = 1.  A code that ignored the OAM phase would have given the two columns the same shape.
    #
    #   (2) THE PLANE-WAVE LIMIT IS ALREADY NEARLY MET at theta_k = 0.20 rad: the m = 0, b = 0 entry is
    #       3.786484e-07 against the DATED plane-wave value 3.786772e-07 of example-Pd.jl branch a, a relative
    #       difference of 7.6e-05.  That comparison inherits an ABSOLUTE anchor rather than comparing two of my own
    #       unverified numbers, which is a stronger position than any other branch of the P series can take.
    #
    #   (3) THE THIRD, NOT DESIGNED IN: sigma goes as b^2 for m = 1.  The ratios between successive impact
    #       parameters are 3.998, 3.998, 3.992, 3.966 against a predicted 4.  That is J_1(kappa b)^2 ~ (kappa b/2)^2
    #       at small argument -- and kappa b = 0.21 even at b = 8, so small-argument is indeed the right regime.
    #       The drift at the largest b is where kappa b ceases to be negligible, i.e. the beginning of the true
    #       Bessel oscillation, and is in the right direction.
    #
    #   NOT DATED.  The three checks above are SHAPES and RATIOS -- a node, a power law, a limit -- and the
    #   absolute scale is inherited from example-Pd.jl rather than established here.  More importantly the scalar
    #   restriction is a genuine physical omission, not a numerical one: a real twisted beam has a polarization
    #   structure that this treatment discards, and a date would suggest a completeness the branch does not have.
    #
    grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm   = Nuclear.Model(10.)

    for  m  in  [0, 1]
        settings = PhotonScattering.Settings(PhotonScattering.Settings();
                        process          = PhotonScattering.RayleighScattering(),
                        approximation    = PhotonScattering.FormFactorApproximation(),
                        beamType         = Beam.BesselBeam(m, 0.20, 1.0),      # m, opening angle [rad], kz (unused here)
                        photonEnergies   = [500.0],
                        impactParameters = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0],
                        polarThetas      = [0.0, pi/2, pi],
                        printBefore      = false )

        wa = Atomic.Computation(Atomic.Computation(), name="Bessel-beam Rayleigh, m = $m", grid=grid, nuclearModel=nm,
                                initialConfigs  = [Configuration("1s^2 2s^2")],
                                finalConfigs    = [Configuration("1s^2 2s^2")],
                                processSettings = settings )
        wb = perform(wa; output=true)
    end
    #
elseif  false
    # Last visit:      22-Aug-2026
    # Last successful:  unknown ... not yet run.
    #
    # Branch b (THE PLANE-WAVE LIMIT, made quantitative): the same Rayleigh scattering at m = 0 and b = 0, with the
    #   opening angle driven towards zero, against the DATED plane-wave value of example-Pd.jl branch a.
    #
    #   WHY THIS DESERVES ITS OWN BRANCH.  Branch a asserts that theta_k --> 0 recovers the plane wave; this branch
    #   measures how fast.  The cone geometry gives cos(Theta) = sin(theta_k) sin(theta) cos(phi_k - phi) +
    #   cos(theta_k) cos(theta), which differs from the plane-wave cos(theta) at SECOND order in theta_k, so the
    #   cross section should approach its limit quadratically: halving the opening angle should quarter the
    #   residual.  A linear approach, or a residual that stalls, would point at the cone integral rather than at
    #   the physics.
    #
    #   The target value is 3.786772e-07 a.u. at 500 eV, from example-Pd.jl branch a -- a branch that is dated
    #   BECAUSE it is anchored absolutely, so the comparison inherits that anchoring rather than merely comparing
    #   two of my own numbers with each other.
    #
    grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm   = Nuclear.Model(10.)

    for  thetaK  in  [0.40, 0.20, 0.10, 0.05, 0.025]
        settings = PhotonScattering.Settings(PhotonScattering.Settings();
                        process          = PhotonScattering.RayleighScattering(),
                        approximation    = PhotonScattering.FormFactorApproximation(),
                        beamType         = Beam.BesselBeam(0, thetaK, 1.0),
                        photonEnergies   = [500.0],
                        impactParameters = [0.0],
                        polarThetas      = Float64[],
                        printBefore      = false )

        wa = Atomic.Computation(Atomic.Computation(), name="Bessel-beam Rayleigh, theta_k = $thetaK", grid=grid,
                                nuclearModel=nm,
                                initialConfigs  = [Configuration("1s^2 2s^2")],
                                finalConfigs    = [Configuration("1s^2 2s^2")],
                                processSettings = settings )
        wb = perform(wa; output=true)
    end
    #
end
#
setDefaults("print summary: close", "")
