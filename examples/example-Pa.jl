
println("Pa) BOUND-FREE PAIR CREATION: gamma + |i(N)> --> |f(N+1)> + e^+, the created electron being captured into a bound orbital.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... NOT VERIFIED; the gauges disagree, see the report below.
    #
    # Branch a (first computation): hydrogen-like uranium capturing the created electron into its own 1s,
    #   gamma + U^91+ (1s ^2S_1/2) --> U^90+ (1s^2 ^1S_0) + e^+, at a photon energy of 40000 a.u. = 1088.5 keV.
    #
    #   WHY THE PROCESS EXISTS, AND WHY IT IS NOT EXOTIC.  A single photon cannot create a pair in vacuum -- momentum
    #   cannot balance -- so free pair creation needs a nucleus to take the recoil.  Here the created ELECTRON does not
    #   merely feel the nucleus, it is captured by it, into a bound orbital of the very ion that provided the field.
    #   At the LHC this is the dominant beam-loss mechanism for fully stripped heavy ions.
    #
    #   THE SIGNATURE IS THE LOWERED THRESHOLD.  Energy conservation gives T_+ = omega - 2 m c^2 + B, with B the binding
    #   energy gained by the captured electron, so the channel opens at omega > 2 m c^2 - B rather than at 1.022 MeV.
    #
    #   WHAT PASSES.
    #     * ENERGY CONSERVATION AND THE THRESHOLD SHIFT, computed rather than put in.  At omega = 1088.455 keV the
    #       positron comes out at 196.62 keV, i.e. B = 196.62 - (1088.455 - 1022.0) = 130.17 keV.  That is the same
    #       1s binding energy of helium-like uranium that the ANNIHILATION module obtained independently the day before
    #       (example-Of.jl branch a), against ~131.8 keV known.  So this channel opens at 891.8 keV, appreciably below
    #       the free-pair threshold, which is exactly the qualitative claim above.
    #     * The magnitude is not obviously wrong, which is more than could be said for the annihilation case.
    #       1.196e-11 a.u. is 3.4e-28 cm^2 = 0.34 mb, and bound-free pair creation by photons on a high-Z ion near
    #       threshold is genuinely of order microbarns to millibarns.
    #
    #   WHAT FAILS.
    #     * GAUGE DISAGREEMENT BY A FACTOR OF 8.5: 1.196150e-11 (Coulomb) against 1.012764e-10 (Babushkin) a.u.  The
    #       annihilation module shows the same disease with a factor of 17.5, and the two share their charge-conjugate
    #       orbital construction and their photon operator, so it is PLAUSIBLE that this is one defect seen twice rather
    #       than two -- but that is a conjecture from shared code, not a measurement, and branch b does not settle it.
    #
    #   A SILENT BUG, fixed, and the SECOND of its kind: the two-step angular construction runs the OPPOSITE way round
    #   from the annihilation case, because the charge-conjugate positron orbital is attached to the INITIAL level here
    #   rather than to the final one.  The intermediate symmetry must therefore come from symf through the multipole and
    #   kappa must couple symi to it.  Carrying over the annihilation ordering raises nothing at all: no CSF survives
    #   the coupling and every amplitude is exactly zero, with a well-formed table of zeros to show for it.  Together
    #   with the kappa-sign error in example-Of.jl this is a pattern worth naming -- an angular-coupling mistake in this
    #   construction is ALWAYS silent, and a table of exact zeros is its signature.
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-7, h=5.0e-3, hp=1.0e-3, rbox=3.0)
    nm   = Nuclear.Model(92.)

    settings = PhotonScattering.Settings(PhotonScattering.Settings();
                    process        = PhotonScattering.BoundFreePairCreation(),
                    approximation  = PhotonScattering.FirstOrderVertex(),
                    photonEnergies = [40000.0],                    # in a.u.; 2 m c^2 = 37558 a.u.
                    multipoles     = [E1, M1, E2],
                    maxKappa       = 4,
                    printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(); name="bound-free pair creation, H-like U",
                            grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s")],
                            finalConfigs    = [Configuration("1s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... the test RAN and gave a clean, quantitative answer, but that answer is a FAILURE,
    #                              so the branch is not dated.  Read the report -- the failure is the useful part.
    #
    # Branch b (DETAILED BALANCE against the annihilation module): the same computation as branch a, but at the photon
    #   energy of example-Of.jl branch a, so that the two runs are exact time-reverses of one another,
    #
    #       e^+ + U^90+ (1s^2)  -->  U^91+ (1s) + gamma        (ParticleScattering, one-photon annihilation)
    #       gamma + U^91+ (1s)  -->  U^90+ (1s^2) + e^+        (this module, bound-free pair creation)
    #
    #   WHY THIS TEST IS WORTH MORE THAN A LITERATURE COMPARISON.  Both cross sections come from JAC, and both carry a
    #   prefactor that was written down rather than derived.  Detailed balance relates them EXACTLY,
    #
    #       g_a p_a^2 sigma_(a->b)  =  g_b p_b^2 sigma_(b->a) ,
    #
    #   so the test needs no measured number and no published calculation.  It cannot say whether either prefactor is
    #   right, but it can say whether they are MUTUALLY consistent -- and if they are not, the violation is a property
    #   of the pair, which is a far smaller thing to search than two independent absolute normalizations.
    #
    #   THE NUMBERS.  With g_a = 2 x 1 (positron spin, J=0 of 1s^2) and p_a = p_+ = 14.1610 a.u. at T_+ = 100 a.u.,
    #   against g_b = 2 x 2 (photon polarizations, J=1/2 of 1s) and p_b = omega/c = 239.896 a.u.:
    #
    #       required        sigma_ann / sigma_pc  =    573.97
    #       observed        Coulomb                =  3.1327e+05        violation factor  545.8
    #                       Babushkin              =  3.2065e+05        violation factor  558.6
    #
    #   A FIRST READING OF THIS AS A CONSTANT 4c = 548.14 WAS WRONG, and is recorded here as withdrawn rather than quietly
    #   removed, because the mistake is instructive: it rested on ONE energy point, where "constant" cannot be measured at
    #   all.  The test was therefore repeated at a SECOND positron energy, T_+ = 400 a.u. against 100, i.e. p_+ doubled
    #   (14.1610 -> 28.4345, ratio 2.008), with omega following from the annihilation run at 902.7196 keV:
    #
    #                        p_+        required      violation (Coulomb)   violation (Babushkin)
    #       T_+ = 100     14.1610        573.97              545.8                 558.6
    #       T_+ = 400     28.4345        144.97              448.3                 462.0
    #
    #   THIS RULES OUT BOTH SIMPLE EXPLANATIONS, which is the whole value of the second point.  A pure PREFACTOR error
    #   would scale as 1/p_+ and the violation would have fallen by 0.498, to 272 / 278.  A CONSTANT error would have left
    #   it at 546 / 559.  It fell by 0.821 (Coulomb) and 0.827 (Babushkin) -- neither, and not a clean power either, the
    #   implied exponent being about p_+^(-0.28).
    #
    #   WHAT THE TWO POINTS TOGETHER DO SAY, and it is a better lead than the withdrawn one.  The energy dependence of the
    #   violation is the SAME in both gauges to within half a percent (0.821 against 0.827), so whatever carries it is
    #   GAUGE-INDEPENDENT -- and therefore is NOT the gauge disagreement that branch a reports, which stays near constant
    #   in energy (17.5 -> 17.2 in the annihilation module across the same factor of four).  So there are TWO faults, not
    #   one, and they are separable by this signature:
    #
    #       (i)  a gauge-INDEPENDENT, energy-DEPENDENT inconsistency between the two directions.  A continuum
    #            normalization convention is the natural suspect -- per unit energy against per unit momentum differ by
    #            exactly a factor that varies with p_+ -- and it would sit in how each module uses the orbital returned by
    #            Continuum.generateOrbitalLocalPotential rather than in the prefactor expressions themselves.
    #       (ii) a gauge-DEPENDENT, energy-INDEPENDENT defect in the amplitude, common to both modules, which the multipole
    #            convergence test of example-Of.jl already showed is not a truncation.
    #
    #   Suspect (i) first: it is the one the two energy points actually measure, and fixing it changes the target that
    #   suspect (ii) must then be measured against.
    #
    #   Note the small mismatch to be tidied when this is next run: omega is set to 32874.62 a.u. against the 32874.334
    #   that the annihilation run's photon energy of 894.5562 keV corresponds to, which puts T_+ at 100.28 a.u. instead
    #   of 100.  It is a 0.28 % effect and does not touch a factor of 548, but it should be made exact before the branch
    #   is ever dated.
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-7, h=5.0e-3, hp=1.0e-3, rbox=3.0)
    nm   = Nuclear.Model(92.)

    settings = PhotonScattering.Settings(PhotonScattering.Settings();
                    process        = PhotonScattering.BoundFreePairCreation(),
                    approximation  = PhotonScattering.FirstOrderVertex(),
                    photonEnergies = [32874.334],                  # = 894.5562 keV, the annihilation run's photon energy
                    multipoles     = [E1, M1, E2],
                    maxKappa       = 5,                            # matches |kappa| <= 5 of example-Of.jl's maxL = 4
                    printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(); name="pair creation at the annihilation energy",
                            grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s")],
                            finalConfigs    = [Configuration("1s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
end
