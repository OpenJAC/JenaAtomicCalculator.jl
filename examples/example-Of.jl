
println("Of) ONE-PHOTON ANNIHILATION of a projectile positron with a bound electron: e^+ + |i(N)> --> |f(N-1)> + photon.")

setDefaults("print summary: open", "zzz-ParticleScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... NOT VERIFIED; two checks fail, see the report below.
    #
    # Branch a (first computation): helium-like uranium, 1s^2 ^1S_0 --> 1s ^2S_1/2 of hydrogen-like uranium, with a
    #   positron of 100 a.u. = 2.72 keV.  Uranium is chosen deliberately and is NOT an arbitrary heavy example: the
    #   annihilation photon carries omega ~ 2 m c^2, i.e. q a_0 ~ 274, so the retardation factor exp(i q.r) oscillates
    #   about q <r> ~ 274/Z times across a K-shell orbital.  That is ~3 for uranium and ~137 for helium, so a LIGHT
    #   target is the worst case for a multipole expansion here, not the easiest -- the opposite of the intuition
    #   carried over from optical transitions.
    #
    #   WHY THE PROCESS EXISTS.  A free electron and positron cannot annihilate into ONE photon: in their common rest
    #   frame the photon would need zero momentum and 2 m c^2 of energy at once.  Binding the electron lifts that,
    #   the nucleus taking up the recoil, so the amplitude is controlled by the part of the electron's momentum
    #   distribution reaching p ~ m c -- which lives close to the nucleus.  Hence a K-shell, high-Z process, with
    #   sigma ~ Z^5 for K-shell annihilation, vanishing in the free limit.  It is the crossing partner of
    #   BOUND-FREE PAIR CREATION, which is computed by PhotonScattering.BoundFreePairCreation() -- see
    #   example-Pa.jl, whose branch a obtains the same 130.17 keV binding energy from the other side of the
    #   crossing.  (It used to point here at module-PairProduction.jl, which never held any computation and was
    #   retired on 25-Aug-2026.)
    #
    #   WHAT PASSES.
    #     * ENERGY CONSERVATION, which is the one thing in this process that has no adjustable content.  The photon
    #       comes out at omega = 894.556 keV against 2 m c^2 + T_+ = 1024.72 keV, i.e. an implied 1s binding energy
    #       of 130.17 keV for helium-like uranium, against ~131.8 keV known.  1.2 % low, which is what a Dirac-Fock
    #       calculation missing correlation should give, and it is computed from the level energies rather than put in.
    #     * The pipeline runs end to end through Atomic.Computation and perform(), the positron partial waves are
    #       generated in the positron potential (reversed sign, no exchange) and the channels are built and summed.
    #
    #   WHAT FAILS, and why this branch carries no date.
    #     * GAUGE DISAGREEMENT BY A FACTOR OF ABOUT 17.  With multipoles [E1, M1, E2] the cross section is
    #       8.152e-22 (Coulomb) against 1.427e-20 (Babushkin) a.u.  Two gauges of the same amplitude must agree to
    #       several digits for a converged calculation; a factor of 17 means the calculation is not converged, or the
    #       amplitude is wrong.  TRUNCATION OF THE MULTIPOLE SERIES HAS BEEN RULED OUT as the cause, which was the
    #       obvious first suspect since q <r> ~ 3 here: repeating the run with [E1, M1, E2, M2, E3, M3, E4] and
    #       maxL = 8 gives 8.305e-22 / 1.433e-20 against 8.152e-22 / 1.427e-20 for [E1, M1, E2] with maxL = 4, i.e.
    #       1.9 % and 0.4 %.  The series is converged and the factor of 17 survives it untouched, so the defect lies
    #       in the amplitude itself or in the charge-conjugation convention, not in how many multipoles are kept.
    #       Worth noting that the insensitivity is itself odd: at q <r> ~ 3 the higher multipoles ought to matter
    #       more than this, which may be the same defect seen from another side.
    #     * THE ABSOLUTE MAGNITUDE IS IMPLAUSIBLE.  8e-22 a.u. is 2e-38 cm^2, where the literature for K-shell
    #       one-photon annihilation on a heavy element at these energies is of order tens of millibarns, i.e. ~1e-26
    #       cm^2.  That is about twelve orders of magnitude, so this is a normalization defect rather than a small
    #       inaccuracy.  ParticleScattering.annihilationCrossSection carries the prefactor and its docstring already
    #       says it was written down rather than derived against a known case; that is where to look first.
    #       UPDATE 21-Aug-2026: the crossing partner now CONSTRAINS it, though it does not yet localize it.  Bound-free pair
    #       creation (example-Pa.jl branch b) is the time-reverse of this process, and detailed balance between the two
    #       requires sigma_ann / sigma_pc = 573.97 where JAC gives 3.1327e+05 (Coulomb) and 3.2065e+05 (Babushkin) --
    #       violation factors of 545.8 and 558.6.  A first reading of that as a constant 4c = 548.14 was WRONG and is
    #       withdrawn: it rested on one energy point, and comparing the two prefactors directly refutes it.  Were the
    #       amplitudes equal in both directions, the coded prefactors alone would give c^2/p_+ = 1326 here -- energy-
    #       dependent, going as 1/p_+, not 548.  Since 546-559 is 0.41 of that, the amplitudes differ between the two
    #       directions as well, and both the prefactor and the amplitude are implicated.  The decisive test is cheap and
    #       unrun: repeat at a second positron energy, a 1/p_+ scaling being the prefactor alone.
    #
    #   A BUG FOUND AND FIXED WHILE WRITING THIS, recorded because it is silent and would recur.  The kappa stored in
    #   an AnnihilationChannel is the kappa of the CHARGE-CONJUGATED orbital, i.e. of the one that enters the CSF and
    #   the angular coupling, whereas the positron partial wave that must be generated is the one with the OPPOSITE
    #   sign, charge conjugation reversing kappa.  Generating the wave with the channel's own kappa produced no error
    #   at all: the conjugated orbital then could not couple to the channel's total symmetry, no CSF survived, and
    #   every amplitude came out exactly zero, with a perfectly well-formed table of zeros to show for it.
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-7, h=5.0e-3, hp=1.0e-3, rbox=3.0)
    nm   = Nuclear.Model(92.)

    settings = ParticleScattering.Settings(ParticleScattering.Settings();
                    projectile     = ParticleScattering.Positron(),
                    process        = ParticleScattering.Annihilation1Photon(),
                    interaction    = ParticleScattering.StaticField(),
                    impactEnergies = [100.0],                      # in a.u., as everywhere in this module
                    multipoles     = [E1, M1, E2],
                    maxL           = 4,
                    printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(); name="1-photon annihilation, He-like U",
                            grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... not yet run; it is blocked on branch a, see below.
    #
    # Branch b (the meaningful validation): the Z-dependence of K-shell one-photon annihilation along the
    #   helium-like sequence, Xe (54), Yb (70), Au (79), U (92), at one fixed positron energy.
    #
    #   WHY THIS IS THE RIGHT FIRST TEST rather than an absolute comparison.  The known result is that the K-shell
    #   one-photon annihilation cross section rises as roughly Z^5 -- the same scaling as the photoelectric effect,
    #   and for the same reason, the amplitude sampling the electron density near the nucleus.  A RATIO of cross
    #   sections is completely insensitive to the overall prefactor, which is exactly the quantity in
    #   ParticleScattering.annihilationCrossSection that has not been derived.  So this branch can confirm or refute
    #   the physics of the amplitude while the normalization is still open, which no absolute number can.
    #
    #   Fitting log sigma against log Z over these four points should give an exponent near 5.  A markedly different
    #   exponent would point at the amplitude or at the charge-conjugation convention rather than at the prefactor.
    #
    #   BLOCKED ON BRANCH a, though no longer for the reason first written down.  The multipole series IS converged,
    #   so the numbers are not merely truncated -- but the surviving factor of 17 between the gauges says the amplitude
    #   or the charge-conjugation convention is wrong, and a Z-exponent fitted to a wrong amplitude measures that error
    #   rather than the physics.  Branch a first.
    #
    #   Note the box must be matched to each Z (Rule 12): a K-shell orbital has <r> ~ a_0/Z, so a box of 3 a.u. is
    #   already generous for uranium and far too small in nothing here -- but the multipole series length grows as
    #   274/Z falls, so the LIGHTER members of this scan are the expensive ones, not the heavier.
    #
    for  Z  in  [54., 70., 79., 92.]
        grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-7, h=5.0e-3, hp=1.0e-3, rbox=3.0)
        nm   = Nuclear.Model(Z)

        settings = ParticleScattering.Settings(ParticleScattering.Settings();
                        projectile     = ParticleScattering.Positron(),
                        process        = ParticleScattering.Annihilation1Photon(),
                        interaction    = ParticleScattering.StaticField(),
                        impactEnergies = [100.0],
                        multipoles     = [E1, M1, E2, M2, E3, M3, E4],
                        maxL           = 8 )

        wa = Atomic.Computation(Atomic.Computation(); name="1-photon annihilation, He-like Z=$Z",
                                grid=grid, nuclearModel=nm,
                                initialConfigs  = [Configuration("1s^2")],
                                finalConfigs    = [Configuration("1s")],
                                processSettings = settings )
        wb = perform(wa; output=true)
    end
    #
end
