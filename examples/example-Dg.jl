
println("Dg) Apply & test the RayleighCompton module with ASF from an internally generated initial- and final-state multiplets.")

setDefaults("print summary: open", "zzz-RayleighCompton.sum")
setDefaults("unit: energy", "eV")

if  false
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... BLOCKED by a defect in src/module-RayleighCompton.jl itself, see below.
    #
    # Branch a: Rayleigh scattering on the ground level of B-like neon, 1s^2 2s^2 2p, at two photon energies.
    #   The branch is SELF-CONTAINED: it generates the Green function of intermediate levels and then uses it,
    #   in that order, in one branch.
    #
    #   THIS BRANCH IS DISABLED AND CANNOT BE MADE TO RUN FROM THE EXAMPLE SIDE.  Six faults in the example were
    #   repaired on 21-Aug-2026 (listed below) and the branch then reached RayleighCompton.computeChannelAmplitude,
    #   where it stops on a defect in the MODULE:
    #
    #       module-RayleighCompton.jl:336-338
    #           amplitude = 0.0im;   lowerNom = lowerDenom = upperNom = upperDenom = 0.
    #           for  (ig, gLevel)  in  enumerate(gChannel.gMultiplet.levels)
    #               if   ig != 1   continue   end          <-- every level except the FIRST is skipped
    #
    #   The comment above that loop reads "Sum over all terms of the Green function channel", and the body runs
    #   for ig == 1 only.  THE SUM OVER INTERMEDIATE STATES IS NOT PERFORMED: however many levels the Green
    #   expansion holds, exactly one term contributes, which is a second-order amplitude reduced to its first
    #   term.  Note the contrast with the pole-DETECTION loop just above it, which correctly reads
    #   `if ig == 1 continue end` -- skipping only the first element so it can compare with ig-1.  Line 338 is
    #   that line with the test inverted, and looks like a debugging leftover.
    #
    #   The MethodError one sees is a consequence, not the disease: all four pole variables are initialised as
    #   Float64 on line 336 and can only be reassigned inside the ig == leftIdx / leftIdx+1 branches, which are
    #   unreachable for any ig != 1.  So upperNom stays Float64 and the call at line 362 fails against the inner
    #   function's upperNom::ComplexF64.  Without a pole the module returns a silent one-term sum; with a pole it
    #   throws.  Its "Last successful: 13May2024" dates cannot have covered this path.
    #
    #   DECISION (maintainer, 21-Aug-2026): the module is NOT repaired.  It was only ever run in a very limited
    #   sense, and getting it right is a new implementation rather than a fix.  Rayleigh scattering will be
    #   implemented independently inside PhotonScattering, with a gMultiplet of just 2-4 intermediate levels,
    #   following the pattern of MultiPhotonTransition's two-photon amplitudes -- which do the same sum, are
    #   validated, and already accept Union{Multiplet, Array{GreenChannel,1}} for the intermediate states.
    #   This branch is kept as the record of what was found and as the input that the new implementation should
    #   reproduce.
    #
    #   WHAT WAS BROKEN, recorded because the file carried a "Last successful: 13May2024" date throughout and
    #   none of it could have run since.  Five faults, each hidden behind the one before:
    #     (1) The active branch used `green`, which was produced only by the FIRST branch -- and that branch was
    #         `if false`.  Any run died on an undefined variable before reaching the physics.
    #     (2) RayleighCompton.Settings(...) was called positionally with TEN arguments against a struct with
    #         ELEVEN fields; `printBefore` was missing.  Now built with the keyword copy-constructor, which is
    #         the JAC convention and cannot go stale by argument count again.
    #     (3) JAC.UseCoulomb / JAC.UseBabushkin -- the module name from before the rename to
    #         JenaAtomicCalculator.
    #     (4) The Green function was generated on Radial.Grid(true) while the scattering ran on a DIFFERENT,
    #         explicitly parametrized grid, so the intermediate-state orbitals and the scattering orbitals did
    #         not share a radial mesh.  Both now use `grid`.
    #     (5) The computation was named "xx" and the Green representation called the target "Be-like", where
    #         1s^2 2s^2 2p is five electrons and therefore B-like.
    #
    #   The energy unit is set to eV explicitly above; it was previously left at whatever the session default
    #   happened to be, while `setDefaults("unit: rate", ...)` was set instead and is not used by this module.
    #
    #   WHAT TO CHECK IN THE OUTPUT, once it runs:
    #     * the two gauges, Coulomb and Babushkin, should agree -- they are printed side by side;
    #     * the elastic cross section should FALL between 10 eV and 100 eV, away from resonance;
    #     * the angle-differential cross section should be largest in the forward direction.
    #   The date stays blank until those are seen and judged, per Rule 7.
    #
    grid             = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    refConfigs       = [Configuration("1s^2 2s^2 2p")]
    # The Green expansion must cover EVERY intermediate symmetry the channels ask for, or
    # RayleighCompton.computeChannelAmplitude dies with "stop a: Green channel not found".  Two E1 photons reach
    # 1/2+, 3/2+, 5/2+ from the 2p_1/2 and 2p_3/2 initial levels, and two M1 photons reach 1/2-, 3/2-, 5/2-.
    # The file previously asked for [1/2+, 3/2-] alone, i.e. a Green function that did not contain the
    # symmetries its own channels required -- fault (6).
    levelSymmetries  = [LevelSymmetry(1//2, Basics.plus),  LevelSymmetry(3//2, Basics.plus),
                        LevelSymmetry(5//2, Basics.plus),  LevelSymmetry(1//2, Basics.minus),
                        LevelSymmetry(3//2, Basics.minus), LevelSymmetry(5//2, Basics.minus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation("Green function for Rayleigh scattering on B-like neon",
                                      Nuclear.Model(10.), grid, refConfigs,
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(),
                                                      levelSymmetries, 3, greenSettings) )
    greenOut         = generate(greenRep, output=true)
    green            = greenOut["Green channels"]

    rayleighSettings = RayleighCompton.Settings(RayleighCompton.Settings();
                            multipoles = [E1, M1], gauges = [Basics.UseCoulomb, Basics.UseBabushkin],
                            photonEnergies = [10., 100.], green = green,
                            calcRayleighRaman = true, calcAngular = true, calcStokes = true, printBefore = true,
                            incidentStokes = ExpStokes(), solidAngles = [SolidAngle(1.0, 0.0)],
                            lineSelection = LineSelection() )

    wa = Atomic.Computation(Atomic.Computation(), name="Rayleigh scattering on B-like neon", grid=grid,
                            nuclearModel=Nuclear.Model(10.),
                            initialConfigs  = [Configuration("1s^2 2s^2 2p")],
                            finalConfigs    = [Configuration("1s^2 2s^2 2p")],
                            processSettings = rayleighSettings )
    wb = perform(wa)
    #
elseif  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... this branch RUNS but computes no scattering; it is the technical half only.
    #
    # Branch b: the Green function of intermediate levels ALONE.  This is the ACTIVE branch, and deliberately so:
    #   branch a is blocked by a module defect that will not be repaired (see there), and a file in which no
    #   branch runs at all is worse than one whose running branch is honest about being partial.  This branch
    #   exercises the whole path up to but not including RayleighCompton -- the grid, the reference
    #   configuration, the six intermediate symmetries and the Green expansion itself -- so a later breakage in
    #   any of those is still caught here.
    #
    #   It is NOT dated, because generating a Green expansion is not a physics result: nothing in its output can
    #   be checked against a known number.  It is a technical check, and the date line says so rather than
    #   implying a verification that did not happen.
    #
    grid             = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    refConfigs       = [Configuration("1s^2 2s^2 2p")]
    # The Green expansion must cover EVERY intermediate symmetry the channels ask for, or
    # RayleighCompton.computeChannelAmplitude dies with "stop a: Green channel not found".  Two E1 photons reach
    # 1/2+, 3/2+, 5/2+ from the 2p_1/2 and 2p_3/2 initial levels, and two M1 photons reach 1/2-, 3/2-, 5/2-.
    # The file previously asked for [1/2+, 3/2-] alone, i.e. a Green function that did not contain the
    # symmetries its own channels required -- fault (6).
    levelSymmetries  = [LevelSymmetry(1//2, Basics.plus),  LevelSymmetry(3//2, Basics.plus),
                        LevelSymmetry(5//2, Basics.plus),  LevelSymmetry(1//2, Basics.minus),
                        LevelSymmetry(3//2, Basics.minus), LevelSymmetry(5//2, Basics.minus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())
    greenRep         = Representation("Green function for Rayleigh scattering on B-like neon",
                                      Nuclear.Model(10.), grid, refConfigs,
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(),
                                                      levelSymmetries, 3, greenSettings) )
    greenOut         = generate(greenRep, output=true)
    green            = greenOut["Green channels"]
    #
end
#
setDefaults("print summary: close", "")
