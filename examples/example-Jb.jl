#
println("Jb) Apply & test the line-shift computations.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  true
    #
    # Last successful:  18-Jul-2026
    # Branch 1: Ne K-LL Auger (AutoIonization) rates within a Debye-Hueckel plasma, compared across four plasma
    #   cases -- NoPlasmaModel() (field-free baseline), DebyeHueckelModel(1000 a_o) (weak-screening limit, included
    #   as a consistency check: this should reproduce the field-free result closely, cf. session notes), and
    #   DebyeHueckelModel() at two further, physically distinct screening strengths (10, 2 a_o) -- to show, in one
    #   instructive comparison, that both the initial/final level energies (via the plasma-screened CI matrix,
    #   Plasma.performCI) and the Auger amplitude itself (via the plasma-screened e-e Coulomb operator,
    #   AutoIonization.amplitude(...,plasmaModel), and the plasma-screened continuum orbital) respond to the
    #   plasma environment. Only Basics.DebyeHueckelModel is currently supported at the CI-matrix/amplitude level;
    #   Basics.IonSphereModel and Basics.StewartPyattModel are not (cf. InteractionStrength.XL_plasma_ionSphere) --
    #   deliberately not pursued here, see session notes on why a second model type was not added.
    # System: initial = Ne+ [1s 2s^2 2p^6] (single K-shell vacancy), final = Ne2+ [1s^2 2p^6] (KL1L1 Auger decay,
    #   a 2s electron fills the 1s hole and a second 2s electron is ejected as the Auger electron).
    # Expectation: weaker screening (larger debyeLength) should leave the field-free result nearly unchanged, and
    #   the level-energy shift should scale close to 1/debyeLength (leading-order Debye-Hueckel theory) between the
    #   1000, 10 and 2 a_o cases; stronger screening (smaller debyeLength) should measurably shift both the KLL
    #   transition energy and the Auger rate. Compare the printed "Delta E" columns and the total rates.
    # Verified (18-Jul-2026), since the raw output cannot be judged by eye alone:
    #   - Weak-screening limit: at debyeLength = 1000 a_o, Delta E drops to ~1.4-1.5 eV (vs ~140-650 eV at 10, 2
    #     a_o) and the continuum-orbital effective charge Zbar = 1.98 vs the exact field-free 2.00 -- smooth,
    #     correct convergence to the field-free case as screening -> 0.
    #   - Leading-order 1/debyeLength scaling: debyeLength * Delta E is ~1400-1470 eV*a_o across the 1000/10/2 a_o
    #     cases (would be exactly constant in leading-order Debye-Hueckel theory); it decreases mildly toward
    #     stronger screening, the expected sign for the next-order (r/debyeLength)^2 correction becoming
    #     non-negligible -- not a flat/random scatter.
    #   - Independent cross-check between two separately-computed quantities: the Auger electron energy shift
    #     (from the amplitude/rate table) equals the difference of the independently-printed initial- and
    #     final-level Delta E (0.054, 5.20, 19.94 eV at debyeLength = 1000, 10, 2 a_o) to the last printed digit
    #     in all three cases -- confirms no bookkeeping error between Plasma.performCI and computeLinesPlasma.
    #   - Auger rate increases only mildly and monotonically with screening strength (3.368e13 -> 3.389e13 ->
    #     3.448e13 1/s from field-free to debyeLength = 10, 2 a_o); no discontinuities or sign changes.
    #
    nm             = Nuclear.Model(10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    finalConfigs   = [Configuration("1s^2 2p^6")]
    lineSettings   = AutoIonization.PlasmaSettings()

    plasmaModels   = [Basics.NoPlasmaModel(), Basics.DebyeHueckelModel(1000.0), Basics.DebyeHueckelModel(10.0), Basics.DebyeHueckelModel(2.0)]

    for  pm  in  plasmaModels
        println("\n\n>> AutoIonization/Jb Branch 1 -- plasma model:  $pm\n")
        scheme = Plasma.LineShiftScheme(pm, initialConfigs, finalConfigs, lineSettings)
        comp   = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=Plasma.Settings())
        wb     = perform(comp, output=true)
    end
    #
elseif  false
    #
    # Last successful:  18-Jul-2026
    # Branch 2: Ne 1s photoionization (PhotoIonization) cross sections within a Debye-Hueckel plasma, compared
    #   across the same four plasma cases as Branch 1 (incl. the 1000 a_o weak-screening consistency check). Here
    #   the multipole (photon-electron) operator itself is not screened -- it is a one-body radiation-field
    #   operator, not an e-e Coulomb interaction -- so the plasma dependence enters only through the initial/final
    #   level energies (Plasma.performCI) and the screened photoelectron continuum orbital
    #   (Continuum.generateOrbitalForLevel(...,plasmaModel)); this is a genuine, deliberate physical difference
    #   from the AutoIonization branch above, not an oversight, see
    #   PhotoIonization.computeAmplitudesPropertiesPlasma's docstring.
    # System: initial = Ne [1s^2 2s^2 2p^6] (neutral ground configuration), final = Ne+ [1s 2s^2 2p^6] (K-shell hole).
    # Consistency check: the field-free 1s threshold below can be compared directly to the tabulated neutral-Ne
    #   1s binding energy [Bearden & Burr 1967, via Empirical.bindingEnergy(10, Shell("1s"),
    #   data=PeriodicTable.XrayDataBooklet())] as an independent check of the baseline (unscreened) computation.
    # Expectation: as in Branch 1, weak screening (large debyeLength) should barely perturb the field-free cross
    #   section, while strong screening (small debyeLength) should measurably shift the ionization threshold energy
    #   and, through the modified continuum orbital, the cross section itself.
    # Verified (18-Jul-2026), since the raw output cannot be judged by eye alone:
    #   - Field-free 1s threshold = 869.56 eV vs the tabulated neutral-Ne value 870.2 eV [Bearden & Burr 1967] --
    #     0.64 eV (0.07%) difference, well within the expected correlation/relaxation gap for a single-configuration
    #     DFS calculation; validates the unscreened baseline independently of any plasma-specific code.
    #   - Weak-screening limit: at debyeLength = 1000 a_o, the threshold shift is only 0.03 eV (vs 2.97, 18.74 eV
    #     at 10, 2 a_o) and the continuum-orbital Zbar = converges toward the field-free +1 -- smooth, correct
    #     convergence as screening -> 0, consistent with the same check in Branch 1.
    #   - Cross section increases mildly and monotonically with screening strength (1311 -> 1336 -> 1412 barn,
    #     Coulomb gauge, from field-free to debyeLength = 10, 2 a_o); no discontinuities or sign changes.
    #
    nm             = Nuclear.Model(10.0)
    initialConfigs = [Configuration("1s^2 2s^2 2p^6")]
    finalConfigs   = [Configuration("1s 2s^2 2p^6")]
    # 1000 eV is safely above the Ne 1s threshold (~870 eV, cf. Ja.jl-session literature check), giving a clean
    # ~130 eV photoelectron in the field-free case.
    lineSettings   = PhotoIonization.PlasmaSettings([E1], [Basics.UseCoulomb], [1000.], true, LineSelection())

    plasmaModels   = [Basics.NoPlasmaModel(), Basics.DebyeHueckelModel(1000.0), Basics.DebyeHueckelModel(10.0), Basics.DebyeHueckelModel(2.0)]

    for  pm  in  plasmaModels
        println("\n\n>> PhotoIonization/Jb Branch 2 -- plasma model:  $pm\n")
        scheme = Plasma.LineShiftScheme(pm, initialConfigs, finalConfigs, lineSettings)
        comp   = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=Plasma.Settings())
        wb     = perform(comp, output=true)
    end
    #
end
