#
println("Jc) Apply & test the Saha-Boltzmann computations for an ionic mixture in LTE.")

setDefaults("print summary: open", "zzz-saha-boltzmann.sum")


if  false
    # Last successful:  18-Jul-2026
    # Compute the Saha-Boltzmann equilibrium densities for a mixture of Carbon ions, restricted to the highest
    # charge states q = 4:6 (bare through He-like), NoExcitations = 1, upperShellNo = 4 -- a light, fast case to
    # exercise the SahaBoltzmannScheme machinery end-to-end.
    # Bugs found and fixed to get this far (in src/module-Plasma-inc-saha-boltzmann-mixture.jl, not this file):
    #   - readEvaluateIonLevelData(filename,...) called JLD2.load(filename) unconditionally, with no isfile check,
    #     contradicting its own docstring ("reads in, if available"); a missing/first-run file crashed instead of
    #     falling back to generation. Fixed with an isfile guard.
    #   - generateIonLevelData called Basics.performCI(...) -- an empty dispatch stub with zero methods -- where
    #     Hamiltonian.performCIwithFrozenOrbitals(...) (matching signature, intended for exactly this per-config,
    #     frozen-mean-field-orbital use) was needed. This is what actually generates each charge state's ionic
    #     levels, so every SahaBoltzmannScheme run needing fresh ion-level data was blocked by this before now.
    # Verified (18-Jul-2026): clean run, no errors/warnings beyond the pre-existing (self-flagged, unrelated)
    #   "Check the grid in Plasma: re-install !!" TODO. Generated 1 level for q=6 (bare nucleus, trivial), 16
    #   levels for q=5 (H-like, n<=4 single excitations), 31 levels for q=4 (He-like). LTE result at T=200 eV,
    #   n_i=2e-5 a.u. (very low density): charge-state distribution q=6:5:4 = 99.53% : 0.466% : 0.00038% within
    #   the allowed q=4:6 range -- directionally sensible (low density favours high ionization at fixed T via the
    #   Saha equation), though qRange=4:6 excludes lower charge states by construction, so this alone does not
    #   confirm the *overall* ionization balance would peak here if wider q were allowed. Independent check: the
    #   code's own built-in pressure cross-check (computed once from the mean charge state, once by explicitly
    #   summing level populations, cf. Plasma.pressure's ">>>> Pressure = ... =!= wp = ..." printout) matched
    #   exactly (0.0010282945050021789 both ways).
    grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.011, 1.0)]  ## , IsotopicFraction(9., 20.2, 0.4)
    scheme      = Plasma.SahaBoltzmannScheme(Basics.NoPlasmaModel(), true, true, 4:6, 10000, 1, 4, ionMixture, String[])

    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=settings)
    wb          = perform(wa, output=true)

elseif  false
    # Last successful:  18-Jul-2026
    # Compute the Saha-Boltzmann equilibrium densities for a mixture of Carbon ions including all charge states
    # q = 2:6 (Be-like through bare), NoExcitations = 2, upperShellNo = 4 (simplified from the original
    # upperShellNo = 8, which was untested and considerably more expensive).
    # Verified (18-Jul-2026): clean run. Generated 1/16/31/67/351 ion-levels for q = 6/5/4/3/2 respectively.
    #   At the same T = 200 eV, n_i = 2e-5 a.u. as Branch 1, but now with the full q = 2:6 range available, the
    #   mean charge state <q> = 5.9953 -- essentially identical to Branch 1's qRange = 4:6 result -- confirming
    #   that Branch 1's near-total ionization was not an artifact of its restricted charge-state range. Same
    #   built-in pressure cross-check matched again (0.0010282945017391558 vs 0.0010282945017391551).
    #   Independent hand-check of the solver's own Saha-Boltzmann formula (Plasma.computeIonLevelNumberDensity):
    #   recomputing n(q=6)/n(q=5,ground) from the converged beta*mu and the level energies gives 910.33, vs the
    #   code's own reported ratio of 910.38 (from 4-sig-fig densities) -- agreement to 0.006%. Note: the log's
    #   "betamu" is literally beta*mu_e (dimensionless), not mu_e itself (module-Plasma-inc-saha-boltzmann-
    #   mixture.jl:807); using it as mu_e directly (my first attempt) gives an answer wrong by a factor of T.
    grid        = Radial.Grid(Radial.Grid(true), rnt = 1.0e-4, h = 5.0e-2, hp = 0., rbox = 100.0)
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.2, 1.0)]
    scheme      = Plasma.SahaBoltzmannScheme(Basics.NoPlasmaModel(), true, false, 2:6, 10000, 2, 4, ionMixture, String[])

    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=settings)
    wb          = perform(wa, output=true)

elseif  true
    # Last successful:  18-Jul-2026
    # Branch 3: same carbon mixture, q = 2:6, but at T = 50 eV instead of 200 eV -- chosen (via NIST/CODATA
    # sequential ionization energies for carbon: 11.26, 24.38, 47.89, 64.49, 392.09, 489.99 eV) to sit between
    # the 3rd (q=2->3, 47.89 eV) and 4th (q=3->4, 64.49 eV) ionization steps, originally expecting q = 2, 3, 4 to
    # all carry real weight instead of one charge state dominating as in Branches 1/2 at T = 200 eV. Reuses the
    # ion-level data generated by Branch 2 (same NoExcitations = 2, upperShellNo = 4, qRange = 2:6) via the cached
    # file below, so only the Saha-Boltzmann solve itself needs to run (no fresh atomic-structure computation).
    # Verified (18-Jul-2026): clean run (9 iterations to converge, vs 4 for Branches 1/2 -- a more genuinely
    #   competitive balance takes more iterations, as expected). Cache hit confirmed
    #   (">>> Ion-level data are found in work/IonicLevelDataZ6A12-q2to6-n4-noexc2.jld for [6, 5, 4, 3, 2]+.").
    #   Result was more interesting than the naive "kT between the 3rd/4th IP" expectation above: the charge-state
    #   distribution is q = 6:5:4:3:2 = 6.05% : 66.37% : 26.83% : 0.734% : 0.0091%, i.e. dominated by q=5 (H-like)
    #   and q=4 (He-like), NOT q=2/3 as first guessed. This is real shell-structure "bottleneck" physics: the huge
    #   392 eV gap from q=4 to q=5 and further to bare (490 eV) is hard to cross at T=50 eV, so population piles up
    #   at the closed-K-shell configurations once it gets there, while the much smaller L-shell gaps (25-65 eV)
    #   below q=4 are comparatively easy to cross -- the same low-density-favours-high-ionization push from
    #   Branches 1/2 still operates, it just gets trapped behind the K-shell bottleneck instead of reaching bare.
    #   Same built-in pressure cross-check matched exactly again (0.00021230960726116257 both ways).
    #   Independent hand-check of two ratios via Plasma.computeIonLevelNumberDensity's own formula, using the
    #   converged beta*mu and the (exact, 1-electron) q=5 ground-state energy: n(q=6)/n(q=5,gs) theory=0.09171
    #   vs code=0.09174 (0.035% off); n(q=5,gs)/n(q=4,gs) theory=2.7289 vs code=2.7307 (0.065% off, using only a
    #   4-significant-figure q=4 ground energy) -- both confirm the solver in this more competitive regime too.
    grid        = Radial.Grid(Radial.Grid(true), rnt = 1.0e-4, h = 5.0e-2, hp = 0., rbox = 100.0)
    nm          = Nuclear.Model(6.0)
    rho         = 2.0e-5      # number density, unchanged from Branches 1/2
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 50.)
    settings    = Plasma.Settings(temp_au, rho, true)
    ionMixture  = [IsotopicFraction(6., 12.2, 1.0)]
    scheme      = Plasma.SahaBoltzmannScheme(Basics.NoPlasmaModel(), true, true, 2:6, 10000, 2, 4, ionMixture,
                                              ["work/IonicLevelDataZ6A12-q2to6-n4-noexc2.jld"])

    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=settings)
    wb          = perform(wa, output=true)

end
#
setDefaults("print summary: close", "")

# Filenames generated
# "IonicLevelDataZ6A12.jld", "IonicLevelDataZ9A20.jld"
