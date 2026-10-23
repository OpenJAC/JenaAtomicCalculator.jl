#
println("Jc) Apply & test the Saha-Boltzmann computations for an ionic mixture in LTE.")

setDefaults("print summary: open", "zzz-saha-boltzmann.sum")


if  false
    # Last visit:  20-Oct-2026
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
    # BOX NOTE (priority item 77, MEASURED 17-Oct-2026) -- the box is adequate and was left alone.
    #   `checkOrbitalBox` flags 4d_3/2, 4d_5/2 (0.995) and 4f_5/2 (0.955) in this 10.0 a.u. box, and item 77
    #   first put this branch at 9.6x oversize, the worst of the eight it listed.  IT IS IN FACT THE SAFEST,
    #   and the reason is a lesson worth keeping: THE WARNING NAMES AN ORBITAL LABEL, BUT THE DAMAGE DEPENDS
    #   ON THE CONFIGURATION IT SITS IN.  The sweep log shows these 4d/4f warnings come from `4d_3/2^1`
    #   (q = 6, BARE carbon -- one electron at Z = 6) and `1s 4d` (q = 5), not from a neutral atom:  this
    #   branch runs q = 4:6, He-like through bare.  Those orbitals are hydrogenic and COMPACT.
    #   MEASURED against the EXACT Dirac energies, which a one-electron system makes available:
    #       4s_1/2 0.006 %   4p_3/2 0.004 %   4d_3/2 0.002 %   4d_5/2 0.002 %   4f_7/2 0.000 %
    #   A neutral-carbon 4d, by contrast, is out by 429 % in this box -- a Z_eff ~ 1 Rydberg orbital is a
    #   different object entirely.  So: adequate here, and read the CONFIGURATION before trusting a ratio.
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
    # Last visit:  18-Oct-2026
    # Last reproduced:  18-Oct-2026 (1 of 6 to 1e-5, 2 more to 1e-3;  3 are NIST ionization energies, not results)
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

elseif  false
    # Last visit:  23-Oct-2026   ## <q> = 5.9548332595 (n_i=2.0e-4) and 5.8020831296 (n_i=1.0e-3) for ALL FOUR
    #   physical models AND NoPlasmaModel;  the cEK = -1.0 mutation gives 5.4624239799 above n_cr and is identical
    #   below it;  n_cr = 2.0316e-03 a.u.;  ten perform runs in 97 s
    # Last successful:  unknown ...
    # Branch d: FIVE IPD MODELS THROUGH perform, AND THE MEASUREMENT THAT THEY ALL GIVE THE SAME ANSWER.
    #   Until 23-Oct-2026 no example exercised ANY ionization-potential-depression model through the perform
    #   route -- every Saha-Boltzmann branch above passes NoPlasmaModel() -- so Basics.DebyeHueckelModel,
    #   Basics.IonSphereModel, Basics.StewartPyattModel and the new Basics.EckerKroellModel were reachable only
    #   by calling Plasma.determineIpShifts directly.  This branch closes that gap, and in doing so it found
    #   something that no direct call to determineIpShifts could have shown.
    #
    #   WHAT THIS BRANCH ACTUALLY ESTABLISHES, and the third item is the point:
    #   (i)   all five models run end-to-end through Plasma.Computation + perform without error;
    #   (ii)  the Ecker-Kroell model CROSSES ITS OWN CRITICAL DENSITY WHILE THE SAHA ITERATION CONVERGES at
    #         n_i = 1.0e-3: iteration 1 starts from the seeded electron density (n = 1.12e-03, weak branch) and
    #         the converged state sits at n = 6.80e-03 (strong branch), with n_cr = 2.0316e-03.  Watch the
    #         printed ">> Ecker-Kroell IPD:" lines, which report n, n_cr, the branch and the shifts themselves
    #         once per iteration.  The branch is decided by the SELF-CONSISTENT density, not by the input one.
    #   (iii) THE MEAN CHARGE STATE IS INSENSITIVE TO THE CHOICE OF IPD MODEL -- all five give an IDENTICAL <q>
    #         to every digit, and identical to NoPlasmaModel, i.e. to no depression at all.
    #
    #   (iii) IS NOT A PROPERTY OF THE MODELS BUT OF HOW THEIR SHIFT IS CONSUMED, and it was calibrated BY
    #   MUTATION rather than argued.  Plasma.solveSahaBoltzmannEquilibrium uses deltaIp only as a bound/unbound
    #   switch -- a level keeps its population if  energy + deltaIp < groundE[q+1]  -- and every physical model
    #   returns deltaIp NEGATIVE, which ADDED to an already negative level energy makes the test EASIER to pass.
    #   So no level is ever removed, however large the depression.  Measured, all at n_i = 1.0e-3:
    #       NoPlasmaModel (deltaIp = -1000, "never cut")          <q> = 5.8020831296
    #       WithoutAutoionizationModel (deltaIp = 0)              <q> = 5.8020831296
    #       DebyeHueckelModel(0.01)  -- deltaIp ~ -600 Hartree    <q> = 5.8020831296
    #       EckerKroellModel with cEK = +1.0 (shipped)            <q> = 5.8020831296
    #       EckerKroellModel with cEK = -1.0 (MUTATED, deltaIp>0) <q> = 5.4624239799
    #   A deliberately enormous NEGATIVE shift changes nothing;  reversing the sign moves <q> by 0.34.  That is
    #   a guard which cannot fail in the shipped direction, and it is what this branch exists to record.
    #   AND THE MUTATION BITES ONLY ABOVE n_cr -- it moves <q> by -0.3397 at n_i = 1.0e-3 and by EXACTLY ZERO at
    #   n_i = 2.0e-4 -- which was not designed and is the sharper test of the two: cEK enters only the strong
    #   branch, so a mutation of cEK that changed anything below n_cr would mean the switch was not working.
    #   The sign is NOT this branch's to fix: it is shared by DebyeHueckel, IonSphere and StewartPyatt, it
    #   predates the Ecker-Kroell model, and the accompanying normalization question (the partition function is
    #   formed from the FULL level list before any level is cut, so a cut does not conserve particles) is a
    #   modelling decision rather than a typo.  It is on the priority list.
    #
    #   NO LITERATURE CLAIM is made here.  What is dated is (i), (ii) and the five numbers of (iii).
    #   COST: ten perform runs, MEASURED at 97 s in total on 23-Oct-2026 -- the first run carries the compilation
    #   and the rest are about ten seconds each.  This is the cheap branch of the file.
    grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    nm          = Nuclear.Model(6.0)
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 200.)
    ionMixture  = [IsotopicFraction(6., 12.011, 1.0)]
    ncr         = 3/(4pi) * (temp_au / 6^2)^3
    println("\n>> Ecker-Kroell critical density for kT = 200 eV, qmax = 6:  n_cr = $ncr a.u.")
    println(">> n_i = 2.0e-4 therefore sits BELOW it (weak branch) and n_i = 1.0e-3 ABOVE it (strong branch).\n")
    #
    # The models are built with zero parameters;  Plasma.updateIpdPlasmaModel fills in the Debye length and the
    # electron density from the converged state at each Saha iteration, exactly as it does for Stewart-Pyatt.
    # cEK = -1.0 is the MUTATION of (iii) above and is included deliberately: without it the four physical models
    # agree and nothing says whether that is physics or a dead switch.
    models = [ ("NoPlasmaModel       ", Basics.NoPlasmaModel()                         ),
               ("DebyeHueckelModel   ", Basics.DebyeHueckelModel(0.0)                  ),
               ("StewartPyattModel   ", Basics.StewartPyattModel(0.0, 0.0, 0.0)        ),
               ("EckerKroell cEK=1.0 ", Basics.EckerKroellModel(0.0, 0.0, 0.0,  1.0)   ),
               ("EckerKroell cEK=-1.0", Basics.EckerKroellModel(0.0, 0.0, 0.0, -1.0)   ) ]
    #
    results = Dict{Tuple{Float64,String},Float64}()
    for  rho  in  [2.0e-4, 1.0e-3]
        for  (label, pm)  in  models
            settings = Plasma.Settings(temp_au, rho, true)
            scheme   = Plasma.SahaBoltzmannScheme(pm, true, false, 4:6, 10000, 1, 4, ionMixture, String[])
            wa       = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                           settings=settings)
            wb       = perform(wa, output=true)
            results[(rho, label)] = wb["mean charge state"]
        end
    end
    #
    println("\n>> Mean charge state <q>, by IPD model and density:\n")
    println("   n_i [a.u.]    " * join([m[1] for m in models], "   "))
    for  rho  in  [2.0e-4, 1.0e-3]
        println("   $rho       " * join([@sprintf("%18.10f", results[(rho, m[1])]) for m in models], " "))
    end
    #
    println("\n>> THE FOUR PHYSICAL MODELS MUST AGREE WITH NoPlasmaModel, AT BOTH DENSITIES -- not because the")
    println("   depression is negligible but because the shipped sign of deltaIp cannot unbind a level:")
    for  rho  in  [2.0e-4, 1.0e-3]
        agree = all(results[(rho, m[1])] == results[(rho, "NoPlasmaModel       ")] for m in models[1:4])
        println("   n_i = $rho :  first four identical to NoPlasmaModel = $agree")
    end
    println("\n>> AND THE MUTATION MUST DIFFER, or the test above is vacuous:")
    for  rho  in  [2.0e-4, 1.0e-3]
        d = results[(rho, "EckerKroell cEK=-1.0")] - results[(rho, "NoPlasmaModel       ")]
        println("   n_i = $rho :  cEK = -1.0 minus NoPlasmaModel = $d")
    end

end
#
setDefaults("print summary: close", "")

# Filenames generated
# "IonicLevelDataZ6A12.jld", "IonicLevelDataZ9A20.jld"
