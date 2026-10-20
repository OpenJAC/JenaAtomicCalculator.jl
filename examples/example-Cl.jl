
println("Cl) Apply & test the StarkShift module with ASF from an internally generated multiplet.")

if  false
    # Last visit:  18-Oct-2026
    # 16-Oct-2026: RUNS (exit 0) after the box was matched, and CANNOT BE DATED.
    #   H(1s) with np perturbers, n = 2..5: alpha_0 = 2.0699 a.u. in both gauges against the EXACT 4.5,
    #   i.e. 46 % of the right answer. alpha_2 = 0.0 exactly, correct for J = 1/2. This is the known
    #   MultipolePolarizibility shortfall, not a fault of this branch: the module is blocked on the
    #   B-spline pseudo-continuum problem, so no branch of this file can be dated until that moves.
    # Last visit:  31-Jul-2026
    # Last successful:  unknown ...
    # Branch a: hydrogen 1s Stark shift in a static field -- reuses the SAME independently-generated
    #   np perturber multiplet (2p..5p) already validated in examples/example-Cj.jl branch a (H(1s)
    #   alpha_0 there: 2.0697 a.u. Babushkin gauge, ~46% of the exact 4.5 a.u. -- the discrete-only
    #   sum, continuum not included; see example-Cj.jl branches a/b for the fuller picture). J=1/2
    #   here, so alpha_2 (and hence the tensor shift) must come out EXACTLY 0. -- this branch is
    #   primarily an internal-consistency/unit-conversion check of the NEW Stark-shift-from-
    #   alpha_0/alpha_2 machinery itself, not a new external-literature test (alpha_0 was already
    #   externally checked in Cj.jl). Field strength E = 100 kV/cm is used throughout this file,
    #   matching the van Leeuwen & Hogervorst (1984) apparatus' own max field (see branches c/d).
    # NOT Radial.Grid(true): its box reaches 614 a.u., and Bsplines.checkGridRepresentation refuses it
    # for H 1s + np perturbers, n = 2..5 -- a box that is much TOO LARGE starves the fixed number of
    # B-splines exactly as badly as one too small (Rule 12). 129.0 a.u. is the box the guard itself
    # names for these subshells.
    grid = Radial.Grid(Radial.Grid(true); rbox = 129.0)
    nm   = Nuclear.Model(1., PointNucleus())

    wa1 = Atomic.Computation(Atomic.Computation(), name="H perturber multiplet (np, n=2..5)", grid=grid,
                            nuclearModel=nm, configs=[Configuration("2p"), Configuration("3p"), Configuration("4p"), Configuration("5p")])
    wb1 = perform(wa1; output=true)
    gMultiplet = wb1[ResultKeys.Multiplet]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="H(1s) Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcome = wb[StarkShift.Settings][1]
    println("\n>> alpha_0 = $(outcome.alpha0) a.u.")
    println(">> alpha_2 = $(outcome.alpha2) a.u.   (expect exactly 0. for J=1/2)")
    for  sub in outcome.Jsublevels
        println(">> M=$(sub.M):  Delta-E = $(sub.energy)  Hartree")
    end
    #
elseif false
    # Last visit:  30-Aug-2026 -- RUNS (exit 0) after its box was matched on 30-Aug; it was the fourth
    #   branch of this file and the one item 76 never named. Li: alpha_0 = 61.26 a.u. Not dated, for the
    #   same reason as the branch above -- the module's absolute scale is not trustworthy yet.
    # Last visit:  31-Jul-2026
    # Last successful:  unknown ...
    # Branch b: Li [He]2s Stark shift -- same internal-consistency idea as branch a, but on a
    #   genuinely different (many-electron, frozen-core) system: reuses the Li perturber multiplet
    #   ([He]2p..[He]10p) already validated in example-Cj.jl branch c against alpha_0(Li,2s) =
    #   164.0740(5) a.u. (Puchalski et al.), where JAC's frozen-core single-CSF treatment landed at
    #   ~37% of that value. J=1/2 again, so alpha_2 must be exactly 0. here too.
    # NOT Radial.Grid(true): 614 a.u. against the 359 the guard names for this branch's subshells, and the
    # exponential family cannot be tuned to 359 (its achievable r_max is quantised). Same fix as the other
    # three branches of this file; hp = rbox/300 is Basics.recommendedGrid's own recipe.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 1.20, rbox = 359.0)
    nm   = Nuclear.Model(3., PointNucleus())
    nMaxLi = 10

    wa1 = Atomic.Computation(Atomic.Computation(), name="Li perturber multiplet ([He]np, n=2..$nMaxLi)", grid=grid,
                            nuclearModel=nm, configs=[Configuration("[He] $(n)p") for n = 2:nMaxLi])
    wb1 = perform(wa1; output=true)
    gMultiplet = wb1[ResultKeys.Multiplet]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Li [He]2s Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[He] 2s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcome = wb[StarkShift.Settings][1]
    println("\n>> alpha_0 = $(outcome.alpha0) a.u.   (example-Cj.jl branch c found 61.20 a.u. Babushkin, ~37% of 164.074)")
    println(">> alpha_2 = $(outcome.alpha2) a.u.   (expect exactly 0. for J=1/2)")
    for  sub in outcome.Jsublevels
        println(">> M=$(sub.M):  Delta-E = $(sub.energy)  Hartree")
    end
    #
elseif false
    # Last visit:  18-Oct-2026
    # Last successful:  unknown ... -- RUNS (exit 0) and is the CLEAREST demonstration of why this file
    #   cannot be dated. Ba 6s8s: alpha_0 = 13098 and 15075 a.u. against the branch's own quoted literature
    #   (van Leeuwen & Hogervorst 1984, Table 1) of 4.68 and 5.68 a.u. -- too large by a factor of about
    #   2600. The two gauges agree with each other EXACTLY, which is the point: the disagreement is with
    #   reality, not between the gauges.
    # NOT Radial.Grid(true) here. Its box is 614 a.u. and Bsplines.checkGridRepresentation refuses it for
    # Ba [Xe]6s + [Xe]6s np, n = 8..10 -- a box much TOO LARGE starves the fixed number of B-splines exactly as
    # badly as one too small (Rule 12). The guard names 359 a.u. as the box these subshells want, and the
    # EXPONENTIAL family cannot be tuned to it: its achievable r_max is quantised (151.4, 214.9, 304.9,
    # 432.7 ...), so no rbox request lands near 359. The non-exponential family honours rbox to 0.1 %.
    # hp = rbox/300 is the recipe Basics.recommendedGrid itself uses.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-06, h = 5.0e-2, hp = 1.20, rbox = 359.0)
    nm   = Nuclear.Model(56., PointNucleus())

    wa1 = Atomic.Computation(Atomic.Computation(), name="Ba perturber multiplet ([Xe]6s np, n=8..10)", grid=grid,
                            nuclearModel=nm, configs=[Configuration("[Xe] 6s 8p"), Configuration("[Xe] 6s 9p"), Configuration("[Xe] 6s 10p")])
    wb1 = perform(wa1; output=true)
    gMultiplet = wb1[ResultKeys.Multiplet]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Ba [Xe]6s8s Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[Xe] 6s 8s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcomes = wb[StarkShift.Settings]

    auToMHzPerKVcm2 = 1.0e-6 * Defaults.convertUnits("energy: from atomic to Hz", 1.0) *
                      0.5 * (1000.0/StarkShift.AU_EFIELD_IN_VCM)^2
    for  outcome in outcomes
        sym = LevelSymmetry(outcome.Jlevel.J, outcome.Jlevel.parity)
        println("\n>> Level $(string(sym)):  alpha_0 = $(outcome.alpha0) a.u.")
        println("   alpha_0 [Coulomb]   = $(outcome.alpha0.Coulomb*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_0 [Babushkin] = $(outcome.alpha0.Babushkin*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_2 [Coulomb]   = $(outcome.alpha2.Coulomb*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_2 [Babushkin] = $(outcome.alpha2.Babushkin*auToMHzPerKVcm2) MHz/(kV/cm)^2")
    end
    println("\n>> Literature (van Leeuwen & Hogervorst 1984, Table 1): 6s8s 1S0: alpha_0=4.68(10);" *
            " 6s8s 3S1: alpha_0=5.68(11), alpha_2=0.000(6)  [MHz/(kV/cm)^2]")
    #
elseif false
    # Last visit:  30-Aug-2026
    # Last successful:  unknown ... -- RUNS (exit 0) but STOPS SHORT TWICE. Ca 3d^2: the SCF reports
    #   "Maximum number of SCF iterations = 24 is reached at accuracy 5.04e-06 ... computations proceed",
    #   i.e. it hands on an UNCONVERGED field -- the AsfSettings() default ceiling of 24, the same trap
    #   that made the first attempt at priority item 50 report a wrong number. Raising maxIterationsScf
    #   is the first thing to try. Beyond that it shares the module-wide absolute-scale problem above.
    # NOT Radial.Grid(true) here. Its box is 614 a.u. and Bsplines.checkGridRepresentation refuses it for
    # Ca [Ar]3d + [Ar]3d np, n = 4..6 -- a box much TOO LARGE starves the fixed number of B-splines exactly as
    # badly as one too small (Rule 12). The guard names 167 a.u. as the box these subshells want, and the
    # EXPONENTIAL family cannot be tuned to it: its achievable r_max is quantised (151.4, 214.9, 304.9,
    # 432.7 ...), so no rbox request lands near 167. The non-exponential family honours rbox to 0.1 %.
    # hp = rbox/300 is the recipe Basics.recommendedGrid itself uses.
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-06, h = 5.0e-2, hp = 0.56, rbox = 167.0)
    nm   = Nuclear.Model(20., PointNucleus())

    wa1 = Atomic.Computation(Atomic.Computation(), name="Ca perturber multiplet ([Ar]3d np, n=4..6)", grid=grid,
                            nuclearModel=nm, configs=[Configuration("[Ar] 3d 4p"), Configuration("[Ar] 3d 5p"), Configuration("[Ar] 3d 6p")])
    wb1 = perform(wa1; output=true)
    gMultiplet = wb1[ResultKeys.Multiplet]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Ca [Ar]3d^2 Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d^2")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcomes = wb[StarkShift.Settings]

    auToMHzPerKVcm2 = 1.0e-6 * Defaults.convertUnits("energy: from atomic to Hz", 1.0) *
                      0.5 * (1000.0/StarkShift.AU_EFIELD_IN_VCM)^2
    for  outcome in outcomes
        sym = LevelSymmetry(outcome.Jlevel.J, outcome.Jlevel.parity)
        println("\n>> Level $(string(sym)):")
        println("   alpha_0 [Coulomb]   = $(outcome.alpha0.Coulomb*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_0 [Babushkin] = $(outcome.alpha0.Babushkin*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_2 [Coulomb]   = $(outcome.alpha2.Coulomb*auToMHzPerKVcm2) MHz/(kV/cm)^2")
        println("   alpha_2 [Babushkin] = $(outcome.alpha2.Babushkin*auToMHzPerKVcm2) MHz/(kV/cm)^2")
    end
    println("\n>> Literature (van Leeuwen & Hogervorst 1984, Table 4): 3d^2 3P0: alpha_0=6.5;" *
            " 3d^2 3P2: alpha_0=-1600(100), alpha_2=1400(100)  [MHz/(kV/cm)^2] -- their OWN best" *
            " calculation lands at 0.30, orders of magnitude off; large disagreement is expected here.")
    #
elseif false
    # Last visit:  19-Oct-2026
    # Last successful:  unknown ... -- RUNS (exit 0) and every number below reproduces, but there is NO external
    #   anchor for it: H 2p's ac polarizability over THIS perturber set (1s, 3s, 4s, 3d, 4d) is not a published
    #   quantity, and the literature n = 2 polarizability is dominated by the degenerate 2s-2p mixing that this
    #   set deliberately excludes. So the branch is internally consistent and externally unchecked, which is not
    #   enough for the date (Rule 7).
    #
    #   CLAIM, 19-Oct-2026:  H 2p, perturbers 1s/3s/4s/3d/4d, EField = 1.0e5 V/cm,
    #       omega -> 0      ac/dc = 0.5000000001  for every one of the six sublevels -- and the residual 1e-10
    #                       is itself right, being of order (omega/dE)^2 = (1e-6/0.0694)^2;
    #       static          alpha_0(2p_1/2) = 160.3008731119,  alpha_0(2p_3/2) = 160.3114279626,
    #                       alpha_2(2p_3/2) = -17.9917340059 a.u.;
    #       1064 nm         alpha_0 = 254.2700089569 / 254.2933778484,  alpha_2 = -29.0167262507 a.u.,
    #                       i.e. ENHANCED by 1.586 over static as omega approaches 0.069444 from below;
    #                       shifts NEGATIVE, so the trap is attractive, and the M = 1/2 to M = 3/2 spread is
    #                       5.49e-9 a.u. on a 2.68e-8 a.u. shift -- a 20 % differential light shift;
    #       omega = 0.30    alpha_0 = -12.6053496371, NEGATIVE, with the pole warning on both levels; the shifts
    #                       turn POSITIVE and the "trap" would repel.  That is the branch's point.
    # Branch e: THE AC STARK SHIFT -- the same machinery at a finite driving frequency, which is what an optical
    #   trap depth, a clock light shift and a magic wavelength all are.  Added 19-Oct-2026 with
    #   StarkShift.Settings.omega.
    #
    #   THE FACTOR OF TWO IS THE WHOLE RISK HERE, so this branch DEMONSTRATES it rather than assuming it.  A
    #   CONSTANT field of strength E and a field oscillating with PEAK amplitude E shift a level by
    #   -(1/2) alpha E^2 and -(1/4) alpha E^2 respectively, because the cycle average of cos^2 is 1/2.  So the
    #   same EField at omega = 0 and at omega -> 0 must give shifts in the ratio EXACTLY 2 : 1, and any other
    #   ratio means the two conventions have been mixed.  The suite asserts this to 1e-9; here it is printed,
    #   because it is the one number a reader of this file should check before using an ac shift.
    #
    #   H 2p GIVES BOTH TENSOR CASES AT ONCE: 2p_1/2 has J = 1/2, where J(2J-1) = 0 forces alpha_2 to vanish
    #   identically and the level shifts as a whole;  2p_3/2 carries a live tensor term, so its sublevels split,
    #   and THAT M-dependence is the differential light shift a clock or a quantum-logic gate has to control.
    #
    #   THE THREE FREQUENCIES, AND THE THIRD IS DELIBERATELY ILLEGITIMATE.  omega = 1e-6 a.u. is the convention
    #   check above.  omega = 0.0428227 a.u. is 1064 nm, the common trap wavelength, and it sits BELOW this
    #   perturber set's lowest excitation energy of 0.069444 a.u. (2p -> 3s/3d), so it is a real ac
    #   polarizability and a real trap depth.  omega = 0.30 a.u. is far ABOVE that threshold, where the true
    #   alpha^(E1,K) is COMPLEX and a real sum over a DISCRETE set cannot represent its imaginary part:
    #   `MultipolePolarizibility.checkFrequency` says exactly that, and what comes back is a
    #   principal-value-like real part with no width.  It is included BECAUSE the warning is the useful output --
    #   alpha_0 even changes sign there, to about -12.6 a.u., and a reader who ignored the warning would quote a
    #   trap that repels instead of attracting.  The ac/dc ratio is meaningless at both finite frequencies and is
    #   printed only for the omega -> 0 case.
    #
    #   NOTE THAT 1s LIES BELOW 2p, and that is deliberate and legitimate here -- unlike a dispersion
    #   coefficient, a Stark shift of an excited level is perfectly well defined, and the 1s perturber simply
    #   enters with a NEGATIVE energy denominator.  It is the single largest contribution and it has the
    #   opposite sign to all the others.
    nm         = Nuclear.Model(1., PointNucleus())
    perturbers = [Configuration("1s"), Configuration("3s"), Configuration("4s"), Configuration("3d"),
                  Configuration("4d")]
    grid       = Basics.recommendedGrid( vcat(perturbers, [Configuration("2p")]), nm )

    wa1 = Atomic.Computation(Atomic.Computation(), name="H perturbers (1s, ns, nd)", grid=grid, nuclearModel=nm,
                             configs=perturbers)
    gMultiplet = perform(wa1; output=true)[ResultKeys.Multiplet]

    EField = 1.0e5          ## V/cm;  the constant field at omega = 0, the PEAK amplitude otherwise
    acOutcomes = function(omega::Float64)
        set = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                  EField=EField, omega=omega)
        wa  = Atomic.Computation(Atomic.Computation(), name="H 2p ac Stark shift", grid=grid, nuclearModel=nm,
                                 configs=[Configuration("2p")], propertySettings=[set])
        return( perform(wa; output=true)[StarkShift.Settings] )
    end

    dcOut = acOutcomes(0.0)
    for  (omega, label, showRatio)  in  ((1.0e-6,    "the convention check:  ac/dc must be exactly 1/2", true),
                                         (0.0428227, "1064 nm -- BELOW the 0.069444 a.u. threshold, a real trap",
                                                     false),
                                         (0.30,      "ABOVE the threshold -- the warning IS the result", false))
        println("\n>> ===== omega = $omega a.u. -- $label =====")
        for  (dc, ac)  in  zip(dcOut, acOutcomes(omega))
            sym = LevelSymmetry(ac.Jlevel.J, ac.Jlevel.parity)
            println(">>   Level $(string(sym)):  alpha_0 = $(ac.alpha0.Babushkin) a.u.,  " *
                    "alpha_2 = $(ac.alpha2.Babushkin) a.u.")
            for  (subDc, subAc)  in  zip(dc.Jsublevels, ac.Jsublevels)
                sa = ">>     M = $(subAc.M):  shift = $(subAc.energy.Babushkin) a.u."
                if  showRatio  &&  subDc.energy.Babushkin != 0.
                    sa = sa * ",   ac/dc = $(subAc.energy.Babushkin/subDc.energy.Babushkin)"
                end
                println(sa)
            end
        end
    end
    println("\n>> -DeltaE at the trap frequency IS the trap depth of that sublevel, and the SPREAD over M is " *
            "the differential light shift.  A magic wavelength, where that shift vanishes for a PAIR of levels, " *
            "is branch g.")
    #
elseif false
    # Last visit:  20-Oct-2026
    # Last successful:  unknown ... (Rule 7: the 831.9 V/m anchor is external and the two scaling laws are exact,
    #   but the shift itself is H 2p over THIS perturber set, which is not a published quantity -- so the physics
    #   judgement, and the date, are the maintainer's)
    #
    #   CLAIM, 20-Oct-2026:  H 2p_1/2, perturbers 1s/3s/4s/3d/4d, alpha_0(0) = 160.3008730798 a.u.,
    #       E_rms(300 K) = 831.9451906740 V/m   against the external 831.9 V/m -- the one anchor here;
    #       shift/leading = 1.000000038 at 1 K, so the T -> 0 collapse onto -(2/15)(alpha_fs pi)^3 T^4 alpha_0 holds;
    #       eta POSITIVE throughout the usable range:  3.776117e-08 (1 K), 3.398507e-07 (3 K), 3.776148e-06 (10 K),
    #           3.398762e-05 (30 K), 3.779303e-04 (100 K), 3.424625e-03 (300 K), 4.160290e-02 (1000 K);
    #       eta ratio / (T2/T1)^2 = 1.000001, 1.000008, 1.000067, 1.000768, 1.006836 from 1 to 300 K -- the T^2 law
    #           to 7e-3 over two and a half decades, which is the assertion a wrong frequency factor cannot pass;
    #       300 K shift = -2.1051401e-16 a.u. (J = 1/2) and -2.1052791e-16 (J = 3/2), against dc Stark shifts of
    #           -3.03e-08 a.u. at 1e5 V/cm -- a ratio of 7e-9, and the J = 3/2 sublevels share ONE bbr shift while
    #           carrying four different dc shifts, which is the rank-0-only structure of isotropic radiation;
    #       and 3000 K / 30000 K are PAST the guard, as intended:  eta 2.201325 and -1.066159, the shift turning
    #           POSITIVE at 30000 K, with the warning naming 6.2 % and 98.5 % of the Planck weight above the lowest
    #           excitation energy of 0.069446 a.u.
    #   NOTE ON READING THE OUTPUT:  the sweep is computed before it is printed, so both Planck-spectrum warnings
    #       appear ABOVE the table rather than beside the 3000 K and 30000 K rows they belong to.
    # Branch f: THE BLACK-BODY RADIATION SHIFT -- the one shift a clock cannot switch off.  Added 19-Oct-2026 as
    #   StarkShift.Settings.temperature, and reached either through `perform` or directly through
    #   StarkShift.computeBlackbodyShift.
    #
    #   IT IS AN AVERAGE OF alpha_0(omega) OVER THE PLANCK SPECTRUM, nothing more:
    #
    #       DeltaE_BBR = -(1/2) <E^2> alpha_eff,    alpha_eff = int x^3 alpha_0(Tx)/(e^x-1) dx / int x^3/(e^x-1) dx
    #
    #   so it carries NO second factor of two -- <E^2> is already a mean square, unlike the ac shift of branch e,
    #   where the peak amplitude costs one.  Only the RANK-0 polarizability enters, because thermal radiation is
    #   isotropic and unpolarized: there is no axis for a rank-2 term to be defined against, and the sublevels of
    #   a J = 3/2 level are therefore shifted ALIKE.  That is the one structural difference from branch e.
    #
    #   THREE THINGS ARE WORTH CHECKING HERE, AND THE FIRST IS THE ONLY EXTERNAL NUMBER IN THIS FILE.
    #     1. <E^2> needs only the fine-structure constant and the temperature, so E_rms at 300 K is a textbook
    #        831.9 V/m -- the value quoted throughout the clock literature.  If that is wrong, nothing else can
    #        be right.  Note what it means for a laboratory field: 8.3 V/cm rms against branch e's 1e5 V/cm, so
    #        the room-temperature BBR shift is some 1e-8 of that dc Stark shift.
    #     2. AS T -> 0 the shift must collapse onto the closed form -(2/15)(alpha_fs pi)^3 T^4 alpha_0(0), which
    #        is the universally quoted leading term.
    #     3. eta = alpha_eff/alpha_0(0) - 1 is the dynamic correction, and it must be POSITIVE and grow as T^2 --
    #        first order in <omega^2>/dE^2.  The sign is not decoration: a quadrature built with the
    #        IMAGINARY-axis factor dE/(dE^2 + omega^2) instead of the real dE/(dE^2 - omega^2) satisfies the
    #        T -> 0 limit to 1e-13 and returns eta of the OPPOSITE sign away from it.  The limit alone would not
    #        have caught it; the scaling did.
    #
    #   EXPECT THE SATURATION WARNING, and do not read past it: five perturbers leave 8.6 % of alpha_0 in the
    #   highest one, so alpha_0 = 160.3 a.u. is a LOWER BOUND of unknown quality and every shift below inherits
    #   that.  It is deliberate -- the set is branch e's, so the two files' numbers can be compared -- and it is
    #   why this branch's claim is the two LAWS rather than the shift: a ratio of shifts at two temperatures, and
    #   eta's T^2 growth, are independent of how good alpha_0 is.
    #
    #   AND AT HIGH T THE WHOLE THING STOPS MEANING ANYTHING, which is why the sweep runs up to 30000 K: once the
    #   Planck spectrum reaches the level's own excitation energies, alpha_0 is being sampled THROUGH its poles
    #   and StarkShift.checkBlackbodyResonance says so.  The warning is the useful output, exactly as in branch e.
    nm         = Nuclear.Model(1., PointNucleus())
    perturbers = [Configuration("1s"), Configuration("3s"), Configuration("4s"), Configuration("3d"),
                  Configuration("4d")]
    grid       = Basics.recommendedGrid( vcat(perturbers, [Configuration("2p")]), nm )

    wa1 = Atomic.Computation(Atomic.Computation(), name="H perturbers (1s, ns, nd)", grid=grid, nuclearModel=nm,
                             configs=perturbers)
    gMultiplet = perform(wa1; output=true)[ResultKeys.Multiplet]

    ## 1. the external anchor, which needs no atom at all
    alphaFs = 1.0 / Defaults.getDefaults("speed of light: c")
    eRms    = sqrt( StarkShift.blackbodyFieldSquare(300.0) ) * StarkShift.AU_EFIELD_IN_VCM * 100.0
    println("\n>> ===== the black-body field itself =====")
    println(">>   E_rms(300 K) = $eRms V/m   against the Stefan-Boltzmann 831.9 V/m")
    println(">>   i.e. $(eRms/100.0) V/cm rms, against the 1.0e5 V/cm dc field of branch e")

    ## the user-facing route: temperature on the Settings, the shift on the Outcome
    set = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                              EField=1.0e5, temperature=300.0)
    wa  = Atomic.Computation(Atomic.Computation(), name="H 2p black-body shift", grid=grid, nuclearModel=nm,
                             configs=[Configuration("2p")], propertySettings=[set])
    outcomes = perform(wa; output=true)[StarkShift.Settings]
    println("\n>> ===== through perform, at 300 K =====")
    for  oc  in  outcomes
        sym = LevelSymmetry(oc.Jlevel.J, oc.Jlevel.parity)
        println(">>   Level $(string(sym)):  bbrShift = $(oc.bbrShift.Babushkin) a.u.,  eta = $(oc.bbrEta)")
        println(">>     dc Stark shifts of the sublevels, for scale:  " *
                join([string(sub.energy.Babushkin) for sub in oc.Jsublevels], ", "))
    end

    ## 2. and 3. the two laws, taken directly so that one sum over states serves every temperature
    level  = outcomes[1].Jlevel
    st, en = StarkShift.polarizabilityTerms(level, gMultiplet, grid)
    alpha0 = StarkShift.rebuiltPolarizability(st, en, 0.0)
    println("\n>> ===== the two scaling laws;  alpha_0(0) = $alpha0 a.u. =====")
    println(">>   T [K]      shift [a.u.]        shift/leading       eta              eta ratio / (T2/T1)^2")
    ## the sweep is collected first and printed afterwards, so that no accumulator is assigned inside a
    ## top-level loop -- Julia's soft scope makes such a variable LOCAL and the read of it then fails
    temps = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 30000.0]
    rows  = [ StarkShift.computeBlackbodyShift(level, gMultiplet, grid, TK)  for TK in temps ]
    for  k = 1:length(temps)
        Tau  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temps[k])
        lead = -2.0/15.0 * (alphaFs*pi)^3 * Tau^4 * alpha0
        sa   = @sprintf(">>   %-8.0f  %-18.6e %-19.9f %-16.6e", temps[k], rows[k][1].Babushkin,
                        rows[k][1].Babushkin/lead, rows[k][2])
        if  k > 1   sa = sa * @sprintf(" %.6f", (rows[k][2]/rows[k-1][2]) / (temps[k]/temps[k-1])^2)   end
        println(sa)
    end
    println(">>   the last column must be 1 while eta is small;  where it is not, the spectrum has reached the " *
            "level's own excitation energies and the resonance warning above says so.")
    #
elseif true
    # Last visit:  20-Oct-2026
    # Last successful:  unknown ... (Rule 7: the root is reproduced from an independent dense scan and from
    #   tools/probe-magicWavelength.jl, but a magic wavelength of H 1s -- 2p is not a measured quantity, so the
    #   date is the maintainer's)
    #
    #   CLAIM, 20-Oct-2026:  H 1s against H 2p_1/2, perturbers 1s/2p/3p/4p/3s/4s/3d/4d,
    #       alpha_0(1s, 0) = 3.4918282527 a.u.,  alpha_0(2p_1/2, 0) = 160.3008730798 a.u.,  ten poles in range in
    #           five near-degenerate PAIRS:  0.069446, 0.069446, 0.093752, 0.093752, 0.375005, 0.375006, 0.444450,
    #           0.444451, 0.468756, 0.468757 a.u.;
    #       ONE magic wavelength in omega = 0.01 ... 0.25 a.u.:  omega = 0.090592606 a.u. = 502.948 nm, where
    #           alpha_common = 3.697721110 a.u. with residual |alpha_A - alpha_B|/alpha = 3.45e-09, slope
    #           -7.7900e+04 a.u./a.u., and aboveResonance TRUE -- it lies above 2p's lowest excitation energy;
    #       2 candidates discarded as POLES rather than zeros;
    #       and the 200 000-point dense scan finds exactly 1 pole-free sign change against 2 that straddle a pole,
    #           so the single returned root is COMPLETE with nothing spurious.
    #       Every figure here also comes out of tools/probe-magicWavelength.jl, which reaches them by a different
    #           driver -- two independent routes to the same root.
    # Branch g: THE MAGIC WAVELENGTH of a PAIR of levels -- the zero of alpha_0(omega, A) - alpha_0(omega, B),
    #   where the trapping light shifts both levels equally and the TRANSITION between them is untouched to first
    #   order in the intensity.  That is what lets an optical trap coexist with a clock or a qubit, and it is the
    #   quantity branch e's differential light shift is engineered to cancel.  Added 19-Oct-2026 as
    #   StarkShift.magicWavelengths, which returns a LIST: there is generically at least one zero between every
    #   pair of consecutive resonances of either level.
    #
    #   WHAT MAKES THIS HARDER THAN A ROOT FIND, and it is the reason the function exists rather than a one-line
    #   bisection: alpha_0 runs to +infinity below every resonance and returns from -infinity above it, so the
    #   DIFFERENCE changes sign at every pole of EITHER level.  A bisection let loose on such a sign change
    #   converges onto the POLE and reports a resonance as a magic wavelength -- the one error that looks entirely
    #   plausible, since what comes back is a sharp, reproducible frequency.  Two things are therefore done, and
    #   the second was found by measurement rather than foreseen:  the range is CUT at every known pole, since the
    #   poles are exactly the excitation energies; and every surviving candidate is verified by its RESIDUAL
    #   |alpha_A - alpha_B|, because the poles come in near-degenerate PAIRS -- 2p -> 3s and 2p -> 3d sit 1e-7
    #   apart -- whose sliver is below numerical resolution, so the sign change inside it is the pair of poles and
    #   not a zero.  The two populations are seven orders apart in residual, so nothing delicate is being decided.
    #
    #   READ THE TWO FLAGS BEFORE USING A ROOT.  `alphaCommon` NEGATIVE is a trap that REPELS rather than holds,
    #   and `aboveResonance` true says the zero sits where the true polarizability is complex and this real sum
    #   over a discrete set has no imaginary part to offer -- a trap there heats the atom by photon scattering,
    #   which no number in this table shows.
    #
    #   WHAT THIS PAIR DOES NOT SHOW, said plainly:  every zero here lies ABOVE the lowest resonance, because
    #   alpha_0(1s) = 3.49 and alpha_0(2p) = 160.3 do not cross below it -- the difference is negative at
    #   omega -> 0 and runs to -infinity at 2p's first pole.  The case a clock actually wants, a magic wavelength
    #   BELOW all resonances, needs a pair whose static polarizabilities are closer together, i.e. an alkali or
    #   alkaline-earth pair rather than two hydrogen levels.  The ALGORITHM is what this branch demonstrates.
    nm    = Nuclear.Model(1., PointNucleus())
    confs = [Configuration("1s"), Configuration("2p"), Configuration("3p"), Configuration("4p"),
             Configuration("3s"), Configuration("4s"), Configuration("3d"), Configuration("4d")]
    grid  = Basics.recommendedGrid(confs, nm)
    ## each level's own selection rules decide which members it reaches, so a set of BOTH parities is right here
    ## and the unreachable members drop out on their own
    one   = function(name::String, cs::Array{Configuration,1})
        perform( Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                                    configs=cs); output=true )[ResultKeys.Multiplet]
    end
    gMultiplet = one("all perturbers", confs)
    levelA     = one("H 1s", [Configuration("1s")]).levels[1]
    levelB     = one("H 2p", [Configuration("2p")]).levels[1]

    sA, eA = StarkShift.polarizabilityTerms(levelA, gMultiplet, grid)
    sB, eB = StarkShift.polarizabilityTerms(levelB, gMultiplet, grid)
    poles  = sort( unique( [dE  for dE in vcat(eA, eB)  if dE > 0.0] ) )
    println("\n>> alpha_0(1s, 0)     = $(StarkShift.rebuiltPolarizability(sA, eA, 0.0)) a.u.")
    println(">> alpha_0(2p_1/2, 0) = $(StarkShift.rebuiltPolarizability(sB, eB, 0.0)) a.u.")
    println(">> the poles of the difference, which are cut out of the search (a.u.):  " *
            "$([round(p, digits=6) for p in poles])")

    mws = StarkShift.magicWavelengths(levelA, levelB, gMultiplet, grid; omegaMin=0.01, omegaMax=0.25)
    println("\n>>   #   omega [a.u.]     lambda [nm]   alpha_common [a.u.]   residual     slope [a.u.]   above?")
    for  (i, mw)  in  enumerate(mws)
        aA = StarkShift.rebuiltPolarizability(sA, eA, mw.omega)
        aB = StarkShift.rebuiltPolarizability(sB, eB, mw.omega)
        println(@sprintf(">>   %-3d %-15.9f %-13.3f %-21.9f %-12.2e %-14.4e %s",
                         i, mw.omega, mw.wavelength, mw.alphaCommon, abs(aA-aB)/max(abs(aA), 1.0e-30),
                         mw.slope, mw.aboveResonance))
    end
    println(">>   a NEGATIVE alpha_common would repel;  `above?` true means the real sum has no width there.")

    ## AN INDEPENDENT CHECK, because a root find that returns nothing and a range with no zero look alike: scan
    ## the range densely and count the sign changes that do NOT straddle a known pole.  That count must equal the
    ## number of roots returned, with nothing extra.
    dAlpha = om -> StarkShift.rebuiltPolarizability(sA, eA, om) - StarkShift.rebuiltPolarizability(sB, eB, om)
    nScan   = 200000;   omegas = range(0.01, 0.25, length=nScan)
    ## collected rather than counted in place:  a counter incremented inside a top-level loop is LOCAL under
    ## Julia's soft scope, and reading it raises
    changes = [ (omegas[k], omegas[k+1])  for k = 1:nScan-1
                                          if dAlpha(omegas[k]) * dAlpha(omegas[k+1]) < 0.0 ]
    nAtPole = count( any(a <= pole <= b  for pole in poles)  for (a, b) in changes; init=0 )
    nFree   = length(changes) - nAtPole
    println("\n>> the dense scan of $nScan points finds $nFree pole-free sign change(s) and $nAtPole across a " *
            "pole;  the function returned $(length(mws)) root(s).  Equal counts mean COMPLETE and nothing " *
            "spurious -- this is the check that the partition and the residual test are both doing their work.")
    #
end
