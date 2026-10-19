
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
elseif true
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
            "the differential light shift.  A magic wavelength is a zero of alpha_0(omega, upper) - " *
            "alpha_0(omega, lower), found by running two levels over a range of omega.")
    #
end
