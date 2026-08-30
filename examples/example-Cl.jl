
println("Cl) Apply & test the StarkShift module with ASF from an internally generated multiplet.")

if  false
    # Last visit:  30-Aug-2026 -- RUNS (exit 0) after the box was matched, and CANNOT BE DATED.
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
    gMultiplet = wb1["multiplet:"]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="H(1s) Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("1s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcome = wb["Stark-shift outcomes:"][1]
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
    gMultiplet = wb1["multiplet:"]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Li [He]2s Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[He] 2s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcome = wb["Stark-shift outcomes:"][1]
    println("\n>> alpha_0 = $(outcome.alpha0) a.u.   (example-Cj.jl branch c found 61.20 a.u. Babushkin, ~37% of 164.074)")
    println(">> alpha_2 = $(outcome.alpha2) a.u.   (expect exactly 0. for J=1/2)")
    for  sub in outcome.Jsublevels
        println(">> M=$(sub.M):  Delta-E = $(sub.energy)  Hartree")
    end
    #
elseif false
    # Last visit:  30-Aug-2026
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
    gMultiplet = wb1["multiplet:"]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Ba [Xe]6s8s Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[Xe] 6s 8s")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcomes = wb["Stark-shift outcomes:"]

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
elseif true
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
    gMultiplet = wb1["multiplet:"]

    starkSettings = StarkShift.Settings(StarkShift.Settings(); calcStarkshifts=true, gMultiplet=gMultiplet,
                                         EField=1.0e5, printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Ca [Ar]3d^2 Stark shift", grid=grid, nuclearModel=nm,
                            configs=[Configuration("[Ar] 3d^2")], propertySettings=[starkSettings] )
    wb = perform(wa; output=true)
    outcomes = wb["Stark-shift outcomes:"]

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
end
