
println("Cl) Apply & test the StarkShift module with ASF from an internally generated multiplet.")

if  false
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
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(1., "point")

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
    # Last visit:  31-Jul-2026
    # Last successful:  unknown ...
    # Branch b: Li [He]2s Stark shift -- same internal-consistency idea as branch a, but on a
    #   genuinely different (many-electron, frozen-core) system: reuses the Li perturber multiplet
    #   ([He]2p..[He]10p) already validated in example-Cj.jl branch c against alpha_0(Li,2s) =
    #   164.0740(5) a.u. (Puchalski et al.), where JAC's frozen-core single-CSF treatment landed at
    #   ~37% of that value. J=1/2 again, so alpha_2 must be exactly 0. here too.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(3., "point")
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
    # Last visit:  31-Jul-2026
    # Last successful:  unknown ...
    # Branch c: barium [Xe]6s8s ^1S_0 and ^3S_1 -- NEW system, NEW external-literature target: van
    #   Leeuwen & Hogervorst, Z. Phys. A 316, 149 (1984) (examples/papers/1984.zpa-vanLeuven-stark-
    #   effect.pdf), Table 1, first two rows: alpha_0^exp(6s8s 1S0) = 4.68(10) MHz/(kV/cm)^2,
    #   alpha_0^exp(6s8s 3S1) = 5.68(11) MHz/(kV/cm)^2, alpha_2^exp(6s8s 3S1) = 0.000(6) MHz/(kV/cm)^2.
    #   J=1 (the 3S1 state) is the FIRST level in this whole example file where a nonzero tensor
    #   alpha_2/shift is even mathematically possible (J<1 forces it to exactly 0.) -- the natural
    #   place to first exercise that part of the formula, even though the paper's own measured
    #   alpha_2 for this state is itself consistent with 0.
    #   SIMPLIFICATION (flagged, not fixed): gMultiplet uses only the 8s->np excitation channel
    #   ([Xe] 6s8p/6s9p/6s10p), treating the 6s electron as a frozen spectator -- the paper's own
    #   treatment also needs the 6s->6p channel plus real CI. Agreement is NOT expected to match Li's
    #   case; this is a genuinely harder, two-valence-electron system.
    #   Unit conversion: alpha[a.u.] -> MHz/(kV/cm)^2 via DeltaE[Hz] = -(1/2) alpha[a.u.] *
    #   (E[V/cm]/StarkShift.AU_EFIELD_IN_VCM)^2 * (Hartree-to-Hz), evaluated at E=1 kV/cm=1000 V/cm
    #   and converted to MHz, then divided by (E in kV/cm)^2 = 1 -- i.e. numerically the same
    #   conversion factor the paper's own units already encode.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(56., "point")

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
    # Last visit:  31-Jul-2026
    # Last successful:  unknown ... (expected large discrepancy -- see comment below, Rule 7 blank on
    #   purpose, not merely "not yet checked").
    # Branch d: calcium [Ar]3d^2 ^3P_0 and ^3P_2 -- NEW external-literature target: van Leeuwen &
    #   Hogervorst (1984), Table 4: alpha_0^exp(3P0) = 6.5, alpha_0^exp(3P2) = -1600(100),
    #   alpha_2^exp(3P2) = 1400(100)  [MHz/(kV/cm)^2]. The paper's OWN sophisticated Coulomb-
    #   approximation calculation lands orders of magnitude BELOW the 3P2 experimental value (their
    #   own Table 4 calc column: 0.30 vs -1600(100) for alpha_0) -- explicitly attributed to admixture
    #   of highly polarizable 4snl Rydberg states into the nominally "pure" 3d^2 states, physics no
    #   simple single-configuration treatment (theirs OR this one) captures. This branch is included
    #   deliberately as an honest large-discrepancy case, not because agreement is expected: JAC's
    #   frozen-core, single-active-electron gMultiplet here is structurally even simpler than the
    #   paper's own already-failing calculation, so an even larger gap is the expected, correct
    #   outcome -- not a sign anything is broken.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(20., "point")

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
