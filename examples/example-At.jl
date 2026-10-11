#
println("At) Apply & test the RandomPhase module: orchestrating one DalgarnoLewis solve per active electron.")
#
# RandomPhase is the many-electron layer built on top of DalgarnoLewis (example-As.jl) -- it is where the
# self-consistent, many-electron physics the name is for belongs, since DalgarnoLewis itself knows nothing about
# more than one electron. See the module's own docstring for the full architectural reasoning.
#
# `settings.couplingEnabled = false` (the default, branch a) is the INDEPENDENT-PARTICLE limit: one
# `DalgarnoLewis.solve` per active electron, each driven by the external operator alone. `couplingEnabled =
# true` (branch b) iterates instead, each electron's source gaining the DIRECT (Hartree, "ring") induced
# potential of every OTHER active electron's current chi (`RadialIntegrals.Yk_ab` at the external operator's
# own multipole rank, projected through `DalgarnoLewis.solve`'s `extraSource` keyword) -- SEE THE MODULE'S OWN
# DOCSTRING for exactly what this has, and has not yet, been checked against before trusting a coupled number.
#
# BRANCH (b) IS DELIBERATELY MODEST ABOUT WHAT IT CLAIMS. The one thing proved RIGOROUSLY there is that a
# SINGLE active electron has no "other" electron to induce anything, so couplingEnabled = true must reproduce
# couplingEnabled = false EXACTLY -- and does. A genuine two-electron case is shown running and returning
# finite numbers, but its CONVERGENCE on this particular toy system (H's own 1s/2s, which has no mean field
# between them to begin with) is explicitly NOT claimed; see the long comment in branch (b) itself.

if  true
    # Last successful:  11-Oct-2026
    # Branch (a): the orchestration adds no physics of its own. Two active electrons of the SAME H-like 1s/2s
    # basis are both driven by an E1 operator at omega = 0.3 a.u. into the SAME 2p_1/2 target channel, through
    # ONE RandomPhase.solve call -- and each returned chi must be bit-identical to what a bare DalgarnoLewis.solve
    # call for that one electron, alone, would have returned.
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 40.0)

    basis = SelfConsistent.performSCF([Configuration("1s 2s")], ni, grid, asfSettings; printout=false).levels[1].basis
    orb_1s = basis.orbitals[Subshell("1s_1/2")];   orb_2s = basis.orbitals[Subshell("2s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    omega = 0.3;   target = Subshell("2p_1/2")

    activeOrbitals = Dict(Subshell("1s_1/2") => orb_1s, Subshell("2s_1/2") => orb_2s)
    targetSubshells = Dict(Subshell("1s_1/2") => target, Subshell("2s_1/2") => target)
    settings = RandomPhase.Settings(RandomPhase.Settings();
                   activeSubshells = [Subshell("1s_1/2"), Subshell("2s_1/2")], targetSubshells = targetSubshells)
    response = RandomPhase.solve(E1, Basics.Coulomb, omega, activeOrbitals, potential, primitives; settings=settings)
    println(">> RandomPhase.solve():  nIterations = $(response.nIterations), converged = $(response.converged)")

    for  sh  in  (Subshell("1s_1/2"), Subshell("2s_1/2"))
        direct = DalgarnoLewis.solve(E1, Basics.Coulomb, omega, targetSubshells[sh], activeOrbitals[sh], potential, primitives)
        viaRP  = response.perturbedOrbitals[sh]
        println(">> $sh:  max|P_direct - P_viaRandomPhase| = $(maximum(abs.(direct.P - viaRP.P))), " *
                "max|Q_direct - Q_viaRandomPhase| = $(maximum(abs.(direct.Q - viaRP.Q)))   (both must be exactly 0.0)")
    end
    #
elseif  true
    # Last successful:  11-Oct-2026
    # Branch (b): THE COUPLED PATH. Part 1 is RIGOROUS: a single active electron has no coupling partner, so
    # couplingEnabled=true must reproduce couplingEnabled=false to the last bit. Part 2 runs the SAME toy
    # two-electron system as branch (a) with coupling switched on, and reports what is ACTUALLY measured, not
    # a hoped-for number: the induced SOURCE starts small (a fraction of a percent of the external one) but
    # the damped fixed-point iteration (Settings.mixing) only creeps toward self-consistency on this case,
    # tens of percent relative change still after many iterations -- not divergent, but not practically
    # convergent in a modest budget either. THE LIKELY READING, not hidden: `Basics.NuclearField` orbitals
    # carry NO mean field between 1s and 2s at all, so this toy case asks the coupling term to supply the
    # ENTIRE electron-electron interaction by itself, with the full (2j_b+1)=2 statistical weight on top -- an
    # unusually stiff fixed point, and quite possibly not representative of a real multi-electron SCF system,
    # where the orbitals already carry a mean field and the induced term is a CORRECTION to it rather than the
    # whole interaction. A genuine multi-electron SCF system with a known RPA/polarizability correction to
    # check against is the natural next validation step -- not attempted here.
    ni   = Nuclear.Model(1.0, Nuclear.PointNucleus())
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 3.0e-2, rbox = 40.0)

    basis = SelfConsistent.performSCF([Configuration("1s 2s")], ni, grid, asfSettings; printout=false).levels[1].basis
    orb_1s = basis.orbitals[Subshell("1s_1/2")];   orb_2s = basis.orbitals[Subshell("2s_1/2")]
    primitives = Bsplines.generatePrimitives(grid)
    potential  = Nuclear.nuclearPotential(ni, grid)
    omega = 0.3;   target = Subshell("2p_1/2")
    activeOrbitals  = Dict(Subshell("1s_1/2") => orb_1s, Subshell("2s_1/2") => orb_2s)
    targetSubshells = Dict(Subshell("1s_1/2") => target, Subshell("2s_1/2") => target)

    # Part 1: the rigorous, decisive check.
    oneActive  = Dict(Subshell("2s_1/2") => target)
    settings1U = RandomPhase.Settings(RandomPhase.Settings(); activeSubshells=[Subshell("2s_1/2")], targetSubshells=oneActive)
    settings1C = RandomPhase.Settings(settings1U; couplingEnabled=true, maxIterations=5, accuracy=1.0e-10)
    r1U = RandomPhase.solve(E1, Basics.Coulomb, omega, activeOrbitals, potential, primitives; settings=settings1U)
    r1C = RandomPhase.solve(E1, Basics.Coulomb, omega, activeOrbitals, potential, primitives; settings=settings1C)
    c1u = r1U.perturbedOrbitals[Subshell("2s_1/2")];   c1c = r1C.perturbedOrbitals[Subshell("2s_1/2")]
    println(">> ONE active electron: couplingEnabled=true converged = $(r1C.converged) in $(r1C.nIterations) " *
            "iteration(s); max|P_coupled - P_uncoupled| = $(maximum(abs.(c1u.P - c1c.P)))   (must be exactly 0.0)")

    # Part 2: the honest, two-electron toy case.
    settings2 = RandomPhase.Settings(RandomPhase.Settings();
                    activeSubshells = [Subshell("1s_1/2"), Subshell("2s_1/2")], targetSubshells = targetSubshells,
                    couplingEnabled = true, maxIterations = 15, printBefore = true)
    respC = RandomPhase.solve(E1, Basics.Coulomb, omega, activeOrbitals, potential, primitives; settings=settings2)
    println(">> TWO active electrons (toy case): converged = $(respC.converged) after $(respC.nIterations) iterations " *
            "-- see this branch's own comment for why that is not yet expected to be true.")
end
