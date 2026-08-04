
println("Fa) Stepwise decay cascades: a fast Mg K-hole reference case, and the larger Si^- 1s-3p case.")

using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


# NOTE on the cascade APPROACH used below (04-Aug-2026). Cascade.AverageSCA() again means what it is defined
# to mean -- every level a single CSF, no configuration mixing of any kind (Fritzsche et al., Symmetry 13, 520
# (2021) Sect. 3.3(i); Eur. Phys. J. D 78, 75 (2024) Sect. 2.3(a)). Between commits 7cc164b and 5893920 it had
# silently acquired configuration mixing within each block, i.e. it behaved as SCA in that one respect while
# keeping AverageSCA's continuum orbitals; any cascade result produced in that window is a hybrid of the two
# approaches and should be re-run. Cascade.RefinedSCA() is now also defined (block CI + per-charge-state
# orbitals + continuum orbitals resolved per fine-structure transition) but its continuum half is not yet
# implemented, so it currently behaves as SCA.


if  false
    # Last visit:      04-Aug-2026
    # Last successful: unknown ... not yet verified against literature
    #
    # Branch a: FAST reference case -- the K-shell hole of Mg^+ (1s^1 2s^2 2p^6 3s^2, 11 electrons), decaying
    #   by a SINGLE electron emission (maxElectronLoss = 1). This is deliberately the same physical system as
    #   the worked example of Symmetry 13, 520 (2021), Sect. 4 ("the stepwise decay cascade of atomic
    #   magnesium, following a 1s inner-shell ionization"), so that the branch doubles as a check against the
    #   published account of the method.
    #
    #   Why this one is fast, and branch b is not: apart from the 1s hole every shell here is closed, so each
    #   charge state contributes only a handful of configurations and each configuration only a few levels.
    #   Restricting to maxElectronLoss = 1 keeps the cascade to the first Auger/radiative generation, which is
    #   where the K-hole physics (the KLL Auger group and the K-alpha line) actually lives. Branch b, by
    #   contrast, has an open 3p^4 valence shell AND maxElectronLoss = 2, which generated 247 cascade steps and
    #   ran for hours -- fine as an occasional benchmark, useless as a working example to iterate on.
    # REPORT (04-Aug-2026): runs in 40 s / 1.4 GB, against ~8-10 h for branch b -- a factor ~800, which is what
    # makes branch-by-branch work possible at all. 10 cascade blocks, 30 steps, 44 radiative and 24 Auger lines.
    #   * K-alpha comes out at 1247-1248 eV against the known Mg K-alpha_1,2 of 1253.6 eV, i.e. 0.5% low --
    #     about what a single-CSF Dirac-Fock treatment without correlation or QED should give -- and it comes
    #     out as a DOUBLET (2 lines), correctly reflecting the 2p_1/2 / 2p_3/2 final states.
    #   * The Auger groups are ordered correctly and consistently ~0.5% low: KL1L1 1083, KL1L2,3 1123-1133,
    #     KL2,3L2,3 1163-1170, KL1M1 1187, KL2,3M1 1231, KM1M1 1271 eV. The ordering is the physical one --
    #     the further out the two final holes sit, the less binding energy is released to the electron.
    #   * Not yet dated "Last successful": the K-alpha agreement is a genuine literature check, but the Auger
    #     energies have only been checked for ordering and internal consistency, not against measured KLL
    #     values, and no rate or yield has been compared with anything.
    #
    # EFFICIENCY FINDING worth acting on: steps 22 and 24 decay the 1s hole to the 2s-hole and 3s-hole
    # configurations, i.e. they are one-electron 2s -> 1s and 3s -> 1s transitions. Both are E1-forbidden
    # (Delta l = 0), and both duly produce 0 lines -- after having been generated, set up and computed. That is
    # 2 of 30 steps here, and the fraction can only grow for cascades with more shells. Symmetry 13, 520 (2021)
    # Sect. 2.3 says precisely that "only steps with at least one non-zero transition amplitude need to be
    # considered ... and that can be read-off from just the occupation of the underlying configurations and,
    # perhaps, by using further selection rules" -- so the filter is part of the concept but is not applied.
    setDefaults("print summary: open", "zzz-Cascade-Mg-Khole.sum")

    name   = "Stepwise decay after Mg 1s ionization"
    # decayShells is passed as Shell[] on purpose: it is CONFIRMED non-functional for StepwiseDecayScheme
    # (Cascade's ForStepwiseDecay path never consults it, and the FromBasis configuration extraction drops
    # zero-occupancy shells), so supplying a list here would only suggest a restriction that is not applied.
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0,
                                         Shell[], Shell[], Shell[])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fa.dat")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      04-Aug-2026
    # Last successful: unknown ... never verified; see the cost note below
    #
    # Branch b: 1s INNER-SHELL IONIZED Si^-, decaying with up to TWO emitted electrons.
    #
    #   Electron bookkeeping, which the previous version of this branch got wrong: Si^- carries 15 electrons
    #   (Z = 14 plus the extra one), so REMOVING a 1s electron leaves 14 -- configuration
    #   1s^1 2s^2 2p^6 3s^2 3p^3. The branch used to carry 1s^1 2s^2 2p^6 3s^2 3p^4, which is 15 electrons,
    #   i.e. Si^- with the 1s electron merely EXCITED into 3p rather than removed; that is a different physical
    #   scenario, and it contradicted the file's own opening line, which announced "Si^+".
    #
    #   With the K hole and an open 3p^3 valence shell this is a genuinely large cascade, unlike branch a.
    #   Should it prove too expensive, the obvious fallback is 1s inner-shell ionized NEUTRAL Si
    #   (1s^1 2s^2 2p^6 3s^2 3p^2, 13 electrons), which removes one more open-shell electron.
    setDefaults("print summary: open", "zzz-Cascade-Si-1s-ionized.sum")

    name   = "Stepwise decay after Si- 1s ionization"
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 2, Dict{Int64,Float64}(), 0,
                                         Shell[], Shell[], Shell[])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2 3p^3")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fa.dat")
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last visit:      04-Aug-2026
    # Last successful: unknown ... not yet run
    #
    # Branch c: the SIMULATION half -- take the amplitudes computed by branch a (or b) and derive the final
    #   ion distribution from them. This is the separation the concept papers insist on: the expensive
    #   computation produces lists of lines once, and any number of cheap simulations post-process them
    #   (Eur. Phys. J. D 78, 75 (2024), Sect. 2.4).
    #
    #   REPAIRED 04-Aug-2026: this branch could not run at all. It called JLD.load while the file imports
    #   JLD2 (there is no JLD module here), and it pointed at
    #   "zzz-cascade-decay-computations-2020-01-25T21.jld" -- a hard-coded artefact from January 2020 that has
    #   long since disappeared. The filename is now a variable that has to be set to the file branch a
    #   actually wrote; a cascade computation prints that name at the end of its run, in the form
    #   zzz-cascade-decay-computations-<YYYY-MM-DDTHH>.jld.
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

    dataFilename = "example-Fa.dat/zzz-cascade-decay-computations-2026-08-04T16.jld"   ## <-- from branch b
    if  !isfile(dataFilename)
        error("Branch c needs the .jld file written by branch a (or b); set dataFilename to the name that " *
              "the cascade computation printed at the end of its run. Looked for: $dataFilename")
    end
    data = [JLD2.load(dataFilename)]
    name = "Simulation of the Si- K-hole decay"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name,
                              properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution()],
                              settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 1.0)]),
                              computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
