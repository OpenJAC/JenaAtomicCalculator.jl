
println("Fb) Cascade SIMULATIONS: one computation, several properties derived from the same cascade data.")

using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


# REDESIGNED 05-Aug-2026. The previous contents could not run at all: it used Cascade.PhotonIonizationScheme
# (no such type -- it is PhotoIonizationScheme), called StepwiseDecayScheme with five arguments where it takes
# seven, loaded hard-coded .jld files from February 2020, and built its simulations with `properties=` (now a
# single `property=`) and a six-field SimulationSettings (now three). Nothing of it was kept.
#
# The new layout follows the separation the concept papers insist upon (Fritzsche et al., Eur. Phys. J. D 78,
# 75 (2024), Sect. 2.4): the expensive COMPUTATION generates the transition amplitudes once, and any number of
# cheap SIMULATIONS post-process that same data into the properties an experiment might actually report.
# Where example-Fa.jl exercises the computation side (one system, three cascade approaches), this file
# exercises the simulation side: branch a computes the cascade once and writes it to disk, and branches b-f
# each derive a DIFFERENT property from exactly that file. Branches b onwards therefore run in seconds, and
# any difference between them is a property of the simulation and never of the underlying atomic data.
#
# The system is the K-shell hole of Mg^+ again, as in example-Fa.jl branch a and in the worked example of
# Symmetry 13, 520 (2021), Sect. 4 -- but with maxElectronLoss = 2 rather than 1, so the decay tree is deep
# enough for the spectra and the ion distribution to be interesting rather than trivial.


if  false
    # Last visit:      unknown ... not yet run
    # Last successful: unknown ... not yet verified
    #
    # Branch a: THE COMPUTATION. Everything below consumes the .jld file this branch writes; its name is
    #   printed at the end of the run and has to be pasted into dataFilename in the branches below.
    setDefaults("print summary: open", "zzz-Cascade-Fb-computation.sum")

    name   = "Mg 1s-hole cascade, two electrons lost"
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 2, Dict{Int64,Float64}(), 0,
                                         Shell[], Shell[], Shell[])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fb.dat")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      unknown ... not yet run
    # Last successful: unknown ... not yet verified
    #
    # Branch b: ION DISTRIBUTION -- how the population is shared between charge states once the cascade has
    #   run its course. The total distributed probability is a real conservation check on the propagation and
    #   should come out as 1.00000.
    setDefaults("print summary: open", "zzz-Cascade-Fb-ionDistribution.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-05T16.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: ion distribution",
                            property=Cascade.IonDistribution([(1, 1.0)], Configuration[]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      unknown ... not yet run
    # Last successful: unknown ... not yet verified
    #
    # Branch c: FINAL LEVEL DISTRIBUTION -- the same propagation as branch b, but resolved into the individual
    #   levels the cascade ends in rather than summed over charge states. Its total must agree with branch b's
    #   1.00000, and that agreement is the cross-check worth making here: the two properties share the
    #   propagation and differ only in how the result is collected.
    setDefaults("print summary: open", "zzz-Cascade-Fb-finalLevels.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-05T16.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: final level distribution",
                            property=Cascade.FinalLevelDistribution([(1, 1.0)], Configuration[], Configuration[]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      unknown ... not yet run
    # Last successful: unknown ... not yet verified
    #
    # Branch d: PHOTON INTENSITIES -- the fluorescence spectrum of the cascade, i.e. the photon lines weighted
    #   by how much population actually passes through their initial levels. This branch has an external
    #   anchor: K-alpha of Mg is measured at 1253.6 eV, and example-Fa.jl branch a puts the computed line at
    #   1247-1248 eV under AverageSCA, so the strongest feature here should sit there. The window is opened
    #   wide (0 ... 2000 eV) so that the L-shell lines at a few tens of eV show up alongside it.
    setDefaults("print summary: open", "zzz-Cascade-Fb-photonIntensities.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-05T16.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: photon intensities",
                            property=Cascade.PhotonIntensities(0., 2000., [(1, 1.0)], Configuration[]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      unknown ... not yet run
    # Last successful: unknown ... not yet verified
    #
    # Branch e: MEAN RELAXATION TIME -- how long the ion takes to reach its ground configuration. This
    #   property was originally designed for HOLLOW IONS, where "how long until the ion is back in its ground
    #   configuration" is the natural question, and it has hardly been used in practice since; a single
    #   K-shell hole is a much simpler case for it, which makes it a reasonable first test. Expect
    #   femtoseconds: a Mg K hole decays predominantly by Auger emission at rates of order 1e14 /s, so an
    #   answer in seconds, or one that is plainly in atomic units, would indicate a units problem rather than
    #   physics. groundConfigs names the configuration the relaxation is measured TOWARDS.
    setDefaults("print summary: open", "zzz-Cascade-Fb-relaxationTime.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-05T16.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: mean relaxation time",
                            property=Cascade.MeanRelaxationTime(1.0e-3, [(1, 1.0)], Configuration[],
                                                                [Configuration("1s^2 2s^2 2p^6 3s^1")]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last visit:      unknown ... NOT YET RUNNABLE, see below
    # Last successful: unknown ... not yet verified
    #
    # Branch f: ELECTRON INTENSITIES -- the Auger-electron spectrum. This is the natural partner of branch d's
    #   photon spectrum and arguably the more informative of the two for a K-hole cascade, since Mg decays
    #   overwhelmingly by electron emission rather than radiatively.
    #
    #   BLOCKED (05-Aug-2026): Cascade.ElectronIntensities is declared, has a constructor and a Base.show, and
    #   is a perfectly ordinary AbstractSimulationProperty -- but Cascade.perform(simulation) never dispatches
    #   on it. That chain of elseif branches handles PhotoAbsorptionSpectrum, IonDistribution,
    #   FinalLevelDistribution, PhotonIntensities, DrRateCoefficients, RrRateCoefficients, MeanRelaxationTime,
    #   ExpansionOpacities and RosselandOpacities, and falls through silently for this one; there is likewise
    #   no simulateElectronIntensities beside simulatePhotonIntensities. It is the same reserved-but-dead
    #   pattern already found in Plasma's aiSettings and in Continuum.Settings.includeExchange -- a type that
    #   advertises a capability which does not exist.
    #   The data it needs is already there: the AutoIonization lines sit in the .jld that branch a writes, so
    #   implementing it means mirroring simulatePhotonIntensities onto those lines.
    setDefaults("print summary: open", "zzz-Cascade-Fb-electronIntensities.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-05T16.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: electron intensities",
                            property=Cascade.ElectronIntensities(0., 2000., [(1, 1.0)]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
