
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
    # Last visit:      06-Aug-2026
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
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch b: ION DISTRIBUTION -- how the population is shared between charge states once the cascade has
    #   run its course. The total distributed probability is a real conservation check on the propagation and
    #   should come out as 1.00000.
    setDefaults("print summary: open", "zzz-Cascade-Fb-ionDistribution.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
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
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch c: FINAL LEVEL DISTRIBUTION -- the same propagation as branch b, but resolved into the individual
    #   levels the cascade ends in rather than summed over charge states. Its total must agree with branch b's
    #   1.00000, and that agreement is the cross-check worth making here: the two properties share the
    #   propagation and differ only in how the result is collected.
    setDefaults("print summary: open", "zzz-Cascade-Fb-finalLevels.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
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
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch d: PHOTON INTENSITIES -- the fluorescence spectrum of the cascade, i.e. the photon lines weighted
    #   by how much population actually passes through their initial levels. This branch has an external
    #   anchor: K-alpha of Mg is measured at 1253.6 eV, and example-Fa.jl branch a puts the computed line at
    #   1247-1248 eV under AverageSCA, so the strongest feature here should sit there. The window is opened
    #   wide (0 ... 2000 eV) so that the L-shell lines at a few tens of eV show up alongside it.
    setDefaults("print summary: open", "zzz-Cascade-Fb-photonIntensities.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
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
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch e: RELAXATION CURVE -- the times by which 70%, 80% and 90% of the initial population has reached
    #   the ground configuration. Renamed from MeanRelaxationTime on 05-Aug-2026: the property never computed
    #   a mean, and a mean over a cascade is not an observable in any case, since the pathways take widely
    #   different times and the distribution is broad and far from exponential. What it computes -- and what
    #   Eur. Phys. J. D 78, 75 (2024) actually describes -- is a set of percentiles of that distribution, i.e.
    #   a curve, and a curve IS comparable to a pump-probe measurement of a final-state yield against delay.
    #
    #   The property was designed for HOLLOW IONS, where "how long until the ion is back in its ground
    #   configuration" is the natural question; a single K-shell hole is a much simpler first test.
    #   Expect femtoseconds: a Mg K hole decays predominantly by Auger emission at rates of order 1e14 /s, so
    #   an answer in seconds would indicate a units problem rather than physics. Note the times are reported in
    #   ATOMIC UNITS of time (1 a.u. = 24.19 as), as is timeStep. groundConfigs names the configuration the
    #   relaxation is measured TOWARDS -- and it must name the ground configuration of EVERY charge state the
    #   cascade can end in, here 11, 10 and 9 electrons. Listing only the 11-electron one makes the target
    #   unreachable: branch b shows just 1.55e-04 of the population ends there, so the propagation can never
    #   reach the 91% it needs and runs without bound. That mistake is easy to make and the code now reports
    #   it rather than spinning.
    setDefaults("print summary: open", "zzz-Cascade-Fb-relaxationTime.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: mean relaxation time",
                            property=Cascade.RelaxationCurve(1.0e-3, [(1, 1.0)], Configuration[],
                                                             [Configuration("1s^2 2s^2 2p^6 3s^1"),
                                                              Configuration("1s^2 2s^2 2p^6"),
                                                              Configuration("1s^2 2s^2 2p^5")]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch f: ELECTRON INTENSITIES -- the Auger-electron spectrum. This is the natural partner of branch d's
    #   photon spectrum and arguably the more informative of the two for a K-hole cascade, since Mg decays
    #   overwhelmingly by electron emission rather than radiatively.
    #
    #   WAS BLOCKED, RESOLVED 05-Aug-2026. Cascade.ElectronIntensities was declared, had a constructor and a
    #   Base.show, and was a perfectly ordinary AbstractSimulationProperty -- but Cascade.perform(simulation)
    #   never dispatched That chain of elseif branches handles PhotoAbsorptionSpectrum, IonDistribution,
    #   FinalLevelDistribution, PhotonIntensities, DrRateCoefficients, RrRateCoefficients, RelaxationCurve,
    #   ExpansionOpacities and RosselandOpacities, and falls through silently for this one; there is likewise
    #   no simulateElectronIntensities beside simulatePhotonIntensities. It is the same reserved-but-dead
    #   pattern already found in Plasma's aiSettings and in Continuum.Settings.includeExchange -- a type that
    #   advertises a capability which does not exist. The propagation had in fact collected the Auger
    #   emissions all along; only the simulate function and the dispatch entry were missing, and both are now
    #   in place. It was the dispatch refactor that made this visible: under the former chain of
    #   `typeof(property) == ...` tests the request simply fell through and returned nothing, silently.
    setDefaults("print summary: open", "zzz-Cascade-Fb-electronIntensities.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
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
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: unknown ... not yet verified
    #
    # Branch g: AUGER-AUGER COINCIDENCE -- the spectrum of the SECOND Auger electron, recorded only for those
    #   cascade events in which the FIRST Auger electron fell into the KLL window. This is the classic cascade
    #   coincidence: gating on the K-shell Auger selects which intermediate (two-hole) state was populated, and
    #   the second spectrum then shows how that state decayed further.
    #
    #   Gates are matched in EMISSION ORDER along each pathway, and all energies are in ATOMIC UNITS, as
    #   everywhere else in these properties -- hence the explicit conversions below, which also keep the
    #   intent readable. The Mg KLL group lies near 1080-1240 eV (see example-Fa.jl branch a), so the gate is
    #   set around it; the second electron is an L-shell Auger of a few tens of eV.
    #
    #   Consistency check worth making: with gates = [] and a wide observation window, this property must
    #   reproduce branch f's ungated electron spectrum exactly.
    setDefaults("print summary: open", "zzz-Cascade-Fb-augerAuger.sum")

    dataFilename = "example-Fb.dat/zzz-cascade-decay-computations-2026-08-06T07.jld"
    if  !isfile(dataFilename)   error("Run branch a first and paste its printed filename here: $dataFilename")   end
    data = [JLD2.load(dataFilename)]

    gateLo = Defaults.convertUnits("energy: to atomic", 1000.);  gateHi = Defaults.convertUnits("energy: to atomic", 1300.)
    obsLo  = Defaults.convertUnits("energy: to atomic",    0.);  obsHi  = Defaults.convertUnits("energy: to atomic",  200.)

    wc = Cascade.Simulation(Cascade.Simulation(), name="Mg K-hole: Auger-Auger coincidence",
                            property=Cascade.ParticleCoincidences([Cascade.ElectronGate(gateLo, gateHi)],
                                                                  Cascade.ElectronGate(obsLo, obsHi), [(1, 1.0)]),
                            method=Cascade.ProbPropagation(),
                            settings=Cascade.SimulationSettings(false, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
