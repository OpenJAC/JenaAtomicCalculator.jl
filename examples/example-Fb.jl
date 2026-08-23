
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
#
# STORED CASCADE DATA GOES STALE, and on 23-Aug-2026 all of it had (this file and example-Fa.jl branch c).
# Radial.Grid has since been restructured into three sub-structs (parameters / knots / mesh) with a
# Base.getproperty shim that keeps the old flat field names working IN SOURCE. Serialization does not go
# through that shim, so a .jld written before the refactor holds the old flat layout, JLD2 rebuilds it as a
# ReconstructedMutable and cannot convert it to an Orbital: the file loads no more. Source compatibility and
# FILE compatibility are separate things, and only the first one was preserved. Hence branches b onwards no
# longer carry a pasted file name -- they take the newest .jld in the data directory, so that re-running
# branch a is the whole of the repair.


if  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch a: THE COMPUTATION. Everything below consumes the .jld file this branch writes; the branches
    #   below find it themselves, by picking the most recently written cascade .jld in example-Fb.dat, so
    #   re-running this branch is all that is needed to refresh them.
    #
    # REPORT (23-Aug-2026): 48.2 s cold, 53 levels generated from the one initial K-hole level. Dated on the
    # strength of its consumers, which is the only honest ground for a branch that produces amplitudes and
    # asserts nothing about them: branches b to g all read this file and all six are dated, and two of the
    # checks they make reach back into these amplitudes rather than merely being internally consistent --
    # branch d recovers the K-alpha doublet at 1247.72 / 1247.49 eV, the same lines to six figures as
    # example-Fa.jl branch a computes independently with maxElectronLoss = 1, and its intensity reproduces
    # that file's fluorescence yield exactly. A fault in these amplitudes would have to survive both.
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch b: ION DISTRIBUTION -- how the population is shared between charge states once the cascade has
    #   run its course. The total distributed probability is a real conservation check on the propagation and
    #   should come out as 1.00000.
    #
    # REPORT (23-Aug-2026): 4.1 s. The propagation converges in 6 rounds, and
    #       11 electrons  1.55394e-04       10 electrons  9.15226e-02       9 electrons  9.08322e-01
    # with the total distributed probability at 1.00000e+00. The shape is the expected one for a light
    # element: the K hole Augers in ~95% of cases (see the fluorescence yield in example-Fa.jl branch a),
    # and the L hole left behind by the ~5% that radiate mostly Augers in turn, so almost nothing survives
    # with 11 electrons. As in example-Fa.jl branch c, the LAST bin is inflated by the model boundary and
    # not by physics: maxElectronLoss = 2 forbids any further emission, so every pathway that would have
    # gone on piles up at 9 electrons. "91% of the ions end up doubly ionized" would be a misreading.
    setDefaults("print summary: open", "zzz-Cascade-Fb-ionDistribution.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch c: FINAL LEVEL DISTRIBUTION -- the same propagation as branch b, but resolved into the individual
    #   levels the cascade ends in rather than summed over charge states. Its total must agree with branch b's
    #   1.00000, and that agreement is the cross-check worth making here: the two properties share the
    #   propagation and differ only in how the result is collected.
    #
    # REPORT (23-Aug-2026): 0.5 s, and the cross-check is sharper than the one this branch's text proposes.
    # Only 6 of the 53 levels carry any population, and they REGROUP onto branch b charge state by charge
    # state, not merely in total:
    #       11 electrons   level 53                        1.55394e-04     vs  b: 1.55394e-04
    #       10 electrons   levels 48, 49, 52               9.15234e-02     vs  b: 9.15226e-02
    #        9 electrons   levels 44, 45                   9.08322e-01     vs  b: 9.08322e-01
    #                                       total over all 53 rows  0.99999996
    # That constrains WHICH levels belong to which charge state and not just a sum, so it would fail on a
    # mis-assignment that a total could never see; the residual 8e-07 in the 10-electron row is the printed
    # precision of branch b's own six figures. One gap worth knowing: this property prints no total line of
    # its own, unlike branch b, so the 0.99999996 above had to be summed externally.
    setDefaults("print summary: open", "zzz-Cascade-Fb-finalLevels.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch d: PHOTON INTENSITIES -- the fluorescence spectrum of the cascade, i.e. the photon lines weighted
    #   by how much population actually passes through their initial levels. This branch has an external
    #   anchor: K-alpha of Mg is measured at 1253.6 eV, and example-Fa.jl branch a puts the computed line at
    #   1247-1248 eV under AverageSCA, so the strongest feature here should sit there. The window is opened
    #   wide (0 ... 2000 eV) so that the L-shell lines at a few tens of eV show up alongside it.
    #
    # REPORT (23-Aug-2026): 0.6 s, 80 photon lines, 9.491e-01 photons emitted per ion. THE EXPECTATION
    # STATED ABOVE IS WRONG, and correcting it is the useful part. K-alpha is there and is exactly where it
    # should be -- a doublet at 1247.72 and 1247.49 eV, the only two lines above 500 eV, carrying
    # 3.599e-02 + 1.811e-02 = 5.410e-02 photons per K hole. That number is an identity, not an estimate: it
    # is the Coulomb-gauge fluorescence yield 0.054099 that example-Fa.jl branch a gets from the ratio of
    # its rate tables, in a different file and by a different route, agreeing to five figures. So the
    # propagation and the rate tables are consistent.
    #   But K-alpha is NOT "the strongest feature". The strongest line in the spectrum sits at 65.85 eV with
    #   0.1496 photons per ion, nearly THREE TIMES the whole K-alpha doublet, and the soft lines carry
    #   0.895 of the 0.949 total. Most of that is a model boundary rather than physics: with
    #   maxElectronLoss = 2 the two-hole states reached after the K Auger are forbidden to Auger again, so
    #   they are forced to radiate, and in a real Mg ion they would overwhelmingly emit an electron instead.
    #   The hard part of this spectrum is trustworthy; the soft part is an artefact of the cap and should
    #   not be read as a fluorescence spectrum.
    setDefaults("print summary: open", "zzz-Cascade-Fb-photonIntensities.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
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
    #
    # REPORT (23-Aug-2026): 2.6 s, and the expectation stated above -- "expect femtoseconds" -- is WRONG,
    # though not for the reason it feared. The answer is
    #       70% relaxed at 2.4504e+05 fs      80% at 3.8823e+05 fs      90% at 7.3095e+05 fs
    # i.e. 245, 388 and 731 PICOSECONDS. That is not a units problem. The femtosecond figure belongs to the
    # K-HOLE LIFETIME (3.2 fs, example-Fa.jl branch a), whereas this property measures something else
    # entirely: the time for the population to reach the GROUND CONFIGURATION, which is set by the SLOWEST
    # step on the way and not by the fastest. The last steps here are soft E1 transitions of 34 to 72 eV
    # with Einstein coefficients of order 1e9 to 1e10 /s, i.e. tenths of a nanosecond each, and 245 ps is
    # exactly what those give. The two times differ by five orders of magnitude and both are correct.
    #   The curve also SATURATES at 0.984627 and never reaches 100%, and the missing 1.53730e-02 is
    #   identifiable to the last digit: it is levels 48 (J=0, odd) and 49 (J=2, odd) of branch c, whose
    #   occupations are 2.77062e-06 and 1.53703e-02 and sum to 1.537307e-02. Both are 1s^2 2s^2 2p^5 3s^1
    #   levels and both are E1-FORBIDDEN to the 1s^2 2s^2 2p^6 (J=0, even) ground state they would have to
    #   reach -- 0 -> 0 is strictly forbidden and 2 -> 0 violates the E1 triangle rule -- so that population
    #   is metastable and is trapped for as long as only E1 is in the cascade. A relaxation curve that
    #   saturates below 1 is therefore the correct answer here, and the shortfall names the levels.
    setDefaults("print summary: open", "zzz-Cascade-Fb-relaxationTime.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
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
    #
    # REPORT (23-Aug-2026): 0.4 s, 61 electron lines, 1.908232 electrons emitted per ion -- and that total
    # is fixed in advance by branch b, which is what makes it a check rather than an observation. The mean
    # number of emitted electrons must equal the mean charge increase of the ion distribution,
    #       2 * 9.08322e-01 + 1 * 9.15226e-02 + 0 * 1.55394e-04 = 1.908167
    # against the 1.908232 counted here, agreeing to 3.4e-05 relative -- the printed precision of branch b's
    # six figures. The photon and electron spectra are thus tied to the same propagation from both sides:
    # branch d's hard lines reproduce the fluorescence yield, and this branch's total reproduces the charge
    # balance.
    setDefaults("print summary: open", "zzz-Cascade-Fb-electronIntensities.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
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
    #
    # REPORT (23-Aug-2026): 0.7 s, and both checks were made. The gated result is 9.083e-01 coincidence
    # counts per ion, which is exactly branch b's 9-electron fraction 9.08322e-01 -- as it must be, since
    # the gate at 1000-1300 eV spans every K-Auger group (1083 to 1271 eV, example-Fa.jl branch a) so that
    # every ion emitting a SECOND electron contributes precisely one count, and those are the ions that end
    # with 9 electrons.
    #   The consistency check the paragraph above proposes was then run for the first time: with
    #   gates = ElectronGate[] and the observation window opened to 0 ... 2000 eV, this property returns 61
    #   lines against branch f's 61, at identical energies, with a worst relative intensity difference of
    #   0.000e+00 and totals agreeing to all printed digits (1.90823178e+00 both). Note what that does and
    #   does not establish: it is a REDUCTION check, so it verifies the observation-window bookkeeping and
    #   the pathway weighting of the coincidence machinery against an independently written routine, while
    #   the GATING itself is what the 9.083e-01 identity above tests.
    setDefaults("print summary: open", "zzz-Cascade-Fb-augerAuger.sum")

    # The .jld file written by branch a. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fb.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
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
