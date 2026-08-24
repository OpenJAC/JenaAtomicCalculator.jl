
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
#
# STORED CASCADE DATA GOES STALE. Branch c reads the .jld that branch a or b wrote, and on 23-Aug-2026 the
# file it had been pointed at (from 04-Aug) could no longer be loaded at all: Radial.Grid has since been
# restructured into three sub-structs with a Base.getproperty shim, which keeps the old field names working
# in SOURCE but not in a serialized file, so JLD2 rebuilds the stored grid as a ReconstructedMutable and
# cannot convert it to an Orbital. Branch c therefore no longer carries a pasted file name; it takes the
# newest .jld in example-Fa.dat, and re-running branch a or b is the whole of the repair.


if  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
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
    #
    # REPORT (23-Aug-2026): the rate check that was missing above has now been made, and it is the one thing
    # here that does NOT come out well. Warm cost 9.1 s; the 40 s above was a COLD first run of a session,
    # which pays the compilation of the whole cascade path and is not comparable with a warm number.
    #   Summing the rates OUT OF the K-hole block gives the K fluorescence yield and the K level width
    #   directly, and both of those are measured quantities:
    #       A_r(K-alpha) = 1.545e13 /s (Babushkin), 1.676e13 /s (Coulomb)
    #       A_a(total)   = 2.931e14 /s over all six Auger groups
    #       omega_K      = 0.050 (Bab), 0.054 (Coul)     against a measured 0.030
    #       Gamma_K      = 0.203 eV, tau = 3.24 fs       against a measured 0.36 eV, 1.83 fs
    #   * The RADIATIVE side is right: Gamma_r = 0.0102 eV against omega_K * Gamma_K = 0.030 * 0.36 =
    #     0.0108 eV, i.e. 6% low. The E1 amplitude and the K-alpha energy are both sound.
    #   * The AUGER side is short by a factor ~1.8, and that one deficit is the whole of the yield error.
    #     The channel set is complete, and the INTERNAL pattern of KL2,3L2,3 is right: its two J=2 lines
    #     carry 6.2e13 while the J=0 and J=1 lines are down by 10^3 and 10^5, which is the 1D dominance one
    #     expects. What is wrong is the SIZE of that group. The computed KLL split is
    #     KL1L1 : KL1L2,3 : KL2,3L2,3 = 0.27 : 0.50 : 0.23, whereas the measured low-Z pattern has
    #     KL2,3L2,3 DOMINANT at about 0.6. Giving that group its experimental share would also close the
    #     factor-1.8 gap in the total, so both discrepancies are consistent with a single cause there.
    #   * It is NOT the continuum normalization. Re-running with "pure Coulomb" in place of "pure sine"
    #     moves the total Auger rate from 2.9308e14 to 2.9339e14 /s -- 0.1%. That candidate is excluded.
    #   * It is NOT the orbital relaxation either. Branches d and e climb the approach ladder and the yield
    #     gets WORSE, not better: 0.050 -> 0.062 -> 0.065 against the measured 0.030.
    #   TWO SUSPECTS REMAIN, and the first version of this note named only one.
    #     (i) MISSING CORRELATION. All three tiers carry at most five CSFs per block and so almost no
    #         configuration mixing -- exactly what a CI-enabled approach would test.
    #     (ii) THE TWO-PARTICLE ANGULAR COEFFICIENTS THEMSELVES. An Auger matrix element moves two bound
    #         electrons, so it is built from spin-angular coefficients for CSF pairs differing in two shells.
    #         Note the shape of what is measured: the ONE-particle quantity here, the radiative width, is right
    #         to 6%, while the TWO-particle one is low by 1.8x with its internal branching pattern wrong.
    #     WHAT HAS BEEN CHECKED, 24-Aug-2026, and what that check cannot do. A reference-free structural test --
    #     for two CSFs of equal J the orbit invariant of the two-particle coefficients must be identical for
    #     (A,B) and (B,A) -- was run over SpinAngular across eight configurations: 9050 equal-occupation, 18580
    #     one-electron-move and 7685 TWO-OR-MORE-ELECTRON-MOVE invariants, zero asymmetric, worst deviation
    #     0.0e+00. That removes one class of defect and NOTHING MORE. Hermiticity is necessary and not
    #     sufficient: it sees only errors that are asymmetric under exchanging bra and ket, and is blind to a
    #     uniform scale error on a family of coefficients or to a family missing from both directions -- which
    #     are precisely the two shapes that would produce what is measured here. "Not asymmetric" is a much
    #     weaker statement than "verified" and must not be recorded as one.
    #     Separating (i) from (ii) needs an INDEPENDENT RECONSTRUCTION of those coefficients, which does not
    #     exist for the two-electron-move case. Until it does, neither suspect is excluded.
    #   Dated on the energies, which are verified against the literature. A YIELD OR A LIFETIME MUST NOT BE
    #   QUOTED FROM THIS BRANCH: omega_K is about twice the measured value and tau is 1.8x too long.
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
    # Last visit:      23-Aug-2026
    # Last successful: unknown ... deliberately not re-run; see the cost note below
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
    #
    # REPORT (04-Aug-2026): 19 minutes -- the fallback proved unnecessary. 233 steps, 12377 radiative and
    # 31806 Auger lines, 946 levels. That is a factor ~30 below the ~8-10 h estimated for the version this
    # branch used to carry, from two causes: the corrected configuration has one electron fewer in the open
    # 3p shell (3p^3 rather than 3p^4), which cuts the number of CSFs sharply, and AverageSCA is again the
    # single-CSF approximation it is defined to be rather than the CI hybrid it had drifted into.
    # The .jld written by this branch feeds branch c.
    #
    # NOT RE-RUN on 23-Aug-2026, deliberately. Two reasons, and the second is the binding one. Its stored
    # output from 04-Aug can no longer be loaded at all (the Radial.Grid refactor; see the head of this
    # file), so the 19 minutes would have to be paid again from scratch. And 233 steps over 946 levels with
    # 44000 lines needs several GB, which this machine does not have free -- a much smaller run was
    # OOM-killed the same afternoon while writing its results. Branch c was therefore run against branch a's
    # data instead, which is what its own text already allows. This branch remains the file's occasional
    # benchmark, and it needs a machine with room.
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
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch c: the SIMULATION half -- take the amplitudes computed by branch a (or b) and derive the final
    #   ion distribution from them. This is the separation the concept papers insist on: the expensive
    #   computation produces lists of lines once, and any number of cheap simulations post-process them
    #   (Eur. Phys. J. D 78, 75 (2024), Sect. 2.4).
    #
    # REPORT (04-Aug-2026): first successful visit to Cascade.Simulation. The probability propagation walks
    # 946 levels and converges cleanly in 9 rounds (residual exactly 0.0), with the total distributed
    # probability coming out at 1.00000 -- a genuine conservation check on the propagation, not just a
    # formatting artefact. Final ion distribution after the Si^- K hole decays:
    #       14 electrons  2.61e-03      (no Auger electron at all: K decayed radiatively AND so did the
    #                                    L hole it left behind -- doubly rare, hence the small number)
    #       13 electrons  8.37e-02      (one Auger electron)
    #       12 electrons  9.14e-01      (two Auger electrons)
    # This is the expected shape for a light element: the K-shell fluorescence yield of Si is only ~0.05, so
    # the K hole overwhelmingly Augers, and whatever does decay radiatively leaves an L hole that Augers in
    # turn. Note the 12-electron figure is inflated by the model boundary, not by physics: maxElectronLoss = 2
    # forbids any further emission, so every pathway that would have continued piles up in the last bin.
    # Reading it as "91% of the ions end up 2-fold ionized" would be wrong.
    #
    #   REPAIRED 04-Aug-2026: this branch could not run at all. It called JLD.load while the file imports
    #   JLD2 (there is no JLD module here), and it pointed at
    #   "zzz-cascade-decay-computations-2020-01-25T21.jld" -- a hard-coded artefact from January 2020 that has
    #   long since disappeared. The filename is now a variable that has to be set to the file branch a
    #   actually wrote; a cascade computation prints that name at the end of its run, in the form
    #   zzz-cascade-decay-computations-<YYYY-MM-DDTHH>.jld.
    #
    # REPORT (23-Aug-2026): 4.3 s. Run against branch a/d/e's Mg data rather than branch b's Si^- data,
    # since the 04-Aug file no longer loads and branch b is too large for this machine; the branch takes
    # whichever computation ran last, and prints the file it chose. On the Mg K hole with
    # maxElectronLoss = 1 (RefinedSCA orbitals) the propagation converges in 6 rounds and gives
    #       11 electrons  7.30642e-05          10 electrons  9.99927e-01      total  1.00000
    # and that first number is worth more than it looks. With only one electron allowed to leave, an ion
    # can still hold 11 electrons ONLY if the K hole radiated AND the 2p hole it left behind radiated too,
    # so the answer must factorise into two branchings that are computable from the rate tables alone.
    # It does, EXACTLY, but only when the fine structure is kept:
    #       2p-hole J=3/2 :  b_K = 0.041655   omega_L = 6.1480e-04   product 2.5609e-05
    #       2p-hole J=1/2 :  b_K = 0.020963   omega_L = 2.2637e-03   product 4.7455e-05
    #                                                          sum   7.306409e-05   vs  7.30642e-05
    # i.e. agreement to 1e-6. The product of the TOTALS gives 6.07e-05 and is simply the wrong quantity:
    # the J=1/2 hole takes only a third of the K-alpha intensity but contributes two thirds of the
    # survivors, because its Auger rate is 3.7x smaller than the J=3/2 one. So this branch checks that the
    # propagation resolves fine structure correctly and not merely that probabilities sum to one -- and it
    # also settles which gauge the simulation uses, since the identity holds in Coulomb and not Babushkin.
    #
    # The Si^- numbers recorded below are from 04-Aug and are kept for reference; they cannot be
    # reproduced without re-running branch b.
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

    # The .jld file written by branch a or b. Its name carries the run time, so it is LOOKED UP here rather
    # than pasted in: a hand-copied name goes stale silently, and a file written before the Radial.Grid
    # refactor cannot be deserialized at all -- see the note at the head of this file.
    dataDir      = "example-Fa.dat"
    dataFiles    = isdir(dataDir) ? filter(f -> startswith(f, "zzz-cascade-decay-computations"), readdir(dataDir)) : String[]
    if  isempty(dataFiles)   error("Run branch a or b first; no cascade .jld file found in $dataDir")   end
    dataFilename = joinpath(dataDir, sort(dataFiles, by = f -> mtime(joinpath(dataDir, f)))[end])
    println(">> reading cascade data from  $dataFilename")
    data = [JLD2.load(dataFilename)]
    name = "Simulation of the K-hole decay computed above"

    # The API had drifted twice since this branch was written, on top of the JLD/JLD2 breakage:
    #   * Cascade.Simulation takes a SINGLE `property`, not a `properties` array;
    #   * Cascade.SimulationSettings lost two fields and now carries only
    #     (printTree, printLongTree, initialPhotonEnergy) -- the initial level occupations, which used to be
    #     its last argument, now belong to the property itself, i.e. to Cascade.IonDistribution.
    # Starting occupation: all of the population in level 1, the K-hole level that the cascade begins from.
    wc   = Cascade.Simulation(Cascade.Simulation(), name=name,
                              property=Cascade.IonDistribution([(1, 1.0)], Configuration[]),
                              method=Cascade.ProbPropagation(),
                              settings=Cascade.SimulationSettings(true, false, 0.),
                              computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
    #
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch d: THE SAME Mg K-hole system as branch a, but one tier up the approach hierarchy --
    #   Basics.SCA() instead of Cascade.AverageSCA(). Everything else is held fixed on purpose: same nucleus,
    #   same initial configuration, same maxElectronLoss, same grid. Every difference between branch a and
    #   branch d is therefore attributable to the approach and to nothing else.
    #
    #   Why this is worth a branch of its own. Until 04-Aug-2026 SCA had never been executed at all: it raised
    #   "Not yet implemented." in three of the eight cascade paths and, in the five where it does exist, it had
    #   simply never been exercised by any example or test. Branch a's K-alpha sits 0.5% below the measured
    #   1253.6 eV, so this branch asks the question the approach hierarchy is really making a claim about --
    #   does the next tier up actually move the answer TOWARDS the measurement?
    #
    #   Caveat on what is being compared. The SCA path in module-Cascade-inc-computations.jl builds its
    #   orbitals with a separate self-consistent field for EVERY block, which is neither the GlobalOrbitals of
    #   AverageSCA nor the per-charge-state set that the published SCA specifies. So this branch measures
    #   "AverageSCA vs SCA-as-implemented", not "vs SCA-as-defined", and it changes two axes at once (level
    #   representation AND bound orbitals) rather than one. That is a property of the current code, not of the
    #   comparison, and it is exactly the inconsistency the bound-orbital work is meant to remove.
    # REPORT (23-Aug-2026): THE 04-AUG REPORT BELOW NO LONGER HOLDS, and the change is real rather than a
    # mis-reading. Warm, in one session, the ladder costs 9.1 s (AverageSCA), 16.1 s (SCA), 25-30 s
    # (RefinedSCA) -- so the 1.7x factor recorded below is CONFIRMED as a ratio (16.1/9.1 = 1.77), and only
    # the absolute seconds were wrong, having been measured cold.
    #   The K-alpha, however, has moved. Today:
    #       AverageSCA  1247.72 / 1247.49 eV      -5.9 eV against the measured 1253.69 / 1253.44
    #       SCA         1264.01 / 1263.72 eV     +10.4 eV
    #       RefinedSCA  1255.28 / 1254.99 eV      +1.6 eV
    #   so the ladder is NOT monotonic, and SCA now OVERSHOOTS by more than AverageSCA undershoots. The
    #   block energies say exactly what changed:
    #       block 10, the initial 1s hole   : -4129 (AverageSCA)   -4120 (SCA)   -4129 (RefinedSCA)
    #       block  2, the 2p-hole final st. : -5376 (AverageSCA)   -5384 (SCA)   -5384 (RefinedSCA)
    #   The 04-Aug note below states that the initial block is "-4129 eV in BOTH branches, unchanged", and
    #   that is no longer true of SCA: it now moves that block up by 9 eV on top of the 8 eV the final state
    #   relaxes, and the two add rather than cancel. RefinedSCA is what leaves the initial block alone and
    #   relaxes only the final one -- so the clean relaxation mechanism described below is now delivered by
    #   branch e, and the credit belongs there. The cause is not in the cascade code: a dozen SCF commits
    #   land between 04-Aug and today, among them a missing chain-rule factor in the one-electron gradient,
    #   Loewdin in place of Gram-Schmidt for the same-kappa orthogonalisation, and a corrected convergence
    #   test -- and SCA is the tier that runs a separate self-consistent field per block, so it is precisely
    #   the tier they act on. Rates do NOT improve up the ladder; see branch a.
    #
    # REPORT (04-Aug-2026, SUPERSEDED): 66.6 s against branch a's 40 s, i.e. SCA costs ~1.7x here. It is the
    # first time SCA has ever been executed in JAC.
    #   * K-alpha moves from 1247-1248 eV (branch a) to 1255 eV, against the measured 1253.6 eV. The error
    #     drops from about -6 eV to +1.4 eV, a factor ~4. So the hierarchy's central claim -- that the next
    #     tier up is genuinely better -- holds here, and is now measured rather than asserted.
    #   * Every Auger group shifts up with it: KL1L1 1083 -> 1088, KL1L2,3 1123-1133 -> 1122-1137,
    #     KL2,3L2,3 1163-1170 -> 1163-1174, KL1M1 1187 -> 1196-1198, KL2,3M1 1231 -> 1238, KM1M1 1271 -> 1290.
    #
    # WHERE THE IMPROVEMENT ACTUALLY COMES FROM -- and it is NOT configuration interaction. The block energies
    # identify the cause exactly:
    #     block 10, the initial 1s hole  : -4129 eV in BOTH branches, unchanged
    #     block  2, the 2p-hole final st.: -5376 eV (AverageSCA)  ->  -5384 eV (SCA)
    # K-alpha is the difference of these two, so the whole 8 eV gain comes from the FINAL state relaxing, and
    # none of it from the initial state. That is exactly what should happen: AverageSCA describes every block
    # with one global orbital set taken from the initial ion's mean field, so the initial block is already
    # optimal there while the decayed blocks are stuck with orbitals that never relax after the 1s hole is
    # filled. SCA gives each block its own self-consistent field, the final state relaxes and drops 8 eV, and
    # K-alpha rises onto the measurement. Standard final-state relaxation, cleanly isolated.
    #   Note the CSF counts per block are IDENTICAL in the two branches (1,2,1,1,4,2,5,4,1,1), and with at most
    #   five CSFs per block spread over several J^P symmetries there is very little for intra-block CI to do.
    #   So of the two axes this branch changes, it is the BOUND-ORBITAL axis that carries the physics here and
    #   the level-representation axis that contributes almost nothing -- worth remembering when weighing which
    #   axis to develop next, and a useful counterweight to the fact that CI is the cheaper of the two.
    setDefaults("print summary: open", "zzz-Cascade-Mg-Khole-SCA.sum")

    name   = "Stepwise decay after Mg 1s ionization (SCA)"
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0,
                                         Shell[], Shell[], Shell[])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid,
                                 approach=Basics.SCA(), scheme=scheme,
                                 initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fa.dat")
    setDefaults("print summary: close", "")
    #
    #
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch e: the SAME Mg K-hole system again, one further tier up -- Cascade.RefinedSCA(). Together with
    #   branches a (AverageSCA) and d (SCA) this gives the whole implemented ladder on identical physics, so
    #   that every difference between the three is attributable to the approach alone.
    #
    #   What RefinedSCA delivers here and what it does not: it relaxes the bound orbitals PER MULTIPLET, which
    #   is implemented and is where its gain comes from. Its second condition, continuum orbitals resolved per
    #   fine-structure transition, is NOT implemented; the run warns and proceeds with one set per step. That
    #   missing half affects Auger RATES only -- no level energy depends on the continuum orbital -- so the
    #   transition energies compared below are unaffected by the omission.
    #
    # REPORT (23-Aug-2026): 25-30 s warm, i.e. ~3x branch a. THIS IS THE BEST OF THE THREE TIERS ON THE ONE
    # quantity that has a measurement: K-alpha comes out at 1255.28 / 1254.99 eV against 1253.69 / 1253.44,
    # i.e. +1.6 eV or 0.13%, where AverageSCA is -5.9 eV and SCA +10.4 eV. The block energies show why, and
    # it is the mechanism branch d's 04-Aug note describes and no longer achieves: RefinedSCA leaves the
    # initial 1s-hole block at -4129 eV, exactly where AverageSCA has it, and relaxes the 2p-hole final
    # state from -5376 to -5384 eV. Final-state relaxation alone, cleanly isolated, with nothing added to
    # the initial state -- which is what a K-alpha energy should be sensitive to.
    #   The RATES do not follow. omega_K comes out at 0.065 (Bab) against a measured 0.030, the WORST of the
    #   three tiers, because the total Auger width falls to 2.56e14 /s while the radiative width rises only
    #   to 1.76e13 /s. So relaxing the bound orbitals buys the transition ENERGY and not the transition
    #   PROBABILITY, and the Auger deficit analysed in branch a survives the whole ladder untouched.
    setDefaults("print summary: open", "zzz-Cascade-Mg-Khole-RefinedSCA.sum")

    name   = "Stepwise decay after Mg 1s ionization (RefinedSCA)"
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], 1, Dict{Int64,Float64}(), 0,
                                         Shell[], Shell[], Shell[])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid,
                                 approach=Cascade.RefinedSCA(), scheme=scheme,
                                 initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fa.dat")
    setDefaults("print summary: close", "")
    #
end
