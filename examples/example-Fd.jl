
println("Fd) Cascade.PhotoIonizationScheme: photo-ionization cross sections of Ne^+, calibrated on neutral Ne.")

using JLD2
#
setDefaults("method: continuum, Galerkin")              ## setDefaults("method: continuum, asymptotic Coulomb")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")                       ## the scheme's photonEnergies are read in THIS unit
setDefaults("unit: cross section", "barn")              ## set explicitly; see the unit note below

grid = Radial.Grid(Radial.Grid(false); rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


# REWRITTEN 06-Aug-2026, second file of the scheme series (see the cascade-schemes plan). Two earlier files
# are superseded by this one and nothing was kept from either:
#   * the previous example-Fd.jl used `using JLD` (the module is JLD2), a Cascade.PhotonExcitationScheme and
#     a Cascade.PhotonIonizationScheme (neither type exists; they are PhotoExcitationScheme and
#     PhotoIonizationScheme), a `NoPoints=` grid keyword that has been replaced by `rbox=`, and three
#     hard-coded .jld files from 2020 and 2021;
#   * example-Fj.jl held the photo-ionization scenario rescued from the old Fb (Si^- at 30 and 80 eV) with
#     the same non-existent PhotonIonizationScheme and a three-field constructor for a nine-field type.
# Fj is thereby free for the HollowIonScheme, as the plan intends.
#
# UNITS -- the one that will bite. Cross sections are printed in the CURRENT cross-section unit, and JAC's
# default is BARN, not Mb. A Ne^+ 2p cross section printed as 5.8e6 is 5.8 Mb. The photon energies of the
# scheme, by contrast, are read in the current ENERGY unit (eV here); the module prints
# ">>> Photon energies must still be given in user-selected units" at every step, which is a warning and not
# a conversion. Both units are therefore set explicitly at the top of this file.
#
# WHAT THIS SCHEME ACTUALLY READS, checked field by field on 06-Aug-2026. PhotoIonizationScheme declares
# nine fields; FOUR are never read anywhere in module-Cascade-inc-photoionization.jl, and a fifth is
# bypassed:
#   + multipoles              used
#   + photonEnergies          used
#   + excitationFromShells    used -- this alone decides which subshells are ionized
#   + lValues                 used -- the partial waves of the photo-electron
#   - electronEnergies        DEAD; never read
#   - electronEnergyShift     DEAD; never read
#   - minCrossSection         DEAD; never read, so no weak lines are suppressed
#   - excitationToShells      DEAD; it occurs only inside a docstring
#   - initialLevelSelection   BYPASSED; it survives only in the commented-out line 32, while line 33
#                             hard-codes LevelSelection(true, indices=[1]).  ONLY the ground level is
#                             ionized, whatever the user asks for.
# The last point matters to anyone building a photo-ionization cascade from an excited or metastable initial
# level: the request is accepted silently and ignored. It is the same class of defect as the hard-wired [E1]
# that example-Fc.jl uncovered in the photo-excitation scheme, and deserves the same fix -- but that is a
# second module, so it is recorded here rather than changed.
#
# Two further observations on the implementation, neither of them changed here:
#   + Cascade.perform(::PhotoIonizationScheme, ...) accepts an outputDirectory argument and never uses it;
#     the .jld file is always written into the working directory.
#   + generateBlocks() for this scheme still builds its orbitals with performSCF per configuration rather
#     than through Cascade.generateBoundOrbitals(comp.approach, ...), so the approach conditions introduced
#     for the stepwise-decay path do not reach it.  AverageSCA is in any case the ONLY usable approach here:
#     generateBlocks() answers SCA() with error("Not yet implemented.").
#
# ONE REAL BUG, FIXED 06-Aug-2026 while writing this file. Base.show(io, ::PhotoIonizationScheme) ended with
#     if  length(photonEnergies) > 0  &&  length(electronEnergies) > 0
# referring to bare variables instead of the fields of the given scheme, so that printing a
# PhotoIonizationScheme -- and hence any Cascade.Computation carrying one -- raised an UndefVarError.  Since
# every example does println(wa), the scheme could not be used in the documented way at all; a smoke test
# that skips the printout runs fine, which is why this survived.  See src/module-Cascade.jl.
#
# HOW GOOD ARE THE NUMBERS?  Answered 23-Aug-2026, and the answer is now split by ENERGY rather than left
# open.  The file was switched to Galerkin above, which is what its own 07-Aug note recommended, and all four
# branches were re-run on the finer grid.  ALL THREE defects recorded in August disappear:
#         gauge ratio      1.5 - 4.5   ->   1.15 - 1.33  at every energy
#         2s/2p at 120 eV  1.65        ->   0.305        (neon data indicate ~0.2)
#         d-wave share     44%         ->   96.8%        (the s wave had carried 56%)
# so all three were artefacts of setDefaults("method: continuum, asymptotic Coulomb"), which gives erratic
# continuum orbitals for valence-shell processes in low-charge ions.  The old numbers are kept in the branch
# texts below, marked SUPERSEDED, because the mistake is worth being able to recognise again.
#
# THE EXTERNAL REFERENCE, which this file previously said it needed.  Branch e ionizes NEUTRAL neon, whose 2p
# cross section is among the best measured in the periodic table, so the calculation can be held against a
# number rather than against an expectation.  Computed (Coulomb / Babushkin, Mb) against the measured curve:
#         25 eV   21.0 / 16.0     vs ~5        40 eV   13.4 /  9.6     vs ~8   (near the measured maximum)
#         30 eV   18.2 / 13.3     vs ~6.5      60 eV    6.87 /  4.92   vs ~6.3
#         50 eV    9.16 /  6.46   vs ~7.5      80 eV    4.14 /  3.06   vs ~4.0
#                                             120 eV    1.82 /  1.44   vs ~2.0
#   + WELL ABOVE THRESHOLD the scheme is now good: within 10-20% at 80 and 120 eV, with the two gauges 1.3
#     apart.  Numbers from this region can be used, with that spread quoted.
#   - NEAR THRESHOLD it is wrong, and wrong in SHAPE and not merely in size.  Neon's measured cross section
#     RISES from about 4 Mb at the 21.56 eV threshold to a delayed maximum near 45 eV; the computed curve
#     falls MONOTONICALLY from 21 Mb at 25 eV.  The overshoot is 1.7x at 40 eV, 2.7x at 30 eV and ~4x at
#     25 eV.  The delayed maximum is a centrifugal-barrier effect -- the 2p -> eps_d channel cannot penetrate
#     the barrier just above threshold -- and its absence says the outgoing d wave carries too much amplitude
#     at low energy.  That is what a frozen-core DFS continuum without polarisation or interchannel coupling
#     would give, and it is a property of AverageSCA rather than of the continuum method: switching to
#     Galerkin fixed the gauge, the subshell ratio and the partial waves, and did NOT fix this.
#   The practical rule: TRUST THIS SCHEME MORE THAN ABOUT 40 eV ABOVE THRESHOLD, AND NOT BELOW.  For Ne^+
#   (threshold ~41 eV) that puts branch b's 50 eV point, only ~9 eV above threshold, in the untrustworthy
#   region while its 120-300 eV points are in the good one.
#
#   RULED OUT as the cause: the continuum normalization.  Repeating branch a with
#   "method: normalization, pure sine" and "pure Coulomb" changes the total by 2.5e-5 relative
#   (0.2062176 vs 0.2062228 a.u.), and both branches of Continuum are genuinely entered (58 orbitals each).
#   At a 40 eV photo-electron energy the Coulomb correction to the asymptotic normalization is simply
#   negligible, so this knob cannot explain anything here.


if  true
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # REPORT (23-Aug-2026, Galerkin): 5 lines, total 18.148 Mb (Coulomb) / 13.686 Mb (Babushkin), i.e. a
    #   gauge ratio of 1.33 where the asymptotic-Coulomb run of 06-Aug gave 3.9 on the same case. 52.6 s as
    #   the cold first run of a session; the smaller branch d runs warm in 14.3 s. At 80 eV the photo-electron
    #   carries ~38 eV, which is just inside the region branch e certifies, so this number is usable with its
    #   33% gauge spread quoted alongside. SUPERSEDED reading of 06-Aug: 5.775 / 1.498 Mb.
    #
    # Branch a: REFERENCE AND SMOKE CASE -- Ne^+ (1s^2 2s^2 2p^5) ionized out of the 2p subshell at a single
    #   photon energy of 80 eV, with s, p and d partial waves for the photo-electron. This is the branch
    #   intended for the per-scheme smoke test in runtests.jl.
    setDefaults("print summary: open", "zzz-Cascade-Fd-reference.sum")

    name   = "Ne^+ 2p photo-ionization at 80 eV"
    scheme = Cascade.PhotoIonizationScheme([E1], [80.0], Float64[], [Shell("2p")], Shell[],
                                           LevelSelection(), [0,1,2], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # REPORT (23-Aug-2026, Galerkin): 21.6 s warm, 25 lines. Totals in Mb (Coulomb / Babushkin):
    #        50 eV  52.46 / 41.22     80 eV  18.15 / 13.69     120 eV  7.348 / 5.667
    #       200 eV   2.100 /  1.725  300 eV   0.705 / 0.611
    #   with gauge ratios of 1.27, 1.33, 1.30, 1.22, 1.15 -- uniform and modest, against the 1.5 to 4.0 of
    #   the superseded run below. The SHAPE has also changed: the old curve peaked at 80 eV, and this one
    #   falls monotonically from 50 eV onwards. That is not automatically the better answer, and branch e is
    #   what decides it: neon really does have a delayed maximum, and the computed neutral-neon curve does
    #   NOT reproduce it either, falling monotonically where the measurement rises. So the 06-Aug curve had
    #   a maximum for the wrong reason and this one lacks a maximum for a real one.
    #   READ THIS TABLE FROM 120 eV UPWARDS. The 50 eV point sits only ~9 eV above the Ne^+ threshold, deep
    #   in the region where branch e shows the calculation overshooting by up to 4x, and 52 Mb should not be
    #   quoted. Dated for the 120-300 eV points and for the gauge behaviour, not for the near-threshold end.
    #   SUPERSEDED (06-Aug, asymptotic Coulomb): 1.148/0.496, 5.775/1.498, 3.409/0.849, 0.965/0.409,
    #   0.314/0.214 Mb at the same five energies.
    #
    # Branch b: PHOTON-ENERGY DEPENDENCE -- the same 2p ionization over a range of photon energies from just
    #   above threshold to well above it. This is the branch that can be checked against something external:
    #   the 2p photo-ionization cross section of neon is among the best measured in the periodic table,
    #   peaking near 8-9 Mb a few eV above threshold and falling to the 1 Mb level by 150-200 eV. Ne^+ has
    #   five 2p electrons rather than six and a higher threshold, so its cross section should be somewhat
    #   smaller and shifted upwards in energy, but of the same shape and order.
    setDefaults("print summary: open", "zzz-Cascade-Fd-energies.sum")

    name   = "Ne^+ 2p photo-ionization, 50 - 300 eV"
    scheme = Cascade.PhotoIonizationScheme([E1], [50.0, 80.0, 120.0, 200.0, 300.0], Float64[],
                                           [Shell("2p")], Shell[], LevelSelection(), [0,1,2], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # REPORT (23-Aug-2026, Galerkin): 23.1 s warm, two steps, 9 lines. At 120 eV
    #        2s block  2.240 / 2.180 Mb        2p block  7.348 / 5.667 Mb        2s/2p = 0.305 (Coulomb)
    #   THE INVERSION IS GONE. The superseded run below put 2s/2p at 1.65, i.e. the 2s subshell ionizing
    #   more readily than the 2p one, which is not what neon does; the ratio is now 0.305 against the ~0.2
    #   the neon data indicate -- the right side of unity and the right order of magnitude.
    #   The internal check also still passes, and passes exactly: the 2p total here, 7.347569e+06 barn,
    #   reproduces branch b's 120 eV entry to all seven printed digits, so adding a second ionized block
    #   leaves the first one numerically untouched. SUPERSEDED (06-Aug): 5.63 Mb and 3.41 Mb, ratio 1.65.
    #
    # Branch c: SUBSHELL DECOMPOSITION -- ionization out of 2s AND 2p in one computation. Two ionized blocks
    #   are generated, so the two subshells appear as separate steps and their cross sections can be compared
    #   directly. Physically the 2s cross section should lie well below the 2p one at these energies, and the
    #   two photo-electron energies should differ by the 2s-2p binding-energy difference.
    setDefaults("print summary: open", "zzz-Cascade-Fd-subshells.sum")

    name   = "Ne^+ 2s and 2p photo-ionization at 120 eV"
    scheme = Cascade.PhotoIonizationScheme([E1], [120.0], Float64[], [Shell("2s"), Shell("2p")], Shell[],
                                           LevelSelection(), [0,1,2], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # REPORT (23-Aug-2026, Galerkin): 14.3 s warm, 5 lines, total 0.585 / 0.504 Mb against branch a's
    #   18.148 / 13.686 Mb with the d wave included. So dropping l = 2 removes 96.8% of the cross section,
    #   which is what an E1 ionization out of 2p must do well above threshold: the photo-electron can only
    #   be s or d, and at a ~38 eV photo-electron energy the d wave carries essentially all of it.
    #   This is the check that changed most. The superseded run below had the s wave alone carrying 56% of
    #   the total, which was reported at the time as a finding about the partial-wave sum; it was in fact a
    #   third symptom of the asymptotic-Coulomb continuum, and it disappears with the rest. lValues is
    #   honoured, the partial-wave sum is not silently truncated, and the hierarchy is now the physical one.
    #   SUPERSEDED (06-Aug): 3.243 Mb against 5.775 Mb, i.e. an s wave carrying 56%.
    #
    # Branch d: PARTIAL-WAVE CONVERGENCE -- branch a repeated with only s and p partial waves. For an E1
    #   ionization out of 2p the photo-electron must be s or d, and the d wave dominates well above
    #   threshold. Dropping l = 2 should therefore remove most of the cross section, which is the cheapest
    #   available test that lValues is honoured and that the partial-wave sum is not silently truncated.
    setDefaults("print summary: open", "zzz-Cascade-Fd-lvalues.sum")

    name   = "Ne^+ 2p photo-ionization at 80 eV, s and p waves only"
    scheme = Cascade.PhotoIonizationScheme([E1], [80.0], Float64[], [Shell("2p")], Shell[],
                                           LevelSelection(), [0,1], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026
    #
    # Branch e: THE EXTERNAL CALIBRATION -- the same scheme applied to NEUTRAL neon, whose 2p photo-ionization
    #   cross section is among the best measured in the periodic table. Everything above computes Ne^+, for
    #   which no comparable measurement exists, so this branch is what turns the file from internally
    #   consistent into externally checked. Seven photon energies from just above the 21.56 eV threshold out
    #   to well beyond it, so that the SHAPE of the curve is tested and not merely one magnitude.
    #
    # REPORT (23-Aug-2026): 51.0 s cold. Two final levels per energy, since Ne^+ 2p^5 is 3/2- and 1/2-.
    #   Computed totals in Mb, Coulomb / Babushkin, against the measured curve:
    #         25 eV   21.0  / 16.0     vs ~5           60 eV    6.87 / 4.92    vs ~6.3
    #         30 eV   18.2  / 13.3     vs ~6.5         80 eV    4.14 / 3.06    vs ~4.0
    #         40 eV   13.4  /  9.6     vs ~8          120 eV    1.82 / 1.44    vs ~2.0
    #         50 eV    9.16 /  6.46    vs ~7.5
    #   The verdict is a BOUNDARY rather than a single number, and it is stated in full at the head of this
    #   file: good to 10-20% above about 60 eV, too large by up to 4x below it, and wrong in SHAPE there --
    #   the measured neon cross section RISES from ~4 Mb at threshold to a delayed maximum near 45 eV, while
    #   the computed one falls monotonically from 21 Mb at 25 eV. The delayed maximum is the centrifugal
    #   barrier keeping the eps_d wave out at low photo-electron energy, and its absence says that wave
    #   carries too much amplitude near threshold -- a frozen-core DFS continuum without polarisation or
    #   interchannel coupling, i.e. a property of AverageSCA and not of the continuum method, since Galerkin
    #   repaired the gauge, the subshell ratio and the partial waves and did not repair this.
    #   Worth keeping methodologically: the gauge ratio is a uniform ~1.3 over the WHOLE range, including
    #   the energies where the result is four times the measurement. The two gauges agree with each other
    #   where both disagree with nature, so gauge consistency is a necessary check and not a sufficient one.
    setDefaults("print summary: open", "zzz-Cascade-Fd-neutralNe.sum")

    name   = "NEUTRAL Ne 2p photo-ionization, calibration against measurement"
    scheme = Cascade.PhotoIonizationScheme([E1], [25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 120.0], Float64[],
                                           [Shell("2p")], Shell[], LevelSelection(), [0,1,2], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^6")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
end
