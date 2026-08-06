
println("Fd) Cascade.PhotoIonizationScheme: photo-ionization cross sections of Ne^+.")

using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")                       ## the scheme's photonEnergies are read in THIS unit
setDefaults("unit: cross section", "barn")              ## set explicitly; see the unit note below

grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)


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
# HOW GOOD ARE THE NUMBERS?  Honest answer, 06-Aug-2026: the thresholds are right, the overall magnitude is
# plausible, and the finer structure is NOT trustworthy.  No branch below is therefore marked
# "Last successful".  What was established:
#   + Thresholds.  The five Ne^2+ (2p^4) levels come out at 39.3 - 44.8 eV against a measured Ne^+ ionization
#     potential of 40.96 eV, and the 2s-hole levels at 65.8 - 73.8 eV.  Good.
#   + Magnitude.  5.8 Mb (Coulomb) for 2p at 80 eV is the right order for a neon 2p subshell.
#   - Gauge.  Coulomb and Babushkin disagree by factors of 1.5 to 4.5 depending on energy and subshell.  That
#     is an internal inconsistency, not a matter of taste, and it alone rules out quoting these numbers.
#   - Subshell ratio.  At 120 eV the computed 2s cross section EXCEEDS the 2p one by 1.65 (Coulomb) or 3.05
#     (Babushkin), where neon data put 2s/2p near 0.2.  The ordering is inverted, not merely inaccurate.
#   - Partial waves.  For 2p ionization the photo-electron must be s or d (E1), and d should dominate far
#     above threshold.  Branch d finds the s wave alone carries 56% of the total.
#   RULED OUT as the cause: the continuum normalization.  Repeating branch a with
#   "method: normalization, pure sine" and "pure Coulomb" changes the total by 2.5e-5 relative
#   (0.2062176 vs 0.2062228 a.u.), and both branches of Continuum are genuinely entered (58 orbitals each).
#   At a 40 eV photo-electron energy the Coulomb correction to the asymptotic normalization is simply
#   negligible, so this knob cannot explain anything here.
#   The remaining suspect is the AverageSCA treatment itself -- single-CSF blocks with no configuration
#   mixing, bound orbitals from a DFS field of the initial configuration, the continuum orbital in the field
#   of the final ion, and no interchannel coupling.  Photoionization is far more sensitive to all four than
#   the decay rates that this approach was designed for.  Establishing that needs an external reference or a
#   better approach, and SCA() is not implemented, so it is left open here.


if  true
    # Last visit:      06-Aug-2026 ... runs; 5 lines, total 5.775 Mb (Coulomb) / 1.498 Mb (Babushkin).
    #                  Warm cost ~7 s, so it fits the smoke budget.  NOT "Last successful": the two gauges
    #                  disagree by 3.9, see the assessment above.
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
    # Last visit:      06-Aug-2026 ... 25 lines, ~24 s.  Totals in Mb (Coulomb / Babushkin):
    #                      50 eV  1.148 / 0.496      80 eV  5.775 / 1.498     120 eV  3.409 / 0.849
    #                     200 eV  0.965 / 0.409     300 eV  0.314 / 0.214
    #                  The curve peaks at 80 eV, i.e. ~40 eV above threshold, and falls by a factor 18 out
    #                  to 300 eV.  A delayed maximum is not by itself wrong for an l -> l+1 channel with a
    #                  centrifugal barrier, so this shape is NOT being called an error -- but neither is it
    #                  confirmed, and the gauge spread over the range (1.5 to 4.0) prevents dating it.
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
    # Last visit:      06-Aug-2026 ... two steps, 9 lines, ~12 s.  At 120 eV the 2s block gives 5.63 Mb and
    #                  the 2p block 3.41 Mb (Coulomb), i.e. 2s/2p = 1.65 where neon data put it near 0.2 --
    #                  the ordering is inverted.  One useful internal check does pass: the 2p part here
    #                  reproduces branch b's 120 eV numbers line by line, so adding the 2s block does not
    #                  disturb the 2p one.  NOT "Last successful".
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
    # Last visit:      06-Aug-2026 ... 5 lines, ~5 s, total 3.243 Mb (Coulomb) against 5.775 Mb for branch a.
    #                  lValues IS honoured -- the number changes -- but the result says the s wave alone
    #                  carries 56% of the 2p cross section, whereas the d wave should dominate at a 40 eV
    #                  photo-electron energy.  Reported as a finding, not as a convergence success.
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
end
