
println("Fc) Cascade.PhotoExcitationScheme: photo-excitation cascades and the configuration pre-filter.")

using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")                       ## set explicitly: the scheme's photon window is read
                                                        ## in the CURRENT user unit, see the note below

grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, rbox = 11.0)


# REWRITTEN 06-Aug-2026, first file of the scheme series (see the cascade-schemes plan). The previous Fc was
# a stepwise-decay cascade of neon 1s^-1 3p, which example-Fa.jl now covers properly, plus a simulation; it
# could not run in any case -- `using JLD` (the module is JLD2), a six-argument StepwiseDecayScheme that takes
# seven, a six-field SimulationSettings that takes three, and a hard-coded .jld file from July 2020. Nothing
# was kept.
#
# UNITS, worth stating because the module is not consistent here. PhotoExcitationScheme.minPhotonEnergy and
# .maxPhotonEnergy are read in the CURRENT USER UNIT and converted internally
# (Defaults.convertUnits("energy: from predefined to atomic unit", ...)), so the numbers below are eV. The
# corresponding window fields of Cascade.PhotonIntensities, by contrast, are taken as raw ATOMIC UNITS. Both
# are input fields of their structs; only the unit convention differs. This file sets the unit explicitly so
# that the window cannot be misread.
#
# HOW THE SCHEME SELECTS CONFIGURATIONS, which is the substance of this file. Before any amplitude is
# computed, Cascade.generateConfigurationsForPhotoexcitation filters the generated configurations twice:
#   (i)  by ENERGY, using an empirical estimate of the total energy (Empirical.totalEnergy with the X-ray
#        Data Booklet), keeping only configurations whose mean energy falls in
#        [minen + 0.2*minPhotonEnergy, maxen + 5*maxPhotonEnergy]. The 0.2 and 5 are deliberate safety
#        margins, because the empirical energies are rough;
#   (ii) by PARITY and multipole: for E1 only configurations of opposite parity to the initial ones survive,
#        while M1, E2 or M2 in the multipole list admit either parity.
# This is exactly the configuration-level filter that Symmetry 13, 520 (2021), Sect. 2.3 argues for -- and
# notably it is NOT applied on the stepwise-decay path, where example-Fa.jl branch a shows two of thirty
# steps being generated, set up and computed only to yield zero lines because they are E1-forbidden.
#
# COST NOTE, because it is counter-intuitive and cost the choice of test case here. A SMALLER system is not
# automatically a cheaper one. Boron-like C^+ (1s^2 2s^2 2p, 5 electrons, a single 2p -> 3s channel, 2 lines)
# was tried first as the cheapest conceivable smoke case and measured 26.7 s warm, repeatably -- against
# 5.0 s for the Ne^+ case below, which has 9 electrons and 27 lines. The C^+ SCF converges cleanly in ten
# iterations, so that is not the cause; the driver was not identified before the case was dropped. Plausible
# suspects are the very low transition energy (~10 eV, i.e. omega ~ 0.37 a.u.) and the diffuse 3s orbital
# (binding 0.19 a.u.) inside an rbox of 11 a.u. Worth remembering when choosing test cases for the other
# schemes: measure, do not assume.
#
# BUG FOUND AND FIXED BY BRANCH c, 06-Aug-2026. Cascade.determineSteps for this scheme built its
# PhotoExcitation.Settings with a hard-wired [E1], so PhotoExcitationScheme.multipoles acted on the
# configuration filter (ii) ONLY. Asking for [E1, M1, E2] therefore enlarged the cascade by same-parity
# blocks that were then set up, diagonalised and computed -- and could not yield a single line, because the
# amplitudes were still E1. Branch c returned 0 lines for both same-parity steps. The settings now take
# scheme.multipoles; see src/module-Cascade-inc-photoexcitation.jl:54.


if  true
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... 2 steps, 27 E1 lines; the parity filter keeps exactly the 2 of 4
    #                  generated configurations that are E1-reachable, and branch c reproduces these 27
    #                  lines identically as its E1 subset.  Warm cost ~6 s.
    #
    # Branch a: REFERENCE AND SMOKE CASE -- Ne^+ (1s^2 2s^2 2p^5) excited from both L subshells into 3s and 3p,
    #   E1 only. Two excited configurations are generated and both are reached, so the cascade has real
    #   branching rather than a single channel. At ~5 s warm (excluding the Julia start-up that a test inside
    #   runtests.jl does not pay) this is also the branch intended for the per-scheme smoke test, comfortably
    #   inside the 10 s budget. See the cost note above for why the ostensibly smaller C^+ case was rejected.
    #
    #   RESULT. Four configurations are generated from the two initial and two final shells; the parity filter
    #   rejects the two same-parity ones (2s -> 3s and 2p -> 3p) and keeps 2s^1 2p^5 3p^1 and 2s^2 2p^4 3s^1.
    #   Two steps follow, giving 16 and 11 lines, 27 in total. This is the filter working exactly as intended,
    #   and is the contrast to the stepwise-decay path noted above, which would have computed all four.
    setDefaults("print summary: open", "zzz-Cascade-Fc-neon.sum")

    name   = "Ne^+ 2s,2p -> 3s,3p photo-excitation"
    scheme = Cascade.PhotoExcitationScheme([E1], 1.0, 200.0, 1, [Shell("2s"), Shell("2p")],
                                           [Shell("3s"), Shell("3p")], LevelSelection(), [0,1], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fc.dat")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... 4 steps, 0 lines -- verified to be the physically correct answer, and
    #                  a documented negative result rather than a useful branch.  Warm cost ~7 s.
    #
    # Branch b: DOUBLE EXCITATION -- branch a with NoExcitations = 2, so that configurations reached by
    #   displacing TWO electrons are generated as well. This tests the generation depth rather than the
    #   filters, and shows how quickly the configuration count grows with that one parameter.
    #
    #   RESULT, and it is a negative one worth keeping. NoExcitations = 2 generates configurations displaced
    #   by EXACTLY two electrons, not by up to two: the four blocks are 2p^5 3s^1 3p^1 (2s emptied),
    #   2s^1 2p^4 3s^2, 2s^1 2p^4 3p^2 and 2p^3 3s^1 3p^1, and the singly-excited configurations that carry
    #   all the physics of branch a are gone. Every one of the four steps then returns 0 lines, which is
    #   correct -- a one-photon, one-electron operator cannot reach a doubly-excited configuration from the
    #   ground configuration. So ~7 s is spent generating, diagonalising (up to 93 CSF) and computing four
    #   steps whose answer is identically zero. The parity filter cannot prevent this, since it tests
    #   configuration parity only and 2p^5 3s^1 3p^1 does have the opposite parity.
    #
    #   TO CONSIDER: either make the generator cumulative (1 .. NoExcitations), which is what the field name
    #   suggests, or reject configurations differing by more than one electron in a one-photon scheme. Not
    #   changed here -- it is a generator question that touches the other schemes too.
    setDefaults("print summary: open", "zzz-Cascade-Fc-double.sum")

    name   = "Ne^+ photo-excitation, up to two displaced electrons"
    scheme = Cascade.PhotoExcitationScheme([E1], 1.0, 200.0, 2, [Shell("2s"), Shell("2p")],
                                           [Shell("3s"), Shell("3p")], LevelSelection(), [0,1], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fc.dat")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... 4 steps, 65 lines, after the multipole bug this branch exposed was
    #                  fixed; multipole hierarchy verified internally, see below.  Warm cost ~6 s.
    #
    # Branch c: THE PARITY FILTER -- branch b repeated with [E1, M1, E2] instead of E1 alone. The energy
    #   pre-filter is unchanged, so any difference in the surviving configurations comes purely from step (ii)
    #   above: with M1 or E2 present, configurations of the SAME parity as the initial one are no longer
    #   rejected. Expect more blocks than branch a, and hence more steps.
    #
    #   RESULT. All four L -> 3l blocks now survive and all four produce lines: 7 (2s -> 3s), 16 (2s -> 3p),
    #   11 (2p -> 3s) and 31 (2p -> 3p), 65 in total. Two independent internal checks:
    #     - the 27 pure-E1 lines are numerically the same 27 lines as branch a, i.e. adding M1 and E2 adds
    #       channels without disturbing the dipole ones;
    #     - the strengths order as they must. Coulomb-gauge oscillator strengths reach 0.132 for E1,
    #       1.5e-6 for E2 and 1.4e-7 for M1. The E2/E1 ratio of ~1e-5 sits at the retardation estimate
    #       (alpha*omega)^2 ~ 6e-5 for omega ~ 30 eV at Z = 10, and M1 is weaker again because its radial
    #       integral is the near-vanishing <3s|2s> overlap, non-zero only through the small components.
    #   BEFORE THE FIX this branch returned 0 lines for steps 1 and 4, which is what exposed the hard-wired
    #   [E1] in determineSteps -- see the note at the head of this file.
    setDefaults("print summary: open", "zzz-Cascade-Fc-multipoles.sum")

    name   = "Ne^+ photo-excitation, E1 + M1 + E2"
    scheme = Cascade.PhotoExcitationScheme([E1, M1, E2], 1.0, 200.0, 1, [Shell("2s"), Shell("2p")],
                                           [Shell("3s"), Shell("3p")], LevelSelection(), [0,1], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fc.dat")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... filter verified to act, but far more loosely than intended; 2 of 3
    #                  configurations excluded, 1 step and 11 lines survive.  Warm cost ~4 s.
    #
    # Branch d: THE ENERGY FILTER -- branch b with the photon window closed down to a few eV, far below what
    #   an L -> 3l excitation of Ne^+ costs. Every generated configuration should now fall outside
    #   [minen + 0.2*minPhotonEnergy, maxen + 5*maxPhotonEnergy] and be rejected, each with a printed
    #   ">>> exclude ... because of energy reasons" line. The branch exists to show the filter working, and to
    #   record what the scheme does when it selects nothing at all -- an empty cascade ought to be reported
    #   clearly rather than fail confusingly further down.
    #
    #   RESULT, which did NOT match that expectation and is the more interesting for it. Only 2s^1 2p^5 3s^1
    #   and 2s^1 2p^5 3p^1 are excluded; 2s^2 2p^4 3s^1 survives and is computed, giving 11 lines at photon
    #   energies near 27 eV -- five times the top of the requested 1-5 eV window. The cause is the deliberate
    #   safety margin: the upper bound is maxen + 5*maxPhotonEnergy, i.e. 25 eV of slack on a 5 eV window, so
    #   the filter cannot be closed tightly. Two consequences worth stating:
    #     - the photon window is a coarse pre-filter on CONFIGURATIONS, not a window on the computed lines.
    #       Anyone wanting only lines inside [minPhotonEnergy, maxPhotonEnergy] must filter afterwards;
    #     - the two excluded configurations are reported with the SAME energy, -69.67304 a.u., because the
    #       empirical estimate (Empirical.totalEnergy) does not resolve 3s from 3p at all. The margins are
    #       generous precisely because that estimate is crude, so the two facts belong together.
    #   The intended empty-cascade case was therefore not reached and remains untested.
    setDefaults("print summary: open", "zzz-Cascade-Fc-window.sum")

    name   = "Ne^+ photo-excitation, window closed to 1-5 eV"
    scheme = Cascade.PhotoExcitationScheme([E1], 1.0, 5.0, 1, [Shell("2s"), Shell("2p")],
                                           [Shell("3s"), Shell("3p")], LevelSelection(), [0,1], 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true, outputDirectory="example-Fc.dat")
    setDefaults("print summary: close", "")
    #
end
