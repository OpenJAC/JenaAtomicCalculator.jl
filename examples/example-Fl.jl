
println("Fl) Cascade.ElectronIonizationScheme, RESONANT channels: resonant electron capture with sequential or")
println("    simultaneous double autoionization, for lithium-like carbon.")
println("    Companion to example-Fi.jl, which computes the OTHER indirect channel of the same scheme --")
println("    impact-excitation with subsequent autoionization -- on the SAME ion, so that the two may be")
println("    compared.  They are kept in separate files because their costs differ by two orders of")
println("    magnitude: a capture is the time reverse of an Auger and needs no partial-wave sum over impact")
println("    energies, so this file runs in a minute where Fi branch a takes 47.")

using JLD2, Printf

grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 14.0)
setDefaults("standard grid", grid)
## Radial.Grid(true) must NOT be used here: it carries hp = 0, and the continuum orbitals of an Auger step
## then cannot be generated at all.

if  true
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had
    #    been resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    #   72 s.  79 resonances between 317 and 324 eV, 8 steps, written to
    #   zzz-cascade-electron-ionization-computations-<date>.jld, which branches b to d read back.
    #
    # Branch a: THE COMPUTATION.  An incident electron is CAPTURED by Li-like carbon into a doubly-excited
    #   resonance, which then sheds TWO electrons, leaving He-like carbon -- one charge state up from where it
    #   started.  Two routes are requested at once, and they differ only in how the resonance decays:
    #     SequentialAuger()     the electrons leave one after the other, through an intermediate that
    #                           autoionizes in its turn   (REDA in the older literature)
    #     SimultaneousAuger()   they leave together, in one double-Auger step   (READI)
    #
    #   THE CHOICE OF SHELLS IS THE WHOLE DIFFICULTY, and it is worth stating because the obvious choice does
    #   not work.  A resonance converges to the threshold it was built on FROM BELOW -- that is what being a
    #   resonance means -- so it can never autoionize into that same threshold.  Capturing into 2s/2p while
    #   exciting 1s -> 2p gives resonances below the 1s2s2p threshold, and the sequential route is then
    #   IDENTICALLY CLOSED: measured, every strength came out zero.  The route needs a resonance on a HIGH
    #   threshold decaying to a LOWER state that is itself autoionizing.  Hence the excitation list spans both
    #   2p and 3s/3p while the capture goes into 3s/3p: the resonances then sit on the 1s2s3l thresholds and can
    #   reach 1s2s2p, which still carries the K hole and so autoionizes again to 1s^2.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-computation.sum")
    #
    scheme = Cascade.ElectronIonizationScheme(Float64[], [Shell("1s")], [Shell("2p"), Shell("3s"), Shell("3p")],
                                              collect(0:5), 1, 0.,
                                              Basics.AbstractProcess[ResonantImpactIonization.SequentialAuger(),
                                                                     ResonantImpactIonization.SimultaneousAuger()],
                                              [Shell("3s"), Shell("3p")], 0.)
    wa = Cascade.Computation(Cascade.Computation(); name="Resonant ionization of Li-like C",
                             nuclearModel=Nuclear.Model(6.), grid=grid, approach=Cascade.AverageSCA(),
                             scheme=scheme, initialConfigs=[Configuration("1s^2 2s")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: as branch a.]  Seconds; no amplitude is recomputed.
    #   79 resonances.  Totals, energy-integrated and in atomic units:
    #       S(sequential) = 4.164005e-04     S(recombination) = 1.876036e-04
    #   and all 158 branching sums came out at exactly 1.00000000.
    #
    # Branch b: THE COMPETING DECAYS, and the branch to read first.  A doubly-excited resonance either RADIATES,
    #   which is dielectronic recombination, or AUTOIONIZES, which is the ionization this file is about.  They
    #   are decays of THE SAME resonance, so one computation gives both and the competition can be read off
    #   directly -- which is the reason to compute resonant ionization inside a cascade rather than on its own.
    #
    #   THE PRIMARY CHECK is the last column: for every resonance the radiative and all autoionization
    #   branchings must sum to 1 to machine precision.  It needs no literature value and no external reference.
    #   BUT READ WHAT IT CHECKS.  Every branching is formed from the total rates of the steps the cascade
    #   actually generated, so a decay route left out of the configuration lists is missing from the denominator
    #   too, and the sum is 1 regardless.  It checks the ARITHMETIC, never the completeness.  A column of exact
    #   ones is therefore necessary and nowhere near sufficient, and it would be easy to mistake for more.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-simulation.sum")
    #
    fn = filter(f -> startswith(f, "zzz-cascade-electron-ionization"), readdir("."))[end]
    println(">> reading $fn")
    sim = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization, strengths",
                             property=Cascade.ResonantIonizationStrengths(1, 0., 0.),
                             method=Cascade.ProbPropagation(), computationData=[JLD2.load(fn)] )
    perform(sim; output=false)
    println("\n     S(sequential) exceeds S(recombination) by about a factor of two here.  That is the expected")
    println("     ordering for a K-hole resonance in a LIGHT ion, where Auger decay dominates radiative decay by")
    println("     orders of magnitude; it reverses at high Z, where the radiative width grows as Z^4.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:  23-Aug-2026
    # Last successful:  unknown ... (this branch QUOTES example-Fi.jl's dated result rather than re-running it;
    #                   a date is not set because nothing here is verified by running it.  To make it a real
    #                   comparison, run Fi branch a on this same grid and put its number in by hand.)
    #
    # Branch c: AGAINST THE OTHER INDIRECT CHANNEL.  Of the three contributions to electron-impact ionization,
    #   this scheme now computes two and the direct one is not implemented at all:
    #
    #     direct                                              Cascade.ImpactIonizationScheme, NOT implemented
    #     impact-excitation with subsequent autoionization     example-Fi.jl, 2836 s at maxKappa = 20
    #     resonant electron capture, sequential/simultaneous   this file, 72 s
    #
    #   WHAT THIS BRANCH DOES NOT DO is put the two on a common footing.  Fi reports a CROSS SECTION at chosen
    #   impact energies; this file reports ENERGY-INTEGRATED RESONANCE STRENGTHS.  The two are different
    #   quantities and cannot be compared without folding the resonances with an electron energy distribution,
    #   which nothing here does.  Setting them side by side would be exactly the kind of comparison that looks
    #   quantitative and is not, so the branch prints the two numbers with their units and stops.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-simulation.sum")
    #
    println("\n  The two indirect channels of electron-impact ionization for Li-like carbon:\n")
    println("    impact-excitation with subsequent autoionization")
    println("       Omega (summed collision strength) = 0.376 / 0.521 / 0.526 at 1000 / 2000 / 4000 eV")
    println("       from example-Fi.jl branch a, dated 08-Aug-2026, at maxKappa = 20; QUOTED, not re-run here.")
    println("    resonant electron capture, sequential")
    println("       S = 4.164e-04 a.u., energy-integrated over 79 resonances between 317 and 324 eV")
    println("       from branch b of this file.")
    println("\n    THESE ARE DIFFERENT QUANTITIES.  A collision strength at an impact energy and an")
    println("    energy-integrated resonance strength are not comparable until the resonances are folded with an")
    println("    electron energy distribution -- a Maxwellian, say -- which would turn both into a rate")
    println("    coefficient.  That fold is not implemented for the resonant channel, so no ratio is quoted here.")
    println("    What CAN be said without it: the resonances lie between 317 and 324 eV, i.e. far below the")
    println("    1000 eV at which the impact-excitation channel was evaluated, so in a plasma cool enough to")
    println("    populate the resonances and not the excitation, the resonant channel is the only one of the two")
    println("    that contributes at all.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: as branch a.]
    #   P(double Auger) = 0.06245 for neon with 2s and 2p passive, and the self-projection limit came out at
    #   exactly 0.0.  At that probability the simultaneous route contributes S = 3.043629e-04 against the
    #   sequential 4.164005e-04, i.e. about three quarters of it.
    #
    # Branch d: THE SIMULTANEOUS ROUTE, AND WHY ITS NUMBER IS THE WEAKEST THING IN THIS FILE.  A double-Auger
    #   width cannot be had without amplitudes, which this module deliberately does not compute.  What it offers
    #   instead is SHAKE-OFF in the sudden approximation: when the first electron leaves, the potential changes
    #   abruptly and a second electron may fail to remain bound.  That probability follows from overlaps of the
    #   orbitals before and after, and JAC has both.
    #
    #   ITS ERROR HAS BOTH SIGNS AND IT IS NOT A BOUND.  A finite set of bound final orbitals lets shake-UP leak
    #   in and INFLATES it; KNOCKOUT -- the departing electron striking a second on its way out -- is a collision,
    #   unreachable from overlaps, comparable to or larger than shake-off, and missing altogether, which DEFLATES
    #   it.  The two do not cancel in any controlled way.
    #
    #   THE CASCADE DOES NOT CHOOSE THE NUMBER FOR YOU, and that is deliberate.  Which subshells count as PASSIVE
    #   is a physics judgement: the two electrons that participate in the Auger transition must be excluded, and
    #   which those are depends on the transition.  So dblAugerProbability defaults to 0., the simultaneous route
    #   then contributes exactly nothing, and a value has to be put in knowingly.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-simulation.sum")
    #
    println("\n  (i)  the shake-off probability, and the limit it must satisfy")
    orbsOf(cfg) = collect(values(SelfConsistent.performSCF([Configuration(cfg)], Nuclear.Model(10.), grid,
                                                            AsfSettings(); printout=false).levels[1].basis.orbitals))
    oi = orbsOf("1s^2 2s^2 2p^6");    of = orbsOf("1s^2 2s^2 2p^5")
    pick(os, s) = os[findfirst(o -> string(o.subshell) == s, os)]
    pass = [pick(oi,"2s_1/2"), pick(oi,"2p_1/2"), pick(oi,"2p_3/2")];    occ = [2, 2, 4]
    println("       projecting the orbitals on THEMSELVES must give exactly zero:  " *
            string(round(ResonantImpactIonization.doubleAugerProbability(pass, occ, oi, grid), sigdigits=3)))
    pDbl = ResonantImpactIonization.doubleAugerProbability(pass, occ, of, grid)
    println("       neon, 2s and 2p passive, one electron fewer after the Auger: P = " *
            string(round(pDbl, digits=5)))
    println("       An Auger raises the IONIC charge at fixed Z.  Comparing orbitals at different NUCLEAR")
    println("       charges is the wrong setup and inflates this several-fold; it gave 0.21 when tried.")
    #
    println("\n  (ii) what that does to the simultaneous channel")
    fn = filter(f -> startswith(f, "zzz-cascade-electron-ionization"), readdir("."))[end]
    for p in [0.0, pDbl]
        println("\n       ---- dblAugerProbability = " * string(round(p, digits=5)) * " ----")
        sim = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization, strengths",
                                 property=Cascade.ResonantIonizationStrengths(1, 0., p),
                                 method=Cascade.ProbPropagation(), computationData=[JLD2.load(fn)] )
        perform(sim; output=false)
    end
    println("\n     At P = 0.062 the simultaneous route reaches about three quarters of the sequential one.")
    println("     The neon probability is used here for want of a carbon one computed for the actual")
    println("     transition; it is an ORDER OF MAGNITUDE standing in for a number, which is all this route")
    println("     supports.  Read it as: the simultaneous channel is not negligible beside the sequential one,")
    println("     and nothing sharper than that.")
    #
    setDefaults("print summary: close", "")
    #
end
