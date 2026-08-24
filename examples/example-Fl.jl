
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


## PICKING THE RIGHT .jld MATTERS HERE, and it did not use to.  Every branch of this file wants the RESONANT data
## that branch a writes -- but example-Fi.jl writes files with the SAME prefix into the SAME directory, and
## readdir(".")[end] sorts ALPHABETICALLY, so as soon as both exist it silently returns whichever name happens to
## sort last.  On 24-Aug-2026 that was Fi's, and branch e duly folded the excitation-autoionization data while
## reporting them as resonant.  Nothing in the output said so.  Choose by CONTENT instead: the resonant data are the
## ones carrying no impact-excitation lines.  Files that no longer load are named and skipped rather than ignored.
function resonantCascadeFile()
    cands = sort(filter(f -> startswith(f, "zzz-cascade-electron-ionization"), readdir(".")), by = f -> mtime(f), rev = true)
    for  f  in  cands
        local d
        try     d = JLD2.load(f)
        catch e println(">> skipping $f -- it cannot be loaded: ", first(split(sprint(showerror, e), "\n")));   continue
        end
        res = d["results"]
        if  !( haskey(res, "impact-excitation lines:")  &&  length(res["impact-excitation lines:"]) > 0 )
            println(">> resonant cascade data from  $f");    return( (f, d) )
        end
    end
    error("No loadable cascade file WITHOUT impact-excitation lines was found in the working directory; " *
          "run branch a of this file first.")
end


if  true
    # Last successful:  23-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had
    #    been resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    #   72 s (87.6 s re-measured 24-Aug-2026).  79 resonances, 8 steps, written to
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
    fn, dat = resonantCascadeFile()
    sim = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization, strengths",
                             property=Cascade.ResonantIonizationStrengths(1, 0., 0.),
                             method=Cascade.ProbPropagation(), computationData=[dat] )
    perform(sim; output=false)
    println("\n     CORRECTED 24-Aug-2026.  The two totals differ by a factor of about two, and reading that as")
    println("     the competition between the two fates of a resonance is WRONG: they are sums over almost")
    println("     entirely DIFFERENT resonances.  The 79 fall into two disjoint groups --")
    println("        31 resonances at 313.2 - 324.0 eV carry ALL the ionization strength, 4.164e-04, and almost")
    println("           no recombination strength, 3.21e-07;")
    println("        48 resonances at 275.4 - 286.8 eV cannot ionize at all -- their sequential route is closed --")
    println("           and carry 99.8% of the recombination strength, 1.873e-04.")
    println("     For the resonances that DO ionize, ionization beats recombination by a factor of 1300, not 2.")
    println("     The Auger-over-radiative ordering expected for a K hole in a light ion is therefore confirmed")
    println("     far more strongly than the ratio of totals suggests; what the ratio of totals measures is which")
    println("     THRESHOLD the strength sits on, not which decay wins.  See branch e, where the two are folded")
    println("     with a Maxwellian and the groups separate by temperature.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:  24-Aug-2026
    # Last successful:  24-Aug-2026 ... seconds, once both computations exist.  THE COMPARISON IS NOW A REAL ONE:
    #   Cascade.EiiRateCoefficients turns both channels into the same quantity, so they can simply be divided.
    #   All in cm^3/s, Babushkin:
    #
    #       T [K]     kT [eV]    alpha^res      alpha^EA      res/EA     alpha^DI (Lotz)   indirect/total
    #       1.0e+06     86.17   2.09908e-13   2.46119e-11   0.008528      1.88454e-08         0.0013
    #       3.0e+06    258.52   4.74611e-13   2.04088e-10   0.002326      2.01413e-08         0.0101
    #       1.0e+07    861.73   1.84744e-13   3.55579e-10   0.000520      1.71275e-08         0.0203
    #       3.0e+07   2585.20   4.54890e-14   2.25738e-10   0.000201      1.31974e-08         0.0168
    #
    #   THE CONCLUSION THIS BRANCH USED TO DRAW IS WRONG, and the reason is worth more than the numbers.  It
    #   argued that the resonances "lie between 275 and 324 eV, i.e. far below the 1000 eV at which the
    #   impact-excitation channel was evaluated, so in a plasma cool enough to populate the resonances and not
    #   the excitation, the resonant channel is the only one of the two that contributes at all."  That confuses
    #   THE ENERGY AT WHICH A CROSS SECTION HAPPENED TO BE SAMPLED with THE THRESHOLD OF THE CHANNEL.  1000 eV
    #   was merely the lowest impact energy example-Fi.jl branch a used to carry; the EA channel's threshold is
    #   the 1s -> 2p excitation energy, 297.3 eV, which is INSIDE the resonance band of 275-324 eV.  The two
    #   channels therefore switch on at essentially the same temperature and there is no window in which only
    #   the resonant one contributes.  Once folded, the resonant channel is 0.85% of the EA channel at
    #   kT = 86 eV and falls to 0.02% by kT = 2585 eV.
    #     What survives of the old intuition is only the TREND: the ratio does fall with temperature by a factor
    #     of 42 across this range, so the resonant channel is relatively more important in a cool plasma.  It is
    #     never the dominant one.  Isolated resonances give alpha ~ T^(-3/2) exp(-E/T), which decays once kT
    #     passes 2E/3, while the EA channel keeps integrating a cross section that extends to high energy -- so
    #     the resonant channel is the one that dies away, and it dies from its own T^(-3/2), not from a threshold.
    #
    #   AND THE DIRECT CHANNEL DWARFS BOTH.  Empirical's Lotz rate is added above purely to give the two
    #   indirect channels a scale, and it settles what "the total is a lower bound" is worth for this ion: the
    #   indirect channels together are 0.13% of the total ionization rate at 1e6 K and at most 2% anywhere in
    #   this range.  For Li-like CARBON the indirect contributions are a correction, not a rival.  That is a
    #   statement about Z and not about the method: the EA share of the Li-like sequence grows with nuclear
    #   charge, which is why Arnaud & Rothenflug fit it at all, and carbon sits at the bottom of that sequence.
    #
    # Branch c: AGAINST THE OTHER INDIRECT CHANNEL, on a common footing at last.  Of the three contributions to
    #   electron-impact ionization, this scheme computes two and the direct one is not implemented:
    #
    #     direct                                              Cascade.ImpactIonizationScheme, NOT implemented
    #     impact-excitation with subsequent autoionization     example-Fi.jl branch a, 4063 s at maxKappa = 20
    #     resonant electron capture, sequential/simultaneous   this file branch a, ~90 s
    #
    #   The branch needs BOTH .jld files present in the working directory -- example-Fi.jl branch a for the EA
    #   channel and branch a of this file for the resonant one -- and sorts out which is which by looking for
    #   impact-excitation lines rather than by trusting a file name.  Files that can no longer be loaded are
    #   named and skipped rather than silently ignored: a cascade .jld written before the Radial.Grid refactor
    #   cannot be deserialized at all, and there are usually several of those lying about.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-channelComparison.sum")
    setDefaults("nuclear: charge", 6.)
    #
    temperatures = [1.0e6, 3.0e6, 1.0e7, 3.0e7]
    cands = sort(filter(f -> startswith(f, "zzz-cascade-electron-ionization"), readdir(".")),
                 by = f -> mtime(f), rev = true)
    ## Collected by PUSHING into an array rather than by assigning to variables declared above the loop: at top
    ## level a `for` body that assigns to an outer name creates a new local instead, so the outer one stays
    ## undefined.  push! mutates the container and never rebinds, so it is safe in that scope.
    loaded = Tuple{String,Any}[]
    for  f  in  cands
        try     push!(loaded, (f, JLD2.load(f)))
        catch e
            println(">> skipping $f -- it cannot be loaded: ", first(split(sprint(showerror, e), "\n")))
        end
    end
    hasExc(d) = haskey(d["results"], "impact-excitation lines:")  &&  length(d["results"]["impact-excitation lines:"]) > 0
    iEA  = findfirst(p ->  hasExc(p[2]), loaded)
    iRes = findfirst(p -> !hasExc(p[2]), loaded)
    if  iEA  === nothing  error("No loadable cascade file with impact-excitation lines: run example-Fi.jl branch a.")  end
    if  iRes === nothing  error("No loadable cascade file with resonances: run branch a of this file.")                end
    dEA  = loaded[iEA];    println(">> EA channel from        $(dEA[1])")
    dRes = loaded[iRes];   println(">> resonant channel from  $(dRes[1])")
    #
    simEA  = Cascade.Simulation(Cascade.Simulation(); name="EA channel",
                                property=Cascade.EiiRateCoefficients(1, temperatures, 0., 0.),
                                method=Cascade.ProbPropagation(), computationData=[dEA[2]] )
    aEA    = perform(simEA;  output=true)["data:"]
    simRes = Cascade.Simulation(Cascade.Simulation(); name="resonant channel",
                                property=Cascade.EiiRateCoefficients(1, temperatures, 0., 0.),
                                method=Cascade.ProbPropagation(), computationData=[dRes[2]] )
    aRes   = perform(simRes; output=true)["data:"]
    #
    fac   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    iConf = Configuration("1s^2 2s");   fConf = Configuration("1s^2")
    println("\n  The three channels of electron-impact ionization for Li-like carbon [cm^3/s]:\n")
    println("      T [K]     kT [eV]     alpha^res       alpha^EA      res/EA     alpha^DI(Lotz)   indirect/total")
    for  i = 1:length(temperatures)
        Tau = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperatures[i])
        ## the Lotz rate announces its Gauss-Legendre grid on every call, which would interleave with the table
        aDI = redirect_stdout(devnull) do
                  Empirical.impactIonizationPlasmaAlpha(Distribution.ElectronMaxwell(Tau), iConf, fConf) * fac
              end
        r   = aRes[i].Babushkin;   e = aEA[i].Babushkin
        println("   ", @sprintf("%9.1e   %8.2f   %.5e   %.5e   %8.6f     %.5e      %8.4f",
                temperatures[i], Defaults.convertUnits("energy: from atomic to eV", Tau), r, e,
                e == 0. ? 0. : r/e, aDI, (r + e)/(r + e + aDI)))
    end
    println("\n    alpha^DI is the DIRECT channel from Empirical (Lotz), which the cascade cannot supply; it is")
    println("    shown only to give the two indirect channels a scale.  Read the truncation table printed above")
    println("    for alpha^EA before quoting any ratio at the highest temperature.")
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
    fn, dat = resonantCascadeFile()
    for p in [0.0, pDbl]
        println("\n       ---- dblAugerProbability = " * string(round(p, digits=5)) * " ----")
        sim = Cascade.Simulation(Cascade.Simulation(); name="Resonant ionization, strengths",
                                 property=Cascade.ResonantIonizationStrengths(1, 0., p),
                                 method=Cascade.ProbPropagation(), computationData=[dat] )
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
elseif  false
    # Last successful:  24-Aug-2026
    #   Seconds on top of branch a; no amplitude is recomputed.
    #
    # Branch e: THE MAXWELLIAN FOLD -- alpha^EII (T), which is what an ionization balance actually needs and what
    #   branch c said did not exist.  A resonance strength is energy-integrated and cannot be compared with a cross
    #   section or with anything a plasma model consumes; folding it with a Maxwellian turns it into a rate
    #   coefficient.  Cascade.EiiRateCoefficients does that fold and adds whichever ionization channels the cascade
    #   data carry, naming the ones that are absent instead of quietly dropping them.
    #
    #   WHAT THE PROPERTY REPORTS, and why each part is there:
    #     + alpha(resonant)        from the resonances of this file, through exactly the Boltzmann factor and unit
    #                              conversion that DielectronicRecombination uses, so the two are comparable.
    #     + alpha(exc-autoion)     structurally ZERO here: this computation requested only the resonant processes,
    #                              so there are no impact-excitation lines to fold.  Said out loud rather than
    #                              printed as a bare 0.
    #     + alpha(direct)          NOT AVAILABLE -- there is no Cascade.perform for ImpactIonizationScheme.  For a
    #                              near-neutral target this is usually the LARGEST channel, so the total is a lower
    #                              bound on the ionization rate rather than the ionization rate.
    #     + alpha^DR of the SAME resonances, free of charge, since recombination is the competing fate of the same
    #                              capture and one computation gives both.
    #
    # REPORT (24-Aug-2026, second visit): the DIRECT channel is now IN the table rather than named as absent -- a
    #   semi-empirical Lotz estimate summed over the occupied subshells of 1s^2 2s at Z = 6.  The resonant column is
    #   unchanged to every digit, so nothing about the computed channel moved:
    #        T [K]      kT [eV]     alpha(res)     alpha(direct)    alpha(TOTAL)     res/total
    #       1.0e+05        8.617   2.498263e-26    1.835123e-09    1.835123e-09      1.4e-17
    #       3.0e+05       25.852   2.307238e-16    1.038287e-08    1.038287e-08      2.2e-08
    #       1.0e+06       86.173   2.099078e-13    1.886584e-08    1.886605e-08      1.1e-05
    #       2.5e+06      215.434   4.876418e-13    2.047767e-08    2.047815e-08      2.4e-05
    #       5.0e+06      430.867   3.610688e-13    1.961617e-08    1.961653e-08      1.8e-05
    #       1.0e+07      861.734   1.847436e-13    1.773149e-08    1.773167e-08      1.0e-05
    #       3.0e+07     2585.203   4.548903e-14    1.392171e-08    1.392176e-08      3.3e-06
    #       1.0e+08     8617.344   8.147902e-15    9.888258e-09    9.888266e-09      8.2e-07
    #   WHY A FIT SITS BESIDE COMPUTED COLUMNS.  Without it the TOTAL read 2.1e-13 where the answer is 1.9e-08 --
    #   wrong by a factor of ninety thousand, not by ninety per cent.  An estimate good to tens of per cent is a
    #   strictly better thing to publish than a total that omits the dominant channel, and the printout says at
    #   length which column is which so that nobody mistakes the Lotz number for cascade output.  Its own shape is
    #   right: it peaks near kT = 215 eV, and the 1s subshell contributes 0.1% of it at kT = 86 eV rising to 5.5% at
    #   2585 eV, as kT approaches the ~490 eV K-shell binding energy.
    #
    # REPORT (24-Aug-2026, first visit).  alpha in cm^3/s, Babushkin, WITHOUT the direct channel:
    #        T [K]      kT [eV]     alpha(res)      alpha^DR      ratio
    #       1.0e+05        8.617    2.498e-26      6.539e-25      0.0382
    #       3.0e+05       25.852    2.307e-16      3.724e-16      0.6195
    #       1.0e+06       86.173    2.099e-13      1.275e-13      1.6463
    #       2.5e+06      215.434    4.876e-13      2.303e-13      2.1171
    #       5.0e+06      430.867    3.611e-13      1.568e-13      2.3022
    #       1.0e+07      861.734    1.847e-13      7.695e-14      2.4008
    #       3.0e+07     2585.203    4.549e-14      1.843e-14      2.4688
    #       1.0e+08     8617.344    8.148e-15      3.268e-15      2.4931
    #
    #   THE SHAPE CHECK PASSES, and it is worth being exact about what it does and does not establish, because the
    #   first version of this note overclaimed it.  alpha ~ T^(-3/2) E exp(-E/T) has d(ln alpha)/dT = 0 at kT = 2E/3
    #   EXACTLY; with a strength-weighted mean resonance energy of 318.503 eV the maximum is predicted at
    #   kT = 212.336 eV, and the tabulated maximum falls at 215.434 eV, one of the two grid points bracketing it.
    #     + It DOES test the fold.  A T^(-1/2) in place of T^(-3/2) puts the true maximum at kT = 2E = 637 eV, two
    #       grid points away, and the bracket test fails; so do a sign error in the exponential and a wrong
    #       Kelvin-to-Hartree conversion.
    #     - It does NOT test that the resonances were correctly identified.  Both the prediction and the curve are
    #       built from the SAME energies, so the level-index collision of 22-Aug-2026 -- which put spurious
    #       resonances at 24 eV where the true ones lie at 320 -- would have given E = 24 eV, a prediction of 16 eV,
    #       and a curve peaking obediently at 16 eV.  It would have passed.  What catches that is a READER seeing
    #       "mean resonance energy 24 eV" printed beside a threshold he knows to be 320.  The property prints E in
    #       eV for exactly that reason, and it is a reporting virtue rather than an automatic check.
    #
    #   THE RATIO COLUMN IS THE PHYSICS, and it is what corrected branch b.  It runs from 0.038 at kT = 8.6 eV to a
    #   plateau of 2.49 at high temperature -- a factor of 65 across the range -- which no single "S(ion)/S(rec) is
    #   about two" can express.  The cause is the two-group structure named in branch b: 48 resonances at 275-287 eV
    #   hold 99.8% of the recombination strength and CANNOT ionize, while the 31 at 313-324 eV hold all of the
    #   ionization strength.  A cool plasma populates the lower group preferentially -- at kT = 8.6 eV the Boltzmann
    #   factor favours it by exp(43/8.6) ~ 150 -- so recombination wins by 26 to 1; a hot plasma populates both, the
    #   exponentials cancel, and the ratio saturates at the energy-weighted ratio of the strengths, 2.50.  The
    #   competition between recombination and ionization is therefore a function of TEMPERATURE, not a number.
    #
    setDefaults("print summary: open", "zzz-Cascade-Fl-rateCoefficients.sum")
    #
    fn, dat = resonantCascadeFile()
    temperatures = [1.0e5, 3.0e5, 1.0e6, 2.5e6, 5.0e6, 1.0e7, 3.0e7, 1.0e8]
    ## directCharge = 6 and directConfig = 1s^2 2s add the SEMI-EMPIRICAL Lotz estimate of the direct channel.  It is
    ## not of the same kind as the two computed channels and the printout says so at length -- but without it the
    ## TOTAL column omits 98% or more of the ionization rate for this ion, and a total wrong by a factor of 100 is a
    ## worse thing to publish than an estimate wrong by 30%.
    sim = Cascade.Simulation(Cascade.Simulation(); name="EII rate coefficients, Li-like C",
                             property=Cascade.EiiRateCoefficients(1, temperatures, 0., 0., 6., Configuration("1s^2 2s")),
                             method=Cascade.ProbPropagation(), computationData=[dat] )
    perform(sim; output=true)
    #
    setDefaults("print summary: close", "")
    #
end
