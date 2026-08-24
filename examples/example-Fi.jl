
println("Fi) Cascade.ElectronIonizationScheme: excitation-autoionization of lithium-like carbon.")

using JLD2, Printf
#
## As in example-Fh.jl, and for the same reasons: Galerkin continuum orbitals, Alok normalization, a box that
## holds several wavelengths of the slowest electron.  See the note in example-Fg.jl on why "asymptotic Coulomb"
## is not usable here.
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
setDefaults("unit: energy", "eV")

grid = Radial.Grid(Radial.Grid(true); rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)


# WRITTEN 07/08-Aug-2026, Stage 4 of the cascade-scheme series.  The previous example-Fi.jl was a second
# electron-impact excitation file (that scheme now lives in example-Fh.jl) and was discarded entirely.
#
# WHAT THIS SCHEME IS.  Electron-impact ionization of an ion proceeds through two quite different channels:
#     DIRECT     e- + A^q+  ->  A^(q+1)+ + 2e-
#     INDIRECT   e- + A^q+  ->  A^q+*  + e-   ->  A^(q+1)+ + 2e-
# where in the indirect case an INNER-shell electron is first excited to a level above the ionization threshold,
# which then autoionizes.  That indirect channel is excitation-autoionization (EA), and for many ions -- along
# the Li-, Na- and Mg-like sequences in particular -- it rivals or exceeds the direct one.  This file computes
# EA, and only EA.
#
# Cascade.ElectronIonizationScheme was written from nothing for this (module-Cascade-inc-electron-ionization.jl
# had been four comment lines).  It is the MIRROR IMAGE of Cascade.DielectronicRecombinationScheme, and the
# implementation follows that file closely:
#       DR :  Auger (capture, by detailed balance)  +  Radiative (stabilization)     ... net recombination
#       EA :  ImpactExc (excitation)                +  Auger (autoionization)        ... net ionization
# Both therefore build two kinds of Cascade.Step and return two kinds of line.
#
# HOW THE AUTOIONIZING LEVELS ARE SELECTED, which is the one point of substance in determineSteps.  An excited
# configuration only contributes to ionization if it lies ABOVE the ionized one, i.e. if it can autoionize at
# all.  Rather than asking the user to know in advance which inner-shell excited configurations qualify, the
# scheme sets up an Auger step only when
#       determineMeanEnergy(excitedBlock) - determineMeanEnergy(ionizedBlock) > 0
# so the energetics do the selecting.  If NO Auger step is generated, that is the physical answer for the shells
# requested, and the subsequent simulation says so explicitly rather than returning a silent zero.
#
# WHAT IS NOT HERE.  The DIRECT channel: Cascade.ImpactIonizationScheme is reserved for it and is deliberately
# NOT implemented.  The obstacle is shape rather than effort -- JAC's ImpactIonization module is semi-empirical
# (its Settings are an AbstractEmpiricalSettings, it is driven from an Empirical.Computation, it takes a Basis
# rather than two multiplets, and it has no Line type at all), so a cascade built on it would store
# ImpactIonization.CrossSection's per (subshell, impact energy) rather than lines, and the cascade approaches
# would influence only the SCF basis behind the binding energies.  That is a workable design but a different
# one, and it is recorded at the struct rather than guessed at here.  The other indirect contributions are the
# two RESONANT-ELECTRON-CAPTURE channels, in which the incident electron is captured into a doubly-excited
# resonance that then sheds two electrons, sequentially or simultaneously (the REDA and READI of the older
# literature).  So of the contributions to electron-impact ionization, this file covers one.
#
# COMPUTATION AND SIMULATION.  The computation returns LINES: ImpactExcitation.Line's for the excitation and
# AutoIonization.Line's for the decay.  Forming a cross section from them is the job of the simulation
# Cascade.EaCrossSections, added with this scheme:
#       sigma^EA(E) = SUM_e sigma^exc(E; i->e) * B_a(e)
# summed over the autoionizing levels e.  Since this scheme computes no radiative rates, B_a = 1 for every level
# that carries an Auger line, so the result is an UPPER BOUND on the EA cross section -- a good one for
# inner-shell excited levels of light ions, where autoionization beats radiative decay by orders of magnitude,
# but an upper bound nonetheless.  The simulation table says so on every printout.
#
# COST, measured 08-Aug-2026.  The excitation half dominates completely, and it scales with the partial-wave cap
# exactly as in example-Fh.jl.  A minimal case -- one impact energy, one target shell, maxKappa = 4 -- takes
# ~156 s and yields 7 excitation and 7 Auger lines.  Anything with several energies and a realistic maxKappa is
# a job of tens of minutes; see branch b before choosing.


if  true
    # Last visit:      24-Aug-2026
    # Last successful: 24-Aug-2026 ... see the REPORT below; the record that follows is the 08-Aug run, whose
    #                  numbers all still stand.
    #
    # Earlier record (08-Aug-2026) ... 2836 s at maxKappa = 20.  4 steps (2 ImpactExc + 2 Auger); 42 excitation
    #                  lines (3 energies x 7 level pairs x 2 excited blocks) and 14 Auger lines, written to
    #                  zzz-cascade-electron-ionization-computations-<date>.jld.  Both excited blocks are admitted
    #                  by the energetic filter, i.e. 1s2s2p and 1s2s3p both autoionize to He-like C.
    #                  Summed collision strength over the autoionizing levels:
    #                      E [eV]           1000       2000       4000
    #                      Omega           0.37567    0.52092    0.52572
    #                      convergence     0.03%      0.58%      2%
    #                  Omega RISES with impact energy, as a dipole-allowed 1s -> 2p excitation must.  Compare
    #                  the first attempt at maxKappa = 10, which gave 0.3411 / 0.3372 / 0.2219 -- FALLING, and
    #                  with a convergence ten to sixty times worse.  Branch b measures that crossover directly.
    #                  CAVEAT, quantified rather than hidden: the 4000 eV point is still roughly 35% low.
    #                  Branch b finds Omega there still climbing by 37% between maxKappa 20 and 40, which is why
    #                  Omega(2000) and Omega(4000) come out almost equal here instead of continuing to rise.
    #                  Converging that point needs maxKappa well beyond 40 and hours of compute; the trend is
    #                  correct and the residual is known, which is what this ansatz is for.
    #
    # REPORT (24-Aug-2026): re-run with six impact energies, 4063 s at maxKappa = 20.  77 excitation lines and
    #   14 Auger lines.  Not 84 excitation lines: the 320 eV point lies BELOW the 1s -> 3p threshold, so that
    #   block contributes only from 400 eV upwards and gives 5 x 7 lines rather than 6 x 7.  The energetic filter
    #   admits both excited blocks as before.  Branch c confirms that the three original energies reproduce their
    #   08-Aug cross sections to five digits, so this is an extension of the dated result and not a replacement
    #   of it.  A pleasant incidental: the computed 1s -> 2p excitation energy is 297.3 eV against the 297.6 eV
    #   that Arnaud & Rothenflug's fitted I_EA gives for this ion -- 0.1%, from two entirely unrelated routes.
    #
    # Branch a: THE COMPUTATION -- excitation-autoionization of Li-like carbon.  A 1s electron is excited into
    #   2p and 3p; the resulting 1s 2s 2p and 1s 2s 3p levels lie above the He-like C limit and autoionize to
    #   1s^2.  The lines are written to a .jld file that branch c reads.  Li-like carbon is chosen because the
    #   semi-empirical formula of Arnaud & Rothenflug covers the Li-like sequence explicitly, which gives
    #   branch d an external reference.
    setDefaults("print summary: open", "zzz-Cascade-Fi-computation.sum")

    name   = "Excitation-autoionization of Li-like C"
    ## maxKappa = 20, which branch b establishes as the MINIMUM at which this cascade is even qualitatively
    ## right: at 10 the collision strength falls with impact energy, which a dipole-allowed 1s -> 2p excitation
    ## cannot do.  At 20 the 1000 eV point is converged to 0.3%, while the 4000 eV point is still ~35% low --
    ## converging it would need maxKappa well beyond 40 and hours of compute.  Semi-converged on purpose: the
    ## trend is correct and the residual is quantified in branch b, which is what this ansatz is for.
    ## SIX energies since 24-Aug-2026, not three.  The three original points all lie at 1000 eV and above, i.e.
    ## at 3.4x the ~297 eV threshold and beyond -- and branch d, which folds these cross sections into a rate
    ## coefficient, then integrates over an unsampled rise.  Worse, the MAXIMUM of sigma^EA turns out to sit at
    ## 400 eV, outside the old range altogether, 32% above the 1000 eV value.  Restoring the three original
    ## energies reproduces the 08-Aug results exactly (see branch c), so nothing is lost by carrying all six.
    scheme = Cascade.ElectronIonizationScheme([320.0, 400.0, 600.0, 1000.0, 2000.0, 4000.0],
                                              [Shell("1s")], [Shell("2p"), Shell("3p")],
                                              collect(0:19), 1, 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      08-Aug-2026
    # Last successful: 08-Aug-2026 ... ~51 min for the three caps.  Summed Omega over the 1s 2s 2p levels, one
    #                  excited block, two impact energies:
    #                      maxKappa   Omega(1000 eV)   Omega(4000 eV)   trend      cost
    #                         10         0.28998          0.18468       FALLING     522 s
    #                         20         0.32136          0.44574       rising      973 s
    #                         40         0.32221          0.61116       rising     1591 s
    #                  THE TREND FLIPS BETWEEN 10 AND 20, and that is the result worth keeping.  At maxKappa = 10
    #                  the collision strength DECREASES with impact energy, which a dipole-allowed 1s -> 2p
    #                  excitation cannot do; at 20 it increases.  A truncated partial-wave sum therefore does not
    #                  merely lower Omega, it inverts its energy dependence -- and it does so silently, since each
    #                  individual number looks perfectly smooth.  maxKappa = 20 is the minimum at which this
    #                  cascade is even qualitatively right.
    #                  CONVERGENCE IS SET BY THE HIGHEST ENERGY, not the typical one.  1000 eV is converged at
    #                  maxKappa = 20 (+0.3% on going to 40), while 4000 eV is still climbing by 37% between 20
    #                  and 40 and is clearly not converged at either.  This is the same behaviour that
    #                  example-Fh.jl branch b measured for H-like carbon, and for the same reason: the Bethe
    #                  logarithm is built from large impact parameters, i.e. from high partial waves, and the
    #                  faster the electron the more of them contribute.
    #                  Dated on the measurement itself: the convergence behaviour and its cost are established
    #                  and internally consistent, whatever the absolute Omega at 4000 eV eventually turns out
    #                  to be.
    #
    # Branch b: COST AND THE PARTIAL-WAVE CAP -- the same cascade at a growing maxKappa, timing each.  The
    #   excitation half of an EA cascade is an ordinary electron-impact excitation and inherits its convergence
    #   behaviour exactly: a truncated partial-wave sum does not merely lower the cross section, it can invert
    #   its energy dependence, silently (see example-Fh.jl branch b).  This branch exists to fix an affordable
    #   cap before anything else in the file is quoted, and to record what that cap costs.
    setDefaults("print summary: open", "zzz-Cascade-Fi-cost.sum")

    ## TWO energies, not one: the decisive question is not whether Omega has settled at a given energy but
    ## whether it RISES between them, as a dipole-allowed excitation must.  A single energy cannot show that,
    ## and a truncated sum looks perfectly smooth at any one energy (see example-Fh.jl branch b).
    for  mk  in  [10, 20, 40]
        println("\n>>> maxKappa = $mk")
        scheme = Cascade.ElectronIonizationScheme([1000.0, 4000.0], [Shell("1s")], [Shell("2p")], collect(0:mk-1), 1, 0.)
        wa = Cascade.Computation(Cascade.Computation(); name="EA of Li-like C, maxKappa = $mk",
                                 nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2 2s")] )
        t  = @elapsed (wb = perform(wa; output=true, outputToFile=false))
        ls = wb["impact-excitation lines:"]
        for  l  in  ls
            println(">>>   E = ", round(Defaults.convertUnits("energy: from atomic", l.initialElectronEnergy), digits=0),
                    " eV   Omega = ", round(l.collisionStrength, sigdigits=5),
                    "   convergence = ", round(l.convergence, sigdigits=2))
        end
        ls2 = wb["impact-excitation lines:"]
        for  en  in  sort(unique([l.initialElectronEnergy for l in ls2]))
            om = sum(l.collisionStrength for l in ls2 if abs(l.initialElectronEnergy-en) < 1.0e-8)
            println(">>>   summed Omega at ", round(Defaults.convertUnits("energy: from atomic", en), digits=0),
                    " eV = ", round(om, sigdigits=5))
        end
        println(">>>   ", length(wb["autoionization lines:"]), " Auger lines;  cost $(round(t, digits=1)) s")
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      24-Aug-2026
    # Last successful: 24-Aug-2026 ... seconds, since no amplitude is recomputed.  With branch a's six energies:
    #                      E [eV]         320         400         600        1000        2000        4000
    #                      sigma^EA    2.8220e5    2.9691e5    2.6031e5    2.2459e5    1.5557e5    7.8347e4   barn
    #                  A REGRESSION CHECK ACROSS SIXTEEN DAYS, and it passes: the three energies that were
    #                  computed on 08-Aug come back as 2.2459e5, 1.5557e5 and 7.8347e4 against the recorded
    #                  2.2460e5, 1.5557e5 and 7.8350e4 -- five digits, after a full re-run of branch a on a
    #                  different day with three extra energies in the list.
    #                  AND THE NEW POINTS CHANGE THE PICTURE.  sigma^EA has its MAXIMUM at 400 eV, about 100 eV
    #                  above the 297 eV threshold, and the peak value is 32% above the 1000 eV one.  Every
    #                  energy in the original list lay on the falling side, so the maximum of the cross section
    #                  was outside the sampled range -- invisible in the old table, and precisely what a
    #                  Maxwellian fold needs most.  The delayed maximum itself is expected: the channel opens at
    #                  threshold with vanishing phase space and the dipole-allowed 1s -> 2p excitation then
    #                  carries it to a peak before the 1/E fall sets in.
    #
    # Superseded record (08-Aug-2026), three energies only:
    #                      E [eV]        1000        2000        4000
    #                      sigma^EA    2.2460e5    1.5557e5    7.8350e4     barn
    #                  TWO checks pass.  (i) The machinery: the simulation identifies 14 distinct autoionizing
    #                  levels from the Auger lines and sums the excitation cross sections over exactly those,
    #                  ignoring excited levels that carry no Auger line -- which is the selection this property
    #                  exists to make.  (ii) Internal consistency: sigma should scale as Omega/E, which predicts
    #                  sigma(1000)/sigma(4000) = (0.37567/0.52572) * 4 = 2.858 against the measured 2.867, i.e.
    #                  agreement to 0.3%.  The cross sections and the collision strengths of branch a are
    #                  therefore mutually consistent, and sigma falls by only 2.87 over a factor 4 in energy --
    #                  more slowly than 1/E, which is the ln(E)/E behaviour a dipole-allowed inner-shell
    #                  excitation requires.  At maxKappa = 10 it fell by 6.2, i.e. faster than 1/E, which is
    #                  impossible; that was the truncation.
    #                  The absolute values remain an UPPER BOUND (branching ratio 1, no radiative decay computed)
    #                  and inherit the ~35% deficit of branch a at 4000 eV.
    #
    # Branch c: THE SIMULATION -- Cascade.EaCrossSections forms sigma^EA(E) from the lines of branch a.  No
    #   amplitude is recomputed; the branch reads the .jld, which is the two-stage pattern of example-Ff.jl,
    #   example-Fg.jl and example-Fh.jl.
    setDefaults("print summary: open", "zzz-Cascade-Fi-simulation.sum")
    setDefaults("unit: cross section", "barn")

    fn   = sort(filter(f -> startswith(f, "zzz-cascade-electron-ionization-"), readdir()), by = f -> stat(f).mtime)[end]
    println(">>> reading the cascade data from  $fn")
    data = [JLD2.load(fn)]
    prop = Cascade.EaCrossSections(1, Float64[])       ## empty list: report at every computed energy
    simu = Cascade.Simulation(Cascade.Simulation(); name="EA cross sections for Li-like C",
                              computationData=data, property=prop,
                              settings=Cascade.SimulationSettings(false, false, 0.) )
    println(simu)
    wd = perform(simu; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      08-Aug-2026 ... the empirical side alone, seconds.  alpha^EA for Li-like C [cm^3/s]:
    #                      1.0e+06 K   (kT =   86 eV)   2.2406e-11
    #                      3.0e+06 K   (kT =  259 eV)   1.5966e-10
    #                      1.0e+07 K   (kT =  862 eV)   2.8943e-10
    #                      3.0e+07 K   (kT = 2585 eV)   3.0973e-10
    #                  Rising steeply and then saturating, which is what an EA rate coefficient must do: the
    #                  Boltzmann factor exp(-E_exc/kT) turns the channel on around kT ~ E_exc/3 and the
    #                  coefficient then flattens.  Kept at "Last visit" because this is only the REFERENCE half
    #                  of the comparison; the cascade side needs branch a, and the two are not yet placed side
    #                  by side (a cross section against a rate coefficient is not a point-by-point comparison).
    #
    # Last visit:      23-Aug-2026 ... 4.8 s, nearly all of it compilation; the empirical evaluation itself is
    #                  sub-second.  alpha^EA = 2.2406e-11, 1.5966e-10, 2.8943e-10, 3.0973e-10 cm^3/s at
    #                  1e6, 3e6, 1e7 and 3e7 K.
    # Last visit:      24-Aug-2026
    # Last successful: 24-Aug-2026 ... THE COMPARISON IS NOW MADE.  Cascade.EiiRateCoefficients supplies the
    #                  Maxwellian fold that was missing, so the cascade's sigma^EA(E) becomes a rate coefficient
    #                  and can be put beside Arnaud & Rothenflug directly.  All in cm^3/s:
    #
    #                      T [K]      kT [eV]    alpha^EA(JAC)   alpha^EA(A&R)   ratio   weight above E_max
    #                      1.0e+06      86.17     2.46119e-11     2.24060e-11    1.098        0.0000
    #                      3.0e+06     258.52     2.04088e-10     1.59657e-10    1.278        0.0000
    #                      1.0e+07     861.73     3.55579e-10     2.89426e-10    1.229        0.0544
    #                      3.0e+07    2585.20     2.25738e-10     3.09731e-10    0.729        0.5421
    #
    #                  AGREEMENT TO 10-28% WHERE THE INTEGRAL IS COMPLETE, and the one point that goes the
    #                  other way is explained by the last column rather than by physics.  That column is the
    #                  share of the Maxwellian weight lying above the largest impact energy computed, which the
    #                  property prints for exactly this purpose: at 3e7 K, 54% of the weight sits beyond 4000 eV
    #                  where sigma^EA is taken as zero, so the cascade value there is a LOWER bound and the
    #                  property says so unprompted.  At 1e6 and 3e6 K the truncation is nil.
    #                  THE SIGN OF THE RESIDUAL IS ALSO RIGHT.  Where the integral is complete the cascade lies
    #                  ABOVE the fit, which is what B_a = 1 predicts: this scheme computes no radiative rates,
    #                  so every autoionizing level is given branching ratio 1 and sigma^EA is an upper bound.
    #                  Working against that, the cascade covers only the 1s -> 2p and 1s -> 3p channels
    #                  requested in branch a while the fit is to the whole EA contribution -- so the two errors
    #                  have opposite signs and 10-28% is the residual of their competition, not of either alone.
    #                  A regression check came free: the four A&R values reproduce those recorded on 23-Aug-2026
    #                  to every printed digit, so Empirical has not drifted under the Cascade work.
    #
    #                  WHAT THE COMPARISON COST.  Branch a had to be re-run with SIX impact energies rather
    #                  than three -- 320, 400, 600, 1000, 2000, 4000 eV -- because a rate coefficient is an
    #                  integral over energy and A&R's own threshold for this ion is I_EA = 297.6 eV, so the
    #                  three original points at 1000 eV and above sat at 3.4x threshold or higher and left the
    #                  entire rise unsampled.  4063 s at maxKappa = 20; 77 excitation lines, not 84, because
    #                  the 320 eV point lies BELOW the 1s -> 3p threshold and that block is closed there.
    #                  Incidentally the computed 1s -> 2p excitation energy comes out at 297.3 eV against
    #                  A&R's fitted I_EA of 297.6 eV, agreeing to 0.1%.
    #
    # Branch d: THE CASCADE AGAINST THE SEMI-EMPIRICAL FIT, on the same quantity and in the same units.
    #   Empirical.excitationAutoionizationPlasmaAlpha with Arnaud1985EA() gives the EA contribution to the
    #   ionization RATE COEFFICIENT for the Li-like sequence, with an explicit correction factor for carbon;
    #   Cascade.EiiRateCoefficients folds this file's own sigma^EA(E) into the same quantity.  Both are
    #   printed side by side, with the truncation diagnostic that says where the fold may be trusted.
    #   READ THE TWO CAVEATS BEFORE THE RATIO.  The cascade value is an UPPER bound, since this scheme
    #   computes no radiative rates and every autoionizing level is therefore given B_a = 1; and it covers
    #   only the 1s -> 2p and 1s -> 3p channels requested in branch a, where the fit is to the whole EA
    #   contribution.  Those two pull in opposite directions, so the residual is the net of them.
    #   Branch a must have been run with the six-energy list before this branch will say anything useful; with
    #   the original three energies at 1000 eV and above, the integral misses everything below 1000 eV.
    setDefaults("print summary: open", "zzz-Cascade-Fi-empirical.sum")
    setDefaults("nuclear: charge", 6.)

    fac   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    temps = [1.0e6, 3.0e6, 1.0e7, 3.0e7]

    ## The cascade side.  perform() wraps what the property returns under the key "data:".
    fn  = filter(f -> startswith(f, "zzz-cascade-electron-ionization"), readdir("."))
    if  isempty(fn)   error("Run branch a first; no cascade .jld file found in the working directory.")   end
    fn  = joinpath(".", sort(fn, by = f -> mtime(f))[end])
    println(">> reading cascade data from  $fn")
    sim = Cascade.Simulation(Cascade.Simulation(); name="EA rate coefficients of Li-like C",
                             property=Cascade.EiiRateCoefficients(1, temps, 0., 0.),
                             method=Cascade.ProbPropagation(), computationData=[JLD2.load(fn)] )
    aJac = perform(sim; output=true)["data:"]

    println("\n  Excitation-autoionization of Li-like C: cascade against Arnaud & Rothenflug [cm^3/s]")
    println("  ----------------------------------------------------------------------------------")
    println("      T [K]      kT [eV]     alpha^EA (cascade)   alpha^EA (A&R 1985)    ratio")
    for  i = 1:length(temps)
        Tau  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temps[i])
        aEmp = Empirical.excitationAutoionizationPlasmaAlpha(Distribution.ElectronMaxwell(Tau),
                                                             Configuration("1s^2 2s")) * fac
        aCas = aJac[i].Babushkin
        println("   ", @sprintf("%9.1e   %8.1f       %12.5e        %12.5e     %7.3f",
                                temps[i], Defaults.convertUnits("energy: from atomic to eV", Tau), aCas, aEmp,
                                aEmp == 0. ? 0. : aCas/aEmp))
    end
    println("  ----------------------------------------------------------------------------------")
    println("    The TRUNCATION table printed by the property above says where this may be read: where the")
    println("    share of the Maxwellian weight above the largest computed impact energy is not small, the")
    println("    cascade column is a lower bound and the ratio must not be quoted.")
    setDefaults("print summary: close", "")
    #
end
