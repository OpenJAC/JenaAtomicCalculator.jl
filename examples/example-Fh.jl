
println("Fh) Cascade.ImpactExcitationScheme: electron-impact excitation of hydrogen-like carbon.")

using JLD2, Printf
#
## The settings below are those of example-Dl.jl, which is the VALIDATED electron-impact example (its branch a
## reproduces the analytic Bethe limit to 0.5%). They are not interchangeable with the defaults used elsewhere
## in the F-series -- see the note on the continuum treatment further down.
setDefaults("method: continuum, Galerkin")              ## NOT "asymptotic Coulomb"
setDefaults("method: normalization, Alok")              ## NOT "pure sine"
setDefaults("unit: energy", "eV")

grid = Radial.Grid(Radial.Grid(true); rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)


# REWRITTEN 07-Aug-2026, Stage 3 of the cascade-scheme series. The previous example-Fh.jl was a second
# photoabsorption file (Fd and Fh were duplicates of each other); photoabsorption now lives in example-Fe.jl,
# so nothing was kept.
#
# THE DIRECT PATH WAS DEAD IN SIX INDEPENDENT WAYS, each hidden behind the one before it. All are fixed:
#   1. ImpactExcitation.computeLinesCascade was CALLED by Cascade.computeSteps and never DEFINED anywhere.
#      Written here, modelled on PhotoExcitation.computeLinesCascade.
#   2. That call omitted comp.nuclearModel, which the continuum orbital of the scattered electron needs.
#   3. Cascade.determineSteps built ImpactExcitation.Settings(false, false, false, false, LineSelection()) --
#      five positional arguments of the wrong types against a ten-field struct, i.e. a MethodError before any
#      physics, and it passed no electron energies at all.
#   4. Basics.extractNonrelativisticConfigurations no longer exists; the FromBasis theme of
#      Basics.extractConfigurations replaces it.
#   5. SelfConsistent.performSCF returns a Multiplet, but generateBlocks still used it as a Basis.
#   6. Cascade.Step's `settings` field is a Union that did not list ImpactExcitation.Settings, so a Step for
#      this scheme could not be constructed however it was built.
# A seventh, harmless one is left alone: perform() stores its result under the key "impact-exited lines:"
# (a typo) while the ElectronExcitation dispatcher reads "impact-excitation lines:" -- they could never match.
#
# COMPUTATION versus SIMULATION -- the point this file had to be corrected on. A cascade COMPUTATION produces
# LINES and nothing else; making physics out of lines is the job of a cascade SIMULATION. That is the split of
# Symmetry 13, 520 (2021) and Eur. Phys. J. D 78, 75 (2024), and it was briefly broken here: the rate
# coefficients were first obtained by switching ImpactExcitation's own calcRateCoefficient on INSIDE the
# computation, which put physics in the wrong stage. That is reverted. The chain is now the same as for
# dielectronic and radiative recombination:
#     branch a   COMPUTATION   -> ImpactExcitation.Line's (sigma, Omega, convergence per level pair and energy)
#                              -> written to a .jld file
#     branch c   SIMULATION    -> Cascade.EieRateCoefficients reads those lines and forms the effective
#                                 collision strengths Upsilon(T) and rate coefficients alpha(T)
# Cascade.EieRateCoefficients was added for this (07-Aug-2026); it mirrors Cascade.DrRateCoefficients and
# Cascade.RrRateCoefficients and simply calls the aggregation that ImpactExcitation already provides. Because
# the effective collision strength is obtained by INTERPOLATING the collision strengths over energy, the
# computation must supply three or more well-spread electron energies; the simulation raises an informative
# error if it is handed fewer.
#
# Basics.Eimex() was renamed to Basics.ImpactExc() at the same time, to pair with the existing
# Basics.ImpactExcAuto() for the resonant channel and to match Auger(), Radiative(), PhotoExc() and the rest.
#
# THE NUMERICAL SETTINGS ARE NOT NEGOTIABLE HERE, and this is the main practical lesson of the file:
#   + CONTINUUM METHOD. Use Galerkin. "asymptotic Coulomb", which the rest of the F-series inherited from the
#     old example headers, gives erratic continuum orbitals for valence transitions in low-charge ions; see
#     the extended note in example-Fg.jl, where it produced a temperature-independent recombination rate.
#   + NORMALIZATION. Use "Alok", as example-Dl.jl does. "pure sine" normalizes on the LAST 300 GRID POINTS,
#     and on a fine grid those span much less than one de Broglie wavelength, so a node there inflates the
#     normalization constant and JAC aborts with "enlarge box-size" even when the box is perfectly adequate.
#     Enlarging rbox from 10 to 22 made it worse, not better, which is the signature of that heuristic.
#   + BOX. rbox must hold several wavelengths of the SLOWEST electron; 20 eV needs ~30 a.u.
#   + THE SYSTEM ITSELF is a numerical choice here, not only a physical one. This file first used the
#     2s -> 2p transition of LI-like C, a valence dipole line in a low-charge ion -- exactly the regime
#     example-Dl.jl records as the slowly-converging and expensive one (its near-neutral Mg II case saturates
#     45% BELOW the Bethe limit even when converged). Combining that with low incident energies, which need a
#     large box and hence many grid points per continuum orbital, made a maxKappa of 60 time out after ten
#     minutes. The file now uses H-like C^5+ 1s -> 2p at 1470 and 5878 eV, i.e. the system and the energy
#     regime of example-Dl.jl branch a, which is validated against the analytic Bethe limit to 0.5%.
#   + PARTIAL WAVES. This is the expensive one. example-Dl.jl records that at maxKappa = 30 the collision
#     strength came out FLAT and 19% low "while looking perfectly smooth and well-behaved", because the Bethe
#     logarithm is built from large impact parameters, i.e. from high partial waves; it uses maxKappa = 120.
#     The first run for this file used maxKappa = 3 and gave a cross section ~500x below a Bethe estimate,
#     with the module's own `convergence` field reporting 10.7%. ALWAYS read that column: JAC flags any line
#     whose convergence exceeds 1e-5, and the failure is silent otherwise.
#
# WHAT IS STILL BLOCKED: Cascade.ElectronExcitationScheme, i.e. the RESONANT half.
# Its perform() is a dispatcher that refers to scheme.calcDirect, scheme.calcResonant and scheme.multipoles --
# and the struct has only two fields, `processes` and `electronEnergies`. It therefore fails on its very first
# line, before it can reach the two further problems the plan already recorded: it constructs a
# Cascade.DielectronicCaptureScheme, a type that does not exist (module-ElectronCapture.jl provides
# ElectronCapture instead), and it references DielectronicCapture.Line, a module that does not exist either.
# Making it work needs a decision on the scheme's own field list AND the electron-capture scheme that the plan
# defers to Stage 6, so it is documented here rather than repaired.


if  true
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... 1015 s at maxKappa = 40; 8 lines (4 energies x 2 final levels), written
    #                  to zzz-cascade-impact-excitation-computations-<date>.jld.  Collision strengths, summed
    #                  over the two 2p levels:
    #                      E_in [eV]    735      1470     2939     5878
    #                      Omega       0.1291   0.2032   0.2868   0.3417     (convergence 2.8e-7 .. 2.3e-3)
    #                  Omega RISES with energy, as a dipole-allowed transition must.  The Bethe slope from these
    #                  four points is (0.3417 - 0.1291)/ln(8) = 0.102 against the analytic 0.1233, i.e. ~17%
    #                  low -- which is precisely what example-Dl.jl reports for an under-converged partial-wave
    #                  sum (19% low at maxKappa = 30).  Branch b reaches 0.8% at maxKappa = 80.  Dated as a
    #                  deliberately SEMI-CONVERGED reference: the trend and the magnitude are right and the
    #                  residual is understood and quantified, which is what this scheme is for.
    #
    # Branch a: REFERENCE CASE -- the 1s -> 2p excitation of H-like C^5+ by electron impact, in the energy
    #   regime of example-Dl.jl branch a. This is a strong dipole-allowed transition, so its collision strength
    #   must rise logarithmically with the incident energy in the Bethe fashion; Dl branch a reproduces that
    #   slope to 0.5% for exactly this system, which is what makes it the right reference here.
    #   This branch WRITES ITS LINES TO A .jld FILE, which branches c and d then read; it is the computation
    #   half of the chain described above and has to be run before them.
    setDefaults("print summary: open", "zzz-Cascade-Fh-reference.sum")

    name   = "Electron-impact 1s -> 2p excitation of H-like C"
    scheme = Cascade.ImpactExcitationScheme([Shell("1s")], [Shell("2p")],
                                            [735.0, 1470.0, 2939.0, 5878.0], collect(0:39), 0, 0., 0.)
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... the Bethe slope reproduced to 0.8%; see below.
    #
    #   RESULT (maxKappa = 10, 20, 40, 80; two 2p levels listed separately; 2 energies):
    #       maxKappa   cost     Omega(1470 eV)      Omega(5878 eV)      convergence
    #          10      180 s    0.0560  0.1121      0.0285  0.0572      3e-2 .. 1.3e-1
    #          20      306 s    0.0672  0.1342      0.0769  0.1538      9e-4 .. 2.4e-2
    #          40      590 s    0.0678  0.1355      0.1139  0.2278      1.6e-6 .. 2.3e-3
    #          80      872 s    0.0678  0.1355      0.1252  0.2503      1.6e-6 .. 7.2e-5
    #   Two things to read off it. First, at maxKappa = 10 the collision strength FALLS with energy, which a
    #   dipole-allowed transition must never do -- truncation does not merely lower Omega, it inverts its
    #   energy dependence, and it does so without any error message. Second, convergence is reached at very
    #   different maxKappa for the two energies: 1470 eV is converged by maxKappa = 40 (1.6e-6) and does not
    #   move thereafter, whereas 5878 eV is still creeping at maxKappa = 80 (7.2e-5, just above the 1e-5
    #   criterion). The faster the electron, the more partial waves -- so maxKappa has to be chosen for the
    #   HIGHEST energy in the run, not the typical one.
    #
    #   THE BETHE CHECK, which is what makes this branch "Last successful". For an optically allowed
    #   transition the collision strength rises logarithmically, Omega -> a*ln(E), with the analytic slope
    #       a = 4 * g_i * f / DeltaE[Ry] = 4 * 2 * 0.4162 / 27 = 0.1233   per unit ln(E)
    #   using f(1s-2p) = 0.4162 and DeltaE = 367.5 eV = 27 Ry for H-like carbon.  Summing the two 2p levels at
    #   maxKappa = 80,
    #       Omega(1470 eV) = 0.0678 + 0.1355 = 0.2033 ,   Omega(5878 eV) = 0.1252 + 0.2503 = 0.3755
    #       slope = (0.3755 - 0.2033) / ln(5878/1470) = 0.1722 / 1.3859 = 0.1243
    #   against the predicted 0.1233 -- agreement to 0.8%, on the same footing as example-Dl.jl branch a
    #   (0.5% by the direct route).  The cascade route therefore reproduces the analytic high-energy limit,
    #   which is the strongest statement available for this scheme.  NB the 5878 eV point is not yet fully
    #   converged, so the true slope is marginally steeper still.
    #
    # Branch b: THE PARTIAL-WAVE SUM -- the branch that decides whether anything else in this file is worth
    #   quoting. The same transition is recomputed with a growing maxKappa, and BOTH the collision strength
    #   and the reported `convergence` are tracked against the cost. Per example-Dl.jl the sum converges only
    #   slowly for a dipole-allowed transition, and a truncated sum fails silently -- it returns a smooth,
    #   plausible, and badly low collision strength.
    setDefaults("print summary: open", "zzz-Cascade-Fh-partialwaves.sum")

    for  mk  in  [10, 20, 40, 80]
        println("\n>>> maxKappa = $mk")
        scheme = Cascade.ImpactExcitationScheme([Shell("1s")], [Shell("2p")],
                                                [1470.0, 5878.0], collect(0:mk-1), 0, 0., 0.)
        wa = Cascade.Computation(Cascade.Computation(); name="EIE, maxKappa = $mk",
                                 nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s")] )
        t  = @elapsed (wb = perform(wa; output=true, outputToFile=false))
        ls = wb["impact-exited lines:"]
        for  l  in  ls
            println(">>>   E = ", round(Defaults.convertUnits("energy: from atomic", l.initialElectronEnergy), digits=1),
                    " eV   Omega = ", round(l.collisionStrength, sigdigits=5),
                    "   convergence = ", round(l.convergence, sigdigits=2))
        end
        println(">>>   cost: $(round(t, digits=1)) s")
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... seconds, since no amplitude is recomputed.  For 1s -> 2p_1/2:
    #                      T [K]      5.0e5      1.0e6      2.0e6      5.0e6      1.0e7
    #                      Upsilon   0.04307    0.04311    0.04374    0.04856    0.05795
    #                      alpha     5.20e-14   2.62e-12   1.58e-11   3.99e-11   5.16e-11   cm^3/s
    #                  THREE internal checks pass:
    #                    (i)   the 2p_3/2 values are EXACTLY twice the 2p_1/2 ones (0.08607 against 0.04307),
    #                          i.e. the (2J+1) statistical ratio, to all digits;
    #                    (ii)  Upsilon is bracketed by the sampled Omega, 0.0431 .. 0.1139 -- a thermal average
    #                          cannot lie outside the range of what is averaged (the free check of Dl branch c);
    #                    (iii) alpha turns on as exp(-DeltaE/kT): at 5e5 K, kT = 43 eV against DeltaE = 367 eV,
    #                          and alpha is suppressed by ~1e-4 relative to its high-temperature value.
    #                  CAVEAT, and it must be read before quoting the low-temperature numbers: Upsilon(5e5 K)
    #                  equals Omega(735 eV) to all digits.  That is not physics, it is the constant
    #                  extrapolation below the sampled energy range (see the note on interpolateCS above).  The
    #                  lowest sampled energy sits 367 eV ABOVE threshold, so at kT = 43 eV almost every
    #                  Gauss-Laguerre node falls below the data and Upsilon is pinned at that floor.  Only the
    #                  two highest temperatures here carry real information about Omega.  Sample from near
    #                  threshold if the low-temperature end matters.
    #
    # Branch c: THE SIMULATION -- Cascade.EieRateCoefficients turns the lines of branch a into effective
    #   collision strengths Upsilon(T) and plasma rate coefficients alpha(T). No amplitude is recomputed here;
    #   the branch reads the .jld written by branch a, which is exactly the two-stage pattern of example-Ff.jl
    #   (dielectronic) and example-Fg.jl (radiative recombination).
    setDefaults("print summary: open", "zzz-Cascade-Fh-rates.sum")

    fn    = sort(filter(f -> startswith(f, "zzz-cascade-impact-excitation-"), readdir()), by = f -> stat(f).mtime)[end]
    println(">>> reading the cascade data from  $fn")
    data  = [JLD2.load(fn)]
    temps = [5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7]
    prop  = Cascade.EieRateCoefficients(1, temps, LevelSelection())
    simu  = Cascade.Simulation(Cascade.Simulation(); name="EIE rate coefficients for H-like C",
                               computationData=data, property=prop,
                               settings=Cascade.SimulationSettings(false, false, 0.) )
    println(simu)
    wd = perform(simu; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... 1s -> 2p summed over both fine-structure components [cm^3/s]:
    #                      T [K]      kT [eV]    cascade      Van Regemorter    ratio
    #                      5.0e+05      43.1     1.5593e-13   5.5414e-13        0.281
    #                      1.0e+06      86.2     7.8459e-12   2.7829e-11        0.282
    #                      2.0e+06     172.3     4.7454e-11   1.6636e-10        0.285
    #                      5.0e+06     430.9     1.1975e-10   4.1234e-10        0.290
    #                      1.0e+07     861.7     1.5480e-10   5.5935e-10        0.277
    #                  The ratio is constant to +-2.3% over 1.3 decades in temperature, so the two agree on the
    #                  SHAPE of alpha(T) and differ by a near-constant factor 3.5, the cascade being lower.
    #                  Two things to keep in mind before reading too much into either number.  The constancy is
    #                  partly structural: both expressions carry the same exp(-DeltaE/kT), so this confirms the
    #                  thermal averaging more strongly than it tests Omega itself.  And the offset has at least
    #                  two plausible sources that are not separated here -- the partial-wave sum of branch a is
    #                  ~17% low at maxKappa = 40, and Van Regemorter's effective Gaunt factor is calibrated on
    #                  neutral and near-neutral atoms and is known to overestimate for highly charged ions
    #                  (JAC also warns that it falls back on Slater-screened hydrogenic binding energies here).
    #                  A factor 3.5 is larger than Van Regemorter's own ~30-50% accuracy claim, so this is a
    #                  consistency check on the trend, not a validation of the absolute rate.
    #
    # Branch d: COMPARISON WITH VAN REGEMORTER -- the external check. Empirical.impactExcitationPlasmaAlpha
    #   with a VanRegemorter1962() approximation gives the rate coefficient of an optically allowed transition
    #   from its oscillator strength and an effective Gaunt factor, which is precisely the regime of the
    #   1s -> 2p line used here. Van Regemorter is itself only a ~30-50% estimate, so agreement within a
    #   factor of about two is the most that should be expected -- but a factor of ten, or a wrong slope in
    #   temperature, would point at the partial-wave truncation of branch b rather than at the formula.
    setDefaults("print summary: open", "zzz-Cascade-Fh-vanregemorter.sum")
    setDefaults("nuclear: charge", 6.)

    fn    = sort(filter(f -> startswith(f, "zzz-cascade-impact-excitation-"), readdir()), by = f -> stat(f).mtime)[end]
    data  = [JLD2.load(fn)]
    temps = [5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7]
    prop  = Cascade.EieRateCoefficients(1, temps, LevelSelection())
    simu  = Cascade.Simulation(Cascade.Simulation(); name="EIE rates", computationData=data, property=prop,
                               settings=Cascade.SimulationSettings(false, false, 0.) )
    rates = redirect_stdout(devnull) do;  perform(simu; output=true)  end
    cascadeAlphas = rates["data:"]        ## one RateCoefficients entry per transition
    #
    ## Van Regemorter gives the rate for the whole 1s -> 2p line, so the two fine-structure components of the
    ## cascade have to be summed before comparing.
    fac = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    println("\n  Electron-impact 1s -> 2p excitation of H-like C: cascade versus Van Regemorter [cm^3/s]")
    println("  ---------------------------------------------------------------------------------")
    println("      T [K]         kT [eV]        cascade        Van Regemorter      ratio")
    for  (i, T)  in  enumerate(temps)
        aCas = sum(r.alphas[i]  for r in cascadeAlphas)
        Tau  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", T)
        aEmp = Empirical.impactExcitationPlasmaAlpha(Distribution.ElectronMaxwell(Tau),
                                                     Configuration("1s"), Configuration("2p");
                                                     approx=Empirical.VanRegemorter1962()) * fac
        println("   ", @sprintf("%9.1e     %8.1f     %12.4e     %12.4e     %7.3f",
                                T, T*8.617333e-5, aCas, aEmp, aCas/aEmp))
    end
    println("  ---------------------------------------------------------------------------------")
    setDefaults("print summary: close", "")
    #
end
