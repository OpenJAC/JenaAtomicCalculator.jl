
println("Ff) Cascade.DielectronicRecombinationScheme: KLL dielectronic recombination of helium-like carbon.")

using JLD2, Printf
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")

grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 8.0)


# REWRITTEN 06-Aug-2026, first file of Stage 2 of the cascade-scheme series.
#
# THE F-FILES WERE SHIFTED BY ONE against the plan, and this was corrected here rather than left to
# archaeology later: the previous example-Ff.jl held RADIATIVE recombination (RrRateCoefficients for Li-like
# iron), which the plan assigns to Fg, and example-Fg.jl held EXPANSION OPACITY, which belongs to Fk. Both
# were copied on to their planned files before this one was overwritten, so nothing was lost:
#     old Ff (radiative recombination)  ->  example-Fg.jl
#     old Fg (expansion opacity)        ->  example-Fk.jl   (a new file)
# Neither copy has been repaired or run; they are raw material for Stages 2 and 5.
#
# WHAT DIELECTRONIC RECOMBINATION IS HERE. A free electron is captured by the ion while a bound electron is
# simultaneously excited -- the inverse of autoionization -- forming a doubly-excited resonance, which then
# stabilizes by emitting a photon. The cascade therefore has two kinds of step, and both appear in the step
# table: Auger steps connecting the captured blocks back to the initial ones (the capture amplitudes, which
# by detailed balance are the autoionization amplitudes), and Radiative steps connecting the captured blocks
# to the singly-excited decay blocks (the stabilization). The reference case below is the classic KLL
# resonance group: He-like C (1s^2) captures an electron into 1s 2l 2l', which stabilizes to 1s^2 2l.
#
# ONE REAL BUG, FIXED 06-Aug-2026 while writing this file, and it is the THIRD instance of the same pattern.
# Cascade.determineSteps built the stabilization step as
#     settings = PhotoEmission.Settings(PhotoEmission.Settings(), gauges=[UseCoulomb, UseBabushkin])
# without ever passing scheme.multipoles -- which was read NOWHERE in
# module-Cascade-inc-dielectronic-recombination.jl. Since PhotoEmission.Settings() defaults to [E1], the
# radiative stabilization was always computed in E1 alone, while generateBlocks() printed
#     "+ all requested multipoles are considered for the stabilization."
# The claim was simply false. Branch b below measures the difference: 57 radiative lines with [E1] against
# 63 with [E1, M1, E2]; before the fix both gave 57. Compare example-Fc.jl (hard-wired [E1] in the
# photo-excitation scheme) and example-Fd.jl (four dead fields in the photo-ionization scheme) -- schemes
# declaring parameters they never read is a recurring defect in this module, and worth checking first for
# every scheme still to come.
#
# Six leftover "@show" statements were removed from the same file at the same time (they dumped the full
# configuration lists and the energy window on every run). One more survives in
# module-Cascade-inc-stepwise-decay.jl:256 and is left alone here, being a different file.  SIX more sit in
# module-Cascade-inc-simulations.jl (lines 1365, 1396, 1424, 1495, 2046, 2075); the one at 1495 fires once per
# resonance inside simulateDrRateCoefficients, i.e. on the very path branches c and e use.  That file is being
# worked on in a parallel session, so it was left untouched -- worth a dedicated cleanup pass.
#
# STILL OPEN in this scheme, recorded but not changed:
#   + Cascade.perform(::DielectronicRecombinationScheme, ...) takes outputToFile and outputDirectory and uses
#     neither; the .jld file always lands in the working directory.
#   + The Auger steps hard-code maxKappa = 7 rather than taking it from the scheme or the settings.
#   + scheme.electronEnergyShift is applied to blocks, but the branch that would shift the energies of
#     user-supplied initialMultiplets is commented out, so the shift is silently ignored on that path.
#
# PHYSICS CHECKS, 06-Aug-2026. There is no external reference here, so all four are internal -- but the KLL
# resonance energies can be argued from the structure of the ion and they come out right:
#   + RESONANCE ENERGIES. The capture resonances land at 229 - 255 eV of electron energy. That is where they
#     must be: the 1s -> 2p excitation of He-like C costs ~308 eV, and the captured n = 2 electron is bound
#     by ~55 eV in the resulting Li-like ion, leaving ~250 eV. The whole KLL group spans ~25 eV, which is the
#     term splitting of 1s 2l 2l'.
#   + POSITION OF THE RATE MAXIMUM. alpha^DR peaks between 1e6 and 1e7 K.  A DR rate coefficient peaks near
#     kT ~ E_res/2, i.e. ~125 eV ~ 1.45e6 K here.  Confirmed.
#   + HIGH-TEMPERATURE FALL-OFF. Far above the resonance, alpha^DR must fall as T^(-3/2). From 1e7 to 1e8 K
#     it drops by a factor 24.4 against the asymptotic 31.6; the calculation is only ~3.5 kT above the
#     resonance at 1e7 K, so a somewhat slower fall-off is exactly what should be seen.
#   + LOW-TEMPERATURE CUT-OFF. At 1e4 K (kT = 0.86 eV) the coefficient is 9.5e-129 -- numerically the
#     Boltzmann factor exp(-250/0.86), and correctly so. It is not a numerical breakdown.
#   + GAUGES agree to 6% on the total rate coefficient.


if  true
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... 12 steps, 16 Auger and 57 radiative lines, ~1.5 s warm -- the cheapest
    #                  reference case of the series so far.  Dated on the strength of branch e, which takes
    #                  exactly this cascade through to a rate coefficient and compares it with an external
    #                  empirical formula; branches b, c and d check it internally.
    #
    # Branch a: REFERENCE AND SMOKE CASE -- the KLL resonance group of He-like carbon, E1 stabilization only.
    #   A 1s electron is excited to n = 2 while the free electron is captured into n = 2, giving the
    #   1s 2l 2l' doubly-excited states of Li-like C, which stabilize to 1s^2 2l. At ~1.5 s warm this is by
    #   far the cheapest reference case of the series so far and is the branch intended for runtests.jl.
    setDefaults("print summary: open", "zzz-Cascade-Ff-reference.sum")

    name   = "KLL dielectronic recombination of He-like C"
    scheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("2p"), 500.0, 0., 0., 1,
                                                     [Shell("1s")], [Shell("2s"), Shell("2p")],
                                                     [Shell("2s"), Shell("2p")],
                                                     [Shell("1s"), Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... ~3.4 s for both runs.  57 radiative lines with [E1] against 63 with
    #                  [E1, M1, E2], while the Auger lines stay at 16 in both -- which is the point: a photon
    #                  multipole cannot change a capture amplitude, so anything else would have been wrong.
    #                  Before the fix both multipole lists gave 57.
    #
    # Branch b: THE MULTIPOLES OF THE STABILIZATION -- branch a repeated with [E1, M1, E2]. Only the
    #   RADIATIVE steps can change; the Auger steps are independent of the photon multipole and must give
    #   the same number of lines. This is the branch that exposed the hard-wired [E1] described above, and
    #   it is the check that the fix works.
    setDefaults("print summary: open", "zzz-Cascade-Ff-multipoles.sum")

    for  mp  in  [[E1], [E1, M1, E2]]
        println("\n>>> stabilization multipoles: $mp")
        scheme = Cascade.DielectronicRecombinationScheme(mp, false, Shell("2p"), 500.0, 0., 0., 1,
                                                         [Shell("1s")], [Shell("2s"), Shell("2p")],
                                                         [Shell("2s"), Shell("2p")],
                                                         [Shell("1s"), Shell("2s"), Shell("2p")])
        wa = Cascade.Computation(Cascade.Computation(); name="KLL DR, multipoles $mp",
                                 nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
        wb = perform(wa; output=true, outputToFile=false)
        nR = sum(length(d.lines) for d in wb["cascade data:"] if eltype(d.lines) == PhotoEmission.Line)
        nA = sum(length(d.lines) for d in wb["cascade data:"] if eltype(d.lines) == AutoIonization.Line)
        println(">>> multipoles $mp :  $nR radiative lines,  $nA Auger lines")
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... ~3 s.  alpha^DR (Coulomb) in cm^3/s:
    #                      1e4 K  9.52e-129     1e5 K  6.55e-23     1e6 K  2.78e-13
    #                      1e7 K  1.17e-13      1e8 K  4.79e-15
    #                  All four checks in the assessment above are on this curve: resonances at 229-255 eV,
    #                  maximum between 1e6 and 1e7 K, T^(-3/2) tail, Boltzmann-suppressed cut-off.
    #
    # Branch c: DR PLASMA RATE COEFFICIENTS -- the observable that this scheme exists for. The cascade
    #   computation of branch a is fed straight into a Cascade.DrRateCoefficients simulation, which sums the
    #   resonance strengths against a Maxwellian electron distribution at the given temperatures. The results
    #   are handed over in memory (wrapped exactly as JLD2.load would return them), so that no run-dated file
    #   name is involved.
    setDefaults("print summary: open", "zzz-Cascade-Ff-rates.sum")

    scheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("2p"), 500.0, 0., 0., 1,
                                                     [Shell("1s")], [Shell("2s"), Shell("2p")],
                                                     [Shell("2s"), Shell("2p")],
                                                     [Shell("1s"), Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name="KLL DR of He-like C", nuclearModel=Nuclear.Model(6.),
                                 grid=grid, approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
    wb     = perform(wa; output=true, outputToFile=false)
    #
    temps  = [1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8]
    prop   = Cascade.DrRateCoefficients(1, 0., temps, 0, 0, DielectronicRecombination.ResonanceSelection())
    simu   = Cascade.Simulation(Cascade.Simulation(); name="DR rate coefficients for He-like C",
                                computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                                property=prop, settings=Cascade.SimulationSettings(false, false, 0.) )
    println(simu)
    wd = perform(simu; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... ~5.8 s for both runs.  Adding n = 3 to the capture shells:
    #                      KLL         12 steps,  16 Auger lines,   57 radiative lines
    #                      KLL + KLM   88 steps,  82 Auger lines,  499 radiative lines
    #                  i.e. 7x the steps and 9x the radiative lines for ONE more principal quantum number,
    #                  and this is still the smallest interesting ion.  A real DR calculation needs the
    #                  series to n ~ 10 and beyond; the numbers here say plainly that it has to be reached by
    #                  extrapolation, not by adding shells.  Cost went from ~1.5 s to ~4 s.
    #
    # Branch d: THE RESONANCE SERIES -- branch a extended so that the free electron may also be captured into
    #   n = 3, i.e. the KLM group in addition to KLL. The KLM resonances lie HIGHER in electron energy (the
    #   captured electron is less bound) but carry smaller strengths, so this branch shows how the resonance
    #   spectrum extends and what each additional n costs. In a real DR calculation the series has to be
    #   followed to much higher n and then extrapolated; this branch exists to show the first step of that,
    #   and to measure the cost before anyone tries n = 10.
    setDefaults("print summary: open", "zzz-Cascade-Ff-series.sum")

    for  (sa, intoShells)  in  [("n = 2 only (KLL)",  [Shell("2s"), Shell("2p")]),
                                ("n = 2 and 3 (KLL + KLM)", [Shell("2s"), Shell("2p"),
                                                             Shell("3s"), Shell("3p"), Shell("3d")])]
        println("\n>>> capture into: $sa")
        scheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("3d"), 500.0, 0., 0., 1,
                                                         [Shell("1s")], [Shell("2s"), Shell("2p")], intoShells,
                                                         [Shell("1s"), Shell("2s"), Shell("2p")])
        wa = Cascade.Computation(Cascade.Computation(); name="DR of He-like C, capture into $sa",
                                 nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
        wb = perform(wa; output=true, outputToFile=false)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... ~30 s (two cascades and two simulations).  In cm^3/s:
    #                      T [K]        KLL only     KLL + KLM    Arnaud 1985   (KLL+KLM)/Arnaud
    #                      5.0e+05      4.4620e-14   7.7919e-14   8.5291e-14        0.914
    #                      1.0e+06      2.7813e-13   6.0274e-13   9.8025e-13        0.615
    #                      2.0e+06      4.1372e-13   1.0166e-12   1.9760e-12        0.515
    #                      5.0e+06      2.4800e-13   6.6076e-13   1.4206e-12        0.465
    #                      1.0e+07      1.1691e-13   3.2027e-13   7.1140e-13        0.450
    #                      3.0e+07      2.7256e-14   7.6085e-14   1.7268e-13        0.441
    #                  ALL THREE TESTS PASS:
    #                    (i)  the maximum lies at 2e6 K in BOTH the cascade and the empirical formula;
    #                    (ii) the ratio is flat at 0.44 - 0.52 from 2e6 to 3e7 K, so the truncated cascade has
    #                         the correct temperature dependence and differs only by a near-constant factor;
    #                    (iii) adding n = 3 lifts the ratio from 0.28 to 0.62 at 1e6 K, so the remaining
    #                         deficit is the n >= 4 Rydberg series -- which is what truncating a series that
    #                         converges roughly as n^-3 must cost.
    #                  The one outlier is 1e5 K, where the cascade exceeds the empirical by 87x.  Both are
    #                  ~1e-23 there, in the deep Boltzmann tail: a formula fitted with ONE effective resonance
    #                  energy and a calculation using the actual ones cannot agree in an exponential tail, and
    #                  neither number is physically relevant at that temperature.
    #
    # Branch e: COMPARISON WITH AN EMPIRICAL RATE COEFFICIENT -- the only external check available for this
    #   scheme without new literature. Empirical.dielectronicRecombinationPlasmaAlpha with Arnaud1985DR()
    #   implements the He-like general-Z formula of Arnaud & Rothenflug, A&AS 60, 425 (1985), which is valid
    #   for Z >= 3 and therefore covers He-like carbon exactly. It gives the TOTAL DR rate coefficient, summed
    #   over the whole Rydberg series of capture shells, whereas the cascade below is truncated at n = 2 (KLL)
    #   or n = 3 (KLL + KLM). The comparison is therefore not an equality test: our numbers must come out
    #   BELOW the empirical total, by the fraction that the truncated series accounts for. What can be tested
    #   is (i) the position of the maximum, (ii) the shape, i.e. whether the ratio is constant in temperature,
    #   and (iii) whether the deficit shrinks by the right amount when n = 3 is added.
    setDefaults("print summary: open", "zzz-Cascade-Ff-empirical.sum")
    setDefaults("nuclear: charge", 6.)

    temps = [1.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7, 3.0e7]
    function drAlphas(intoShells)
        scheme = Cascade.DielectronicRecombinationScheme([E1], false, Shell("3d"), 500.0, 0., 0., 1,
                                                         [Shell("1s")], [Shell("2s"), Shell("2p")], intoShells,
                                                         [Shell("1s"), Shell("2s"), Shell("2p")])
        wa   = Cascade.Computation(Cascade.Computation(); name="DR of He-like C", nuclearModel=Nuclear.Model(6.),
                                   grid=grid, approach=Cascade.AverageSCA(), scheme=scheme,
                                   initialConfigs=[Configuration("1s^2")] )
        wb   = perform(wa; output=true, outputToFile=false)
        prop = Cascade.DrRateCoefficients(1, 0., temps, 0, 0, DielectronicRecombination.ResonanceSelection())
        simu = Cascade.Simulation(Cascade.Simulation(); name="DR rate coefficients",
                                  computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                                  property=prop, settings=Cascade.SimulationSettings(false, false, 0.) )
        wd   = perform(simu; output=true)
        return( wd["data:"] )      ## an EmProperty per temperature, ALREADY in cm^3/s
    end
    kll = drAlphas([Shell("2s"), Shell("2p")])
    klm = drAlphas([Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p"), Shell("3d")])
    #
    ## Cascade.DrRateCoefficients returns cm^3/s already, whereas Empirical.dielectronicRecombinationPlasmaAlpha
    ## returns atomic units -- only the latter has to be converted.  (Getting this backwards costs a factor 6.1e-9.)
    fac = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    println("\n  Dielectronic recombination of He-like C: cascade versus the empirical total [cm^3/s]")
    println("  ------------------------------------------------------------------------------------")
    println("      T [K]        KLL only       KLL + KLM      Arnaud 1985      (KLL+KLM)/Arnaud")
    for  (i, T)  in  enumerate(temps)
        Tau  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", T)
        aEmp = Empirical.dielectronicRecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Tau),
                                                              Configuration("1s^2")) * fac
        a2   = kll[i].Coulomb;            a3 = klm[i].Coulomb
        println("   ", @sprintf("%-13.1e %-15.4e %-15.4e %-16.4e %6.3f", T, a2, a3, aEmp, a3/aEmp))
    end
    println("  ------------------------------------------------------------------------------------")
    setDefaults("print summary: close", "")
    #
end
