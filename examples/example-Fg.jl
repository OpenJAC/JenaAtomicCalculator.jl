
println("Fg) Cascade.RadiativeRecombinationScheme: radiative recombination of helium-like carbon.")

using JLD2, Printf
#
setDefaults("method: continuum, Galerkin")              ## NOT "asymptotic Coulomb" -- see the note below
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")

## The grid is that of example-Dd.jl branch b, which is the validated radiative-recombination case: hp must be
## small enough that 15*hp stays below the de Broglie wavelength at the highest electron energy used.
grid = Radial.Grid(Radial.Grid(false); rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


# REWRITTEN 06-Aug-2026, second file of Stage 2. The previous content of this file was an expansion-opacity
# script, which belongs to Fk and was copied there first (see the header of example-Ff.jl for the one-place
# shift of the whole F-series). The radiative-recombination material that used to sit in Ff -- a
# RrRateCoefficients run for Li-like iron with capture into n = 2..8 -- was the right SCHEME in the wrong
# file; it is superseded here by a case small enough to iterate on.
#
# RADIATIVE RECOMBINATION is the simplest of the recombination processes: a free electron is captured and a
# photon is emitted in one step, with no intermediate resonance. It is therefore non-resonant and varies
# smoothly with the electron energy -- the exact opposite of the dielectronic recombination of example-Ff.jl,
# and the two together make up the total recombination rate of an ion.
#
# TWO REAL BUGS, both fixed 06-Aug-2026 while writing this file, and the RR rate coefficients could not be
# computed at all before them:
#   (1) Cascade.simulateRrRateCoefficients read the line list as
#           simulation.computationData[1]["results"]["photo-recombination line data:"].linesR
#       but Cascade.Data carries its lines in the field `lines`.  `.linesR` has not existed since that struct
#       was generalised, so EVERY call raised a FieldError.  The path was simply dead.
#   (2) With that repaired the simulation ran and returned exactly zero at every temperature.  The summation
#       filtered the lines with a HARD-WIRED initial symmetry,
#           if  LevelSymmetry(line.initialLevel.J, line.initialLevel.parity) != LevelSymmetry(AngularJ64(1), Basics.plus)
#               continue    end
#       i.e. it kept only J^P = 1+ initial levels -- a leftover from one particular test case.  Every ion with
#       a closed-shell (0+) ground state, which includes the He-like ion here, had all of its lines discarded
#       and produced a rate coefficient of exactly 0.0000e+00 without any warning.  The filter now uses the
#       symmetry taken from the data, which is also the one the table header announces.
# Both live in module-Cascade-inc-simulations.jl.  Note that this scheme does NOT share the hard-wired-
# multipole defect of the photo-excitation and dielectronic-recombination schemes: it passes
# scheme.multipoles to PhotoRecombination.Settings correctly.  Its scheme.minPhotonEnergy, however, is read
# nowhere, and seven "@show" statements remain in module-Cascade-inc-photorecombination.jl.
#
# THE NUMERICAL POINT OF THIS FILE, which turned out to matter more than any of the above. The free-electron
# energies are placed on a Gauss-Legendre grid over [0, maxFreeElectronEnergy] with NoFreeElectronEnergies
# points, and the SAME grid is then used for every temperature. A rate coefficient is an integral of
# sigma(E) E exp(-E/kT) over that grid, so the grid has to resolve the Maxwellian peak at ~kT AND extend
# well beyond it. One grid cannot do both over a wide temperature range, and branch c shows what happens
# when it does not: at 1e7 K the coefficient is wrong by a factor 3 when the grid stops at 500 eV, and at
# 1e5 K no grid tried here converges at all. Treat maxFreeElectronEnergy and NoFreeElectronEnergies as
# physical parameters to be chosen per temperature range, not as technical defaults.
#
# ASYMPTOTIC COULOMB versus GALERKIN -- the lesson of this file, and it applies to the whole F-series.
# This file was first written with setDefaults("method: continuum, asymptotic Coulomb"), copied from the
# headers of the older example files. With that setting alpha^RR came out almost INDEPENDENT of temperature
# (2.60e-15, 2.61e-15, 2.14e-15 at 1e6, 3e6, 1e7 K), which radiative recombination must never be. That sent a
# long hunt through the energy-grid units, the quadrature range and finally a suspected missing 1/E factor in
# PhotoRecombination.computeCrossSectionForMultipoles -- all of them wrong turns. The cause was the continuum
# method alone. Measured on one and the same carbon case, sigma [barn] at E = 20, 100, 500, 2000 eV:
#       asymptotic Coulomb, coarse grid :  1.73    0.52    1.45    0.47      erratic, essentially flat
#       asymptotic Coulomb, fine grid   :  1.46    0.56    1.52    0.48      still erratic
#       Galerkin, fine grid             : 176.6   27.9     3.57    0.364     clean, ~ E^-1.33
# Only the Galerkin form reproduces the E^-1.3 fall-off that example-Dd.jl validated against Stobbe and
# against Ichihara & Eichler, ADNDT (2000). Asymptotic Coulomb matches the orbital to an analytic Coulomb
# form outside a matching radius, which is defensible when the residual potential really is Coulombic where
# the matrix element has weight -- a bare-nucleus K-shell case. For capture into a VALENCE shell of a
# low-charge ion the matrix element draws its weight from the core region, where that assumption fails.
# USE GALERKIN for any cascade scheme with a continuum electron, together with a grid whose hp satisfies
# 15*hp < the de Broglie wavelength at the highest electron energy of the run.
#
# VERDICT, 07-Aug-2026, after switching to Galerkin: the temperature dependence is now RIGHT. Branch d gives
#       T [K]        cascade        empirical (ScaledHydrogenic)    ratio
#       1.0e+06      6.0825e-14     4.4704e-15                      13.6
#       3.0e+06      2.4644e-14     1.3690e-15                      18.0
#       1.0e+07      7.6587e-15     3.2296e-16                      23.7
# i.e. alpha^RR falls as T^-0.90 against the empirical T^-1.14 -- the same shape, where before it was flat.
# The remaining factor of 13-24 is NOT a new defect of this cascade: example-Dd.jl branch a already records
# that JAC's ab-initio radiative recombination exceeds Empirical.ScaledHydrogenic by ~25-32x and
# Empirical.UsingJAC by ~4.3-7.2x for a comparable capture channel, and calls that an open, unresolved case.
# The cascade path now reproduces the behaviour of the validated module path, which is what this file can
# reasonably establish; the offset against the empirical estimates remains that same OPEN question.


if  true
    # Last visit:      07-Aug-2026 ... 2 steps, 12 Rec lines.  With Galerkin orbitals on the finer grid this
    #                  costs more than the ~0.6 s it took with asymptotic Coulomb, but it is the setting that
    #                  gives correct cross sections; see the note above.
    #
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026 ... 3.9 s WARM (40 s cold, the first cascade of a session paying the
    #                  compilation of the whole path).  Two steps, 12 Rec lines, as claimed.  The cross
    #                  sections fall with electron energy as radiative recombination must:
    #                      capture into 2s   94.8 -> 15.4 -> 6.22 -> 3.96  at 99, 229, 399, 530 eV
    #
    #                  TWO CORRECTIONS TO WHAT THIS BRANCH USED TO CLAIM, both found on dating it.
    #                  (i)  "~0.6 s warm" is wrong by a factor of six; the measured warm cost is 3.9 s.  It is
    #                       still the cheapest branch of the file, so the conclusion drawn from it stands, but
    #                       the number did not.
    #                  (ii) THE GAUGE AGREEMENT IS NOT UNIFORM AND WAS NOT MENTIONED AT ALL.  Coulomb against
    #                       Babushkin, over the twelve lines:
    #                           capture into 2s      2.1, 6.5, 10.2, 10.4 %
    #                           capture into 2p     31.9, 20.5, 32.5, 35.0 %   and  39.8, 3.3, 14.1, 15.7 %
    #                       The 2s channel is converged to ten per cent or better; the 2p channels are not,
    #                       reaching 40 %.  That is a property of the continuum orbital rather than of this
    #                       branch, and branch c is the one that measures whether the free-electron grid is
    #                       responsible.  Dated with the residual quantified rather than left blank, in the
    #                       manner of branch a of example-Fi.jl -- but a number taken from the 2p channels
    #                       here carries a 40 % gauge uncertainty and should not be quoted without it.
    #
    # Branch a: REFERENCE AND SMOKE CASE -- radiative recombination of He-like C (1s^2) into the n = 2 shells,
    #   giving Li-like C. Two steps and 12 Rec lines; the cheapest branch of the whole series and an obvious
    #   smoke-test candidate.
    setDefaults("print summary: open", "zzz-Cascade-Fg-reference.sum")

    name   = "Radiative recombination of He-like C into n = 2"
    scheme = Cascade.RadiativeRecombinationScheme([E1], [0,1], 4, 500.0, 0., 0., [Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=false)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... alpha^RR = 6.08e-14, 2.46e-14, 7.66e-15 cm^3/s at 1e6, 3e6, 1e7 K, i.e.
    #                  falling as T^-0.90 against the empirical T^-1.14.  Before the two bug fixes this path
    #                  raised a FieldError, then returned exactly zero; with asymptotic Coulomb it returned
    #                  finite but temperature-independent numbers.  Dated on the SHAPE, not the magnitude --
    #                  see the verdict above for the open offset against the empirical estimates.
    #
    # Branch b: RR PLASMA RATE COEFFICIENTS -- the observable, and the path that both bugs above had blocked.
    #   The grid is chosen deliberately here rather than left at the branch-a default: 32 points up to
    #   5000 eV, which is ~6 kT at the highest temperature requested. See branch c for why.
    setDefaults("print summary: open", "zzz-Cascade-Fg-rates.sum")

    scheme = Cascade.RadiativeRecombinationScheme([E1], [0,1], 32, 5000.0, 0., 0., [Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name="RR of He-like C", nuclearModel=Nuclear.Model(6.),
                                 grid=grid, approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
    wb     = perform(wa; output=true, outputToFile=false)
    #
    temps  = [1.0e6, 3.0e6, 1.0e7]
    prop   = Cascade.RrRateCoefficients(1, temps, [E1], LevelSelection(), Configuration[])
    simu   = Cascade.Simulation(Cascade.Simulation(); name="RR rate coefficients for He-like C",
                                computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                                property=prop, settings=Cascade.SimulationSettings(false, false, 0.) )
    println(simu)
    wd = perform(simu; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026 ... ~4.3 s.  alpha^RR [cm^3/s] against the grid:
    #                      points  E_max      T = 1e5      T = 1e6      T = 1e7
    #                         4     500     1.7000e-15   3.1590e-15   6.4530e-16
    #                        16     500     6.3400e-16   2.4780e-15   6.3240e-16
    #                        24    2000     8.3100e-16   2.4400e-15   1.9050e-15
    #                        32    5000     7.4790e-16   2.6030e-15   2.1440e-15
    #                        48    5000     5.0730e-16   2.4850e-15   2.1400e-15
    #                  kT is 8.6, 86 and 862 eV at the three temperatures.  Reading across: 1e7 K needs
    #                  E_max well above 500 eV (it changes by a factor 3.4), while 1e5 K never settles at all
    #                  -- it wanders over 5.1e-16 to 8.3e-16 with no trend, because a Gauss-Legendre grid on
    #                  [0, 5000 eV] simply has too few nodes near 8.6 eV.  Only the 1e6 K column looks
    #                  converged, and branch d shows even that is suspect.
    #
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026 ... 92 s for the five grids.  alpha^RR [cm^3/s], kT = 8.6, 86 and 862 eV:
    #                      points  E_max     T = 1e5       T = 1e6       T = 1e7
    #                         4     500     3.7307e-14    5.7168e-14    5.5986e-15
    #                        16     500     3.3324e-13    6.2644e-14    5.7797e-15
    #                        24    2000     3.0283e-13    6.1291e-14    7.5374e-15
    #                        32    5000     2.8333e-13    6.0825e-14    7.6587e-15
    #                        48    5000     3.2262e-13    6.1772e-14    7.6880e-15
    #
    #                  BOTH PREDICTED FAILURES ARE VISIBLE, and they are of quite different size.
    #                    * TOO FEW POINTS: at T = 1e5 the four-point grid gives 3.73e-14 against ~3e-13 for
    #                      every larger one -- LOW BY A FACTOR OF NINE.  Four Gauss-Legendre points spread
    #                      over 500 eV cannot resolve a Maxwellian peaked at 8.6 eV; the integrand is missed
    #                      rather than approximated.
    #                    * TOO SHORT A REACH: at T = 1e7 the two 500 eV grids give 5.6-5.8e-15 against
    #                      7.5-7.7e-15 once E_max reaches 2000 eV -- 30 % low, the truncated exponential tail.
    #                  They pull opposite ways, and no single grid here serves all three temperatures: at
    #                  48 points and 5000 eV the two higher temperatures are converged to about 1 %, while
    #                  T = 1e5 still scatters by +-8 % between the last three rows.
    #
    #                  IT ALSO VALIDATES BRANCH b RETROSPECTIVELY, which is the point of running it.  Branch b
    #                  quotes 6.08e-14 at 1e6 K and 7.66e-15 at 1e7 K; those are the 32-point, 5000 eV row to
    #                  the last digit printed.  So branch b's grid was an adequate choice for the temperatures
    #                  it reports, and it is this branch rather than that one which establishes it.
    #
    # Branch c: THE FREE-ELECTRON ENERGY GRID -- the branch that decides whether any number above is worth
    #   quoting. The same rate coefficients are recomputed on a sequence of Gauss-Legendre grids of growing
    #   size and reach. Two separate effects are visible and they pull in opposite directions: too few points
    #   fails to resolve the Maxwellian peak near kT, while too small a maxFreeElectronEnergy truncates the
    #   exponential tail. A single grid has to satisfy both at every temperature requested, which is exactly
    #   what it cannot do over a wide range.
    setDefaults("print summary: open", "zzz-Cascade-Fg-grid.sum")

    temps = [1.0e5, 1.0e6, 1.0e7]
    function rrAlphas(nE, eMax)
        scheme = Cascade.RadiativeRecombinationScheme([E1], [0,1], nE, eMax, 0., 0., [Shell("2s"), Shell("2p")])
        wa = Cascade.Computation(Cascade.Computation(); name="RR", nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
        wb = redirect_stdout(devnull) do;  perform(wa; output=true, outputToFile=false)  end
        prop = Cascade.RrRateCoefficients(1, temps, [E1], LevelSelection(), Configuration[])
        simu = Cascade.Simulation(Cascade.Simulation(); name="RR", property=prop,
                                  computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                                  settings=Cascade.SimulationSettings(false, false, 0.) )
        wd = redirect_stdout(devnull) do;  perform(simu; output=true)  end
        return( [x.Coulomb for x in wd["data:"]] )
    end
    println("\n  RR rate coefficients [cm^3/s] against the free-electron energy grid")
    println("  ---------------------------------------------------------------------------")
    println("    points   E_max [eV]      T = 1e5         T = 1e6         T = 1e7")
    for  (nE, eMax)  in  [(4, 500.0), (16, 500.0), (24, 2000.0), (32, 5000.0), (48, 5000.0)]
        a = rrAlphas(nE, eMax)
        println("   ", @sprintf("%5d   %9.1f    %12.4e    %12.4e    %12.4e", nE, eMax, a[1], a[2], a[3]))
    end
    println("  ---------------------------------------------------------------------------")
    println("    kT is 8.6 eV, 86 eV and 862 eV at these three temperatures.")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      07-Aug-2026
    # Last successful: 07-Aug-2026 ... ~38 s.  Ratios 13.6, 18.0 and 23.7 at 1e6, 3e6 and 1e7 K.  This is the
    #                  branch that first exposed the flat temperature dependence and then confirmed the
    #                  Galerkin fix; the residual offset matches the one already documented in example-Dd.jl.
    #
    # Branch d: COMPARISON WITH AN EMPIRICAL RATE COEFFICIENT -- the counterpart of branch e of example-Ff.jl.
    #   Empirical.photorecombinationPlasmaAlpha with a ScaledHydrogenic() approximation gives the RR rate
    #   coefficient for capture into a given shell, so it can be compared shell by shell with the cascade.
    #   As in Ff the cascade here is truncated -- capture into n = 2 only -- so it must come out below any
    #   total, and the comparison tests the shape and the magnitude rather than an equality.
    setDefaults("print summary: open", "zzz-Cascade-Fg-empirical.sum")
    setDefaults("nuclear: charge", 6.)

    temps  = [1.0e6, 3.0e6, 1.0e7]
    scheme = Cascade.RadiativeRecombinationScheme([E1], [0,1], 32, 5000.0, 0., 0., [Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name="RR of He-like C", nuclearModel=Nuclear.Model(6.),
                                 grid=grid, approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^2")] )
    wb     = redirect_stdout(devnull) do;  perform(wa; output=true, outputToFile=false)  end
    prop   = Cascade.RrRateCoefficients(1, temps, [E1], LevelSelection(), Configuration[])
    simu   = Cascade.Simulation(Cascade.Simulation(); name="RR", property=prop,
                                computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ],
                                settings=Cascade.SimulationSettings(false, false, 0.) )
    wd     = redirect_stdout(devnull) do;  perform(simu; output=true)  end
    #
    fac = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    println("\n  Radiative recombination of He-like C into n = 2 [cm^3/s]")
    println("  ------------------------------------------------------------------")
    println("      T [K]         cascade         empirical        ratio")
    for  (i, T)  in  enumerate(temps)
        Tau  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", T)
        aEmp = Empirical.photorecombinationPlasmaAlpha(Distribution.ElectronMaxwell(Tau),
                                                       Distribution.PhotonVacuumField(0.),
                                                       Configuration("1s^2"), Configuration("1s^2 2s")) * fac
        aCas = wd["data:"][i].Coulomb
        println("   ", @sprintf("%9.1e     %12.4e     %12.4e     %7.3f", T, aCas, aEmp, aCas/aEmp))
    end
    println("  ------------------------------------------------------------------")
    setDefaults("print summary: close", "")
    #
end
