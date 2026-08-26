println("Dvnew) Apply & test the PhotoRecombinationInterference module: the RADIATIVE (RR) and the")
println("    DIELECTRONIC (DR) recombination amplitudes added COHERENTLY, and the interference between")
println("    them in the cross section, the anisotropy beta_2 and the linear polarization of the emitted")
println("    photon.  Both processes end in the same final state -- a recombined ion plus one photon --")
println("    so the amplitudes must be summed before, not after, they are squared.  JAC could not see")
println("    this before: module-PhotoRecombination.jl gives RR alone, and module-DielectronicRecombination.jl")
println("    is rate-based (S = C Gamma_r/(Gamma_a+Gamma_r)) and so discards the phase that interference")
println("    lives on -- it also has no angular distribution and no polarization at all.")
println("    Physics reference: X.-M. Tong, Phys. Rev. A 107 (2023) 052801, Li-like target and Be-like")
println("    product through the KLL resonances; and Phys. Rev. Lett. 74 (1995) 54 for the original")
println("    observation of the Fano asymmetry in uranium.")

if  true
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    # Branch a: THE RR LIMIT IS EXACT.  With includeDR = false the module must reproduce PhotoRecombination digit for
    #   digit, and not merely closely: `toPhotoRecombinationLine` packages the very same amplitudes into a
    #   PhotoRecombination.Line and hands them to the very same PhotoRecombination.computeCrossSectionForMultipoles and
    #   PhotoRecombination.computeAnisotropyParameter.  So this is a bookkeeping test, and it is the one that fails
    #   first and loudest if a partial wave or a total symmetry has been dropped or double-counted.
    #
    #   RESULT, 20-Aug-2026, Li-like Fe23+ -> Be-like Fe22+ 1s^2 2s^2:
    #      E [eV]     sigma PhotoRec     sigma PRI(RR)      ratio     beta_2 PhotoRec   beta_2 PRI     ratio
    #      2000       3.044576835e-6     3.044576835e-6     1.0       -0.49986513       -0.49986513    1.0
    #      4000       1.334206023e-6     1.334206023e-6     1.0       -0.49983597       -0.49983597    1.0
    #   and sigma_DR = interference = 0 identically.  beta_2 = -0.4999 is itself the textbook value for E1 radiative
    #   recombination into an s shell, so the reference is not merely self-consistent but right.
    setDefaults("print summary: open", "zzz-PhotoRecombinationInterference.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 4.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(26.)
    energies = [2000.0, 4000.0]
    #
    # ---- THE INPUT, spelled out ---------------------------------------------------------------------------------
    #   INITIAL       Li-like Fe23+, ground configuration     1s^2 2s      N = 3 electrons, one level, J^P = 1/2+
    #   INTERMEDIATE  Be-like Fe22+, KLL doubly excited       1s 2s^2 2p   N = 4, a 1s hole plus a 2p electron;
    #                                                                     four levels, J^P = 0-, 1-, 2-, 1-
    #   FINAL         Be-like Fe22+, ground configuration     1s^2 2s^2    N = 4, one level, J^P = 0+
    #   The captured electron is the fourth: Li-like + e --> Be-like, and the KLL label says that a K-shell electron
    #   is promoted to L while the free electron is captured into L as well.
    initialConfigs      = [Configuration("1s^2 2s")]
    intermediateConfigs = [Configuration("1s 2s^2 2p")]
    finalConfigs        = [Configuration("1s^2 2s^2")]
    Basics.displayConfigurations(26., initialConfigs;      sa="INITIAL, Li-like Fe23+:        ")
    Basics.displayConfigurations(26., intermediateConfigs; sa="INTERMEDIATE, Be-like KLL:     ")
    Basics.displayConfigurations(26., finalConfigs;        sa="FINAL, Be-like Fe22+:          ")
    #
    prSet = PhotoRecombination.Settings(PhotoRecombination.Settings(); multipoles=[E1], gauges=[UseCoulomb, UseBabushkin],
                electronEnergies=energies, maxKappa=2, calcAnisotropy=true, printBefore=false)
    cRR = Atomic.Computation(Atomic.Computation(), name="RR reference", grid=grid, nuclearModel=nModel,
              initialConfigs = initialConfigs, finalConfigs = finalConfigs,
              processSettings = prSet)
    rrLines = perform(cRR; output=true)["photo recombination lines:"]
    #
    priSet = PhotoRecombinationInterference.Settings(PhotoRecombinationInterference.Settings();
                multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], electronEnergies=energies, maxKappa=2,
                includeRR=true, includeDR=false, calcAnisotropy=true, calcPolarization=true, printBefore=false)
    cPRI = Atomic.Computation(Atomic.Computation(), name="RR limit of the interference module", grid=grid, nuclearModel=nModel,
               initialConfigs      = initialConfigs,
               intermediateConfigs = intermediateConfigs,
               finalConfigs        = finalConfigs,
               processSettings = priSet)
    priPaths = perform(cPRI; output=true)["photorecombination-interference pathways:"]
    #
    println("\n\n  RR limit:  the interference module against PhotoRecombination\n")
    println("   E_elec [a.u.]     sigma PhotoRec      sigma PRI (RR)      ratio          beta_2 PhotoRec   beta_2 PRI        ratio")
    for  p in priPaths
        for  l in rrLines
            if  abs(l.electronEnergy - p.electronEnergy) < 1.0e-10  &&  l.finalLevel.index == p.finalLevel.index
                b1 = PhotoRecombination.computeAnisotropyParameter(2, l)
                b2 = length(p.betaNu) > 0 ? p.betaNu[1] : EmPropertyC(0.0im)
                println("  " * rpad(round(p.electronEnergy, sigdigits=8), 17) *
                        rpad(round(l.crossSection.Coulomb, sigdigits=10), 20) *
                        rpad(round(p.crossSectionRR.Coulomb, sigdigits=10), 20) *
                        rpad(round(p.crossSectionRR.Coulomb/l.crossSection.Coulomb, digits=12), 15) *
                        rpad(round(real(b1.Coulomb), sigdigits=8), 18) *
                        rpad(round(real(b2.Coulomb), sigdigits=8), 17) *
                        string(round(real(b2.Coulomb)/real(b1.Coulomb), digits=12)))
            end
        end
    end
    println("\n  and with includeDR = false the resonant term must vanish identically:")
    for  p in priPaths
        println("     E = $(round(p.electronEnergy, sigdigits=6))   sigma_DR = $(p.crossSectionDR.Coulomb)   " *
                "interference = $(p.interference.Coulomb)")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch b: THE DR LIMIT, AND WITH IT THE CALIBRATION OF N.  The resonant amplitude carries one overall constant that
    #   the direct amplitude does not, and getting it wrong would leave every interference term wrong by that factor.
    #   `PhotoRecombinationInterference.dielectronicNormalization` DERIVES it -- by requiring the resonant term alone,
    #   integrated over an isolated resonance, to reproduce the standard dielectronic resonance strength, and carrying
    #   JAC's own three normalizations (Gamma_a = 2 pi Sum|A_Auger|^2, A_r = 8 pi alpha omega/(2J_d+1) Sum|A_rad|^2,
    #   sigma_RR = (1/beta gamma^2) 8 pi^3 alpha^3 omega/(2J_f+1) Sum|A_RR|^2) through that condition, whereupon
    #   everything cancels except N = sqrt((2J_f+1)/(2J_i+1)).  This branch does NOT trust that derivation; it measures it.
    #
    #   The test is done at EXACT resonance, which is equivalent to testing the energy integral and costs one point
    #   instead of a scan: for a Lorentzian, integral = peak * pi Gamma/2, so matching the peak matches the integral.
    #   Expected peak, built from JAC's OWN rates for the level in question:
    #        sigma_peak = (4 pi / k^2) * (2J_d+1)/(2(2J_i+1)) * Gamma_a * A_r / Gamma_d^2 ,   k^2 = 2 E_res .
    #
    #   RESULT, 20-Aug-2026.  Li-like Fe23+, intermediate 1s 2s^2 2p level 4 (J = 1-) at E_res = 172.464625 a.u.
    #   = 4693.00 eV, with Gamma_a = 1.392415e-4, A_r = 1.082293e-2, Gamma_d = 1.096217e-2 a.u.:
    #        sigma_DR (computed) = 3.4109213e-4 a.u.
    #        sigma_DR (expected) = 3.4265846e-4 a.u.
    #        ratio               = 0.9954289
    #   The residual 0.46 % is NOT numerical slop and NOT partial-wave truncation -- the ratio is bit-identical for
    #   maxKappa = 3, 4, 5 and 6.  It is the RELATIVISTIC KINEMATICS of the incoming electron: the module carries
    #   1/(beta^2 gamma^2), whereas the textbook Breit-Wigner reference above uses the non-relativistic k^2 = 2E.  The
    #   ratio between the two is k^2/(c^2 beta^2 gamma^2) = 0.9954290 at this energy, against the measured 0.9954289:
    #   dividing one by the other leaves 0.99999989, i.e. agreement to ONE PART IN 10^7.  N is therefore confirmed
    #   exactly, and it is the reference formula, not the module, that is the approximate one.
    setDefaults("print summary: open", "zzz-PhotoRecombinationInterference.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 4.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(26.)
    #
    # ---- THE INPUT, spelled out ---------------------------------------------------------------------------------
    #   INITIAL       Li-like Fe23+, ground configuration     1s^2 2s      N = 3 electrons, one level, J^P = 1/2+
    #   INTERMEDIATE  Be-like Fe22+, KLL doubly excited       1s 2s^2 2p   N = 4, a 1s hole plus a 2p electron;
    #                                                                     four levels, J^P = 0-, 1-, 2-, 1-
    #   FINAL         Be-like Fe22+, ground configuration     1s^2 2s^2    N = 4, one level, J^P = 0+
    #   The captured electron is the fourth: Li-like + e --> Be-like, and the KLL label says that a K-shell electron
    #   is promoted to L while the free electron is captured into L as well.
    initialConfigs      = [Configuration("1s^2 2s")]
    intermediateConfigs = [Configuration("1s 2s^2 2p")]
    finalConfigs        = [Configuration("1s^2 2s^2")]
    Basics.displayConfigurations(26., initialConfigs;      sa="INITIAL, Li-like Fe23+:        ")
    Basics.displayConfigurations(26., intermediateConfigs; sa="INTERMEDIATE, Be-like KLL:     ")
    Basics.displayConfigurations(26., finalConfigs;        sa="FINAL, Be-like Fe22+:          ")
    #
    dIndex = 4;    E_res = 172.464625                       # a.u., from the level structure printed by this same run
    gammaA = 1.392415e-4;   A_r = 1.082293e-2               # a.u., Coulomb gauge
    gammaD = gammaA + A_r
    sigmaPeakExpected = (4pi/(2*E_res)) * (3.0/4.0) * gammaA * A_r / gammaD^2
    wc     = Defaults.getDefaults("speed of light: c")
    gam    = 1.0 + E_res/wc^2;   bg2 = (1.0 - 1.0/gam^2) * gam^2
    relCorr = (2*E_res) / (wc^2 * bg2)                      # k^2 / (c^2 beta^2 gamma^2)
    #
    priSet = PhotoRecombinationInterference.Settings(PhotoRecombinationInterference.Settings();
                multipoles=[E1], gauges=[UseCoulomb, UseBabushkin],
                electronEnergies=[Defaults.convertUnits("energy: from atomic", E_res)], maxKappa=3,
                includeRR=false, includeDR=true, calcAnisotropy=true, calcPolarization=false, printBefore=true,
                pathwaySelection = PathwaySelection(true; indexTriples=[(1, dIndex, 1)]))
    comp = Atomic.Computation(Atomic.Computation(), name="DR limit, calibration of N", grid=grid, nuclearModel=nModel,
               initialConfigs      = initialConfigs,
               intermediateConfigs = intermediateConfigs,
               finalConfigs        = finalConfigs,
               processSettings = priSet)
    paths = perform(comp; output=true)["photorecombination-interference pathways:"]
    #
    println("\n\n  DR limit:  calibration of the resonant normalization N\n")
    for  p in paths
        ratio = p.crossSectionDR.Coulomb / sigmaPeakExpected
        println("     sigma_DR (computed)                     = $(p.crossSectionDR.Coulomb)  a.u.")
        println("     sigma_DR (non-relativistic Breit-Wigner) = $sigmaPeakExpected  a.u.")
        println("     ratio                                   = $ratio")
        println("     relativistic factor k^2/(c^2 b^2 g^2)   = $relCorr")
        println("     ratio / relativistic factor             = $(ratio/relCorr)   <-- this is the number that must be 1")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch c: THE INTERFERENCE ITSELF, Li-like Fe23+ across one KLL resonance.  Both terms on, scanning the electron
    #   energy through the 1s 2s^2 2p level 4 (J = 1-) resonance in units of its own width Gamma_d = 0.2983 eV.  The
    #   pathway selection keeps that ONE resonance, so the profile is not contaminated by its neighbours: this is an
    #   interference study, not a spectrum.
    #
    #   RESULT, 20-Aug-2026, Coulomb gauge.  The interference term behaves exactly as a Fano profile must:
    #      detuning/Gamma     sigma_total     sigma_RR       sigma_DR      interference    beta_2     lin.pol.
    #        -40             9.394374e-7    1.090520e-6    5.315738e-8   -2.042398e-7   -0.463072    0.564018
    #         -4             4.304246e-6    1.087279e-6    5.246290e-6   -2.029322e-6    0.245356   -0.419498
    #          0             3.421791e-4    1.086920e-6    3.410921e-4   +6.7e-11        0.105394   -0.166886
    #         +4             8.368723e-6    1.086560e-6    5.248856e-6   +2.033307e-6   -0.121671    0.172040
    #        +40             1.345043e-6    1.083336e-6    5.341896e-8   +2.082879e-7   -0.477384    0.578091
    #
    #   Three things are worth reading off, and each is a check:
    #   (i)  The interference term passes through ZERO at exact resonance (6.4e-11 against a total of 3.4e-4, i.e. 2e-7
    #        relative) and changes sign there.  That is the signature of the resonant amplitude's phase rotating through
    #        90 degrees, and it is not something that was put in by hand.
    #   (ii) It falls off as 1/(E - E_d), not as 1/(E - E_d)^2: 2.042e-7 at 40 widths against 4.103e-7 at 20 widths is a
    #        ratio of 2.009, while sigma_DR over the same interval falls by 4.003.  So the total returns to sigma_RR only
    #        SLOWLY -- at 40 widths the interference is still 19 % of sigma_RR.  An expectation that it should have died
    #        away by then would be wrong, and it is the amplitude, not the cross section, that decides.
    #   (iii) beta_2 swings from -0.463 (essentially the pure-RR value) through +0.105 at resonance and back, and the
    #        linear polarization with it, from +0.564 to -0.167 and back to +0.578.  That sign-changing swing across the
    #        resonance IS the effect Tong (2023) reports, here in a Li-like target.
    #
    #   ON THE SIGN OF THE ASYMMETRY, which is the one thing here that is NOT independently verified.  The MAGNITUDE of
    #   every column above rests on two calibrations that each reproduce an independent reference (branch a to machine
    #   precision, branch b to one part in 10^7), and neither of those can test a phase, since both square the amplitude.
    #   The sign instead follows from a convention question that had to be settled by hand: `PhotoRecombination` and
    #   `AutoIonization` attach DIFFERENT phase factors to the continuum state, i^l exp(+i phi) against i^l exp(-i phi),
    #   so the capture amplitude belonging to the former is AutoIonization.amplitude * exp(2 i phi) and NOT its conjugate,
    #   the two differing by (-1)^l.  `ElectronCapture.amplitude` does take the conjugate -- correctly for itself, since it
    #   only squares the result -- and following it here would mirror every profile above.  In each of the three cases run
    #   on 20-Aug-2026 the channel was J_t = 1- built on a 1/2+ ion, so only ODD l contributes and the difference happened
    #   to act as a global sign; in a mixed-parity channel it would instead reweight the partial waves against each other.
    #   The sign as printed has NOT been checked against the measured Fano profile of PRL 74 (1995) 54, and that check is
    #   the natural next step for this branch.
    setDefaults("print summary: open", "zzz-PhotoRecombinationInterference.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 4.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(26.)
    #
    # ---- THE INPUT, spelled out ---------------------------------------------------------------------------------
    #   INITIAL       Li-like Fe23+, ground configuration     1s^2 2s      N = 3 electrons, one level, J^P = 1/2+
    #   INTERMEDIATE  Be-like Fe22+, KLL doubly excited       1s 2s^2 2p   N = 4, a 1s hole plus a 2p electron;
    #                                                                     four levels, J^P = 0-, 1-, 2-, 1-
    #   FINAL         Be-like Fe22+, ground configuration     1s^2 2s^2    N = 4, one level, J^P = 0+
    #   The captured electron is the fourth: Li-like + e --> Be-like, and the KLL label says that a K-shell electron
    #   is promoted to L while the free electron is captured into L as well.
    initialConfigs      = [Configuration("1s^2 2s")]
    intermediateConfigs = [Configuration("1s 2s^2 2p")]
    finalConfigs        = [Configuration("1s^2 2s^2")]
    Basics.displayConfigurations(26., initialConfigs;      sa="INITIAL, Li-like Fe23+:        ")
    Basics.displayConfigurations(26., intermediateConfigs; sa="INTERMEDIATE, Be-like KLL:     ")
    Basics.displayConfigurations(26., finalConfigs;        sa="FINAL, Be-like Fe22+:          ")
    #
    E_res_eV  = 4693.0014857
    gammaD_eV = Defaults.convertUnits("energy: from atomic", 1.096217e-2)
    detunings = [-40.0, -20.0, -8.0, -4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 20.0, 40.0]
    energies  = [E_res_eV + d*gammaD_eV for d in detunings]
    #
    priSet = PhotoRecombinationInterference.Settings(PhotoRecombinationInterference.Settings();
                multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], electronEnergies=energies, maxKappa=3,
                includeRR=true, includeDR=true, calcAnisotropy=true, calcPolarization=true, printBefore=false,
                pathwaySelection = PathwaySelection(true; indexTriples=[(1, 4, 1)]))
    comp = Atomic.Computation(Atomic.Computation(), name="RR + DR interference, Fe23+ KLL", grid=grid, nuclearModel=nModel,
               initialConfigs      = initialConfigs,
               intermediateConfigs = intermediateConfigs,
               finalConfigs        = finalConfigs,
               processSettings = priSet)
    paths = perform(comp; output=true)["photorecombination-interference pathways:"]
    #
    println("\n\n  Interference across the Fe23+ KLL resonance, Gamma_d = $(round(gammaD_eV,digits=4)) eV\n")
    println("   detuning/Gamma   sigma_total     sigma_RR        sigma_DR       interference   cos(phase)   beta_2     lin.pol.")
    for  p in paths
        d  = (Defaults.convertUnits("energy: from atomic", p.electronEnergy) - E_res_eV)/gammaD_eV
        st = p.crossSection.Coulomb;   sr = p.crossSectionRR.Coulomb
        sd = p.crossSectionDR.Coulomb; si = p.interference.Coulomb
        cp = (sr > 0. && sd > 0.) ? si/(2*sqrt(sr*sd)) : NaN
        b2 = length(p.betaNu) > 0 ? real(p.betaNu[1].Coulomb) : NaN
        println("  " * rpad(round(d, digits=2), 16) * rpad(round(st, sigdigits=7), 16) * rpad(round(sr, sigdigits=7), 16) *
                rpad(round(sd, sigdigits=7), 15) * rpad(round(si, sigdigits=7), 15) * rpad(round(cp, digits=5), 13) *
                rpad(round(b2, digits=6), 11) * string(round(real(p.linearPolarization.Coulomb), digits=6)))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch d: HIGH Z, where the two processes become comparable.  Li-like Xe51+ (Z = 54), the same KLL structure, now
    #   through the 1s 2s^2 2p level 2 (J = 1-) resonance at 20551.5 eV.  The radiative width has grown by roughly Z^4 --
    #   Gamma_r = 6.92e-2 a.u. here against 8.74e-4 for the corresponding Fe level -- which is exactly the regime in
    #   which DR ceases to be a small resonant decoration on RR.
    #
    #   Rule 12 and the continuum grid pull in OPPOSITE directions here, and the grid is chosen for the harder of the
    #   two: the 2p turning point scales as 1/(Z-1) and wants a small box, but a 20 keV continuum electron has a de
    #   Broglie wavelength of 0.15 a.u. and needs hp small enough to sample it.  Continuum.gridConsistency() enforces
    #   the latter and will stop the run rather than return plausible-looking numbers from an under-sampled orbital --
    #   it did stop an earlier attempt at Z = 92, which is why this branch is Xe and not uranium.
    #
    #   RESULT, 20-Aug-2026, Coulomb gauge:
    #      detuning/Gamma     sigma_total     sigma_RR       sigma_DR      interference    beta_2     lin.pol.
    #        -40             1.192834e-6    1.156404e-6    2.710938e-8   +9.320773e-9   -0.485518    0.586017
    #         -4             3.897887e-6    1.150922e-6    2.659961e-6   +8.700351e-8   -0.004109    0.006151
    #          0             1.741576e-4    1.150315e-6    1.730046e-4   +2.7e-9         0.247659   -0.423991
    #         +4             3.727003e-6    1.149708e-6    2.663334e-6   -8.603918e-8    0.058686   -0.090690
    #        +40             1.162957e-6    1.144270e-6    2.694755e-8   -8.260864e-9   -0.472466    0.573273
    #   beta_2 swings from -0.486 to +0.248 and the polarization from +0.586 to -0.424, a larger swing than at Fe.  The
    #   normalized interference cos(phase) is only ~0.026 here against ~0.43 at Fe, i.e. the two amplitudes are far more
    #   nearly orthogonal in phase at high Z even though they are comparable in magnitude -- which is why the cross
    #   section barely departs from sigma_RR + sigma_DR while the POLARIZATION is transformed.  That dissociation
    #   between the two observables is the practical argument for measuring the polarization rather than the yield.
    setDefaults("print summary: open", "zzz-PhotoRecombinationInterference.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 5.0e-7, h = 3.0e-2, hp = 5.0e-3, rbox = 4.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(54.)
    #
    # ---- THE INPUT, spelled out ---------------------------------------------------------------------------------
    #   INITIAL       Li-like Xe51+, ground configuration     1s^2 2s      N = 3 electrons, one level, J^P = 1/2+
    #   INTERMEDIATE  Be-like Xe50+, KLL doubly excited       1s 2s^2 2p   N = 4, a 1s hole plus a 2p electron;
    #                                                                     four levels, J^P = 0-, 1-, 2-, 1-
    #   FINAL         Be-like Xe50+, ground configuration     1s^2 2s^2    N = 4, one level, J^P = 0+
    #   The captured electron is the fourth: Li-like + e --> Be-like, and the KLL label says that a K-shell electron
    #   is promoted to L while the free electron is captured into L as well.
    initialConfigs      = [Configuration("1s^2 2s")]
    intermediateConfigs = [Configuration("1s 2s^2 2p")]
    finalConfigs        = [Configuration("1s^2 2s^2")]
    Basics.displayConfigurations(54., initialConfigs;      sa="INITIAL, Li-like Xe51+:        ")
    Basics.displayConfigurations(54., intermediateConfigs; sa="INTERMEDIATE, Be-like KLL:     ")
    Basics.displayConfigurations(54., finalConfigs;        sa="FINAL, Be-like Xe50+:          ")
    #
    E_res_eV  = 20551.48
    gammaD_eV = Defaults.convertUnits("energy: from atomic", 7.13171e-2)
    detunings = [-40.0, -20.0, -8.0, -4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 20.0, 40.0]
    energies  = [E_res_eV + d*gammaD_eV for d in detunings]
    #
    priSet = PhotoRecombinationInterference.Settings(PhotoRecombinationInterference.Settings();
                multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], electronEnergies=energies, maxKappa=3,
                includeRR=true, includeDR=true, calcAnisotropy=true, calcPolarization=true, printBefore=false,
                pathwaySelection = PathwaySelection(true; indexTriples=[(1, 2, 1)]))
    comp = Atomic.Computation(Atomic.Computation(), name="RR + DR interference, Xe51+ KLL", grid=grid, nuclearModel=nModel,
               initialConfigs      = initialConfigs,
               intermediateConfigs = intermediateConfigs,
               finalConfigs        = finalConfigs,
               processSettings = priSet)
    paths = perform(comp; output=true)["photorecombination-interference pathways:"]
    #
    println("\n\n  Interference across the Xe51+ KLL resonance, Gamma_d = $(round(gammaD_eV,digits=4)) eV\n")
    println("   detuning/Gamma   sigma_total     sigma_RR        sigma_DR       interference   cos(phase)   beta_2     lin.pol.")
    for  p in paths
        d  = (Defaults.convertUnits("energy: from atomic", p.electronEnergy) - E_res_eV)/gammaD_eV
        st = p.crossSection.Coulomb;   sr = p.crossSectionRR.Coulomb
        sd = p.crossSectionDR.Coulomb; si = p.interference.Coulomb
        cp = (sr > 0. && sd > 0.) ? si/(2*sqrt(sr*sd)) : NaN
        b2 = length(p.betaNu) > 0 ? real(p.betaNu[1].Coulomb) : NaN
        println("  " * rpad(round(d, digits=2), 16) * rpad(round(st, sigdigits=7), 16) * rpad(round(sr, sigdigits=7), 16) *
                rpad(round(sd, sigdigits=7), 15) * rpad(round(si, sigdigits=7), 15) * rpad(round(cp, digits=5), 13) *
                rpad(round(b2, digits=6), 11) * string(round(real(p.linearPolarization.Coulomb), digits=6)))
    end
    #
    setDefaults("print summary: close", "")
    #
end
