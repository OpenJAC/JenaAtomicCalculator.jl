
println("Od) Elastic electron scattering: the transport cross section of helium against MEASURED values " *
        "(Jablonski, Salvat & Powell 2004, Table 3).")

setDefaults("print summary: open", "zzz-ParticleScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:  17-Aug-2026
    # Last successful:  unknown ...
    # THE FIRST COMPARISON OF THIS MODULE WITH MEASURED DATA. Jablonski, Salvat & Powell, J. Phys. Chem. Ref. Data 33,
    # 409 (2004) -- examples/papers/2004.NIST-jablonski-DCS.pdf -- give in their Table 3 the transport cross section of
    # helium derived from MEASURED differential cross sections, together with their own values from two atomic
    # potentials: Thomas-Fermi-Dirac (TFD) and Dirac-Hartree-Fock (DHF). The transport cross section is their sigma_tr,
    # our sigma_1, and the two are the same integral -- their Eq. (15) of NSRDS-64 is
    #
    #        sigma_tr = 2 pi INT (1 - cos th) (dsigma/dOmega) sin th dth
    #
    # so this is a like-for-like comparison and not a rescaled one. All values in units of a0^2.
    #
    #   REPORT:
    #
    #       E [eV]    JAC (DFS)     measured        TFD        DHF
    #        50       3.403         2.51 / 3.32     2.45       2.77
    #       100       1.157         0.754           0.834      0.928
    #
    #   HOW TO READ IT, and the two halves point in different directions.
    #
    #   The ENERGY DEPENDENCE is right to better than 1.5 %. The ratio sigma_1(50 eV)/sigma_1(100 eV) comes out as 2.94
    #   here, against 2.94 for TFD and 2.98 for DHF. Whatever is wrong is not the shape of the cross section, and not
    #   the partial-wave machinery that produces it.
    #
    #   The ABSOLUTE SCALE is high: 25 % above DHF and 39 % above TFD at 100 eV, and 53 % above the measured value.
    #   The likeliest reason is a MISSING PIECE OF PHYSICS rather than a defect, and the paper itself points at it: its
    #   whole subject is how much the choice of atomic potential moves these cross sections, and TFD and DHF already
    #   differ from each other by 11 % at 100 eV. This module uses neither -- it uses JAC's own DFS field with Slater
    #   exchange -- and, more importantly, it has NO correlation-polarization potential at all. At 50-100 eV a helium
    #   target is strongly polarised by the passing electron, the induced dipole attracts it, and omitting that is
    #   known to leave the static-exchange cross section too large. ELSEPA carries such a term as an option; adding one
    #   is the natural next step for this module and is listed as such in the plan.
    #
    #   Two further cautions from the same table, worth keeping in view. The two measured values at 50 eV are 2.51 and
    #   3.32 -- they differ by 32 % from each other, so "the measurement" is not a single number here. And for mercury
    #   the same table has measured 4.29 against calculated 8.70 (TFD) and 7.70 (DHF) at 100 eV, i.e. the reference
    #   calculations are themselves a factor ~1.8 above experiment for a heavy target. Agreement to tens of percent is
    #   the state of this comparison, not a slack tolerance we have chosen.
    #
    #   NO Last successful date is claimed. The energy dependence is verified; the absolute value is not, and will not
    #   be until a polarization potential exists.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 20.0)
    literature = Dict( 50.0 => ("2.51 / 3.32", "2.45", "2.77"),  100.0 => ("0.754", "0.834", "0.928") )
    println("\n   E [eV]    JAC (DFS)     measured        TFD        DHF        [a0^2]")
    for  en in [50.0, 100.0]
        psSettings = ParticleScattering.Settings(ParticleScattering.Settings(),
                                                 impactEnergies = [en], polarThetas = [1.0], polarPhis = [0.0],
                                                 printBefore = false, epsPartialWave = 1.0e-7, maxL = 120)
        wc = Atomic.Computation(Atomic.Computation(), name="Transport cross section", grid=grid,
                                nuclearModel    = Nuclear.Model(2.0),
                                initialConfigs  = [Configuration("1s^2")],
                                finalConfigs    = [Configuration("1s^2")],
                                processSettings = psSettings )
        wd    = perform(wc; output=true)
        event = wd["particle-scattering events:"][1]
        m, t, d = literature[en]
        println("   ", rpad(en, 10), rpad(round(event.integrated.sigmaMomentumTransfer, sigdigits=4), 14),
                rpad(m, 16), rpad(t, 11), d)
    end
    #
elseif  false
    # Last successful:  unknown ...
    #
end
#
setDefaults("print summary: close", "")
