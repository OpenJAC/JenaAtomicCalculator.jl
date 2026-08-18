
println("Od) Elastic electron scattering: the transport cross section of helium against MEASURED values " *
        "(Jablonski, Salvat & Powell 2004, Table 3).")

setDefaults("print summary: open", "zzz-ParticleScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:  17-Aug-2026
    # Last successful:  17-Aug-2026
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
    #   REPORT (rbox = 20 a0), our three interaction models against their two potentials, all in a0^2:
    #
    #       E [eV]   static   Slater   Furness-McCarthy  |  measured     TFD     DHF
    #        50      2.346    3.403    2.814             |  2.51/3.32    2.45    2.77
    #       100      0.8261   1.157    0.9432            |  0.754        0.834   0.928
    #
    #   THE MODELS PAIR OFF, and that is the result. Our Furness-McCarthy exchange -- the term ELSEPA uses, NSRDS-64
    #   Eq. (10) -- reproduces their DHF column to 1.6 % at BOTH energies, and our plain static field reproduces their
    #   cruder TFD column to 1-4 %. Each of our models lands on its counterpart, which is a quantitative agreement with
    #   a reference calculation rather than an order-of-magnitude one, and it validates the amplitudes, the phase-shift
    #   convention, the transport integral and now the exchange term together.
    #
    #   THE SLATER COLUMN IS THE CAUTIONARY ONE. Using JAC's DFS field for the projectile -- which is what this module
    #   did until 17-Aug-2026 -- puts sigma_1 25 % above DHF. That term is built for a BOUND electron and does not know
    #   the impact energy; the Furness-McCarthy term does, and weakens as the projectile becomes fast. It is kept as
    #   ParticleScattering.StaticFieldSlaterExchange() for comparison, but it is no longer the default.
    #
    #   WHAT REMAINS, and it is NOT ours. Against the MEASURED value at 100 eV we are 25 % high -- but so is the
    #   reference: DHF gives 0.928 against 0.754, i.e. 23 % high. Our residual disagreement with experiment is the
    #   reference calculation's own, and its cause is known: neither treatment carries a correlation-polarization
    #   potential, and at 50-100 eV a helium target is strongly polarised by the passing electron. Adding such a term
    #   is the next ELSEPA feature to implement, and it should move BOTH us and them towards the measurement.
    #
    #   Two cautions from the same table, so that the tolerance is not mistaken for slackness. The two measured He
    #   values at 50 eV are 2.51 and 3.32 -- 32 % apart -- so "the measurement" is not one number here. And for mercury
    #   at 100 eV the same table has measured 4.29 against 8.70 (TFD) and 7.70 (DHF), the reference calculations being
    #   a factor ~1.8 above experiment for a heavy target.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 20.0)
    literature = Dict( 50.0 => ("2.51/3.32", "2.45", "2.77"),  100.0 => ("0.754", "0.834", "0.928") )
    models     = [ ParticleScattering.StaticField(), ParticleScattering.StaticFieldSlaterExchange(),
                   ParticleScattering.StaticFieldFurnessMcCarthy() ]
    println("\n   E [eV]   static   Slater   Furness-McCarthy  |  measured     TFD     DHF     [a0^2]")
    for  en in [50.0, 100.0]
        sigmas = Float64[]
        for  model in models
            psSettings = ParticleScattering.Settings(ParticleScattering.Settings(), interaction = model,
                                                     impactEnergies = [en], polarThetas = [1.0], polarPhis = [0.0],
                                                     printBefore = false, epsPartialWave = 1.0e-7, maxL = 120)
            wc = Atomic.Computation(Atomic.Computation(), name="Transport cross section", grid=grid,
                                    nuclearModel    = Nuclear.Model(2.0),
                                    initialConfigs  = [Configuration("1s^2")],
                                    finalConfigs    = [Configuration("1s^2")],
                                    processSettings = psSettings )
            event = perform(wc; output=true)["particle-scattering events:"][1]
            push!( sigmas, event.integrated.sigmaMomentumTransfer )
        end
        m, t, d = literature[en]
        println("   ", rpad(en, 9), rpad(round(sigmas[1], sigdigits=4), 9), rpad(round(sigmas[2], sigdigits=4), 9),
                rpad(round(sigmas[3], sigdigits=4), 18), "|  ", rpad(m, 12), rpad(t, 8), d)
    end
    #
elseif  false
    # Last successful:  unknown ...
    #
end
#
setDefaults("print summary: close", "")
