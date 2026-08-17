
println("Oc) Elastic electron scattering: the Sherman function and its growth with the nuclear charge -- the " *
        "observable that a single-kappa treatment cannot produce at all.")

setDefaults("print summary: open", "zzz-ParticleScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:  17-Aug-2026
    # Last successful:  17-Aug-2026
    # THE PHYSICS. The Sherman function S(theta) is the analysing power of Mott polarimetry: the left-right asymmetry
    # with which a transversely polarised electron beam scatters. It is built from the INTERFERENCE of the direct and
    # spin-flip amplitudes,
    #
    #        S = i (f g* - f* g) / (|f|^2 + |g|^2) ,
    #
    # and g(theta) is a difference of the two spin-orbit partners,  exp(2i delta^-) - exp(2i delta^+). A treatment that
    # computes only kappa = -l-1, as this module did before 17-Aug-2026, gives the two partners the SAME phase shift,
    # so g vanishes identically and S is zero everywhere. Nothing about the spin can be said in that case. The numbers
    # below are therefore not an improved result but a new observable.
    #
    # WHAT IS COMPARED WITH THE LITERATURE. The amplitudes are those of the NIST Electron Elastic-Scattering
    # Cross-Section Database (NSRDS-64, Jablonski, Salvat & Powell 2016; examples/papers/2026.NIST-NSRDS-64.pdf), whose
    # Eqs. (1)-(3) and (13) this module reproduces expression for expression:
    #
    #     NIST Eq. (1)   dsigma/dOmega = |f|^2 + |g|^2
    #     NIST Eq. (2)   f(th) = 1/(2iK) SUM_l {(l+1)[exp(2i d_l^+) - 1] + l[exp(2i d_l^-) - 1]} P_l(cos th)
    #     NIST Eq. (3)   g(th) = 1/(2iK) SUM_l [exp(2i d_l^-) - exp(2i d_l^+)] P_l^1(cos th)
    #     NIST Eq. (13)  K = sqrt(E (E + 2 m c^2)) / c
    #     NIST Eq. (15)  sigma_tr = 2 pi INT (1 - cos th) (dsigma/dOmega) sin th dth      [our sigma_1]
    #
    # with d_l^+ = delta(kappa = -l-1), j = l+1/2, and d_l^- = delta(kappa = +l), j = l-1/2.
    #
    # ONE DELIBERATE DIFFERENCE, and it is the point of doing this in JAC: NSRDS-64 builds its potential from Desclaux's
    # tabulated Dirac-Hartree-Fock density with the Furness-McCarthy local exchange, whereas this module uses JAC's own
    # DFS field (Slater exchange) for the same atom. Absolute cross sections may therefore differ by the amount that
    # choice of potential is worth; the guide notes that transport cross sections are far less sensitive to it than the
    # elastic total, which our own grid study reproduces (see below).
    #
    #   REPORT (500 eV electrons, plane wave, static field with exchange, rbox = 20 a0):
    #
    #       target   Z     max |S|      at theta      sigma_el [a.u.]
    #       He        2    4.157e-4     116.9 deg      1.6868
    #       Ne       10    4.427e-3     110.9 deg      6.2805
    #       Ar       18    1.319e-2     113.9 deg     14.807
    #       Kr       36    1.647e-1     119.9 deg     19.051
    #
    #   Two things to read off. The Sherman function grows by nearly THREE ORDERS OF MAGNITUDE from He to Kr while the
    #   elastic cross section grows only by one: S measures the spin-orbit term of the scattering potential, which rises
    #   far faster with Z than the potential itself. This is exactly why Mott polarimeters are built with gold foils and
    #   never with light gases. And the maximum sits near 110-120 deg for every target, which is where Mott detectors
    #   are in fact placed.
    #
    #   GRID, calibrated against the reference. NSRDS-64 states the ranges rmax at which it truncates the potential:
    #   12.27 a0 for Z = 1 and 17-21 a0 up to Z = 96. Our own convergence study for e + He at 500 eV gives
    #
    #       rbox      sigma_el        sigma_1         DCS(180 deg)
    #       10        1.682109        0.07436136      8.410941e-4
    #       15        1.686657        0.07436101      8.387165e-4
    #       20        1.686769        0.074361        8.463693e-4
    #       25        1.686772        0.074361        8.463840e-4
    #
    #   so rbox = 20 is converged and rbox = 10 costs 0.27 % in sigma_el. It also confirms a statement of the guide:
    #   sigma_1 is converged to SEVEN digits already at rbox = 10, i.e. it is about a hundred times less sensitive to
    #   the potential than sigma_el -- Jablonski and Powell say the transport cross section "depend[s] very weakly on
    #   the potential", and that is what we measure.
    #
    #   NOT claimed: agreement with a measured Sherman function or DCS. NSRDS-64 is a user's guide and carries no
    #   numerical tables; the database itself would be needed for that.
    thetas = [t for t in 0.05 : (pi/60) : pi]
    println("\n   target   Z     max |S|        at theta [deg]    sigma_el [a.u.]")
    for  (name, Z, conf) in (("He",  2.0, "1s^2"),
                             ("Ne", 10.0, "1s^2 2s^2 2p^6"),
                             ("Ar", 18.0, "1s^2 2s^2 2p^6 3s^2 3p^6"),
                             ("Kr", 36.0, "1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6"))
        grid       = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 20.0)
        psSettings = ParticleScattering.Settings(ParticleScattering.Settings(),
                                                 impactEnergies = [500.0], polarThetas = thetas, polarPhis = [0.0],
                                                 printBefore = false, epsPartialWave = 1.0e-7, maxL = 120)
        wc = Atomic.Computation(Atomic.Computation(), name="Sherman function", grid=grid, nuclearModel=Nuclear.Model(Z),
                                initialConfigs  = [Configuration(conf)],
                                finalConfigs    = [Configuration(conf)],
                                processSettings = psSettings )
        wd    = perform(wc; output=true)
        event = wd["particle-scattering events:"][1]
        k     = argmax( [abs(o.sherman) for o in event.angular] )
        println("   ", rpad(name, 9), rpad(Int(Z), 6), rpad(round(abs(event.angular[k].sherman), sigdigits=4), 15),
                rpad(round(event.angular[k].theta * 180/pi, digits=1), 18),
                round(event.integrated.sigmaElastic, sigdigits=5))
    end
    #
elseif  false
    # Last successful:  unknown ...
    #
end
#
setDefaults("print summary: close", "")
