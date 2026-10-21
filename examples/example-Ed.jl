
println("Ed) Test of the Photoionization module, field-free and Debye-screened (PlasmaShift was retired) with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoIonization-Plasma.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

if  true
    # Last visit:  20-Oct-2026
    # Last successful:  15-Oct-2026
    # REWRITTEN 31-Aug-2026. It drove the retired PlasmaShift module -- PlasmaShift.PhotoSettings passed to an
    #   Atomic.Computation -- and PlasmaShift no longer exists anywhere in src/: it was folded into Plasma on
    #   09-Aug-2026, and the route through Atomic.Computation is gone with it. The current route is
    #   Plasma.LineShiftScheme inside a Plasma.Computation, which is what examples/example-Jb.jl branch b uses.
    #
    # AND IT NOW COMPUTES THE FIELD-FREE CASE TOO, WHICH IS THE POINT OF THE REWRITE. Priority item 64 (closed
    #   31-Aug-2026) was a normalisation in PhotoIonization.computeAmplitudesPropertiesPlasma that was wrong by
    #   4 pi/(alpha omega)^2 -- 174.7 at 1 keV -- and it survived because the ONLY caller in the package,
    #   example-Jb.jl branch b, ran all four of its plasma models THROUGH THE SAME PLASMA FUNCTION, including the
    #   one it called "field-free". Its checks were monotonicity, absence of discontinuities and convergence to
    #   the field-free limit; every one of those is SCALE-INVARIANT, so a factor common to all four was invisible.
    #   This branch therefore computes the field-free cross section through the OTHER route --
    #   PhotoIonization.computeLines, whose normalisation is validated against Stobbe -- and prints the two side
    #   by side. The screened value must approach the field-free one as the Debye length grows, and now that
    #   comparison crosses the two routes instead of comparing one route with itself.
    #
    # System: Kr^26+ (Ne-like), 1s^2 2s^2 2p^6 -> 2p^-1 and 2s^-1, at 3000 and 3200 eV, E1, Coulomb gauge.
    #
    # MEASURED 31-Aug-2026, total cross section summed over final levels, in barn:
    #                                   omega = 3000 eV        omega = 3200 eV
    #     field-free  computeLines        2.212e5                1.734e5
    #     screened    lambda = 1000 a_o   2.097e5  (0.948)       1.693e5  (0.976)
    #     screened    lambda =   10 a_o   1.343e5                1.400e5
    #     screened    lambda =    2 a_o   1.559e5                1.316e5
    #
    # THE CHECK THAT MATTERS is the first pair of rows: at weak screening the plasma route lands within 5 % and
    #   2.4 % of the field-free route, and those two numbers come from DIFFERENT FUNCTIONS with independently
    #   written normalisations. Before item 64 was fixed the same comparison would have shown a factor of ~175,
    #   which is precisely what no check inside example-Jb.jl could see.
    # THE NON-MONOTONICITY AT 3000 eV IS A CHANNEL OPENING, NOT AN ARTIFACT -- explained 15-Oct-2026, and the
    #   control is already in this branch's own output. The cross section falls 2.212 -> 2.097 -> 1.343 and then
    #   RISES to 1.559e5 barn, and the "No lines" column goes 2, 2, 2 -> 3 across exactly those four rows:
    #   screening destabilises the bound orbitals and so lowers the ionisation threshold, and the 2s-hole channel
    #   that is CLOSED at 3000 eV under weak screening is OPEN at lambda = 2. The extra channel adds cross
    #   section, which is why the sequence stops falling.
    #   AT 3200 eV THE COUNT IS 3 IN EVERY ROW AND THE SEQUENCE IS MONOTONIC (1.734 -> 1.693 -> 1.400 -> 1.316),
    #   which is the control: same ion, same models, same code, and the only difference is that no channel opens.
    #   So nothing here is unexplained, and the earlier "recorded, not explained" note is withdrawn.
    #   RE-RUN 15-Oct-2026: all EIGHT numbers above reproduce to every recorded digit -- 2.212102 / 1.734163 /
    #   2.096794 / 1.693182 / 1.343052 / 1.400447 / 1.558702 / 1.316345 e5 barn -- and the cross-route ratios
    #   come out 0.9479 and 0.9764 against the 0.948 and 0.976 recorded. Nothing drifted in six weeks.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    nm   = Nuclear.Model(36.)
    iC   = [Configuration("1s^2 2s^2 2p^6")]
    fC   = [Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6")]
    omegas = [3000., 3200.]

    println("\n>> (a) FIELD-FREE, through PhotoIonization.computeLines (Stobbe-validated normalisation):")
    piSettings = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=[E1], gauges=[UseCoulomb],
                                          photonEnergies=omegas, printBefore=false)
    wa = Atomic.Computation(Atomic.Computation(), name="Ed field-free", grid=grid, nuclearModel=nm,
                            initialConfigs=iC, finalConfigs=fC, processSettings=piSettings)
    wb = perform(wa; output=true)

    println("\n>> (b) SCREENED, through Plasma.LineShiftScheme -> PhotoIonization.computeLinesPlasma:")
    lineSettings = PhotoIonization.PlasmaSettings([E1], [Basics.UseCoulomb], omegas, false, LineSelection())
    for  lambda  in  [1000.0, 10.0, 2.0]
        println("\n>> Debye length = $lambda a_o")
        scheme = Plasma.LineShiftScheme(Basics.DebyeHueckelModel(lambda), iC, fC, lineSettings)
        comp   = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                    settings=Plasma.Settings())
        perform(comp, output=true)
    end
    #
end
#
setDefaults("print summary: close", "")

