
println("Dq) Apply & test the InternalConversion module (bound-electron-to-continuum ejection driven")
println("    by a nuclear gamma-ray transition of given energy/multipolarity) with initial- and")
println("    final-state multiplets declared directly via Atomic.Computation's initialConfigs/")
println("    finalConfigs/processSettings, exactly as for any other JAC process (e.g. AutoIonization,")
println("    see example-Dc.jl). Internal conversion coefficients (ICC) depend only on atomic")
println("    structure (Z, shell, continuum energy, multipolarity) -- never on nuclear input -- since")
println("    the nuclear matrix element cancels in the alpha = (IC rate)/(gamma rate) ratio; the")
println("    nuclear transition energy (settings.gammaEnergy) is the only externally-given quantity.")

if  true
    # Last successful:  29-Jul-2026 -- verified: ICC identical (to all printed digits) across
    #   maxKappa=2,3,5; ICC strictly monotonically decreasing over gammaEnergy=35,40,60,100 Hartree.
    # Branch a: internal-consistency checks -- run FIRST, before any literature comparison, mirroring
    #   how CrystalField's own textbook/group-theory checks preceded its literature comparisons. No
    #   external table is used here; what is checked is that the new code behaves the way any correct
    #   ICC implementation MUST behave, independent of any specific literature value:
    #   (i) maxKappa convergence -- summing over ever more (unphysical, angular-momentum-forbidden)
    #       high-|kappa| partial waves must NOT change the result once the allowed-channel list is
    #       exhausted (a real bug -- e.g. a broken selection rule silently admitting extra channels --
    #       would instead show ICC drifting upward without bound as maxKappa grows).
    #   (ii) monotonic, steep energy-scaling -- at fixed multipolarity and shell, ICC must decrease
    #       monotonically, and quite steeply, as the electron kinetic energy (gammaEnergy - binding
    #       energy) increases -- the qualitative, model-independent hallmark of every ICC table ever
    #       published (a sign of a wrong radial integral would be a flat, oscillating, or increasing
    #       trend instead).
    #   System: Ne(1s^2 2s^2 2p^6) -> Ne+(1s 2s^2 2p^6), M1, K-shell.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    nm   = Nuclear.Model(10.)
    initialConfigs = [Configuration("1s^2 2s^2 2p^6")];   finalConfigs = [Configuration("1s 2s^2 2p^6")]
    #
    println(">> (i) maxKappa convergence, M1, gammaEnergy = 40. Hartree:")
    iccByMaxKappa = Float64[]
    for  mk in [1, 2, 3, 5]
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M1], gammaEnergy=40.0, maxKappa=mk)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-a-mk$mk", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        icc = wb["internal conversion lines:"][1].ICC
        push!(iccByMaxKappa, icc)
        println("   maxKappa=$mk   ICC=$icc")
    end
    println(">> ICC must be identical for maxKappa=2,3,5 (channel list is exhausted at maxKappa=2 for")
    println("   M1 from a 1s_1/2 hole): $(iccByMaxKappa[2]) , $(iccByMaxKappa[3]) , $(iccByMaxKappa[4])")
    #
    println("")
    println(">> (ii) energy scaling, M1, maxKappa = 3:")
    iccByEnergy = Float64[]
    for  ge in [35.0, 40.0, 60.0, 100.0]
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M1], gammaEnergy=ge, maxKappa=3)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-a-ge$ge", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        line = wb["internal conversion lines:"][1]
        push!(iccByEnergy, line.ICC)
        println("   gammaEnergy=$ge Ha   E_k=$(line.electronEnergy) Ha   ICC=$(line.ICC)")
    end
    println(">> ICC must decrease monotonically as gammaEnergy increases: $(issorted(iccByEnergy, rev=true))")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  unknown ...
    # Branch b: 207Pb, K-shell M4 alpha_K vs. Raman et al., Phys. Rev. C 66, 044312 (2002), Table VII /
    #   Fig. 1 (gamma energy 1063.6 keV): experimental alpha_K = 0.0945 +/- 0.0022. This is a
    #   heavy, multi-shell (Z=82) SCF plus a highly relativistic (~975 keV kinetic energy) continuum
    #   orbital -- computationally substantial (many minutes), so run separately/in background before
    #   trusting or dating this branch.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    gammaEnergy_Ha = Defaults.convertUnits("energy: to atomic", 1063.6e3)
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=2.0e-3, rbox=6.0)
    nm   = Nuclear.Model(82.)
    #
    icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M4], gammaEnergy=gammaEnergy_Ha, maxKappa=6, printBefore=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dq-b-Pb207-K", grid=grid, nuclearModel=nm,
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 4f^14 5s^2 5p^6 5d^10 6s^2 6p^2")],
                               finalConfigs  =[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 4f^14 5s^2 5p^6 5d^10 6s^2 6p^2")],
                               processSettings=icSettings)
    wb = perform(comp; output=true)
    line = wb["internal conversion lines:"][1]
    println(">> 207Pb K-shell M4 alpha_K = $(line.ICC)   (exp: 0.0945 +/- 0.0022, Raman Table VII)")
    #
    setDefaults("print summary: close", "")
    #
end
