
println("Dq) Apply & test the InternalConversion module (bound-electron-to-continuum ejection driven")
println("    by a nuclear gamma-ray transition of given energy/multipolarity) with initial- and")
println("    final-state multiplets declared directly via Atomic.Computation's initialConfigs/")
println("    finalConfigs/processSettings, exactly as for any other JAC process (e.g. AutoIonization,")
println("    see example-Dc.jl). Internal conversion coefficients (ICC) depend only on atomic")
println("    structure (Z, shell, continuum energy, multipolarity) -- never on nuclear input -- since")
println("    the nuclear matrix element cancels in the alpha = (IC rate)/(gamma rate) ratio; the")
println("    nuclear transition energy (settings.gammaEnergy) is the only externally-given quantity.")

if  true
    # Last visit:  22-Oct-2026
    # Last reproduced:  22-Oct-2026   ## both claims of the dated line below: ICC = 333.60223553412624
    #   identical across maxKappa=2,3,5, and 500.3079 > 333.6022 > 99.1099 > 21.6239 over 35,40,60,100 Ha
    # Last successful:  15-Aug-2026 -- verified: ICC identical (to all printed digits) across
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
    #   acceptUnvalidated=true now only silences the module's printed note on validated scope; it is no
    #   longer a refusal (22-Oct-2026, see the module's STATUS block). Nothing here compares against a
    #   table: both checks are internal, and both remain meaningful whatever overall factor the absolute
    #   values carry -- which is why this branch's claim survived the 22-Oct repair of those values.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    nm   = Nuclear.Model(10.)
    initialConfigs = [Configuration("1s^2 2s^2 2p^6")];   finalConfigs = [Configuration("1s 2s^2 2p^6")]
    #
    println(">> (i) maxKappa convergence, M1, gammaEnergy = 40. Hartree:")
    iccByMaxKappa = Float64[]
    for  mk in [1, 2, 3, 5]
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M1], gammaEnergy=40.0, maxKappa=mk,
                                                 acceptUnvalidated=true)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-a-mk$mk", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        icc = wb[InternalConversion.Settings][1].ICC
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
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M1], gammaEnergy=ge, maxKappa=3,
                                                 acceptUnvalidated=true)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-a-ge$ge", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        line = wb[InternalConversion.Settings][1]
        push!(iccByEnergy, line.ICC)
        println("   gammaEnergy=$ge Ha   E_k=$(line.electronEnergy) Ha   ICC=$(line.ICC)")
    end
    println(">> ICC must decrease monotonically as gammaEnergy increases: $(issorted(iccByEnergy, rev=true))")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:  22-Oct-2026   ## alpha_K = 9.602494e-02 against the measured 0.0945(22), ratio 1.016
    # Last successful:  unknown ...
    # Branch b: 207Pb, K-shell M4 alpha_K vs. Raman et al., Phys. Rev. C 66, 044312 (2002), Table VII /
    #   Fig. 1 (gamma energy 1063.6 keV): experimental alpha_K = 0.0945 +/- 0.0022. This is a
    #   heavy, multi-shell (Z=82) SCF plus a highly relativistic (~975 keV kinetic energy) continuum
    #   orbital -- computationally substantial, so run separately/in background before trusting or
    #   dating this branch. Two Z=82 SCF passes dominate the cost at roughly 800 s each; WITHOUT the
    #   lineSelection below the Cartesian product of the 6s^2 6p^2 initial levels with the 1s-hole
    #   final ones took 5268 s on 22-Oct-2026 to compute lines nobody reads.
    #   This is the RELATIVISTIC anchor and branch c is its counterpart: see branch c for why the two
    #   are a pair and not a duplication.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    gammaEnergy_Ha = Defaults.convertUnits("energy: to atomic", 1063.6e3)
    # hp = 5.0e-4, not the 1.5e-3 this branch carried until 22-Oct-2026. Continuum.gridConsistency, which
    # set that value, is an A PRIORI check and passes 1.5e-3 at 15.65 points per oscillation -- and the
    # orbital it returns there is nevertheless unusable: Continuum.checkContinuumAmplitude, a POSTERIOR
    # check on the normalized orbital, now REFUSES it outright. Both figures are correct; they measure
    # different things, and only the second one looks at what came back.
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=5.0e-4, rbox=6.0)
    nm   = Nuclear.Model(82.)
    #
    icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M4], gammaEnergy=gammaEnergy_Ha, maxKappa=6, printBefore=true,
                                             lineSelection=LineSelection(true, indexPairs=[(1,1)]), acceptUnvalidated=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dq-b-Pb207-K", grid=grid, nuclearModel=nm,
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 4f^14 5s^2 5p^6 5d^10 6s^2 6p^2")],
                               finalConfigs  =[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 4f^14 5s^2 5p^6 5d^10 6s^2 6p^2")],
                               processSettings=icSettings)
    wb = perform(comp; output=true)
    line = wb[InternalConversion.Settings][1]
    println(">> 207Pb K-shell M4 alpha_K = $(line.ICC)   (exp: 0.0945 +/- 0.0022, Raman Table VII)")
    #
    setDefaults("print summary: close", "")
    #

elseif  false
    # Last visit:  22-Oct-2026   ## alpha_K = 1.530727e+03, against the no-vacancy 1544(1) of Table IV
    # Last successful:  unknown ...
    # Branch c: 119Sn, K-shell M4 alpha_K vs. Nica et al., Phys. Rev. C 89, 014303 (2014), Table IV
    #   (gamma energy 65.660 keV). That table is why this branch exists: it quotes THREE calculations of
    #   the same coefficient which differ ONLY in how the K-shell vacancy is treated --
    #       no vacancy                 1544(1)
    #       vacancy, frozen orbitals   1618(1)
    #       vacancy, SCF of the ion    1603(1)
    #   against a measurement of 1621(25). A number here can therefore be placed against a MODEL and not
    #   merely against an experiment, which is what makes this anchor diagnostic rather than only a check.
    #   JAC gives 1530.7, which sits 0.9 % from the no-vacancy value and 4.5 % below the SCF-of-ion one.
    #   That is an OPEN QUESTION and the most interesting thing in this file: the continuum orbital here
    #   IS generated in the potential of the final (1s-hole) level, so the vacancy is present, yet its
    #   effect comes out under a per cent where the table puts it at +3.8 %. Isolating it needs one run
    #   with the ion potential against one with the neutral potential at the SAME electron energy, which
    #   the perform route cannot express today -- the choice sits inside Continuum. That is a src/ gap.
    #   WHY THIS BRANCH AND NOT ONLY BRANCH b: same multipole, same shell, OPPOSITE kernel regime. Here
    #   the ejected electron has 36 keV and k_photon*r = 0.35 at the 1s radius, so the Hankel kernel h_4
    #   is still in its r^(-5) limit exactly where the bound orbital sits; in branch b the electron has
    #   975 keV and k_photon*r = 3.5, so the kernel is oscillatory there instead. The inner round-off
    #   floor that InternalConversion.icInnermostReliablePoint excludes was located by precisely this
    #   contrast: a cut at the NUCLEAR radius lies outside the floor for 207Pb and inside it for 119Sn,
    #   so it looked stable on the one case while moving alpha_K by a factor 1.94 on the other.
    #   lineSelection restricts the Cartesian product: Sn's 5s^2 5p^2 ground configuration and the 1s-hole
    #   ion both carry several levels, so the unrestricted run generated 770 continuum orbitals on
    #   22-Oct-2026 where the one physically meaningful pair needs 2. It buys far less WALL time than that
    #   ratio suggests -- both runs came in near 14 minutes, because two Z=50 SCF passes dominate and the
    #   continuum work does not. Neither figure is a clean timing: three Julia processes shared the machine.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    gammaEnergy_Ha = Defaults.convertUnits("energy: to atomic", 65.660e3)
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=3.0e-3, rbox=20.0)
    nm   = Nuclear.Model(50.)
    #
    icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M4], gammaEnergy=gammaEnergy_Ha,
                                             maxKappa=6, lineSelection=LineSelection(true, indexPairs=[(1,1)]),
                                             acceptUnvalidated=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dq-c-Sn119-K", grid=grid, nuclearModel=nm,
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^2")],
                               finalConfigs  =[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^2")],
                               processSettings=icSettings)
    wb   = perform(comp; output=true)
    line = wb[InternalConversion.Settings][1]
    println(">> 119Sn K-shell M4 alpha_K = $(line.ICC)")
    println(">>   against Nica et al. Table IV:  no vacancy 1544(1), SCF of ion 1603(1), frozen orbitals 1618(1),")
    println(">>   measured 1621(25).   Ratio to the no-vacancy value: $(line.ICC/1544.0)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:  22-Oct-2026   ## E1 3.2658e+03, M1 3.3360e+02, E2 4.3358e+05, M2 2.8184e+05,
    #   E3 3.2065e+07, M3 4.0524e+07 -- structure only, no literature claim
    # Last successful:  unknown ...
    # Branch d: the MULTIPOLE dependence at fixed Z, shell and transition energy -- Ne K-shell, 40 Ha --
    #   and the one trap in this module's interface.
    #   (i) alpha falls steeply with L and, at equal L, the electric and magnetic coefficients differ in
    #       structure and not merely in size: the magnetic radial factor is the single integral of Eq. (2),
    #       while the electric one is the two-term combination of Eq. (3), carrying h_(L-1) beside h_L.
    #       This branch prints them side by side so that the two forms can be seen to behave differently.
    #   (ii) ONE MULTIPOLE PER RUN. The Settings field is an array and its docstring once advertised
    #       "mixed multipolarity supported, e.g. [M1,E2]", but the module sums |amplitude|^2 over every
    #       channel of every multipole it is given. For a genuinely mixed transition the coefficient is
    #       (alpha(M1) + delta^2 alpha(E2)) / (1 + delta^2) with the NUCLEAR mixing ratio delta, which the
    #       module neither receives nor could know. An array therefore yields a sum with no physical
    #       meaning. The multipoles are computed separately here and combined by whoever knows delta.
    #   NO LITERATURE CLAIM: the electric multipoles are not validated against any table (both anchors in
    #   this file are magnetic), so this branch documents structure and interface, not accuracy.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    nm   = Nuclear.Model(10.)
    initialConfigs = [Configuration("1s^2 2s^2 2p^6")];   finalConfigs = [Configuration("1s 2s^2 2p^6")]
    #
    println(">> Ne K-shell, gammaEnergy = 40 Hartree, ONE multipole per run:")
    for  mp in [E1, M1, E2, M2, E3, M3]
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[mp], gammaEnergy=40.0,
                                                 maxKappa=6, acceptUnvalidated=true)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-d-$mp", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        println("   $mp  (L=$(mp.L), $(mp.electric ? "electric, Eq. 3, two terms" : "magnetic, Eq. 2, one term"))" *
                "   alpha = $(wb[InternalConversion.Settings][1].ICC)")
    end
    println(">> These must NOT be added together to describe a mixed transition; see (ii) above.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:  23-Oct-2026   ## alpha_K(119Sn M4) = 1530.73, 1529.73, 1519.26 for hp=3.0e-3, 1.5e-3,
    #   1.0e-3 at rbox=20; 1558.37, 1530.73, 1524.03 for rbox=15, 20, 30 at hp=3.0e-3; guard test raised
    # Last successful:  unknown ...
    # Branch e: the NUMERICAL approximation -- is alpha_K a property of the physics or of the grid?
    #   This is the branch that demonstrates the 22-Oct-2026 repair, and it does so on the perform route
    #   rather than by calling the module's internals.
    #   (i) alpha_K under refinement of hp and of rbox. Refining the grid also MOVES the inner match point
    #       that InternalConversion.icInnermostReliablePoint selects, since that point is read off the
    #       orbital rather than fixed at a radius -- so a stable alpha_K across these runs tests the match
    #       point and the quadrature together, which is the honest test of the repair.
    #   (ii) a grid too coarse for the ejected electron is REFUSED rather than quietly returning a number.
    #       The last run below is a GUARD TEST: raising is the PASS, and completing would be the defect.
    #       Which guard fires is itself informative, and MEASURED here rather than assumed: at hp=5.0e-2
    #       the A PRIORI Continuum.gridConsistency refuses first, at 2.43 points per oscillation against a
    #       required 15. The POSTERIOR Continuum.checkContinuumAmplitude is therefore NOT exercised by
    #       this branch and cannot be -- it exists precisely for the grids that pass the a-priori test and
    #       still carry an unusable orbital. Judge the branch by whether the error names a grid or
    #       amplitude guard; an error from anywhere else would be an earlier fault masking the guard.
    #   Sn rather than Ne, because the grid question is sharpest where the ejected electron is fast, and
    #   lineSelection again, for the reason given in branch c. THIS IS THE EXPENSIVE BRANCH OF THE FILE:
    #   six independent runs, each with its own pair of Z=50 SCF passes on a grid that the refinement
    #   makes progressively finer. Measured 23-Oct-2026 at 3 h 10 min wall, sharing the machine with two
    #   other Julia processes -- so budget hours and run it in the background, not interactively.
    setDefaults("print summary: open", "zzz-InternalConversion.sum")
    #
    gammaEnergy_Ha = Defaults.convertUnits("energy: to atomic", 65.660e3)
    nm   = Nuclear.Model(50.)
    initialConfigs = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^2")]
    finalConfigs   = [Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^2")]
    #
    println(">> (i) 119Sn K-shell M4 alpha_K under grid refinement:")
    for  (hp, rbox) in [(3.0e-3, 20.0), (1.5e-3, 20.0), (1.0e-3, 20.0), (3.0e-3, 15.0), (3.0e-3, 30.0)]
        grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=hp, rbox=rbox)
        icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M4], gammaEnergy=gammaEnergy_Ha,
                                                 maxKappa=6, lineSelection=LineSelection(true, indexPairs=[(1,1)]),
                                                 acceptUnvalidated=true)
        comp = Atomic.Computation(Atomic.Computation(), name="Dq-e-hp$hp-rbox$rbox", grid=grid, nuclearModel=nm,
                                   initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
        wb = perform(comp; output=true)
        println("   hp=$hp  rbox=$rbox   alpha_K=$(wb[InternalConversion.Settings][1].ICC)")
    end
    #
    println("")
    println(">> (ii) GUARD TEST -- a grid too coarse for a 36 keV electron must RAISE, not return a number:")
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=5.0e-2, rbox=20.0)
    icSettings = InternalConversion.Settings(InternalConversion.Settings(); multipoles=[M4], gammaEnergy=gammaEnergy_Ha,
                                             maxKappa=6, lineSelection=LineSelection(true, indexPairs=[(1,1)]),
                                             acceptUnvalidated=true)
    comp = Atomic.Computation(Atomic.Computation(), name="Dq-e-guard", grid=grid, nuclearModel=nm,
                               initialConfigs=initialConfigs, finalConfigs=finalConfigs, processSettings=icSettings)
    try
        wb = perform(comp; output=true)
        println(">> GUARD TEST FAILED: it completed and returned $(wb[InternalConversion.Settings][1].ICC)")
    catch err
        println(">> GUARD TEST PASSED (raising is the expected outcome):")
        println("   ", sprint(showerror, err)[1:min(end,600)])
    end
    #
    setDefaults("print summary: close", "")
    #
end
