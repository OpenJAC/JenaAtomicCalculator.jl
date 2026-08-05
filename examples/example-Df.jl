
println("Df) Apply & test the DielectronicRecombination module with ASF from internally generated initial-, intermediate and final-state multiplets.")

setDefaults("print summary: open", "zzz-DielectronicRecombination.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")
setDefaults("unit: strength", "cm^2 eV")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, asymptotic Coulomb")
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")

## NOTE (04-Aug-2026): this file was rewritten from scratch. Its seven previous branches were ALL dead code:
## every one of them called `Dielectronic.Settings(...)`, but there is no module or alias `Dielectronic`
## anywhere in src/ -- the module is `DielectronicRecombination`. They also passed 8 or 9 positional arguments
## (two branches used mutually inconsistent signatures) to a Settings struct that has THIRTEEN fields, used the
## removed `process = Dierec()` keyword, and relied on a `JAC.` binding that is not exported. Two branches were
## simultaneously `true`, so the second could never run. Five carried "Last successful: unknown", two claimed
## "13May2024". The physics INTENT is preserved below -- the Tu et al. and Xu et al. comparisons are kept as
## targets -- but none of the old code was salvageable.
##
##   a)  He-like C4+ + e- --> Li-like C3+, K-LL          -- the small anchor; every check is internal
##   b)  the same system, plus the rate coefficient alpha(T)
##   c)  branch a)'s system through the RE-STRUCTURED route -- equivalence + the DR satellite spectrum
##   d)  branch b)'s Xe53+ case through the re-structured route -- the Harman regression
##
## THE TWO STRATEGIES OF THIS MODULE, for orientation:
##   * FINE-STRUCTURE resolved (module-DielectronicRecombination-inc-FS-resolved.jl) -- the working path. Two
##     sub-modes: full pathway computation (default), and `calcOnlyPassages`, which collects each resonance's
##     contribution directly and is meant for high-n Rydberg intermediates and semi-empirical corrections.
##   * HYPERFINE resolved (module-DielectronicRecombination-inc-HF-resolved.jl) -- reachable ONLY as
##     `calcHyperfineResolved && calcOnlyPassages`, i.e. a sub-mode of passages, not a peer of the FS path.
##     Physically it gives the INITIAL (recombining) ion the full nuclear moments and sets mu = Q = 0 for the
##     intermediate and final levels, so that F is a good quantum number for the capture step.
##     As of 04-Aug-2026 this path CANNOT RUN: computeHyperfineAmplitudes, displayHyperfineResults and
##     displayHyperfineRateCoefficients are each called but never defined (the amplitude function is defined
##     under the name `computecomputeHyperfineAmplitudes`), displayHyperfinePassages is handed an undefined
##     variable `passages`, and the non-distributed path passes a Multiplet where an Hfs.HfMultiplet is required.
##
## Note that Settings() defaults to gauges = UseGauge[] -- EMPTY. Set the gauges explicitly or no radiative
## channel is computed at all.

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


if  false
    # Last visit:  04-Aug-2026
    # --- Branch a: He-like C4+ (Z=6) + e- --> Li-like C3+, K-LL dielectronic recombination.
    #
    # WHY THIS FIRST. It is the smallest system in which DR is a genuine three-step process, and small enough
    # that every quantity can be checked by hand from the module's own printout, with no literature at all:
    #
    #   1. ENERGY CONSISTENCY. electronEnergy must equal E(intermediate) - E(initial), and photonEnergy must
    #      equal E(intermediate) - E(final). Both are readable off the level tables that precede the results.
    #
    #   2. THE RESONANCE-STRENGTH FORMULA. The Resonance table prints resonanceEnergy, captureRate, augerRate
    #      and photonRate SEPARATELY from resonanceStrength, so the strength can be recomputed by hand from
    #
    #           S = (h^2 / (2 m_e E_res)) * (g_d / (2 g_i)) * A_a * A_r / (A_a + A_r)
    #
    #      and compared. This is the same style of check that confirmed the Maxwellian average in example-Dl.jl.
    #
    #   3. THE SUM RULE. Summing the per-pathway reducedStrength over all final levels of one intermediate level
    #      must reproduce that intermediate level's resonanceStrength.
    #
    #   4. GAUGE AGREEMENT on everything radiative (photonRate, resonanceStrength). The capture side is
    #      gauge-independent by construction, so a gauge spread here localises to the radiative step.
    #
    # The K-LL resonances are the 1s2l2l' doubly-excited levels of Li-like carbon: a free electron is captured
    # into a 2l orbital while a 1s electron is simultaneously promoted to 2l', and the resulting autoionizing
    # level then stabilises by emitting a photon down to 1s^2 2l. Literature target, not yet compared here:
    # Xu et al., PRA 2016.
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles  = [E1],
                                                    gauges      = [UseCoulomb, UseBabushkin],
                                                    printBefore = false )
    wa = Atomic.Computation(Atomic.Computation(), name="Df-a: K-LL DR into He-like C4+", grid=grid,
                            nuclearModel        = Nuclear.Model(6.),
                            initialConfigs      = [Configuration("1s^2")],
                            intermediateConfigs = [Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            processSettings     = drSettings )
    perform(wa)
    #
elseif  false
    # Last visit:  04-Aug-2026
    # --- Branch b: K-LL DR of H-like Xe53+ (Z=54). THE LITERATURE ANCHOR of this file.
    #
    # This is the one case for which independent reference data sit inside this repository:
    #   apps/b24-apps-schippers-dr-breit-xenon/b24.01-harmann-Xe53DRtheo4.dat  -- Harman's theoretical resonance
    #   energies, strengths and widths, and apps-dr-xenon.report-01..03, which hold earlier JAC output for the
    #   same system. (The app's own script b24.01-Xe53KLLDR.jl is itself dead under the current API: it calls
    #   `Dielectronic.Settings` with 10 positional arguments.)
    #
    # It is also STRUCTURALLY SIMPLER than branch a): the recombining ion is H-like, so ONE electron is present
    # before capture and the resonance is a doubly-excited He-like 2l2l' state with an EMPTY K shell. The K-LL
    # resonances therefore sit near 20.6-21.6 keV, which is what makes the branch expensive rather than difficult.
    #
    # DELIBERATE DIFFERENCE FROM THE APP SCRIPT. Its finalConfigs were {1s^2, 1s2s, 2s2p} -- which omits 1s2p and
    # instead lists 2s2p, an autoionizing configuration identical to one of the intermediates. The stabilizing
    # transition is 2p --> 1s, so 2p^2 --> 1s2p is the main radiative channel of the 2p^2 resonances and must be
    # present; 2s^2 has no E1 channel at all (2s --> 1s is forbidden), and 1s^2 cannot be reached from 2l2l' by a
    # single E1 at all. The final set used here is therefore {1s^2, 1s2s, 1s2p}: 1s^2 is kept only so that its
    # absence of channels is visible rather than assumed.
    #
    # WHAT THIS BRANCH IS FOR. Branch a) verified ratios, a sum rule and the width relation -- all of which are
    # blind to an overall scale factor. The module's reducedStrength carries (2J_m+1)/(2J_i+1) with a factor 2
    # explicitly commented out and marked "factor 2 is not really clear" in
    # module-DielectronicRecombination-inc-FS-resolved.jl. Only an ABSOLUTE comparison can decide that, and this
    # is it.
    grid54     = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    # QED matters at Z = 54: it shifts the K-shell binding by tens of eV and hence the resonance energies. The
    # standard self-consistent field is used; no special SCF field is needed or wanted here.
    asfSet54   = AsfSettings(AsfSettings(), qedModel=QedPetersburg(), eeInteraction=CoulombBreit(1.0))
    #
    # THE BREIT INTERACTION IS NOT OPTIONAL HERE. Run first with the default augerOperator = CoulombInteraction()
    # and the RADIATIVE rates come out within 0.1-1.6% of the 2024 reference, while the AUGER rates fall short by
    # 19%, 39%, 73% and 93% for m = 1, 5, 2 and 6 respectively -- the shortfall growing as the Auger rate itself
    # gets smaller. That is the signature of Breit opening capture channels that the Coulomb interaction alone
    # barely supports, and at Z = 54 it dominates the weak resonances. The application directory this benchmark
    # comes from is named "schippers-dr-breit-xenon" precisely because that is the effect it was built to study.
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles    = [E1],
                                                    gauges        = [UseCoulomb, UseBabushkin],
                                                    augerOperator = CoulombBreit(1.0),
                                                    printBefore   = false )
    wb = Atomic.Computation(Atomic.Computation(), name="Df-b: K-LL DR of H-like Xe53+", grid=grid54,
                            nuclearModel        = Nuclear.Model(54., "Fermi"),
                            # NOTE the explicit ZERO occupations. DielectronicRecombination.checkConsistentMultiplets
                            # requires the initial, intermediate and final states to be built over the SAME ordered
                            # subshell list (the initial one may be shorter), a restriction coming from the angular
                            # coefficients. Writing Configuration("2s^2") would drop 1s from the intermediate subshell
                            # list and abort the run; "1s^0 2s^2 2p^0" keeps it. This is why the app script uses this
                            # notation throughout -- it is load-bearing, not decoration.
                            initialConfigs      = [Configuration("1s^1 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^0 2s^2 2p^0"), Configuration("1s^0 2s^1 2p^1"),
                                                   Configuration("1s^0 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^0 2p^0"), Configuration("1s^1 2s^1 2p^0"),
                                                   Configuration("1s^1 2s^0 2p^1")],
                            initialAsfSettings  = asfSet54, intermediateAsfSettings = asfSet54, finalAsfSettings = asfSet54,
                            processSettings     = drSettings )
    perform(wb)
    #
    # ============================== REPORT (04-Aug-2026) ==============================
    #
    # Wall time 1:50, peak RSS 2.5 GB. Ten K-LL resonances, 20.6 - 21.6 keV.
    #
    # COMPARISON WITH HARMAN (apps/b24-apps-schippers-dr-breit-xenon/b24.01-harmann-Xe53DRtheo4.dat):
    #
    #   Harman E [eV]   Harman S      this run (m)   S [cm^2 eV], Coulomb   ratio
    #     20632.4       2.28e-20        m1, 20625      2.2381e-20           0.982
    #     20631.7       7.6e-21         m2, 20628      7.7012e-21           1.013
    #     20655.6       1.75e-20        m3, 20653      1.7102e-20           0.977
    #     21051         7.36e-21        m5, 21048      7.3938e-21           1.005
    #     21095.8       1.05e-21        m6, 21093      1.0826e-21           1.031
    #     21112.3       6.44e-20        m7, 21113      6.3310e-20           0.983
    #     20726.8       1.16e-21        m4, 20721      8.3063e-22           0.716   <-- outlier
    #
    # SUCCESSFUL. Six of seven resonances agree with independent theory to 0.5 - 3.1%.
    #
    #  1. THE FACTOR-2 QUESTION IS SETTLED. reducedStrength carries (2J_m+1)/(2J_i+1) with a factor 2 commented
    #     out and marked "factor 2 is not really clear" in module-DielectronicRecombination-inc-FS-resolved.jl.
    #     Agreement at the 1-3% level across six independent resonances cannot coexist with a missing factor of
    #     two. The factor must stay commented out. Do not reopen this without an equally strong comparison.
    #
    #  2. BREIT IS ESSENTIAL AT Z = 54, and the evidence is clean because this branch was first run WITHOUT it.
    #     With augerOperator = CoulombInteraction() the radiative rates were already within 0.07-1.6% of the
    #     reference while the Auger rates fell short by 19% (m1), 39% (m5), 73% (m2) and 93% (m6) -- the deficit
    #     growing as the Auger rate itself gets smaller. Breit opens capture channels that Coulomb alone barely
    #     supports. Switching to CoulombBreit(1.0) moved the strengths from 18-93% off to 0.5-3.1%.
    #
    #  3. THE OUTLIER, m4 at 20721 eV, 28% low. m4 is the 2p^2 3P_0 level: J = 0 and even parity, so
    #     autoionization into 1s (J = 1/2) admits a SINGLE channel, eps s_1/2, and as a pure triplet its Auger
    #     amplitude is small (1.42e13 against 4.49e14 for m1). Its width is radiative-dominated, so S is very
    #     nearly proportional to A_a and a 28% deficit in the strength is a 28% deficit in that one
    #     near-cancelling single-channel amplitude. Adding 1s2p to the final configurations -- absent from the
    #     2024 app run -- improved it from 32% to 28%, so the missing configuration was NOT the cause. This is a
    #     plausible mechanism consistent with everything else seen here, NOT a demonstration; it remains open.
    #
    # WHY "Last visit" AND NOT "Last successful" (Rule 7). The agreement is real, but the reference is a single
    # theoretical dataset, not experiment, and the one outlier is unexplained. A "Last successful" date should
    # wait until the m4 amplitude is understood.
    #
elseif  false
    # Last visit:  05-Aug-2026
    # --- Branch c: branch a)'s system through the RE-STRUCTURED route (CaptureLine + PhotonLine).
    #
    # Same He-like C4+ K-LL case as branch a), but with useFsNewRoute = true, so that
    # module-DielectronicRecombination-inc-FSnew-resolved.jl is exercised instead of the Pathway/Passage/Resonance
    # code. The point is EQUIVALENCE: the "Total Auger rates, radiative rates and resonance strengths" table must
    # come out bit-identical to branch a). Any difference is a defect in the new route, since the physics, the
    # grid and the settings are otherwise the same.
    #
    # calcPhotonSpectrum = true additionally prints the DR satellite spectrum, which the old route cannot produce:
    # it computes the per-final-level photon energies and then discards them. Two identities are checkable there
    # without any reference data, and they are the unit test for setTotalRates():
    #     sum_f A_r(m,f)      must equal the totalPhotonRate column of every capture line sharing that m
    #     sum_f intensity     must equal that capture line's resonance strength
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    printBefore        = false,
                                                    useFsNewRoute      = true,
                                                    calcPhotonSpectrum = true )
    wc = Atomic.Computation(Atomic.Computation(), name="Df-c: K-LL DR into He-like C4+, FSnew route", grid=grid,
                            nuclearModel        = Nuclear.Model(6.),
                            initialConfigs      = [Configuration("1s^2 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^1 2s^2 2p^0"), Configuration("1s^1 2s^1 2p^1"),
                                                   Configuration("1s^1 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^1 2p^0"), Configuration("1s^2 2s^0 2p^1")],
                            processSettings     = drSettings )
    perform(wc)
    #
elseif  true
    # Last visit:  05-Aug-2026
    # --- Branch d: branch b)'s Xe53+ literature case through the RE-STRUCTURED route.
    #
    # The regression that matters most: the six resonances matched against Harman in branch b) must stay within
    # 0.5-3.1%, and the m4 outlier must stay at 28% -- neither improved nor worsened by the rewrite. A rewrite
    # that changed the agreement would be changing physics, which it must not.
    grid54     = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    asfSet54   = AsfSettings(AsfSettings(), qedModel=QedPetersburg(), eeInteraction=CoulombBreit(1.0))
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                                                    multipoles         = [E1],
                                                    gauges             = [UseCoulomb, UseBabushkin],
                                                    augerOperator      = CoulombBreit(1.0),
                                                    printBefore        = false,
                                                    useFsNewRoute      = true,
                                                    calcPhotonSpectrum = true )
    wd = Atomic.Computation(Atomic.Computation(), name="Df-d: K-LL DR of H-like Xe53+, FSnew route", grid=grid54,
                            nuclearModel        = Nuclear.Model(54., "Fermi"),
                            initialConfigs      = [Configuration("1s^1 2s^0 2p^0")],
                            intermediateConfigs = [Configuration("1s^0 2s^2 2p^0"), Configuration("1s^0 2s^1 2p^1"),
                                                   Configuration("1s^0 2s^0 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^0 2p^0"), Configuration("1s^1 2s^1 2p^0"),
                                                   Configuration("1s^1 2s^0 2p^1")],
                            initialAsfSettings  = asfSet54, intermediateAsfSettings = asfSet54, finalAsfSettings = asfSet54,
                            processSettings     = drSettings )
    perform(wd)
    #
end
#
setDefaults("print summary: close", "")
