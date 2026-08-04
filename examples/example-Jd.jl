#
println("Jd) Apply & test the Plasma.SatelliteDiagnosticScheme (dielectronic-satellite-to-parent-line")
println("    intensity-ratio Te diagnostic).")

# Original design sketch (18-Jul-2026, condensed): proposed a dielectronic-satellite-to-parent-line
# intensity-ratio Te diagnostic for an open-shell recombining ion (Cl-like 3s^2 3p^5, 2P_3/2+2P_1/2),
# built around fields recombiningIonRefConfig/captureFromShells/captureToShells/captureExcitationDegree/
# referenceInitialConfigs/referenceFinalConfigs/lineSelection/plasmaModel. None of these names survived
# into the real implementation below: Plasma.Computation already supplies refConfigs/nuclearModel/grid/
# asfSettings for every scheme, so they are never duplicated on the scheme itself; the capture/reference
# fields collapsed into fromShells/toShells/intoShells/decayShells, reusing Basics.ForDielectronicRecombination
# directly; and DielectronicRecombination.Settings/ImpactExcitation.Settings are embedded as-is
# (drSettings/ieSettings) rather than re-flattened. See branch b's header for the full account of what
# changed and why, including the resolution of the sketch's open "electron-impact excitation" gate.


if  true
    # Last visit:  02-Aug-2026
    # Last successful:  unknown ...
    # Branch a: the Gabriel (1972) case -- He-like Fe (Fe XXV), the classic, most extensively tabulated
    #   dielectronic-satellite Te diagnostic in the literature, and the natural PRIMARY validation case
    #   for this scheme: simpler than branch b's F-like example (He-like ground = 1 level, no fine
    #   structure at all; "1s2p" excited manifold = only 4 levels, vs. F-like's 28), and -- unlike
    #   branch b -- the "branching fraction = 1 (pure radiative, no autoionization)" assumption the
    #   driver makes is not just a reasonable default here, it is STRUCTURALLY GUARANTEED: autoionizing
    #   a singly-excited configuration needs a THIRD electron to eject while the promoted one falls back
    #   into the hole; a 2-electron system never has one, regardless of which shell is promoted. This is
    #   exactly why the w/x/y/z lines are ordinary, sharp, bound-bound transitions in every textbook
    #   treatment, not an assumption Gabriel's own analysis has to make.
    #   fromShells=[1s], toShells=[2p]: the K-shell promotion behind the w/x/y lines (from 1s2p's
    #   4 levels; z, from 1s2s ^3S_1, is a DIFFERENT core promotion not attempted here).
    #   intoShells=[2s]: the classic KLL (n'=2) Li-like satellite capture, "1s2s2p", giving the
    #   satellites closest to and most commonly discussed alongside w/x/y (Gabriel's j/k-type lines);
    #   not the same n'=3 regime as branch b.
    #   decayShells=[1s]: the promoted 2p electron decays back into the 1s hole, exactly mirroring
    #   fromShells, as in branch b.
    #   lineSelection targets level 4 of "1s2p" (J=1, the highest-energy of the 4 sublevels): verified
    #   directly before running the full scheme -- SCF/CI alone (no DR/IE yet) gives level ordering
    #   J=0,1,2,1 by INCREASING energy, matching the well-known high-Z He-like pattern
    #   ^3P_0 < ^3P_1 < ^3P_2 < ^1P_1; the w-line photon energy from level 4 comes out as 246.567 Ha =
    #   6710.5 eV, within 0.15% of the literature value 6700 eV (Bely-Dubau, Gabriel & Volonte; via
    #   web search, 01-Aug-2026) -- a strong check before the DR/IE machinery is even engaged. Level 1
    #   (J=0) is E1-forbidden to the J=0 ground (a J=0-->J=0 electric-dipole transition is forbidden
    #   outright); level 3 (J=2) is E1-forbidden by DeltaJ<=1. Only levels 2 (^3P_1, the y line) and 4
    #   (^1P_1, the w line) are E1-allowed; level 4 (w) is used here as the classic, strongest case.
    #
    #   STATUS (02-Aug-2026), PAUSED HERE FOR THIS SESSION -- SatelliteDiagnosticScheme itself is NOT
    #   being pushed on further for now:
    #     - Satellite/DR side (steps 1-3 of the driver) IS validated: alpha_DR comes out real, finite
    #       and physically sensible (~3e-14 cm^3/s, the right order of magnitude for a DR rate
    #       coefficient) once two real bugs were fixed (see branch b's header for the gauges=[] default,
    #       and the driver's satelliteAlpha loop for a per-resonance NaN guard needed specifically here:
    #       one resonance had BOTH augerRate==0 and photonRate==0 -- a genuine 0/0, unlike branch b's
    #       NaN case which had a real nonzero radiative rate -- and would otherwise have silently
    #       poisoned the whole sum via NaN-propagation through addition).
    #     - Parent/IE side (step 5) is BLOCKED: alpha_exc comes out catastrophically wrong (~1e26 to
    #       1e30 cm^3/s -- no real atomic rate coefficient exceeds ~1e-7 cm^3/s, so this is off by
    #       ~33-37 orders of magnitude), even after fixing the temperature range to properly straddle
    #       the ~7.79e7 K threshold (an initial attempt at 5e5-1e7 K, copied from branch b's much lower
    #       ~33 eV/~3.8e5 K threshold without rescaling, was ALSO wrong, but rescaling alone did not fix
    #       the blow-up -- it is a genuine ImpactExcitation numerical bug in the high-threshold/high-Z
    #       rate-coefficient integration, not a temperature-choice mistake).
    #   Per explicit instruction: module-ImpactExcitation.jl needs to be investigated and hardened
    #   INDEPENDENTLY, in its own dedicated session (real bugs already on record there: the stale
    #   RateSettings reference in example-Dl.jl, the printBefore=true displayLines() crash, and now this
    #   high-threshold rate-coefficient blow-up) -- not as a side effect of exercising this scheme. Once
    #   that is done, re-run this branch and, if alpha_exc comes out sane, this is likely close to a
    #   genuine, literature-comparable R(Te) curve (the satellite side and the w-line energy check both
    #   already agree well with the literature).
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(26.)                               # Fe XXV, He-like

    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                    calcRateAlpha = true,
                    temperatures  = [2e7, 5e7, 1e8, 2e8, 5e8] )               # [K] -- straddles the w-line threshold
                    # (~6710 eV = 7.79e7 K); an earlier attempt at 5e5-1e7 K (matching F-like Ne+'s much
                    # lower ~33 eV threshold, ~3.8e5 K) was far below threshold and gave both meaningless
                    # near-zero DR/IE rates at the low end and a real numerical blow-up in the IE rate
                    # coefficient at the high end (parentAlpha ~1e16-1e23 cm^3/s -- not a physical value
                    # at any temperature) -- not yet understood, and side-stepped here rather than chased,
                    # by using temperatures that are actually appropriate for this threshold.

    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                    calcRateCoefficient  = true,
                    temperatures         = Float64[],                         # left empty -- driver reconciles from drSettings
                    lineSelection        = LineSelection(true, indexPairs=[(1,4)]),
                    maxKappa             = 15,                                # convergence ceiling, not a fixed loop count
                    numElectronEnergies  = 3,
                    maxEnergyMultiplier  = 4.0 )

    scheme = Plasma.SatelliteDiagnosticScheme(Plasma.SatelliteDiagnosticScheme();
                 fromShells  = [Shell("1s")],
                 toShells    = [Shell("2p")],
                 intoShells  = [Shell("2s")],                                 # KLL (n'=2) Li-like satellite
                 decayShells = [Shell("1s")],
                 drSettings  = drSettings,
                 ieSettings  = ieSettings )

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                     refConfigs=[Configuration("1s^2")], asfSettings=AsfSettings(),
                                     settings=Plasma.Settings() )

    wb = perform(computation, output=true)
    #
elseif  true
    # Last visit:  02-Aug-2026
    # Last successful:  unknown ...
    # Branch b: F-like Ne+, the more complex, OPEN-SHELL recombining-ion case (2P_3/2+2P_1/2 ground,
    #   28-level "2p^4 3d" excited manifold) worked through and debugged FIRST, before branch a's
    #   simpler Gabriel case was added -- kept as the "real bugs found and fixed here" record, and as a
    #   genuine second validation once branch a is unblocked. This branch's own R(Te) came out real,
    #   finite and PLAUSIBLE-looking (0.409 at 50000 K down to 0.012 at 1.0e6 K, Coulomb/Babushkin
    #   agreeing to ~11%, decreasing monotonically with Te as physically expected) -- but branch a later
    #   found module-ImpactExcitation.jl to have a severe rate-coefficient numerical bug (off by ~33-37
    #   orders of magnitude) in a different (high-threshold/high-Z) regime, which this branch's own
    #   ~33 eV/Z=10 case did not trigger but was never independently cross-checked against either.
    #   Given that, "plausible-looking" is deliberately NOT being treated as "validated" here -- no
    #   "Last successful" date until module-ImpactExcitation.jl is hardened independently and this
    #   branch is re-run to confirm the same numbers still come out.
    #   What changed from the original design sketch, and why:
    #     - recombiningIonRefConfig/nuclearModel/grid/asfSettings dropped from the scheme entirely --
    #       Plasma.Computation already carries refConfigs/nuclearModel/grid/asfSettings for every scheme
    #       (confirmed by reading its struct directly), so LineShiftScheme-style schemes never duplicate
    #       them, and neither does this one.
    #     - captureFromShells/captureToShells renamed fromShells/toShells (shared by BOTH the DR capture
    #       and the parent-line electron-impact excitation, since it's the SAME core promotion populated
    #       two different ways); captureExcitationDegree dropped in favour of the already-existing
    #       Basics.ForDielectronicRecombination(fromShells,toShells,intoShells,decayShells) theme, reused
    #       directly rather than re-implemented.
    #     - referenceInitialConfigs/referenceFinalConfigs dropped -- the "reference/parent" excited
    #       configuration is just fromShells-->toShells applied to refConfigs, generated on the fly.
    #     - lineSelection (satellite-window) and plasmaModel not carried over: DielectronicRecombination
    #       .Settings and ImpactExcitation.Settings are embedded directly (drSettings/ieSettings) and
    #       reused AS-IS (their own pathwaySelection/lineSelection/temperatures/maxKappa/
    #       numElectronEnergies/maxEnergyMultiplier apply), rather than re-flattening copies of their
    #       fields onto the scheme.
    #     - item (5)'s gate (electron-impact excitation) is CLOSED: ImpactExcitation.computeLines works
    #       for this, once routed around two real bugs found today (module-ImpactExcitation.jl: a stale
    #       RateSettings type referenced in example-Dl.jl that no longer exists -- Settings itself already
    #       carries temperatures/calcRateCoefficient/maxKappa/etc. directly; and displayLines() crashing
    #       with printBefore=true because it runs before channels is populated -- worked around by always
    #       forcing printBefore=false in the driver). Fixing module-ImpactExcitation.jl itself is deferred
    #       to a separate session.
    #     - The parent line's "branching ratio" question from the sketch's design notes is answered
    #       DIFFERENTLY than assumed: NOT a competing-autoionization computation. Autoionization needs a
    #       final state with one FEWER bound electron (the ion after ejecting one); the bare excited level
    #       (fromShells-->toShells, no captured spectator) has the SAME electron count as refConfigs -- an
    #       ordinary valence excitation, not an ionizing process. A genuine autoionization channel only
    #       opens up when EI excites an INNER-SHELL electron of a system with enough electrons for a
    #       third one to be ejected (branch a's He-like case rules this out structurally; here it is a
    #       simplifying assumption, not verified) -- not attempted; branching fraction is fixed at 1
    #       (pure radiative), with an explicit @warn() in the driver documenting this limitation.
    #   Also debugged here (both fixes now live in the driver, so they apply equally to branch a):
    #     - DielectronicRecombination.Settings()'s own default has gauges=UseGauge[] (EMPTY) -- with no
    #       gauge selected, every DR photon rate silently comes out EmProperty(0.,0.), which then makes
    #       the whole R(Te) ratio NaN (0/0) once summed. The driver now defaults gauges to
    #       [UseCoulomb,UseBabushkin] whenever the caller leaves scheme.drSettings.gauges empty.
    #     - scheme.ieSettings.lineSelection pairs are (ground-index, excited-index), matching IE's own
    #       initialConfigs=refConfigs / finalConfigs=bareExcitedConfs direction; the (now-removed) Auger
    #       step went the OPPOSITE direction and needed each pair transposed -- reusing IE's pairs
    #       verbatim there silently targeted the wrong excited level.
    #   System: F-like Ne+ (Z=10, q=1), 2p --> 3d core excitation, captured shell 3s only (n'=3, the
    #   cheapest of three n'=3 captures explored today: 3s=220s, 3d=471s, 3p=665s standalone). Scoped to
    #   ONE captured shell, not the full n'=3 family (3s+3p+3d combined is known to exceed 41 min and was
    #   abandoned; the full family means separate top-level calls, one per shell, not attempted here).
    #   lineSelection targets level 2 of 2p^4 3d (J=5/2+, 0.0104 eV above level 1), NOT level 1 (J=7/2+):
    #   level 1 is E1-FORBIDDEN to both F-like ground levels (J=3/2-,1/2-; DeltaJ=2 exceeds E1's DeltaJ<=1
    #   rule). Of 2p^4 3d's 28 levels (all + parity), 21 are E1-allowed (J=1/2,3/2,5/2); level 2 is the
    #   lowest-energy allowed one -- same cost as level 1, physically valid.
    #   Grid: matches ForPedestrians.computeResonanceStrength's own default fallback (rbox=20), used in
    #   every standalone n'=3 test today, deliberately NOT rbox=90. rbox=90 was only ever needed for
    #   n'>=4 captured shells ("orbital not bound" at the default box); for n'=3 alone (this branch) the
    #   default box already worked cleanly. Using rbox=90 anyway was tried once and gave a real,
    #   reproducible numerical artifact: EVERY radiative (PhotoEmission) rate in the DR resonance table
    #   came back EXACTLY zero (Auger rates stayed nonzero), even though the identical physics on the
    #   default box gives real, nonzero radiative rates -- an unexplained bound-bound dipole-integral
    #   degradation on the much larger grid, not a driver bug. Root cause not investigated further; noted
    #   here as a real constraint for later (e.g. if intoShells is ever extended to n'>=4).
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    nm   = Nuclear.Model(10.)                              # F-like Ne+ (Z=10, q=1)

    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings();
                    calcRateAlpha = true,
                    temperatures  = [5e4, 1e5, 2e5, 5e5, 1e6] )                # [K]

    ieSettings = ImpactExcitation.Settings(ImpactExcitation.Settings();
                    calcRateCoefficient  = true,
                    temperatures         = Float64[],                         # left empty -- driver reconciles from drSettings
                    lineSelection        = LineSelection(true, indexPairs=[(1,2)]),
                    maxKappa             = 15,                                # convergence ceiling, not a fixed loop count
                    numElectronEnergies  = 3,
                    maxEnergyMultiplier  = 4.0 )

    scheme = Plasma.SatelliteDiagnosticScheme(Plasma.SatelliteDiagnosticScheme();
                 fromShells  = [Shell("2p")],
                 toShells    = [Shell("3d")],
                 intoShells  = [Shell("3s")],                                 # single shell for this first test
                 decayShells = [Shell("2p")],
                 drSettings  = drSettings,
                 ieSettings  = ieSettings )

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid,
                                     refConfigs=[Configuration("1s^2 2s^2 2p^5")], asfSettings=AsfSettings(),
                                     settings=Plasma.Settings() )

    wb = perform(computation, output=true)
    #
end
