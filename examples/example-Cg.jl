
println("Cg) Apply & test the DecayYield module with ASF from an internally generated multiplet.")
println("    Branches follow Koziol, arXiv:2204.13428 (2022), 'Relativistic calculated K-shell level")
println("    widths and fluorescence yields for atoms with 20<=Z<=30' (MCDF/RATIP, same lineage as JAC),")
println("    itself cross-checked against Krause (1979), Bambynek (1972), and recent experiments;")
println("    see project memory project_decayyield_fluorescence.md.")

setDefaults("print summary: open", "zzz-DecayYield.sum")

if  false
    # Last successful:  31-Jul-2026
    # Branch a: Mg (Z=12) K-shell hole -- basic functionality test, NOT a literature comparison.
    #   BUG FOUND (found this session, fixed here): the original stub used configs=["1s 2s^2 2p^6"], only
    #   9 electrons -- MISSING Mg's own valence 3s^2 shell entirely (neutral Mg has 12 electrons; a K-hole
    #   ion should have 11: "1s 2s^2 2p^6 3s^2"). The original config silently modeled a stripped, mostly
    #   fictitious 9-electron ion, not a real Mg K-hole. Also, the original call
    #   DecayYield.Settings("SCA", true, LevelSelection()) (3 args) no longer matched the (already-drifted)
    #   4-field Settings struct.
    #   SETTINGS REFACTORED (this session): `DecayYield.Settings.approach` was a raw String
    #   ("AverageSCA"/"SCA") that got re-decoded into the real cascade-approach singleton types deep inside
    #   computeOutcomes. Moving to `approach::Basics.AbstractCascadeApproach` directly surfaced a genuine
    #   CIRCULAR MODULE DEPENDENCY: Cascade-inc-computations.jl already needs DecayYield.Outcome/Settings as
    #   compile-time function-argument types, so DecayYield cannot compile-time-depend on Cascade.
    #   AbstractCascadeApproach and its three singletons (AverageSCA, SCA, UserMCA) were therefore moved OUT
    #   of module-Cascade.jl and INTO module-Basics-inc-abstract.jl (Basics.AbstractCascadeApproach etc.) --
    #   the same early-loaded-shared-supertype pattern JAC already uses for AbstractPropertySettings/
    #   AbstractProcessSettings. Cascade.jl re-exposes them unchanged via its existing `using ..Basics`
    #   (both bare AverageSCA()/SCA() and self-qualified Cascade.AverageSCA()/Cascade.SCA() usage throughout
    #   Cascade's own included files keep working with no further edits there). Also added a new
    #   `decayShells::Array{Shell,1}` field, UNIONed with the shells auto-extracted from the initial
    #   configuration(s) -- needed to explicitly communicate a spectator shell that is entirely ABSENT from
    #   a double-core-hole configuration (e.g. "1s^2 2p^5" needs 2s added explicitly, since it never appears
    #   in that configuration string at all, yet can become populated during the decay cascade). Empty by
    #   default (Shell[]), so existing behavior is unchanged unless used -- see branch e.
    #   VERIFIED bit-identical to the pre-refactor run: omega_R=3.1704e-2/3.2604e-2 (Cou/Bab),
    #   omega_A=9.6830e-1/9.6740e-1, omega_R+omega_A=1.0 exactly (both gauges) -- confirms this was a purely
    #   structural/type-safety change with zero physics impact.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Cg-a-Mg-Khole", grid=grid, nuclearModel=Nuclear.Model(12.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection()) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch b: Ca (Z=20) K-shell hole, "1s 2s^2 2p^6 3s^2 3p^6 4s^2" (19 electrons; neutral Ca [Ar]4s^2
    #   minus one 1s electron) -- matches Koziol's own valence configuration "4s^2" for Ca exactly.
    #   Target omega_K (Koziol Table VI/VII, four MCDF/RATIP model variants C+F/B+F/C+R/B+R):
    #   0.163 (C+F), 0.170 (B+F), 0.154 (C+R), 0.161 (B+R). No experimental point at Z=20 in Koziol's
    #   Table VII.
    #   VERIFIED: omega_R=0.1703/0.1720 (Cou/Bab), omega_A=0.8297/0.8280 -- omega_R+omega_A=1.0 exactly, both
    #   gauges. JAC's values (0.170-0.172) sit right at/just above the top of Koziol's 0.154-0.170 spread --
    #   good agreement (~1-10% depending on which of the 4 Koziol variants is compared), consistent with JAC
    #   using its own SCF/Auger-rate treatment rather than Koziol's exact MCDF/Fac-Ratip combination.
    #   NOTE: Radial.Grid(true) sets hp=0, breaking Auger continuum-orbital generation -- fixed via a custom
    #   grid, same style as branch a/e's working Mg/Ne grids.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Cg-b-Ca-Khole", grid=grid, nuclearModel=Nuclear.Model(20.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^6 4s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection()) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  unknown ... ABANDONED (see REMARK below), not a dating candidate.
    # Branch c: Fe (Z=26) K-shell hole, "1s 2s^2 2p^6 3s^2 3p^6 3d^6 4s^2" (25 electrons; neutral Fe
    #   [Ar]3d^6 4s^2 minus one 1s electron) -- matches Koziol's own valence configuration "3d^6 4s^2" for
    #   Fe exactly. This is the most literature-rich comparison point available: Koziol's own MCDF/RATIP
    #   models give omega_K = 0.335 (B+R) - 0.356 (C+F); independent recent experiments in the same table
    #   give Durak&Ozdemir 0.330+/-0.005, Simsek 0.331+/-0.012, Sogut 0.366+/-0.033, Han 0.358+/-0.029 --
    #   all clustering tightly around 0.33-0.36.
    #   NOTE: Radial.Grid(true) sets hp=0, breaking Auger continuum-orbital generation -- fixed via a custom
    #   grid. Fe's 3d^6 is a genuinely open shell (unlike Ca/Zn), so this K-hole configuration produces
    #   MULTIPLE levels -- was going to also be used to test levelSelection explicitly (all levels vs.
    #   LevelSelection(true, indices=[1])) but never got that far -- see REMARK.
    #   REMARK (31-Jul-2026): user killed this run after 140 minutes, no result -- an open 3d^6 shell right
    #   at the K-hole level (Fe's own valence, not just a spectator elsewhere) creates an intractably large
    #   CSF/Auger-final-state space for this stepwise-decay approach, consistent with Koziol's own Table I
    #   showing Fe's SECONDARY (post-Auger) hole-state counts already in the thousands (L^-1M^-1=2172,
    #   M^-1N^-1=1160, etc.) even before JAC's own CSF generation/angular-coefficient machinery multiplies
    #   that further. NOT a good example for this kind of quick literature-comparison branch -- an
    #   open-valence-shell K-hole case is simply too expensive here, independent of any bug. Abandoned;
    #   branch d (Zn, CLOSED 3d^10) tried next as a likely more tractable substitute. If a Fe-like case is
    #   wanted later, a MUCH more restrictive decayShells/decay-configuration space would likely be needed
    #   first -- natural follow-up once the Cascade module gets its planned longer revisit.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Cg-c-Fe-Khole-allLevels", grid=grid, nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^6 4s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection()) ] )
    perform(wa)

    wa2  = Atomic.Computation(Atomic.Computation(), name="Cg-c-Fe-Khole-level1only", grid=grid, nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^6 4s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(true, indices=[1])) ] )
    perform(wa2)
    #
elseif  false
    # Last successful:  31-Jul-2026
    # Branch d: Zn (Z=30) K-shell hole, "1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2" (29 electrons; neutral Zn
    #   [Ar]3d^10 4s^2 minus one 1s electron; 3d^10 is a CLOSED subshell, so this system's CSF space is much
    #   smaller than Fe's despite the higher Z -- CONFIRMED tractable, unlike branch c's abandoned open-3d^6
    #   Fe attempt) -- matches Koziol's valence configuration "3d^10 4s^2" for Zn exactly.
    #   Target omega_K (Koziol Table VI/VII): 0.470 (C+R) - 0.487 (B+F); experiments in the same table:
    #   Yashoda 0.471+/-0.018, Durak 0.482+/-0.032, Gudennavar 0.464+/-0.010, Simsek 0.482+/-0.022,
    #   Sogut 0.525+/-0.050, Han 0.477+/-0.038 -- all clustering around 0.46-0.49 (Sogut a bit higher,
    #   largest error bar).
    #   VERIFIED: omega_R=0.4911/0.4909 (Cou/Bab), omega_A=0.5089/0.5091 -- omega_R+omega_A=1.0 exactly, both
    #   gauges. JAC's value sits just above Koziol's own theoretical band (0.470-0.487) but comfortably
    #   within the experimental scatter, close to Durak (0.482+/-0.032) and Simsek (0.482+/-0.022) within
    #   their error bars -- good agreement.
    #   NOTE: Radial.Grid(true) sets hp=0, breaking Auger continuum-orbital generation -- fixed via a custom
    #   grid, same style as branch a/b/e's working grids.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="Cg-d-Zn-Khole", grid=grid, nuclearModel=Nuclear.Model(30.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection()) ] )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  31-Jul-2026 (all 4 sub-computations ran cleanly; decayShells' effect is CONFIRMED
    #   ABSENT for StepwiseDecayScheme, root cause found and documented below, not a fixed positive result)
    # Branch e: Mg (Z=12) double-core-hole "1s 2p^6 3s^2" -- explicit decayShells test, per user request/
    #   suggestion. This configuration has a 1s HOLE (1 electron, the K-hole to be decayed), 2s COMPLETELY
    #   ABSENT (0 electrons, not even written in the configuration string -- a genuine second, spectator
    #   core hole), 2p FULL, and 3s FULL (a real electron source/donor for the 1s hole's radiative/Auger
    #   decay, unlike the earlier "1s^2 2p^5" attempt which had no donor shell above the hole at all and
    #   gave an uninformative null result regardless of decayShells).
    #   ROOT CAUSE FOUND (this session, before rerunning): Cascade.StepwiseDecayScheme.decayShells is NOT
    #   consulted by the currently-active decay-configuration generator (Basics.generateConfigurations(
    #   Basics.ForStepwiseDecay(N), initialConfigs), module-Cascade-inc-stepwise-decay.jl:350) -- it decides
    #   which shells may participate purely via haskey(conf.shells, shell) (Basics.ValenceShells(),
    #   module-BasicsAZ-inc-configurations.jl:1120), i.e. whether the shell is an explicit KEY in the
    #   Configuration's own shell dict, regardless of occupation value. A shell absent from the configuration
    #   string entirely (like 2s here) fails haskey and is silently skipped no matter what decayShells says.
    #   The OLDER, now-disabled (Oct-2025) Cascade.generateConfigurationsForStepwiseDecay DID insert
    #   zero-occupation decayShells entries into the configuration before generating decay configs; its
    #   replacement dropped this step. FIXED in module-DecayYield.jl (NOT in Cascade.jl, confined to this one
    #   module, per Rule 6): DecayYield.computeOutcomes now pre-inserts zero-occupation entries for every
    #   settings.decayShells shell directly into the Configuration objects passed to Cascade.Computation,
    #   restoring the old behavior locally without touching Cascade internals (a longer, separate revisit of
    #   the Cascade module is planned for later).
    #   Two settings compared on the SAME initial configuration:
    #   (1) decayShells=Shell[] (default) -- 2s never becomes a valence shell, cannot participate at all.
    #   (2) decayShells=[Shell("2s")] -- 2s is injected as an explicit zero-occupation entry beforehand.
    # NOTE: Radial.Grid(true) (the plain default) sets hp=0, which breaks continuum-orbital generation --
    # required for the Auger channels DecayYield always computes -- with "Improper grid for continuum
    # processes with grid.hp = 0." (also hit by branch b, Ca, before this was diagnosed). A custom grid with
    # an explicit nonzero hp is required for ANY DecayYield/StepwiseDecayScheme computation; reusing branch
    # a's grid here since it is the same Z=12 Mg.
    #
    # RESULT (Mg sub-test): IDENTICAL for both settings (No_R=2, No_A=10, omega_R=4.2330e-2/4.3723e-2,
    # omega_A=9.5767e-1/9.5628e-1, Cou/Bab). Understood, not a failure: Basics.ValenceShells() always
    # identifies the FIRST under-full shell in canonical order as "the hole to fill" (AddElectrons targets
    # ONLY valenceShells[1]); everything else is only ever a SOURCE (RemoveElectrons). Here 1s is already
    # correctly identified as the hole (full-shell scan reaches it first), so injecting empty 2s only adds it
    # to the source pool -- and you cannot remove an electron from an already-empty shell, so it does
    # nothing. decayShells can only matter when the injected shell would itself become the FIRST identified
    # hole (i.e. it must sit BEFORE the real hole in canonical shell order) -- which is exactly the ORIGINAL
    # Ne "1s^2 2p^5" case (2s precedes 2p canonically), re-tested below with the same grid fix.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    nm12 = Nuclear.Model(12.)
    wa1 = Atomic.Computation(Atomic.Computation(), name="Cg-e-Mg-Khole-noShell2s", grid=grid, nuclearModel=nm12,
                            configs=[Configuration("1s 2p^6 3s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(), decayShells=Shell[]) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cg-e-Mg-Khole-withShell2s", grid=grid, nuclearModel=nm12,
                            configs=[Configuration("1s 2p^6 3s^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(), decayShells=[Shell("2s")]) ] )
    perform(wa2)

    # Ne (Z=10) "1s^2 2p^5" re-test, same grid style as the TestFrames Cascade-StepwiseDecay test for Ne.
    # Here 2s (canonically BEFORE 2p) should, once injected, be identified as the hole ITSELF (displacing
    # 2p), not merely added as an inert spectator -- a fundamentally different calculation, not just a
    # bigger valence-shell list.
    # RESULT: STILL IDENTICAL (both zero) -- deeper root cause found: Cascade.perform(scheme::
    # StepwiseDecayScheme, comp) does NOT use comp.initialConfigs (the one my DecayYield.jl fix modifies)
    # directly for decay-configuration generation. It first runs SCF on it, then RE-EXTRACTS a configuration
    # from the resulting basis via Basics.extractConfigurations(Basics.FromMultiplet(), multiplets), which
    # bottoms out in Basics.extractConfigurations(Basics.FromBasis(), basis)
    # (module-BasicsAZ-inc-configurations.jl:655): `if occ > 0  newShells = merge(...)`  -- a shell is kept
    # ONLY if some CSF gives it NONZERO occupation. My injected "2s=>0" is faithfully passed into SCF, but
    # since no CSF in this single-configuration basis ever populates 2s, this re-extraction step drops it
    # again right before Basics.ForStepwiseDecay ever sees it. So the fix IS necessary (ensures the initial
    # SCF/level computation itself is aware of the spectator hole) but NOT sufficient to make decayShells
    # influence which decay configurations get generated -- that gap lives in the SHARED
    # module-BasicsAZ-inc-configurations.jl extraction machinery, not in DecayYield.jl or even Cascade.jl
    # specifically. Per explicit user decision, NOT pursued further this session -- deferred to a later,
    # dedicated revisit of the Cascade module. decayShells IS confirmed to work correctly for OTHER cascade
    # schemes that consult it directly (DielectronicRecombinationScheme, HollowIonScheme) -- only
    # StepwiseDecayScheme/DecayYield is affected by this gap.
    gridNe = Radial.Grid(Radial.Grid(false), rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, rbox = 9.5)
    nm10   = Nuclear.Model(10.)
    wa3 = Atomic.Computation(Atomic.Computation(), name="Cg-e-Ne-2p5-noShell2s", grid=gridNe, nuclearModel=nm10,
                            configs=[Configuration("1s^2 2p^5")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(), decayShells=Shell[]) ] )
    perform(wa3)

    wa4 = Atomic.Computation(Atomic.Computation(), name="Cg-e-Ne-2p5-withShell2s", grid=gridNe, nuclearModel=nm10,
                            configs=[Configuration("1s^2 2p^5")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(), decayShells=[Shell("2s")]) ] )
    perform(wa4)
    #
elseif  true
    # Last successful:  31-Jul-2026
    # Branch f: Si (Z=14) K-shell hole, "1s 2s^2 2p^6 3s^2 3p^2" (13 electrons; neutral Si [Ne]3s^2 3p^2
    #   minus one 1s electron) -- explicit levelSelection test, per user request. Si's open 3p^2 shell gives
    #   MULTIPLE levels (3P0,1,2-derived terms coupled with the K-hole) with only 2 valence electrons in an
    #   open shell, vs. Fe's 6 in an open d-shell (branch c, abandoned as computationally intractable after
    #   a 140-minute timeout) -- chosen specifically to be a much smaller, tractable multi-level case. No
    #   literature target needed; this is purely a mechanics test: run once with the default (all levels),
    #   once restricted to LevelSelection(true, indices=[1]) (just the lowest level), and confirm the
    #   restricted run reproduces exactly the level-1 row of the full run.
    #   VERIFIED: full run gives 6 levels; level 1 (J=1/2+, E=-6.02721096e+03 eV): No_R=16, No_A=180,
    #   omega_R=5.3512e-02/5.4721e-02, omega_A=9.4649e-01/9.4528e-01 (Cou/Bab, sum=1.0 exactly). The
    #   level1-only run reproduces this row BIT-IDENTICALLY -- confirms levelSelection both correctly
    #   restricts output and reproduces exactly the same physics for the selected level, no side effects.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    wa1 = Atomic.Computation(Atomic.Computation(), name="Cg-f-Si-Khole-allLevels", grid=grid, nuclearModel=Nuclear.Model(14.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection()) ] )
    perform(wa1)

    wa2 = Atomic.Computation(Atomic.Computation(), name="Cg-f-Si-Khole-level1only", grid=grid, nuclearModel=Nuclear.Model(14.),
                            configs=[Configuration("1s 2s^2 2p^6 3s^2 3p^2")],
                            propertySettings=[ DecayYield.Settings(DecayYield.Settings(); approach=Basics.SCA(), printBefore=true,
                                                geant4=false, levelSelection=LevelSelection(true, indices=[1])) ] )
    perform(wa2)
    #
end
#
setDefaults("print summary: close", "")


