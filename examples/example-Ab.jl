#
println("Ab) Apply & test the  a SCF field: B-spline primitives and one-particle spectra in a local potential.")
#

if  true
    # Last successful:  29-Jul-2026 (re-verified bit-identical after the ALField promotion/rename, see below)
    # Branch 0 (19-Jul-2026, migrated to ALFieldClaude2 28-Jul-2026, PROMOTED to plain ALField 29-Jul-2026):
    # the simplest possible non-trivial case for the full (non-local-exchange) average-level (AL) Dirac-Fock
    # SCF -- helium, 1s^2, a single CSF, a single subshell, a single "same-subshell doubled" two-electron
    # coefficient. Originally ran on Basics.ALFieldClaude() (the FIRST kink-aware code line); ALFieldClaude
    # was retired, superseded by the bVector-native line developed under the working name ALFieldClaude2
    # (ALFieldClaude's own SCF loop, solveAverageLevelFieldClaude, was Claude1-exclusive and deleted).
    # 29-Jul-2026: the ORIGINAL, buggy Basics.ALField() (the pre-Claude2 implementation shown for comparison
    # in the elseif branch just below) was retired/deleted entirely, and ALFieldClaude2 was renamed to take
    # over the plain Basics.ALField()/solveAverageLevelField name -- there is now only ONE ALField, and it is
    # the validated one. The underlying kink-aware radial-integral machinery (XL_CoulombClaude/
    # XL_CoulombTensorClaude/computeFunctionalClaude) is unchanged by this rename, just no longer
    # Claude2-suffixed at the scField/SCF-driver level. The new,
    # isolated code line (module-Basics-inc-abstract.jl, module-RadialIntegrals.jl,
    # module-InteractionStrength.jl, module-SelfConsistent.jl -- all additive, the original DFS/ALField code
    # paths are completely untouched) replaces the naive tensor-product Gauss-Legendre double sum used for
    # the two-electron Slater integral R^k(abcd) by a kink-aware construction: the "screened potential"
    # V_k(r) from the fixed orbital pair is built via a cubic spline (Dierckx) of the pair density plus an
    # explicit split of the inner integral at r=s, evaluated by adaptive Gauss-Kronrod quadrature (QuadGK) --
    # see RadialIntegrals.buildScreenedPotentialClaude for the full derivation. This was diagnosed as the
    # actual root cause of the "wrong by <1%" DF total-energy bug: the r_</r_>^(k+1) kernel has a kink
    # (continuous value, discontinuous first derivative) at r=s, which the original single, fixed-rule
    # quadrature cannot handle correctly whenever that kink falls inside a break-point cell --
    # disproportionately large for compact same-shell integrals like He's F^0(1s,1s).
    #
    # Reference values (helium ground state, 1s^2 1S0):
    #   - Non-relativistic Hartree-Fock limit:  E = -2.861680 Hartree (Clementi & Roetti / standard textbook
    #     value; essentially exact for a single-configuration HF calculation of He).
    #   - The exact (fully correlated) non-relativistic ground-state energy is -2.903724377 Hartree (Pekeris
    #     1959) -- NOT the target here (AL is single-configuration, no correlation), but useful as a sanity
    #     bound: the AL/DF result must lie ABOVE -2.9037 (variational bound) and close to -2.8617.
    #   - Relativistic (Dirac-Fock) corrections for Z=2 are tiny -- of order (Z*alpha)^2 ~ 2e-4 relative,
    #     i.e. a few times 1e-4 Hartree in absolute terms -- so the converged AL total energy is expected to
    #     sit just barely below -2.861680 Hartree, NOT differ from it at the 1e-2 or 1e-1 Hartree level.
    #   - The un-fixed ALField() result on this same (62-spline) grid was found to be WRONG-SIGN, i.e.
    #     slightly ABOVE -2.861680 (see the elseif branch just below for the reference run and its number);
    #     an energy-evaluation-only recompute (not a full re-optimized SCF) with the fixed integral flipped
    #     this to the physically-correct, more-negative side.
    #
    # RESULT (19-Jul-2026, confirmed on ALFieldClaude; RE-CONFIRMED 28-Jul-2026 on ALFieldClaude2 after the
    # migration; RE-CONFIRMED AGAIN 29-Jul-2026 after the ALField promotion/rename, bit-identical to 12
    # decimal places, exactly as expected for a pure rename with no logic change): the full SCF converges
    # cleanly and monotonically in 8 iterations (orbital-conv -> 1 to 9.6e-10). Final total energy, from the
    # "Level energies" table (which now also goes through the kink-aware integral for the CI step -- see
    # Hamiltonian.performCIClaude/setupMatrixClaude, added once this branch first showed a mismatch between
    # the SCF-internal energy and the CI-reported one; performCI itself is not touched, so DFS/EOL runs are
    # unaffected):
    #   E(He-AL) = -2.861813307613 Hartree  (29-Jul-2026, post-promotion, bit-identical to the 28-Jul-2026
    #              pre-rename number)
    #   Delta vs. non-rel. HF reference (-2.861680)     = -1.33e-4 Hartree  (correct sign: below HF limit)
    #   Delta vs. Pekeris exact non-rel. ground state    = well above -2.903724377, as required
    # This is the first fully self-consistent, correctly-signed Dirac-Fock total energy obtained for helium
    # in this code line -- the kink-quadrature diagnosis and fix are confirmed end to end.
    #
    wa = Atomic.Computation(Atomic.Computation(), name="He-AL", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(2.),
                            configs=[Configuration("1s^2")],
                            asfSettings=AsfSettings(AsfSettings(), scField=Basics.ALField(), accuracyScf = 1.0e-8, maxIterationsScf = 60)  )

    wb = perform(wa)

elseif  false
    # Last successful:  unknown -- HISTORICAL RECORD, no longer reproducible as written (see note below)
    # Branch 0-reference: originally the same He 1s^2 case as Branch 0 above, but with the ORIGINAL, un-fixed
    # Basics.ALField() -- kept here for direct comparison. Result on the 62-spline baseline grid (19-Jul-2026):
    # total energy = -2.861635343768 Hartree, i.e. +4.4656e-5 Hartree ABOVE the non-relativistic HF reference
    # of -2.861680 -- the wrong sign (a converged Dirac-Fock result must lie slightly BELOW the non-relativistic
    # HF limit, not above it), diagnosed as the two-electron kink-quadrature bug described in Branch 0.
    # RETIRED 29-Jul-2026: the buggy pre-fix Basics.ALField()/solveAverageLevelField implementation shown
    # here was deleted entirely and ALFieldClaude2 was promoted to take over the plain ALField name -- the
    # code below, if run today, would silently execute the NEW, validated implementation (same as Branch 0),
    # NOT reproduce the wrong number shown above. Kept only as a historical record of the bug that was fixed;
    # do not re-enable this branch expecting the old (wrong) result.
    wa = Atomic.Computation(Atomic.Computation(), name="He-AL", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(2.),
                            configs=[Configuration("1s^2")],
                            asfSettings=AsfSettings(AsfSettings(), scField=Basics.ALField(), accuracyScf = 1.0e-8, maxIterationsScf = 60)  )

    wb = perform(wa)

elseif  false
    # Last successful:  unknown
    # Branch 1 (was: `if true`): the original neon/ALField case that started this debugging session.
    # State as of 19-Jul-2026, before switching to the simpler He case above:
    #   - Fixed two missing `using` imports in module-SelfConsistent.jl (InteractionStrength, RadialIntegrals)
    #     that made this branch crash outright before even one SCF iteration could run -- these were real,
    #     but are very unlikely to be THE decades-old "wrong by <1%" bug, since they blocked execution
    #     entirely rather than giving a slightly-wrong number. Most likely a newer regression from the
    #     SelfConsistent module restructuring mentioned earlier this session (cf. project memory).
    #   - With those two fixed, the SCF total energy OSCILLATES rather than converges cleanly: -124.667 ->
    #     -122.021 -> -123.130 -> -122.707 -> -122.877 Hartree over 5 iterations, while the loop's own
    #     accuracy metric (orbital-conv, an orbital-overlap product) drops below accuracyScf and the loop
    #     stops anyway -- i.e. the stopping criterion does not track whether the TOTAL ENERGY has actually
    #     stabilized, only whether individual orbitals stopped changing overlap-wise between iterations.
    #   - The FINAL reported "Level energies" total (-122.7067 Hartree) matches iteration 4's functional
    #     value, not iteration 5's (-122.8771, the last one actually computed) -- a separate, additional
    #     stale-orbitals/off-by-one issue in how the converged orbitals are handed to the final CI step.
    #   - Neither number is within "less than 1%" of the expected ~-128.6 to -128.9 Hartree range for Ne
    #     (this is roughly 4.5% off) -- substantially bigger than what the user remembers, suggesting the
    #     code has degraded further since the original "1%" bug was last characterized, and that oscillation
    #     issue needs fixing before Ne can be used to isolate the original, subtler discrepancy again.
    # Resume here once the He 1s^2 case above is fully understood.
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne]")],
                            asfSettings=AsfSettings(AsfSettings(), scField=Basics.ALField(), accuracyScf = 1.0e-4)  )

    wb = perform(wa)

elseif  true
    # Last successful:  unknown
    # Compute different direct potentials for the charge density of a given level 
    # (a Radial.Grid(...; NoPoints = 900) call stood here until 29-Aug-2026: the keyword became `rbox`, and the
    #  line was dead in any case -- the next line reassigns `grid` before it is ever used.)
    grid = Radial.Grid(true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(16., PointNucleus()), 
                            configs=[Configuration("[Ne] 3s^2 3p^2")],  ## , Basics.DFSField()
                            asfSettings=AsfSettings(AsfSettings(), scField=Basics.DFSField())  )

    wb = perform(wa)

elseif false
    # Last successful:  unknown
    # Test for Björn, 7. Mai 2020
    ## grid=JAC.Radial.Grid(true)
    # `NoPoints` was retired as a Radial.Grid keyword; `rbox` replaced it. rbox = 85.0 reproduces the
    # ~2000 points this line asked for (1911 at 80, 2114 at 90) with rnt/h/hp unchanged. This is the
    # SECOND occurrence in this file -- item 79 repaired the one at line 115 on 29-Aug and missed this one.
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5, h = 5.0e-2, hp = 5.0e-2, rbox = 85.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(8.0, FermiNucleus()), 
                            configs=[Configuration("1s^2 2s"), Configuration("1s 2s^2"), Configuration("1s 2s 3p"), Configuration("1s 2s 4p")], 
                            ## configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s 2s 3s"), Configuration("1s 2p 3s"), 
                            ##          Configuration("1s 2s 2p")],  ## 
                            # WAS a fully POSITIONAL AsfSettings of the retired shape -- EIGHTEEN values into a struct
                            # that now has FOURTEEN fields, so it could not be remapped position by position. Only the
                            # values that are unambiguous are carried over: DFSField, hydrogenic start, 40 iterations,
                            # 1.0e-6, NoneQed, LSjjSettings(false). DELIBERATELY DROPPED because the old positions
                            # cannot be identified with confidence: two level lists, [1] and [1,2,3,4], of which one
                            # was presumably levelSelectionCI; and the four bare Bools. The branch is undated, so
                            # nothing rests on them. The same file already uses this keyword form at four other calls.
                            asfSettings=AsfSettings(AsfSettings(); scField = Basics.DFSField(),
                                                    startScfFrom = StartFromHydrogenic(), maxIterationsScf = 40,
                                                    accuracyScf = 1.0e-6, qedModel = NoneQed(),
                                                    jjLS = LSjjSettings(false))  )

    wb = perform(wa)

end
