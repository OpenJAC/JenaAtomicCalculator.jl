#
println("Ar) Apply & test the second-order treatment of the Q space of a RAS step.")
#
# A RAS step admits a correlation class by putting every one of its CSFs into the CI, or not at all.  This file
# is about the third possibility: rank the configurations a step would add, promote the strongly coupled ones
# into the CI, fold the weakly coupled remainder in through second-order perturbation theory, and drop the rest.
# The two thresholds of Basics.SecondOrder(promoteAbove, discardBelow) decide which, and the step REFUSES to
# stay silent when its own precondition fails -- see branch (b).
#
# WHAT P AND Q ARE HERE, because the whole scheme rests on it.  P is the space the previous VARIATIONAL step
# built: correlation decided by physics stays decided, and no threshold re-litigates it.  Q is what the
# perturbative step adds OVER THE SAME ORBITALS -- it introduces no new shells, so nothing is optimized and no
# SCF runs.  A perturbative step is therefore a LEAF: it says what the excitations left out of the previous step
# are still worth, and is never inherited by a later variational step.
#
# AND Q IS FINITE.  It is the set of CSFs the current orbital set can build but the previous step excluded, not
# the complete second-order sum over an infinite virtual space.  A correction from here is a LOWER BOUND on the
# true second-order correction; a new layer ENLARGES Q rather than shrinking it.

if  false
    # Last visit:  12-Sep-2026
    # Last successful:  12-Sep-2026 -- step 1 -459.68354038, step 2 -459.74254091, step 3 -459.93473441 Ha;
    #                   splittings 22628.6/37880.4, 20295.7/31601.4, 22764.4/56203.3 cm^-1.
    #                   P = 66 CSF (5 configurations from step 2), Q = 842 CSF in 11 configurations; 3 of the
    #                   66 levels in P carry the reference and set the weights; 8 configurations promoted,
    #                   3 folded, 0 discarded, so the CI uses 799 of 908 CSF.  Sum |c|^2 over the folded part
    #                   = 1.13e-06 / 1.20e-06 / 1.10e-06 -- a weak perturbation, and the step says so.
    #
    #    VERIFIED AGAINST GROUND TRUTH, which is what makes this branch worth keeping.  Re-run with
    #    SecondOrder(1.0e-30, 1.0e-40) every configuration is promoted, the CI uses all 908 CSF and the scheme
    #    reduces to an ordinary CI: -459.93426153 Ha, splittings 22773.1 / 56286.4 cm^-1.  Against that,
    #    FOLDING COSTS 0.47 mHa ON THE TOTAL AND -8.7 / -83.1 cm^-1 ON THE SPLITTINGS (0.04 % and 0.15 %),
    #    for 799 CSF instead of 908.  The folded answer lies BELOW the exact one: second order over-binds,
    #    as it must.
    #
    #    THE SAVING HERE IS ONLY 12 %, AND THAT IS THE HONEST RESULT FOR THIS ION rather than a disappointment.
    #    2s,2p is the shell ADJACENT to the 3s,3p valence, so its correlation is not weak and the threshold
    #    rightly promotes most of it -- which is the criterion rediscovering GRASP's own rule that the upper
    #    one or two core layers should be opened rather than folded.  The gain belongs to DEEP cores, where
    #    |c|^2 ~ 1/D^2 is small: see branch (c).
    #
    # a) Cl III 1s^2 2s^2 2p^6 3s^2 3p^3, J = 3/2 odd.  Three steps:
    #      step 1   the reference alone
    #      step 2   VALENCE correlation 3s,3p -> 3d, variational.  The orbitals are optimized HERE.
    #      step 3   the CORE excitations 2s,2p -> 3d, folded.  NO new shells: the same 3s,3p,3d, everything
    #               frozen, and only the excitation pattern is widened.  The delta against step 2 IS the Q space.
    #
    #    NOTE WHY step 3 NAMES BOTH THE CORE AND THE VALENCE in its from-list.  A core-VALENCE double needs one
    #    electron from each, and generateBasis makes doubles by applying the single-excitation generator twice;
    #    naming only the core would therefore generate core-core alone and CV could never appear.  The driver
    #    subtracts what step 2 already built, so the widened list costs nothing.
    #
    Z      = 17.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]
    sym    = LevelSymmetry(AngularJ64(3//2), Basics.minus)
    val    = [Shell("3s"), Shell("3p")];    core = [Shell("2s"), Shell("2p")]
    to     = [Shell("3s"), Shell("3p"), Shell("3d")]
    frozen1 = [Shell("1s"), Shell("2s"), Shell("2p")]
    frozen2 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p"), Shell("3d")]
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to, frozen=frozen1),
               RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                  frozen=frozen2, treatment = Basics.SecondOrder(1.0e-5, 1.0e-12)) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 20.)
    rasSettings = RasSettings(Int64[], 24, 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Cl III -- valence variational, core folded on the same orbitals",
                            Nuclear.Model(Z), grid, refs, RasExpansion([sym], 15, steps, rasSettings) )
    println("wa = $wa")
    wb = generate(wa, output=true)
    for i = 1:length(steps)
        k = "step" * string(i)
        if  haskey(wb, k)
            ls = sort(wb[k].levels, by = l -> l.energy)
            println(">> step $i : $(length(ls)) levels, lowest = $(ls[1].energy) Ha" *
                    (length(ls) > 2 ? ";  splittings " *
                     string(round((ls[2].energy-ls[1].energy)*219474.6313702, digits=1)) * " / " *
                     string(round((ls[3].energy-ls[1].energy)*219474.6313702, digits=1)) * " cm^-1" : ""))
        end
    end

elseif  false
    # Last visit:  12-Sep-2026
    # Last successful:  12-Sep-2026 -- the CI uses all 908 CSF; -459.93426153 Ha, 22773.1 / 56286.4 cm^-1,
    #                   identical to a plain RAS over the same space.
    #
    # b) THE GUARD.  The same case as (a) with the thresholds set so low that EVERY configuration is promoted.
    #    The scheme must then reduce EXACTLY to an ordinary CI -- nothing folded, nothing discarded -- so this
    #    branch has an exact expected answer rather than a plausible one, and it is what would catch a
    #    regression in the ranking, the promotion or the effective Hamiltonian.  Run it after any change to
    #    Hamiltonian.performCIwithSecondOrderQ.
    #
    Z      = 17.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]
    sym    = LevelSymmetry(AngularJ64(3//2), Basics.minus)
    val    = [Shell("3s"), Shell("3p")];    core = [Shell("2s"), Shell("2p")]
    to     = [Shell("3s"), Shell("3p"), Shell("3d")]
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to,
                                  frozen=[Shell("1s"), Shell("2s"), Shell("2p")]),
               RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                  frozen=[Shell("1s"),Shell("2s"),Shell("2p"),Shell("3s"),Shell("3p"),Shell("3d")],
                                  treatment = Basics.SecondOrder(1.0e-30, 1.0e-40)) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 20.)
    rasSettings = RasSettings(Int64[], 24, 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Cl III -- the guard: promote everything, reproduce the exact CI",
                            Nuclear.Model(Z), grid, refs, RasExpansion([sym], 15, steps, rasSettings) )
    println("wa = $wa")
    wb = generate(wa, output=true)
    ls = sort(wb["step3"].levels, by = l -> l.energy)
    println(">> step 3 : $(length(ls)) levels, lowest = $(ls[1].energy) Ha  -- EXPECTED -459.93426153 Ha")

elseif  true
    # Last visit:  12-Sep-2026
    # Last successful:  unknown -- AND DELIBERATELY SO; this branch is NOT yet trustworthy.  It runs, the gate
    #                   correctly REFUSES (Sum |c|^2 = 6.26e-05 against a 1.0e-5 threshold) and names what to
    #                   promote, and the reference level is sound: -1183.042412 / -1183.132879 / -1184.092930 Ha
    #                   with w_ref 1.000 / 0.974 / 0.939, i.e. -0.091 Ha of valence and -0.960 Ha of core
    #                   correlation, both the right order for Fe.
    #                   BUT A SPURIOUS INTRUDER APPEARS: a 2p-hole configuration lands 10 Ha BELOW the 3s^2
    #                   reference, which is impossible -- the hole must cost ~28 Ha.  The cause is that step 2's
    #                   correlation orbitals are CONTRACTED pseudo-orbitals ("4s" at <r> = 0.698 a.u., on top of
    #                   the 3s valence) and the perturbative step reuses them to build CORE-HOLE CSFs, where they
    #                   are unphysical.  The gate cannot see this: it tests the FOLDED remainder, not whether the
    #                   PROMOTED configurations are sound.  See priority item 21.  Do not quote numbers from this
    #                   branch until that is resolved.
    #
    # c) Fe XV 1s^2 2s^2 2p^6 3s^2, J = 0 even -- THE CASE THE METHOD IS FOR.
    #
    #    Branch (a) saved only 12 %, and that was the honest answer for Cl III: its "core" 2s,2p sits directly
    #    below the 3s,3p valence, so the coupling is not weak and the threshold rightly promotes most of it.
    #    The weight goes as |c|^2 ~ 1/D^2, so the fold-versus-promote boundary is set by how FAR the shell sits
    #    from the valence -- and the gain therefore belongs to DEEP cores, which is to say to heavy ions, which
    #    is exactly where the variational alternative is least affordable.
    #
    #    Fe XV is Mg-like: the 3s valence sits over a 1s,2s,2p core that is tens of a.u. away, against ~10 a.u.
    #    for Cl III.  And the arithmetic of the space is stark -- valence correlation into n=3 and n=4 gives
    #    17 CSF, the same with the core opened gives 334.  The core excitations are 95 % of the space.
    #    If they fold, that 95 % is had for the cost of a ranking pass.
    #
    Z      = 26.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2")]
    sym    = LevelSymmetry(0, Basics.plus)
    val    = [Shell("3s")];   core = [Shell("2s"), Shell("2p")]
    to     = [Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to,
                                  frozen=[Shell("1s"), Shell("2s"), Shell("2p")]),
               RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                  frozen=vcat([Shell("1s"),Shell("2s"),Shell("2p")], to),
                                  treatment = Basics.SecondOrder(1.0e-5, 1.0e-12)) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 10.)
    rasSettings = RasSettings(Int64[], 24, 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Fe XV 3s^2 -- a DEEP core folded on the valence orbitals",
                            Nuclear.Model(Z), grid, refs, RasExpansion([sym], 12, steps, rasSettings) )
    println("wa = $wa")
    wb = generate(wa, output=true)
    # REPORT THE REFERENCE LEVEL, NOT THE LOWEST ONE.  Once the core is opened, a core-excited configuration
    # can sink below the reference in the energy ORDER, and then "lowest" names the wrong state -- the same
    # trap SelfConsistent.selectTargetLevelsEOL exists to avoid.  The step diagnostic prints w_ref per level;
    # here the level is picked by that weight.
    for i = 1:length(steps)
        k = "step" * string(i)
        haskey(wb, k)  ||  continue
        mp   = wb[k];   b = mp.levels[1].basis
        rIdx = [ r for r = 1:length(b.csfs)
                     if Basics.extractConfiguration(Basics.FromBasis(), b, b.csfs[r]) in refs ]
        best = nothing;  bw = -1.0
        for lev in mp.levels
            w = sum(lev.mc[r]^2  for r in rIdx; init=0.0)
            if  w > bw    bw = w;   best = lev    end
        end
        println(">> step $i : $(length(mp.levels)) levels;  the REFERENCE level is at $(best.energy) Ha " *
                "with w_ref = $(round(bw, digits=5))")
    end
end
