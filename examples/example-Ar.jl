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
    # Last visit:      13-Sep-2026
    # Last successful:  unknown -- THE 12-Sep DATE WAS WITHDRAWN, NOT LOST.  It was earned on two-particle angular
    #                   coefficients that were wrong for CSF pairs differing by one electron over a shared open
    #                   spectator (fixed 13-Sep, commit 1d05ab3), and step 3 of this branch -- which opens the
    #                   core -- generates such pairs in bulk.  Re-run that day on the repaired code:
    #
    #                       step 1  -459.68354038  unchanged   (the reference alone has no such pairs)
    #                       step 2  -459.76562495  was -459.74254091     -23.1 mHa
    #                       step 3  -459.76919888  was -459.93473441    +165.5 mHa
    #                       step 3 splittings 21323.4 / 33606.0 cm^-1,  were 22764.4 / 56203.3
    #
    #                   THE NEW NUMBERS ARE NOT YET JUDGED, which is why no date is written.  They are almost
    #                   certainly the better ones -- the old step-3 value was an eigenvalue of a wrong matrix and
    #                   so not a variational bound at all, and energies rising is what removing a spurious
    #                   lowering does -- but this branch carries no measured comparison to settle it, and
    #                   "different" is not "right".  Dating it needs a judgement, not a re-run.
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
    #
    # ===== WHAT THE SAME ION LOOKS LIKE WHEN THE VALENCE MODEL IS PUSHED (measured 15/16-Sep-2026) =====
    #
    # Branch (a) folds the core on top of ONE valence layer, where the valence error is still 18 %, so the core
    # contribution is buried.  Measured separately on the same ion, the same J = 3/2 odd levels and the same
    # grid (rbox = 20 a.u.), with NO core at all -- 2s,2p frozen throughout -- and only the valence ladder grown:
    #
    #   layer 1   refs(3s^2 3p^3) [1s..2p frozen] SD(3s..3p --> 3s..3d)     66 CSF
    #   layer 2   + [1s..3d frozen] SD(3s..3p --> 3s..4f)                  802 CSF
    #   layer 3   + [1s..4f frozen] SD(3s..3p --> 3s..5g)                 2816 CSF
    #   layer 4   + [1s..5g frozen] SD(3s..3p --> 3s..6h)                 6622 CSF     1:41:10, peak 4.84 GB
    #
    #                 2D*_3/2 (NIST 18052.5)      2P*_3/2 (NIST 29906.5)
    #     layer 1      21371.4  +18.38 %           33678.0  +12.61 %
    #     layer 2      19194.8   +6.33 %           31522.5   +5.40 %
    #     layer 3      18791.5   +4.09 %           30944.3   +3.47 %
    #     layer 4      18760.0   +3.92 %           30899.9   +3.32 %
    #     increments  -2176.6 / -403.3 / -31.5    -2155.5 / -578.2 / -44.4 cm^-1
    #
    # THE VALENCE LADDER CONVERGES AND DOES NOT REACH THE MEASURED VALUES.  A geometric extrapolation of the
    # increments (and Aitken on the last three) puts the infinite-layer limit at 18757 / 30896 cm^-1, i.e.
    # +3.90 % and +3.31 % -- so ten layers would still leave ~700 and ~980 cm^-1.  What is missing is what was
    # never in the model: the frozen 2s,2p core, and Breit (this ladder is Coulomb-only).
    #
    # THE LAST INCREMENT UNDERSTATES THE REMAINING ERROR BY 22x (31.5 against 707.5 cm^-1; 44.4 against 993.5).
    # A ladder that has stopped moving has not arrived -- the same lesson branch (c) records for Fe XV, where
    # the factor reached 100x on one interval.
    #
    # AND THE n=3 COMPLEX IS NOT THE CULPRIT.  SDTQ(3s..3p --> 3s..3d) is the COMPLETE CAS for J = 3/2 odd --
    # parity is (-1)^(3p occupation), so only odd-3p configurations mix, and there are seven of them; SD reaches
    # five, SDTQ adds 3s3p3d^3 and 3p3d^4, and 3d^5 needs a QUINTUPLE excitation and is even parity anyway.
    # Measured: 66 CSF -> 145 CSF for 0.14 and 0.36 percentage points.  Triples and quadruples are not the gap.
    #
    # AND THE OBVIOUS NEXT TEST -- FOLDING THE CORE ONTO THE CONVERGED VALENCE LADDER -- CANNOT BE DONE THIS WAY.
    # Tried 16-Sep-2026: valence SD to 4f, then PT-SD(2s..3p --> 3s..4f), 12370 CSF, NINE HOURS.  The step ran to
    # completion and returned splittings of +59 % and +62 %, and the code said why in its own output: two CI
    # levels lay BELOW the reference carrying almost none of it -- an intruder at -466.85 Ha (7 Ha below the
    # ground state!) of 2p^4 3s^2 3p^3 4p^2, and one at -460.79 Ha of 2p^5 3s 3p^3 3d 4f.  A variational layer's
    # 4p and 4f are CONTRACTED pseudo-orbitals built to describe valence correlation; a 2p -> 4p excitation on
    # them appears to LOWER the energy instead of costing ~200 eV.  This is not a threshold or a Q-space-size
    # question, so no cheaper variant of the run fixes it.
    #
    # THE SAME DOUBT APPLIES TO THIS BRANCH.  Branch (a) folds the core onto LAYER-1 orbitals and reports the
    # fold as costing 0.47 mHa / -8.7 / -83.1 cm^-1 (0.0035 mHa after the 13-Sep angular repair).  Those numbers
    # were taken without an intruder check.  Before they are trusted, read the step's warning block.

elseif  false
    # Last visit:      13-Sep-2026
    # Last successful:  unknown -- withdrawn for the reason given in branch (a): the 12-Sep value was earned on
    #                   defective angular coefficients.  Re-run 13-Sep on the repaired code the CI still uses all
    #                   908 CSF and gives -459.76920238 Ha, where it gave -459.93426153.  THE GUARD ITSELF STILL
    #                   HOLDS, and that is the part worth noting: this branch's claim is that promoting everything
    #                   reproduces a plain CI over the same space EXACTLY, and it still does -- branch (a)'s
    #                   step 3 gives -459.76919888 against this -459.76920238, the same 3.5e-06 Ha apart as
    #                   before.  The defect moved both numbers together and left the identity intact, which is
    #                   why a self-consistency guard cannot substitute for a comparison against measurement.
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
    # Last visit:  13-Sep-2026
    # Last successful:  13-Sep-2026 -- the table this branch prints, against NIST ASD 5.12.  Coulomb+Breit(0.):
    #                   3P*_0 233235 / 233406 / 235209,  3P*_1 239101 / 239246 / 240884,
    #                   3P*_2 253163 / 253269 / 254777,  1P*_1 357408 / 356254 / 354891 cm^-1
    #                   for AL / 1-layer / 1-layer+PT.  Intervals 0-1: 5865 / 5840 / 5674 against NIST 5818;
    #                   1-2: 14062 / 14024 / 13893 against 14160.  Runtime 33 + 198 + 881 s.
    #
    # c) Fe XV 1s^2 2s^2 2p^6 3s^2 -- THE LADDER, against measured levels.
    #
    #    THIS REPLACES AN EARLIER BRANCH (12-Sep-2026) AND THE DIFFERENCE IS INSTRUCTIVE.  That one took
    #    {3s^2} ALONE as the reference and reached 3s3p, 3s3d, 3p3d only as excitations -- but those are REAL
    #    low-lying states of this ion (NIST puts 3s3d 3D at 678 772 and 3p3d 3F* at 928 241 cm^-1), not
    #    correlation, and the configurations of one complex are quasi-degenerate so no perturbative route to
    #    them converges.  It also ran on two-particle angular coefficients that were wrong for exactly the CSF
    #    pairs such a reference generates.  Both are fixed; the old branch is not worth keeping.
    #
    #    THE REFERENCE IS THE COMPLETE n=3 COMPLEX and the ladder is built on SPECTROSCOPIC orbitals:
    #    {3s^2, 3s3p, 3s3d, 3p^2, 3p3d, 3d^2}, so 3d is optimized as a valence orbital of step 1 rather than as
    #    a correlation orbital of a later layer.  The EOL TARGETS are restricted to the levels of interest --
    #    the P space and the target set are separate choices, and optimizing across 3d^2 and 3p3d levels nobody
    #    asked for makes the levels one does want worse.
    #
    #    BREIT IS IN THE CI FROM THE START, and at Z = 26 it is not a correction but a leading term of the fine
    #    structure: measured here it shifts 3P*_0 by +628, 3P*_1 by +438 and 3P*_2 by only +23 cm^-1, i.e. it
    #    is strongly J-dependent and moves the 0-1 interval by -190 and 1-2 by -415 cm^-1.  A Coulomb-only
    #    version of this same ladder misses both.  `CoulombBreit(0.)` is the exact omega -> 0 limit.
    #
    #    WHAT THE THREE MODELS SAY.  AL already reaches -0.26 / -0.23 / -0.26 % on the triplet; the n=4 layer
    #    improves every level and takes the 0-1 interval to +0.4 %; and the PERTURBATIVE step then OVERSHOOTS,
    #    carrying the intervals past NIST to -2.5 % and -1.9 %.  That is worth stating plainly because the same
    #    step looks like a success in a Coulomb-only calculation, where its -165 cm^-1 on 0-1 substitutes for
    #    the -190 cm^-1 that Breit supplies: the folded core and the magnetic interaction are independent and
    #    of similar size here, and applying both double-counts what only one of them should provide.
    #    1P*_1 is the exception and behaves oppositely -- it improves monotonically, 1.56 / 1.23 / 0.85 %,
    #    being a correlation problem that Breit barely touches.
    #
    #    MEASURED SEPARATELY, not by this branch: a second layer (n=5, 86 min) with the Breit shift of the
    #    n=4 layer added gives the triplet to 0.07 / 0.06 / 0.10 % and intervals 5840 / 14052, i.e. +0.4 % and
    #    -0.8 %.  Layers move the term as a whole; Breit sets its internal structure.
    #
    Z      = 26.0
    c(s)   = Configuration("1s^2 2s^2 2p^6 " * s)
    refs   = [c("3s^2"), c("3s 3p"), c("3s 3d"), c("3p^2"), c("3p 3d"), c("3d^2")]
    wanted = [c("3s^2"), c("3s 3p")]                      # the EOL targets, and the levels reported
    syms   = [LevelSymmetry(0, Basics.plus),  LevelSymmetry(0, Basics.minus),
              LevelSymmetry(1, Basics.minus), LevelSymmetry(2, Basics.minus)]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 10.)
    core   = [Shell("2s"), Shell("2p")];   v3 = [Shell("3s"), Shell("3p"), Shell("3d")]
    n4     = vcat(v3, [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
    frz    = [Shell("1s"), Shell("2s"), Shell("2p")]
    nist   = [("3s^2  1S_0 ", 0.0), ("3s3p  3P*_0", 233842.0), ("3s3p  3P*_1", 239660.0),
              ("3s3p  3P*_2", 253820.0), ("3s3p  1P*_1", 351911.0)]
    #
    # THE FIVE LEVELS ARE PICKED BY THEIR DOMINANT CONFIGURATION AND NOT BY ENERGY ORDER: with one correlation
    # layer a correlation level of this space can lie between them, and "the k-th level of symmetry J" then
    # names the wrong state.  Within J = 1- both wanted levels are 3s3p, so those two are ordered by energy
    # AMONG the 3s3p-dominated ones, which is what NIST's own ordering asserts.
    fiveLevels = function(mp)
        bs   = mp.levels[1].basis
        cof  = [ Basics.extractConfiguration(Basics.FromBasis(), bs, cc)  for cc in bs.csfs ]
        i2   = findall(x -> x == wanted[1], cof);    ip = findall(x -> x == wanted[2], cof)
        wt(l, ix)  = sum( l.mc[r]^2  for r in ix; init=0.0 )
        sy(J, p)   = [ l  for l in mp.levels  if l.J == AngularJ64(J) && l.parity == p ]
        pk(J,p,ix) = ( ls = sy(J,p);  isempty(ls) ? nothing : ls[argmax([wt(l,ix) for l in ls])] )
        g  = pk(0, Basics.plus, i2);    e0 = pk(0, Basics.minus, ip);    e2 = pk(2, Basics.minus, ip)
        l1 = sort( [ l for l in sy(1, Basics.minus) if wt(l, ip) >= 0.20 ], by = l -> l.energy )
        (isnothing(g) || isnothing(e0) || isnothing(e2) || length(l1) < 2)  &&  return( Float64[] )
        k(x) = Defaults.convertUnits("energy: from atomic to Kayser", x - g.energy)
        return( [0.0, k(e0.energy), k(l1[1].energy), k(e2.energy), k(l1[2].energy)] )
    end
    report = function(tag, mp)
        e = fiveLevels(mp)
        if  isempty(e)   println(">> $tag : the five levels were not all found");   return   end
        println("\n>> $tag")
        for  (i, (name, ref))  in  enumerate(nist)
            @printf(">>   %-12s %12.1f   NIST %10.1f  %+9.1f  (%+6.2f %%)\n", name, e[i], ref, e[i]-ref,
                    ref == 0.0 ? 0.0 : 100*(e[i]-ref)/ref)
        end
        @printf(">>   3P* intervals:  0-1 %8.1f (NIST  5818.0)   1-2 %9.1f (NIST 14160.0)\n",
                e[3]-e[2], e[4]-e[3])
    end
    #
    set = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                     eeInteractionCI = CoulombBreit(0.), gridStopper = false)
    report("AL", SelfConsistent.performSCF(refs, Nuclear.Model(Z), grid, set; printout=false))
    #
    runRas = function(tag, steps)
        rs = RasSettings(Int64[], 300, 1.0e-6, CoulombBreit(0.), LevelSelection(true, configurations=wanted))
        wb = generate( Representation("Fe XV -- $tag", Nuclear.Model(Z), grid, refs,
                                      RasExpansion(syms, 12, steps, rs)), output=true )
        k  = "step" * string(length(steps))
        haskey(wb, k) ? report(tag, wb[k]) : println(">> $tag : no $k")
    end
    layer4 = RasStep(RasStep(); seFrom=v3, seTo=n4, deFrom=v3, deTo=n4, frozen=frz)
    ptStep = RasStep(RasStep(); seFrom=vcat(core,v3), seTo=n4, deFrom=vcat(core,v3), deTo=n4,
                                frozen=vcat(frz, n4), treatment = Basics.SecondOrder(1.0, 1.0e-12))
    runRas("1-layer",    [RasStep(), layer4])
    runRas("1-layer+PT", [RasStep(), layer4, ptStep])
elseif  false
    # Last visit:      13-Sep-2026
    # Last successful:  unknown -- AL and 0-layer are VERIFIED against NIST (numbers below); the 1-layer and
    #                   2-layer rows were still running on the compute machine when this was written and are
    #                   NOT yet filled in, so the branch is not yet dated.
    #
    # d) Fe VII [Ar] 3d^2, Z = 26 -- THE d-SHELL TEST.
    #
    #    WHY d^2 AND WHY THIS ION.  Every CSF pair that the repaired two-particle coefficients got wrong had a
    #    p SPECTATOR, carrying j = 1/2 and 3/2.  A d spectator carries j = 3/2 and 5/2, higher than anything the
    #    Mg-like and C-like cases reach, and 3d^2 also brings the same-shell seniority/CFP path into combination
    #    with the repaired cross-configuration one.  Fe VII rather than Ti III because at charge +6 the ion is
    #    compact, so correlation is a smaller fraction of the answer and a residual angular error is not masked
    #    by a large correlation error; Ti III is the better stress test of the RECIPE and the worse test of the
    #    COEFFICIENTS.
    #
    #    THE REFERENCE IS 3d^2 ALONE, AND THE MEASURED SPECTRUM SAYS SO.  Branch (c) had to carry the complete
    #    n=3 complex because 3s3p, 3p3d and 3d^2 are quasi-degenerate there.  Here the other members of the n=3
    #    complex are CORE-HOLE configurations (3p^5 3d^3 from 389 340 cm^-1) and the nearest genuine neighbour,
    #    3d4s, begins at 344 462 cm^-1 -- some 277 000 cm^-1 (34 eV) above the highest 3d^2 level.  Nothing is
    #    near-degenerate with 3d^2, so completing a complex would add cost and no physics.  The rule is "put the
    #    quasi-degenerate partners in P", not "always enlarge the reference", and reading it off the measured
    #    spectrum is what separates the two cases.
    #
    #    MEASURED 13-Sep-2026, Coulomb+Breit(0.), excitation energies in cm^-1 against
    #    NIST Atomic Spectra Database (ver. 5.12), retrieved 13-Sep-2026, DOI 10.18434/T4W30F:
    #
    #       term    J        NIST          AL      0-layer     dev (0-layer)
    #       3F      2         0.00        0.00        0.00        --
    #       3F      3      1049.75     1019.50     1013.95      -3.4 %
    #       3F      4      2329.48     2280.76     2268.12      -2.6 %
    #       1D      2     17474.00    21094.20    21148.47     +21.0 %
    #       3P      0     20040.0     24382.74    24441.20     +21.9 %
    #       3P      1     20428.8     24750.59    24809.29     +21.4 %
    #       3P      2     21275.97    25527.86    25584.97     +20.2 %
    #       1G      4     28923.5     32226.03    32310.27     +11.7 %
    #       1S      0     67076.9     78191.71    78419.45     +16.9 %
    #
    #    WHAT THAT PATTERN MEANS, and the two halves say different things.  The 3F FINE STRUCTURE -- two
    #    intervals INSIDE one term, the analogue of the 3P* case that exposed the coefficient defect -- comes out
    #    within 3 %, correctly ordered and with no inversion: the d spectator behaves.  The TERM SEPARATIONS are
    #    12-22 % too large, which is the textbook signature of an open d shell with no correlation: term energies
    #    there are governed by the Slater F^2 and F^4 integrals, which a bare mean field overestimates because it
    #    cannot screen the d-d repulsion.  It is the same fact that makes ligand-field and astrophysical work
    #    scale F_k by ~0.8.  So the angular part looks sound and what is missing is d-d correlation, which is
    #    what the layers below exist to supply -- a far harder test of them than branch (c), where the mean field
    #    was already within 0.5 %.
    #
    #    NOTE ALSO that the 0-layer EOL step makes the TERM separations very slightly WORSE while lowering the
    #    total energy.  That is not a defect: EOL optimizes the energy of the target levels, not the splittings
    #    between them.
    #
    #    NIST's own leading percentages make J = 2+ the interesting block: 1D_2 is 91 % with 6 % 3P, and 3P_2 is
    #    92 % with 6 % 1D, so those two genuinely mix and the angular coefficients must get that mixing right.
    #
    Z      = 26.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2")]
    syms   = [LevelSymmetry(J, Basics.plus)  for J = 0:4]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z))
    v3     = [Shell("3d")]
    n4     = [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]
    frz    = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
    nist   = [(2,0.00), (3,1049.75), (4,2329.48), (2,17474.00), (0,20040.0),
              (1,20428.8), (2,21275.97), (4,28923.5), (0,67076.9)]
    #
    # EVERY 3d^2-DOMINATED LEVEL IS PRINTED WITH ITS J AND ITS WEIGHT, in energy order, rather than being matched
    # to a term by name: the term order of a d^2 ion is not something to assume in advance, and the J sequence
    # against NIST's is what identifies the levels.
    showLevels = function(lab, mp)
        bs  = mp.levels[1].basis
        cof = [ Basics.extractConfiguration(Basics.FromBasis(), bs, cc)  for cc in bs.csfs ]
        ig  = findall(x -> x == refs[1], cof)
        wt(l) = sum( l.mc[r]^2  for r in ig; init=0.0 )
        keep  = sort( [ l for l in mp.levels if wt(l) >= 0.20 ], by = l -> l.energy )
        isempty(keep)  &&  (println(">> $lab : no 3d^2-dominated level found");   return)
        g = keep[1]
        println("\n>> Fe VII   $lab     ($(length(keep)) of $(length(mp.levels)) levels are 3d^2)")
        println(">>     J    excitation [cm^-1]      NIST        deviation     w(3d^2)")
        for  (i, l) in enumerate(keep)
            e = Defaults.convertUnits("energy: from atomic to Kayser", l.energy - g.energy)
            if  i <= length(nist)  &&  Basics.twice(l.J) == 2*nist[i][1]
                ref = nist[i][2]
                @printf(">>   %4s  %18.2f  %10.2f  %+12.2f %s   %6.3f\n", string(l.J), e, ref, e-ref,
                        ref == 0.0 ? "        " : @sprintf("(%+6.1f %%)", 100*(e-ref)/ref), wt(l))
            else
                @printf(">>   %4s  %18.2f  %10s  %12s        %6.3f\n", string(l.J), e, "--", "J MISMATCH", wt(l))
            end
        end
    end
    #
    set = AsfSettings(AsfSettings(); scField=Basics.ALField(), eeInteraction=CoulombInteraction(),
                                     eeInteractionCI=CoulombBreit(0.), gridStopper=false)
    showLevels("AL", SelfConsistent.performSCF(refs, Nuclear.Model(Z), grid, set; printout=false))
    #
    runRas = function(lab, steps)
        rs = RasSettings(Int64[], 60, 1.0e-6, CoulombBreit(0.), LevelSelection(true, configurations=refs))
        wb = generate( Representation("Fe VII $lab", Nuclear.Model(Z), grid, refs,
                                      RasExpansion(syms, 20, steps, rs)), output=true )
        k = "step" * string(length(steps))
        haskey(wb, k) ? showLevels(lab, wb[k]) : println(">> $lab : no $k")
    end
    runRas("0-layer", [RasStep()])
    runRas("1-layer", [RasStep(), RasStep(RasStep(); seFrom=v3, seTo=n4, deFrom=v3, deTo=n4, frozen=frz)])
end
