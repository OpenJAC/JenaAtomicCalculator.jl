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
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(),
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
    # Last visit:  07-Nov-2026
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
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Cl III -- the guard: promote everything, reproduce the exact CI",
                            Nuclear.Model(Z), grid, refs, RasExpansion([sym], 15, steps, rasSettings) )
    println("wa = $wa")
    wb = generate(wa, output=true)
    ls = sort(wb["step3"].levels, by = l -> l.energy)
    println(">> step 3 : $(length(ls)) levels, lowest = $(ls[1].energy) Ha  -- EXPECTED -459.93426153 Ha")

elseif  true
    # Last visit:  08-Nov-2026 (TIMED OUT at 500 s)
    # Last visit:  07-Nov-2026   (TIMED OUT at 500 s)
    # Last visit:  20-Oct-2026 (TIMED OUT at 500 s)
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
        rs = RasSettings(Int64[], Basics.RotationRoute(300), 1.0e-6, CoulombBreit(0.), LevelSelection(true, configurations=wanted))
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
    # Last visit:      16-Oct-2026
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
        rs = RasSettings(Int64[], Basics.RotationRoute(60), 1.0e-6, CoulombBreit(0.), LevelSelection(true, configurations=refs))
        wb = generate( Representation("Fe VII $lab", Nuclear.Model(Z), grid, refs,
                                      RasExpansion(syms, 20, steps, rs)), output=true )
        k = "step" * string(length(steps))
        haskey(wb, k) ? showLevels(lab, wb[k]) : println(">> $lab : no $k")
    end
    runRas("0-layer", [RasStep()])
    runRas("1-layer", [RasStep(), RasStep(RasStep(); seFrom=v3, seTo=n4, deFrom=v3, deTo=n4, frozen=frz)])
elseif  false
    # Last visit:  07-Nov-2026
    # Last visit:  02-Nov-2026
    #
    # e) THE TWO Q-SPACE PERTURBATION ROUTES, ON BRANCH (a)'s OWN Cl III CASE (competition item 11).
    #    Basics.AbstractQSpacePerturbation gives secondOrderBlock a second way to fold a Q-space
    #    configuration: DiagonalCSFPerturbation (every CSF isolated, the default and what this file has
    #    always used) or BlockDiagonalizedPerturbation (a folded configuration's own CSFs diagonalized
    #    together first, then folded using ITS eigenstates). This branch runs EXACTLY branch (a)'s own
    #    three-step ladder -- reference, valence 3s,3p->3d, core 2s,2p+3s,3p fold -- through BOTH routes
    #    and compares them, rather than inventing a new case. tools/probe-secondOrderQ.jl already measured
    #    that an effective Hamiltonian is safe in exactly this (core, weakly-coupled) regime; this branch
    #    is the production-path confirmation of that, not a new finding.
    Z      = 17.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]
    sym    = LevelSymmetry(AngularJ64(3//2), Basics.minus)
    val    = [Shell("3s"), Shell("3p")];    core = [Shell("2s"), Shell("2p")]
    to     = [Shell("3s"), Shell("3p"), Shell("3d")]
    frozen1 = [Shell("1s"), Shell("2s"), Shell("2p")]
    frozen2 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p"), Shell("3d")]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 20.)
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    results = Tuple{String,Float64,Float64}[]
    for  (lab, pert)  in  [ ("DiagonalCSFPerturbation (default)", Basics.DiagonalCSFPerturbation()),
                            ("BlockDiagonalizedPerturbation",     Basics.BlockDiagonalizedPerturbation()) ]
        steps = [ RasStep(),
                  RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to, frozen=frozen1),
                  RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                     frozen=frozen2,
                                     treatment = Basics.SecondOrder(1.0e-5, 1.0e-12, 1.0e-2, pert)) ]
        wa = Representation("Cl III -- $lab", Nuclear.Model(Z), grid, refs, RasExpansion([sym], 15, steps, rasSettings) )
        wb = generate(wa, output=true)
        e2 = sort(wb["step2"].levels, by = l -> l.energy)[1].energy
        e3 = sort(wb["step3"].levels, by = l -> l.energy)[1].energy
        push!(results, (lab, e2, e3))
        println("\n>> $lab :  step 2 (valence) = $e2 Ha,  step 3 (core fold) = $e3 Ha,  increment = " *
                "$((e3-e2)*1000) mHa")
    end
    println("\n  The two routes agree on step 2 (unaffected by the new code) and on step 3 to:")
    println("  " * string(abs(results[1][3] - results[2][3])) * " Ha " *
            "(tools/probe-secondOrderQ.jl: the two forms are indistinguishable in this weakly-coupled, " *
            "core regime -- +-0.2 cm^-1 there, which is the same order as this difference).")
    abs(results[1][2] - results[2][2]) < 1.0e-12  ||
        error("example-Ar.jl branch e: the two perturbation routes disagree on step 2, which neither " *
              "route touches -- something is leaking between them.")
    #
elseif  false
    # Last visit:  02-Nov-2026
    #
    # f) THE DUAL-ROLE SHELL GUARD (competition item 11), built small rather than reproducing the nine
    #    hours branch (a)'s own header records. A core-valence step (seFrom/deFrom naming a shell this
    #    step itself freezes) is not allowed to send its excitation into a shell an EARLIER VARIATIONAL
    #    layer already built a correlation orbital for -- exactly the pathology that gave branch (a)'s
    #    16-Sep-2026 failure two intruder levels, one 7 Ha below the ground state, from a 2p -> 4p
    #    "excitation" that should cost ~200 eV instead looking free. On Be, cheap: layer 2 introduces 2p
    #    as a valence correlation shell; layer 3 then (wrongly) tries to reuse 2p as a core-excitation
    #    target from 1s,2s. Check `jac-warn.report` after running this branch for the physics-worded
    #    message this is supposed to produce.
    name        = "Beryllium -- dual-role shell guard (deliberately reproduces the failure pattern, cheaply)"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1]) )
    steps = [ RasStep(),
              RasStep(RasStep(); seFrom=[Shell("2s")], seTo=[Shell("2p")], deFrom=[Shell("2s")], deTo=[Shell("2p")],
                                 frozen=[Shell("1s")]),
              RasStep(RasStep(); seFrom=[Shell("1s"), Shell("2s")], seTo=[Shell("2p")],
                                 deFrom=[Shell("1s"), Shell("2s")], deTo=[Shell("2p")],
                                 frozen=[Shell("1s"), Shell("2s"), Shell("2p")],
                                 treatment=Basics.SecondOrder(1.0e-5, 1.0e-12)) ]
    wa = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs,
                        RasExpansion([LevelSymmetry(0, Basics.plus)], 4, steps, rasSettings) )
    wb = generate(wa, output=true)
    println("\n>> Check jac-warn.report now: it must name step 3, the shells 1s, 2s, and the shell 2p, in")
    println(">> the physics-worded message (\"... which is already a valence correlation orbital in this")
    println(">> ladder ... add the core-valence excitation into an unused shell instead\").")
    #
elseif  false
    # Last visit:  04-Nov-2026
    #
    # g) Cl III 1s^2 2s^2 2p^6 3s^2 3p^3 -- THE FULL GROUND TERM, valence correlation plus the NEW
    #    shape-crossing core-valence correction (AtomicState.RasStep.shapedTo).
    #
    #    THREE SYMMETRIES, ONE RUN: J = 1/2, 3/2, 5/2 odd together cover all five fine-structure levels of
    #    the three LS terms of p^3 (4S, 2D, 2P).  Three steps, same pattern as branch (a):
    #      step 1   the reference alone
    #      step 2   VALENCE correlation 3s,3p -> 3d, variational.  3d is optimized HERE, as an ordinary
    #               SCFShape correlation orbital.
    #      step 3   shapedTo = [3d]: the core excitations 2s,2p -> 3d get their OWN, independently
    #               contracted CoreValenceCorrelationShape realization of 3d, built and refined on its own
    #               small Q space (BasicsAZ.retagShapedSubshells / generateOrbitalsForShapedBasis), then
    #               coupled back to step 2's own P via Hamiltonian.secondOrderShapeCrossingCorrection --
    #               BiOrthogonal.computeTransformation(ShapeCrossingJunction(), ...), summed over P's CI
    #               mixing coefficients and Q's genuinely-new CSFs only.  Unlike branch (a)'s step 3, which
    #               reuses step 2's 3d UNCHANGED and is correspondingly cheap but shape-blind, this is the
    #               case the dual-role-shell guard WOULD refuse outside this sanctioned path -- see
    #               BasicsAZ.runRasStep!'s own exemption for shapedTo.
    #
    #    MEASURED, 04-Nov-2026, Coulomb, rbox = 20 a.u., RotationRoute(24), ~5.5 min total for all three
    #    symmetries in the ONE combined RasExpansion the code below actually runs:
    #
    #      #  term J     E(step2), no MBPT    MBPT          E(step2+MBPT)     cm^-1 no-MBPT  cm^-1 +MBPT   NIST cm^-1
    #      1  4S  3/2   -459.77080847 Ha    -3.9346 mHa   -459.77474306 Ha        0.00           0.00          0.000
    #      2  2D  3/2   -459.67457008 Ha    -4.1435 mHa   -459.67871357 Ha    21121.88       21076.04      18052.460
    #      3  2D  5/2   -459.67423974 Ha    -4.1487 mHa   -459.67838841 Ha    21194.38       21147.40      18118.430
    #      4  2P  1/2   -459.61915232 Ha    -4.1627 mHa   -459.62331503 Ha    33284.68       33234.61      29813.920
    #      5  2P  3/2   -459.61877330 Ha    -4.1747 mHa   -459.62294802 Ha    33367.86       33315.16      29906.450
    #
    #    Term labels by energy order (unambiguous here: 4S lies far below 2D, which lies far below 2P, exactly
    #    as expected for p^3) and NIST ASD observed energies (physics.nist.gov, via the CHIANTI cl_3.elvlc table,
    #    which cites NIST ASD directly) are both NEW this visit.  THE EOL TARGET SET MATTERS: running J = 3/2
    #    alone gives step 2 energies and a level 4/5 ORDER DIFFERENT from running all three symmetries together,
    #    since EOL optimizes orbitals for a functional averaged over however many levels are named as targets --
    #    three symmetries at once is a five-level average, one symmetry alone is a smaller, different average.
    #    The numbers above are what THIS branch, as written, actually produces; a single-symmetry run is not a
    #    cheaper way to get the same numbers.
    #
    #    THE CORRECTION IS OVERWHELMINGLY A COMMON SHIFT -- -3.93 to -4.17 mHa, a spread of only 0.24 mHa
    #    across every level of all three terms -- which is what reusing the SAME orbital for the SAME
    #    channel on every level should give.  What survives into the SPLITTINGS is much smaller and, here,
    #    NARROWS both fine-structure intervals slightly: level 3 - level 2 by -1.1 cm^-1 (72.50 -> 71.36,
    #    against NIST's 65.97), level 5 - level 4 by -2.6 cm^-1 (83.18 -> 80.55, against NIST's 92.53).  Against
    #    run2.jl's own "clean but shell-starved" capture of -0.2478 mHa (a fresh, never-reused n5 shell, branch
    #    (a)'s sibling scratch script), reshaping the ALREADY-CONTRACTED 3d instead of using a poorly-matched
    #    fresh one captures roughly SIXTEEN TIMES as much core-valence correlation energy on the common shift
    #    alone.
    #
    #    AGAINST NIST, THE TERM SEPARATIONS THEMSELVES ARE STILL TOO WIDE BY ~3000-3500 cm^-1 (about 15-17 %) --
    #    4S-2D is 3069-3076 cm^-1 too high, 4S-2P is 3409-3471 cm^-1 too high, and MBPT closes only a small part
    #    of that (45-52 cm^-1, under 2 %).  This is NOT a defect: the gap is the one-layer valence correlation
    #    (3s,3p -> 3d only) being too small on its own to fix the term separations, which is exactly why the
    #    user's own framing of this round asked for "1-2 correlation layers" rather than a converged result --
    #    the shape-crossing core-valence piece corrects a DIFFERENT, smaller effect (the fine-structure splitting
    #    within a term) than the missing valence layers would.  A converged comparison against NIST needs more
    #    layers and is not what this branch sets out to show.
    #
    Z      = 17.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]
    syms   = [LevelSymmetry(AngularJ64(1//2), Basics.minus), LevelSymmetry(AngularJ64(3//2), Basics.minus),
              LevelSymmetry(AngularJ64(5//2), Basics.minus)]
    val    = [Shell("3s"), Shell("3p")];    core = [Shell("2s"), Shell("2p")]
    to     = [Shell("3s"), Shell("3p"), Shell("3d")]
    frozen1 = [Shell("1s"), Shell("2s"), Shell("2p")]
    frozen2 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to, frozen=frozen1),
               RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                  shapedTo=[Shell("3d")], frozen=frozen2,
                                  treatment = Basics.SecondOrder(1.0e-5, 1.0e-12)) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 20.)
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Cl III -- full ground term, valence + shape-crossing core-valence MBPT",
                            Nuclear.Model(Z), grid, refs, RasExpansion(syms, 15, steps, rasSettings) )
    wb     = generate(wa, output=true)
    refKeys = Set( Hamiltonian.configurationKey(c)  for c in refs )
    refLevels = function(mp)
        ls = Level[]
        for  lv in mp.levels
            bs = lv.basis
            w  = sum( lv.mc[r]^2  for r = 1:length(bs.csfs)
                      if  Hamiltonian.configurationKey(
                              Basics.extractConfiguration(Basics.FromBasis(), bs, bs.csfs[r])) in refKeys )
            if  w >= 0.5    push!(ls, lv)    end
        end
        return( sort(ls, by = l -> l.energy) )
    end
    println("\n>> step 2 (valence correlation only), reference-dominated levels, lowest first:")
    for  (i, lv)  in  enumerate(refLevels(wb["step2"]))
        println(">>   level $i   J = $(lv.J)   E = $(lv.energy) Ha")
    end
    println(">> step 3's own printout above carries the MBPT correction for each of these, in the same ")
    println(">> energy order within each J -- see 'ShapeCrossingPartition' and the 'deltaE' lines.")
    #
elseif  false
    # Last visit:  04-Nov-2026
    #
    # h) Cl III -- DEEPENING THE VALENCE LADDER ALONE, n3 -> n3+n4, NO shape-crossing MBPT step.  Answers a
    #    narrower, separate question from branch (g): how much of the ~3000-3500 cm^-1 NIST gap left after ONE
    #    valence layer (3s,3p -> 3d) closes once a SECOND valence layer (-> 4s,4p,4d,4f) is added -- "more
    #    valence layers, a separate well-understood scope question", independent of how close the result
    #    lands to experiment.
    #
    #    THREE STEPS: step 1 reference, step 2 the SAME n3 layer as branch (g) (3s,3p -> 3d, variational),
    #    step 3 a SECOND layer, 3s,3p -> 3d,4s,4p,4d,4f, freezing the n3 layer's own 3d (so only the NEW n4
    #    shells vary here).  No shapedTo step is attached: STACKING THE SHAPE-CROSSING MBPT CORRECTION ON TOP
    #    OF A DEEPENED LADDER HIT A NEW, GENUINE LIMITATION, measured this visit, NOT fixed -- once the PRIOR
    #    step carries more subshells than the MBPT step's own 2s,2p -> 3d excitation space touches (9 here,
    #    against this ladder's 16), Hamiltonian.secondOrderShapeCrossingCorrection's own safety assertion
    #    (P and Q must agree (n,kappa) for (n,kappa), position for position) correctly REFUSES rather than
    #    silently relabelling the wrong orbitals -- branch (g) never hit this because there P WAS step 2's own
    #    9-subshell basis.  Composing the shape-crossing step with a deeper ladder is an open question for
    #    `src/` (priority item 5), not attempted here.
    #
    #    MEASURED, 04-Nov-2026, Coulomb, rbox = 28 a.u., RotationRoute(32), combined 3 symmetries, 3m18s total:
    #    the n4 layer's own doubling test: -52.1556 mHa at 32 iterations against -48.3675 mHa at 16, a 7.3 %
    #    change -- CONVERGED in its budget, unlike an earlier n4-alone probe at RotationRoute(24) (16.1 %,
    #    not yet converged).
    #
    #      #  term J     E(n3) [Ha]          E(n3+n4) [Ha]       NIST cm^-1   gap(n3) cm^-1  gap(n3+n4) cm^-1  closes
    #      1  4S  3/2   -459.77079156 Ha    -459.80175423 Ha         0.000          0.00           0.00         --
    #      2  2D  3/2   -459.67457044 Ha    -459.71494778 Ha     18052.460       3065.63         999.35       67.4 %
    #      3  2D  5/2   -459.67423456 Ha    -459.71460272 Ha     18118.430       3073.38        1009.12       67.2 %
    #      4  2P  1/2   -459.61915520 Ha    -459.65991407 Ha     29813.920       3466.41        1316.40       62.0 %
    #      5  2P  3/2   -459.61876353 Ha    -459.65949487 Ha     29906.450       3459.85        1315.87       62.0 %
    #
    #    ONE LAYER OF n4 CLOSES ROUGHLY TWO THIRDS OF THE GAP left after n3 alone -- a clean, LAYER-BY-LAYER
    #    convergence (Rule 21: compare increments, not just endpoints), the 2D and 2P pairs each closing by
    #    near-identical amounts, which is itself evidence the residual ~1000-1350 cm^-1 is further CORRELATION
    #    rather than noise.  Whether a THIRD layer (n5) would close most of the rest, or hit diminishing
    #    returns, is left open: n4's own CSF count already grew 12.5x over n3 alone (172 -> 2151, J=3/2 block),
    #    and extrapolating that growth rate to n5 risks a run that would not finish, or not converge, inside a
    #    reasonable time budget -- STOPPED HERE deliberately rather than attempted and left unconverged.
    #
    Z      = 17.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]
    syms   = [LevelSymmetry(AngularJ64(1//2), Basics.minus), LevelSymmetry(AngularJ64(3//2), Basics.minus),
              LevelSymmetry(AngularJ64(5//2), Basics.minus)]
    val    = [Shell("3s"), Shell("3p")]
    to3    = [Shell("3s"), Shell("3p"), Shell("3d")]
    to4    = [Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]
    frozen1 = [Shell("1s"), Shell("2s"), Shell("2p")]
    frozen2 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
    frozen3 = vcat(frozen2, [Shell("3d")])
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to3, deFrom=val, deTo=to3, frozen=frozen1),
               RasStep(RasStep(); seFrom=val, seTo=to4, deFrom=val, deTo=to4, frozen=frozen3) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 28.)
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(32), 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Cl III -- two valence layers, n3+n4, no MBPT",
                            Nuclear.Model(Z), grid, refs, RasExpansion(syms, 15, steps, rasSettings) )
    wb     = generate(wa, output=true)
    refKeys = Set( Hamiltonian.configurationKey(c)  for c in refs )
    refLevels = function(mp)
        ls = Level[]
        for  lv in mp.levels
            bs = lv.basis
            w  = sum( lv.mc[r]^2  for r = 1:length(bs.csfs)
                      if  Hamiltonian.configurationKey(
                              Basics.extractConfiguration(Basics.FromBasis(), bs, bs.csfs[r])) in refKeys )
            if  w >= 0.5    push!(ls, lv)    end
        end
        return( sort(ls, by = l -> l.energy) )
    end
    println("\n>> step 2 (n3 valence only), reference-dominated levels, lowest first:")
    for  (i, lv)  in  enumerate(refLevels(wb["step2"]))
        println(">>   level $i   J = $(lv.J)   E = $(lv.energy) Ha")
    end
    println("\n>> step 3 (n3+n4 valence), reference-dominated levels, lowest first:")
    for  (i, lv)  in  enumerate(refLevels(wb["step3"]))
        println(">>   level $i   J = $(lv.J)   E = $(lv.energy) Ha")
    end
    #

elseif  false
    # Last visit:  08-Nov-2026 (TIMED OUT at 500 s -- an overnight sweep's own cap; this branch's own
    #                            comment below already measures it at 980 s by hand, so this is expected,
    #                            not a regression)
    # Last visit:  04-Nov-2026
    #
    # i) Ti III 1s^2 2s^2 2p^6 3s^2 3p^6 3d^2 -- THE FULL GROUND TERM on a d^2 CONFIGURATION, same
    #    shape-crossing core-valence architecture as branch (g), tried on a DIFFERENT angular case.
    #
    #    FIVE SYMMETRIES, ONE RUN: J = 0,1,2,3,4 even together cover all nine fine-structure levels of
    #    the five LS terms of d^2 (3F, 1D, 3P, 1G, 1S).  Three steps, same pattern as branch (g), one
    #    shell up: 3d plays the role 3s,3p played for Cl III (the correlated VALENCE shell), 3s,3p play
    #    the role 2s,2p played (the CORE being correlated into the valence correlation orbital):
    #      step 1   the reference alone
    #      step 2   VALENCE correlation 3d -> 4d, variational.  4d is optimized HERE, an ordinary
    #               SCFShape correlation orbital.
    #      step 3   shapedTo = [4d]: the core excitations 3s,3p -> 4d get their OWN, independently
    #               contracted CoreValenceCorrelationShape realization of 4d, coupled back to step 2's
    #               own P via Hamiltonian.secondOrderShapeCrossingCorrection, exactly as in branch (g).
    #
    #    MEASURED, 04-Nov-2026, Coulomb, rbox = 24 a.u. (Rule 12: r_+ ~ 9.55 a.u. for 4d at Z_eff = 3,
    #    2.5x that is 23.9), RotationRoute(24), 980 s (16.3 min) for all five symmetries in the ONE
    #    combined RasExpansion the code below actually runs:
    #
    #      #  term J     E(step2), no MBPT    MBPT          E(step2+MBPT)     cm^-1 no-MBPT  cm^-1 +MBPT   NIST cm^-1
    #      1  3F  2     -852.18031725 Ha   -386.1516 mHa  -852.56646885 Ha        0.00           0.00          0.000
    #      2  3F  3     -852.17947578 Ha   -386.7165 mHa  -852.56619228 Ha      184.68          60.70        184.900
    #      3  3F  4     -852.17838493 Ha   -388.3614 mHa  -852.56674633 Ha      424.10         -60.90        420.400
    #      4  1D  2     -852.13107813 Ha   -431.4523 mHa  -852.56253043 Ha    10806.74         864.38       8473.500
    #      5  3P  0     -852.12108194 Ha   -414.9108 mHa  -852.53599274 Ha    13000.65        6688.73      10538.400
    #      6  3P  1     -852.12080630 Ha   -434.1391 mHa  -852.55494540 Ha    13061.15        2529.11      10603.600
    #      7  3P  2     -852.12024151 Ha   -434.0010 mHa  -852.55424251 Ha    13185.10        2683.37      10721.200
    #      8  1G  4     -852.10401746 Ha   -400.1674 mHa  -852.50418486 Ha    16745.87       13669.76      14397.600
    #      9  1S  0     -851.99825540 Ha   -537.3595 mHa  -852.53561490 Ha    39957.96        6771.66      32475.500
    #
    #    (NIST: Sugar & Corliss, "Energy Levels of Titanium, Ti I through Ti XXII", J. Phys. Chem. Ref.
    #    Data, the NIST compiled Ti III table.)
    #
    #    STEP 2 ALONE IS A GOOD, REPRODUCIBLE RESULT: term separations come out 16-28 % too wide against
    #    NIST -- the SAME ballpark and direction as branch (g)'s Cl III result (15-17 % too wide), one
    #    valence layer behaving consistently across a p^3 AND a d^2 case.
    #
    #    THE MBPT LAYER DOES NOT BEHAVE THE SAME WAY HERE, AND THE HONEST READING IS THAT IT IS NOT
    #    SPECTROSCOPICALLY USEFUL ON THIS CASE AS IT STANDS.  Unlike Cl III's -3.93 to -4.17 mHa common
    #    shift (spread 0.24 mHa), the Ti III correction ranges -386 to -537 mHa -- a 151 mHa SPREAD, i.e.
    #    ~33200 cm^-1, LARGER than the entire NIST term span (32476 cm^-1 to 1S0).  It does not merely
    #    fail to help: it DESTROYS the ordering step 2 already had right -- level 3 (3F4) drops BELOW
    #    the nominal ground level (-60.90 cm^-1), and 1S0 collapses from 39958 to 6772 cm^-1, landing
    #    inside the 3P cluster instead of far above it.
    #
    #    MOST LIKELY CAUSE, recorded for whoever next touches this (priority item 5): Q here holds 5371
    #    genuinely-new CSFs (printed by ShapeCrossingPartition) against Cl III's far smaller p^3-derived
    #    Q space, because a d^2 -> (3d,4d) manifold is intrinsically much denser.  With 9 target levels
    #    each weighed against that whole pool, some terms plausibly sit at SMALL-BUT-NOT-TINY denominators
    #    (hqq[s] - E_level) -- small enough to inflate the 1/Delta sum, but not caught by
    #    secondOrderShapeCrossingCorrection's existing `abs(dd) < 1.0e-10` floor, which only guards the
    #    literally-degenerate case.  This is NOT a box-size artifact (rbox matched Rule 12's own formula to
    #    within 0.1 a.u., and the run produced no box or convergence warnings at all).  A genuine repair --
    #    a denominator floor wide enough to exclude near-degenerate Q-CSFs, or restricting Q the way
    #    condition (a) for explicit second-order amplitudes already asks for elsewhere in JAC -- is
    #    separate, deliberate work, not attempted here.
    #
    Z      = 22.0
    refs   = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2")]
    syms   = [ LevelSymmetry(AngularJ64(j), Basics.plus)  for j = 0:4 ]
    val    = [Shell("3d")];    core = [Shell("3s"), Shell("3p")]
    to     = [Shell("3d"), Shell("4d")]
    frozen1 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
    frozen2 = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p"), Shell("3d")]
    steps  = [ RasStep(),
               RasStep(RasStep(); seFrom=val, seTo=to, deFrom=val, deTo=to, frozen=frozen1),
               RasStep(RasStep(); seFrom=vcat(core,val), seTo=to, deFrom=vcat(core,val), deTo=to,
                                  shapedTo=[Shell("4d")], frozen=frozen2,
                                  treatment = Basics.SecondOrder(1.0e-5, 1.0e-12)) ]
    grid   = Basics.recommendedGrid(refs, Nuclear.Model(Z); rbox = 24.)
    rasSettings = RasSettings(Int64[], Basics.RotationRoute(24), 1.0e-6, CoulombInteraction(),
                              LevelSelection(true, configurations=refs))
    wa     = Representation("Ti III -- full 3d^2 ground term, valence + shape-crossing core-valence MBPT",
                            Nuclear.Model(Z), grid, refs, RasExpansion(syms, 15, steps, rasSettings) )
    wb     = generate(wa, output=true)
    refKeys = Set( Hamiltonian.configurationKey(c)  for c in refs )
    refLevels = function(mp)
        ls = Level[]
        for  lv in mp.levels
            bs = lv.basis
            w  = sum( lv.mc[r]^2  for r = 1:length(bs.csfs)
                      if  Hamiltonian.configurationKey(
                              Basics.extractConfiguration(Basics.FromBasis(), bs, bs.csfs[r])) in refKeys )
            if  w >= 0.5    push!(ls, lv)    end
        end
        return( sort(ls, by = l -> l.energy) )
    end
    println("\n>> step 2 (valence correlation only), reference-dominated levels, lowest first:")
    for  (i, lv)  in  enumerate(refLevels(wb["step2"]))
        println(">>   level $i   J = $(lv.J)   E = $(lv.energy) Ha")
    end
    println(">> step 3's own printout above carries the MBPT correction for each of these, in the same ")
    println(">> energy order within each J -- see 'ShapeCrossingPartition' and the 'deltaE' lines.")
    #

end
