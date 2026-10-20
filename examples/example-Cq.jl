println("Cq) Multipole polarizabilities alpha^(Mp,K) at real and imaginary frequency, and the long-range " *
        "dispersion coefficients C_6, C_8, C_10 and C_3 built from them.")

## THE NOTATION USED THROUGHOUT THIS FILE IS alpha^(Mp,K):  ONE parenthesised superscript, the multipole first
## and the rank second, and no subscripts anywhere.  The literature's `alpha_2` means the rank-2 DIPOLE
## polarizability in Stark work and the SCALAR QUADRUPOLE one in dispersion work, which are different
## quantities;  here they are alpha^(E1,2) and alpha^(E2,0) and cannot be confused.  alpha^(E1,0) and
## alpha^(E1,2) are the alpha_0 and alpha_2 of Angel & Sandars.
##
## OF THE SIX BRANCHES BELOW, ONLY (c) CARRIES A `Last successful` DATE, and the rest say `unknown` on
## purpose (Rule 7): every branch's claim line quotes numbers that the branch reproduces, but the physical
## judgement that earns the date is the maintainer's and is not machine-writable. A blank date here is a
## statement that nobody has judged the physics yet, not that the branch failed.

if  false
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (2 of 2 numbers on the claim line, after the branch was converted to the
    #   standard perform route -- which is what says the conversion changed no physics)
    # Last successful:  unknown ... (the numbers below reproduce; the date is the maintainer's to set -- Rule 7)
    # Branch a: alpha^(Mp,K) FOR ANY MULTIPOLE AND ANY RANK A LEVEL CAN CARRY.  H 2p_3/2 with nd perturbers is
    #   the smallest case in which the rank-2 polarizability is NOT ZERO, so it is the one that actually tests
    #   the rank decomposition rather than just the scalar part.
    #
    #   CLAIM, 19-Oct-2026:   alpha^(E1,0) = +1.639094759924e+02 a.u.,   alpha^(E1,2) = -1.639099954070e+01 a.u.,
    #                         and alpha^(E2,0) = alpha^(E2,2) = 0 EXACTLY -- see the note on the E2 zeros below.
    #
    #   The three selection rules are carried BY the angular coefficient and need not be imposed: K is EVEN
    #   (the two coupled operators are identical), K <= 2L (they are 2^L-poles) and K <= 2J (Wigner-Eckart
    #   between the substates of one level).  So a J = 1/2 or J = 0 level has only K = 0 however high the
    #   multipole, and this 2p_3/2 level carries K = 0 and K = 2 and nothing else.
    #   MultipolePolarizibility.allowedRanks returns that list for any (Mp, J).
    #
    #   THE E2 ZEROS BELOW ARE CORRECT PHYSICS AND NOT A FAILURE, and they are kept because a zero one can
    #   account for is worth more than a non-zero one taken on trust.  An E2 operator acting on a p level (l = 1)
    #   reaches l' = 1 and l' = 3, i.e. p and f perturbers;  this set contains nd only, so no perturber is
    #   reached and every E2 rank vanishes identically.  `MultipolePolarizibility.checkPerturberSet` says so
    #   loudly rather than letting a bare 0.0 be read as a defect.  Note also that the rank LIST is a property of
    #   (Mp, J) alone, so E2 and E3 still report [0, 2] here even though their values are zero -- what the
    #   perturber set can reach is a separate question from what the level can carry.
    #
    #   The MEANING of alpha^(Mp,K) is fixed on an observable rather than on a convenient factor: it is the
    #   coefficient of the angular shape A_K(J,M) = C^(JM)_(JM,K0)/C^(JJ)_(JJ,K0) in the substate polarizability
    #   alpha(J,M) = sum_K A_K(J,M) alpha^(Mp,K).  Hence A_0 = 1, A_K(J,J) = 1, and for K = 2 the shape is the
    #   familiar [3M^2 - J(J+1)]/[J(2J-1)].
    grid = Basics.recommendedGrid([Configuration("2p"), Configuration("3d"), Configuration("4d"),
                                   Configuration("5d"), Configuration("6d")], Nuclear.Model(1., PointNucleus());
                                   printout=false)
    nm   = Nuclear.Model(1., PointNucleus())

    wa   = Atomic.Computation(Atomic.Computation(), name="H 2p", grid=grid, nuclearModel=nm,
                              configs=[Configuration("2p")])
    mp2p = perform(wa; output=true)[ResultKeys.Multiplet]
    ## NOT a `for` loop with an assignment: at TOP LEVEL Julia's soft scope would make that assignment
    ## local to the loop body, leaving p32 at levels[1] -- which is the 2p_1/2 level and carries only K = 0.
    p32  = mp2p.levels[ findfirst(lev -> Basics.twice(lev.J) == 3, mp2p.levels) ]
    wb   = Atomic.Computation(Atomic.Computation(), name="H nd perturbers (n=3..6)", grid=grid, nuclearModel=nm,
                              configs=[Configuration("3d"), Configuration("4d"), Configuration("5d"),
                                       Configuration("6d")])
    gMultiplet = perform(wb; output=true)[ResultKeys.Multiplet]

    println("\n>> ranks a J = 3/2 level can carry:  E1 $(MultipolePolarizibility.allowedRanks(E1, p32.J)), " *
            "E2 $(MultipolePolarizibility.allowedRanks(E2, p32.J)), E3 $(MultipolePolarizibility.allowedRanks(E3, p32.J))")
    ## THROUGH THE STANDARD ROUTE:  the multipoles asked for go on the Settings and the whole alpha^(Mp,K) table
    ## comes back on the Outcome, one entry per (multipole, rank, frequency)
    set = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings(); multipoles=[E1, E2],
                                           gMultiplet=gMultiplet)
    wc  = Atomic.Computation(Atomic.Computation(), name="H 2p alpha^(Mp,K)", grid=grid, nuclearModel=nm,
                             configs=[Configuration("2p")], propertySettings=[set])
    for  oc  in  perform(wc; output=true)[MultipolePolarizibility.Settings]
        Basics.twice(oc.level.J) == 3   ||   continue
        for  (mp, K, omega, alpha)  in  oc.alphas
            println(">> alpha^($mp,$K) = $(alpha.Babushkin) a.u.")
        end
    end
    println(">> (both E2 ranks vanish because an E2 operator reaches p and f from a p level, and this perturber " *
            "set is nd only)")
    #
elseif false
    # Last visit:  19-Oct-2026
    # Last successful:  unknown ... (the numbers below reproduce; the date is the maintainer's to set -- Rule 7)
    # Branch b: THE DYNAMIC POLARIZABILITY ON THE REAL AXIS, through the user-facing `Settings.omegas`.  This is
    #   the ac Stark / trap-depth / magic-wavelength quantity, and it costs nothing beyond the static one.
    #
    #   CLAIM, 19-Oct-2026:  H 1s with np n = 2..5 on a 129 a.u. box,  alpha^(E1,0)(omega) =
    #                        3.5525962681, 3.6134773051, 3.8095506404 a.u. at omega = 0, 0.05, 0.10 a.u.,
    #                        and -2.6987673676 a.u. at omega = 0.60, which is ABOVE the first resonance and is
    #                        therefore the real part of a complex quantity -- see the warning it triggers.
    #
    #   THE FREQUENCY ENTERS ONLY THROUGH THE DENOMINATOR: with dE = W_J' - W_J the static 1/dE becomes
    #   dE/(dE^2 - omega^2), and NOTHING else changes -- the same reduced matrix elements, the same angular
    #   coefficients, the same perturber set.  omega = 0 is always computed as well, so the Outcome's alpha0 and
    #   alpha2 keep their meaning;  the full (Mp, K, omega, alpha) list is in `outcome.alphas`.
    #
    #   The rise with omega is correct BELOW the first resonance, which here sits at dE(1s->2p) = 0.375 a.u.
    #   Above it the true alpha(omega) is COMPLEX -- the atom absorbs -- and a sum over a DISCRETE perturber set
    #   can only return a principal-value-like real part with no width;  `MultipolePolarizibility.checkFrequency`
    #   says so loudly, and the last frequency below is there to make it fire.
    grid = Radial.Grid(Radial.Grid(true); rbox = 129.0)
    nm   = Nuclear.Model(1., PointNucleus())

    wa1 = Atomic.Computation(Atomic.Computation(), name="H perturbers (np, n=2..5)", grid=grid, nuclearModel=nm,
                             configs=[Configuration("2p"), Configuration("3p"), Configuration("4p"),
                                      Configuration("5p")])
    gMultiplet = perform(wa1; output=true)[ResultKeys.Multiplet]

    polSettings = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings(), multipoles=[E1],
                        gMultiplet=gMultiplet, omegas=[0.05, 0.10, 0.60], printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="H(1s) dynamic dipole polarizability", grid=grid,
                            nuclearModel=nm, configs=[Configuration("1s")], propertySettings=[polSettings])
    outcome = perform(wa; output=true)[MultipolePolarizibility.Settings][1]

    println("\n>> alpha^(Mp,K)(omega), as (multipole, rank, omega, value):")
    for  (mp, K, omega, alpha) in outcome.alphas
        println(">>    alpha^($mp,$K) at omega = $omega  ->  $(alpha.Babushkin) a.u.")
    end
    println(">> the omega = 0.60 entry lies ABOVE the first resonance at 0.375 a.u. and is therefore the real")
    println(">> part of a complex quantity;  the warning printed above is the point of including it.")
    #
elseif false
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (3 of 3 numbers on the claim line, after the branch was converted to the
    #   standard perform route:  4.499751433446605 a.u., and the ratios 0.9999995158 and 0.9999545283 -- every
    #   quoted digit unchanged, which is the evidence that the conversion was faithful)
    # Last successful:  19-Oct-2026   ## alpha^(E1,0) 4.4997514334 and C_6 6.4983890995 a.u., at 1.4e-8 and
    #   4.8e-7 of the Dirac reference;  C_3 0.2499886321 a.u., which is 0.854 alpha^2 BELOW the
    #   NONRELATIVISTIC exact 1/4 -- the gap is the relativistic correction, not an unconverged sum.
    # Branch c: THE IMAGINARY AXIS, AND THE TWO COEFFICIENTS THAT NEED ONLY THE DIPOLE RESPONSE -- the
    #   van der Waals C_6 of H-H and the atom-wall C_3 of H.  Both have exact comparison values, and C_3's
    #   needs no literature at all.
    #
    #   CLAIM, 19-Oct-2026:  H 1s with np n = 2..16 on a 64 a.u. box,
    #       alpha^(E1,0)(0) = 4.4997514334 a.u.   (reference 4.499751495178,  0.999999986)
    #       C_6(H-H)        = 6.4983890995 a.u.   (Dirac    6.498392245754,   0.999999516)
    #       C_3(H, perfect conductor) = 0.2499886321 a.u.  (the nonrel. exact 1/4,  0.999954528)
    #
    #   WHY THE THIRD RATIO IS THE INTERESTING ONE, AND WHY IT IS NOT A CONVERGENCE FAILURE.  The first two
    #   comparisons are Dirac against Dirac and agree to 1.4e-8 and 4.8e-7, so they say the machinery is right
    #   but nothing about the physics.  The third compares a DIRAC C_3 with the NONRELATIVISTIC exact 1/4, which
    #   follows from <r^2> = 3 a.u. for H 1s -- two different quantities, and the 4.5e-5 gap is therefore
    #   EXPECTED rather than a defect.  Its size settles which it is:  4.5472e-05 is 0.854 alpha^2, i.e. exactly
    #   the order a relativistic correction to <r^2> must have.  Read as incompleteness of the 15-member np sum
    #   it would instead be an unexplained coincidence with alpha^2.
    #
    #   AND THE SAME TEST TIES THE THREE NUMBERS TOGETHER.  The Dirac reference's own alpha is 1.037 alpha^2
    #   below the exact nonrelativistic 9/2, and its C_6 is 1.768 TIMES that relative shift below the exact
    #   nonrelativistic 6.49902670540.  A factor between 1 and 2 is what a double sum over a response shifted
    #   NON-uniformly in frequency must give -- 2 exactly would mean a uniform rescaling of alpha(i omega).
    #   So the three numbers are consistent with one another through a single correction, which no amount of
    #   reproducing them separately could establish.
    #
    #   ON THE IMAGINARY AXIS THERE ARE NO POLES AT ALL -- the denominator is dE^2 + omega^2 -- which is why a
    #   dispersion coefficient is built there and why none of the real-axis warnings apply.
    #
    #   AND NO NUMERICAL omega QUADRATURE IS NEEDED.  On a discrete perturber spectrum each polarizability is a
    #   sum of single poles and int_0^inf [a/(a^2+w^2)][b/(b^2+w^2)] dw = pi/(2(a+b)) exactly, so
    #   C_6 = (3/2) sum_ab A_a A_b/(dE_a + dE_b) -- Dalgarno's oscillator-strength form -- and
    #   C_3 = (1/8) sum_a A_a.  The second is checkable by hand: with A_a = (2/3)|<a|d|0>|^2 the dipole sum rule
    #   gives C_3 = <d^2>/12, which is also what an electrostatic image dipole yields, and <r^2> = 3 for H 1s.
    #
    #   THE RADIAL BOX IS 64 a.u. ON PURPOSE, AND THIS IS THE OPPOSITE OF RULE 12'S USUAL ADVICE.  A sum over
    #   states needs the UPPER members of its series to stand in for the CONTINUUM.  A box matched to the highest
    #   orbital leaves them as genuine bound states, and the set is then bound-only: measured on this very series,
    #   `Basics.recommendedGrid` gives 775 a.u. and alpha^(E1,0) drops to 0.861 of the reference while C_3 drops
    #   to 0.770 -- see branch e, which makes the guard fire on exactly that.
    grid = Radial.Grid(Radial.Grid(true); rbox = 64.0)
    nm   = Nuclear.Model(1., PointNucleus())
    asf  = AsfSettings(AsfSettings(); gridStopper=false)       ## the high-n members are pseudo-continuum, not orbitals

    confs = Configuration[];    for  n = 2:16    push!(confs, Configuration("$(n)p"))    end
    wa1  = Atomic.Computation(Atomic.Computation(), name="H 1s", grid=grid, nuclearModel=nm,
                              configs=[Configuration("1s")], asfSettings=asf)
    level = perform(wa1; output=true)[ResultKeys.Multiplet].levels[1]
    wa2  = Atomic.Computation(Atomic.Computation(), name="H perturbers (np, n=2..16)", grid=grid, nuclearModel=nm,
                              configs=confs, asfSettings=asf)
    gMultiplet = perform(wa2; output=true)[ResultKeys.Multiplet]

    ## THROUGH THE STANDARD ROUTE:  the perturber multiplet, the orders and the wall go on the Settings, and
    ## computeOutcomes reports the saturation of the set and displays the results itself
    set = DispersionCoefficient.Settings(DispersionCoefficient.Settings(); orders=[6], multipoles=[E1],
                              calcAtomWall=true, gMultiplet=gMultiplet)
    wa3 = Atomic.Computation(Atomic.Computation(), name="H(1s) C_6 and C_3", grid=grid, nuclearModel=nm,
                              configs=[Configuration("1s")], propertySettings=[set], asfSettings=asf)
    oc  = perform(wa3; output=true)[DispersionCoefficient.Settings][1]

    ## and the scan below is ARITHMETIC on the terms that one computation stored -- a sum of single poles
    ## evaluated at a new frequency.  There is no Settings field for a list of imaginary frequencies and none is
    ## wanted:  a second perform would repeat the SCF to learn nothing new.
    terms = DispersionCoefficient.selectTerms(oc.partnerA, 1)
    println("\n>> alpha^(E1,0)(i omega) on the imaginary axis -- positive, monotone, and no poles:")
    for  omega in [0.0, 0.1, 0.3, 1.0, 3.0, 10.0]
        println(">>    omega = $omega  ->  $(DispersionCoefficient.polarizability(terms, omega)) a.u.")
    end
    println(">> C_6 / published Dirac 6.498392245754  =  $(oc.coefficients[1].value/6.498392245754)")
    println(">> C_3 / the exact 1/4                   =  $(oc.atomWallC3/0.25)")
    #
elseif false
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (3 of 3 numbers on the claim line, after the branch was converted to the
    #   standard perform route with ONE combined perturber computation in place of three:  0.9999995158,
    #   1.0000009257, 1.0000234525 -- every quoted digit unchanged, and C_6 bit-identical to branch c's)
    # Last successful:  unknown ... (the numbers below reproduce; the date is the maintainer's to set -- Rule 7)
    # Branch d: C_6, C_8 AND C_10 OF H(1s)-H(1s) TOGETHER, which is the full Casimir-Polder machinery and needs
    #   the E1, E2 and E3 responses at once.  This is the branch that anchors the module.
    #
    #   CLAIM, 19-Oct-2026, each series on its own 64 a.u. box (np n = 2..16, nd n = 3..16, nf n = 4..22):
    #       C_6  =    6.4983890995 a.u.   against the Dirac    6.498392245754   ratio 1.000000
    #       C_8  =  124.3841457941 a.u.   against             124.384030658     ratio 1.000001
    #       C_10 = 3285.4478384364 a.u.   against            3285.370788289     ratio 1.000023
    #   (Tang, Zhang, Zhang, Jiang & Mitroy, arXiv:1205.6883 -- the same paper that anchors the multipole
    #   polarizabilities themselves, so the comparison is internally consistent.)
    #
    #   THE WHOLE RELATION IS ONE BINOMIAL, not a table of special cases.  With n = l1 + l2 + 1,
    #       C_2n = (1/2 pi) sum over ORDERED (l1,l2) of binomial(2 l1 + 2 l2, 2 l1) int alpha_l1^A alpha_l2^B dw
    #   which reproduces the four prefactors the literature writes out by hand: 3/pi for C_6, 15/2 pi for each
    #   of the two C_8 terms, and 28/pi (the (1,3) and (3,1) pair) plus 35/pi ((2,2)) for C_10.  The printed
    #   table shows those pieces, which is what says WHICH response has to be converged for WHICH coefficient.
    #
    #   This branch runs three deep perturber series and takes a few minutes.
    nm   = Nuclear.Model(1., PointNucleus())
    asf  = AsfSettings(AsfSettings(); gridStopper=false)
    grid = Radial.Grid(Radial.Grid(true); rbox = 64.0)

    ## THROUGH THE STANDARD ROUTE, AND NOTE WHAT THAT COSTS HERE:  `Settings` carries ONE gMultiplet for all
    ## multipoles, so the three series go in together and each multipole picks the members its own selection
    ## rules allow -- E1 the np, E2 the nd, E3 the nf -- with the rest dropping out.  That is one perturber
    ## computation of 48 configurations instead of three of 15, 14 and 19;  the orbitals are the same, since for
    ## a one-electron system each subshell is solved in the nuclear field independently of the others, and the
    ## claim line above is what checks that -- measured 20-Oct-2026, all three ratios unchanged to every quoted
    ## digit.  FOR A MANY-ELECTRON ATOM THE COMBINATION IS NOT FREE:  one SCF over all three series is a
    ## different mean field from three separate ones, so there the perturber sets would have to be chosen with
    ## that in mind rather than simply concatenated.
    confs = vcat( [Configuration("$(n)p")  for n = 2:16], [Configuration("$(n)d")  for n = 3:16],
                  [Configuration("$(n)f")  for n = 4:22] )
    gm  = perform( Atomic.Computation(Atomic.Computation(), name="H np + nd + nf perturbers", grid=grid,
                   nuclearModel=nm, configs=confs, asfSettings=asf); output=true )[ResultKeys.Multiplet]
    set = DispersionCoefficient.Settings(DispersionCoefficient.Settings(); orders=[6, 8, 10],
                   multipoles=[E1, E2, E3], calcAtomWall=true, gMultiplet=gm)
    wa  = Atomic.Computation(Atomic.Computation(), name="H(1s) C_6, C_8, C_10", grid=grid, nuclearModel=nm,
                   configs=[Configuration("1s")], propertySettings=[set], asfSettings=asf)
    oc  = perform(wa; output=true)[DispersionCoefficient.Settings][1]
    for  (c, ref) in zip(oc.coefficients, [6.498392245754, 124.384030658, 3285.370788289])
        println(">> C_$(c.n) / published = $(c.value/ref)")
    end
    #
elseif false
    # Last visit:  19-Oct-2026 (RAISED AS DESIGNED -- this branch is a guard test)
    # Last successful:  unknown ... (a guard test; the date is the maintainer's to set -- Rule 7)
    # Branch e: THE TWO GUARDS, AND THIS BRANCH IS SUPPOSED TO RAISE.  Completing it would be the defect.
    #
    #   (1) AN EXCITED PARTNER IS REFUSED, not warned about.  The Casimir-Polder construction integrates along
    #       the imaginary axis because every excitation energy then carries the SAME sign.  If a partner is
    #       excited at least one dE is negative, the contour argument fails, and the true interaction gains an
    #       additional REAL-frequency (resonant) term that no imaginary-axis integral contains (Lepers & Dulieu,
    #       arXiv:1703.02833, Eq. 75).  Since excited levels are exactly what JAC is for, this is the COMMON case
    #       here and not an edge one -- and the number that would otherwise come back is finite, smooth and
    #       wrong, which is why `DispersionCoefficient.checkGroundLevel` raises.
    #
    #   (2) A BOUND-ONLY PERTURBER SET IS WARNED ABOUT, by the span of its excitation energies.  A bound Rydberg
    #       series accumulates at the ionization threshold, so max(dE)/min(dE) stays near 1; a series running
    #       into the continuum spreads over orders of magnitude.  Measured on np n = 2..16: span 54.6 in a 64 a.u.
    #       box (alpha 1.00000 of reference, C_3 0.99995 of exact) against span 1.38 in a 775 a.u. box
    #       (alpha 0.86054, C_3 0.77026).  The second half of this branch makes that fire.
    nm  = Nuclear.Model(1., PointNucleus())
    asf = AsfSettings(AsfSettings(); gridStopper=false)

    println("\n>> (2) first the saturation guard, on np n = 2..16 in a 775 a.u. box:")
    gridBad  = Radial.Grid(Radial.Grid(true); rbox = 775.0)
    confsBad = Configuration[];   for  n = 2:16    push!(confsBad, Configuration("$(n)p"))    end
    levBad = perform(Atomic.Computation(Atomic.Computation(), name="H 1s", grid=gridBad, nuclearModel=nm,
                     configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]
    gmBad  = perform(Atomic.Computation(Atomic.Computation(), name="np, too large a box", grid=gridBad,
                     nuclearModel=nm, configs=confsBad, asfSettings=asf); output=true)[ResultKeys.Multiplet]
    bad    = DispersionCoefficient.Partner("H(1s), 775 a.u. box", levBad,
                     [DispersionCoefficient.multipoleTerms(E1, levBad, gmBad, gridBad)])
    DispersionCoefficient.checkSaturation(bad)
    println(">> alpha^(E1,0) = $(DispersionCoefficient.staticPolarizability(DispersionCoefficient.selectTerms(bad,1)))" *
            "  against the reference 4.499751495178")
    println(">> C_3          = $(DispersionCoefficient.computeAtomWallCoefficient(bad, Inf))  against the exact 0.25")

    println("\n>> (1) now the ground-level guard: H 2p_3/2 with 1s among its perturbers MUST raise.")
    grid  = Basics.recommendedGrid([Configuration("2p"), Configuration("1s"), Configuration("3d")], nm;
                                    printout=false)
    mp2p  = perform(Atomic.Computation(Atomic.Computation(), name="H 2p", grid=grid, nuclearModel=nm,
                    configs=[Configuration("2p")], asfSettings=asf); output=true)[ResultKeys.Multiplet]
    p32   = mp2p.levels[ findfirst(lev -> Basics.twice(lev.J) == 3, mp2p.levels) ]   ## soft scope, see branch a
    gmx   = perform(Atomic.Computation(Atomic.Computation(), name="H 1s + 3d perturbers", grid=grid,
                    nuclearModel=nm, configs=[Configuration("1s"), Configuration("3d")], asfSettings=asf);
                    output=true)[ResultKeys.Multiplet]
    excited = DispersionCoefficient.Partner("H(2p_3/2)", p32,
                    [DispersionCoefficient.multipoleTerms(E1, p32, gmx, grid)])
    DispersionCoefficient.computeCoefficient(6, excited, excited)
    println(">> IF THIS LINE IS REACHED THE GUARD IS BROKEN -- the call above must have raised.")
    #
elseif true
    # Last visit:  20-Oct-2026
    # Last successful:  unknown ... (Rule 7: the conductor limit is checked against an EXACT 1/4 that needs no
    #   literature, and the three limits and the per-pole factor are algebra, but the silica-like number itself
    #   has no published counterpart -- so the judgement, and the date, are the maintainer's)
    #
    #   CLAIM, 20-Oct-2026:  H(1s) with np n = 2..16 on a 64 a.u. box, eps(0) = 3.9 throughout.
    #       THROUGH `perform`:  C_6 = 6.49838910e+00 a.u., the same number branch d anchors against the published
    #           Dirac value, and against a silica-like wall (w_0 = 10 eV = 0.3675 a.u., hence Omega = 0.5752 a.u.)
    #           C_3 = 8.18504949e-02 a.u.;
    #       THE THREE LIMITS:  w_0 -> Inf gives 0.14795245571963425, BIT-IDENTICAL to the constant-eps value, so
    #           the new path does not disturb the old one;  eps(0) -> Inf gives 0.2499886320780027 against the
    #           EXACT 1/4, a ratio of 0.999955 which is the perturber set and not the wall;  and w_0 = 1e-6 gives
    #           5.2105e-07, falling LINEARLY in w_0 toward zero;
    #       WHAT THE RESONANCE COSTS, C_3(w_0)/C_3(const eps):  0.203400 (2 eV), 0.386166 (5), 0.553222 (10),
    #           0.708593 (20), 0.856002 (50), 0.958810 (200 eV) -- so at the silica-like 10 eV the static
    #           dielectric constant alone OVERSTATES C_3 by 81 %, and only by 200 eV is the limit nearly reached;
    #       and eps(0) = 0.5 RAISES.
    #       Every figure agrees with tools/probe-dielectricWall.jl to better than 1e-8 by an independent route --
    #           that probe calls the coefficient directly, this branch arrives through perform.
    # Branch f: A WALL THAT HAS A RESONANCE, and the point of the branch is that it still needs NO QUADRATURE.
    #   Added 20-Oct-2026 as DispersionCoefficient.Settings.wallResonance;  this is also the only branch in this
    #   file that goes through `perform`, so it exercises the whole user-facing path -- computeOutcomes, the
    #   ground-level guard, the new wall guard and displayResults -- where branch d calls computeOutcome directly.
    #
    #   WHY THERE IS NO INTEGRAL TO DO.  The atom-wall coefficient is
    #
    #       C_3 = (1/4 pi) int_0^inf alpha^(E1,0)(i w) r(i w) dw
    #
    #   and ONE Lorentz oscillator, eps(i w) = 1 + w_p^2/(w_0^2 + w^2), has
    #
    #       r(i w) = (eps-1)/(eps+1) = R/(Omega^2 + w^2),    R = w_p^2/2,   Omega^2 = w_0^2 + w_p^2/2
    #
    #   which is the SAME single-pole shape as one polarizability term.  Since
    #   int_0^inf dw/((a^2+w^2)(b^2+w^2)) = pi/(2ab(a+b)), the whole thing collapses to a double sum,
    #
    #       C_3 = (R/8) sum_a A_a / [Omega (dE_a + Omega)],   Omega = w_0 sqrt((eps(0)+1)/2),
    #                                                         R     = (eps(0)-1) w_0^2/2
    #
    #   so the frequency dependence of a real dielectric costs exactly nothing beyond the static case.  That is
    #   the same reason C_6, C_8 and C_10 need no quadrature in branch d, applied to the wall instead of to the
    #   second atom.  `DispersionCoefficient.wallReflectionPole` is the one line that turns (eps(0), w_0) into
    #   (R, Omega).
    #
    #   AND IT SHOWS THE CONSTANT-eps OPTION TO BE A LIMIT RATHER THAN A MODEL.  Scaling C_3 by (eps-1)/(eps+1),
    #   which is what this module shipped first, is the Omega >> dE case -- a wall whose charges respond
    #   INSTANTANEOUSLY on the atom's own frequency scale -- and its error is one factor per pole, exactly
    #   1 + dE_a/Omega.  For a silica-like wall that is not a small correction, as the scan below shows.
    #
    #   THREE LIMITS, EACH A PHYSICAL WALL, and the third is the one worth pausing on:
    #       w_0 -> Inf      the instantaneous-response limit, i.e. the constant-eps answer;
    #       eps(0) -> Inf   a perfect conductor, reflecting fully at EVERY frequency, so no resonance can enter;
    #       w_0 -> 0        ZERO -- a wall whose charges cannot follow the atom does not reflect at the
    #                       frequencies that carry C_3, which is why a wall is not characterised by eps(0) alone.
    #
    #   THE BOX IS 64 a.u. AND THAT IS DELIBERATELY SMALL -- see branch c.  C_3 is the quantity a poor perturber
    #   set damages FIRST, because its sum carries no energy denominator at all, so the saturation report printed
    #   below is the first thing to read.
    nm   = Nuclear.Model(1., PointNucleus())
    asf  = AsfSettings(AsfSettings(); gridStopper=false)
    grid = Radial.Grid(Radial.Grid(true); rbox = 64.0)
    eV   = Defaults.convertUnits("energy: from eV to atomic", 1.0)

    confs = [Configuration("$(n)p")  for n = 2:16]
    gm    = perform(Atomic.Computation(Atomic.Computation(), name="H np perturbers", grid=grid, nuclearModel=nm,
                    configs=confs, asfSettings=asf); output=true)[ResultKeys.Multiplet]

    ## the user-facing route:  the wall is two numbers on the Settings and C_3 comes back on the Outcome
    set = DispersionCoefficient.Settings(DispersionCoefficient.Settings(); orders=[6], multipoles=[E1],
                     calcAtomWall=true, wallDielectric=3.9, wallResonance=10.0*eV, gMultiplet=gm)
    wa  = Atomic.Computation(Atomic.Computation(), name="H(1s) against a silica-like wall", grid=grid,
                     nuclearModel=nm, configs=[Configuration("1s")], propertySettings=[set], asfSettings=asf)
    oc  = perform(wa; output=true)[DispersionCoefficient.Settings][1]

    ## everything else is arithmetic on the terms that one computation already produced -- no further SCF
    partner = oc.partnerA
    cConst  = DispersionCoefficient.computeAtomWallCoefficient(partner, 3.9)
    println("\n>> ===== the three limits =====")
    println(">>   w_0 -> Inf   (instantaneous response) = " *
            "$(DispersionCoefficient.computeAtomWallCoefficient(partner, 3.9, Inf))" *
            "   against the constant-eps $cConst")
    println(">>   eps(0) -> Inf (perfect conductor)     = " *
            "$(DispersionCoefficient.computeAtomWallCoefficient(partner, Inf, 10.0*eV))" *
            "   against the EXACT 0.25")
    println(">>   w_0 -> 0     (charges cannot follow)  = " *
            "$(DispersionCoefficient.computeAtomWallCoefficient(partner, 3.9, 1.0e-6))" *
            "   which must fall to zero like w_0")

    println("\n>> ===== what the resonance costs:  eps(0) = 3.9 throughout =====")
    println(">>   w_0 [eV]    Omega [a.u.]      C_3(w_0)             C_3(const eps)       ratio")
    for  w0eV  in  (2.0, 5.0, 10.0, 20.0, 50.0, 200.0)
        c3 = DispersionCoefficient.computeAtomWallCoefficient(partner, 3.9, w0eV*eV)
        Om = DispersionCoefficient.wallReflectionPole(3.9, w0eV*eV)[2]
        println(@sprintf(">>   %-11.1f %-17.6f %-20.12f %-20.12f %.6f", w0eV, Om, c3, cConst, c3/cConst))
    end
    println(">>   the ratio is the per-pole factor 1/(1 + dE_a/Omega) averaged over the set;  H(1s) -> 2p sits " *
            "at dE = 0.375 a.u., so a wall whose Omega is comparable to that is nowhere near the limit.")

    ## and the guard, which now runs BEFORE any perturber set is built
    println("\n>> ===== the guard =====")
    try
        DispersionCoefficient.checkWall( DispersionCoefficient.Settings(DispersionCoefficient.Settings();
                                         calcAtomWall=true, wallDielectric=0.5) )
        println(">>   eps(0) = 0.5 was ACCEPTED -- the guard is broken;  (eps-1)/(eps+1) is negative there.")
    catch
        println(">>   eps(0) = 0.5 raised, as it must:  no dielectric has a negative reflection factor, and the " *
                "refusal arrives before the first perturber set is built rather than after it.")
    end
    #
end
