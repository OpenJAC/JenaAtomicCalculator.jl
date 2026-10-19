#
# Validation of module-DispersionCoefficient.jl, added 19-Oct-2026 (priority item 80).  Five cases, each
# testing something the previous one cannot:
#
#   CASE A  THE COMBINATORICS ALONE, with no atomic input at all.  `casimirPolderPrefactor` summed over
#           ORDERED (l1,l2) pairs must reproduce the four prefactors the literature writes by hand:
#           3/pi (C_6), 15/2pi (each C_8 term), 28/pi and 35/pi (the two C_10 terms).
#   CASE B  THE CLOSED-FORM OVERLAP INTEGRAL against a numerical quadrature of the same integrand.  This
#           is what justifies having no omega quadrature in `src/`.
#   CASE C  THE ATOM-WALL C_3 against an EXACT value that needs no literature:  with A_a = (2/3)|<a|d|0>|^2
#           the dipole sum rule gives C_3 = <d^2>/12, which is also what an electrostatic image dipole
#           yields, so C_3(H 1s) = 3/12 = 1/4 EXACTLY.  The measured shortfall of a bound-only set is the
#           calibration of `checkSaturation` and is the reason that guard exists.
#   CASE D  END TO END: H(1s)-H(1s) C_6, C_8 and C_10 against Tang, Zhang, Zhang, Jiang & Mitroy,
#           arXiv:1205.6883 -- 6.498392245754, 124.384030658, 3285.370788289 (Dirac).  Each alpha^(E l,0)
#           that feeds them is printed against its own reference from the same paper, so a coefficient that
#           misses can be attributed to the response under it rather than to the module.
#   CASE E  THE GUARDS.  The ground-level guard must RAISE for an excited partner (it is refused, not
#           warned about, per item 80), a C_n missing a multipole must be reported INCOMPLETE, and the
#           spectral-span guard must fire on a bound-only perturber set and stay silent on a good one.
#
# RESULTS, 19-Oct-2026.  H 1s with np n = 2..16, nd n = 3..16, nf n = 4..22, each on a 64 a.u. box:
#   alpha^(E1,0) 4.4997514334 (1.00000 of ref), alpha^(E2,0) 14.998853851 (1.00000), alpha^(E3,0)
#   131.24328699 (1.00004).
#
#   CASE A  all six ordered-pair prefactors reproduce the literature to 0.00e+00.
#   CASE B  the closed-form overlap integral against a 40000-point quadrature of the same integrand:
#           1.6e-10 worst over the (1,1), (1,2), (2,2) and (1,3) pairs, which is the quadrature's own
#           error and not the module's.  So no omega quadrature is needed in `src/`.
#   CASE C  C_3 = 0.2499886321 against the EXACT 1/4 -- ratio 0.999955, and sum_a A_a = 1.9999090566
#           against the exact 2 = (2/3)<r^2>.  Against a wall of eps = 2, C_3 = 0.0833295440 = 1/3 of it.
#   CASE D  C_6  =    6.4983890995   against the Dirac    6.498392245754   ratio 1.000000
#           C_8  =  124.3841457941   against             124.384030658     ratio 1.000001
#           C_10 = 3285.4478384364   against            3285.370788289     ratio 1.000023
#           -- so the whole chain reproduces three published dispersion coefficients, the worst of them
#           to 2.3e-05, and item 80's DONE WHEN is met.
#   CASE E  the ground-level guard RAISES on H 2p_3/2 with 1s among the perturbers (one dE is then
#           -0.375 a.u.);  a C_10 built from the E1 response alone reports complete = false and keeps 0
#           of its 3 ordered pairs;  and the span guard fires on a 775 a.u. box and is silent on 64 a.u.
#
# THE BOX IS THE ONE THING TO GET RIGHT HERE, AND IT IS THE OPPOSITE WAY ROUND FROM RULE 12'S USUAL
# ADVICE.  A sum over states needs the upper members of its series to stand in for the CONTINUUM;  a box
# matched to the HIGHEST orbital leaves them as genuine bound states, and the set is then bound-only.
# Measured, one np n = 2..16 / nd n = 3..16 / nf n = 4..22 set per multipole, alpha/reference:
#
#        box      40      64     100     129     200     400     775(recommendedGrid)
#        E1   1.00000 1.00000 1.00000 1.00000 0.99992 0.98888 0.86054
#        E2   1.00000 1.00000 1.00000 0.99999 0.99952 0.91838 0.44283
#        E3   1.00004 1.00004 1.00001 0.99999 1.00059 1.02861 1.27518
#
# E1 and E2 fail DOWNWARD (the continuum is missing) and E3 UPWARD (the r^L weighting of a too-large box),
# so there is an optimum rather than a direction, and 40-100 a.u. is the plateau for all three.
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra

# arXiv:1205.6883, Dirac;  the same paper that anchored the polarizabilities themselves
const ALPHA_REF = Dict(1 => 4.499751495177639, 2 => 14.998829822856442, 3 => 131.237821447844662)
const CN_DIRAC  = Dict(6 => 6.498392245754, 8 => 124.384030658, 10 => 3285.370788289)
const CN_NONREL = Dict(6 => 6.499026705406, 8 => 124.399084,    10 => 3285.828415)
const C3_EXACT  = 0.25                      # <r^2>/12 with <r^2> = 3 for H 1s

function terms1s(multipole::EmMultipole, letter::String, nmax::Int64, nm::Nuclear.Model)
    confs = [Configuration("1s")]
    for  n = (multipole.L+1):nmax    push!(confs, Configuration("$n$letter"))    end
    # THE BOX IS DELIBERATELY SMALL, and the box scan in the header above is why:  a box matched to the
    # HIGHEST orbital of the series leaves every member a genuine bound state, and the sum then misses the
    # continuum entirely.  64 a.u. is on the plateau for all three multipoles.  The grid guard is overridden
    # because the upper members are not meant to be accurate orbitals -- they stand in for the continuum.
    grid = Radial.Grid(Radial.Grid(true); rbox = 64.0)
    asf  = AsfSettings(AsfSettings(); gridStopper=false)
    lev  = perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=grid, nuclearModel=nm,
                   configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]
    gm   = perform(Atomic.Computation(Atomic.Computation(); name="$letter series", grid=grid, nuclearModel=nm,
                   configs=confs[2:end], asfSettings=asf); output=true)[ResultKeys.Multiplet]
    return( (DispersionCoefficient.multipoleTerms(multipole, lev, gm, grid), lev, grid) )
end

function main()
    nm = Nuclear.Model(1., PointNucleus())

    println("\n", "="^120)
    println(" CASE A -- the Casimir-Polder combinatorics alone:  sum over ORDERED (l1,l2) of the prefactor")
    println("="^120)
    @printf("   %-6s %-14s %22s %22s %12s\n", "C_n", "(l1,l2)", "prefactor", "literature", "rel.diff")
    for  (n, pairs, lits, names) in [ (6,  [(1,1)],         [3/pi],          ["3/pi"]),
                                      (8,  [(1,2),(2,1)],   [15/(2pi), 15/(2pi)], ["15/2pi", "15/2pi"]),
                                      (10, [(1,3),(3,1)],   [14/pi, 14/pi],  ["14/pi", "14/pi"]),
                                      (10, [(2,2)],         [35/pi],         ["35/pi"]) ]
        for  (i, (l1,l2)) in enumerate(pairs)
            pf = DispersionCoefficient.casimirPolderPrefactor(l1, l2)
            @printf("   %-6s %-14s %22.14e %22.14e %12.2e\n", "C_$n", "($l1,$l2)", pf, lits[i],
                    abs(pf - lits[i])/lits[i])
        end
    end
    println("   >> (1,3)+(3,1) together give 28/pi, which is how the literature quotes that C_10 term.")

    println("\n", "="^120)
    println(" CASE C/D -- building H 1s:  alpha^(E1,0) from np, alpha^(E2,0) from nd, alpha^(E3,0) from nf,")
    println("             each on ITS OWN matched box, which is what the Partner constructor is for")
    println("="^120)
    t1, lev1, g1 = terms1s(E1, "p", 16, nm)
    t2, _,    g2 = terms1s(E2, "d", 16, nm)
    t3, _,    g3 = terms1s(E3, "f", 22, nm)
    @printf("\n   %-6s %-14s %8s %7s %20s %22s %10s\n", "Mp", "perturbers", "box", "poles", "alpha^(Mp,0)",
            "Dirac reference", "of ref")
    for  (t, lab, gr) in [(t1, "2p..16p", g1), (t2, "3d..16d", g2), (t3, "4f..22f", g3)]
        a0 = DispersionCoefficient.staticPolarizability(t);     ref = ALPHA_REF[t.multipole.L]
        @printf("   %-6s %-14s %8.1f %7d %20.10e %22.12e %10.5f\n", string(t.multipole), lab, gr.r[end],
                length(t.strengths), a0, ref, a0/ref)
    end
    hydrogen = DispersionCoefficient.Partner("H(1s)", lev1, [t1, t2, t3])

    println("\n", "="^120)
    println(" CASE B -- the closed-form overlap integral against a numerical quadrature of the SAME integrand")
    println("="^120)
    @printf("   %-12s %22s %22s %12s\n", "(l1,l2)", "closed form", "quadrature", "rel.diff")
    for  (ta, tb, lab) in [(t1,t1,"(1,1)"), (t1,t2,"(1,2)"), (t2,t2,"(2,2)"), (t1,t3,"(1,3)")]
        closed = DispersionCoefficient.overlapIntegral(ta, tb)
        nq = 40000;   acc = 0.0
        for  k = 1:nq
            tt = (k - 0.5)/nq;    w = tt/(1-tt);    jac = 1/(1-tt)^2
            acc = acc + DispersionCoefficient.polarizability(ta, w) *
                        DispersionCoefficient.polarizability(tb, w) * jac / nq
        end
        @printf("   %-12s %22.14e %22.14e %12.2e\n", lab, closed, acc, abs(acc-closed)/closed)
    end

    println("\n", "="^120)
    println(" CASE C -- the atom-wall C_3 against the EXACT 1/4 (= <r^2>/12, no literature needed)")
    println("="^120)
    c3 = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, Inf)
    @printf("   C_3 (perfect conductor) = %.10f   against the exact %.10f   ratio %.6f\n", c3, C3_EXACT,
            c3/C3_EXACT)
    @printf("   sum_a A_a = %.10f   against the exact 2 = (2/3)<r^2>   ratio %.6f\n", sum(t1.strengths),
            sum(t1.strengths)/2.0)
    println("   >> THE SHORTFALL IS THE POINT:  the SAME perturber set gives alpha^(E1,0) to 5 digits (above)")
    println("      and C_6 to 5 digits (below), because those weight each pole by 1/dE and 1/(dE_a+dE_b);")
    println("      sum_a A_a weights it by NOTHING, so C_3 is the least saturated quantity in the module.")
    @printf("   C_3 against a wall of eps = 2 (reflection 1/3) = %.10f\n",
            DispersionCoefficient.computeAtomWallCoefficient(hydrogen, 2.0))

    println("\n", "="^120)
    println(" CASE D -- END TO END:  H(1s)-H(1s) atom-atom coefficients against arXiv:1205.6883")
    println("="^120)
    set = DispersionCoefficient.Settings(DispersionCoefficient.Settings(); calcAtomWall=true)
    oc  = DispersionCoefficient.computeOutcome(hydrogen, hydrogen, set)
    DispersionCoefficient.displayResults(stdout, oc)
    @printf("\n   %-6s %22s %22s %10s %22s %10s\n", "C_n", "this module", "Dirac", "of Dirac",
            "non-relativistic", "of n.r.")
    for  c in oc.coefficients
        @printf("   %-6s %22.10f %22.10f %10.6f %22.10f %10.6f\n", "C_$(c.n)", c.value, CN_DIRAC[c.n],
                c.value/CN_DIRAC[c.n], CN_NONREL[c.n], c.value/CN_NONREL[c.n])
    end

    println("\n", "="^120)
    println(" CASE E -- THE GUARDS")
    println("="^120)
    println("   (E1) an EXCITED partner must be REFUSED, not warned about:")
    # H 2p_3/2 with 1s among the perturbers:  one excitation energy is then NEGATIVE
    confs = [Configuration("2p"), Configuration("1s"), Configuration("3d")]
    gr    = Basics.recommendedGrid(confs, nm; printout=false)
    asf   = AsfSettings(AsfSettings(); gridStopper=false)
    mp2p  = perform(Atomic.Computation(Atomic.Computation(); name="2p", grid=gr, nuclearModel=nm,
                    configs=[Configuration("2p")], asfSettings=asf); output=true)[ResultKeys.Multiplet]
    p32   = mp2p.levels[1];    for l in mp2p.levels    Basics.twice(l.J) == 3 && (p32 = l)    end
    gmx   = perform(Atomic.Computation(Atomic.Computation(); name="1s+3d", grid=gr, nuclearModel=nm,
                    configs=[Configuration("1s"), Configuration("3d")], asfSettings=asf);
                    output=true)[ResultKeys.Multiplet]
    excited = DispersionCoefficient.Partner("H(2p_3/2)", p32, [DispersionCoefficient.multipoleTerms(E1, p32, gmx, gr)])
    try
        DispersionCoefficient.computeCoefficient(6, excited, excited)
        println("   >> FAILED: no error was raised, which means the guard is not working.")
    catch err
        println("   >> RAISED AS DESIGNED:");   println("      ", sprint(showerror, err))
    end
    println("\n   (E2) a C_10 asked for with the E1 response ALONE must be reported INCOMPLETE:")
    thin = DispersionCoefficient.Partner("H(1s), E1 only", lev1, [t1])
    c10  = DispersionCoefficient.computeCoefficient(10, thin, thin)
    @printf("   >> C_10 = %.6e   complete = %s   pieces kept: %d of the 3 ordered pairs\n", c10.value,
            c10.complete, length(c10.pieces))
    println("\n   (E3) the SPAN guard:  it must fire on a bound-only set and stay silent on the good one.")
    println("        First the good 64 a.u. set -- NOTHING may be printed between here and the next label:")
    DispersionCoefficient.checkSaturation(hydrogen)
    println("        Now the same np n = 2..16 series in a 775 a.u. box, where every member is a genuine")
    println("        bound state and the continuum is missing -- the guard MUST fire:")
    gridBad = Radial.Grid(Radial.Grid(true); rbox = 775.0)
    confsBad = [Configuration("1s")];   for n = 2:16    push!(confsBad, Configuration("$(n)p"))    end
    levBad  = perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=gridBad, nuclearModel=nm,
              configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]
    gmBad   = perform(Atomic.Computation(Atomic.Computation(); name="np, big box", grid=gridBad,
              nuclearModel=nm, configs=confsBad[2:end], asfSettings=asf); output=true)[ResultKeys.Multiplet]
    tBad    = DispersionCoefficient.multipoleTerms(E1, levBad, gmBad, gridBad)
    bad     = DispersionCoefficient.Partner("H(1s), 775 a.u. box", levBad, [tBad])
    DispersionCoefficient.checkSaturation(bad)
    aBad = DispersionCoefficient.staticPolarizability(tBad)
    cBad = DispersionCoefficient.computeAtomWallCoefficient(bad, Inf)
    @printf("   >> the damage it warns about:  alpha^(E1,0) = %.6f (%.5f of ref),  C_3 = %.6f (%.5f of exact)\n",
            aBad, aBad/ALPHA_REF[1], cBad, cBad/C3_EXACT)
end

main()
