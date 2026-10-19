#
# Validation of the DYNAMIC multipole polarizability alpha^(Mp,K)(omega) added to
# module-MultipolePolarizibility.jl on 19-Oct-2026 (step 2; priority item 69's remaining half).
# Five cases, each testing something the previous one cannot:
#
#   CASE A  THE STATIC LIMIT.  alpha(omega) through the DYNAMIC branch as omega -> 0 must return the static
#           value;  omega = 0 itself takes the static branch deliberately, to stay bit-identical.
#   CASE B  THE SMALL-omega LAW, an EXACT internal check with no literature in it.  Expanding the denominator,
#               alpha(omega)   = alpha(0) + omega^2 S3 + O(omega^4)
#               alpha(i omega) = alpha(0) - omega^2 S3 + O(omega^4),     S3 = sum_a A_a / dE_a^3
#           so the two axes must depart from the static value by the SAME amount in OPPOSITE directions, with a
#           slope S3 that is computable in closed form from the stored matrix elements.  This also checks that
#           `Contribution.reducedMe` round-trips, which step 3 depends on.
#   CASE C  SHAPE ON THE IMAGINARY AXIS: positive, monotone decreasing, below alpha(0), and -> 0.
#   CASE D  THE REAL-AXIS GUARDS fire on a pole and above the first excitation energy.
#   CASE E  END TO END: C_6(H-H) from alpha^(E1,0)(i omega), BOTH as the closed-form double sum and by numerical
#           quadrature over frequencyFactor, against the published 6.498392245754 (Dirac) / 6.499026705406
#           (non-relativistic), Tang, Zhang, Zhang, Jiang & Mitroy, arXiv:1205.6883 Table VIII.
#           This is the dry run for the DispersionCoefficient module (priority item 80).
#
# RESULTS, 19-Oct-2026, H 1s with np n = 2..16 on a 151 a.u. box, 30 perturbers, lowest dE = 0.375005 a.u.:
#
#   alpha^(E1,0)(0) = 4.4997326653 against the exact 4.5, and REBUILT FROM `Contribution.reducedMe` to
#   rel.diff 0.00e+00 -- so the stored matrix element round-trips exactly, which step 3 depends on.
#
#   CASE A  the dynamic branch returns the static value as omega -> 0, departing as omega^2 exactly:
#           5.91e-12, 5.91e-08, 5.91e-06 at omega = 1e-6, 1e-4, 1e-3.
#   CASE B  S3 = 2.658116e+01 in closed form -- quoted to six figures ON PURPOSE: two rebuilds of the same
#           case gave 2.658116118e+01 and 2.658115926e+01, a 7e-08 spread, which is this project's known
#           build sensitivity and not a property of S3.  The measured slopes are +2.658116e+01 (real) and
#           -2.658116e+01 (imaginary) -- THE SAME DEPARTURE IN OPPOSITE DIRECTIONS, as the expansion
#           demands -- to 7.3e-08 and 7.5e-08 at omega = 1e-4, and the residual itself scales as omega^2
#           (6.5e-06 at 1e-3, 6.5e-04 at 1e-2), i.e. it IS the O(omega^4) remainder and nothing else.
#   CASE C  alpha(i omega)/alpha(0) = 1, 0.944512, 0.661854, 0.164986, 0.023316, 0.002203, 0.000022 at
#           omega = 0, 0.1, 0.3, 1, 3, 10, 100;  monotone decreasing, and the 1/omega^2 tail is right
#           (9.98e-05 at omega = 100 against alpha(0)/omega^2 scaling).
#   CASE D  all three behave: ON A POLE fires at 0.50 % from the 2p perturber; at omega = 0.60 BOTH the
#           complex-alpha and the pole warning fire; and on the IMAGINARY axis at the same omega NOTHING
#           fires, which is the point -- there are no poles there.
#   CASE E  C_6(H-H) closed-form double sum  6.4982808622
#           C_6(H-H) numerical quadrature    6.4982808662   -- the two agree to 6.2e-10, which validates
#           `frequencyFactor` on the imaginary axis by two independent routes.
#           against the published Dirac 6.498392245754  ->  ratio 0.999983
#           against the non-relativistic 6.499026705406 ->  ratio 0.999885
#           So the WHOLE CHAIN -- reduced matrix elements, rank coefficient, imaginary-frequency factor,
#           Casimir-Polder double sum -- reproduces a published dispersion coefficient to 1.7e-05.
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra

const C6_DIRAC  = 6.498392245754
const C6_NONREL = 6.499026705406

# alpha^(E1,0) and the per-perturber pieces A_a, dE_a that every later case is built from
function pieces(level::Level, gMultiplet::Multiplet, grid::Radial.Grid, mp::EmMultipole)
    alphas, cs = MultipolePolarizibility.computeMultipolePolarizabilities(mp, level, gMultiplet, grid)
    a0 = 0.0;   for (K,a) in alphas    K == 0 && (a0 = a.Babushkin)   end
    As = Float64[];   dEs = Float64[]
    for  c in cs
        dE = -c.deltaE                                     # = W_J' - W_J
        # A_a from the STORED reduced matrix element, which is the quantity step 3 will re-sum
        A  = MultipolePolarizibility.rankCoefficient(mp.L, 0, level.J, c.nuLevel.J) * abs2(c.reducedMe)
        push!(As, A);   push!(dEs, dE)
    end
    return( (a0, As, dEs, cs) )
end

alphaFrom(As, dEs, omega, imag) =
    sum( As[i] * MultipolePolarizibility.frequencyFactor(dEs[i], omega, imag)  for i = 1:length(As) )

function main()
    nm   = Nuclear.Model(1., PointNucleus())
    grid = Radial.Grid(Radial.Grid(true); rbox = 129.0)
    asf  = AsfSettings(AsfSettings(); gridStopper=false)
    confs = Configuration[];   for n = 2:16    push!(confs, Configuration("$(n)p"))   end
    lev = perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=grid, nuclearModel=nm,
                  configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]
    gm  = perform(Atomic.Computation(Atomic.Computation(); name="np 2..16", grid=grid, nuclearModel=nm,
                  configs=confs, asfSettings=asf); output=true)[ResultKeys.Multiplet]

    a0, As, dEs, cs = pieces(lev, gm, grid, E1)
    # the SAME alpha_0 rebuilt from the stored matrix elements: if these disagree, reducedMe is not round-tripping
    a0rebuilt = alphaFrom(As, dEs, 0.0, false)
    println("\n", "="^118)
    @printf(" H 1s, np n = 2..16 on a %.1f a.u. box:  %d perturbers,  lowest dE = %.6f a.u.\n",
            grid.r[end], length(As), minimum(dEs))
    @printf(" alpha^(E1,0)(0) = %.10f   (exact 4.5)      rebuilt from reducedMe %.10f   rel.diff %.2e\n",
            a0, a0rebuilt, abs(a0rebuilt-a0)/a0)
    println("="^118)

    println("\n CASE A -- the static limit through the DYNAMIC branch")
    @printf("   %-14s %22s %22s %12s\n", "omega", "alpha(omega)", "alpha(i omega)", "rel to a0")
    for w in [1.0e-6, 1.0e-4, 1.0e-3]
        ar = alphaFrom(As, dEs, w, false);   ai = alphaFrom(As, dEs, w, true)
        @printf("   %-14.1e %22.14e %22.14e %12.2e\n", w, ar, ai, abs(ar-a0)/a0)
    end

    println("\n CASE B -- the small-omega law:  (alpha(w) - a0)/w^2 -> +S3  and  (alpha(iw) - a0)/w^2 -> -S3")
    S3 = sum( As[i] / dEs[i]^3  for i = 1:length(As) )
    @printf("   S3 in closed form = %.12e\n", S3)
    @printf("   %-10s %20s %20s %14s %14s\n", "omega", "slope real", "slope imaginary", "rel err real", "rel err imag")
    for w in [1.0e-4, 1.0e-3, 1.0e-2]
        sr = (alphaFrom(As, dEs, w, false) - a0)/w^2
        si = (alphaFrom(As, dEs, w, true)  - a0)/w^2
        @printf("   %-10.1e %20.12e %20.12e %14.2e %14.2e\n", w, sr, si, abs(sr-S3)/S3, abs(si+S3)/S3)
    end

    println("\n CASE C -- shape on the imaginary axis (no poles there at all)")
    @printf("   %-10s %22s %12s\n", "omega", "alpha(i omega)", "of alpha(0)")
    prev = a0;  mono = true
    for w in [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0]
        ai = alphaFrom(As, dEs, w, true)
        ai > prev + 1.0e-14 && (mono = false)
        @printf("   %-10.2f %22.14e %12.6f\n", w, ai, ai/a0)
        prev = ai
    end
    println("   >> monotone decreasing: $mono    (and alpha(i omega) must never exceed alpha(0))")

    println("\n CASE D -- the real-axis guards")
    println("   (D1) omega just below the 2p perturber -- ON A POLE must fire:")
    MultipolePolarizibility.computeMultipolePolarizabilities(E1, lev, gm, grid; omega = 0.995*minimum(dEs))
    println("\n   (D2) omega well above the lowest excitation energy -- COMPLEX-alpha warning must fire:")
    MultipolePolarizibility.computeMultipolePolarizabilities(E1, lev, gm, grid; omega = 0.60)
    println("\n   (D3) the imaginary axis at the same omega -- NOTHING may fire:")
    MultipolePolarizibility.computeMultipolePolarizabilities(E1, lev, gm, grid; omega = 0.60, imaginaryOmega = true)
    println("   (if no red text appeared between D3 and here, the guard is correctly silent)")

    println("\n CASE E -- END TO END: C_6(H-H) from alpha^(E1,0)(i omega)")
    # closed form:  int_0^inf [a/(a^2+w^2)][b/(b^2+w^2)] dw = pi/(2(a+b))  =>  C_6 = (3/2) sum_ab A_a A_b/(dE_a+dE_b)
    c6sum = 1.5 * sum( As[i]*As[j]/(dEs[i]+dEs[j])  for i = 1:length(As), j = 1:length(As) )
    # and by NUMERICAL quadrature over frequencyFactor, on the substitution w = t/(1-t) to reach infinity
    nq = 20000;  acc = 0.0
    for k = 1:nq
        tt = (k - 0.5)/nq;   w = tt/(1-tt);   jac = 1/(1-tt)^2
        acc = acc + alphaFrom(As, dEs, w, true)^2 * jac / nq
    end
    c6quad = 3.0/pi * acc
    @printf("   C_6 closed-form double sum  = %.10f\n", c6sum)
    @printf("   C_6 numerical quadrature    = %.10f   (agree to %.2e -- validates frequencyFactor on the i axis)\n",
            c6quad, abs(c6quad-c6sum)/c6sum)
    @printf("   published Dirac             = %.10f   ratio %.6f\n", C6_DIRAC,  c6sum/C6_DIRAC)
    @printf("   published non-relativistic  = %.10f   ratio %.6f\n", C6_NONREL, c6sum/C6_NONREL)
end

main()
