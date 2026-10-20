#
# Validation of the FREQUENCY-DEPENDENT atom-wall coefficient added to module-DispersionCoefficient.jl,
# 20-Oct-2026.  Five cases.  The question the whole file answers:  the module shipped a wall of CONSTANT
# dielectric constant, scaling C_3 by (eps-1)/(eps+1).  Is that a model of a dielectric wall, or a LIMIT of
# one -- and if it is a limit, how far from it does a real wall sit?
#
#   CASE W1  THE CONTROL.  The new three-argument method with wallResonance = Inf must reproduce the
#            two-argument one BIT FOR BIT, including the committed C_3(H 1s) = 0.2499886321 and the eps = 2
#            value 0.0833295440.  Nothing that was published may move.
#   CASE W2  THE THREE LIMITS, each a physical wall:  Omega >> dE returns the constant-eps factor;
#            eps(0) -> Inf returns a perfect conductor; omega_0 -> 0 returns ZERO, because a wall whose
#            charges cannot follow the atom does not reflect at the frequencies that carry C_3.
#   CASE W3  THE EXACT PER-POLE FACTOR, on a hand-built ONE-POLE atom that needs no SCF at all:  the
#            closed form says C_3(omega_0) = C_3(const eps) / (1 + dE/Omega) exactly, so this is an
#            algebra check at machine precision and not a convergence statement.
#   CASE W4  A REAL WALL.  H(1s) against a silica-like oscillator, and a scan in omega_0, to put a number
#            on what the constant-eps shortcut costs where it matters.
#   CASE W5  THE GUARDS.  eps(0) < 1 and omega_0 <= 0 are not dielectrics and must RAISE.
#
# RESULTS, 20-Oct-2026.  All five cases pass, and the answer to the question above is that the shipped
# constant-eps wall is a LIMIT and that a real wall sits FAR from it.
#
#   W1  0.0e+00 on all three walls -- the published numbers are untouched; conductor/exact 1/4 = 0.999955,
#       which is the perturber set (np n = 2..16, 64 a.u. box) and not the new code.
#   W2  Omega >> dE agrees with the constant-eps factor to 3-4e-07 at omega_0 = 1e6;  eps(0) = 1e8 reaches
#       the conductor to 1.8e-04;  and omega_0 -> 0 falls LINEARLY in omega_0 (5.02e-03, 5.21e-05, 5.21e-07
#       for omega_0 = 1e-02, 1e-04, 1e-06), i.e. to zero as it must.
#   W3  worst relative difference 3.24e-16 over 12 (eps(0), omega_0) pairs.  An algebra check, so machine
#       precision is the only acceptable result;  the silica-like row gives Omega = 0.575228, f = 1.651915.
#   W4  eps(0) = 3.9, C_3(omega_0)/C_3(const eps) against omega_0:  0.2034 (2 eV), 0.3862 (5), 0.5532 (10),
#       0.7086 (20), 0.8560 (50), 0.9588 (200 eV).  AT THE SILICA-LIKE 10 eV THE CONSTANT-eps NUMBER IS
#       1.808x TOO LARGE -- 81 %, and that is the set-averaged figure now quoted in the module docstring;
#       the 65 % of the single 2p pole (f = 1.6519) is the per-pole value and is the smaller of the two,
#       since the higher np poles carry larger dE_a.
#   W5  all three guards raise.
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra

const C3_EXACT = 0.25            # <r^2>/12 with <r^2> = 3 for H 1s;  needs no literature
const EV       = 1.0 / 27.211386245988

function terms1s(nmax::Int64, nm::Nuclear.Model)
    confs = [Configuration("$(n)p")  for n = 2:nmax]
    grid  = Radial.Grid(Radial.Grid(true); rbox = 64.0)
    asf   = AsfSettings(AsfSettings(); gridStopper=false)
    lev   = perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=grid, nuclearModel=nm,
                    configs=[Configuration("1s")], asfSettings=asf); output=true)[ResultKeys.Multiplet].levels[1]
    gm    = perform(Atomic.Computation(Atomic.Computation(); name="np series", grid=grid, nuclearModel=nm,
                    configs=confs, asfSettings=asf); output=true)[ResultKeys.Multiplet]
    return( (DispersionCoefficient.Partner("H(1s)", lev, [DispersionCoefficient.multipoleTerms(E1, lev, gm, grid)]),
             lev, grid) )
end

function main()
    nm = Nuclear.Model(1., PointNucleus())

    ## CASE W3 and W5 first: they need no atomic structure whatever
    println("\n", "="^120)
    println(" CASE W3 -- the exact per-pole factor, on a ONE-POLE atom (A = 2, dE = 0.375):  the closed form")
    println("            says C_3(omega_0) = C_3(const eps)/(1 + dE/Omega) with no approximation anywhere")
    println("="^120)
    onePole = DispersionCoefficient.Partner("one pole", Level(),
                  [DispersionCoefficient.MultipoleTerms(E1, [2.0], [0.375])])
    @printf("   %-10s %-10s %12s %18s %18s %18s %12s\n", "eps(0)", "omega_0", "Omega", "C_3(omega_0)",
            "C_3(const eps)/f", "f = 1+dE/Omega", "rel.diff")
    worst = 0.0
    for  eps0 in [1.5, 3.9, 12.0],  om0 in [0.1, 0.3675, 1.0, 10.0]
        R, Om = DispersionCoefficient.wallReflectionPole(eps0, om0)
        c3w   = DispersionCoefficient.computeAtomWallCoefficient(onePole, eps0, om0)
        fac   = 1.0 + 0.375/Om
        c3c   = DispersionCoefficient.computeAtomWallCoefficient(onePole, eps0, Inf) / fac
        rd    = abs(c3w - c3c) / abs(c3c);   worst = max(worst, rd)
        @printf("   %-10.4f %-10.4f %12.6f %18.12f %18.12f %18.12f %12.2e\n", eps0, om0, Om, c3w, c3c, fac, rd)
    end
    @printf("   >> worst relative difference over 12 (eps,omega_0) pairs: %.2e  -- an algebra check, so this\n", worst)
    println("   >> must be machine precision and nothing else.")

    println("\n", "="^120)
    println(" CASE W5 -- the guards:  eps(0) < 1 and omega_0 <= 0 are not dielectrics")
    println("="^120)
    for  (eps0, om0) in [(0.5, 0.4), (3.9, 0.0), (3.9, -1.0)]
        try
            DispersionCoefficient.wallReflectionPole(eps0, om0)
            println("   eps(0) = $eps0, omega_0 = $om0  ->  NO ERROR RAISED  <<< DEFECT")
        catch e
            println("   eps(0) = $eps0, omega_0 = $om0  ->  raised as it should")
        end
    end

    ## now the atom
    println("\n", "="^120)
    println(" building H(1s) with np n = 2..16 on a 64 a.u. box (the plateau; see probe-dispersionCoefficient.jl)")
    println("="^120)
    hydrogen, lev, grid = terms1s(16, nm)
    terms = DispersionCoefficient.selectTerms(hydrogen, 1)
    @printf("   poles %d,  sum_a A_a = %.10f (exact 2),  alpha^(E1,0) = %.10f\n", length(terms.strengths),
            sum(terms.strengths), DispersionCoefficient.staticPolarizability(terms))

    println("\n", "="^120)
    println(" CASE W1 -- THE CONTROL:  wallResonance = Inf must reproduce the two-argument method bit for bit")
    println("="^120)
    @printf("   %-16s %22s %22s %14s\n", "wall", "2-argument", "3-arg, omega_0 = Inf", "difference")
    for  (lab, eps0) in [("conductor", Inf), ("eps = 2", 2.0), ("eps = 3.9", 3.9)]
        c2 = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, eps0)
        c3 = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, eps0, Inf)
        @printf("   %-16s %22.12f %22.12f %14.1e\n", lab, c2, c3, abs(c2-c3))
    end
    @printf("   >> conductor against the EXACT 1/4:  ratio %.6f\n",
            DispersionCoefficient.computeAtomWallCoefficient(hydrogen, Inf) / C3_EXACT)

    println("\n", "="^120)
    println(" CASE W2 -- the three limits")
    println("="^120)
    c3cond = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, Inf)
    for  eps0 in [2.0, 3.9]
        c3const = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, eps0)
        c3big   = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, eps0, 1.0e6)
        println(@sprintf("   Omega >> dE:   eps(0) = %-6.2f  omega_0 = 1e6  ->  %18.12f  against the constant-eps ",
                         eps0, c3big) * @sprintf("%18.12f   rel.diff %.2e", c3const, abs(c3big-c3const)/c3const))
    end
    c3hiEps = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, 1.0e8, 0.4)
    println(@sprintf("   eps(0) -> Inf: eps(0) = 1e8   omega_0 = 0.4  ->  %18.12f  against the conductor ", c3hiEps) *
            @sprintf("%18.12f   rel.diff %.2e", c3cond, abs(c3hiEps-c3cond)/c3cond))
    for  om0 in [1.0e-2, 1.0e-4, 1.0e-6]
        @printf("   omega_0 -> 0:  eps(0) = 3.90   omega_0 = %.0e  ->  %18.12e   (must fall to zero like omega_0)\n",
                om0, DispersionCoefficient.computeAtomWallCoefficient(hydrogen, 3.9, om0))
    end

    println("\n", "="^120)
    println(" CASE W4 -- a REAL wall:  what the constant-eps shortcut costs.  eps(0) = 3.9 throughout;")
    println("            omega_0 = 10 eV is a silica-like UV resonance, and H(1s) -> 2p sits at dE = 0.375 a.u.")
    println("="^120)
    c3const = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, 3.9)
    @printf("   %-14s %12s %20s %20s %14s\n", "omega_0 [eV]", "Omega [a.u.]", "C_3(omega_0)", "C_3(const eps)",
            "ratio")
    for  omEv in [2.0, 5.0, 10.0, 20.0, 50.0, 200.0]
        om0   = omEv * EV
        _, Om = DispersionCoefficient.wallReflectionPole(3.9, om0)
        c3w   = DispersionCoefficient.computeAtomWallCoefficient(hydrogen, 3.9, om0)
        @printf("   %-14.1f %12.6f %20.12f %20.12f %14.6f\n", omEv, Om, c3w, c3const, c3w/c3const)
    end
    println("   >> the ratio is the per-pole factor 1/(1 + dE_a/Omega) averaged over the set, and it is the")
    println("   >> answer to the question this file asks:  the constant-eps number is a LIMIT, and a real")
    println("   >> wall sits far from it.")

    println("\n   >> end of probe-dielectricWall.jl")
end

main()
