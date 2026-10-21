#
# probe-dynamicCorePolarizability.jl  --  21-Oct-2026,  priority item 92 (was challenges S19)
#
# THE ONE QUESTION OF STEP 1:  does the COMPLETE-BASIS chi route, which is the only route whose core
# polarizability is box-independent, still work at omega != 0 -- and does its dynamic assembly reduce
# EXACTLY to the validated static number as omega -> 0?
#
# WHY THE QUESTION IS NOT TRIVIAL.  A static polarizability needs ONE solve; a dynamic one needs TWO,
# because the two denominators (E_n - E_a - omega) and (E_n - E_a + omega) no longer coincide:
#
#     alpha(omega) = (1/3) sum_a sum_t CL(a,1,t)^2 [ <phi_a|r|chi_at(+omega)> + <phi_a|r|chi_at(-omega)> ]
#
# `DalgarnoLewis.solve` already takes omega and solves at targetEnergy = E_a + omega, so nothing has to be
# built -- but the PREFACTOR SPLIT is new: at omega = 0 the two solves coincide and 1/3 + 1/3 must give back
# the static 2/3 to machine precision.  That is the calibration, and it is where such a thing goes wrong.
#
#   CASE A  THE omega = 0 REDUCTION, to machine precision, against the static assembly on the same grid.
#           H 1s, where the static chi route gives 4.49975 against the EXACT 4.5.
#   CASE B  THE QUADRATIC ONSET, an exact internal law with no literature in it:
#               alpha(omega) = alpha(0) + omega^2 S3 + O(omega^4),
#           so (alpha(omega) - alpha(0))/omega^2 must approach a CONSTANT, and the residual must itself fall
#           as omega^2.  This is the same law `probe-dynamicPolarizability.jl` CASE B used on the
#           sum-over-states route, applied here to a route that shares none of its code.
#   CASE C  AN INDEPENDENT CROSS-ROUTE CHECK.  On a compact box H 1s is reached by BOTH routes (chi: 4.49975;
#           sum-over-states with np n = 2..16: 4.4997), so at an omega far from every pole the two must agree
#           to about that same accuracy.  The routes share no code path -- a B-spline linear solve against a
#           stored perturber list -- so agreement is real evidence and disagreement localises the error.
#
# Run:  julia --project=. tools/probe-dynamicCorePolarizability.jl  [a] [b] [c]   (no argument runs all)
#
# ===================================================================================================

using JenaAtomicCalculator, Printf, LinearAlgebra
const DL = JenaAtomicCalculator.DalgarnoLewis
const RP = JenaAtomicCalculator.RandomPhase
const IS = InteractionStrength
const RI = RadialIntegrals
const AM = AngularMomentum

const dls    = DL.Settings(DL.Settings(); printBefore=false)
const runAll = !any(a -> a in ARGS, ["a", "b", "c"])

# 1064 nm, the workhorse trap wavelength, in a.u. -- the frequency this item is ultimately about.
const OMEGA_1064 = 2.0 * pi * 137.035999 / (1064.0e-9 / 5.29177210903e-11)


# Projects the occupied orbitals of the target channel out of chi, so the Pauli principle is respected.
# Copied unchanged from tools/probe-rydbergDeliverable.jl, where the assembly is validated statically.
function projectPauli(chi, t, orbs, shells, grid)
    P = copy(chi.P);   Q = copy(chi.Q)
    for  o  in  shells
        (o.kappa == t.kappa && haskey(orbs, o) && orbs[o].energy < 0.0)  ||  continue
        ov = RI.overlap(orbs[o].P, P, grid) + RI.overlap(orbs[o].Q, Q, grid)
        nn = RI.overlap(orbs[o].P, orbs[o].P, grid) + RI.overlap(orbs[o].Q, orbs[o].Q, grid)
        nn <= 0.0 && continue
        for i = 1:min(length(P), length(orbs[o].P))   P[i] -= ov/nn * orbs[o].P[i]   end
        for i = 1:min(length(Q), length(orbs[o].Q))   Q[i] -= ov/nn * orbs[o].Q[i]   end
    end
    return( Orbital(chi.subshell, false, true, chi.targetEnergy, P, Q,
                    zeros(length(P)), zeros(length(Q)), Radial.Grid()) )
end


# The extraSource that turns DalgarnoLewis.solve's own drive into the PURE radial drive r, i.e. R_i - b_i.
function unitFix(aSub, vSub, orbs, prim, grid)
    nsL = grid.nsL;  nsS = grid.nsS;  R = zeros(nsL+nsS);  b = zeros(nsL+nsS)
    for i = 1:nsL
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesL[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, full, zeros(grid.NoPoints),
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[i] = RI.rkDiagonal(1, orbs[aSub], syn, grid);   b[i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    for i = 1:nsS
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesS[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, zeros(grid.NoPoints), full,
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[nsL+i] = RI.rkDiagonal(1, orbs[aSub], syn, grid)
        b[nsL+i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    return( R - b )
end


"""
One channel's contraction <phi_a| r |chi_at(omega)>, with chi solved at targetEnergy = E_a + omega.
`omega` may be of either sign;  that is the whole of what step 1 adds.
"""
function channelTerm(aSub, vSub, omega, orbs, pot, prim, grid, shells)
    chi = redirect_stdout(devnull) do
        DL.solve(E1, Basics.Babushkin, omega, vSub, orbs[aSub], pot, prim; settings=dls,
                 staticDrive=true, extraSource = unitFix(aSub, vSub, orbs, prim, grid)) end
    chiP = projectPauli(chi, vSub, orbs, shells, grid)

    return( RI.rkDiagonal(1, orbs[aSub], chiP, grid) )
end


"""
THE OCCUPATION FACTOR, and it is why hydrogen is not simply the formula's own output.  Both assemblies below
carry N_a/(2j_a+1) per subshell -- exactly 1 for every CLOSED subshell, so it is invisible for the ion cores
this item is about, and 1/2 for hydrogen's single 1s electron.  The assembly copied from
tools/probe-rydbergDeliverable.jl has that factor set to 1 implicitly, because every system it was validated on
was closed-shell;  applied to H 1s unchanged it returns 8.99950, which is the 1s^2 value and NOT a defect.
Measured 21-Oct-2026: 8.999502620270 = 2 x 4.499751310135.

The STATIC assembly exactly as validated (one solve per channel, prefactor 2/3).  Kept separate from the
dynamic one on purpose: CASE A compares them, so they must not share the code being tested.
"""
function alphaStatic(shells, orbs, pot, prim, grid, occ)
    tot = 0.0
    for  aSub  in  shells
        (!haskey(orbs, aSub) || orbs[aSub].energy >= 0.0) && continue
        fOcc = occ[aSub] / (Basics.subshell_2j(aSub) + 1)
        for  vSub  in  RP.reachableChannels(1, aSub)
            tot += fOcc * 2.0/3.0 * AM.CL_reduced_me(aSub, 1, vSub)^2 *
                   channelTerm(aSub, vSub, 0.0, orbs, pot, prim, grid, shells)
        end
    end

    return( tot )
end


"""
The DYNAMIC assembly: two solves per channel, at +omega and -omega, each carrying 1/3.  At omega = 0 the two
solves are identical and this must return `alphaStatic` to machine precision -- CASE A.
"""
function alphaDynamic(shells, orbs, pot, prim, grid, omega, occ)
    tot = 0.0
    for  aSub  in  shells
        (!haskey(orbs, aSub) || orbs[aSub].energy >= 0.0) && continue
        fOcc = occ[aSub] / (Basics.subshell_2j(aSub) + 1)
        for  vSub  in  RP.reachableChannels(1, aSub)
            cl2  = fOcc * AM.CL_reduced_me(aSub, 1, vSub)^2
            tot += 1.0/3.0 * cl2 * channelTerm(aSub, vSub, +omega, orbs, pot, prim, grid, shells)
            tot += 1.0/3.0 * cl2 * channelTerm(aSub, vSub, -omega, orbs, pot, prim, grid, shells)
        end
    end

    return( tot )
end


# Hydrogen on one box, with the 1s orbital and the two p channels it reaches.
function hydrogen(rbox)
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = rbox/120, rbox = rbox)
    prim = Bsplines.generatePrimitives(grid)
    nm   = Nuclear.Model(1.0, PointNucleus())
    pot  = Nuclear.nuclearPotential(nm, grid)
    shells = [Subshell("1s_1/2")]
    orbs = redirect_stdout(devnull) do
        Bsplines.generateOrbitals(shells, pot, nm, prim; printout=false) end
    occ  = Dict(Subshell("1s_1/2") => 1.0)          ## ONE electron, not a closed 1s^2

    return( grid, prim, pot, orbs, shells, occ )
end


println("\n", "="^110)
@printf("ITEM 92, STEP 1 [DYN-CHI-CALIB].  omega(1064 nm) = %.8f a.u.\n", OMEGA_1064)
println("="^110)

grid, prim, pot, orbs, shells, occ = hydrogen(80.0)
@printf("\nH 1s on rbox = 80 a.u.,  %d+%d splines,  eps(1s) = %+.10f  (exact -0.5000066565)\n",
        grid.nsL, grid.nsS, orbs[Subshell("1s_1/2")].energy)

# ---------------------------------------------------------------------------------------------------
if  runAll || "a" in ARGS
    println("\n CASE A -- the omega = 0 reduction:  the DYNAMIC assembly must return the STATIC one exactly")
    aS = alphaStatic( shells, orbs, pot, prim, grid, occ)
    aD = alphaDynamic(shells, orbs, pot, prim, grid, 0.0, occ)
    @printf("   static  assembly (2/3, one solve)    alpha = %.12f   against the exact 4.5:  %.2e\n",
            aS, abs(aS/4.5 - 1))
    @printf("   dynamic assembly (1/3+1/3, two)      alpha = %.12f\n", aD)
    @printf("   difference                                   %.3e  %s\n", abs(aD - aS),
            abs(aD - aS) < 1.0e-12 ? "<-- the prefactor split is RIGHT" :
                                     "<-- THE SPLIT IS WRONG;  stop and fix it before CASE B")
end

# ---------------------------------------------------------------------------------------------------
if  runAll || "b" in ARGS
    println("\n CASE B -- the quadratic onset:  (alpha(w) - alpha(0))/w^2 -> S3, and the residual falls as w^2")
    a0 = alphaDynamic(shells, orbs, pot, prim, grid, 0.0, occ)
    @printf("   %12s %18s %16s %14s\n", "omega", "alpha(omega)", "(a(w)-a0)/w^2", "alpha/alpha(0)")
    prev = nothing
    for w in [1.0e-4, 1.0e-3, 1.0e-2, 3.0e-2, OMEGA_1064, 1.0e-1]
        aw = alphaDynamic(shells, orbs, pot, prim, grid, w, occ)
        @printf("   %12.8f %18.10f %16.6f %14.8f%s\n", w, aw, (aw - a0)/w^2, aw/a0,
                abs(w - OMEGA_1064) < 1e-12 ? "   <-- 1064 nm" : "")
    end
    println("   >> a FLAT third column near small omega is the law;  its drift at larger omega is O(w^4).")
end

# ---------------------------------------------------------------------------------------------------
if  runAll || "c" in ARGS
    println("\n CASE C -- cross-route:  S3 from this route against S3 from the SUM-OVER-STATES route")
    a0 = alphaDynamic(shells, orbs, pot, prim, grid, 0.0, occ)
    s1 = (alphaDynamic(shells, orbs, pot, prim, grid, 1.0e-4, occ) - a0) / 1.0e-8
    s2 = (alphaDynamic(shells, orbs, pot, prim, grid, 1.0e-3, occ) - a0) / 1.0e-6
    s0 = s1 - (s2 - s1)/99.0                      ## Richardson, since the residual is O(omega^2)
    ref = 26.581161                               ## tools/probe-dynamicPolarizability.jl CASE B, 19-Oct-2026
    @printf("   chi route (complete basis, B-spline solve)   S3 = %.6f   [extrapolated omega -> 0]\n", s0)
    @printf("   sum over states (closed form, stored me's)   S3 = %.6f   [recorded 19-Oct-2026]\n", ref)
    @printf("   relative difference                               %.2e\n", abs(s0/ref - 1))
    println("   >> the two routes share NO code path -- a linear solve against a stored perturber list -- so")
    println("   >> this is independent evidence and not a self-check.")
end


# ===================================================================================================
# CASE D -- STEP 2 [DYN-CORE-BOX].  The item's own question: is alpha_core(omega) box-independent at a
# TRAP frequency, as the static one is (2.0e-4 over a fivefold box change)?  One variable: the box.
# ===================================================================================================
if  "d" in ARGS
    for (label, Z, conf) in (("Na+ (Ne-like)", 11.0, "1s^2 2s^2 2p^6"),)
        @printf("\n=== CASE D:  %s, chi route;  static reference RRPA 0.9457 a.u. ===\n", label)
        @printf(" %7s %5s %13s %13s %11s %11s\n", "rbox", "nsL", "alpha(0)", "alpha(1064nm)", "enhance %", "/RRPA")
        cShells = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
        cOcc    = Dict(Subshell("1s_1/2") => 2.0, Subshell("2s_1/2") => 2.0,
                       Subshell("2p_1/2") => 2.0, Subshell("2p_3/2") => 4.0)
        for rbox in [40.0, 80.0, 120.0, 200.0]
            local g, pr, nm, asf, mp, potNuc, pot, os, a0, aw
            try
                g  = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = rbox/120, rbox = rbox)
                pr = Bsplines.generatePrimitives(g)
                nm = Nuclear.Model(Z, PointNucleus())
                asf = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
                mp = redirect_stdout(devnull) do
                    SelfConsistent.performSCF([Configuration(conf)], nm, g, asf; printout=false) end
                potNuc = Nuclear.nuclearPotential(nm, g)
                pot    = Basics.add(potNuc, Basics.computePotential(Basics.DFSField(), g, mp.levels[1].basis))
                os = redirect_stdout(devnull) do
                    Bsplines.generateOrbitals(cShells, pot, nm, pr; printout=false) end
                a0 = alphaDynamic(cShells, os, pot, pr, g, 0.0,        cOcc)
                aw = alphaDynamic(cShells, os, pot, pr, g, OMEGA_1064, cOcc)
                @printf(" %7.1f %5d %13.6f %13.6f %11.4f %11.5f\n",
                        rbox, g.nsL, a0, aw, 100*(aw/a0 - 1), a0/0.9457)
            catch e
                @printf(" %7.1f   RAISED  %s\n", rbox, string(typeof(e)))
            end
            flush(stdout)
        end
    end
    println("\n>> TWO questions, and they are independent: is the alpha(1064nm) COLUMN flat (box-independence")
    println(">> survives at omega != 0), and is the ENHANCEMENT small (a trap paper may neglect the frequency")
    println(">> dependence of the core term)?  The /RRPA column is item 75/S30's exchange, not this item's.")
end

# ===================================================================================================
# CASE E -- is 0.067 % a NUMBER or a LAW?  The same measurement across the Ne-like sequence and one
# heavier shell.  If the enhancement is S3 omega^2 / alpha(0) with S3/alpha(0) = 1/dE_eff^2, then the
# extracted dE_eff must track the system's real excitation energy and rise with the charge.
# One box (120 a.u.), since CASE D settled the box question.
# ===================================================================================================
if  "e" in ARGS
    @printf("\n=== CASE E:  the enhancement across species, rbox = 120 a.u., omega(1064nm) = %.6f ===\n", OMEGA_1064)
    @printf(" %-16s %8s %13s %13s %11s %12s\n", "system", "Z", "alpha(0)", "alpha(1064nm)", "enhance %", "dE_eff a.u.")
    for (label, Z, conf, shellStr, occs) in
            (("Ne   (neutral)",  10.0, "1s^2 2s^2 2p^6", ["1s_1/2","2s_1/2","2p_1/2","2p_3/2"], [2.,2.,2.,4.]),
             ("Na+  (Ne-like)",  11.0, "1s^2 2s^2 2p^6", ["1s_1/2","2s_1/2","2p_1/2","2p_3/2"], [2.,2.,2.,4.]),
             ("Mg2+ (Ne-like)",  12.0, "1s^2 2s^2 2p^6", ["1s_1/2","2s_1/2","2p_1/2","2p_3/2"], [2.,2.,2.,4.]),
             ("Ar   (neutral)",  18.0, "1s^2 2s^2 2p^6 3s^2 3p^6",
              ["1s_1/2","2s_1/2","2p_1/2","2p_3/2","3s_1/2","3p_1/2","3p_3/2"], [2.,2.,2.,4.,2.,2.,4.]),
             ("K+   (Ar-like)",  19.0, "1s^2 2s^2 2p^6 3s^2 3p^6",
              ["1s_1/2","2s_1/2","2p_1/2","2p_3/2","3s_1/2","3p_1/2","3p_3/2"], [2.,2.,2.,4.,2.,2.,4.]))
        local g, pr, nm, asf, mp, pot, os, a0, aw, sh, oc
        try
            sh = [Subshell(x) for x in shellStr];   oc = Dict(zip(sh, occs))
            g  = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = 1.0, rbox = 120.0)
            pr = Bsplines.generatePrimitives(g)
            nm = Nuclear.Model(Z, PointNucleus())
            asf = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
            mp = redirect_stdout(devnull) do
                SelfConsistent.performSCF([Configuration(conf)], nm, g, asf; printout=false) end
            pot = Basics.add(Nuclear.nuclearPotential(nm, g),
                             Basics.computePotential(Basics.DFSField(), g, mp.levels[1].basis))
            os = redirect_stdout(devnull) do
                Bsplines.generateOrbitals(sh, pot, nm, pr; printout=false) end
            a0 = alphaDynamic(sh, os, pot, pr, g, 0.0,        oc)
            aw = alphaDynamic(sh, os, pot, pr, g, OMEGA_1064, oc)
            @printf(" %-16s %8.0f %13.6f %13.6f %11.4f %12.4f\n", label, Z, a0, aw,
                    100*(aw/a0 - 1), OMEGA_1064 / sqrt(max(aw/a0 - 1, 1.0e-30)))
        catch e
            @printf(" %-16s %8.0f   RAISED  %s\n", label, Z, string(typeof(e)))
        end
        flush(stdout)
    end
    println("\n>> dE_eff = omega / sqrt(enhancement) is the single oscillator the response behaves as.  If it")
    println(">> RISES along Ne -> Na+ -> Mg2+ and is SMALLER for the looser Ar/K+ 3p shell, the 0.067 % is a law")
    println(">> with a known scaling and not a coincidence of one system.")
end

println()
