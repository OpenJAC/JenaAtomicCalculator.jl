#
# probe-rydbergDeliverable.jl  --  20-Oct-2026
#
# THE ONE QUESTION:  of the numbers the Rydberg / optical-tweezer community actually asks for, WHICH ONES can
# JAC supply, and what sets the error on each?  Three cases, each against a reference that owes nothing to JAC.
#
#   CASE 1  the CORE polarizability by the complete-basis chi route, on a box ladder -- is it box-independent,
#           where a sum over states is not?   Na+ against the RRPA 0.9457 a.u., Ne as the control.
#   CASE 2  a Rydberg QUANTUM DEFECT in a frozen ion core, with and without the core-polarization term
#           alpha_core/(2 r^3) -- how accurately must alpha_core be known?   Na against the measured defects.
#   CASE 3  a Rydberg-Rydberg radial DIPOLE MATRIX ELEMENT against the closed form
#           <n,l+1| r |n,l> = (3/2) n sqrt(n^2 - (l+1)^2),  which needs no literature.   H at n = 16..60.
#
# Run:  julia --project=. tools/probe-rydbergDeliverable.jl  [c1] [c2] [c3]      (no argument runs all three)
#
# ===================================================================================================
#  RESULTS, 20-Oct-2026.  JAC'S STRENGTH HERE IS THE MATRIX ELEMENT, NOT THE ENERGY -- AND THAT IS THE
#  OPPOSITE OF WHAT THE COMMUNITY'S OWN TOOLS ARE GOOD AT, WHICH IS WHY IT IS WORTH SAYING.
#
#  CASE 1  THE BOX STOPS BEING A TUNING KNOB.  alpha(Na+) over a FIVEFOLD box change:
#
#            rbox      40      80     120     200  a.u.
#            alpha  1.277625 1.277776 1.277831 1.277876      spread 2.0e-4 relative
#
#          Against a SUM OVER STATES for the same system (2p^5 ns + 2p^5 nd, n <= 8), which gives
#          1.1213 / 0.8719 / 0.6861 / 0.6135 at rbox 12 / 20 / 32 / 50 -- a monotone slide through the
#          reference, so the right answer is reachable only by CHOOSING the box.  For a closed-shell core
#          the continuum carries about 35 % (hydrogen 1s: 18.6 %), and the chi route includes it implicitly.
#
#          WHAT IS LEFT IS PHYSICS, NOT NUMERICS:  1.2778 is 1.351 x the RRPA 0.9457, and the Ne control
#          gives 3.7663 against the measured 2.669, i.e. 1.411 x.  Both are the UNSCREENED response.  The
#          missing induced field is priority item 75 (exchange in the induced potential), so THIS CASE IS
#          WHAT GIVES ITEM 75 AN OUTSIDE CONSUMER.
#
#  CASE 2  AND alpha_core IS NOT THE LIMITING ERROR -- THE MEAN FIELD IS.  Na quantum defects, mean over
#          n = 8..14, in a frozen Na+ DFS core:
#
#                                 delta(s)  delta(p)  delta(d)  delta(f)
#            bare frozen core      1.47052   0.99503   0.08961   0.00598
#            + pol, alpha=0.9457   1.49199   1.01684   0.10555   0.00754
#            + pol, alpha=1.2778   1.49947   1.02459   0.11156   0.00810
#            MEASURED              1.34796   0.85538   0.01554   0.00145
#
#          Read the rows in order.  A 35 % error in alpha_core (0.332 absolute) moves delta by 0.0075,
#          while the BARE core is already 0.123 away from the measured value -- SIXTEEN TIMES as much.  And
#          the polarization term moves delta the WRONG WAY here, because a frozen ION core without
#          correlation already over-binds the valence electron, so making it more attractive is a
#          correction applied to an error of the opposite sign.  (Note this is the mirror image of
#          example-Cj.jl branch d, where the NEUTRAL's self-consistent DFS UNDER-binds Li 2s: the frozen-ion
#          and self-consistent-neutral local fields fail in opposite directions.)
#
#          A Rydberg level moves by d(delta)/n^3, so 0.001 of delta is 30 MHz at n = 60.  The alpha_core
#          error is then ~230 MHz and the BARE MEAN FIELD's ~3.7 GHz -- and Foerster defects that matter for
#          blockade are tens of MHz to about a GHz.  Since C_6 ~ 1/Delta, an ALKALI blockade C_6 built on
#          JAC's OWN energies cannot be competitive with one built on MEASURED defects.  That is an honest
#          negative and it is the single most useful thing this probe establishes.
#
#  CASE 3  BUT GIVEN THE ENERGIES, THE MATRIX ELEMENT IS RIGHT.  H, one box of 10580 a.u. (= 2.5 r_+ for
#          n = 46) with hp = rbox/3000, 502 splines;  ratio to the exact closed form:
#
#               n      s_1/2 -> p_3/2     p_1/2 -> d_3/2
#              16        0.9999964          0.9999988
#              30        1.0004044          1.0000858
#              46        1.0005039          1.0010421
#              60        0.9244425          1.0005247
#
#          So 1e-6 at n = 16 and still 1e-3 at n = 46 -- for elements of ~3000 a.u., which is the regime
#          where a sum over states has nothing to say.  THE n = 60 s-CHANNEL FAILURE IS RULE 12 AND IT
#          ANNOUNCES ITSELF:  60s wants about 2.5 x 2 x 60^2 = 18000 a.u. and this box is 10580, and the
#          same orbital's ENERGY is simultaneously off (ratio 0.99918 against 0.9999975 for 60p_1/2).  So
#          the energy is a free self-check on whether the box still carries the matrix element.
#
#  WHAT THE THREE CASES SAY TOGETHER.  The community's standard tools take MEASURED quantum defects and
#  compute matrix elements from a one-electron model potential.  JAC is strong exactly where that is weak --
#  the radial integral, the response, the core -- and weak exactly where that is strong, the energy.  So the
#  useful contributions are (a) alpha_core, and especially alpha_core(omega) at trap frequencies, which is a
#  genuine INPUT to everyone's model potential and far less available dynamically than statically;  and
#  (b) systems where no measured defect exists to fit -- divalent (Sr, Yb) and ionic Rydberg series -- where
#  a model potential has no footing and a many-electron code does.  Competing on Rb nS defects does not.
#
#  TWO ROUTES MEASURED AND REJECTED -- DO NOT REBUILD THEM:
#   - a sum over states for a closed-shell CORE polarizability (CASE 1's second table);  the box becomes a
#     fitted parameter.
#   - assembling alpha from `RandomPhase.solve`'s coupled chi with either
#     (2/3) CL^2 rkDiagonal  or  (2/3)(2j_a+1) eMultipole.  Calibrated on H 1s with coupling ON (one active
#     electron has no OTHER electron, so the induced term vanishes and the answer must be the exact 4.5):
#     the two give -5.4395 and +7.3828.  The coupled path's chi carries `eMultipole`'s normalisation and a
#     chi-convention sign, and the VALIDATED coupled assembly already exists in
#     `tools/mResolvedPolarizability.jl` (`runCoupledCase`, which self-checks against an uncoupled control).
#     Use that one;  this probe therefore reports only the UNSCREENED total in CASE 1.
# ===================================================================================================

using JenaAtomicCalculator, Printf, LinearAlgebra
const DL = JenaAtomicCalculator.DalgarnoLewis
const RP = JenaAtomicCalculator.RandomPhase
const IS = InteractionStrength
const RI = RadialIntegrals
const AM = AngularMomentum

const dls = DL.Settings(DL.Settings(); printBefore=false)
const runAll = !any(a -> a in ARGS, ["c1", "c2", "c3"])


# Projects the occupied orbitals of the target channel out of chi, so that the Pauli principle is respected;
# copied from tools/mResolvedPolarizability.jl, where the whole assembly is validated on alpha(H 1s) = 4.5.
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


# The unscreened static dipole polarizability of a CLOSED-SHELL system, assembled as
# alpha = sum_a sum_t (2/3) CL(a,1,t)^2 <phi_a| r |chi_at>, with chi driven by the pure radial drive.
function coreAlpha(shells, orbs, pot, prim, grid)
    tot = 0.0
    for  aSub  in  shells
        (!haskey(orbs, aSub) || orbs[aSub].energy >= 0.0) && continue
        for  vSub  in  RP.reachableChannels(1, aSub)
            chi = redirect_stdout(devnull) do
                DL.solve(E1, Basics.Babushkin, 0.0, vSub, orbs[aSub], pot, prim; settings=dls,
                         staticDrive=true, extraSource = unitFix(aSub, vSub, orbs, prim, grid)) end
            chiP = projectPauli(chi, vSub, orbs, shells, grid)
            tot += 2.0/3.0 * AM.CL_reduced_me(aSub, 1, vSub)^2 * RI.rkDiagonal(1, orbs[aSub], chiP, grid)
        end
    end

    return( tot )
end


# Builds the frozen DFS core of a closed-shell system and the bound orbitals of `shells` within it.
function frozenCore(Z, confString, coreShells, shells, grid, prim)
    nm  = Nuclear.Model(Z, PointNucleus())
    asf = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
    mp  = redirect_stdout(devnull) do
        SelfConsistent.performSCF([Configuration(confString)], nm, grid, asf; printout=false) end
    potNuc = Nuclear.nuclearPotential(nm, grid)
    pot    = Basics.add( potNuc, Basics.computePotential(Basics.DFSField(), grid, mp.levels[1].basis) )
    orbs   = redirect_stdout(devnull) do
        Bsplines.generateOrbitals(shells, pot, nm, prim; printout=false) end

    return( nm, pot, orbs )
end


# ===================================================================================================
# CASE 1 -- is the chi route's core polarizability BOX-INDEPENDENT?
# ===================================================================================================
if  runAll || "c1" in ARGS
    for (label, Z, conf, ref, refName, boxes) in
            (("Na+ (Ne-like)", 11.0, "1s^2 2s^2 2p^6", 0.9457, "RRPA alpha(Na+)",      [40.0, 80.0, 120.0, 200.0]),
             ("Ne  (CONTROL)", 10.0, "1s^2 2s^2 2p^6", 2.669,  "measured alpha(Ne)",   [80.0, 120.0]))
        @printf("\n=== CASE 1:  %s, unscreened chi route;  reference %s = %.4f a.u. ===\n", label, refName, ref)
        @printf(" %8s %6s %7s %15s %12s %10s\n", "rbox", "nsL", "points", "eps(outer)", "alpha", "/ref")
        shells = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
        for rbox in boxes
            local grid, prim, nm, pot, orbs, al
            grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = rbox/120, rbox = rbox)
            prim = Bsplines.generatePrimitives(grid)
            try
                nm, pot, orbs = frozenCore(Z, conf, shells, shells, grid, prim)
                al = coreAlpha(shells, orbs, pot, prim, grid)
                @printf(" %8.1f %6d %7d %+15.8f %12.6f %10.5f\n", rbox, grid.nsL, grid.NoPoints,
                        orbs[shells[end]].energy, al, al/ref)
            catch e
                @printf(" %8.1f   RAISED  %s\n", rbox, string(typeof(e)))
            end
        end
    end
    println("\n>> If the alpha column is flat, the box is no longer a fitted parameter and what remains is")
    println(">> the SCREENING, i.e. priority item 75.")
end


# ===================================================================================================
# CASE 2 -- how accurately must alpha_core be known for a Rydberg quantum defect?
# ===================================================================================================
if  runAll || "c2" in ARGS
    measured = Dict(0 => 1.34796, 1 => 0.85538, 2 => 0.015543, 3 => 0.001453)
    labels   = Dict(0 => "s", 1 => "p", 2 => "d", 3 => "f")
    rbox = 1000.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = rbox/3000, rbox = rbox)
    prim = Bsplines.generatePrimitives(grid)
    vShells = Subshell[]
    for n = 8:14, kap in (-1, 1, 2, 3)   push!(vShells, Subshell(n, kap))   end
    coreShells = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
    nm, potCore, _ = frozenCore(11.0, "1s^2 2s^2 2p^6", coreShells, coreShells, grid, prim)
    @printf("\n=== CASE 2:  Na in a frozen Na+ DFS core;  r_max = %.1f, %d+%d splines, Z(r_max) = %.4f ===\n",
            grid.r[end], grid.nsL, grid.nsS, potCore.Zr[end])

    # Z(r) + alpha/(2 r^3) W(r),  W = 1 - exp(-(r/rc)^6):  the standard cut-off core-polarization term,
    # since Z(r) = -r V(r) and V_pol = -alpha/(2 r^4) W(r).
    function withPol(pot, alpha, rc)
        Zr = copy(pot.Zr)
        for i = 1:length(Zr)
            r = grid.r[i];   r <= 0.0 && continue
            Zr[i] += alpha/(2*r^3) * (1.0 - exp(-min(50.0, (r/rc)^6)))
        end
        return( Radial.Potential(@sprintf("core + pol(%.4f, %.2f)", alpha, rc), Zr, grid) )
    end

    function defects(pot, label)
        orbs = redirect_stdout(devnull) do
            Bsplines.generateOrbitals(vShells, pot, nm, prim; printout=false) end
        @printf("\n--- %s ---\n%4s %10s %10s %10s %10s\n", label, "n", "delta(s)", "delta(p)", "delta(d)", "delta(f)")
        acc = Dict(l => Float64[] for l = 0:3)
        for n = 8:14
            row = Float64[]
            for (l, kap) in ((0,-1), (1,1), (2,2), (3,3))
                s = Subshell(n, kap)
                if  haskey(orbs, s) && orbs[s].energy < 0.0
                    d = n - 1.0/sqrt(-2*orbs[s].energy);   push!(row, d);   push!(acc[l], d)
                else    push!(row, NaN)
                end
            end
            @printf("%4d %10.5f %10.5f %10.5f %10.5f\n", n, row[1], row[2], row[3], row[4])
        end
        out = Dict( l => isempty(acc[l]) ? NaN : sum(acc[l])/length(acc[l])  for l = 0:3 )
        @printf("%4s", "mean");   for l = 0:3   @printf(" %10.5f", out[l])               end
        @printf("\n%4s", "meas");  for l = 0:3   @printf(" %10.5f", measured[l])          end
        @printf("\n%4s", "diff");  for l = 0:3   @printf(" %10.5f", out[l]-measured[l])   end
        println()

        return( out )
    end

    d0 = defects(potCore,              "bare frozen DFS core, no core-polarization term")
    dA = defects(withPol(potCore, 0.9457, 1.5), "+ pol, alpha = 0.9457 (the RRPA reference), rc = 1.5")
    dB = defects(withPol(potCore, 1.2778, 1.5), "+ pol, alpha = 1.2778 (CASE 1's unscreened value), rc = 1.5")
    @printf("\n>> SENSITIVITY to d(alpha) = %+.4f, against the BARE core's own error:\n", 1.2778-0.9457)
    @printf("%6s %14s %16s %12s\n", "l", "d(delta)", "bare - measured", "ratio")
    for l = 0:3
        @printf("%6s %+14.5f %+16.5f %12.1f\n", labels[l], dB[l]-dA[l], d0[l]-measured[l],
                (d0[l]-measured[l])/(dB[l]-dA[l]))
    end
    # A Rydberg level sits at E = -1/(2(n-delta)^2), so dE = d(delta)/(n-delta)^3 ~ d(delta)/n^3.
    @printf("\n>> 0.001 of delta is %.3e Ha = %.1f MHz at n = 60, so the last column is how many times\n",
            0.001/60.0^3, 0.001/60.0^3 * 6.5796839e9)
    @printf(">> larger the mean-field error is than the alpha_core error it is being corrected with --\n")
    @printf(">> here %.2f GHz against %.2f GHz, where a blockade-relevant Foerster defect is tens of MHz\n",
            (d0[0]-measured[0])/60.0^3 * 6.5796839e6, (dB[0]-dA[0])/60.0^3 * 6.5796839e6)
    @printf(">> to about a GHz.\n")
end


# ===================================================================================================
# CASE 3 -- given the energies, is the Rydberg-Rydberg dipole matrix element right?
# ===================================================================================================
if  runAll || "c3" in ARGS
    rbox = 10580.0                     ## = 2.5 r_+ for n = 46;  see tools/probe-rydbergBoxConflict.jl
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = rbox/3000, rbox = rbox)
    prim = Bsplines.generatePrimitives(grid)
    nm   = Nuclear.Model(1.0, PointNucleus())
    pot  = Nuclear.nuclearPotential(nm, grid)
    ns   = [16, 30, 46, 60]
    shells = Subshell[]
    for n in ns, kap in (-1, 1, -2, 2, -3)   push!(shells, Subshell(n, kap))   end
    orbs = redirect_stdout(devnull) do
        Bsplines.generateOrbitals(shells, pot, nm, prim; printout=false) end
    @printf("\n=== CASE 3:  H, r_max = %.1f, %d points, %d+%d splines ===\n", grid.r[end], grid.NoPoints,
            grid.nsL, grid.nsS)
    @printf("\n%4s %-9s %-9s %16s %16s %12s\n", "n", "from", "to", "JAC <a|r|b>", "exact", "ratio")
    for n in ns
        for (l, kFrom, kTo, nameF, nameT) in ((0, -1, -2, "s_1/2", "p_3/2"), (1, 1, -3, "p_1/2", "d_3/2"))
            local a, b, jac, exact
            a = Subshell(n, kFrom);   b = Subshell(n, kTo)
            if  !(haskey(orbs,a) && haskey(orbs,b)) || orbs[a].energy >= 0.0 || orbs[b].energy >= 0.0
                @printf("%4d %-9s %-9s %16s\n", n, nameF, nameT, "MISSING");   continue
            end
            jac   = RI.rkDiagonal(1, orbs[a], orbs[b], grid)
            exact = 1.5 * n * sqrt(n^2 - (l+1)^2)
            @printf("%4d %-9s %-9s %16.4f %16.4f %12.7f\n", n, nameF, nameT, abs(jac), exact, abs(jac)/exact)
        end
    end
    @printf("\n%4s %-9s %18s %18s %12s\n", "n", "subshell", "JAC energy", "exact Dirac", "ratio")
    for n in ns
        for kap in (-1, 1, -2)
            local s, ex
            s = Subshell(n, kap)
            haskey(orbs, s) && orbs[s].energy < 0.0 || continue
            ex = Basics.computeDiracEnergy(s, 1.0)
            @printf("%4d %-9s %18.10f %18.10f %12.7f\n", n, string(s), orbs[s].energy, ex, orbs[s].energy/ex)
        end
    end
    println("\n>> The ENERGY column is the free self-check: where the box has stopped carrying an orbital, its")
    println(">> energy says so before the matrix element does.")
end
