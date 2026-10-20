#
# probe-rydbergBoxConflict.jl  --  19-Oct-2026
#
# THE ONE QUESTION:  for a magic-wavelength calculation between a COMPACT state and a RYDBERG state, is the
# radial-box conflict a matter of PRINCIPLE or of COST?  Three cases, each with an EXACT reference.
#
#   CASE 0  how the B-spline count depends on (h, hp, rbox).  No atomic input at all.
#   CASE 1  on a ladder of boxes x spline counts: the orbital energies of 1s, 16s and 46s against
#           Basics.computeDiracEnergy (closed form), and alpha(H 1s) from np n = 2..16 against the EXACT 4.5.
#   CASE 2  on a wide box with the basis held fixed: does extending n_max restore alpha?
#
# ===================================================================================================
#  RESULTS, 19-Oct-2026.  THE BOX WAS NEVER THE OBSTACLE;  THE PERTURBER SET IS.
#
#  CASE 0  THE SPLINE COUNT FOLLOWS hp, NOT rbox.  hp = rbox/300 -- which is what Basics.recommendedGrid
#          uses -- pins nsL near 110 at EVERY box from 64 to 10580 a.u., so that recipe spends the same
#          ~110 splines over 165x more space.  hp = rbox/3000 buys ~500 at any box for ~5x the points.
#
#  CASE 1  exact Dirac:  1s = -0.5000066566   16s = -0.0019531312   46s = -0.0002362952 Ha
#
#            rbox  hp/rbox  nsL |  E(1s) err   E(16s)/ex    E(46s)/ex |  alpha(1s)     /4.5
#              64    1/300   99 |   1.3e-12   -115.85     -10445.24   |  4.3077437   0.95728
#              64   1/3000  485 |   1.3e-12   -117.80     -10343.53   |  4.3116131   0.95814
#             775    1/300  106 |   1.8e-10      1.000000    -45.57   |  3.6514892   0.81144
#             775   1/3000  492 |   1.3e-12      1.000000    -46.03   |  3.6514938   0.81144
#            3000   1/3000  496 |   1.4e-11      1.000000      0.5005 |  3.6514940   0.81144
#           10580    1/300  113 |   2.5e-09      0.877699      0.5611 |  3.6621223   0.81380
#           10580   1/1000  213 |   1.4e-09      1.113205      0.9914 |  3.6489926   0.81089
#           10580   1/3000  499 |   3.5e-10      1.000001      1.000001|  3.6514926   0.81144
#
#          (a) THE BOX BLOCKER IS GONE.  rbox = 10580 (= 2.5 r_+ for 46s) with hp = rbox/3000 gives 46s to
#              1.000001 AND 1s to 3.5e-10 ON THE SAME GRID.  A compact and a Rydberg orbital coexist.
#          (b) Rule 12's STARVATION is real and visible here -- 46s/exact runs 0.561 / 0.991 / 1.000001 at
#              113 / 213 / 499 splines on the 10580 box -- but it is a DIFFERENT failure from (c).
#          (c) alpha(1s) = 3.6515 on boxes 775, 3000 AND 10580 and at 106, 206, 210, 492, 496 and 499
#              splines -- IDENTICAL TO FIVE DIGITS across a five-fold change in resolution.  So it is not
#              a basis problem.  **IT IS THE EXACT BOUND-STATE LIMIT:** the hand sum of the exact
#              hydrogenic 1s->np oscillator strengths over the SAME n = 2..16 gives 3.6516.  The complete
#              bound limit is 3.6634 and THE CONTINUUM CARRIES 0.8366 OF THE EXACT 4.5, i.e. 18.6 %.
#              The narrow 64 a.u. box gets the TOTAL right only by mis-representing its upper members as
#              pseudo-continuum states that stand in for that 18.6 %.
#
#  CASE 2  MORE BOUND STATES CANNOT REPLACE THE CONTINUUM.  dE_max stays pinned at the 0.5 Ha ionization
#          threshold, so every added configuration is one more bound Rydberg state accumulating there with
#          an n^-3 oscillator strength:
#
#            box 3000 a.u.:  n_max 16 -> 90 (30 -> 178 poles)  alpha/4.5 0.81144 -> 0.81782   dE_max 0.4981 -> 0.5027
#            box  775 a.u.:  n_max 16 -> 90 (30 -> 178 poles)  alpha/4.5 0.81144 -> 0.87753   dE_max 0.4981 -> 0.5582
#
#          Six times the configurations buy 0.6 % on the 3000 a.u. box while 18.6 % is missing.  The 775
#          a.u. row improves only because dE_max finally creeps PAST threshold and a few pseudo-continuum
#          states appear -- the narrow box's mechanism, arriving late.
#
#  CONCLUSION:  one basis cannot serve both, AND IT DOES NOT NEED TO.  A magic wavelength is a crossing of
#  TWO INDEPENDENT polarizabilities and nothing requires them to share a grid.  Moreover the RYDBERG side
#  needs no sum at all at optical frequency: 46S -> 46P sits at 67.6 GHz = 1.03e-05 a.u. while 1064 nm is
#  omega = 4.28e-02, a factor 4168 ABOVE every nearby-n spacing, so the response goes to the ponderomotive
#  limit alpha -> -1/omega^2 = -545.3 a.u., NEGATIVE -- which is why a red-detuned tweezer that traps the
#  ground state anti-traps the Rydberg state.  The real work is all COMPACT: the ground state's
#  alpha_0(omega), the ION CORE's dynamic polarizability, and the correction to that analytic limit.
# ===================================================================================================
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

## ===================================================================================================
##  RUN S19-1.  THE ONE QUESTION:
##     Is the Rydberg/compact box conflict a matter of PRINCIPLE or of COST?
##
##  S19-0 established that the spline count follows hp, not rbox:  hp = rbox/300 -- which is what
##  `Basics.recommendedGrid` uses -- pins nsL near 110 at EVERY box, so a 10580 a.u. box spends the
##  same 110 splines over 165x more space than a 64 a.u. box.  hp = rbox/3000 buys ~500 splines at
##  any box, for ~5x the points.
##
##  So: on a box large enough for a diffuse Rydberg orbital, does raising the spline count restore a
##  COMPACT state's polarizability?  If yes, S19 is a cost problem.  If no, one basis cannot serve
##  both and a magic wavelength needs two.
##
##  ONE VARIABLE PER ROW.  Reference values are EXACT, not literature:
##     alpha(H 1s) = 4.5 a.u. exactly (nonrelativistic);  the Dirac orbital energies from
##     Basics.computeDiracEnergy, which is a closed form.
## ===================================================================================================

nm   = Nuclear.Model(1.0, PointNucleus())
scf  = Basics.NuclearField()
exact1s  = Basics.computeDiracEnergy(Subshell("1s_1/2"),  1.0)
exact16s = Basics.computeDiracEnergy(Subshell("16s_1/2"), 1.0)
exact46s = Basics.computeDiracEnergy(Subshell("46s_1/2"), 1.0)
@printf("\nexact Dirac:  1s = %.10f   16s = %.10f   46s = %.10f  Ha\n", exact1s, exact16s, exact46s)

oneOrbital = function(conf::String, grid::Radial.Grid)
    try
        asf = AsfSettings(AsfSettings(); scField=scf, gridStopper=false)
        mp  = redirect_stdout(devnull) do
            perform(Atomic.Computation(Atomic.Computation(); name=conf, grid=grid, nuclearModel=nm,
                     configs=[Configuration(conf)], asfSettings=asf); output=true)[ResultKeys.Multiplet]
        end
        return( mp.levels[1].energy )
    catch e
        return( NaN )
    end
end

alpha1s = function(grid::Radial.Grid)
    try
        asf = AsfSettings(AsfSettings(); scField=scf, gridStopper=false)
        gm  = redirect_stdout(devnull) do
            perform(Atomic.Computation(Atomic.Computation(); name="np", grid=grid, nuclearModel=nm,
                     configs=[Configuration("$(n)p") for n = 2:16], asfSettings=asf);
                     output=true)[ResultKeys.Multiplet]
        end
        set = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings();
                  multipoles=[E1], gMultiplet=gm, omegas=[0.], printBefore=false)
        oc  = redirect_stdout(devnull) do
            perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=grid, nuclearModel=nm,
                     configs=[Configuration("1s")], asfSettings=asf, propertySettings=[set]);
                     output=true)[MultipolePolarizibility.Settings][1]
        end
        return( oc.alpha0.Babushkin )
    catch e
        return( NaN )
    end
end

@printf("\n%8s %7s %6s | %12s %12s %12s | %12s %9s\n",
        "rbox", "hp/rbox", "nsL", "E(1s) err", "E(16s)/ex", "E(46s)/ex", "alpha(1s)", "/4.5")
@printf("%s\n", "-"^95)
for  rbox in (64.0, 775.0, 3000.0, 10580.0)
    for  (tag, div) in (("1/300", 300.0), ("1/1000", 1000.0), ("1/3000", 3000.0))
        g  = Radial.Grid(Radial.Grid(false), rnt=2.0e-6, h=5.0e-2, hp=rbox/div, rbox=rbox)
        t0 = time()
        e1  = oneOrbital("1s",  g)
        e16 = oneOrbital("16s", g)
        e46 = oneOrbital("46s", g)
        a   = alpha1s(g)
        @printf("%8.0f %7s %6d | %12.3e %12.6f %12.6f | %12.7f %9.5f   [%.0f s]\n",
                rbox, tag, g.nsL, abs(e1-exact1s), e16/exact16s, e46/exact46s, a, a/4.5, time()-t0)
        flush(stdout)
    end
end

## ===================================================================================================
##  RUN S19-2.  THE ONE QUESTION, forced by S19-1:
##     On a LARGE box, can the perturber SET be extended until alpha(H 1s) is recovered?
##
##  S19-1 showed that raising the spline count from 99 to 492 leaves alpha(1s)/4.5 at 0.81144 --
##  IDENTICAL to five digits -- on a 775 a.u. box.  So that failure is NOT a starved basis and more
##  resolution cannot cure it.  The remaining explanation is COMPLETENESS OF THE SET: on a wide box
##  the np n = 2..16 states are all genuine BOUND orbitals, so the sum has no continuum in it at all,
##  whereas on a narrow box the upper members are squeezed into pseudo-continuum states that stand in
##  for it.
##
##  If so, extending n_max on the wide box must restore alpha -- and S19 is a COST problem after all.
##  ONE VARIABLE: n_max.  Box and spline count held fixed.  Reference: the EXACT 4.5 a.u.
## ===================================================================================================

nm   = Nuclear.Model(1.0, PointNucleus())
scf  = Basics.NuclearField()
asf  = AsfSettings(AsfSettings(); scField=scf, gridStopper=false)

alphaOf = function(grid::Radial.Grid, nmax::Int64)
    try
        gm = redirect_stdout(devnull) do
            perform(Atomic.Computation(Atomic.Computation(); name="np", grid=grid, nuclearModel=nm,
                     configs=[Configuration("$(n)p") for n = 2:nmax], asfSettings=asf);
                     output=true)[ResultKeys.Multiplet]
        end
        set = MultipolePolarizibility.Settings(MultipolePolarizibility.Settings();
                  multipoles=[E1], gMultiplet=gm, omegas=[0.], printBefore=false)
        oc  = redirect_stdout(devnull) do
            perform(Atomic.Computation(Atomic.Computation(); name="H 1s", grid=grid, nuclearModel=nm,
                     configs=[Configuration("1s")], asfSettings=asf, propertySettings=[set]);
                     output=true)[MultipolePolarizibility.Settings][1]
        end
        dEmax = maximum(c.nuLevel.energy for c in oc.contributions) -  oc.level.energy
        dEmin = minimum(c.nuLevel.energy for c in oc.contributions) -  oc.level.energy
        return( (oc.alpha0.Babushkin, length(oc.contributions), dEmin, dEmax) )
    catch e
        return( (NaN, 0, NaN, NaN) )
    end
end

for  (rbox, div) in ((775.0, 1000.0), (3000.0, 1000.0))
    g = Radial.Grid(Radial.Grid(false), rnt=2.0e-6, h=5.0e-2, hp=rbox/div, rbox=rbox)
    @printf("\n=== box %.0f a.u.,  nsL = %d  (exact alpha = 4.5) ===\n", rbox, g.nsL)
    @printf("%7s %8s %11s %13s %13s %11s %8s\n",
            "n_max", "poles", "alpha(1s)", "dE_min", "dE_max", "span", "/4.5")
    for  nmax in (16, 25, 40, 60, 90)
        t0 = time()
        a, np, dmin, dmax = alphaOf(g, nmax)
        @printf("%7d %8d %11.6f %13.6e %13.6e %11.2f %8.5f   [%.0f s]\n",
                nmax, np, a, dmin, dmax, dmax/dmin, a/4.5, time()-t0)
        flush(stdout)
    end
end
