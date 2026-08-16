#
println("Ao) Apply & test the average-level (AL) and optimized-level (EOL) self-consistent fields, and where they differ.")

using Printf

# WRITTEN 10-Aug-2026; REWRITTEN 16-Aug-2026, when EOLField was wired to the orbital-rotation solver.
# The AL field optimizes a configuration AVERAGE; the EOL field optimizes the (2J+1)-weighted energy of the
# SELECTED target level(s).  For those levels EOL must therefore come out at or below AL, and these branches
# check exactly that.
#
# WHY THE SOLVER WAS CHANGED.  Until 16-Aug-2026 EOLField reached SelfConsistent.solveOptimizedLevelField,
# which folds the off-diagonal CSF-pair terms into the same (1/occ)-scaled Fock matrix.  That converges to a
# DEGENERATE stationary point whenever two near-degenerate CSFs compete for one correlation channel: the
# correlating weight runs away to zero, the correlation orbital then no longer enters the energy at all, and
# its gradient vanishes for a trivial reason -- it looks converged because it IS stationary, at the wrong
# critical point.  Measured with that solver, and the reason it was replaced:
#
#     case                          AL            EOL          EOL - AL     correlating weight
#     Be-like C^2+  2s^2 + 2p^2   -36.49373010  -36.45536243   +0.038       0.2258 -> 0.0002
#     Li-like C^3+  2s + 3s + 3d  -34.74190558  -34.74194657   -0.000041    no competition
#     Si^+  3s^2 3p + 3p^3       -289.22410301 -289.22455590   -0.00045     no competition
#     Si^2+ 3s^2 + 3p^2          -288.66107449 -288.64011958   +0.021       0.1825 -> -0.0000
#
# It failed precisely where TWO NEAR-DEGENERATE CSFs COMPETE FOR ONE CORRELATION CHANNEL (2s^2 vs 2p^2,
# 3s^2 vs 3p^2) and succeeded where a single CSF dominates -- not a matter of system size, nuclear charge or
# the presence of a core: the light two-electron-valence case and the Ne-core one behaved identically, and
# the two charge states of the SAME element disagreed with each other.  That solver remains callable
# directly, so the comparison can still be made.
#
# Minimizing over orbital ROTATIONS escapes the degenerate point, and that is what EOLField now does
# (SelfConsistent.solveOptimizedLevelFieldByRotation, started from an average-level basis, since a rotation
# is a LOCAL optimizer).  On Be 1s^2 2s^2 + 1s^2 2p^2 it reaches 5.3 mHa BELOW AL where the old solver sat
# 19.4 mHa above it; on Ne^2+ 2p^4 + 2p^2 3s^2 + 2p^2 3p^2, optimizing the lowest four levels, it gains
# 178 mHa over AL and -- the sharper test -- puts the 2p^4 3P multiplet in its correct INVERTED order
# (3P_2 lowest, as a more-than-half-filled shell requires), where AL has it the wrong way round.
#
# ORTHONORMALITY.  Each branch also reports the worst same-kappa overlap, which the CSF expansion requires
# to vanish.  Until 10-Aug-2026 it did not: the SCF damping step, mixed = 0.5*old + 0.5*raw, destroyed the
# orthogonality that Hamiltonian.projectHamiltonian had just enforced, and nothing restored it -- converged
# Li reached <2s|3s> = 1.1e-03 while all 44 approved tests passed.  Both drivers now re-orthonormalize;
# the values below should all be <= 1e-08, and TestFrames.testMethod_OrbitalOrthonormality asserts it.
#
# COST.  Branches a and b run in well under a minute and are the ones to re-run routinely; c and d carry a
# real Ne core and take two to four minutes per field.  Ions rather than neutrals throughout: a higher
# charge contracts the valence orbitals, so a much smaller box suffices and the cost falls with it.


"""
`compareFields(label, Z, configs, grid)`  ... runs the same configuration list through the AL and the EOL
    field, and reports for each: the target-level energy, the leading mixing coefficients, and the worst
    same-kappa deviation from orthonormality. Nothing is returned.
"""
function compareFields(label, Z, configs, grid)
    println("\n  ", label)
    energies = Dict{String,Float64}()
    for  (tag, field)  in  [("AL ", Basics.ALField()), ("EOL", Basics.EOLField())]
        ## redirect_stdout: the driver prints a per-orbital "overlap = ..." line for every subshell of
        ## every iteration, and those printlns are NOT gated by the printout keyword it accepts.
        local multiplet
        t = @elapsed redirect_stdout(devnull) do
                multiplet = SelfConsistent.performSCF(configs, Nuclear.Model(Z), grid,
                                AsfSettings(AsfSettings(); scField=field, maxIterationsScf=20); printout=false)
            end
        basis = multiplet.levels[1].basis;    orbitals = basis.orbitals;    worst = 0.
        for  (i, sha)  in  enumerate(basis.subshells),  (j, shb)  in  enumerate(basis.subshells)
            if  sha.kappa != shb.kappa   ||   j < i    continue    end
            oa  = orbitals[sha];   ob = orbitals[shb]
            mtp = min( size(oa.P,1), size(ob.P,1), length(grid.wr) )
            ov  = sum( grid.wr[k] * (oa.P[k]*ob.P[k] + oa.Q[k]*ob.Q[k])  for k = 1:mtp )
            worst = max( worst, abs( ov - (i == j ? 1.0 : 0.0) ) )
        end
        energies[strip(tag)] = multiplet.levels[1].energy
        @printf("     %s  E = %15.8f   mix %s   worst same-kappa overlap %.1e   (%.0f s)\n", tag,
                multiplet.levels[1].energy,
                join([@sprintf("%+.4f", x) for x in multiplet.levels[1].mc[1:min(2,end)]], " "), worst, t)
    end
    d = energies["EOL"] - energies["AL"]
    @printf("     EOL - AL = %+.8f Ha   %s\n", d,
            d <= 0. ? "EOL binds its own target level better, as it must" :
                      "<== EOL is WORSE on the level it optimizes: the degenerate stationary point")
    return( nothing )
end


if  true
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ... absolute energies not yet checked against literature
    #
    # a) Be-like C^2+, 1s^2 2s^2 + 1s^2 2p^2.  The canonical failure: two CSFs of the same J = 0+ block
    #    competing for one correlation channel.  EOL drives the 2p^2 weight from AL's 0.2258 down to 0.0002
    #    and lands 0.038 Ha ABOVE AL on the very level it is optimizing.  Cheapest branch here; use it when
    #    re-checking anything about EOL.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    compareFields("Be-like C^2+   1s^2 2s^2 + 1s^2 2p^2", 6.,
                  [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")], grid)
    #
elseif  false
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ... absolute energies not yet checked against literature
    #
    # b) Li-like C^3+, 1s^2 2s + 1s^2 3s + 1s^2 3d.  A single valence electron, so no two CSFs compete and
    #    EOL behaves correctly -- 4.1e-05 Ha BELOW AL.  Its value here is the kappa = -1 block, which holds
    #    THREE orbitals (1s, 2s, 3s): that is what exposes the orthonormality defect, since two orbitals
    #    only reach ~5e-05 while three reach 1.1e-03.  The reported overlap should now be ~1e-09.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 10.0)
    compareFields("Li-like C^3+   1s^2 2s + 1s^2 3s + 1s^2 3d", 6.,
                  [Configuration("1s^2 2s"), Configuration("1s^2 3s"), Configuration("1s^2 3d")], grid)
    #
elseif  false
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ... absolute energies not yet checked against literature
    #
    # c) Al-like Si^+, [Ne] 3s^2 3p + [Ne] 3p^3.  An OPEN valence shell with J != 0, a real Ne core, and one
    #    dominant CSF -- so again no competition, and EOL comes out 4.5e-04 Ha below AL.  Together with (d),
    #    which uses the same element at a different charge, this shows the defect is about the CSF structure
    #    and not about the atom.  ~2-3 minutes per field.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    compareFields("Al-like Si^+   [Ne] 3s^2 3p + [Ne] 3p^3", 14.,
                  [Configuration("1s^2 2s^2 2p^6 3s^2 3p"), Configuration("1s^2 2s^2 2p^6 3p^3")], grid)
    #
elseif  false
    # Last visit:  16-Aug-2026
    # Last successful:  unknown ... absolute energies not yet checked against literature
    #
    # d) Mg-like Si^2+, [Ne] 3s^2 + [Ne] 3p^2.  Branch (a) again, now with a real Ne core beneath it: the
    #    same closed-valence competition, the same collapse (0.1825 -> -0.0000), the same sign and a similar
    #    size of error (+0.021 Ha).  Compare with (c), the SAME element one charge state lower, where EOL is
    #    correct.  ~2-4 minutes per field.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    compareFields("Mg-like Si^2+  [Ne] 3s^2 + [Ne] 3p^2", 14.,
                  [Configuration("1s^2 2s^2 2p^6 3s^2"), Configuration("1s^2 2s^2 2p^6 3p^2")], grid)
    #
end
