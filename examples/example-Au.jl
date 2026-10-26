#
println("Au) Apply & test that the SCF+RAS multi-threading (competition item 17, steps 1-6) changes no result, only wall time.")

using JenaAtomicCalculator, Printf, LinearAlgebra

# WRITTEN 26-Oct-2026, after closing competition item 17 (SCF+RAS multi-threading; commits b2efae13, 76c1f673,
# d10e7336) on the small Be 6-config toy case used while building it. This file puts the SAME claim -- bit-
# identical energy at any thread count, only the wall time changes -- on a mid-size, physically realistic
# case: branch (c) of example-Ao.jl, Al-like Si^+ ([Ne] 3s^2 3p + [Ne] 3p^3), a real Ne core with one open
# valence shell and the EOL rotation route.
#
# MEASURED, 26-Oct-2026, machine quiet (ps aux clean, peer tree unchanged throughout every run), BLAS pinned
# to 1 thread throughout so only Threads.jl-level parallelism can vary:
#
#     code                             -t   wall     alloc    gc%    E (Ha)
#     OLD (9eaea2c5, before steps 4-5)  1   17.65 s  8.71 GB  15.6%  -289.2276140694
#     OLD (9eaea2c5, before steps 4-5) 16   16.24 s  8.71 GB   9.7%  -289.2276140694
#     NEW (current HEAD, steps 1-6)     1   18.59 s  9.03 GB  13.0%  -289.2276140694
#     NEW (current HEAD, steps 1-6)     8   13.37 s  9.03 GB  10.2%  -289.2276140694
#     NEW (current HEAD, steps 1-6)    16   13.34 s  9.03 GB  10.7%  -289.2276140694
#
# ENERGY IS BIT-IDENTICAL ON EVERY ROW -- the correctness claim this file exists to make, now confirmed on a
# case bigger than the toy one. The real, user-facing gain is NEW -t16 against OLD -t1: 24.4% less wall time
# (17.65s -> 13.34s), with -t8 already capturing essentially all of it (Amdahl-capped). There is also a
# small, real cost: the new code is 5.3% SLOWER at -t1 (18.59s vs 17.65s) -- the two-pass bookkeeping step 4
# added to avoid racing on vkCache/grad. Both numbers are honest; neither should be hidden behind the other.
#
# TO REPRODUCE: run this branch twice, as  julia --project=. -t 1  and  julia --project=. -t 16  (or any
# other count), and compare the printed energy -- it must match to every printed digit -- and the printed
# wall time, which will differ. A mismatched energy is the regression this file watches for; a changed wall
# time by itself is not, and should not be treated as one without re-checking the machine was quiet.


"""
`checkThreadingInvariance(label, Z, configs, grid)`  ... runs one SCF (EOL field) and reports the thread
    count it ran under, the resulting energy, and the wall time -- the numbers to compare ACROSS separate
    runs at different -t values to check that threading changes only the wall time. Nothing is returned.
"""
function checkThreadingInvariance(label, Z, configs, grid)
    println("\n  ", label)
    set = AsfSettings(AsfSettings(); scField=Basics.EOLField(), scfRoute=Basics.standardRoute(Basics.EOLField(), 20))
    local multiplet
    LinearAlgebra.BLAS.set_num_threads(1)     ## isolate Threads.jl-level parallelism from BLAS's own pool
    t = @elapsed redirect_stdout(devnull) do
            multiplet = SelfConsistent.performSCF(configs, Nuclear.Model(Z), grid, set; printout=false)
        end
    @printf("     -t %2d   E = %.10f Ha   wall = %6.2f s\n", Threads.nthreads(), multiplet.levels[1].energy, t)

    return( nothing )
end


if  true
    # Last visit:  26-Oct-2026
    # Last successful:  26-Oct-2026 -- E = -289.2276140694 Ha bit-identical at -t1/-t8/-t16 and against the
    #                   OLD pre-item code; NEW -t16 24.4% faster than OLD -t1 (17.65s -> 13.34s), NEW -t1
    #                   5.3% slower than OLD -t1 (bookkeeping overhead, see the table above).
    #
    # Al-like Si^+, [Ne] 3s^2 3p + [Ne] 3p^3 -- the same case as example-Ao.jl branch (c), chosen because it
    # is mid-size (a real core, one open valence shell) rather than the tiny toy case the threading work was
    # built against. ~15-20 s per run on this machine.
    grid = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    checkThreadingInvariance("Al-like Si^+   [Ne] 3s^2 3p + [Ne] 3p^3", 14.,
                  [Configuration("1s^2 2s^2 2p^6 3s^2 3p"), Configuration("1s^2 2s^2 2p^6 3p^3")], grid)
end
