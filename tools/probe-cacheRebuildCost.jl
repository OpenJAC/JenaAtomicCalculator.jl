#
# probe-cacheRebuildCost.jl   --   WHAT IT WOULD COST IN TIME to hold one symmetry block's angular coefficients
#                                  at a time instead of all of them.
#
# THE TRADE BEING PRICED.  An EOL/RAS solve currently builds the per-CSF-pair angular coefficient store for EVERY
# symmetry block and keeps them ALL in memory for the whole SCF run, because each outer iteration revisits them.
# Measured 14-Sep-2026 on the Ti III balanced row (25 085 CSFs over five blocks), the largest block carries 31 %
# of the pair cost, so keeping only one block at a time would take the peak from about 26 GB to 9-10 GB -- an
# ordinary desktop instead of the compute machine.  The price is that each block's store must then be REBUILT
# every iteration instead of being kept.
#
# SO THE QUESTION THIS FILE ANSWERS IS: what fraction of one SCF iteration is the cache build?  If the build is a
# few per cent of an iteration, the trade is nearly free and should simply be taken.  If it is comparable to the
# iteration, the trade roughly doubles the run and becomes a choice rather than an improvement.
#
# WHAT IS TIMED, AND WHY IT IS NOT THE WHOLE SOLVE.  `cacheCsfPairCoefficientsEOL` is timed per block, directly,
# on a bare `Basics.generateBasis` -- the coefficients are orbital-independent, so no SCF is needed to build them
# and the number is clean.  The iteration cost is taken from a real solve, divided by the iterations it ran, so
# the ratio compares like with like.  The rebuild would happen once per block per iteration, which is exactly
# what is reported.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

Z     = 22.0
refs  = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s"),
         Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")]
syms  = [LevelSymmetry(J, Basics.plus)  for J = 0:4]
grid  = Basics.recommendedGrid(refs, Nuclear.Model(Z); printout=false)
core  = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
valp  = [Shell("3p"), Shell("3d"), Shell("4s")]

for (tag, to) in [("3p-v",  [Shell("3d"), Shell("4s")]),
                  ("3p-4p", [Shell("3d"), Shell("4s"), Shell("4p")])]
    step  = RasStep(RasStep(); seFrom=valp, seTo=to, deFrom=valp, deTo=to, frozen=core)
    basis = Basics.generateBasis(refs, syms, step)
    blocks = unique( [ LevelSymmetry(c.J, c.parity)  for c in basis.csfs ] )
    println("\n", "="^92);   println("$tag :  $(length(basis.csfs)) CSFs over $(length(blocks)) blocks");  println("="^92)
    tBuild = 0.0
    for sym in blocks
        n = count(c -> LevelSymmetry(c.J, c.parity) == sym, basis.csfs)
        t = @elapsed SC.cacheCsfPairCoefficientsEOL(sym, basis)
        tBuild += t
        @printf("   block %-6s n %6d   build %8.2f s\n", string(sym), n, t);   flush(stdout)
    end
    @printf("   ALL BLOCKS BUILT ONCE: %.2f s\n", tBuild);   flush(stdout)

    # A real solve, so the per-iteration cost is measured rather than assumed.
    set = AsfSettings(AsfSettings(); scField=Basics.EOLField(), eeInteraction=CoulombInteraction(),
                                     eeInteractionCI=CoulombInteraction(), gridStopper=false,
                                     scfRoute=Basics.RotationRoute(6),
                                     levelSelectionCI=LevelSelection(true, configurations=refs))
    confs = unique( vcat(refs, Basics.generateConfigurations(
                              Basics.generateConfigurations(refs, valp, to), valp, to)) )
    tSolve = @elapsed (mp = redirect_stdout(devnull) do
                                SelfConsistent.performSCF(confs, Nuclear.Model(Z), grid, set; printout=false)
                            end)
    nIter = 6
    @printf("   solve %.1f s over <= %d iterations  ->  ~%.1f s/iteration\n", tSolve, nIter, tSolve/nIter)
    @printf("   REBUILD PER ITERATION would add %.2f s, i.e. %.0f %% of one iteration\n",
            tBuild, 100*tBuild/(tSolve/nIter));   flush(stdout)
end
