#
# verify-alphavariation-old-behavior.jl -- the "before" half of challenge S6's closure evidence (08-Oct-2026).
# Companion to verify-alphavariation-overlap-matching.jl (read that one first).
#
# This script must be run against the PRE-FIX source (`git stash` or `git checkout -- src/module-
# AlphaVariation.jl` with the fix patch saved first) to reproduce the wrong number positional-index matching
# gave through Cf17+'s real level crossing at the module's own default step. Running it against the FIXED
# source reproduces the companion script's x=0.125 row instead and proves nothing about the old behavior.
#

using JenaAtomicCalculator, Printf
const B = JenaAtomicCalculator.Basics

Zcf  = 98.0;  core = "[Xe] 4f^14 5d^10"
refs = [Configuration("$core 6s^2 5f^1"), Configuration("$core 6s^2 6p^1")]
nm   = Nuclear.Model(Zcf)
grid = Basics.recommendedGrid(refs, nm; rnt = 2.0e-7)

asfAL = AsfSettings(AsfSettings(); scField=Basics.ALField(), scfRoute=Basics.AverageLevelRoute(60),
                     eeInteractionCI=Basics.CoulombBreit(0.0))
waAL  = Atomic.Computation(Atomic.Computation(), name="Cf17-AL", grid=grid, nuclearModel=nm,
                           configs=refs, asfSettings=asfAL)
wbAL  = perform(waAL; output=true)
multiplet = wbAL[ResultKeys.Multiplet]

x = 0.125
settings = AlphaVariation.Settings(true, x, true, LevelSelection())
outcomes = AlphaVariation.computeOutcomes(multiplet, nm, grid, refs, asfAL, settings; output=true)
for o in outcomes
    @printf(">>> OLD CODE  J=%-5s  omega = %14.6f Ha  q = %14.3f cm^-1\n",
            string(o.level.J), o.omega, Defaults.convertUnits("energy: from atomic to Kayser", o.q))
end
