#
# verify-alphavariation-overlap-matching.jl -- the harness behind challenge S6's closure (08-Oct-2026).
#
# S6 was: AlphaVariation.computeOutcomes matched a level across the two shifted-alpha calculations BY ITS
# POSITION in the energy-sorted multiplet, which is unsafe on an ion whose levels can re-order under the
# shift -- exactly the near-degeneracy these alpha-variation-search candidates are chosen FOR. Fixed by
# AlphaVariation.matchLevelByOverlap, which matches by the mixing-coefficient OVERLAP between the two
# calculations instead (same CSF basis throughout, since configs and asfSettings are unchanged).
#
# This script reruns the exact case the challenge was measured on (Cf17+'s own reference levels,
# apps-very-schaefer-californium) and prints the clock-transition q-factor -- q(J=1/2) - q(J=5/2) -- at a
# safe step (x=0.01, where the levels do not cross) and at the module's own dangerous DEFAULT step
# (x=0.125, where they do). Keep this script; it is the regression evidence for the fix, not scratch.
#
# MEASURED 08-Oct-2026, for the record (run tools/verify-alphavariation-old-behavior.jl to reproduce the
# "before" row; that script reverts to the pre-fix code for one run and must not be run against the fixed
# source):
#
#     method                  x       clock q [cm^-1]      against published (-450 000 cm^-1)
#     ---------------------------------------------------------------------------------------
#     fixed (overlap match)  0.01     -448 974              0.2 %
#     fixed (overlap match)  0.125    -450 473              0.1 %   (through the real crossing)
#     OLD   (index match)    0.125    +99 470               WRONG SIGN, wrong by a factor of ~4.5
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
println(">>> reference multiplet energies (Ha), lowest to highest:")
for l in sort(multiplet.levels, by=x->x.energy)
    @printf(">>>   J=%-5s parity=%s  E = %18.10f\n", string(l.J), string(l.parity), l.energy)
end

for x in (0.01, 0.125)
    println("\n>>> ===== variationX = $x =====")
    settings = AlphaVariation.Settings(true, x, true, LevelSelection())
    outcomes = AlphaVariation.computeOutcomes(multiplet, nm, grid, refs, asfAL, settings; output=true)
    byJ = Dict(string(o.level.J) => o for o in outcomes)
    qClock = Defaults.convertUnits("energy: from atomic to Kayser", byJ["1/2"].q - byJ["5/2"].q)
    @printf(">>> clock-transition q = q(J=1/2) - q(J=5/2) = %12.1f cm^-1  (published: -450000 cm^-1)\n", qClock)
end
