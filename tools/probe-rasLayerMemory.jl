#
# probe-rasLayerMemory.jl   --   WHICH SPACE THE EOL COST MESSAGE IS TALKING ABOUT, and what the RAS layer
#                                around it actually costs.
#
# THE DISCREPANCY THIS FILE EXISTS TO SETTLE.  On 14-Sep-2026 the EOL cost law was re-fitted from peak RSS over
# eleven solves driven through `SelfConsistent.performSCF` (tools/probe-eolPeakRss.jl): 14 to 1231 CSFs, all of
# them between 1.42 and 1.67 GB.  At the same time a Ti III RAS ladder on the compute machine held 26.9 GB while
# its `[EOL-C3] cost estimate` lines named spaces of 14 and 70 CSFs and predicted 1.6 and 1.7 GB.
#
# THE SUSPECTED CAUSE IS A SCOPE MISMATCH, NOT A WRONG LAW.  A RAS ladder is driven by
# `generate(Representation(..., RasExpansion(...)))`, which is a DIFFERENT code path from a direct
# `performSCF` on a configuration list.  The EOL solver is invoked inside it, and the message it prints
# describes THE SPACE THAT SOLVER WAS HANDED -- which for a reference step is small -- while the layer built
# around it may be orders of magnitude larger.  If that is right, the law is not out by 20x; it is answering a
# question about a 14-CSF sub-problem and being read as a prediction for a 25 085-CSF job, which is a defect in
# WHAT IS REPORTED rather than in the arithmetic.
#
# THE TEST.  Run one RAS ladder whose reference is small and whose layer is genuinely large, in its own process,
# and collect BOTH: every `[EOL-C3] cost estimate` line with the CSF count it names, and the process's peak RSS.
# The diagnosis is confirmed if the messages keep naming the reference space while the peak tracks the layer.
#
using JenaAtomicCalculator, Printf

Z     = 22.0
c3d2  = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2")
c3d4s = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s")
c4s2  = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")
refs  = [c3d2, c3d4s, c4s2]
syms  = [LevelSymmetry(J, Basics.plus)  for J = 0:4]
grid  = Basics.recommendedGrid(refs, Nuclear.Model(Z); printout=false)
core  = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
val   = [Shell("3d"), Shell("4s")]

# THE LAYER IS THE ARGUMENT, so the same script gives a small and a large layer without changing anything else.
# OPENING 3p IS WHAT MAKES THE SPACE EXPLODE, and it is what the real ladder does: 3p holds six electrons, so
# allowing excitations out of it multiplies the CSF count by orders of magnitude while the SUBSHELL count barely
# moves.  A layer built from the two valence electrons alone reaches only tens of CSFs -- measured 19 and 49 here
# on 14-Sep-2026 -- and is no test of the regime this file is about.  The 3p-opened rows are.
layer = length(ARGS) >= 1 ? ARGS[1] : "4p"
valp  = [Shell("3p"), Shell("3d"), Shell("4s")]
(from, to) =
  layer == "none"  ? (Shell[],  Shell[]) :
  layer == "4p"    ? (val,  [Shell("3d"), Shell("4s"), Shell("4p")]) :
  layer == "4d"    ? (val,  [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d")]) :
  layer == "3p-v"  ? (valp, [Shell("3d"), Shell("4s")]) :
  layer == "3p-4p" ? (valp, [Shell("3d"), Shell("4s"), Shell("4p")]) :
                     (valp, [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d")])

ref0  = RasStep(RasStep(); frozen=core)
steps = isempty(to) ? [ref0] :
        [ref0, RasStep(RasStep(); seFrom=from, seTo=to, deFrom=from, deTo=to, frozen=core)]
rs    = RasSettings(Int64[], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, configurations=refs))

t  = @elapsed (wb = generate( Representation("Ti III probe", Nuclear.Model(Z), grid, refs,
                                             RasExpansion(syms, 20, steps, rs)), output=true ))
k  = "step" * string(length(steps))
n  = haskey(wb, k) ? length(wb[k].levels[1].basis.csfs) : -1
@printf("\nRESULT  layer %-5s  steps %d  finalCsf %7d  wall %8.1f\n", layer, length(steps), n, t)
