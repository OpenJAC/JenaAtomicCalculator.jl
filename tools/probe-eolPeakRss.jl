#
# probe-eolPeakRss.jl   --   ONE EOL solve, run in its OWN process, so its peak memory can be measured.
#
# THE TWO-POINT MEASUREMENT ITEM 29 OWES, so the EOL cost law can be re-fitted.  The law at
# module-SelfConsistent-inc-optimizedlevel.jl (~line 1149) reads
#
#     memEol = 1.58 + 1.0e-3 * nCsf * (1.55 + 0.097*(nLev-1))          GB
#
# and was fitted on 04-Sep-2026 to peak RSS of real runs.  Item 29 then replaced the pair cache -- the term that
# grew as n^2 within a block -- and measured it 15x smaller, which makes 1.55 an UPPER BOUND rather than a
# prediction.  The law was left in place with a note that re-fitting it needs that measurement redone.
#
# THIS SCRIPT TAKES ONE CASE INDEX AND SOLVES ONE CASE, AND THAT IS THE WHOLE POINT.  Measuring several solves
# inside one process does NOT work, and the way it fails is quiet: RSS never falls, so once the first solve has
# grown Julia's heap, every later solve allocates inside that heap and its peak reads as almost nothing.
# Measured 14-Sep-2026 in one process: 0.74 GB for the first case and then 0.03 / 0.06 / 0.10 GB for cases with
# 2.5x, 5.6x and 10.5x its CSF count -- which is heap reuse, not a cost that fell.  Resetting VmHWM through
# /proc/self/clear_refs does not help, because the reset floor is the RSS the process already holds.
#
# SO THE DRIVER RUNS THIS SCRIPT ONCE PER CASE UNDER `/usr/bin/time -v`, and the number that means something is
# that process's Maximum resident set size.  A baseline run with case index 0 loads JAC and solves nothing, so
# the JAC-load constant can be subtracted rather than assumed.
#
using JenaAtomicCalculator, Printf

const CORE = "1s^2 2s^2 2p^6 3s^2 3p^6 "
# THE CASES COME IN TWO SERIES, AND THE SECOND IS THE ONE THAT DECIDES THE SLOPE.
#
# Series A grows the VALENCE SHELLS, so the CSF count and the SUBSHELL count rise together (10 -> 16).  A fit
# over it cannot say which of the two the cost follows, and that matters: the law is written per CSF, while a
# good deal of the machinery -- B-spline matrices per kappa, the kink-aware tensor cache per rank -- follows the
# subshell count instead.  Series A is kept because it is what a user actually builds.
#
# Series B HOLDS THE SUBSHELL SET FIXED at series A's largest (16 subshells, through 4f) and grows the CSF count
# alone, by opening the 3p core into those same shells.  Every space in series B spans the same subshells, so
# any rise across it is a CSF effect and nothing else.  That is the measurement the re-fit needs.
const BASE4F = [CORE*"3d^2", CORE*"3d 4s", CORE*"4s^2", CORE*"3d 4p", CORE*"4s 4p", CORE*"4p^2",
                CORE*"3d 4d", CORE*"4s 4d", CORE*"4p 4d", CORE*"4d^2", CORE*"3d 4f", CORE*"4s 4f",
                CORE*"4p 4f", CORE*"4d 4f", CORE*"4f^2"]
const OPEN3P = "1s^2 2s^2 2p^6 3s^2 3p^5 "
const CASES = [
  ("baseline: JAC loaded, nothing solved",  22.0, String[]),
  ("A  3d^2 + 3d4s + 4s^2",                 22.0, BASE4F[1:3]),
  ("A  + 4p",                               22.0, BASE4F[1:6]),
  ("A  + 4d",                               22.0, BASE4F[1:10]),
  ("A  + 4f   (nSub 16, series B base)",    22.0, BASE4F),
  ("B  + 3p^5 3d^2 4p",                     22.0, vcat(BASE4F, [OPEN3P*"3d^2 4p"])),
  ("B  + 3p^5 3d^2 4f",                     22.0, vcat(BASE4F, [OPEN3P*"3d^2 4p", OPEN3P*"3d^2 4f"])),
  ("B  + 3p^5 3d 4s 4p, 3d 4s 4f",          22.0, vcat(BASE4F, [OPEN3P*"3d^2 4p", OPEN3P*"3d^2 4f",
                                                                OPEN3P*"3d 4s 4p", OPEN3P*"3d 4s 4f"])) ]

idx  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 0
# THE TARGET-LEVEL COUNT IS AN ARGUMENT because the law carries a per-level term and that term has never been
# re-measured.  Holding nLev at one value and then quoting a law that varies with it is how the 04-Sep fit came
# to over-predict; series A and B above both run at nLev = 3, and a third series varies nLev at fixed space.
nLev = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
(tag, Z, cstrs) = CASES[idx + 1]

if  isempty(cstrs)
    @printf("RESULT  %d  %-34s  nCsf %6d  nSub %4d  nLev %3d  wall %8.1f\n", idx, tag, 0, 0, 0, 0.0)
else
    confs = [Configuration(c) for c in cstrs]
    nm    = Nuclear.Model(Z)
    grid  = Basics.recommendedGrid(confs, nm; printout=false)
    set   = AsfSettings(AsfSettings(); scField = Basics.EOLField(), eeInteraction = CoulombInteraction(),
                                       eeInteractionCI = CoulombInteraction(), gridStopper = false,
                                       scfRoute = Basics.RotationRoute(12),
                                       levelSelectionCI = LevelSelection(true, indices=collect(1:nLev)))
    t  = @elapsed (mp = redirect_stdout(devnull) do
                            SelfConsistent.performSCF(confs, nm, grid, set; printout=false)
                        end)
    b  = mp.levels[1].basis
    @printf("RESULT  %d  %-34s  nCsf %6d  nSub %4d  nLev %3d  wall %8.1f\n",
            idx, tag, length(b.csfs), length(b.subshells), nLev, t)
end
