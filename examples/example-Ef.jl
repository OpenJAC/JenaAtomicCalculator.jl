
println("Ef) RETIRED -- the plasma-shifted Auger path moved to Plasma.LineShiftScheme; see example-Jb.jl.")

# RETIRED 18-Aug-2026.  This file computed Auger rates in a screened (Debye-Hueckel) plasma by handing an
# AutoIonization.PlasmaSettings to an Atomic.Computation.  That route no longer exists, and the file had never
# earned a date -- its header still read "Last successful:  unknown".  It was broken in THREE independent ways,
# each hidden behind the one before, which is why the first of them was all anyone ever saw:
#
#   1. Basics.DebyeHueckel() was renamed Basics.DebyeHueckelModel(), so the very first line raised.
#   2. AutoIonization.PlasmaSettings was called with SIX arguments against a struct that has TWO,
#      `printBefore` and `lineSelection`.  The plasma model, the Debye length, the ion-sphere radius and the
#      bound-electron count are no longer settings at all: the model is passed SEPARATELY, as its own argument.
#   3. Even repaired, perform(::Atomic.Computation) has no branch for AutoIonization.PlasmaSettings.  The
#      dispatch chain over processSettings ends in error("stop b"), so the computation could not have run.
#
# NOTHING IS LOST, and that is why this file is retired rather than repaired.  The physics has a better home:
# a plasma is a property of the COMPUTATION, not of one process's settings, so it belongs to Plasma.Computation
# with a Plasma.LineShiftScheme, which carries the model once and applies it to the CI, the continuum orbital
# and the amplitude alike.  That path is live and dispatches on AutoIonization.PlasmaSettings itself, in
# Plasma.perform(::LineShiftScheme) -> AutoIonization.computeLinesPlasma(..., plasmaModel).
#
# WHERE TO LOOK:  examples/example-Jb.jl, branch 1, Last successful 18-Jul-2026.  It runs the SAME ion this
# file used, Ne 1s 2s^2 2p^6 -> 1s^2 2p^6, over Basics.NoPlasmaModel() and DebyeHueckelModel at 1000, 10 and
# 2 a_o, and it is verified rather than merely run: the Auger electron energy shift taken from the rate table
# equals the difference of the independently printed initial- and final-level shifts (0.054, 5.20, 19.94 eV) to
# the last printed digit, and the rate rises monotonically 3.368e13 -> 3.448e13 1/s with screening.
