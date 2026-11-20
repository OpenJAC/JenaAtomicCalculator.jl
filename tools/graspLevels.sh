#!/bin/bash
#
# graspLevels.sh -- drive GRASP2018 for ONE Dirac-Hartree-Fock case and print its level energies.
#
#   usage:  tools/graspLevels.sh <Z> "<configuration>" "<active set>" "<2J lower,2J higher>" [workdir]
#   e.g.    tools/graspLevels.sh 92 "1s(1,*)"            "1s"      "1,1"
#           tools/graspLevels.sh 54 "1s(2,*)"            "1s"      "0,0"
#           tools/graspLevels.sh 26 "1s(2,*)2s(2,*)2p(6,*)" "1s,2s,2p" "0,0"
#
# It is the GRASP half of the cross-code benchmark (priority item 23).  The JAC half is
# tools/crossCheckBenchmark.jl, which also holds the design and the calibrations.
#
# WHAT IS FIXED HERE, AND WHY EACH CHOICE IS FORCED.  A cross-code comparison is only sharp if both
# codes solve the SAME problem, so everything that could differ is pinned:
#
#   * A POINT NUCLEUS (mass number 0 to rnucleus).  Two codes' finite-nucleus models are not
#     identical, and a nuclear-model difference would be read as a code difference.
#   * NO Breit, NO QED -- rmcdhf alone, rci is NOT run.  Those are the terms a closed-form Dirac
#     energy does not contain, and the terms whose implementations legitimately differ.
#   * A SINGLE reference configuration, so the CSF list is trivially identical in both codes and
#     the comparison is Dirac-Hartree-Fock against Dirac-Hartree-Fock with no correlation model in
#     between.  This is what makes a disagreement a BUG rather than a modelling difference, and it
#     is the one thing no partner other than GRASP allows.
#   * GRASP'S OWN DEFAULT GRID, left untouched.  For a point nucleus that is
#     RNT = exp(-65/16)/Z, H = 0.0625, HP = 0 (PURELY EXPONENTIAL) and N = 220
#     -- see src/appl/rmcdhf90/getscd.f90.  The JAC side is then given these same four numbers,
#     rather than GRASP being driven off its default path through the interactive grid prompts.
#     NOTE that this is NOT JAC's own default grid (rnt = 2e-6, h = 0.05, hp != 0), so a run of
#     either code on its own recommended grid is NOT a matched comparison.
#
# The binaries are expected at $GRASP_BINDIR, or at ~/grasp-oracle/bin.  The commit of the oracle
# is printed, because a number whose code state is unknown is not evidence.
set -e

Z="$1";  CONF="$2";  ACTIVE="$3";  JRANGE="$4";  WORK="${5:-}"
if [ -z "$JRANGE" ]; then
    >&2 echo "usage: $0 <Z> \"<configuration>\" \"<active set>\" \"<2J lower,2J higher>\" [workdir]"
    exit 2
fi
B="${GRASP_BINDIR:-$HOME/grasp-oracle/bin}"
for exe in rnucleus rcsfgenerate rwfnestimate rangular rmcdhf; do
    [ -x "$B/$exe" ] || { >&2 echo "ERROR: $B/$exe not found or not executable"; exit 127; }
done

if [ -z "$WORK" ]; then
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/graspLevels.XXXXXXXX")
    trap 'rm -rf "$WORK"' EXIT
fi
mkdir -p "$WORK";   cd "$WORK"

# rnucleus:  Z, mass number 0 = POINT nucleus, then neutral mass / spin / dipole / quadrupole,
# none of which enters a point-nucleus Dirac-Hartree-Fock energy.
"$B/rnucleus" > rnucleus.log 2>&1 <<EOF
$Z
0
1.0
0
0
0
EOF

# rcsfgenerate:  default ordering, no predefined core, one reference configuration, the active set,
# the 2J window, and ZERO excitations -- the single-CSF condition above.
"$B/rcsfgenerate" > rcsfgenerate.log 2>&1 <<EOF
*
0
$CONF

$ACTIVE
$JRANGE
0
n
EOF
mv -f rcsf.out rcsf.inp

"$B/rwfnestimate" > rwfnestimate.log 2>&1 <<EOF
y
2
*
EOF

"$B/rangular" > rangular.log 2>&1 <<EOF
y
EOF

# rmcdhf with DEFAULT settings (the first `y`), so the grid is GRASP's own, as documented above.
"$B/rmcdhf" > rmcdhf.log 2>&1 <<EOF
y
1
1
5
*
*
100
EOF

NCSF=$(grep -c "" rcsf.inp 2>/dev/null || echo "?")
echo "# graspLevels.sh   Z = $Z   conf = $CONF   active = $ACTIVE   2J = $JRANGE"
echo "# oracle commit $(cd "$(dirname "$B")" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "# grid: GRASP default for a POINT nucleus -- RNT = exp(-65/16)/Z, H = 0.0625, HP = 0, N = 220"
echo "# level energies (Hartree), from rmcdhf.sum:"
awk '/^ *[0-9]+ +[0-9]+\/?[0-9]* *[+-]/ { print "  " $0 }' rmcdhf.sum | head -20
