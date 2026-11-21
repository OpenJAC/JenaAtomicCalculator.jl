#!/bin/bash
#
# graspTransitions.sh -- drive GRASP2018 for ONE Dirac-Hartree-Fock case and print its E1 lines.
#
#   usage:  tools/graspTransitions.sh <Z> "<conf1>;<conf2>[;...]" "<active set>" "<2J lo,2J hi>"
#   e.g.    tools/graspTransitions.sh 18 "1s(2,*);1s(1,*)2p(1,*)" "1s,2p" "0,4"
#
# The transition half of the cross-code benchmark;  the level
# half is tools/graspLevels.sh and the JAC half is tools/crossCheckBenchmark.jl.
#
# WHAT IS PINNED, and why each choice is forced -- the same discipline as graspLevels.sh:
#   * a POINT nucleus, Coulomb only, no Breit and no QED (rci is NOT run);
#   * the SAME configuration list on both sides, so the CSF space is identical;
#   * ONE orbital set for the initial and the final state.  That last point matters more here
#     than for a level:  with two different sets, GRASP must go through `rbiotransform` and JAC
#     through its own biorthogonal transformation, and the comparison would then also be testing
#     two implementations of that transformation.  With one set, GRASP ITSELF says the
#     transformation is unnecessary -- rbiotransform refuses with "Initial and final states are
#     identical ... just copy name.w to name.bw" -- and that is exactly what this script does.
#
# THE TRANSITIONS MUST CONNECT DIFFERENT CONFIGURATIONS.  The velocity form of an electric
# multipole is structurally zero between two orbitals of the SAME subshell, so an
# intra-configuration line has no meaningful Coulomb-gauge rate to compare.  Give at least two
# configurations of opposite parity.
#
# OUTPUT:  the lines of GRASP's name.name.t table, one row per gauge, as
#     E1  <upper J/P index>  <lower J/P index>  <E in Kays>  <A_Coulomb>  <A_Babushkin>
set -e

Z="$1";  CONFS="$2";  ACTIVE="$3";  JRANGE="$4";  KEEP="${5:-}"
[ -n "$JRANGE" ] || { >&2 echo "usage: $0 <Z> \"<conf1>;<conf2>\" \"<active>\" \"<2J lo,2J hi>\" [keepdir]"; exit 2; }
B="${GRASP_BINDIR:-$HOME/grasp-oracle/bin}"
for exe in rnucleus rcsfgenerate rwfnestimate rangular rmcdhf rtransition; do
    [ -x "$B/$exe" ] || { >&2 echo "ERROR: $B/$exe missing"; exit 127; }
done
# A KEEPDIR IS OFFERED BECAUSE THE LOGS ARE THE ONLY RECORD OF WHAT GRASP WAS ACTUALLY ASKED.
# The prompt answers in scf.in below are positional, so the only way to confirm which levels entered
# the average -- and under which weighting -- is to read rmcdhf.log afterwards.  Without this the
# directory is deleted and that question cannot be answered at all.
if [ -n "$KEEP" ]; then
    mkdir -p "$KEEP";  WORK=$(cd "$KEEP" && pwd)
else
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/graspTrans.XXXXXXXX");  trap 'rm -rf "$WORK"' EXIT
fi
cd "$WORK"

"$B/rnucleus" > rnucleus.log 2>&1 <<EOF
$Z
0
1.0
0
0
0
EOF

# rcsfgenerate takes the reference configurations one per line, then a blank line.
{ echo "*"; echo "0"
  echo "$CONFS" | tr ';' '\n'
  echo ""
  echo "$ACTIVE"; echo "$JRANGE"; echo "0"; echo "n"; } > gen.in
"$B/rcsfgenerate" > rcsfgenerate.log 2>&1 < gen.in
mv -f rcsf.out rcsf.inp

# THE PER-BLOCK LEVEL LINES ARE DERIVED, NOT TYPED.  rmcdhf asks for the ASF serial numbers of
# each block in turn, and the number of blocks and their sizes are whatever rcsfgenerate decided.
# Reading them back from its own report is the only way this stays correct when the configuration
# list changes;  a hard-coded list of "1" lines silently mis-assigns the moment a block appears.
awk '/block  *J\/P/{f=1;next} f&&NF==3{print $3} f&&NF<3{f=0}' rcsfgenerate.log > blocks.txt
NB=$(grep -c . blocks.txt)
[ "$NB" -ge 1 ] || { >&2 echo "ERROR: could not read the block structure"; exit 1; }
{ echo "y"; while read -r n; do echo "1-$n"; done < blocks.txt; echo "1"; echo "5"; echo "*"
  echo "*"; echo "100"; } > scf.in

"$B/rwfnestimate" > rwfnestimate.log 2>&1 <<EOF
y
2
*
EOF
"$B/rangular" > rangular.log 2>&1 <<< "y"
"$B/rmcdhf"   > rmcdhf.log    2>&1 < scf.in

# GRASP's own instruction for one orbital set:  copy rather than transform.
cp rcsf.inp xx.c;  cp rwfn.out xx.w;  cp rmix.out xx.m
cp xx.w xx.bw;     cp xx.m xx.bm
printf 'y\nn\nxx\nxx\nE1\n' | "$B/rtransition" > rtransition.log 2>&1

echo "# graspTransitions.sh   Z = $Z   configs = $CONFS   active = $ACTIVE   2J = $JRANGE"
echo "# oracle commit $(cd "$(dirname "$B")" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "# one orbital set;  point nucleus;  Coulomb only;  no Breit;  no QED"
echo "# upperJP upperIdx lowerJP lowerIdx   E[Kays]        A_Coulomb       A_Babushkin"
awk '
  /^ f  +[0-9]/ {
      # GRASP row:  f <ulev> <uJ> <uP>  f <llev> <lJ> <lP>  <E/Kays>  C <A>  <gf>  <S>
      ui=$2; uj=$3; up=$4; li=$6; lj=$7; lp=$8; ee=$9;
      if ($10=="C") { ac=$11 }
      getline; if ($1=="B") { ab=$2 }
      gsub("D","E",ac); gsub("D","E",ab); gsub("D","E",ee);
      printf "  %s%s %-3s %s%s %-3s %18.4f %16.6e %16.6e\n", uj,up,ui, lj,lp,li, ee+0, ac+0, ab+0
  }' xx.xx.t
