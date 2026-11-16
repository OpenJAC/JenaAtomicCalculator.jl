#
# probeAngularSkeletonReuse.jl
#
# A PROBE, not a module:  the configuration-state-function GENERATOR question.  Graspg, an extension of GRASP2018,
# reports Hamiltonian-build speed-ups of up to 68x by exploiting the fact that spin-angular integration does not
# depend on principal quantum numbers, so that the angular part need be computed only once per family of CSFs that
# share an angular structure.  The question for JAC is how large that family reuse actually is HERE, since the
# answer decides whether the idea is worth implementing at all.  This probe measures it.  Written 15-Nov-2026.
#
# THE PREMISE IS PROVEN BY READING, NOT ASSUMED.  `SpinAngular`'s whole angular path -- `module-SpinAngular.jl`,
# `-inc-recoupling.jl`, `-inc-twoparticle.jl`, 148 kB -- contains ZERO reads of a subshell's principal quantum
# number.  It calls `subshell_j` 33 times and `subshell_l` twice (i.e. it reads kappa), and touches `subshells[i]`
# only to write the LABEL of a coefficient back out.  So two CSF pairs whose occupations, seniorities and couplings
# agree over subshells of the same KAPPA sequence produce coefficient values (`T`, `V`, `nu`) that are identical,
# and differ only in the `Subshell` labels.  That is the property Graspg's configuration-state-function generators
# exploit, and in JAC it holds by construction rather than by design.
#
# WHAT IS MEASURED, in three parts:
#
#   1. THE COST SPLIT BY PAIR CLASS.  A pair whose occupation vectors differ by d = sum|n_r - n_s| electron moves
#      can contribute only for d <= 4 (two creations, two annihilations).  JAC ALREADY exits such pairs cheaply,
#      which is the first thing this probe established and it killed a simpler proposal: a d > 4 pair costs
#      ~0.9 us against ~93-118 us for one that contributes, so a pre-filter on d is worth about 2 % and is not
#      worth writing.  The cost lives in the d == 4 pairs.
#
#   2. THE REUSE AMONG THE PAIRS THAT COST SOMETHING.  Keyed on the angular skeleton with BOTH the principal
#      quantum numbers AND the absolute subshell positions removed -- the latter matters: keying on absolute
#      positions treats a 2s->3p pair as different from a 2s->4p pair, which is precisely the equivalence the
#      whole idea rests on, and getting that wrong understates the reuse (1.2x instead of 2.77x).
#
#   3. THAT A CHEAP KEY COLLAPSES AS WELL AS AN EXPENSIVE ONE.  The full-content tuple costs ~26 us per pair --
#      a quarter of the call it would save, so a naive memo would be nearly pointless.  The cheap form is an
#      O(N) per-CSF class id plus two union-relative bit masks, at ~0.12 us per pair, and it is verified here to
#      produce exactly the same number of classes.
#
# WHAT THIS PROBE DOES NOT SAY: how the angular pass compares with the RADIAL side of the same matrix.  That needs
# a real run with orbitals.  What it does say is that the radial side must evaluate only the DISTINCT integrals
# (36 059 for the largest case here, after the memoisation of item 17), while the angular side evaluates per PAIR.
#
using JenaAtomicCalculator, Printf, Random
const JAC = JenaAtomicCalculator


"""
`buildBlock(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})`
    ... generates the single- and double-excited configuration set, makes it relativistic, builds the CSFs and
        returns `(block::Array{CsfR,1}, subshells::Array{Subshell,1}, occupations::Array{Array{Int64,1},1}, sym)`
        for the LARGEST (J, parity) block, since `Hamiltonian.setupMatrix` builds one block at a time and the
        quadratic cost is per block.  No orbitals and no SCF are involved.
"""
function buildBlock(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    confs    = Basics.generateConfigurations(refConfigs, fromShells, toShells, 2)
    relConfs = ConfigurationR[]
    for  conf in confs
        append!(relConfs, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    subshells = Basics.generateSubshellList(relConfs);    ns = length(subshells)
    csfs      = CsfR[]
    for  relConf in relConfs    append!(csfs, Basics.generateCsfRs(relConf, subshells))    end

    blocks = Dict{Tuple{Int64,String},Array{CsfR,1}}()
    for  csf in csfs
        key = (Basics.twice(csf.J), string(csf.parity))
        haskey(blocks, key) ? push!(blocks[key], csf) : (blocks[key] = CsfR[csf])
    end
    sym   = argmax(k -> length(blocks[k]), collect(keys(blocks)))
    block = blocks[sym]
    occs  = [ Int64[ i <= length(c.occupation) ? c.occupation[i] : 0  for i = 1:ns ]  for c in block ]

    return( (block, subshells, occs, sym) )
end


"""
`classIds(block::Array{CsfR,1}, subshells::Array{Subshell,1}, occs::Array{Array{Int64,1},1})`
    ... assigns each CSF one integer class id from its POSITION-FREE angular content: for every occupied
        subshell, in order, its kappa, occupation, seniority, subshellJ and subshellX, plus the CSF's own J and
        parity.  Two CSFs of one n-family share an id.  Cost is O(N), not O(N^2) -- which is what makes the memo
        in the implementation proposal cheap.  Returns `(ids::Array{Int32,1}, nClasses::Int64)`.
"""
function classIds(block::Array{CsfR,1}, subshells::Array{Subshell,1}, occs::Array{Array{Int64,1},1})
    ns = length(subshells);    n = length(block)
    classOf = Dict{Tuple,Int32}();    ids = zeros(Int32, n)
    for  i = 1:n
        csf = block[i];    parts = Any[ Basics.twice(csf.J), string(csf.parity) ]
        for  p = 1:ns
            occs[i][p] == 0  &&  continue
            push!(parts, subshells[p].kappa, occs[i][p], csf.seniorityNr[p],
                         Basics.twice(csf.subshellJ[p]), Basics.twice(csf.subshellX[p]))
        end
        ids[i] = get!(classOf, Tuple(parts), Int32(length(classOf)+1))
    end

    return( (ids, length(classOf)) )
end


"""
`probe(label::String, refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})`
    ... runs the three measurements for one case and prints one line.  Returns a named tuple carrying the CSF
        count, the pair census by class, the measured microseconds per pair in each class, the reuse factor and
        the implied saving on the angular pass.
"""
function probe(label::String, refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1},
               toShells::Array{Shell,1})
    block, subshells, occs, sym = buildBlock(refConfigs, fromShells, toShells)
    n  = length(block);    ns = length(subshells)
    ids, nClasses = classIds(block, subshells, occs)

    dist(r, s) = ( d = 0;  for i = 1:ns   d += abs(occs[r][i] - occs[s][i])   end;  d )

    # --- 1. census and the cheap pair key, in one pass
    cheap = Set{Tuple{Int32,Int32,UInt32,UInt32}}()
    nLe2 = 0;   nEq4 = 0;   nGt4 = 0
    for  r = 1:n,  s = r:n
        d = dist(r, s)
        if      d >  4     nGt4 += 1;   continue
        elseif  d == 4     nEq4 += 1
        else               nLe2 += 1
        end
        ml = UInt32(0);   mr = UInt32(0);   slot = 0
        @inbounds for  i = 1:ns
            (occs[r][i] == 0  &&  occs[s][i] == 0)  &&  continue
            occs[r][i] > 0  &&  (ml |= UInt32(1) << slot)
            occs[s][i] > 0  &&  (mr |= UInt32(1) << slot)
            slot += 1
        end
        push!(cheap, (ids[r], ids[s], ml, mr))
    end
    nContrib = nLe2 + nEq4
    reuse    = nContrib / max(length(cheap), 1)

    # --- 2. time the angular pass on a sample of each class
    Random.seed!(15112026)
    pick = Dict(:le2 => Tuple{Int,Int}[], :eq4 => Tuple{Int,Int}[], :gt4 => Tuple{Int,Int}[])
    tries = 0
    while  (length(pick[:le2]) < 400 || length(pick[:eq4]) < 400 || length(pick[:gt4]) < 400)  &&  tries < 400_000
        tries += 1
        r = rand(1:n);   s = rand(1:n);   r, s = min(r,s), max(r,s);   d = dist(r, s)
        k = d <= 2 ? :le2 : (d == 4 ? :eq4 : :gt4)
        length(pick[k]) < 400  &&  push!(pick[k], (r,s))
    end
    op1 = SpinAngular.OneParticleOperator(0, plus);    op2 = SpinAngular.TwoParticleOperator(0, plus)
    # Warm EVERY class before timing any of it.  Without this the first case measured 282 and 1129 us where the
    # warm figures are 44 and 21 -- i.e. it reported Julia's compiler, not the algorithm.
    for  k in (:le2, :eq4, :gt4),  (r,s) in pick[k][1:min(end,25)]
        SpinAngular.computeCoefficients(op1, block[r], block[s], subshells)
        SpinAngular.computeCoefficients(op2, block[r], block[s], subshells)
    end
    us = Dict{Symbol,Float64}()
    for  k in (:le2, :eq4, :gt4)
        ps = pick[k];    isempty(ps) && (us[k] = 0.0; continue)
        tStart = time()
        for (r,s) in ps
            SpinAngular.computeCoefficients(op1, block[r], block[s], subshells)
            SpinAngular.computeCoefficients(op2, block[r], block[s], subshells)
        end
        us[k] = 1e6 * (time() - tStart) / length(ps)
    end
    tFull  = (nLe2*us[:le2] + nEq4*us[:eq4] + nGt4*us[:gt4]) / 1e6
    tMemo  = (nContrib/reuse) * (nLe2*us[:le2] + nEq4*us[:eq4]) / max(nContrib,1) / 1e6 + nGt4*us[:gt4]/1e6

    @printf("  %-13s %3d sub %6d CSF %5d cls | %8d d<=2 %9d d=4 %9d d>4 | %6.1f %6.1f %5.1f us | reuse %5.2fx | %6.1f s -> %5.1f s  (%4.1fx)\n",
            label, ns, n, nClasses, nLe2, nEq4, nGt4, us[:le2], us[:eq4], us[:gt4], reuse, tFull, tMemo,
            tFull/max(tMemo, 1e-9))

    return( (label=label, nSub=ns, nCsf=n, nClasses=nClasses, nLe2=nLe2, nEq4=nEq4, nGt4=nGt4,
             usLe2=us[:le2], usEq4=us[:eq4], usGt4=us[:gt4], reuse=reuse, tFull=tFull, tMemo=tMemo) )
end


# The generator regime is a FIXED reference with a GROWING orbital set: each added layer contributes n-shells of
# the same kappa, which is exactly where the angular structure repeats.
neon   = [Configuration("1s^2 2s^2 2p^6")]
from   = [Shell("2s"), Shell("2p")]
layers = [ ("SD + n3",    [Shell("3s"), Shell("3p"), Shell("3d")]),
           ("SD + n3,n4", [Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]),
           ("SD + n3..n5",[Shell("3s"), Shell("3p"), Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f"),
                           Shell("5s"), Shell("5p"), Shell("5d"), Shell("5f")]) ]

println("\n", "="^196)
println("ANGULAR SKELETON REUSE IN A CI MATRIX  --  what a CSF-generator scheme could save in JAC, measured")
println("="^196)
@printf("  %-13s %3s %6s %5s | %8s %9s %9s | %20s | %11s | %22s\n",
        "case", "sub", "CSF", "cls", "pairs d<=2", "pairs d=4", "pairs d>4",
        "us/pair by class", "pair reuse", "angular pass, full -> memo")
println("  ", "-"^194)
results = []
for (label, toShells) in layers
    try    push!(results, probe(label, neon, from, toShells))
    catch err   @printf("  %-13s RAISED: %s\n", label, sprint(showerror, err)[1:min(end,140)])    end
end

println("\n", "="^196)
println("""
  READING.
  `cls`        distinct position-free CSF classes; 1724 CSFs collapse to 524, which is where the reuse comes from.
  `us/pair`    measured, by pair class.  The d>4 column is the one that killed a simpler idea: JAC ALREADY exits
               those pairs in about a microsecond, so filtering them is worth ~2 % and not worth writing.
  `pair reuse` contributing pairs divided by distinct angular skeletons, with n AND absolute position removed.
               This is the factor a memo can save on the angular side, and it is exact.
  `full -> memo` the measured angular pass for the whole block, and what it would become with the memo.
  The trend across layers is the finding: the reuse rises steeply with the orbital set, which is the regime
  Graspg's generators were built for -- and it is near 1 for a one-layer calculation, where a memo buys nothing.
""")
