#
# probe-eolMemoryTerms.jl   --   WHERE THE MEMORY OF AN EOL/RAS LAYER ACTUALLY SITS, measured term by term.
#
# THIS IS THE MEASUREMENT ITEM 29 STILL OWES.  The pair cache was re-formed on 14-Sep-2026 and measured 15x
# smaller (125.24 -> 8.30 MB over the five symmetry blocks of an 887-CSF Ti III space).  That proves the CACHE
# shrank.  It does NOT prove the JOB fits, and the cost law at
# `module-SelfConsistent-inc-optimizedlevel.jl` is deliberately left as an UPPER BOUND until the terms are
# re-measured -- which is what this file does.
#
# THE TERMS OF AN EOL SOLVE THAT LIVE FOR THE WHOLE RUN:
#   (A) blockCaches        the PairCoefficientCache of every relevant symmetry block          -- the item-29 term
#   (B) coeffs1p/coeffs2p  the COMBINED angular coefficients returned by combineAngularCoefficientsEOL
#   (C) the UNCONDENSED intermediate that (B) builds before condensing -- never measured, never freed early
#   (D) tensorCaches, primitives, storage                                                     -- orbital machinery
#
# (C) IS THE ONE TO WATCH.  `combineAngularCoefficientsEOL` pushes EVERY coefficient of EVERY non-zero CSF pair
# into one flat vector and only then condenses it by an O(N^2) double loop.  That is the same explicit per-pair
# storage item 29 removed from the cache, rebuilt at full size in a temporary -- and its condensation is
# quadratic in the entry count, so it is a time term as well as a memory one.
#
# THE ELEMENT TYPE IS ALSO MEASURED, and deliberately: those vectors are built as `Coefficient2p[]`, where
# `Coefficient2p` is a UnionAll (it carries a kind parameter).  A Vector of a non-concrete element type stores
# POINTERS to individually heap-boxed structs rather than the structs inline, which costs both bytes and a GC
# object per entry.  The probe prints the measured bytes/entry beside the inline size so the difference is a
# number rather than a claim.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent
const SA = JenaAtomicCalculator.SpinAngular
const B  = JenaAtomicCalculator.Basics

mb(x) = Base.summarysize(x) / 1024^2

"Build an AL basis + multiplet for `confs`, so the CSFs, orbitals and mixing vectors are all real."
function buildCase(confs, Z)
    nm   = Nuclear.Model(Z)
    grid = Basics.recommendedGrid(confs, nm; printout=false)
    set  = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                      eeInteractionCI = CoulombInteraction(), gridStopper = false)
    mp   = redirect_stdout(devnull) do
               SelfConsistent.performSCF(confs, nm, grid, set; printout=false)
           end
    return( (mp.levels[1].basis, mp, grid, nm) )
end

"Count, WITHOUT building it, the uncondensed entries combineAngularCoefficientsEOL would push."
function uncondensedEntries(blockCaches, targetLevels)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumW    = sum( twiceJp1(level.J)  for level in targetLevels )
    weights = [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    n1 = 0;   n2 = 0;   nPairs = 0
    for  (_, cache)  in  blockCaches
        idx = cache.idxCsf;   n = length(idx)
        for r = 1:n, s = 1:n
            drs = 0.
            for (i, level) in enumerate(targetLevels)   drs += weights[i]*level.mc[idx[r]]*level.mc[idx[s]]   end
            drs == 0.  &&  continue
            nPairs += 1
            k  = (s-1)*n + r
            n1 += cache.ptr1p[k+1] - cache.ptr1p[k]
            n2 += cache.ptr2p[k+1] - cache.ptr2p[k]
        end
    end
    return( (n1, n2, nPairs) )
end

println("="^118)
println("WHERE AN EOL LAYER'S MEMORY SITS -- term by term.  Item 29's owed measurement.")
println("="^118)

cases = [ ("Mg-like  n=3 complex",  12.0,
             ["1s^2 2s^2 2p^6 3s^2", "1s^2 2s^2 2p^6 3s 3p", "1s^2 2s^2 2p^6 3s 3d",
              "1s^2 2s^2 2p^6 3p^2", "1s^2 2s^2 2p^6 3p 3d", "1s^2 2s^2 2p^6 3d^2"]),
          ("Si-like  n=3 complex",  14.0,
             ["1s^2 2s^2 2p^6 3s^2 3p^2", "1s^2 2s^2 2p^6 3s 3p^3", "1s^2 2s^2 2p^6 3p^4",
              "1s^2 2s^2 2p^6 3s^2 3p 3d", "1s^2 2s^2 2p^6 3s^2 3d^2", "1s^2 2s^2 2p^6 3s 3p^2 3d",
              "1s^2 2s^2 2p^6 3p^3 3d"]),
          ("Ti III  3d^2 complex",  22.0,
             ["1s^2 2s^2 2p^6 3s^2 3p^6 3d^2", "1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s",
              "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2", "1s^2 2s^2 2p^6 3s^2 3p^6 3d 4d",
              "1s^2 2s^2 2p^6 3s^2 3p^6 4p^2", "1s^2 2s^2 2p^6 3s^2 3p^6 3d 4p",
              "1s^2 2s^2 2p^6 3s^2 3p^6 4s 4d", "1s^2 2s^2 2p^6 3s^2 3p^6 4d^2"]) ]

@printf("\n%-24s %7s %5s %9s %10s %10s %10s %9s %9s %8s\n",
        "case", "nCsf", "nSub", "pairs!=0", "cache MB", "allocMB", "uncondMB", "cond MB", "B/entry", "comb s")
println("-"^118)

for (tag, Z, cstrs) in cases
    confs = [Configuration(c) for c in cstrs]
    (basis, mp, grid, nm) = buildCase(confs, Z)
    syms  = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    caches = Dict( sym => SC.cacheCsfPairCoefficientsEOL(sym, basis)  for sym in syms )
    # the target levels an EOL run would carry: the lowest few, with their real mixing vectors
    tgt   = mp.levels[1:min(3, length(mp.levels))]
    (n1, n2, nPairs) = uncondensedEntries(caches, tgt)
    allocB = @allocated ((c1, c2) = SC.combineAngularCoefficientsEOL(caches, tgt))
    tComb  = @elapsed  SC.combineAngularCoefficientsEOL(caches, tgt)

    cacheMB  = sum( mb(c) for (_, c) in caches )
    condMB   = mb(c1) + mb(c2)
    # measured bytes per stored entry of the CONDENSED vector -- compare with the inline struct size
    bPer     = (length(c1) + length(c2)) == 0 ? 0.0 :
               (Base.summarysize(c1) + Base.summarysize(c2)) / (length(c1) + length(c2))
    uncondMB = (n1 * (Base.summarysize(c1) / max(length(c1),1)) +
                n2 * (Base.summarysize(c2) / max(length(c2),1))) / 1024^2
    cmp      = (n1^2 + n2^2) / 2

    @printf("%-24s %7d %5d %9d %10.3f %10.3f %10.3f %10.3f %9.1f %8.2f\n",
            tag, length(basis.csfs), length(basis.subshells), nPairs, cacheMB,
            allocB/1024^2, uncondMB, condMB, bPer, tComb)
    @printf("%-24s   uncondensed entries: 1p %d, 2p %d  ->  condensed: 1p %d, 2p %d  (%.1fx);  O(N^2) cmp %.2e\n",
            "", n1, n2, length(c1), length(c2),
            (n1+n2) / max(length(c1)+length(c2), 1), cmp)
    flush(stdout)
end

println("\nINLINE struct sizes for reference:  Coefficient1p ", sizeof(SA.Coefficient1p{SA.OrdinaryKind}),
        " B,  Coefficient2p ", sizeof(SA.Coefficient2p{SA.EffectiveStrengthKind}), " B.")
println("The B/entry column is what the CONDENSED vector actually costs; a value well above the inline size")
println("means the vector's element type is not concrete and every entry is a separately boxed heap object.")
