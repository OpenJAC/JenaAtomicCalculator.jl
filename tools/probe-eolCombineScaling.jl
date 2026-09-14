#
# probe-eolCombineScaling.jl   --   how the EOL angular-coefficient TEMPORARY scales with the CSF count.
#
# THE COMPANION TO probe-eolMemoryTerms.jl, which established the RANKING of the terms on small spaces.  This
# one establishes the SCALING, which is what the cost law in module-SelfConsistent-inc-optimizedlevel.jl needs
# and has never had: the law of 04-Sep-2026 was fitted to the OLD pair cache, and item 29 replaced that cache
# without re-fitting it.
#
# WHY THE ORBITALS DO NOT MATTER HERE, and why this probe is therefore cheap.  The per-CSF-pair angular
# coefficients are ORBITAL-INDEPENDENT BY CONSTRUCTION -- that is the whole reason they are cached once and
# reused across SCF iterations.  So a converged AL basis and a one-step DFS basis give byte-identical caches,
# and the SCF can be cut to its cheapest setting without touching what is measured.  Only the mixing vectors
# (level.mc) depend on the orbitals, and they enter only through which pairs have d_rs != 0 -- reported as
# `pairs!=0` so that its effect is visible rather than assumed.
#
# WHAT IS BEING FITTED.  Per SCF ITERATION (this call sits inside the iteration loop, at line ~2275, and runs
# TWICE per iteration when the off-diagonal split is on, which is the default for the rotation route):
#
#     N        the number of UNCONDENSED coefficient entries pushed into one flat vector
#     N x B    the bytes that temporary costs, B measured at 105-178 B/entry because the vector's element type
#              is not concrete and every entry is separately boxed
#     N^2 / 2  the comparisons its condensation tail performs, since that tail is a nested double loop
#
# and the CONDENSED result, which is bounded by the number of distinct (nu,a,b,c,d) labels and therefore by the
# SUBSHELL count rather than the CSF count -- i.e. it saturates.  If that holds over a decade of CSF counts,
# then the whole cost of this routine is a temporary whose answer is small, which is a fixable shape.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

mb(x) = Base.summarysize(x) / 1024^2

"A basis for `confs` at the cheapest SCF setting; see the note above on why this does not affect the measurement."
function cheapBasis(confs, Z)
    nm   = Nuclear.Model(Z)
    grid = Basics.recommendedGrid(confs, nm; printout=false)
    set  = AsfSettings(AsfSettings(); scField = Basics.DFSField(), eeInteraction = CoulombInteraction(),
                                      eeInteractionCI = CoulombInteraction(), gridStopper = false)
    mp   = redirect_stdout(devnull) do
               SelfConsistent.performSCF(confs, nm, grid, set; printout=false)
           end
    return( (mp.levels[1].basis, mp) )
end

function entryCounts(blockCaches, targetLevels)
    twiceJp1(J) = ( J.den == 1 ? 2*J.num : J.num ) + 1
    sumW    = sum( twiceJp1(level.J)  for level in targetLevels )
    weights = [ twiceJp1(level.J) / sumW  for level in targetLevels ]
    n1 = 0;  n2 = 0;  nPairs = 0
    for  (_, cache)  in  blockCaches
        idx = cache.idxCsf;   n = length(idx)
        for r = 1:n, s = 1:n
            drs = 0.
            for (i, level) in enumerate(targetLevels)   drs += weights[i]*level.mc[idx[r]]*level.mc[idx[s]]   end
            drs == 0.  &&  continue
            nPairs += 1;   k = (s-1)*n + r
            n1 += cache.ptr1p[k+1] - cache.ptr1p[k]
            n2 += cache.ptr2p[k+1] - cache.ptr2p[k]
        end
    end
    return( (n1, n2, nPairs) )
end

# Ti III, the ion the 25 085-CSF layer belongs to, in spaces of growing size.  The core stays closed, as the
# maintainer's own recipe asks; the growth is in the valence shells the layers open.
const CORE = "1s^2 2s^2 2p^6 3s^2 3p^6 "
spaces = [
  ("3d^2 + 3d4s + 4s^2",                    [CORE*"3d^2", CORE*"3d 4s", CORE*"4s^2"]),
  ("+ 4p",                                  [CORE*"3d^2", CORE*"3d 4s", CORE*"4s^2", CORE*"3d 4p",
                                             CORE*"4s 4p", CORE*"4p^2"]),
  ("+ 4d",                                  [CORE*"3d^2", CORE*"3d 4s", CORE*"4s^2", CORE*"3d 4p",
                                             CORE*"4s 4p", CORE*"4p^2", CORE*"3d 4d", CORE*"4s 4d",
                                             CORE*"4p 4d", CORE*"4d^2"]),
  ("+ 4f",                                  [CORE*"3d^2", CORE*"3d 4s", CORE*"4s^2", CORE*"3d 4p",
                                             CORE*"4s 4p", CORE*"4p^2", CORE*"3d 4d", CORE*"4s 4d",
                                             CORE*"4p 4d", CORE*"4d^2", CORE*"3d 4f", CORE*"4s 4f",
                                             CORE*"4p 4f", CORE*"4d 4f", CORE*"4f^2"]) ]

println("="^124)
println("HOW THE EOL ANGULAR-COEFFICIENT TEMPORARY SCALES -- Ti III, closed core, growing valence space.")
println("="^124)
@printf("\n%-22s %7s %5s %10s %10s %11s %10s %10s %9s %10s\n",
        "space", "nCsf", "nSub", "pairs!=0", "cache MB", "entries 2p", "temp MB", "cond 2p", "B/entry", "N^2/2")
println("-"^124)

for (tag, cstrs) in spaces
    confs = [Configuration(c) for c in cstrs]
    local basis, mp
    try
        (basis, mp) = cheapBasis(confs, 22.0)
    catch e
        @printf("%-22s  SKIPPED: %s\n", tag, sprint(showerror, e)[1:min(70,end)]);   flush(stdout);   continue
    end
    syms   = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
    caches = Dict( sym => SC.cacheCsfPairCoefficientsEOL(sym, basis)  for sym in syms )
    tgt    = mp.levels[1:min(3, length(mp.levels))]
    (n1, n2, nPairs) = entryCounts(caches, tgt)
    allocB = @allocated ((c1, c2) = SC.combineAngularCoefficientsEOL(caches, tgt))
    cacheMB = sum( mb(c) for (_, c) in caches )
    bPer    = (length(c1)+length(c2)) == 0 ? 0.0 :
              (Base.summarysize(c1)+Base.summarysize(c2)) / (length(c1)+length(c2))
    @printf("%-22s %7d %5d %10d %10.3f %11d %10.3f %10d %9.1f %10.2e\n",
            tag, length(basis.csfs), length(basis.subshells), nPairs, cacheMB, n2,
            allocB/1024^2, length(c2), bPer, (n1+n2)^2/2)
    flush(stdout)
end

println("\nREAD THE LAST TWO COLUMNS TOGETHER.  `cond 2p` is the ANSWER the routine returns and `N^2/2` is what it")
println("pays to get there.  If the answer saturates while the price keeps climbing, the cost is entirely in a")
println("temporary -- which is the shape item 29 set out to remove from the cache and did not remove from here.")
