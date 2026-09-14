#
# probe-coefficientsPerPair.jl   --   HOW MANY ANGULAR COEFFICIENTS A CSF PAIR CARRIES, and what that depends on.
#
# WHY THIS EXISTS.  An EOL/RAS layer's memory is dominated by what is stored PER CSF PAIR within a symmetry block,
# so the natural predictor is  SUM_b n_b^2  x  (coefficients per pair), with b running over the J^P blocks.  On
# 14-Sep-2026 four Ti III anchors were measured -- 887 -> 1.647 GB, 3167 -> 2.149, 11 298 -> 12.784, 25 085 ->
# 29.73 -- and NO law in the CSF count alone fits them: the space-dependent cost per CSF is 1.7e-4, 2.1e-4, 1.0e-3,
# 1.1e-3 GB, a FIVEFOLD JUMP between the second and third anchor and then flat.
#
# THE MAINTAINER'S READING, WHICH THIS FILE TESTS.  Those four are not four sizes of one thing: they are layers
# that add 4p, then 4d, then 4f, and the jump lands exactly where 4d enters.  The number of two-particle
# coefficients a pair carries depends on which RANKS k are allowed, and k runs to 2l -- so p contributes k <= 2,
# d <= 4, f <= 6.  A law in the CSF count must then fail across p-, d- and f-layers, which is what it did.
#
# BUT l ALONE IS NOT ENOUGH, AND THAT IS THE SECOND QUESTION HERE.  The rank RANGE is set by l, but the NUMBER of
# contributing terms is not: a diagonal element sums over the electron pairs of the CSF, and d^5 has ten
# open-shell pairs where d^3 has three.  So the count should grow with OCCUPATION as well, and `c(l_max)` is too
# crude a predictor.
#
# AND THE THIRD QUESTION SEPARATES THE TWO EFFECTS: d^5 against p^3 d^2 carries the same five open-shell electrons,
# but spreads them over two shells of lower individual l (p-p gives k <= 2, p-d k <= 3, d-d k <= 4) while offering
# MORE distinct subshell quadruples.  Which effect wins is not something to reason out; it is measured below.
#
# NO SCF IS NEEDED -- the pair coefficients are orbital-independent by construction, which is why they are cached
# once and reused across iterations -- so every row here costs only the angular work.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

const CORE = "1s^2 2s^2 2p^6 "
cases = [ ("p^2      3p^2",            CORE*"3p^2"),
          ("p^3      3p^3",            CORE*"3p^3"),
          ("p^4      3p^4",            CORE*"3p^4"),
          ("d^2      3d^2",            CORE*"3s^2 3p^6 3d^2"),
          ("d^3      3d^3",            CORE*"3s^2 3p^6 3d^3"),
          ("d^4      3d^4",            CORE*"3s^2 3p^6 3d^4"),
          ("d^5      3d^5",            CORE*"3s^2 3p^6 3d^5"),
          ("mixed    3p^3 3d^2",       CORE*"3s^2 3p^3 3d^2"),
          ("mixed    3p^4 3d^1",       CORE*"3s^2 3p^4 3d"),
          ("f^2      4f^2",            CORE*"3s^2 3p^6 3d^10 4s^2 4p^6 4f^2"),
          ("f^3      4f^3",            CORE*"3s^2 3p^6 3d^10 4s^2 4p^6 4f^3") ]

@printf("\n%-22s %7s %6s %11s %12s %12s %11s %11s\n",
        "configuration", "nCsf", "blocks", "sum n_b^2", "pairs != 0", "2p entries", "per pair", "fill %")
println("-"^104)

for (tag, cstr) in cases
    conf  = Configuration(cstr)
    rconfs = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
    subsh  = Basics.generateSubshellList(rconfs)
    csfs   = CsfR[]
    for rc in rconfs   append!(csfs, Basics.generateCsfRs(rc, subsh))   end
    nEl    = sum(csfs[1].occupation)
    basis  = Basis(true, nEl, subsh, csfs, Subshell[], Dict{Subshell,Orbital}())
    blocks = unique( [ LevelSymmetry(c.J, c.parity) for c in csfs ] )
    tot2p = 0;   totPairs = 0;   totFilled = 0;   s2 = 0.0
    for sym in blocks
        n = count(c -> LevelSymmetry(c.J, c.parity) == sym, csfs)
        n == 0 && continue
        s2 += Float64(n)^2
        c  = SC.cacheCsfPairCoefficientsEOL(sym, basis)
        np = length(c.ptr2p) - 1
        filled = count(i -> c.ptr2p[i+1] > c.ptr2p[i], 1:np)
        tot2p += length(c.lab2p);   totPairs += np;   totFilled += filled
    end
    @printf("%-22s %7d %6d %11.4g %12d %12d %11.2f %11.1f\n",
            tag, length(csfs), length(blocks), s2, totFilled, tot2p,
            tot2p/max(totFilled,1), 100*totFilled/max(totPairs,1))
    flush(stdout)
end
println("\nREAD THE `per pair` COLUMN.  If it is set by l alone, p^2/p^3/p^4 agree with each other and d^n with each")
println("other;  if it also grows with OCCUPATION, it rises down each family.  And d^5 against 3p^3 3d^2 says whether")
println("spreading the same five electrons over two lower-l shells costs more or less than keeping them in one.")
