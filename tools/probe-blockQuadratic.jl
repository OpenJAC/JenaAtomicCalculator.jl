#
# probe-blockQuadratic.jl   --   IS THE EOL PAIR CACHE QUADRATIC IN THE SYMMETRY-BLOCK SIZE?
#
# WHY THIS MATTERS MORE THAN THE COST LAW IT CORRECTS.  Two linear fits of EOL peak memory disagree by an order
# of magnitude: the 04-Sep-2026 law (1.58 + 1.74e-3*nCsf GB) was fitted on spaces of thousands of CSFs, and a
# re-fit of 14-Sep over 14 to 1231 CSFs gave a slope 24x smaller.  Measured peaks: 887 CSFs cost 1.6 GB and
# 25 085 CSFs cost about 26 GB -- 28x the CSFs for 16x the memory overall, but against a ~1.5 GB fixed baseline
# the SPACE-DEPENDENT part goes from ~0.1 to ~25 GB, a factor of 250.  That is not a slope either fit got wrong;
# it is the wrong FUNCTIONAL FORM.  Two linear fits to different parts of a parabola will always disagree.
#
# THE MECHANISM TO TEST.  `PairCoefficientCache` is built per SYMMETRY BLOCK and indexed by ORDERED PAIRS, so
# its ptr arrays alone are n^2 entries for a block of n CSFs, and the interned label/value arrays grow with the
# number of non-empty pairs, which is also O(n^2).  Item 29 reduced the CONSTANT in front of that n^2 -- the
# struct's own docstring says the gain grows with block size -- but nothing in it changed the EXPONENT.  If the
# cost is quadratic in the LARGEST BLOCK rather than linear in the TOTAL CSF count, then a law written per CSF
# cannot be right at any size, and the quantity to warn on is the block.
#
# WHY NO SCF IS NEEDED, AND WHY ONE SPACE IS ENOUGH.  The pair coefficients are orbital-independent, so the
# cache can be built from a bare `Basics.generateBasis` with no solve at all.  And a single RAS space already
# contains several blocks of DIFFERENT sizes -- one per J -- so the scaling curve comes out of one run, with
# every point sharing the same ion, core, layer and subshell set.  That removes the confounding that made the
# earlier series ambiguous: here nothing varies but the block size.
#
# READ THE LAST COLUMN.  MB/n^2 constant across the blocks means quadratic; MB/n constant would mean linear.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

Z     = 22.0
refs  = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s"),
         Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")]
syms  = [LevelSymmetry(J, Basics.plus)  for J = 0:6]
grid  = Basics.recommendedGrid(refs, Nuclear.Model(Z); printout=false)
core  = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
valp  = [Shell("3p"), Shell("3d"), Shell("4s")]

layer = length(ARGS) >= 1 ? ARGS[1] : "3p-v"
to    = layer == "3p-v"  ? [Shell("3d"), Shell("4s")] :
        layer == "3p-4p" ? [Shell("3d"), Shell("4s"), Shell("4p")] :
                           [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d")]

step  = RasStep(RasStep(); seFrom=valp, seTo=to, deFrom=valp, deTo=to, frozen=core)
basis = Basics.generateBasis(refs, syms, step)
println("\nlayer $layer :  $(length(basis.csfs)) CSFs over $(length(basis.subshells)) subshells")

# The orbitals are irrelevant to the angular coefficients but the struct wants them; any set will do.
blocks = unique( [ LevelSymmetry(csf.J, csf.parity)  for csf in basis.csfs ] )
sort!(blocks, by = s -> -count(c -> LevelSymmetry(c.J, c.parity) == s, basis.csfs))

# THE SPARSITY COLUMNS ARE THE POINT OF THE SECOND HALF OF THIS PROBE.  A CI matrix element vanishes unless the
# two CSFs differ by at most TWO orbitals (Slater-Condon), so in a layer of single AND double excitations two
# CSFs that are each doubly excited can differ by four and contribute strictly nothing.  The cache already sees
# this -- 61 % of ordered pairs carried nothing on one 239-CSF block -- but its `ptr` array allocates a slot for
# every ordered pair regardless, i.e. 4n^2 bytes per block whether the pairs exist or not.
#
# SO `fill %` IS THE DIAGNOSTIC.  If it FALLS as the block grows, the non-empty pairs grow slower than n^2 and
# the measured quadratic is partly the dense INDEX rather than the physics -- in which case making the pair index
# sparse (item 29 interned the values but left the index dense) is the fix, and it gets cheaper the larger the
# block, which is where it is needed.  If `fill %` is flat, the pairs really are quadratic and only hermiticity
# and recomputation are left.
@printf("\n%-10s %8s %14s %14s %14s %10s %12s %12s\n",
        "block", "n", "cache MB", "MB / n", "MB / n^2", "fill %", "ptr MB", "entries")
println("-"^104)
tot = 0.0
for  sym  in  blocks
    n = count(c -> LevelSymmetry(c.J, c.parity) == sym, basis.csfs)
    n == 0  &&  continue
    c  = SC.cacheCsfPairCoefficientsEOL(sym, basis)
    m  = Base.summarysize(c) / 1024^2
    global tot += m
    nPair  = length(c.ptr2p) - 1
    filled = count( i -> c.ptr2p[i+1] > c.ptr2p[i], 1:nPair ) +
             count( i -> c.ptr1p[i+1] > c.ptr1p[i], 1:(length(c.ptr1p)-1) )
    nEntry = (length(c.lab1p) + length(c.lab2p))
    ptrMB  = 4.0 * (length(c.ptr1p) + length(c.ptr2p)) / 1024^2
    @printf("%-10s %8d %14.4f %14.3e %14.3e %10.1f %12.4f %12d\n",
            string(sym), n, m, m/n, m/n^2, 100*filled/(2*nPair), ptrMB, nEntry)
    flush(stdout)
end
@printf("\ntotal cache over all blocks: %.3f MB\n", tot)
println("MB/n^2 roughly constant down the column => QUADRATIC in the block size, and a per-CSF law is the wrong")
println("shape.  MB/n roughly constant instead => linear, and the disagreement between the two fits is elsewhere.")
println("\nTHEN READ `fill %` AGAINST n.  Falling => the non-empty pairs grow slower than n^2 and the dense pair")
println("index is paying for pairs the physics does not have;  flat => the pairs are genuinely quadratic.")
