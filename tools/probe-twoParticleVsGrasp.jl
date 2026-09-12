#
# probe-twoParticleVsGrasp.jl   --   the JAC half of the GRASP2018 comparison behind priority item 27.
#
# Dumps JAC's TWO-PARTICLE spin-angular coefficients, converted to GRASP's plain-Slater (R^k) convention by
# `SpinAngular.toGraspCoulomb`, for every CSF pair of {3s3p, 3p3d} at J = 0-, 1-, 2- with no core -- the same
# 12 CSFs GRASP's `rcsfgenerate` produces from the same input.  The GRASP half is
# `tools/grasp-drivers/angdump2.f90`, which calls RKCO_GG and prints the identical eight columns.
#
# MEASURED 12-Sep-2026: 44 of the 52 CSF pairs agree exactly; 8 DISAGREE, and every one of those is a CROSS pair
# between a 3s3p CSF and a 3p3d CSF -- the interaction that inverts the 3P* fine structure (item 27, and
# `tools/probe-fineStructureSpectator.jl` for the spectroscopic symptom).  Every pair WITHIN one configuration
# agrees, k = 1 included, so the disagreement is not a convention mismatch in the conversion.  Some cross pairs
# differ by an exact sign; others omit a k = 1 coefficient outright and carry a k = 2 coefficient wrong by
# sqrt(2/7) or sqrt(14)/3.
#
# THE CSF ORDERS DIFFER between the two codes and must be matched before comparing -- each CSF here is identified
# by (occupied relativistic subshells, 2J), which is unique for these two-electron configurations.  The mapping
# measured on that day was JAC -> GRASP = 3, 8, 1, 4, 5, 9, 2, 6, 10, 11, 7, 12; re-derive it rather than trust
# it, since it depends on the order `rcsfgenerate` happens to emit.
#
# Building the oracle: source `/home/fritzsch/fri/grasp/grasp-master` (read-only, Rule 6; copy elsewhere to
# build).  Compile each library with `-std=legacy -fallow-argument-mismatch` in repeated passes until no new
# object appears, in the order libmod, lib9290, libmcp90, librang90; link without iniest2, maneig, dvdson and
# spodmv, which need a LAPACK this machine does not have.
#
using JenaAtomicCalculator, Printf
const SA = JenaAtomicCalculator.SpinAngular

subsh = [Subshell(3,-1), Subshell(3,1), Subshell(3,-2), Subshell(3,2), Subshell(3,-3)]   # 3s 3p- 3p 3d- 3d
csfs  = CsfR[]
for  cs in ["3s 3p", "3p 3d"]
    for rc in Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(cs))
        append!(csfs, Basics.generateCsfRs(rc, subsh))
    end
end
keep = [ c for c in csfs if c.parity == Basics.minus && Basics.twice(c.J) in [0,2,4] ]
# label each CSF by (occupied subshell indices, 2J) so it can be matched to GRASP's list
function tag(c)
    occ = [ i for i = 1:length(subsh) if c.occupation[i] > 0 ]
    return( (occ, Basics.twice(c.J)) )
end
println("# JAC CSFs")
for (i, c) in enumerate(keep)
    (o, tj) = tag(c)
    @printf("# %2d  subshells %-8s  2J = %d\n", i, string(o), tj)
end
println("# ICSF JCSF   a    b    c    d   K  COEFF     (GRASP Coulomb convention)")
op = SA.TwoParticleOperator(0, Basics.plus)
for  (i, ci) in enumerate(keep), (j, cj) in enumerate(keep)
    tag(ci)[2] == tag(cj)[2]  ||  continue
    for cf in SA.computeCoefficients(op, ci, cj, subsh)
        g = SA.toGraspCoulomb(cf)
        ia = findfirst(==(g.a), subsh);  ib = findfirst(==(g.b), subsh)
        ic = findfirst(==(g.c), subsh);  id = findfirst(==(g.d), subsh)
        abs(g.V) > 1.0e-14 &&
            @printf("%5d%5d %4d %4d %4d %4d %3d %26.17e\n", i, j, ia, ib, ic, id, g.nu, g.V)
    end
end
