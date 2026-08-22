#
println("Aq) Compare the spin-angular coefficients of SpinAngularNew against SpinAngular and against GRASP2018.")

#
# SpinAngularNew is deliberately NOT included from JenaAtomicCalculator.jl while it is under development, so that a broken
# intermediate state cannot break the package. It is included here directly.
#
include("../src/module-SpinAngularNew.jl")

configs = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]

relconfList = ConfigurationR[]
for  conf in configs
    wax = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf);   append!(relconfList, wax)
end
subshellList = Basics.generateSubshellList(relconfList)
Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)

csfList = CsfR[]
for  relconf in relconfList
    newCsfs = Basics.generateCsfRs(relconf, subshellList);   append!(csfList, newCsfs)
end

if  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    # Branch a (AGAINST GRASP2018, the independent oracle): rank-0 one-particle coefficients for 1s^2 2s and 1s^2 2p,
    #   compared with the coefficients that GRASP2018's librang90 -- G. Gaigalas's library, the direct ancestor of
    #   SpinAngular -- produces for the same three CSFs.
    #
    #   WHY GRASP AND NOT THE PRESENT MODULE.  A comparison against SpinAngular can only establish that the two agree;
    #   it cannot establish that either is right, and SpinAngular has a documented normalization defect that shipped
    #   (module-Hfs.jl:370-379).  GRASP2018 is a genuinely independent implementation, and since SpinAngularNew adopts
    #   GRASP's convention this branch is an EQUALITY test rather than a conversion test.
    #
    #   HOW THE REFERENCE WAS OBTAINED, so that it can be reproduced.  GRASP2018 source (read-only, outside this
    #   repository) copied to a scratch directory and built there; the libraries must be made in the order
    #   libmod -> lib9290 -> libmcp90 -> librang90, and gfortran 13.3 needs -std=legacy -fallow-argument-mismatch.
    #   Only SETMC requires LAPACK, for four DLAMCH calls ('L','O','U','E'), which a ~30-line shim supplies.  A driver
    #   calls SETMC; SETCON; FACTT; SETCSLA(FNAME,NCORE) and then loops CSF pairs over ONESCALAR.  Two traps: FACTT is
    #   required and easy to omit, and SETCSLA takes CHARACTER(LEN=24) and builds its file name via INDEX(NAME,' '),
    #   so a length-4 literal 'rcsf' yields an EMPTY name and it silently reads fort.21 instead.
    #
    #   REPORT (22-Aug-2026).  EXACT AGREEMENT on every coefficient, and the new module is the more accurate of the two.
    #
    #       CSF pair          subshell      GRASP2018                SpinAngularNew
    #        (1,1)             1s_1/2       2.00000000000000044      2.0
    #        (1,1)             2s_1/2       1.00000000000000000      1.0
    #        (2,2)             1s_1/2       2.00000000000000044      2.0
    #        (2,2)             2p_3/2       1.00000000000000000      1.0
    #        (3,3)             1s_1/2       2.00000000000000044      2.0
    #        (3,3)             2p_1/2       1.00000000000000000      1.0
    #
    #   (JAC's CSF ORDER differs from GRASP's -- JAC 2 <-> 3 here -- so the rows are matched by CSF content, not index.)
    #
    #   THE 4.4e-16 IS GRASP'S, NOT OURS.  GRASP reaches 2.0 through a chain of Float64 multiplications and square roots
    #   and lands 4.4e-16 above it; SpinAngularNew reaches it as an integer occupation and is exact.  So this branch also
    #   measures the accumulated rounding error of the reference, which is a free by-product of computing the quantity
    #   the exact way rather than the transcribed way.
    #
    #   RANK > 0 IS NOT COVERED HERE.  It is not yet implemented in SpinAngularNew and the call refuses rather than
    #   returning a number; the control that JAC and GRASP already agree to ~1e-15 at rank > 0 was measured separately
    #   and is what establishes that the comparison method itself is sound.
    #
    graspReference = Dict( (1,"1s_1/2") => 2.00000000000000044, (1,"2s_1/2") => 1.0,
                           (2,"1s_1/2") => 2.00000000000000044, (2,"2p_3/2") => 1.0,
                           (3,"1s_1/2") => 2.00000000000000044, (3,"2p_1/2") => 1.0 )

    op = SpinAngularNew.OneParticleOperator(0, Basics.plus)
    println("\n  CSF pair   subshell     GRASP2018              SpinAngularNew          difference")
    for  (ic, csf) in enumerate(csfList)
        coeffs = SpinAngularNew.computeCoefficients(op, csf, csf, subshellList)
        for  c in coeffs
            key = (ic, string(c.a))
            if  haskey(graspReference, key)
                gr = graspReference[key]
                @printf("   (%d,%d)     %-9s  %22.17f  %22.17f  %10.2e\n", ic, ic, string(c.a), gr, c.T, abs(gr - c.T))
            end
        end
    end
    #
elseif  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    # Branch b (THE OCCUPATION IDENTITY): the one exact identity a scalar one-particle operator must satisfy, asserted on
    #   every diagonal matrix element.
    #
    #   WHAT IT TESTS AND WHY IT CAN BE TRUSTED.  For any CSF and any coupling whatever,
    #
    #       <Psi| sum_i f(i) |Psi>  =  sum_a  N_a <a| f |a>
    #
    #   so in GRASP convention the rank-0 diagonal coefficient of each subshell is exactly its occupation number.  This is
    #   an IDENTITY, not a tolerance: the deviation is zero or the module is wrong.  It is also the practical argument for
    #   adopting GRASP's convention over a "plain" coefficient of N_a/sqrt(2j+1) -- the latter is not recognizably anything,
    #   so a caller with the wrong convention sees a plausible number, whereas here it sees 1.414 where it expects 2.
    #
    #   The predecessor module has NO direct test coverage at all -- not one of the 51 tests in test/runtests.jl mentions
    #   SpinAngular -- so a check that runs on real work rather than on a fixture is worth having.
    #
    #   REPORT (22-Aug-2026): deviation 0.0 for all three CSFs, i.e. zero exactly and not merely to rounding, because the
    #   coefficient is produced as an integer occupation rather than accumulated in floating point.
    #
    op = SpinAngularNew.OneParticleOperator(0, Basics.plus)
    maxDeviation = 0.0
    println("")
    for  (ic, csf) in enumerate(csfList)
        coeffs    = SpinAngularNew.computeCoefficients(op, csf, csf, subshellList)
        deviation = SpinAngularNew.checkOccupationIdentity(coeffs, csf, subshellList)
        global maxDeviation = max(maxDeviation, deviation)
        println("  CSF $ic  (occ = $(csf.occupation)):   deviation = $deviation")
    end
    println("  maximum deviation over all CSFs = $maxDeviation      (must be 0)")
    #
elseif  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    # Branch c (AGAINST THE PRESENT MODULE): the same rank-0 coefficients from SpinAngular, so that the convention change
    #   is recorded as a number and so that a later change to either module shows up here.
    #
    #   TWO DIFFERENCES ARE EXPECTED, AND BOTH APPEAR.
    #
    #   (1) THE NORMALIZATION, sqrt(2j_a+1).  SpinAngular's computeCoefficientsScalar has the line
    #       `# wa = wa * sqrt(Basics.twice(ji) + 1.)` COMMENTED OUT, while the identical line is ACTIVE in
    #       computeCoefficientsNonScalar under the comment "GRASP like".  So its rank-0 coefficients are smaller than
    #       GRASP's by exactly sqrt(2j_a+1).  MEASURED, and it is exact on both j values present:
    #
    #           subshell    SpinAngular      SpinAngularNew    ratio      sqrt(2j+1)
    #           1s_1/2      1.4142135624     2.0               1.414214   1.414214
    #           2s_1/2      0.7071067812     1.0               1.414214   1.414214
    #           2p_1/2      0.7071067812     1.0               1.414214   1.414214
    #           2p_3/2      0.5              1.0               2.000000   2.000000
    #
    #       NEITHER MODULE IS WRONG ON ITS OWN.  module-Hamiltonian.jl:277 re-applies sqrt(jj+1) at the call site, and
    #       IsotopeShift.amplitude does the same, so JAC is self-consistent; Hfs.amplitude instead DIVIDES the factor out
    #       for rank > 0.  The same physics is expressed two opposite ways and every caller must know which -- which is
    #       why four modules have guessed wrong, one of them shipping a hyperfine constant too large by sqrt(2).
    #
    #   (2) TWO COEFFICIENTS THAT SHOULD NOT EXIST.  SpinAngular emits rank-0 coefficients for the CSF pairs (1,3) and
    #       (3,1), connecting 2s_1/2 to 2p_1/2 -- that is kappa = -1 to kappa = +1, and across OPPOSITE PARITY.  A scalar
    #       operator cannot do either.  GRASP2018 emits nothing for those pairs, and neither does SpinAngularNew, which
    #       decides it in SpinAngularNew.isAllowed1p from the triangle and parity conditions before any float is computed.
    #
    #       This is LATENT rather than harmful today: the Hamiltonian only pairs CSFs within one symmetry block, so those
    #       pairs never arise there, and the radial integral would vanish in any case.  It is recorded because it is
    #       exactly what a magnitude threshold cannot catch and a selection rule can -- `abs(wa) >= 2.0e-10` asks whether
    #       a number came out small, where the question is whether the quantum numbers permit it at all.
    #
    op    = SpinAngularNew.OneParticleOperator(0, Basics.plus)
    opOld = SpinAngular.OneParticleOperator(0, Basics.plus, true)
    println("\n  rank ICSF JCSF  a         b          SpinAngular        SpinAngularNew       ratio")
    for  (ic, l) in enumerate(csfList),  (ir, r) in enumerate(csfList)
        oldCoeffs = SpinAngular.computeCoefficients(opOld, l, r, subshellList)
        newCoeffs = SpinAngularNew.computeCoefficients(op, l, r, subshellList)
        for  oc in oldCoeffs
            if  abs(oc.T) < 1.0e-14    continue    end
            idx = findfirst(nc -> nc.a == oc.a  &&  nc.b == oc.b, newCoeffs)
            if  idx === nothing
                @printf("     0   %2d   %2d  %-9s %-9s %18.10f   %18s\n", ic, ir, string(oc.a), string(oc.b), oc.T,
                        "(none -- forbidden)")
            else
                nc = newCoeffs[idx]
                @printf("     0   %2d   %2d  %-9s %-9s %18.10f   %18.10f  %10.6f\n", ic, ir, string(oc.a), string(oc.b),
                        oc.T, nc.T, nc.T/oc.T)
            end
        end
    end
    #
end
#
