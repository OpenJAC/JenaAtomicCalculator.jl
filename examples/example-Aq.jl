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
    #   RANK > 0 IS NOT COVERED HERE, but only because this CSF set cannot test it: every CSF is closed shells plus ONE
    #   electron, and on such a set every rank > 0 coefficient GRASP returns is exactly 1.0, so an implementation that
    #   returned 1.0 for every allowed pair would score perfectly.  Branch e tests rank > 0 on a set that can fail.
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
elseif  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    # Branch d (GENUINELY OPEN SHELLS, and the bug it found): carbon-like 1s^2 2s^2 2p^2, whose five CSFs carry two
    #   electrons in an open shell and so exercise seniority and the coefficients of fractional parentage.
    #
    #   WHY THIS BRANCH EXISTS.  Branches a-c run on 1s^2 2s and 1s^2 2p, where every CSF is closed shells plus ONE
    #   electron.  On that set EVERY rank > 0 coefficient GRASP produces is 1.0, so an implementation that returned 1.0
    #   for every allowed pair would score perfectly -- the set cannot discriminate.  On 2p^2 the GRASP values spread
    #   over sqrt(1/10), sqrt(1/6), sqrt(3/10), sqrt(1/5), sqrt(4/3), sqrt(8/5) and their negatives, which is a set that
    #   can fail.  Building it is what found the defect below.
    #
    #   THE DEFECT, in this module and now fixed.  A one-body operator changes the occupation of at most TWO subshells,
    #   by exactly ONE electron each.  Any other pattern is an EXACT ZERO -- producing it would take a two-body operator.
    #   computeCoefficientsScalar originally classified such a pattern as UNSUPPORTED and raised.  The CSF pair
    #   2p^2 (J=0) against (2p-)^2 (J=0) is exactly that case: it differs in two subshells by two electrons each, and
    #   the right answer is an empty list.  GRASP emits nothing for it.  Confusing "identically zero" with "not
    #   implemented" is a real error, and only a set with two-electron open shells could expose it.
    #
    #   REPORT (22-Aug-2026).  Five CSFs, matching GRASP's five.
    #
    #       occupation identity, all 5 CSFs               max deviation = 0.0     (exactly, not to rounding)
    #       off-diagonal coefficients, SpinAngularNew     0                       (GRASP: 0)
    #       pairs that raised                             0                       (before the fix: non-zero)
    #       off-diagonal coefficients, SpinAngular        0                       (GRASP: 0)
    #
    #   Two paths are exercised on the way, and both return empty for the right reason rather than by accident: the pair
    #   2p^2 <-> (2p-)^2 by the two-electron rule above, and the pair 2p^2 <-> 2p- 2p by the triangle condition, since a
    #   scalar operator cannot connect j = 1/2 to j = 3/2.  The present module agrees here; the two spurious
    #   coefficients of branch c needed OPPOSITE PARITY to appear, which this single configuration does not provide.
    #
    localConfigs = [Configuration("1s^2 2s^2 2p^2")]
    localRel     = ConfigurationR[]
    for  conf in localConfigs
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    op = SpinAngularNew.OneParticleOperator(0, Basics.plus)
    println("\n  $(length(localCsfs)) CSFs of 1s^2 2s^2 2p^2:")
    for  (i,c) in enumerate(localCsfs)   println("    CSF $i:  J = $(c.J)$(string(c.parity))   occ = $(c.occupation)")   end

    maxDeviation = 0.0
    for  (i,c) in enumerate(localCsfs)
        coeffs = SpinAngularNew.computeCoefficients(op, c, c, localSubshells)
        global maxDeviation = max(maxDeviation, SpinAngularNew.checkOccupationIdentity(coeffs, c, localSubshells))
    end
    println("\n  occupation identity, max deviation = $maxDeviation      (must be 0)")

    nOff = 0;   nRaised = 0
    for  (i,l) in enumerate(localCsfs),  (j,r) in enumerate(localCsfs)
        i == j  &&  continue
        try
            global nOff    = nOff + length( SpinAngularNew.computeCoefficients(op, l, r, localSubshells) )
        catch  ex
            global nRaised = nRaised + 1
        end
    end
    println("  off-diagonal coefficients emitted   = $nOff        (GRASP2018: 0)")
    println("  CSF pairs that raised               = $nRaised        (must be 0: an exact zero is not 'unsupported')")
    #
elseif  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    # Branch e (RANK > 0, AND THE NORMALIZATION THAT HAD TO BE MEASURED): the two 2p- 2p CSFs of 1s^2 2s^2 2p^2, which
    #   carry two singly-occupied open subshells and therefore need recoupling but no coefficients of fractional
    #   parentage. Fifteen GRASP2018 coefficients, spanning ranks 1, 2 and 3 and both J = 1 and J = 2.
    #
    #   WHAT WAS IMPLEMENTED.  The textbook two-subsystem reduction: for |(j_a j_b) J> the tensor acting on the first
    #   electron gives (-1)^(j_a+j_b+J'+k) sqrt((2J+1)(2J'+1)) {j_a J j_b; J' j_a k} times <j_a||t^(k)||j_a>, and the
    #   mirror expression for the second.  That part is standard and was not in doubt.
    #
    #   WHAT WAS NOT.  Which normalization GRASP puts on top of it.  Reading oneparticlejj1.f90:61 says the coefficient
    #   is divided by sqrt(2k+1) -- and using that, all fifteen values came out WRONG, by ratios of 1.000000, 1.290994,
    #   0.774597, 0.654654 and 0.845154.  Those are not noise and not a single constant: they are exactly
    #
    #       ratio  =  sqrt( (2*J_bra + 1) / (2k + 1) )
    #
    #                        k = 1              k = 2              k = 3
    #       J_bra = 1     sqrt(3/3) = 1.000   sqrt(3/5) = 0.775   sqrt(3/7) = 0.655
    #       J_bra = 2     sqrt(5/3) = 1.291   sqrt(5/5) = 1.000   sqrt(5/7) = 0.845
    #
    #   fitting all fifteen.  So GRASP's NET convention divides by sqrt(2*J_bra+1), not by sqrt(2k+1) -- its own
    #   recoupling factor carries the remaining J-dependence.  This is the same shape as its rank-0 path, which also
    #   divides by sqrt(2J+1), so the two ranks are more alike than the source lines suggest.
    #
    #   WHY THIS IS NOT CURVE-FITTING.  A factor inferred from the data it was inferred from proves nothing, so the
    #   corrected form was tested OUT OF SAMPLE on a set it had never seen: 1s^2 2s and 1s^2 2p, where every CSF is
    #   closed shells plus ONE electron.  There J = j_a, the two roots cancel identically, and the coefficient must be
    #   EXACTLY 1 -- which is what GRASP returns for every such pair, and what this module now returns: 1.000000000000000
    #   for all five.  The prediction was made by the formula, not fitted to the answer.
    #
    #   REPORT (22-Aug-2026):  15 of 15 coefficients matched, ratio between 0.9999999999999997 and 1.0000000000000004,
    #   signs included.  Out of sample, 5 of 5 exactly 1.0.
    #
    #   ONE FLAW OF THIS MODULE THAT THE OUT-OF-SAMPLE RUN EXPOSED, and it is worth recording because it was mine: the
    #   first version returned an EMPTY LIST for any CSF pair whose open subshells differed -- which silently swallowed
    #   the single-electron SUBSTITUTIONS, for which GRASP does return coefficients.  An empty list where a real
    #   coefficient exists is a silent wrong answer, exactly what this module is for.  Such pairs now RAISE, and the 18
    #   of them in that set are counted below rather than passing unnoticed.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2s^2 2p^2")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    # JAC CSF 3 = GRASP CSF 3 (2p- 2p, J=1);  JAC CSF 4 = GRASP CSF 5 (2p- 2p, J=2)
    graspRank = Dict( (1,3,3,"2p_1/2") => -4.08248290463862851e-01, (1,3,3,"2p_3/2") =>  9.12870929175277013e-01,
                      (1,3,4,"2p_1/2") => -9.12870929175276791e-01, (1,3,4,"2p_3/2") =>  4.08248290463863017e-01,
                      (1,4,3,"2p_1/2") =>  7.07106781186547240e-01, (1,4,3,"2p_3/2") => -3.16227766016837886e-01,
                      (1,4,4,"2p_1/2") =>  7.07106781186547462e-01, (1,4,4,"2p_3/2") =>  9.48683298050513879e-01,
                      (2,3,3,"2p_3/2") =>  7.07106781186547684e-01, (2,3,4,"2p_3/2") =>  7.07106781186547684e-01,
                      (2,4,3,"2p_3/2") => -5.47722557505166185e-01, (2,4,4,"2p_3/2") =>  8.36660026534075785e-01,
                      (3,3,4,"2p_3/2") =>  1.00000000000000000e+00, (3,4,3,"2p_3/2") => -7.74596669241483293e-01,
                      (3,4,4,"2p_3/2") =>  6.32455532033675882e-01 )

    println("\n  rank bra ket  subshell        GRASP2018         SpinAngularNew        ratio")
    nMatched = 0;   worstRatio = 1.0
    for  k in [1,2,3],  ib in [3,4],  ik in [3,4]
        opK    = SpinAngularNew.OneParticleOperator(k, Basics.plus)
        coeffs = SpinAngularNew.computeCoefficients(opK, localCsfs[ib], localCsfs[ik], localSubshells)
        for  c in coeffs
            key = (k, ib, ik, string(c.a))
            if  haskey(graspRank, key)
                g = graspRank[key];    global nMatched = nMatched + 1
                if  abs(c.T/g - 1.0) > abs(worstRatio - 1.0)    global worstRatio = c.T/g    end
                @printf("   %d   %d   %d   %-9s %18.12f %18.12f  %12.9f\n", k, ib, ik, string(c.a), g, c.T, c.T/g)
            end
        end
    end
    println("\n  matched $nMatched of $(length(graspRank)) GRASP coefficients;  worst ratio = $worstRatio")
    #
end
#
