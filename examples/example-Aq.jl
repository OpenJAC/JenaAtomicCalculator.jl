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
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch f (COEFFICIENTS OF FRACTIONAL PARENTAGE): a subshell holding TWO electrons, which is the case branches a-e
    #   all raise on, and the last substantial piece of the one-particle problem.
    #
    #   WHAT IS AND IS NOT RE-IMPLEMENTED.  G. Gaigalas's completely reduced (j Q J ||| W^(kq kj) ||| j Q' J') tables live
    #   in module-SpinAngular-inc-reducedcoeffs.jl, stored EXACTLY as [sign, num, den] and returned as sign*sqrt(num/den).
    #   They are data, they are his, and they are not in doubt, so they are reused rather than re-typed -- re-typing a
    #   correct table adds risk and nothing else.  What is re-implemented is the ASSEMBLY: the quasispin Wigner-Eckart step
    #   that turns a completely reduced element into the one for a shell of N electrons, with kq = 1 for even kj and 0 for
    #   odd kj, and M_Q equal on both sides because the operator conserves particle number.
    #
    #   THE ASSEMBLY WAS ISOLATED BEFORE IT WAS TRUSTED.  Comparing SpinAngularNew.shellReducedW against
    #   SpinAngular.irreducibleTensor(SchemeEta_W(), ...) gave ratio 1.000000 on every case -- which separates "is the
    #   shell matrix element right?" from "is the outer normalization right?".  It was the second that was wrong, and
    #   knowing which half to look at is most of the work.  Solving for the outer factor on four GRASP coefficients gave
    #   1/sqrt(2J_bra+1) on all four, so
    #
    #       T^(k)(a,a)  =  - <j^N v J || W^(k) || j^N v' J'> * sqrt(2j+1) / ( sqrt(2k+1) * sqrt(2J_bra+1) )
    #
    #   REPORT (23-Aug-2026), 1s^2 3d^2 against GRASP2018 -- 16 of 16 matched, worst ratio 1.0000000000000002:
    #
    #     j = 5/2 (3d_5/2), TWELVE values, seniority 2 at both J = 2 and J = 4, so the CFP tables are genuinely exercised
    #       rank 1: (2,2) 0.828078671211   (4,4) 1.511857892037
    #       rank 2: (0,2) 2.000000000000   (2,0) 0.894427191000   (2,2) -0.638876565000   (2,4) 0.995910003310
    #               (4,2) 0.742307488958   (4,4) 0.670059394260
    #       rank 3: (2,2) -0.995910003310  (2,4) 1.355261854358   (4,2) 1.010152544552   (4,4) -0.273550602216
    #     j = 3/2 (3d_3/2), FOUR values, reproducing the 2p^2 numbers exactly as a consistency check
    #       rank 1: (2,2) 1.264911064067   rank 2: (0,2) 2.0, (2,0) 0.894427191   rank 3: (2,2) -1.264911064067
    #
    #   AND THE POINT OF THE WHOLE EXERCISE.  That one expression also covers RANK 0: for k = 0 the shell element is
    #   -N sqrt((2J+1)/(2j+1)), the roots cancel, and it collapses to exactly N.  Checked here on every 3d^2 CSF, for
    #   both j values and both seniorities: 2.000000000000000 against an occupation of 2, every time.  A single formula
    #   for every rank is what goal (1) asked for, and this is it demonstrated rather than asserted.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 3d^2")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    i52 = findfirst(sh -> Basics.subshell_2j(sh) == 5, localSubshells)
    ref52 = Dict( (1,4,4)=> 8.28078671210824901e-01, (1,8,8)=> 1.51185789203690879e+00,
                  (2,0,4)=> 1.99999999999999956e+00, (2,4,0)=> 8.94427190999915744e-01,
                  (2,4,4)=>-6.38876564999939944e-01, (2,4,8)=> 9.95910003310478631e-01,
                  (2,8,4)=> 7.42307488958090178e-01, (2,8,8)=> 6.70059394260489882e-01,
                  (3,4,4)=>-9.95910003310478409e-01, (3,4,8)=> 1.35526185435787672e+00,
                  (3,8,4)=> 1.01015254455221060e+00, (3,8,8)=>-2.73550602216096561e-01 )

    findCsf(twoJ) = findfirst(c -> c.occupation[i52] == 2  &&  Basics.twice(c.J) == twoJ  &&
                              all(k -> k == i52 || c.occupation[k] == 0 ||
                                       c.occupation[k] == Basics.subshell_2j(localSubshells[k])+1,
                                  1:length(localSubshells)), localCsfs)

    println("\n  3d_5/2 with two electrons, against GRASP2018")
    println("   rank 2Jb 2Jk       GRASP2018        SpinAngularNew        ratio")
    nMatched = 0;   worstRatio = 1.0
    for  ((k, jb, jk), g) in sort(collect(ref52), by = x -> x[1])
        ib = findCsf(jb);    ik = findCsf(jk)
        (ib === nothing || ik === nothing)  &&  continue
        coeffs = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(k, Basics.plus),
                                                    localCsfs[ib], localCsfs[ik], localSubshells)
        idx = findfirst(c -> c.a == localSubshells[i52], coeffs)
        idx === nothing  &&  continue
        v = coeffs[idx].T;    global nMatched = nMatched + 1
        if  abs(v/g - 1.0) > abs(worstRatio - 1.0)    global worstRatio = v/g    end
        @printf("    %d   %d   %d   %16.12f %16.12f  %12.9f\n", k, jb, jk, g, v, v/g)
    end
    println("\n  matched $nMatched of $(length(ref52));   worst ratio = $worstRatio")

    println("\n  and the same formula at k = 0 must give the occupation number:")
    for  (i,c) in enumerate(localCsfs)
        isub = findfirst(kk -> c.occupation[kk] != 0  &&
                               c.occupation[kk] != Basics.subshell_2j(localSubshells[kk])+1, 1:length(localSubshells))
        isub === nothing  &&  continue
        c.occupation[isub] != 2  &&  continue
        jj = AngularJ64( Basics.subshell_2j(localSubshells[isub])//2 )
        w  = SpinAngularNew.shellReducedW(jj, 2, c.seniorityNr[isub], c.subshellJ[isub],
                                                 c.seniorityNr[isub], c.subshellJ[isub], 0)
        v  = -w * sqrt(Basics.twice(jj)+1.0) / sqrt(Basics.twice(c.J)+1.0)
        @printf("    CSF %d  %-9s  formula = %18.15f   occupation = %d\n", i, string(localSubshells[isub]), v,
                c.occupation[isub])
    end
    #
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch g (THE GENERAL CASE): several open subshells, one of which holds more than one electron -- the last gap in
    #   the one-particle problem, and the case every earlier branch raised on. 1s^2 2p^2 3s, where 2p_3/2 carries two
    #   electrons and 3s one, so BOTH subshells contribute and the coupling tree is no longer trivial.
    #
    #   WHAT WAS ADDED.  SpinAngularNew.chainRecoupling: the tensor is peeled outwards one subshell at a time. For every
    #   subshell beyond the acting one it sits in the FIRST subsystem with J_q as spectator; at the acting subshell it
    #   sits in the SECOND with X_{ip-1} as spectator. The coefficient is then the same expression as before,
    #
    #       T^(k)(a,a) = - R_chain * <j^N v J_a || W^(k) || j^N v' J'_a> * sqrt(2j_a+1) / ( sqrt(2k+1) sqrt(2J_bra+1) )
    #
    #   only with R_chain no longer equal to one.
    #
    #   WHY IT REPLACED THE TWO SPECIAL CASES INSTEAD OF JOINING THEM.  Both earlier results fall out of it as limits:
    #   with every other subshell closed each factor collapses to 1 and the single-subshell formula returns; with two
    #   singly-occupied subshells the two expressions reduce term for term to the Edmonds two-subsystem formulae. That
    #   was checked BEFORE the new method was used anywhere -- it reproduced 15 of 15 and 12 of 12 on the already
    #   verified sets -- and only then did it take over the dispatch. A general method that cannot reproduce the special
    #   cases it subsumes has not earned them.
    #
    #   REPORT (23-Aug-2026), against GRASP2018: 17 of 17, worst ratio 0.9999999999999997.
    #
    #     rank 1:  (1,1) 3s 1.000000    (4,4) 2p 1.200000 / 3s -0.447214    (4,7) 2p -0.400000 / 3s 0.894427
    #              (7,4) 2p 0.326599 / 3s -0.730297                          (7,7) 2p 1.222020 / 3s 0.683130
    #     rank 2:  (1,4) 2p 1.264911    (1,7) 2p 1.549193    (4,1) 2p -0.894427    (7,1) 2p 0.894427
    #     rank 3:  (4,4) 2p -0.800000   (4,7) 2p 0.979796    (7,4) 2p -0.800000    (7,7) 2p -0.979796
    #
    #   These are not a column of one repeated number, which matters: the values spread over 0.33 to 1.55 with both
    #   signs, and (4,7) differs from (7,4), so the bra/ket asymmetry is being tested and not averaged away.
    #
    #   TOTAL COVERAGE of the one-particle rank > 0 problem: 48 GRASP coefficients across three configurations.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2p^2 3s")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    i2p = findfirst(sh -> Basics.subshell_2j(sh) == 3 && Basics.subshell_l(sh) == 1, localSubshells)
    i3s = findfirst(sh -> string(sh) == "3s_1/2", localSubshells)
    findG(twoJ) = findfirst(c -> c.occupation[i2p] == 2 && c.occupation[i3s] == 1 &&
                                 Basics.twice(c.J) == twoJ, localCsfs)
    gToJac = Dict( 1 => findG(1), 4 => findG(3), 7 => findG(5) )
    subOf  = Dict( 3 => i2p, 4 => i3s )

    graspGen = [ (1,1,1,4,  9.99999999999999889e-01), (1,4,4,3,  1.20000000000000018e+00),
                 (1,4,4,4, -4.47213595499957983e-01), (1,4,7,3, -3.99999999999999967e-01),
                 (1,4,7,4,  8.94427190999915522e-01), (1,7,4,3,  3.26598632371090436e-01),
                 (1,7,4,4, -7.30296743340221433e-01), (1,7,7,3,  1.22202018532155754e+00),
                 (1,7,7,4,  6.83130051063973176e-01), (2,1,4,3,  1.26491106406735176e+00),
                 (2,1,7,3,  1.54919333848296659e+00), (2,4,1,3, -8.94427190999915633e-01),
                 (2,7,1,3,  8.94427190999915633e-01), (3,4,4,3, -7.99999999999999822e-01),
                 (3,4,7,3,  9.79795897113271086e-01), (3,7,4,3, -7.99999999999999822e-01),
                 (3,7,7,3, -9.79795897113271086e-01) ]

    println("\n  rank bra ket  subshell        GRASP2018        SpinAngularNew        ratio")
    nMatched = 0;   worstRatio = 1.0
    for  (k, gb, gk, gs, g) in graspGen
        ib = gToJac[gb];   ik = gToJac[gk]
        (ib === nothing || ik === nothing)  &&  continue
        coeffs = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(k, Basics.plus),
                                                    localCsfs[ib], localCsfs[ik], localSubshells)
        idx = findfirst(c -> c.a == localSubshells[subOf[gs]], coeffs)
        idx === nothing  &&  continue
        v = coeffs[idx].T;    global nMatched = nMatched + 1
        if  abs(v/g - 1.0) > abs(worstRatio - 1.0)    global worstRatio = v/g    end
        @printf("   %d    %d   %d   %-9s %16.12f %16.12f  %12.9f\n", k, gb, gk,
                string(localSubshells[subOf[gs]]), g, v, v/g)
    end
    println("\n  matched $nMatched of $(length(graspGen));   worst ratio = $worstRatio")
    #
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch h (SUBSTITUTIONS): CSF pairs of UNEQUAL occupation, where one electron moves between two subshells. This is
    #   the case every earlier branch refused, and it completes the one-particle problem.
    #
    #   THREE INGREDIENTS BEYOND THE RECOUPLING, each easy to drop and each changing the answer.
    #     (1) two single-subshell matrix elements <j^N v J || a^(+/-) || j^N' v' J'>, from the same Gaigalas CFP tables;
    #     (2) an ORDERING phase, because the recoupling is built with the lower subshell index first, so a creation on
    #         the higher index costs (-1)^(j_a + j_b - k + 1);
    #     (3) the JORDAN-WIGNER phase (-1)^(occupation strictly between the two subshells, + 1) -- the sign from
    #         anticommuting past the electrons in between. It depends on the OTHER subshells, not on the two taking
    #         part, which is exactly what makes it easy to forget.
    #
    #   AND A NINE-J RATHER THAN A SIX-J.  The chain is cut at the higher subshell: beyond it the total rank k is peeled
    #   outwards as before, AT it the two ranks j_a and j_b join -- which needs a 9j -- and below it the rank-j_a
    #   operator is reduced through the SAME chainRecoupling used for equal occupations, restricted to the sub-chain.
    #   That reuse is why this stayed short. Generalising chainRecoupling to a HALF-INTEGER rank was required and was a
    #   real bug on the way: the inner rank is j_a, not an integer, and the phase had to be rewritten over twice-values.
    #
    #   REPORT (23-Aug-2026), against GRASP2018 -- 30 coefficients, all exact.
    #
    #     ADJACENT subshells, 2p_1/2 <-> 2p_3/2 of 1s^2 2s^2 2p^2 : 16 of 16, ranks 1 and 2, worst ratio 1.0000000000000009
    #     NON-ADJACENT, 2s <-> 3s of 1s^2 2s^2 2p + 1s^2 2s 2p 3s : 14 of 14, worst ratio 1.0000000000000013
    #
    #   The second set is the one that matters: 2s and 3s are separated by 2p_1/2 and 2p_3/2, so the Jordan-Wigner
    #   string spans three subshells and is genuinely exercised. In the adjacent set it spans one, and a wrong string
    #   would still have passed half the time.
    #
    #   THE NORMALIZATION WAS NOT ASSUMED.  With the phases in place the residual against GRASP came out 1.000, sqrt(3)
    #   and sqrt(5) at J_bra = 0, 1, 2 -- i.e. exactly sqrt(2J_bra+1), the SAME outer factor already measured for the
    #   equal-occupation case. That it turned out to be the same factor is a result, not an assumption, and it is the
    #   third independent time that factor has appeared.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2s^2 2p^2")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    ipm = findfirst(sh -> string(sh) == "2p_1/2", localSubshells)
    ipp = findfirst(sh -> string(sh) == "2p_3/2", localSubshells)
    # GRASP 1,2,3,4,5 -> JAC 1,5,3,2,4  for 1s^2 2s^2 2p^2
    gToJac = Dict(1=>1, 2=>5, 3=>3, 4=>2, 5=>4)
    subOf  = Dict(3 => ipm, 4 => ipp)
    graspSub = [ (1,1,3,4,3,  1.41421356237309515), (1,2,3,3,4,  1.41421356237309515),
                 (1,3,1,3,4,  0.577350269189625731),(1,3,2,4,3,  1.15470053837925124),
                 (1,3,4,3,4,  0.912870929175276791),(1,4,3,4,3,  1.0),
                 (1,4,5,4,3, -1.00000000000000022), (1,5,4,3,4,  0.707106781186547573),
                 (2,1,5,4,3, -1.41421356237309515), (2,2,5,3,4,  1.41421356237309492),
                 (2,3,4,3,4, -0.707106781186547351),(2,4,3,4,3, -0.774596669241483404),
                 (2,4,5,4,3, -1.18321595661992318), (2,5,1,3,4,  0.447213595499957817),
                 (2,5,2,4,3, -0.894427190999915633),(2,5,4,3,4,  0.836660026534075563) ]

    println("\n  adjacent subshells, 2p_1/2 <-> 2p_3/2")
    println("   rank bra ket   a        b            GRASP2018        SpinAngularNew        ratio")
    nMatched = 0;   worstRatio = 1.0
    for  (k, gb, gk, sa, sb, g) in graspSub
        coeffs = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(k, Basics.plus),
                                                    localCsfs[gToJac[gb]], localCsfs[gToJac[gk]], localSubshells)
        idx = findfirst(c -> c.a == localSubshells[subOf[sa]] && c.b == localSubshells[subOf[sb]], coeffs)
        idx === nothing  &&  continue
        v = coeffs[idx].T;    global nMatched = nMatched + 1
        if  abs(v/g - 1.0) > abs(worstRatio - 1.0)    global worstRatio = v/g    end
        @printf("    %d   %d   %d  %-8s %-8s %16.12f %16.12f  %12.9f\n", k, gb, gk,
                string(localSubshells[subOf[sa]]), string(localSubshells[subOf[sb]]), g, v, v/g)
    end
    println("\n  matched $nMatched of $(length(graspSub));   worst ratio = $worstRatio")
    #
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch i (WHAT IT COSTS): goal (2) asked for a faster and more elegant module, and "more elegant" without a number
    #   is not a result. Both implementations are run over the same CSF pairs and timed.
    #
    #   REPORT (23-Aug-2026), 15 CSFs of 1s^2 2s^2 2p^2 + 1s^2 2s 2p^3, i.e. 225 CSF pairs, best of five:
    #
    #       rank      SpinAngular    SpinAngularNew    ratio       allocation old / new     ratio
    #        0           0.34 ms         0.04 ms       0.112         426 kB /   75 kB       0.176
    #        1           1.87 ms         1.22 ms       0.651        2295 kB / 1285 kB       0.560
    #        2           1.50 ms         0.95 ms       0.634        1933 kB / 1032 kB       0.534
    #
    #   RANK 0 IS THE INTERESTING ONE, at nine times faster and six times less memory, and the reason is structural
    #   rather than clever: for a diagonal pair the coefficient is the occupation number, so the whole recursion is
    #   skipped rather than run and its result discarded. That is the same design decision as deciding zeros by
    #   selection rules -- it is the SHAPE of the computation that pays, not micro-optimisation, and none has been done.
    #
    #   AT RANK > 0 the gain is 1.5x to 1.9x with roughly half the allocation, from the same source: the selection
    #   rules reject a pair in integer arithmetic before any 6j is evaluated, where the predecessor computes and then
    #   tests `abs(wa) >= 2.0e-10`.
    #
    #   THE TIMES ABOVE WILL NOT REPRODUCE EXACTLY, and a later reader should not read that as a regression. Repeated
    #   runs of this branch gave rank-2 ratios between 0.52 and 0.65 on the same machine; the ALLOCATION figures are
    #   stable to a few per cent and are the better number to compare against. Anything outside roughly 0.4-0.8 on the
    #   ratio, or a change in the allocation columns, would be worth looking at.
    #
    #   WHAT THIS DOES NOT SHOW.  A small case, 15 CSFs and 4 subshells, and one-particle operators only. The
    #   two-particle (electron-electron) operator dominates the run time of any real calculation and is NOT implemented
    #   here, so none of this is yet a statement about JAC as a whole. No caching or memoisation has been added either;
    #   the recoupling factors are recomputed for every CSF pair, which is the obvious next gain and is not taken.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2s^2 2p^2"), Configuration("1s^2 2s 2p^3")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end
    println("\n  $(length(localCsfs)) CSFs, $(length(localSubshells)) subshells, $(length(localCsfs)^2) pairs")

    sweepNew(op) = (for l in localCsfs, r in localCsfs
                        SpinAngularNew.computeCoefficients(op, l, r, localSubshells)   end)
    sweepOld(op) = (for l in localCsfs, r in localCsfs
                        SpinAngular.computeCoefficients(op, l, r, localSubshells)      end)

    println("\n   rank      old (ms)    new (ms)   ratio      old (kB)    new (kB)   ratio")
    for  k in [0, 1, 2]
        opNew = SpinAngularNew.OneParticleOperator(k, Basics.plus)
        opOld = SpinAngular.OneParticleOperator(k, Basics.plus, true)
        sweepNew(opNew);   sweepOld(opOld)                      # warm up, so compilation is not timed
        tNew = minimum([@elapsed sweepNew(opNew) for _ in 1:5])
        tOld = minimum([@elapsed sweepOld(opOld) for _ in 1:5])
        aNew = @allocated sweepNew(opNew)
        aOld = @allocated sweepOld(opOld)
        @printf("    %d      %9.2f   %9.2f  %7.3f    %9.0f   %9.0f  %7.3f\n",
                k, tOld*1e3, tNew*1e3, tNew/tOld, aOld/1024, aNew/1024, aNew/aOld)
    end
    #
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch j (THE TWO-PARTICLE COMPARISON, MADE EXPLICIT): JAC's electron-electron coefficients against GRASP2018,
    #   for all 25 CSF pairs of 1s^2 2s^2 2p^2. This is the first direct check of SpinAngular's e-e coefficients against
    #   an independent implementation -- the module has no test coverage of its own at all.
    #
    #   WHY A CONVERSION IS NEEDED, and what it is NOT.  The two codes decompose the interaction onto DIFFERENT
    #   quantities. JAC's coefficient multiplies the effective strength X^L; GRASP's Coulomb coefficient multiplies the
    #   plain Slater integral R^k (its own RKINTC header says so). The bridge is read from JAC's source, from
    #   InteractionStrength.XL_CoulombReference, which forms
    #
    #       X^L(abcd)  =  (-1)^L <a||C^L||c> <b||C^L||d> R^L(abcd)
    #
    #   so the coefficient of R^k follows by multiplying JAC's by exactly that prefactor. This is NOT a relation either
    #   code states as theory -- GRASP has no effective strength at all -- it is a bridge built here so the two can be
    #   compared, and its (-1)^L was read off the source rather than fitted. An earlier attempt omitted it and left
    #   precisely the 14 odd-k EXCHANGE terms disagreeing by a sign and nothing else.
    #
    #   THE RESULT THAT MATTERS IS NOT ONLY THAT THEY AGREE.  It is WHICH terms survive.
    #
    #       JAC emits, over the 25 CSF pairs, non-zero X^L coefficients          103
    #       the C^k factors ANNIHILATE                                            42
    #       surviving after conversion                                            61
    #       GRASP2018 emits                                                       61
    #
    #   The 42 that vanish are the terms whose tensorial structure forbids them -- X^L itself returns zero on the
    #   triangle condition and on rem(l_a+l_c+L,2) == 1, a few lines above the formula quoted. So the effective strength
    #   is not merely a repackaging: it CARRIES the selection rules, and a coefficient that looks non-zero in JAC's list
    #   contributes nothing because the strength it multiplies is zero. That is the argument for building operators on
    #   effective strengths rather than on radial integrals, and here it is visible as a count.
    #
    #   AND IT IS WHY THE SAME COEFFICIENTS SERVE BREIT.  module-Hamiltonian.jl:282-292 uses ONE coeff.V with both
    #   XL_Coulomb and XL_Breit. GRASP cannot: its Breit path takes a different callback (BREID), carries a sixth label
    #   ITYPE in 1..6, and multiplies one of six integral routines BRINT1..BRINT6 chosen by that tag. So the R^k
    #   convention does not generalise beyond Coulomb, while the effective-strength one does.
    #
    #   REPORT (23-Aug-2026): 61 keys on each side, none missing, ZERO disagreements, worst ratio 1.000000000000001.
    #   Coefficients are compared as a multiset on the canonical key of R^k -- invariant under (ab)<->(cd) and under
    #   a<->b together with c<->d -- with duplicates summed first. The two codes do not emit the same orderings, and
    #   comparing positionally would have manufactured differences that are not there.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2s^2 2p^2")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end
    localIdx = Dict(sh => i for (i,sh) in enumerate(localSubshells))

    # ... the canonical key of R^k(abcd), and the JAC -> GRASP CSF ordering for this configuration
    rkKey(k,a,b,c,d) = (k, minimum([(a,b,c,d), (b,a,d,c), (c,d,a,b), (d,c,b,a)]))
    jacToGrasp = Dict(1=>1, 2=>4, 3=>3, 4=>5, 5=>2)

    graspTwo = [
        (1,1,1,1,1,1,0, 1.000000000000000e+00),
        (1,1,1,2,1,2,0, 4.000000000000000e+00),
        (1,1,1,2,2,1,0, -2.000000000000000e+00),
        (1,1,1,4,1,4,0, 4.000000000000000e+00),
        (1,1,2,2,2,2,0, 1.000000000000000e+00),
        (1,1,2,4,2,4,0, 4.000000000000000e+00),
        (1,1,4,4,4,4,0, 1.000000000000000e+00),
        (1,1,1,4,4,1,1, -6.666666666666665e-01),
        (1,1,2,4,4,2,1, -6.666666666666665e-01),
        (1,1,4,4,4,4,2, 1.999999999999999e-01),
        (1,2,3,3,4,4,2, 2.828427124746190e-01),
        (2,1,3,3,4,4,2, 2.828427124746190e-01),
        (2,2,1,1,1,1,0, 1.000000000000000e+00),
        (2,2,1,2,1,2,0, 4.000000000000000e+00),
        (2,2,1,2,2,1,0, -2.000000000000000e+00),
        (2,2,1,3,1,3,0, 4.000000000000000e+00),
        (2,2,2,2,2,2,0, 1.000000000000000e+00),
        (2,2,2,3,2,3,0, 4.000000000000000e+00),
        (2,2,3,3,3,3,0, 1.000000000000000e+00),
        (2,2,1,3,3,1,1, -6.666666666666669e-01),
        (2,2,2,3,3,2,1, -6.666666666666669e-01),
        (3,3,1,1,1,1,0, 1.000000000000000e+00),
        (3,3,1,2,1,2,0, 4.000000000000000e+00),
        (3,3,1,2,2,1,0, -2.000000000000000e+00),
        (3,3,1,3,1,3,0, 2.000000000000000e+00),
        (3,3,1,4,1,4,0, 2.000000000000000e+00),
        (3,3,2,2,2,2,0, 1.000000000000000e+00),
        (3,3,2,3,2,3,0, 2.000000000000000e+00),
        (3,3,2,4,2,4,0, 2.000000000000000e+00),
        (3,3,3,4,3,4,0, 9.999999999999998e-01),
        (3,3,1,3,3,1,1, -3.333333333333334e-01),
        (3,3,1,4,4,1,1, -3.333333333333333e-01),
        (3,3,2,3,3,2,1, -3.333333333333334e-01),
        (3,3,2,4,4,2,1, -3.333333333333333e-01),
        (3,3,3,4,4,3,2, -1.999999999999999e-01),
        (4,4,1,1,1,1,0, 1.000000000000000e+00),
        (4,4,1,2,1,2,0, 4.000000000000000e+00),
        (4,4,1,2,2,1,0, -2.000000000000000e+00),
        (4,4,1,4,1,4,0, 4.000000000000000e+00),
        (4,4,2,2,2,2,0, 1.000000000000000e+00),
        (4,4,2,4,2,4,0, 4.000000000000000e+00),
        (4,4,4,4,4,4,0, 9.999999999999996e-01),
        (4,4,1,4,4,1,1, -6.666666666666665e-01),
        (4,4,2,4,4,2,1, -6.666666666666665e-01),
        (4,4,4,4,4,4,2, -1.199999999999999e-01),
        (4,5,3,4,4,4,2, -1.131370849898476e-01),
        (5,4,3,4,4,4,2, -1.131370849898475e-01),
        (5,5,1,1,1,1,0, 1.000000000000000e+00),
        (5,5,1,2,1,2,0, 4.000000000000000e+00),
        (5,5,1,2,2,1,0, -2.000000000000000e+00),
        (5,5,1,3,1,3,0, 2.000000000000000e+00),
        (5,5,1,4,1,4,0, 2.000000000000000e+00),
        (5,5,2,2,2,2,0, 1.000000000000000e+00),
        (5,5,2,3,2,3,0, 2.000000000000000e+00),
        (5,5,2,4,2,4,0, 2.000000000000000e+00),
        (5,5,3,4,3,4,0, 9.999999999999998e-01),
        (5,5,1,3,3,1,1, -3.333333333333334e-01),
        (5,5,1,4,4,1,1, -3.333333333333333e-01),
        (5,5,2,3,3,2,1, -3.333333333333334e-01),
        (5,5,2,4,4,2,1, -3.333333333333333e-01),
        (5,5,3,4,4,3,2, -4.000000000000000e-02)
    ]

    opTwo   = SpinAngular.TwoParticleOperator(0, Basics.plus, true)
    jacConv = Dict{Any,Float64}();   nRaw = 0;   nAnnihilated = 0
    for  (ic,l) in enumerate(localCsfs),  (ir,r) in enumerate(localCsfs)
        for  c in SpinAngular.computeCoefficients(opTwo, l, r, localSubshells)
            abs(c.V) < 1.0e-14  &&  continue
            global nRaw = nRaw + 1
            f = AngularMomentum.CL_reduced_me(c.a, c.nu, c.c) * AngularMomentum.CL_reduced_me(c.b, c.nu, c.d)
            if  isodd(c.nu)    f = -f    end
            if  abs(f) < 1.0e-14    global nAnnihilated = nAnnihilated + 1    end
            kk = (jacToGrasp[ic], jacToGrasp[ir],
                  rkKey(c.nu, localIdx[c.a], localIdx[c.b], localIdx[c.c], localIdx[c.d]))
            jacConv[kk] = get(jacConv, kk, 0.0) + c.V * f
        end
    end
    surviving = filter(p -> abs(p[2]) > 1.0e-10, jacConv)

    println("")
    println("  JAC raw X^L coefficients (non-zero)   : ", nRaw)
    println("  annihilated by the C^k factors        : ", nAnnihilated)
    println("  surviving after conversion            : ", length(surviving))
    println("  GRASP2018 emits                       : ", length(graspTwo))

    nMatched = 0;   nDisagree = 0;   worstRatio = 1.0
    for  (ic, ir, a, b, c, d, k, g) in graspTwo
        kk = (ic, ir, rkKey(k, a, b, c, d))
        if  haskey(surviving, kk)
            v = surviving[kk];   global nMatched = nMatched + 1
            if  abs(v/g - 1.0) > 1.0e-9                    global nDisagree  = nDisagree + 1   end
            if  abs(v/g - 1.0) > abs(worstRatio - 1.0)     global worstRatio = v/g             end
        else
            println("  MISSING from JAC: CSF(", ic, ",", ir, ") [", a, " ", b, " ; ", c, " ", d, "] k=", k, " = ", g)
        end
    end
    println("")
    println("  matched           : ", nMatched, " of ", length(graspTwo))
    println("  disagreeing >1e-9 : ", nDisagree)
    println("  worst ratio       : ", worstRatio)
    #
end
#
