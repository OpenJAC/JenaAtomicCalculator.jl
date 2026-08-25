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
    #   is not a result. Both implementations are run over the same CSF pairs and timed, warmed up first so that
    #   compilation is not measured, best of five.
    #
    #   REPORT (25-Aug-2026).  Set A is 15 CSFs of 1s^2 2s^2 2p^2 + 1s^2 2s 2p^3; the two occupation-changing rows use
    #   larger sets, 441 CSFs of 3d^4 4p / 3d^3 4p^2 and 241 of 2p^2 3d^2 / 2p 3d 4s 4p, because those cases barely occur
    #   in a small one.
    #
    #       case                        old (ms)   new (ms)   ratio     old (kB)   new (kB)   ratio
    #       1-particle rank 0              0.26       0.04     0.145        427         69     0.162
    #       1-particle rank 1              1.74       0.60     0.343       2296        494     0.215
    #       1-particle rank 2              1.35       0.41     0.305       1934        429     0.222
    #       2-particle equal occupation    2.09       1.20     0.575       3361       1256     0.374
    #       2-particle one electron moved  604        391      0.646     420638     248737     0.591
    #       2-particle two electrons moved 290        305      1.054     186769     183119     0.980
    #
    #   WHERE THE TIME ACTUALLY GOES, MEASURED RATHER THAN GUESSED.  Three rounds of plausible optimisation -- hoisting
    #   matrix elements out of loops they did not depend on, forming whole rank vectors at once, caching the Racah
    #   transform as a matrix -- moved the one-electron-move row only from 1.63 to 1.29. A sampling profile then showed
    #   the sweep was dominated by neither recoupling nor bookkeeping but by EXACT RATIONAL ARITHMETIC: BigInt allocation
    #   inside the Wigner-symbol package, reached through the Clebsch-Gordan in `shellReducedA`. Memoising that one
    #   function took the row from 1.29 to 0.65 in a single step, and the two-electron row from 1.46 to 1.05.
    #     THE LESSON IS THE ORDER OF OPERATIONS: the three hoists were right, worth keeping, and together worth less than
    #   a quarter of what the profile found in one measurement. They should have come second.
    #
    #   THE CACHES ARE PURE MEMOISATION and change no number: the arguments are a handful of small quantum numbers, the
    #   key is the complete argument list, and every verification in this file was re-run after adding them --
    #   45 952 one-electron-move coefficients, 386 678 two-electron-move coefficients and the whole one-particle
    #   inventory, all unchanged. `SpinAngularNew.clearCaches()` empties them and returns how many entries were held.
    #
    #   AGAINST GRASP2018, per coefficient on the same CSF list, this module is now about 4x slower rather than 9x. The
    #   remaining gap is structural rather than algorithmic: RKCO_GG writes into preallocated module-level buffers and
    #   allocates nothing per CSF pair, where both Julia modules build a fresh vector for every pair. An in-place API
    #   filling a caller-supplied buffer is the next gain and is NOT taken here.
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
elseif  true
    # Last visit:      23-Aug-2026
    # Last successful:  23-Aug-2026
    #
    # Branch k (A DEFECT IN GRASP2018, NOT IN JAC): the one place where JAC's two-particle coefficients and
    #   GRASP2018's disagree, and the argument that JAC is the one that is right.
    #
    #   WHERE IT SHOWS.  With at most TWO open subshells the two codes agree exactly -- 2p^2 61/61, 3d^2 82/82,
    #   4f^2 152/152, 3d^3 309/309, i.e. 604 coefficients.  With THREE singly-occupied open subshells (1s^2 2p^2 3d)
    #   six CSF pairs disagree, and not marginally: assembling the matrix element with ARBITRARY radial integrals --
    #   legitimate because it is LINEAR in them, so agreement decides equivalence of the coefficient sets whatever the
    #   split -- gives JAC 0.056399383 against GRASP -4.874986331 on the worst pair.
    #
    #   THE PAIR.  GRASP CSFs 19 and 20 of that configuration are
    #
    #       CSF 19:  1s ( 2)  2p-( 1)  2p ( 1)  3d-( 1)      1/2  3/2  3/2      coupling 1  ->  5/2+
    #       CSF 20:  1s ( 2)  2p-( 1)  2p ( 1)  3d-( 1)      1/2  3/2  3/2      coupling 2  ->  5/2+
    #
    #   Same subshells, same J in every subshell, same total J and parity; they differ ONLY in the intermediate
    #   coupling of 2p- with 2p.  Distinct, orthogonal CSFs of one symmetry.
    #
    #   THREE INDEPENDENT ARGUMENTS, none of which requires trusting either code.
    #
    #   (1) AN EXACT IDENTITY.  A k = 0 DIRECT term (a,b,a,b) carries C^0, which is the identity, so the operator is
    #       n_a n_b.  Number operators are diagonal in the occupations, CSFs are orthonormal eigenstates of the
    #       occupations, so  <i| n_a n_b |j>  =  N_a N_b delta_ij  and MUST vanish for i /= j, whatever the coupling
    #       differs by.  GRASP emits, for the pair (20,19),
    #
    #           (2p- , 3d- ; 2p- , 3d-)  k=0   0.99999999999999978        = N(2p-) * N(3d-) = 1
    #           (2p  , 3d- ; 2p  , 3d-)  k=0   0.99999999999999956        = N(2p)  * N(3d-) = 1
    #           (1s  , 3d- ; 1s  , 3d-)  k=0   2.00000000000000000        = N(1s)  * N(3d-) = 2
    #
    #       -- the DIAGONAL occupation products, on an OFF-DIAGONAL pair.  JAC emits none of them, and produces a
    #       strict subset of GRASP's terms with no extras anywhere.
    #
    #   (2) GRASP CONTRADICTS ITSELF, with JAC nowhere in the argument.  For that same pair its own ONESCALAR returns
    #       NOTHING, i.e. <20| n_a |19> = 0, exactly as orthogonality requires -- while its RKCO_GG returns n_a n_b.
    #       The one-particle and two-particle outputs of a single code disagree.
    #
    #   (3) THE AUTHOR DOCUMENTED THIS BUG CLASS IN 2009.  src/lib/librang90/ReadMe, in Lithuanian, opens
    #       "Orginalioje programoje Grasp2K yra klaidu !!!!" -- there are errors in the original Grasp2K -- prescribes
    #       inserting IF(JA.NE.JB) RETURN in rkco_gg.f, and closes: "The program computes wrongly when the calculation
    #       contained two or more configurations whose difference reduced only to the value of the SENIORITY quantum
    #       number."  Dated 2009.11.22, G.G.  CSFs 19 and 20 differ only in a coupling quantum number: exactly that.
    #
    #   WHY THE FIX DOES NOT COVER IT.  In rkco_gg.f90 the guard sits at about line 356, AFTER "IF (INCOR .LT. 1)
    #   RETURN" and "IF (NCORE .EQ. 0) RETURN", so it protects only the CORE section.  The offending terms survive with
    #   INCOR = 0, so they come from the MAIN path above that line, which the guard never reaches.  The onescalar half
    #   of the same 2009 fix IS effective, which is precisely why argument (2) works.
    #
    #   WHY THIS SURVIVED DECADES OF USE.  The defect is narrow: it bites only CSF pairs differing solely in a
    #   coupling or seniority label.  That is why it does not corrupt ordinary calculations, why it took a systematic
    #   comparison to surface, and why the author had to leave a note about it at all.
    #
    #   WHAT THIS BRANCH CHECKS, since the GRASP numbers cannot be recomputed without the harness in tools/: that JAC
    #   obeys the identity.  Note what the identity does and does not forbid -- and the first version of this branch
    #   got it wrong, which is why the check is stated carefully now.  A two-body operator MAY connect two CSFs that
    #   differ in their coupling, through exchange and through higher multipoles; JAC returns ten such coefficients
    #   for this pair, and that is correct.  What must vanish is only the k = 0 DIRECT term, because that one alone
    #   reduces to n_a n_b.  So the branch counts k = 0 direct terms specifically.
    #
    #   REPORT (23-Aug-2026): JAC returns 0 rank-0 one-particle coefficients and 0 k=0 DIRECT two-particle
    #   coefficients for the pair, in both orderings, out of ten two-particle coefficients in total.  GRASP returns
    #   three k=0 direct terms, equal to the occupation products.  The harness that produced the GRASP side is
    #   tools/diag-grasp-angular.jl.
    #
    localRel = ConfigurationR[]
    for  conf in [Configuration("1s^2 2p^2 3d")]
        append!(localRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf))
    end
    localSubshells = Basics.generateSubshellList(localRel)
    Defaults.setDefaults("relativistic subshell list", localSubshells; printout=false)
    localCsfs = CsfR[]
    for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSubshells))    end

    i2pm = findfirst(sh -> string(sh) == "2p_1/2", localSubshells)
    i2p  = findfirst(sh -> string(sh) == "2p_3/2", localSubshells)
    i3dm = findfirst(sh -> string(sh) == "3d_3/2", localSubshells)

    # ... the two CSFs that differ ONLY in the intermediate coupling of 2p- with 2p
    cand = [i for (i,c) in enumerate(localCsfs)
                if c.occupation[i2pm] == 1 && c.occupation[i2p] == 1 && c.occupation[i3dm] == 1 &&
                   Basics.twice(c.J) == 5]
    println("\n  CSFs of 1s^2 2p_1/2 2p_3/2 3d_3/2 with J = 5/2:")
    for  i in cand
        c = localCsfs[i]
        println("    JAC CSF $i:  subshellX = ", Basics.twice.(c.subshellX) .// 2, "   J = $(c.J)$(string(c.parity))")
    end

    if  length(cand) >= 2
        i, j = cand[1], cand[2]
        n1a  = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(0, Basics.plus),
                                                  localCsfs[i], localCsfs[j], localSubshells)
        n1b  = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(0, Basics.plus),
                                                  localCsfs[j], localCsfs[i], localSubshells)
        o2a  = SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus, true),
                                               localCsfs[i], localCsfs[j], localSubshells)
        o2b  = SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus, true),
                                               localCsfs[j], localCsfs[i], localSubshells)
        nz(v)     = count(c -> abs(c.V) > 1.0e-14, v)
        direct0(v)= count(c -> c.nu == 0 && c.a == c.c && c.b == c.d && abs(c.V) > 1.0e-14, v)
        println("\n  JAC, for the pair (", i, ",", j, ") and its transpose:")
        println("    rank-0 one-particle coefficients        : ", length(n1a), " and ", length(n1b),
                "      (must be 0 -- orthogonality)")
        println("    two-particle coefficients, all told     : ", nz(o2a), " and ", nz(o2b),
                "    (may be non-zero: exchange and higher multipoles)")
        println("    of those, k=0 DIRECT terms (a,b;a,b)    : ", direct0(o2a), " and ", direct0(o2b),
                "      (must be 0 -- these reduce to n_a n_b)")
        println("\n  GRASP2018, same pair, returns THREE k=0 direct terms: 1.0, 1.0, 2.0 -- the occupation products.")
    end
    #
elseif  true
    # Last visit:      24-Aug-2026
    # Last successful:  24-Aug-2026
    #
    # Branch l (THE TWO-PARTICLE DIAGONAL CASE, COMPLETE): the electron-electron coefficients from SpinAngularNew
    #   against SpinAngular's, compared as complete LISTS rather than term by term -- so that a coefficient present on
    #   one side and absent on the other counts as a failure rather than going unnoticed.
    #
    #   THREE TERMS, EACH SETTLED SEPARATELY AND EACH THE SUBJECT OF A WRONG TURN WORTH RECORDING.
    #
    #   (1) DIRECT, between two distinct subshells. The scalar product of two rank-k tensors, one on each subshell, so
    #       it reuses shellReducedW and substitutionRecoupling with both ranks equal to k coupled to zero. Its reach was
    #       MEASURED rather than assumed -- exact in all four occupation classes, closed/closed, closed/open,
    #       closed/single, single/single.
    #
    #       THE WRONG TURN: a closed form fitted to CLOSED-SHELL data, which gave this term right at k = 0 and wrong
    #       everywhere else. A closed shell forces J = 0 and hides the J-dependence entirely, so the set it was fitted
    #       on could not discriminate. It was withdrawn rather than shipped.
    #
    #   (2) EXCHANGE, between two distinct subshells. Not the direct term relabelled: it expands over the direct
    #       channel at intermediate ranks by the Racah transformation,
    #
    #           V^k(a,b;b,a) = sum_K (2K+1) { j_a j_b k ; j_b j_a K } V^K(a,b;a,b)
    #
    #       arrived at BY ELIMINATION. Of five candidate 6j orderings three VANISH at ranks where the coefficient does
    #       not, which refutes them outright; of the two survivors only the phase (-1)^(j_a+j_b+k) holds. The general
    #       sum then had to reproduce that collapsed case before replacing it.
    #
    #   (3) SAME-SUBSHELL. The two-body quasispin object, LESS a normal-ordering correction:
    #
    #           V^k(a,a;a,a) = 0.5 * ( WW(k)/sqrt(2k+1) - (-1)^(2j+k) W(0)/sqrt(2j+1) ) / sqrt(2J+1)
    #
    #       THE WRONG TURN HERE IS THE INSTRUCTIVE ONE. Two attempts computed WW and compared it DIRECTLY against the
    #       coefficient. It disagrees everywhere -- for j = 3/2, N = 2, J = 0 the coefficient is 0.25 at every rank
    #       while WW gives 1.0, 0, 2.236, 0 -- and I recorded that as refuting the closure ROUTE. That was wrong twice
    #       over: the route is what SpinAngular itself uses, and the object was already correct. What was missing was
    #       the factor of one half and the subtraction. A route wrongly marked dead does not get revisited, which is
    #       why the withdrawal was worth its own commit.
    #
    #       Physically the subtraction is not decoration: both operators act on one shell, so the two-body object
    #       counts each electron with itself, and that self-interaction is removed by the rank-0 one-body element --
    #       not by anything rank-dependent.
    #
    #   REPORT (24-Aug-2026), complete lists over nine configurations:
    #
    #       configuration        CSFs   matched   missing   extra
    #       1s^2 2s^2 2p^2          5        95         0       0
    #       1s^2 3d^2               9       115         0       0
    #       1s^2 3d^3              19       345         0       0
    #       1s^2 2s^2 2p^4          5       127         0       0
    #       1s^2 2s 2p              4        48         0       0
    #       1s^2 2p^2 3s            8       151         0       0
    #       1s^2 3d^2 4s           16       350         0       0
    #       1s^2 2s^2 2p 3d        12       308         0       0
    #       1s^2 4f^2              13       205         0       0
    #       TOTAL                          1744         0       0     worst ratio 1.000000000000
    #
    #   OFF-DIAGONAL CSF PAIRS STILL RAISE, and that is the remaining gap. Such a pair needs recoupling this code does
    #   not do, and returning the diagonal answer for it would be a wrong number rather than a missing one.
    #
    localCfgs = ["1s^2 2s^2 2p^2", "1s^2 3d^2", "1s^2 3d^3", "1s^2 2s^2 2p^4", "1s^2 2s 2p",
                 "1s^2 2p^2 3s", "1s^2 3d^2 4s", "1s^2 2s^2 2p 3d", "1s^2 4f^2"]
    tMatched = 0;   tMissing = 0;   tExtra = 0;   tWorst = 1.0
    println("\n  configuration        CSFs  matched  missing  extra   worst ratio")
    for  cfg in localCfgs
        localRel = Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(cfg))
        localSub = Basics.generateSubshellList(localRel)
        Defaults.setDefaults("relativistic subshell list", localSub; printout=false)
        local localCsfs = CsfR[]
        for  relconf in localRel    append!(localCsfs, Basics.generateCsfRs(relconf, localSub))    end
        local nm = 0;  local nn = 0;  local nx = 0;  local wr = 1.0
        for  c in localCsfs
            old = [x for x in SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0, Basics.plus, true),
                                                              c, c, localSub)   if abs(x.V) > 1.0e-14]
            new = SpinAngularNew.computeCoefficients(SpinAngularNew.TwoParticleOperator(), c, c, localSub)
            kk(x) = (x.nu, string(x.a), string(x.b), string(x.c), string(x.d))
            dOld = Dict{Any,Float64}();   for x in old   dOld[kk(x)] = get(dOld,kk(x),0.0) + x.V   end
            dNew = Dict{Any,Float64}();   for x in new   dNew[kk(x)] = get(dNew,kk(x),0.0) + x.V   end
            for  (k,v) in dOld
                if  haskey(dNew,k)
                    nm += 1;   r = v/dNew[k];   abs(r-1) > abs(wr-1) && (wr = r)
                else    nn += 1
                end
            end
            for  k in keys(dNew)    haskey(dOld,k) || (nx += 1)    end
        end
        @printf("  %-18s %4d %8d %8d %6d   %.12f\n", cfg, length(localCsfs), nm, nn, nx, wr)
        global tMatched += nm;  global tMissing += nn;  global tExtra += nx
        abs(wr-1) > abs(tWorst-1) && (global tWorst = wr)
    end
    println("")
    @printf("  TOTAL matched %d   missing %d   extra %d   worst ratio %.12f\n", tMatched, tMissing, tExtra, tWorst)
    #
elseif  true
    # Last visit:      24-Aug-2026
    # Last successful:  24-Aug-2026
    #
    # Branch m (OFF-DIAGONAL IN THE COUPLING): CSF pairs of equal occupation but different coupling, which branch l
    #   refused. They turn out to need almost nothing new -- and finding that out exposed a defect of exactly the kind
    #   this file criticises GRASP2018 for in branch k.
    #
    #   WHY ALMOST NOTHING WAS NEEDED.  The direct and exchange terms already go through substitutionRecoupling, which
    #   takes bra and ket separately and returns zero when their couplings cannot be connected. So 56 of 56 coefficients
    #   on coupling-off-diagonal pairs were already exactly right the first time they were asked for.
    #
    #   THE DEFECT.  Twenty coefficients were also produced that SpinAngular does not produce, and every one of them was
    #   the same-subshell term of the CLOSED 1s^2 shell, valued 0.5 -- its DIAGONAL value -- on a pair that differs only
    #   in the coupling of some OTHER subshell.  That is the same shape as the eighteen spurious terms GRASP emits in
    #   branch k: a diagonal formula applied to an off-diagonal pair.  Having spent the morning arguing that GRASP was
    #   wrong to do it, I had it in my own same-subshell routine within the hour.
    #
    #   WHY IT HAPPENED, AND THE FIX.  twoParticleSameShell read only its own subshell's quantum numbers. For a closed
    #   1s^2 the bra and ket shell terms are identical whatever the rest of the CSF does, so it returned the diagonal
    #   value.  But (W^(k) x W^(k))^(0) is a SCALAR in the subshell it acts in: it cannot change any coupling, so every
    #   other subshell must be in an identical state and the running couplings X must agree throughout, or the two CSFs
    #   are orthogonal and the element vanishes.  The routine now checks that.  The direct and exchange terms never had
    #   the problem because their recoupling factor enforces it for them -- the same-subshell term was the one place
    #   where nothing did.
    #
    #   THE LESSON IS ABOUT WHERE TO LOOK, not about the sign of a 6j.  A routine that consults only the object it acts
    #   on cannot know whether the rest of the state permits it to act at all, and it will happily return the diagonal
    #   answer for an off-diagonal pair. Both codes made that mistake in the same week, in the same term.
    #
    #   REPORT (24-Aug-2026), ALL equal-occupation pairs, diagonal and coupling-off-diagonal, as complete lists:
    #
    #       configuration        pairs   matched   missing   extra
    #       1s^2 2s^2 2p^2           9        95         0       0
    #       1s^2 3d^2               29       115         0       0
    #       1s^2 3d^3              127       389         0       0
    #       1s^2 2s^2 2p^4           9       127         0       0
    #       1s^2 2s 2p               8        48         0       0
    #       1s^2 2p^2 3s            26       163         0       0
    #       1s^2 3d^2 4s            98       386         0       0
    #       1s^2 2s^2 2p 3d         40       308         0       0
    #       1s^2 4f^2               61       205         0       0
    #       TOTAL                  407      1836         0       0    worst ratio 1.000000000000
    #
    #   WHAT STILL RAISES: pairs that MOVE electrons between subshells. Those need creation and annihilation on two
    #   shells at once -- SpinAngular's twoParticleDiffOcc2, twoParticle6, twoParticle15to18 and twoParticle19to42 --
    #   and none of it is written here. Returning the equal-occupation answer for such a pair would be a wrong number.
    #
    localCfgs = ["1s^2 2s^2 2p^2","1s^2 3d^2","1s^2 3d^3","1s^2 2s^2 2p^4","1s^2 2s 2p",
                 "1s^2 2p^2 3s","1s^2 3d^2 4s","1s^2 2s^2 2p 3d","1s^2 4f^2"]
    tP = 0;  tM = 0;  tN = 0;  tX = 0;  tW = 1.0
    println("\n  configuration        pairs   matched  missing  extra   worst ratio")
    for  cfg in localCfgs
        local lRel = Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(cfg))
        local lSub = Basics.generateSubshellList(lRel)
        Defaults.setDefaults("relativistic subshell list", lSub; printout=false)
        local lCsfs = CsfR[]
        for  rc in lRel    append!(lCsfs, Basics.generateCsfRs(rc, lSub))    end
        local np = 0;  local nm = 0;  local nn = 0;  local nx = 0;  local wr = 1.0
        for  (i,l) in enumerate(lCsfs),  (j,r) in enumerate(lCsfs)
            l.occupation == r.occupation  ||  continue
            np += 1
            old = [x for x in SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0,Basics.plus,true),
                                                              l, r, lSub)   if abs(x.V) > 1.0e-14]
            new = SpinAngularNew.computeCoefficients(SpinAngularNew.TwoParticleOperator(), l, r, lSub)
            kk(x) = (x.nu, string(x.a), string(x.b), string(x.c), string(x.d))
            dO = Dict{Any,Float64}();   for x in old   dO[kk(x)] = get(dO,kk(x),0.0) + x.V   end
            dN = Dict{Any,Float64}();   for x in new   dN[kk(x)] = get(dN,kk(x),0.0) + x.V   end
            for  (k,v) in dO
                if  haskey(dN,k);   nm += 1;  rr = v/dN[k];  abs(rr-1) > abs(wr-1) && (wr = rr)
                else   nn += 1
                end
            end
            for  k in keys(dN)    haskey(dO,k) || (nx += 1)    end
        end
        @printf("  %-18s %5d %9d %8d %6d   %.12f\n", cfg, np, nm, nn, nx, wr)
        global tP += np;  global tM += nm;  global tN += nn;  global tX += nx
        abs(wr-1) > abs(tW-1) && (global tW = wr)
    end
    println("")
    @printf("  TOTAL %d pairs: matched %d, missing %d, extra %d, worst ratio %.12f\n", tP, tM, tN, tX, tW)
    #
elseif  true
    # Last visit:      24-Aug-2026
    # Last successful:  24-Aug-2026
    #
    # Branch n (MOVING ONE ELECTRON BETWEEN SUBSHELLS): the first two-particle case in which the two CSFs do NOT have the
    #   same occupations. Branch m closed the equal-occupation problem; this one opens the occupation-changing one, where
    #   the two-body operator creates an electron in one subshell, annihilates one in another, and acts a second time on a
    #   spectator subshell. SpinAngular spends four routines and a forty-branch index tree on this case
    #   (twoParticleDiffOcc2 -> twoParticle7to14 -> twoParticle7to8 / 9to10 / 11to14).
    #
    #   ONE RECOUPLING ROUTINE REPLACES THE CASE TREE, and that is the substantial part. An operator acting on any number
    #   of subshells, its shell tensors coupled along the subshell chain, needs only three kinds of factor: a peel between
    #   acting subshells, a nine-j junction at each acting subshell above the lowest, and the acting factor at the lowest
    #   one. `treeRecoupling` is those three in a loop, and it is not an alternative to the two routines already validated
    #   here -- it CONTAINS them. Checked rather than claimed: on 33 CSFs over five configurations it reproduces
    #   `chainRecoupling` on 1476 comparisons and `substitutionRecoupling` on 164, worst difference 2.2e-16, one ulp.
    #
    #   THREE PHASES DECIDE THE SIGN, and each was established on data rather than assumed:
    #     (1) the Jordan-Wigner string over the occupations between the two changing subshells, taken UNCONDITIONALLY;
    #     (2) (-1)^(j_a + j_d - k + 1) when the creation sits on the higher subshell index -- the same phase the
    #         one-particle substitution already carried;
    #     (3) (-1)^(j_a + j_d + k + 1) when the spectator lies BETWEEN the two changing subshells, because the coupling
    #         tree then joins the acceptor with the spectator instead of with the donor.
    #
    #   HOW (3) WAS FOUND, AND WHY THE FIRST TWO ATTEMPTS AT IT WERE WRONG.  With the spectator OUTSIDE, the assembly was
    #   exact on 1136 coefficients the first time it was run. With it BETWEEN, every magnitude was still exact -- 518 of
    #   518 -- and only the sign was wrong, so a phase was all that was missing. Two guesses failed and are recorded
    #   because each looked convincing: excluding the spectator from the Jordan-Wigner string (partially right, hence
    #   misleading), and the sign (-1)^(sum of the nine arguments) of the middle junction's nine-j, which matched 234 of
    #   518, i.e. nothing. A parity search then returned an exact rule -- but only after the search BASIS was corrected:
    #   the first search ranged over 2j_a, 2j_d, 2j_s and found nothing, because a re-pairing phase such as
    #   (-1)^(j_a+j_d-k) is not a parity of the doubled values. The rule it found, (-1)^(k + j_a + j_d + N_s + 1), then
    #   simplified: N_s only undid the wrong Jordan-Wigner exclusion, leaving a phase in the three RANKS alone -- which is
    #   what a re-pairing of three tensors coupled to zero must be, so the rule is not merely a fit.
    #   IT WAS THEN TESTED OUT OF SAMPLE: 8820 further coefficients from six configuration sets not used to find it, all
    #   exact.
    #
    #   THE PARTNER TERM RUNS IN THE OPPOSITE DIRECTION.  Each spectator contributes two coefficients per rank, one for
    #   each of the two pairings, and the second follows from the first by the same Racah sum the equal-occupation
    #   exchange term uses. But the sum has a direction: with the spectator outside, the assembly yields the DIRECT
    #   pairing and the transform produces the crossed one; with it between, the assembly yields the CROSSED pairing and
    #   the transform must run back. Reading the six-j's bottom row on the primary quadruple's own c and d gets this
    #   right; using one direction for both leaves 214 of 872 coefficients wrong.
    #
    #   REPORT (24-Aug-2026).  Term by term against SpinAngular over ten configuration sets, split by arrangement:
    #
    #       spectator outside     61 470 coefficients    0 differing   0 missing   0 extra
    #       spectator between     20 718 coefficients    0 differing   0 missing   0 extra
    #       TOTAL                 82 188 coefficients    0 differing   0 missing   0 extra
    #
    #   and, below, as COMPLETE lists -- every coefficient SpinAngular emits for the pair, no topology set aside and
    #   NOTHING REFUSED: every one-electron-move pair is now handled. Measured over eleven configuration sets,
    #   8574 pairs and 46 020 coefficients agree, nothing missing, nothing extra, worst ratio 1.000000000000.
    #
    #   THE SAME-SUBSHELL SPECTATOR, ADDED 24-Aug-2026, COMPLETES THIS CASE.  When the spectator coincides with the
    #   acceptor or the donor, three of the four one-electron operators fall on one subshell and the shell matrix element
    #   becomes a coupled tensor, built by CLOSURE over intermediate subshell terms from pieces already here. Both
    #   orderings were needed and each reproduces the predecessor's own tensor exactly, 8897 of 8897 values apiece over
    #   every j up to 9/2.
    #
    #   WHICH ORDERING GOES WHERE IS PHYSICS, NOT A PHASE, and getting it wrong is invisible almost everywhere. W belongs
    #   on the side holding FEWER electrons: (a x W) when the spectator is the acceptor, (W x a) when it is the donor.
    #   Taking (a x W) for both agrees with the predecessor on 2621 of 2621 coefficients of 3d^4 4p and on every case in
    #   five other configurations -- and is WRONG, because a and W on one subshell do not commute, so the two orderings
    #   are different operators rather than one with a sign. The error surfaces only when the donor subshell is CLOSED in
    #   the ket: W^(k>0) on a closed shell vanishes, so a rank-1 coefficient that physically exists is silently dropped.
    #   ONE configuration of the thirteen below, 1s^2 2s^2 3s against 1s 2s^2 3s^2, showed it -- two coefficients wrong by
    #   a factor of two and two missing. A calibration over thousands of values on the wrong side of a structural
    #   distinction is not evidence, and this is the clearest example of it in this file.
    #
    localCfgs = [(["1s^2 2s 3d","1s^2 2p 3d"],        "1s^2 2s 3d / 2p 3d"),
                 (["2p^2 3d","2p 3d^2"],              "2p^2 3d / 2p 3d^2"),
                 (["4d^2 5p","4d 5p^2"],              "4d^2 5p / 4d 5p^2"),
                 (["3p^3 3d","3p^2 3d^2"],            "3p^3 3d / 3p^2 3d^2"),
                 (["1s^2 2s 4f","1s^2 2p 4f"],        "1s^2 2s 4f / 2p 4f"),
                 (["1s^2 2s 3d 4s","1s^2 2p 3d 4s"],  "two spectators at once"),
                 (["3d^3 4s","3d^2 4s^2"],            "3d^3 4s / 3d^2 4s^2"),
                 (["3d^4 4p","3d^3 4p^2"],            "3d^4 4p / 3d^3 4p^2"),
                 (["1s^2 2s^2 3s","1s 2s^2 3s^2"],    "closed 2s BETWEEN"),
                 (["1s^2 2s^2 2p^4","1s^2 2s 2p^5"],  "closed 2p^4"),
                 (["4f^2 5s","4f 5s^2"],              "4f^2 5s / 4f 5s^2")]
    hT = 0;  rT = 0;  mT = 0;  nT = 0;  xT = 0;  wT = 1.0
    println("\n  configuration            handled  raised   matched  missing  extra   worst ratio")
    for  (cfgs, tag) in localCfgs
        local lRel = ConfigurationR[]
        for  c in cfgs    append!(lRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(),
                                                                     Configuration(c)))    end
        local lSub = Basics.generateSubshellList(lRel)
        Defaults.setDefaults("relativistic subshell list", lSub; printout=false)
        local lCsfs = CsfR[]
        for  rc in lRel    append!(lCsfs, Basics.generateCsfRs(rc, lSub))    end
        local nh = 0;  local nr = 0;  local nm = 0;  local nn = 0;  local nx = 0;  local wr = 1.0
        for  l in lCsfs,  r in lCsfs
            l.J == r.J  &&  l.parity == r.parity  ||  continue
            local dd = l.occupation - r.occupation
            count(!=(0), dd) == 2  &&  sum(abs, dd) == 2   ||  continue
            local mv = findall(!=(0), dd)
            local iC = dd[mv[1]] > 0 ? mv[1] : mv[2];     local iA = dd[mv[1]] > 0 ? mv[2] : mv[1]
            nh += 1
            old = [x for x in SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0,Basics.plus,true),
                                                             l, r, lSub)   if abs(x.V) > 1.0e-14]
            new = SpinAngularNew.computeCoefficients(SpinAngularNew.TwoParticleOperator(), l, r, lSub)
            kk(x) = (x.nu, string(x.a), string(x.b), string(x.c), string(x.d))
            dO = Dict{Any,Float64}();   for x in old   dO[kk(x)] = get(dO,kk(x),0.0) + x.V   end
            dN = Dict{Any,Float64}();   for x in new   dN[kk(x)] = get(dN,kk(x),0.0) + x.V   end
            for  (k,v) in dO
                if  haskey(dN,k);   nm += 1;  rr = v/dN[k];  abs(rr-1) > abs(wr-1) && (wr = rr)
                else   nn += 1
                end
            end
            for  k in keys(dN)    haskey(dO,k) || (nx += 1)    end
        end
        @printf("  %-22s %7d %7d %9d %8d %6d   %.12f\n", tag, nh, nr, nm, nn, nx, wr)
        global hT += nh;  global rT += nr;  global mT += nm;  global nT += nn;  global xT += nx
        abs(wr-1) > abs(wT-1) && (global wT = wr)
    end
    println("")
    @printf("  TOTAL: %d pairs handled, %d raised | matched %d, missing %d, extra %d, worst ratio %.12f\n",
            hT, rT, mT, nT, xT, wT)
    #

elseif  true
    # Last visit:      24-Aug-2026
    # Last successful:  24-Aug-2026
    #
    # Branch o (THE COMPLETE DIFFERENCE INVENTORY): the question this branch answers is the one that decides whether
    #   SpinAngularNew may ever replace SpinAngular -- do the two return the SAME set of coefficients, and where they do
    #   not, exactly why? Every other branch compares selected cases and reports agreement. This one takes EVERY key
    #   either module emits, over twelve configurations and ranks 0 to 3, and forces each into a named class. A key that
    #   fell into no class would be an unexplained difference, and the count of those is printed.
    #
    #   THERE ARE EXACTLY TWO CLASSES OF DIFFERENCE, and both are deliberate:
    #
    #   (1) THE RANK-0 NORMALIZATION, ratio exactly sqrt(2 j_a + 1). This is the convention change branch c documents:
    #       SpinAngularNew carries the factor at EVERY rank, SpinAngular has it commented out at rank 0 and its callers
    #       compensate (module-Hamiltonian.jl:277 re-applies it, Hfs.amplitude divides it out). A migration must delete
    #       those compensating lines, and that is the one edit that carries real risk.
    #
    #   (2) THE PARITY SELECTION RULE. SpinAngularNew's `isAllowed1p` refuses a pair whose l_a + l_b does not match the
    #       operator's parity; SpinAngular does not, and returns everything the triangle condition allows. That is not a
    #       guess about its intent: given the SAME CSFs and ranks, SpinAngular returns BIT-IDENTICAL lists for
    #       Basics.plus and Basics.minus, so it ignores the parity it is handed and leaves it to the caller's radial
    #       integral, which vanishes anyway. Every one of the 5648 keys in this class was checked individually, and all
    #       5648 have odd l_a + l_b against an even-parity operator -- none unexplained.
    #
    #       THIS ONE IS BENIGN IN THE HAMILTONIAN AND NOT BENIGN IN GENERAL. The products are zero either way, so no
    #       matrix element changes. But a caller that passes the WRONG parity gets a silently shorter list from the new
    #       module and the full list from the old, and a caller that counts terms sees different numbers. Anyone
    #       migrating a module must check which parity it passes.
    #
    #   AND ONE REAL BUG, FOUND BY THIS BRANCH AND FIXED (24-Aug-2026), which is why the branch exists rather than being
    #   a formality. Before the fix the inventory showed 949 keys EXTRA in the new module at ranks 1 to 3. They were
    #   spurious: for two CSFs of 1s^2 3d^2 4s with identical occupations, identical subshellJ and identical seniority,
    #   differing ONLY in the intermediate coupling X_3 (1 against 2), the new module returned a 4s_1/2 coefficient where
    #   the old correctly returned none. A rank-1 operator acting on subshell 4 cannot change X_3.
    #     THE GIVE-AWAY WAS A SYMMETRY, NOT A REFERENCE VALUE: the pairs (7,8) and (8,7) have equal J, so a Hermitian
    #   operator must give equal magnitudes, and it gave +0.745 against -0.447.
    #     THE CAUSE: `nonScalarGeneral` guarded the OTHER subshells' subshellJ and seniority but not their intermediate
    #   couplings, and `chainRecoupling` takes X_{ip-1} from the BRA alone -- so for a pair differing below the acting
    #   subshell it returned the diagonal answer instead of zero. The guard now sits in `actingFactor`, where the
    #   operator is the identity below the lowest acting subshell, and `chainRecoupling` was made to delegate to
    #   `peelRange` and `actingFactor` so that the guard cannot be half-applied.
    #     THIS IS THE FOURTH TIME THIS EXACT DEFECT HAS APPEARED HERE: in GRASP2018 (branch k), in this module's
    #   same-subshell two-particle term (branch m), and now in its one-particle term. Each time the shape is identical --
    #   a routine that consults only the subshell it acts on, returning the diagonal answer for an off-diagonal pair.
    #   A new method in this module should be assumed to have it until a coupling-off-diagonal pair says otherwise.
    #
    #   REPORT (24-Aug-2026), twelve configurations, 12055 CSF pairs at each rank:
    #
    #       operator          identical   sqrt(2j+1)   parity-only-in-old   num.zero-new   EXTRA   DISAGREE   unexplained
    #       1-particle k=0            0          684                  156              0       0          0             0
    #       1-particle k=1         3112            0                 1324              3       0          0             0
    #       1-particle k=2         4073            0                 2308              1       0          0             0
    #       1-particle k=3         3626            0                 1860              1       0          0             0
    #       2-particle (scalar)    6603            0                    0              0       0          0             0
    #
    #   (3) A THIRD CLASS, five keys in all and harmless, but listed rather than folded into another column. The new
    #       module decides zeros by SELECTION RULE and not by magnitude, so it has no `abs(wa) >= 2.0e-10` cutoff and
    #       will occasionally emit a coefficient of order 1e-17 that the old module's cutoff removes. Nothing downstream
    #       can notice a 1e-17, but a comparison that counted it as a disagreement would be wrong and one that hid it in
    #       another column would be worse.
    #
    #   Neither module ever emitted a duplicate key, so no summing is needed before comparing -- which is worth knowing,
    #   since example-Ak.jl compared totals precisely because it assumed otherwise.
    #
    localSets = [["1s^2 2s"], ["1s^2 2p"], ["1s^2 2s^2 2p^2"], ["1s^2 3d^2"], ["1s^2 3d^3"], ["1s^2 4f^2"],
                 ["1s^2 2s^2 2p 3d"], ["1s^2 3d^2 4s"], ["1s^2 2s 3d","1s^2 2p 3d"], ["2p^2 3d","2p 3d^2"],
                 ["1s^2 2s^2 2p","1s^2 2s 2p^2"], ["4d^2 5p","4d 5p^2"]]
    println("\n  operator            identical  sqrt(2j+1)  parity-old-only  num.zero-new   EXTRA  DISAGREE  unexplained  dup.keys")
    for  rank in [0, 1, 2, 3, -1]
        nId = 0; nCv = 0; nPa = 0; nEx = 0; nDs = 0; nUn = 0; nDp = 0; nZr = 0
        for  cfgSet in localSets
            local lRel = ConfigurationR[]
            for  c in cfgSet    append!(lRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(),
                                                                           Configuration(c)))    end
            local lSub = Basics.generateSubshellList(lRel)
            Defaults.setDefaults("relativistic subshell list", lSub; printout=false)
            local lCsfs = CsfR[]
            for  rc in lRel    append!(lCsfs, Basics.generateCsfRs(rc, lSub))    end
            for  l in lCsfs,  r in lCsfs
                local old = [];   local new = []
                if  rank >= 0
                    old = [x for x in SpinAngular.computeCoefficients(
                              SpinAngular.OneParticleOperator(rank, Basics.plus, true), l, r, lSub) if abs(x.T) > 1.0e-14]
                    new = SpinAngularNew.computeCoefficients(SpinAngularNew.OneParticleOperator(rank, Basics.plus),
                                                             l, r, lSub)
                else
                    local dd = l.occupation - r.occupation
                    local reach = false
                    if      all(==(0), dd)                          reach = (l.J == r.J && l.parity == r.parity)
                    elseif  count(!=(0), dd) == 2 && sum(abs, dd) == 2
                        local mv = findall(!=(0), dd)
                        local iC = dd[mv[1]] > 0 ? mv[1] : mv[2];   local iAn = dd[mv[1]] > 0 ? mv[2] : mv[1]
                        reach = !(r.occupation[iC] >= 1 || r.occupation[iAn] >= 2) && l.J == r.J && l.parity == r.parity
                    end
                    reach  ||  continue
                    old = [x for x in SpinAngular.computeCoefficients(
                              SpinAngular.TwoParticleOperator(0, Basics.plus, true), l, r, lSub) if abs(x.V) > 1.0e-14]
                    new = SpinAngularNew.computeCoefficients(SpinAngularNew.TwoParticleOperator(), l, r, lSub)
                end
                kk(x)  = rank >= 0 ? (x.nu, string(x.a), string(x.b)) :
                                     (x.nu, string(x.a), string(x.b), string(x.c), string(x.d))
                vv(x)  = rank >= 0 ? x.T : x.V
                dO = Dict{Any,Float64}();  for x in old   haskey(dO,kk(x)) && (nDp += 1);  dO[kk(x)] = vv(x)   end
                dN = Dict{Any,Float64}();  for x in new   haskey(dN,kk(x)) && (nDp += 1);  dN[kk(x)] = vv(x)   end
                for  (k,vO) in dO
                    if  haskey(dN,k)
                        local rr = dN[k]/vO
                        local cv = rank == 0 ? sqrt(Basics.subshell_2j(Subshell(k[2])) + 1.0) : 1.0
                        if      abs(rr - 1.0) < 1.0e-10    nId += 1
                        elseif  abs(rr - cv)  < 1.0e-9     nCv += 1
                        else                               nDs += 1
                        end
                    else
                        # ... only in the old module: the parity rule the new one applies is the whole explanation
                        local la = Basics.subshell_l(Subshell(k[2]));   local lb = Basics.subshell_l(Subshell(k[3]))
                        iseven(la + lb) ? (nUn += 1) : (nPa += 1)
                    end
                end
                for  (k,vN) in dN
                    if  !haskey(dO,k)
                        abs(vN) < 1.0e-14 ? (nZr += 1) : (nEx += 1)
                    end
                end
            end
        end
        @printf("  %-18s %9d %11d %16d %13d %7d %9d %12d %9d\n",
                rank >= 0 ? "1-particle k=$rank" : "2-particle scalar", nId, nCv, nPa, nZr, nEx, nDs, nUn, nDp)
    end
    #

elseif  true
    # Last visit:      25-Aug-2026
    # Last successful:  25-Aug-2026
    #
    # Branch p (MOVING TWO ELECTRONS BETWEEN FOUR DISTINCT SUBSHELLS): the last large topology, and the one whose shape
    #   differs from everything before it. SpinAngular spends twenty-four numbered sub-cases on it
    #   (twoParticle19to42 -> 19to26 / 27to34 / 35to42, each twice).
    #
    #   WHY A SUM APPEARS HERE AND NOWHERE ELSE.  The two-body operator pairs each creation with its own annihilation,
    #   and with four subshells that pairing generally CROSSES the subshell chain -- the two pairs straddle each other.
    #   The CSF's coupling tree runs along the chain, so the operator has to be re-expressed in it, and that costs a SUM
    #   over the chain's one free intermediate rank rather than the phase that sufficed for one electron. Only three
    #   pairing patterns are possible, and only two of them occur for any one family:
    #
    #       :chain    pairs (1,2) and (3,4)   already the chain's own pairing, R forced to k, weight 1/sqrt(2k+1)
    #       :cross13  pairs (1,3) and (2,4)   (-1)^(j2+j3+R+k+1) sqrt(2R+1) { j1 j2 R ; j4 j3 k }
    #       :cross14  pairs (1,4) and (2,3)   (-1)^(j2+j3+k+1)   sqrt(2R+1) { j1 j2 R ; j3 j4 k }
    #
    #   THE THIRD ONE HAS NO R IN ITS PHASE, and that is not a fitted quirk: exchanging j3 and j4 in the :cross13 form
    #   brings a factor (-1)^(j3+j4-R) whose R CANCELS the R already there. Every candidate tried with an R in that
    #   phase failed, including the one symmetry alone would suggest; the algebra says why, and the measurement then
    #   agreed on 43180 of 43180 non-trivial cases.
    #
    #   THE TWO FAMILIES ARE INDEPENDENT CONTRACTIONS, NOT AN EXCHANGE PAIR, and assuming otherwise cost a full cycle.
    #   For the equal-occupation and one-electron cases the second family follows from the first by a Racah sum, because
    #   there both belong to the same pair of subshells. Here R^k(c1,c2,a1,a2) and R^k(c1,c2,a2,a1) are DIFFERENT
    #   integrals over four distinct orbitals, related by no symmetry, and each must be assembled with its own pairing.
    #   Taking the Racah transform left thousands of coefficients wrong and thousands more missing -- and the give-away
    #   was the MISSING count, not the wrong values: the partner has its own rank range, bounded by different j's.
    #
    #   THE JORDAN-WIGNER STRING IS A PREFIX COUNT for four operators -- the occupations lying BEFORE each acting
    #   subshell -- and not the count of gaps between them that the one-electron move uses. The two differ exactly when
    #   the creations and annihilations interleave, which they do in half the arrangements. Carrying the gap form
    #   instead leaves whole classes of pair wrong while leaving others untouched, which is the hardest kind of error to
    #   read from a total.
    #
    #   REPORT (25-Aug-2026), against SpinAngular over eight configuration sets, split by family and pattern:
    #
    #       family    pattern     coefficients   differing   missing   extra
    #       direct    :chain            96 950           0         0       0
    #       direct    :cross13          81 796           0         0       0
    #       crossed   :cross13          61 936           0         0       0
    #       crossed   :cross14         145 996           0         0       0
    #       TOTAL                      386 678           0         0       0
    #
    #   and as complete lists below: 118 728 pairs, 383 812 coefficients, worst ratio 1.000000000000.
    #
    #   NOTHING RAISES ANY MORE, as of 25-Aug-2026: the two remaining topologies -- two electrons moved with ONE
    #   subshell doubled, and with BOTH doubled -- were finished the same day and the module now covers every occupation
    #   pattern a two-body operator can connect. Measured over sixteen configuration sets, all five classes together:
    #
    #       class                          coefficients   differing   missing   extra
    #       equal occupation                    105 283           0         0       0
    #       one electron moved                  127 872           0         0       0
    #       two moved, both doubled               2 828           0         0       0
    #       two moved, one doubled                  526           0         0       0
    #       two moved, four subshells            64 494           0         0       0
    #       TOTAL                               301 003           0         0       0
    #
    localSets = [(["1s^2 2s 2p","1s^2 3s 3p"],        "1s^2 2s2p / 3s3p"),
                 (["2p^2 3d^2","2p 3d 4s 4p"],        "2p^2 3d^2 / 2p3d4s4p"),
                 (["2s^2 3d^2","2s 2p 3d 4f"],        "2s^2 3d^2 / 2s2p3d4f"),
                 (["3s^2 3d^2","3s 3p 3d 4s"],        "3s^2 3d^2 / 3s3p3d4s"),
                 (["3d^2 4f^2","3d 4f 5s 5p"],        "3d^2 4f^2 / 3d4f5s5p"),
                 (["3p^2 4d^2","3p 4d 4f 5s"],        "3p^2 4d^2 / 3p4d4f5s"),
                 (["4d^2 5p^2","4d 5p 5d 5f"],        "4d^2 5p^2 / 4d5p5d5f")]
    tP = 0;  tM = 0;  tN = 0;  tX = 0;  tW = 1.0
    println("\n  configuration            pairs   matched  missing  extra   worst ratio")
    for  (cfgs, tag) in localSets
        local lRel = ConfigurationR[]
        for  c in cfgs   append!(lRel, Basics.generateConfigurations(Basics.RelativisticConfigurations(),
                                                                    Configuration(c)))   end
        local lSub = Basics.generateSubshellList(lRel)
        Defaults.setDefaults("relativistic subshell list", lSub; printout=false)
        local lCsfs = CsfR[]
        for  rc in lRel    append!(lCsfs, Basics.generateCsfRs(rc, lSub))    end
        local np = 0;  local nm = 0;  local nn = 0;  local nx = 0;  local wr = 1.0
        for  l in lCsfs,  r in lCsfs
            l.J == r.J  &&  l.parity == r.parity   ||  continue
            local dd = l.occupation - r.occupation
            count(!=(0), dd) == 4  &&  sum(abs, dd) == 4   ||  continue
            np += 1
            old = [x for x in SpinAngular.computeCoefficients(SpinAngular.TwoParticleOperator(0,Basics.plus,true),
                                                              l, r, lSub)   if abs(x.V) > 1.0e-14]
            new = SpinAngularNew.computeCoefficients(SpinAngularNew.TwoParticleOperator(), l, r, lSub)
            kk(x) = (x.nu, string(x.a), string(x.b), string(x.c), string(x.d))
            dO = Dict{Any,Float64}();   for x in old   dO[kk(x)] = get(dO,kk(x),0.0) + x.V   end
            dN = Dict{Any,Float64}();   for x in new   dN[kk(x)] = get(dN,kk(x),0.0) + x.V   end
            for  (k,v) in dO
                if  haskey(dN,k);   nm += 1;  rr = v/dN[k];  abs(rr-1) > abs(wr-1) && (wr = rr)
                else   nn += 1
                end
            end
            for  k in keys(dN)    haskey(dO,k) || (nx += 1)    end
        end
        @printf("  %-22s %7d %9d %8d %6d   %.12f\n", tag, np, nm, nn, nx, wr)
        global tP += np;  global tM += nm;  global tN += nn;  global tX += nx
        abs(wr-1) > abs(tW-1) && (global tW = wr)
    end
    println("")
    @printf("  TOTAL %d pairs: matched %d, missing %d, extra %d, worst ratio %.12f\n", tP, tM, tN, tX, tW)
    #

end
#
