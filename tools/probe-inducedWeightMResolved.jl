using JenaAtomicCalculator, Printf
const RP = JenaAtomicCalculator.RandomPhase
const AM = AngularMomentum
# =========================================================================================================
# ITEM 75, STEP 2a.  THE PER-PARTNER WEIGHT OF THE DIRECT INDUCED POTENTIAL, two ways.
#
# The module builds the perturbed density as   rho += weightB * 2 * (P_b P_chi + Q_b Q_chi)   with
#     weightB = sqrt(2j_b+1) * CL_reduced_me(t_b, L, b) / (2L+1)
# and then takes -Yk_ab(L, r, rho).  The m-RESOLVED count builds the SAME radial density and the same Yk_ab,
# so every shared convention -- the 4pi/(2L+1), the r factor in Zr, the sign -- cancels in the RATIO, and
# what is left is purely angular:
#     W = sum_{m_b}  <b m|C^L_0|t_b m> * <t_b m|C^L_0|b m>
# because the m-dependence of chi^m is exactly <t_b m|C^L_0|b m> (step 1), and the density pairs it with
# <b m|C^L_0|t_b m>.  No 6j anywhere;  Wigner-Eckart only.
#
# WHAT WOULD CONFIRM THE ITEM'S SUSPICION: weightB/W = sqrt(2j_b+1) on every partner.
# WHAT WOULD REFUTE IT:              weightB/W = 1.
# A third outcome -- some other clean factor -- is also informative, so the ratio is printed, not asserted.
# =========================================================================================================
function meC(aSub::Subshell, bSub::Subshell, L::Int, m::Rational)
    ja = Basics.subshell_2j(aSub)//2;   jb = Basics.subshell_2j(bSub)//2
    (abs(m) > ja || abs(m) > jb) && return 0.0
    w3 = AM.Wigner_3j(AngularJ64(ja), AngularJ64(L), AngularJ64(jb),
                      AngularM64(-m),  AngularM64(0), AngularM64(m))
    return( (-1.0)^Float64(ja - m) * w3 * AM.CL_reduced_me(aSub, L, bSub) )
end

L = 1
occ = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"),
       Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2"),
       Subshell("4f_5/2"), Subshell("4f_7/2")]
@printf("\n%-10s %-10s %16s %16s %12s %10s\n",
        "b", "t_b", "module x chi_unit", "W(m-resolved)", "ratio", "verdict")
println("-"^92)
nSqrt = 0;  nOne = 0;  nOther = 0
for  b  in  occ
    jb2 = Basics.subshell_2j(b)
    for  tb  in  RP.reachableChannels(L, b)
        # THE MODULE'S WEIGHT MULTIPLIES ITS OWN chi, WHICH ALREADY CARRIES ONE CL FACTOR.  `DL.solve` builds
        # the source as eMultipole(1, phi_b, B_i) = CL(b,L,t_b)/sqrt(2j_b+1) * R_i, so
        #     chi_module = [CL(b,L,t_b)/sqrt(2j_b+1)] * chi_unit,
        # while the m-resolved count below multiplies chi_unit, which has NO angular factor at all.  Comparing
        # weightB against W directly is apples to oranges and was the author's first error here;  the
        # like-for-like quantity is the weight AS IT MULTIPLIES chi_unit.
        wRaw = sqrt(jb2 + 1.0) * RP.angularCoefficient(L, tb, b) / (2L + 1)
        wMod = wRaw * AM.CL_reduced_me(b, L, tb) / sqrt(jb2 + 1.0)
        jb = jb2//2
        W = 0.0
        for  m  in  -jb:1//1:jb    W += meC(b, tb, L, m) * meC(tb, b, L, m)    end
        (wMod == 0.0 || W == 0.0) && continue
        ratio = wMod / W;   sq = sqrt(jb2 + 1.0)
        tag = abs(ratio - 1.0) < 1e-10 ? "+1 AGREE" :
              abs(ratio + 1.0) < 1e-10 ? "-1 SIGN"  : "other"
        tag == "-1 SIGN" ? (global nSqrt += 1) : tag == "+1 AGREE" ? (global nOne += 1) : (global nOther += 1)
        @printf("%-10s %-10s %+16.10f %+16.10f %+12.8f %10s\n", string(b), string(tb), wMod, W, ratio, tag)
    end
end
println("-"^92)
@printf(">> %d partners differ by exactly -1 (SIGN),  %d agree exactly (+1),  %d neither\n", nSqrt, nOne, nOther)
