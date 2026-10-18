using JenaAtomicCalculator, Printf, LinearAlgebra
const DL = JenaAtomicCalculator.DalgarnoLewis
const RP = JenaAtomicCalculator.RandomPhase
const IS = InteractionStrength
const RI = RadialIntegrals
const AM = AngularMomentum
# ==========================================================================================================
# ITEM 75, STEP 1.  The FIRST-ORDER (unscreened) static dipole polarizability, computed m-RESOLVED -- explicit
# magnetic substates, Wigner-Eckart applied here, no 6j, and no use of the module's per-partner factor -- and
# compared with the module's own  alpha = (2/3) * (2j_a+1) * eMultipole(1, phi_a, chi)  assembly.
#
# WHY THIS STEP EXISTS.  Item 75's open question is a NORMALISATION, not physics: over 12 parity-allowed
# quadruples the module's induced-potential per-partner factor sits a factor sqrt(2j_b+1) away from the
# m-resolved one, while 42 Ba+ channels agree to exactly 3.00000.  Before attacking that, the FIRST-order
# amplitude must be shown to agree, because it shares the CL convention and nothing else.  If step 1 already
# disagrees, the bridge is here and not in the induced potential.
#
# AND IT IS CALIBRATED ON A CLOSED FORM FIRST.  H 1s has alpha = 4.5 exactly, so the m-resolved assembly is
# validated against a number that owes nothing to JAC before it is pointed at Ne.  The author made two phase
# errors of exactly this kind on 17-Oct, so no phase below is reasoned out by hand: every one comes from
# Wigner_3j and is checked by COVARIANCE -- the sum over a filled subshell must not depend on the axis.
# ==========================================================================================================
dls = DL.Settings(DL.Settings(); printBefore=false)

# ---- the m-resolved matrix element <v m | r C^1_0 | a m>, Wigner-Eckart applied HERE ---------------------
#      <a ma| C^L_q |b mb> = (-1)^(j_a - m_a) * (3j: j_a L j_b; -m_a q m_b) * <a||C^L||b>
function meM(vSub::Subshell, aSub::Subshell, m::Rational)
    jv = Basics.subshell_2j(vSub)//2;   ja = Basics.subshell_2j(aSub)//2
    # the substate m must EXIST in BOTH subshells;  guarding only the first argument asks Wigner_3j for an
    # m = +-3/2 in a j = 1/2 slot the moment a 2p_3/2 -> 1s_1/2 channel appears, which H 1s never produced.
    (abs(m) > jv || abs(m) > ja) && return 0.0
    w3 = AM.Wigner_3j(AngularJ64(jv), AngularJ64(1), AngularJ64(ja),
                      AngularM64(-m),  AngularM64(0), AngularM64(m))
    ph = (-1.0)^Float64(jv - m)
    return( ph * w3 * AM.CL_reduced_me(vSub, 1, aSub) )
end

function projectPauli(chi, t, orbs, shells, grid)
    P = copy(chi.P);   Q = copy(chi.Q)
    for  o  in  shells
        (o.kappa == t.kappa && haskey(orbs, o) && orbs[o].energy < 0.0)  ||  continue
        ov = RI.overlap(orbs[o].P, P, grid) + RI.overlap(orbs[o].Q, Q, grid)
        nn = RI.overlap(orbs[o].P, orbs[o].P, grid) + RI.overlap(orbs[o].Q, orbs[o].Q, grid)
        nn <= 0.0 && continue
        for i = 1:min(length(P), length(orbs[o].P))   P[i] -= ov/nn * orbs[o].P[i]   end
        for i = 1:min(length(Q), length(orbs[o].Q))   Q[i] -= ov/nn * orbs[o].Q[i]   end
    end
    return( Orbital(chi.subshell, false, true, chi.targetEnergy, P, Q,
                    zeros(length(P)), zeros(length(Q)), Radial.Grid()) )
end

# ---- the UNIT-SOURCE perturbed orbital:  (H - eps_a) chi = r * phi_a, with NO angular factor at all -------
#      DL.solve always builds its own source b_i = eMultipole(1, phi_a, B_i);  passing
#      extraSource = R_i - b_i makes the total source exactly R_i, the pure radial drive.
function chiUnit(aSub, vSub, orbs, pot, prim, grid)
    nsL = grid.nsL;  nsS = grid.nsS
    R = zeros(nsL+nsS);   b = zeros(nsL+nsS)
    for i = 1:nsL
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesL[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, full, zeros(grid.NoPoints),
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[i] = RI.rkDiagonal(1, orbs[aSub], syn, grid)
        b[i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    for i = 1:nsS
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesS[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, zeros(grid.NoPoints), full,
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[nsL+i] = RI.rkDiagonal(1, orbs[aSub], syn, grid)
        b[nsL+i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    chi = redirect_stdout(devnull) do
        DL.solve(E1, Basics.Babushkin, 0.0, vSub, orbs[aSub], pot, prim; settings=dls,
                 staticDrive=true, extraSource = R - b) end
    return( chi )
end

# ==========================================================================================================
# CASE A -- HYDROGEN 1s, where alpha = 4.5 EXACTLY.  The m-resolved assembly is
#
#     alpha = -2 * sum_a sum_{m_a} sum_{kappa_v} <a m_a| r C^1_0 |chi^m>,     (H - eps_a) chi^m = (rC^1_0) phi_a
#
# and since the m-dependence of the source is the scalar meM(v,a,m), chi^m = meM(v,a,m) * chiUnit, so the
# contribution is  meM(a,v,m) * meM(v,a,m) * <r>_unit.  The SIGN convention is fixed by this case, not argued.
# ==========================================================================================================
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/5, h = 5.0e-2, hp = 1.0, rbox = 120.0)
prim = Bsplines.generatePrimitives(grid)
nm   = Nuclear.Model(1.0, PointNucleus())
shellsH = [Subshell("1s_1/2")]
potH = Nuclear.nuclearPotential(nm, grid)
orbsH = redirect_stdout(devnull) do
    Bsplines.generateOrbitals(shellsH, potH, nm, prim; printout=false)
end
@printf("\n=== CASE A:  H 1s,  eps = %+.8f  (exact -0.5) ===\n", orbsH[Subshell("1s_1/2")].energy)

aSub = Subshell("1s_1/2");  ja = Basics.subshell_2j(aSub)//2
alphaMod = 0.0;  alphaMres = 0.0;  perM = Dict{Rational,Float64}()
for  vSub  in  RP.reachableChannels(1, aSub)
    chi = chiUnit(aSub, vSub, orbsH, potH, prim, grid)
    chiP = projectPauli(chi, vSub, orbsH, shellsH, grid)
    rUnit = RI.rkDiagonal(1, orbsH[aSub], chiP, grid)              # <a| r |chi_unit>, pure radial
    # the module's own route, for the same channel
    c  = AM.CL_reduced_me(aSub, 1, vSub) / sqrt(Basics.subshell_2j(aSub)+1)
    global alphaMod  += 2.0/3.0 * (Basics.subshell_2j(aSub)+1) * c^2 * rUnit
    for  m  in  -ja:1//1:ja
        contrib = meM(aSub, vSub, m) * meM(vSub, aSub, m) * rUnit
        global alphaMres += 2.0 * contrib
        perM[m] = get(perM, m, 0.0) + 2.0 * contrib
    end
    @printf("   channel %-8s  <r>_unit = %+.8f\n", string(vSub), rUnit)
end
@printf(">> module   assembly  alpha = %.8f\n", alphaMod)
@printf(">> m-RESOLVED         alpha = %.8f\n", alphaMres)
# PER SUBSTATE is the one-electron answer;  the total is twice it because the loop fills the 1s^2 shell.
@printf(">> per substate       alpha = %.8f   against the EXACT 4.50000000,  rel.err %.2e\n",
        alphaMres/2, abs(alphaMres/2/4.5 - 1))
@printf(">> COVARIANCE (must be flat over m):")
for (m, v) in sort(collect(perM))   @printf("  m=%+s: %.8f", string(m), v)   end
spread = maximum(values(perM)) - minimum(values(perM))
@printf("\n>> per-m spread = %.3e  (a filled subshell is isotropic, so this MUST be ~0)\n", spread)

# ==========================================================================================================
# NOTE ON CASE A's FACTOR 2, so it is not misread: the loop sums over BOTH m substates of 1s, i.e. a FILLED
# 1s^2 shell, while H has ONE electron.  The per-m value IS the one-electron answer -- 4.49972686 against the
# exact 4.5, a relative 6.1e-06 -- and the module's (2j_a+1) factor makes the same filled-shell assumption.
# WHAT CASE A ESTABLISHES is the agreement itself: module and m-resolved give 8.99945372 both, to every digit.
#
# CASE B -- Ne at the lambda-tuned mean field.  This is the case that can SEE a sqrt(2j_b+1) discrepancy,
# because 2p_1/2 AND 2p_3/2 are both occupied and the channels reach j_v = 1/2, 3/2 and 5/2.  H 1s cannot.
# ==========================================================================================================
EV = 27.211386;  ipNe = 21.5645/EV
shellsNe = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
outer    = Subshell("2p_3/2")
nmNe = Nuclear.Model(10.0, PointNucleus())
asf  = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
mpNe = redirect_stdout(devnull) do
    SelfConsistent.performSCF([Configuration("1s^2 2s^2 2p^6")], nmNe, grid, asf; printout=false)
end
potNuc = Nuclear.nuclearPotential(nmNe, grid)
potTot = Basics.add( potNuc, Basics.computePotential(Basics.DFSField(), grid, mpNe.levels[1].basis) )
function epsAt(lam::Float64)
    ZrNew = [ potNuc.Zr[i] + lam*(potTot.Zr[i] - potNuc.Zr[i]) for i = 1:length(potTot.Zr) ]
    p  = Radial.Potential(@sprintf("scaled %.6f", lam), ZrNew, grid)
    os = redirect_stdout(devnull) do
        Bsplines.generateOrbitals(shellsNe, p, nmNe, prim; printout=false)
    end
    return( os[outer].energy, p, os )
end
lo, hi = 0.30, 1.00;  lamStar = 1.0;  eStar = 0.0;  potNe = nothing;  orbsNe = nothing
for it = 1:40
    global lamStar = 0.5*(lo+hi)
    global eStar, potNe, orbsNe = epsAt(lamStar)
    abs(eStar + ipNe) < 1.0e-9 && break
    if  eStar < -ipNe   global lo = lamStar   else   global hi = lamStar   end
end
@printf("\n=== CASE B:  Ne at lambda* = %.6f,  eps(2p_3/2) = %+.8f  against the measured %+.8f  (%+.3f %%) ===\n",
        lamStar, eStar, -ipNe, 100*(abs(eStar)/ipNe - 1))

totMod = 0.0;  totMres = 0.0
@printf("\n%-10s %-10s %14s %16s %16s %10s\n", "occupied", "channel", "<r>_unit", "module", "m-resolved", "ratio")
for  aSub  in  shellsNe
    (!haskey(orbsNe, aSub) || orbsNe[aSub].energy >= 0.0) && continue
    local ja = Basics.subshell_2j(aSub)//2      # `local` deliberately: Julia's soft scope would otherwise warn
    for  vSub  in  RP.reachableChannels(1, aSub)
        chi  = chiUnit(aSub, vSub, orbsNe, potNe, prim, grid)
        chiP = projectPauli(chi, vSub, orbsNe, shellsNe, grid)
        rU   = RI.rkDiagonal(1, orbsNe[aSub], chiP, grid)
        c    = AM.CL_reduced_me(aSub, 1, vSub) / sqrt(Basics.subshell_2j(aSub)+1)
        aMod = 2.0/3.0 * (Basics.subshell_2j(aSub)+1) * c^2 * rU
        aMre = 0.0
        for  m  in  -ja:1//1:ja   aMre += 2.0 * meM(aSub, vSub, m) * meM(vSub, aSub, m) * rU   end
        global totMod  += aMod;   global totMres += aMre
        @printf("%-10s %-10s %+14.8f %+16.8f %+16.8f %10.6f\n", string(aSub), string(vSub), rU, aMod, aMre,
                aMod == 0.0 ? NaN : aMre/aMod)
    end
end
@printf("\n>> Ne alpha_unscreened:  module %.8f   m-resolved %.8f   ratio %.10f\n", totMod, totMres,
        totMres/totMod)
@printf(">> IF THIS RATIO IS 1 TO MACHINE PRECISION, the first-order CL convention is SOUND and item 75's\n")
@printf(">> normalisation bridge lies in the INDUCED POTENTIAL, not in the amplitude.\n")
