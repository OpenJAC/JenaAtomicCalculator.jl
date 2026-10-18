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

# the extraSource that turns `DalgarnoLewis.solve`'s own drive into the PURE radial drive r, i.e. R_i - b_i
function unitFix(aSub, vSub, orbs, prim, grid)
    nsL = grid.nsL;  nsS = grid.nsS;  R = zeros(nsL+nsS);  b = zeros(nsL+nsS)
    for i = 1:nsL
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesL[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, full, zeros(grid.NoPoints),
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[i] = RI.rkDiagonal(1, orbs[aSub], syn, grid);   b[i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    for i = 1:nsS
        full = zeros(grid.NoPoints);  bsp = prim.bsplinesS[i];  full[bsp.lower:bsp.upper] = bsp.bs
        syn  = Orbital(vSub, false, true, orbs[aSub].energy, zeros(grid.NoPoints), full,
                       zeros(grid.NoPoints), zeros(grid.NoPoints), Radial.Grid())
        R[nsL+i] = RI.rkDiagonal(1, orbs[aSub], syn, grid)
        b[nsL+i] = IS.eMultipole(1, orbs[aSub], syn, grid)
    end
    return( R - b )
end

# ---- the UNIT-SOURCE perturbed orbital:  (H - eps_a) chi = r * phi_a, with NO angular factor at all -------
#      DL.solve always builds its own source b_i = eMultipole(1, phi_a, B_i);  passing
#      extraSource = R_i - b_i makes the total source exactly R_i, the pure radial drive.
function chiUnit(aSub, vSub, orbs, pot, prim, grid)
    chi = redirect_stdout(devnull) do
        DL.solve(E1, Basics.Babushkin, 0.0, vSub, orbs[aSub], pot, prim; settings=dls,
                 staticDrive=true, extraSource = unitFix(aSub, vSub, orbs, prim, grid)) end
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


# ==========================================================================================================
# CASE C -- ITEM 75, STEP 2b.  The SCREENED alpha of Ne, computed m-RESOLVED, against the module's own.
#
# WHY THIS IS THE STEP THAT CAN DECIDE ANYTHING.  Step 2a compared the module's per-partner WEIGHT with the
# m-resolved one and found the magnitudes identical on all 28 partners and the signs split by Delta j.  A
# weight comparison cannot arbitrate that, because the two routes carry the Wigner-Eckart phase in different
# places and a phase convention is not an observable.  The SCREENED alpha is one.
#
# WHAT IS REPLACED AND WHAT IS REUSED, so that the comparison isolates the angular bookkeeping and nothing
# else.  REUSED, verbatim and deliberately:  `RadialIntegrals.Yk_ab` for the Poisson step,
# `RandomPhase.inducedSourceVector` for the projection onto the B-spline basis, `RandomPhase.gmres` for the
# coupled solve, `DalgarnoLewis.solve` for every channel, and the same grid, orbitals and lambda* as CASE B.
# REPLACED:  every angular factor.  The module carries one `CL_reduced_me` at each end plus a weight of
# sqrt(2j_b+1) and a (-1)^L;  this route carries the explicit substate sum
#     sum_{m_b} <b m_b|C^1_0|t_b m_b> <t_b m_b|C^1_0|b m_b>
# at the density end and NOTHING at the receiving end, because the receiving end's own matrix element
# factors out of the radial equation exactly as it does in the drive.
#
# THREE THINGS ARE FIXED BY DERIVATION HERE, NOT BY A FIT, and each is stated so it can be disputed:
#
#   (1) THE m-SUM SUPPLIES THE 1/(2L+1).  3j orthogonality at fixed middle projection gives
#       sum_m (j_t 1 j_b; -m 0 m)^2 = 1/3 exactly, so the module's /(2L+1) -- which its own docstring records
#       as DERIVED from the Ba+ moment ratio of 3.00000 -- is this orthogonality factor.  The two routes
#       agree about it for a reason, which is worth knowing before any disagreement elsewhere is read.
#
#   (2) THE 4*pi CANCELS COMPLETELY.  The L = 1 component of the density carries (2L+1)/(4pi) and the Poisson
#       solution carries 4pi/(2L+1), so the m-resolved route needs no spherical-harmonic normalisation at all.
#
#   (3) THE DENSITY CARRIES A MINUS SIGN RELATIVE TO THE chi THAT `DalgarnoLewis.solve` RETURNS, and this is
#       fixed by CASE A rather than by hand.  CASE A assembles alpha as +2 <a|r|chi> and reproduces the exact
#       4.5, so JAC's chi is the NEGATIVE of the textbook first-order correction;  the physical density
#       perturbation is therefore -2 sum phi chi, and with that sign the induced field OPPOSES the drive.  A
#       route that screened only because its sign had been chosen to screen would prove nothing.
# ==========================================================================================================

# the m-resolved density weight of the channel (bSub -> tSub): the explicit substate sum over a FILLED
# subshell, with no closed form inserted, so that the 3j orthogonality is measured rather than assumed
function densityWeight(bSub::Subshell, tSub::Subshell)
    jb = Basics.subshell_2j(bSub)//2;   w = 0.0
    for  m  in  -jb:1//1:jb    w = w + meM(bSub, tSub, m) * meM(tSub, bSub, m)    end
    return( w )
end

# THE m-RESOLVED PAIR OF ANGULAR FACTORS.  The density weight is the explicit substate sum;  the MINUS is
# item (3) above, since JAC's chi is the negative of the first-order correction, and with it the induced field
# opposes the drive.  The receiving end carries NOTHING, because its own matrix element factors out of the
# radial equation exactly as the drive's does -- which is the whole content of the m-resolved formulation.
wMres(b, tb) = -densityWeight(b, tb)
aMres(a, ta) = 1.0

# THE MODULE'S PAIR, written in the same radial normalisation so that the two are comparable.  `RandomPhase`
# solves for chi_mod = [CL(a,1,t)/sqrt(2j_a+1)] g, because `DalgarnoLewis.solve` builds its drive from
# `InteractionStrength.eMultipole`;  carrying that scalar through the module's own density weight
# sqrt(2j_b+1)*CL(t_b,1,b)/(2L+1), and through its receiving coefficient aCoeff = -CL(t_a,1,a), gives the two
# factors below.  Note what the receiving one contains and the m-resolved one does not: a sqrt(2j_a+1) that
# comes from `eMultipole`'s normalisation being applied to the drive and NOT to the induced term, although
# both are pieces of the same one-body perturbation in the same equation.
wMod(b, tb) = AM.CL_reduced_me(tb, 1, b) * AM.CL_reduced_me(b, 1, tb) / 3.0
aMod(a, ta) = -sqrt(Basics.subshell_2j(a) + 1.0) * AM.CL_reduced_me(ta, 1, a) / AM.CL_reduced_me(a, 1, ta)


# ==========================================================================================================
# THE COUPLED CASE, PARAMETERISED BY THE SYSTEM IT IS HANDED.  One body of angular code serves every system,
# deliberately: Rule 18's lesson is that a second copy of this algebra drifts, and the Ne numbers below are
# the control that the parameterisation changed nothing.
#
# `uncRef` is the uncoupled alpha this system must reproduce with the induced term switched off (0.0 to skip
# the control), and `refValue`/`refName` the external number both routes are finally compared against.
# A tuple (alpha_uncoupled, alpha_mResolved, alpha_moduleConvention) is returned.
# ==========================================================================================================
function runCoupledCase(label::String, grid, prim, orbs, pot, shells, occ, uncRef::Float64,
                        refValue::Float64, refName::String; showOperator::Bool=true)
    chans = Tuple{Subshell,Subshell}[]
    for  a  in  shells
        (!haskey(orbs, a) || orbs[a].energy >= 0.0)   &&   continue
        for  t  in  RP.reachableChannels(1, a)   push!(chans, (a, t))   end
    end
    msz  = grid.nsL + grid.nsS
    fixV = Dict( ch => unitFix(ch[1], ch[2], orbs, prim, grid)   for ch in chans )
    @printf("\n=== %s:  %d channels, %d splines per channel, %d unknowns ===\n", label, length(chans), msz,
            length(chans)*msz)

    function solveAll(s::Array{Float64,1})
        gs = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()
        for  (k, ch)  in  enumerate(chans)
            gs[ch] = redirect_stdout(devnull) do
                DL.solve(E1, Basics.Babushkin, 0.0, ch[2], orbs[ch[1]], pot, prim; settings=dls,
                         staticDrive=true, extraSource = fixV[ch] + s[(k-1)*msz+1 : k*msz])
            end
        end
        return( gs )
    end

    # ONE APPLICATION OF THE COUPLED MAP, WITH THE ANGULAR BOOKKEEPING AS AN ARGUMENT.  `wFun` is the weight
    # the induced DENSITY gives the channel (b -> t_b) and `aFun` the coefficient the RECEIVING channel
    # (a -> t_a) applies to the projected potential.  Everything else -- the Poisson step, the B-spline
    # projection, the Pauli projection, the grid, the orbitals, the GMRES -- is shared, so a difference
    # between two (wFun, aFun) pairs is a difference in the angular algebra and in nothing else.
    function rebuild(gs, wFun::Function, aFun::Function)
        out = zeros(length(chans)*msz)
        gP  = Dict( ch => projectPauli(gs[ch], ch[2], orbs, shells, grid)   for ch in chans )
        for  (k, ch)  in  enumerate(chans)
            aSub = ch[1];   rho = zeros(grid.NoPoints);   mtp = 1
            for  ch2  in  chans
                b = ch2[1];   tb = ch2[2];   nb = occ[b];   jb2 = Basics.subshell_2j(b)
                w = wFun(b, tb) * nb / (jb2 + 1.0)
                # ONE electron is removed from its own subshell.  This is the m-AVERAGED self-exclusion,
                # which is what (N-1)/N is:  the substate's own share of its subshell's L = 1 density is
                # exactly 1/N of it only after averaging over m, so both routes make the same approximation
                # here and the comparison between them is unaffected by it.
                b == aSub   &&   (w = w * (nb - 1.0) / nb)
                w == 0.0   &&   continue
                g = gP[ch2];   ob = orbs[b]
                mtpB = min(length(ob.P), length(g.P), length(ob.Q), length(g.Q))
                for  i = 1:mtpB   rho[i] = rho[i] + w * 2.0 * (ob.P[i]*g.P[i] + ob.Q[i]*g.Q[i])   end
                mtp = max(mtp, mtpB)
            end
            VindZr = [ -RI.Yk_ab(1, grid.r[i], rho, mtp, grid)   for i = 1:grid.NoPoints ]
            out[(k-1)*msz+1 : k*msz] = aFun(aSub, ch[2]) * RP.inducedSourceVector(VindZr, orbs[aSub], prim)
        end
        return( out )
    end

    function alphaFrom(gs)
        tot = 0.0;   per = Dict{Tuple{Subshell,Subshell},Float64}()
        for  ch  in  chans
            aSub = ch[1];   tSub = ch[2];   ja = Basics.subshell_2j(aSub)//2
            gPp  = projectPauli(gs[ch], tSub, orbs, shells, grid)
            rU   = RI.rkDiagonal(1, orbs[aSub], gPp, grid)
            a    = 0.0
            for  m  in  -ja:1//1:ja   a = a + 2.0 * meM(aSub, tSub, m) * meM(tSub, aSub, m) * rU   end
            per[ch] = a * occ[aSub] / (Basics.subshell_2j(aSub) + 1.0);   tot = tot + per[ch]
        end
        return( tot, per )
    end

    # (i) the UNCOUPLED limit -- the control on the whole plumbing, where an independent value exists
    gs0 = solveAll(zeros(length(chans)*msz))
    a0, per0 = alphaFrom(gs0)
    if  uncRef != 0.0
        @printf(">> CONTROL, induced term OFF:  alpha = %.8f   against %.8f   (diff %.2e)\n", a0, uncRef,
                abs(a0 - uncRef))
    else
        @printf(">> induced term OFF:  alpha = %.8f\n", a0)
    end

    # (ii) the two coupled solves.  The map is affine, F(s) = b + M s, so b = F(0) and (I - M) s = b;  the
    #      iteration's spectral radius is above one here, which is why this is a GMRES solve and not a sweep.
    function coupled(wFun::Function, aFun::Function)
        b0 = rebuild(gs0, wFun, aFun)
        applyA(s) = s - (rebuild(solveAll(s), wFun, aFun) - b0)
        sStar, nKry, relRes = RP.gmres(applyA, b0, 20, 1.0e-8)
        return( solveAll(sStar), nKry, relRes )
    end
    gsMres, nkM, rrM = coupled(wMres, aMres)
    aScr, perScr = alphaFrom(gsMres)
    @printf(">> m-RESOLVED  GMRES:  %d Krylov steps, relative residual %.2e\n", nkM, rrM)
    gsMod, nkD, rrD = coupled(wMod, aMod)
    aModTot, perMod = alphaFrom(gsMod)
    @printf(">> MODULE      GMRES:  %d Krylov steps, relative residual %.2e\n", nkD, rrD)

    # (iii) what the two routes do to ONE element of the coupled operator, so that the bridge is visible as a
    #       number per (receiving, source) pair rather than inferred from the totals alone.
    if  showOperator
        @printf("\n%-21s %-21s %12s %12s %10s\n", "receiving (a -> t_a)", "source (b -> t_b)", "module",
                "m-resolved", "ratio")
        for  cha  in  chans[max(1,length(chans)-2):end]
            for  chb  in  chans[max(1,length(chans)-2):end]
                mM = aMod(cha[1], cha[2]) * wMod(chb[1], chb[2])
                mR = aMres(cha[1], cha[2]) * wMres(chb[1], chb[2])
                @printf("%-21s %-21s %+12.6f %+12.6f %10.6f\n", "$(cha[1]) -> $(cha[2])",
                        "$(chb[1]) -> $(chb[2])", mM, mR, mM/mR)
            end
        end
    end

    @printf("\n%-10s %-10s %14s %14s %14s %10s %10s\n", "a", "t_a", "uncoupled", "m-res scr.", "module scr.",
            "ratio", "sqrt(2ja+1)")
    for  ch  in  chans
        @printf("%-10s %-10s %+14.8f %+14.8f %+14.8f %10.6f %10.6f\n", string(ch[1]), string(ch[2]), per0[ch],
                perScr[ch], perMod[ch], perMod[ch]/perScr[ch], sqrt(Basics.subshell_2j(ch[1]) + 1.0))
    end
    @printf("\n>> %s alpha:  unscreened %.8f   m-resolved screened %.8f   module screened %.8f\n", label, a0,
            aScr, aModTot)
    @printf(">> screening k:  m-resolved %.6f   module %.6f      (k = 1 - alpha_screened/alpha_unscreened)\n",
            1.0 - aScr/a0, 1.0 - aModTot/a0)
    @printf(">> against %s = %.4f:  m-resolved %+.1f %%,  module %+.1f %%\n", refName, refValue,
            100*(aScr/refValue - 1), 100*(aModTot/refValue - 1))

    return( (a0, aScr, aModTot) )
end


# ==========================================================================================================
# CASE C -- Ne at lam*, which is the CONTROL on everything above: its three numbers were measured before the
# case was parameterised and must come back unchanged.
# ==========================================================================================================
occNe = Dict( Subshell("1s_1/2") => 2.0, Subshell("2s_1/2") => 2.0,
              Subshell("2p_1/2") => 2.0, Subshell("2p_3/2") => 4.0 )
runCoupledCase("CASE C -- Ne at lam*", grid, prim, orbsNe, potNe, shellsNe, occNe, totMres, 2.669,
               "the measured alpha(Ne)")

# ==========================================================================================================
# CASE D -- ITEM 75, STEP 2c.  THE CORE POLARIZABILITY OF Ba+ [Xe], WHICH IS WHAT ARBITRATES.
#
# WHY HERE AND NOT ON Ne.  Both routes sit 34 % below Ne's measured 2.669, for reasons that have nothing to
# do with the angular algebra -- the mean field and the missing exchange -- so Ne cannot say which convention
# is right.  Ba+'s core has an external target that CAN: the published RPA alpha_core = 10.6(5), which item
# 73 reached at 9.944 with the module's convention.  Two things are therefore checked at once: whether this
# harness reproduces item 73's numbers (unscreened ~18.16, screened ~9.94) with the module's angular factors,
# which validates the emulation, and where the m-resolved route lands against the same published value.
#
# THE GRID IS `Basics.recommendedGrid`'s, NOT a hand-chosen box -- Rule 12, and the reason it matters here is
# the span: a box sized for the diffuse 5p must still carry the 1s of Z = 56.  The perturbed orbital chi of a
# STATIC drive is bound (it decays as the source does), so the function's continuum caveat does not apply.
# ==========================================================================================================
if  "ba"  in  ARGS
occBaShell = Dict( Shell("1s") => 2, Shell("2s") => 2, Shell("2p") => 6, Shell("3s") => 2, Shell("3p") => 6,
                   Shell("3d") => 10, Shell("4s") => 2, Shell("4p") => 6, Shell("4d") => 10, Shell("5s") => 2,
                   Shell("5p") => 6, Shell("6s") => 1 )
gridBa = Basics.recommendedGrid(occBaShell, 56.0; printout=true)
primBa = Bsplines.generatePrimitives(gridBa)
nmBa   = Nuclear.Model(56.0)
coreBa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"),
           Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2"),
           Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2"), Subshell("4d_3/2"), Subshell("4d_5/2"),
           Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2") ]
occBa  = Dict( sh => Float64(Basics.subshell_2j(sh) + 1)   for sh in coreBa )
asfBa  = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
mpBa   = redirect_stdout(devnull) do
    SelfConsistent.performSCF([Configuration("[Xe] 6s")], nmBa, gridBa, asfBa; printout=false)
end
potBa  = Basics.add( Nuclear.nuclearPotential(nmBa, gridBa),
                     Basics.computePotential(Basics.DFSField(), gridBa, mpBa.levels[1].basis) )
orbsBa = redirect_stdout(devnull) do
    Bsplines.generateOrbitals(coreBa, potBa, nmBa, primBa; printout=false)
end
@printf("\n=== CASE D:  Ba+ [Xe] core on recommendedGrid, rbox = %.2f, %d points, %d+%d splines ===\n",
        gridBa.r[gridBa.NoPoints], gridBa.NoPoints, gridBa.nsL, gridBa.nsS)
# the two spin-orbit partners of the outermost core subshell must agree in mean radius -- Rule 12's second
# guard, and the one that catches a WRONG STATE rather than an inaccurate one
for  sh  in  [Subshell("1s_1/2"), Subshell("4d_3/2"), Subshell("4d_5/2"), Subshell("5p_1/2"),
              Subshell("5p_3/2")]
    @printf("   %-9s  eps = %+14.6f   <r> = %8.4f\n", string(sh), orbsBa[sh].energy,
            RI.rkDiagonal(1, orbsBa[sh], orbsBa[sh], gridBa))
end
runCoupledCase("CASE D -- Ba+ [Xe] core", gridBa, primBa, orbsBa, potBa, coreBa, occBa, 0.0, 10.6,
               "the published RPA alpha_core(Ba+)")
else
    println("\n(CASE D -- Ba+ [Xe] core -- skipped;  pass \"ba\" as an argument to run it, ~15 min.)")
end


# ==========================================================================================================
# CASE E -- ITEM 75, STEP 2d.  THE SIMILARITY INVARIANT, WHICH IS WHAT DECIDES WHETHER THE TWO ROUTES ARE
# TWO CONVENTIONS OR TWO DIFFERENT OPERATORS.
#
# WHY RECIPROCITY WAS THE WRONG TEST, recorded because it was proposed and then killed by its own derivation.
# The coupled operator is `diag(aFun) o K o diag(wFun)` with K -- the radial Coulomb coupling, the Poisson
# step and the B-spline projection -- INDEPENDENT of the angular convention.  Writing the bilinear Coulomb
# form as C[B,A], every block is aFun(B) * C[B,A] * wFun(A), and C is symmetric because the Coulomb kernel is.
# So a route is symmetrizable whenever `wFun(A)/aFun(A)` has the SAME SIGN for every channel, and BOTH routes
# satisfy that:  m-resolved gives -densityWeight(A), negative-definite since densityWeight is a sum of
# squares;  the module gives -CL(a,1,t_a)^2 / (3 sqrt(2j_a+1)), also negative for every channel.  Reciprocity
# therefore holds for both and distinguishes nothing.
#
# WHAT DOES DISTINGUISH THEM.  Matching `diag(aMod) K diag(wMod) = D [K diag(wMres)] D^-1` block by block and
# setting B = A kills the d's completely and leaves
#
#     aFun(A) * wFun(A)   must be EQUAL in the two routes, for every channel A,
#
# and that product is therefore invariant under ANY rescaling of the channel variables -- which is exactly
# what a normalisation convention is free to change.  Measured below.  If the two differ, no choice of
# normalisation connects them and one of the two is wrong as a matter of fact rather than of convention.
# ==========================================================================================================
@printf("\n=== CASE E -- the similarity invariant  aFun(A) * wFun(A)  over the Ne channels ===\n")
@printf("%-10s %-10s %16s %16s %10s %12s\n", "a", "t_a", "module", "m-resolved", "ratio", "sqrt(2ja+1)")
worst = 0.0
for  a  in  shellsNe
    for  t  in  RP.reachableChannels(1, a)
        pMod = aMod(a, t) * wMod(a, t);   pMre = aMres(a, t) * wMres(a, t)
        @printf("%-10s %-10s %+16.10f %+16.10f %10.6f %12.6f\n", string(a), string(t), pMod, pMre,
                pMod/pMre, sqrt(Basics.subshell_2j(a) + 1.0))
        global worst = max(worst, abs(pMod/pMre / sqrt(Basics.subshell_2j(a) + 1.0) - 1.0))
    end
end
@printf(">> the ratio is sqrt(2j_a+1) on every channel, worst deviation %.2e\n", worst)
@printf(">> AND IT IS THE SAME SIGN IN BOTH, so each route is separately symmetrizable and RECIPROCITY CANNOT\n")
@printf(">> DECIDE between them;  the invariant above can, and it says the two are DIFFERENT OPERATORS.\n")
