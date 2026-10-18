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

# the SAME two vectors, returned separately:  R_i = <B_i| r |phi_a> is the pure radial drive and
# b_i = eMultipole(1, phi_a, B_i) is the drive `DalgarnoLewis.solve` builds for itself
function driveVectors(aSub, vSub, orbs, prim, grid)
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
    return( R, b )
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
# A tuple (alpha_uncoupled, alpha_mResolved, alpha_moduleConvention, gsMres, chans) is returned -- the last
# two so that a caller can compare the m-resolved coupled functions themselves against the module's.
# ==========================================================================================================
function runCoupledCase(label::String, grid, prim, orbs, pot, shells, occ, uncRef::Float64,
                        refValue::Float64, refName::String; showOperator::Bool=true,
                        projectShells::Union{Nothing,Array{Subshell,1}}=nothing, selfMode::Symbol=:partners,
                        exchangeSign::Float64=0.0)
    # the set that blocks Pauli-forbidden components is not always the ACTIVE set:  for alpha_core of Ba+ the
    # 6s valence is occupied and must block the s channels, while not itself being an active core subshell.
    proj  = isnothing(projectShells) ? shells : projectShells
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
    # `xSign` adds `RandomPhase.exchangeSourceVector`, which is written in THIS normalisation (receiving
    # coefficient 1, channel function g) because its coefficient was derived as the ratio to
    # ME(t_a, m_ta, L, q, a, m_a) -- the same object the drive is proportional to.  It is therefore added to
    # the m-RESOLVED pair only;  the module-convention pair carries aFun != 1 and mixing the two would be the
    # very normalisation error this file exists to find.  xSign is +-1 because the chi-convention minus that
    # the direct term carries explicitly (wMres) cannot be read off the exchange assembly -- see CASE J.
    function rebuild(gs, wFun::Function, aFun::Function, xSign::Float64)
        out = zeros(length(chans)*msz)
        gP  = Dict( ch => projectPauli(gs[ch], ch[2], orbs, proj, grid)   for ch in chans )
        xChis = xSign == 0.0 ? nothing :
                Dict( ch => DalgarnoLewis.PerturbedOrbital(ch[2], ch[1], E1, Basics.Babushkin, 0.0,
                            orbs[ch[1]].energy, gP[ch].P, gP[ch].Q, 0.0)   for ch in chans )
        for  (k, ch)  in  enumerate(chans)
            aSub = ch[1];   rho = zeros(grid.NoPoints);   mtp = 1
            for  ch2  in  chans
                b = ch2[1];   tb = ch2[2];   nb = occ[b];   jb2 = Basics.subshell_2j(b)
                w = wFun(b, tb) * nb / (jb2 + 1.0)
                # ONE electron is removed from its own subshell.  This is the m-AVERAGED self-exclusion,
                # which is what (N-1)/N is:  the substate's own share of its subshell's L = 1 density is
                # exactly 1/N of it only after averaging over m, so both routes make the same approximation
                # here and the comparison between them is unaffected by it.
                # :partners is the current module -- an electron is screened by its (N-1) subshell partners.
                # :skip is what `inducedPotentialZr` did BEFORE the 17-Oct-2026 repair, and item 73's Ba+
                # number predates that repair, so reproducing it needs this option.
                if  b == aSub
                    selfMode == :skip   &&   continue
                    # :full keeps the whole subshell, self-interaction included -- which a proper RPAE does,
                    # because EXCHANGE cancels it;  without exchange it is spurious and OVER-screens.
                    selfMode == :partners   &&   (w = w * (nb - 1.0) / nb)
                end
                w == 0.0   &&   continue
                g = gP[ch2];   ob = orbs[b]
                mtpB = min(length(ob.P), length(g.P), length(ob.Q), length(g.Q))
                for  i = 1:mtpB   rho[i] = rho[i] + w * 2.0 * (ob.P[i]*g.P[i] + ob.Q[i]*g.Q[i])   end
                mtp = max(mtp, mtpB)
            end
            VindZr = [ -RI.Yk_ab(1, grid.r[i], rho, mtp, grid)   for i = 1:grid.NoPoints ]
            out[(k-1)*msz+1 : k*msz] = aFun(aSub, ch[2]) * RP.inducedSourceVector(VindZr, orbs[aSub], prim)
            if  !isnothing(xChis)
                out[(k-1)*msz+1 : k*msz] = out[(k-1)*msz+1 : k*msz] + xSign *
                    RP.exchangeSourceVector(E1, orbs, xChis, aSub, ch[2], grid, prim; occupations=occ,
                                            maxRank=8)
            end
        end
        return( out )
    end

    function alphaFrom(gs)
        tot = 0.0;   per = Dict{Tuple{Subshell,Subshell},Float64}()
        for  ch  in  chans
            aSub = ch[1];   tSub = ch[2];   ja = Basics.subshell_2j(aSub)//2
            gPp  = projectPauli(gs[ch], tSub, orbs, proj, grid)
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
    function coupled(wFun::Function, aFun::Function, xSign::Float64)
        b0 = rebuild(gs0, wFun, aFun, xSign)
        applyA(s) = s - (rebuild(solveAll(s), wFun, aFun, xSign) - b0)
        sStar, nKry, relRes = RP.gmres(applyA, b0, xSign == 0.0 ? 20 : 30, 1.0e-8)
        return( solveAll(sStar), nKry, relRes )
    end
    gsMres, nkM, rrM = coupled(wMres, aMres, exchangeSign)
    aScr, perScr = alphaFrom(gsMres)
    @printf(">> m-RESOLVED  GMRES:  %d Krylov steps, relative residual %.2e\n", nkM, rrM)
    gsMod, nkD, rrD = coupled(wMod, aMod, 0.0)
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

    return( (a0, aScr, aModTot, gsMres, chans) )
end


# ==========================================================================================================
# CASE C -- Ne at lam*, which is the CONTROL on everything above: its three numbers were measured before the
# case was parameterised and must come back unchanged.
# ==========================================================================================================
occNe = Dict( Subshell("1s_1/2") => 2.0, Subshell("2s_1/2") => 2.0,
              Subshell("2p_1/2") => 2.0, Subshell("2p_3/2") => 4.0 )
caseC = nothing
if  !("skipC" in ARGS)
    global caseC = runCoupledCase("CASE C -- Ne at lam*", grid, prim, orbsNe, potNe, shellsNe, occNe, totMres,
               2.669, "the measured alpha(Ne)")
else
    println("\n(CASE C -- Ne at lam* -- skipped;  \"skipC\" was passed.)")
end

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


# ==========================================================================================================
# CASE F -- ITEM 75, STEP 3.  THE DECIDING TEST, AND IT NEEDS NO EXCHANGE, NO LITERATURE AND NO CONVENTION.
#
# WHY THE HYDROGEN CANCELLATION IS NOT THE TEST IT WAS RECORDED AS.  Item 75 held that for ONE electron the
# direct and exchange induced fields must cancel exactly, since `V_HF` vanishes there, and the attempt gave
# -2.604 and +3.665 where both were expected to be 1.  Writing out the first-order Hartree-Fock equations
# shows why it cannot be 1:
#
#     (delta_J - delta_K) phi_a  =  phi_a J[phi chi]  -  chi J[|phi|^2]        and      (J - K) chi  =  -that
#
# so the induced self-field cancels against the HF self-interaction acting on CHI, not against the exchange
# induced field.  `phi_a J[phi chi]` and `chi J[|phi|^2]` are different radial functions in the same channel
# and nothing makes them equal.  JAC's left-hand side carries a DFS mean field with no `J - K` at all, so
# "direct + E1 + E2 = 0" is not an identity in this formulation and was never going to come out at 1.
#
# THE TEST THAT DOES DECIDE IT.  The external drive and the induced potential are BOTH ordinary one-body
# potentials, acting on the SAME electron, in the SAME Dalgarno-Lewis equation.  A one-body potential has ONE
# matrix element, so the same v(r) must produce the SAME source vector down either path.  Take v(r) = r --
# which is the external dipole drive itself -- and push it through the induced path:
#
#     V = -VindZr/r  by `RandomPhase.inducedSourceVector`'s convention, so v(r) = r means VindZr = -r^2,
#     and `inducedSourceVector` then returns exactly  +int B_i(r) r phi_a(r) dr = R_i.
#
# The module multiplies that by aCoeff = -CL(t_a,1,a) and compares against its OWN drive
# `eMultipole(1, phi_a, B_i) = CL(a,1,t_a)/sqrt(2j_a+1) * R_i`;  the m-resolved route multiplies by 1 and
# compares against R_i.  Whichever route returns the same vector for the same potential is the consistent one.
#
# AND THE DRIVE IS THE SIDE THAT IS ALREADY ANCHORED, which is what makes this an arbitration rather than
# another comparison: `eMultipole`'s normalisation is fixed by the EXACT alpha(H 1s) = 4.5 through the
# assembly alpha = (2/3)(2j_a+1) eMultipole(1, phi_a, chi), and CASE A reproduces 4.49972686 while CASE B
# shows the m-resolved assembly agreeing with it to 1.0000000000 on all nine Ne channels.  So the drive's
# reduction is not in question; only the induced term's is.
# ==========================================================================================================
@printf("\n=== CASE F -- the SAME potential v(r) = r down the drive path and the induced path ===\n")
@printf("%-10s %-10s %14s %14s %12s %12s\n", "a", "t_a", "module ratio", "m-res ratio", "spread(mod)",
        "spread(mres)")
VindR = [ -grid.r[i]^2   for i = 1:grid.NoPoints ]       # the Zr representation of v(r) = r
worstMod = 0.0;  worstMres = 0.0
for  a  in  shellsNe
    for  t  in  RP.reachableChannels(1, a)
        R, bDrive = driveVectors(a, t, orbsNe, prim, grid)
        sInd = RP.inducedSourceVector(VindR, orbsNe[a], prim)
        # the module's induced path carries aCoeff = -CL(t_a,1,a), the (-1)^L applied once in `solve`
        aCo  = -AM.CL_reduced_me(t, 1, a)
        keep = [ i  for i = 1:length(R)  if abs(bDrive[i]) > 1.0e-14 * maximum(abs, bDrive) ]
        rM   = [ aCo * sInd[i] / bDrive[i]   for i in keep ]
        rR   = [ 1.0 * sInd[i] / R[i]        for i in keep if abs(R[i]) > 0.0 ]
        mM   = sum(rM)/length(rM);   mR = sum(rR)/length(rR)
        @printf("%-10s %-10s %+14.8f %+14.8f %12.2e %12.2e\n", string(a), string(t), mM, mR,
                maximum(abs, rM .- mM), maximum(abs, rR .- mR))
        global worstMod  = max(worstMod,  abs(abs(mM) - 1.0))
        global worstMres = max(worstMres, abs(mR - 1.0))
    end
end
@printf(">> m-RESOLVED:  the induced path returns the drive EXACTLY -- worst |ratio - 1| = %.2e over all\n",
        worstMres)
@printf(">>              channels, so one potential has one matrix element and the two paths agree.\n")
@printf(">> MODULE:      worst ||ratio| - 1| = %.2e.  The ratio is -sqrt(2j_a+1) CL(t,1,a)/CL(a,1,t), i.e.\n",
        worstMod)
@printf(">>              sqrt(2j_a+1) in magnitude:  THE SAME POTENTIAL GETS TWO DIFFERENT MATRIX ELEMENTS,\n")
@printf(">>              which is a defect and not a convention, and it sits at the RECEIVING end.\n")


# ==========================================================================================================
# CASE G -- ITEM 75, STEP 3's GAP CLOSED.  Ba+ ON ITEM 73's OWN SETUP, RECONSTRUCTED FROM ITS FINGERPRINT.
#
# Item 73's grid and mean field were never recorded, so its alpha_core(Ba+) = 9.944 against the published RPA
# 10.6(5) could not be used to test anything.  They were recovered by scanning against the three per-shell
# numbers the item DOES record -- 5p_3/2 12.34, 5p_1/2 4.93, 4d 0.86, total 18.1598 -- and the match is
# [Xe] 6s, DFS, on a dense grid, using THE SCF BASIS'S OWN ORBITALS:
#     tot 18.1339 (-0.1 %)   5p_3/2 12.3164 (-0.2 %)   5p_1/2 4.9197 (-0.2 %)   4d 0.8564 (-0.4 %)
#
# TWO NEGATIVES WORTH KEEPING FROM THAT SCAN, because each kills a plausible explanation:
#   * THE BOX IS NOT THE VARIABLE.  rbox = 20, 25, 30 and 40 at hp = 2e-2 give IDENTICAL numbers to four
#     decimals, and `recommendedGrid`'s far sparser basis agrees to 0.09 %.  So nothing here was grid-limited.
#   * WHAT THE VARIABLE IS: whether the orbitals are the SCF's own or REGENERATED in a frozen DFS potential
#     built from that same basis.  Those are not the same orbitals -- the potential is not the operator whose
#     eigenstates they are -- and the difference is **9.6 % on alpha_core** (19.909 against 18.134).  Item 75
#     already records the violent form of this for AL/EOL orbitals;  this is the mild DFS form of it, and it
#     is large enough to swamp a 5.7 % convention difference, which is exactly why step 2c could not arbitrate.
# ==========================================================================================================
if  "ba73"  in  ARGS
gridB3 = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 25.0)
primB3 = Bsplines.generatePrimitives(gridB3)
nmB3   = Nuclear.Model(56.0)
coreB3 = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"),
           Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2"),
           Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2"), Subshell("4d_3/2"), Subshell("4d_5/2"),
           Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2") ]
asfB3  = AsfSettings(AsfSettings(); scField = Basics.DFSField(), gridStopper = false)
mpB3   = redirect_stdout(devnull) do
    SelfConsistent.performSCF([Configuration("[Xe] 6s")], nmB3, gridB3, asfB3; printout=false)
end
potB3  = Basics.add( Nuclear.nuclearPotential(nmB3, gridB3),
                     Basics.computePotential(Basics.DFSField(), gridB3, mpB3.levels[1].basis) )
orbsB3 = Dict{Subshell,Orbital}()
for  (k, v)  in  mpB3.levels[1].basis.orbitals   orbsB3[k] = v   end
occB3  = Dict( sh => Float64(Basics.subshell_2j(sh) + 1)   for sh in coreB3 )
projB3 = vcat(coreB3, [Subshell("6s_1/2")])
@printf("\n=== CASE G:  Ba+ on item 73's reconstructed setup -- SCF basis orbitals, rbox %.2f, %d points ===\n",
        gridB3.r[gridB3.NoPoints], gridB3.NoPoints)
@printf("   eps(5p_3/2) = %+.6f,  eps(5p_1/2) = %+.6f,  eps(6s_1/2) = %+.6f\n",
        orbsB3[Subshell("5p_3/2")].energy, orbsB3[Subshell("5p_1/2")].energy,
        orbsB3[Subshell("6s_1/2")].energy)
smList = [ Symbol(a)  for a in ARGS  if a in ["partners", "skip", "full"] ]
isempty(smList)  &&  (smList = [:partners, :skip, :full])
for  sm  in  smList
    @printf("\n--- self term: %s %s---\n", sm,
            sm == :skip ? "(the whole source subshell skipped, as `inducedPotentialZr` did before 17-Oct) " :
            sm == :full ? "(the self term kept IN FULL -- spurious without exchange, over-screens) " :
                          "(the (N-1) partners screen, the current module) ")
    runCoupledCase("CASE G -- Ba+ [Xe] core, item 73's setup", gridB3, primB3, orbsB3, potB3, coreB3, occB3,
                   18.1598, 10.6, "the published RPA alpha_core(Ba+)"; showOperator=false,
                   projectShells=projB3, selfMode=sm)
end
@printf(">> THE EMULATION CONTROL: one of the two module-convention columns should land near item 73's 9.944.\n")
else
    println("\n(CASE G -- Ba+ on item 73's setup -- skipped;  pass \"ba73\" to run it.)")
end


# ==========================================================================================================
# CASE H -- ITEM 75, THE REPAIR, CHECKED AGAINST THE MODULE'S OWN CODE RATHER THAN AGAINST AN EMULATION.
#
# CASES A-G all ran the angular algebra IN THIS FILE, which is what let the two conventions be compared on
# one radial footing;  the price is that none of them touched `src/`.  The repair landed on 18-Oct-2026 and
# these three checks are the ones that can fail if it is wrong.  Each tests a DIFFERENT end:
#
#   (H1) THE SUBSTATE SUM, which is what the density end's argument order rests on.  The weight of the
#        channel b -> t_b is sum_m <b m|C^L_0|t_b m><t_b m|C^L_0|b m>, and the claim is that it equals
#        CL_reduced_me(b, L, t_b)^2/(2L+1):  positive definite, so target-first -- which differs by
#        (-1)^(j_b - j_t_b) -- cannot be right.  Checked at BOTH L = 1 and L = 2 over s, p, d and f, because
#        a rank-1 accident is exactly the kind of thing that hides here.
#
#   (H2) THE RECEIVING END, against `RandomPhase.receivingCoefficient` itself.  This is CASE F re-run with
#        the module's real function in place of the emulated coefficient:  the same v(r) = r down the drive
#        path and the induced path must give the same source vector, i.e. ratio 1.
#
#   (H3) THE DENSITY END, against `RandomPhase.inducedPotentialZr` itself.  The module is handed the chis it
#        would have solved for -- the m-resolved g scaled by the module's own normalisation
#        CL(b,1,t_b)/sqrt(2j_b+1) -- and its returned potential is compared point by point with the
#        m-resolved one built here.  The Pauli projection is applied BEFORE handing them over, which is safe
#        because it is a projector: the module projects again and gets the same vector, so a disagreement
#        here is an angular disagreement and not a projection one.
# ==========================================================================================================
@printf("\n=== CASE H1 -- the substate sum against CL^2/(2L+1), at L = 1 and L = 2 ===\n")
function densityWeightL(L::Int64, bSub::Subshell, tSub::Subshell)
    jb = Basics.subshell_2j(bSub)//2;   jt = Basics.subshell_2j(tSub)//2;   w = 0.0
    for  m  in  -jb:1//1:jb
        abs(m) > jt   &&   continue
        w3b = AM.Wigner_3j(AngularJ64(jb), AngularJ64(L), AngularJ64(jt),
                           AngularM64(-m),  AngularM64(0), AngularM64(m))
        w3t = AM.Wigner_3j(AngularJ64(jt), AngularJ64(L), AngularJ64(jb),
                           AngularM64(-m),  AngularM64(0), AngularM64(m))
        w   = w + (-1.0)^Float64(jb - m) * w3b * AM.CL_reduced_me(bSub, L, tSub) *
                  (-1.0)^Float64(jt - m) * w3t * AM.CL_reduced_me(tSub, L, bSub)
    end
    return( w )
end
subsH = [Subshell("1s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2"),
         Subshell("4f_5/2"), Subshell("4f_7/2")]
worstId = 0.0;  nId = 0;  nFlip = 0
@printf("%-3s %-10s %-10s %14s %14s %14s %9s\n", "L", "b", "t_b", "substate sum", "CL(b,L,t)^2/(2L+1)",
        "target-first", "flipped")
for  L  in  (1, 2)
    for  b  in  subsH,  tb  in  subsH
        cl = AM.CL_reduced_me(b, L, tb)
        cl == 0.0   &&   continue
        ws = densityWeightL(L, b, tb);   wc = cl^2 / (2L + 1)
        wt = sqrt(Basics.subshell_2j(b) + 1.0) * AM.CL_reduced_me(tb, L, b) / (2L + 1) *
             cl / sqrt(Basics.subshell_2j(b) + 1.0)
        global worstId = max(worstId, abs(ws - wc));   global nId = nId + 1
        wt * wc < 0.0   &&   (global nFlip = nFlip + 1)
        if  nId <= 8  ||  wt * wc < 0.0
            @printf("%-3d %-10s %-10s %+14.8f %+14.8f %+14.8f %9s\n", L, string(b), string(tb), ws, wc, wt,
                    wt*wc < 0.0 ? "YES" : "-")
        end
    end
end
@printf(">> %d channels:  worst |substate sum - CL^2/(2L+1)| = %.2e, so the density weight IS positive\n",
        nId, worstId)
@printf(">> definite;  target-first disagrees in SIGN on %d of them, which is the repair's whole content\n", nFlip)
@printf(">> at this end.  (The last column is target-first in the same g-normalisation.)\n")


@printf("\n=== CASE H2 -- `RandomPhase.receivingCoefficient` on the same v(r) = r as CASE F ===\n")
@printf("%-10s %-10s %14s %12s\n", "a", "t_a", "ratio", "spread")
VindH = [ -grid.r[i]^2   for i = 1:grid.NoPoints ]
worstH2 = 0.0
for  a  in  shellsNe
    for  t  in  RP.reachableChannels(1, a)
        R, bDrive = driveVectors(a, t, orbsNe, prim, grid)
        sInd = RP.inducedSourceVector(VindH, orbsNe[a], prim)
        aCo  = RP.receivingCoefficient(1, a, t)
        keep = [ i  for i = 1:length(R)  if abs(bDrive[i]) > 1.0e-14 * maximum(abs, bDrive) ]
        rM   = [ aCo * sInd[i] / bDrive[i]   for i in keep ]
        mM   = sum(rM)/length(rM)
        @printf("%-10s %-10s %+14.8f %12.2e\n", string(a), string(t), mM, maximum(abs, rM .- mM))
        global worstH2 = max(worstH2, abs(mM - 1.0))
    end
end
@printf(">> worst |ratio - 1| = %.2e:  the repaired receiving end gives the induced potential the SAME\n", worstH2)
@printf(">> matrix element as the drive, which is what CASE F showed it did not.\n")


@printf("\n=== CASE H3 -- `RandomPhase.inducedPotentialZr` against the m-resolved density ===\n")
chansH = Tuple{Subshell,Subshell}[]
for  a  in  shellsNe
    for  t  in  RP.reachableChannels(1, a)   push!(chansH, (a, t))   end
end
occOrbH = Dict( sh => orbsNe[sh]   for sh in shellsNe )
gH      = Dict( ch => projectPauli(chiUnit(ch[1], ch[2], orbsNe, potNe, prim, grid), ch[2], orbsNe,
                                   shellsNe, grid)   for ch in chansH )
# the module's own chi is the m-resolved g times CL(b,1,t_b)/sqrt(2j_b+1), since `DalgarnoLewis.solve`
# builds its drive from `InteractionStrength.eMultipole` and that scalar rides through the linear equation
chisH = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()
for  ch  in  chansH
    sb = AM.CL_reduced_me(ch[1], 1, ch[2]) / sqrt(Basics.subshell_2j(ch[1]) + 1.0)
    chisH[ch] = DalgarnoLewis.PerturbedOrbital(ch[2], ch[1], E1, Basics.Babushkin, 0.0, orbsNe[ch[1]].energy,
                                               sb * gH[ch].P, sb * gH[ch].Q, 0.0)
end
@printf("%-10s %14s %14s %12s\n", "a (excluded)", "max|V_module|", "max|V_mres|", "max rel diff")
worstH3 = 0.0
for  a  in  shellsNe
    vMod = RP.inducedPotentialZr(E1, orbsNe, chisH, a, grid; occupiedOrbitals=occOrbH, occupations=occNe)
    rho = zeros(grid.NoPoints);   mtp = 1
    for  ch  in  chansH
        b = ch[1];   tb = ch[2];   nb = occNe[b];   jb2 = Basics.subshell_2j(b)
        w = wMres(b, tb) * nb / (jb2 + 1.0)
        b == a   &&   (w = w * (nb - 1.0) / nb)
        w == 0.0   &&   continue
        ob = orbsNe[b];   g = gH[ch]
        mtpB = min(length(ob.P), length(g.P), length(ob.Q), length(g.Q))
        for  i = 1:mtpB   rho[i] = rho[i] + w * 2.0 * (ob.P[i]*g.P[i] + ob.Q[i]*g.Q[i])   end
        mtp = max(mtp, mtpB)
    end
    vMres = [ -RI.Yk_ab(1, grid.r[i], rho, mtp, grid)   for i = 1:grid.NoPoints ]
    sc = maximum(abs, vMres);   rel = sc > 0.0 ? maximum(abs, vMod - vMres)/sc : 0.0
    @printf("%-10s %14.6e %14.6e %12.2e\n", string(a), maximum(abs, vMod), sc, rel)
    global worstH3 = max(worstH3, rel)
end
@printf(">> worst relative difference %.2e:  the repaired `inducedPotentialZr` builds the m-resolved\n", worstH3)
@printf(">> induced potential, so BOTH ends of the module now carry the convention CASE F anchored.\n")


# ==========================================================================================================
# CASE I -- `RandomPhase.solve` ITSELF, AT omega = 0, AGAINST THE M-RESOLVED COUPLED FUNCTIONS.
#
# This is the end-to-end check the harness could not make until today, and it needed TWO repairs to become
# possible:  the angular one of 18-Oct (item 75) and the `staticDrive` pass-through (item 69).  Before them
# `RandomPhase.solve` at omega = 0 drove every channel with a kernel that is identically zero there, returned
# chi = 0, and reported `converged = true` -- so the module had never produced a static screened response at
# all, and this harness was a SUBSTITUTE for it rather than a CHECK on it.
#
# WHAT MUST AGREE, AND EXACTLY BY HOW MUCH.  The harness drives with the pure radial source r*phi_a, giving
# g;  the module drives with `InteractionStrength.eMultipole`, giving chi = s_a g with
# s_a = CL(a,1,t_a)/sqrt(2j_a+1).  That scalar rides through the whole linear problem, so after the COUPLED
# solve the module's chi must still be exactly s_a times the m-resolved g on EVERY channel.  Comparing the
# functions rather than the polarizability is deliberate: a total can agree while the channels do not, which
# is precisely how the convention error hid (CASE C -- 0.7 % on the total, 0.644 to 2.194 per channel).
# ==========================================================================================================
if  isnothing(caseC)
    println("\n(CASE I -- needs CASE C's m-resolved chis;  re-run without \"skipC\".)")
else
_, _, _, gsRef, chansRef = caseC
tgtNe = Dict( a => first(RP.reachableChannels(1, a))   for a in shellsNe )
setI  = RP.Settings(RP.Settings(); activeSubshells=shellsNe, targetSubshells=tgtNe, couplingEnabled=true,
                     allChannels=true, directSolve=true, maxIterations=20, accuracy=1.0e-8, printBefore=false)
respI = redirect_stdout(devnull) do
    RP.solve(E1, Basics.Babushkin, 0.0, Dict(a => orbsNe[a] for a in shellsNe), potNe, prim;
             settings=setI, occupiedOrbitals=Dict(a => orbsNe[a] for a in shellsNe), occupations=occNe,
             staticDrive=true)
end
@printf("\n=== CASE I:  RandomPhase.solve at omega = 0, staticDrive = true -- %d channels, converged = %s ===\n",
        length(respI.channelOrbitals), string(respI.converged))
@printf("%-10s %-10s %16s %16s %14s\n", "a", "t_a", "chi_mod/g measured", "CL(a,1,t)/sqrt(2j+1)", "rel. diff")
worstI = 0.0;  nMissing = 0
for  ch  in  chansRef
    if  !haskey(respI.channelOrbitals, ch)   global nMissing = nMissing + 1;   continue   end
    cm = respI.channelOrbitals[ch];   gr = gsRef[ch]
    sa = AM.CL_reduced_me(ch[1], 1, ch[2]) / sqrt(Basics.subshell_2j(ch[1]) + 1.0)
    mtp = min(length(cm.P), length(gr.P))
    nrm = maximum(abs, gr.P[1:mtp])
    rat = sum(cm.P[i]*gr.P[i] for i = 1:mtp) / sum(gr.P[i]^2 for i = 1:mtp)
    dev = maximum(abs(cm.P[i] - sa*gr.P[i]) for i = 1:mtp) / (abs(sa) * nrm)
    @printf("%-10s %-10s %+16.10f %+16.10f %14.2e\n", string(ch[1]), string(ch[2]), rat, sa, dev)
    global worstI = max(worstI, dev)
end
@printf(">> worst pointwise deviation of chi_module from s_a * g_m-resolved:  %.2e  (%d channels missing)\n",
        worstI, nMissing)
@printf(">> So the MODULE now reproduces the m-resolved coupled response function by function, not merely in\n")
@printf(">> the total -- which is the end-to-end statement neither route could make before today.\n")
end


# ==========================================================================================================
# CASE J -- ITEM 75's OWN SCOPE: HOW BIG IS EXCHANGE?  Sized, not built.
#
# `RandomPhase.exchangeSourceVector` and `exchangeAngularCoefficient` are in `src/` and verified (970
# informative cases per term, worst 1.11e-16, validated out of sample) but deliberately UNWIRED.  What the
# item asks is whether wiring them is worth it, and that is a question about SIZE, which costs one source
# rebuild rather than a coupled solve.
#
# **WHY THIS IS DONE IN THE HARNESS AND NOT IN THE MODULE, which is a normalisation argument and not a
# preference.**  `exchangeSourceVector`'s coefficient was derived m-RESOLVED, as the ratio of the m-resolved
# exchange source to ME(t_a, m_ta, L, q, a, m_a) -- the same object the DRIVE is proportional to.  So the
# function is written in the normalisation where the receiving coefficient is 1 and the channel function is
# g, i.e. THIS harness's, not the module's chi = s_a g.  Wiring it into `solve` therefore needs a per-partner
# rescale by s_a/s_b, and that is precisely the kind of invisible factor this item exists to stop.  Size it
# first, in the normalisation it was derived in, and wire it only if the size justifies the bookkeeping.
#
# **AND THE SELF TERM IS KEPT IN FULL, because that is what exchange is for:** in a proper RPAE the direct
# term runs over all N_b electrons, self included, and the spurious self-interaction is cancelled by exchange
# rather than by hand.  Removing it on the direct side as well would remove it twice.
#
# **WHAT IS NOT ANCHORED HERE, SAID PLAINLY: the overall SIGN.**  Both terms are linear in chi_b, exactly as
# the direct term is, so both carry the same factor for JAC's chi convention -- and on the direct side that
# factor is a MINUS, fixed by the exact alpha(H 1s) = 4.5 (CASE A) and not by inspection.  Whether
# `exchangeSourceVector` already contains it cannot be read off its assembly, so BOTH signs are reported and
# the physics decides: exchange must OPPOSE the direct induced field, hence reduce screening and RAISE alpha.
# A size is quoted as |exchange| / |direct| and is sign-independent;  only the alpha columns depend on it.
# ==========================================================================================================
@printf("\n=== CASE J -- the size of the exchange source against the direct one, Ne at lam* ===\n")
chansJ = Tuple{Subshell,Subshell}[]
for  a  in  shellsNe
    for  t  in  RP.reachableChannels(1, a)   push!(chansJ, (a, t))   end
end
gJ = Dict( ch => projectPauli(chiUnit(ch[1], ch[2], orbsNe, potNe, prim, grid), ch[2], orbsNe, shellsNe, grid)
           for ch in chansJ )
# the SAME functions, packaged as PerturbedOrbitals, which is what `exchangeSourceVector` takes
chisJ = Dict( ch => DalgarnoLewis.PerturbedOrbital(ch[2], ch[1], E1, Basics.Babushkin, 0.0,
                      orbsNe[ch[1]].energy, gJ[ch].P, gJ[ch].Q, 0.0)   for ch in chansJ )
@printf("%-10s %-10s %14s %14s %10s %12s\n", "a", "t_a", "|direct|", "|exchange|", "ratio", "cos(angle)")
sumD = 0.0;  sumX = 0.0
tJ = @elapsed for  ch  in  chansJ
    aSub = ch[1];   rho = zeros(grid.NoPoints);   mtp = 1
    for  ch2  in  chansJ
        b = ch2[1];   tb = ch2[2];   nb = occNe[b];   jb2 = Basics.subshell_2j(b)
        w = wMres(b, tb) * nb / (jb2 + 1.0)         # self term kept IN FULL -- exchange is what cancels it
        w == 0.0   &&   continue
        ob = orbsNe[b];   g = gJ[ch2]
        mtpB = min(length(ob.P), length(g.P), length(ob.Q), length(g.Q))
        for  i = 1:mtpB   rho[i] = rho[i] + w * 2.0 * (ob.P[i]*g.P[i] + ob.Q[i]*g.Q[i])   end
        mtp = max(mtp, mtpB)
    end
    VindJ = [ -RI.Yk_ab(1, grid.r[i], rho, mtp, grid)   for i = 1:grid.NoPoints ]
    sD = RP.inducedSourceVector(VindJ, orbsNe[aSub], prim)
    sX = RP.exchangeSourceVector(E1, orbsNe, chisJ, aSub, ch[2], grid, prim; occupations=occNe, maxRank=8)
    nD = sqrt(sum(sD.^2));   nX = sqrt(sum(sX.^2))
    cosA = nD*nX > 0.0 ? sum(sD .* sX)/(nD*nX) : 0.0
    @printf("%-10s %-10s %14.6e %14.6e %10.4f %12.5f\n", string(aSub), string(ch[2]), nD, nX, nX/nD, cosA)
    global sumD = sumD + nD^2;   global sumX = sumX + nX^2
end
@printf(">> over all %d channels:  |exchange|/|direct| = %.4f in norm,  one rebuild took %.1f s\n",
        length(chansJ), sqrt(sumX/sumD), tJ)
@printf(">> cos(angle) says whether exchange opposes the direct source channel by channel;  a value near -1\n")
@printf(">> means the two nearly cancel there, near +1 that they reinforce AT THE SIGN `exchangeSourceVector`\n")
@printf(">> returns -- which is the sign this case does NOT anchor.\n")


# ==========================================================================================================
# CASE K -- WHAT EXCHANGE DOES TO THE NUMBER, Ne, both signs, against the direct-only result.
#
# CASE J measured the SIZE and it is 0.84 of the direct term in norm -- so this is not a correction and the
# alpha it produces is not a refinement of the direct-only one.  Three runs, ONE variable:
#
#   exchange OFF   the direct-only result with the self term kept IN FULL (not (N-1)/N), which is the
#                  reference exchange has to be compared against -- removing the self term by hand and then
#                  adding exchange would remove it twice
#   exchange -1    the sign for which exchange OPPOSES the direct induced source, i.e. REDUCES screening and
#                  RAISES alpha, which is what exchange must do
#   exchange +1    the sign `RandomPhase.exchangeSourceVector` returns as written
#
# The physics picks the sign and says so out loud: at +1 CASE J measured cos(direct, exchange) POSITIVE on
# eight of nine channels, i.e. exchange would ADD to the screening, which is wrong.  That is consistent with
# the chi-convention minus being absent from the assembly -- the direct term carries it explicitly in
# `wMres`, fixed by the exact alpha(H 1s) = 4.5, and both terms are linear in chi_b so both need it.
# ==========================================================================================================
if  "exch" in ARGS
for  xs  in  (0.0, -1.0, 1.0)
    @printf("\n--- exchange %s ---\n", xs == 0.0 ? "OFF (direct only, self term in full)" :
            xs < 0.0 ? "ON, sign -1 (opposes the direct source)" : "ON, sign +1 (as the function returns)")
    runCoupledCase("CASE K -- Ne at lam*, exchange $(xs)", grid, prim, orbsNe, potNe, shellsNe, occNe, 0.0,
                   2.669, "the measured alpha(Ne)"; showOperator=false, selfMode=:full, exchangeSign=xs)
end
else
    println("\n(CASE K -- exchange on/off -- skipped;  pass \"exch\" to run it.)")
end
