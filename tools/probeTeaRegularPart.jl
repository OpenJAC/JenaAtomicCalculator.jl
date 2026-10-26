#
# tools/probeTeaRegularPart.jl
#
# STEP 1 of building R(eps), the regular part of the continuum branch of the two-electron Auger amplitude
# (priority items 93 and 98).  Written 26-Oct-2026.  It is a diagnostic, not a test, and it answers ONE
# question before any machinery is built:
#
#     CAN JAC FORM  <f_core, eps_f kappa_f | V^(e-e) | g_core, eps kappa>,  i.e. a two-body matrix element
#     with a CONTINUUM ORBITAL ON BOTH SIDES at DIFFERENT energies -- and can the answer be trusted?
#
# WHY THIS IS THE WHOLE OF R(eps).  The second matrix element of the continuum branch splits exactly,
#
#     <f | V | g_core, eps>  =  delta(eps - eps_f) <f_core| H - E |g_core>  +  R(eps),
#
# and the delta piece -- the continuum electron a pure spectator of a junction acting on the core alone -- is
# already built and committed (AutoIonization.computeTeaContinuumAmplitude, 0288c867).  R(eps) is everything
# else: the part in which the free electron is NOT a spectator but PARTICIPATES.  Reading g_core -> f_core,
# a 3p electron fills the 2p hole and the free electron carries off the ~240 eV released, going eps -> eps_f.
# That is an inelastic scattering of the outgoing Auger electron off the ion it just left, and it is the
# mechanism the published comparison calculation consists of.
#
# ITS SHAPE IS FIXED BY ARITHMETIC, NOT BY CHOICE.  f = (f_core, eps_f) and g = (g_core, eps) each carry 15
# bound electrons plus one free, so V connects them, and they differ in TWO orbitals: the core's 3p -> 2p and
# the free electron's eps -> eps_f.  Two orbital differences are exactly where H - E reduces to V^(e-e)
# identically (item 98), so the operator here is the plain two-body interaction and no junction reasoning is
# needed.  Energy bookkeeping confirms it: eps_f = eps + (E_g - E_f), and with E_g - E_f ~ +240 eV and the
# pole at eps* = E_i - E_g ~ 210 eV that gives eps_f ~ 450 eV, which is the electron energy branch f prints.
#
# THE TRAP THIS PROBE EXISTS TO CATCH.  Two continuum subshells that differ only in their LABEL count as
# different orbitals, and two that share a label are collapsed by Basics.merge into one -- silently.  JAC
# happens to label the observed electron Subshell(101,kappa) in AutoIonization and the intermediate one
# Subshell(103,kappa) in PhotoDoubleIonization, so the labels are already distinct; but relying on that is
# luck, and a collapsed label would not raise -- it would return a plausible wrong number.
#
# SO THE CALIBRATION IS A MUTATION, NOT A PLAUSIBILITY READING.  If the intermediate radial orbital genuinely
# enters the integral, R MUST MOVE when eps is changed.  A machinery that collapsed the two labels, or looked
# up the wrong radial function, would return a value INDEPENDENT of eps -- which looks entirely reasonable on
# a single line of output.  Three checks are therefore made:
#
#     (A) EXISTS      R is finite and non-zero for an allowed partial-wave pair.
#     (B) HERMITIAN   <f|V|n> = <n|V|f>.  The operator is real symmetric and JAC's continuum orbitals are
#                     real, so these must agree; they will NOT agree if the subshell merging or the orbital
#                     lookup is wrong on one side.
#     (C) MOVES       R varies with eps.  This is the mutation test, and it is the one that cannot be passed
#                     by accident.
#
# The phase factor is stripped deliberately: AutoIonization.amplitude multiplies by i^l exp(-i*phase) for the
# kappa it is handed, and passing kappa = -1 with phase = 0 makes that factor i^0 exp(0) = 1 exactly.  The
# observed electron's factor belongs to the TOTAL amplitude, applied once, and must NOT be applied to each
# matrix element separately -- the intermediate electron is an internal line of a sum over a complete set of
# REAL basis states and carries no boundary-condition phase of its own.
#
# Run:  julia --project=<JAC-root> tools/probeTeaRegularPart.jl
#
# WHAT IT ESTABLISHED, 26-Oct-2026.  ALL THREE CHECKS PASS, on 246 matrix elements over 42 cores, 2 waves and
# 3 energies each:
#
#     (A) EXISTS      138 of 246 non-zero.
#     (B) HERMITIAN   <f|V|n> = 5.045200e-03 against <n|V|f> = 5.045200e-03, relative difference 3.8e-10
#                     (strongest pair: core 14, kappa = -4, eps = 317.07 eV).
#     (C) MOVES       46 of 46 live groups depend on eps, 23 of them CHANGING SIGN across the three energies;
#                     relative spreads 0.13 to 0.98.  36 groups are identically zero.
#
# So the two-sided continuum matrix element can be formed by AutoIonization.amplitude exactly as it stands, and
# the intermediate radial orbital demonstrably enters the integral.
#
# AND THE 36 ZERO GROUPS ARE THE USEFUL RESULT, because they are a SELECTION RULE the machinery derived rather
# than noise.  Once the free electron's own eps -> eps_f is counted it consumes ONE of V's two slots, so R(eps)
# reaches only cores ONE orbital away from f_core = 2p^6 3s^2 3p^3 -- while the delta piece, built on H - E,
# reached all three.  Of the three Ar(3+) cores only 2p^5 3s^2 3p^4 qualifies; 2p^5 3s 3p^5 and 2p^5 3p^6 are
# THREE changes away and vanish.  Level 39 (J=1/2 +) is pure even parity -- only 2p^5 3s 3p^5 carries even
# parity -- hence exactly zero, which is why the zeros track parity.  THE TWO BRANCHES THEREFORE SAMPLE
# DIFFERENT CORE SPACES and must be converged separately; the delta piece's reasoning does not carry over.
#
# TWO CONSEQUENCES FOR THE COST, which is what makes the full run affordable.  A naive loop is 42 cores x ~4
# kappas x (mesh + 1 pole) intermediate levels, each needing a continuum orbital and two matrix elements, PER
# observed partial wave and PER final level.  (1) A (core, kappa) block that vanishes at one energy vanishes at
# all of them -- verified here -- so the block is skipped after a single test rather than evaluated at every
# node.  (2) <n|V|i> does not depend on kappa_f, only <f|V|n> does, so the Auger factor is computed once per
# intermediate level and reused across the observed partial waves.
#
# AND R IS NOT SMALL: its largest value here, 5.0e-3, is the size of the delta piece's own junctions
# (0.005 to 0.027), so this is the larger half of the continuum branch rather than a correction to it.
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

# [AR-TEA-RPART-STEP1]  ONE question: does the two-sided continuum matrix element exist and move with eps?
# The setup is the branch-f setup verbatim (examples/example-De.jl), including both continuum defaults.

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
asfSettings = AsfSettings(AsfSettings(), scField=Basics.HSField())
nm          = Nuclear.Model(18.)
nucPot      = Nuclear.nuclearPotential(nm, grid)
cSettings   = Continuum.Settings(false, grid.NoPoints-50)
kind        = CoulombInteraction()

toEv(x) = Defaults.convertUnits("energy: from atomic to eV", x)

spaces = [ ("i  initial  2p^-2", [Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")]),
           ("g  cores    2p^-1", [Configuration("1s^2 2s^2 2p^5 3s^2 3p^4"),
                                  Configuration("1s^2 2s^2 2p^5 3s 3p^5"),
                                  Configuration("1s^2 2s^2 2p^5 3p^6")]),
           ("f  final    3p^-3", [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]) ]
mults = Multiplet[]
for (label, configs) in spaces
    comp = Atomic.Computation(Atomic.Computation(), name=label, grid=grid, nuclearModel=nm,
                              configs=configs, asfSettings=asfSettings)
    push!(mults, perform(comp; output=true)[ResultKeys.Multiplet])
end
iMult, gMult, fMult = mults

iLevel = iMult.levels[1];    symi = LevelSymmetry(iLevel.J, iLevel.parity)
fLevel = fMult.levels[1];    symf = LevelSymmetry(fLevel.J, fLevel.parity)
epsf   = iLevel.energy - fLevel.energy

println("\n" * "="^120)
@printf("  AR TWO-ELECTRON AUGER -- R(eps) STEP 1:  can  <f, eps_f kappa_f | V | g, eps kappa>  be formed?\n")
println("="^120)
@printf("  initial  2p^-2    J^P = %-6s  E_i = %12.4f a.u.\n", string(symi), iLevel.energy)
@printf("  final    3p^-3    J^P = %-6s  E_f = %12.4f a.u.   eps_f = %10.4f a.u. = %9.2f eV\n",
        string(symf), fLevel.energy, epsf, toEv(epsf))

# The observed electron: the partial waves that couple the final level to the initial total symmetry.
kappasF = [k for k in AngularMomentum.allowedKappaSymmetries(symf, symi) if abs(k) <= 4]
@printf("  observed partial waves kappa_f (|kappa| <= 4): %s\n", string(kappasF))

# The final SCATTERING state, built exactly as AutoIonization.computeAmplitudesProperties builds it.
redFLevel = Basics.generateLevelWithSymmetryReducedBasis(fLevel, fLevel.basis.subshells)
kappaF    = kappasF[1]
cOrbital, cPhase = Continuum.generateOrbitalForLevel(epsf, Subshell(101, kappaF), redFLevel, nm, grid, cSettings;
                                                     nuclearPot=nucPot)
fScatt = Basics.generateLevelWithExtraElectron(cOrbital, symi, redFLevel)
@printf("  final scattering state: kappa_f = %d, phase = %.6f, %d CSFs, %d subshells\n",
        kappaF, cPhase, length(fScatt.basis.csfs), length(fScatt.basis.subshells))

# The RAW matrix element: kappa = -1 and phase = 0 make AutoIonization.amplitude's own factor
# i^l exp(-i phase) = i^0 exp(0) = 1 exactly, so what comes back is the bare reduced matrix element.
rawME(bra, ket) = AutoIonization.amplitude(kind, -1, 0.0, bra, ket, grid; printout=false)

println("\n  " * "-"^116)
@printf("  R(eps) FOR EVERY CORE AND EVERY ALLOWED INTERMEDIATE WAVE.  pole at eps* = E_i - E_g.\n")
println("  " * "-"^116)
@printf("  %-4s %-8s %10s %10s | %-7s %12s %14s %14s\n",
        "g", "J^P_g", "E_g[a.u.]", "eps*[eV]", "kappa", "eps [eV]", "Re R(eps)", "|R(eps)|")
println("  " * "-"^116)

records = []
for (ig, gLevel) in enumerate(gMult.levels)
    symg     = LevelSymmetry(gLevel.J, gLevel.parity)
    epsStar  = iLevel.energy - gLevel.energy
    epsStar <= 0.  &&  continue
    redG     = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
    kappasN  = [k for k in AngularMomentum.allowedKappaSymmetries(symg, symi) if abs(k) <= 4]
    isempty(kappasN)  &&  continue
    # Three energies: well below the pole, AT the pole, and above it.  The mutation test lives in this column.
    for kappaN in kappasN[1:min(2,length(kappasN))]
        for frac in [0.5, 1.0, 1.5]
            en = frac * epsStar
            nOrb, nPh = Continuum.generateOrbitalForLevel(en, Subshell(103, kappaN), redG, nm, grid, cSettings;
                                                          nuclearPot=nucPot)
            nScatt = Basics.generateLevelWithExtraElectron(nOrb, symi, redG)
            R      = rawME(fScatt, nScatt)
            push!(records, (ig, kappaN, en, R))
            @printf("  %-4d %-8s %10.4f %10.2f | %-7d %12.2f %14.6e %14.6e\n",
                    gLevel.index, string(symg), gLevel.energy, toEv(epsStar), kappaN, toEv(en), real(R), abs(R))
        end
    end
end

println("\n  " * "="^116)
println("  THE THREE CHECKS")
println("  " * "="^116)

# (A) EXISTS
nz = count(r -> abs(r[4]) > 0., records)
@printf("  (A) EXISTS     %d of %d matrix elements are non-zero.  %s\n", nz, length(records),
        nz > 0 ? "PASS" : "FAIL -- nothing was formed at all")

# (C) MOVES -- the mutation test.  Group by (core, kappa) and ask whether R depends on eps.
#     A group that is IDENTICALLY ZERO is NOT a suspect constant but a selection-rule zero, and is counted
#     separately: reaching an EVEN-parity core from f_core = 2p^6 3s^2 3p^3 costs three orbital changes once
#     the continuum electron's own eps -> eps_f is counted, and V is a TWO-body operator.  Lumping those in
#     with the live groups would report a failure where the machinery is right.
println("\n  (C) MOVES      R must DEPEND on eps.  A collapsed subshell label would give a constant column.")
let nLive = 0, nMoved = 0, nZero = 0
    for key in unique([(r[1], r[2]) for r in records])
        grp = [r for r in records if (r[1], r[2]) == key]
        length(grp) < 2  &&  continue
        vals = [abs(r[4]) for r in grp]
        if  maximum(vals) == 0.
            nZero = nZero + 1
            continue
        end
        nLive  = nLive + 1
        spread = (maximum(vals) - minimum(vals)) / maximum(vals)
        moved  = spread > 1.0e-6
        moved  &&  (nMoved = nMoved + 1)
        # A SIGN CHANGE across the three energies is stronger evidence still: an oscillatory radial overlap
        # cannot be faked by a mislookup that returns one fixed orbital.
        sgn    = length(unique(sign.([real(r[4]) for r in grp]))) > 1 ? "  SIGN CHANGES" : ""
        @printf("      core %2d, kappa %+d :  |R| = %s   rel. spread %8.4f   %s%s\n",
                key[1], key[2], join([@sprintf("%.4e", v) for v in vals], "  "), spread,
                moved ? "MOVES" : "CONSTANT -- SUSPECT", sgn)
    end
    @printf("\n      %d live groups, %d of them move;  %d groups identically zero (selection rule, see above)\n",
            nLive, nMoved, nZero)
    @printf("      verdict: %s\n", (nLive > 0 && nMoved == nLive) ?
            "PASS -- every live group depends on eps, so the intermediate radial orbital DOES enter" :
            "FAIL -- a live group is independent of eps")
end

# (B) HERMITIAN
println("\n  (B) HERMITIAN  <f|V|n> against <n|V|f>, same pair, same orbitals.")
if !isempty(records)
    ig, kappaN, en, Rfn = records[argmax([abs(r[4]) for r in records])]
    gLevel = gMult.levels[ig];    redG = Basics.generateLevelWithSymmetryReducedBasis(gLevel, gLevel.basis.subshells)
    nOrb, nPh = Continuum.generateOrbitalForLevel(en, Subshell(103, kappaN), redG, nm, grid, cSettings;
                                                  nuclearPot=nucPot)
    nScatt = Basics.generateLevelWithExtraElectron(nOrb, symi, redG)
    Rnf    = rawME(nScatt, fScatt)
    rel    = abs(Rfn) == 0. ? 0. : abs(Rfn - Rnf) / abs(Rfn)
    @printf("      strongest pair: core %d, kappa %+d, eps = %.2f eV\n", ig, kappaN, toEv(en))
    @printf("      <f|V|n> = %14.6e      <n|V|f> = %14.6e      relative difference = %.3e   %s\n",
            real(Rfn), real(Rnf), rel, rel < 1.0e-8 ? "PASS" : "FAIL -- the two sides are NOT equivalent")
end

println("\n=== TEA R(eps) STEP-1 PROBE FINISHED ===")
