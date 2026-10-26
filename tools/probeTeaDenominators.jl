#
# tools/probeTeaDenominators.jl
#
# Measures the TWO second-order denominators of the argon two-electron Auger case -- the one the BOUND
# intermediates carry and the one the CONTINUUM intermediates carry -- and asks whether they are the same
# quantity with opposite sign.  Written 26-Oct-2026 for priority items 93 and 98.  It is a diagnostic, not a
# test, and it is cheap: four BOUND multiplets and no continuum orbital at all.
#
# WHY THIS IS THE FIRST THING TO MEASURE.  The sum over intermediate states of the second-order amplitude runs
# over a COMPLETE set, which partitions by bound-electron count: 16 bound electrons (what JAC computes today)
# and 15 bound + 1 free (what Zitnik et al. compute, and what JAC has none of).  The two are disjoint and add
# COHERENTLY, so their relative SIGN decides whether adding the missing half raises our number or lowers it.
# The step-2 matrix element of the continuum branch splits exactly,
#
#     <f| V |n(eps)>  =  delta(eps - eps_f) * C  +  R(eps),
#
# because the continuum electron is a pure SPECTATOR of the core junction C and therefore contributes only its
# own overlap, which for per-energy normalization is a delta function.  Integrating the delta piece puts the
# intermediate electron at eps_f and leaves the denominator
#
#     E_i - E_g - eps_f  =  E_f - E_g       (using E_i = E_f + eps_f, i.e. energy conservation)
#
# while the bound branch carries  E_i - E_n.  Both differences are "move one 3p electron into the 2p hole",
# measured once in the 16-electron ion and once in the 15-electron one -- so they should come out with the SAME
# MAGNITUDE AND OPPOSITE SIGN, and the delta piece of the continuum branch would then interfere DESTRUCTIVELY
# with everything JAC computes today.  That is a prediction about a 2.25x excess, and this probe tests only its
# energetic half; the numerators are a separate question and need the amplitudes.
#
# WHAT IT ESTABLISHED, 26-Oct-2026.  THE PREDICTION HOLDS, and closer than it had any right to:
#
#     bound, 2p^-1 3p^-1  (one 2p hole filled)      E_i - E_n  =  +252.7, +255.2, +256.7 eV
#     bound, 3p^-2        (both 2p holes filled)    E_i - E_n  =  +493.6, +495.7 eV
#     continuum, 2p^-1 3p^-2 eps, delta piece       E_f - E_g  =  -234.0 ... -246.8 eV
#
# The dominant bound family and the WHOLE continuum family are the same physical quantity -- one 3p electron
# moves into a 2p hole -- measured once in the 16-electron ion and once in the 15-electron ion plus a free
# electron, and they agree in magnitude to about 6 % with opposite sign.  The eps_f column came out at 450.6 to
# 455.8 eV, matching the electron energies branch f prints, which is the internal check that the two charge
# states are on one energy scale.
#
# TWO THINGS THIS DOES NOT SHOW, and they are the reason it is only half the argument.  The numerators carry
# the two operators in OPPOSITE ORDER on the two paths -- V^(e-e) then H - E going through a bound
# intermediate, H - E then V^(e-e) going through a continuum one -- with different angular couplings, so the
# cancellation need not be as complete as the propagators suggest.  And the deep 3p^-2 bound family at +495 eV
# has NO delta-piece partner at all: its counterpart core would be the final state itself, where the
# denominator vanishes and the junction is diagonal rather than a transition.
#
# RE-RUN AND CERTIFIED the same morning, bit-identical, because a peer session had touched src/ six seconds
# before the first run exited (module-SelfConsistent-inc-optimizedlevel.jl, which this HSField path never
# reaches).  The two runs agree across two different states of that file, which is itself the control.
#
# Run:  julia --project=<JAC-root> tools/probeTeaDenominators.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

# [AR-TEA-DENOMINATORS]  ONE question: is  E_i - E_n(bound)  equal and opposite to  E_f - E_g(core)?
# Everything is the branch-f setup verbatim (examples/example-De.jl): same grid, same HSField, same Z.

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
asfSettings = AsfSettings(AsfSettings(), scField=Basics.HSField())
nm          = Nuclear.Model(18.)

# The four bound spaces.  (i) and (n) carry 16 electrons, (g) and (f) carry 15.
spaces = [ ("i  initial   2p^-2",      [Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")]),
           ("n  bound     2p^-1 3p^-1", [Configuration("1s^2 2s^2 2p^5 3s^2 3p^5"),
                                         Configuration("1s^2 2s^2 2p^6 3s^2 3p^4")]),
           ("g  core      2p^-1 (..)",  [Configuration("1s^2 2s^2 2p^5 3s^2 3p^4"),
                                         Configuration("1s^2 2s^2 2p^5 3s 3p^5"),
                                         Configuration("1s^2 2s^2 2p^5 3p^6")]),
           ("f  final     3p^-3",       [Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")]) ]

mults = Multiplet[]
for (label, configs) in spaces
    comp = Atomic.Computation(Atomic.Computation(), name=label, grid=grid, nuclearModel=nm,
                              configs=configs, asfSettings=asfSettings)
    push!(mults, perform(comp; output=true)[ResultKeys.Multiplet])
end
iMult, nMult, gMult, fMult = mults

toEv(x) = Defaults.convertUnits("energy: from atomic to eV", x)

iLevel = iMult.levels[1];    symi = LevelSymmetry(iLevel.J, iLevel.parity)
println("\n" * "="^112)
@printf("  INITIAL  2p^-2  level 1:  J^P = %s,  energy %14.6f a.u. = %12.3f eV\n", string(symi), iLevel.energy, toEv(iLevel.energy))
println("="^112)

# (1) The bound denominators actually used today: same total symmetry as the initial level, self-term excluded.
println("\n  (1) BOUND intermediates -- the denominator  E_i - E_n  that example-De.jl branch f divides by")
println("  " * "-"^108)
@printf("  %-5s %-10s %16s %16s\n", "n", "J^P", "E_n [eV]", "E_i - E_n [eV]")
for nLevel in nMult.levels
    symn = LevelSymmetry(nLevel.J, nLevel.parity)
    if  symn != symi                          continue    end
    if  abs(nLevel.energy - iLevel.energy) <= 1.0e-8   continue    end
    @printf("  %-5d %-10s %16.3f %16.3f\n", nLevel.index, string(symn), toEv(nLevel.energy), toEv(iLevel.energy - nLevel.energy))
end

# (2) The continuum denominators of the delta piece.  The junction is a SCALAR on the core, and the continuum
#     electron is a spectator with a fixed kappa, so the total symmetry matches only for J_g = J_f and equal
#     parity -- hence the pairing below is per final level.
println("\n  (2) CONTINUUM intermediates, delta piece -- the denominator  E_f - E_g  it divides by")
println("  " * "-"^108)
@printf("  %-5s %-10s %12s %10s | %-5s %-10s %12s %14s\n",
        "f", "J^P", "eps_f [eV]", "E_f [eV]", "g", "J^P", "E_g [eV]", "E_f - E_g [eV]")
for fLevel in fMult.levels
    symf  = LevelSymmetry(fLevel.J, fLevel.parity)
    epsf  = iLevel.energy - fLevel.energy
    first = true
    for gLevel in gMult.levels
        symg = LevelSymmetry(gLevel.J, gLevel.parity)
        if  symg != symf    continue    end
        if  first
            @printf("  %-5d %-10s %12.3f %10.3f | %-5d %-10s %12.3f %14.3f\n", fLevel.index, string(symf),
                    toEv(epsf), toEv(fLevel.energy), gLevel.index, string(symg), toEv(gLevel.energy),
                    toEv(fLevel.energy - gLevel.energy))
            first = false
        else
            @printf("  %-5s %-10s %12s %10s | %-5d %-10s %12.3f %14.3f\n", "", "", "", "",
                    gLevel.index, string(symg), toEv(gLevel.energy), toEv(fLevel.energy - gLevel.energy))
        end
    end
    if  first   @printf("  %-5d %-10s %12.3f %10.3f | %s\n", fLevel.index, string(symf), toEv(epsf),
                        toEv(fLevel.energy), "NO core level of this J^P -- the delta piece does not reach it")    end
end

println("\n  eps_f above must match the 450-456 eV electron energies branch f prints; that is the internal check.")
println("  THE TEST: the  E_i - E_n  column of (1) against the  E_f - E_g  column of (2), magnitude and sign.")
println("\n=== TEA DENOMINATOR PROBE FINISHED ===")
