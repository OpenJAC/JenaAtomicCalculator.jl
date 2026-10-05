# ANSWERED 04-Oct-2026, AND THE ANSWER IS THAT THERE IS NO ANOMALY.  The -21.49 this probe sets out to explain
# is a GRID artefact: Hg(2+) [Xe] 4f^14 5d^10 needs some 27 B-splines per a.u. before gamma_inf settles, and its
# recommended grid supplies 6.4.  Refined, gamma_inf(Hg2+) = -62.00, which sits where the neighbouring ions put
# it and needs no special explanation.  The outer d and f shells are what make this ion demand the finest mesh of
# the four measured -- so the configuration WAS the discriminator, just not in the way asked here.
# KEPT as the record of how the question was put, and because its channel decomposition is still the right tool.
#
# probe-nuclearShieldingHgAnomaly.jl -- WHY IS THE UNCOUPLED gamma OF Hg(2+) SO SMALL?
#
# The 14-ion scan of 04-Oct-2026 gives, UNCOUPLED (no screening, the module's default):
#       Ba(2+)  Z=56  [Xe]                      -101.83
#       Lu(3+)  Z=71  [Xe] 4f^14                 -84.99
#       Hg(2+)  Z=80  [Xe] 4f^14 5d^10           -21.49     <-- four to five times smaller
#       Pb(4+)  Z=82  [Xe] 4f^14 5d^10           -46.32     <-- same configuration, 2.2x larger
#       Bi(3+)  Z=83  [Xe] 4f^14 5d^10 6s^2      -53.37
# gamma should vary smoothly along an isoelectronic sequence, so Hg(2+) against Pb(4+) is already wrong before
# any screening enters.  Two things can do that, and they are distinguishable:
#   (a) A CHANNEL IS BEING LOST.  gamma is a sum over (subshell, kappa') channels; if one large channel is
#       dropped -- most likely by the Pauli-blocking test, which excludes intermediate states whose energy
#       matches an occupied one -- the total collapses.  The per-channel breakdown shows it immediately.
#   (b) THE GRID.  A neutral-ish heavy ion is diffuse; `recommendedGrid` sizes the box for the BOUND orbitals,
#       and Rule 12 of CLAUDE.md is explicit that a box mismatched to the orbitals produces a wrong NUMBER
#       rather than an obvious failure.  Hg(2+) is the least charged of the family (2+ against 4+ for Pb), so
#       it is the MOST diffuse and the most exposed to this.  Varying the box answers it.
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
const NS  = JenaAtomicCalculator.NuclearShielding

function run(name, Z, confs; rbox=0.0)
    conf = Configuration(confs);   nm = Nuclear.Model(Z)
    grid = rbox > 0 ? Basics.recommendedGrid([conf], nm; rbox=rbox, printout=false) :
                      Basics.recommendedGrid([conf], nm; printout=false)
    setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    tmp  = tempname()
    out  = open(tmp,"w") do io;  redirect_stdout(io) do
               mp = perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                            configs=[conf], asfSettings=asf); output=true)["multiplet:"]
               NS.computeOutcomes(mp, nm, grid, NS.Settings(); output=true)  end  end
    rm(tmp, force=true)
    return( out[1], grid )
end

println("\n  IS THE BOX DEPENDENCE CONFINED TO THE HEAVY IONS, OR IS IT EVERYWHERE?")
println("  " * "="^100)
println("  The heavy ions were found on 04-Oct to grow without limit with the box.  If the VALIDATED light")
println("  ions do the same, then the agreement with Sternheimer was obtained at one particular box and the")
println("  quantity is not converged for anybody -- which is a different and much larger statement.")
println("  " * "-"^100)
@printf("  %-7s %4s %-22s %9s %12s %14s\n", "ion", "Z", "configuration", "r_box", "gamma", "vs default")
for  (name, Z, confs)  in  (("Y^3+", 39.0, "[Kr]"), ("Ba^2+", 56.0, "[Xe]"), ("Th^4+", 90.0, "[Rn]"))
    ref = NaN
    for  rb  in  (0.0, 30.0, 60.0, 90.0)
        oc, grid = run(name, Z, confs; rbox=rb)
        isnan(ref)  &&  (ref = oc.gammaE2)
        @printf("  %-7s %4s %-22s %9.1f %12.3f %13.1f %%\n", rb == 0.0 ? name : "", rb == 0.0 ? string(Int(Z)) : "",
                rb == 0.0 ? confs : "", grid.tL[end], oc.gammaE2, 100*(oc.gammaE2-ref)/abs(ref))
        flush(stdout)
    end
    println()
end
println("  " * "-"^100)
println("  Sternheimer (raw, for comparison):  Y(3+) -34.1,  Ba(2+) --,  Th(4+) -- .")
