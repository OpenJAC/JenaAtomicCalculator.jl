# probe-nuclearShielding.jl -- THE ACCEPTANCE TEST FOR module-NuclearShielding.jl
#
# The quadrupole antishielding factor has a validation set that exists independently of us: Sternheimer,
# Phys. Rev. 159, 266 (1967), Table I, computed by a different method (direct solution of the inhomogeneous
# equation) on different orbitals (Herman-Skillman, NEUTRAL atom) with no experiment in it.  Three of his four
# ions are closed-shell and therefore in scope here; Am(3+) has an open 5f shell and this module refuses it.
#
# TWO THINGS MUST BE KEPT IN MIND WHEN READING THE COMPARISON, or agreement will be over- or under-claimed.
#   (a) His PUBLISHED numbers carry a hand-applied factor 1/1.1, which corrects for his use of neutral-atom
#       wave functions.  We use proper ION orbitals, so the honest comparison is against his RAW values,
#       i.e. the published ones times 1.1.
#   (b) His l -> l+-2 channels ("angular modes") were ESTIMATED from a Thomas-Fermi model and are 5-10 % of
#       his total.  We compute them properly, so a difference of that order is expected and is ours to keep.
# Together these put the meaningful agreement threshold at roughly 15 %, not at 1 %.
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
setDefaults("unit: energy", "eV")

# ion, Z, configuration, Sternheimer's PUBLISHED gamma (after his 1/1.1)
cases = [("Y^3+",  39.0, "[Kr]",        -31.0),
         ("In^3+", 49.0, "[Kr] 4d^10",  -23.0),
         ("Bi^3+", 83.0, "[Xe] 4f^14 5d^10 6s^2", -39.0),
         ("Th^4+", 90.0, "[Rn]",         NaN )]

println("="^104)
println("  ACCEPTANCE TEST:  relativistic quadrupole antishielding, uncoupled, uniform external field")
println("="^104)
@printf("  %-8s %5s %-26s %12s %12s %12s %9s\n",
        "ion", "Z", "configuration", "this work", "Sternh.pub", "Sternh.raw", "vs raw")
for (name, Z, confs, gPub) in cases
    conf = Configuration(confs)
    nm   = Nuclear.Model(Z)
    grid = Basics.recommendedGrid([conf], nm; printout=false)
    setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    tmp  = tempname()
    out  = try
        open(tmp,"w") do io;  redirect_stdout(io) do
            mp = perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                         configs=[conf], asfSettings=asf); output=true)["multiplet:"]
            NuclearShielding.computeOutcomes(mp, nm, grid, NuclearShielding.Settings(); output=true)
        end;  end
    catch err
        println("  $name FAILED: ", first(split(sprint(showerror,err),"\n")));  nothing
    end
    rm(tmp, force=true)
    out === nothing && continue
    g = out[1].gammaE2
    if  isnan(gPub)
        @printf("  %-8s %5.0f %-26s %12.2f %12s %12s %9s\n", name, Z, confs, g, "--", "--", "--")
    else
        raw = gPub * 1.1
        @printf("  %-8s %5.0f %-26s %12.2f %12.1f %12.1f %8.2f\n", name, Z, confs, g, gPub, raw, g/raw)
    end
    flush(stdout)
end
println("")
