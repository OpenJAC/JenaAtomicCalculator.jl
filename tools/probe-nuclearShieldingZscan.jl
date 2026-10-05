# ITS INPUT IS UNRELIABLE, 04-Oct-2026.  This probe asks whether the self-consistent shift follows Z or the
# shells present, and answers from 14 ions computed on their recommended grids.  Those grids are 0.1-2.5 % out
# for ions with no outer d or f shell and a FACTOR OF 2.9 out for Hg(2+) -- which is one of the two ions the
# conclusion rests on.  So the scan must be repeated with hp halved until each UNCOUPLED gamma has settled before
# anything is read from it; see priority item 43 and `tools/probe-nuclearShieldingDensity.jl`.
# KEPT because the isoelectronic-pair design is right and is what refuted the first Z hypothesis.
#
# probe-nuclearShieldingZscan.jl -- IS THE SIGN REVERSAL OF THE SELF-CONSISTENT SHIFT A Z EFFECT OR A SHELL EFFECT?
#
# Measured on 04-Oct-2026, |gamma| FALLS under self-consistency for Y(3+) and In(3+) and RISES for Th(4+).  Three
# ions cannot tell which variable is responsible, and the two candidates make different predictions:
#
#   IF IT IS Z,     the change should vary smoothly along the whole ladder, and ALSO within an isoelectronic
#                   family where the configuration is held fixed and only the charge moves.
#   IF IT IS THE SHELLS PRESENT, it should jump when a new shell appears -- and the obvious suspect is the 4f,
#                   which Y(3+) and In(3+) do not have and Th(4+) does -- while staying flat inside a family.
#
# The set below is built to separate them: FOUR isoelectronic pairs or triples (same configuration, different Z)
# spread along a ladder that crosses the points where 4d, 4f and 5d first appear.  Bi(3+) is the decisive single
# point: it has the filled 4f and sits at Z = 83, so if the effect follows the shells it should rise like
# thorium, and if it follows Z it should sit between In(3+) and Th(4+).
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
const NS  = JenaAtomicCalculator.NuclearShielding

cases = [("Y^3+",  39.0, "[Kr]"),                      ("Zr^4+", 40.0, "[Kr]"),
         ("In^3+", 49.0, "[Kr] 4d^10"),                ("Sn^4+", 50.0, "[Kr] 4d^10"),
         ("Ba^2+", 56.0, "[Xe]"),                      ("La^3+", 57.0, "[Xe]"),
         ("Lu^3+", 71.0, "[Xe] 4f^14"),                ("Hf^4+", 72.0, "[Xe] 4f^14"),
         ("Hg^2+", 80.0, "[Xe] 4f^14 5d^10"),          ("Pb^4+", 82.0, "[Xe] 4f^14 5d^10"),
         ("Bi^3+", 83.0, "[Xe] 4f^14 5d^10 6s^2"),
         ("Ra^2+", 88.0, "[Rn]"),                      ("Th^4+", 90.0, "[Rn]"),
         ("U^6+",  92.0, "[Rn]")]

setsU = NS.Settings()
setsS = NS.Settings(NS.Settings(); selfConsistent=true, scfIterations=16, scfMixing=0.5)

println("\n  THE SELF-CONSISTENT SHIFT ACROSS Z, WITH ISOELECTRONIC FAMILIES HELD TOGETHER")
println("  " * "="^104)
@printf("  %-7s %4s %-24s %12s %16s %10s   %s\n",
        "ion", "Z", "configuration", "uncoupled", "self-consistent", "change", "4f?")
println("  " * "-"^104)
for  (name, Z, confs)  in  cases
    conf = Configuration(confs);   nm = Nuclear.Model(Z)
    grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    tmp  = tempname()
    out  = try
        mp = open(tmp,"w") do io;  redirect_stdout(io) do
                 perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                         configs=[conf], asfSettings=asf); output=true)["multiplet:"]  end  end
        gU = open(tmp,"w") do io;  redirect_stdout(io) do
                 NS.computeOutcomes(mp, nm, grid, setsU; output=true)[1].gammaE2  end  end
        gS = open(tmp,"w") do io;  redirect_stdout(io) do
                 NS.computeOutcomes(mp, nm, grid, setsS; output=true)[1].gammaE2  end  end
        (gU, gS)
    catch err
        println("  $name FAILED: ", first(split(sprint(showerror,err),"\n")));   nothing
    end
    rm(tmp, force=true)
    out === nothing && continue
    gU, gS = out
    @printf("  %-7s %4.0f %-24s %12.2f %16.2f %9.1f %%   %s\n", name, Z, confs, gU, gS,
            100*(abs(gS)-abs(gU))/abs(gU), occursin("4f", confs) ? "yes" : "no")
    flush(stdout)
end
println("  " * "-"^104)
println("  Read DOWN the column for a Z trend; read WITHIN a family (same configuration) for a pure charge")
println("  effect; and compare the 'no 4f' rows against the 'yes' ones for a shell effect.")
