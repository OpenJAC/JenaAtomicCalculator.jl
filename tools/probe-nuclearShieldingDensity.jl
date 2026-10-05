#
# tools/probe-nuclearShieldingDensity.jl
#
# THE DENSITY REFINEMENT TEST FOR A SHIELDING FACTOR -- the only check of the three tried that caught every bad
# case measured on 04-Oct-2026.  It keeps the box that Basics.recommendedGrid chooses and HALVES hp repeatedly,
# so the B-spline basis grows while the box does not, and reports gamma_inf at each density.  Two successive
# densities that agree mean the number is converged in the mesh -- and, because gamma_inf turns out not to
# depend on the box at fixed density, converged in the box as well, so no separate box study is needed.
#
# WHY A FIXED THRESHOLD WILL NOT DO, which is what this tool exists to replace.  The density an ion needs is a
# property of the ion: Th(4+) [Rn] is converged at 4.0 splines per a.u., while Hg(2+) [Xe] 4f^14 5d^10 at 4.70
# per a.u. in a 30 a.u. box returns -2.6e+07 -- seven orders of magnitude wrong, with both of the module's own
# guards silent.  Hg(2+) needs about 11 per a.u.  A threshold set high enough for Hg(2+) would reject Th(4+)'s
# perfectly sound default grid, so no single number serves and the refinement must actually be run.
#
# Usage:   julia --project=. tools/probe-nuclearShieldingDensity.jl
#
using JenaAtomicCalculator, Printf
const NS = JenaAtomicCalculator.NuclearShielding


"""
`gammaOnGrid(name::String, Z::Float64, confString::String, rbox::Float64, hp::Union{Nothing,Float64})`  
    ... computes gamma_inf for the given closed-shell ion on a grid of the given box and mesh parameter, with all
        of JAC's own printout suppressed; a tuple (rbox::Float64, nsL::Int64, density::Float64, gamma::Float64,
        a1::Float64) is returned, where density = nsL/rbox and a1 is the cancellation ratio of the channel sum.
"""
function gammaOnGrid(name::String, Z::Float64, confString::String, rbox::Float64, hp::Union{Nothing,Float64})
    conf = Configuration(confString);    nm = Nuclear.Model(Z)
    grid = isnothing(hp)  ?  Basics.recommendedGrid([conf], nm; rbox=rbox, printout=false)  :
                             Basics.recommendedGrid([conf], nm; rbox=rbox, hp=hp, printout=false)
    setDefaults("standard grid", grid)
    asfSettings = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    tmpFile = tempname()
    outcome = open(tmpFile, "w") do io
        redirect_stdout(io) do
            multiplet = perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                                configs=[conf], asfSettings=asfSettings); output=true)["multiplet:"]
            NS.computeOutcomes(multiplet, nm, grid, NS.Settings(); output=true)[1]
        end
    end
    rm(tmpFile, force=true)
    a1 = maximum(abs(c.value) for c in outcome.contributions) / abs(outcome.gammaE2)

    return( (grid.tL[end], grid.nsL, grid.nsL/grid.tL[end], outcome.gammaE2, a1) )
end


"""
`refine(name::String, Z::Float64, confString::String; nRefine::Int64=4, tolerance::Float64=0.005)`  
    ... runs the refinement test for one ion: gamma_inf is computed on the recommended grid and then with hp
        halved nRefine times, and the sequence is printed with the relative change at each step.  The test is
        passed at the first step whose relative change is below tolerance; a converged gamma::Float64 is
        returned, or NaN if the sequence has not settled within nRefine refinements.
"""
function refine(name::String, Z::Float64, confString::String; nRefine::Int64=4, tolerance::Float64=0.005)
    conf = Configuration(confString);    nm = Nuclear.Model(Z)
    grid0 = Basics.recommendedGrid([conf], nm; printout=false)
    rbox  = grid0.tL[end];    hp = grid0.hp
    @printf("\n  %s   (Z = %.0f, %s):  box kept at %.2f a.u., hp halved from %.4f\n", name, Z, confString, rbox, hp)
    println("  " * "-"^84)
    @printf("  %10s %6s %10s %14s %12s %8s\n", "hp", "nsL", "per a.u.", "gamma_inf", "rel. change", "A1")
    gammaPrev = NaN;    gammaConverged = NaN
    for  i = 0:nRefine
        hpi = hp / 2.0^i
        rb, nsL, density, gamma, a1 = gammaOnGrid(name, Z, confString, rbox, hpi)
        change = isnan(gammaPrev) ? NaN : abs(gamma - gammaPrev) / abs(gamma)
        @printf("  %10.5f %6d %10.2f %14.4f %12s %8.2f\n", hpi, nsL, density, gamma,
                isnan(change) ? "--" : @sprintf("%.2f %%", 100change), a1)
        flush(stdout)
        if  !isnan(change)  &&  change < tolerance  &&  isnan(gammaConverged)     gammaConverged = gamma    end
        gammaPrev = gamma
    end
    if  isnan(gammaConverged)
        println("  NOT CONVERGED within $nRefine refinements -- gamma_inf must not be quoted for this ion yet.")
    else
        @printf("  CONVERGED:  gamma_inf = %.2f  (to better than %.1f %%)\n", gammaConverged, 100tolerance)
    end

    return( gammaConverged )
end


println("\n  THE DENSITY REFINEMENT TEST  (box fixed, B-spline basis doubled at each step)")
println("  " * "="^84)
refine("Y^3+",  39.0, "[Kr]")
refine("Ba^2+", 56.0, "[Xe]")
refine("Th^4+", 90.0, "[Rn]")
refine("Hg^2+", 80.0, "[Xe] 4f^14 5d^10")
println("\n  " * "="^84)
