# Controlled comparison for ImpactExcitation's caches, run BEFORE any infrastructure is replaced.
#
#   mode :none   ... cacheAngularCoefficients = false, cacheRadialIntegrals = false  (the pre-cache code path)
#   mode :cached ... both caches on (per-run AngularCoefficientCache, per-line InteractionStrength.XLCache)
#
# Every channel amplitude and sigma are compared against :none, and wall time and allocations are measured on a second
# (compiled) run. The DCS of the :cached run is integrated over the sphere and compared with sigma.
#
# Run from <JAC-root>:   julia --project=. -t 1 tools/diag-ie-radialcache.jl     (and again with -t 2)

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation, InteractionStrength

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))

const GRID   = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)
const ENERGY = 50.0      # eV above the 2s -> 2p threshold region of C3+ (about 8 eV); comfortably above it
const KMAX   = 10

# RESULT OF THE FIRST RUN (24-Sep-2026, pre-refactor code, 1 thread): the three modes :none / :local / :xl gave
# bit-identical amplitudes and sigma = 1.492167e+00, 2.985159e+00 a.u.; times 27.4 / 25.3 / 24.8 s. On that evidence the
# local Dict was replaced by InteractionStrength.XLCache. The :xl override that made the comparison is therefore gone, and
# the harness now compares :none against :cached (both caches on) in the refactored code, and against those numbers.
const PRE_REFACTOR_SIGMA = [1.492167e+00, 2.985159e+00]

function runMode(mode::Symbol; dcs::Bool=false)
    cache = mode != :none
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [ENERGY], calcRateCoefficient = false,
                                    maxKappa = KMAX, printBefore = false, operator = CoulombInteraction(),
                                    calcDCS = dcs, angles = collect(0.0:0.5:180.0),
                                    cacheAngularCoefficients = cache, cacheRadialIntegrals = cache)
    wa = Atomic.Computation(Atomic.Computation(), name="diag: Li-like C3+ 2s -> 2p", grid=GRID,
                            nuclearModel = Nuclear.Model(6.),
                            initialConfigs = [Configuration("1s^2 2s")], finalConfigs = [Configuration("1s^2 2p")],
                            processSettings = set)
    local res
    out = redirect_stdout(devnull) do
        stats = @timed (res = perform(wa, output=true))
        stats
    end
    return res[ImpactExcitation.Settings], out.time, out.bytes
end

function amplitudes(lines)
    [ [ch.amplitude for ch in ln.channels] for ln in lines ]
end

# First pass compiles every path; second pass is timed.
if  get(ENV, "DCS_ONLY", "") == ""
for m in (:none, :cached)   runMode(m)   end
results = Dict{Symbol,Any}()
for m in (:none, :cached)   results[m] = runMode(m)   end

ref = amplitudes(results[:none][1])
println("\n=== ImpactExcitation cache comparison: C3+ 1s2 2s -> 1s2 2p,  E = $ENERGY eV,  maxKappa = $KMAX ===")
@printf("%-7s  %12s  %12s  %14s  %14s  %22s\n", "mode", "time [s]", "alloc [MB]", "max|dA|/max|A|", "max|dsig|/sig", "sigma per line [a.u.]")
for m in (:none, :cached)
    lines, t, b = results[m]
    amps = amplitudes(lines)
    dA   = maximum( maximum(abs.(amps[i] .- ref[i])) / maximum(abs.(ref[i])) for i in eachindex(ref) )
    dS   = maximum( abs(lines[i].crossSection - results[:none][1][i].crossSection) / results[:none][1][i].crossSection
                    for i in eachindex(lines) )
    @printf("%-7s  %12.2f  %12.1f  %14.3e  %14.3e  %s\n", m, t, b/2^20, dA, dS,
            join([@sprintf("%.6e", ln.crossSection) for ln in lines], "  "))
end
println("threads = $(Threads.nthreads());  pre-refactor sigma = $(PRE_REFACTOR_SIGMA)")

end

# DCS consistency: integral over the sphere against sigma, for the :cached path.
lines, _, _ = runMode(:cached, dcs=true)
println("\n=== DCS consistency: 2 pi * int dcs sin(theta) dtheta  against  sigma ===")
for ln in lines
    th  = deg2rad.(ln.angles);   f = ln.dcs .* sin.(th)
    int = 2pi * sum( (th[i+1]-th[i]) * (f[i]+f[i+1]) / 2 for i in 1:length(th)-1 )
    @printf("  J_i = %s -> J_f = %s :  integral = %.6e   sigma = %.6e   ratio = %.6f\n",
            string(ln.initialLevel.J), string(ln.finalLevel.J), int, ln.crossSection, int/ln.crossSection)
end
