# Exact non-local exchange in the continuum (continuumExchange = :nonlocal) against the local models, and its cost.
#
#   (1) Be-like C2+ 2s2 -> 2s2p at 20 and 50 eV: the 3P levels are reached by exchange, so the models spread by up to 1.9x
#       (tools/diag-ie-exchange.jl); :nonlocal is the one without an approximation to the exchange operator.
#   (2) Mg+ 3s -> 3p at 10 and 20 eV, maxKappa = 20, with the tail: static against non-local distorted waves, to be put
#       beside the relativistic distorted-wave code run with and without its own non-local exchange (work/rdw-mg).
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-nonlocalExchange.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

function run(grid, nm, ic, fc, energies, exch, kmax; tail=false)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = energies, calcRateCoefficient = false,
                                    maxKappa = kmax, printBefore = false, continuumExchange = exch, asymptoticTail = tail)
    wa  = Atomic.Computation(Atomic.Computation(), name="diag", grid=grid, nuclearModel=nm,
                             initialConfigs=[ic], finalConfigs=[fc], processSettings=set)
    st  = redirect_stdout(devnull) do;   @timed perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( st.value, st.time )
end

# (1)
g1 = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)
models = (:Slater, :FurnessMcCarthy, :none, :nonlocal)
res = Dict(m => run(g1, Nuclear.Model(6.), Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"), [20.0, 50.0], m, 10) for m in models)
println("\n=== (1) Be-like C2+ 2s2 -> 2s2p, maxKappa = 10: Omega, and wall time of the whole run ===")
@printf("  %6s %-8s %12s %12s %12s %12s\n", "E [eV]", "J_f (#)", "Slater", "FM", "none", "nonlocal")
for  i in eachindex(res[:Slater][1])
    l = res[:Slater][1][i]
    @printf("  %6.1f %-8s %s\n", Defaults.convertUnits("energy: from atomic to eV", l.initialElectronEnergy),
            string(l.finalLevel.J) * " (" * string(l.finalLevel.index) * ")",
            join([@sprintf("%12.5e", res[m][1][i].collisionStrength) for m in models], " "))
end
println("  wall time [s]:  ", join([@sprintf("%s %.1f", m, res[m][2]) for m in models], ",  "))

# (2)
g2 = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 25.0)
res2 = Dict(m => run(g2, Nuclear.Model(12.), Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), [10.0, 20.0], m, 20; tail=true)
            for m in (:none, :nonlocal))
println("\n=== (2) Mg+ 3s -> 3p, maxKappa = 20, tail on: sigma [a0^2] ===")
@printf("  %6s %-6s %12s %12s %10s\n", "E [eV]", "J_f", "static", "nonlocal", "nl/static")
for  i in eachindex(res2[:none][1])
    a = res2[:none][1][i];   b = res2[:nonlocal][1][i]
    @printf("  %6.1f %-6s %12.4f %12.4f %10.4f\n", Defaults.convertUnits("energy: from atomic to eV", a.initialElectronEnergy),
            string(a.finalLevel.J), a.crossSection, b.crossSection, b.crossSection/a.crossSection)
end
println("  wall time [s]:  static $(round(res2[:none][2], digits=1)),  nonlocal $(round(res2[:nonlocal][2], digits=1))")
