# The Coulomb-Bethe continuation of the partial-wave sum (coulombBetheKappa) against brute-force references and against the
# geometric top-up, for the two dipole lines where the geometric top-up was measured to be a lower bound.
#
#   (1) H-like C5+ 1s -> 2p at 2939 eV; reference example-Dl branch a (maxKappa = 120): 9.638804e-02, 1.927089e-01.
#   (2) Na-like Mg II 3s -> 3p at 200 eV, rbox = 25; reference example-Dl branch d (rbox = 120, maxKappa = 300):
#       Omega_total = 112.83.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-coulombBethe.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

function omegas(grid, Z, ic, fc, E, kmax; topUp=false, cb=0, tail=true)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [E], calcRateCoefficient = false,
                                    maxKappa = kmax, printBefore = false, asymptoticTail = tail, topUp = topUp,
                                    coulombBetheKappa = cb, continuumExchange = :Slater)
    wa  = Atomic.Computation(Atomic.Computation(), name="cb", grid=grid, nuclearModel=Nuclear.Model(Z),
                             initialConfigs=[ic], finalConfigs=[fc], processSettings=set)
    st  = redirect_stdout(devnull) do;   @timed perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( [l.collisionStrength for l in st.value], st.time )
end

PART = get(ENV, "PART", "12")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))
g1  = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)
ref = [9.638804e-02, 1.927089e-01]
if  occursin("1", PART)
println("\n=== (1) C5+ 1s -> 2p, 2939 eV: Omega / Omega(maxKappa = 120) ===")
for  (lab, kw) in [("maxKappa 20, geometric top-up", (topUp=true,)), ("maxKappa 20, Coulomb-Bethe to 150", (cb=150,)),
                   ("maxKappa 20, CB to 150 + top-up", (cb=150, topUp=true))]
    om, t = omegas(g1, 6., Configuration("1s"), Configuration("2p"), 2939., 20; kw...)
    @printf("  %-36s  %.4f  %.4f    %6.1f s\n", lab, om[1]/ref[1], om[2]/ref[2], t)
end
end

if  occursin("2", PART)
g2 = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 25.0)
println("\n=== (2) Mg II 3s -> 3p, 200 eV, rbox = 25: Omega_total / 112.83 (rbox = 120, maxKappa = 300) ===")
for  (lab, kw) in [("maxKappa 40, geometric top-up", (topUp=true,)), ("maxKappa 40, Coulomb-Bethe to 300", (cb=300,)),
                   ("maxKappa 40, CB to 300 + top-up", (cb=300, topUp=true))]
    om, t = omegas(g2, 12., Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), 200., 40; kw...)
    @printf("  %-36s  %.4f    %6.1f s\n", lab, sum(om)/112.83, t)
end
end

# RESULT, 25-Sep-2026:
#   (1) C5+ 1s -> 2p, 2939 eV, against maxKappa = 120 computed with TODAY's code (9.216e-02 / 1.843e-01; the value recorded in
#       example-Dl of 04-Aug-2026 is 4.4 % higher and no longer reproduces -- tools/diag-ie-branchA-drift.jl):
#         maxKappa 20 + geometric top-up 0.997,  maxKappa 20 + Coulomb-Bethe to 150 0.998  (plain maxKappa 20: 0.878).
#       Coulomb-Bethe waves 21-40 reproduce the distorted waves 21-40 to 98.5 %; single dipole integrals agree to 2 % at
#       l = 30 and 0.06 % at l = 80.
#   (2) Mg II 3s -> 3p, 200 eV, rbox = 25, against rbox = 120 / maxKappa = 300 (112.83):
#         maxKappa 40 + geometric top-up 0.923 (39 s),  maxKappa 40 + Coulomb-Bethe to 300 0.993 (56 s).
