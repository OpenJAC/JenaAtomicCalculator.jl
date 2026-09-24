# How much of the partial-wave truncation does the geometric top-up recover, for a DIPOLE transition?
#
# H-like C5+ 1s -> 2p at 2939 eV, the grid of example-Dl.jl branch a. Reference: that branch's REPORT (04-Aug-2026,
# maxKappa = 120, dated Last successful 13-Aug-2026): Omega(2p_1/2) = 9.638804e-02, Omega(2p_3/2) = 1.927089e-01.
# A dipole transition is the hard case -- its partial strengths fall off slowly (Bethe) -- so this bounds from below
# what the top-up does for any line.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-topUp.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))

const GRID = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)
const REF  = [9.638804e-02, 1.927089e-01]

function omegas(kmax, topUp)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [2939.], calcRateCoefficient = false,
                                    maxKappa = kmax, printBefore = false, topUp = topUp, asymptoticTail = true)
    wa  = Atomic.Computation(Atomic.Computation(), name="C5+", grid=GRID, nuclearModel=Nuclear.Model(6.),
                             initialConfigs=[Configuration("1s")], finalConfigs=[Configuration("2p")], processSettings=set)
    lines = redirect_stdout(devnull) do;   perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( [l.collisionStrength for l in lines] )
end

println("\n=== C5+ 1s -> 2p at 2939 eV: Omega / Omega(maxKappa = 120, example-Dl) ===")
@printf("  %8s   %-24s   %-24s\n", "maxKappa", "plain  (2p_1/2, 2p_3/2)", "with top-up")
for  kmax in (20, 30, 40)
    a = omegas(kmax, false);   b = omegas(kmax, true)
    @printf("  %8d   %.4f  %.4f            %.4f  %.4f\n", kmax, a[1]/REF[1], a[2]/REF[2], b[1]/REF[1], b[2]/REF[2])
end
