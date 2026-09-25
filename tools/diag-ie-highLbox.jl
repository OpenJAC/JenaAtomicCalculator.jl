# Is the brute-force maxKappa = 120, rbox = 30 reference of example-Dl branch a inflated at high l by the box?
# At 2939 eV the direct dipole integral of a high partial wave is a near-cancellation between the part inside the box and
# the part beyond it (l = 80: grid -3.8e-5, tail +2.8e-5), so a calculation WITHOUT the tail overstates high-l waves.
# Distorted waves only, maxKappa = 60, C5+ 1s -> 2p at 2939 eV; Omega divided by the maxKappa = 120 reference.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-highLbox.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))

const REF = [9.638804e-02, 1.927089e-01]
function om(rbox, tail)
    g   = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = rbox)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [2939.], calcRateCoefficient = false,
                                    maxKappa = 60, printBefore = false, asymptoticTail = tail)
    wa  = Atomic.Computation(Atomic.Computation(), name="C5+", grid=g, nuclearModel=Nuclear.Model(6.),
                             initialConfigs=[Configuration("1s")], finalConfigs=[Configuration("2p")], processSettings=set)
    l = redirect_stdout(devnull) do;  perform(wa, output=true)["impact-excitation lines:"]  end
    return( [x.collisionStrength for x in l] ./ REF )
end

println("\n=== C5+ 1s -> 2p, 2939 eV, distorted waves, maxKappa = 60: Omega / Omega(maxKappa = 120, rbox = 30, no tail) ===")
for  rbox in (30.0, 60.0)
    a = om(rbox, false);   b = om(rbox, true)
    @printf("  rbox %5.1f   no tail %.4f %.4f    with tail %.4f %.4f\n", rbox, a[1], a[2], b[1], b[2])
end
