# Is the SHAPE of the DCS right?  The integral 2 pi int dcs dOmega = sigma cannot tell, since it is blind to the relative
# phases of the partial waves. The exact reference is the first Born approximation for H 1s -> 2p at high energy:
#     dsigma/dOmega = (k_f/k_i) (4/q^4) 3 |I(q)|^2,   I(q) = int P_2p P_1s j_1(qr) dr,
# with the analytic hydrogen orbitals -- independent of every part of JAC's scattering code. At 1000 eV the distorted-wave
# DCS must approach it at small and intermediate angles.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-dcsBorn.jl

using JenaAtomicCalculator, Printf, SpecialFunctions
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

const E   = 1000.0
const ANG = [0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 90.0]

function bornDCS(theta)
    ki = sqrt(2E/27.211386);   kf = sqrt(2(E - 10.2043)/27.211386)
    q  = sqrt(ki^2 + kf^2 - 2ki*kf*cosd(theta))
    j1(x) = x < 1e-6 ? x/3 : sin(x)/x^2 - cos(x)/x
    h = 1e-3;   I = sum( (r = i*h;   (r^2*exp(-r/2)/(2sqrt(6))) * (2r*exp(-r)) * j1(q*r) * h) for i = 1:60000 )
    return( kf/ki * 4/q^4 * 3 * I^2 )
end

grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 5.0e-3, rbox = 20.0)
set  = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [E], calcRateCoefficient = false,
                                 maxKappa = 30, printBefore = false, asymptoticTail = true, coulombBetheKappa = 400,
                                 continuumExchange = :none, calcDCS = true, angles = ANG)
wa   = Atomic.Computation(Atomic.Computation(), name="H", grid=grid, nuclearModel=Nuclear.Model(1., UniformNucleus(), 1., 0.8797, AngularJ64(1//2), 0.0, 0.0, 0.0),
                          initialConfigs=[Configuration("1s")], finalConfigs=[Configuration("2p")], processSettings=set)
lines = redirect_stdout(devnull) do;   perform(wa, output=true)[ImpactExcitation.Settings]   end
jac   = sum(l.dcs for l in lines)                 # 2p_1/2 + 2p_3/2
println("\n=== H 1s -> 2p at $E eV: DCS [a0^2/sr], JAC (2p_1/2 + 2p_3/2) against first Born ===")
@printf("  %6s %14s %14s %10s\n", "theta", "Born", "JAC", "JAC/Born")
for  (i, th) in enumerate(ANG)
    b = bornDCS(th);   @printf("  %6.1f %14.5e %14.5e %10.4f\n", th, b, jac[i], jac[i]/b)
end
@printf("  sigma JAC = %.5f a0^2\n", sum(l.crossSection for l in lines))
