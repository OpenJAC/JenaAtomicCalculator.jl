# Do the Coulomb-Bethe channels belong in the DCS?  Mg+ 3s -> 3p_1/2 at 50 eV, maxKappa = 40, tail on, static waves:
# without and with coulombBetheKappa = 200. Two tests:
#   (1) 2 pi int dcs sin(theta) dtheta must reproduce sigma in both runs (then the CB channels enter coherently and correctly
#       normalized);
#   (2) the DCS against the relativistic distorted-wave code (Fortran_Codes/magnesium/.../dcs2p1-50ev.dat: angle [deg],
#       DCS [a0^2/sr]), which converges the small-angle peak by Born closure.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. -t 2 tools/diag-ie-dcsCoulombBethe.jl

using JenaAtomicCalculator, Printf, DelimitedFiles
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

const GRID   = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 25.0)
const ANGLES = collect(0.0:0.5:180.0)

function line12(cb)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [50.0], calcRateCoefficient = false,
                                    maxKappa = 40, printBefore = false, asymptoticTail = true, continuumExchange = :none,
                                    coulombBetheKappa = cb, calcDCS = true, angles = ANGLES)
    wa  = Atomic.Computation(Atomic.Computation(), name="Mg+", grid=GRID, nuclearModel=Nuclear.Model(12.),
                             initialConfigs=[Configuration("[Ne] 3s")], finalConfigs=[Configuration("[Ne] 3p")], processSettings=set)
    lines = redirect_stdout(devnull) do;   perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( lines[findfirst(l -> l.finalLevel.J == AngularJ64(1//2), lines)] )
end

integral(l) = (th = deg2rad.(l.angles);   f = l.dcs .* sin.(th);   2pi * sum((th[i+1]-th[i]) * (f[i]+f[i+1]) / 2 for i in 1:length(th)-1))

a = line12(0);   b = line12(200)
println("\n=== (1) Mg+ 3s -> 3p_1/2, 50 eV: 2 pi int dcs dOmega / sigma ===")
@printf("  maxKappa 40            : %.6f   (sigma = %.4f a0^2)\n", integral(a)/a.crossSection, a.crossSection)
@printf("  maxKappa 40 + CB to 200: %.6f   (sigma = %.4f a0^2)\n", integral(b)/b.crossSection, b.crossSection)

ref = readdlm("Fortran_Codes/magnesium/magnesium/dcs2p1-50ev.dat")
println("\n=== (2) DCS [a0^2/sr] against the RDW code ===")
@printf("  %6s %12s %12s %12s %10s\n", "theta", "RDW", "JAC k40", "JAC k40+CB", "CB/RDW")
for  th in (0, 1, 2, 3, 5, 10, 20, 40, 60, 90, 120, 150, 180)
    k = findfirst(==(Float64(th)), ANGLES);   r = ref[findfirst(==(th), Int.(ref[:,1])), 2]
    @printf("  %6d %12.4f %12.4f %12.4f %10.4f\n", th, r, a.dcs[k], b.dcs[k], b.dcs[k]/r)
end
