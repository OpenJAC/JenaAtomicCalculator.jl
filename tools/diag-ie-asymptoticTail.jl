# Does the asymptotic tail make the collision strength independent of the radial box?
#
# H-like C5+ 1s -> 2p, a dipole transition: its direct integrals carry the long-range 1/r^2 tail, and it is the case
# where example-Dl.jl found the Bethe limit. The 2p orbital has its outer turning point at ~1.1 a.u., so every box
# below is already much larger than the bound orbitals (Rule 12); only the continuum-continuum part should move.
#
#   (1) the quality of the Coulomb continuation of a few continuum orbitals (fitError, relative to the amplitude);
#   (2) sigma for rbox = 12, 30, 60 a.u., with asymptoticTail = false / true.
#
# Run from <JAC-root>:   julia --project=. -t 2 tools/diag-ie-asymptoticTail.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation, Continuum

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

const NM     = Nuclear.Model(6.)
const ENERGY = 735.0
const KMAX   = 30
grid(rbox)   = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = rbox)

# (1) fit quality
g   = grid(30.0)
mp  = Atomic.Computation(Atomic.Computation(), name="C5+ 1s", grid=g, nuclearModel=NM, configs=[Configuration("1s")])
lev = redirect_stdout(devnull) do;   perform(mp, output=true)[JAC.ResultKeys.Multiplet].levels[1]   end
println("\n=== (1) Coulomb continuation of the continuum orbitals, C5+ 1s, rbox = 30, E = $ENERGY eV, Z_ion = 5 ===")
for  kappa in (-1, 1, -2, 5, -10, 20)
    c, _ = redirect_stdout(devnull) do
        Continuum.generateOrbitalForLevel(ENERGY/27.211386, Subshell(101, kappa), lev, NM, g, Continuum.Settings(false, g.NoPoints-200))
    end
    t = ImpactExcitation.continuumTail(c, g, 5.0)
    @printf("  kappa = %4d   l = %3d   r0 = %7.3f   q = %7.4f   eta = % .4f   amplitude = %.4e   fitError = %.2e\n",
            kappa, t.l, t.r0, t.q, t.eta, sqrt(t.A^2 + t.B^2), t.fitError)
end

# (2) box study
function sigma(rbox, tail)
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [ENERGY], calcRateCoefficient = false,
                                    maxKappa = KMAX, printBefore = false, asymptoticTail = tail)
    wa  = Atomic.Computation(Atomic.Computation(), name="C5+", grid=grid(rbox), nuclearModel=NM,
                             initialConfigs=[Configuration("1s")], finalConfigs=[Configuration("2p")], processSettings=set)
    lines = redirect_stdout(devnull) do;   perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( [l.collisionStrength for l in lines] )
end

println("\n=== (2) Omega(1s -> 2p_1/2), Omega(1s -> 2p_3/2) against the box, E = $ENERGY eV, maxKappa = $KMAX ===")
@printf("  %6s   %-30s   %-30s\n", "rbox", "no tail", "with tail")
for  rbox in (12.0, 30.0, 60.0)
    a = sigma(rbox, false);   b = sigma(rbox, true)
    @printf("  %6.1f   %.6e  %.6e   %.6e  %.6e\n", rbox, a[1], a[2], b[1], b[2])
end
