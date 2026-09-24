# What orthogonalizing the continuum to the bound orbitals of the same kappa does to electron-impact excitation.
#
#   (1) the overlaps <b|eps kappa> before and after Continuum.orthogonalizeToBasis, for the eps s and eps p waves of C3+;
#   (2) sigma with orthogonalizeContinuum = false / true, with bound orbitals from the DFS field (the default, the SAME local
#       operator as the continuum, hence already orthogonal) and from a Hartree-Slater field (not orthogonal), for
#         Li-like C3+  1s2 2s -> 1s2 2p     one valence electron, direct-dominated
#         Be-like C2+  1s2 2s2 -> 1s2 2s2p  two valence electrons; the J=0,1,2 levels of 3P are reached by exchange
#
# Run from <JAC-root>:   julia --project=. -t 2 tools/diag-ie-orthogonality.jl

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation, Continuum, RadialIntegrals

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))

const GRID = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 30.0)
const NM   = Nuclear.Model(6.)

# (1) overlaps
mp  = Atomic.Computation(Atomic.Computation(), name="C3+ 2s", grid=GRID, nuclearModel=NM,
                         configs=[Configuration("1s^2 2s")], asfSettings=AsfSettings())
mlt = redirect_stdout(devnull) do;   perform(mp, output=true)[JAC.ResultKeys.Multiplet]   end
lev = mlt.levels[1]
println("\n=== (1) overlaps of the continuum with the bound orbitals of the same kappa, C3+ 1s2 2s, E = 50 eV ===")
for  sh in [Subshell(101, -1), Subshell(101, 1), Subshell(101, -2)]
    for  orth in (false, true)
        c, _ = redirect_stdout(devnull) do
            Continuum.generateOrbitalForLevel(50.0/27.211386, sh, lev, NM, GRID, Continuum.Settings(false, GRID.NoPoints-50);
                                              orthogonalize=orth)
        end
        ovs = [ @sprintf("<%s|c> = % .3e", string(s), RadialIntegrals.overlap(b, c, GRID))
                for (s, b) in lev.basis.orbitals if b.isBound && s.kappa == sh.kappa ]
        println("  kappa = $(sh.kappa)  orthogonalize = $(rpad(orth,5))  ", isempty(ovs) ? "(no bound orbital of this kappa)" : join(ovs, "   "))
    end
end

# (2) cross sections
function sigmas(iConf, fConf, energy, orth; field=Basics.DFSField())
    set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = [energy], calcRateCoefficient = false,
                                    maxKappa = 10, printBefore = false, orthogonalizeContinuum = orth)
    wa  = Atomic.Computation(Atomic.Computation(), name="diag", grid=GRID, nuclearModel=NM,
                             initialConfigs=[iConf], finalConfigs=[fConf], processSettings=set,
                             initialAsfSettings=AsfSettings(AsfSettings(); scField=field),
                             finalAsfSettings=AsfSettings(AsfSettings(); scField=field))
    lines = redirect_stdout(devnull) do;   perform(wa, output=true)[ImpactExcitation.Settings]   end
    return( lines )
end

for  (label, ic, fc) in [ ("Li-like C3+ 1s2 2s -> 1s2 2p",   Configuration("1s^2 2s"),   Configuration("1s^2 2p")),
                          ("Be-like C2+ 1s2 2s2 -> 1s2 2s2p", Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p")) ]
  for  field in (Basics.DFSField(), Basics.HSField())
    a = sigmas(ic, fc, 50.0, false; field=field);   b = sigmas(ic, fc, 50.0, true; field=field)
    println("\n=== (2) $label,  bound orbitals from $field,  E = 50 eV,  maxKappa = 10 ===")
    @printf("  %-8s %-8s %14s %14s %10s\n", "J_i", "J_f (#)", "sigma plain", "sigma orth", "orth/plain")
    for  i in eachindex(a)
        @printf("  %-8s %-8s %14.6e %14.6e %10.4f\n", string(a[i].initialLevel.J), string(a[i].finalLevel.J) * " (" *
                string(a[i].finalLevel.index) * ")", a[i].crossSection, b[i].crossSection, b[i].crossSection/a[i].crossSection)
    end
  end
end
