#
# probe-cachingCost.jl
#
# What does caching actually buy in the CSF-pair matrix build, and how does the cost scale with the
# number of CSFs?  This decides whether the EXPLICIT route (cost ~ |P_CSF| x |Q_CSF|) can ever compete
# with Gaigalas' IMPLICIT route (cost ~ |P_conf| x |Q_conf|) for an open d- or f-shell.
#
# Two caches are separable:  SpinAngular's shell-level memo (angular) and InteractionStrength.XLCache
# (radial Slater integrals).  Timed cold and warm, independently.
#
using JenaAtomicCalculator, LinearAlgebra, Printf
const JAC = JenaAtomicCalculator

nm       = Nuclear.Model(17.)
refConf  = Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")
sym      = LevelSymmetry(AngularJ64(3//2), Basics.minus)
coreSh   = [Shell("2s"), Shell("2p")];   valSh = [Shell("3s"), Shell("3p")]
grid     = Basics.recommendedGrid([refConf], nm; rbox = 20.0, printout=false)
settings = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                      eeInteractionCI = CoulombInteraction(), gridStopper = false)
primitives = Bsplines.generatePrimitives(grid)
mpRef      = SelfConsistent.performSCF([refConf], nm, grid, settings, printout=false)
pot        = Basics.add( Nuclear.nuclearPotential(nm, grid),
                         Basics.computePotential(Basics.DFSField(), grid, mpRef.levels[1].basis) )

println("="^104)
@printf("%-26s %7s %11s %11s %11s %11s %11s\n", "case", "N CSF",
        "cold [s]", "warm ang", "warm rad", "warm both", "us/pair")
println("="^104)

for (tag, fromSh, toSh) in [("C  -> 3d",        coreSh,               [Shell("3d")]),
                            ("C  -> 3d,4s",     coreSh,               [Shell("3d"), Shell("4s")]),
                            ("VV -> 3d,4s",     valSh,                [Shell("3d"), Shell("4s")]),
                            ("ALL-> 3d,4s",     vcat(coreSh,valSh),   [Shell("3d"), Shell("4s")])]
    confS    = Basics.generateConfigurations([refConf], fromSh, toSh)
    confD    = Basics.generateConfigurations(confS, fromSh, toSh)
    confList = Configuration[]
    for confa in vcat([refConf], confS, confD)
        addTo = true;  for confb in confList   if confa == confb   addTo = false; break  end  end
        addTo && push!(confList, confa)
    end
    basis0  = redirect_stdout(devnull) do;  Basics.generateBasis(confList, [sym])  end
    allOrbs = Bsplines.generateOrbitals(basis0.subshells, pot, nm, primitives; printout=false)
    basis   = Basis(true, basis0.NoElectrons, basis0.subshells, basis0.csfs, basis0.coreSubshells, allOrbs)
    n       = length(basis.csfs);   npair = n*(n+1)/2

    # cold: both caches empty
    SpinAngular.clearCaches();   c1 = InteractionStrength.XLCache()
    tCold = @elapsed Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, c1, printout=false)
    # warm angular only (fresh radial cache)
    c2 = InteractionStrength.XLCache()
    tWarmA = @elapsed Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, c2, printout=false)
    # warm radial only (fresh angular)
    SpinAngular.clearCaches()
    tWarmR = @elapsed Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, c2, printout=false)
    # warm both
    tWarmB = @elapsed Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, c2, printout=false)

    @printf("%-26s %7d %11.3f %11.3f %11.3f %11.3f %11.2f\n",
            tag, n, tCold, tWarmA, tWarmR, tWarmB, 1e6*tWarmB/npair)
end
println("="^104)
println("cold      = both caches empty      warm ang = angular memo full, radial fresh")
println("warm rad  = radial full, angular cleared   warm both = steady state")
println("us/pair is the STEADY-STATE cost of one CSF-pair matrix element -- the number that decides")
println("whether |P_CSF| x |Q_CSF| is affordable for an open d/f shell.")
