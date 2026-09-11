#
# probe-wignerCacheValue.jl
#
# What are SpinAngular's memos actually worth, now that AngularMomentum's default Wigner method is
# FloatingWigner (25-Aug-2026) rather than ExactWigner?  The memos' own rationale is written in terms of
# BigInt rational arithmetic, which that change removed -- so the rationale needs re-measuring, not re-reading.
#
# Two measurements, because they answer different questions:
#   (A) END TO END.  setupMatrix with the memos cleared before EVERY call, against the memos left warm.
#       This is the FIRST-FILL cost only: within one call the memo fills and serves the rest of that call.
#   (B) PER CALL.  shellReducedW / shellReducedA against their own *Uncached bodies, on keys drawn from a
#       real run.  Times the number of repeat lookups, this is the memo's TOTAL value.
#
using JenaAtomicCalculator, LinearAlgebra, Printf, Statistics
const JAC = JenaAtomicCalculator
const SA  = JenaAtomicCalculator.SpinAngular

function buildBasis(refConf, sym, fromSh, toSh, nm, grid, pot, primitives)
    confS    = Basics.generateConfigurations([refConf], fromSh, toSh)
    confD    = Basics.generateConfigurations(confS, fromSh, toSh)
    confList = Configuration[]
    for confa in vcat([refConf], confS, confD)
        addTo = true;  for confb in confList  if confa == confb  addTo = false; break end end
        addTo && push!(confList, confa)
    end
    basis0  = redirect_stdout(devnull) do; Basics.generateBasis(confList, [sym]) end
    orbs    = Bsplines.generateOrbitals(basis0.subshells, pot, nm, primitives; printout=false)
    return( Basis(true, basis0.NoElectrons, basis0.subshells, basis0.csfs, basis0.coreSubshells, orbs) )
end

cases = [("Cl III 3s^2 3p^3  (short trees)", 17.0, "1s^2 2s^2 2p^6 3s^2 3p^3", AngularJ64(3//2), Basics.minus,
          [Shell("3s"),Shell("3p")], [Shell("3d"),Shell("4s")]),
         ("Ti III 3d^2       (long trees)",  22.0, "1s^2 2s^2 2p^6 3s^2 3p^6 3d^2", AngularJ64(2), Basics.plus,
          [Shell("3d")], [Shell("4s"),Shell("4p"),Shell("4d")])]

println("="^110)
@printf("%-34s %6s %11s %11s %9s %10s %10s\n", "case","N CSF","cold [s]","warm [s]","memo %","W keys","A keys")
println("="^110)

for (tag, Z, cstr, J, par, fromSh, toSh) in cases
    nm   = Nuclear.Model(Z);   refConf = Configuration(cstr);   sym = LevelSymmetry(J, par)
    grid = Basics.recommendedGrid([refConf], nm; rbox = 20.0, printout=false)
    set  = AsfSettings(AsfSettings(); scField=Basics.ALField(), eeInteraction=CoulombInteraction(),
                                      eeInteractionCI=CoulombInteraction(), gridStopper=false)
    prim = Bsplines.generatePrimitives(grid)
    mpR  = SelfConsistent.performSCF([refConf], nm, grid, set, printout=false)
    pot  = Basics.add(Nuclear.nuclearPotential(nm, grid),
                      Basics.computePotential(Basics.DFSField(), grid, mpR.levels[1].basis))
    bas  = buildBasis(refConf, sym, fromSh, toSh, nm, grid, pot, prim)

    Hamiltonian.setupMatrix(sym, bas, nm, grid, set, InteractionStrength.XLCache(), printout=false)  # compile
    cold = Float64[];  warm = Float64[]
    for r = 1:3
        SA.clearCaches()
        push!(cold, @elapsed Hamiltonian.setupMatrix(sym, bas, nm, grid, set, InteractionStrength.XLCache(), printout=false))
    end
    c = InteractionStrength.XLCache()
    Hamiltonian.setupMatrix(sym, bas, nm, grid, set, c, printout=false)
    for r = 1:3
        push!(warm, @elapsed Hamiltonian.setupMatrix(sym, bas, nm, grid, set, c, printout=false))
    end
    @printf("%-34s %6d %11.3f %11.3f %8.1f%% %10d %10d\n", tag, length(bas.csfs),
            median(cold), median(warm), 100*(median(cold)-median(warm))/median(cold),
            length(SA.SHELL_W_CACHE), length(SA.SHELL_A_CACHE))
end
println("="^110)

# (B) what one memoised call saves, measured directly against its own uncached body
println("\nPER-CALL cost, cached lookup vs the uncached body (1e5 calls each):")
SA.clearCaches()
jj = AngularJ64(5//2)
SA.shellReducedW(jj, 3, 1, AngularJ64(5//2), 3, AngularJ64(9//2), 2)       # prime the memo
tC = @elapsed for i=1:100000  SA.shellReducedW(jj, 3, 1, AngularJ64(5//2), 3, AngularJ64(9//2), 2)          end
tU = @elapsed for i=1:100000  SA.shellReducedWUncached(jj, 3, 1, AngularJ64(5//2), 3, AngularJ64(9//2), 2)  end
@printf("   shellReducedW    cached %8.4f s   uncached %8.4f s   speed-up %6.1fx\n", tC, tU, tU/tC)
SA.shellReducedA(jj, 3, 1, AngularJ64(5//2), 2, 2, AngularJ64(2), AngularM64(1//2))
tC = @elapsed for i=1:100000  SA.shellReducedA(jj, 3, 1, AngularJ64(5//2), 2, 2, AngularJ64(2), AngularM64(1//2))          end
tU = @elapsed for i=1:100000  SA.shellReducedAUncached(jj, 3, 1, AngularJ64(5//2), 2, 2, AngularJ64(2), AngularM64(1//2))  end
@printf("   shellReducedA    cached %8.4f s   uncached %8.4f s   speed-up %6.1fx\n", tC, tU, tU/tC)
println("\nDONE.")
