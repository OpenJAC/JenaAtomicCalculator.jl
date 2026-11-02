# An independent check of the non-local exchange operator Continuum.exchangeMatrix: for a CLOSED-shell atom the Dirac-Fock
# orbitals satisfy (h + V_Hartree + X) phi = eps phi exactly, with X the very operator K represents (its b = a term removes the
# Hartree self-interaction). So the orbital energies of JAC's average-level Dirac-Fock SCF (AL; for a closed shell it is plain
# Dirac-Fock), a code that shares nothing with exchangeMatrix, must re-appear as eigenvalues of H_static + K in the same
# B-spline basis. Without K they must not.
#
# Run from <JAC-root>:   OPENBLAS_NUM_THREADS=2 julia --project=. tools/diag-exchange-vs-DiracFock.jl

using JenaAtomicCalculator, Printf, LinearAlgebra
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: Continuum, Bsplines

function check(Z, config, subshells; rbox=15.0)
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = rbox)
    nm   = Nuclear.Model(Float64(Z))
    set  = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction())
    lev  = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="DF", grid=grid, nuclearModel=nm, configs=[Configuration(config)],
                                   asfSettings=set), output=true)[JAC.ResultKeys.Multiplet].levels[1]
    end
    prim = Bsplines.generatePrimitives(grid);   np = Nuclear.nuclearPotential(nm, grid)
    pot  = Continuum.continuumPotential(0.0, lev, np, grid, :none)
    nL   = length(prim.bsplinesL);   nS = length(prim.bsplinesS)
    st   = Dict{String,Array{Float64,2}}()
    S    = zeros(nL+nS, nL+nS)
    S[1:nL,1:nL]             = Bsplines.generateTTpMatrix!("LL-overlap", 0, prim, st)
    S[nL+1:end,nL+1:end]     = Bsplines.generateTTpMatrix!("SS-overlap", 0, prim, st)
    println("\n=== $config (Z = $Z): AL Dirac-Fock orbital energies against eigenvalues of H_static (+ K) [Hartree] ===")
    @printf("  %-8s %16s %16s %12s %16s\n", "orbital", "Dirac-Fock", "H_static + K", "rel. diff", "H_static alone")
    for  s in subshells
        sh  = Subshell(s);   eDF = lev.basis.orbitals[sh].energy
        H0  = redirect_stdout(devnull) do;   Bsplines.generateGalerkinMatrix(sh, 0.0, pot, prim)   end
        K   = Continuum.exchangeMatrix(sh, lev, prim)
        near(M) = (ev = real.(eigvals(M, S));   ev[argmin(abs.(ev .- eDF))])
        eK  = near(H0 + K);   e0 = near(H0)
        @printf("  %-8s %16.8f %16.8f %12.2e %16.8f\n", s, eDF, eK, abs(eK - eDF)/abs(eDF), e0)
    end
end

check(10, "1s^2 2s^2 2p^6", ["1s_1/2", "2s_1/2", "2p_1/2", "2p_3/2"])
check(18, "1s^2 2s^2 2p^6 3s^2 3p^6", ["1s_1/2", "2s_1/2", "2p_3/2", "3s_1/2", "3p_1/2", "3p_3/2"]; rbox = 20.0)
