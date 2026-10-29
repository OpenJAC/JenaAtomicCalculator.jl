#
#  [CONT-XC]  Is exchange counted TWICE when a continuum orbital is generated with a nonlocal exchange operator
#             on top of a Dirac-Fock-SLATER local potential?
#
#  WHY THIS EXISTS.  `Continuum.generateOrbitalForLevel` builds its local part as
#  `Basics.computePotential(Basics.DFSField(1.0), grid, level)`, and the DFS potential is BY DEFINITION
#      V_DFS(r) = int rho_t(r')/r_>  -  [ 3/(4 pi^2 r^2) rho_t(r) ]^(1/3) * strength ,
#  whose second term IS Slater's LOCAL exchange;  `strength` multiplies that term and nothing else (verified in
#  `Basics.computePotential(::DFSField, ...)`).  Adding `Bsplines.setupExchangeMatrix` -- the full NONLOCAL Fock
#  exchange -- on top of it therefore counts exchange twice, once locally and once nonlocally.
#
#  THE DISCRIMINATOR IS REFERENCE-FREE, which is why it is worth having: a bound orbital and a continuum orbital
#  of the SAME kappa are eigenfunctions of the same self-adjoint operator at different eigenvalues ONLY IF the
#  continuum side carries the same exchange the bound side did -- and then their overlap must VANISH.  So
#  <2p_1/2 | eps p_1/2> names which operator the bound side actually used, with no literature needed.
#
#  MEASURED 29-Oct-2026, Ne, E = 2.0 Ha, normalized overlap:
#
#      slater   nonlocal   <2p_1/2|eps p_1/2>
#        0.0     false       -8.15e-02      ... NO exchange at all, and the worst of the four
#        0.0     true        +1.08e-03      ... pure direct + nonlocal Fock: 76x BETTER than no exchange, so the
#                                               nonlocal operator is doing real and correct work
#        1.0     true        -2.52e-02      ... the Slater term ON TOP of it gives back most of that gain: this is
#                                               the DOUBLE-COUNTING signature, and it is visible from both sides
#        1.0     false       -1.46e-08      ... the CONSISTENT pair: this level's bound orbitals ARE DFS orbitals,
#                                               so a DFS continuum orbital is orthogonal to them to round-off
#
#  SO THE FINDING IS NOT "exchange is wrong" BUT "exchange is right and the Slater term must come off with it".
#  Read rows 1 and 2 for the operator (76x), rows 2 and 3 for the double counting (23x).
#
#  WHAT THE RESIDUAL 1.08e-03 MEANS, and it is NOT a defect: it cannot go to zero here, because the bound
#  orbitals of this level are DFS orbitals and no Dirac-Fock continuum operator can be orthogonal to those.
#  CONCLUSION: a nonlocal exchange term is only meaningful when (a) the local part is PURE DIRECT, i.e.
#  DFSField(0.0), and (b) the level's own bound orbitals came from a Dirac-Fock field.  THE DECISIVE FOLLOW-UP
#  is to repeat row 2 with bound orbitals from an actual Fock field and check the overlap falls to ~1e-8.
#
#  RUN:  julia --project=. tools/probeContinuumExchangeConsistency.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

##  THIS PROBE NEEDS A CAPABILITY THAT NEED NOT BE IN THIS VERSION.  The nonlocal-exchange rows call
##  Continuum.generateOrbitalGalerkin with `coreOrbitals`/`occupations` keywords;  where those do not exist, say so
##  and run the local-only reference row alone, rather than dying on a MethodError that reads like a defect.
const hasExchangeKeywords = any(m -> :coreOrbitals in Base.kwarg_decl(m),
                                methods(Continuum.generateOrbitalGalerkin))

grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
nm   = Nuclear.Model(10.)
wa   = perform(Atomic.Computation(Atomic.Computation(), name="Ne", nuclearModel=nm, grid=grid,
                                  configs=[Configuration("1s^2 2s^2 2p^6")]); output=true)
lev  = wa[ResultKeys.Multiplet].levels[1]
setDefaults("method: continuum, Galerkin")
sh    = Subshell(100, 1);   occSh = Subshell(2, 1);   E = 2.0
mtp   = Continuum.gridConsistency(E, grid)
prim  = Bsplines.generatePrimitives(grid)
nucP  = Nuclear.nuclearPotential(nm, grid)
coreO = lev.basis.orbitals;   occs = Basics.extractMeanOccupation(lev.basis)
println("\n   slater   nonlocal   <2p_1/2|eps p_1/2>      interpretation")
rows = hasExchangeKeywords ? [(1.0, false), (1.0, true), (0.0, true), (0.0, false)] : [(1.0, false), (0.0, false)]
if  !hasExchangeKeywords
    println("   (no `coreOrbitals` keyword in this version -- the nonlocal rows are skipped; the table in this")
    println("    file's header is the record of what they showed when they could be run.)")
end
for  (slater, nonlocal)  in  rows
    pot = Basics.add(nucP, Basics.computePotential(Basics.DFSField(slater), grid, lev))
    co  = nonlocal ?
          Continuum.generateOrbitalGalerkin(E, sh, pot, Continuum.Settings(true, mtp); primitives=prim,
                                            coreOrbitals=coreO, occupations=occs) :
          Continuum.generateOrbitalGalerkin(E, sh, pot, Continuum.Settings(false, mtp); primitives=prim)
    ov = JAC.RadialIntegrals.overlap(lev.basis.orbitals[occSh], co, grid)
    nrm = sqrt(JAC.RadialIntegrals.overlap(co, co, grid))
    @printf("   %5.1f    %-8s   %+.6e  (norm %.3e)\n", slater, nonlocal, ov/nrm, nrm)
end
println("\n   The bound orbitals of this level came from JAC's DEFAULT field;  whichever row reproduces ~0 names")
println("   the operator the bound side actually used.")
