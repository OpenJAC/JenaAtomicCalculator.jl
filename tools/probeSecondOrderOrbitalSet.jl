#
# probeSecondOrderOrbitalSet.jl
#
# A PROBE, not a module:  priority item 1's OPEN question.  Both validity conditions of a second-order
# amplitude are settled (reference: the Green space must be disjoint from the two CI spaces, and all three
# multiplets must share ONE orbital set).  What has never been measured is WHICH one set.
#
# The 0.996749 that `testModule_TwoElectronOnePhoton` asserts was obtained with an SCF over the UNION of all
# three spaces.  The item's OPEN paragraph recommends the INITIAL state's set instead, on the argument that
# the series is referenced to E_i.  Nobody has run the alternatives, so the recommendation rests on an
# argument rather than on a number.
#
# ONE VARIABLE IS CHANGED per run (nothing else moves):  the configuration list the single SCF is performed
# over.  Everything after that line -- grid, nuclear model, AsfSettings, the five CI spaces, both routes and
# the line selection -- is copied verbatim from the suite test, so the UNION row below must reproduce
# 0.996749 and that is this probe's own calibration.
#
# TWO numbers are reported per variant, because they answer different questions:
#   ratio   = route(ii)/route(i) on THAT set -- does the second-order construction reproduce the full CI
#             on the basis it was built on?  This is the validity question.
#   A_CI    = the full-CI rate itself -- which set gives the better-converged physics?  A ratio of 1.000 on
#             a poor basis means only that both routes are equally wrong.
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.2e-2, rbox = 10.0)
nm   = Nuclear.Model(10.);              asf  = AsfSettings()
cRef = [Configuration("2s 2p")]
cG   = [Configuration("1s 2p"), Configuration("1s 2s")]
c2nd = [Configuration("1s^2"), Configuration("2p^2")]
cCIi = [Configuration("2s 2p"), Configuration("1s 2p")]
cCIf = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("2p^2")]

# The variants.  Each must span the subshells 1s, 2s, 2p or the CI spaces cannot be built on it; the initial
# configuration 2s2p alone does NOT span 1s, which is itself part of the answer and is noted in the report.
variants = [ ("UNION  (the suite test)",      vcat(cRef, cG, cCIf)),
             ("INITIAL + final (no Green)",   vcat(cRef, cCIf)),
             ("FINAL CI space alone",         cCIf),
             ("GREEN space alone",            cG),
             ("INITIAL + 1s^2 only",          vcat(cRef, [Configuration("1s^2")])) ]

function runVariant(name, scfConfigs)
    orbitals = SelfConsistent.performSCF(scfConfigs, nm, grid, asf; printout=false).levels[1].basis.orbitals
    iMp  = Hamiltonian.performCIwithFrozenOrbitals(cRef, orbitals, nm, grid, asf; printout=false)
    gMp  = Hamiltonian.performCIwithFrozenOrbitals(cG,   orbitals, nm, grid, asf; printout=false)
    fMp  = Hamiltonian.performCIwithFrozenOrbitals(c2nd, orbitals, nm, grid, asf; printout=false)
    ciMp = Hamiltonian.performCIwithFrozenOrbitals(cCIi, orbitals, nm, grid, asf; printout=false)
    cfMp = Hamiltonian.performCIwithFrozenOrbitals(cCIf, orbitals, nm, grid, asf; printout=false)

    teopSettings = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,UseBabushkin], false,
                                                 LineSelection(true, indexPairs=[(i,1) for i = 1:4]), 0.,
                                                 CoulombInteraction(), gMp)
    teopLines = TwoElectronOnePhoton.computeLines(fMp, iMp, nm, grid, teopSettings; output=true)
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1],
                        gauges=[UseCoulomb,UseBabushkin], printBefore=false, calcBiorthogonal=false,
                        lineSelection=LineSelection(true, indexPairs=[(i,1) for i = 1:8]) )
    photoLines = PhotoEmission.computeLines(cfMp, ciMp, grid, photoSettings; output=true)
    conditions = TwoElectronOnePhoton.checkSecondOrderConditions(fMp, iMp, gMp, grid)

    teopLine  = teopLines[argmax([line.teopRate.Coulomb  for line in teopLines])]
    toGround  = [line  for line in photoLines  if line.finalLevel.index == 1]
    photoLine = toGround[argmin([abs(line.omega - teopLine.omega)  for line in toGround])]
    aTeop  = teopLine.teopRate.Coulomb
    aPhoto = photoLine.photonRate.Coulomb
    return( (name = name, ratio = aTeop/aPhoto, aTeop = aTeop, aPhoto = aPhoto,
             dOmega = abs(photoLine.omega - teopLine.omega)/teopLine.omega,
             omega = teopLine.omega, conditions = conditions) )
end

results = []
for (name, configs) in variants
    @printf("\n>>>> VARIANT: %s\n", name)
    try
        push!(results, runVariant(name, configs))
    catch err
        @printf("     RAISED: %s\n", sprint(showerror, err)[1:min(end,200)])
        push!(results, (name = name, ratio = NaN, aTeop = NaN, aPhoto = NaN, dOmega = NaN,
                        omega = NaN, conditions = ["RAISED"]))
    end
end

println("\n", "="^118)
println("PRIORITY ITEM 1, THE OPEN QUESTION:  which orbital set does a second-order amplitude live on?")
println("He-like Ne, 2s2p -> 1s^2 two-electron one-photon, Coulomb gauge; one SCF per row, everything else identical.")
println("="^118)
@printf("  %-28s %12s %14s %14s %10s  %s\n", "SCF performed over", "ratio ii/i", "A_2nd [a.u.]", "A_CI [a.u.]",
        "d(omega)", "conditions")
println("  ", "-"^116)
for r in results
    @printf("  %-28s %12.6f %14.5e %14.5e %10.2e  %s\n", r.name, r.ratio, r.aTeop, r.aPhoto, r.dOmega,
            isempty(r.conditions) ? "met" : string(r.conditions))
end
println("  ", "-"^116)
base = findfirst(r -> startswith(r.name, "UNION"), results)
if base !== nothing  &&  !isnan(results[base].aPhoto)
    println("\n  A_CI relative to the UNION set (which basis is the better physics):")
    for r in results
        isnan(r.aPhoto) && continue
        @printf("     %-28s  A_CI / A_CI(union) = %8.4f\n", r.name, r.aPhoto/results[base].aPhoto)
    end
end
println()
