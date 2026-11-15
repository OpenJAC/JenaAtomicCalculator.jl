#
# probe-pooledPseudostateTiming.jl -- localizes where time and memory go in the pooled, two-kappa, resonance-
# matched Fe K-LL pseudostate CI (priority-list item 17 / the sparse-solver wiring plan's own open blocker).
#
# Exactly example-Dz.jl branch b's system (Fe K-LL), but replicates PhotoIonizationPseudostates.computeSpectrum's
# own body step by step -- ladder build per kappa, CSF/basis assembly, Hamiltonian.performCI, the amplitude loop
# -- with an explicit, FLUSHED timing line after each, so a background run redirected to a file still shows
# where the wall-clock actually goes instead of staying silent until completion or death.
#
# Run from the JAC root:   julia --project=. -t 8 tools/probe-pooledPseudostateTiming.jl
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra, Distributed, Profile

setDefaults("unit: energy", "eV")

windowLo = parse(Float64, get(ENV, "PROBE_WINDOW_LO", "171.0"))
windowHi = parse(Float64, get(ENV, "PROBE_WINDOW_HI", "171.6"))
doProfile = get(ENV, "PROBE_PROFILE", "0") == "1"

function  stage(f::Function, label::String)
    print(">>> [", label, "] starting at ", Dates.format(Dates.now(), "HH:MM:SS"), " ... ")
    flush(stdout)
    t0  = time()
    out = f()
    dt  = time() - t0
    @printf("done in %.1f s\n", dt)
    flush(stdout)

    return( out )
end

using Dates

println("probe-pooledPseudostateTiming.jl: Fe K-LL, pooled kappas -2/+1, rbox=17000, window=(171.0,171.6)")
flush(stdout)

grid   = Radial.Grid( Radial.Grid(false); hp=0.02, rbox=4.0 )
nModel = Nuclear.Model(26.)

iMult = stage("initial multiplet SCF") do
    SelfConsistent.performSCF([Configuration("1s^2 2s^2"), Configuration("1s 2s^2 2p")], nModel, grid,
                               AsfSettings(); printout=false)
end
initialLevel = iMult.levels[1]

fMult = stage("final-ion frozen-orbital CI") do
    Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s^2 2s")], initialLevel.basis.orbitals,
                                             nModel, grid, AsfSettings(); printout=false)
end
finalIonLevel = fMult.levels[1]

thrEv = Defaults.convertUnits("energy: from atomic to eV", finalIonLevel.energy - initialLevel.energy)
@printf("\n>>> ionization threshold = %.4f eV\n\n", thrEv)

println(">>> pseudostateWindow = ($windowLo, $windowHi),  PROBE_PROFILE=$doProfile")
flush(stdout)
ppSet = PhotoIonizationPseudostates.Settings(PhotoIonizationPseudostates.Settings(); multipoles=[E1],
            gauges=[UseCoulomb, UseBabushkin], pseudostateWindow=(windowLo, windowHi), pseudostateRbox=17000.0,
            printBefore=false)

kappas   = [-2, 1]
symmetry = LevelSymmetry(AngularJ64(1), Basics.minus)

ladderGrid = stage("ladderGridFor") do
    PhotoIonizationPseudostates.ladderGridFor(grid, ppSet)
end

stage("checkFrozenOrbitalConsistency") do
    PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(initialLevel, finalIonLevel, ladderGrid)
end

laddersByKappa = Dict{Int64,PhotoIonizationPseudostates.KappaLadder}()
for  kappa  in  kappas
    laddersByKappa[kappa] = stage("buildPseudostateLadder, kappa=$kappa") do
        PhotoIonizationPseudostates.buildPseudostateLadder(kappa, finalIonLevel, nModel, ladderGrid, grid, ppSet)
    end
    println("      -> ", length(laddersByKappa[kappa].states), " pseudostate(s) retained for kappa = ", kappa)
    flush(stdout)
end

ladderSubshells = vcat( [collect(keys(laddersByKappa[kappa].states))  for kappa in kappas]... )
println(">>> total pseudostate subshells across both kappas: ", length(ladderSubshells))
flush(stdout)

dominantIdx    = argmax( abs.(finalIonLevel.mc) )
coreCsf        = finalIonLevel.basis.csfs[dominantIdx]
coreOccupation = Dict{Subshell,Int64}( finalIonLevel.basis.subshells[i] => coreCsf.occupation[i]
                                        for i = 1:length(finalIonLevel.basis.subshells)
                                        if  coreCsf.occupation[i] > 0 )

extendedSubshells, paddedBasis, boundCsfs, combinedCsfs, kappaRanges = stage("CSF generation + basis assembly") do
    extendedSubshells = Base.sort( vcat(initialLevel.basis.subshells, ladderSubshells), lt=Base.isless )
    paddedBasis       = ManyElectron.Basis(initialLevel.basis, extendedSubshells)
    boundCsfs         = CsfR[];   for  csf in paddedBasis.csfs
                                      LevelSymmetry(csf.J, csf.parity) == symmetry   &&   push!(boundCsfs, csf)   end
    combinedCsfs = copy(boundCsfs)
    kappaRanges  = Dict{Int64,UnitRange{Int64}}()
    for  kappa  in  kappas
        pseudoRelConfs = ConfigurationR[]
        for  sh  in  keys(laddersByKappa[kappa].states)
            subshells = merge( coreOccupation, Dict(sh => 1) )
            push!( pseudoRelConfs, ConfigurationR(subshells, finalIonLevel.basis.NoElectrons + 1) )
        end
        pseudoCsfs = CsfR[]
        for  relconf in pseudoRelConfs
            newCsfs = Basics.generateCsfRs(relconf, extendedSubshells)
            for  csf in newCsfs    LevelSymmetry(csf.J, csf.parity) == symmetry   &&   push!(pseudoCsfs, csf)    end
        end
        startIdx = length(combinedCsfs) + 1
        append!( combinedCsfs, pseudoCsfs )
        kappaRanges[kappa] = startIdx:length(combinedCsfs)
    end

    return( extendedSubshells, paddedBasis, boundCsfs, combinedCsfs, kappaRanges )
end
println(">>> combined basis: ", length(boundCsfs), " bound CSF(s) + ", length(combinedCsfs) - length(boundCsfs),
        " pseudostate CSF(s) = ", length(combinedCsfs), " total")
flush(stdout)

cacheSize = parse(Int64, get(ENV, "PROBE_MAXCACHE", "64"))
println(">>> LazyLadderOrbitals maxCache = ", cacheSize)
flush(stdout)
combinedOrbitals = PhotoIonizationPseudostates.LazyLadderOrbitals(initialLevel.basis.orbitals, laddersByKappa;
                                                                   maxCache=cacheSize)
coreSubshellList = Subshell[]
for  (k, sh)  in  enumerate(extendedSubshells)
    mocc = Basics.subshell_2j(sh) + 1;    isFilled = true
    for  csf in combinedCsfs    if  csf.occupation[k] != mocc    isFilled = false;   break   end   end
    isFilled   &&   push!(coreSubshellList, sh)
end
combinedBasis = Basis(true, initialLevel.basis.NoElectrons, extendedSubshells, combinedCsfs, coreSubshellList,
                      combinedOrbitals)

gcBefore = Base.gc_num()
mp = stage("Hamiltonian.performCI over the combined basis") do
    if  doProfile
        Profile.clear()
        @profile Hamiltonian.performCI(combinedBasis, nModel, ladderGrid, AsfSettings(); printout=true, writeSummary=false)
    else
        Hamiltonian.performCI(combinedBasis, nModel, ladderGrid, AsfSettings(); printout=true, writeSummary=false)
    end
end
gcDiff = Base.GC_Diff(Base.gc_num(), gcBefore)
@printf(">>> performCI returned %d level(s);  allocations = %d,  bytes = %.3f GB,  GC time = %.1f s\n",
        length(mp.levels), gcDiff.malloc + gcDiff.realloc + gcDiff.poolalloc + gcDiff.bigalloc,
        gcDiff.allocd / 1024^3, gcDiff.total_time / 1e9)
flush(stdout)

if  doProfile
    println("\n>>> CPU PROFILE (flat, self-time, top frames):")
    io = IOBuffer()
    Profile.print(io; format=:flat, sortedby=:count, mincount=200, noisefloor=2.0)
    print(String(take!(io)))
    flush(stdout)
end

sortedLevels = sort(mp.levels, by = lev -> lev.energy)

stage("amplitude loop over all $(length(sortedLevels)) level(s), both gauges") do
    for  (n, lev)  in  enumerate(sortedLevels)
        omega = lev.energy - initialLevel.energy
        for  mp_  in  ppSet.multipoles
            PhotoEmission.amplitude(Absorption(), mp_, Basics.Coulomb,   omega, lev, initialLevel, ladderGrid; printout=false)
            PhotoEmission.amplitude(Absorption(), mp_, Basics.Babushkin, omega, lev, initialLevel, ladderGrid; printout=false)
        end
        if  n % 20 == 0  ||  n == length(sortedLevels)
            @printf("      -> amplitude(s) done for %d of %d levels at %s\n", n, length(sortedLevels),
                    Dates.format(Dates.now(), "HH:MM:SS"))
            flush(stdout)
        end
    end
end

println("\n>>> PROBE COMPLETE.")
