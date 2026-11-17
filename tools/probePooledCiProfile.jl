#
# probePooledCiProfile.jl -- a STATISTICAL profile of the pooled 342-CSF Fe K-LL `performCI` (item 17), taken
# WHILE it runs rather than after it finishes.
#
# The 16-Nov-2026 retry's kill-time backtrace showed the signalled thread inside `RadialIntegrals.GrantIab` at
# the moment `timeout` fired -- but that is ONE stack sample, and `GrantIab` already got a subshell-pair-keyed
# memo on 09-Nov (`8ff45d2d`, `i1Cache` in `Hamiltonian.setupMatrix`). A single sample cannot say whether
# radial-integral MISSES still dominate the wall time or whether `SpinAngular.computeCoefficients` -- run once
# per CSF pair, ~342*341/2 = 58311 times, with NO memo at all -- is the real cost. This script gets the answer
# the honest way: `@profile` the real call and print a snapshot from a SEPARATE thread after a short, fixed
# wall-clock window, without waiting for completion.
#
# Run from the JAC root:   julia --project=. -t 8 tools/probePooledCiProfile.jl
# Env: PROBE_SNAPSHOT_S (default 420) -- seconds to let performCI run before the snapshot is taken and the
#      process exits. Keep this well under the 3 h that killed the full run; the profile's shape should already
#      be stable long before that.
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra, Distributed, Profile, Dates

setDefaults("unit: energy", "eV")

snapshotS = parse(Float64, get(ENV, "PROBE_SNAPSHOT_S", "420"))

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

println("probePooledCiProfile.jl: same basis as probe-pooledPseudostateTiming.jl, profiled mid-flight")
println(">>> snapshot after $(snapshotS) s of performCI, Threads.nthreads() = $(Threads.nthreads())")
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

ppSet = PhotoIonizationPseudostates.Settings(PhotoIonizationPseudostates.Settings(); multipoles=[E1],
            gauges=[UseCoulomb, UseBabushkin], pseudostateWindow=(171.0, 171.6), pseudostateRbox=17000.0,
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

println(">>> launching performCI under @profile on a background task, snapshot in $(snapshotS) s ...")
flush(stdout)

Profile.clear()
Profile.init(n = 10_000_000, delay = 0.005)

ciTask = Threads.@spawn begin
    @profile Hamiltonian.performCI(combinedBasis, nModel, ladderGrid, AsfSettings(); printout=false, writeSummary=false)
end

sleep(snapshotS)

nsamples = length(Profile.fetch())
println("\n>>> SNAPSHOT at t=$(snapshotS)s, $(nsamples) raw sample slots collected, ciTask done = $(istaskdone(ciTask))")
flush(stdout)

println("\n>>> FLAT PROFILE, self-time, top frames (mincount filters noise):")
io = IOBuffer()
Profile.print(io; format=:flat, sortedby=:count, mincount=50, noisefloor=2.0)
print(String(take!(io)))
flush(stdout)

println("\n>>> PROBE COMPLETE -- exiting without waiting for performCI to finish.")
exit(0)
