#
# probe-eolMemoryComposition.jl   --   WHERE AN EOL RUN'S PEAK MEMORY ACTUALLY GOES, phase by phase.
#
# THE QUESTION, AND WHY IT IS THE ONLY ONE LEFT WORTH ASKING ABOUT EOL MEMORY.  Two pieces of work aimed at the
# per-CSF-pair angular store: item 29 made it 15x smaller, and item 30 would have held one symmetry block at a
# time instead of all of them.  Item 29 was decisive at 25 085 CSFs -- that row was SIGKILLed before and runs at
# 29.73 GB after -- but item 30 was measured on 14-Sep-2026 to save **1.4 % of the peak for 2.4x the time**
# (3.023 GB against 2.980 GB, 208.8 s against 500.4 s, on 7 062 CSFs, energies bit-identical) and was closed as
# not worth doing.  The two results together say the store is NO LONGER what sets the peak: about 1.4 GB of a
# 29.7 GB peak.  So roughly 95 % of the memory of a large EOL run is UNACCOUNTED FOR, and every further attempt
# to economise is aimed at the wrong term until that is measured.
#
# WHAT THIS FILE DOES.  It walks the phases of an EOL solve in order, in ONE process, reading the kernel's own
# figures at each step:  VmRSS (what is held now) and VmHWM (the high-water mark so far).  The DIFFERENCE between
# successive VmHWM readings attributes the peak to the phase that caused it, which a sum of `summarysize` over
# live structures cannot do -- it misses every temporary, and temporaries are exactly what is suspected here.
#
# READ THE `dHWM` COLUMN.  It is the memory each phase added to the high-water mark;  the phase with the largest
# entry is what a machine must be sized for, and it is NOT necessarily the phase holding the most afterwards.
#
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

function rss()
    r = h = 0.0
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:")  &&  (r = parse(Float64, split(line)[2]) / 1024^2)
        startswith(line, "VmHWM:")  &&  (h = parse(Float64, split(line)[2]) / 1024^2)
    end
    return( (r, h) )
end

lastH = 0.0
function mark(tag)
    GC.gc();   (r, h) = rss()
    @printf("  %-46s  RSS %7.3f GB   HWM %7.3f GB   dHWM %+7.3f GB\n", tag, r, h, h - lastH)
    global lastH = h;   flush(stdout)
end

Z    = 22.0
refs = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s"),
        Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")]
valp = [Shell("3p"), Shell("3d"), Shell("4s")]
to   = length(ARGS) >= 1 && ARGS[1] == "big" ? [Shell("3d"),Shell("4s"),Shell("4p"),Shell("4d")] :
                                               [Shell("3d"),Shell("4s"),Shell("4p")]

println("\n", "="^118);   println("WHERE THE MEMORY GOES -- phase by phase, one EOL solve");   println("="^118)
mark("baseline: JAC loaded")

nm    = Nuclear.Model(Z)
grid  = Basics.recommendedGrid(refs, nm; printout=false)
confs = unique( vcat(refs, Basics.generateConfigurations(
                          Basics.generateConfigurations(refs, valp, to), valp, to)) )
mark("configuration list ($(length(confs)) configurations)")

relconfs = ConfigurationR[]
for c in confs   append!(relconfs, Basics.generateConfigurations(Basics.RelativisticConfigurations(), c))   end
mark("relativistic configurations ($(length(relconfs)))")

subsh = Basics.generateSubshellList(relconfs)
csfs  = CsfR[]
for rc in relconfs   append!(csfs, Basics.generateCsfRs(rc, subsh))   end
mark("CSF list ($(length(csfs)) CSFs, $(length(subsh)) subshells)")

nEl   = sum(csfs[1].occupation)
basis = Basis(true, nEl, subsh, csfs, Subshell[], Dict{Subshell,Orbital}())
blocks = unique( [ LevelSymmetry(c.J, c.parity) for c in csfs ] )
@printf("  ... %d symmetry blocks, largest %d CSFs\n",
        length(blocks), maximum( count(c -> LevelSymmetry(c.J,c.parity) == s, csfs) for s in blocks ))

prim = Bsplines.generatePrimitives(grid)
mark("B-spline primitives")

caches = Dict{LevelSymmetry,Any}()
for sym in blocks   caches[sym] = SC.cacheCsfPairCoefficientsEOL(sym, basis)   end
storeMB = sum( Base.summarysize(c) for (_,c) in caches ) / 1024^2
mark(@sprintf("angular store, ALL blocks (%.0f MB live)", storeMB))

# THE SOLVER'S OWN SETUP, replicated in the order `solveOptimizedLevelField` builds it.  The screened-potential
# tensor caches are the leading suspect: one is built PER RANK that occurs, three per rank (LL, LS, SS), and each
# is a tensor over B-spline pairs -- so an f-shell layer, which reaches rank 6, builds twenty-one of them.
storage = Dict{String,Array{Float64,2}}()
nsL = prim.grid.nsL;   nsS = prim.grid.nsS
matrixB = zeros( nsL+nsS, nsL+nsS )
matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, prim, storage)
matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, prim, storage)
mark("B-spline overlap matrix + TTp storage")

nucPot = Nuclear.nuclearPotential(nm, grid)
posSpec = SC.positiveBranchSpectrum(subsh, prim, nucPot, matrixB, storage)
mark("positive-branch spectrum")

lMax = maximum( Basics.subshell_l(sh) for sh in subsh )
tensorCaches = Dict{Int64,Any}()
for L = 0:(2*lMax)
    cLL = RadialIntegrals.buildScreenedPotentialCache(L, prim.bsplinesL, prim.bsplinesL, grid; rtol=1.0e-6)
    cLS = RadialIntegrals.buildScreenedPotentialCache(L, prim.bsplinesL, prim.bsplinesS, grid; rtol=1.0e-6)
    cSS = RadialIntegrals.buildScreenedPotentialCache(L, prim.bsplinesS, prim.bsplinesS, grid; rtol=1.0e-6)
    tensorCaches[L] = (cLL, cLS, cSS)
    mark(@sprintf("screened-potential tensor cache, rank %d (of 0..%d)", L, 2*lMax))
end
@printf("  ... tensor caches live: %.0f MB\n", Base.summarysize(tensorCaches)/1024^2)

# AND NOW THE PHASES THAT RUN EVERY ITERATION.  All of the setup above reaches only ~1.08 GB where a real solve
# of this same space peaks at 3.02 GB, so about two thirds of the peak is raised INSIDE the iteration -- and the
# leading suspect is `combineAngularCoefficientsEOL`, which pushes every coefficient of every contributing CSF
# pair into one flat vector and condenses it with a nested double loop.  Its elements are `Coefficient2p`, a
# UnionAll, so the vector is NOT concretely typed and every entry is a separately boxed heap object at 105-178
# bytes measured.
alSet = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                   eeInteractionCI = CoulombInteraction(), gridStopper = false)
mpRef = redirect_stdout(devnull) do
            SelfConsistent.performSCF(refs, nm, grid, alSet; printout=false)
        end
mark("reference multiplet (AL on the reference configurations)")
# The correlation subshells have no orbital yet, so a spectrum is generated for them in the screened potential --
# the same construction `Basics.generate(RasExpansion,...)` uses, and a phase of the real run in its own right.
elecPot = Basics.computePotential(Basics.DFSField(1.0), grid, mpRef.levels[1].basis)
meanPot = Basics.add(nucPot, elecPot)
startOrbs = redirect_stdout(devnull) do
                Bsplines.generateOrbitals(subsh, meanPot, nm, prim, printout=false)
            end
orbs  = merge(startOrbs, mpRef.levels[1].basis.orbitals)
basis = Basis(true, nEl, subsh, csfs, Subshell[], orbs)
mark("start orbitals for the correlation subshells")

r1 = Dict{Tuple{Subshell,Subshell},Float64}();  r2 = Dict{Tuple{Int64,Subshell,Subshell,Subshell,Subshell},Float64}()
levels = Level[]
for sym in blocks
    cache = caches[sym]
    m = SC.buildCIMatrixEOL(cache, orbs, grid, nucPot, r1, r2)
    append!(levels, SC.diagonalizeBlockEOL(sym, cache.idxCsf, m, basis))
end
mark("CI matrices built and diagonalised, ALL blocks")

mp   = Basics.sortByEnergy( Multiplet("EOL", levels) )
tgt  = mp.levels[1:min(3,length(mp.levels))]
mark("multiplet sorted")

(c1, c2) = SC.combineAngularCoefficientsEOL(caches, tgt)
mark(@sprintf("combineAngularCoefficientsEOL (returns %d + %d)", length(c1), length(c2)))

genOcc = SC.computeGeneralizedOccupationEOL(caches, tgt, basis)
mark("computeGeneralizedOccupationEOL")

println("\n>> THE PHASE WITH THE LARGEST dHWM IS WHAT THE MACHINE MUST BE SIZED FOR.  A phase can also add nothing")
println(">> to the high-water mark while holding a great deal -- that means an EARLIER phase already reached")
println(">> higher, and the memory it needed was transient.  Those transients are what a `summarysize` of live")
println(">> structures cannot see, and are the leading suspect for the ~95 % that the angular store does not explain.")
