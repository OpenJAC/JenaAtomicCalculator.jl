#
# diag-grasp-angular.jl
#
# A BILATERAL REGRESSION NET between JAC's spin-angular coefficients and GRASP2018's.
#
# The point is not a single comparison. Angular coefficients are notoriously hard for general open-shell
# configurations and general operators, and a defect may sit on EITHER side; if one of the two codes is later
# corrected for a special case, we want to SEE it rather than rediscover it. So this file turns "compare a new
# configuration" from an afternoon into one call, and prints reference tables ready to paste into an example.
#
# USAGE, from the JAC root:
#     include("work/diag-grasp-angular.jl")
#     g = GraspAngular.build()                                  # copies + builds GRASP2018 and the two drivers
#     GraspAngular.compare1p(g, ["1s^2 2s^2 2p^2"])             # one-particle, all ranks
#     GraspAngular.compare2p(g, ["1s^2 2s^2 2p^2"])             # two-particle (Coulomb)
#     GraspAngular.reference2p(g, ["1s^2 3d^2"])                # emit a paste-ready reference table
#
# WHAT IT NEEDS: gfortran, and a readable GRASP2018 source tree. Nothing is written into that tree -- it is copied
# first -- so a read-only reference installation is fine and Rule 6 is respected.
#
# TRAPS THIS FILE ALREADY HANDLES, each of which cost time to find once:
#   * the libraries must be built in the order libmod -> lib9290 -> libmcp90 -> librang90;
#   * gfortran 13 needs -std=legacy -fallow-argument-mismatch, or librang90 will not compile;
#   * only SETMC wants LAPACK, for four DLAMCH calls, which dlamch_shim.f90 supplies;
#   * rcsfgenerate's dialogue puts the 2*J range BEFORE the excitation count, and 2*J must have the parity of the
#     electron number -- odd 2*J for an odd number of electrons, or the generator dies with an EOF;
#   * a CSF list with fixed occupations yields NO substitutions; two reference configurations are needed for those.
#
module GraspAngular

using JenaAtomicCalculator, Printf

const DEFAULT_SOURCE = "/home/fritzsch/fri/grasp/grasp2018/grasp-master"
const FFLAGS = ["-O2", "-fno-automatic", "-std=legacy", "-fallow-argument-mismatch"]


"""
`GraspAngular.build(; source, workDir)`
    ... to copy a GRASP2018 source tree into `workDir` and build the libraries, `rcsfgenerate` and the two dump
        drivers. Nothing is written into `source`. Returns the run directory, which every other call takes.
"""
function build(; source::String = DEFAULT_SOURCE,
                 workDir::String = joinpath(tempdir(), "grasp-angular"))
    isdir(source) || error("GraspAngular.build: no GRASP source at $source")
    root = joinpath(workDir, "grasp-master")
    if !isdir(root)
        mkpath(workDir);   run(`cp -r $source $workDir/`)
    end
    lib = joinpath(root, "lib");    bin = joinpath(root, "bin")
    mkpath(lib);   mkpath(bin)
    env = copy(ENV)
    env["FC"] = "gfortran";  env["FC_FLAGS"] = join(FFLAGS, " ");  env["FC_LD"] = " "
    env["GRASP"] = root;     env["LAPACK_LIBS"] = ""
    for d in ["libmod", "lib9290", "libmcp90", "librang90"]          # ORDER MATTERS
        isfile(joinpath(lib, "lib" * (d == "libmod" ? "mod" : d[4:end]) * ".a")) && continue
        # stdout is noise, stderr is NOT suppressed: a build tool that hides its own errors is useless
        # exactly when it breaks, which cost time here once already.
        cd(joinpath(root, "src", "lib", d)) do
            run(pipeline(setenv(`make`, env), stdout = devnull))
        end
    end
    if !isfile(joinpath(bin, "rcsfgenerate"))
        cd(joinpath(root, "src", "appl", "rcsfgenerate90")) do
            run(pipeline(setenv(`make`, env), stdout = devnull))
        end
    end
    run_ = joinpath(workDir, "run");   mkpath(run_)
    drivers = joinpath(@__DIR__, "grasp-drivers")
    incs = ["-I" * joinpath(root, "src", "lib", d) for d in ["libmod","lib9290","librang90","libmcp90"]]
    for (src, exe) in [("angdump1.f90","angdump1"), ("angdump2.f90","angdump2")]
        isfile(joinpath(run_, exe)) && continue
        cd(run_) do
            run(`gfortran $FFLAGS -o $exe $(joinpath(drivers,src)) $(joinpath(drivers,"dlamch_shim.f90"))
                 $incs -L$lib -lrang90 -l9290 -lmcp90 -lmod`)
        end
    end

    return( run_ )
end


"""
`GraspAngular.csfs(runDir, configs; twoJ, excitations)`
    ... to drive `rcsfgenerate` for the given spectroscopic configurations and leave `rcsf.c` in `runDir`. The
        configurations are GRASP strings such as "1s(2,i)2s(2,i)2p(2,i)". Returns the number of CSFs generated.
"""
function csfs(runDir::String, configs::Array{String,1}; active::String = "", twoJ::Tuple{Int,Int} = (0,4),
              excitations::Int = 0)
    root = joinpath(dirname(runDir), "grasp-master")
    act  = active == "" ? "2s,2p" : active
    dial = "*\n0\n" * join(configs, "\n") * "\n*\n" * act * "\n" *
           "$(twoJ[1]),$(twoJ[2])\n$(excitations)\nn\n"
    cd(runDir) do
        rm("rcsf.out", force = true);   rm("rcsf.c", force = true)
        open("dialogue.txt","w") do io   write(io, dial)   end
        try  run(pipeline(`$(joinpath(root,"bin","rcsfgenerate"))`, stdin = "dialogue.txt",
                          stdout = devnull, stderr = devnull))
        catch;  end
        isfile("rcsf.out") || error("GraspAngular.csfs: rcsfgenerate produced nothing -- check the 2*J parity " *
                                    "(odd 2*J for an odd electron number) and the active set")
        cp("rcsf.out", "rcsf.c", force = true)
    end
    # ... the CSF count is read from the driver rather than guessed from the CSL layout, which is easy to
    #     miscount: a block spans several lines and the leading subshell differs from case to case.
    n = 0
    out = cd(() -> read(`./angdump1`, String), runDir)
    for line in split(out, '\n')
        if occursin("# NCF =", line)
            n = parse(Int, split(split(line, "NCF =")[2])[1]);   break
        end
    end

    return( n )
end


"""
`GraspAngular.dump1p(runDir)`
    ... to run the one-particle driver on the CSF list in `runDir`; returns a Dict keyed by
        (rank, icsf, jcsf, ia, ib) with the GRASP coefficient as value.
"""
function dump1p(runDir::String)
    out = cd(() -> read(`./angdump1`, String), runDir)
    d = Dict{NTuple{5,Int},Float64}()
    for line in split(out, '\n')
        p = split(strip(line))
        length(p) == 8 && all(!isnothing, tryparse.(Int, p[1:5])) || continue
        d[(parse(Int,p[1]),parse(Int,p[2]),parse(Int,p[3]),parse(Int,p[4]),parse(Int,p[5]))] = parse(Float64,p[8])
    end

    return( d )
end


"""
`GraspAngular.dump2p(runDir; incor)`
    ... to run the two-particle driver on the CSF list in `runDir`; returns a Dict keyed by (icsf, jcsf, canonical
        R^k key) with the GRASP Coulomb coefficient as value, duplicates summed. The canonical key folds the
        symmetries of R^k(abcd), without which the two codes' differing orderings manufacture false differences.
"""
function dump2p(runDir::String; incor::Int = 1)
    out = cd(() -> read(`./angdump2 $incor`, String), runDir)
    d = Dict{Any,Float64}()
    for line in split(out, '\n')
        p = split(strip(line))
        length(p) == 8 && all(!isnothing, tryparse.(Int, p[1:7])) || continue
        ic,ir,a,b,c,dd,k = (parse(Int,p[i]) for i in 1:7)
        key = (ic, ir, rkKey(k,a,b,c,dd))
        d[key] = get(d,key,0.0) + parse(Float64,p[8])
    end

    return( d )
end


"""
`GraspAngular.rkKey(k, a, b, c, d)`
    ... the canonical key of R^k(abcd). Since a and c both belong to electron 1 and b and d both to electron 2,
        the integral is invariant under a <-> c, under b <-> d, and under exchanging the two electrons -- an orbit
        of EIGHT index patterns, not four.

        GETTING THIS WRONG MANUFACTURES A DISAGREEMENT. With only the four patterns of the electron swap, terms
        that GRASP writes as (c,b,a,d) and JAC as (a,b,c,d) land on different keys, and the same physical integral
        is then counted as "missing on one side" and "disagreeing on the other". Returns a tuple usable as a key.
"""
function rkKey(k::Int, a::Int, b::Int, c::Int, d::Int)
    # a and c both sit on electron 1, b and d both on electron 2, so R^k is symmetric under a<->c, under
    # b<->d, AND under the electron swap (ab)<->(cd). That is an EIGHT-element orbit, not four.
    orbit = [(a,b,c,d), (c,b,a,d), (a,d,c,b), (c,d,a,b),
             (b,a,d,c), (d,a,b,c), (b,c,d,a), (d,c,b,a)]

    return( (k, minimum(orbit)) )
end



"""
`GraspAngular.signatures(runDir)`
    ... to read, from the one-particle driver's output, a signature for every GRASP CSF: its occupation vector, its
        2J and its parity. The occupations come from the RANK-0 DIAGONAL coefficients, which in GRASP's convention
        ARE the occupation numbers, so no parsing of the CSL file's coupling tree is needed. Returns
        (subshells::Array{Tuple{Int,Int},1} as (n, kappa), sigs::Array{Tuple{Array{Int,1},Int,Int},1}).
"""
function signatures(runDir::String)
    out = cd(() -> read(`./angdump1`, String), runDir)
    ncf = 0;   nw = 0;   sub = Tuple{Int,Int}[];   twoJ = Int[];   par = Int[]
    occ = Dict{Int,Array{Int,1}}()
    for line in split(out, '\n')
        if      occursin("# NCF =", line)
            ncf = parse(Int, split(split(line,"NCF =")[2])[1]);   nw = parse(Int, split(split(line,"NW =")[2])[1])
        elseif  occursin("# subshell", line)
            p = split(line);   push!(sub, (parse(Int,p[6]), parse(Int,p[9])))
        elseif  occursin("# csf", line)
            p = split(line);   push!(twoJ, parse(Int,p[6]) - 1);   push!(par, parse(Int,p[9]))
        end
    end
    for i = 1:ncf    occ[i] = zeros(Int, nw)    end
    for line in split(out, '\n')
        p = split(strip(line))
        length(p) == 8 && all(!isnothing, tryparse.(Int, p[1:5])) || continue
        parse(Int,p[1]) == 0 || continue                       # rank 0 only
        ic = parse(Int,p[2]);  ir = parse(Int,p[3]);  ia = parse(Int,p[4]);  ib = parse(Int,p[5])
        ic == ir && ia == ib || continue                        # diagonal only
        occ[ic][ia] = round(Int, parse(Float64,p[8]))
    end

    return( sub, [(occ[i], twoJ[i], par[i]) for i = 1:ncf] )
end


"""
`GraspAngular.fingerprints(runDir)`
    ... to give every GRASP CSF a fingerprint that distinguishes CSFs of the SAME occupation, 2J and parity but
        DIFFERENT coupling. The fingerprint is the sorted vector of that CSF's diagonal rank-1..3 one-particle
        coefficients, which depend on the coupling tree and so separate what (occupation, 2J, parity) cannot.

        THIS BOOTSTRAPS ON THE ONE-PARTICLE RESULT and that should be stated: JAC and GRASP agree on rank > 0
        one-particle coefficients to ~1e-15 over 78 checked values, so the same quantity computed on both sides is
        a legitimate label. It is an assumption, not a theorem; if that agreement ever broke, this matching would
        break with it. Returns Array{Array{Float64,1},1}, indexed by GRASP CSF.
"""
function fingerprints(runDir::String)
    out = cd(() -> read(`./angdump1`, String), runDir)
    ncf = 0
    for line in split(out, '\n')
        occursin("# NCF =", line) && (ncf = parse(Int, split(split(line,"NCF =")[2])[1]))
    end
    fp = [Float64[] for _ = 1:ncf]
    for line in split(out, '\n')
        p = split(strip(line))
        length(p) == 8 && all(!isnothing, tryparse.(Int, p[1:5])) || continue
        k = parse(Int,p[1]);  ic = parse(Int,p[2]);  ir = parse(Int,p[3])
        k > 0 && ic == ir || continue
        push!(fp[ic], round(parse(Float64,p[8]), digits = 9))
    end

    return( [sort(v) for v in fp] )
end


"""
`GraspAngular.compare2p(runDir, jacConfig; incor)`
    ... to compare JAC's two-particle coefficients with GRASP's for the given JAC configuration, matching the two
        CSF orderings automatically and converting JAC's X^L coefficients into GRASP's Coulomb convention. Prints
        a report and returns a NamedTuple of the counts.

        THE REPORT SAYS WHAT WAS COMPARED AND WHAT WAS NOT, which matters: a GRASP coefficient belonging to a CSF
        pair that could not be matched is NOT compared, and reporting "0 disagreements" over the matched subset
        while quoting GRASP's full total beside it reads as a far stronger claim than was made. Both are printed,
        and a partial run says so on its own line.

        The conversion is  V_grasp = V_jac * (-1)^k * <a||C^k||c> * <b||C^k||d>, the prefactor of
        `InteractionStrength.XL_CoulombReference`. Coefficients are compared as a multiset on the canonical R^k key.
"""
function compare2p(runDir::String, jacConfig::String; incor::Int = 1)
    gsub, gsig = signatures(runDir)
    gfp        = fingerprints(runDir)
    grasp      = dump2p(runDir, incor = incor)

    rel = Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(jacConfig))
    sub = Basics.generateSubshellList(rel)
    Defaults.setDefaults("relativistic subshell list", sub; printout = false)
    cs  = CsfR[];   for r in rel   append!(cs, Basics.generateCsfRs(r, sub))   end

    jsub = [(sh.n, sh.kappa) for sh in sub]
    if  jsub != gsub
        println("  SUBSHELL ORDERINGS DIFFER -- JAC $jsub  vs GRASP $gsub");   return( nothing )
    end

    function jacFingerprint(c)
        v = Float64[]
        for k = 1:3
            for x in SpinAngular.computeCoefficients(SpinAngular.OneParticleOperator(k, Basics.plus, true), c, c, sub)
                abs(x.T) > 1.0e-14 && push!(v, round(x.T, digits = 9))
            end
        end
        sort(v)
    end
    jacSig(c) = (collect(c.occupation), Basics.twice(c.J), c.parity == Basics.plus ? 1 : -1)

    jacToGrasp = Dict{Int,Int}();   unmatched = Int[]
    for (i,c) in enumerate(cs)
        hits = findall(==(jacSig(c)), gsig)
        if      length(hits) == 1     jacToGrasp[i] = hits[1]
        elseif  length(hits) >  1
            fp = jacFingerprint(c)
            refined = [h for h in hits if gfp[h] == fp]
            length(refined) == 1 ? (jacToGrasp[i] = refined[1]) : push!(unmatched, i)
        else    push!(unmatched, i)
        end
    end
    matchedG = Set(values(jacToGrasp))

    idx  = Dict(sh => i for (i,sh) in enumerate(sub))
    op   = SpinAngular.TwoParticleOperator(0, Basics.plus, true)
    conv = Dict{Any,Float64}();   nRaw = 0;   nAnn = 0
    for (ic,l) in enumerate(cs), (ir,r) in enumerate(cs)
        (haskey(jacToGrasp,ic) && haskey(jacToGrasp,ir)) || continue
        for c in SpinAngular.computeCoefficients(op, l, r, sub)
            abs(c.V) < 1.0e-14 && continue
            nRaw += 1
            f = AngularMomentum.CL_reduced_me(c.a, c.nu, c.c) * AngularMomentum.CL_reduced_me(c.b, c.nu, c.d)
            isodd(c.nu) && (f = -f)
            abs(f) < 1.0e-14 && (nAnn += 1)
            kk = (jacToGrasp[ic], jacToGrasp[ir], rkKey(c.nu, idx[c.a], idx[c.b], idx[c.c], idx[c.d]))
            conv[kk] = get(conv,kk,0.0) + c.V*f
        end
    end
    surv = filter(p -> abs(p[2]) > 1.0e-10, conv)

    inScope  = filter(p -> p[1][1] in matchedG && p[1][2] in matchedG, grasp)
    outScope = length(grasp) - length(inScope)
    nMatched = 0;  nDisagree = 0;  worst = 1.0;  nMissing = 0
    for (kk,g) in inScope
        if  haskey(surv, kk)
            v = surv[kk];  nMatched += 1
            abs(v/g - 1.0) > 1.0e-9 && (nDisagree += 1)
            abs(v/g - 1.0) > abs(worst - 1.0) && (worst = v/g)
        else
            nMissing += 1
        end
    end
    nExtra = length(surv) - nMatched

    println("  ", jacConfig)
    @printf("    CSFs                : JAC %d, GRASP %d;  matched %d;  UNMATCHED %d\n",
            length(cs), length(gsig), length(jacToGrasp), length(unmatched))
    @printf("    CSF pairs compared  : %d of %d\n", length(jacToGrasp)^2, length(cs)^2)
    @printf("    GRASP coefficients  : %d inside those pairs, %d OUTSIDE and NOT compared\n",
            length(inScope), outScope)
    @printf("    JAC coefficients    : %d raw, %d annihilated by C^k, %d surviving\n", nRaw, nAnn, length(surv))
    @printf("    VERDICT             : %d matched, %d disagreeing, %d missing, %d extra;  worst ratio %.15f\n",
            nMatched, nDisagree, nMissing, nExtra, worst)
    outScope > 0 && @printf("    NOTE: PARTIAL RUN -- %d GRASP coefficients were never looked at.\n", outScope)

    return( (csfs = length(cs), matchedCsfs = length(jacToGrasp), unmatchedCsfs = length(unmatched),
             inScope = length(inScope), outScope = outScope, matched = nMatched, disagree = nDisagree,
             missings = nMissing, extra = nExtra, worst = worst) )
end


"""
`GraspAngular.provenance(runDir)`
    ... to describe WHICH code produced the numbers in `runDir`: the source tree it was copied from, its git
        revision if it has one, the compiler, and whether the angular entry points this file depends on are the
        ones it expects. Returns a Dict, and `GraspAngular.stamp` renders it as comment lines for an example file.

        THIS EXISTS BECAUSE A REFERENCE TABLE WITHOUT A PROVENANCE IS AN ORPHAN. GRASP is not one thing: 92, 2K,
        2013, 2018 and GRASPG differ, and librang90's own ReadMe opens by saying there are errors in the original
        Grasp2K which it corrects. If the angular generation changes on either side -- ours or theirs -- the only
        way to see it is to compare apples with apples, and that requires knowing which apple. A number recorded
        without its source cannot later be re-attributed.

        WHAT HAPPENS IF THE OTHER SIDE CHANGES. Two cases, and they fail differently. If an ENTRY POINT changes
        name or signature, the drivers stop compiling -- loud, immediate, unmistakable. If the NUMBERS change while
        the interface does not, nothing breaks: the comparison simply starts disagreeing, which is the case this
        whole file exists to catch. The `entryPoints` check below covers the first; the comparison covers the
        second; the provenance is what lets you say which version each belongs to.
"""
function provenance(runDir::String)
    root = joinpath(dirname(runDir), "grasp-master")
    d    = Dict{String,Any}()
    d["graspRoot"]  = root
    d["builtFrom"]  = get(ENV, "GRASP_SOURCE", DEFAULT_SOURCE)
    rev = "none"
    # ... a reference installation is usually not a git checkout; that is normal, not an error
    try   rev = strip(read(pipeline(`git -C $root rev-parse --short HEAD`, stderr = devnull), String))
    catch;  end
    d["graspRevision"] = rev
    ver = "unknown"
    try   ver = split(read(`gfortran --version`, String), '\n')[1]   catch;  end
    d["compiler"] = ver
    d["flags"]    = join(FFLAGS, " ")
    # ... the entry points this file depends on; a rename here is what a new GRASP generation would most likely change
    wanted = ["onescalar.f90", "oneparticlejj.f90", "rkco_gg.f90", "speak.f90", "buffer_C.f90"]
    found  = String[];   missing = String[]
    for w in wanted
        hits = String[]
        for (rt, _, fs) in walkdir(joinpath(root, "src"))
            w in fs && push!(hits, joinpath(rt, w))
        end
        isempty(hits) ? push!(missing, w) : push!(found, w)
    end
    d["entryPointsFound"]   = found
    d["entryPointsMissing"] = missing

    return( d )
end


"""
`GraspAngular.stamp(runDir)`
    ... to render the provenance as comment lines, for pasting above a reference table in an example file so that
        the numbers carry the identity of the code that produced them. Returns a String.
"""
function stamp(runDir::String)
    d = provenance(runDir)
    io = IOBuffer()
    println(io, "    #   REFERENCE PROVENANCE -- which code produced the numbers below:")
    println(io, "    #     source tree : ", d["builtFrom"])
    println(io, "    #     git revision: ", d["graspRevision"])
    println(io, "    #     compiler    : ", d["compiler"])
    println(io, "    #     flags       : ", d["flags"])
    if  isempty(d["entryPointsMissing"])
        println(io, "    #     entry points: all expected present (onescalar, oneparticlejj, rkco_gg, speak, buffer_C)")
    else
        println(io, "    #     entry points: MISSING ", join(d["entryPointsMissing"], ", "),
                    "  -- a different GRASP generation?")
    end

    return( String(take!(io)) )
end

end # module
