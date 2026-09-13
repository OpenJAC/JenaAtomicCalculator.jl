#
# probe-twoParticleSweepVsGrasp.jl  --  a SWEEP of JAC's two-particle coefficients against GRASP2018, over
#                                       configuration pairs that differ by ONE electron.  Step 1 of priority item 27.
#
# WHY THIS EXISTS.  The only direct GRASP two-particle comparison in the repository (`example-Aq.jl` branch j) runs
# on ONE configuration, 1s^2 2s^2 2p^2, where every CSF pair has the SAME occupations.  The occupation-changing case
# -- one electron moved between subshells -- was validated instead against `SpinAngularGaigalas`, the parked JAC
# predecessor, at more than 128 000 coefficients.  Two implementations with a common ancestor can agree perfectly
# and share an error, and they do: see item 27.  This file closes that gap by putting GRASP on the other side.
#
# SCOPE, deliberately narrow.  Only configuration pairs of TWO electrons outside closed shells, because a CSF is then
# labelled uniquely by (its occupied relativistic subshells, 2J) and the two codes' CSF orders can be matched without
# parsing GRASP's coupling tree.  That is enough to characterise the defect; widening it means a real parser.
#
# USAGE.  Build the oracle as described in `tools/probe-twoParticleVsGrasp.jl`, then
#     JAC_GRASP_BIN=<dir holding angdump2 and gen/rcsfgenerate>  julia --project=. tools/probe-twoParticleSweepVsGrasp.jl
#
using JenaAtomicCalculator, Printf
const SA = JenaAtomicCalculator.SpinAngular

binDir = get(ENV, "JAC_GRASP_BIN", "")
isempty(binDir) && error("Set JAC_GRASP_BIN to the directory holding angdump2 and gen/rcsfgenerate.")
work = mktempdir()

# subshell label  ->  Subshell,  e.g. "3p-" is j = 1/2 (kappa +1), "3p" is j = 3/2 (kappa -2)
function labelToSubshell(lab::String)
    n  = parse(Int64, lab[1:1])
    l  = findfirst(isequal(lab[2]), "spdfghi") - 1
    return( endswith(lab, "-") ? Subshell(n, l) : Subshell(n, -(l+1)) )
end

"""
`graspCsfs(fname::String)`
    ... parses a GRASP `rcsf` file and returns, per CSF in GRASP's own order, the tuple
        `(subshells::Array{Subshell,1}, twoJ::Int64)`. Only CSFs with two occupied subshells are described
        completely by that pair, which is why this file restricts itself to them.
"""
function graspCsfs(fname::String)
    lines = readlines(fname);    i = findfirst(l -> startswith(l, "CSF(s):"), lines) + 1
    out   = Tuple{Array{Subshell,1},Int64}[]
    while  i <= length(lines)
        l = lines[i]
        if  strip(l) == "*"    i += 1;   continue    end
        if  isempty(strip(l))  i += 1;   continue    end
        subs = Subshell[]
        for  m in eachmatch(r"([1-9][spdfghi]-?)\s*\(\s*(\d+)\)", l)
            for  _ = 1:parse(Int64, m.captures[2])   push!(subs, labelToSubshell(String(m.captures[1])))   end
        end
        # the total J and parity sit on the last line of the CSF's block
        j = i + 1
        while  j <= length(lines) && !occursin(r"[0-9]\s*[+-]\s*$", lines[j])    j += 1    end
        j > length(lines)  &&  break
        mm = match(r"(\d+)(?:/(\d+))?\s*([+-])\s*$", strip(lines[j]))
        tw = isnothing(mm.captures[2]) ? 2*parse(Int64, mm.captures[1]) : parse(Int64, mm.captures[1])
        push!(out, (subs, tw))
        i = j + 1
    end
    return( out )
end

function runGrasp(confA::String, confB::String, active::String, tag::String)
    dir = joinpath(work, tag);   mkpath(dir)
    inp = "*\n0\n$(confA)\n$(confB)\n*\n$(active)\n0,6\n0\nn\n"
    open(joinpath(dir, "in.txt"), "w") do io    write(io, inp)    end
    cd(dir) do
        run(pipeline(`$(joinpath(binDir,"gen","rcsfgenerate"))`, stdin=joinpath(dir,"in.txt"),
                     stdout=joinpath(dir,"gen.log"), stderr=devnull))
        cp(joinpath(dir,"rcsf.out"), joinpath(dir,"rcsf.c"), force=true)
        run(pipeline(`$(joinpath(binDir,"angdump2")) 1`, stdout=joinpath(dir,"g.log"), stderr=devnull))
    end
    coeffs = Dict{Tuple{Int,Int}, Vector{NTuple{6,Any}}}()
    for  line in eachline(joinpath(dir,"g.log"))
        p = split(line);    length(p) == 8 || continue
        all(x -> occursin(r"^-?\d+$", x), p[1:7]) || continue
        ic, ir = parse(Int64,p[1]), parse(Int64,p[2])
        push!(get!(coeffs,(ic,ir),[]), (parse(Int64,p[7]), parse(Int64,p[3]), parse(Int64,p[4]),
                                        parse(Int64,p[5]), parse(Int64,p[6]), parse(Float64,p[8])))
    end
    return( (graspCsfs(joinpath(dir,"rcsf.out")), coeffs) )
end

canon(k,a,b,c,d) = (k, minimum([(a,b,c,d),(b,a,d,c),(c,d,a,b),(d,c,b,a)])...)
function bag(list)
    m = Dict{NTuple{5,Int64},Float64}()
    for (k,a,b,c,d,v) in list    key = canon(k,a,b,c,d);   m[key] = get(m,key,0.0) + v    end
    return( Dict(kk => vv  for (kk,vv) in m  if abs(vv) > 1.0e-12) )
end

cases = [ ("3s(1,i)3p(1,i)", "3p(1,i)3d(1,i)", "3s,3p,3d",    "3s3p / 3p3d"),
          ("3s(1,i)3p(1,i)", "3s(1,i)4p(1,i)", "3s,3p,4p",    "3s3p / 3s4p"),
          ("3s(1,i)3p(1,i)", "3p(1,i)4s(1,i)", "3s,3p,4s",    "3s3p / 3p4s"),
          ("2s(1,i)2p(1,i)", "2p(1,i)3d(1,i)", "2s,2p,3d",    "2s2p / 2p3d"),
          ("2s(1,i)2p(1,i)", "2s(1,i)3p(1,i)", "2s,2p,3p",    "2s2p / 2s3p"),
          ("3p(1,i)3d(1,i)", "3d(1,i)4f(1,i)", "3p,3d,4f",    "3p3d / 3d4f"),
          ("3s(1,i)4f(1,i)", "4f(1,i)5g(1,i)", "3s,4f,5g",    "3s4f / 4f5g") ]

@printf("\n%-16s %7s %8s %9s %8s %7s   %s\n", "configurations", "pairs", "agreeing", "differing",
        "missing", "extra", "verdict")
for  (cA, cB, act, tag)  in  cases
    # A case whose CSF list GRASP declines to build is SKIPPED and SAID SO, never silently dropped: the sweep is
    # evidence, and a case missing from the table without a reason looks like a case that passed.
    local gCsfs, gCoef
    try
        (gCsfs, gCoef) = runGrasp(cA, cB, act, replace(tag, r"[ /]" => "_"))
    catch err
        @printf("%-16s %7s %8s %9s %8s %7s   rcsfgenerate declined this input\n", tag, "-","-","-","-","-")
        continue
    end
    # JAC's CSFs for the same two configurations, labelled the same way
    relconfs = ConfigurationR[]
    for  c in [replace(cA, r"\((\d+),i\)" => s" "), replace(cB, r"\((\d+),i\)" => s" ")]
        append!(relconfs, Basics.generateConfigurations(Basics.RelativisticConfigurations(), Configuration(c)))
    end
    subsh = Basics.generateSubshellList(relconfs)
    jCsfs = CsfR[];    for rc in relconfs   append!(jCsfs, Basics.generateCsfRs(rc, subsh))   end
    idOf  = c -> (sort([ subsh[i] for i = 1:length(subsh) if c.occupation[i] > 0 ], by = s -> (s.n, s.kappa)),
                  Basics.twice(c.J))
    gId   = [ (sort(g[1], by = s -> (s.n, s.kappa)), g[2])  for g in gCsfs ]
    jToG  = Dict{Int,Int}()
    for  (ij, c) in enumerate(jCsfs)
        k = findfirst(isequal(idOf(c)), gId);    isnothing(k) || (jToG[ij] = k)
    end
    op = SA.TwoParticleOperator(0, Basics.plus)
    six = Dict(sh => i for (i,sh) in enumerate(subsh))
    jCoef = Dict{Tuple{Int,Int}, Vector{NTuple{6,Any}}}()
    for  (ii,ci) in enumerate(jCsfs), (jj,cj) in enumerate(jCsfs)
        (haskey(jToG,ii) && haskey(jToG,jj) && ci.J == cj.J && ci.parity == cj.parity) || continue
        for  cf in SA.computeCoefficients(op, ci, cj, subsh)
            g = SA.toGraspCoulomb(cf);    abs(g.V) > 1.0e-14 || continue
            push!(get!(jCoef,(jToG[ii],jToG[jj]),[]), (g.nu, six[g.a], six[g.b], six[g.c], six[g.d], g.V))
        end
    end
    nPair = 0; nOk = 0; nDif = 0; nMiss = 0; nExtra = 0
    for  key in union(keys(gCoef), keys(jCoef))
        bg = bag(get(gCoef,key,[]));    bj = bag(get(jCoef,key,[]))
        nPair += 1;    bad = false
        for  (kk,v) in bg
            if      !haskey(bj,kk)                   nMiss  += 1;  bad = true
            elseif  abs(v - bj[kk]) > 1.0e-9         nDif   += 1;  bad = true
            end
        end
        for  kk in keys(bj)    haskey(bg,kk) || (nExtra += 1;  bad = true)    end
        bad ? nothing : (nOk += 1)
    end
    @printf("%-16s %7d %8d %9d %8d %7d   %s\n", tag, nPair, nOk, nDif, nMiss, nExtra,
            (nDif+nMiss+nExtra) == 0 ? "agree" : "DISAGREE")
    flush(stdout)
end
