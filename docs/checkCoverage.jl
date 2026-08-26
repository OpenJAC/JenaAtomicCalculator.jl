#
# checkCoverage.jl -- does the published documentation still cover what Rule 16 of CLAUDE.md says it should?
#
# Run from the JAC root:   julia --project=. docs/checkCoverage.jl
#
# This exists because the failure it detects is SILENT. The `@autodocs` blocks in docs/src/api-*.md select
# content either by module (the whole module) or by a `Pages` allowlist, and Documenter says nothing about a
# source file that no block asked for. Until 18-Aug-2026 whole modules were missing from the published API and
# nothing reported it. The script re-derives the sets rather than trusting a list, so it cannot go stale.
#
# It PRINTS and returns a status; it changes nothing.

using Printf

const ROOT     = normpath(joinpath(@__DIR__, ".."))
const SRC      = joinpath(ROOT, "src")
const DOCSRC   = joinpath(ROOT, "docs", "src")
const EXAMPLES = joinpath(ROOT, "examples")

"""
`stripComments(text::String)`
    ... removes everything after a `#` on each line, so that a module named only in a comment does not count as
        used. A code::String is returned.
"""
function stripComments(text::String)
    return( join([split(l, "#")[1] for l in split(text, "\n")], "\n") )
end


"""
`datedBranchCount(path::String)`
    ... counts the `# Last successful:` markers of one example file that carry a real date, i.e. excluding the
        `unknown ...` placeholder. A count::Int64 is returned.
"""
function datedBranchCount(path::String)
    n = 0
    for m in eachmatch(r"#\s*Last successful:\s*([^\n]*)", read(path, String))
        d = strip(m.captures[1])
        if  !isempty(d) && !startswith(lowercase(d), "unknown")    n += 1    end
    end

    return( n )
end


"""
`qualifyingModules()`
    ... determines which physics modules Rule 16 would publish: those USED IN THE CODE of an example file that
        carries at least two dated branches. A set::Set{String} of module names is returned.
"""
function qualifyingModules()
    mods = Set{String}()
    for f in readdir(SRC)
        if startswith(f, "module-") && endswith(f, ".jl") && !occursin("-inc-", f)
            push!(mods, f[8:end-3])
        end
    end
    good = Set{String}()
    for f in readdir(EXAMPLES)
        (startswith(f, "example-") && endswith(f, ".jl"))  ||  continue
        path = joinpath(EXAMPLES, f)
        datedBranchCount(path) >= 2  ||  continue
        code = stripComments(read(path, String))
        for m in mods
            occursin(Regex("\\b" * m * "\\."), code)  &&  push!(good, m)
        end
    end

    return( good )
end


"""
`documentedModules()`
    ... collects every module named in an `@autodocs` block of `docs/src/*.md`, together with the union of the
        source files those blocks select. A tuple (modules::Set{String}, pages::Set{String}) is returned.
"""
function documentedModules()
    mods = Set{String}();    pages = Set{String}()
    for f in readdir(DOCSRC)
        endswith(f, ".md")  ||  continue
        text = read(joinpath(DOCSRC, f), String)
        for blk in eachmatch(r"```@autodocs(.*?)```"s, text)
            b = blk.captures[1]
            for m in eachmatch(r"Modules\s*=\s*\[(.*?)\]"s, b),  name in split(m.captures[1], ",")
                nm = strip(String(name))
                isempty(nm)  ||  push!(mods, String(last(split(nm, "."))))
            end
            for m in eachmatch(r"Pages\s*=\s*\[(.*?)\]"s, b),  q in eachmatch(r"\"([^\"]+)\"", m.captures[1])
                push!(pages, q.captures[1])
            end
        end
    end

    return( (mods, pages) )
end


## Modules the MAINTAINER has decided not to publish yet, each with its reason and the condition that releases it.
## This is NOT a way of silencing the check: a withheld module is PRINTED on every run, and the reason is carried here
## rather than in anyone's memory, so the decision stays visible and can be revisited. Adding an entry to make a message
## go away, rather than because the maintainer decided it, is the misuse this comment exists to name.
const WITHHELD = Dict(
    "SpinAngularNew" => "maintainer, 26-Aug-2026: not published before the sqrt(2j+1) convention migration (Rule 19) " *
                        "is done, and then only the functions other modules actually call -- today nothing in src/ " *
                        "calls it, so the published set would be empty in any case." )

qual              = qualifyingModules()
(documented, pgs) = documentedModules()
withheld          = sort(collect(intersect(qual, keys(WITHHELD))))
missing           = sort(collect(setdiff(qual, documented, keys(WITHHELD))))

println("\n", "="^96)
println("  Documentation coverage against Rule 16 of CLAUDE.md")
println("="^96)
@printf("  physics modules qualifying (>=2 dated example branches) : %3d\n", length(qual))
@printf("  modules referenced by docs/src/api-*.md                 : %3d\n", length(documented))
@printf("  QUALIFYING BUT NOT DOCUMENTED                           : %3d\n", length(missing))
for m in missing    println("       * ", m)    end
@printf("  deliberately WITHHELD by the maintainer                  : %3d\n", length(withheld))
for m in withheld   println("       * ", m, "\n           ", WITHHELD[m])   end

## A module documented through a Pages allowlist must have ALL of its files listed, or the rest vanish silently.
holes = String[]
for m in sort(collect(intersect(documented, qual)))
    files = filter(f -> startswith(f, "module-" * m * ".jl") || startswith(f, "module-" * m * "-inc-"), readdir(SRC))
    isempty(files)  &&  continue
    ## only a Pages-limited module can lose files; one selected whole is complete by construction
    if  any(f -> f in pgs, files)
        lost = filter(f -> !(f in pgs), files)
        isempty(lost)  ||  append!(holes, lost)
    end
end
@printf("\n  files of a Pages-limited module NOT listed anywhere      : %3d\n", length(holes))
for f in sort(holes)    println("       * ", f)    end

ok = isempty(missing) && isempty(holes)
println("\n  ", ok ? "OK -- the published set matches Rule 16." :
                     "MISMATCH -- see above; re-derive the sets, do not hand-edit an allowlist.")
println("="^96, "\n")
