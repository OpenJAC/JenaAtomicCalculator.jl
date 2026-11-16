#
# measureDevelopmentRecord.jl -- re-derive the AI-assistance development-record numbers quoted in
# docs/src/development.md.
#
#   Run from the JAC root:   julia --project=. tools/measureDevelopmentRecord.jl
#
# WHY THIS EXISTS.  Priority item 24: the development-record numbers (how many commits carry an AI
# co-author trailer, and since when) were typed by hand on 15-Nov-2026 from a one-off count, and nothing
# re-derives them -- this project's most repeated failure mode (Rule 7 exists because of it, and three
# separate hand-typed figures went stale in one week back in August for the same reason). This script
# re-derives the numbers from `git log` on the current checkout alone -- no stored count, no manual tally --
# and rewrites ONLY the region between the markers in the page:
#
#     <!-- BEGIN generated: development record -->  ...  <!-- END generated: development record -->
#
# Everything outside the markers is prose and stays under human control. The script PRINTS the table and,
# with --write, replaces the region; it changes nothing else. It also stamps the commit and the date into
# the page, because a measured figure whose code state is unknown is not evidence (Rule 21).
#
# A commit counts as AI-assisted if its body carries a `Co-Authored-By:` trailer naming Claude or Anthropic --
# the convention this project's own commits follow, and the only machine-readable signal available. A commit
# whose assistant forgot the trailer undercounts silently; that is a real limitation of the method, not of
# this script, and is why the page states the method rather than only the number.
#
using Printf, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))


# ===== 1. Walk the full commit history once: the author date and whether the body carries the trailer ========

function  commitRecords()

    raw = read(`git -C $(ROOT) log --format=%ad%x1f%B%x1e --date=format:%Y-%m-%d`, String)

    recs = Tuple{Date,Bool}[]
    for  chunk  in  split(raw, "\x1e"; keepempty=false)
        parts = split(chunk, "\x1f"; limit=2)
        length(parts) == 2  ||  continue
        d = tryparse(Date, strip(parts[1]))
        d === nothing  &&  continue
        hasTrailer = occursin(r"^Co-Authored-By:.*\b(Claude|Anthropic)\b"mi, parts[2])
        push!(recs, (d, hasTrailer))
    end

    return( recs )

end


# ===== 2. Build the year-by-year table and the whole-history totals ============================================

function  buildTable(recs::Vector{Tuple{Date,Bool}})

    byYear = Dict{Int, Tuple{Int,Int}}()
    for  (d, trailer)  in  recs
        y = Dates.year(d)
        (tot, tr) = get(byYear, y, (0, 0))
        byYear[y] = (tot + 1, tr + (trailer ? 1 : 0))
    end

    years = sort(collect(keys(byYear)))
    rows  = [ @sprintf("| %d | %d | %d |", y, byYear[y][1], byYear[y][2])  for y in years ]

    totalAll     = length(recs)
    totalTrailer = count(r -> r[2], recs)

    trailerDates = [ d  for  (d, t)  in  recs  if  t ]
    firstTrailer = isempty(trailerDates)  ?  nothing  :  minimum(trailerDates)

    return( (rows = rows, totalAll = totalAll, totalTrailer = totalTrailer, firstTrailer = firstTrailer) )

end


# ===== 3. Render the generated block ============================================================================

function  renderBlock(t)

    sinceLine = t.firstTrailer === nothing  ?  "no commit in this history carries the trailer" :
                "the earliest carries author date **" * Dates.format(t.firstTrailer, "d U yyyy") * "**"

    commit = try
        strip(read(`git -C $(ROOT) rev-parse --short HEAD`, String))
    catch
        "unknown"
    end

    block = "<!-- BEGIN generated: development record -->\n" *
            "| | |\n|---|---|\n" *
            "| commits in the whole history | $(t.totalAll) |\n" *
            "| of those, carrying an AI co-author trailer | **$(t.totalTrailer)** |\n\n" *
            "| year | commits | carrying the trailer |\n|---|---|---|\n" *
            join(t.rows, "\n") * "\n\n" *
            "Of the commits carrying the trailer, $sinceLine.\n\n" *
            "*Measured " * Dates.format(Dates.now(), "d U yyyy") * " at commit `" * commit *
            "`; regenerate with `julia --project=. tools/measureDevelopmentRecord.jl --write`.*\n" *
            "<!-- END generated: development record -->"

    return( block )

end

if  abspath(PROGRAM_FILE) == @__FILE__

    recs  = commitRecords()
    table = buildTable(recs)
    block = renderBlock(table)
    println("\n", block)

    if  "--write" in ARGS
        p = joinpath(ROOT, "docs", "src", "development.md");   s = read(p, String)
        a = findfirst("<!-- BEGIN generated: development record -->", s)
        b = findfirst("<!-- END generated: development record -->", s)
        if  a === nothing  ||  b === nothing
            error("measureDevelopmentRecord: the markers are missing from docs/src/development.md; add them " *
                  "around the development-record table before running with --write.")
        end
        write(p, s[1:first(a)-1] * block * s[last(b)+1:end])
        println("\n>> docs/src/development.md updated between the markers.")
    end

end
