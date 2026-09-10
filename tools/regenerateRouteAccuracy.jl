#
# regenerateRouteAccuracy.jl -- re-derive the route-accuracy figures quoted in docs/src/scf-routes.md.
#
#   Run from the JAC root:   julia --project=. tools/regenerateRouteAccuracy.jl
#
# WHY THIS EXISTS.  The page tells a user how far to trust each SCF route, and until now those figures were
# typed in by hand.  Three times in one week a hand-recorded number in this project turned out to have been
# measured on a code state that had since moved -- the 28-Aug Stobbe comparison, the example-Ai tables, and the
# route accuracies themselves.  A figure a user relies on to decide whether to publish a number must be
# re-derivable, so this script re-derives it and rewrites ONLY the region between the markers in the page:
#
#     <!-- BEGIN generated: route accuracy -->  ...  <!-- END generated: route accuracy -->
#
# Everything outside the markers is prose and stays under human control.  The script PRINTS the table and, with
# --write, replaces the region; it changes nothing else.  It also stamps the commit and the date into the page,
# because a measured figure whose code state is unknown is not evidence (Rule 21).
using JenaAtomicCalculator, Printf
const SC = JenaAtomicCalculator.SelfConsistent

cases = [ ("Be-like, Z = 4",   4.0, [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]),
          ("Be-like, Z = 10", 10.0, [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]),
          ("Be-like, Z = 26", 26.0, [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]),
          ("Be-like, Z = 92", 92.0, [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]),
          ("C-like, Z = 6",    6.0, [Configuration("1s^2 2s^2 2p^2"), Configuration("1s^2 2s^1 2p^3"),
                                     Configuration("1s^2 2p^4")]) ]
quiet(f) = (t = tempname(); r = open(t,"w") do io; redirect_stdout(io) do; f() end end; rm(t,force=true); r)

rows = String[]
for (label, Z, cfs)  in  cases
    nm   = Nuclear.Model(Z);   grid = Basics.recommendedGrid(cfs, nm)
    sel  = LevelSelection(true, indices=[1])
    mkS(route) = AsfSettings(AsfSettings(); scField=Basics.EOLField(), scfRoute=route,
                             accuracyScf=1.0e-8, levelSelectionCI=sel)
    eRot = quiet(() -> sort(SC.performSCF(cfs, nm, grid, mkS(Basics.RotationRoute(2000));
                                          printout=false).levels, by=l->l.energy)[1].energy)
    eFck = quiet(() -> sort(SC.performSCF(cfs, nm, grid, mkS(Basics.FockRoute(60));
                                          printout=false).levels, by=l->l.energy)[1].energy)
    dAbs = (eFck - eRot) * 1000.0                      # mHa, positive = Fock lies ABOVE
    dRel = abs(eFck - eRot) / abs(eRot)
    push!(rows, @sprintf("| %-16s | %.1f mHa | %.0e |", label, dAbs, dRel))
    @printf("  %-16s  rotation %18.9f   Fock %18.9f   %+8.2f mHa   %.1e\n", label, eRot, eFck, dAbs, dRel)
end

commit = strip(read(`git rev-parse --short HEAD`, String))
block  = "<!-- BEGIN generated: route accuracy -->\n" *
         "| system | absolute | relative to the total energy |\n|---|---|---|\n" *
         join(rows, "\n") * "\n\n" *
         "*Measured " * Libc.strftime("%d %B %Y", time()) * " at commit `" * commit *
         "`; regenerate with `julia --project=. tools/regenerateRouteAccuracy.jl --write`.*\n" *
         "<!-- END generated: route accuracy -->"
println("\n", block)

if  "--write" in ARGS
    p = joinpath(@__DIR__, "..", "docs", "src", "scf-routes.md");   s = read(p, String)
    a = findfirst("<!-- BEGIN generated: route accuracy -->", s)
    b = findfirst("<!-- END generated: route accuracy -->", s)
    if  a === nothing  ||  b === nothing
        error("regenerateRouteAccuracy: the markers are missing from docs/src/scf-routes.md; add them around " *
              "the accuracy table before running with --write.")
    end
    write(p, s[1:first(a)-1] * block * s[last(b)+1:end])
    println("\n>> docs/src/scf-routes.md updated between the markers.")
end
