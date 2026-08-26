# push!(LOAD_PATH,"../src/")

# Build the plots WITHOUT opening a window. Six @example blocks in docs/src/examples.md call `plot(...)`, and Documenter
# EVALUATES them, so on a machine with a display the GR backend pops up a window for each one on every build. CI never
# noticed because a build server has no display. GKSwstype = "100" is GR's headless mode: savefig still writes the .svg
# files the pages embed, nothing appears on screen. Must be set BEFORE Plots/GR is first loaded.
ENV["GKSwstype"] = "100"

using Documenter, JenaAtomicCalculator

makedocs(;
    modules=[JenaAtomicCalculator],
    format=Documenter.HTML(prettyurls=false, repolink = "https://github.com/OpenJAC/JenaAtomicCalculator.jl",
                           size_threshold = 1000 * 1024,),
    pages=[
        "Home"                          => "index.md",
        "Getting Started"               => "getting-started.md", 
        "For pedestrians"               => "for-pedestrians.md", 
        "Demos"                         => "demos.md",
        "Examples"                      => "examples.md",
        "News"                          => "news.md",
        "API Atomic computations"       => "api-atomic.md",
        "API Atomic processes"          => "api-processes.md",
        "API Atomic properties"         => "api-properties.md",
        "API Basics"                    => "api-basics.md",
        "API Infrastructure"            => "api-infrastructure.md",
        "API Cascade computations"      => "api-cascades.md",
        "API Empirical computations"    => "api-empirical.md",
        "API Plasma computations"       => "api-plasma.md",
        "API Racah algebra"             => "api-racah.md",
        "Bibliography to JAC"           => "reference.md",
        "Getting involved"              => "getting-involved.md",
        "Contributors"                  => "contributors.md",
        "License"                       => "license.md",
    ],
    repo     = "https://github.com/OpenJAC/JenaAtomicCalculator.jl",
    sitename = "JenaAtomicCalculator.jl",
    authors  = "Stephan Fritzsche",
    # `warnonly = Documenter.except()` downgraded EVERY Documenter check to a warning, so the build could not
    # fail whatever was missing -- which is why whole modules were absent from the published API and nothing
    # said so.  The checks below are the ones worth having; the remaining classes stay warn-only so that an
    # incomplete cross-reference does not block a release.
    checkdocs = :exports,
    warnonly  = [:cross_references, :linkcheck, :missing_docs],
)

deploydocs(;
    repo     = "github.com/OpenJAC/JenaAtomicCalculator.jl"
)
