# push!(LOAD_PATH,"../src/")

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
    warnonly = Documenter.except(),
)

deploydocs(;
    repo     = "github.com/OpenJAC/JenaAtomicCalculator.jl"
)
