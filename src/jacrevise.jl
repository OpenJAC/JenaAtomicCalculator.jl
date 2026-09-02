
#  Revise is NOT a dependency of JenaAtomicCalculator -- it was removed from Project.toml on 02-Sep-2026, since no
#  code in src/ uses it and every user paid for it in install and precompile time.  This developer helper still
#  finds it through Julia's default LOAD_PATH entry @v#.#, i.e. from the shared environment where Revise is
#  normally installed; `Pkg.add("Revise")` there if `using Revise` below fails.
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using JenaAtomicCalculator

# ] pkg> activate
#   pkg> dev ..   # development environment
#   Pkg.recompile()
#   using JAC


