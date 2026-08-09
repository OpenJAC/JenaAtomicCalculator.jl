
"""
`module  JAC.BasicsAZ`  
 	   ... a submodel of JAC that contains methods from the Basics module but with reference to many other modules.
"""
module BasicsAZ

using  Printf,  LinearAlgebra, GaussQuadrature, JenaAtomicCalculator, ..Basics, ..TableStrings
       ## using JenaAtomicCalculator ... since otherwise almost all other modules must be included explicitly
       


include("module-BasicsAZ-inc-AG.jl")
include("module-BasicsAZ-inc-compute.jl")
include("module-BasicsAZ-inc-configurations.jl")
include("module-BasicsAZ-inc-generate.jl")
include("module-BasicsAZ-inc-HP.jl")
include("module-BasicsAZ-inc-perform.jl")
include("module-BasicsAZ-inc-QZ.jl")

end # module
