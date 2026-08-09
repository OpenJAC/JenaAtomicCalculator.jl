
"""
`module  JAC.Tools`  
... a submodel of JAC that contains all methods for providing a toolbox for -- more o less --
    simple tools for the JAC program.
"""
module Tools


## using Interact

"""
`struct  Tools.Settings`  
    ... defines a type for the settings in estimating double-Auger and autoionization rates.

    + printBefore  ::Bool   ... True, if all energies and lines are printed before their evaluation.
"""
struct Settings
    printBefore    ::Bool
end 





"""
`Tools.taskGridCalculatorResults(rnt::Float64, h::Float64, hp::Float64, rmax::Float64)`  ... prints the results of the grid computations
"""
function taskGridCalculatorResults(rnt::Float64, h::Float64, hp::Float64, rmax::Float64)
    println("** $rnt  $h  $hp  $rmax")
end

end # module

