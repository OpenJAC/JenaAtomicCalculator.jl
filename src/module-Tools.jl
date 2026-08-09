
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



#==
"""
`Tools.perform(obs::Observable{Any})`  ... performs different tasks due to the given Observable.
"""
function perform(obs::Observable{Any})
    #
    if      obs[] == "Start task"          println("Nothing is implemented here.")
    elseif  obs[] == "Grid calculator"     Tools.taskGridCalculator()
    else                                   error("stop a")
    end
    return( nothing )
end


"""
`Tools.select()`  ... select different tools from a menu for which no further input is required.
"""
function select()
    t1 = "Simple tools JAC without input: "
    b1 = dropdown(["Start task", "Grid calculator"])
    update = button("Update")
    ui = vbox( hbox( pad(1em, t1) ),
                hbox( pad(1em, b1), pad(1em, update) )
                )
    Interact.display(ui)
    output = Interact.@map (&update;  observe(b1)[] ) 
    return( output )
end


"""
`Tools.select(dict::Dict)`  
    ... select different tools from a menu if proper results::dict are given in terms of a dictionary as obtained
        usually from an Atomic.Computation.
"""
function select(dict::Dict)
    println("Tools.select(): A menu if a dict comes in; $dict ")
end


"""
`Tools.taskGridCalculator()`  ... grid calculator to determine useful grid parameters
"""
function taskGridCalculator()
    t1 = "Grid calculator: Determine 'one' of the grid parameters if all others are given explicitly by the user. " *
            "There are two forms of the radial grid: "
    t2 = "   (i)    logarithmic grid   r(i)  = rnt * (Base.MathConstants.e^{(i-1)*h} - 1)            ... for bound-state electrons  \n" *
            "   (ii)   log-lin grid       ln(r(i)/rnt+1) +  (h/hp)*r(i) = (i-1)*h    ... for free electrons.  \n\n"
    t3 = "The logarithmic grid is used always if hp = 0. \n" *
            "The use of this 'calculator' has been found useful in Auger and photoionization studies as well as if " *
            "Rydberg-type electrons need to be taken into account in the bound-state density."
    #
    b1 = spinbox(label="rnt  ";    value=2.0e-6)
    b2 = spinbox(label="h     ";   value=5.0e-2)
    b3 = spinbox(label="hp    ";   value=0.0)
    b4 = spinbox(label="r_max ";   value=5.6e+2)
    update = button("Update")
    ui = vbox( hbox( pad(0em, t1) ), 
                hbox( pad(0em, t2) ),
                hbox( pad(0em, t3) ),
                hbox( pad(1em, b1), pad(1em, b2), pad(1em, b3), pad(1em, b4), pad(1em, update) )
                )
    Interact.display(ui) 
    Interact.@map (&update; taskGridCalculatorResults(observe(b1)[], observe(b2)[], observe(b3)[], observe(b4)[]) ) 
    return( nothing )
end
==#


"""
`Tools.taskGridCalculatorResults(rnt::Float64, h::Float64, hp::Float64, rmax::Float64)`  ... prints the results of the grid computations
"""
function taskGridCalculatorResults(rnt::Float64, h::Float64, hp::Float64, rmax::Float64)
    println("** $rnt  $h  $hp  $rmax")
end

end # module

