#
println("Gb) Symbolic simplification of Racah expressions by means of sum rules.")
#
using SymEngine

if  true
    # Last successful:  13May2024
    # Test for sum rules; the 'simple form' (without additional phases) has now been tested for all (1..47) 
    # selectRacahExpression(1..45); with additional phases for selectRacahExpression(1..30, 32..47)
    # test simplification of deltas with 31, 38, 39
    rex = RacahAlgebra.selectRacahExpression(41);          println("rex-original        = $rex") 
    ## rex = RacahAlgebra.equivalentForm(rex);                println("rex-equivalent      = $rex")
    rex = RacahAlgebra.evaluate(rex);
    #
end
