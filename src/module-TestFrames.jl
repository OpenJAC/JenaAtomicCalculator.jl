
"""
`module  JAC.TestFrames`
... a submodel of JAC that contains all methods for testing individual methods, modules, ....
"""
module TestFrames


using  Printf, SymEngine, JLD2, JenaAtomicCalculator,
       ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

export testDummy


"""
`TestFrames.testCompareLines(oldLine::String, newLine::String)`
    ... compares two output lines token by token: tokens that are parseable as Float64 are compared
        with relative tolerance rtol; all other tokens are compared exactly as strings. Returns true
        if every token agrees within tolerance, false otherwise.
"""
function testCompareLines(oldLine::String, newLine::String; rtol::Float64=1.0e-6)
    oldTokens = split(strip(oldLine))
    newTokens = split(strip(newLine))
    length(oldTokens) != length(newTokens)  &&  return( false )
    for (ot, nt) in zip(oldTokens, newTokens)
        ov = tryparse(Float64, ot)
        nv = tryparse(Float64, nt)
        if  ov !== nothing  &&  nv !== nothing
            denom = max(abs(ov), abs(nv), 1.0e-100)
            abs(ov - nv) / denom > rtol  &&  return( false )
        else
            ot != nt  &&  return( false )
        end
    end
    return( true )
end


function testCompareFiles(fold::String, fnew::String, sa::String, noLines::Int64; rtol::Float64=1.0e-6)
    success = true
    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # Compare the test computations with previous results
    oldLines = readlines(fold)
    newLines = readlines(fnew)
    iold = 0;  for i=1:length(oldLines)   line = oldLines[i];   if  occursin(sa, line)  iold = i;   break   end   end
    inew = 0;  for i=1:length(newLines)   line = newLines[i];   if  occursin(sa, line)  inew = i;   break   end   end
    if  iold == 0   ||   inew == 0    success = false
        println(stdout, "Tries to compare two inappropriate files fold = $(fold); fnew =$(fnew) on string $sa  ($iold, $inew)")
        if printTest   println(iostream, "Tries to compare two inappropriate files fold = $(fold); fnew =$(fnew) on string $sa  ($iold, $inew)")  end
        return( success )
    end
    #
    ii = inew + 1
    for  i = iold+2:iold+noLines
        ii = ii + 1
        if   length(oldLines[i]) < 5    continue    end
        if   !testCompareLines(oldLines[i], newLines[ii]; rtol=rtol)    success = false
            if  printTest  println(iostream, "    *** Old::  " * oldLines[i])
                           println(iostream, "    *** New::  " * newLines[ii])
            end
        end
    end

    return( success )
end


function testPrint(sa::String, success::Bool)
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    ok(succ) =  succ ? "[OK]" : "[Fail]"
    sb = sa * TableStrings.hBlank(110);   sb = sb[1:100] * ok(success);    println(iostream, sb)
    return( nothing )
end


"""
`TestFrames.testDummy(; short::Bool=true)`
    ... placeholder test that always returns true; used to stub out sections not yet implemented.
"""
function testDummy(; short::Bool=true)
    return( true )
end


include("module-TestFrames-inc-general.jl")
include("module-TestFrames-inc-properties.jl")
include("module-TestFrames-inc-processes.jl")
include("module-TestFrames-inc-cascades.jl")
include("module-TestFrames-inc-empirical.jl")
include("module-TestFrames-inc-plasma.jl")
include("module-TestFrames-inc-strongfield.jl")
include("module-TestFrames-inc-liouville.jl")
include("module-TestFrames-inc-deep-learning.jl")


end # module
