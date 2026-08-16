
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


"""
`TestFrames.testCompareValues(oldValues::Array{Float64,1}, newValues::Array{Float64,1})`
    ... compares two lists of numbers with the same relative-tolerance rule that TestFrames.testCompareLines
        applies to its numeric tokens. Returns true if every entry agrees within rtol, false otherwise.
"""
function testCompareValues(oldValues::Array{Float64,1}, newValues::Array{Float64,1}; rtol::Float64=1.0e-6)
    length(oldValues) != length(newValues)  &&  return( false )
    for (ov, nv) in zip(oldValues, newValues)
        denom = max(abs(ov), abs(nv), 1.0e-100)
        abs(ov - nv) / denom > rtol  &&  return( false )
    end
    return( true )
end


"""
`TestFrames.testParseSummaryRow(line::String)`
    ... parses one data row of the two table families that JAC's summary files are built from:

        transition row      3 --   1     1 - --> 0 +    3.512378e+01   E1   Babushkin   <values>
        level row           1    1/2 +   -9.711504e+02  <further energies>
        continuation row                               Coulomb        <values>

        Returns a tuple (identity, labels, indices, values). The identity is what fixes the row PHYSICALLY:
        the angular symmetries. The labels are the non-numeric fields that follow them -- multipole and
        gauge, where present -- and join the identity to form the key under which rows are grouped. The
        printed level numbers are returned separately as `indices` and are deliberately kept OUT of both,
        since they are precisely what a reordering of (near-)degenerate levels changes.

        A continuation row carries no identity of its own -- JAC prints the second gauge of a transition
        this way -- and is returned with an empty identity for the caller to fill in from the row above.
        Returns nothing if the row belongs to none of these shapes.
"""
function testParseSummaryRow(line::String)
    tokens  = string.(split(strip(line)))
    isJ(sa) = occursin(r"^[0-9]+(/2)?$", sa)
    isP(sa) = sa in ["+", "-"]
    numbersOf(tokenList) = Float64[ tryparse(Float64, ta)  for ta in tokenList  if !isnothing(tryparse(Float64, ta)) ]
    labelsOf(tokenList)  = String[ ta  for ta in tokenList  if isnothing(tryparse(Float64, ta)) ]
    ##
    ## A transition row:  i  --  f   J^P  -->  J^P   <energy> [<multipole> <gauge>] <values>
    if  length(tokens) > 8   &&   tokens[2] == "--"   &&   tokens[6] == "-->"   &&
        isJ(tokens[4])  &&  isP(tokens[5])  &&  isJ(tokens[7])  &&  isP(tokens[8])
        iIndex = tryparse(Int64, tokens[1]);    fIndex = tryparse(Int64, tokens[3])
        (isnothing(iIndex)  ||  isnothing(fIndex))   &&   return( nothing )
        rest   = tokens[9:end];    values = numbersOf(rest)
        length(values) == 0   &&   return( nothing )
        identity = "transition  " * tokens[4] * tokens[5] * " --> " * tokens[7] * tokens[8]
        return( (identity=identity, labels=labelsOf(rest), indices=[iIndex, fIndex], values=values) )
    end
    ##
    ## A level row:  level   J  parity   <energies> [<gauge>] <values>
    if  length(tokens) > 3   &&   isJ(tokens[2])   &&   isP(tokens[3])
        lIndex = tryparse(Int64, tokens[1])
        isnothing(lIndex)   &&   return( nothing )
        rest   = tokens[4:end];    values = numbersOf(rest)
        length(values) == 0   &&   return( nothing )
        return( (identity="level  " * tokens[2] * tokens[3], labels=labelsOf(rest), indices=[lIndex], values=values) )
    end
    ##
    ## A continuation row:  <labels> <values>, with the identity inherited from the row above
    if  length(tokens) > 1   &&   isnothing(tryparse(Float64, tokens[1]))
        ifirst = findfirst(ta -> !isnothing(tryparse(Float64, ta)), tokens)
        (isnothing(ifirst)  ||  ifirst == 1)                                        &&   return( nothing )
        all(ta -> !isnothing(tryparse(Float64, ta)), tokens[ifirst:end])            ||   return( nothing )
        return( (identity="", labels=tokens[1:ifirst-1], indices=Int64[], values=numbersOf(tokens[ifirst:end])) )
    end
    ##
    return( nothing )
end


"""
`TestFrames.testCollectSummaryRows(lines::Array{String,1}, ianchor::Int64, noLines::Int64, fname::String)`
    ... collects the parsed data rows of the noLines-block that follows the anchor line ianchor, and returns
        them as tuples (key, indices, values). Blank lines, ruler lines and header/unit lines -- recognized by
        carrying no number at all -- are skipped silently.

        A continuation row inherits from the row above: its identity, its level numbers, and those leading
        values that it does not print itself, so that the inherited transition energy remains available for
        the ordering. It is accepted only if the row above has the same number of labels and more values,
        which keeps an ordinary line of prose from being mistaken for one.

        A line that does carry numbers but belongs to none of the supported shapes is REFUSED loudly and makes
        the method return nothing; a comparator that quietly does something plausible on a table it does not
        understand would be worse than none.
"""
function testCollectSummaryRows(lines::Array{String,1}, ianchor::Int64, noLines::Int64, fname::String)
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    keyOf(identity, labels) = identity * (isempty(labels) ? "" : "  " * join(labels, " "))
    rows = Any[];    previous = nothing
    for  i = ianchor+2:min(ianchor+noLines, length(lines))
        line = lines[i]
        length(strip(line)) < 5                                      &&   continue
        occursin(r"^[\s\-=]+$", line)                                &&   continue
        all(ta -> isnothing(tryparse(Float64, ta)), split(line))     &&   continue
        row  = TestFrames.testParseSummaryRow(line)
        ##
        ## Resolve a continuation row against the row above, or refuse the line
        if  !isnothing(row)  &&  row.identity == ""
            if  isnothing(previous)                             ||   isempty(previous.labels)  ||
                length(row.labels) != length(previous.labels)   ||   length(row.values) >= length(previous.values)
                row = nothing
            else
                nInherited = length(previous.values) - length(row.values)
                row = (identity=previous.identity, labels=row.labels, indices=previous.indices,
                       values=vcat(previous.values[1:nInherited], row.values))
            end
        end
        if  isnothing(row)
            sa = "*** Unsupported table shape in $(fname), line $i;  no comparison by key is attempted:"
            println(stdout, sa);   println(stdout, "    " * line)
            if printTest   println(iostream, sa);   println(iostream, "    " * line)   end
            return( nothing )
        end
        previous = row
        push!(rows, (key=keyOf(row.identity, row.labels), indices=row.indices, values=row.values))
    end
    return( rows )
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
