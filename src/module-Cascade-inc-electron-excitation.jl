
# Functions and methods for scheme::Cascade.ElectronExcitationSchem computations

"""
`Cascade.perform(scheme::ElectronExcitationScheme, comp::Cascade.Computation)`  
    ... to set-up and perform an electron-excitation cascade computation that accounts for the direct and resonant part. 
        Apart from the direct electron-impact excitation of the atom, it enables one to treat the (resonant) dielectronic
        capture with subsequent electron emission (re-autoionization) in the cascade. Typical electron-excitation
        properties are energy-dependent collision strengths, impact-excitation cross sections, effective collision strengths,
        and several others.

`Cascade.perform(scheme::ElectronExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in 
        subsequent cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::ElectronExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the two cascade computations for the electron-impact and dielectronic-capture with re-autoionzation separately, 
    # even if this requires to run the computation of the initial multiplet twice.
    if  scheme.calcDirect
        println("\n> Direct impact-excitation part of the electron-excitation cascade computations")
        println(  "  =============================================================================  \n")
        ieScheme = Cascade.ImpactExcitationScheme(scheme.multipoles)
        ieComp   = Cascade.Computation(comp, scheme=ieScheme)
        ieOut    = Cascade.perform(ieScheme, ieComp; output=true, outputToFile=false)
    end
    #
    #
    if  scheme.calcResonant
        println("\n> Resonant dielectronic capture part of the electron-excitation cascade computations")
        println(  "  ==================================================================================  \n")
        dcScheme = Cascade.DielectronicCaptureScheme(scheme.multipoles)
        dcComp   = Cascade.Computation(comp, scheme=dcScheme)
        dcOut    = Cascade.perform(dcScheme, dcComp; output=true, outputToFile=false)
    end
    #
    # Collect the output from the computations above if required
    data = Cascade.Data[]
    if output   &&   scheme.calcDirect 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => ieOut["initial multiplets:"]) )    
        results = Base.merge( results, Dict("impact-excited multiplets:"    => ieOut["generated multiplets:"]) )    
        results = Base.merge( results, Dict("impact-excitation lines:"      => ieOut["impact-excitation lines:"]) )
        push!(data, Cascade.Data{ImpactExcitation.Line}(ieOut["impact-excitation lines:"]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end
    #
    if output   &&   scheme.calcResonant 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => dcOut["initial multiplets:"]) )    
        results = Base.merge( results, Dict("dielectronic multiplets:"      => dcOut["generated multiplets:"]) )    
        results = Base.merge( results, Dict("dielectronic-capture lines:"   => dcOut["dielectronic-capture lines:"]) )
        push!(data, Cascade.Data{DielectronicCapture.Line}(dcOut["dielectronic-capture lines:"]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end
    #
    #  Write out the result to file to later continue with simulations on the cascade data
    if  outputToFile
        filename = "zzz-cascade-photoabsorption-computations-" * string(Dates.now())[1:13] * ".jld"
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        JLD2.@save filename results
    end
    
    return( results )
end
