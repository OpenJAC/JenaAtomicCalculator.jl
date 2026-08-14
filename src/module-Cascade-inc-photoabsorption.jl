
# Functions and methods for scheme::Cascade.PhotoAbsorptionSchem computations


"""
`Cascade.perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation)`  
    ... to set-up and perform a photoabsorption cascade computation that that accounts for the direct and resonant 
        of photoabsorption. Apart from the direct photoionization of the atom, it enables one to treat the photoexcitation 
        and subsequent electron emission (autoionization) in the cascade. Typical photoabsorption properties are 
        energy-dependent photoionization cross sections, photoabsorption spectra, PI plasma rate coefficients, and others.
        The results of a photoabsorption cascade computation are comprised into (output) data::ExcitationData, while 
        nothing is returned otherwise.

`Cascade.perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary that is written to disk and can be used in subsequent
        cascade simulation. The particular output depends on the specifications of the cascade.
"""
function perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    # Perform the two cascade computations for the photoionization and photoexcitation separately, even if this
    # requires to run the computation of the initial multiplet twice.
    if  scheme.calcDirect
        println("\n> Direct photoionization part of the photoabsorption cascade computations")
        println(  "  =======================================================================  \n")
        iScheme = Cascade.PhotoIonizationScheme(scheme.multipoles, scheme.photonEnergies, scheme.electronEnergies, 
                            scheme.excitationFromShells, scheme.excitationToShells, scheme.initialLevelSelection, 
                            scheme.lValues, scheme.electronEnergyShift, scheme.minCrossSection)
        iComp   = Cascade.Computation(comp, scheme=iScheme)
        iOut    = Cascade.perform(iScheme, iComp; output=true, outputToFile=false)
    end
    #
    #
    if  scheme.calcResonant
        println("\n> Resonant photoexcitation part of the photoabsorption cascade computations")
        println(  "  =========================================================================  \n")
        eScheme = Cascade.PhotoExcitationScheme(scheme.multipoles, scheme.photonEnergies[1], scheme.photonEnergies[end], 1, 
                            scheme.excitationFromShells, scheme.excitationToShells, scheme.initialLevelSelection, 
                            scheme.lValues, scheme.electronEnergyShift, scheme.minCrossSection)
        eComp   = Cascade.Computation(comp, scheme=eScheme)
        eOut    = Cascade.perform(eScheme, eComp; output=true, outputToFile=false)
    end
    #
    # Collect the output from the computations above if required
    data = Cascade.Data[]
    if output   &&   scheme.calcDirect 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => iOut["initial multiplets:"]) )    
        results = Base.merge( results, Dict("photoionized multiplets:"      => iOut["generated multiplets:"]) )    
        results = Base.merge( results, Dict("photoionization lines:"        => iOut["photoionization lines:"]) )
        push!(data, Cascade.Data{PhotoIonization.Line}(iOut["photoionization lines:"]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end
    #
    if output   &&   scheme.calcResonant 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => eOut["initial multiplets:"]) )    
        results = Base.merge( results, Dict("photoexcited multiplets:"      => eOut["generated multiplets:"]) )    
        results = Base.merge( results, Dict("photoexcitation lines:"        => eOut["photoexcitation lines:"]) )
        push!(data, Cascade.Data{PhotoExcitation.Line}(eOut["photoexcitation lines:"]) )
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
        Cascade.writeDataFile(filename, results)
    end
    
    return( results )
end
