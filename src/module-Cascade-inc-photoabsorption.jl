
# Functions and methods for scheme::Cascade.PhotoAbsorptionSchem computations


"""
`Cascade.perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true,
                    outputDirectory::String="")`
    ... to set-up and perform a photoabsorption cascade computation that accounts for both the direct and the resonant part of
        photoabsorption. Apart from the direct photoionization of the atom, it enables one to treat the photoexcitation and the
        subsequent electron emission (autoionization) in the cascade. Typical photoabsorption properties are energy-dependent
        photoionization cross sections, photoabsorption spectra, PI plasma rate coefficients and others. The method runs the
        photoionization and the photoexcitation part in turn and merges their results.
        With output = false, nothing is returned. With output = true, a results::Basics.PerformResults is returned -- an
        AbstractDict{String,Any} -- that carries the initial and the generated multiplets, the individual line lists, and
        under the key "cascade data:" the Cascade.Data of the two parts that a subsequent Cascade.Simulation reads. With
        outputToFile = true the same data is written to a zzz-cascade-photoabsorption-computations-<date>.jld file, into outputDirectory when one is
        named and into the working directory otherwise.
"""
function perform(scheme::PhotoAbsorptionScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")

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


    if  scheme.calcResonant
        println("\n> Resonant photoexcitation part of the photoabsorption cascade computations")
        println(  "  =========================================================================  \n")
        eScheme = Cascade.PhotoExcitationScheme(scheme.multipoles, scheme.photonEnergies[1], scheme.photonEnergies[end], 1, 
                            scheme.excitationFromShells, scheme.excitationToShells, scheme.initialLevelSelection, 
                            scheme.lValues, scheme.electronEnergyShift, scheme.minCrossSection)
        eComp   = Cascade.Computation(comp, scheme=eScheme)
        eOut    = Cascade.perform(eScheme, eComp; output=true, outputToFile=false)
    end

    # Collect the output from the computations above if required
    data = Cascade.Data[]
    if output   &&   scheme.calcDirect 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => iOut[ResultKeys.InitialMultiplets]) )    
        results = Base.merge( results, Dict("photoionized multiplets:"      => iOut[ResultKeys.GeneratedMultiplets]) )    
        results = Base.merge( results, Dict("photoionization lines:"        => iOut[PhotoIonization.Settings]) )
        push!(data, Cascade.Data{PhotoIonization.Line}(iOut[PhotoIonization.Settings]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end

    if output   &&   scheme.calcResonant 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => eOut[ResultKeys.InitialMultiplets]) )    
        results = Base.merge( results, Dict("photoexcited multiplets:"      => eOut[ResultKeys.GeneratedMultiplets]) )    
        results = Base.merge( results, Dict("photoexcitation lines:"        => eOut[PhotoExcitation.Settings]) )
        push!(data, Cascade.Data{PhotoExcitation.Line}(eOut[PhotoExcitation.Settings]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end

    #  Write out the result to file to later continue with simulations on the cascade data
    if  outputToFile
        if  outputDirectory != ""    mkpath(outputDirectory)    end
        filename = joinpath(outputDirectory, "zzz-cascade-photoabsorption-computations-" * string(Dates.now())[1:13] * ".jld")
        println("\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file.")
        if  printSummary   println(iostream, "\n* Write all results to disk; use:\n   JLD2.save(''$filename'', results) \n   using JLD2 " *
                                                "\n   results = JLD2.load(''$filename'')    ... to load the results back from file." )      end      
        Cascade.writeDataFile(filename, results)
    end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end
