
# Functions and methods for scheme::Cascade.ElectronExcitationSchem computations

"""
`Cascade.perform(scheme::ElectronExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true,
                 outputDirectory::String="")`
    ... sets up and performs an electron-excitation cascade computation. The scheme covers two channels and selects them by
        its own processes list: Basics.ImpactExc() requests the DIRECT electron-impact excitation, and Basics.ImpactExcAuto()
        the RESONANT part, i.e. dielectronic capture followed by re-autoionization. Typical properties are energy-dependent
        collision strengths, impact-excitation cross sections and effective collision strengths. A results::Dict{String,Any}
        is returned if output=true, and nothing otherwise; with outputToFile=true the same dictionary is written to a .jld
        file for use in a subsequent cascade simulation.

        The direct channel hands the scheme's own configuration -- shells, l-values, free-electron grid and energy shift -- to a
        Cascade.ImpactExcitationScheme, and raises if the shell lists are empty rather than failing further in.

        ONLY THE DIRECT CHANNEL IS AVAILABLE. The resonant one raises, because Cascade.DielectronicCaptureScheme has no
        perform method in JAC; see the error text for what is missing. Until 16-Aug-2026 this function could not run at all:
        it read scheme.calcDirect, scheme.calcResonant and scheme.multipoles, none of which are fields of
        Cascade.ElectronExcitationScheme, so it raised on the first line that touched the scheme.
"""
function perform(scheme::ElectronExcitationScheme, comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true, outputDirectory::String="")
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # The two channels are selected by the scheme's own processes list; there are no calcDirect/calcResonant fields.
    calcDirect   = Basics.ImpactExc()      in  scheme.processes
    calcResonant = Basics.ImpactExcAuto()  in  scheme.processes
    if  !calcDirect  &&  !calcResonant
        error("Cascade.perform(::ElectronExcitationScheme): the scheme requests neither channel.  Put " *
              "Basics.ImpactExc() into its processes for the direct electron-impact excitation, and/or " *
              "Basics.ImpactExcAuto() for the resonant (dielectronic capture with re-autoionization) part; " *
              "the scheme currently carries $(scheme.processes).")
    end
    #
    # Perform the two cascade computations for the electron-impact and dielectronic-capture with re-autoionzation separately, 
    # even if this requires to run the computation of the initial multiplet twice.
    if  calcDirect
        println("\n> Direct impact-excitation part of the electron-excitation cascade computations")
        println(  "  =============================================================================  \n")
        # The direct channel IS an impact excitation, so the scheme hands its configuration straight over.
        if  isempty(scheme.fromShells)  ||  isempty(scheme.toShells)
            error("Cascade.perform(::ElectronExcitationScheme): the direct channel needs both fromShells and " *
                  "toShells; the scheme carries fromShells = $(scheme.fromShells) and toShells = " *
                  "$(scheme.toShells).  Without them the impact-excitation cascade has nothing to excite " *
                  "between and fails further in with an unrelated message.")
        end
        ieScheme = Cascade.ImpactExcitationScheme(scheme.fromShells, scheme.toShells, scheme.electronEnergies,
                                                  scheme.lValues, scheme.NoFreeElectronEnergies,
                                                  scheme.maxFreeElectronEnergy, scheme.electronEnergyShift)
        ieComp   = Cascade.Computation(comp, scheme=ieScheme)
        ieOut    = Cascade.perform(ieScheme, ieComp; output=true, outputToFile=false)
    end
    #
    #
    if  calcResonant
        error("Cascade.perform(::ElectronExcitationScheme): the RESONANT channel is not available.  It would run a " *
              "Cascade.DielectronicCaptureScheme, for which JAC has no perform method, and would collect its lines as " *
              "Cascade.Data{DielectronicCapture.Line} -- a module that does not exist either.  Both have to be written " *
              "before Basics.ImpactExcAuto() can be requested; the DIRECT channel, Basics.ImpactExc(), works and can be " *
              "used on its own.  Raising here rather than failing later on a missing name is deliberate (Rule 13).")
    end
    #
    # Collect the output from the computations above if required
    data = Cascade.Data[]
    if output   &&   calcDirect 
        results = Base.merge( results, Dict("name"                          => comp.name) ) 
        results = Base.merge( results, Dict("cascade scheme"                => comp.scheme) ) 
        results = Base.merge( results, Dict("initial multiplets:"           => ieOut["initial multiplets:"]) )    
        results = Base.merge( results, Dict("impact-excited multiplets:"    => ieOut["generated multiplets:"]) )    
        results = Base.merge( results, Dict("impact-excitation lines:"      => ieOut["impact-excitation lines:"]) )
        push!(data, Cascade.Data{ImpactExcitation.Line}(ieOut["impact-excitation lines:"]) )
        results = Base.merge( results, Dict("cascade data:"                 => data ) )
    end
    #
    ## The resonant output block is unreachable while the resonant channel raises above, and it referenced
    ## DielectronicCapture.Line -- a module JAC does not have.  It is removed rather than left as dead code
    ## that names something undefined; the error above records what has to be built.

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
