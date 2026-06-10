

using Printf, ..Basics, ..Defaults, ..Nuclear, ..ManyElectron, ..Radial, ..TableStrings

#################################################################################################################################
#################################################################################################################################


"""
`Plasma.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    asfSettings::AsfSettings, settings::Plasma.Settings; output=true)` 
    ... to compute a new multiplet from the plasma-modified Hamiltonian matrix as specified by the given 
        settings. The diagonalization of the (plasma-modified) Hamiltonian matrix follows the asfSettings that
        were applied before for the computation of the multiplet; warnings are issued if features were selected
        that are not supported for plasma calculations, such as the Breit interaction, QED shifts and others.
        The energies and plasma shifts (with regard to the unperturbed multiplet) are printed to screen 
        and into the summary file. The generated (plasma) multiplet is returned if output=true and nothing otherwise,
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                            asfSettings::AsfSettings, settings::Plasma.Settings; output=true)
    println("")
    printstyled("Plasma.computeOutcomes(): The computation of plasma-shifts starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------- \n", color=:light_green)
    #
    basis        = multiplet.levels[1].basis
    newMultiplet = Basics.perform("computation: CI for plasma", basis, nm, grid, asfSettings, settings)
    # Print all results to screen
    Plasma.displayResults(stdout, multiplet, newMultiplet, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Plasma.displayResults(iostream, multiplet, newMultiplet, settings)   end
    #
    if    output    return( newMultiplet )
    else            return( nothing )
    end
end


"""
`Plasma.displayResults(stream::IO, multiplet::Multiplet, newMultiplet::Multiplet, settings::Plasma.Settings)`  
    ... to display the energies, M_ms and F-parameters, etc. for the  selected levels. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, multiplet::Multiplet, pMultiplet::Multiplet, settings::Plasma.Settings)
    nx = 64
    println(stream, " ")
    println(stream, " ")
    println(stream, "  Plasma shifts for $(settings.plasmaModel):")
    
    # Write out the relevant plasma parameters
    if settings.plasmaModel == Plasma.DebyeHueckel()
        println(stream, "\n     + lambda = $(settings.lambdaDebye)")
        println(stream,   "     + Plasma screening included perturbatively in CI matrix but not in SCF field.")
    else  error("Unsupported plasma model = $(plasmaSettings.plasmaModel)")
    end
    
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(18, "Energy w/o plasma"; na=4)              
    sb = sb * TableStrings.center(18, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Delta E";     na=4)              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  i  in  1:length(multiplet.levels)
        sa  = "  ";    sym = LevelSymmetry( multiplet.levels[i].J, multiplet.levels[i].parity);    
                    newsym = LevelSymmetry( pMultiplet.levels[i].J, pMultiplet.levels[i].parity)
        sa = sa * TableStrings.center(10, TableStrings.level(multiplet.levels[i].index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=6)
        energy = multiplet.levels[i].energy
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))              * "     "
        if  sym == newsym   deltaE = pMultiplet.levels[i].energy - energy;    
                            sa     = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", deltaE))   
        else                sa     = sa * "Level crossing."
        end
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end


    
"""
`Plasma.perform(scheme::Plasma.LineShiftScheme, computation::Plasma.Computation; output::Bool=true)`  
    ... to perform a Saha-Boltzmann equilibrium computation for a given ion mixture. For output=true, a dictionary 
        is returned from which the relevant results can be can easily accessed by proper keys.
"""
function  perform(scheme::Plasma.LineShiftScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
        
    nm = computation.nuclearModel
        
    # Use computation.grid, computation.nuclearModel, scheme.initialConfigs, scheme.finalConfigs, scheme.settings, 
    initialMultiplet = SelfConsistent.performSCF(scheme.initialConfigs, nm, computation.grid, computation.asfSettings)
    finalMultiplet   = SelfConsistent.performSCF(scheme.finalConfigs, nm, computation.grid, computation.asfSettings)
    #
    if      typeof(scheme.settings)  == AutoIonization.PlasmaSettings
        outcome = AutoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nm, computation.grid, scheme.settings) 
        if output    results = Base.merge( results, Dict("AutoIonization lines in plasma:" => outcome) )         end
    elseif  typeof(scheme.settings)  == PhotoIonization.PlasmaSettings
        outcome = PhotoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nm, computation.grid, scheme.settingss) 
        if output    results = Base.merge( results, Dict("Photo lines in plasma:" => outcome) )                  end
    else
        error("stop a")
    end
    
    println("Line-Shift computation complete ...")
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end

