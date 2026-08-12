

using Printf, ..Basics, ..Defaults, ..Nuclear, ..ManyElectron, ..Radial, ..TableStrings

#################################################################################################################################
#################################################################################################################################


"""
`Plasma.performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                  plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=false)`
    ... computes and diagonalizes the CI matrix for all CSF in the given basis under the given plasma model,
        mirroring Hamiltonian.performCI but replacing the field-free Hamiltonian.setupMatrix by the plasma-screened
        Basics.compute(JP, basis, nm, grid, settings, plasmaModel). For plasmaModel = Basics.NoPlasmaModel(), this
        delegates directly to Hamiltonian.performCI, so that the field-free computation of the standard processes
        is never touched by the plasma-specific code path. A multiplet::Multiplet is returned.
"""
function performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                   plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=false)
    if  typeof(plasmaModel) == Basics.NoPlasmaModel
        return( Hamiltonian.performCI(basis, nm, grid, settings; printout=printout) )
    end

    # First determine the number of CSF in each J^P symmetry block
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else                              symmetries[sym] = 1
        end
    end

    NoCsf = 0;   for (sym,v) in symmetries   NoCsf = NoCsf + v   end
    if  NoCsf != length(basis.csfs)   error("stop b; NoCsf = $NoCsf ")   end

    # Calculate for each symmetry block the corresponding plasma-screened CI matrix, diagonalize it and append
    # a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        if  !Basics.selectSymmetry(sym, settings.levelSelectionCI)     continue    end
        matrix = Basics.compute(sym, basis, nm, grid, settings, plasmaModel; printout=printout)
        eigen  = Basics.fixEigenvectorPhase!( Basics.diagonalize(MatrixWithLinearAlgebra(), matrix) )

        levels = Level[]
        for  ev = 1:length(eigen.values)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs)
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector )
            push!( levels, newlevel)
        end
        wa = Multiplet(string(sym) * "+ (plasma)", levels)
        push!( multiplets, wa)
    end

    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)

    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)
        if  Basics.selectLevel(level, settings.levelSelectionCI)    push!(levelNos, ilev)    end
    end

    if  printout    Basics.tabulate(stdout, mp, levelNos)    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Basics.tabulate(iostream, mp, levelNos)    end

    return( mp )
end


"""
`Plasma.displayResults(stream::IO, multiplet::Multiplet, pMultiplet::Multiplet, plasmaModel::Basics.AbstractPlasmaModel)`
    ... to display the energies, M_ms and F-parameters, etc. for the  selected levels. A neat table is printed but nothing
        is returned otherwise.
"""
function  displayResults(stream::IO, multiplet::Multiplet, pMultiplet::Multiplet, plasmaModel::Basics.AbstractPlasmaModel)
    nx = 64
    println(stream, " ")
    println(stream, " ")
    println(stream, "  Plasma shifts for $(plasmaModel):")
    println(stream,   "     + Plasma screening included perturbatively in the CI matrix but not in the SCF field.")

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
    pm = scheme.plasmaModel

    # Determine self-consistent (field-free) orbitals for the initial and final configurations; the orbitals
    # themselves are not re-optimized in the plasma potential (only the CI matrix and, further downstream, the
    # continuum orbital and transition amplitude are screened -- see AutoIonization/PhotoIonization). The CI
    # step is then redone under the given plasma model to obtain the actual plasma-shifted multiplets.
    initialFieldFree = SelfConsistent.performSCF(scheme.initialConfigs, nm, computation.grid, computation.asfSettings)
    finalFieldFree   = SelfConsistent.performSCF(scheme.finalConfigs,   nm, computation.grid, computation.asfSettings)
    initialMultiplet = Plasma.performCI(initialFieldFree.levels[1].basis, nm, computation.grid, computation.asfSettings, pm)
    finalMultiplet   = Plasma.performCI(finalFieldFree.levels[1].basis,   nm, computation.grid, computation.asfSettings, pm)
    #
    # Display the plasma-induced level shifts of the initial and final multiplets, before the process-specific rates
    println("\n  Initial-state levels:")
    Plasma.displayResults(stdout, initialFieldFree, initialMultiplet, pm)
    println("  Final-state levels:")
    Plasma.displayResults(stdout, finalFieldFree, finalMultiplet, pm)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        Plasma.displayResults(iostream, initialFieldFree, initialMultiplet, pm)
        Plasma.displayResults(iostream, finalFieldFree, finalMultiplet, pm)
    end
    #
    if      typeof(scheme.settings)  == AutoIonization.PlasmaSettings
        outcome = AutoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nm, computation.grid, scheme.settings, pm)
        if output    results = Base.merge( results, Dict("AutoIonization lines in plasma:" => outcome) )         end
    elseif  typeof(scheme.settings)  == PhotoIonization.PlasmaSettings
        outcome = PhotoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nm, computation.grid, scheme.settings, pm)
        if output    results = Base.merge( results, Dict("Photo lines in plasma:" => outcome) )                  end
    else
        error("Unsupported line-shift settings type = $(typeof(scheme.settings)).")
    end

    println("Line-Shift computation complete ...")
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end

