
"""
`module  JAC.ReducedDensityMatrix`  
... a submodel of JAC that contains all methods for computing reduced density matrices, natural orbitals, radial density
    distributions and related information for some level(s).
"""
module ReducedDensityMatrix


using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
                ..SpinAngular, ..TableStrings

"""
`struct  ReducedDensityMatrix.Settings  <:  AbstractPropertySettings`  
    ... defines a type for providing and parameters for reduced density matrices and other parameters.

    + calcNatural              ::Bool             ... True if natural orbitals need to be calculated, and false otherwise.
    + calcDensity              ::Bool             ... True if the radial density need to be calculated, and false otherwise.
    + calcIpq                  ::Bool             ... True if orbital interaction I_pq need to be calculated, and false otherwise.
    + calc2pRDM                ::Bool             ... True if the 2p RDM need to be calculated, and false otherwise;
                                                        the 1pRDM is calculated by default.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    calcNatural                ::Bool 
    calcDensity                ::Bool 
    calcIpq                    ::Bool 
    calc2pRDM                  ::Bool 
    printBefore                ::Bool 
    levelSelection             ::LevelSelection
end 


"""
`ReducedDensityMatrix.Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false, 
                                    printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection()) 
    ... keyword constructor to overwrite selected value of reduced-density matrix computations.
"""
function Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false,
                    calc2pRDM::Bool=false, printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection())
    Settings(calcNatural, calcDensity, calcIpq, calc2pRDM, printBefore, levelSelection)
end

# `Base.show(io::IO, settings::ReducedDensityMatrix.Settings)`  
#       ... prepares a proper printout of the variable settings::ReducedDensityMatrix.Settings.
function Base.show(io::IO, settings::ReducedDensityMatrix.Settings) 
    println(io, "calcNatural:              $(settings.calcNatural)  ")
    println(io, "calcDensity:              $(settings.calcDensity)  ")
    println(io, "calcIpq:                  $(settings.calcIpq)  ")
    println(io, "calc2pRDM:                $(settings.calc2pRDM)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  ReducedDensityMatrix.Outcome`  
    ... defines a type to keep the outcome of a reduced density-matrix and natural-orbital computation.

    + level                     ::Level                 ... Atomic level to which the outcome refers to.
    + naturalSubshells          ::Array{Subshell,1}     ... List of natural orbitals (subshells).
    + naturalOccupation         ::Array{Float64,1}      ... Occupation numbers of natural orbitals.
    + naturalOrbitalExpansion   ::Dict{Subshell, Array{Float64,1}} 
        ... Dictionary of the expansion of (one-electron) natural orbitals in terms of the standard orbitals
            from the given basis.
    + naturalOrbitals           ::Dict{Subshell, Orbital} . Dictionary of (one-electron) natural orbitals.
    + orbitalInteraction        ::Array{Float64,2}      ... Orbital interaction I_pq = I[ip,iq].
    + rho1p                     ::Array{Float64,2}      ... One-particle RDM rho^(1p) [ip,iq].
    + rho2p                     ::Array{Float64,4}      ... Two-particle RDM rho^(1p) [ip,iq].
    + electronDensity           ::Radial.Density        ... Radial density distribution.
"""
struct Outcome 
    level                       ::Level 
    naturalSubshells            ::Array{Subshell,1}
    naturalOccupation           ::Array{Float64,1}
    naturalOrbitalExpansion     ::Dict{Subshell, Array{Float64,1}}
    naturalOrbitals             ::Dict{Subshell, Orbital}
    orbitalInteraction          ::Array{Float64,2}
    rho1p                       ::Array{Float64,2}
    rho2p                       ::Array{Float64,4}
    electronDensity             ::Radial.Density
end 


"""
`ReducedDensityMatrix.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of isotope-shift properties.
"""
function Outcome()
    Outcome(Level(), Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}, Dict{Subshell, Orbital}(), 
            zeros(1,1), zeros(1,1), zeros(1,1,1,1), Radial.Density() )
end


# `Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome)`  ... prepares a proper printout of the variable outcome::ReducedDensityMatrix.Outcome.
function Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome) 
    println(io, "level:                   $(outcome.level)  ")
    println(io, "naturalSubshells:        $(outcome.naturalSubshells)  ")
    println(io, "naturalOccupation:       $(outcome.naturalOccupation)  ")
    println(io, "naturalOrbitalExpansion: $(outcome.naturalOrbitalExpansion)  ")
    println(io, "naturalOrbitals:         $(outcome.naturalOrbitals)  ")
    println(io, "orbitalInteraction:      $(outcome.orbitalInteraction)  ")
    println(io, "rho1p:                   $(outcome.rho1p)  ")
    println(io, "rho2p:                   $(outcome.rho2p)  ")
    println(io, "electronDensity:         $(outcome.electronDensity)  ")
end


"""
`ReducedDensityMatrix.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::ReducedDensityMatrix.Settings; output=true)` 
    ... to compute (as selected) the 1-particle and 2-particle reduced density matrices, natural orbitals or other
        requested information for the levels of interest.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                            settings::ReducedDensityMatrix.Settings; output=true)
    println("")
    printstyled("ReducedDensityMatrix.computeOutcomes(): The computation of the reduced density matrices starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    outcomes = ReducedDensityMatrix.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    ReducedDensityMatrix.displayOutcomes(outcomes, settings)    end
    # Calculate all requested properties for each outcome
    newOutcomes = ReducedDensityMatrix.Outcome[]
    for  outcome in outcomes
        newOutcome = ReducedDensityMatrix.computeProperties(outcome, nm, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    ReducedDensityMatrix.displayResults(stdout, newOutcomes, grid, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    ReducedDensityMatrix.displayResults(iostream, newOutcomes, grid, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`ReducedDensityMatrix.computeNaturalOrbitalExpansion(rho1p::Array{Float64,2}, level::Level)`  
    ... to perform the expansion of the natural orbitals for the given level in terms of the standard orbitals
        as defined by the basis. A tuple (naturalOcc::Array{Float64,1}, naturalExp::Dict{Subshell, Array{Float64,1}})
        is returned which provides the natural occpuation numbers and expansion coefficients with regard
        to the given list of natural orbitals. This method makes use of the 1-particle RDM to extract and diagonalize the 
        expansion matrix.
"""
function  computeNaturalOrbitalExpansion(rho1p::Array{Float64,2}, level::Level)
    naturalOcc = Float64[];      naturalExp = Dict{Subshell, Array{Float64,1}}();       lenNO = length(level.basis.subshells)
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", rho1p)
    orbIndices = Float64[]
    for  vector in eigen.vectors    wx = findmax(vector);   push!(orbIndices, wx[2])    end
    @show "  "
    @show orbIndices
    for i = 1:lenNO     
        wx = findfirst(x->x==i, orbIndices)
    end
    
    return( naturalOcc, naturalExp )
end


"""
`ReducedDensityMatrix.computeNaturalOrbitals(naturalSubshells::Array{Subshell,1}, naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)`  
    ... to compute the natural orbitals as superposition of the standard orbitals; for each subshell in naturalExp,
        an orbital is computed, and naturalOrbitals::Dict{Subshell, Orbital} returned
"""
function  computeNaturalOrbitals(naturalSubshells::Array{Subshell,1}, naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)
    naturalOrbitals = Dict{Subshell, Orbital}()
    println(">> computeNaturalOrbitals() ... not yet implemented !!")

    return( naturalOrbitals )
end


"""
`ReducedDensityMatrix.computeOrbitalInteractions(naturalSubshells::Array{Subshell,1}, naturalOrbitals::Dict{Subshell, Orbital}, level::Level)`  
    ... to compute the natural orbitals as superposition of the standard orbitals; for each subshell in naturalExp,
        an orbital is computed, and naturalOrbitals::Dict{Subshell, Orbital} returned
"""
function  computeOrbitalInteractions(naturalSubshells::Array{Subshell,1}, naturalOrbitals::Dict{Subshell, Orbital}, level::Level)
    ns = length(naturalSubshells)
    orbitalInteraction = zeros(ns,ns)    
    println(">> computeOrbitalInteractions() ... not yet implemented !!")

    return( orbitalInteraction )
end


"""
`ReducedDensityMatrix.computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::ReducedDensityMatrix.Settings) 
    ... to compute all properties for a given level; an outcome::ReducedDensityMatrix.Outcome is returned for 
        which the properties are now evaluated explicitly.
"""
function  computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::ReducedDensityMatrix.Settings)
    naturalOccupation = outcome.naturalOccupation;        naturalOrbitalExp  = outcome.naturalOrbitalExpansion;  
    naturalOrbitals   = outcome.naturalOrbitals;          orbitalInteraction = outcome.orbitalInteraction;       
    electronDensity   = outcome.electronDensity;          rho1p              = outcome.rho1p;             rho2p = outcome.rho2p
    
    naturalSubshells   = copy(outcome.level.basis.subshells)
    
    if  settings.calcNatural
        rho1p                                = compute1pRDM(outcome.level)
        naturalOccupation, naturalOrbitalExp = computeNaturalOrbitalExpansion(rho1p, outcome.level)
        naturalOrbitals                      = computeNaturalOrbitals(naturalSubshells, naturalOrbitalExp, outcome.level)
        # For test purposes
        naturalOccupation                    = ones( length(naturalSubshells) )
        naturalOrbitalExp[Subshell(1,-1)]    = zeros( length(naturalSubshells) )
        naturalOrbitalExp[Subshell(2,-1)]    = ones(  length(naturalSubshells) )
    end

    if  settings.calcDensity
        electronDensity = computeRadialDistribution(naturalOrbitalExp, grid, outcome.level)
        # For test purposes
        electronDensity = Radial.Density("Radial dist", ones(grid.NoPoints), grid)
    end

    if  settings.calcIpq
        orbitalInteraction = computeOrbitalInteractions(naturalSubshells, naturalOrbitals, outcome.level)
    end

    if  settings.calc2pRDM
        rho2p = compute2pRDM(outcome.level)
    end

    newOutcome = ReducedDensityMatrix.Outcome( outcome.level, naturalSubshells,  naturalOccupation, naturalOrbitalExp,
                                                naturalOrbitals, orbitalInteraction, rho1p, rho2p, electronDensity )
    return( newOutcome )
end


"""
`ReducedDensityMatrix.computeRadialDistribution(naturalExp::Dict{Subshell, Array{Float64,1}}, grid::Radial.Grid, level::Level)`  
    ... to compute the radial distribution function from the given standard orbitals and natural expansion
        coefficients. An radial disribution D::Array{Float64,1} is returned that specifies the radial distribution
        w.r.t. the given grid.
"""
function  computeRadialDistribution(naturalExp::Dict{Subshell, Array{Float64,1}}, grid::Radial.Grid, level::Level)
    D = Radial.Density()
    println(">> computeRadialDistribution() ... not yet implemented !!")

    return( D )
end


"""
`ReducedDensityMatrix.compute1pRDM(level::Level)`  
    ... to compute 1p RDM; a rdm::Array{Float64,2} is returned.
"""
function  compute1pRDM(level::Level)
    subshellList = level.basis.subshells
    lenNO = length(subshellList);    rho_pq     = zeros(lenNO,lenNO)
    # Cycle over all matrix elements of the CSF basis
    opa          = SpinAngular.OneParticleOperator(0, plus, true)
    for  (ir, rcsf) in enumerate(level.basis.csfs)
        for  (is, scsf) in enumerate(level.basis.csfs)
            # Calculate angular coefficient for rank-0 operator
            wa = SpinAngular.computeCoefficients(opa, rcsf, scsf, subshellList) 
            # Cycle over the pair of natural subshells in rho_pq
            for (ip,p)  in  enumerate(subshellList)
                for (iq,q)  in  enumerate(subshellList)
                    for  coeff in wa
                        if  (p == coeff.a   &&  q == coeff.b)  ||  (p == coeff.b   &&  q == coeff.a)
                            jj = Basics.subshell_2j(level.basis.orbitals[coeff.a].subshell)
                            rho_pq[ip,iq] = rho_pq[iq,ip] = rho_pq[ip,iq] + level.mc[ir] * coeff.T * sqrt( jj + 1) * level.mc[is]
                        end
                    end
                end
            end
        end
    end

    return( rho_pq )
end


"""
`ReducedDensityMatrix.compute2pRDM(level::Level)`  
    ... to compute 2p RDM for all pairs (p,q;r,s) of subshells; a rdm::Array{Float64,4} is returned.
"""
function  compute2pRDM(level::Level)
    subshellList = level.basis.subshells
    lenNO = length(subshellList);    rdm     = zeros(lenNO,lenNO,lenNO,lenNO)
    # Cycle over all matrix elements of the CSF basis
    opa = SpinAngular.TwoParticleOperator(0, plus, true)
    for  (i, icsf) in enumerate(level.basis.csfs)
        for  (j, jcsf) in enumerate(level.basis.csfs)
            # Calculate angular coefficient for rank-0 operator
            wa = SpinAngular.computeCoefficients(opa, icsf, jcsf, subshellList) 
            # Cycle over the pair of natural subshells in rho_pqrs
            for (ip,p)  in  enumerate(subshellList)
                for (iq,q)  in  enumerate(subshellList)
                    for (ir,r)  in  enumerate(subshellList)
                        for (is,s)  in  enumerate(subshellList)
                            for  coeff in wa
                                if  (p == coeff.a   &&  q == coeff.b   &&  r == coeff.c   &&  s == coeff.d)  ||  
                                    (p == coeff.c   &&  q == coeff.d   &&  r == coeff.a   &&  s == coeff.b)
                                    rdm[ip,iq,ir,is] = rdm[ir,is,ip,iq] = rdm[ip,iq,ir,is] + level.mc[i] * coeff.V * level.mc[j]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return( rdm )
end


"""
`ReducedDensityMatrix.determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings)`  
    ... to determine a list of Outcomes's for the computation of reduced density matrices and other information 
        for the given multiplet. It takes into account the particular selections and settings. 
        An Array{ReducedDensityMatrix.Outcome,1} is returned. Apart from the level specification, all physical 
        properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings) 
    outcomes = ReducedDensityMatrix.Outcome[]
    @show settings.levelSelection
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, ReducedDensityMatrix.Outcome(level, Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}(), 
                                    Dict{Subshell, Orbital}(), zeros(1,1), zeros(1,1), zeros(1,1,1,1), Radial.Density()) )
        end
    end
    return( outcomes )
end


"""
`ReducedDensityMatrix.displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise. Moreover, the
        selected properties and information is printed as well, though not yet calculated.
"""
function  displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)
    nx = 43
    println(" ")
    println("  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    println("\n\n  The following properties are calculated: \n")
    println("  + Natural occupations numbers")
    println("  + Natural orbital expansion in terms of the standard orbitals.")
    if  settings.calcNatural   println("  + Representation of the natural orbitals.")  end
    if  settings.calcDensity   println("  + Radial electron density distribution.")    end
    if  settings.calcIpq       println("  + Orbital interaction I_pq.")                end
    
    return( nothing )
end


"""
`ReducedDensityMatrix.displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, grid::Radial.Grid, 
                                        settings::ReducedDensityMatrix.Settings)`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise. Moreover, the
        selected properties and information is printed as well, though not yet calculated.
"""
function  displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, grid::Radial.Grid, 
                            settings::ReducedDensityMatrix.Settings)
    nx = 43
    println(stream, " ")
    println(stream, "  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    # Now print all selected result for each outcome
    for  outcome in outcomes
        nx = 120;       sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity) 
        println(stream, " ")
        println(stream, " ")
        println(stream, "  =============================")
        println(stream, "  Level: $(outcome.level.index) with symmetry $sym  ")
        println(stream, "  =============================")
        #
        #
        println(stream, " ")
        println(stream, "  1-particle RDM:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
        for  (i, sa)  in  enumerate(wa)
            if  i == 1      sb = "     Orb | Orb  " * sa
            else            sb = "                " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  (i, subsh)  in  enumerate(outcome.level.basis.subshells)
            sb = string(subsh)
            wa = TableStrings.floatList(10, outcome.rho1p[i,:])
            for  (j, sa)  in  enumerate(wa)
                if       j == 1  sb = "     " * sb * "     " * sa
                else             sb = "                    " * sa      end
                println(stream,  sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        #
        println(stream, " ")
        println(stream, "  Natural occupation numbers:  ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.naturalSubshells)
        for  (i, sa)  in  enumerate(wa)
            if       i == 1  sb = "     NO   " * sa
            else             sb = "          " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        wa = TableStrings.floatList(10, outcome.naturalOccupation)
        for  (i, sa)  in  enumerate(wa)
            println(stream,  "          " * sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        #
        nx = 120
        println(stream, " ")
        println(stream, "  Natural orbital expansion:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.naturalSubshells)
        for  (i, sa)  in  enumerate(wa)
            if       i == 1  sb = "     Orbitals: " * sa
            else             sb = "               " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  (i, subsh)  in  enumerate(outcome.naturalSubshells)
            sb = string(subsh)
            if  i > 2   println(stream, "     " * sb * "    test ...");    break     end
            wa = TableStrings.floatList(10, outcome.naturalOrbitalExpansion[subsh])
            for  (j, sa)  in  enumerate(wa)
                if       j == 1  sb = "     " * sb * "   " * sa
                else             sb = "                  " * sa      end
                println(stream,  sb)
            end
            end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        #
        if  settings.calcNatural
            println(stream, "\n  Natural orbitals are calculated for the following subshells and are available by outcome.naturalOrbitals:\n")
            wa = TableStrings.subshellList(10, outcome.naturalSubshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     NO:  " * sa
                else             sb = "          " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        #
        if  settings.calcDensity
            println(stream, "\n  Radial electron distribution:\n")
            nx = 42;    sa = "         i)         r[i]       density[i]"
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream,  sa)
            println(stream, "  ", TableStrings.hLine(nx))
            for  i = 1:30:outcome.electronDensity.grid.NoPoints
                println(stream, "   " * ("      "*string(i))[end-6:end] * ")      " *
                                @sprintf("%.4e", grid.r[i]) * "    " * @sprintf("%.4e", outcome.electronDensity.Dr[i])  )
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        #
        if  settings.calcIpq
            nx = 120
            println(stream, "\n  Orbital interactions I_pq:\n")
            println(stream, "  ", TableStrings.hLine(nx))
            wa = TableStrings.subshellList(10, outcome.naturalSubshells)
            for  (i, sa)  in  enumerate(wa)
                if  i == 1      sb = "     Orb | Orb  " * sa
                else            sb = "                " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))
            #
            for  (i, subsh)  in  enumerate(outcome.naturalSubshells)
                sb = string(subsh)
                wa = TableStrings.floatList(10, outcome.orbitalInteraction[i,:])
                for  (j, sa)  in  enumerate(wa)
                    if       j == 1  sb = "     " * sb * "     " * sa
                    else             sb = "                    " * sa      end
                    println(stream,  sb)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        #
        if  settings.calc2pRDM
            println(stream, "\n  Reduced 2-particle RDM were calculated for the following subshells and are available by outcome.rho2p:\n")
            wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     Orb:  " * sa
                else             sb = "           " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))
            #
            # Compute the diagonal part of rho2p
            lenNO = length(outcome.level.basis.subshells);   rho_pprr = zeros(lenNO,lenNO)
            for  ip = 1:lenNO   for  ir = 1:lenNO     rho_pprr[ip,ir] = outcome.rho2p[ip,ip, ir,ir]     end   end
            #
            println(stream, "\n  The diagonal part of the 2-particle RDM refers the mean (product) of occuaption numbers " *
                            "rho_pp,rr = (n_p * n_r)_av:\n")
            wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     Orb:  " * sa
                else             sb = "           " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))
            #
            for  (i, subsh)  in  enumerate(outcome.level.basis.subshells)
                sb = string(subsh)
                wa = TableStrings.floatList(10, outcome.orbitalInteraction[i,:])
                for  (j, sa)  in  enumerate(wa)
                    if       j == 1  sb = "     " * sb * "     " * sa
                    else             sb = "                    " * sa      end
                    println(stream,  sb)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
    end # outcomes
    #
    return( nothing )
end

end # module
