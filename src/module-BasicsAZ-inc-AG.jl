
export recastAG

"""
`Basics.add(pota::Radial.Potential, potb::Radial.Potential)`  
    ... to add two radial potentials together that are defined on the same grid. A potential::RadialPotential is returned 
        that inherits its radial size from the potential that is defined in a larger range of r-values.
"""
function  Basics.add(pota::Radial.Potential, potb::Radial.Potential)
    if  pota.grid.NoPoints != potb.grid.NoPoints  ||   pota.grid.rnt != potb.grid.rnt  ||  pota.grid.h != potb.grid.h  ||  
        pota.grid.hp != potb.grid.hp                   error("stop a")   end
    name = pota.name * "+" * potb.name;   nmax = max(length(pota.Zr), length(potb.Zr));   nmin = min(length(pota.Zr), length(potb.Zr))
    Zr   = zeros(nmax);   
    nx = length(pota.Zr);    Zr[1:nx] = Zr[1:nx] + pota.Zr[1:nx] 
    nx = length(potb.Zr);    Zr[1:nx] = Zr[1:nx] + potb.Zr[1:nx]
    Zr[nmin+1:nmax] .= Zr[nmin] 

    potential = Radial.Potential(name, Zr, pota.grid)
    return( potential )
end



"""
`Basics.addZerosToCsfR(nz::Int64, csf::CsfR)`  
    ... 'adds' a number of zeros to the occupation, seniorityNr, etc. of a given CsfR.
"""
function  Basics.addZerosToCsfR(nz::Int64, csf::CsfR)
    !csf.useStandardSubshells   &&   error("Zeros can only be added to CSF with standard subshell order.")
    
    occupation = csf.occupation;   seniorityNr = csf.seniorityNr;    subshellJ = csf.subshellJ;    subshellX = csf.subshellX 
    
    for  i = 1:nz
        push!(occupation, 0);    push!(seniorityNr, 0);    push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellJ[end])
    end
    
    newCsf = CsfR( csf.useStandardSubshells, csf.J, csf.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[])
    return(newCsf)
end



"""
`Basics.analyze(::LevelDecompositionOfNRconfigurations, level::Level)`
    ... analyzes and lists the non-relativistic configurations with a weight larger than 5%; a list of NR
        configurations with their corresponding weights is printed, but nothing is returned.
"""
function Basics.analyze(::LevelDecompositionOfNRconfigurations, level::Level)
    !(level.hasBasis)   &&   error("Levels without a basis cannot be analyzed.")

    confList   = generate("configuration list: NR, from basis", level.basis)
    percentage = zeros( length(confList) )

    # Now generate for each CSF in this basis a Configuration, compare with confList and add contribution to the corresponding percentage
    for  k = 1:length(level.basis.csfs)
        csf = level.basis.csfs[k]
        shellList = Dict{Shell,Int64}[]
        for  i = 1:length(csf.occupation)
            n = csf.subshell[i].n;    l = Basics.subshell_l( csf.subshell[i] );    sh = Shell(n, l)
            if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
            else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
            end
        end
        confNew = Configuration( shellList, NoElectrons )
        # Compare and add the weight to the right configuration
        for  i = 1:length(confList)
            if   confList[i] == confNew    percentage[i] = percentage[i] + level.mc[k]^2;    break    end
        end
    end

    # Check that all weights sum up properly
    wa = sum( percentage );   abs(1. - wa) > 1.0e-3  &&   error("Total percentage $wa must add to 1.")
    # Sort the percentage and print both, configurations and percentage
    wb = sortperm( percentage )
    for i in wb
        println("Configuration: $(confList[i])  ... with  $(percentage[i]) %")
    end

    return( nothing )
end


"""
`Basics.analyze(::LevelDecompositionOfCsfR, level::Level, N::Int64)`
    ... analyzes and lists (up to) N relativistic CSF together with their weight |c_n|^2 in the expansion of
        the given level; a list of CSF and their corresponding weights is printed, but nothing is returned.
        **Not yet implemented.**
"""
function Basics.analyze(::LevelDecompositionOfCsfR, level::Level, N::Int64)
    error("Not yet implemented !")

    return( nothing )
end



"""
`Basics.analyze(orbitals::Dict{Subshell, Orbital}; printout::Bool=false)`  
    ... to anaylze the phase of the orbitals and to ensure that all orbitals start `positve' near r --> 0; 
        a set of newOrbitals::Dict{Subshell, Orbital} is returned.
"""
function Basics.analyze(orbitals::Dict{Subshell, Orbital}; printout::Bool=false)
    newOrbitals = Dict{Subshell, Orbital}();    mOrbital = Orbital(Subshell("1s_1/2"), -1.0)
    # Analyze and modify the orbital if necessary
    for (k,v) in orbitals
        found = false
        for  P in v.P   
            if      abs(P) > 1.0e-5  &&  abs(P) == P
                mOrbital = v;   found = true;   break
            elseif  abs(P) > 1.0e-5  &&  abs(P) == -P
                mOrbital = Orbital( v.subshell, v.isBound, v.useStandardGrid, v.energy, -1.0 .* v.P,  -1.0 .* v.Q,  
                                    -1.0 .* v.Pprime,  -1.0 .* v.Qprime, v.grid)
                if  printout    println(">>> Sign changed for orbital $k")   end                    
                found = true;   break
            end
        end
        if  !found  error("stop a")     end
        newOrbitals[k] = mOrbital
    end

    return( newOrbitals )
end


    
"""
`Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)`  
    ... to analyze the convergence for two blocks of eigenvalues und eigenvectors; the value > 0 indicates the relative maximum 
        difference abs( a-b/a+b ) between two eigenvalues. The value is set arbitrarely set to 1.0e+2 if the number of eigenvalues
        differ in the two blocks.
"""
function Basics.analyzeConvergence(blocka::Basics.Eigen, blockb::Basics.Eigen)
    va = blocka.values;   vb = blockb.values
    if       length(va) !=  length(vb)    wa = 1.0e+2   
    else     wa = 0.;       for  i = 1:length(va)   wa = max(wa, abs( (va[i]-vb[i]) / (va[i]+vb[i]) ))   end
    end
    return( wa )
end


    
"""
`Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)`  
    ... to analyze the convergence for two radial orbitals; the value > 0 indicates the relative maximum difference abs( a-b/a+b ) 
        between the energies or (large or small) components of the orbitals. The function assumes that the orbitals are represented 
        on the same grid and that the components are zero if only one of the orbitals is defined beyond some maximum grid point.
"""
function Basics.analyzeConvergence(a::Radial.Orbital, b::Radial.Orbital)
    if       a.subshell !=  b.subshell                                             error("stop a")   
    elseif   a.useStandardGrid !=  b.useStandardGrid   ||  !(a.useStandardGrid)    @show a.useStandardGrid, b.useStandardGrid;   error("stop b")   
    else     wa = 0.;   wa = max(wa, abs( (a.energy-b.energy) / (a.energy+b.energy) ));   nx = min(length(a.P), length(b.P))
            #=== for  i = 1:nx 
                wp = abs(a.P[i])+abs(b.P[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.P[i]-b.P[i]) / wp ))   end 
                wp = abs(a.Q[i])+abs(b.Q[i]);    if  wp > 1.0e-10   wa = max(wa, abs( (a.Q[i]-b.Q[i]) / wp ))   end   
            end ===#
    end
    return( wa )
end



    
"""
`Basics.analyzeGrid(shells::Array{Shell,1}, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to analyze the representation of hydrogenic orbitals from shells for the given nuclear charge and grid.
        The procedure orders the shells, tabulates their (hydrogenic) energy in the given grid and compares them
        with the exact solution. Nothing is returned otherwise..
"""
function Basics.analyzeGrid(shells::Array{Shell,1}, nm::Nuclear.Model, grid::Radial.Grid)
    subshells = Basics.generateSubshellList(shells::Array{Shell,1})
    ns = length(subshells);    sortedSubshells = sort(subshells);    isCalculated = falses(ns)
    sortedLines = String[];    for n = 1:ns   push!(sortedLines, "aa")   end
    #
    # Calculate the hydrogenic orbitals for all requested shells and take over the results for this shell
    # as well as all other shells of the same symmetry
    wa = Bsplines.generatePrimitives(grid)
    wb = Bsplines.generateOrbitalsHydrogenic(wa, nm, subshells, printout=true) 
    #
    # Print the table in a neat format
    nx = 67
    println("  ", TableStrings.hLine(nx))
    sa = "  ";      
    sa = sa * TableStrings.center(10, "Subshell";         na=3)
    sa = sa * TableStrings.center(17, "Energies [a.u.]";  na=2)
    sa = sa * TableStrings.center(17, "Dirac-E  [a.u.]";  na=2)
    sa = sa * TableStrings.center(17, "Delta-E / |E|";    na=2)
    println(sa)
    #
    for subsh in subshells
        sa = "" *  TableStrings.flushright(10, string(subsh); na=6)
        en = Basics.computeDiracEnergy(subsh, nm.Z)
        sa = sa * @sprintf("%.8e", wb[subsh].energy)                 * "    "
        sa = sa * @sprintf("%.8e", en)                               * "    "
        sa = sa * @sprintf("%.8e", (wb[subsh].energy-en)/abs(en))
        println(sa)
    end
    println("  ", TableStrings.hLine(nx))
    
    return( nothing )
end

    

"""
`Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)`  
    ... to determine the NoEnergySharings sharings by using the Gauss-Legendre integration points in the interval
        [0., energy]; an Array{Tuple{Float64,Float64,Float64},1} is returned. This methods assumes that the package 
        GaussQuadrature is 'used'.
"""
function Basics.determineEnergySharings(energy::Float64, NoEnergySharings::Int64)
    xx, ww = GaussQuadrature.legendre(NoEnergySharings)
    sharings = Tuple{Float64,Float64,Float64}[]
    for  (ix, x) in enumerate(xx)
        xrel = (x+1.) / 2
        push!(sharings, (xrel*energy, energy - xrel*energy, ww[ix]*energy) )
    end
    return( sharings )
end



"""
`Basics.determineHoleShells(conf::Configuration)`  
    ... to determine the hole-shells of a given non-relativistic configuration; a shellList::Array{Shell,1} is returned.
"""
function Basics.determineHoleShells(conf::Configuration)

    shellList = Shell[]
    for  (k,v) in conf.shells
    if   conf.shells[k] < 2*(2k.l+1)   push!(shellList, k)    end 
    end

    return( shellList )
end



"""
`Basics.determineMeanEnergy(conf::Configuration, orbitals::Dict{Subshell, Orbital}, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to determine the mean energy of a given non-relativistic configuration; this 'mean' energy is calculated as the mean
        energy of all single-CSF levels and by using the the given set of orbitals. A value::Float64 is returned.
"""
function Basics.determineMeanEnergy(conf::Configuration, orbitals::Dict{Subshell, Orbital}, nm::Nuclear.Model, grid::Radial.Grid)
    asfSettings = AsfSettings();    meanEnergy = 0.;    nlev = 0
    
    multiplet   = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, nm, grid, asfSettings, printout=false)
    for  level in multiplet.levels
        nl = Basics.twice(level.J) + 1;     nlev = nlev + nl;     meanEnergy = meanEnergy + nl * level.energy
    end  
    meanEnergy = meanEnergy / nlev
    
    if true  println(">>> mean energy of $conf  is:  $meanEnergy  [a.u.]")   end
    
    return( meanEnergy )
end


"""
`Basics.determineMeanEnergy(multiplet::Multiplet)`  
    ... to determine the mean energy of a given multiplet. A value::Float64 is returned.
"""
function Basics.determineMeanEnergy(multiplet::Multiplet)
    meanEnergy = 0.;    nlev = 0
    for  level in multiplet.levels
        nl = Basics.twice(level.J) + 1;     nlev = nlev + nl;     meanEnergy = meanEnergy + nl * level.energy
    end  
    meanEnergy = meanEnergy / nlev
    
    if false  println(">>> mean energy of $(nlev)-level multiplet is:  $meanEnergy  [a.u.]")   end
    
    return( meanEnergy )
end


"""
`Basics.determineNearestPoints(x0::Float64, n::Int64, values::Array{Float64,1})`  
    ... to determine the n points in values that are nearest to x0;  
        two (short) arrays ws::Array{Float64,1}, diffs::Array{Float64,1} of length  n   are returned that contain the 
        nearest values and differences (value - x0) from x0. An error message is issued of n > length(values).
"""
function Basics.determineNearestPoints(x0::Float64, n::Int64, values::Array{Float64,1})
    valuesx = copy(values)
    @show n, values
    if  n > length(values)
        println(">>> Warning in Basics.determineNearestPoints():  n=$n > No of $values ")
        for  i=1:100   push!(valuesx, values[end] + i*0.1 );   if length(valuesx) == n   break   end     end
    end
    ws = zeros(n);      diffs = zeros(n)
    dabs = Float64[];     for  v in valuesx   push!(dabs, abs(v-x0))   end
    for  i=1:n 
        ix = findmin(dabs)[2];    ws[i] = valuesx[ix];    diffs[i] = valuesx[ix] - x0
        dabs[ix] = 1.0e99
    end
    
    if  length(unique(ws)) != n   error("stop a")    end 
    
    return( ws, diffs )
end


"""
`Basics..determineNonorthogonalShellOverlap(finalSubshells::Array{Subshell,1}, initialSubshells::Array{Subshell,1}, 
                                            grid::Radial.Grid)`  
    ... determines the overlap of those shells that do not agree in both sets but have the same symmetry; a zero
        overlap and Subshell(99,-1) are returned in cases that no such pair can be found.
        A triple (shella, shellb, overlap) is returned.
"""
function Basics.determineNonorthogonalShellOverlap(finalSubshells::Array{Subshell,1}, initialSubshells::Array{Subshell,1}, 
                                                    grid::Radial.Grid)
    neSubshells = setdiff(union(finalSubshells, initialSubshells),  intersect(finalSubshells, initialSubshells))
    ret = (Subshell(99,-1), Subshell(99,-1), 0.);     ne = length(neSubshells)
    for  a = 1:ne
        for  b = a+1:ne
            if   neSubshells[a].kappa == neSubshells[b].kappa   
                println(">>>> non-zero overlap integral < $(neSubshells[a]) | $(neSubshells[b]) > = 0.1  ... but incorrect !!!")
                ret = (neSubshells[a], neSubshells[b], 0.1)
            end
        end
    end
    
    if  ret[3] == 0.  println(">>>> zero overlap integral for $neSubshells ")  end
    
    return( ret )
end


"""
`Basics.determinePolarizationVector(q::Int64, kind::Basics.AbstractPolarization; star::Bool)`  
    ... determines the q-th component of the (complex and unit) polarization vector u for the given kind of polarization;
        for star = true, the complex-conjugate component is returned.  a comp::Float64 is returned.
"""
#---OLD VERSIONS:---
#function Basics.determinePolarizationVector(q::Int64, kind::Basics.LinearPolarization; star::Bool=false)
#    if      q == 1   return( -1/sqrt(2.) )
#    elseif  q == -1  return( 1/sqrt(2.) )
#    else   return( 0 )
#    end
#end
#
#function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftCircular; star::Bool=false)
#    if  q == 1   return( -1. )   else   return( 0. )                end
#end

#function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightCircular; star::Bool=false)
#    if  q == -1  return( 1. )   else   return( 0. )                end
#end

#---NEW VERSIONS:---
function Basics.determinePolarizationVector(q::Int64, kind::Basics.LinearPolarization; star::Bool)
    if  star == false
        if      q == 1   return( -1/sqrt(2.) )
        elseif  q == -1  return( 1/sqrt(2.) )
        else   return( 0 )
        end
    else
        if      q == -1   return( 1/sqrt(2.) )
        elseif  q == 1  return( -1/sqrt(2.) )
        else   return( 0 )
        end
    end
end
#
function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftCircular; star::Bool)
    if  star == false
        if  q == 1   return( -1. )   else   return( 0. )                end
    else
        if  q == -1   return( 1. )   else   return( 0. )                end
    end
end
#
function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightCircular; star::Bool)
    if  star == false
        if  q == -1  return( 1. )   else   return( 0. )                end
    else
        if  q == 1  return( -1. )   else   return( 0. )                end
    end
end
#
function Basics.determinePolarizationVector(q::Int64, kind::Basics.LeftElliptical; star::Bool)
    epsilon = kind.ellipticity
    if  star == false
        if  q == 1   return( -1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( 1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )                end
    else
        if  q == -1   return( 1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( -1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )               end
    end
end
#
function Basics.determinePolarizationVector(q::Int64, kind::Basics.RightElliptical; star::Bool)
    epsilon = kind.ellipticity
    if  star == false
        if  q == -1  return( 1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( -1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )               end
    else
        if  q == 1  return( -1/sqrt(2*(1. + epsilon^2)) * (1. - epsilon) )   else   return( 1/sqrt(2*(1. + epsilon^2)) * (1. + epsilon) )                end
    end
end

"""
`Basics.determinePolarizationLambda(kind::Basics.AbstractPolarization)`  
    ... determines the (standard) lambda value for the given kind of polarization; a lambda::Int64 is returned.
"""
function Basics.determinePolarizationLambda(kind::Basics.LeftCircular)   return(  1 )    end
function Basics.determinePolarizationLambda(kind::Basics.RightCircular)  return( -1 )    end
function Basics.determinePolarizationLambda(kind::Basics.LeftElliptical)   return(  1 )    end
function Basics.determinePolarizationLambda(kind::Basics.RightElliptical)  return( -1 )    end


"""
`Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)`  
    ... to determine the specified lines as tupels of initial- and final levels. A level index is 0 in lineList always refers to
        all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64},1} is returned.
"""
function Basics.determineSelectedLines(lineList::Array{Tuple{Int64,Int64},1}, initialMultiplet::Multiplet, finalMultiplet::Multiplet)
    newList = Tuple{Int64,Int64}[]
    for tupl in lineList
        if      tupl[1] >  length(initialMultiplet.levels)   error("Initial level index $(tupl[1]) must be <= No. initial levels.")
        elseif  tupl[2] >  length(finalMultiplet.levels)     error("Final level index $(tupl[2]) must be <= No. final levels.")
        elseif  tupl[1] == 0  &&   tupl[2] == 0    error("Tupel (0,0) is not allowed.")
        elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)   push!(newList, (i, tupl[2]))   end
        elseif  tupl[2] == 0       for i = 1:length(finalMultiplet.levels)     push!(newList, (tupl[1], i))   end
        else                                                                   push!(newList, (tupl[1], tupl[2]))
        end
    end
    newList = unique(newList)
    return( newList)
end     


"""
`Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                    finalMultiplet::Multiplet)`  
    ... to determine the specified pathways as tupels of initial-, intermediate- and final levels. A level index 0 in pathwayList always 
        refers to all levels from the corresponding multiplet. A (unique) newList::Array{Tuple{Int64,Int64,Int64},1} is returned.
"""
function Basics.determineSelectedPathways(pathwayList::Array{Tuple{Int64,Int64,Int64},1}, initialMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                finalMultiplet::Multiplet)
    newList = Tuple{Int64,Int64,Int64}[]
    for tupl in pathwayList
        if      tupl[1] >  length(initialMultiplet.levels)       error("Initial level index $(tupl[1]) must be <= No. initial levels.")
        elseif  tupl[2] >  length(intermediateMultiplet.levels)  error("Intermediate level index $(tupl[2]) must be <= No. intermediate levels.")
        elseif  tupl[3] >  length(finalMultiplet.levels)         error("Final level index $(tupl[3]) must be <= No. final levels.")
        elseif  (tupl[1] + tupl[2] == 0  || tupl[1] + tupl[3] == 0  ||  tupl[2] + tupl[3] == 0)   &&  error("Only tupels with one 0-index are allowed.")
        elseif  tupl[1] == 0       for i = 1:length(initialMultiplet.levels)          push!(newList, (i, tupl[2], tupl[3]))   end
        elseif  tupl[2] == 0       for i = 1:length(intermediateMultiplet.levels)     push!(newList, (tupl[1], i, tupl[3]))   end
        elseif  tupl[3] == 0       for i = 1:length(finalMultiplet.levels)            push!(newList, (tupl[1], tupl[2], i))   end
        else                                                                          push!(newList, (tupl[1], tupl[2], tupl[3]))
        end
    end
    newList = unique(newList)
    return( newList)
end     



"""
`Basics.diagonalize(::MatrixWithLinearAlgebra, matrixA::Array{Float64,2}; range=(0:0)::UnitRange{Int64})`
    ... applies LinearAlgebra.eigen() to the symmetric matrix matrixA; only the upper-triangular part is used
        by an explicit symmetrization; an eigen::Basics.Eigen is returned.
"""
function Basics.diagonalize(::MatrixWithLinearAlgebra, matrixA::Array{Float64,2}; range=(0:0)::UnitRange{Int64})
    mA = LinearAlgebra.Symmetric(matrixA)
    if  range == 0:0
        wa = LinearAlgebra.eigen( mA )
        vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
        for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end
    else
        wa = LinearAlgebra.eigen( mA, range )
        vectors = Vector{Float64}[];    for  i = range   push!( vectors, wa.vectors[:,1] )   end
    end

    return( Basics.Eigen( wa.values, vectors ) )
end


"""
`Basics.diagonalize(::GeneralizedEigenvaluesWithLinearAlgebra, matrixA::Array{Float64,2}, matrixB::Array{Float64,2})`
    ... applies LinearAlgebra.eigen() to the generalized eigenvalue problem (matrixA, matrixB); only the
        upper-triangular parts of both matrices are used by an explicit symmetrization;
        an eigen::Basics.Eigen is returned.
"""
function Basics.diagonalize(::GeneralizedEigenvaluesWithLinearAlgebra, matrixA::Array{Float64,2}, matrixB::Array{Float64,2})
    mA = LinearAlgebra.Symmetric(matrixA)
    mB = LinearAlgebra.Symmetric(matrixB)
    wa = LinearAlgebra.eigen(mA, mB)
    vectors = Vector{Float64}[];    wb = wa.vectors;    d = size(wb)[1]
    for  i = 0:d-1    push!(vectors, wb[i*d+1:i*d+d])    end

    return( Basics.Eigen( wa.values, vectors ) )
end


"""
`Basics.diracDelta(x::Float64, dx::Float64)`
    ... evaluates Dirac's function  delta(x) = 0  for abs(x) > dx/2   and   delta(x) = 1/dx    for abs(x) <= dx/2;
        a values::Float64 is returned.
"""
function Basics.diracDelta(x::Float64, dx::Float64)    
    if  abs(x) > dx/2   return( 0. )    else    return( 1/dx )  end
end


"""
`Basics.display(::PhysicalConstants)`
    ... to display (all) currently defined physical constants; nothing is returned. Cf. Defaults.setDefaults().
"""
function Basics.display(::PhysicalConstants)
    println("Physical constants are defines as follows:  \n",
            "------------------------------------------  \n")
    sb = "  + Fine-structure constant:    " * string( Defaults.getDefaults("alpha") );                println(sb)
    sb = "  + Electron mass [kg]:         " * string( Defaults.getDefaults("electron mass: kg") );    println(sb)
    sb = "  + Electron mass [amu]:        " * string( Defaults.getDefaults("electron mass: amu") );   println(sb)
    println()

    return( nothing )
end


"""
`Basics.display(::CurrentSettings)`
    ... to display (all) currently defined settings of the JAC module; nothing is returned.
"""
function Basics.display(::CurrentSettings)
    println("Current settings of the JAC module:  \n",
            "-----------------------------------  \n")
    sb = "  + Framework:                              " * string( Defaults.getDefaults("framework") );                println(sb)
    sb = "  + Energy unit:                            " * string( Defaults.getDefaults("unit: energy") );             println(sb)
    sb = "  + Rate and transition probability unit:   " * string( Defaults.getDefaults("unit: rate") );               println(sb)
    sb = "  + Cross section unit:                     " * string( Defaults.getDefaults("unit: cross section") );      println(sb)
    sb = "  + Time unit:                              " * string( Defaults.getDefaults("unit: time") );               println(sb)
    println()
    if      Defaults.getDefaults("standard grid") != false    println("  + A standard grid has been defined; cf. Defaults.getDefaults()" )
    elseif  Defaults.getDefaults("standard grid") == false    println("  + No standard grid has yet been defined; cf. Defaults.setDefaults()" )
    end
    println()

    return( nothing )
end


"""
`Basics.display(stream::IO, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; longTable::Bool=false)`  
    ... displays the (generated)  orbitals in a neat and compact form; the default (longTable==false) just prints the subshell,
        isbound and energies. If longTable==true, it prints in addition the r_max, <r>, ... values.  nothing is returned.
"""
function  Basics.display(stream::IO, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid; longTable::Bool=false)
    println(stream, " ")
    println(stream, "  Relativistic orbitals:")
    println(stream, " ")
    if  longTable
        nx = 108
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "   Subshell   isBound   energy [a.u.]     energy " * TableStrings.inUnits("energy") * 
                        "   st-grid   r_max [a.u.]    <r> [a.u.]    <r^2> [a.u.]  ")
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  sh  in Defaults.GBL_STANDARD_SUBSHELL_LIST
            if haskey(orbitals, sh)
                v = orbitals[sh]
                sen_au = @sprintf("%.8e", v.energy)
                sen    = @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", v.energy))
                if !v.useStandardGrid   error("Non-standard grids not yet implemented.")    end
                p2 = zeros( size(v.P, 1) );    p2 = p2 .+ (v.P .* v.P);    p2 = p2 .+ (v.Q .* v.Q)
                imax    = findmax(p2)
                rmax    = grid.r[imax[2]]
                rmean   = RadialIntegrals.rkDiagonal(1, v, v, grid)
                r2mean  = RadialIntegrals.rkDiagonal(2, v, v, grid)
                srmax   = @sprintf("%.5e", rmax)
                srmean  = @sprintf("%.5e", rmean)
                sr2mean = @sprintf("%.5e", r2mean)
                sa      = "     $sh    $(v.isBound)   " * sen_au * "   " * sen * "   $(v.useStandardGrid)     " 
                sa      = sa * srmax * "    " * srmean * "    " * sr2mean
                println(sa)
            end
        end 
        println("  ", TableStrings.hLine(nx))
    else
        nx = 62
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "   Subshell   isBound   energy [a.u.]     energy " * TableStrings.inUnits("energy") * "   st-grid ")
        println(stream, "  ", TableStrings.hLine(nx))
        #
        for  sh  in Defaults.GBL_STANDARD_SUBSHELL_LIST
            if haskey(orbitals, sh)
                v = orbitals[sh]
                sen_au = @sprintf("%.8e", v.energy)
                sen    = @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", v.energy))
                sa     = "     $sh    $(v.isBound)   " * sen_au * "   " * sen * "   $(v.useStandardGrid)"
                println(sa)
            end
        end 
        println("  ", TableStrings.hLine(nx))
    end
    
    return( nothing )
end



"""
`Basics.display(stream::IO, channels::Array{AtomicState.GreenChannel,1})`  
    ... displays the (generated Green function)  channels in a neat and compact form; nothing is returned in this case.
"""
function  Basics.display(stream::IO, channels::Array{AtomicState.GreenChannel,1})
    println(stream, " ")
    println(stream, "  Green function channels:")
    println(stream, " ")
    for  channel in channels
        energies = Float64[]
        for  level in channel.gMultiplet.levels   push!(energies, level.energy)  end
        minenergy = minimum(energies);    maxenergy = maximum(energies)
        println(stream, "  Channel with $(channel.symmetry) symmetry, $(length(channel.gMultiplet.levels)) levels and energies [Hartree]:" *
                        " $minenergy ... $maxenergy")
    end
    
    return( nothing )
end


"""
`Basics.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10, E0::Float64=0.)`  
    ... display on stream the level energies of the N lowest levels in ascending order, and together with the configuration 
        of most leading term in the level expansion. All level energies refer to the lowest + E0 [Hartree], which enables 
        one to refer to any level as the (total) 0. energy.  A neat table is printed but nothing is returned otherwise.
"""
function Basics.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10, E0::Float64=0.)
    allLevels = Level[];    nx = 104
    for  multiplet  in multiplets
        for  level in  multiplet.levels   push!(allLevels, level)      end
    end
    println(stream, ">>> Total number of levels is $(length(allLevels))  in all multiplets together. \n´")
    sortedLevels = Base.sort( allLevels, lt=Base.isless)
    energy0      = sortedLevels[1].energy + E0
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    J Parity   No. electrons    Energy " * TableStrings.inUnits("energy") * "        Leading configuration") 
    println(stream, "  ", TableStrings.hLine(nx))
    for  (n, level)  in  enumerate(sortedLevels)
        if  n > N   break   end
        mc2  = level.mc .* level.mc;   index = findmax(mc2)[2]
        conf = Basics.extractConfiguration(Basics.FromBasis(), level.basis, level.basis.csfs[index])
        #
        sa   = TableStrings.flushright(12, string(LevelSymmetry(level.J, level.parity)))  * 
                TableStrings.flushright(12, string(conf.NoElectrons))                      * 
                TableStrings.flushright(20, @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy - energy0))) * "       " *
                TableStrings.flushleft(64, string(conf)) 
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    
    return( nothing )
end


"""
`   + (stream::IO, channels::Array{AtomicState.GreenChannel,1}; N::Int64=10)`  
    ... display on stream the lowest N level energies of all levels from the (Green) channles in ascending order, and together with the 
        configuration of most leading term in the level expansion. A neat table is printed but nothing is returned otherwise.
"""
function Basics.displayLevels(stream::IO, channels::Array{AtomicState.GreenChannel,1}; N::Int64=10)
    multiplets = Multiplet[]
    for  channel in channels    push!(multiplets, channel.gMultiplet)   end
    Basics.displayLevels(stream, multiplets)
    
    return( nothing )
end


"""
`Basics.displayMeanEnergies(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10, E0::Float64=0.)`  
    ... display on stream the mean energy as well as the minimum and maximum energy of all levels from the 
        same (leading) configuration. All these energies are taken with regard to the lowest level energy
        or are shifted by lowest + E0 [Hartree], which enables one to shift the zero energy to any position.
        These mean, min and max energies are listed together with the corresponding configuration.
        A neat table is printed but nothing is returned otherwise.
"""
function Basics.displayMeanEnergies(stream::IO, multiplets::Array{Multiplet,1}; N::Int64=10, E0::Float64=0.)
    # Combine all levels into a single list
    allLevels = Level[]
    for  multiplet  in multiplets
        for  level in  multiplet.levels   push!(allLevels, level)      end
    end
    println(stream, "\n>> Total number of levels is $(length(allLevels))  in all multiplets together.")
    sortedLevels = Base.sort( allLevels, lt=Base.isless)
    energy0      = sortedLevels[1].energy + E0
    # Now determine and compare the leading configuration for all levels
    nall = length(sortedLevels);    notConsideredYet = trues(nall)
    leadingConfigs = Configuration[]
    for  (n, level)  in  enumerate(sortedLevels)
        mc2  = level.mc .* level.mc;   index = findmax(mc2)[2]
        conf = Basics.extractConfiguration(Basics.FromBasis(), level.basis, level.basis.csfs[index])
        push!(leadingConfigs, conf)
    end
    #
    configs = Configuration[];   meanEnergies = Float64[];   minEnergies = Float64[];   maxEnergies = Float64[];   noLevels = Int64[]
    emin = emax = emean = 0.
    for  (na, aLevel)  in  enumerate(sortedLevels)
        if  notConsideredYet[na]     
            notConsideredYet[na] = false;    selectedConf = leadingConfigs[na];   
            emin = 1.0e10;   emax = -1.0e10;    ne = 1;   emean = aLevel.energy
            # Now compare with those levels, which were not yet considered
            for  (nb, bLevel)  in  enumerate(sortedLevels)
                if  notConsideredYet[nb]   &&   selectedConf == leadingConfigs[nb]
                    notConsideredYet[nb] = false;      emean = emean + bLevel.energy
                    emin = min(emin, bLevel.energy);   emax  = max(emax, bLevel.energy);    ne = ne + 1
                end
            end
            push!(minEnergies, emin);      push!(maxEnergies, emax);    push!(meanEnergies, emean / ne)
            push!(configs, selectedConf);  push!(noLevels, ne)
        end 
    end
    #
    nx = 150
    println(stream, ">> All energies are shown in:  " * TableStrings.inUnits("energy") * "\n") 
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "    No. e's  P      mean En        min En    ...    max En     No. Levels    Leading configuration") 
    println(stream, "  ", TableStrings.hLine(nx))
    for  (n, nconf)  in  enumerate(configs)
        if  n > N   break   end
        #
        sa = " " * 
        TableStrings.flushright( 8, string(nconf.NoElectrons)) * "   " *
        TableStrings.flushright( 2, string( Basics.extractFromConfiguration(Basics.GetParity(), nconf) ))                      * "   " *
        TableStrings.flushright(12, @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", meanEnergies[n] - energy0))) * "   " *
        TableStrings.flushright(12, @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", minEnergies[n]  - energy0))) * " ... " *
        TableStrings.flushright(12, @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", maxEnergies[n]  - energy0))) * "  " *
        TableStrings.flushright( 8, string( noLevels[n] )) * "      " *  
        TableStrings.flushleft(64, string(nconf)) 
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    
    return( nothing )
end


"""
`Basics.displayOrbitalOverlap(stream::IO, leftOrbitals::Dict{Subshell, Orbital}, rightOrbitals::Dict{Subshell, Orbital}, 
                                grid::Radial.Grid; kappas::Array{Int64,1} = [-1,1,-2])`  
    ... display on stream the overlap of the orbitals from the set of leftOrbs and rightOrbs, and for the given
        kappa-symmetries. A neat table is printed but nothing is returned otherwise.
"""
function Basics.displayOrbitalOverlap(stream::IO, leftOrbitals::Dict{Subshell, Orbital}, rightOrbitals::Dict{Subshell, Orbital}, 
                                        grid::Radial.Grid; kappas::Array{Int64,1} = [-1,1,-2])
    nx = 150
    println(stream, "\n  Orbital overlaps") 
    println(stream,   "  ================") 
    #
    for kappa in kappas
        sa = "\n  "
        for  (kl,vl) in leftOrbitals
            if   kl.kappa != kappa          continue    end
            for  (kr,vr) in rightOrbitals
                if   kr.kappa != kappa      continue    end
                overlap = RadialIntegrals.overlap(vl, vr, grid)
                if     length(sa) < nx - 20   sa =    sa * "    |<$kl|$kr>|^2 = " * @sprintf("%.4e", overlap*overlap)
                else   println(stream, sa);   sa =  "  " * "    |<$kl|$kr>|^2 = " * @sprintf("%.4e", overlap*overlap) 
                end
            end
        end
        println(stream, sa)
    end
    ## println(stream, "  ", TableStrings.hLine(nx))
    
    return( nothing )
end


"""
`Basics.displayOrbitalProperties(stream::IO, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64, radiusWR::Float64, 
                                    nm::Nuclear.Model, grid::Radial.Grid)`  
    ... displays the orbital properties:
    
            subshell, energy, Delta energy (wrt Dirac), f(energy; mu, T), occupation No, <r>, <1/r>
            
        for all orbitals with occNo > 0.1 in descending order.
        A neat table is printed and nothing is returned otherwise.
"""
function Basics.displayOrbitalProperties(stream::IO, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64, radiusWR::Float64, 
                                            nm::Nuclear.Model, grid::Radial.Grid)
    subshells = Subshell[];     ens  = Float64[];     rExps    = Float64[];     rInvs = Float64[];     
    occNos    = Float64[];      fFDs = Float64[];     deltaEns = Float64[];     totalOcc = 0.
    # Prepare lists of the individual properties for later sorting and printing
    for  (k,v) in orbitals
        fFD      = Basics.FermiDirac(v.energy, chemMu, temp)
        occNo    = (Basics.twice(Basics.subshell_j(k)) + 1) * fFD
        rExp     = RadialIntegrals.rkDiagonal( 1, v, v, grid) 
        rInv     = RadialIntegrals.rkDiagonal(-1, v, v, grid)
        enDirac  = HydrogenicIon.energy(k, nm.Z)
        totalOcc = totalOcc + occNo
        push!(subshells,  k);    push!(ens, v.energy);   push!(rExps, rExp);    push!(rInvs, rInv)     
        push!(occNos, occNo);    push!(fFDs, fFD);       push!(deltaEns, v.energy - enDirac)
    end
    ps = sortperm(ens, lt = isless)
    
    nx = 99
    println(stream, " ")
    println(stream, "  Orbital properties ordered by occupation numbers for chemical potential" * 
                    "  mu = $chemMu [Hartree], \n  temperature = $temp [Hartree] and R^WS = $radiusWR [a_o]")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(12, "Subshell"   ; na=0);                sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(14, "Energy   "  ; na=0);                sb = sb * TableStrings.center(14, "[Hartree] " ; na=0)
    sa = sa * TableStrings.center(14, "Delta en"   ; na=0);                sb = sb * TableStrings.center(14, "[Hartree]"  ; na=0)
    sa = sa * TableStrings.center(16, "occupation No"    ; na=0);          sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(16, "f^FD (en, mu, T)" ; na=0);          sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(12, "< r >"      ; na=2);                sb = sb * TableStrings.center(12, "  [a_o]"    ; na=2)
    sa = sa * TableStrings.center(12, "< 1/r >"    ; na=0);                sb = sb * TableStrings.center(12, "  [1/a_o]"  ; na=0)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    for  (ip,p)  in  enumerate(ps)
        sa = "    " * string(subshells[p]) * "   ";                   sa = sa[1:12]
        sa = sa * "    " * @sprintf("%.4e", ens[p])      * "    ";    sa = sa[1:28]
        sa = sa * "    " * @sprintf("%.4e", deltaEns[p]) * "    ";    sa = sa[1:44]
        sa = sa * "  "   * @sprintf("%.4e", occNos[p])   * "    "
        sa = sa * "  "   * @sprintf("%.4e", fFDs[p])     * "   "
        sa = sa * "  "   * @sprintf("%.4e", rExps[p])    * " "
        sa = sa * "  "   * @sprintf("%.4e", rInvs[p])    * "    "
        println(stream, sa)
    end
    println(stream, "    total                                     " * @sprintf("%.4e", totalOcc))
    println(stream, "  ", TableStrings.hLine(nx))
    
    return ( nothing )
end

    
"""
`Basics.expandOrbital(orbital::Orbital, boundOrbitals::Dict{Subshell, Orbital}, threshold::Float64, grid::Radial.Grid; 
                        printout::Bool=false)`  
    ... expand the orbital in the given orbitalList; only those orbitals with  <orb | orbital>^2  >= threshold are
        taken into account. If printout=true, this function 'reports' about the expansion by showing the 
        orbitals, their overlap as well as the total contribution. In general, the given orbital is not normalized
        and the expansion is usually incomplete. A list of (non-normalized) orbitals is returned in which
        the large and small components are mutiplied by the corresponding overlap.. 
"""
function Basics.expandOrbital(orbital::Orbital, boundOrbitals::Dict{Subshell, Orbital}, threshold::Float64, grid::Radial.Grid; 
                                printout::Bool=false)
    expansion = Orbital[];    n = orbital.subshell.n;    kappa = orbital.subshell.kappa;    tovl2 = 0.
    overlap   = RadialIntegrals.overlap(orbital, orbital, grid)   # | <orbital|overlap> |^2 of given orbital
    for  (subsh, orb) in boundOrbitals
        if  subsh.kappa == kappa  
            ovl = RadialIntegrals.overlap(orb, orbital, grid)
            if  ovl^2 >= threshold
                tovl2  = tovl2 + ovl^2
                newOrb = Orbital(orb.subshell, orb.isBound, orb.useStandardGrid, orb.energy, ovl * orb.P, ovl * orb.Q, 
                                    ovl * orb.Pprime, ovl * orb.Qprime, orb.grid)
                push!(expansion, newOrb)
                ## println(">>>> $subsh [$(newOrb.subshell)] with overlap $ovl appended to expansion of $(orbital.subshell).")   
            end
        end
    end
    println(">>> Expansion of $(orbital.subshell) with |<..|..>|^2 = $overlap into $(length(expansion)) partial waves with " *
            "total overlap |<new|given>|^2 = $tovl2 ")
    
    return( expansion )
end

    
"""
`Basics.extractMeanEnergy(pqn::Int64, basis::Basis)`  
    ... extracts the mean orbital energy of all orbitals with principal quantum number pqn; an en::Float64 [Hartree] is returned.
        This value has a physical meaning only for sufficiently large pqn (Rydberg orbitals).
"""
function Basics.extractMeanEnergy(pqn::Int64, basis::Basis)
    en = 0.;  nn = 0
    for  (k,v)  in basis.orbitals
        if  k.n == pqn   en = en + v.energy;   nn = nn + 1    end 
    end 
    return( en/nn )
end

    
"""
`Basics.extractMeanOccupation(basis::Basis)`  
    ... extracts the mean occupation numbers for all subshells of the given basis; an 
        meanOccs::Dict{Subshell, Float64} is returned.
"""
function Basics.extractMeanOccupation(basis::Basis)
    meanOccs = Dict{Subshell, Float64}();   ncsf = length(basis.csfs)
    for  (is, subsh)  in  enumerate(basis.subshells)
        occ = 0.;   for  csf in basis.csfs   occ = occ + csf.occupation[is]   end   
        meanOccs[subsh] = occ / ncsf
    end
    return( meanOccs )
end


"""
`Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1})  
    ... extract a nonrelativistic shell list from the given subshell list; the shells appear in the same order
        as the corresponding subshells. A shellList::Array{Subshell,1} is returned.
"""
function Basics.extractNonrelativisticShellList(subshells::Array{Subshell,1}) 
    shellList = Shell[]

    for  subsh  in subshells
        l = Basics.subshell_l(subsh);   shell = Shell( subsh.n, l)
        if  shell in shellList      else  push!(shellList, shell)   end
    end 
    
    return( shellList )
end


"""
`Basics.extractNonrelativisticShellList(confs::Array{Configuration,1})  
    ... extract a nonrelativistic shell list from the given configuration list; the shells are not ordered. 
        A shellList::Array{Subshell,1} is returned.
"""
function Basics.extractNonrelativisticShellList(confs::Array{Configuration,1}) 
    shellList = Shell[]
    
    for conf in confs
        for (k,v) in conf.shells
            if  k in shellList      else    push!(shellList, k)     end
        end
    end
    
    return( shellList )
end


"""
`Basics.extractNonrelativisticShellList(multiplet::Multiplet)  
    ... extract a nonrelativistic shell list from the given multiplet; however, the shells are not ordered. 
        A shellList::Array{Shell,1} is returned.
"""
function Basics.extractNonrelativisticShellList(multiplet::Multiplet) 
    shellList = Shell[]
    confs     = Basics.extractConfigurations(Basics.FromBasis(), multiplet.levels[1].basis)
    
    for conf in confs
        for (k,v) in conf.shells
            if  k in shellList      else    push!(shellList, k)     end
        end
    end
    
    return( shellList )
end


"""
`Basics.extractNonrelativisticShellList(multiplets::Array{Multiplet,1})  
    ... extract a nonrelativistic shell list from the given list of multiplets; however, the shells are not ordered. 
        A shellList::Array{Shell,1} is returned.
"""
function Basics.extractNonrelativisticShellList(multiplets::Array{Multiplet,1}) 
    shellList = Shell[]
    for  mp in multiplets
        shellList = Basics.merge(shellList, Basics.extractNonrelativisticShellList(mp))
    end
    
    return( shellList )
end


"""
`Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)`  
    ... extracts the quantum numbers of the nonrelativistic open shells from the given csfr. Here, we assume that csfR is 
        defined in the given basis, and it is checked that the subshells of each nonrelativistic shell occur subsequently
        in basis and in the right sequence j = l - 1/2, j = l + 1/2, respectively. 
        An rQN::Dict{Shell, Tuple{Tuple{Int64,Int64,Int64,Int64,Int64}, Tuple{Int64,Int64,Int64,Int64,Int64}}}
                = Dict(shell => (rQNm, rQNp) ) is returned with rQNm = (idxm, Nm, Qm, Jm, Jmx) and 
        rQNp = (idxp, Np, Qp, Jp, Jpx). The two tuples rQNm, rQNp contains all the subshell quantum numbers and their
        coupling in csfR.
"""
function Basics.extractOpenShellQNfromCsfR(csfR::CsfR, basis::Basis)
    wa = Dict{Shell,Tuple{NTuple{5,Int64},NTuple{5,Int64}}}()
    for  s = 1:length(basis.subshells)
        sx = s + 1;     subsh = basis.subshells[s];     n = subsh.n;        occ   = csfR.occupation[s]
        l  = Basics.subshell_l(subsh);    j2 = Basics.subshell_2j(subsh);   shell = Shell(n,l)
        # Get l,j from next subshell ... if still defined
        if  sx <= length(basis.subshells)       subshx = basis.subshells[sx]
            nx = subshx.n;    lx = Basics.subshell_l(subshx);   j2x = Basics.subshell_2j(subshx);   occx   = csfR.occupation[sx]
        end
        # Now collect the quantum numbers into tuples
        if  l == 0   &&   0 < occ < 2
            ## wa = Base.merge( wa, Dict( shell => ( (s, csfR.occupation[s], csfR.seniorityNr[s], 
            ##                                        Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
            ##                                        (-9, 0, -9, 0, Basics.twice(csfR.subshellX[s]) ) ) ) )
            wa = Base.merge( wa, Dict( shell => ( (-9, 0, -9, 0, 0 ),
                                                (s, csfR.occupation[s], csfR.seniorityNr[s], 
                                                Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ) ) ) )
        elseif  2l + 1 == j2    continue
        elseif  2l - 1 == j2    &&    (n != nx  ||  l != lx  ||  j2x != j2 + 2)
                error("stop a")
        elseif  0 < occ + occx < 2* (2l + 1)
            wa = Base.merge( wa, Dict(shell => ( (s, csfR.occupation[s], csfR.seniorityNr[s], 
                                                Basics.twice(csfR.subshellJ[s]), Basics.twice(csfR.subshellX[s]) ), 
                                                (sx, csfR.occupation[sx], csfR.seniorityNr[sx], 
                                                Basics.twice(csfR.subshellJ[sx]), Basics.twice(csfR.subshellX[sx]) ) ) ) )
        end
    end
    
    return( wa )
end


"""
`Basics.extractRelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)`  
    ... extract the nonrelativistic configuration from the given CSF that is defined in basis.csfs. 
        A conf::Configuration is returned.
"""
function  Basics.extractRelativisticConfigurationFromCsfR(csf::CsfR,  basis::Basis)
    subshells = basis.subshells
    
    newSubshells = Dict{Subshell,Int64}();    NoElectrons = 0
    for  (i, subsh) in enumerate(subshells)
        occ = csf.occupation[i]
        if   occ > 0  newSubshells = Base.merge( newSubshells, Dict( subsh => occ));     NoElectrons = NoElectrons + occ   end
    end
    
    if      basis.NoElectrons != NoElectrons    error("stop a")
    else    conf = ConfigurationR(newSubshells, NoElectrons)
    end
    
    return( conf )
end


"""
`Basics.extractRelativisticSubshellList(confList::Array{Configuration,1})`  
    ... extract all (relativistic) subshells that (will) contribute to the given configuration list, and for which (start) 
        orbitals will be needed in course of the computation. A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractRelativisticSubshellList(confList::Array{Configuration,1})
    # First extract all shells, then unique this list and, finally, convert in corresponding subshells
    shellList = Shell[]
    
    for  conf in confList
        for  (k,v) in conf.shells   push!( shellList, k)    end
    end
    
    shellList = unique(shellList)
    subshellList = Subshell[]
    
    for  sh in shellList    
        if  sh.l == 0    push!( subshellList, Subshell(sh.n, -1))
        else             push!( subshellList, Subshell(sh.n, sh.l));    push!( subshellList, Subshell(sh.n, -sh.l -1))
        end
    end
    subshellList = Base.sort(subshellList)
    println(">> Generated subshell list:  $subshellList ")
    return( subshellList )
end


"""
`Basics.extractRelativisticSubshellList(csf::CsfR, subshells::Array{Subshell,1})`  
    ... extract all (relativistic) subshells that are occupied in the given CSF with regard to the given subshells. 
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractRelativisticSubshellList(csf::CsfR, subshells::Array{Subshell,1})
    reducedSubshells = Subshell[]
    
    for (i, subsh) in enumerate(subshells)
        if  csf.occupation[i]  != 0     push!(reducedSubshells, subsh)  end
    end
    
    return( reducedSubshells )
end


"""
`Basics.extractRelativisticSubshellList(level::Level)`  
    ... extract all (relativistic) subshells that are occupied in the given level in CSF with weigth > 0.1. 
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractRelativisticSubshellList(level::Level)
    # Only consider CSF in the basis with non-zero mixing coefficients
    reducedSubshells = Subshell[]

    for (i, coeff) in enumerate(level.mc)
        if  coeff*coeff  >  0.1   
            redSubshells     = Basics.extractRelativisticSubshellList(level.basis.csfs[i], level.basis.subshells)
            reducedSubshells = append!(reducedSubshells, redSubshells)
        end
    end
    reducedSubshells     = unique(reducedSubshells)
    
    return( reducedSubshells )
end


"""
`Basics.extractRelativisticSubshellList(rep::Representation)`  
    ... extract all (relativistic) subshells that (will) contribute to the given configuration list, and for which (start) 
        orbitals will be needed in course of the computation. A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractRelativisticSubshellList(rep::Representation)
    # First extract all shells, then unique this list and, finally, convert in corresponding subshells
    shellList = Shell[]
    
    for  conf in rep.refConfigs
        for  (k,v) in conf.shells   push!( shellList, k)    end
    end
    
    for step in rep.repType.steps
        append!( shellList, step.seFrom);    append!( shellList, step.seTo)    
        append!( shellList, step.deFrom);    append!( shellList, step.deTo)   
        append!( shellList, step.teFrom);    append!( shellList, step.teTo)    
        append!( shellList, step.qeFrom);    append!( shellList, step.qeTo)   
    end
    
    shellList = unique(shellList)
    subshellList = Subshell[]
    
    for  sh in shellList    
        if  sh.l == 0    push!( subshellList, Subshell(sh.n, -1))
        else             push!( subshellList, Subshell(sh.n, sh.l));    push!( subshellList, Subshell(sh.n, -sh.l -1))
        end
    end
    
    return( subshellList )
end


"""
`Basics.extractRydbergSubshellList(level::Level, n0::Int64, wx::Float64)`  
    ... extract all (relativistic) Rydberg subshells with principal quantum number n > n0, which are occupied 
        in the given level in CSF with weigth > wx. A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractRydbergSubshellList(level::Level, n0::Int64, wx::Float64)
    # Only consider CSF in the basis with non-zero mixing coefficients
    reducedSubshells = Subshell[]

    for (i, coeff) in enumerate(level.mc)
        if  coeff*coeff  >  wx
            redSubshells     = Basics.extractRelativisticSubshellList(level.basis.csfs[i], level.basis.subshells)
            reducedSubshells = append!(reducedSubshells, redSubshells)
        end
    end
    reducedSubshells = unique(reducedSubshells)
    rydbergSubshells = Subshell[]
    for  subsh  in  reducedSubshells   if  subsh.n > n0   push!(rydbergSubshells, subsh)   end    end
    
    return( rydbergSubshells )
end


"""
`Basics.extractShellList(confList::Array{Configuration,1})`  
    ... extract all (nonrelativistic) shells that (will) contribute to the given configuration list. 
        A shellList::Array{Shell,1} is returned.
"""
function Basics.extractShellList(confList::Array{Configuration,1})
    # First extract all shells, then unique this list.
    shellList = Shell[]
    
    for  conf in confList
        for  (k,v) in conf.shells   push!( shellList, k)    end
    end
    
    shellList = unique(shellList)
    shellList = sort(shellList)
    ## println(">> Generated shell list:  $shellList ")
    
    return( shellList )
end


"""
`Basics.extractSubshellList(confList::Array{ConfigurationR,1})`  
    ... extract all (relativistic) subshells that (will) contribute to the given configuration list. 
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractSubshellList(confList::Array{ConfigurationR,1})
    # First extract all shells, then unique this list.
    subshList = Subshell[]
    
    for  conf in confList
        for  (k,v) in conf.subshells   push!( subshList, k)    end
    end
    
    subshList = unique(subshList)
    subshList = sort(subshList)
    println(">> Generated subshell list:  $subshList ")
    
    return( subshList )
end


"""
`Basics.extractSubshellList(conf::Configuration, orbitals::Dict{Subshell, Orbital})`  
    ... extract all (relativistic) subshells that do NPT have an orbital in orbitals. 
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.extractSubshellList(conf::Configuration, orbitals::Dict{Subshell, Orbital})
    shells    = Basics.extractShellList([conf]);    allSubshells = Subshell[];   notSubshells = Subshell[]
    for  shell in shells    append!(allSubshells, Basics.shellSplitIntoSubshells(shell))   end
    
    for subsh in allSubshells
        if haskey(orbitals, subsh)    else   push!(notSubshells, subsh)   end
    end
    
    
    return( notSubshells )
end


"""
`Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)`  
    ... extracts the nonrelativistic shell list as well as the (nonrelativistic) occupation of the given CsfR. 
        Here, we assume that csfR is defined in the given basis and that all subshells occur in standard sequence
        for a jjJ-LSJ transformation. A tuple tp::Tuple{Array{Shell,1},Array{Int64,1}}  of two arrays of equal 
        length is returned that contains the nonrelativistic shells and their occupation in the specification of 
        csfR. 
"""
function Basics.extractShellOccupationFromCsfR(csfR::CsfR, basis::Basis)
    shellList = Basics.extractNonrelativisticShellList(basis.subshells) 
    occList   = Int64[]
    
    for  s = 1:length(basis.subshells)
        l = Basics.subshell_l(basis.subshells[s]);   j2 = Basics.subshell_2j(basis.subshells[s])
        if      l == 0          push!( occList, csfR.occupation[s])
        elseif  2*l + 1 == j2   continue    # contributions from j = l + 1/2 are collected together with j = l - 1/2 below
        else    push!( occList, csfR.occupation[s] + csfR.occupation[s+1])
        end
    end 
    
    return( (shellList, occList) )
end


"""
`Basics.extractValenceShell(basis::Basis)`  
    ... extract the valence shell from basis, ie. the Shell with the (in its magnitude) lowest binding energy.
        A shell::Shell is returned.
"""
function Basics.extractValenceShell(basis::Basis)
    # Find the subshell with the lowest binding energy
    subShell = Subshell(100,-1);   energy = 1.0e6
    for (k,v) in basis.orbitals
        if  abs(v.energy) < energy    energy = abs(v.energy);    subShell = k   end 
    end
    shell = Shell(subShell.n, Basics.subshell_l(subShell) )
    
    return( shell )
end


"""
`Basics.FermiDirac(epsilon::Float64, mu::Float64, temp::Float64)`  
    ... computes the Fermi-Dirac function f(epsilon, mu; temp) if all parameters are given in atomic
        (Hartree) units.
"""
function Basics.FermiDirac(epsilon::Float64, mu::Float64, temp::Float64)
    wa = (epsilon - mu) / temp
    wa = 1.0 / (exp(wa) + 1.)
    return( wa )
end
