
export  isless

"""
`Basics.hasSubshell(pqn::Int64, subshells::Array{Subshell,1})`  
    ... returns true if one (or more) of subshells has the principal quantum number pqn, and false otherwise.
        A value::Bool is returned
"""
function Basics.hasSubshell(pqn::Int64, subshells::Array{Subshell,1})
    for  subsh in subshells   if  subsh.n == pqn   return( true)  end   end
    return( false )    
end


"""
`Basics.interpolateOnGridGrasp92((F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid) )`  
    ... to interpolate the (radial) function F from the oldgrid to newgrid, by a call to the Grasp92 Fortran procedures; 
        an function G::Array{Float64,1} is returned. **Not yet implemented !**    
"""
function Basics.interpolateOnGridGrasp92(from::Tuple{Array{Float64,1},Radial.Grid}, to::Tuple{DataType,Radial.Grid})
    # First prepare all fields and settings to call quad_grasp92()
    error("Not yet implemented !")
end


"""
`Basics.interpolateOnGridTrapezRule((F::Array{Float64,1}, oldGrid::Radial.Grid), (G::Void, newGrid::Radial.Grid)` 
    ... to integrate the (radial) function F from the oldgrid to newgrid, by using a trapez rule from the Numerical Recipies. 
        A function G::Array{Float64,1} is returned. **Not yet implemented !**.
"""
function Basics.interpolateOnGridTrapezRule(from::Tuple{Array{Float64,1},Radial.Grid}, to::Tuple{DataType,Radial.Grid})
    error("Not yet implemented !")
end


"""
`Basics.integrate(::NewtonCotes, F::Array{Float64,1}, grid::Radial.Grid)`
    ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r), by using
        a 5-point Newton-Cotes integration formula; a value::Float64 is returned.
"""
function Basics.integrate(::NewtonCotes, F::Array{Float64,1}, grid::Radial.Grid)
    return integrateOnGridNewtonCotes(F, grid)
end


"""
` + Basics.integrate(::SimpsonRule, F::Array{Float64,1}, grid::Radial.Grid)`
    ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r), by using
        Simpson's rule; a value::Float64 is returned.
"""
function Basics.integrate(::SimpsonRule, F::Array{Float64,1}, grid::Radial.Grid)
    return integrateOnGridSimpsonRule(F, grid)
end


"""
` + Basics.integrate(::TrapezRule, F::Array{Float64,1}, grid::Radial.Grid)`
    ... to integrate the (radial) function F over the given (radial) grid, int_0^infinity dr F(r), by using
        a simple trapezoid rule; a value::Float64 is returned.
"""
function Basics.integrate(::TrapezRule, F::Array{Float64,1}, grid::Radial.Grid)
    return integrateOnGridTrapezRule(F, grid)
end


"""
`Basics.integrateOnGridNewtonCotes()`  ... integrates by using a 5-point Newton-Cotes integration formula; see Basic.integrate().
"""
function Basics.integrateOnGridNewtonCotes(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45 * grid.h
    n = size(coefficients, 1)

    result = 0.;    i0     = 1

    while  abs(F[i0]) == 0
        i0 += 1
    end
    gamma  = log(F[i0]/F[i0+1] * grid.rp[i0+1]/grid.rp[i0]) / log(grid.r[i0]/grid.r[i0+1])
    result = 1/(gamma + 1) * grid.r[i0] * F[i0] / grid.rp[i0]

    result += F[i0] * coefficients[n]
    for i = i0+1:size(F, 1)
        result += coefficients[(i - i0) % (n - 1) + 1] * F[i]
    end

    if (size(F, 1) - i0) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end

    return( result )
end


"""
`Basics.integrateOnGridSimpsonRule()`  ... integrates by using Simpson's rule; see Basics.integrate().
"""
function Basics.integrateOnGridSimpsonRule(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 1, 4, 1]) / 3 * grid.h
    n = size(coefficients, 1)

    result = 0.
    result += F[1] * coefficients[n]
    for i = 2:size(F, 1)
        result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
    end

    if (size(F, 1) - 1) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end

    return( result )
end


"""
`Basics.integrateOnGridTrapezRule()`  ... integrates by using a simple trapez rule; see Basics.integrate().
"""
function Basics.integrateOnGridTrapezRule(F::Array{Float64,1}, grid::Radial.Grid)
    coefficients = Array{Float64}([2 * 1, 1]) / 2 * grid.h
    n = size(coefficients, 1)

    result = 0.

    result += F[1] * coefficients[n]
    for i = 2:size(F, 1)
        result += coefficients[(i - 1) % (n - 1) + 1] * F[i]
    end

    if (size(F, 1) - 1) % (n - 1) + 1 == 1
        result -= coefficients[n] * F[size(F, 1)]
    end

    return( result )
end


function Base.isless(x::Basics.Shell, y::Basics.Shell)
    if  x.n == y.n      return( x.l < y.l )
    else                return( x.n < y.n )
    end
end


function Base.isless(x::Basics.Subshell, y::Basics.Subshell)
    if      x.n == y.n   &&  abs(x.kappa) == abs(y.kappa)     return( x.kappa < y.kappa )
    elseif  x.n == y.n                                        return( abs(x.kappa) < abs(y.kappa) )
    else                                                      return( x.n < y.n )
    end
end


function Base.isless(x::ManyElectron.Level, y::ManyElectron.Level)
    return x.energy < y.energy
end


function Base.isless(x::Basics.AngularJ64, y::Basics.AngularJ64)
    return x.num < y.num
end

    
"""
`Basics.isSymmetric(matrix::Array{Float64,2})`  
    ... returns true if the matrix is symmetric and false otherwise; if not symmetric, the function prints the 
        first 10 pairs  M_ij, M_ji  which violate the symmetry condition.
"""
function Basics.isSymmetric(matrix::Array{Float64,2})
    symmetric = true;    na = 0
    if  size(matrix,1) !=  size(matrix,2)  error("not quadratic")   end
    #
    for j = 1:size(matrix,2)
        for i = j+1:size(matrix,1)
            if abs( matrix[i,j] - matrix[j,i]) > 1.0e-10   
                na = na + 1;   symmetric = false;   println("*** $i  $j  $(matrix[i,j])  $(matrix[j,i]) ")
                if  na >= 10   return( symmetric )  end
            end
        end
    end
    
    return( symmetric )
end

    
"""
`Basics.isZero(matrix::Array{Float64,2})`  
    ... returns true if the matrix is zero and false otherwise; it sums up the |M[i,j]^2 and prints the results. 
        If not zero, the function prints the first 10 matrix elements  M_ij  which violate this condition.
"""
function Basics.isZero(matrix::Array{Float64,2})
    isz = true;    na = 0;   wa = 0.
    #
    for j = 1:size(matrix,2)
        for i = 1:size(matrix,1)
            wa = wa + abs( matrix[i,j])^2
            #
            if abs( matrix[i,j]) > 1.0e-10   
                na = na + 1;   isz = false;   println("*** $i  $j  $(matrix[i,j]) ")
                if  na >= 10   return( isz )    end
            end
        end
    end
    println(">> isZero():   sum^2 = $(wa)")
    
    return( isz )
end

    
"""
`Basics.isViolated()`  ... returns true if some rule/limitation is violated, and false otherwise.

+ `(conf::Configuration, restriction::AbstractConfigurationRestriction)`  
    ... returns true if restriction is violated by the given conf, and false otherwise.
"""
function Basics.isViolated(conf::Configuration, restriction::AbstractConfigurationRestriction)
    wa = false;  ne = 0
    if      typeof(restriction) == RestrictNoElectronsTo
        for (sh,v) in  conf.shells
            if  sh.n >= restriction.nmin    ||    sh.l >= restriction.lmin    ne = ne + v   end
        end
        if  ne > restriction.ne                                                                wa = true   end
    elseif  typeof(restriction) == RestrictMaximumDisplacements
        shells = Basics.extractShellList([conf, restriction.conf])
        dis = 0;   
        for shell in shells
            if      haskey(conf.shells, shell)  &&  haskey(restriction.conf.shells, shell)  
                dis = dis + abs( conf.shells[shell] - restriction.conf.shells[shell] )
            elseif  haskey(conf.shells, shell)              dis = dis + abs( conf.shells[shell])
            elseif  haskey(restriction.conf.shells, shell)  dis = dis + abs( restriction.conf.shells[shell])
            else    error("stop a")
            end
        end
        if  dis > restriction.maxDisplace                                                      wa = true   end
    elseif  typeof(restriction) == RestrictParity
        if  Basics.extractFromConfiguration(Basics.GetParity(), conf) != restriction.parity    wa = true   end
    elseif  typeof(restriction) == RestrictToShellDoubles
        for (sh,v) in  conf.shells
            if  sh.n >= restriction.nmin    
                if  sh.l >= restriction.lmin 
                    println("restrictions:  $sh  $(sh.n)  $v  $(restriction.nmin)")
                    if !(v in [0,2])                              wa = true   end
                end
            end 
        end
    elseif  typeof(restriction) == RequestMinimumOccupation
        for (sh,v) in  conf.shells
            if  sh in restriction.shells   ne = ne + v     end
        end
        if  ne < restriction.ne                                   wa = true   end
    elseif  typeof(restriction) == RequestMaximumOccupation
        for (sh,v) in  conf.shells
            if  sh in restriction.shells   ne = ne + v     end
        end
        if  ne > restriction.ne                                   wa = true   end
    else    error("stop a")
    end
    
    return( wa )
end

    
"""
`Basics.isStandardSubshellList(basis::Basis)`  
    ... returns true if the subshell list basis.subshells is in standard order, and false otherwise.
        It requests especially that both subshells of the same shell (n,l) occur in the sequence j = l-1/2, j = l+1/2
        and that l increases before n increases.
"""
function Basics.isStandardSubshellList(basis::Basis)
    function  decimal(n,l,jnum)
        return( 10000n + 100l + jnum)
    end
    subshells = basis.subshells
    na = subshells[1].n;    la = Basics.subshell_l(subshells[1]);    ja = Basics.subshell_j(subshells[1])
    
    for  i = 2:length(subshells)
        nb = subshells[i].n;    lb = Basics.subshell_l(subshells[i]);    jb = Basics.subshell_j(subshells[1])
        if  decimal(nb,lb,jb.num) <= decimal(na,la,ja.num)    return( false )   end
    end
    
    return( true )
end

    
"""
`Basics.lastPoint(wa::Array{Float64,1}, eps::Float64)`  
    ... returns the last point p of wa[1:end]  for which |wa[p]| >= eps; an index p::Int64 is returned.
"""
function Basics.lastPoint(wa::Array{Float64,1}, eps::Float64)
    if  eps < 0.    error("Improper eps = $eps")    end
    for  p = length(wa):-1:1
        if  abs(wa[p]) >= eps       return(p)       end
    end
    error("No array element abs(wa[p]) > eps = $eps")
end


"""
`Basics.LevelSymmetry(subsh::Subshell)`  ... constructor for a given (Subshell).
"""
function  Basics.LevelSymmetry(subsh::Subshell)
    if  rem(Basics.subshell_l(subsh), 2) == 0   sa = "+"    else    sa = "-"    end
    LevelSymmetry( Basics.subshell_j(subsh), Parity(sa) )    
end


"""
Basics.merge(bases::Array{Basis,1})`  
    ... to merge two (or more) atomic bases into a single basis::Basis. This method assumes the same number of electrons in all basis and 
        that the subshell lists are the same or can be made 'consistent' to each other. Two bases have a consistent subshell list if all 
        subshells, that appear in any of the two lists appear always in the same sequence (if they are not missing at all). In the merged 
        basis, the radial orbitals are taken from the basis (in the bases-array}, from where they are found first.
"""
function Basics.merge(bases::Array{Basis,1})

    if  length(bases) > 2 
        println("Number of bases = $(length(bases)) ...")   
        bs       = Basics.merge([ bases[1], bases[2] ])
        basesNew = [bs]
        for i = 3:length(bases)    push!( basesNew, bases[i])    end
        Basics.merge(basesNew)
    elseif  length(bases) == 1    return( basis[1] )
    elseif  length(bases) == 2
        # This is the essential step here; first create a consistent subshell list
        basisA = bases[1];    basisB = bases[2]
        wa = deepcopy( basisA.subshells );   wb = [ 100i  for i = 1:length(wa) ]
        # Add the subshells of basisB
        for  i = 1:length( basisB.subshells )
            if    basisB.subshells[i] in wa
            else  push!( wa, basisB.subshells[i] );    push!( wb, wb[i-1]+1 )   
            end
        end
        # Now sort wb and create a new subshell list from wa due to this sorting
        wc = sortperm( wb )
        newSubshells = Subshell[]
        for  i in wc   push!( newSubshells, wa[i] )    end
        println("Sorted subshell list: $(newSubshells) ")
        error("Not yet fully implemented, see source code ")
        # Transform all CSF from basisA and basisB to newSubshells ... and apply CsfRExcludeDouble()
        # return new basis
        #
    else    error("stop a")
    end
end


"""
`Basics.merge(multiplets::Array{Multiplet,1})`  
    ... to merge two (or more) atomic multiplets into a single multiplet::Multiplet. This method assumes (and checks) that all 
        levels have level.hasStateRep = true and that all levels refer to the same basis.
"""
function Basics.merge(multiplets::Array{Multiplet,1})

if      length(multiplets) == 1    return( multiplets[1] )
    elseif  length(multiplets) == 2
        levels = copy(multiplets[1].levels);   basis = multiplets[1].levels[1].basis
        for  lev  in  multiplets[2].levels
            if  !(lev.hasStateRep)   ||   lev.basis != basis    error("Levels of multiplets do not refer to the same basis.")   end
            push!(levels, lev)
        end
        mp = Multiplet( multiplets[1].name * "+" * multiplets[2].name, levels )
        return( mp )
    elseif  length(multiplets) > 2
        mp = Basics.merge([ multiplets[1], multiplets[2] ])
        for  k = 3:length(multiplets)
            mp = Basics.merge([ mp, multiplets[k] ])
        end
        return( mp )
    else    error("stop a")
    end
end


"""
`Basics.merge(aList::Array{Shell,1}, bList::Array{Shell,1}, ...)`  
    ... to merge two (or more) shell lists into a single list, to unify and to order them. 
        A cList::Array{Shell,1} is returned
"""
function Basics.merge(aList::Array{Shell,1}, bList::Array{Shell,1})
    cList = copy(aList);      append!(cList, bList);      cList = Base.unique(cList)
    cList = Basics.sort(cList)
    return( cList )
end
function Basics.merge(aList::Array{Shell,1}, bList::Array{Shell,1}, cList::Array{Shell,1})
    dList = Basics.merge(aList, bList);     dList = Basics.merge(dList, cList)
    dList = Basics.sort(dList)
    return( dList )
end


"""
`Basics.merge(aList::Array{Subshell,1}, bList::Array{Subshell,1}, ...)`  
    ... to merge two (or more) subshell lists into a single list, to unify and to order them. 
        A cList::Array{Subshell,1} is returned
"""
function Basics.merge(aList::Array{Subshell,1}, bList::Array{Subshell,1})
    cList = copy(aList);      append!(cList, bList);      cList = Base.unique(cList)
    cList = Basics.sort(cList)
    return( cList )
end
function Basics.merge(aList::Array{Subshell,1}, bList::Array{Subshell,1}, cList::Array{Subshell,1})
    dList = Basics.merge(aList, bList);     dList = Basics.merge(dList, cList)
    dList = Basics.sort(dList)
    return( dList )
end


"""
`Basics.modifyLevelEnergies(multiplet::Multiplet)`  
    ... to shift the total energies of the atomic levels in multiplet interactively; a (new) multiplet::Multiplet is 
        returned in which the total energies are just modified accordingly. **Not yet implemented !**
"""
function Basics.modifyLevelEnergies(multiplet::Multiplet)
    error("Not yet implemented !")

    return( nothing )
end


"""
`Basics.modifyLevelMixing(level::Level, enhancementFaktor::Float64)`  
    ... to "enhance" the mixing of all CSF with |c_ik|^2 < 1.0e-3 by the given enhancementFaktor. 
        The weights all CSF are renormalized eventually. A (normalized) newLevel::Level is returned.
"""
function Basics.modifyLevelMixing(level::Level, enhancementFaktor::Float64)
    mcx = Float64[];   
    for  mc in level.mc   if  abs(mc)^2 < 1.0e-3   push!(mcx, enhancementFaktor * mc)   else  push!(mcx, mc)   end    end
    
    # Renormalize the mixing coefficients
    wx  = 0.;          for  mc  in  mcx   wx = wx + abs(mc)^2         end
    mcy = Float64[];   for  mc  in  mcx   push!(mcy, mc / sqrt(wx))   end
    wy  = 0.;          for  mc  in  mcy   wy = wy + abs(mc)^2         end
    
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, 
                     level.hasStateRep, level.basis, mcy)
    return( newLevel )
end


# using Plots
# pyplot()

using RecipesBase

"""
`Basics.plot(::RadialPotentials, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)`
    ... to plot one or more radial potentials, where N::Int64 is the number of grid points to be considered.
        call:  using Plots; pyplot()    ... to access this method by plot(...)
"""
function Basics.plot(::RadialPotentials, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end


@recipe function f(::RadialPotentials, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N = 0)
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(potentials)
    for  pot in potentials
        wb = wc
        nx = min(length(pot.Zr), N);   wb[1:nx] = pot.Zr[1:nx]
        append!(wa, wb)
        push!(labels, pot.name)
    end
    x = grid.r[1:N];     y = reshape(wa, (N, np))
    label --> reshape(labels, (1,np))
    x, y
end


"""
`Basics.plot(::RadialOrbitalsLarge, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`
    ... to plot the large component of one or more radial orbitals, where N::Int64 is the number of grid points.
        call:  using Plots; pyplot()    ... to access this method by plot(...)

`Basics.plot(::RadialOrbitalsSmall, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`
    ... to plot the small component of one or more radial orbitals.

`Basics.plot(::RadialOrbitalsBoth, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`
    ... to plot both the large and small component of one or more radial orbitals.
"""
function Basics.plot(::RadialOrbitalsLarge, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end

function Basics.plot(::RadialOrbitalsSmall, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end

function Basics.plot(::RadialOrbitalsBoth, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end


@recipe function f(::RadialOrbitalsLarge, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(orbitals)
    for  orb in orbitals
        wb = wc
        nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx]
        append!(wa, wb)
        push!(labels, "$(orb.subshell):large")
    end
    x = grid.r[1:N];     y = reshape(wa, (N, np))
    label --> permutedims(labels)
    x, y
end

@recipe function f(::RadialOrbitalsSmall, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(orbitals)
    for  orb in orbitals
        wb = wc
        nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx]
        append!(wa, wb)
        push!(labels, "$(orb.subshell):small")
    end
    x = grid.r[1:N];     y = reshape(wa, (N, np))
    label --> permutedims(labels)
    x, y
end

@recipe function f(::RadialOrbitalsBoth, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(orbitals)
    for  orb in orbitals
        wb = wc;    nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx];    append!(wa, wb)
        wb = wc;    nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx];    append!(wa, wb)
        push!(labels, "$(orb.subshell):large");   push!(labels, "$(orb.subshell):small")
    end
    x = grid.r[1:N];     y = reshape(wa, (N, 2np))
    label --> permutedims(labels)
    x, y
end


"""
`Basics.plot(theme::AbstractPlotTheme, lines::Array{PhotoEmission.Line,1})`
    ... to plot transition rates or oscillator strengths as function of transition energies. **Not yet implemented.**
"""
function Basics.plot(theme::AbstractPlotTheme, lines::Array{PhotoEmission.Line,1})
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
end


"""
`Basics.plot(theme::AbstractPlotTheme, lines::Array{PhotoEmission.Line,1}, widths::Float64)`
    ... to plot transition rates with Gaussian or Lorentzian broadening. **Not yet implemented.**
"""
function Basics.plot(theme::AbstractPlotTheme, lines::Array{PhotoEmission.Line,1}, widths::Float64)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
end
