
using  ..Basics,  ..Defaults, ..Pulse


"""
`struct  Liouville.RamanLevel`  
    ... defines a struct to comprise the level information for a Liouville time evolution in the stimulated-Raman scheme.

    + leadingConfig    ::Configuration       ... leading configuration of the level.
    + leadingNotation  ::String              ... leading notation, need to be added manually.
    + level            ::Level               ... atomic level as calculated self-consistently.
"""
struct RamanLevel 
    leadingConfig      ::Configuration 
    leadingNotation    ::String  
    level              ::Level  
end 


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Liouville.RamanLevel()`  ... constructor for the default values of Liouville.RamanLevel set
"""
function RamanLevel()
    RamanLevel( Configuration(), "xx", Level() )
end


# `Base.show(io::IO, lioulevel::Liouville.RamanLevel)`  ... prepares a proper printout of the variable settings::Liouville.RamanLevel.
function Base.show(io::IO, data::Liouville.RamanLevel) 
    println(io, "leadingConfig:          $(lioulevel.leadingConfig)  ")
    println(io, "leadingNotation:        $(lioulevel.leadingNotation)  ")
    println(io, "level:                  $(lioulevel.level)  ")
end


"""
`Liouville.determineLightField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)`  
    ... determines the field amplitude A(t) at the given time t at the (central) position of the atom or
        atomic cloud. All pulses must be of computational type. An fieldAmplitude::Float64 is returned.
"""
function determineLightField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)
    function gaussianEnvelope(t::Float64, sigma::Float64)
        # Determine the value of a Gaussian distribution, centered at 0., at time t and sigma
        wa = t^2/ (2*sigma)^2
        return( exp(-wa) / (sigma * sqrt(2pi)) )
    end 
    
    fieldA = 0.
    for pulse in pulses
        if    typeof(pulse) == Pulse.GaussianSimplified
            if  abs(t - pulse.timeDelay)  < 5 * pulse.fwhm    
                fieldA = fieldA  +  pulse.A0 * gaussianEnvelope(t - pulse.timeDelay, pulse.fwhm/2.0)
            end
        else  error("stop a; pulse = $pulse ")
        end
    end
    
    @show fieldA
    
    return( fieldA )
end


"""
`Liouville.displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, densityM::Matrix{ComplexF64})`  
    ... Display the current density matrix together with the ; nothing is returned.
"""
function displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, densityM::Matrix{ComplexF64})

    println(stream, " ")
    println(stream, "  Selected Liouville Raman levels and current density matrix:")
    println(stream, " ")
    for  (idx, liouLevel)  in  enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  ";   sa = sa[end-10:end]
        sa = sa * string(liouLevel.leadingConfig)    * "   "
        sa = sa * string(liouLevel.leadingNotation)  * "                          "
        sa = sa[1:70]
        row = densityM[idx, :]
        for z in row   sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))    end
        println(sa)
    end 
    println(stream, " ")
    
    return( nothing )
end

    
"""
`Liouville.displayGenericHamiltonian(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, atomicHM::Matrix{ComplexF64},
                                     couplingHM::Array{Function, 2})`  
    ... Display the current density matrix together with the ; nothing is returned.
"""
function displayGenericHamiltonian(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, atomicHM::Matrix{ComplexF64},
                                   couplingHM::Array{Function, 2})

    # Add the two hamiltonian matrix first
    noLevels         = length(liouvilleLevels)
    totalHamiltonian = [ t -> couplingHM[i,j](t) + atomicHM[i,j]  for i in 1:noLevels, j in 1:noLevels ]
        
    t = 0.1
    totalH = [f(t) for f in totalHamiltonian]
        
    
    println(stream, " ")
    println(stream, "  Selected Liouville Raman levels and generic Hamiltonian matrix, evaluated for t=0.1:")
    println(stream, " ")
    for  (idx, liouLevel)  in  enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  ";   sa = sa[end-10:end]
        sa = sa * string(liouLevel.leadingConfig)    * "   "
        sa = sa * string(liouLevel.leadingNotation)  * "                          "
        sa = sa[1:70]
        row = totalH[idx, :]
        for z in row   sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))    end
        println(sa)
    end 
    println(stream, " ")
    
    return( nothing )
end


"""
`Liouville.displayLevels(stream, liouvilleLevels::Array{Liouville.RamanLevel,1})`  
    ... Displays the list of Liouville Raman levels in a near format; nothing is returned.
"""
function displayLevels(stream, liouvilleLevels::Array{Liouville.RamanLevel,1})

    println(stream, " ")
    println(stream, "  Selected Liouville Raman levels:")
    println(stream, " ")
    for  (idx, liouLevel)  in  enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  ";   sa = sa[end-10:end]
        sa = sa * string(liouLevel.leadingConfig)    * "   "
        sa = sa * string(liouLevel.leadingNotation)  * "   "
        println(sa)
    end 
    println(stream, " ")
    
    return( nothing )
end


"""
`Liouville.initializeLevels(scheme::Liouville.StimulatedRamanScheme, multiplet::Multiplet)`  
    ... initialize the levels::Vector{Liouville.RamanLevel} on which the Liouville time evolution is based on.
        A list levels::Vector{Liouville.RamanLevel}
"""
function initializeLevels(scheme::Liouville.StimulatedRamanScheme, multiplet::Multiplet)
    liouvilleLevels = Liouville.RamanLevel[];   noLevels = length(scheme.levelSelection.indices);   foundIndices = falses(noLevels)
    
    # Proper level notations must be provided via the scheme
    if  length(scheme.levelNotations) - 1 != noLevels  
        error("Expect $noLevels strings for level notations, got $(scheme.levelNotations) ")
    end 
        
    for  (idx, index)  in  enumerate(scheme.levelSelection.indices)
        for  level in multiplet.levels
            if  index == level.index   foundIndices[idx] = true  
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                liouvLevel  = Liouville.RamanLevel(leadingConf, scheme.levelNotations[idx], level);  
                push!(liouvilleLevels, liouvLevel)
            end
        end 
    end
    if  all(foundIndices) == false  error("Not all level indices are found in multiplet; foundIndices = $foundIndices ")    end
    
    # Add a autoionizing/decay level to model the loss from the system
    push!(liouvilleLevels, Liouville.RamanLevel(Configuration("[He]"), scheme.levelNotations[end], Level() ) )
    
    Liouville.displayLevels(stdout, liouvilleLevels)
    
    return( liouvilleLevels )
end


"""
`Liouville.initializeDensityMatrix(liouvilleLevels::Array{Liouville.RamanLevel,1})`  
    ... initialize the density matrix that applies before the light-field start to act upon the atom.
        A dm::Matrix{Float64} is returned
"""
function initializeDensityMatrix(liouvilleLevels::Array{Liouville.RamanLevel,1})
    noLevels = length(liouvilleLevels);     energies = Float64[]
    densityM = zeros(ComplexF64, noLevels,noLevels)
    
    # Determine level with lowest energy and set population to 1.
    for  lioulevel  in  liouvilleLevels   push!(energies, lioulevel.level.energy)   end
    lowestEn = minimum(energies);         idx = findfirst(==(lowestEn), energies)
    densityM[idx, idx] = 1.0
    
    return ( densityM )
end


"""
`Liouville.initializeAtomicHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, liouvilleLevels::Array{Liouville.RamanLevel,1})`  
    ... initialize the Hamiltonian matrix ... in such a form that it can later be readily "updated"
        for some given light field. A hamiltonian::... is returned. A hamiltonian::Matrix{ComplexF64} is returned
"""
function initializeAtomicHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, liouvilleLevels::Array{Liouville.RamanLevel,1})
    noLevels    = length(liouvilleLevels);     energies = Float64[]
    hamiltonian = zeros(ComplexF64, noLevels,noLevels)
    
    # Determine level with lowest energy and set diagonal matrix elements to the excitation energie
    for  lioulevel  in  liouvilleLevels   push!(energies, lioulevel.level.energy)   end
    lowestEn = minimum(energies);         
    for n = 1:noLevels   hamiltonian[n,n] = energies[n] - lowestEn   end
    
    return ( hamiltonian )
end


"""
`Liouville.initializeCouplingHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, liouvilleLevels::Array{Liouville.RamanLevel,1})`  
    ... initialize the Hamiltonian matrix ... in such a form that it can later be readily "updated"
        for some given light field. A hamiltonian::... is returned. A hamiltonian::Matrix{Function} is returned
"""
function initializeCouplingHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, liouvilleLevels::Array{Liouville.RamanLevel,1})
    println("==> initializeCouplingHamiltonianMatrix() ... not yet implemented properly.")
    noLevels    = length(liouvilleLevels)
    hamiltonian = Matrix{Function}(undef, noLevels, noLevels)
    for  i in 1:noLevels,  j in 1:noLevels   hamiltonian[i,j] = t -> 0.0   end        # Assign a function f(t) = 0. to all matrix elements 
    
    @warn("Interaction matrix elements must be set manually at present.")
    
    # Assign manually the individual matrix elements
    hamiltonian[1,2]  =  hamiltonian[2,1]  =  t -> sin(t)
       
    return ( hamiltonian )
end


"""
`Liouville.evaluateHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, t::Float64, pulses::Vector{Pulse.AbstractPulse},
                                     atomicHM::Matrix{ComplexF64}, couplingHM::Array{Function, 2})`  
    ... to evaluate the Hamiltonian matrix for the light-field at given time t. 
"""
function evaluateHamiltonianMatrix(scheme::Liouville.StimulatedRamanScheme, t::Float64, pulses::Vector{Pulse.AbstractPulse})
    println("==> Liouville.evaluateHamiltonianMatrix() ... not yet implemented properly.")

    # Add the two hamiltonian matrix first
    noLevels         = size(atomicHM,1)
    totalHamiltonian = [ t -> couplingHM[i,j](t) + atomicHM[i,j]  for i in 1:noLevels, j in 1:noLevels ]
        
    t = 0.1
    totalH = [f(t) for f in totalHamiltonian]

    
    return ( totalH )
end


"""
`Liouville.evolveDensityMatrix(scheme::Liouville.StimulatedRamanScheme, approach::Liouville.AbstractEvolutionApproach, 
                             densityM::Array{Float64,2}, hamiltonianM::Array{Float64,2}, lightField::Float64).jena
                             `  
    ... to evolve the density matrix for a certain time-interval ... for a given density matrix and some generic Hamiltonian matrix. 
"""
function evolveDensityMatrix(scheme::Liouville.StimulatedRamanScheme, approach::Liouville.AbstractEvolutionApproach, 
                             densityM::Array{Float64,2}, hamiltonianM::Array{Float64,2}, lightField::Float64)
    println("Liouville.evolveDensityMatrix() ... not yet implemented properly.")
    dm = zeros(5,5)
    
    return ( dm )
end


"""
`Liouville.analyzeDensityMatrix(scheme::Liouville.StimulatedRamanScheme, ...)`  
    ... to display the density matrix 
"""
function analyzeDensityMatrix(scheme::Liouville.StimulatedRamanScheme, comp::Liouville.Computation)
    println("Liouville.analyzeDensityMatrix() ... not yet implemented properly.")
    dm = zeros(5,5)
    
    return ( dm )
end



"""
`Liouville.perform(scheme::Liouville.StimulatedRamanScheme, computation::Liouville.Computation; output::Bool=true)`  
    ... to perform a Liouville time-evolution computation for a stimulated Raman scheme. For output=true, a dictionary 
        is returned from which the relevant results can be can easily accessed by proper keys.
"""
function  perform(scheme::Liouville.StimulatedRamanScheme, computation::Liouville.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    
    println("")
    printstyled("Liouville.perform(): The Liouville computation of stimulated Raman scattering starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------- \n", color=:light_green)
    
    # Convert the pulses into a computational useful form
    pulses = Pulse.AbstractPulse[]
    for  pulse in computation.pulses
        if      typeof(pulse) == Pulse.GaussianSimplified     push!(pulses, pulse)
        elseif  typeof(pulse) == Pulse.FelPulse               push!(pulses, Pulse.convertPulse(pulse) )
        else    error("Unknown pulse = $pulse ")
        end
    end

    # Compute the level structure for all reference configurations
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel, computation.grid, computation.asfSettings)
    #
    # Initialize the Liouville levels as well as the density and Hamiltonian matrix
    levels      = Liouville.initializeLevels(scheme, multiplet)
    densityM    = Liouville.initializeDensityMatrix(levels)
    atomicHM    = Liouville.initializeAtomicHamiltonianMatrix(scheme, levels)
    couplingHM  = Liouville.initializeCouplingHamiltonianMatrix(scheme, levels)
    
    if  computation.settings.printBefore
        Liouville.displayDensityMatrix(stdout, levels, densityM)
        Liouville.displayGenericHamiltonian(stdout, levels, atomicHM, couplingHM)
    end 
    
    # Generate the (complete) light field due to the chosen pulse sequence
    lightField = Liouville.determineLightField(pulses, -10.)
    
    # Evolve the density matrix for the given light field and by using the proper evolution approach
    # Where do we need to evaluate the Hamiltonian with the parameters of the field: Liouville.evaluateHamiltonian(.. with current parameters)
    ## finalDM = Liouville.evolveDensityMatrix(scheme, computation.approach, densityM, hamiltonianM, lightField)
    
    # Analyze the final density matrix for its stimulated-Raman spectrum
    ## Liouville.analyzeDensityMatrix()
        
    # Return results if required   ... collect all relevant information into this dictionary
    if  output   
        results["pulses:"]  = computation.pulses               
    end
    println(" ")

    
    println("\n> Stimulated Raman scattering computation complete ...")
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    
    return( results )
end
