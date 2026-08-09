
"""
`module  JAC.AtomicCompass`  
	... a submodel of JAC that contains all methods to specify, set-up and simulate an "atomic compass" in terms 
        of some time-envolved density matrix for a (quantum-optical) few-level systems, interacting with external 
        fields.
"""
module AtomicCompass


using  .. AngularMomentum, ..Basics, ..Beam, ..Defaults, ..Nuclear

export  AbstractTimeEvoluationMethod


"""
`abstract type AtomicCompass.AbstractObservableType` 
    ... defines an abstract type to specify different observables that can be derived from the density matrix
        at some given time.
        
        !!! I can imagine quite a number of different observables, though I didn't look much into the details 
            of yours or Theo's experimental scheme. These observables (types) will likely be developed as the needs 
            arise from case studies or experiments. From my viewpoint, it is helpful important to specify different 
            types of observables with proper names and parameters in order to ensure a good classification of 
            different experimental scenarios !!!!
    
    + struct NoCompassObservable   ... No observable defined for atomic-compass simulations.
    + struct AbsorptionPattern2D   ... to specify the calculation of an absorption pattern for atoms distributed over
                                       the cross section of the beam.
 """
abstract type  AbstractObservableType                             end
struct     NoCompassObservable     <:  AbstractObservableType     end


# `Base.show(io::IO, obs::AbstractObservableType)`  ... prepares a proper printout of the variable obs::AbstractObservableType.
function Base.show(io::IO, obs::AbstractObservableType) 
    print(io, string(obs) )
end


# `Base.string(obs::AbstractObservableType)`  ... provides a proper printout of the variable obs::AbstractObservableType.
function Base.string(obs::AbstractObservableType) 
    if       obs == NoCompassObservable()     return("No observable defined for atomic-compass simulations.")
    else     error("stop a")
    end
end 


"""
`struct  AtomicCompass.AbsorptionPattern2D   <:  AbstractObservableType`  
    ... to compute the absorption pattern for a given 2D atomic target, i.e. target atoms at selected positions in 2D.

    + nVector      ::Basics.Cartesian3DFieldVector
        ... 3D unit vector that specifies the normal direction of the absorption pattern (2D target).
    + times        ::Array{Float64,1}     ... Lists of times at which the absorption pattern need to be determined.
"""
struct  AbsorptionPattern2D                 <:  AbstractObservableType
    nVector        ::Basics.Cartesian3DFieldVector
    times          ::Array{Float64,1} 
end 


# `Base.show(io::IO, obs::AbsorptionPattern2D)`  ... prepares a proper printout of the variable obs::AbsorptionPattern2D.
function Base.show(io::IO, obs::AbsorptionPattern2D) 
    sa = "Absorption pattern are calculated at times = $(obs.times) for a 2D target with normal nVector = $(obs.nVector)."
    print(io, sa)
end


"""
`abstract type AtomicCompass.AbstractTimeEvoluationMethod` 
    ... defines an abstract type to specify the numerical method for time-evolution of the density matrix.
        !!! Apart from certain predictor-corrector methods, there are likely Julia packages that help
            with this time evolution; in this case, the associated concrete method just need to provide 
            a valid interface to data types to make use of these solvers. !!!!
    
    + struct OneTimeStepMethod    ... to integrate the differential equation with a simple (but inaccurate)
                                      one-step method
 """
abstract type  AbstractTimeEvoluationMethode                                           end
struct     OneTimeStepMethod       <:  AtomicCompass.AbstractTimeEvoluationMethode     end


"""
`struct  AtomicCompass.BasisState`  
    ... to specify a IJFM-coupled basis state of the (quantum-optical) few-level system.

    + name        ::Char             ... Character to name the level ('a', 'e', 'g').
    + totalJ      ::AngularJ64       ... total J of the electronic level.
    + F           ::AngularJ64       ... total F of the hyperfine level.
    + M           ::AngularM64       ... total M of the hyperfine level.
"""
struct  BasisState
    name          ::Char  
    totalJ        ::AngularJ64 
    F             ::AngularJ64
    M             ::AngularM64
end 


# `Base.show(io::IO, state::AtomicCompass.BasisState)`  ... prepares a proper printout of the variable state::AtomicCompass.BasisState.
function Base.show(io::IO, state::AtomicCompass.BasisState) 
    sa = "Basis state (name, J, F, M) = ($(state.name), $(state.totalJ), $(state.F), $(state.M)).";   print(io, sa)
end


"""
`struct  AtomicCompass.Level`  
    ... to specify an electronic level of the (quantum-optical) few-level system.

    + name        ::Char                   ... Character to name the level ('a', 'e', 'g').
    + totalJ      ::AngularJ64             ... total J of the electronic level.
    + Fs          ::Array{AngularJ64,1}    ... Vector of total F values to be considered for this level.
"""
struct  Level
    name          ::Char  
    totalJ        ::AngularJ64 
    Fs            ::Array{AngularJ64,1} 
end 


# `Base.show(io::IO, level::AtomicCompass.Level)`  ... prepares a proper printout of the variable level::AtomicCompass.Level.
function Base.show(io::IO, level::AtomicCompass.Level) 
    sa = "Electronic level $(level.name) with J = $(level.totalJ) and total hyperfine momenta F's = $(level.Fs).";   print(io, sa)
end


"""
`struct  AtomicCompass.Observable`  
    ... to specify an observable for the atomic-compass simulation, based on a given (quantum-optical) few-level system.

    + obsType          ::AtomicCompass.AbstractObservableType   ... type of the observable to be simulated.
    + time             ::Float64                                ... Time at which the observable is determined.
    + distribution1D   ::Array{Float64,1}                       
        ... Vector of values that can be displayed and interpreted in line with the given type of observable.
"""
struct  Observable
    obsType            ::AtomicCompass.AbstractObservableType
    time               ::Float64
    distribution1D     ::Array{Float64,1}
end 


# `Base.show(io::IO, obs::AtomicCompass.Observable)`  ... prepares a proper printout of the variable obs::AtomicCompass.Observable.
function Base.show(io::IO, obs::AtomicCompass.Observable) 
    println(io, "obsType:           $(obs.obsType)  ")
    println(io, "time:              $(obs.time)  ")
    println(io, "distribution1D:    $(obs.distribution1D)  ")
end


"""
`struct  AtomicCompass.Settings`  
    ... to specify the settings for a particular atomic-compass simulation, based on a given (quantum-optical) few-level system.

    + qAxis            ::Basics.AbstractQuantizationAxis    ... Specify a useful quantization axis.
    + timeRange        ::StepRange{Int64, Int64}            
        ... time range for integrating the density matrix [user-specified units].
    + multipoles       ::Array{EmMultipole,1}               ... Specifies the (radiat. field) multipoles to be included.
    + mpAmplitudes     ::Array{Float64,1}                   
        ... Associated list of multipole amplitudes; it is presently assumed that the "plane-wave" multipole amplitudes
            are given as external parameters.
    + laserOmegas      ::Array{Float64,1}                   ..  List of laser frequencies at which lasers are available. 
    + decayGamma       ::Float64                            ... Gamma decay parameter in the Liouville equations.
    + impactVectors    ::Array{Cartesian2DFieldVector,1}    ... Cartesian field vectors at which target atoms are considered.
    + previousDM       ::Array{ComplexF64,2}                
        ... A (previously-generated) density matrix from which the time-evolution can be continued; this idea is not 
            yet worked out and will likely require several matrices at different times and together with the specification 
            of the basis states to be able to decide, whether something useful is "coming" to the code.
            This could be anything which is expensive to compute ... so that a restart from scratch is usually not desirable.
"""
struct  Settings
    qAxis              ::Basics.AbstractQuantizationAxis
    timeRange          ::StepRange{Int64, Int64}            
    multipoles         ::Array{EmMultipole,1} 
    mpAmplitudes       ::Array{Float64,1}                   
    laserOmegas        ::Array{Float64,1}                   
    decayGamma         ::Float64    
    impactVectors      ::Array{Cartesian2DFieldVector,1} 
    previousDM         ::Array{ComplexF64,2}                
end 


"""
`AtomicCompass.Settings()`  ... constructor for setting the default values.
"""
function Settings()
    Settings(Basics.NoQuantizationAxis(), 0:1:1, [E1], [0.], Float64[], 0., 
             Cartesian2DFieldVector[], Matrix{ComplexF64}(undef, 0,0) )
end


# `Base.show(io::IO, settings::AtomicCompass.Settings)`  ... prepares a proper printout of the variable settings::AtomicCompass.Settings.
function Base.show(io::IO, settings::AtomicCompass.Settings) 
    println(io, "qAxis:              $(settings.qAxis)  ")
    println(io, "timeRange:          $(settings.timeRange)  ")
    println(io, "multipoles:         $(settings.multipoles)  ")
    println(io, "mpAmplitudes:       $(settings.mpAmplitudes)  ")
    println(io, "laserOmegas:        $(settings.laserOmegas)  ")
    println(io, "decayGamma:         $(settings.decayGamma)  ")
    println(io, "impactVectors:      $(settings.impactVectors)  ")
    println(io, "previousDM:         $(settings.previousDM)  ")
end


"""
`struct  AtomicCompass.Simulation`  
    ... to define and specify the "scenario" for an atomic-compass simulation, based on a given (quantum-optical) few-level system.

    + levels           ::Array{AtomicCompass.Level,1}       ... List of atomic levels of the few-level system.
    + photonBeam       ::Beam.PhotonBeam                    ... Photon beam that interacts with the few-level system.
    + nuclearModel     ::Nuclear.Model                      ... Nuclear model of the atomic few-level system.
    + magneticFields   ::Array{Basics.AbstractEmField,1}    ... List of magnetic fields that affect the few-level system.
    + observables      ::Array{AtomicCompass.AbstractObservableType,1}   
        ... List of observables to be calculated for the given few-level system.
    + settings         ::AtomicCompass.Settings             ... Settings to provided detailed physical parameters for the simulation.   
"""
struct  Simulation
    levels             ::Array{AtomicCompass.Level,1} 
    photonBeam         ::Beam.PhotonBeam  
    nuclearModel       ::Nuclear.Model  
    magneticFields     ::Array{Basics.AbstractEmField,1}  
    observables        ::Array{AtomicCompass.AbstractObservableType,1}   
    settings           ::AtomicCompass.Settings  
end 


"""
`AtomicCompass.Simulation()`  ... constructor for setting the default values.
"""
function Simulation()
    Simulation( AtomicCompass.Level[], Beam.PhotonBeam(), Nuclear.Model(1.0), Basics.AbstractEmField[],
                AtomicCompass.AbstractObservableType[], AtomicCompass.Settings() )
end


# `Base.show(io::IO, sim::AtomicCompass.Simulation)`  ... prepares a proper printout of the variable sim::AtomicCompass.Simulation.
function Base.show(io::IO, sim::AtomicCompass.Simulation) 
    println(io, "levels:              $(sim.levels)  ")
    println(io, "photonBeam:          $(sim.photonBeam)  ")
    println(io, "nuclearModel:        $(sim.nuclearModel)  ")
    println(io, "magneticFields:      $(sim.magneticFields)  ")
    println(io, "magneticFields:      $(sim.magneticFields)  ")
    println(io, "settings:            $(sim.settings)  ")
end


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`AtomicCompass.checkDataConsistency(simulation::AtomicCompass.Simulation)`  
    ... Check the consistency of the input data and terminate with some proper error message,
        if certain basic restrictions/limitations are not fullfilled. Nothing is returned.
        
        I suggest to check explicitly for a few such basic restrictions/limitations to quickly recognize whether 
        your program can or cannot deal with certain input values; in Julia, this is NOT absolutely necessary, 
        since these limitations should be recognized also by the available code itself, as it will not find the right code
        for rather arbitrary input. Still, I find this check helpful here ... as we do (often) not well understand
        right from the beginning, how complex such simulations might be "advanced" eventually. A termination of the 
        code by such a check "confirms" that you have thought reasonable careful about the "allowed input".
        
        However, don't make too strong restrictions by this check ... as we wish to handle a good class of different 
        atomic-compass systems with this code.
"""
function checkDataConsistency(simulation::AtomicCompass.Simulation)
    # Check for something silly
    if      3 == 4
        error("Input requests 3 == 4; improve value ...")
    elseif  3 == 5
        error("Input requests 3 == 5; improve value ...")
    end
    
    return( nothing )
end    


"""
`AtomicCompass.computeObservables(observables::Array{AtomicCompass.Observable,1}, 
                                  bStates::Array{AtomicCompass.BasisState,1}, settings::AtomicCompass.Settings)`  
    ... computes all requested observables in a -- more or less -- efficient way. The number and presentation of these
        observables have beed derived/determined before from the given input. Here, it need to be organized that
        these observables are either computed in turn or by some other sequence/organization so that time and memory
        can be utilized more efficient. While, typically, several observables are computed by a single simulation,
        the total number (1..20) will be likely moderate. We might also consider a splitting between the solution
        of a set of density matrices and the (independent) evaluation of the requested observables.
        
        At the beginning, I suggest to start with typically a moderate number and kind of observables that are returned
        as results of the simulation, and from which later different figures can be constructed externally.
        
        The computation of the observables includes the initialization and set-up of the Hamiltonian and density matrices, 
        the selection/specification of proper atomic targets, the time-evolution of the density matrix as well as 
        the extraction of the relevant results that are communicated to the user.
        
        
        I'm aware that this is the "central" part of the atomic-compass simulations ... and I'm aware that you need
        to think rather careful in which sequence different (parts of these) observables are generate. However, I wish to 
        push you into the "computational way" (of thinking), namely that parts of the computations can typically be re-used 
        for different observables if you organize your code and data types properly. --- Let's talk a bit more what are 
        standard simulations ... and what might come in the near future to make your work efficient, please.
        Or, in short, how to make the super-duper atomic-compass code !!!
        
        A list of (new) observables::::Array{AtomicCompass.Observable,1} is returned which now contains all requested results.
"""
function computeObservables(observables::Array{AtomicCompass.Observable,1}, 
                            bStates::Array{AtomicCompass.BasisState,1}, settings::AtomicCompass.Settings)
    newObservables = AtomicCompass.Observable[]
    
    @warn("computeObservables() ... not yet implemented.")
    
    # Here comes some major organization of the code !! Loop through all observables and generate the requested data.
    
    # Push the results to newObservables
    # obs = AtomicCompass.Observable(...)
    # push(newObservables, obs)

    return( newObservables )
end


"""
`AtomicCompass.computeRelexationMatrix(time::Float64, aState::BasisState, bState::BasisState, dm::Array{ComplexF64,2})`  
    ... computes the relaxation matrix element R_ab (...)
        A (complex) transition me::ComplexF64 is returned.
"""
function computeRelexationMatrix(time::Float64, aState::BasisState, bState::BasisState, dm::Array{ComplexF64,2})
    me = ComplexF64(0.)
    
    @warn("computeRelexationMatrix() ... not yet implemented.")
    
    # I do not understand well how "complex" this relaxation matrix might become in different models ... though it is
    # generally a matrix (and not a constant); for the beginning, simply start with your simple model.
    # Perhaps, you will need the settings::AtomicCompass.Settings as well as input parameter to control how this
    # relaxation matrix is generated and utilized..

    return( me )
end


"""
`AtomicCompass.computeTransitionME(beam::Beam.PhotonBeam, aState::BasisState, bState::BasisState,  nm::Nuclear.Model,
                                   mp::EmMultipole)`  
    ... computes the transition matrix element V_ab^(beam) (mp, polarization, phiB, I, ...)
        A (complex) transition me::ComplexF64 is returned.
"""
function computeTransitionME(beam::Beam.PhotonBeam, aState::BasisState, bState::BasisState,  nm::Nuclear.Model,
                             mp::EmMultipole)
    me = ComplexF64(0.)
    
    @warn("computeTransitionM() ... not yet implemented.")
    
    # Most (all) of the parameters should be accessible from the input parameters
    # Please, implement/search for Julia package to get the Wigner matrices.
    # Bessel functions have been used at several places in StrongField, ParticleScattering, ... but, please,
    # compare with Mathematica also.
    # Please, use functions from AngularMomentum module to compute phases and [..] brackets with AngularJ64 values.
    # Let's use the physical data types as long as useful ... to better recognize possible bugs and to ensure a 
    # consistency of the code.

    return( me )
end


"""
`AtomicCompass.computeTotalMagneticField(magneticFields::Array{Basics.AbstractEmField,1}, time::Float64)`  
    ... compute the total magnetic field vector at the given time.
        A vector::Basics.Cartesian3DFieldVector is returned.
"""
function computeTotalMagneticField(magneticFields::Array{Basics.AbstractEmField,1}, time::Float64)
    totalB = Basics.Cartesian3DFieldVector(1., 0., 0.)
    
    @warn("computeTotalMagneticField() ... not yet implemented.")
    
    # Make use of the time and the field definitions ... to construct a total B-field.
    # Typically write some error("...message ...") if you will get parameters which you do not expect
    # or cannot be handled properly ... specific details can always be thought later, if we shall ever come
    # to this point. It is a key of good physics programming to find a good balance between having reasonable
    # general tools ... without getting lost in details.

    return( totalB )
end


"""
`AtomicCompass.convertTimeRange(userRange::StepRange{Int64, Int64})`  
    ... Convert the time-range from user-defined units into atomic times; since the atomic time-unit
        is rather short, you can make a useful conversion to have not to many time steps.
        We could also use a NoTimeStep parameter in the settings ... to get some simple control
        at how many places you wish to compute the density matrix.
        
        Perhaps, we should "allow" and indicate that other units (like msec) are applied through in this module
        as the atomic time unit is vey ry short ... although this would contradict the basic philosophy of JAC. 
        Let's discuss this openly.
        
        A timeRange::StepRange{Int64, Int64} is returned.
"""
function convertTimeRange(userRange::StepRange{Int64, Int64})
    timeRange = 0:10000:10000000  ## This could be [a.u.] here

    @warn("convertTimeRange() ... not yet implemented.")
    
    # Make the conversion explicit; use functions from Defaults.getUnits(...) to understand of what comes "in".
   
    return( timeRange )
end    


"""
`AtomicCompass.determineObservables(simulation::AtomicCompass.Simulation)`  
    ... determine all the observables that need to be simulated. This means that the requested observables (types) and 
        the given fields and targets are brought together to understand, which detailed computations/simulations need 
        to be carried out. It is expected that every simulation results in about 1--20 individual observables which
        give rise to some large but still moderate amount of data each. Here, we make only the initialization of the 
        observables, while their explicit computation starts later
        
        A list of observables::::Array{AtomicCompass.Observable,1} but without the final data; this list only specifies
        which data are requested in this simulation.
"""
function determineObservables(simulation::AtomicCompass.Simulation)
    observables = AtomicCompass.Observable[]
    
    @warn("determineObservables() ... not yet implemented.")
    
    # Here comes some major organization of the code !!
    
    # Push the results to newObservables
    # obs = AtomicCompass.Observable(...)
    # push(newObservables, obs)

    return( observables )
end



"""
`AtomicCompass.generateIJFMbasis(levels::Array{AtomicCompass.Level,1}, nm::Nuclear.Model)`  
    ... generates the IJFM basis states for the evolution of the density matrix; the sequence of these basis states
        defines the order of raws/columns in the density matrix. A list of basisStates::Array{AtomicCompass.BasisState,1}
        is returned.
"""
function generateIJFMbasis(levels::Array{AtomicCompass.Level,1}, nm::Nuclear.Model)
    basisStates = AtomicCompass.BasisState[]
    
    @warn("generateIJFMbasis() ... not yet implemented.")
    
    # Loop through all levels and level.Fs, M ... to generate the basisStates
    # push!(basisStates, AtomicCompass.BasisState(level.name, ...))

    return( basisStates )
end


"""
`Basics.perform(simulation::AtomicCompass.Simulation; output=true)`  
    ... performs the (full) simulation with the given physical and technical parameters; 
        all intermediate and final results are printed to screen. 
        A results::Dict{String,Any} is returned, if output=true, and nothing otherwise.
"""
function Basics.perform(simulation::AtomicCompass.Simulation; output=true)
    results = Dict{String,Any}()
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    
    println("")
    printstyled("The atomic-compass simulation starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------- \n", color=:light_green)
    println("")
    
    # Check for data consistency
    AtomicCompass.checkDataConsistency(simulation)
    
    # Generate the basis states, ...
    bStates = AtomicCompass.generateIJFMbasis(simulation.levels, simulation.nuclearModel)
    
    # Determine and display the observables (though without explicit computations)
    observables = AtomicCompass.determineObservables(simulation)
    ## AtomicCompass.displayInitializedObservables(observables, simulation.settings) # this often help with coding.

    # Compute the observables; this includes the initialization of the density matrix, its time evolution
    # as well as a proper loop about this evolution process to extract all desired information
    results = AtomicCompass.computeObservables(observables, bStates, simulation.settings)
    ## AtomicCompass.displayFinalObservables(observables, simulation.settings)
    
    if  output    return( results )
    else          return( nothing )
    end
end

end # module



