
"""
`module  JAC.READI`  
... a submodel of JAC that contains all methods for computing resonant-excitation Auger double-ionization cross sections and rates.
"""
module READI 


using ..ManyElectron, ..DoubleAuger


"""
`struct  READI.Settings`  ... defines a type for the details and parameters in computing resonant-excitation Auger double-ionization 
                                pathways |i(N)>  --> |m(N+1)>  --> |n(N-1)>.

    + electronEnergies   ::Array{Float64,1}                   ... List of impact-energies of the incoming elecgtrons.
    + selectPathways     ::Bool                               ... True if particular pathways are selected for the computations.
    + selectedPathways   ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
    + maxKappa           ::Int64                              ... Maximum kappa value of partial waves to be included.
"""
struct Settings
    electronEnergies     ::Array{Float64,1}
    selectPathways       ::Bool
    selectedPathways     ::Array{Tuple{Int64,Int64,Int64},1}
    maxKappa             ::Int64 
end 


"""
`READI.Settings()`  ... constructor for the default values of resonant-excitation Auger double-ionization settings.
"""
function Settings()
    Settings( Float64[], false, Tuple{Int64,Int64,Int64}[], 0)
end


# `Base.show(io::IO, settings::READI.Settings)`  ... prepares a proper printout of the variable settings::READI.Settings.  
function Base.show(io::IO, settings::READI.Settings) 
    println(io, "electronEnergies:     $(settings.electronEnergies)  ")
    println(io, "selectPathways:       $(settings.selectPathways)  ")
    println(io, "selectedPathways:     $(settings.selectedPathways)  ")
    println(io, "maxKappa:             $(settings.maxKappa)  ")
end


"""
`struct  READI.Channel`  ... defines a type for a resonant-excitation Auger double-ionization channel that specifies 
                                all quantum numbers, phases and amplitudes.

    + excitationChannel  ::ImpactExcitation.Channel      ... Channel that describes the electron-impact excitation process.
    + augerChannel1      ::AutoIonization.PartialWave        ... Channel that describes the first subsequent Auger/autoionization process.
    + augerChannel2      ::AutoIonization.PartialWave        ... Channel that describes the second Auger/autoionization process.
"""
struct  Channel
    excitationChannel    ::ImpactExcitation.Channel
    channel              ::DoubleAuger.Channel
end 


"""
`struct  READI.Pathway`  
    ... defines a type for a resonant-excitation Auger double-ionization pathway that may include the definition of different 
        excitation and double Auger channels and their corresponding amplitudes.

    + initialLevel           ::Level         ... initial-(state, N-electron) level
    + intermediateLevel      ::Level         ... intermediate-(state, N+1 electron) level m
    + finalLevel             ::Level         ... final-(state,  N-1 electron) level
    + electronInEnergy       ::Float64       ... energy of the (incoming) electron
    + electronOutEnergy      ::Float64       ... energy of the (outgoing, scattered) electron
    + electronAugerEnergy1   ::Float64       ... energy of the (first emitted Auger) electron
    + electronAugerEnergy2   ::Float64       ... energy of the (second emitted Auger) electron
    + crossSection           ::Float64       ... total cross section of this pathway
    + channels               ::Array{READI.Channel,1} ... List of channels of this pathway.
"""
struct  Pathway
    initialLevel             ::Level
    intermediateLevel        ::Level
    finalLevel               ::Level
    electronInEnergy         ::Float64
    electronOutEnergy        ::Float64
    electronAugerEnergy1     ::Float64
    electronAugerEnergy2     ::Float64
    crossSection             ::Float64 
    channels                 ::Array{READI.Channel,1}
end 


"""
`READI.Pathway()`  
    ... constructor for an electron-impact excitation Auger double ionization pathway between a specified 
        initial, intermediate and final level.
"""
function Pathway()
    Pathway(Level(), Level(), Level(), 0., 0., 0., 0., 0., false, READI.Channel[] )
end


# `Base.show(io::IO, pathway::READI.Pathway)`  ... prepares a proper printout of the variable pathway::READI.Pathway.
function Base.show(io::IO, pathway::READI.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "electronInEnergy:           $(pathway.electronInEnergy)  ")
    println(io, "electronOutEnergy:          $(pathway.electronOutEnergy)  ")
    println(io, "electronAugerEnergy1:       $(pathway.electronAugerEnergy1)  ")
    println(io, "electronAugerEnergy2:       $(pathway.electronAugerEnergy2)  ")
    println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "channels:                   $(pathway.channels)  ")
end

end # module
