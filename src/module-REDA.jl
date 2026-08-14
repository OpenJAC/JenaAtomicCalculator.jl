
"""
`module  JAC.REDA`  
... a submodel of JAC that contains all methods for computing resonant-excitation (sequential) double-autoionization 
    cross sections and rates.
"""
module REDA 

using ..ManyElectron, ..AutoIonization, ..ImpactExcitation


"""
`struct  REDA.Settings`  ... defines a type for the details and parameters in computing resonant-excitation (sequential) 
                                double-autoionization pathways |i(N)>  --> |m(N+1)>  --> |n(N)>  --> |f(N-1)>.

    + electronEnergies   ::Array{Float64,1}          ... List of impact-energies of the incoming elecgtrons.
    + selectPathways     ::Bool                      ... True if particular pathways are selected for the computations.
    + selectedPathways   ::Array{Tuple{Int64,Int64,Int64,Int64},1}  
        ... List of list of pathways, given by tupels (inital, inmediate1, inmediate2, final).
    + maxKappa           ::Int64                     ... Maximum kappa value of partial waves to be included.
"""
struct Settings
    electronEnergies     ::Array{Float64,1}
    selectPathways       ::Bool
    selectedPathways     ::Array{Tuple{Int64,Int64,Int64,Int64},1}
    maxKappa             ::Int64 
end 


"""
`REDA.Settings()`  ... constructor for the default values of resonant-excitation (sequential) double-autoionization settings.
"""
function Settings()
    Settings( Float64[], false, Tuple{Int64,Int64,Int64,Int64}[], 0)
end


# `Base.show(io::IO, settings::REDA.Settings)`  ... prepares a proper printout of the variable settings::REDA.Settings.  
function Base.show(io::IO, settings::REDA.Settings) 
    println(io, "electronEnergies:     $(settings.electronEnergies)  ")
    println(io, "selectPathways:       $(settings.selectPathways)  ")
    println(io, "selectedPathways:     $(settings.selectedPathways)  ")
    println(io, "maxKappa:             $(settings.maxKappa)  ")
end


"""
`struct  REDA.Channel`  ... defines a type for a resonant-excitation (sequential) double-autoionization channel that specifies 
            all quantum numbers, phases and amplitudes.

    + excitationChannel  ::ImpactExcitation.Channel      ... Channel that describes the electron-impact excitation process.
    + augerChannel1      ::AutoIonization.PartialWave        ... Channel that describes the first subsequent Auger/autoionization process.
    + augerChannel2      ::AutoIonization.PartialWave        ... Channel that describes the second Auger/autoionization process.
"""
struct  Channel
    excitationChannel    ::ImpactExcitation.Channel
    augerChannel1        ::AutoIonization.PartialWave
    augerChannel2        ::AutoIonization.PartialWave
end 


"""
`struct  REDA.Pathway`  
    ... defines a type for a electron-impact excitation (sequential) double-autoionization pathway that may include the definition of 
        different excitation and autoionization channels and their corresponding amplitudes.

    + initialLevel           ::Level           ... initial-(state, N-electron) level
    + intermediateLevelm     ::Level           ... intermediate-(state, N+1 electron) level m
    + intermediateLeveln     ::Level           ... intermediate-(state, N-electron) level n
    + finalLevel             ::Level           ... final-(state,  N-1 electron) level
    + electronInEnergy       ::Float64         ... energy of the (incoming) electron
    + electronOutEnergy      ::Float64         ... energy of the (outgoing, scattered) electron
    + electronAugerEnergy1   ::Float64         ... energy of the (first emitted Auger) electron
    + electronAugerEnergy2   ::Float64         ... energy of the (second emitted Auger) electron
    + crossSection           ::Float64         ... total cross section of this pathway
    + channels               ::Array{REDA.Channel,1}  ... List of channels of this pathway.
"""
struct  Pathway
    initialLevel             ::Level
    intermediateLevelm       ::Level
    intermediateLeveln       ::Level
    finalLevel               ::Level
    electronInEnergy         ::Float64
    electronOutEnergy        ::Float64
    electronAugerEnergy1     ::Float64
    electronAugerEnergy2     ::Float64
    crossSection             ::Float64 
    channels                 ::Array{REDA.Channel,1}
end 


"""
`REDA.Pathway()`  
    ... constructor for an electron-impact excitation (sequential) double-autoionization pathway between a specified initial, two 
        intermediate and final level.
"""
function Pathway()
    Pathway(Level(), Level(), Level(), Level(), 0., 0., 0., 0., 0., false, REDA.Channel[] )
end


# `Base.show(io::IO, pathway::REDA.Pathway)`  ... prepares a proper printout of the variable pathway::REDA.Pathway.
function Base.show(io::IO, pathway::REDA.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevelm:         $(pathway.intermediateLevelm)  ")
    println(io, "intermediateLeveln:         $(pathway.intermediateLeveln)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "electronInEnergy:           $(pathway.electronInEnergy)  ")
    println(io, "electronOutEnergy:          $(pathway.electronOutEnergy)  ")
    println(io, "electronAugerEnergy1:       $(pathway.electronAugerEnergy1)  ")
    println(io, "electronAugerEnergy2:       $(pathway.electronAugerEnergy2)  ")
    println(io, "crossSection:               $(pathway.crossSection)  ")
    println(io, "channels:                   $(pathway.channels)  ")
end

end # module
