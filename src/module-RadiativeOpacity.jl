
"""
`module  JAC.RadiativeOpacity`  
... a submodel of JAC that contains all methods for computing radiative -- spectral & mean -- opacities.
"""
module RadiativeOpacity


using Printf, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..TableStrings

"""
`struct  RadiativeOpacity.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing radiative opacities.

    + approach           ::AbstractOpacityApproach         
        ... Determines the approach for calculating radiative opacities: {SincleConf(), ...}.
    + NoExcitations      ::Int64          
        ... Number of excitations with regard to the given (reference) configurations.
    + NoLevels           ::Int64    
        ... maximum number of levels (to terminate long computations in a regular fashion).
    + maxExcitation      ::Float64    
        ... maximum excitation energy w.r.t. the ground level for which levels are considered.
    + deltaLambda        ::Float64          ... Delta Lambda in the expansion opacity.
    + time               ::Float64          ... Time after merger in the expansion opacity.
    + densities          ::Array{Float64,1} ... densities [g/cm^3] in the expansion opacity.
    + temperatures       ::Array{Float64,1} ... temperatures [K] in the expansion opacity.
"""
struct Settings  <:  AbstractPropertySettings 
    approach             ::Int64        
    NoExcitations        ::Int64          
    NoLevels             ::Int64    
    maxExcitation        ::Float64    
    DeltaLambda          ::Float64 
    time                 ::Float64 
    densities            ::Array{Float64,1}
    temperatures         ::Array{Float64,1}
end 


"""
`RadiativeOpacity.Settings()`  
    ... constructor for an `empty` instance of RadiativeOpacity.Settings for the computation of radiative opacities.
"""
function Settings()
    Settings(1, 0, 0, 0., 0., 0., Float64[], Float64[])
end


# `Base.show(io::IO, settings::RadiativeOpacity.Settings)`  ... prepares a proper printout of settings::RadiativeOpacity.Settings.
function Base.show(io::IO, settings::RadiativeOpacity.Settings) 
    println(io, "approach:              $(settings.approach)  ")
    println(io, "NoExcitations:         $(settings.NoExcitations)  ")
    println(io, "NoLevels:              $(settings.NoLevels)  ")
    println(io, "maxExcitation:         $(settings.maxExcitation)  ")
    println(io, "DeltaLambda:           $(settings.DeltaLambda)  ")
    println(io, "time:                  $(settings.time)  ")
    println(io, "densities:             $(settings.densities)  ")
    println(io, "temperatures:          $(settings.temperatures)  ")
end


"""
`struct  RadiativeOpacity.Outcome`  
    ... defines a type to keep the outcome of a opacity computation as well other results.

    + spectralOpacity       ::Float64            ... Spectral opacity.
    + expansionOpacity      ::Float64            ... Expansion opacity.
    + rosselandOpacity      ::Float64            ... Mean Rosseland opacity.
"""
struct Outcome 
    spectralOpacity         ::Float64 
    expansionOpacity        ::Float64 
    rosselandOpacity        ::Float64
end 


"""
`RadiativeOpacity.Outcome()`  
    ... constructor for an `empty` instance of RadiativeOpacity.Outcome for the computation of opacities.
"""
function Outcome()
    Outcome(0., 0., 0.)
end


# `Base.show(io::IO, outcome::RadiativeOpacity.Outcome)`  ... prepares a proper printout of outcome::RadiativeOpacity.Outcome.
function Base.show(io::IO, outcome::RadiativeOpacity.Outcome) 
    println(io, "spectralOpacity:         $(outcome.spectralOpacity)  ")
    println(io, "expansionOpacity:        $(outcome.expansionOpacity)  ")
    println(io, "rosselandOpacity:        $(outcome.rosselandOpacity)  ")
end

## Introduce some sort Opacity.Line to keep the data E_i, E_f, level_i, level_f, J_i, J_f, gf ... for all calculated
## lines ... but not the full representation of amplitudes and levels
## Perhaps, levels::Array{Level,1}  ... to keep a number of levels
## calculate multiplet of configurations
## generate a list of configurations for which levels are taken into account.

end # module
