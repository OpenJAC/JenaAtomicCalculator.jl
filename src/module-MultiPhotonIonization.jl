
"""
`module  JAC.MultiPhotonIonization`  
... a submodel of JAC that contains all methods for computing multi-photon (single-electron) ionization cross sections.
"""
module MultiPhotonIonization


using  ..AngularMomentum, ..Basics, ..Continuum, ..InteractionStrength, ..ManyElectron, ..Radial

"""
`struct  MultiPhotonIonization.Settings`  
    ... defines a type for the settings in estimating multi-photon (single-electron) ionization cross sections

    + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
    + photonEnergies          ::Array{Float64,1}             ... List of photon energies.
    + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.

"""
struct Settings  <:  Basics.AbstractProcessSettings 
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    photonEnergies            ::Array{Float64,1} 
    printBefore               ::Bool  
    lineSelection             ::LineSelection  
end 


"""
`MultiPhotonIonization.Settings()`  ... constructor for the default values of multi-photon (single-electron) ionization estimates.
"""
function Settings()
    Settings(EmMultipole[], UseGauge[], Float64[], false, LineSelection() )
end


# `Base.show(io::IO, settings::MultiPhotonIonization.Settings)`  
# 	... prepares a proper printout of the variable settings::MultiPhotonIonization.Settings.
function Base.show(io::IO, settings::MultiPhotonIonization.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end


"""
`struct  MultiPhotonIonization.Channel`  
    ... defines a type for a multi-photon channel to help characterize multi-photon (single-electron)ionization with 
        well-defined multipolarities.

    + NoPhotons      ::EmMultipole            ... Number of photons in the multi-photon ionization
    + multipoles     ::Array{EmMultipole,1}   ... Multipoles of all N incoming photons.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + kappa          ::Int64                  ... partial-wave of the free electron
    + symmetry       ::LevelSymmetry          ... total angular momentum and parity of the scattering state
    + phase          ::Float64                ... phase of the partial wave
    + amplitude      ::Complex{Float64}       ... multi-photon ionization amplitude associated with the given channel.
"""
struct  Channel
    NoPhotons        ::EmMultipole 
    multipoles       ::Array{EmMultipole,1}
    gauge            ::EmGauge 
    kappa            ::Int64
    symmetry         ::LevelSymmetry 
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  MultiPhotonIonization.Line`  ... defines a type for a multi-photon ionization line that may include the definition of channels.

    + initialLevel     ::Level                  ... initial-(state) level
    + finalLevel       ::Level                  ... final-(state) level
    + photonEnergy     ::Float64                ... Energy of the incoming photons; all photons are assumed to have equal energy.
    + crossSection     ::EmProperty             ... Cross section for this multi-photon ionization.
    + channels         ::Array{MultiPhotonIonization.Channel,1}  ... List of MultiPhotonIonization.Channels of this line.
"""
struct  Line
    initialLevel       ::Level
    finalLevel         ::Level
    photonEnergy       ::Float64
    crossSection       ::EmProperty
    channels           ::Array{MultiPhotonIonization.Channel,1}
end


# `Base.show(io::IO, line::MultiPhotonIonization.Line)`  ... prepares a proper printout of the variable line::MultiPhotonIonization.Line.
function Base.show(io::IO, line::MultiPhotonIonization.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "channels:          $(line.channels)  ")
end



"""
`MultiPhotonIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                            settings::MultiPhotonIonization.Settings; output=true)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{MultiPhotonIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                        settings::MultiPhotonIonization.Settings; output=true)
    println("")
    printstyled("MultiPhotonIonization.computeLines(): The computation of multi-photon ionization amplitudes starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------- ----------------------------------------------------- \n", color=:light_green)
    println("")
    #
    #
    pathways = "Not yet implemented !"
    #
    if    output    return( pathways )
    else            return( nothing )
    end
end

# Inlcude one-electron code
include("module-MultiPhotonIonization-inc-one-electron.jl")

end # module
