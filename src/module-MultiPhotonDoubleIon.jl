
"""
`module  JAC.MultiPhotonDoubleIon`  
... a submodel of JAC that contains all methods for computing multi-photon two-electron (double) ionization cross sections.
"""
module MultiPhotonDoubleIon


using ..Basics, ..ManyElectron, ..Radial


"""
`struct  MultiPhotonDoubleIon.Settings`  
    ... defines a type for the settings in estimating multi-photon two-electron (double) ionization cross sections.

    + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
    + NoPhotons               ::Int64                        ... Number of photons in the multi-photon ionization
    + photonEnergies          ::Array{Float64,1}             ... List of photon energies.
    + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.

"""
struct Settings  <:  Basics.AbstractProcessSettings 
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    NoPhotons                 ::Int64
    photonEnergies            ::Array{Float64,1} 
    printBefore               ::Bool  
    lineSelection             ::LineSelection
end 


"""
`MultiPhotonDoubleIon.Settings()`  ... constructor for the default values of multi-photon two-electron (double) ionization estimates.
"""
function Settings()
    Settings(EmMultipole[], UseGauge[], 0, Float64[], false, LineSelection() )
end


# `Base.show(io::IO, settings::MultiPhotonDoubleIon.Settings)`  
#		... prepares a proper printout of the variable settings::MultiPhotonDoubleIon.Settings.
function Base.show(io::IO, settings::MultiPhotonDoubleIon.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "NoPhotons:                $(settings.NoPhotons)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end


"""
`struct  MultiPhotonDoubleIon.Channel`  
    ... defines a type for a multi-photon channel to help characterize multi-photon (single-electron) ionization with 
        well-defined multipolarities.

    + NoPhotons      ::EmMultipole            ... Number of photons in the multi-photon ionization
    + multipoles     ::Array{EmMultipole,1}   ... Multipoles of all N incoming photons.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + kappa1         ::Int64                  ... partial-wave of the first free electron
    + kappa2         ::Int64                  ... partial-wave of the second free electron
    + phase1         ::Float64                ... phase of the first partial wave
    + phase2         ::Float64                ... phase of the second partial wave
    + symmetry1      ::LevelSymmetry          ... total angular momentum and parity of scattering state after coupling of first electron
    + symmetry2      ::LevelSymmetry          ... total angular momentum and parity of scattering state after coupling of second electron
    + amplitude      ::Complex{Float64}       ... multi-photon ionization amplitude associated with the given channel.
"""
struct  Channel
    NoPhotons        ::EmMultipole 
    multipoles       ::Array{EmMultipole,1}
    gauge            ::EmGauge 
    kappa1           ::Int64 
    kappa2           ::Int64   
    phase1           ::Float64 
    phase2           ::Float64 
    symmetry1        ::LevelSymmetry 
    symmetry2        ::LevelSymmetry
    amplitude        ::Complex{Float64}
end




# MultiPhotonDoubleIon line between initial and final (bound-state) levels
"""
`struct  MultiPhotonDoubleIon.Line`  ... defines a type for a multi-photon ionization line that may include the definition of channels.

    + initialLevel     ::Level           ... initial-(state) level
    + finalLevel       ::Level           ... final-(state) level
    + NoPhotons        ::Int64           ... Number of photons in the multi-photon process
    + photonEnergy     ::Float64         ... Energy of the incoming photons; all photons are assumed to have equal energy.
    + crossSection     ::EmProperty      ... Cross section for this multi-photon ionization.
    + channels         ::Array{MultiPhotonDoubleIon.Channel,1}  ... List of MultiPhotonDoubleIon.Channels of this line.
"""
struct  Line
    initialLevel       ::Level
    finalLevel         ::Level
    NoPhotons          ::Int64 
    photonEnergy       ::Float64
    crossSection       ::EmProperty
    channels           ::Array{MultiPhotonDoubleIon.Channel,1}
end


# `Base.show(io::IO, line::MultiPhotonDoubleIon.Line)`  ... prepares a proper printout of the variable line::MultiPhotonDoubleIon.Line.
function Base.show(io::IO, line::MultiPhotonDoubleIon.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "NoPhotons:         $(line.NoPhotons)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "channels:          $(line.channels)  ")
end



"""
`MultiPhotonDoubleIon.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                    settings::MultiPhotonDoubleIon.Settings; output=true)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{MultiPhotonDoubleIon.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                        settings::MultiPhotonDoubleIon.Settings; output=true)
    println("")
    printstyled("MultiPhotonDoubleIon.computeLines(): The computation of multi-photon double ionization amplitudes starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------- ---------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    #
    pathways = "Not yet implemented !"
    #
    if    output    return( pathways )
    else            return( nothing )
    end
end

end # module
