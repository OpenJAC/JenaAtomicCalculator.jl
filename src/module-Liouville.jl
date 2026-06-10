
"""
`module  JAC.Liouville`  
... a submodel of JAC that contains all methods to set-up and process Liouville-evolution computations 
    and simulations. In these computations, the atomic ensemble is prepared at t0 (= 0) in some initial
    density matrix and later interacts with one or several pulses until conclusion can be drawn about 
    some intermediate of final density matrix of the ensemble. Obviously, the time t0 does not occur
    explicitly in the formalism. Instead, the interaction of the atom with the pulses is controlled by
    the delayTime of each pulse, and the time evolution need to be done long enough that all the pulses
    have interacted with the atoms.
"""
module Liouville

using  Dates, JLD2, Printf
using  ..AtomicState, ..Basics, ..BsplinesN, ..Defaults, ..ManyElectron, ..Nuclear, ..Pulse, ..Radial, ..RadialIntegrals, 
       ..SelfConsistent, ..TableStrings


"""
`abstract type Liouville.AbstractLiouvilleScheme` 
    ... defines an abstract type to distinguish different kinds of Liouville computations; see also:
    
    + struct StimulatedRamanScheme    
        ... to perform a stimulate Raman scattering computation based on the time-evolution of the
            density matrix in light field of one or several pulses.
"""
abstract type  AbstractLiouvilleScheme       end


"""
`struct  Liouville.StimulatedRamanFelScheme   <:  Liouville.AbstractLiouvilleScheme`  
    ... a struct a stimulate Raman scattering computation at an FEL. In a stimulated Raman scheme at an FEL, 
        atoms interaction with a (stimulated Raman-FEL) light field, which is based on one or several 
        linearly-polarized pulses with well-defined amplitudes A0, FWHM in time as well as selected energy spread 
        deltaE. These spreads are assumed to be much larger than those due to the energy-time uncertainty. 
        Other properties of the pulses, such as the CEP or special polarizations, are not taken into account.
        In the simplest case, such a stimulate Raman scattering computation at an FEL is based on five levels
        (g -- ground, e -- inner-shell excited, s -- stimulated decay level, r- spontaneus decay level,
        a -- all autoionizing and loss levels). 

    + levelSelection        ::LevelSelection       
        ... to specify the levels (in the multiplet of refConfigs) that are to be involved in the computation.
    + levelNotations        ::Array{String,1} 
        ... vector of level notations to introduce a proper name beyond the leading configuration;
            this should include a notation of the spontanous and loss channel.
    + gammaR                ::Float64              
        ... Loss rate gammaR(e --> r) of the inner-shell excited level due to the spontaneous decay.
    + gammaA                ::Float64 
        ... Loss rate gammaA(e --> a) of the inner-shell excited level due to autoionization processes.
    + calcPopulations       ::Bool                 ... True, if stimulated Raman lines are to be calculated.
"""
struct   StimulatedRamanScheme   <:  Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection 
    levelNotations          ::Array{String,1} 
    gammaR                  ::Float64
    gammaA                  ::Float64
    calcPopulations         ::Bool
end  


"""
`Liouville.StimulatedRamanScheme()`  ... constructor for an 'default' instance of a Liouville.StimulatedRamanScheme.
"""
function StimulatedRamanScheme()
    StimulatedRamanScheme( LevelSelection(false), String[], 0., 0., false )
end


# `Base.string(scheme::StimulatedRamanScheme)`  ... provides a String notation for the variable scheme::StimulatedRamanScheme.
function Base.string(scheme::StimulatedRamanScheme)
    sa = "Stimulated Raman scattering at FEL computation:\n"
    return( sa )
end


# `Base.show(io::IO, scheme::StimulatedRamanScheme)`  ... prepares a proper printout of the scheme::StimulatedRamanScheme.
function Base.show(io::IO, scheme::StimulatedRamanScheme)
    sa = Base.string(scheme);                print(io, sa)
    println(io, "levelSelection:          $(scheme.levelSelection)  ")
    println(io, "levelNotations:          $(scheme.levelNotations)  ")
    println(io, "gammaR:                  $(scheme.gammaR)  ")
    println(io, "gammaA:                  $(scheme.gammaA)  ")
    println(io, "calcPopulations:         $(scheme.calcPopulations)  ")
end


"""
`abstract type Liouville.AbstractEvolutionApproach` 
    ... defines an abstract type to distinguish different approaches of Liouville time-evolutions; see also:
    
    + struct FirstOrderTimeEvolution    
        ... to perform a simple first-order time-evolution of the density matrix; this mainly serves for tests.
"""
abstract type  AbstractEvolutionApproach       end
struct   FirstOrderTimeEvolution   <:  Liouville.AbstractEvolutionApproach   end  



"""
`struct  Liouville.Settings`  ... defines a type for the details and parameters of Liouville computation; not yet worked out.

    + printBefore               ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
"""
struct Settings 
    printBefore                 ::Bool     
end 


"""
`Liouville.Settings()`  ... constructor for the default values of Liouville computations
"""
function Settings()
    Settings(false)
end


# `Base.show(io::IO, settings::Liouville.Settings)`  ... prepares a proper printout of the variable settings::Liouville.Settings.
function Base.show(io::IO, settings::Liouville.Settings) 
    println(io, "printBefore:         $(settings.printBefore)  ")
end


"""
`struct  Liouville.Computation`  
    ... defines a type for defining  Liouville evolution computations for the interaction of atoms with a light field.

    + scheme                         ::Liouville.AbstractLiouvilleScheme    ... Scheme (kind) of Liouville evolution computation.
    + approach                       ::Liouville.AbstractEvolutionApproach  ... Approach used to solve the time evolution of the density matrix.
    + freeTime                       ::Float64                              
        ... time in [a.u.] that is chosen for free propagation before the first pulse arrives and after the last pulse has left the atoms.
            In practice, each pulse has an internally chosen final length, and the time before the first and after the last pulse are cosidered
            to follow a trivial time evolutions just due to the atomic Hamiltonian.
    + pulses                         ::Vector{Pulse.AbstractPulse}          ... List of pulses that act upon the atom.
    + nuclearModel                   ::Nuclear.Model                        ... Model, charge and parameters of the nucleus.
    + grid                           ::Radial.Grid                          ... The radial grid to be used for the computation.
    + refConfigs                     ::Array{Configuration,1}               ... A list of non-relativistic configurations.
    + asfSettings                    ::AsfSettings                     
        ... Provides the settings for the SCF process (under Liouville-evolution conditions) and the associated CI calculations.
    + settings                       ::Liouville.Settings              ... communicates general setting of these Liouville computations.
"""
struct  Computation
    scheme                           ::AbstractLiouvilleScheme 
    approach                         ::Liouville.AbstractEvolutionApproach
    freeTime                         ::Float64 
    pulses                           ::Vector{Pulse.AbstractPulse}
    nuclearModel                     ::Nuclear.Model
    grid                             ::Radial.Grid
    refConfigs                       ::Array{Configuration,1}
    asfSettings                      ::AsfSettings                     
    settings                         ::Liouville.Settings
end 


"""
`Liouville.Computation()`  ... constructor for an 'empty' instance::Liouville.Computation.
"""
function Computation()
    Computation(StimulatedRamanScheme(), Liouville.FirstOrderTimeEvolution(), 0., Pulse.AbstractPulse[], Nuclear.Model(1.), Radial.Grid(), 
                Configuration[], AsfSettings(), Liouville.Settings() )
end


"""
`Liouville.Computation(comp::Liouville.Computation;`

    scheme=..,                  approach=..,                  freeTime=..,              pulses=..,                  
    nuclearModel=..,            grid=..,                      refConfigs=..,            asfSettings=..,               
    settings=..,                printout::Bool=false)
                    
    ... constructor for modifying the given Liouville.Computation by 'overwriting' the previously selected parameters.
"""
function Computation(comp::Liouville.Computation;
    scheme::Union{Nothing,Liouville.AbstractLiouvilleScheme}=nothing,           approach::Union{Nothing,Liouville.AbstractEvolutionApproach}=nothing, 
    freeTime::Union{Nothing,Float64}=nothing,                                   pulses::Union{Nothing,Vector{Pulse.AbstractPulse}}=nothing,           
    nuclearModel::Union{Nothing,Nuclear.Model}=nothing,                         grid::Union{Nothing,Radial.Grid}=nothing,      
    refConfigs::Union{Nothing,Array{Configuration,1}}=nothing,                  asfSettings::Union{Nothing,AsfSettings}=nothing, 
    settings::Union{Nothing,Liouville.Settings}=nothing, 
    printout::Bool=false)
    
    if  isnothing(scheme)            schemex            = comp.scheme            else  schemex                  = scheme                   end 
    if  isnothing(approach)          approachx          = comp.approach          else  approachx                = approach                 end 
    if  isnothing(freeTime)          freeTimex          = comp.freeTime          else  freeTimex                = freeTime                 end 
    if  isnothing(pulses)            pulsesx            = comp.pulses            else  pulsesx                  = pulses                   end 
    if  isnothing(nuclearModel)      nuclearModelx      = comp.nuclearModel      else  nuclearModelx            = nuclearModel             end 
    if  isnothing(grid)              gridx              = comp.grid              else  gridx                    = grid                     end 
    if  isnothing(refConfigs)        refConfigsx        = comp.refConfigs        else  refConfigsx              = refConfigs               end 
    if  isnothing(asfSettings)       asfSettingsx       = comp.asfSettings       else  asfSettingsx             = asfSettings              end 
    if  isnothing(settings)          settingsx          = comp.settings          else  settingsx                = settings                 end 
    
    
    cp = Computation(schemex, approachx, freeTimex, pulsesx, nuclearModelx, gridx, refConfigsx, asfSettingsx, settingsx) 
                        
    if printout  Base.show(cp)      end
    return( cp )
end


"""
`Liouville.Computation( ... example for a Liouville evolution computations)`  

        grid     = Radial.Grid(true)
        nuclearM = Nuclear.Model(18., "Fermi")
        ...
        refConfigs  = [Configuration("[Ne]", Configuration("1s 2s^2 2p^6 3p")]
        Liouville.Computation(Liouville.Computation(), grid=grid, nuclearModel=nuclearM, refConfigs=refConfigs, asfSettings=... )
    
    ... These simple examples can be further improved by overwriting the corresponding parameters.
"""
function Computation(wa::Bool)    
    Liouville.Computation()    
end


# `Base.string(comp::Liouville.Computation)`  ... provides a String notation for the variable comp::Liouville.Computation.
function Base.string(comp::Liouville.Computation)
    sa = "Liouville evolution computation:  for Z = $(comp.nuclearModel.Z), "
    return( sa )
end


# `Base.show(io::IO, comp::Liouville.Computation)`  ... prepares a printout of comp::Liouville.Computation.
function Base.show(io::IO, comp::Liouville.Computation)
    sa = Base.string(comp);             print(io, sa, "\n")
    println(io, "scheme:                $(comp.scheme)  ")
    println(io, "approach:              $(comp.approach)  ")
    println(io, "pulses:                $(comp.pulses)  ")
    println(io, "nuclearModel:          $(comp.nuclearModel)  ")
    println(io, "grid:                  $(comp.grid)  ")
    println(io, "refConfigs:            $(comp.refConfigs)  ")
    println(io, "asfSettings:           $(comp.asfSettings)  ")
    println(io, "settings:              $(comp.settings)  ")
end


"""
`struct  Liouville.DensityMatrixData`  
    ... defines a type to collect density matrices at different times.

    + time     ::Float64               ... time
    + dm       ::Array{Float64,2}      ... density matrix.
"""
struct DensityMatrixData 
    time       ::Float64
    dm         ::Array{Float64,2}   
end 


"""
`Liouville.DensityMatrixData()`  ... constructor for the default values of Liouville.DensityMatrixData set
"""
function DensityMatrixData()
    DensityMatrixData( 0., zeros(2,2) )
end


# `Base.show(io::IO, data::Liouville.DensityMatrixData)`  ... prepares a proper printout of the variable settings::Liouville.DensityMatrixData.
function Base.show(io::IO, data::Liouville.DensityMatrixData) 
    println(io, "time:          $(data.time)  ")
    println(io, "dm:            $(data.dm)  ")
end



include("module-Liouville-inc-stimulated-raman.jl")

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Basics.perform(comp::Liouville.Computation)`  
    ... to set-up and perform a Liouville evolution computation that starts from a given set of reference configurations 
        and supports the evolution of a (time-dependent) density matrix due to the interaction with one or several 
        pulses. The results of all individual steps of the computations are printed to screen but nothing is returned 
        otherwise.

`Basics.perform(comp::Liouville.Computation; output::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the Liouville evolution computation but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Liouville.Computation; output::Bool=false)
    Liouville.perform(comp.scheme, comp::Liouville.Computation, output=output)
end


end # module
