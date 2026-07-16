
"""
`module  JAC.Empirical`  
... a submodel of JAC that contains all methods to set-up and process (simple) empirical computations such as 
    electron-impact ionization, charge exchange, etc.
"""
module Empirical


using  ..AtomicState, ..Basics, ..Continuum, ..Defaults, ..Distribution, ..Radial, ..ManyElectron, ..Nuclear,
       ..InteractionStrength, ..ImpactIonization, ..PeriodicTable, ..Semiempirical, ..SelfConsistent

"""
`abstract type Empirical.AbstractEmpiricalApproximation` 
    ... defines an abstract and a number of singleton types to deal with different line profiles, for instance,
        in a plasma.

    + GaussianProfile       ... assumes a Gaussian line profile L^Gaussian (omega)
    + LorentzianProfile     ... assumes a Lorentzian line profile L^Lorentzian (omega)
"""
abstract type  AbstractEmpiricalApproximation                                end
struct     ScaledHydrogenic              <:  AbstractEmpiricalApproximation  end
struct     UsingJAC                      <:  AbstractEmpiricalApproximation  end

struct     Bohr1913                      <:  AbstractEmpiricalApproximation  end
struct     Bethe1931                     <:  AbstractEmpiricalApproximation  end
struct     Axelrod1980                   <:  AbstractEmpiricalApproximation  end
struct     KozmaFranson1992              <:  AbstractEmpiricalApproximation  end

struct     Lotz1967                      <:  AbstractEmpiricalApproximation  end
struct     VanRegemorter1962             <:  AbstractEmpiricalApproximation  end


export  AbstractEmpiricalApproximation, GivenEinsteinA, ScaledHydrogenic, UsingJAC,
        Bohr1913, Bethe1931, Axelrod1980, KozmaFranson1992, Lotz1967, VanRegemorter1962
   

"""
`struct  GivenEinsteinA  <:  Empirical.AbstractEmpiricalApproximation`  
    ... to communicate an externally given Einstein-A value, together with the multipole and energy
        of the transition i --> f.

    + multipole   ::EmMultipole     ... Multipole of the transition.
    + energy      ::Float64         ... Energy of transition.
    + rate        ::Float64         ... Einstein-A value.
"""
struct  GivenEinsteinA  <:  Empirical.AbstractEmpiricalApproximation
    multipole     ::EmMultipole
    energy        ::Float64
    rate          ::Float64
end


# `Base.show(io::IO, tripleA::GivenEinsteinA)`  ... provides a String notation for the variable tripleA::GivenEinsteinA.
function Base.show(io::IO, tripleA::GivenEinsteinA)
    sa = "Given Einstein-A value for $(tripleA.multipole) transition with E_if [Hartree] = $(tripleA.energy) " *
         "and A_if [a.u.] = $(tripleA.rate) "
    print(io, sa)
end
   
        

#################################################################################################################################
#################################################################################################################################


"""
`struct  Empirical.Computation`  
    ... defines a type for an empirical computation of various ionization and charge exchange processes.

    + name                   ::String                           ... A name associated to the computation.
    + nuclearModel           ::Nuclear.Model                    ... Model, charge and parameters of the nucleus.
    + grid                   ::Radial.Grid                      ... The radial grid to be used for the computation.
    + configs                ::Array{Configuration,1}           ... A list of non-relativistic configurations.
    + settings               ::Basics.AbstractEmpiricalSettings ... Provides the settings for the selected computations.
"""
struct  Computation
    name                     ::String
    nuclearModel             ::Nuclear.Model
    grid                     ::Radial.Grid
    configs                  ::Array{Configuration,1}
    settings                 ::Basics.AbstractEmpiricalSettings
end 


"""
`Empirical.Computation()`  ... constructor for an 'empty' instance::Empirical.Computation.
"""
function Computation()
    Computation("", Nuclear.Model(1.), Radial.Grid(), Configuration[], ImpactIonization.Settings() )
end


"""
`Empirical.Computation(comp::Empirical.Computation;`

    name=..,                nuclearModel=..,            grid=..,                    configs=..,                   settings=..,  
    printout::Bool=false)
                    
    ... constructor for modifying the given Empirical.Computation by 'overwriting' the previously selected parameters.
"""
function Computation(comp::Empirical.Computation;
    name::Union{Nothing,String}=nothing,               nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
    grid::Union{Nothing,Radial.Grid}=nothing,          configs::Union{Nothing,Array{Configuration,1}}=nothing,       settings::Union{Nothing,Any}=nothing,            
    printout::Bool=false)
    
    if  isnothing(name)                     namex                    = comp.name                    else  namex                    = name                     end 
    if  isnothing(nuclearModel)             nuclearModelx            = comp.nuclearModel            else  nuclearModelx            = nuclearModel             end 
    if  isnothing(grid)                     gridx                    = comp.grid                    else  gridx                    = grid                     end 
    if  isnothing(configs)                  configsx                 = comp.configs                 else  configsx                 = configs                  end 
    if  isnothing(settings)                 settingsx                = comp.settings                else  settingsx                = settings                 end 
    
    
    cp = Computation(namex, nuclearModelx, gridx, configsx, settingsx)
                        
    if printout  Base.show(cp)      end
    return( cp )
end


# `Base.string(comp::Empirical.Computation)`  ... provides a String notation for the variable comp::Empirical.Computation.
function Base.string(comp::Empirical.Computation)
    sa = "Empirical computation:    $(comp.name) for Z = $(comp.nuclearModel.Z), "
    sa = sa * " with the \nconfigurations:        "
    for  config  in  comp.configs   sa = sa * string(config) * ",  "                end
    sa = sa * "\n and for the process/settings: \n $(comp.settings) "
    return( sa )
end


# `Base.show(io::IO, comp::Empirical.Computation)`  ... prepares a printout of comp::Empirical.Computation.
function Base.show(io::IO, comp::Empirical.Computation)
    sa = Base.string(comp);             print(io, sa, "\n")
    println(io, "nuclearModel:          $(comp.nuclearModel)  ")
    println(io, "grid:                  $(comp.grid)  ")
end



"""
`Basics.perform(computation::Empirical.Computation)`  
    ... to set-up and perform an empirical computation that starts from a given nuclear model and set of configurations,
        and which is mainly controlled by its settings. The results are printed to screen but nothing is returned otherwise.

`Basics.perform(computation::Empirical.Computation; output=true)`  
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the kind
        and specifications of the empirical computation but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(computation::Empirical.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nm          = computation.nuclearModel
    asfSettings = AsfSettings()
    
    if typeof(computation.settings)  in  [ImpactIonization.Settings]
        # Generate an SCF basis for the given configurations to extract the one-particle energies for all shells
        multiplet  = SelfConsistent.performSCF(computation.configs, nm, computation.grid, asfSettings)
        basis      = multiplet.levels[1].basis
        
    else
        error("stop a")
    end
    
    if typeof(computation.settings) == ImpactIonization.Settings
        outcome = ImpactIonization.computeCrossSections(basis, computation.grid, nm, computation.settings)        
        if output    results = Base.merge( results, Dict("EII cross sections:" => outcome) )           end
        #
    else
        error("stop b")
    end
    
    return( results )
end

include("module-Empirical-inc-energies.jl")
include("module-Empirical-inc-plasma-rates.jl")
include("module-Empirical-inc-stopping-powers.jl")

end # module
