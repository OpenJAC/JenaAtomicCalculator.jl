
"""
`module  JAC.Atomic`  
	... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
	    atomic representations.
"""
module Atomic


## using Interact
using   ..Basics, ..Defaults, ..LSjj, ..Radial, ..ManyElectron, ..Nuclear, ..PhotoEmission, ..SelfConsistent
export  perform



"""
`struct  Computation`  
    ... defines a type for defining  (the model of simple) atomic computation of a single multiplet, 
        including the SCF and CI as well as level properties and transition property calculations.

    + name                           ::String                          ... A name associated to the computation.
    + nuclearModel                   ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
    + grid                           ::Radial.Grid                     ... The radial grid to be used for the computation.
    + propertySettings               ::Array{Basics.AbstractPropertySettings,1}  ... List of atomic properties to be calculated.
    + configs                        ::Array{Configuration,1}          ... A list of non-relativistic configurations.
    + asfSettings                    ::AsfSettings                     
        ... Provides the settings for the SCF process and for the CI and QED calculations.
    + initialConfigs                 ::Array{Configuration,1}        
        ... A list of initial-state configurations for some transition property calculation, such as radiative transition, Auger, etc. 
    + initialAsfSettings             ::AsfSettings                     ... Provides the SCF and CI settings for the initial-state multiplet.
    + intermediateConfigs            ::Array{Configuration,1}          ... A list of initial-state configurations.
    + intermediateAsfSettings        ::AsfSettings                     ... Provides the SCF settings for the intermediate-state multiplet.
    + finalConfigs                   ::Array{Configuration,1}          ... A list of final-state configurations.
    + finalAsfSettings               ::AsfSettings                     ... Provides the SCF and CI settings for the final-state multiplet.
    + processSettings                ::Basics.AbstractProcessSettings  ... Provides the settings for the selected process.
"""
struct  Computation
    name                           ::String
    nuclearModel                   ::Nuclear.Model
    grid                           ::Radial.Grid
    propertySettings               ::Array{Basics.AbstractPropertySettings,1}
    configs                        ::Array{Configuration,1}
    asfSettings                    ::AsfSettings
    initialConfigs                 ::Array{Configuration,1} 
    initialAsfSettings             ::AsfSettings
    intermediateConfigs            ::Array{Configuration,1} 
    intermediateAsfSettings        ::AsfSettings
    finalConfigs                   ::Array{Configuration,1}
    finalAsfSettings               ::AsfSettings
    processSettings                ::Basics.AbstractProcessSettings
end 


"""
`Atomic.Computation()`  ... constructor for an 'empty' instance::Atomic.Computation.
"""
function Computation()
    Computation("", Nuclear.Model(1.), Radial.Grid(), Basics.AbstractPropertySettings[], 
                Configuration[], AsfSettings(),
                Configuration[], AsfSettings(),
                Configuration[], AsfSettings(),
                Configuration[], AsfSettings(), PhotoEmission.Settings() )
end


"""
`Atomic.Computation(comp::Atomic.Computation;`

    name=..,                nuclearModel=..,            grid=..,                    configs=..,                   asfSettings=..,
    initialConfigs=..,      initialAsfSettings=..,      intermediateConfigs=..,     intermediateAsfSettings=..,
    finalConfigs=..,        finalAsfSettings=..,        propertySettings=..,
    process=..,             processSettings=..,         printout::Bool=false)
                    
    ... constructor for modifying the given Atomic.Computation by 'overwriting' the previously selected parameters.
"""
function Computation(comp::Atomic.Computation;
    name::Union{Nothing,String}=nothing,                                        nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
    grid::Union{Nothing,Radial.Grid}=nothing,                                   propertySettings::Union{Nothing,Array{Basics.AbstractPropertySettings,1},Any}=nothing,   
    configs::Union{Nothing,Array{Configuration,1}}=nothing,                     asfSettings::Union{Nothing,AsfSettings}=nothing, 
    initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,              initialAsfSettings::Union{Nothing,AsfSettings}=nothing, 
    intermediateConfigs::Union{Nothing,Array{Configuration,1}}=nothing,         intermediateAsfSettings::Union{Nothing,AsfSettings}=nothing, 
    finalConfigs::Union{Nothing,Array{Configuration,1}}=nothing,                finalAsfSettings::Union{Nothing,AsfSettings}=nothing, 
    processSettings::Union{Nothing,Any}=nothing,            
    printout::Bool=false)
    
    if  isnothing(name)                     namex                    = comp.name                    else  namex                    = name                     end 
    if  isnothing(nuclearModel)             nuclearModelx            = comp.nuclearModel            else  nuclearModelx            = nuclearModel             end 
    if  isnothing(grid)                     gridx                    = comp.grid                    else  gridx                    = grid                     end 
    if  isnothing(propertySettings)         propertySettingsx        = comp.propertySettings        else  propertySettingsx        = propertySettings         end 
    if  isnothing(configs)                  configsx                 = comp.configs                 else  configsx                 = configs                  end 
    if  isnothing(asfSettings)              asfSettingsx             = comp.asfSettings             else  asfSettingsx             = asfSettings              end 
    if  isnothing(initialConfigs)           initialConfigsx          = comp.initialConfigs          else  initialConfigsx          = initialConfigs           end 
    if  isnothing(initialAsfSettings)       initialAsfSettingsx      = comp.initialAsfSettings      else  initialAsfSettingsx      = initialAsfSettings       end 
    if  isnothing(intermediateConfigs)      intermediateConfigsx     = comp.intermediateConfigs     else  intermediateConfigsx     = intermediateConfigs      end 
    if  isnothing(intermediateAsfSettings)  intermediateAsfSettingsx = comp.intermediateAsfSettings else  intermediateAsfSettingsx = intermediateAsfSettings  end 
    if  isnothing(finalConfigs)             finalConfigsx            = comp.finalConfigs            else  finalConfigsx            = finalConfigs             end 
    if  isnothing(finalAsfSettings)         finalAsfSettingsx        = comp.finalAsfSettings        else  finalAsfSettingsx        = finalAsfSettings         end 
    if  isnothing(processSettings)          prsx                     = comp.processSettings         else  prsx                     = processSettings          end 
    
    
    cp = Computation(namex, nuclearModelx, gridx, propertySettingsx, configsx, asfSettingsx, initialConfigsx, initialAsfSettingsx,     
                        intermediateConfigsx,  intermediateAsfSettingsx, finalConfigsx, finalAsfSettingsx,        
                        prsx) 
                        
    if printout  Base.show(cp)      end
    return( cp )
end


"""
`Atomic.Computation( ... example for SCF computations)`  

        grid     = Radial.Grid(true)
        nuclearM = Nuclear.Model(18., FermiNucleus())
        settings = AsfSettings(AsfSettings(), selectLevelsCI = true, selectedLevelsCI = [1,2, 4,5, 7,8], jjLS = LSjjSettings(false) )
        configs  = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")]
        Atomic.Computation(Atomic.Computation(), name="Example", grid=grid, nuclearModel=nuclearM, configs=configs, asfSettings=settings )

`Atomic.Computation( ... example for the computation of one atomic process)`

        grid           = Radial.Grid(true)
        initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")]
        finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")] 
        photoSettings  = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb], printBefore=true)
        Atomic.Computation(Atomic.Computation(), name="Example", grid=grid, nuclearModel=nuclearM;
                            initialConfigs=initialConfigs, finalConfigs=finalConfigs, 
                            process = Radiative(), processSettings=photoSettings ); 
    
    ... These simple examples can be further improved by overwriting the corresponding parameters.
"""
function Computation(wa::Bool)    
    Atomic.Computation()    
end


# `Base.string(comp::Atomic.Computation)`  ... provides a String notation for the variable comp::Atomic.Computation.
function Base.string(comp::Atomic.Computation)
    sa = "Atomic computation:    $(comp.name) for Z = $(comp.nuclearModel.Z), "
    if  comp.processSettings != Nothing   sa = sa * "for the process (comp.process) and with the \n\ninitial configurations:       "
        for  config  in  comp.initialConfigs   sa = sa * string(config) * ",  "         end
        if  length(comp.intermediateConfigs) > 0
            sa = sa * "\nintermediate configurations:  "
            for  config  in  comp.intermediateConfigs     sa = sa * string(config) * ",  "     end
        end
        sa = sa * "\nfinal configurations:         "
        for  config  in  comp.finalConfigs     sa = sa * string(config) * ",  "         end
    else                                sa = sa * "for the properties $(comp.propertySettings) and with the \nconfigurations:        "
        for  config  in  comp.configs   sa = sa * string(config) * ",  "                end
    end
    return( sa )
end


# `Base.show(io::IO, comp::Atomic.Computation)`  ... prepares a printout of comp::Atomic.Computation.
function Base.show(io::IO, comp::Atomic.Computation)
    sa = Base.string(comp);             print(io, sa, "\n\n")
    println(io, "nuclearModel:          $(comp.nuclearModel)  ")
    println(io, "grid:                  $(comp.grid)  ")
    #
    # For the computation of some given atomic process
    if  comp.processSettings != Nothing
        println(io, "processSettings:              \n$(comp.processSettings)  ")
    if  comp.initialAsfSettings != AsfSettings()     
        println(io, "initialAsfSettings:           \n$(comp.initialAsfSettings)  ")         end
    if  comp.intermediateAsfSettings != AsfSettings()     
        println(io, "intermediateAsfSettings:      \n$(comp.intermediateAsfSettings)  ")    end
    if  comp.finalAsfSettings != AsfSettings()     
        println(io, "finalAsfSettings:             \n$(comp.finalAsfSettings)  ")           end
    #
    else
    # For the computation of no or several atomic properties
    if  comp.asfSettings != AsfSettings()     
        println(io, "asfSettings:                  \n$(comp.asfSettings)  ")                end
    if  Isotope in comp.properties    &&  comp.isotopeSettings  != IsotopeShift.Settings()
        println(io, "isotopeSettings:              \n$(comp.isotopeSettings)  ")            end
    end
end


"""
`Atomic.warnAboutUnusedSettings(comp::Atomic.Computation)`
    ... warns about every AsfSettings field of the given computation that the path it will actually take can NEVER
        read, and whose value differs from the default. An Atomic.Computation carries FOUR such fields -- asfSettings
        for the `configs` path, and initialAsfSettings / intermediateAsfSettings / finalAsfSettings for the
        initial/final path -- deliberately, so that each multiplet may carry its own correlation model. Only the
        fields of the chosen path are consulted, and a value placed in one of the others is silently ignored: someone
        who sets eeInteractionCI = CoulombBreit(...) in `asfSettings` while giving initialConfigs gets a Coulomb-only
        CI with no indication why.

        THE COMPARISON AGAINST THE DEFAULT IS WHAT MAKES THIS USABLE. Every computation carries all four fields
        whether or not the user touched them, so a check that merely asked "is this path reading the field" would
        fire on every run and be ignored within a week. AsfSettings defines Base.:(==), so `!= AsfSettings()`
        separates "the user set this" from "untouched".

        IT IS A WARNING AND NEVER AN ERROR: one Computation may legitimately be built once and reused across paths.
        And it catches only a field that is NEVER consulted -- it cannot catch someone who fills initialAsfSettings
        where finalAsfSettings was meant, since both are read. Nothing is returned.
"""
function warnAboutUnusedSettings(comp::Atomic.Computation)
    usesConfigs = length(comp.configs) != 0
    unused      = usesConfigs ? [ ("initialAsfSettings",      comp.initialAsfSettings),
                                  ("intermediateAsfSettings", comp.intermediateAsfSettings),
                                  ("finalAsfSettings",        comp.finalAsfSettings) ] :
                                [ ("asfSettings",             comp.asfSettings) ]
    named       = usesConfigs ? "configs"       : "initialConfigs/finalConfigs"
    reads       = usesConfigs ? "asfSettings"   : "initialAsfSettings, intermediateAsfSettings and finalAsfSettings"
    for  (name, set)  in  unused
        if  set != AsfSettings()
            println(">> WARNING: this computation is built from $named, so it reads $reads;  the $name you gave " *
                    "is NOT used and its contents -- an eeInteractionCI or scfRoute set there, for instance -- " *
                    "will not reach the calculation.")
        end
    end

    return( nothing )
end

end # module
