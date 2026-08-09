
"""
`module  JAC.Plasma`  
... a submodel of JAC that contains all methods to set-up and process (simple) plasma computations 
    and simulations.
"""
module Plasma

using  Dates, JLD2, LinearAlgebra, Printf
using  ..Atomic, ..AtomicState, ..Basics, ..Bsplines, ..Defaults, ..DielectronicRecombination, ..Hamiltonian, ..ImpactExcitation,
       ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
       ..Semiempirical, ..TableStrings, ..FormFactor, ..PhotoEmission, ..PhotoIonization, ..AutoIonization, ..SelfConsistent


"""
`abstract type Plasma.AbstractPlasmaScheme`
    ... defines an abstract type to distinguish different kinds of plasma computations; see also:

    + struct AverageAtomScheme
        ... to perform an average-atom computation.
    + struct LineShiftScheme
        ... to compute the energy shifts and properties of atomic/ionic lines in some selected plasma model.
    + struct SahaBoltzmannScheme
        ... to compute thermodynamic properties of a Saha-Boltzmann LTE mixture.
    + struct SatelliteDiagnosticScheme
        ... to compute a dielectronic-satellite-to-parent-line intensity-ratio Te diagnostic.
    + struct CollisionalRadiativeScheme
        ... to compute the (relative) collisional-radiative population balance among a small set of
            levels of one ion.
"""
abstract type  AbstractPlasmaScheme       end


"""
`struct  Plasma.AverageAtomScheme  <:  Plasma.AbstractPlasmaScheme`  
    ... a struct to perform an average-atom computation.

    + nMax                  ::Int64                ... maximum principal quantum number for the subshells of the AA model.
    + lMax                  ::Int64                ... maximum orbital quantum number for the subshells of the AA model.
    + scField               ::AbstractScField      ... maximum orbital quantum number for the subshells of the AA model.
    + calcPhotoionizationCs ::Bool                 ... True, if photoionization cross sections are to be calculated.
    + calcFormFactor        ::Bool                 ... True, if the form factor need to be calculated.
    + calcScatteringFactor  ::Bool                 ... True, if scattering factor need to be calculated.
    + piSubshells           ::Array{Subshell,1}    ... Bound subshells to be included into the photoionization cross sections.
    + omegas                ::Array{Float64,1}     ... energies for the photoionization & scattering factors [in a.u.].
    + qValues               ::Array{Float64,1}     ... q-values for calculating the form factor [in a.u.].
"""
struct   AverageAtomScheme  <:  Plasma.AbstractPlasmaScheme
    nMax                    ::Int64
    lMax                    ::Int64
    scField                 ::AbstractScField
    calcPhotoionizationCs   ::Bool          
    calcFormFactor          ::Bool             
    calcScatteringFactor    ::Bool              
    piSubshells             ::Array{Subshell,1} 
    omegas                  ::Array{Float64,1} 
    qValues                 ::Array{Float64,1} 
end  


"""
`Plasma.AverageAtomScheme()`  ... constructor for an 'default' instance of a Plasma.AverageAtomScheme.
"""
function AverageAtomScheme()
    AverageAtomScheme( 1, 0, Basics.AaHSField(), false, false, false, Subshell[], Float64[], Float64[] )
end


# `Base.string(scheme::AverageAtomScheme)`  ... provides a String notation for the variable scheme::AverageAtomScheme.
function Base.string(scheme::AverageAtomScheme)
    sa = "Average-atom computation with nMax =$(scheme.nMax) and lMax =$(scheme.lMax):"
    return( sa )
end


# `Base.show(io::IO, scheme::AverageAtomScheme)`  ... prepares a proper printout of the scheme::AverageAtomScheme.
function Base.show(io::IO, scheme::AverageAtomScheme)
    sa = Base.string(scheme);                print(io, sa)
    println(io, "\nscField:                  $(scheme.scField)  ")
    println(io, "calcPhotoionizationCs:      $(scheme.calcPhotoionizationCs)  ")
    println(io, "calcFormFactor:             $(scheme.calcFormFactor)  ")
    println(io, "calcScatteringFactor:       $(scheme.calcScatteringFactor)  ")
    println(io, "piSubshells:                $(scheme.piSubshells)  ")
    println(io, "omegas:                     $(scheme.scField)  ")
    println(io, "qValues:                    $(scheme.qValues)  ")
end


"""
`struct  Plasma.LineShiftScheme  <:  Plasma.AbstractPlasmaScheme`  
    ... defines a type for the details and parameters of computing level energies with plasma interactions.

    + plasmaModel      ::AbstractPlasmaModel          ... Specify a particular plasma model, e.g. ion-sphere, Debye.
    + initialConfigs   ::Array{Configuration,1}       ... List of one or several configurations that define the initial-state multiplet.
    + finalConfigs     ::Array{Configuration,1}       ... List of one or several configurations that define the final-state multiplet.
    + settings         ::AbstractLineShiftSettings    ... Specify the process and settings for which line-shifts need to be computed.
    ## + NoBoundElectrons ::Int64                      ... Effective number of bound electrons.
"""
struct LineShiftScheme  <:  Plasma.AbstractPlasmaScheme
    plasmaModel        ::AbstractPlasmaModel
    initialConfigs     ::Array{Configuration,1}
    finalConfigs       ::Array{Configuration,1}
    settings           ::AbstractLineShiftSettings    
    ## NoBoundElectrons   ::Int64
end 


"""
`Plasma.LineShiftScheme()`  ... constructor for a standard instance of Plasma.LineShiftScheme.
"""
function LineShiftScheme()
    LineShiftScheme( Basics.NoPlasmaModel(), Configuration[], Configuration[], Basics.NoLineShiftSettings() )
end


# `Base.show(io::IO, scheme::Plasma.LineShiftScheme)`  ... prepares a proper printout of the scheme::Plasma.LineShiftScheme.
function Base.show(io::IO, scheme::Plasma.LineShiftScheme)
    println(io, "plasmaModel:            $(scheme.plasmaModel)  ")
    println(io, "initialConfigs:         $(scheme.initialConfigs)  ")
    println(io, "finalConfigs:           $(scheme.finalConfigs)  ")
    println(io, "settings:               $(scheme.settings)  ")
end


"""
`struct  Plasma.SahaBoltzmannScheme  <:  Plasma.AbstractPlasmaScheme`  
    ... a struct to thermodynamic properties of a Saha-Boltzmann LTE mixture..

    + plasmaModel           ::Basics.AbstractPlasmaModel                 
        ... A plasma model that "restricts" the Saha-Boltzmann equilibrium densities by some additional parameters, for instance,
            due to ionization-potential-depression (IPD) or others.
    + calcLTE               ::Bool                 ... True, if the Saha-Boltzmann equilibrium densities should be calculated.         
    + printIonLevels        ::Bool                 ... True, for printing detailed information about all ionic levels.         
    + qRange                ::UnitRange{Int64}     ... Range of charge states q for which ions are taken into account. 
    + maxNoIonLevels        ::Int64                ... (maximum) No of ionic levels for any charge state of the ions in the mixture.
    + NoExcitations         ::Int64                
        ... No of excitations (S, D, T) that are taken into account with regard to the reference (ground) configuration of the ions.
            This number is taken as a second qualifier to characterize the quality of the ionic-level data. Usually, NoExcitations = 1,2.
    + upperShellNo          ::Int64                
        ... upper-most princicpal quantum number n for which orbitals are taken into account into the ionic-level computations;
            this is often chosen upperShellNo= 4..10; if ionic-level data are given, this number decides which of the levels
            are taken into account.
    + isotopicMixture       ::Array{Basics.IsotopicFraction,1}   
        ... List of (non-normlized) isotopic fractions that form the requested mixture; the fractions will first be renormalized
            to 1 in course of the Saha-Boltzmann LTE computations.
    + isotopeFilenames      ::Array{String,1}     
        ... set of files names from which the ionic-level data can be read in for the different isotopes (Z,A) in the mixture.
    
"""
struct  SahaBoltzmannScheme  <:  Plasma.AbstractPlasmaScheme
    plasmaModel             ::Basics.AbstractPlasmaModel                 
    calcLTE                 ::Bool        
    printIonLevels          ::Bool         
    qRange                  ::UnitRange{Int64} 
    maxNoIonLevels          ::Int64 
    NoExcitations           ::Int64                
    upperShellNo            ::Int64                
    isotopicMixture         ::Array{Basics.IsotopicFraction,1}   
    isotopeFilenames        ::Array{String,1}     
end  


"""
`Plasma.SahaBoltzmannScheme()`  ... constructor for an 'default' instance of a Plasma.SahaBoltzmannScheme.
"""
function SahaBoltzmannScheme()
    SahaBoltzmannScheme( Basics.NoPlasmaModel(), false, false, 0:0, 0., 0., 0., 0., IsotopicFraction[], String[] )
end


# `Base.string(scheme::SahaBoltzmannScheme)`  ... provides a String notation for the variable scheme::SahaBoltzmannScheme.
function Base.string(scheme::SahaBoltzmannScheme)
    sa = "Saha-Boltzmann computation with of an ionic mixture:"
    return( sa )
end


# `Base.show(io::IO, scheme::SahaBoltzmannScheme)`  ... prepares a proper printout of the scheme::SahaBoltzmannScheme.
function Base.show(io::IO, scheme::SahaBoltzmannScheme)
    sa = Base.string(scheme);             println(io, sa)
    println(io, "plasmaModel:       $(scheme.plasmaModel)  ")
    println(io, "calcLTE:           $(scheme.calcLTE)  ")
    println(io, "printIonLevels:    $(scheme.printIonLevels)  ")
    println(io, "qRange:            $(scheme.qRange)  ")
    println(io, "maxNoIonLevels:    $(scheme.maxNoIonLevels)  ")
    println(io, "NoExcitations:     $(scheme.NoExcitations)  ")
    println(io, "upperShellNo:      $(scheme.upperShellNo)  ")
    println(io, "isotopicMixture:   $(scheme.isotopicMixture)  ")
    println(io, "isotopeFilenames:  $(scheme.isotopeFilenames)  ")
end


"""
`struct  Plasma.SatelliteDiagnosticScheme  <:  Plasma.AbstractPlasmaScheme`
    ... dielectronic-satellite-to-parent-line intensity-ratio Te diagnostic for one recombining ion
        (computation.refConfigs). The core promotion fromShells-->toShells drives BOTH the parent
        line (electron-impact excited, via ieSettings) and the satellites (dielectronically captured
        into intoShells, via drSettings) -- the same physical excitation, populated two different ways.

    + fromShells   ::Array{Shell,1}   ... Core shell(s) excited; shared by the DR and impact-excitation steps.
    + toShells     ::Array{Shell,1}   ... Core-excitation target shell(s); shared by the DR and impact-excitation steps.
    + intoShells   ::Array{Shell,1}
        ... Captured Rydberg shell(s) for DR; an EXPLICIT list of individually-resolved shells (e.g. n'=3: 3s,3p,3d),
            not a generated high-n range -- only spectrally-resolvable low-n' shells belong here. Same naming AND same
            semantics as Basics.ForDielectronicRecombination.intoShells: every shell given here is treated TOGETHER in
            one combined doubly-excited representation/CI, exactly as ForDielectronicRecombination already works. Cost
            grows sharply with the number of shells combined (measured: 220-665s for a single n=3 shell of F-like Ne+;
            combining all three n=3 shells together exceeded 41 minutes). If cost becomes prohibitive, use fewer shells
            here, or issue several separate Plasma.Computation calls (one per shell or small subset) -- the scheme does
            not make this scoping choice automatically.
    + decayShells  ::Array{Shell,1}   ... (Radiative) stabilization target shell(s) for DR, e.g. [Shell("2p")]. Same
                                          role as Basics.ForDielectronicRecombination.decayShells.
    + drSettings   ::DielectronicRecombination.Settings
        ... Reused as-is (not re-flattened): temperatures, pathwaySelection, calcRateAlpha, corrections,
            augerOperator, multipoles, gauges all come from here.
    + ieSettings   ::ImpactExcitation.Settings
        ... Reused as-is: temperatures, lineSelection, calcRateCoefficient, maxKappa, numElectronEnergies,
            maxEnergyMultiplier, operator all come from here. printBefore is forced to false internally by
            the driver regardless of what is set here (module-ImpactExcitation.jl displayLines() bug workaround).
"""
struct  SatelliteDiagnosticScheme  <:  Plasma.AbstractPlasmaScheme
    fromShells    ::Array{Shell,1}
    toShells      ::Array{Shell,1}
    intoShells    ::Array{Shell,1}
    decayShells   ::Array{Shell,1}
    drSettings    ::DielectronicRecombination.Settings
    ieSettings    ::ImpactExcitation.Settings
end


"""
`Plasma.SatelliteDiagnosticScheme()`  ... constructor for an 'default' instance of a Plasma.SatelliteDiagnosticScheme.
"""
function SatelliteDiagnosticScheme()
    SatelliteDiagnosticScheme( Shell[], Shell[], Shell[], Shell[],
                               DielectronicRecombination.Settings(), ImpactExcitation.Settings() )
end


"""
`Plasma.SatelliteDiagnosticScheme(scheme::Plasma.SatelliteDiagnosticScheme;`

        fromShells=..,        toShells=..,        intoShells=..,      decayShells=..,
        drSettings=..,        ieSettings=..)

    ... constructor for modifying the given Plasma.SatelliteDiagnosticScheme by 'overwriting' the previously
        selected parameters.
"""
function SatelliteDiagnosticScheme(scheme::Plasma.SatelliteDiagnosticScheme;
    fromShells::Union{Nothing,Array{Shell,1}}=nothing,     toShells::Union{Nothing,Array{Shell,1}}=nothing,
    intoShells::Union{Nothing,Array{Shell,1}}=nothing,     decayShells::Union{Nothing,Array{Shell,1}}=nothing,
    drSettings::Union{Nothing,DielectronicRecombination.Settings}=nothing,
    ieSettings::Union{Nothing,ImpactExcitation.Settings}=nothing)

    if  isnothing(fromShells)    fromShellsx  = scheme.fromShells    else  fromShellsx  = fromShells    end
    if  isnothing(toShells)      toShellsx    = scheme.toShells      else  toShellsx    = toShells      end
    if  isnothing(intoShells)    intoShellsx  = scheme.intoShells    else  intoShellsx  = intoShells    end
    if  isnothing(decayShells)   decayShellsx = scheme.decayShells   else  decayShellsx = decayShells   end
    if  isnothing(drSettings)    drSettingsx  = scheme.drSettings    else  drSettingsx  = drSettings     end
    if  isnothing(ieSettings)    ieSettingsx  = scheme.ieSettings    else  ieSettingsx  = ieSettings     end

    SatelliteDiagnosticScheme( fromShellsx, toShellsx, intoShellsx, decayShellsx, drSettingsx, ieSettingsx )
end


# `Base.show(io::IO, scheme::SatelliteDiagnosticScheme)`  ... prepares a proper printout of the scheme::SatelliteDiagnosticScheme.
function Base.show(io::IO, scheme::SatelliteDiagnosticScheme)
    println(io, "fromShells:        $(scheme.fromShells)  ")
    println(io, "toShells:          $(scheme.toShells)  ")
    println(io, "intoShells:        $(scheme.intoShells)  ")
    println(io, "decayShells:       $(scheme.decayShells)  ")
    println(io, "drSettings:        $(scheme.drSettings)  ")
    println(io, "ieSettings:        $(scheme.ieSettings)  ")
end


"""
`struct  Plasma.CollisionalRadiativeScheme  <:  Plasma.AbstractPlasmaScheme`
    ... a small, closed-form (non-iterative) collisional-radiative population balance for one ion: for a
        generated set of its levels, computes how its population is distributed among them (a set of
        fractions summing to 1). Not a full multi-hundred-level CR code -- for N levels, the steady-state
        balance (dn_i/dt = 0) is an N x N linear system, solved directly, once per (Te, Ne) point.
        Answers "how is this ion's population distributed among its levels", not "how much of this ion is
        present" -- the latter (the ionization balance across charge states) is a separate question,
        answered separately (e.g. by Plasma.SahaBoltzmannScheme) if absolute densities are wanted.
        The ion itself (element and charge state) is not a field of this scheme: it is set entirely by
        the enclosing Plasma.Computation's nuclearModel (the element) and refConfigs (whose electron
        count against nuclearModel.Z fixes the charge state), exactly as for every other scheme.

    + NoExcitations     ::Int64
        ... how many electrons may be excited at once, relative to the ion's own ground configuration,
            when building its level set (1 = singles only, 2 = up to doubles, etc.) -- controls how rich
            a level set is generated.
    + upperShellNo      ::Int64
        ... the highest principal quantum number n allowed for an excited electron when building that
            level set -- controls how far up in energy the generated levels reach.
    + ieSettings        ::ImpactExcitation.Settings
        ... controls the electron-collisional (upward, excitation) rates between levels; the reverse
            (de-excitation) direction is obtained from these by detailed balance, not computed separately.
    + peSettings        ::PhotoEmission.Settings
        ... controls the radiative decay (A-value) rates between levels.
    + aiSettings        ::AutoIonization.Settings
        ... RESERVED, CURRENTLY UNUSED (verified 04-Aug-2026): this field is accepted and stored, but the
            driver never reads it and NO autoionization is computed -- every generated level is treated as
            purely bound (radiative + collisional only), as the @warn in Plasma.perform states. Setting it
            therefore changes nothing today. It is kept because the competing autoionization (electron-
            ejecting) decay is genuinely needed for any generated level that lies above the next charge
            state's ionization threshold -- a real possibility once NoExcitations/upperShellNo reach far
            enough up in energy, and the reason example-Je.jl branch b (boron-like Ne5+) cannot be trusted
            for its inner-shell block, which sits ~720 eV above threshold.
    + levelsFilenames   ::Array{String,1}
        ... candidate files (tried in order) that may contain a previously-generated, still-matching level
            set; [] (the default) means always regenerate. A fresh file is auto-written (name printed to
            screen) whenever no candidate matches, for adoption into a later run -- same manual-reuse
            workflow as Plasma.SahaBoltzmannScheme.isotopeFilenames.
    + ratesFilenames    ::Array{String,1}
        ... candidate files (tried in order) that may contain previously-computed, still-matching
            ImpactExcitation collision data (pre-thermal-average, so any scheme.ieSettings.temperatures
            list can reuse them); [] (the default) means always recompute. Same auto-write/manual-reuse
            workflow as levelsFilenames.
    + cacheDirectory    ::String
        ... subdirectory that freshly-generated levelsFilenames/ratesFilenames files are written into
            (created if missing); "" (the default) writes into the current directory. Only governs where
            FRESH files go -- candidates in levelsFilenames/ratesFilenames are looked up as given, so a
            cache file from elsewhere can always be reused directly.
"""
struct  CollisionalRadiativeScheme  <:  Plasma.AbstractPlasmaScheme
    NoExcitations       ::Int64
    upperShellNo        ::Int64
    ieSettings          ::ImpactExcitation.Settings
    peSettings          ::PhotoEmission.Settings
    aiSettings          ::AutoIonization.Settings
    levelsFilenames     ::Array{String,1}
    ratesFilenames      ::Array{String,1}
    cacheDirectory      ::String
end


"""
`Plasma.CollisionalRadiativeScheme()`  ... constructor for an 'default' instance of a Plasma.CollisionalRadiativeScheme.
"""
function CollisionalRadiativeScheme()
    CollisionalRadiativeScheme( 1, 4, ImpactExcitation.Settings(), PhotoEmission.Settings(), AutoIonization.Settings(),
                                String[], String[], "" )
end


"""
`Plasma.CollisionalRadiativeScheme(scheme::Plasma.CollisionalRadiativeScheme;`

        NoExcitations=..,        upperShellNo=..,
        ieSettings=..,           peSettings=..,           aiSettings=..,
        levelsFilenames=..,      ratesFilenames=..,       cacheDirectory=..)

    ... constructor for modifying the given Plasma.CollisionalRadiativeScheme by 'overwriting' the previously
        selected parameters.
"""
function CollisionalRadiativeScheme(scheme::Plasma.CollisionalRadiativeScheme;
    NoExcitations::Union{Nothing,Int64}=nothing,            upperShellNo::Union{Nothing,Int64}=nothing,
    ieSettings::Union{Nothing,ImpactExcitation.Settings}=nothing,
    peSettings::Union{Nothing,PhotoEmission.Settings}=nothing,
    aiSettings::Union{Nothing,AutoIonization.Settings}=nothing,
    levelsFilenames::Union{Nothing,Array{String,1}}=nothing,
    ratesFilenames::Union{Nothing,Array{String,1}}=nothing,
    cacheDirectory::Union{Nothing,String}=nothing)

    if  isnothing(NoExcitations)     NoExcitationsx   = scheme.NoExcitations     else   NoExcitationsx   = NoExcitations     end
    if  isnothing(upperShellNo)      upperShellNox    = scheme.upperShellNo      else   upperShellNox    = upperShellNo      end
    if  isnothing(ieSettings)        ieSettingsx      = scheme.ieSettings        else   ieSettingsx      = ieSettings        end
    if  isnothing(peSettings)        peSettingsx      = scheme.peSettings        else   peSettingsx      = peSettings        end
    if  isnothing(aiSettings)        aiSettingsx      = scheme.aiSettings        else   aiSettingsx      = aiSettings        end
    if  isnothing(levelsFilenames)   levelsFilenamesx = scheme.levelsFilenames   else   levelsFilenamesx = levelsFilenames   end
    if  isnothing(ratesFilenames)    ratesFilenamesx  = scheme.ratesFilenames    else   ratesFilenamesx  = ratesFilenames    end
    if  isnothing(cacheDirectory)    cacheDirectoryx  = scheme.cacheDirectory    else   cacheDirectoryx  = cacheDirectory    end

    CollisionalRadiativeScheme( NoExcitationsx, upperShellNox, ieSettingsx, peSettingsx, aiSettingsx,
                                levelsFilenamesx, ratesFilenamesx, cacheDirectoryx )
end


# `Base.show(io::IO, scheme::CollisionalRadiativeScheme)`  ... prepares a proper printout of the scheme::CollisionalRadiativeScheme.
function Base.show(io::IO, scheme::CollisionalRadiativeScheme)
    println(io, "NoExcitations:     $(scheme.NoExcitations)  ")
    println(io, "upperShellNo:      $(scheme.upperShellNo)  ")
    println(io, "ieSettings:        $(scheme.ieSettings)  ")
    println(io, "peSettings:        $(scheme.peSettings)  ")
    println(io, "aiSettings:        $(scheme.aiSettings)  ")
    println(io, "levelsFilenames:   $(scheme.levelsFilenames)  ")
    println(io, "ratesFilenames:    $(scheme.ratesFilenames)  ")
    println(io, "cacheDirectory:    $(scheme.cacheDirectory)  ")
end


"""
`struct  Plasma.Settings`  ... defines a type for the details and parameters of computing photoionization lines.

    + temperature               ::Float64
        ... Plasma temperature. Unit convention is scheme-dependent: Plasma.perform(scheme::AverageAtomScheme, ...)
            expects temperature in [K] and converts it internally; Plasma.perform(scheme::SahaBoltzmannScheme, ...)
            expects temperature already in atomic (Hartree) units, cf. example-Jc.jl. Watch this when re-using
            the same Settings across different schemes.
    + density                   ::Float64     ... Plasma density in [g/cm^3].
    + useNumberDensity          ::Bool
        ... true, if the density above is taken as (total ion) number density ni, and false otherwise.
"""
struct Settings 
    temperature                 ::Float64     
    density                     ::Float64 
    useNumberDensity            ::Bool     
end 


"""
`Plasma.Settings()`  ... constructor for the default values of plasma computations
"""
function Settings()
    Settings(0., 0., false)
end


"""
`Plasma.Settings(set::Plasma.Settings;`

        temperature=..,         density=..,         useNumberDensity =..)
                    
    ... constructor for modifying the given Plasma.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::Plasma.Settings;    
    temperature::Union{Nothing,Float64}=nothing,                            density::Union{Nothing,Float64}=nothing,
    useNumberDensity::Union{Nothing,Bool}=nothing)  
    
    if  isnothing(temperature)        temperaturex      = set.temperature        else  temperaturex      = temperature       end 
    if  isnothing(density)            densityx          = set.density            else  densityx          = density           end 
    if  isnothing(useNumberDensity)   useNumberDensityx = set.useNumberDensity   else  useNumberDensityx = useNumberDensity  end 

    Settings( temperaturex, densityx, useNumberDensityx )
end


# `Base.show(io::IO, settings::Plasma.Settings)`  ... prepares a proper printout of the variable settings::Plasma.Settings.
function Base.show(io::IO, settings::Plasma.Settings) 
    println(io, "temperature:               $(settings.temperature)  ")
    println(io, "density:                   $(settings.density)  ")
    println(io, "useNumberDensity:          $(settings.useNumberDensity)  ")
end


"""
`struct  Computation`  
    ... defines a type for defining  (simple) plasma computation for atoms and ions in a given set of
        reference configurations. It also support different atomic processes under plasma conditions.
        The plasma environment is typical described in terms of its temperature, density, etc.

    + scheme                         ::AbstractPlasmaScheme            ... Scheme (kind) of plasma computation.
    + nuclearModel                   ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
    + grid                           ::Radial.Grid                     ... The radial grid to be used for the computation.
    + refConfigs                     ::Array{Configuration,1}          ... A list of non-relativistic configurations.
    + asfSettings                    ::AsfSettings                     
        ... Provides the settings for the SCF process (under plasma conditions) and the associated CI calculations.
    + settings                       ::Plasma.Settings                 ... communicates the properties of the plasma
"""
struct  Computation
    scheme                           ::AbstractPlasmaScheme 
    nuclearModel                     ::Nuclear.Model
    grid                             ::Radial.Grid
    refConfigs                       ::Array{Configuration,1}
    asfSettings                      ::AsfSettings                     
    settings                         ::Plasma.Settings
end 


"""
`Plasma.Computation()`  ... constructor for an 'empty' instance::Plasma.Computation.
"""
function Computation()
    Computation(AverageAtomScheme(), Nuclear.Model(1.), Radial.Grid(), Configuration[], AsfSettings(), Plasma.Settings() )
end


"""
`Plasma.Computation(comp::Plasma.Computation;`

    scheme=..,                  nuclearModel=..,            grid=..,                refConfigs=..,              asfSettings=..,     
    settings=..,
    printout::Bool=false)
                    
    ... constructor for modifying the given Plasma.Computation by 'overwriting' the previously selected parameters.
"""
function Computation(comp::Plasma.Computation;
    scheme::Union{Nothing,Plasma.AbstractPlasmaScheme}=nothing,                  
    nuclearModel::Union{Nothing,Nuclear.Model}=nothing,                         grid::Union{Nothing,Radial.Grid}=nothing,      
    refConfigs::Union{Nothing,Array{Configuration,1}}=nothing,                  asfSettings::Union{Nothing,AsfSettings}=nothing, 
    settings::Union{Nothing,Plasma.Settings}=nothing, 
    printout::Bool=false)
    
    if  isnothing(scheme)            schemex            = comp.scheme            else  schemex                  = scheme                   end 
    if  isnothing(nuclearModel)      nuclearModelx      = comp.nuclearModel      else  nuclearModelx            = nuclearModel             end 
    if  isnothing(grid)              gridx              = comp.grid              else  gridx                    = grid                     end 
    if  isnothing(refConfigs)        refConfigsx        = comp.refConfigs        else  refConfigsx              = refConfigs               end 
    if  isnothing(asfSettings)       asfSettingsx       = comp.asfSettings       else  asfSettingsx             = asfSettings              end 
    if  isnothing(settings)          settingsx          = comp.settings          else  settingsx                = settings                 end 
    
    
    cp = Computation(schemex, nuclearModelx, gridx, refConfigsx, asfSettingsx, settingsx) 
                        
    if printout  Base.show(cp)      end
    return( cp )
end


"""
`Plasma.Computation( wa::Bool)`  

        grid     = Radial.Grid(true)
        nuclearM = Nuclear.Model(18., "Fermi")
        ...
        refConfigs  = [Configuration("[Ne] 3s^2 3p^5")]
        Plasma.Computation(Plasma.Computation(), grid=grid, nuclearModel=nuclearM, refConfigs=refConfigs, asfSettings=... )
    
    ... These simple examples can be further improved by overwriting the corresponding parameters.
"""
function Computation(wa::Bool)    
    Plasma.Computation()    
end


# `Base.string(comp::Plasma.Computation)`  ... provides a String notation for the variable comp::Plasma.Computation.
function Base.string(comp::Plasma.Computation)
    sa = "Plasma computation:  for Z = $(comp.nuclearModel.Z), "
    return( sa )
end


# `Base.show(io::IO, comp::Plasma.Computation)`  ... prepares a printout of comp::Plasma.Computation.
function Base.show(io::IO, comp::Plasma.Computation)
    sa = Base.string(comp);             print(io, sa, "\n")
    println(io, "scheme:                $(comp.scheme)  ")
    println(io, "nuclearModel:          $(comp.nuclearModel)  ")
    println(io, "grid:                  $(comp.grid)  ")
    println(io, "refConfigs:            $(comp.refConfigs)  ")
    println(io, "asfSettings:           $(comp.asfSettings)  ")
    println(io, "settings:              $(comp.settings)  ")
end




"""
`struct  Plasma.PartialWaveData`  
    ... defines a type to collect for partial-waves (kappa) cross sections, rates, etc. at different energies

    + kappa    ::Int64                            ... kappa
    + pairs    ::Vector{Tuple{Float64, Float64}}  
        ... vector of pairs, for instance (energy, cs), to later combine data in a proper manner.
"""
struct PartialWaveData 
    kappa     ::Int64      
    pairs     ::Vector{Tuple{Float64, Float64}}  
end 


"""
`Plasma.PartialWaveData()`  ... constructor for the default values of Plasma.PartialWaveData set
"""
function PartialWaveData()
    PartialWaveData( -1, Tuple{Float64, Float64}[] )
end


# `Base.show(io::IO, data::Plasma.PartialWaveData)`  ... prepares a proper printout of the variable settings::Plasma.PartialWaveData.
function Base.show(io::IO, data::Plasma.PartialWaveData) 
    println(io, "kappa:          $(data.kappa)  ")
    println(io, "pairs:          $(data.pairs)  ")
end



include("module-Plasma-inc-average-atom.jl")
include("module-Plasma-inc-line-shifts.jl")
include("module-Plasma-inc-saha-boltzmann-mixture.jl")
include("module-Plasma-inc-satellite-diagnostic.jl")
include("module-Plasma-inc-collisional-radiative.jl")

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Basics.perform(comp::Plasma.Computation)`  
    ... to set-up and perform a plasma computation that starts from a given set of reference configurations and 
        support both, an atomic-average SCF procedure and the computation of various plasma properties and processe.
        The results of all individual steps are printed to screen but nothing is returned otherwise.

`Basics.perform(comp::Plasma.Computation; output::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the plasma computation but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Plasma.Computation; output::Bool=false)
    Plasma.perform(comp.scheme, comp::Plasma.Computation, output=output)
end



end # module
