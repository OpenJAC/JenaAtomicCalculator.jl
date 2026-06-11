

"""
`abstract type Basics.AbstractAngularMomentum` 
    ... defines an abstract and a number of concrete types for dealing with angular momentum variables.

    + AngularJ64        ... to deal with (total) angular momenta J >= 0
    + AngularM64        ... to deal with magnetic quantum numbers -J <= M <= J
    + LevelSymmetry     ... to deal with an overall J^P symmetry of a level.
"""
abstract type  AbstractAngularMomentum                         end

export  AbstractAngularMomentum, AngularJ64, AngularM64, LevelSymmetry, oplus


"""
`struct  AngularJ64 <: AbstractAngularMomentum`  ... defines a type for angular momenta j.

    + num  ::Int64              ... numerator
    + den  ::Int64              ... denominator, must be 1 or 2
"""
struct  AngularJ64 <: AbstractAngularMomentum
    num      ::Int64
    den      ::Int64
end


"""
`Basics.AngularJ64(m::Integer)`  ... constructor for a given integer (numerator).
"""
function AngularJ64(j::Integer)
    j < 0   &&   error("j must be positive.")
    AngularJ64(j, 1)
end


"""
`Basics.AngularJ64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
"""
function AngularJ64(rational::Rational{Int64})
    !(rational.den in [1,2])   &&   error("Denominator must be 1 or 2.")
    rational.num < 0           &&   error("j must be positive.")
    AngularJ64(rational.num, rational.den)
end


# `Base.show(io::IO, j::AngularJ64)`  ... prepares a proper printout of the variable j::AngularJ64.
function Base.show(io::IO, j::AngularJ64) 
    if      j.den == 1    print(io, j.num)
    elseif  j.den == 2    print(io, j.num, "/2")
    else    error("stop a")
    end
end


"""
`Basics.oplus(ja::AngularJ64, jb::AngularJ64)`  
    ... adds the angular momenta ja `oplus` jb and returns a list::Array{AngularJ64,1} of j-valuescin the interval 
        |ja - jb| <= j <= ja + jb.
"""
function  oplus(ja::AngularJ64, jb::AngularJ64)
    if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
    if  jb.den == 1   jb2 = 2jb.num   else   jb2 = jb.num   end
    jList = AngularJ64[]
    for  j = abs(ja2 - jb2):2:ja2+jb2    push!(jList, AngularJ64(j//2) )    end
    return( jList )
end

oplus(ja::AngularJ64, jb::Int64) = oplus(ja, AngularJ64(jb))
oplus(ja::Int64, jb::AngularJ64) = oplus(AngularJ64(ja), jb)
oplus(ja::Int64, jb::Int64)      = oplus(AngularJ64(ja), AngularJ64(jb))

Base.Float64(ja::AngularJ64) = ja.num / ja.den


"""
`Basics.projections(ja::AngularJ64)`  
    ... returns all allowed projections of a given angular momenta ja as a list::Array{AngularM64,1} of m-values, i.e. -ja, -ja+1, ..., j.
"""
function  projections(ja::AngularJ64)
    if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
    mList = AngularM64[]
    for  m = -ja2:2:ja2    push!(mList, AngularM64(m//2) )    end
    return( mList )
end


"""
`struct  Basics.AngularM64 <: AbstractAngularMomentum`  ... defines a type for angular momentum projections m.

    + num  ::Int64              ... numerator
    + den  ::Int64              ... denominator, must be 1 or 2
"""
struct  AngularM64 <: AbstractAngularMomentum
    num      ::Int64
    den      ::Int64
end


"""
`Basics.AngularM64(m::Integer)`  ... constructor for a given integer (numerator).
"""
function AngularM64(m::Integer)
    AngularM64(m, 1)
end


"""
`Basics.AngularM64(m::Integer, j::AngularJ64)`  
    ... constructor for a given integer (numerator) that is consistent with a given j-value.
"""
function AngularM64(m::Integer, j::AngularJ64)
    !(j.den == 1)      &&  error("m must be integer for j = $(j).")
    j.num < abs(m)     &&  error("abs(m) must be <= j = $(j).")
    AngularM64(m, 1)
end


"""
`Basics.AngularM64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
"""
function AngularM64(rational::Rational{Int64})
    !(rational.den in [1,2])  &&  error("Denominator must be 1 or 2.")
    AngularM64(rational.num, rational.den)
end


"""
`Basics.AngularM64(rational::Rational{Int64}, j::AngularJ64)`  
    ... constructor for a given  Rational{Int64} that is consistent with a given m-value.
"""
function AngularM64(rational::Rational{Int64}, j::AngularJ64)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    !(j.den == rational.den)      &&   error("j,m must be both integer or half-integer.")
    j.num < abs(rational.num)     &&   error("abs(m) must be <= j = $(j).")
    AngularM64(rational.num, rational.den)
end


"""
`Basics.AngularM64(j::AngularJ64)`  ... constructor for a given j::AngularJ64 to define m = j.
"""
function AngularM64(j::AngularJ64)
    AngularM64(j.num, j.den)
end


# `Base.show(io::IO, m::AngularM64)`  ... prepares a proper printout of the variable m::AngularM64.
function Base.show(io::IO, m::AngularM64) 
    if      m.den == 1    print(io, m.num)
    elseif  m.den == 2    print(io, m.num, "/2")
    else    error("stop a")
    end
end


"""
`Basics.add(ma::AngularM64, mb::AngularM64)`  
    ... adds the projections of the angular momenta ma + mb and returns a mc::AngularM64.
"""
function  add(ma::AngularM64, mb::AngularM64)
    if  ma.den == 1   ma2 = 2ma.num   else   ma2 = ma.num   end
    if  mb.den == 1   mb2 = 2mb.num   else   mb2 = mb.num   end
    return( AngularM64( (ma2+mb2)//2 ) )
end


# Conversion between HalfInt and AngularJ64, AngularM64

twice(x::Union{AngularJ64,AngularM64}) = ifelse(isone(x.den), twice(x.num), x.num)
twice(x) = x + x


"""
`struct  Basics.LevelSymmetry  <  AbstractAngularMomentum`  ... defines a struct for defining the overall J^P symmetry of a level.

    + J          ::AngularJ64  ... total angular momentum of a level
    + parity     ::Parity      ... total parity of the level
"""
struct  LevelSymmetry <: AbstractAngularMomentum
    J            ::AngularJ64
    parity       ::Parity  
end


"""
`Basics.LevelSymmetry(rational::Rational{Int64}, parity::Parity)`  ... constructor for a given (Rational,Parity).
"""
function  LevelSymmetry(rational::Rational{Int64}, parity::Parity)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    LevelSymmetry( AngularJ64(rational.num, rational.den), parity )    
end


"""
`Basics.LevelSymmetry(i::Int64, parity::Parity)`  ... constructor for a given (Int64,Parity).
"""
function  LevelSymmetry(i::Int64, parity::Parity)
    LevelSymmetry( AngularJ64(i), parity )    
end


"""
`Basics.LevelSymmetry(rational::Rational{Int64}, sa::String)`  ... constructor for a given (Rational,String).
"""
function  LevelSymmetry(rational::Rational{Int64},sa::String)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    LevelSymmetry( AngularJ64(rational.num, rational.den), Parity(sa) )    
end


"""
`Basics.LevelSymmetry(i::Int64, sa::String)`  ... constructor for a given (Int64,String).
"""
function  LevelSymmetry(i::Int64,sa::String)
    LevelSymmetry( AngularJ64(i), Parity(sa) )    
end


# `Base.show(io::IO, sym::LevelSymmetry)`  ... prepares a proper printout of the variable sym::LevelSymmetry.
function Base.show(io::IO, sym::LevelSymmetry) 
    print(io, string(sym) )
end


# `Base.string(sym::LevelSymmetry)`  ... provides a proper printout of the variable sym::LevelSymmetry.
function Base.string(sym::LevelSymmetry) 
    if      sym.parity == plus     return( "$(sym.J) +" )
    elseif  sym.parity == minus    return( "$(sym.J) -" )  
    else    error("stop a")
    end
end

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractContinuumNormalization` 
    ... defines an abstract and a number of singleton types for dealing with the normalization of continuum orbitals.

    + PureSineNorm       ... normalize with regard to an (asymtotic) pure sine funtion, sin(kr).
    + CoulombSineNorm    ... normalize with regard to an (asymtotic) Coulombic-sine funtion, sin(kr + ...).
    + OngRussekNorm      ... normalize by following Ong & Russek (1973).
    + AlokNorm           ... normalize following Salvats Radial code
"""
abstract type  ContinuumNormalization                          end
struct     PureSineNorm         <:  ContinuumNormalization     end
struct     CoulombSineNorm      <:  ContinuumNormalization     end
struct     OngRussekNorm        <:  ContinuumNormalization     end
struct     AlokNorm             <:  ContinuumNormalization     end

export  AbstractContinuumNormalization,   PureSineNorm,   CoulombSineNorm,   OngRussekNorm,   AlokNorm

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractConfigurationRestriction` 
    ... defines an abstract types for dealing with restrictions that need to be applied to a list of configurations.
        Typically, a loop through is made through all given restrictions and all configurations are tested to obey all
        these restrictions. Two contradicting restrictions, for instance RestrictParity(plus) & RestrictParity(minus),
        therefore leads zero configurations in all cases. It remains the reponsibility of the user to make sure that
        the given restrictions are consistent with what is to be achieved. The given set of restrictions can be easily
        extended if this need arises by the users.

    + RestrictMaximumDisplacements(..)  ... to restrict the maximum replacements wrt a (2nd) configuration.
    + RestrictNoElectronsTo(..)         ... to restrict the total number of electron in high subshells.
    + RestrictParity(..)                ... to restrict to configurations with a given parity.
    + RestrictToShellDoubles(..)        ... to allow only double occupations in high subshells.
    + RequestMinimumOccupation(..)      ... to request a minimum occupation in a given set of shells.
    + RequestMaximumOccupation(..)      ... to request a maximum occupation in a given set of shells.
"""
abstract type  AbstractConfigurationRestriction      end

export  AbstractConfigurationRestriction, RestrictMaximumDisplacements, RestrictNoElectronsTo, RestrictParity, 
        RestrictToShellDoubles, RequestMinimumOccupation, RequestMaximumOccupation

"""
`struct  Basics.RestrictMaximumDisplacements  <: AbstractConfigurationRestriction`   
    ... restrict the maximum replacements w.r.t. (another) given configuration, which can have a different number of electrons.
        A "displacement" is simply defined as the difference of occuation numbers. This restriction is useful to determine allowed
        configurations in a second-order treatment of atomic processes. An odd number of displacement naturally arise for all 
        configurations which differ by one or three electrons from the reference configurations. Several restrictions of this
        type can be formulated but are treated separately.

    + conf          ::Configuration   ... configuration w.r.t. which displacements are taken.
    + maxDisplace   ::Int64           ... maximum number of displacements >= 0.
"""
struct   RestrictMaximumDisplacements  <: AbstractConfigurationRestriction
    conf            ::Configuration
    maxDisplace     ::Int64
end


function Base.string(res::RestrictMaximumDisplacements)
    sa = "Restrict to configurations with maximum $(res.maxDisplace) displacements w.r.t. $(res.conf)."
    return( sa )
end

function Base.show(io::IO, res::RestrictMaximumDisplacements)
    sa = string(res);       print(io, sa)
end


"""
`struct  Basics.RestrictNoElectronsTo  <: AbstractConfigurationRestriction`   
    ... restrict the number of electron in all shells with principal quantum number n >= nmin or orbital angular momentum l >= lmin 
        to a total of ne electrons.

    + ne            ::Int64     ... maximum number of (allowed) electrons in the specicied higher subshells.
    + nmin          ::Int64     ... principal quantum number nmin.
    + lmin          ::Int64     ... orbital angular momentum lmin.
"""
struct   RestrictNoElectronsTo  <: AbstractConfigurationRestriction
    ne              ::Int64
    nmin            ::Int64
    lmin            ::Int64
end


function Base.string(res::RestrictNoElectronsTo)
    sa = "Restrict to configurations with a maximum of $(res.ne) electrons in shells with n >= $(res.nmin) & l >= $(res.lmin)."
    return( sa )
end

function Base.show(io::IO, res::RestrictNoElectronsTo)
    sa = string(res);       print(io, sa)
end


"""
`struct  Basics.RestrictParity  <: AbstractConfigurationRestriction`   
    ... restrict to configurations with a given parity.

    + parity        ::Basics.Parity   ... given parity.
"""
struct   RestrictParity  <: AbstractConfigurationRestriction
    parity          ::Basics.Parity
end


function Base.string(res::RestrictParity)
    sa = "Restrict to configurations with parity $(res.parity)."
    return( sa )
end

function Base.show(io::IO, res::RestrictParity)
    sa = string(res);       print(io, sa)
end


"""
`struct  Basics.RestrictToShellDoubles  <: AbstractConfigurationRestriction`   
    ... restrict to a double electron occupation in all shells with principal quantum number n >= nmin or orbital angular momentum l >= lmin. 

    + nmin          ::Int64     ... principal quantum number nmin.
    + lmin          ::Int64     ... orbital angular momentum lmin.
"""
struct   RestrictToShellDoubles  <: AbstractConfigurationRestriction
    nmin            ::Int64
    lmin            ::Int64
end


function Base.string(res::RestrictToShellDoubles)
    sa = "Restrict to configurations with a double electron occupation in shells with n >= $(res.nmin) & l >= $(res.lmin)."
    return( sa )
end

function Base.show(io::IO, res::RestrictToShellDoubles)
    sa = string(res);       print(io, sa)
end


"""
`struct  Basics.RequestMinimumOccupation  <: AbstractConfigurationRestriction`   
    ... request a minimum occupation ne in the given (list of) shells. 

    + ne            ::Int64           ... minimum electron occupation.
    + shells        ::Array{Shell,1}  ... list of shells.
"""
struct   RequestMinimumOccupation  <: AbstractConfigurationRestriction
    ne              ::Int64 
    shells          ::Array{Shell,1}
end


function Base.string(res::RequestMinimumOccupation)
    sa = "Request a minimum occupation of $(res.ne) electron in the (list of) shells $(res.shells)."
    return( sa )
end

function Base.show(io::IO, res::RequestMinimumOccupation)
    sa = string(res);       print(io, sa)
end


"""
`struct  Basics.RequestMaximumOccupation  <: AbstractConfigurationRestriction`   
    ... request a maximum occupation ne in the given (list of) shells. 

    + ne            ::Int64           ... maximum electron occupation.
    + shells        ::Array{Shell,1}  ... list of shells.
"""
struct   RequestMaximumOccupation  <: AbstractConfigurationRestriction
    ne              ::Int64 
    shells          ::Array{Shell,1}
end


function Base.string(res::RequestMaximumOccupation)
    sa = "Request a maximum occupation of $(res.ne) electron in the (list of) shells $(res.shells)."
    return( sa )
end

function Base.show(io::IO, res::RequestMaximumOccupation)
    sa = string(res);       print(io, sa)
end


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractConfigurationTheme` 
    ... defines an abstract and a number of concrete (detailed and singleton) data types to distinguish a good number of concrete
        themes for manipulating (list of) configurations. These themes are listed and briefly explained below; they typically provide
        the data central to the particular theme, while other input for applying the theme is handled by multiple dispatch.

    + AddElectrons            ... to add to the given configurations one or several electrons into specified shells.
    + ExciteElectrons         ... to excite for the given configurations one or several electrons into specified shells.
    + RemoveElectrons         ... to remove from the given configurations one or several electrons from specified shells.
    + RestrictExcitations     ... to restrict the given configurations due to one or several configuration restrictions.
    
    + ForAutoIonization       ... to generate configurations that are related by autoionization (deexcitation + single remove).
    + ForDielectronicCapture  ... to generate configurations that are related by dielectronic capture (excitation + single capture).
    + ForDielectronicRecombination     ... to compute the DR resonance strength (for pedestrians only).
    + ForElectronCapture      ... to generate configurations that are related by the capture of an additional electron.
    + ForImpactIonization     ... to estimate electron impact-ionization cross sections (for pedestrians only)
    + ForHollowIons           ... to generate configurations that are related by multiple capture into high-n shells.
    + ForPhotoEmission        ... to generate configurations that are related by photoemission (single replacement of electrons).
    + ForPhotoIonization      ... to generate configurations that are related by photoionization (single removement of electrons).
    + ForPhotoRecombination   ... to generate configurations that are related by radiative capture (just single-electron capture).
    + ForStepwiseDecay        ... to generate configurations that are related by photoemission and autoionization.

    + ForGivenConfigs         ... to perform computations for given configurations.
    
    + GroundConfiguration     ... to generate the ground configuration for a given number of electrons.
    + MeanConfiguration       ... to generate the mean configuration, i.e. a configuration with mean occupation numbers.
    + RelativisticConfigurations  to generate/deal with relativistic configurations.
    + SuperConfiguration      ... to generate configurations from a given super-configuration.

    + AllShells               ... to extract all occupied shells from a set of configurations.
    + ByMultipoles            ... to extract the configurations due to multipole selection themes.
    + ByParity                ... to extract the configurations due to parity.
    + ClosedCore              ... to extract the closed core from the given configuration.
    + ClosedShells            ... to extract the closed shells from given configurations.
    + ClosedSubshells         ... to extract the closed subshells from given configurations.
    + ContractShells          ... to extract the configurations without empty shell occupation (contracted shells).
    + ExcitationLevel         ... to determine the excitation level of a configuration.
    + ExpandShells            ... to extract the configurations with empty shell occupation (expanded shells).
    + FromBasis               ... to extract the configuration from the basis.
    + FromMultiplet           ... to extract the configuration from the multiplet(s).
    + GeneralizedConfiguration .. to extract the generalized configuration.
    + GetParity               ... to extract the parity of a configurations.
    + IsOccupied              ... to extract where a shell or subshell is occupied in the given configuration.
    + LeadingConfiguration    ... to extract the leading configuration.
    + LeadingConfigurationR   ... to extract the leading relativistic configuration.
    + MeanOccupation          ... to extract occupation numbers from given configurations.
    + Multiplicity            ... to extract the multiplicity of a configuration.
    + NonrelativisticBasis    ... to extract the configurations from a non-relativistic basis.
    + NumberOfElectrons       ... to extract the number of electrons from a set of configurations.
    + OccupationDifference    ... to extract differences of occupation numbers from given configurations.
    + OpenShells              ... to extract the open shells from the given configurations.
    + OpenSubshells           ... to extract the open subshells from the given configurations.
    + TotalAM                 ... to extract the total angular momenta J that are associated with the configuration.
    + ValenceOccupation       ... to extract the remaining configuration beyond a given (closed) core configuration.
    
    + FineStructure           ... to display the total J fine-structure levels of a configuration (without energies).
    + FineStructureLS         ... to display the total LSJ fine-structure levels of a configuration (without energies).
    + HyperfineStructure      ... to display the total F hyperfine-structure levels of a configuration (without energies).
    + HundRules               ... to display the total LSJ fine-structure levels, ordered by Hund's themes (not yet).
        
"""
abstract type  AbstractConfigurationTheme                                  end
#
#
struct   ForAutoIonization              <:  AbstractConfigurationTheme     end
struct   ForImpactIonization            <:  AbstractConfigurationTheme     end
struct   ForPhotoEmission               <:  AbstractConfigurationTheme     end
struct   ForPhotoIonization             <:  AbstractConfigurationTheme     end

struct   ForGivenConfigs                <:  AbstractConfigurationTheme     end
#
struct   MeanConfiguration              <:  AbstractConfigurationTheme     end
struct   RelativisticConfigurations     <:  AbstractConfigurationTheme     end
struct   SplitByEnerby                  <:  AbstractConfigurationTheme     end
struct   SuperConfiguration             <:  AbstractConfigurationTheme     end
#
struct   AllShells                      <:  AbstractConfigurationTheme     end
struct   ByMultipoles                   <:  AbstractConfigurationTheme     end
struct   ClosedCore                     <:  AbstractConfigurationTheme     end
struct   ClosedShells                   <:  AbstractConfigurationTheme     end
struct   ClosedSubshells                <:  AbstractConfigurationTheme     end
struct   ContractShells                 <:  AbstractConfigurationTheme     end
struct   ExcitationLevel                <:  AbstractConfigurationTheme     end
struct   FromBasis                      <:  AbstractConfigurationTheme     end
struct   FromMultiplet                  <:  AbstractConfigurationTheme     end
struct   GeneralizedConfigurations      <:  AbstractConfigurationTheme     end
struct   GetParity                      <:  AbstractConfigurationTheme     end 
struct   IsOccupied                     <:  AbstractConfigurationTheme     end 
struct   LeadingConfiguration           <:  AbstractConfigurationTheme     end
struct   LeadingConfigurationR          <:  AbstractConfigurationTheme     end
struct   MeanOccupation                 <:  AbstractConfigurationTheme     end
struct   Multiplicity                   <:  AbstractConfigurationTheme     end
struct   NonrelativisticBasis           <:  AbstractConfigurationTheme     end
struct   NumberOfElectrons              <:  AbstractConfigurationTheme     end
struct   OccupationDifference           <:  AbstractConfigurationTheme     end
struct   OpenShellNumber                <:  AbstractConfigurationTheme     end
struct   OpenShells                     <:  AbstractConfigurationTheme     end
struct   OpenSubshells                  <:  AbstractConfigurationTheme     end
struct   ValenceOccupation              <:  AbstractConfigurationTheme     end
struct   ValenceShells                  <:  AbstractConfigurationTheme     end
#
struct   FineStructure                  <:  AbstractConfigurationTheme     end
struct   FineStructureLS                <:  AbstractConfigurationTheme     end
struct   HundsRules                     <:  AbstractConfigurationTheme     end

export  AbstractConfigurationTheme, AddElectrons, ExciteElectrons, RemoveElectrons, RestrictExcitations,
        ForAutoIonization, ForElectronCapture, ForDielectronicCapture, ForDielectronicRecombination, ForHollowIons, 
        ForPhotoEmission, ForPhotoIonization,  ForPhotoRecombination, ForRasExcitations, ForStepwiseDecay, 
        GeneralizedConfigurations, GroundConfiguration, MeanConfiguration, RelativisticConfigurations, 
        SuperConfiguration,
        AllShells, ByNumber, ByParity, ClosedCore, ClosedShells, ClosedSubshells, ContractShells, ExcitationLevel, 
        ExpandShells, FromBasis, FromMultiplet, GetParity, IsOccupied, LeadingConfiguration, LeadingConfigurationR, 
        MeanOccupation, Multiplicity, NonrelativisticBasis, NumberOfElectrons, OccupationDifference, OpenShellNumber, 
        OpenShells, OpenSubshells, TotalAM, ValenceOccupation, ValenceShells,
        FineStructure, FineStructureLS, HundsRules, HyperfineStructure

        
"""
`struct  Basics.AddElectrons          <:  AbstractConfigurationTheme`   
    ... to add ne electrons in the intoshells for each of the given configurations. 

    + ne             ::Int64            ... number of electrons to be added.
    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are added.
"""
struct   AddElectrons                 <:  AbstractConfigurationTheme
    ne               ::Int64           
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.AddElectrons)
    sa = "AddElectrons theme with number of electrons ne=$(theme.ne) added into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.AddElectrons)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ByNumber              <:  AbstractConfigurationTheme`   
    ... to select configurations due to given numbers of electrons. 

    + NoElectrons    ::Array{Int64,1}   ... List of selected electron numbers
"""
struct   ByNumber                     <:  AbstractConfigurationTheme
    NoElectrons      ::Array{Int64,1}           
end


function Base.string(theme::Basics.ByNumber)
    sa = "ByNumber theme with selected numbers of electrons $(theme.NoElectrons)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ByNumber)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ByParity              <:  AbstractConfigurationTheme`   
    ... select configurations due to a given parity. 

    + P              ::Basics.Parity  ... Selected parity.
"""
struct   ByParity                     <:  AbstractConfigurationTheme
    P                ::Basics.Parity           
end


function Base.string(theme::Basics.ByParity)
    sa = "ByParity theme with selected parity P=$(theme.P)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ByParity)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ExciteElectrons       <:  AbstractConfigurationTheme`   
    ... to excite ne electrons fromShells to the intoshells for each of the given configurations. 

    + ne             ::Int64            ... number of electrons to be added.
    + fromShells     ::Array{Shell,1}   ... The shells from which the electrons are excited.
    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are excited.
"""
struct   ExciteElectrons              <:  AbstractConfigurationTheme
    ne               ::Int64           
    fromShells       ::Array{Shell,1}
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.ExciteElectrons)
    sa = "ExciteElectrons theme with number of electrons ne=$(theme.ne) from $(theme.fromShells) into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ExciteElectrons)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ExpandShells          <:  AbstractConfigurationTheme`   
    ... expand a configuration with empty shells due to a given number of shells

    + shells         ::Array{Shell,1}  ... Shells that will occur in the expanded form.
"""
struct   ExpandShells                 <:  AbstractConfigurationTheme
    shells           ::Array{Shell,1}
end


function Base.string(theme::Basics.ExpandShells)
    sa = "ExpandShells theme with shells=$(theme.shells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ExpandShells)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ForDielectronicCapture  <:  AbstractConfigurationTheme`   
    ... to excite one electron fromShells to the toshells and add another electron into the intoShells
        for each of the given configurations. 

    + fromShells     ::Array{Shell,1}   ... The shells from which the electrons are excited from --> to.
    + toShells       ::Array{Shell,1}   ... The shells to which the electrons are excited from --> to.
    + intoShells     ::Array{Shell,1}   ... The shells into which the additional electron is captured.
"""
struct   ForDielectronicCapture         <:  AbstractConfigurationTheme
    fromShells       ::Array{Shell,1}
    toShells         ::Array{Shell,1}
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.ForDielectronicCapture)
    sa = "ForDielectronicCapture theme with electron excitaton from $(theme.fromShells) to $(theme.toShells) " *
         "as well as the capture of an additional electron into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForDielectronicCapture)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ForDielectronicRecombination     <:  AbstractConfigurationTheme`   
    ... to excite one electron fromShells to the toshells and add another electron into the intoShells
        for each of the given configurations. The subsequent stabiliztion is considered into the 
            decayShells.

    + fromShells     ::Array{Shell,1}   ... The shells from which the electrons are excited from --> to.
    + toShells       ::Array{Shell,1}   ... The shells to which the electrons are excited from --> to.
    + intoShells     ::Array{Shell,1}   ... The shells into which the additional electron is captured.
    + decayShells    ::Array{Shell,1}   ... The shells into which (radiative) stabilization occurs.
"""
struct   ForDielectronicRecombination           <:  AbstractConfigurationTheme
    fromShells       ::Array{Shell,1}
    toShells         ::Array{Shell,1}
    intoShells       ::Array{Shell,1}
    decayShells      ::Array{Shell,1}
end


function Base.string(theme::Basics.ForDielectronicRecombination)
    sa = "ForDielectronicRecombination theme with electron excitaton from $(theme.fromShells) to $(theme.toShells), " *
         "together with the capture of an additional electron into $(theme.intoShells) " *
        "\nand the (raditive) stabilization in the $(theme.decayShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForDielectronicRecombination)
    sa = string(theme);       print(io, sa)
end  

  
"""
`struct  Basics.ForElectronCapture      <:  AbstractConfigurationTheme`   
    ... to add ne electrons in the intoshells for each of the given configurations. 

    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are added.
"""
struct   ForElectronCapture             <:  AbstractConfigurationTheme
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.ForElectronCapture)
    sa = "ForElectronCapture theme with one additional electron into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForElectronCapture)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ForHollowIons           <:  AbstractConfigurationTheme`   
    ... to excite ne electrons fromShells to the intoshells for each of the given configurations. 

    + ne             ::Int64            ... number of electrons to be captured.
    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are captured
    + decayShells    ::Array{Shell,1}   ... The shells which need to be considered between the shells of the ions and 
                                            the intoShells in order to model the decay of hollow ions.
"""
struct   ForHollowIons                  <:  AbstractConfigurationTheme
    ne               ::Int64           
    intoShells       ::Array{Shell,1}
    decayShells      ::Array{Shell,1}
end


function Base.string(theme::Basics.ForHollowIons)
    sa = "ForHollowIons theme with number of electrons ne=$(theme.ne), the capture into into $(theme.intoShells) and " *
         "the subsequent decay into $(theme.decayShells)"
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForHollowIons)
    sa = string(theme);       print(io, sa)
end
  
  
"""
`struct  Basics.ForPhotoRecombination   <:  AbstractConfigurationTheme`   
    ... to add ne electrons in the intoshells for each of the given configurations. 

    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are added.
"""
struct   ForPhotoRecombination          <:  AbstractConfigurationTheme
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.ForPhotoRecombination)
    sa = "ForPhotoRecombination theme with one additional electron into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForPhotoRecombination)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.ForRasExcitations       <:  AbstractConfigurationTheme`   
    ... to excite ne electrons fromShells to the intoshells for each of the given configurations. 

    + se             ::Bool             ... True if single excitations to be included, and false otherwise.
    + de             ::Bool             ... True if double excitations to be included, and false otherwise.
    + te             ::Bool             ... True if triple excitations to be included, and false otherwise.
    + qe             ::Bool             ... True if quadruple excitations to be included, and false otherwise.
    + fromShells     ::Array{Shell,1}   ... The shells from which the electrons are excited.
    + intoShells     ::Array{Shell,1}   ... The shells into which the electrons are excited.
"""
struct   ForRasExcitations               <:  AbstractConfigurationTheme
    se               ::Bool   
    de               ::Bool   
    te               ::Bool   
    qe               ::Bool   
    fromShells       ::Array{Shell,1}
    intoShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.ForRasExcitations)
    sa = "ForRasExcitations theme with SDTQ = ($(theme.se), $(theme.de), $(theme.te), $(theme.qe)) and excitations " *
         "from $(theme.fromShells) into $(theme.intoShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForRasExcitations)
    sa = string(theme);       print(io, sa)
end
  
  
"""
`struct  Basics.ForStepwiseDecay        <:  AbstractConfigurationTheme`   
    ... to generate all those configurations that occur due to the stepwise photoemission 
        autoionization of configurations with some inner-shell hole.

    + maximallyReleased   ::Int64   ... Maximum number of electrons that can be released from the 
                                        given configurations.
"""
struct   ForStepwiseDecay               <:  AbstractConfigurationTheme
    maximallyReleased     ::Int64
end


function Base.string(theme::Basics.ForStepwiseDecay)
    sa = "ForStepwiseDecay theme with the maximum number of released electrons $(theme.maximallyReleased)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.ForStepwiseDecay)
    sa = string(theme);       print(io, sa)
end
  
  
"""
`struct  Basics.HyperfineStructure      <:  AbstractConfigurationTheme`   
    ... to display the total F hyperfine-structure levels of a configuration (without energies).

    + spinI        ::AngularJ64   ... Nuclear spin I.
"""
struct   HyperfineStructure             <:  AbstractConfigurationTheme
    spinI          ::AngularJ64 
end


function Base.string(theme::Basics.HyperfineStructure)
    sa = "HyperfineStructure theme with nuclear spin $(theme.spinI)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.HyperfineStructure)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.RemoveElectrons       <:  AbstractConfigurationTheme`   
    ... to remove ne electrons fromShells for each of the given configurations. 

    + ne             ::Int64            ... number of electrons to be added.
    + fromShells     ::Array{Shell,1}   ... The shells from which the electrons are removed.
"""
struct   RemoveElectrons              <:  AbstractConfigurationTheme
    ne               ::Int64           
    fromShells       ::Array{Shell,1}
end


function Base.string(theme::Basics.RemoveElectrons)
    sa = "RemoveElectrons theme with number of electrons ne=$(theme.ne) from $(theme.fromShells)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.RemoveElectrons)
    sa = string(theme);       print(io, sa)
end

        
"""
`struct  Basics.RestrictExcitations   <:  AbstractConfigurationTheme`   
    ... to restrict (reduce the number of) the given configurations by applying one or several configuration
        restrictions. 

    + ne             ::Int64                                       ... number of electrons to be added.
    + fromShells     ::Array{Shell,1}                              ... The shells from which the electrons are excited.
    + toShells       ::Array{Shell,1}                              ... The shells into which the electrons are excited.
    + restrictions   ::Array{AbstractConfigurationRestriction,1}   ... set of restrictions that will be applied.
"""
struct   RestrictExcitations          <:  AbstractConfigurationTheme
    ne               ::Int64        
    fromShells       ::Array{Shell,1} 
    toShells         ::Array{Shell,1} 
    restrictions     ::Array{AbstractConfigurationRestriction,1}
end


"""
`Basics.RestrictExcitations(restrictions::Array{AbstractConfigurationRestriction,1})`  
    ... constructor to just provide a set of restriction.
"""
function Basics.RestrictExcitations(restrictions::Array{AbstractConfigurationRestriction,1})
    Basics.RestrictExcitations(0, Shell[], Shell[], restrictions)
end


"""
`Basics.RestrictExcitations(theme::Basics.ExciteElectrons, restrictions::Array{AbstractConfigurationRestriction,1})`  
    ... constructor to specify the paraemterss by the ExciteElectron() theme.
"""
function Basics.RestrictExcitations(theme::Basics.ExciteElectrons, restrictions::Array{AbstractConfigurationRestriction,1})
    Basics.RestrictExcitations(theme.ne, theme.fromShells, theme.intoShells, restrictions)
end


function Base.string(theme::Basics.RestrictExcitations)
    sa = "RestrictExcitations theme with the excitation of $(theme.ne) electrons from $(theme.fromShells) " *
         "to $(theme.toShells) and with the following $(length(theme.restrictions)) restriction:"
    return( sa )
end

function Base.show(io::IO, theme::Basics.RestrictExcitations)
    sa = string(theme);       print(io, sa)
    for  restriction  in  theme.restrictions
        print(io, restriction)
    end
end


"""
`struct  Basics.GroundConfiguration   <:  AbstractConfigurationTheme`   
    ... to generate the ground configuration for a given (Z,N), i.e. the (nuclearCharge, NoElectrons). 

    + Z             ::Float64         ... nuclear charge Z
    + NoElectrons   ::Int64           ... number of electrons.
"""
struct   GroundConfiguration          <:  AbstractConfigurationTheme
    Z               ::Float64
    NoElectrons     ::Int64
end


function Base.string(theme::Basics.GroundConfiguration)
    sa = "GroundConfiguration theme with charge Z=$(theme.Z) and NoElectrons=$(theme.NoElectrons)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.GroundConfiguration)
    sa = string(theme);       print(io, sa)
end


"""
`struct  Basics.TotalAM      <:  AbstractConfigurationTheme`   
    ... to extract the total angular momenta (AM) to which the CSF of a configuration can couple. 

    + allJ          ::Bool                  ... True, if all J-values (including multiple couplings) are to be returned, 
                                                and false otherwise.
    + totalJs       ::Array{AngularJ64,1}   ... Selected total angular momenta J.
"""
struct   TotalAM             <:  AbstractConfigurationTheme
    allJ            ::Bool
    totalJs         ::Array{AngularJ64,1}
end


function Base.string(theme::Basics.TotalAM)
    sa = "Total angular momenta with allJ=$(theme.allJ) and totalJs=$(theme.totalJs)."
    return( sa )
end

function Base.show(io::IO, theme::Basics.TotalAM)
    sa = string(theme);       print(io, sa)
end


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractContinuumSolutions` 
    ... defines an abstract and a number of singleton types for solving the continuum orbitals in a given potential.

    + ContBessel              ... generate a pure Bessel function for the large component together with kinetic balance.
    + ContSine                ... generate a pure Sine function for the large component together with kinetic balance.
    + NonrelativisticCoulomb  ... generate a non-relativistic Coulomb function for the large component together with kinetic balance.
    + AsymptoticCoulomb       ... generate a pure (asymptotic) Coulombic function for both components.
    + BsplineGalerkin         ... generate a continuum orbital with the Galerkin method.dealing with warnings that are made during a run or REPL session.
"""
abstract type  AbstractContinuumSolutions                            end
struct     ContBessel             <:  AbstractContinuumSolutions     end
struct     ContSine               <:  AbstractContinuumSolutions     end
struct     NonrelativisticCoulomb <:  AbstractContinuumSolutions     end
struct     AsymptoticCoulomb      <:  AbstractContinuumSolutions     end
struct     BsplineGalerkin        <:  AbstractContinuumSolutions     end

export  AbstractContinuumSolutions, ContBessel, ContSine, NonrelativisticCoulomb, AsymptoticCoulomb, BsplineGalerkin 

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractEeInteraction` 
    ... defines an abstract and a number of singleton types for specifying the electron-electron interaction.

    + struct DiagonalCoulomb                   ... to represent the Coulomb part of the e-e interaction for just diagonal ME.        
    + struct CoulombInteraction                ... to represent the Coulomb part of the electron-electron interaction.        
    + struct CoulombGaunt                      ... to represent the Coulomb part of the electron-electron interaction.        
    + struct BreitInteraction(factor::Float64) 
        ... to represent the (frequency-dependent) Breit part of the electron-electron interaction with factor*omega,
            including the zero-Breit approximation (if factor=0.).
    + struct CoulombBreit(factor::Float64)     ... to represent the Coulomb+Breit part of the electron-electron interaction.        
"""
abstract type  AbstractEeInteraction                          end
struct     DiagonalCoulomb      <:  AbstractEeInteraction     end
struct     CoulombInteraction   <:  AbstractEeInteraction     end
struct     CoulombGaunt         <:  AbstractEeInteraction     end
struct     BreitInteraction     <:  AbstractEeInteraction     
    factor ::Float64
end
struct     CoulombBreit         <:  AbstractEeInteraction     
    factor ::Float64
end

export  AbstractEeInteraction, DiagonalCoulomb, CoulombInteraction, CoulombGaunt, BreitInteraction, CoulombBreit

function Base.show(io::IO, eeint::Union{BreitInteraction,CoulombBreit}) 
    sa = "$(typeof(eeint)) [factor=$(eeint.factor)]";                print(io, sa)
end

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractEmissionKind`
    ... defines an abstract and two singleton types to distinguish between the emission and
        absorption direction of a radiative amplitude.

    + struct Emission    ... to compute the photon-emission amplitude  <f || O^(Mp) || i>.
    + struct Absorption  ... to compute the photon-absorption amplitude <f || O^(Mp) || i>;
                            equal to the conjugate of the emission amplitude with initial
                            and final levels interchanged.
"""
abstract type  AbstractEmissionKind                      end
struct         Emission    <:  AbstractEmissionKind      end
struct         Absorption  <:  AbstractEmissionKind      end

export  AbstractEmissionKind, Emission, Absorption


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractExcitationScheme` 
    ... defines an abstract and a number of singleton types to distinguish between different schemes for
        generating configuration lists as they frequently occur in Green function and cascade computations.

    + struct NoExcitationScheme        
    ... dummy scheme for (unsupported) initialization of this abstract tpye.

    + struct DeExciteSingleElectron        
    ... generates all excitations and de-excitations of a single electron from a given list of bound
        electron configurations. The number of electrons of the generated configurations is the same as 
        for the given bound configurations.

    + struct DeExciteTwoElectrons       
    ... generates all excitations and de-excitations of one or two electrons from a given list of bound
        electron configurations. The number of electrons of the generated configurations is the same as 
        for the given bound configurations.
        
    + struct AddSingleElectron             
    ... generates configurations by just adding a single electrons to a given list of bound
        electron configurations. The number of electrons of the generated configurations is N+1.
        
    + struct ExciteByCapture             
    ... generates all excitations and de-excitations of one or more electron from a given list of bound
        electron configurations, together with an capture of an additional electron. 
        The number of electrons of the generated configurations is N+1.
"""
abstract type  AbstractExcitationScheme                               end
struct         NoExcitationScheme      <:  AbstractExcitationScheme   end
struct         DeExciteSingleElectron  <:  AbstractExcitationScheme   end
struct         DeExciteTwoElectrons    <:  AbstractExcitationScheme   end
struct         AddSingleElectron       <:  AbstractExcitationScheme   end
struct         ExciteByCapture         <:  AbstractExcitationScheme   end
    
export  AbstractExcitationScheme, NoExcitationScheme, DeExciteSingleElectron, DeExciteTwoElectrons, AddSingleElectron, ExciteByCapture

#################################################################################################################################
#################################################################################################################################



"""
`abstract type Basics.AbstractEmpiricalSettings` 
    ... defines an abstract type to distinguish between different settings of empirical processes/computations.
"""
abstract type  AbstractEmpiricalSettings                end

export  AbstractEmpiricalSettings

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractFieldValue` 
    ... to specify and deal with the function values of different physical fields, such as f(x), f(x,y,z), 
        f(rho,phi), f(r,theta,phi), and for both, scalar and vector fields. Often, such field values are the outcome
        of some computation which can be used for integration, display, etc.
        
    + struct Cartesian2DFieldValue{Type}     ... to specify a field value of type T in terms of x, y.
    + struct Cartesian3DFieldValue{Type}     ... to specify a field value of type T in terms of x, y, z.
    + struct PolarFieldValue{Type}           ... to specify a field value of type T in terms of rho, phi.
    + struct SphericalFieldValue{Type}       ... to specify a field value of type T in terms of r, theta, phi.
    
"""
abstract type  AbstractFieldValue               end


"""
`struct  Basics.Cartesian2DFieldValue{Type}   <: Basics.AbstractFieldValue`  
    ... to specify a scalar field value of type Type in terms of x, y.

    + x       ::Float64       ... x-coordinate.
    + y       ::Float64       ... y-coordinate.
    + val     ::Type          ... field value of type Type.
"""
struct  Cartesian2DFieldValue{Type}        <: Basics.AbstractFieldValue
    x         ::Float64
    y         ::Float64
    val       ::Type 
end 


# `Base.show(io::IO, value::Cartesian2DFieldValue{Type})`  ... prepares a proper printout of the variable value::Cartesian2DFieldValue{Type}.
function Base.show(io::IO, value::Cartesian2DFieldValue{Type}) 
    sa = "Cartesian 2D field value f(x,y) = f($(value.x),$(value.y)) = $(value.val)";                print(io, sa)
end


"""
`struct  Basics.Cartesian3DFieldValue{Type}   <: Basics.AbstractFieldValue`  
    ... to specify a scalar field value of type Type in terms of x, y, z.

    + x       ::Float64       ... x-coordinate.
    + y       ::Float64       ... y-coordinate.
    + z       ::Float64       ... z-coordinate.
    + val     ::Type          ... field value of type Type.
"""
struct  Cartesian3DFieldValue{Type}        <: Basics.AbstractFieldValue
    x         ::Float64
    y         ::Float64
    z         ::Float64
    val       ::Type 
end 


# `Base.show(io::IO, value::Cartesian3DFieldValue{Type})`  ... prepares a proper printout of the variable value::Cartesian3DFieldValue{Type}.
function Base.show(io::IO, value::Cartesian3DFieldValue{Type}) 
    sa = "Cartesian 3D field value f(x,y,z) = f($(value.x),$(value.y),$(value.z)) = $(value.val)";                print(io, sa)
end


"""
`struct  Basics.PolarFieldValue{Type}     <: Basics.AbstractFieldValue`  
    ... to specify a field value of type Type in terms of rho, phi.

    + rho     ::Float64       ... rho-coordinate.
    + phi     ::Float64       ... phi-coordinate.
    + val     ::Type          ... field value of type Type.
"""
struct  PolarFieldValue{Type}             <: Basics.AbstractFieldValue
    rho       ::Float64
    phi       ::Float64
    val       ::Type 
end 


# `Base.show(io::IO, value::PolarFieldValue{Type})`  ... prepares a proper printout of the variable value::PolarFieldValue{Type}.
function Base.show(io::IO, value::PolarFieldValue{Type}) 
    sa = "Polar field value f(rho,phi) = f($(value.rho),$(value.phi)) = $(value.val)";                print(io, sa)
end    


"""
`struct  Basics.SphericalFieldValue{Type}     <: Basics.AbstractFieldValue`  
    ... to specify a field value of type Type in terms of r, theta, phi.

    + r       ::Float64       ... r-coordinate.
    + theta   ::Float64       ... theta-coordinate.
    + phi     ::Float64       ... phi-coordinate.
    + val     ::Type          ... field value of type Type.
"""
struct  SphericalFieldValue{Type}             <: Basics.AbstractFieldValue
    r         ::Float64
    theta     ::Float64
    phi       ::Float64
    val       ::Type 
end 


# `Base.show(io::IO, value::SphericalFieldValue{Type})`  ... prepares a proper printout of the variable value::SphericalFieldValue{Type}.
function Base.show(io::IO, value::SphericalFieldValue{Type}) 
    sa = "Spherical field value f(r,theta,phi) = f($(value.r),$(value.theta),$(value.phi)) = $(value.val)";       print(io, sa)
end

export  AbstractFieldValue, Cartesian2DFieldValue, Cartesian3DFieldValue, PolarFieldValue, SphericalFieldValue

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractFieldVector` 
    ... to specify and deal with different physical field vectors, such as (Ax, Ay),  (Ax, Ay, Az),  (A_1, A_0, A_-1).
        Often, such field values are used to charaterize electric or magnetic fields (vector potentials).
        
    + struct Cartesian2DFieldVector{Type}    ... to specify a field vector of type T in terms of (Ax, Ay).
    + struct Cartesian3DFieldVector{Type}    ... to specify a field vector of type T in terms of (Ax, Ay, Az).
    + struct SphericalFieldVector{Type}      ... to specify a field vector of type T in terms of (A_1, A_0, A_-1).
    
"""
abstract type  AbstractFieldVector               end


"""
`struct  Basics.Cartesian2DFieldVector{Type}   <: Basics.AbstractFieldVector`  
    ... to specify a scalar field vector of type T in terms of (Ax, Ay).

    + x       ::Type       ... x-component.
    + y       ::Type       ... y-component.
"""
struct  Cartesian2DFieldVector{Type}        <: Basics.AbstractFieldVector
    x         ::Type
    y         ::Type 
end 


# `Base.show(io::IO, vector::Cartesian2DFieldVector{Type})`  ... prepares a proper printout of the variable vector::Cartesian2DFieldVector{Type}.
function Base.show(io::IO, vector::Cartesian2DFieldVector{Type}) 
    sa = "Cartesian 2D field vector (Ax, Ay) = ($(vector.x),$(vector.y)).";                print(io, sa)
end


"""
`struct  Basics.Cartesian3DFieldVector{Type}   <: Basics.AbstractFieldVector`  
    ... to specify a scalar field vector of type T in terms of (Ax, Ay, Az).

    + x       ::Type       ... x-component.
    + y       ::Type       ... y-component.
    + z       ::Type       ... z-component.
"""
struct  Cartesian3DFieldVector{Type}        <: Basics.AbstractFieldVector
    x         ::Type
    y         ::Type 
    z         ::Type 
end 


# `Base.show(io::IO, vector::Cartesian3DFieldVector{Type})`  ... prepares a proper printout of the variable vector::Cartesian3DFieldVector{Type}.
function Base.show(io::IO, vector::Cartesian3DFieldVector{Type}) 
    sa = "Cartesian 3D field vector (Ax, Ay, Az) = ($(vector.x),$(vector.y),$(vector.z))";                print(io, sa)
end

export  AbstractFieldVector, Cartesian2DFieldVector, Cartesian3DFieldVector, SphericalFieldVector


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractEmField` 
    ... defines an abstract type to distinguish between different static and time-dependent electric and magnetic fields;
        this is useful for atomic compass simulations and for the computation of Stark shifts.
        
    + NoEmField               ... No electric and magnetic field is defined.
    + StaticField             ... A static field that is characterized by its (real) amplitude into a given direction.
    + TimeHarmonicField       ... Define a time-harmonic field that is characterized by its (real) amplitude into a 
                                  given direction and a (harmonic) frequency.
"""
abstract type  AbstractEmField                    end
struct         NoEmField    <:  AbstractEmField   end


# `Base.show(io::IO, field::AbstractEmField)`  ... prepares a proper printout of the variable field::AbstractEmField.
function Base.show(io::IO, field::AbstractEmField) 
    print(io, string(field) )
end


# `Base.string(field::AbstractEmField)`  ... provides a proper printout of the variable field::AbstractEmField.
function Base.string(field::AbstractEmField) 
    if       field == NoEmField()     return("No electric and magnetic field is defined.")
    else     error("stop a")
    end
end 


"""
`struct  Basics.StaticField   <:  AbstractEmField`  
    ... to specify a static -- electric or magnetic -- field that is characterized by its (real) amplitude into a 
        given direction.

    + amplitude    ::Basics.Cartesian3DFieldVector{Float64}  ... 3D vector that represents the amplitude A_o of the field.
"""
struct  StaticField           <:  AbstractEmField
    amplitude      ::Basics.Cartesian3DFieldVector{Float64}
end 


# `Base.show(io::IO, field::StaticField)`  ... prepares a proper printout of the variable field::StaticField.
function Base.show(io::IO, field::StaticField) 
    sa = "Static field with amplitude A_o = $(field.amplitude)."
    print(io, sa)
end


"""
`struct  Basics.TimeHarmonicField   <:  AbstractEmField`  
    ... to specify a time-harmonic -- electric or magnetic -- field that is characterized by its (real) amplitude into a 
        given direction and a (harmonic) frequency.

    + amplitude    ::Basics.Cartesian3DFieldVector{Float64}  ... 3D vector that represents the amplitude A_o of the field.
    + omega        ::Float64                                ... Frequency of the time-harmonic field.
                                                    
"""
struct  TimeHarmonicField           <:  AbstractEmField
    amplitude      ::Basics.Cartesian3DFieldVector{Float64}
    omega          ::Float64
end 


# `Base.show(io::IO, field::TimeHarmonicField)`  ... prepares a proper printout of the variable field::TimeHarmonicField.
function Base.show(io::IO, field::TimeHarmonicField) 
    sa = "Time-harmonic field with frequency $(field.omega) and amplitude A_o = $(field.amplitude)."
    print(io, sa)
end

export AbstractEmField, NoEmField, StaticField, TimeHarmonicField


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractLevelPopulation` 
    ... defines an abstract and a number of singleton types to distinguish between different level population (models).

    + struct BoltzmannLevelPopulation     ... to represent a Boltzmann level population.       
    + struct SahaLevelPopulation          ... to represent a Saha level population.        
"""
abstract type  AbstractLevelPopulation                                end
struct    BoltzmannLevelPopulation      <:  AbstractLevelPopulation   end
struct    SahaLevelPopulation           <:  AbstractLevelPopulation   end

export  AbstractLevelPopulation, BoltzmannLevelPopulation, SahaLevelPopulation

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractMesh` 
    ... to specify different mesh types in terms of their basic parameters; these meshes can be used, for example,
        for integration or for defining the representation of observables. 
        
        
    + struct Cartesian2DMesh    ... to specify a 2D Cartesian mesh in terms of x and y.
    + struct GLegenreMesh       ... to specify a Gauss-Legendre mesh in terms of [a,b] and number of zeros.
    + struct LinearMesh         ... to specify a linear mesh in terms of [a,b] and number of points.
    + struct PolarMesh          ... to specify a 2D mesh for rho and phi.
    + struct SphercialMesh      ... to specify a 3D mesh in terms of r, theta, phi.
    
"""
abstract type  AbstractMesh                     end


"""
`struct  Basics.Cartesian2DMesh   <: Basics.AbstractMesh`  
    ... to specify a 2D Cartesian mesh in terms of x and y.

    + xMesh       ::Basics.AbstractMesh       ... mesh for the x-coordinate.
    + yMesh       ::Basics.AbstractMesh       ... mesh for the y-coordinate.
"""
struct  Cartesian2DMesh   <: Basics.AbstractMesh
    xMesh         ::Basics.AbstractMesh
    yMesh         ::Basics.AbstractMesh
end 


# `Base.show(io::IO, cMesh::Cartesian2DMesh)`  ... prepares a proper printout of the variable cMesh::Cartesian2DMesh.
function Base.show(io::IO, cMesh::Cartesian2DMesh) 
    sa = "2D Cartesian mesh with xMesh = $(cMesh.xMesh) and yMesh = $(cMesh.yMesh)";   print(io, sa)
end


"""
`struct  Basics.GLegenreMesh   <: Basics.AbstractMesh`  
    ... to specify a Gauss-Legendre mesh in terms of [a,b] and number of zeros.

    + a           ::Float64       ... mesh as defined in the interval [a,b].
    + b           ::Float64
    + NoZeros     ::Int64         ... Number of GL zeros of the mesh.
"""
struct  GLegenreMesh   <: Basics.AbstractMesh
    a             ::Float64
    b             ::Float64
    NoZeros       ::Int64 
end 


"""
`Basics.GLegenreMesh()`  ... constructor for the default settings of GLegenreMesh.
"""
function GLegenreMesh()
    GLegenreMesh( 0., 1.0, 4)
end


# `Base.show(io::IO, glMesh::GLegenreMesh)`  ... prepares a proper printout of the variable glMesh::GLegenreMesh.
function Base.show(io::IO, glMesh::GLegenreMesh) 
    sa = "Gauss-Legendre mesh for the interval  [a,b] = [$(glMesh.a),$(glMesh.b)]  and with $(glMesh.NoZeros) zeros."
    print(io, sa)
end


"""
`struct  Basics.LinearMesh   <: Basics.AbstractMesh`  
    ... to specify a linear mesh in terms of [a,b] and number of mesh points.

    + a           ::Float64       ... mesh as defined in the interval [a,b].
    + b           ::Float64
    + NoPoints    ::Int64         ... Number of mesh points, including a, b.
"""
struct  LinearMesh   <: Basics.AbstractMesh
    a             ::Float64
    b             ::Float64
    NoPoints      ::Int64 
end 


"""
`Basics.LinearMesh()`  ... constructor for the default settings of LinearMesh.
"""
function LinearMesh()
    GLegenreMesh( 0., 1.0, 10)
end


# `Base.show(io::IO, lMesh::LinearMesh)`  ... prepares a proper printout of the variable lMesh::LinearMesh.
function Base.show(io::IO, lMesh::LinearMesh) 
    sa = "Linear mesh for the interval  [a,b] = [$(glMesh.a),$(glMesh.b)]  and with $(lMesh.NoPoints) points."
    print(io, sa)
end


"""
`struct  Basics.PolarMesh   <: Basics.AbstractMesh`  
    ... to specify a 2D polar mesh in terms of rho and phi.

    + rhoMesh     ::Basics.AbstractMesh       ... mesh for the rho-coordinate.
    + phiMesh     ::Basics.AbstractMesh       ... mesh for the phi-coordinate.
"""
struct  PolarMesh   <: Basics.AbstractMesh
    rhoMesh       ::Basics.AbstractMesh
    phiMesh       ::Basics.AbstractMesh
end 


"""
`Basics.PolarMesh()`  ... constructor for the default settings of PolarMesh.
"""
function PolarMesh()
    PolarMesh( Basics.GLegenreMesh(0.0, 1.0, 12), Basics.GLegenreMesh(0.0, 2pi, 12))
end


# `Base.show(io::IO, pMesh::PolarMesh)`  ... prepares a proper printout of the variable pMesh::PolarMesh.
function Base.show(io::IO, pMesh::PolarMesh) 
    sa = "2D polar mesh with  rhoMesh::$(typeof(pMesh.rhoMesh))  and  phiMesh::$(typeof(pMesh.phiMesh))";     print(io, sa)
end


"""
`struct  Basics.SphericalMesh   <: Basics.AbstractMesh`  
    ... to specify a 3D polar mesh in terms of r, theta and phi.

    + rMesh       ::Basics.AbstractMesh         ... mesh for the r-coordinate.
    + thetaMesh   ::Basics.AbstractMesh         ... mesh for the theta-coordinate.
    + phiMesh     ::Basics.AbstractMesh         ... mesh for the phi-coordinate.
"""
struct  SphericalMesh   <: Basics.AbstractMesh
    rMesh         ::Basics.AbstractMesh
    thetaMesh     ::Basics.AbstractMesh
    phiMesh       ::Basics.AbstractMesh
end 


"""
`Basics.SphericalMesh()`  ... constructor for the default settings of Spherical.
"""
function SphericalMesh()
    SphericalMesh( Basics.GLegenreMesh(0., 1.0, 12), Basics.GLegenreMesh(0., pi, 8), Basics.GLegenreMesh(0., 2pi, 12) )
end


# `Base.show(io::IO, sMesh::SphericalMesh)`  ... prepares a proper printout of the variable sMesh::SphericalMesh.
function Base.show(io::IO, sMesh::SphericalMesh) 
    sa = "3D spherical mesh with  rMesh::$(typeof(sMesh.rMesh)),  thetaMesh::$(typeof(sMesh.thetaMesh))  and  " * 
            "phiMesh::$(typeof(sMesh.phiMesh))"
    print(io, sa)
end

export  AbstractMesh, Cartesian2DMesh, GLegenreMesh, LinearMesh, PolarMesh, SphercialMesh

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractPlasmaModel` 
    ... defines an abstract and a number of singleton types for the the (allowed) plasma models.

    + NoPlasmaModel                 ... No plasma model defined.
    + DebyeHueckelModel             ... Debye-Hueckel plasma model.
    + IonSphereModel                ... Ion-sphere (not yet supported).
    + StewartPyattModel             ... Stewart-Pyatt model (not yet supported).
    + WithoutAutoionizationModel    ... Just excludes all autoionizing levels; no original plasma model.
"""
abstract type  AbstractPlasmaModel                                    end
struct         NoPlasmaModel                <:  AbstractPlasmaModel   end
struct         DebyeBox                     <:  AbstractPlasmaModel   end
struct         WithoutAutoionizationModel   <:  AbstractPlasmaModel   end


# `Base.show(io::IO, model::AbstractPlasmaModel)`  ... prepares a proper printout of the variable model::AbstractPlasmaModel.
function Base.show(io::IO, model::AbstractPlasmaModel) 
    print(io, string(model) )
end


# `Base.string(model::AbstractPlasmaModel)`  ... provides a proper printout of the variable model::AbstractPlasmaModel.
function Base.string(model::AbstractPlasmaModel) 
    if       model == NoPlasmaModel()               return("No plasma model.")
    elseif   model == DebyeBox()                    return("Debye-box model.")
    elseif   model == WithoutAutoionizationModel()  return("Without autoionization model (just exclude all autoionizing levels).")
    else     error("stop a")
    end
end 


"""
`struct  Basics.DebyeHueckelModel   <:  AbstractPlasmaModel`  
    ... to specify (the parameters of) a Debye-Hückel potential.

    + debyeLength  ::Float64               ... the Debye length D [a_o].
    + radius       ::Float64               ... the Debye radius R_D [a_o].
"""
struct  DebyeHueckelModel   <:  AbstractPlasmaModel
    debyeLength   ::Float64
    debyeRadius   ::Float64
end 


"""
`Basics.DebyeHueckelModel()`  ... constructor for the default settings of Basics.DebyeHueckelModel().
"""
function DebyeHueckelModel()
    DebyeHueckelModel( 0.1, 0. )
end


# `Base.show(io::IO, model::DebyeHueckelModel)`  ... prepares a proper printout of the variable model::DebyeHueckelModel.
function Base.show(io::IO, model::DebyeHueckelModel) 
    sa = "Debye-Hueckel plasma model with Debye length D = $(model.debyeLength) a_o and radius R_D = $(model.debyeRadius) a_o."
    print(io, sa)
end


"""
`struct  Basics.IonSphereModel   <:  AbstractPlasmaModel`  
    ... to specify (the parameters of) a ion-sphere potential.

    + radius          ::Float64               ... the ion-sphere radius R [a_o].
    + electronDensity ::Float64               ... electron density n_e (T).
"""
struct  IonSphereModel           <:  AbstractPlasmaModel
    radius            ::Float64
    electronDensity   ::Float64 
end 


"""
`Basics.IonSphereModel()`  ... constructor for the default settings of Basics.IonSphereModel().
"""
function IonSphereModel()
    IonSphereModel( 0.9, 0. )
end


# `Base.show(io::IO, model::IonSphereModel)`  ... prepares a proper printout of the variable model::IonSphereModel.
function Base.show(io::IO, model::IonSphereModel) 
    sa = "Ion-sphere plasma model with (Wigner-Seitz radius) R = $(model.radius) a_o and electron density n_e (T) = " *
         "$(model.electronDensity)"
    print(io, sa)
end


"""
`struct  Basics.StewartPyattModel   <:  AbstractPlasmaModel`  
    ... to specify (the parameters of) a Stewart-Pyatt plasma model.

    + radius          ::Float64               ... the Stewart-Pyatt radius R [a_o].
    + electronDensity ::Float64               ... electron density n_e (T).
    + lambda          ::Float64               ... lambda^(ST) (T, ne, ni)
"""
struct  StewartPyattModel           <:  AbstractPlasmaModel
    radius            ::Float64
    electronDensity   ::Float64 
    lambda            ::Float64 
end 


"""
`Basics.StewartPyattModel()`  ... constructor for the default settings of Basics.StewartPyattModel().
"""
function StewartPyattModel()
    StewartPyattModel( 0., 0., 0. )
end


# `Base.show(io::IO, model::StewartPyattModel)`  ... prepares a proper printout of the variable model::StewartPyattModel.
function Base.show(io::IO, model::StewartPyattModel) 
    sa = "Stewart-PyattModel plasma model with (Wigner-Seitz radius) R = $(model.radius) a_o, electron density n_e (T) = " *
         "$(model.electronDensity), and lambda = $(model.lambda)"
    print(io, sa)
end

export  AbstractPlasmaModel, NoPlasmaModel, DebyeHueckelModel, DebeyBox, IonSphereModel, StewartPyattModel

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractLineShiftSettings` 
    ... defines an abstract and a number of singleton types for the the (allowed) plasma models.

    + Basics.NoLineShiftSettings         ... No line-shift settings are defined.
    + AutoIonization.PlasmaSettings      ... Settings for Auger-line computations.
    + PhotoIonization.PlasmaSettings     ... Settings for photoionization-line computations.
"""
abstract type  AbstractLineShiftSettings                end

# `Base.show(io::IO, settings::AbstractLineShiftSettings)`  ... prepares a proper printout of the variable settings::AbstractLineShiftSettings.
function Base.show(io::IO, settings::AbstractLineShiftSettings) 
    print(io, string(settings) )
end


# `Base.string(settings::AbstractLineShiftSettings)`  ... provides a proper printout of the variable settings::AbstractLineShiftSettings.
function Base.string(settings::AbstractLineShiftSettings) 
    if       typeof(settings) == Basics.NoLineShiftSettings        return("No line-shift settings.")
    else     return("Plasma settings for $(typeof(settings)) computations.")
    end
end 

export  AbstractLineShiftSettings, NoLineShiftSettings

#################################################################################################################################
#################################################################################################################################

"""
`abstract type Basics.AbstractPolarization` 
    ... defines an abstract type to comprise various polarizations of light and electron beams.

    + LinearPolarization        ... to specify a linearly-polarized pulse/beam.
    + LeftCircular              ... to specify a left-circularly polarized pulse/beam.
    + RightCircular             ... to specify a right-circularly polarized pulse/beam.
    + LeftElliptical            ... to specify an elliptically polarized pulse/beam.
    + RightElliptical           ... to specify an elliptically polarized pulse/beam.
    + NonePolarization          ... to specify an upolarized pulse/beam.
    + DensityMatrixPolarization ... to specify the polarization of a pulse/beam by its (2x2) 
                                    density matrix (not yet).
"""
abstract type  AbstractPolarization  end

struct         LinearPolarization      <:  AbstractPolarization   end
struct         LeftCircular            <:  AbstractPolarization   end
struct         RightCircular           <:  AbstractPolarization   end


"""
`struct     Basics.LeftElliptical          <:  Basics.AbstractPolarization`   

        + ellipticity      ::Float64     ... Ellipticity of the beam in the range 0...1.
"""
struct         LeftElliptical          <:  AbstractPolarization
        ellipticity      ::Float64
end


"""
`struct     Basics.RightElliptical          <:  Basics.AbstractPolarization`   

        + ellipticity      ::Float64     ... Ellipticity of the beam in the range 0...1.
"""
struct         RightElliptical          <:  AbstractPolarization
        ellipticity      ::Float64
end


struct         NonePolarization        <:  AbstractPolarization   end

function Base.string(pol::LinearPolarization)   return( "linearly-polarized" )            end
function Base.string(pol::LeftCircular)         return( "left-circularly polarized" )     end
function Base.string(pol::RightCircular)        return( "right-circularly polarized" )    end
function Base.string(pol::LeftElliptical)       return( "left-elliptically polarized with ellipticity $(pol.ellipticity)" )   end
function Base.string(pol::RightElliptical)      return( "right-elliptically polarized with ellipticity $(pol.ellipticity)" )  end
function Base.string(pol::NonePolarization)     return( "unpolarized" )                   end
    
export  AbstractPolarization, LinearPolarization, LeftCircular, RightCircular, LeftElliptical, RightElliptical, NonePolarization, 
        DensityMatrixPolarization

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractPotential` 
    ... defines an abstract and a number of singleton types to distinguish between different (electronic) atomic potentials.

    + struct DFSpotential     ... to represent a Dirac-Fock-Slater potential.        
    + struct CoreHartree      ... to represent a core-Hartree potential.        
    + struct KohnSham         ... to represent a Kohn-Sham potential.        
    + struct HartreeSlater    ... to represent a Hartree-Slater potential.        
"""
abstract type  AbstractPotential                        end
struct    DFSpotential          <:  AbstractPotential   end
struct    CoreHartree           <:  AbstractPotential   end
struct    KohnSham              <:  AbstractPotential   end
struct    HartreeSlater         <:  AbstractPotential   end

export  AbstractPotential, DFSpotential, CoreHartree, KohnSham, HartreeSlater

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractProcessSettings` 
    ... defines an abstract type to distinguish between different settings of atomic processes.
"""
abstract type  AbstractProcessSettings                  end
struct   NoProcessSettings  <: AbstractProcessSettings  end
    
export   AbstractProcessSettings

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractProcess` 
    ... defines an abstract and a number of singleton types to distinguish different atomic processes.

    + struct Auger            ... Auger transitions, i.e. single autoionization or the emission of a single free electron into the continuum.
    + struct AugerInPlasma    ... Auger transitions but calculated for a specified plasma model.
    + struct Compton          ... Rayleigh-Compton scattering cross sections.
    + struct Coulex           ... Coulomb-excitation of target or projeticle electrons by fast, heavy ions.
    + struct Coulion          ... Coulomb-ionization of target or projeticle electrons by fast, heavy ions.
    + struct Dierec           ... di-electronic recombination, i.e. the dielectronic capture of a free electron and the subsequent emission of a photon.
    + struct DoubleAuger      ... Double Auger rates.
    + struct ImpactExcAuto    ... di-electronic recombination, i.e. the dielectronic capture of a free electron and the subsequent emission of a photon.
    + struct MultiPhotonDE    ... multi-photon excitation and decay rates, including 2-photon, etc. processes.
    + struct MultiPI          ... multi-photon (single-electron) ionization.
    + struct MultiPDI         ... multi-photon (single-electron) double ionization.
    + struct Photo            ... Photoionization processes, i.e. the emission of a single free electron into the continuum due to an external light field.
    + struct PhotoDouble      ... Photo-double ionization rates.
    + struct PhotoExc         ... Photoexcitation rates.
    + struct PhotoExcFluor    ... photoexcitation fluorescence rates and cross sections.
    + struct PhotoExcAuto     ... photoexcitation autoionization cross sections and collision strengths.
    + struct PhotoInPlasma    ... Photoionization processes but calculated for a specified plasma model.
    + struct PhotoIonFluor    ... photoionization fluorescence rates and cross sections.
    + struct PhotoIonAuto     ... photoionization autoionization cross sections and collision strengths.
    + struct Radiative        ... Radiative (multipole) transitions between bound-state levels of the same charge state.
    + struct Rec              ... radiative electron capture, i.e. the capture of a free electron with the simultaneous emission of a photon.
    + struct Eimex            ... electron-impact excitation cross sections and collision strengths.
    + struct RAuger           ... Radiative Auger rates.
"""
abstract type  AbstractProcess                          end
struct    NoProcess             <:  AbstractProcess     end
struct    Auger                 <:  AbstractProcess     end
struct    AugerInPlasma         <:  AbstractProcess     end
struct    Compton               <:  AbstractProcess     end
struct    Coulex                <:  AbstractProcess     end
struct    Coulion               <:  AbstractProcess     end
struct    Dierec                <:  AbstractProcess     end
struct    DoubleAuger           <:  AbstractProcess     end
struct    ElecCapture           <:  AbstractProcess     end
struct    ImpactExcAuto         <:  AbstractProcess     end
struct    InternalConv          <:  AbstractProcess     end
struct    MultiPhotonDE         <:  AbstractProcess     end
struct    MultiPI               <:  AbstractProcess     end
struct    MultiPDI              <:  AbstractProcess     end
struct    Photo                 <:  AbstractProcess     end
struct    PhotoDouble           <:  AbstractProcess     end
struct    PhotoExc              <:  AbstractProcess     end
struct    PhotoExcFluor         <:  AbstractProcess     end
struct    PhotoExcAuto          <:  AbstractProcess     end
struct    PhotoInPlasma         <:  AbstractProcess     end
struct    PhotoIonFluor         <:  AbstractProcess     end
struct    PhotoIonAuto          <:  AbstractProcess     end
struct    Radiative             <:  AbstractProcess     end
struct    Rec                   <:  AbstractProcess     end
struct    Eimex                 <:  AbstractProcess     end
struct    RAuger                <:  AbstractProcess     end
struct    PairA1P               <:  AbstractProcess     end
    

export  AbstractProcess, NoProcess, Auger, AugerInPlasma, Compton, Coulex, Coulion, Dierec, DoubleAuger, Eimex, ElecCapture, 
        ImpactExcAuto, InternalConv, MultiPhotonDE, MultiPI, MultiPDI, Photo, PhotoDouble, PhotoExc, PhotoExcAuto, PhotoExcFluor, 
        PhotoInPlasma, PhotoIonAuto, PhotoIonFluor, Radiative, RAuger, Rec, PairA1P, Coulion

function Base.string(propc::NoProcess)          return( "no process" )                         end
function Base.string(propc::Auger)              return( "Auger" )                              end
function Base.string(propc::AugerInPlasma)      return( "Auger in plasma" )                    end
function Base.string(propc::Compton)            return( "Rayleigh-Compton" )                   end
function Base.string(propc::Dierec)             return( "Dielectronic recombination" )         end
function Base.string(propc::DoubleAuger)        return( "Double Auger" )                       end
function Base.string(propc::ElecCapture)        return( "Electron capture" )                   end
function Base.string(propc::ImpactExcAuto)      return( "ImpactExcAuto" )                      end
function Base.string(propc::InternalConv)       return( "InternalConv" )                       end
function Base.string(propc::MultiPhotonDE)      return( "multi-photon excitation & decay" )    end
function Base.string(propc::Photo)              return( "Photo-Ionization" )                   end
function Base.string(propc::PhotoExc)           return( "Photo-Excitation" )                   end
function Base.string(propc::PhotoDouble)        return( "single-photon double ionization" )    end
function Base.string(propc::PhotoExcFluor)      return( "Photo-Excitation-Fluoresence" )       end
function Base.string(propc::PhotoExcAuto)       return( "Photo-Excitation-Autoionization" )    end
function Base.string(propc::PhotoInPlasma)      return( "Photo in plasma" )                    end
function Base.string(propc::PhotoIonFluor)      return( "Photo-Ionization-Fluoresence" )       end
function Base.string(propc::PhotoIonAuto)       return( "Photo-Ionization-Autoionization" )    end
function Base.string(propc::Radiative)          return( "Radiative" )                          end
function Base.string(propc::RAuger)             return( "Radiative Auger" )                    end
function Base.string(propc::Rec)                return( "Rec" )                                end

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractPropertySettings` 
    ... defines an abstract type to distinguish between different settings of atomic level properties.
"""
abstract type  AbstractPropertySettings                 end

export AbstractPropertySettings 

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractQuantizationAxis` 
    ... defines an abstract type to distinguish between different choices of the quantization axis in light-atom 
        interactions.
        
    + DefaultQuantizationAxis      ... Use the (default) z-axis for quantization of atomic levels.
    + StaticQuantizationAxis       ... Define a static quantization axis in terms of a unit vector.
    + HarmonicQuantizationAxis     ... Define a time-harmonic quanization axis in terms of a unit vector and a frequency.
"""
abstract type  AbstractQuantizationAxis                 end
struct         DefaultQuantizationAxis  <:  AbstractQuantizationAxis   end


# `Base.show(io::IO, axis::AbstractQuantizationAxis)`  ... prepares a proper printout of the variable axis::AbstractQuantizationAxis.
function Base.show(io::IO, axis::AbstractQuantizationAxis) 
    print(io, string(axis) )
end


# `Base.string(axis::AbstractQuantizationAxis)`  ... provides a proper printout of the variable axis::AbstractQuantizationAxis.
function Base.string(axis::AbstractQuantizationAxis) 
    if       axis == DefaultQuantizationAxis()     return("Default (z-) quantization axis.")
    else     error("stop a")
    end
end 


"""
`struct  Basics.StaticQuantizationAxis   <:  AbstractQuantizationAxis`  
    ... to specify a static quantization axis in terms of a unit vector for its direction.

    + nVector      ::Basics.Cartesian3DFieldValue{Float64}
        ... 3D unit vector that specifies the quantization axis.
"""
struct  StaticQuantizationAxis   <:  AbstractQuantizationAxis
    nVector        ::Basics.Cartesian3DFieldValue{Float64}
end 


# `Base.show(io::IO, axis::StaticQuantizationAxis)`  ... prepares a proper printout of the variable axis::StaticQuantizationAxi.
function Base.show(io::IO, axis::StaticQuantizationAxis) 
    sa = "Static quantization axis along nVector = $(axis.nVector)."
    print(io, sa)
end



"""
`struct  Basics.HarmonicQuantizationAxis   <:  AbstractQuantizationAxis`  
    ... to specify a time-harmonic quantization axis in terms of a unit vector for its direction and a frequency omega.
        !!! It need to be explained how omega is related to the components of nVector; perhaps, some further
            further specification is required to make this axis unique.

    + nVector  ::Basics.Cartesian3DFieldValue{Float64} ... 3D unit vector that specifies the quantization axis.
    + omega    ::Float64                               ... Frequency of the time-harmonic motion.
                                                    
"""
struct  HarmonicQuantizationAxis   <:  AbstractQuantizationAxis
    nVector    ::Basics.Cartesian3DFieldValue{Float64}
    omega      ::Float64
end 


# `Base.show(io::IO, axis::HarmonicQuantizationAxis)`  ... prepares a proper printout of the variable axis::HarmonicQuantizationAxis.
function Base.show(io::IO, axis::HarmonicQuantizationAxis) 
    sa = "Time-harmonic quantization axis with frequency $(axis.omega) along nVector = $(axis.nVector)."
    print(io, sa)
end

export AbstractQuantizationAxis, DefaultQuantizationAxis, StaticQuantizationAxis, HarmonicQuantizationAxis


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractScField` 
    ... defines an abstract and a number of singleton types to distinguish between different self-consistent fields

    + struct ALField          ... to represent an average-level field.        
    + struct EOLField         ... to represent an (extended) optimized-level field.        
    + struct DFSField         ... to represent an mean Dirac-Fock-Slater field.        
    + struct DFSwCPField      ... to represent an mean Dirac-Fock-Slater with core-polarization field.        
    + struct HSField          ... to represent an mean Hartree-Slater field.        
    + struct EHField          ... to represent an mean extended-Hartree field.        
    + struct KSField          ... to represent an mean Kohn-Sham field.        
    + struct HartreeField     ... to represent an mean Hartree field.        
    + struct CHField          ... to represent an mean core-Hartree field.        
    + struct NuclearField     ... to represent a pure nuclear (potential) field.        
"""
abstract type  AbstractScField                          end
struct     ALField              <:  AbstractScField     end
struct     EOLField             <:  AbstractScField     end
struct     HSField              <:  AbstractScField     end
struct     EHField              <:  AbstractScField     end
struct     KSField              <:  AbstractScField     end
struct     HartreeField         <:  AbstractScField     end
struct     CHField              <:  AbstractScField     end 
struct     NuclearField         <:  AbstractScField     end
struct     AaDFSField           <:  AbstractScField     end   
struct     AaHSField            <:  AbstractScField     end   

    
"""
`struct  Basics.DFSField        <:  AbstractScField`  
    ... defines a type to describe a mean Dirac-Fock-Slater field.

    + strength           ::Float64   ... Strength factor of the DFS potential (default: strength=1.0)
"""
struct     DFSField             <:  AbstractScField     
    strength             ::Float64 
end

# `Basics.DFSField()`  ... defines the default strength=1.0
function DFSField()
    DFSField(1.0)
end

    
"""
`struct  Basics.DFSwCPField          <:  AbstractScField`  
    ... defines a type to describe a mean Dirac-Fock-Slater field with core-polarization.

    + corePolarization  ::CorePolarization   ... Parametrization of the core-polarization potential/contribution.
"""
struct     DFSwCPField          <:  AbstractScField   
    corePolarization    ::CorePolarization
end

export  AbstractScField, AaDFSField, AaHSField, ALField, EOLField, DFSField, DFSwCPField, HSField, NuclearField

#################################################################################################################################
#################################################################################################################################

"""
`abstract type Basics.AbstractSelection` 
    ... defines an abstract and a number of concrete types to distinguish between level- and line-selectors

    + struct LevelSelection   ... to specify a list of levels by means of their (level) indices or level symmetries.         
    + struct LineSelection    ... to specify a list of lines by means of their (level) indices or level symmetries.      
    + struct PathwaySelection ... to specify a list of lines by means of their (level) indices or level symmetries.      
    + struct ShellSelection   ... to specify a list of lines by means of their (level) indices or level symmetries.      
"""
abstract type  AbstractSelection      end


"""
`struct  Basics.LevelSelection  <  Basics.AbstractSelection`  
    ... defines a struct to specify a list of levels by means of their (level) indices or level symmetries.

    + active       ::Bool                     ... true, if some selection has been made.
    + indices      ::Array{Int64,1}           ... List of selected indices.
    + symmetries   ::Array{LevelSymmetry,1}   ... List of selected symmetries
"""
struct  LevelSelection  <:  AbstractSelection
    active         ::Bool  
    indices        ::Array{Int64,1}
    symmetries     ::Array{LevelSymmetry,1}
end


"""
`Basics.LevelSelection()`  ... constructor for an inactive LevelSelection.
"""
function  LevelSelection()
    LevelSelection( false, Int64[], LevelSymmetry[])    
end


"""
`Basics.LevelSelection(active::Bool; indices::Array{Int64,1}=Int64[], symmetries::Array{LevelSymmetry,1}=LevelSymmetry[])`  
    ... constructor for specifying the details of a LevelSelection.
"""
function  LevelSelection(active::Bool; indices::Array{Int64,1}=Int64[], symmetries::Array{LevelSymmetry,1}=LevelSymmetry[])
    if  active   LevelSelection( true, indices, symmetries)  
    else         LevelSelection()
    end
end

function Base.show(io::IO, selection::LevelSelection) 
    print(io, string(selection) )
end

function Base.string(selection::LevelSelection) 
    if  selection.active   sa = "LevelSelection:  indices = $(selection.indices);    symmetries = $(selection.symmetries)."
    else                   sa = "Inactive LevelSelection."
    end
    return( sa )
end


"""
`struct  Basics.LineSelection  <  Basics.AbstractSelection`  
    ... defines a struct to specify a list of level pair by means of their (level) indices or level symmetries.

    + active        ::Bool                                          ... true, if some selection has been made.
    + indexPairs    ::Array{Tuple{Int64,Int64},1}                   ... List of selected index pairs.
    + symmetryPairs ::Array{Tuple{LevelSymmetry,LevelSymmetry},1}   ... List of selected symmetry pairs.
"""
struct  LineSelection  <: Basics.AbstractSelection
    active          ::Bool  
    indexPairs      ::Array{Tuple{Int64,Int64},1}
    symmetryPairs   ::Array{Tuple{LevelSymmetry,LevelSymmetry},1}
end


"""
`Basics.LineSelection()`  ... constructor for an inactive LineSelection.
"""
function  LineSelection()
    LineSelection( false, Tuple{Int64,Int64}[], Tuple{LevelSymmetry,LevelSymmetry}[])    
end


"""
`Basics.LineSelection(active::Bool; indexPairs::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                                    symmetryPairs::Array{Tuple{LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry}[])`  
    ... constructor for specifying the details of a LineSelection.
"""
function  LineSelection(active::Bool; indexPairs::Array{Tuple{Int64,Int64},1}=Tuple{Int64,Int64}[],
                                        symmetryPairs::Array{Tuple{LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry}[])
    if  active   LineSelection( true, indexPairs, symmetryPairs)  
    else         LineSelection()
    end
end

function Base.show(io::IO, selection::LineSelection) 
    print(io, string(selection) )
end

function Base.string(selection::LineSelection) 
    if  selection.active   sa = "LineSelection:  indexPairs = $(selection.indexPairs);    symmetryPairs = $(selection.symmetryPairs)."
    else                   sa = "Inactive LineSelection."
    end
    return( sa )
end


"""
`struct  Basics.PathwaySelection  <  Basics.AbstractSelection`  
    ... defines a struct to specify a list of level triple (pathways) by means of their (level) indices or level symmetries.

    + active          ::Bool                                          ... true, if some selection has been made.
    + indexTriples    ::Array{Tuple{Int64,Int64,Int64},1}             ... List of selected index triples.
    + symmetryTriples ::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}  ... List of selected symmetry triples.
"""
struct  PathwaySelection  <:  Basics.AbstractSelection
    active            ::Bool  
    indexTriples      ::Array{Tuple{Int64,Int64,Int64},1}  
    symmetryTriples   ::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}
end


"""
`Basics.PathwaySelection()`  ... constructor for an inactive PathwaySelection.
"""
function  PathwaySelection()
    PathwaySelection( false, Tuple{Int64,Int64,Int64}[], Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])    
end


"""
`Basics.PathwaySelection(active::Bool; indexTriples::Array{Tuple{Int64,Int64,Int64},1}=Tuple{Int64,Int64,Int64}[],
                symmetryTriples::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])`  
    ... constructor for specifying the details of a PathwaySelection.
"""
function  PathwaySelection(active::Bool; indexTriples::Array{Tuple{Int64,Int64,Int64},1}=Tuple{Int64,Int64,Int64}[],
            symmetryTriples::Array{Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry},1}=Tuple{LevelSymmetry,LevelSymmetry,LevelSymmetry}[])
    if  active   PathwaySelection( true, indexTriples, symmetryTriples)  
    else         PathwaySelection()
    end
end

function Base.show(io::IO, selection::PathwaySelection) 
    print(io, string(selection) )
end

function Base.string(selection::PathwaySelection) 
    if  selection.active   sa = "PathwaySelection:  indexTriples = $(selection.indexTriples);    symmetryTriples = $(selection.symmetryTriples)."
    else                   sa = "Inactive PathwaySelection."
    end
    return( sa )
end


"""
`struct  Basics.ShellSelection  <  Basics.AbstractSelection`  
    ... defines a struct to specify a list of shells by means of different constructors.

    + active       ::Bool                     ... true, if some selection has been made.
    + shells       ::Array{Shell,1}           ... List of explicitly selected shells.
    + lSymmetries  ::Array{Int64,1}           ... List of selected l-symmetries
"""
struct  ShellSelection  <: Basics.AbstractSelection
    active         ::Bool  
    shells         ::Array{Shell,1}  
    lSymmetries    ::Array{Int64,1}
end


"""
`Basics.ShellSelection()`  ... constructor for an inactive ShellSelection.
"""
function  ShellSelection()
    ShellSelection( false, Shell[], Int64[])    
end


"""
`Basics.ShellSelection(active::Bool; shells::Array{Shell,1}=Shell[], lSymmetries::Array{Int64,1}=Int64[])`  
    ... constructor for specifying the details of a ShellSelection.
"""
function  ShellSelection(active::Bool; shells::Array{Shell,1}=Shell[], lSymmetries::Array{Int64,1}=Int64[])
    if  active   ShellSelection( true, shells, lSymmetries)  
    else         ShellSelection()
    end
end

function Base.show(io::IO, selection::ShellSelection) 
    print(io, string(selection) )
end

function Base.string(selection::ShellSelection) 
    if  selection.active   sa = "ShellSelection:  shells = $(selection.shells);    symmetries = $(selection.lSymmetries)."
    else                   sa = "Inactive LevelSelection."
    end
    return( sa )
end

export  AbstractSelection, LevelSelection, LineSelection, PathwaySelection, ShellSelection

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractSpectrumKind` 
    ... defines an abstract and a number of singleton types for specifying different kinds of spectra that need to be displayed

    + struct BarIntensities                         
        ... to display just intensity-bars at given x-positions.  
    + struct DiscreteLines                         
        ... to display lines (for guiding the eyes) but which are only defined at discrete x-points.
    + struct DiscretePoints                         
        ... to display discrete points that are defined at discrete x-points.
    + struct LorentzianIntensitiesConstantWidths    
        ... to display the superposition of Lorentzians with given position and intensity but constant widths.        
    + struct LorentzianIntensities  
        ... to display the superposition of Lorentzians with given position, intensity and individual widths.     
"""
abstract type  AbstractSpectrumKind                                          end
struct     BarIntensities                       <:  AbstractSpectrumKind     end
struct     DiscreteLines                        <:  AbstractSpectrumKind     end
struct     DiscretePoints                       <:  AbstractSpectrumKind     end
struct     LorentzianIntensitiesConstantWidths  <:  AbstractSpectrumKind     end
struct     LorentzianIntensities                <:  AbstractSpectrumKind     end

export  AbstractSpectrumKind, BarIntensities, DiscreteLines, DiscretePoints, LorentzianIntensitiesConstantWidths, LorentzianIntensities

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractWarning` 
    ... defines an abstract and a number of singleton types for dealing with warnings that are made during a run or REPL session.
        Cf. Defaults.warn().

    + AddWarning        ... add a Warning to a warningList.
    + PrintWarnings     ... print all warnings into a jac-warn.report file.
    + ResetWarnings     ... reset (empty) the warningList, usually at the beginning of a new run.to distinguish between different warnings
"""
abstract type  AbstractWarning                          end
struct     AddWarning           <:  AbstractWarning     end
struct     PrintWarnings        <:  AbstractWarning     end
struct     ResetWarnings        <:  AbstractWarning     end

export  AbstractWarning, AddWarning, PrintWarnings, ResetWarnings


"""
`abstract type Basics.AbstractIntegrationRule`
    ... defines an abstract and three singleton types to select the numerical integration rule
        used by Basics.integrate() on a radial grid.

    + NewtonCotes  ... 5-point Newton-Cotes formula.
    + SimpsonRule  ... Simpson's rule.
    + TrapezRule   ... simple trapezoid rule.
"""
abstract type  AbstractIntegrationRule                          end
struct         NewtonCotes  <:  AbstractIntegrationRule         end
struct         SimpsonRule  <:  AbstractIntegrationRule         end
struct         TrapezRule   <:  AbstractIntegrationRule         end

export  AbstractIntegrationRule, NewtonCotes, SimpsonRule, TrapezRule


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractGenerateTheme`
    ... defines an abstract and a number of singleton types to select the generation theme for Basics.generate(),
        replacing the former string-key dispatch.

    + CondensedMultiplet                   ... condense/reduce the CSF basis of a multiplet by a single weight.
    + ConfigurationListNRFromBasis         ... generate the NR configuration list from a given basis.
    + ConfigurationListNRFromConfiguration ... generate an NR configuration list from a reference configuration with excitations.
    + CsfList                              ... construct the CSF list from a single relativistic configuration.
    + OrderedShellList                     ... generate an ordered NR shell list from a set of configurations.
    + OrderedSubshellList                  ... generate an ordered relativistic subshell list from configurations or two bases.
    + SlaterTypeSpectrum                   ... generate a complete single-electron STO spectrum (positive and negative states).
    + SlaterTypeSpectrumPositive           ... generate the same but return only the positive states.
"""
abstract type  AbstractGenerateTheme                                              end
struct         CondensedMultiplet              <:  AbstractGenerateTheme          end
struct         ConfigurationListNRFromBasis    <:  AbstractGenerateTheme          end
struct         ConfigurationListNRFromConfiguration  <:  AbstractGenerateTheme    end
struct         CsfList                         <:  AbstractGenerateTheme          end
struct         OrderedShellList                <:  AbstractGenerateTheme          end
struct         OrderedSubshellList             <:  AbstractGenerateTheme          end
struct         SlaterTypeSpectrum              <:  AbstractGenerateTheme          end
struct         SlaterTypeSpectrumPositive      <:  AbstractGenerateTheme          end

export  AbstractGenerateTheme, CondensedMultiplet, ConfigurationListNRFromBasis, ConfigurationListNRFromConfiguration,
        CsfList, OrderedShellList, OrderedSubshellList, SlaterTypeSpectrum, SlaterTypeSpectrumPositive


#################################################################################################################################
#################################################################################################################################


"""
`abstract type Basics.AbstractComputeTheme`
    ... defines an abstract and a number of singleton types to select the computation theme for Basics.compute(),
        replacing the former string-key dispatch.

    + AngularCoeffsEeRatip2013    ... compute electron-electron angular coefficients via the Ratip2013 interface.
    + AngularCoeffs1pRatip2013    ... compute single-particle angular coefficients via the Ratip2013 interface.
    + AngularCoeffs1pGrasp92      ... compute single-particle angular coefficients via the Grasp92 interface.
    + CImatrixWithSymmetryJP      ... compute the CI Hamiltonian matrix for a given J^P symmetry block.
    + RadialOrbitalBunge1993      ... generate a start orbital from Bunge (1993) Roothaan-Hartree-Fock data.
    + RadialOrbitalMcLean1981     ... generate a start orbital from McLean (1981) Roothaan-Hartree-Fock data.
    + RadialOrbitalHydrogenic     ... generate a hydrogenic start orbital.
    + RadialOrbitalThomasFermi    ... generate a Thomas-Fermi start orbital.
"""
abstract type  AbstractComputeTheme                                              end
struct         AngularCoeffsEeRatip2013   <:  AbstractComputeTheme              end
struct         AngularCoeffs1pRatip2013   <:  AbstractComputeTheme              end
struct         AngularCoeffs1pGrasp92     <:  AbstractComputeTheme              end
struct         CImatrixWithSymmetryJP     <:  AbstractComputeTheme              end
struct         RadialOrbitalBunge1993     <:  AbstractComputeTheme              end
struct         RadialOrbitalMcLean1981    <:  AbstractComputeTheme              end
struct         RadialOrbitalHydrogenic    <:  AbstractComputeTheme              end
struct         RadialOrbitalThomasFermi   <:  AbstractComputeTheme              end

export  AbstractComputeTheme, AngularCoeffsEeRatip2013, AngularCoeffs1pRatip2013, AngularCoeffs1pGrasp92,
        CImatrixWithSymmetryJP, RadialOrbitalBunge1993, RadialOrbitalMcLean1981, RadialOrbitalHydrogenic, RadialOrbitalThomasFermi


"""
`abstract type Basics.AbstractDisplayTheme`
    ... defines an abstract and a number of singleton types to select the display theme for Basics.display(),
        replacing the former string-key dispatch.

    + PhysicalConstants   ... display all currently defined physical constants.
    + CurrentSettings     ... display all currently defined settings of the JAC module.
"""
abstract type  AbstractDisplayTheme                                end
struct         PhysicalConstants   <:  AbstractDisplayTheme       end
struct         CurrentSettings     <:  AbstractDisplayTheme       end

export  AbstractDisplayTheme, PhysicalConstants, CurrentSettings


"""
`abstract type Basics.AbstractPlotTheme`
    ... defines an abstract and a number of singleton types to select the plot theme for Basics.plot(),
        replacing the former string-key dispatch.

    + RadialPotentials    ... plot one or more radial potentials.
    + RadialOrbitalsLarge ... plot the large component of one or more radial orbitals.
    + RadialOrbitalsSmall ... plot the small component of one or more radial orbitals.
    + RadialOrbitalsBoth  ... plot both components of one or more radial orbitals.
"""
abstract type  AbstractPlotTheme                                   end
struct         RadialPotentials    <:  AbstractPlotTheme           end
struct         RadialOrbitalsLarge <:  AbstractPlotTheme           end
struct         RadialOrbitalsSmall <:  AbstractPlotTheme           end
struct         RadialOrbitalsBoth  <:  AbstractPlotTheme           end

export  AbstractPlotTheme, RadialPotentials, RadialOrbitalsLarge, RadialOrbitalsSmall, RadialOrbitalsBoth
