
#== September 2025

This code/module has been written without any tests so far; apart from syntax errors, a number of
logical errors need to be expected. Nonetheless, it provides a first set of functions in order to deal with the 
autoionization and (x-ray) photon emission decay of muonic atoms. Not all the functions have been realized in 
full details yet.

At present, this module is not (yet) included into the JAC code.
We need to think and discuss more carefully which of the properties and spectra can be predicted
for muonic atoms, and how they should be implemented as a new cascade scheme.

We also need to implement the two cascade schemes: Cascade.MuonicDecayScheme
together with the associated parameters and -inc- files. In the first instance, these computations are based
on pre-calculated muonic and electronic orbitals, from which all other orbitals are obtained. --- In this
scheme allow different decay processes (Auger(), Photo()) and the maximum number of released electrons ... follow
StepwiseDecay. 

Report properly about the basic assumptions/limitations in this scheme;
split Cascade-inc-Simulations  into   
Cascade-inc-Simulations-ions
Cascade-inc-Simulations-electrons
Cascade-inc-Simulations-photons
Cascade-inc-Simulations-coincidences   ... to better separate different topics, together with the associated
                                        data types.
                                        
New SimulationProperties:  MuonicXrays, ElectronicXraysAfterMuonCapture, IonDistributionAfterMuonCapture.

Likely, we need to introduce a Cascade.MuonLevel (energy, configuration, muonSubshell) to deal with full cascades
of the muon-decay.

Introduce a function  Basics.interpolate(orbital::Orbital, newGrid::Radial.Grid) to interpolate an orbital to 
a new grid.

==#


"""
`module  JAC.MuonicAtom`  
    ... a submodel of JAC that contains all methods for modeling the levels and decay of muonic atoms and ions.
    In contrast to the (purely) electronic case, MuonicLevel's are characterized by a full electronic configuration
    and a muon ins an subshell::MuonSubshell. Moreover, the electronic processes associated with the decay
    of such muonic levels are described by MuonicAtom.PhotoEmissionLine's and MuonicAtom.AutoIonizationLine's,
    which will contain a free electron in a particular partial-wave (subshell).
    
    This modules mainly provides data structures and functions that are needed to deal with muonic atoms
    in a Cascade.MuonicXrayScheme and Cascade.MuonicAugerScheme. The decay of muonic atoms assumes that muon has
    a binding energy similar or larger than the most inner (1s) electrons. The hydrogenic muon subshells, for which 
    this condition is fullfilled, can be readily estimated as is indicated during the modeling of the photon
    and electron spectra.
"""
module MuonicAtom


using Base.Threads, Distributed, Printf, ProgressMeter, SpecialFunctions,
        ..AngularMomentum, ..Defaults, ..TableStrings



"""
`struct  MuonicAtom.MuonSubshell`  ... defines a type for the allowed values of relativistic muon subhells.  

    + n        ::Int64  ... principal quantum number 
    + kappa    ::Int64  ... relativistic angular quantum number
"""
struct  MuonSubshell
    n          ::Int64 
    kappa      ::Int64 
end


"""
`MuonicAtom.MuonSubshell(sa::String)`  ... constructor for a given String, such as 1s_1/2, 2s_1/2, 2p_3/2, ... .
"""
function MuonSubshell(sa::String) 
    wa = strip( sa );   
    wb = findnext("_", wa, 1);    if  wb == nothing   error("No underscore in sa = $sa")           else   wb = wb[1]   end
    wc = findnext("/", wa, wb);   if  wc == nothing   error("No / in wa[] = $(wa[wb[1]+1:end])")   else   wc = wc[1]   end
    wd = string( wa[wb+1:wc-1] )
    n  = parse(Int64, wa[1:wb-2]);                  l  = shellNotation( string(wa[ wb-1 ]) )
    j2 = parse(Int64, wd );                         !(string(wa[end]) == "2")    &&  error("Unrecognized subshell string sa = $sa") 
    if  l - j2/2 > 0     kappa = Int64( (j2+1)/2 )    else    kappa = Int64( -(j2+1)/2 )   end
    MuonSubshell(n, kappa)   
end


# `Base.show(io::IO, sh::MuonSubshell)`  ... prepares a proper printout of the variable sh::MuonSubshell.
function Base.show(io::IO, sh::MuonSubshell) 
    print(io, string(sh) )
end


# `Base.string(sh::MuonSubshell)`  ... provides a proper printout of the variable sh::MuonSubshell.
function Base.string(sh::MuonSubshell) 
    sa = "muon:" * string(sh.n) * shellNotation( Basics.subshell_l(sh) ) * "_" * string( Basics.subshell_j(sh) )
    return( sa )
end

#######################################################################################################################
#######################################################################################################################


"""
`struct  MuonicAtom.Level`  
    ... defines a Level type for muonic atoms and ions; such levels are simply modeled by the subshell::MuonicSubshell
        of the muon and the configuration for the occupation of the electrons.

    + muonSubshell   ::MuonSubshell                 ... Subshell of the muon.
    + muonOrbital    ::Orbital                      ... Orbital of the muon but with a grid, scaled by the muon mass.
    + config         ::Configuration                ... Configuration to model the electronic shell structure.
    + orbitals       ::Dict{Subshell, Orbital}      ... Orbitals of the bound electrons
    + energy         ::Float64                      ... Total mean energy of the electron configuration.
    + totalenergy    ::Float64                      ... (Estimated) total energy of the muonic atom.
"""
struct  Level
    muonSubshell   ::MuonSubshell 
    muonOrbital    ::Orbital 
    config         ::Configuration 
    orbitals       ::Dict{Subshell, Orbital} 
    energy         ::Float64 
    totalenergy    ::Float64       
end 


"""
`MuonicAtom.Level()`  ... constructor for an 'empty' instance of a MuonicAtom.Level
"""
function  Level()
    Level( MuonSubshell(1,-1), Orbital(), Configuration("[He]"), Dict{Subshell, Orbital}(), 0., 0. )
end


# `Base.show(io::IO, level::MuonicAtom.Level)`  ... prepares a proper printout of level::MuonicAtom.Level.
function Base.show(io::IO, level::MonteCarlo.Level) 
    println(io, "muonSubshell:   $(level.muonSubshellJ)  ")
    println(io, "muonOrbitaly:   $(level.muonOrbital)  ")
    println(io, "config:         $(level.config)  ")
    println(io, "orbitals:       $(level.orbitals)  ")
    println(io, "energy:         $(level.energy)  ")
    println(io, "totalenergy:    $(level.totalenergy)  ")
end


#######################################################################################################################
#######################################################################################################################

"""
`struct  MuonicAtom.AutoIonizationChannel`  
    ... defines a  Channel type to deal with the (electronic) autoionization of muonic atoms. Since all electrons are 
        treated only in terms of their (configuration-) averaged density, the partial waves alone characterize the 
        individual autoionization channels.

    + kappa          ::Int64                 ... partial wave of the emitted electron.
    + amplitude      ::Complex{Float64}      ... associated two-particle matrix element/amplitude of the channel.
"""
struct  AutoIonizationChannel
    kappa            ::Int64  
    amplitude        ::Complex{Float64} 
end 


"""
`MuonicAtom.AutoIonizationChannel()`  ... constructor for an 'empty' instance of a MuonicAtom.AutoIonizationChannel
"""
function  AutoIonizationChannel()
    AutoIonizationChannel( 0., ComplexF64(0.) )
end


# `Base.show(io::IO, channel::AutoIonizationChannel)`  ... prepares a proper printout of channel::AutoIonizationChannel.
function Base.show(io::IO, channel::AutoIonizationChanne) 
    println(io, "kappa:          $(channel.kappa)  ")
    println(io, "amplitude:      $(channel.amplitude)  ")
end


"""
`struct  MuonicAtom.AutoIonizationLine`  
    ... defines a  Line type to deal with the (electronic) autoionization of muonic atoms, once the muon occupies
        a shell "inside" of the bound electron density. While the muon can then be described by Dirac solutions,
        the electrons are treated in terms of configurations.

    + initialLevel   ::MuonicAtom.Level                 ... initial level of the muonic atom.
    + finalLevel     ::MuonicAtom.Level                 ... final level of the muonic atom, with one electron less.
    + energy         ::Float64                          ... (Mean) free energy of the electrons.
    + channels       ::MuonicAtom.AutoIonizationChannel 
        ... Partial-wave channels for the autoionization of muonic atoms.
"""
struct  AutoIonizationLine
    initialLevel     ::MuonicAtom.Level
    finalLevel       ::MuonicAtom.Level
    energy           ::Float64  
    channels         ::MuonicAtom.AutoIonizationChannel 
end 


"""
`MuonicAtom.AutoIonizationLine()`  ... constructor for an 'empty' instance of a MuonicAtom.AutoIonizationLine
"""
function  AutoIonizationLine()
    AutoIonizationLine( MuonicAtom.Level(), MuonicAtom.Level(), 0., AutoIonizationChannel[] )
end


# `Base.show(io::IO, line::AutoIonizationLine)`  ... prepares a proper printout of line::AutoIonizationLine.
function Base.show(io::IO, line::AutoIonizationLine) 
    println(io, "initialLevel:   $(line.initialLevel)  ")
    println(io, "finalLevel:     $(line.finalLevel)  ")
    println(io, "energyl:        $(line.energy)  ")
    println(io, "channels:       $(line.channels)  ")
end


"""
`struct  MuonicAtom.PhotoEmissionChannel`  
    ... defines a  Channel type to deal with the (x-ray) photo-emission of muonic atoms. Since all electrons are 
        treated only in terms of their (configuration-) averaged density, the channels are characterized by the 
        multipoles alone.

    + multipole      ::EmMultipole           ... multipole of the emitted photon.
    + amplitude      ::Complex{Float64}      ... associated two-particle matrix element/amplitude of the channel.
"""
struct  PhotoEmissionChannel
    multipole        ::EmMultipole 
    amplitude        ::Complex{Float64} 
end 


"""
`MuonicAtom.PhotoEmissionChannel()`  ... constructor for an 'empty' instance of a MuonicAtom.PhotoEmissionChannel
"""
function  PhotoEmissionChannel()
    PhotoEmission( E1, ComplexF64(0.) )
end


# `Base.show(io::IO, channel::PhotoEmissionChannel)`  ... prepares a proper printout of channel::PhotoEmissionChannel.
function Base.show(io::IO, channel::PhotoEmissionChannel) 
    println(io, "multipole:      $(channel.multipole)  ")
    println(io, "amplitude:      $(channel.amplitude)  ")
end


"""
`struct  MuonicAtom.PhotoEmissionLine`  
    ... defines a  Line type to deal with the (x-ray) photoemission of muonic atoms, once the muon occupies
        a shell "inside" of the bound electron density. While the muon can then be described by Dirac solutions,
        the electrons are treated in terms of configurations only.

    + initialLevel   ::MuonicAtom.Level                 ... initial level of the muonic atom.
    + finalLevel     ::MuonicAtom.Level                 ... final level of the muonic atom, with one electron less.
    + energy         ::Float64                          ... Energy of the emitted photon = transition energy.
    + channels       ::MuonicAtom.PhotoEmissionChannel  ... Multipole channels of the photoemission from muonic atoms.
"""
struct  PhotoEmissionLine
    initialLevel     ::MuonicAtom.Level
    finalLevel       ::MuonicAtom.Level
    energy           ::Float64  
    channels         ::MuonicAtom.PhotoEmissionChannel 
end 


"""
`MuonicAtom.PhotoEmissionLine()`  ... constructor for an 'empty' instance of a MuonicAtom.PhotoEmissionLine
"""
function  PhotoEmissionLine()
    PhotoEmissionLine( MuonicAtom.Level(), MuonicAtom.Level(), 0., PhotoEmissionChannel[]   )
end


# `Base.show(io::IO, line::PhotoEmissionLine)`  ... prepares a proper printout of line::PhotoEmissionLine.
function Base.show(io::IO, line::PhotoEmissionLine) 
    println(io, "initialLevel:   $(line.initialLevel)  ")
    println(io, "finalLevel:     $(line.finalLevel)  ")
    println(io, "energyl:        $(line.energy)  ")
    println(io, "channels:       $(line.channels)  ")
end


#######################################################################################################################
#######################################################################################################################


"""
`MuonicAtom.computeAutoionizationAmplitude(... provide orbitals for the matrix element )`  
    ... to compute the autoionization amplitude from the interacting two muonic and two electron orbitals by taking 
        their Coulomb interaction into account. One of the electron describe the free (outgoing) electron.
"""
function computeAutoionizationAmplitude()
    amplitude = ComplexF64(0.)  
    
    return( amplitude )
end


"""
`MuonicAtom.computeAutoionizationnLine(line::MuonicAtom.AutoionizationLine)`  
    ... to compute all properties of a MuonicAtom.AutoionizationLine. A newLine::MuonicAtom.AutoionizationLine is
        returned.
"""
function computeAutoionizationLine(line::MuonicAtom.AutoionizationLine)
    newLine = MuonicAtom.AutoionizationLine()
    
    return( newLine )
end


"""
`MuonicAtom.computePhotoEmissionAmplitude(... provide orbitals for the matrix element )`  
    ... to compute the photoemission amplitude from the interacting two muonic or electron orbitals. Likely, we need 
        to distinguish between two cases, a photo emission from the muonic transitions and from electronic transitions.
        Perhaps, this can be treated together if the muonic orbitals are properly interpolated to some useful
        electronic grid.
"""
function computePhotoEmissionnAmplitude()
    amplitude = ComplexF64(0.)  
    
    return( amplitude )
end


"""
`MuonicAtom.computePhotoEmissionLine(line::MuonicAtom.PhotoEmissionLine)`  
    ... to compute all properties of a MuonicAtom.PhotoEmissionLine. A newLine::MuonicAtom.PhotoEmissionLine is
        returned.
"""
function computePhotoEmissionLine(line::MuonicAtom.PhotoEmissionLine)
    newLine = MuonicAtom.PhotoEmissionLine()
    
    return( newLine )
end


"""
`MuonicAtom.estimateEnergy(subshell::MuonicAtom.MuonSubshell)`  
    ... to estimate the binding energy of a muon in the subshell::MuonicAtom.MuonSubshell. An energy::Float64 
        [in Hartree] is returned. Dirac's energy formula is applied with the proper mass of muons.
"""
function estimateEnergy(subshell::MuonicAtom.MuonSubshell)
    energy = 0.
                
    return( energy )
end


"""
`MuonicAtom.estimateRkExpectation(subshell::MuonicAtom.MuonSubshell, k::Int64)`  
    ... to estimate the r^k expectation values of a muon in the subshell::MuonicAtom.MuonSubshell. An value::Float64 
        is returned. Dirac's expectation values is applied with the proper mass of muons.
"""
function estimateRkExpectation(subshell::MuonicAtom.MuonSubshell, k::Int64)
    value = 0.
                
    return( value )
end


"""
`MuonicAtom.generateDiracOrbitals(subshells::Array{MuonSubshell,1}, nm::Nuclear.Model, 
                                  primitives::BsplinesN.Primitives; printout::Bool=true)`
    ... generates all muonic (Dirac) orbitals from subshells for the nuclear potential as specified by nm. A set of 
        orbitals::Dict{Subshell, Orbital} is returned, for which the proper muon mass is taken into account for the 
        given nuclear model. Note that a properly chosen grid is essential to model the muonic solutions, which 
        can later be interpolated upon a suitable (electronic) grid.
"""
function generateDiracOrbitals(subshells::Array{MuonSubshell,1}, nm::Nuclear.Model, 
                               primitives::BsplinesN.Primitives; printout::Bool=true)
                
    return( nothing )
end


"""
`MuonicAtom.generateMuonSpektrum(...)`
    ... generates all muonic (Dirac) orbitals from subshells ... which are kept troughout the computations.
        ... there is no need to make these orbitals too good.
"""
function generateMuonSpektrum()
                
    return( nothing )
end


"""
`MuonicAtom.runx()`  
    ... to run 
"""
function runx()
    events = Float64[]
                
    return( events )
end



end # module
