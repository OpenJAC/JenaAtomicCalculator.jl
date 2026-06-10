
"""
`module  JAC.Cascade
... a submodel of JAC that contains all data types and methods to set-up and process (simple) 
    atomic cascade computations of various kind.
"""
module Cascade


using Dates, JLD2, Printf, FastGaussQuadrature, Distributed, ProgressMeter,
        ..AngularMomentum, ..AutoIonization, ..Basics, ..BsplinesN, ..Continuum, ..Defaults, 
        ..DecayYield, ..DielectronicRecombination, ..ElectronCapture, ..Empirical, ..Hamiltonian, ..ImpactExcitation, ..Radial, ..ManyElectron, ..Nuclear, 
        ..PeriodicTable, ..PhotoEmission, ..PhotoExcitation, ..PhotoIonization, ..PhotoRecombination, 
        ..SelfConsistent, ..Semiempirical, ..TableStrings


"""
`abstract type Cascade.AbstractCascadeScheme` 
    ... defines an abstract type to distinguish different excitation, ionization and decay schemes of an atomic cascade; see also:
    
    + struct DielectronicCaptureScheme  
        ... to model just the (dielectronic) capture and the formation of doubly-excited levels up to a maximum excitation 
            energy and for a given list of subshells; the excitation energy refers to the lowest level of the reference 
            configurations, and the cascade blocks are built only by means of the given subshells.
            NOTE: The original DielectronicCaptureScheme --> DielectronicRecombinationScheme has been renamed in 
            August 2023 in order to enlarge the consistency of the notations and code !!
    + struct DielectronicRecombinationScheme  
        ... to model the (dielectronic) recombination of an electron up to a maximum excitation energy and for a given 
            list of subshells; the excitation energy refers to the lowest level of the reference configurations, and 
            cascade blocks are built only by means of the given subshells.
    + struct ElectronExcitationScheme  
        ... to model the electron excitation spectra in terms of the direct (EIE) and resonant contributions, i.e. the 
            dielectronic capture of an electron with subsequent re-autoionization. Typical electron-excitation properties are 
            energy-dependent EIE cross sections, effective collision strengths, EIE plasma rate coefficients, and several
            others (not yet).
    + struct ElectronIonizationScheme  
        ... to model the electron ionization spectra including the direct (EII) and resonant contributions, i.e. the 
            dielectronic capture of an electron with subsequent double-autoionization. For this double autoionization, a 
            branching factor will just be estimated. Typical electron-ionization properties are energy-dependent EII cross 
            sections, effective collision strengths for impact-ionization, EII plasma rate coefficients, and several others 
            (not yet).
    + struct ExpansionOpacityScheme  
        ... to model the expansion opacity of an ion in its ground or some low-lying state; this scheme takes a maximum photon
            (transition) energy and the excitation from the fromShells to the toShells in order to select the relevant 
            configurations. These shell lists refer to the given set of reference configurations.
    + struct HollowIonScheme    
        ... to model the capture of one or several electrons into a list of subshells; various distributions of the
            electron among these shells are supported. For the subsequent decay, the list of decay shells need to be specified
            as well.
    + struct ImpactExcitationScheme    
        ... to model the (direct) electron-impact excitation  (collision strength) of atoms from some initial to final 
            fine-structure level, and for a list of impact energies (not yet).
    + struct ImpactIonizationScheme    
        ... to model the (direct) electron-impact ionization  (collision strength) of atoms from some initial to final 
            fine-structure level, and for a list of impact energies. It typically applies some empirical cross sections 
            (not yet).
    + struct PhotoAbsorptionScheme    
        ... to model photoabsortion spectra, including the direct and resonant contributions, i.e. the photoexcitation
            of an inner-shell electron with subsequent electron emission. Typical photoabsorption properties are the
            energy-dependent photoionization cross sections, photoabsorption spectra, PI plasma rate coefficients, and 
            several others.
    + struct PhotoExcitationScheme    
        ... to model the (prior) photo-excitation part of an overall photoabsorption process; it considers a set of 
            inner-shell excitations with regard to the list of reference configurations. This cascade scheme is mainly
            used to compute the resonant contributions to photoabsorption or the initial excitations for a subsequent
            decay cascade.
    + struct PhotoIonizationScheme    
        ... to model the photoionization a basic part of photoabsortion. Typical photoionization properties are the 
            energy-dependent partial and total photoionization cross sections for a range of photon energies
            and several others.
    + struct RadiativeRecombinationScheme  
        ... to model the radiative recombination (capture) of an electron up to a maximum free-electron energy as well as
            for a given list of shells (intoShells), into which the capture is considered; the energy of the emitted photons then
            refer to the levels of the reference configurations, and all cascade blocks are built only with the given intoShells.
    + struct StepwiseDecayScheme       
        ... to model a standard decay scheme in terms of radiative and non-radiative transitions by starting from 
            the levels of one or several initial multiplets. Typical decay properties are ion distributions as well as
            photon and electron spectra from such cascades.
"""
abstract type  AbstractCascadeScheme       end


"""
`struct  Cascade.DielectronicCaptureScheme  <:  Cascade.AbstractCascadeScheme`  
    ... to model just the (dielectronic) capture and the formation of doubly-excited levels up to a maximum excitation 
        energy and for a given list of subshells; the excitation energy refers to the lowest level of the reference 
        configurations, and the cascade blocks are built only by means of the given subshells.
        NOTE: The original DielectronicCaptureScheme --> DielectronicRecombinationScheme has been renamed in 
        August 2023 in order to enlarge the consistency of the notations and code !!

    + maxExcitationEnergy   ::Float64                 
        ... Maximum excitation energy [in a.u.] with regard to the reference configurations/levels that restrict the number 
            of excited configurations to be taken into accout. This maximum excitation energy has to be derived from the maximum 
            temperature for which DR coefficients need to be derived and is typically set to 5x T_e,max.
    + electronEnergyShift   ::Float64                 
        ... Energy shift for all resonance energies; this is realized by shifting the initial level energies by the negative amount.
            The shift is taken in the user-defined units.
    + NoExcitations         ::Int64                 
        ... (Maximum) Number of electron replacements in the doubly-excited configuration with regard to the initial 
            configurations/multiplets, apart from one additional electron due to the electron capture itself.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which excitations are to be considered.
    + excitationToShells  ::Array{Shell,1}    
        ... List of shells to which (core-shell) excitations are to be considered.
    + intoShells            ::Array{Shell,1}
        ... List of shells into which electrons are initially placed (captured).
"""
struct   DielectronicCaptureScheme  <:  Cascade.AbstractCascadeScheme
    maxExcitationEnergy     ::Float64   
    electronEnergyShift     ::Float64 
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    intoShells              ::Array{Shell,1}
end


"""
`Cascade.DielectronicCaptureScheme()`  ... constructor for an 'default' instance of a Cascade.DielectronicCaptureScheme.
"""
function DielectronicCaptureScheme()
    DielectronicCaptureScheme(1.0, 0., 1, Shell[], Shell[], Shell[], Shell[] )
end


# `Base.string(scheme::DielectronicCaptureScheme)`  ... provides a String notation for the variable scheme::DielectronicCaptureScheme.
function Base.string(scheme::DielectronicCaptureScheme)
    sa = "Dielectronic capture (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::DielectronicCaptureScheme)`  ... prepares a proper printout of the scheme::DielectronicCaptureScheme.
function Base.show(io::IO, scheme::DielectronicCaptureScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "maxExcitationEnergy:        $(scheme.maxExcitationEnergy)  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "NoExcitations:              $(scheme.NoExcitations)  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "intoShells:                 $(scheme.intoShells)  ")
end


"""
`struct  Cascade.DielectronicRecombinationScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe the dielectronic recombination of electrons for an atom in some initial state/configuration;
        for such a scheme, the doubly-excited configurations due to the electron capture are generated automatically due to
        given maximal numbers of the (into-) shells (nl) as well as the maximum displacement with regard to the initial configuration.
        An additional maximum excitation energy need to be provided due to the maximum temperatures for which DR plasma rate coefficients
        are to be determined, cf. Basics.convert().

    + multipoles            ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the radiative stabilization processes.
    + calcWithoutIntoShells ::Bool
        ... If true, the intoShells are not considered explicitly but deterined automatically from the given 
            initial configurations, the from shells as well as the total excitation energy. This features 
            should be set to true, whenever the intoShells are dificult to determine in advance.
    + maxIntoShell          ::Shell   
        ... maximum shell into which an electron is captured, if no intoShells are given explicitly (calcWithoutIntoShells=true).
            From this maxIntoShell, the maximum n_max and l_max is derived for generating the doubly-excited configurations
    + maxExcitationEnergy   ::Float64                 
        ... Maximum excitation energy [in a.u.] with regard to the reference configurations/levels that restrict the number 
            of excited configurations to be taken into accout. This maximum excitation energy has to be derived from the maximum 
            temperature for which DR coefficients need to be derived and is typically set to 5x T_e,max.
    + electronEnergyShift   ::Float64                 
        ... Energy shift for all resonance energies; this is realized by shifting the initial level energies by the negative amount.
            The shift is taken in the user-defined units.
    + minPhotonEnergy       ::Float64                 
        ... Minimum (mean) photon energy [in a.u.] in the radiative stabilization of the doubly-excited configurations; 
            If cascade blocks are separated be less than this energy, the radiative stabilization is neglected.
    + NoExcitations         ::Int64                 
        ... (Maximum) Number of electron replacements in the doubly-excited configuration with regard to the initial 
            configurations/multiplets, apart from one additional electron due to the electron capture itself.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which excitations are to be considered.
    + excitationToShells  ::Array{Shell,1}    
        ... List of shells to which (core-shell) excitations are to be considered.
    + intoShells            ::Array{Shell,1}
        ... List of shells into which electrons are initially placed (captured).
    + decayShells           ::Array{Shell,1}
        ... List of shells into which electrons the electrons can decay (apart from the core shells).
"""
struct   DielectronicRecombinationScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    calcWithoutIntoShells   ::Bool
    maxIntoShell            ::Shell       
    maxExcitationEnergy     ::Float64   
    electronEnergyShift     ::Float64 
    minPhotonEnergy         ::Float64                 
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    intoShells              ::Array{Shell,1}
    decayShells             ::Array{Shell,1}
end


"""
`Cascade.DielectronicRecombinationScheme()`  ... constructor for an 'default' instance of a Cascade.DielectronicRecombinationScheme.
"""
function DielectronicRecombinationScheme()
    DielectronicRecombinationScheme([E1], false, Shell(1,0), 1.0, 0., 0., 1, Shell[], Shell[], Shell[], Shell[] )
end


# `Base.string(scheme::DielectronicRecombinationScheme)`  ... provides a String notation for the variable scheme::DielectronicRecombinationScheme.
function Base.string(scheme::DielectronicRecombinationScheme)
    sa = "Dielectronic capture & stabilization (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::DielectronicRecombinationScheme)`  ... prepares a proper printout of the scheme::DielectronicRecombinationScheme.
function Base.show(io::IO, scheme::DielectronicRecombinationScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "calcWithoutIntoShells:      $(scheme.calcWithoutIntoShells)  ")
    println(io, "maxIntoShell:               $(scheme.maxIntoShell)  ")
    println(io, "maxExcitationEnergy:        $(scheme.maxExcitationEnergy)  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
    println(io, "NoExcitations:              $(scheme.NoExcitations)  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "intoShells:                 $(scheme.intoShells)  ")
    println(io, "decayShells:                $(scheme.decayShells)  ")
end


"""
`struct  Cascade.ElectronExcitationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... to compute electron excitation spectra including the direct (EIE) and resonant contributions, i.e. the dielectronic
            capture with subsequent re-autoionization. Typical electron-excitation properties are energy-dependent 
            EIE cross sections, effective collision strengths, EIE plasma rate coefficients, and others.

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + electronEnergies      ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
"""
struct   ElectronExcitationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    electronEnergies        ::Array{Float64,1}
end


"""
`Cascade.ElectronExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.ElectronExcitationScheme.
"""
function ElectronExcitationScheme()
    ElectronExcitationScheme([Eimex], Float64[] )
end


# `Base.string(scheme::ElectronExcitationScheme)`  ... provides a String notation for the variable scheme::ElectronExcitationScheme.
function Base.string(scheme::ElectronExcitationScheme)
    sa = "Electron-impact excitation (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ElectronExcitationScheme)`  ... prepares a proper printout of the scheme::ElectronExcitationScheme.
function Base.show(io::IO, scheme::ElectronExcitationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
    println(io, "processes:                   $(scheme.processes)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
end


"""
`struct  Cascade.ElectronIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... to compute electron ionization spectra including the direct (EII) and resonant contributions, i.e. the dielectronic
            capture with subsequent double-autoionization. For this double autoionization, a branching factor is estimated.
            Typical electron-ionization properties are energy-dependent EII cross sections, effective collision strengths, 
            EII plasma rate coefficients, and others (not yet).

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + electronEnergies      ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
"""
struct   ElectronIonizationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    electronEnergies        ::Array{Float64,1}
end


"""
`Cascade.ElectronIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.ElectronIonizationScheme.
"""
function ElectronIonizationScheme()
    ElectronIonizationScheme([Eimex], Float64[] )
end


# `Base.string(scheme::ElectronIonizationScheme)`  ... provides a String notation for the variable scheme::ElectronIonizationScheme.
function Base.string(scheme::ElectronIonizationScheme)
    sa = "Electron-impact ionization (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ElectronIonizationScheme)`  ... prepares a proper printout of the scheme::ElectronIonizationScheme.
function Base.show(io::IO, scheme::ElectronIonizationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
    println(io, "processes:                   $(scheme.processes)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
end


"""
`struct  Cascade.ExpansionOpacityScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe the expansion opacity of ions in some initial state/configuration; for this scheme,
        the excited (even- and odd-parity) configurations due to the photoabsorption and emission in a plasma are generated automatically 
        in terms of the chosen excitation scheme and a maximum photon (transition) energy that is taken into account.
        
    + multipoles            ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included for the radiative transitions in the plasma.
    + minPhotonEnergy       ::Float64                 
        ... Minimum photon (transition) energy [in a.u.] that are taken into account for all absorption lines; 
            this transition energy refers to the longest wavelength for which transition amplitudes are calculated.
    + maxPhotonEnergy       ::Float64                 
        ... Maximum photon (transition) energy [in a.u.] that are taken into account for all absorption lines; 
            this transition energy refers to the shortest wavelength for which the opacity is needed.
    + meanEnergyShift       ::Float64                 
        ... Energy shift for all excited configurations [in a.u.]; this allows to correct for missing correlation 
            contributions.
    + NoExcitations         ::Int64                 
        ... (Maximum) Number of electron replacements in the excited configuration with regard to the initial configurations/multiplets.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which excitations are to be considered.
    + excitationToShells    ::Array{Shell,1}    
        ... List of shells to which excitations are to be considered.
    + printTransitions      ::Bool      ... Print transition data for comparison, if true.
"""
struct   ExpansionOpacityScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    minPhotonEnergy         ::Float64 
    maxPhotonEnergy         ::Float64   
    meanEnergyShift         ::Float64                 
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    printTransitions        ::Bool
end


"""
`Cascade.ExpansionOpacityScheme()`  ... constructor for an 'default' instance of a Cascade.ExpansionOpacityScheme.
"""
function ExpansionOpacityScheme()
    ExpansionOpacityScheme([E1], 0., 1.0, 0., 1, Shell[], Shell[], false)
end


# `Base.string(scheme::ExpansionOpacityScheme)`  ... provides a String notation for the variable scheme::ExpansionOpacityScheme.
function Base.string(scheme::ExpansionOpacityScheme)
    sa = "Expansion opacity calculation (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ExpansionOpacityScheme)`  ... prepares a proper printout of the scheme::ExpansionOpacityScheme.
function Base.show(io::IO, scheme::ExpansionOpacityScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
    println(io, "maxPhotonEnergy:            $(scheme.maxPhotonEnergy)  ")
    println(io, "meanEnergyShift:            $(scheme.meanEnergyShift)  ")
    println(io, "NoExcitations:              $(scheme.NoExcitations)  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "printTransitions:           $(scheme.printTransitions)  ")
end


"""
`struct  Cascade.HollowIonScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe the formation and decay of a hollow ion, e.g. an electronic core configuration 
        into which one or several additional electrons are captured into (high) nl shells. Both the shell (lists) for 
        the initial capture (intoShells) and the subsequent decay (decayShells) need to be specified explicitly to 
        readily control the size of the computations.

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the decay scheme.  
    + multipoles            ::Array{EmMultipole,1}           
        ... Multipoles of the radiation field that are to be included into the radiative stabilization processes.
    + NoCapturedElectrons   ::Int64   
        ... Number of captured electrons, e.g. placed in the intoShells.
    + intoShells            ::Array{Shell,1}
        ... List of shells into which electrons are initially placed (captured).
    + decayShells           ::Array{Shell,1}
        ... List of shells into which electrons the electrons can decay (apart from the core shells).
"""
struct   HollowIonScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}    
    multipoles              ::Array{EmMultipole}  
    NoCapturedElectrons     ::Int64
    intoShells              ::Array{Shell,1}
    decayShells             ::Array{Shell,1}
end


"""
`Cascade.HollowIonScheme()`  ... constructor for an 'default' instance of a Cascade.HollowIonScheme.
"""
function HollowIonScheme()
    HollowIonScheme([Auger(), Radiative()], [E1], 0, Shell[], Shell[] )
end


# `Base.string(scheme::HollowIonScheme)`  ... provides a String notation for the variable scheme::HollowIonScheme.
function Base.string(scheme::HollowIonScheme)
    sa = "Hollow ion (scheme) due to processes:"
    return( sa )
end


# `Base.show(io::IO, scheme::HollowIonScheme)`  ... prepares a proper printout of the scheme::HollowIonScheme.
function Base.show(io::IO, scheme::HollowIonScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "processes:                  $(scheme.processes)  ")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "NoCapturedElectrons:        $(scheme.NoCapturedElectrons)  ")
    println(io, "intoShells:                 $(scheme.intoShells)  ")
    println(io, "decayShells:                $(scheme.decayShells)  ")
end


"""
`struct  Cascade.ImpactExcitationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... to compute the (direct) electron-impact excitation spectrum for a list of impact energies (not yet).

    + processes              ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + fromShells             ::Array{Shell,1}    
        ... List of shells from which impact-excitations are to be considered.
    + toShells               ::Array{Shell,1}    
        ... List of shells into which impact-excitations are to be considered, including possibly already occupied shells.
    + electronEnergies       ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
    + lValues                ::Array{Int64,1}
        ... Orbital angular momentum values of the free-electrons, for which partial waves are considered for the RR.
    + NoFreeElectronEnergies ::Int64             
        ... Number of free-electron energies that a chosen for a Gauss-Laguerre integration.
    + maxFreeElectronEnergy  ::Float64             
        ... Maximum free-electron energies [in a.u.] that restrict the energy of free-electron orbitals; this maximum energy has to 
            be derived from the maximum temperature for which RR plasma coefficients need to be obtained and is typically set to 
            about 5x T_e,max.
    + electronEnergyShift    ::Float64                 
        ... Energy shift for all bound-state energies relative to the levels from the reference configuration; this is realized by 
            shifting the initial level energies by the negative amount. The shift is taken in the user-defined units.
            
    Either a list of electronEnergies or the NoFreeElectronEnergies can be specified; the program terminates, if "non-zero"
    entries appears for these two subfields.
"""
struct   ImpactExcitationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    fromShells              ::Array{Shell,1}
    toShells                ::Array{Shell,1}
    electronEnergies        ::Array{Float64,1}
    lValues                 ::Array{Int64,1}
    NoFreeElectronEnergies  ::Int64 
    maxFreeElectronEnergy   ::Float64 
    electronEnergyShift     ::Float64    
end


"""
`Cascade.ImpactExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.ImpactExcitationScheme.
"""
function ImpactExcitationScheme()
    ImpactExcitationScheme([Eimex], Shell[], Shell[], Float64[], Int64[], 0, 0., 0. )
end


# `Base.string(scheme::ImpactExcitationScheme)`  ... provides a String notation for the variable scheme::ImpactExcitationScheme.
function Base.string(scheme::ImpactExcitationScheme)
    sa = "Direct electron-impact excitation (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ImpactExcitationScheme)`  ... prepares a proper printout of the scheme::ImpactExcitationScheme.
function Base.show(io::IO, scheme::ImpactExcitationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
    println(io, "processes:                   $(scheme.processes)  ")
    println(io, "fromShells:                  $(scheme.fromShells)  ")
    println(io, "toShells:                    $(scheme.toShells)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
    println(io, "lValues:                     $(scheme.lValues)  ")
    println(io, "NoFreeElectronEnergies:      $(scheme.NoFreeElectronEnergies)  ")
    println(io, "maxFreeElectronEnergy:       $(scheme.maxFreeElectronEnergy)  ")
    println(io, "electronEnergyShift:         $(scheme.electronEnergyShift)  ")
    
    if  scheme.NoFreeElectronEnergies != 0   &&   length(scheme.electronEnergies) != 0   
        error("Either electronEnergies  <OR>  NoFreeElectronEnergies can be specified explicitly.")
    end
end


"""
`struct  Cascade.ImpactIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... to compute the (direct) electron-impact excitation spectrum for a list of impact energies (not yet).

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + electronEnergies      ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
"""
struct   ImpactIonizationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    electronEnergies        ::Array{Float64,1}
end


"""
`Cascade.ImpactIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.ImpactIonizationScheme.
"""
function ImpactIonizationScheme()
    ImpactIonizationScheme([Eimex], Float64[] )
end


# `Base.string(scheme::ImpactIonizationScheme)`  ... provides a String notation for the variable scheme::ImpactIonizationScheme.
function Base.string(scheme::ImpactIonizationScheme)
    sa = "Direct electron-impact ionization (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ImpactIonizationScheme)`  ... prepares a proper printout of the scheme::ImpactIonizationScheme.
function Base.show(io::IO, scheme::ImpactIonizationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
    println(io, "processes:                   $(scheme.processes)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
end


"""
`struct  Cascade.PhotoAbsorptionScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe a photo-absorption calculation for an atom in some initial state/configuration
        and for a given range of photon energies, processes and multipoles, etc.

    + multipoles            ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the excitation/ionization processes.
    + photonEnergies        ::Array{Float64,1}
        ... List of photon energies (in user-selected units) for which absorption cross sections/spectra are to be
            calculated; this describes the list, distribution and resolution of energies. It is checked that either
            photonEnergies or electronEnergies are given only
    + electronEnergies       ::Array{Float64,1}
        ... List of electron energies (in user-selected units) for which absorption cross sections/spectra are to be
            calculated; this describes the list, distribution and resolution of energies.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which photo-excitations are to be considered.
    + excitationToShells    ::Array{Shell,1}    
        ... List of shells into which photo-excitations are to be considered, including possibly already occupied shells.
    + initialLevelSelection ::LevelSelection    
        ... Specifies the selected initial levels of some given initial-state configurations; these initial level numbers/
            symmetries always refer to the set of initial configurations.
    + lValues               ::Array{Int64,1}
        ... Orbital angular momentum values of the free-electrons, for which partial waves are considered for the PI.
    + calcDirect            ::Bool      ... True, if the direct contributions need to be calculated.                
    + calcResonant          ::Bool      ... True, if the resonant contributions need to be calculated.                
    + electronEnergyShift   ::Float64                 
        ... Energy shift for all bound-state energies relative to the levels from the reference configuration; this is realized by 
            shifting the initial level energies by the negative amount. The shift is taken in the user-defined units.
    + minCrossSection       ::Float64                 
        ... minimum cross section (in user-selected units) for which contributions are accounted for in the list of
            photoionization lines.
"""
struct   PhotoAbsorptionScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    photonEnergies          ::Array{Float64,1}                 
    electronEnergies        ::Array{Float64,1}                 
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    initialLevelSelection   ::LevelSelection 
    lValues                 ::Array{Int64,1}
    calcDirect              ::Bool               
    calcResonant            ::Bool               
    electronEnergyShift     ::Float64
    minCrossSection         ::Float64
end


"""
`Cascade.PhotoAbsorptionScheme()`  ... constructor for an 'default' instance of a Cascade.PhotoAbsorptionScheme.
"""
function PhotoAbsorptionScheme()
    PhotoAbsorptionScheme([E1], Float64[], Float64[], Shell[], Shell[], LevelSelection(false), [0], true, false, 0., 0.)
end


# `Base.string(scheme::PhotoAbsorptionScheme)`  ... provides a String notation for the variable scheme::PhotoAbsorptionScheme.
function Base.string(scheme::PhotoAbsorptionScheme)
    sa = "Photoabsorption (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::PhotoAbsorptionScheme)`  ... prepares a proper printout of the scheme::PhotoAbsorptionScheme.
function Base.show(io::IO, scheme::PhotoAbsorptionScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "photonEnergies:             $(scheme.photonEnergies)  ")
    println(io, "electronEnergies:           $(scheme.electronEnergies)  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "initialLevelSelection:      $(scheme.initialLevelSelection)  ")
    println(io, "lValues:                    $(scheme.lValues )  ")
    println(io, "calcDirect:                 $(scheme.calcDirect )  ")
    println(io, "calcResonant:               $(scheme.calcResonant )  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "minCrossSection:            $(scheme.minCrossSection)  ")
    #
    if  length(scheme.photonEnergies) > 0  &&  length(scheme.electronEnergies) > 0    
        error("Only photon or electron energies can be specified.")
    end
end


"""
`struct  Cascade.PhotoExcitationScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe a photo-excitation calculation for an atom in some initial state/configuration
        and for a given set of shell-excitations

    + multipoles            ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the excitation processes.
    + minPhotonEnergy       ::Float64                 
        ... Minimum photon energy [in a.u.] that restrict the number of excited configurations to be taken into accout.
    + maxPhotonEnergy       ::Float64                 
        ... Maximum photon energy [in a.u.] that restrict the number of excited configurations to be taken into accout.
    + NoExcitations         ::Int64                 
        ... (Maximum) Number of electron replacements with regard to the initial configurations/multiplets.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which photo-excitations are to be considered.
    + excitationToShells    ::Array{Shell,1}    
        ... List of shells into which photo-excitations are to be considered, including possibly already occupied shells.
    + initialLevelSelection ::LevelSelection    
        ... Specifies the selected initial levels of some given initial-state configurations; these initial level numbers/
            symmetries always refer to the set of initial configurations.
    + lValues               ::Array{Int64,1}
        ... Orbital angular momentum values of the free-electrons, for which partial waves are considered for the PI.
    + electronEnergyShift   ::Float64                 
        ... Energy shift for all bound-state energies relative to the levels from the reference configuration; this is realized by 
            shifting the initial level energies by the negative amount. The shift is taken in the user-defined units.
    + minCrossSection       ::Float64                 
        ... minimum cross section (in user-selected units) for which contributions are accounted for in the list of
            photoionization lines. This may seriously restrict the amount of data that is prepared for the subsequent simulation 
            of photoabsorption spectra.
"""
struct   PhotoExcitationScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    minPhotonEnergy         ::Float64                 
    maxPhotonEnergy         ::Float64                 
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    initialLevelSelection   ::LevelSelection 
    lValues                 ::Array{Int64,1}
    electronEnergyShift     ::Float64
    minCrossSection         ::Float64
end


"""
`Cascade.PhotoExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.PhotoExcitationScheme.
"""
function PhotoExcitationScheme()
    PhotoExcitationScheme([E1], 1.0, 100., 0, Shell[], Shell[], LevelSelection(false), [0], 0., 0. )
end


# `Base.string(scheme::PhotoExcitationScheme)`  ... provides a String notation for the variable scheme::PhotoExcitationScheme.
function Base.string(scheme::PhotoExcitationScheme)
    sa = "Photon-excitation (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::PhotoExcitationScheme)`  ... prepares a proper printout of the scheme::PhotoExcitationScheme.
function Base.show(io::IO, scheme::PhotoExcitationScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
    println(io, "maxPhotonEnergy:            $(scheme.maxPhotonEnergy)  ")
    println(io, "NoExcitations :             $(scheme.NoExcitations )  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "initialLevelSelection:      $(scheme.initialLevelSelection)  ")
    println(io, "lValues:                    $(scheme.lValues)  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "minCrossSection:            $(scheme.minCrossSection)  ")
end


"""
`struct  Cascade.PhotoIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe a photo-absorption calculation for an atom in some initial state/configuration
        and for a given range of photon energies and multipoles, etc.

    + multipoles            ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the excitation/ionization processes.
    + photonEnergies        ::Array{Float64,1}
        ... List of photon energies (in user-selected units) for which absorption cross sections/spectra are to be
            calculated; this describes the list, distribution and resolution of energies. It is checked that either
            photonEnergies or electronEnergies are given only.
    + electronEnergies       ::Array{Float64,1}
        ... List of electron energies (in user-selected units) for which absorption cross sections/spectra are to be
            calculated; this describes the list, distribution and resolution of energies.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which photo-excitations are to be considered.
    + excitationToShells    ::Array{Shell,1}    
        ... List of shells into which photo-excitations are to be considered, including possibly already occupied shells.
    + initialLevelSelection ::LevelSelection    
        ... Specifies the selected initial levels of some given initial-state configurations; these initial level numbers/
            symmetries always refer to the set of initial configurations.
    + lValues               ::Array{Int64,1}
        ... Orbital angular momentum values of the free-electrons, for which partial waves are considered for the PI.
    + electronEnergyShift   ::Float64                 
        ... Energy shift for all bound-state energies relative to the levels from the reference configuration; this is realized by 
            shifting the initial level energies by the negative amount. The shift is taken in the user-defined units.
    + minCrossSection       ::Float64                 
        ... minimum cross section (in user-selected units) for which contributions are accounted for in the list of
            photoionization lines.
"""
struct   PhotoIonizationScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    photonEnergies          ::Array{Float64,1}                 
    electronEnergies        ::Array{Float64,1}                 
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    initialLevelSelection   ::LevelSelection 
    lValues                 ::Array{Int64,1}
    electronEnergyShift     ::Float64
    minCrossSection         ::Float64
end


"""
`Cascade.PhotoIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.PhotoIonizationScheme.
"""
function PhotoIonizationScheme()
    PhotoIonizationScheme([E1], Float64[], Float64[], Shell[], Shell[], LevelSelection(false), [0], 0., 0.)
end


# `Base.string(scheme::PhotoIonizationScheme)`  ... provides a String notation for the variable scheme::PhotoIonizationScheme.
function Base.string(scheme::PhotoIonizationScheme)
    sa = "Photoionization (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::PhotoIonizationScheme)`  ... prepares a proper printout of the scheme::PhotoIonizationScheme.
function Base.show(io::IO, scheme::PhotoIonizationScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "photonEnergies:             $(scheme.photonEnergies)  ")
    println(io, "electronEnergies:           $(scheme.electronEnergies)  ")
    println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    println(io, "initialLevelSelection:      $(scheme.initialLevelSelection)  ")
    println(io, "lValues:                    $(scheme.lValues )  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "minCrossSection:            $(scheme.minCrossSection)  ")
    #
    if  length(photonEnergies) > 0  &&  length(electronEnergies) > 0    
        error("Only photon or electron energies can be specified.")
    end
end


#==
"""
`struct  Cascade.PhotoIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to represent (and generate) a mean-field orbital basis.

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + maxPhotoElectrons     ::Int64                 
        ... (Maximum) Number of photo-electrons in which the initial-and (photo-ionized) final-state configurations can differ 
            from each other; this affects the number of ionized multiplets in the cascade.
    + photonEnergies        ::Array{Float64,1}      
        ... List of photon energies for which this photo-ionization scheme is to be calculated.
"""
struct   PhotoIonizationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    maxPhotoElectrons       ::Int64
    photonEnergies          ::Array{Float64,1}
end


"""
`Cascade.PhotoIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.PhotoIonizationScheme.
"""
function PhotoIonizationScheme()
    PhotoIonizationScheme([PhotoIonization], 0, Float64[] )
end


# `Base.string(scheme::PhotoIonizationScheme)`  ... provides a String notation for the variable scheme::PhotoIonizationScheme.
function Base.string(scheme::PhotoIonizationScheme)
    sa = "Photo-ionization (scheme) as prior part of an atomic cascade:"
    return( sa )
end


# `Base.show(io::IO, scheme::PhotoIonizationScheme)`  ... prepares a proper printout of the scheme::PhotoIonizationScheme.
function Base.show(io::IO, scheme::PhotoIonizationScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "processes:                  $(scheme.processes)  ")
    println(io, "maxPhotoElectrons:          $(scheme.maxPhotoElectrons)  ")
    println(io, "photonEnergies:             $(scheme.photonEnergies)  ")
end
==#


"""
`struct  Cascade.RadiativeRecombinationScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to define and describe the radiative recombination of electrons for an atom in some initial state/configuration;
        for such a scheme, the configurations due to the electron capture are generated automatically due to given maximal numbers 
        of the (into-) shells (nl) as well as the maximum displacement with regard to the initial configuration.
        An additional maximum excitation energy need to be provided due to the maximum temperatures for which RR plasma rate coefficients
        are to be determined, cf. Basics.convert().

    + multipoles             ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the radiative stabilization processes.
    + lValues                ::Array{Int64,1}
        ... Orbital angular momentum values of the free-electrons, for which partial waves are considered for the RR.
    + NoFreeElectronEnergies ::Int64             
        ... Number of free-electron energies that a chosen for a Gauss-Laguerre integration.
    + maxFreeElectronEnergy  ::Float64             
        ... Maximum free-electron energies [in a.u.] that restrict the energy of free-electron orbitals; this maximum energy has to 
            be derived from the maximum temperature for which RR plasma coefficients need to be obtained and is typically set to 
            about 5x T_e,max.
    + electronEnergyShift    ::Float64                 
        ... Energy shift for all bound-state energies relative to the levels from the reference configuration; this is realized by 
            shifting the initial level energies by the negative amount. The shift is taken in the user-defined units.
    + minPhotonEnergy       ::Float64                 
        ... Minimum (mean) photon energy [in a.u.] for which the radiative decay is taken into account.
    + intoShells            ::Array{Shell,1}
        ... List of shells into which electrons are initially placed (captured).
"""
struct   RadiativeRecombinationScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    lValues                 ::Array{Int64,1}
    NoFreeElectronEnergies  ::Int64  
    maxFreeElectronEnergy   ::Float64
    electronEnergyShift     ::Float64 
    minPhotonEnergy         ::Float64                 
    intoShells              ::Array{Shell,1}
end


"""
`Cascade.RadiativeRecombinationScheme()`  ... constructor for an 'default' instance of a Cascade.RadiativeRecombinationScheme.
"""
function RadiativeRecombinationScheme()
    RadiativeRecombinationScheme([E1], Int64[], 0, 0., 0., 0., Shell[] )
end


# `Base.string(scheme::RadiativeRecombinationScheme)`  ... provides a String notation for the variable scheme::RadiativeRecombinationScheme.
function Base.string(scheme::RadiativeRecombinationScheme)
    sa = "Radiative recombination (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::RadiativeRecombinationScheme)`  ... prepares a proper printout of the scheme::RadiativeRecombinationScheme.
function Base.show(io::IO, scheme::RadiativeRecombinationScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "multipoles:                 $(scheme.multipoles)  ")
    println(io, "lValues:                    $(scheme.lValues)  ")
    println(io, "NoFreeElectronEnergies:     $(scheme.NoFreeElectronEnergies)  ")
    println(io, "maxFreeElectronEnergy:      $(scheme.maxFreeElectronEnergy)  ")
    println(io, "electronEnergyShift:        $(scheme.electronEnergyShift)  ")
    println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
    println(io, "intoShells:                 $(scheme.intoShells)  ")
end


"""
`struct  Cascade.StepwiseDecayScheme  <:  Cascade.AbstractCascadeScheme`  
    ... a struct to represent (and generate) a mean-field orbital basis.

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + maxElectronLoss       ::Int64             
        ... (Maximum) Number of electrons in which the initial- and final-state configurations can differ from each other; 
            this also determines the maximal steps of any particular decay path.
    + chargeStateShifts     ::Dict{Int64,Float64} 
        ... (N => en) total energy shifts of all levels with N electrons; these shifts [in a.u.] help open/close decay 
            channels by simply shifting the total energies of all levels.
    + NoShakeDisplacements  ::Int64             
        ... Maximum number of electron displacements due to shake-up  or shake-down processes in any individual step of cascade.
    + decayShells           ::Array{Shell,1}        ... List of shells that may occur during the decay.
    + shakeFromShells       ::Array{Shell,1}        ... List of shells from which shake transitions may occur.
    + shakeToShells         ::Array{Shell,1}        ... List of shells into which shake transitions may occur.
"""
struct   StepwiseDecayScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    maxElectronLoss         ::Int64
    chargeStateShifts       ::Dict{Int64,Float64}
    NoShakeDisplacements    ::Int64
    decayShells             ::Array{Shell,1}
    shakeFromShells         ::Array{Shell,1}
    shakeToShells           ::Array{Shell,1}
end


"""
`Cascade.StepwiseDecayScheme()`  ... constructor for an 'default' instance of a Cascade.StepwiseDecayScheme.
"""
function StepwiseDecayScheme()
    StepwiseDecayScheme([Radiative()], 0, Dict{Int64,Float64}(), 0, Shell[], Shell[], Shell[] )
end


# `Base.string(scheme::StepwiseDecayScheme)`  ... provides a String notation for the variable scheme::StepwiseDecayScheme.
function Base.string(scheme::StepwiseDecayScheme)
    sa = "Stepwise decay (scheme) of an atomic cascade with:"
    return( sa )
end


# `Base.show(io::IO, scheme::StepwiseDecayScheme)`  ... prepares a proper printout of the scheme::StepwiseDecayScheme.
function Base.show(io::IO, scheme::StepwiseDecayScheme)
    sa = Base.string(scheme);                print(io, sa, "\n")
    println(io, "processes:                  $(scheme.processes)  ")
    println(io, "maxElectronLoss:            $(scheme.maxElectronLoss)  ")
    println(io, "chargeStateShifts:          $(scheme.chargeStateShifts)  ")
    println(io, "NoShakeDisplacements:       $(scheme.NoShakeDisplacements)  ")
    println(io, "decayShells:                $(scheme.decayShells)  ")
    println(io, "shakeFromShells:            $(scheme.shakeFromShells)  ")
    println(io, "shakeToShells:              $(scheme.shakeToShells)  ")
end

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
`abstract type Cascade.AbstractCascadeApproach` 
    ... defines an abstract and a number of singleton types for the computational approach/model that is applied in order to 
        generate and evaluate all many-electron amplitudes of a given cascade.

    + struct AverageSCA         
    ... all levels in the cascade are described in single-configuration and single-CSF approximation; this (rather crude) approach 
        neglects all configuration-interactions and also applies just a single set of one-electron orbitals (from the least-ionized charge
        state) for all considered charge states.
        
    + struct SCA                
    ... all levels in the cascade are described in single-configuration approximation but with 'mixtures' within the configuration;
        an individual mean-field is generated for each charge state and all continuum orbitals are generated for the correct transition
        energy in the field of the remaining ion. Moreover, all the fine-structure transitions are calculated individually.
"""
abstract type  AbstractCascadeApproach                   end
struct         AverageSCA  <:  AbstractCascadeApproach   end
struct         SCA         <:  AbstractCascadeApproach   end
struct         UserMCA     <:  AbstractCascadeApproach   end



"""
`struct  Cascade.Block`  
    ... defines a type for an individual block of configurations that are treatet together within the cascade. Such an block is given 
        by a list of configurations that may occur as initial- and/or final-state configurations in some step of the canscade and that 
        give rise to a common multiplet in order to allow for configuration interactions but to avoid 'double counting' of individual 
        levels in the cascade.

    + NoElectrons     ::Int64                     ... Number of electrons in this block.
    + confs           ::Array{Configuration,1}    ... List of one or several configurations that define the multiplet.
    + hasMultiplet    ::Bool                      
        ... true if the (level representation in the) multiplet has already been computed and false otherwise.
    + multiplet       ::Multiplet                 ... Multiplet of the this block.
"""
struct  Block
    NoElectrons       ::Int64
    confs             ::Array{Configuration,1} 
    hasMultiplet      ::Bool
    multiplet         ::Multiplet  
end 


"""
`Cascade.Block()`  ... constructor for an 'empty' instance of a Cascade.Block.
"""
function Block()
    Block( )
end


# `Base.show(io::IO, block::Cascade.Block)`  ... prepares a proper printout of the variable block::Cascade.Block.
function Base.show(io::IO, block::Cascade.Block) 
    println(io, "NoElectrons:        $(block.NoElectrons)  ")
    println(io, "confs :             $(block.confs )  ")
    println(io, "hasMultiplet:       $(block.hasMultiplet)  ")
    println(io, "multiplet:          $(block.multiplet)  ")
end


"""
`struct  Cascade.Step`  
    ... defines a type for an individual step of an excitation and/or decay cascade. Such a step is determined by the two lists 
        of initial- and final-state configuration as well as by the atomic process, such as Auger, PhotoEmission, or others, which related
        the initial- and final-state levels to each other. Since the (lists of) initial- and final-state configurations treated (each)
        by a single multiplet (for parities and total angular momenta), a cascade step supports full configuration interaction within
        the multiplet but also help  avoid 'double counting' of individual levels. Indeed, each electron configuration may occur only in
        one cascade block. In contrast, each list of initial- and final-state multiplets (cascade blocks) can occur in quite different 
        steps due to the considered processes and parallel decay pathes in a cascade.

    + process          ::JBasics.AbstractProcess   ... Atomic process that 'acts' in this step of the cascade.
    + settings         ::Union{PhotoEmission.Settings, AutoIonization.Settings, PhotoIonization.Settings, PhotoExcitation.Settings,
                                DielectronicCapture.Settings}        
                                                    ... Settings for this step of the cascade.
    + initialConfigs   ::Array{Configuration,1}    ... List of one or several configurations that define the initial-state multiplet.
    + finalConfigs     ::Array{Configuration,1}    ... List of one or several configurations that define the final-state multiplet.
    + initialMultiplet ::Multiplet                 ... Multiplet of the initial-state levels of this step of the cascade.
    + finalMultiplet   ::Multiplet                 ... Multiplet of the final-state levels of this step of the cascade.
"""
struct  Step
    process            ::Basics.AbstractProcess
    settings           ::Union{PhotoEmission.Settings, AutoIonization.Settings, PhotoIonization.Settings, PhotoExcitation.Settings, 
                                PhotoRecombination.Settings}
    initialConfigs     ::Array{Configuration,1}
    finalConfigs       ::Array{Configuration,1}
    initialMultiplet   ::Multiplet
    finalMultiplet     ::Multiplet
end 


"""
`Cascade.Step()`  ... constructor for an 'empty' instance of a Cascade.Step.
"""
function Step()
    Step( Basics.NoProcess, PhotoEmission.Settings, Configuration[], Configuration[], Multiplet(), Multiplet())
end


# `Base.show(io::IO, step::Cascade.Step)`  ... prepares a proper printout of the variable step::Cascade.Step.
function Base.show(io::IO, step::Cascade.Step) 
    println(io, "process:                $(step.process)  ")
    println(io, "settings:               $(step.settings)  ")
    println(io, "initialConfigs:         $(step.initialConfigs)  ")
    println(io, "finalConfigs:           $(step.finalConfigs)  ")
    println(io, "initialMultiplet :      $(step.initialMultiplet )  ")
    println(io, "finalMultiplet:         $(step.finalMultiplet)  ")
end


"""
`struct  Cascade.Computation`  
    ... defines a type for a cascade computation, i.e. for the computation of a whole photon excitation, photon ionization and/or 
        decay cascade. The -- input and control -- data from this computation can be modified, adapted and refined to the practical needs, 
        and before the actual computations are carried out explictly. Initially, this struct just contains the physical meta-data about the 
        cascade, that is to be calculated, but a new instance of the same Cascade.Computation gets later enlarged in course of the 
        computation in order to keep also wave functions, level multiplets, etc.

    + name               ::String                          ... A name for the cascade
    + nuclearModel       ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
    + grid               ::Radial.Grid                     ... The radial grid to be used for the computation.
    + asfSettings        ::AsfSettings                     ... Provides the settings for the SCF process.
    + scheme             ::Cascade.AbstractCascadeScheme   ... Scheme of the atomic cascade (photoionization, decay, ...)
    + approach           ::Cascade.AbstractCascadeApproach 
        ... Computational approach/model that is applied to generate and evaluate the cascade; possible approaches are: 
            {AverageSCA(), SCA(), ...}
    + initialConfs       ::Array{Configuration,1}          
        ... List of one or several configurations that contain the level(s) from which the cascade starts.
    + initialMultiplets  ::Array{Multiplet,1}              
        ... List of one or several (initial) multiplets; either initialConfs 'xor' initialMultiplets can  be specified 
            for a given cascade computation.
"""
struct  Computation
    name                 ::String
    nuclearModel         ::Nuclear.Model
    grid                 ::Radial.Grid
    asfSettings          ::AsfSettings
    scheme               ::Cascade.AbstractCascadeScheme 
    approach             ::Cascade.AbstractCascadeApproach
    initialConfigs       ::Array{Configuration,1}
    initialMultiplets    ::Array{Multiplet,1} 
end 


"""
`Cascade.Computation()`  ... constructor for an 'default' instance of a Cascade.Computation.
"""
function Computation()
    Computation("Default cascade computation",  Nuclear.Model(10.), Radial.Grid(), AsfSettings(), Cascade.StepwiseDecayScheme(), 
                Cascade.AverageSCA(), Configuration[], Multiplet[] )
end


"""
`Cascade.Computation(comp::Cascade.Computation;`
    
            name=..,               nuclearModel=..,             grid=..,              asfSettings=..,     
            scheme=..,             approach=..,                 initialConfigs=..,    initialMultiplets=..)
            
    ... constructor for re-defining the computation::Cascade.Computation.
"""
function Computation(comp::Cascade.Computation;                              
    name::Union{Nothing,String}=nothing,                                  nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
    grid::Union{Nothing,Radial.Grid}=nothing,                             asfSettings::Union{Nothing,AsfSettings}=nothing,  
    scheme::Union{Nothing,Cascade.AbstractCascadeScheme}=nothing,         approach::Union{Nothing,Cascade.AbstractCascadeApproach}=nothing, 
    initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,        initialMultiplets::Union{Nothing,Array{Multiplet,1}}=nothing)
    
    if  initialConfigs != nothing   &&   initialMultiplets != nothing   
        error("Only initialConfigs=..  'xor'  initialMultiplets=..  can be specified for a Cascade.Computation().")
    end

    if  name                 == nothing   namex                 = comp.name                    else  namex = name                              end 
    if  nuclearModel         == nothing   nuclearModelx         = comp.nuclearModel            else  nuclearModelx = nuclearModel              end 
    if  grid                 == nothing   gridx                 = comp.grid                    else  gridx = grid                              end 
    if  asfSettings          == nothing   asfSettingsx          = comp.asfSettings             else  asfSettingsx = asfSettings                end 
    if  scheme               == nothing   schemex               = comp.scheme                  else  schemex = scheme                          end 
    if  approach             == nothing   approachx             = comp.approach                else  approachx = approach                      end 
    if  initialConfigs       == nothing   initialConfigsx       = comp.initialConfigs          else  initialConfigsx = initialConfigs          end 
    if  initialMultiplets    == nothing   initialMultipletsx    = comp.initialMultiplets       else  initialMultipletsx = initialMultiplets    end
    
    Computation(namex, nuclearModelx, gridx, asfSettingsx, schemex, approachx, initialConfigsx, initialMultipletsx)
end


# `Base.string(comp::Cascade.Computation)`  ... provides a String notation for the variable comp::Cascade.Computation.
function Base.string(comp::Cascade.Computation)
    if        typeof(comp.scheme) == Cascade.StepwiseDecayScheme              sb = "stepwise decay scheme"
    elseif    typeof(comp.scheme) == Cascade.DielectronicCaptureScheme        sb = "(di-) electronic capture scheme"
    elseif    typeof(comp.scheme) == Cascade.DielectronicRecombinationScheme  sb = "dielectronic recombination scheme"
    elseif    typeof(comp.scheme) == Cascade.ElectronExcitationScheme         sb = "electron-excitation scheme"
    elseif    typeof(comp.scheme) == Cascade.ElectronIonizationScheme         sb = "electron-ionization scheme"
    elseif    typeof(comp.scheme) == Cascade.ExpansionOpacityScheme           sb = "expansion opacity scheme"
    elseif    typeof(comp.scheme) == Cascade.HollowIonScheme                  sb = "hollow ion scheme"
    elseif    typeof(comp.scheme) == Cascade.ImpactExcitationScheme           sb = "electron-impact excitation scheme"
    elseif    typeof(comp.scheme) == Cascade.ImpactIonizationScheme           sb = "electron-impact ionization scheme"
    elseif    typeof(comp.scheme) == Cascade.PhotoAbsorptionScheme            sb = "photoabsorption scheme"
    elseif    typeof(comp.scheme) == Cascade.PhotoExcitationScheme            sb = "photo-excitation scheme"
    elseif    typeof(comp.scheme) == Cascade.PhotoIonizationScheme            sb = "photoionization scheme"
    elseif    typeof(comp.scheme) == Cascade.RadiativeRecombinationScheme     sb = "radiative recombination capture scheme"
    else      error("unknown typeof(comp.scheme)")
    end
    
    sa = "Cascade computation   $(comp.name)  for a $sb,  in $(comp.approach) approach as well as "
    sa = sa * "for Z = $(comp.nuclearModel.Z) and initial configurations: \n "
    if  comp.initialConfigs    != Configuration[]  for  config  in  comp.initialConfigs     sa = sa * string(config) * ",  "     end     end
    if  comp.initialMultiplets != Multiplet[]      
        for  mp      in  comp.initialMultiplets  sa = sa * mp.name * " with " * string(length(mp.levels)) * " levels,  "         end     end
    return( sa )
end


# `Base.show(io::IO, comp::Cascade.Computation)`  ... prepares a proper printout comp::Cascade.Computation.
function Base.show(io::IO, comp::Cascade.Computation)
    sa = Base.string(comp)
    sa = sa * "\n ... in addition, the following parameters/settings are defined: ";       print(io, sa, "\n")
    println(io, "> cascade scheme:           $(comp.scheme)  ")
    println(io, "> nuclearModel:             $(comp.nuclearModel)  ")
    println(io, "> grid:                     $(comp.grid)  ")
    println(io, "> asfSettings:            \n$(comp.asfSettings) ")
end

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

include("module-Cascade-inc-simulations-structs.jl")
include("module-Cascade-inc-computations.jl")
include("module-Cascade-inc-dielectronic-recombination.jl")
include("module-Cascade-inc-electron-excitation.jl")
include("module-Cascade-inc-electron-ionization.jl")
include("module-Cascade-inc-expansion-opacity.jl")
include("module-Cascade-inc-hollow-ion.jl")
include("module-Cascade-inc-impact-excitation.jl")
include("module-Cascade-inc-impact-ionization.jl")
include("module-Cascade-inc-photoabsorption.jl")
include("module-Cascade-inc-photoexcitation.jl")
include("module-Cascade-inc-photoionization.jl")
include("module-Cascade-inc-photorecombination.jl")
include("module-Cascade-inc-stepwise-decay.jl")
include("module-Cascade-inc-simulations.jl")


end # module


