
"""
`module  JAC.Cascade
... a submodel of JAC that contains all data types and methods to set-up and process (simple) 
    atomic cascade computations of various kind.
"""
module Cascade


using Dates, JLD2, Printf, FastGaussQuadrature, Distributed, ProgressMeter,
        ..AngularMomentum, ..AtomicState, ..AutoIonization, ..Basics, ..Bsplines, ..Continuum, ..Defaults,
        ..DecayYield, ..DielectronicRecombination, ..ElectronCapture, ..Empirical, ..Hamiltonian, ..ImpactExcitation, ..Radial, ..ManyElectron, ..Nuclear, 
        ..PeriodicTable, ..PhotoEmission, ..PhotoExcitation, ..PhotoIonization, ..PhotoRecombination, ..ResonantImpactIonization,
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
`Cascade.DielectronicCaptureScheme()`  ... constructor for a 'default' instance of a Cascade.DielectronicCaptureScheme.
"""
function DielectronicCaptureScheme()
    DielectronicCaptureScheme(0., 0., 0, Shell[], Shell[], Shell[] )
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
        ... List of the atomic processes that are supported and should be included into the cascade; Basics.ImpactExc()
            selects the DIRECT channel and Basics.ImpactExcAuto() the RESONANT one.
    + fromShells            ::Array{Shell,1}
        ... List of shells from which the impact excitation proceeds.
    + toShells              ::Array{Shell,1}
        ... List of shells into which the impact excitation proceeds.
    + electronEnergies      ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
    + lValues               ::Array{Int64,1}
        ... Orbital angular momenta of the incoming/outgoing free electron to be included.
    + NoFreeElectronEnergies ::Int64
        ... Number of free-electron energies of the internal (Gauss-Legendre) grid.
    + maxFreeElectronEnergy ::Float64
        ... Maximum free-electron energy [in a.u.] of that grid.
    + electronEnergyShift   ::Float64
        ... Energy shift of all electron energies, in the user-selected units.

        These are the same configuration data that a Cascade.ImpactExcitationScheme carries, because the direct channel IS
        an impact excitation and delegates to that scheme. Until 16-Aug-2026 this struct held only `processes` and
        `electronEnergies`, which left the direct channel with no shells to excite between and therefore unusable.
"""
struct   ElectronExcitationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    fromShells              ::Array{Shell,1}
    toShells                ::Array{Shell,1}
    electronEnergies        ::Array{Float64,1}
    lValues                 ::Array{Int64,1}
    NoFreeElectronEnergies  ::Int64
    maxFreeElectronEnergy   ::Float64
    electronEnergyShift     ::Float64
    maxExcitationEnergy     ::Float64
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    intoShells              ::Array{Shell,1}
end


"""
`Cascade.ElectronExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.ElectronExcitationScheme.
"""
function ElectronExcitationScheme()
    ElectronExcitationScheme([ImpactExc()], Shell[], Shell[], Float64[], Int64[], 0, 0., 0., 0., 0, Shell[], Shell[], Shell[] )
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
    println(io, "fromShells:                  $(scheme.fromShells)  ")
    println(io, "toShells:                    $(scheme.toShells)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
    println(io, "lValues:                     $(scheme.lValues)  ")
    println(io, "NoFreeElectronEnergies:      $(scheme.NoFreeElectronEnergies)  ")
    println(io, "maxFreeElectronEnergy:       $(scheme.maxFreeElectronEnergy)  ")
    println(io, "electronEnergyShift:         $(scheme.electronEnergyShift)  ")
end


"""
`struct  Cascade.ElectronIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... to compute the INDIRECT contribution to electron-impact ionization, i.e. excitation-autoionization (EA):
            a free electron excites an inner-shell electron into a level that lies above the ionization threshold of
            the ion, and that level then autoionizes.  The net effect is ionization, and for many ions -- notably
            along the Li-, Na- and Mg-like sequences -- EA rivals or exceeds the direct channel.

            The scheme is built exactly like Cascade.DielectronicRecombinationScheme, with which it is the mirror
            image: DR captures an electron and stabilizes radiatively, EA excites an electron and stabilizes by
            emitting one.  Both therefore generate two kinds of Cascade.Step,

                DR :  Auger (capture, by detailed balance)   +  Radiative (stabilization)
                EA :  ImpactExc (excitation)                 +  Auger (autoionization)

            A cascade computation returns the ImpactExcitation.Line's and AutoIonization.Line's; combining them
            into an EA cross section or plasma rate coefficient is the task of a Cascade.Simulation.

            NOT INCLUDED, and deliberately so:
              + the DIRECT electron-impact ionization channel.  Cascade.ImpactIonizationScheme is reserved for it
                and is not implemented; see the note there.
              + the two RESONANT-ELECTRON-CAPTURE channels, in which the incident electron is CAPTURED into a
                doubly-excited resonance which then sheds two electrons, either one after the other
                (resonant-electron-capture-with-sequential-double-autoionization, the REDA of the older
                literature) or both at once (...-with-simultaneous-double-autoionization, READI).  Together with
                the impact-excitation channel above they make up the resonant part of the ionization; only the
                impact-excitation one is available here.

                THE OLD NOTE HERE SAID THESE "require the electron-capture scheme that is still outstanding".
                That has not been true since 16-Aug-2026: Cascade.perform(::DielectronicCaptureScheme) exists and
                runs, and Cascade.perform(::ElectronExcitationScheme) already delegates to it for its own resonant
                channel.  What is still missing is the second Auger generation and the strengths, not the capture.

    + electronEnergies      ::Array{Float64,1}                
        ... List of impact energies of the incoming electron, in the user-selected units.  The collision strengths
            are computed AT these energies, so three or more well-spread values are needed if a rate coefficient
            is to be formed from them afterwards.
    + excitationFromShells  ::Array{Shell,1}
        ... List of (inner) shells out of which the electron is excited.
    + excitationToShells    ::Array{Shell,1}
        ... List of shells into which it is excited.  Only those excited configurations are kept whose mean energy
            lies ABOVE the ionization threshold of the ion, since only those can autoionize.
    + lValues               ::Array{Int64,1}
        ... Orbital angular momenta of the free electron, i.e. the partial waves summed over in the excitation.
            maxKappa is taken as maximum(lValues)+1.  A truncated sum does not merely lower the cross section, it
            can invert its energy dependence, and it does so silently; read the `convergence` column.
    + NoExcitations         ::Int64
        ... Number of displaced electrons for the generation of the excited configurations; 1 in nearly all cases.
    + electronEnergyShift   ::Float64
        ... Energy shift for all bound-state energies relative to the levels of the reference configuration, taken
            in the user-defined units.
    + processes             ::Array{Basics.AbstractProcess,1}
        ... The channels to be computed, selected exactly as Cascade.ElectronExcitationScheme selects its own.  Any
            combination may be given, and the cost of the selection varies by three orders of magnitude:

              Basics.ImpactExcAuto()                            impact-excitation with subsequent autoionization
                  ... the inner-shell excitation channel described above.  EXPENSIVE: it needs a partial-wave sum
                      over the free electron at every impact energy; ~2800 s for the Li-like carbon of example-Fi.jl.
              ResonantImpactIonization.SequentialAuger()        resonant-electron-capture-with-sequential-
                                                                double-autoionization  (REDA in older literature)
                  ... CHEAP, seconds rather than hours: a capture is the time reverse of an Auger and needs no
                      partial-wave sum over impact energies at all.
              ResonantImpactIonization.SimultaneousAuger()      resonant-electron-capture-with-simultaneous-
                                                                double-autoionization  (READI)
                  ... rides on the resonances of the sequential channel and adds almost nothing to the cost.

    + intoShells            ::Array{Shell,1}
        ... Shells into which the incident electron is CAPTURED, for the two resonant channels only; ignored by the
            impact-excitation channel.  The doubly-excited resonances are built by capturing into these shells while
            exciting from excitationFromShells into excitationToShells.
    + dblAugerOverride      ::Float64
        ... Optional replacement for the computed double-Auger probability of the simultaneous channel.  0. means
            "not set", in which case the shake-off estimate of ResonantImpactIonization is used; a positive value
            replaces it and is reported as having done so.

    APPROXIMATIONS, so that it is clear what a computation of this scheme actually delivers:

    + the DIRECT channel is absent
        ... Cascade.ImpactIonizationScheme is reserved for it and is not implemented, so a total ionization cross
            section CANNOT be obtained from this scheme alone.  What it gives is the indirect and resonant parts.
    + single-CSF blocks, no configuration mixing
        ... in Cascade.AverageSCA() every block is generated from single-CSF levels with orbitals from a
            Dirac-Fock-Slater potential, computed independently per block.  Level energies are therefore of
            cascade quality, not of spectroscopic quality.
    + only the Coulomb interaction in the autoionization
        ... the Breit interaction is not included in any Auger rate here.
    + a partial-wave truncation that is silent
        ... lValues fixes maxKappa for BOTH the impact excitation and the Auger rates.  Too small a value does not
            merely lower a cross section, it can invert its energy dependence; example-Fi.jl branch b measures that
            crossover explicitly.  The capture rates of the resonant channels inherit the same truncation.
    + isolated resonances, for the resonant channels
        ... resonance strengths are ADDED, with no interference between resonances and no overlap of their widths.
            Sound while the resonances are narrow against their spacing; it degrades for high capture shells.
    + branchings that sum to one BY CONSTRUCTION
        ... every branching ratio is formed from the total rates of the decay steps the cascade actually generated.
            A decay route left out of the configuration lists is absent from the denominator too, so the branchings
            still sum to one.  That sum checks the arithmetic; it can never reveal a missing channel.
    + the simultaneous channel is an ESTIMATE, not a computed rate
        ... its double-Auger width comes from shake-off in the sudden approximation, whose error has BOTH SIGNS.
            See ResonantImpactIonization.shakeProbability before quoting anything from it.
"""
struct   ElectronIonizationScheme  <:  Cascade.AbstractCascadeScheme
    electronEnergies        ::Array{Float64,1}
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    lValues                 ::Array{Int64,1}
    NoExcitations           ::Int64
    electronEnergyShift     ::Float64
    processes               ::Array{Basics.AbstractProcess,1}
    intoShells              ::Array{Shell,1}
    dblAugerOverride        ::Float64
end


"""
`Cascade.ElectronIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.ElectronIonizationScheme.
"""
function ElectronIonizationScheme()
    ElectronIonizationScheme(Float64[], Shell[], Shell[], Int64[], 1, 0., Basics.AbstractProcess[Basics.ImpactExcAuto()], Shell[], 0. )
end


"""
`Cascade.ElectronIonizationScheme(electronEnergies::Array{Float64,1}, excitationFromShells::Array{Shell,1},`
                                  `excitationToShells::Array{Shell,1}, lValues::Array{Int64,1}, NoExcitations::Int64,`
                                  `electronEnergyShift::Float64)`
    ... constructor for the impact-excitation channel alone, i.e. the six arguments this scheme had before the two
        resonant channels were added; processes defaults to [Basics.ImpactExcAuto()], intoShells to empty and
        dblAugerOverride to 0.  It exists so that computations written against the earlier six-field struct -- such
        as the dated branches of example-Fi.jl -- keep running unchanged.
"""
function ElectronIonizationScheme(electronEnergies::Array{Float64,1}, excitationFromShells::Array{Shell,1},
                                  excitationToShells::Array{Shell,1}, lValues::Array{Int64,1}, NoExcitations::Int64,
                                  electronEnergyShift::Float64)
    ElectronIonizationScheme(electronEnergies, excitationFromShells, excitationToShells, lValues, NoExcitations,
                             electronEnergyShift, Basics.AbstractProcess[Basics.ImpactExcAuto()], Shell[], 0. )
end


# `Base.string(scheme::ElectronIonizationScheme)`  ... provides a String notation for the variable scheme::ElectronIonizationScheme.
function Base.string(scheme::ElectronIonizationScheme)
    sa = "Electron-impact ionization (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ElectronIonizationScheme)`  ... prepares a proper printout of the scheme::ElectronIonizationScheme.
function Base.show(io::IO, scheme::ElectronIonizationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
    println(io, "excitationFromShells:        $(scheme.excitationFromShells)  ")
    println(io, "excitationToShells:          $(scheme.excitationToShells)  ")
    println(io, "lValues:                     $(scheme.lValues)  ")
    println(io, "NoExcitations:               $(scheme.NoExcitations)  ")
    println(io, "electronEnergyShift:         $(scheme.electronEnergyShift)  ")
    println(io, "processes:                   $(scheme.processes)  ")
    println(io, "intoShells:                  $(scheme.intoShells)  ")
    println(io, "dblAugerOverride:            $(scheme.dblAugerOverride)  ")
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
    + NoExcitations         ::Int64                 
        ... (Maximum) Number of electron replacements in the excited configuration with regard to the initial configurations/multiplets.
    + excitationFromShells  ::Array{Shell,1}    
        ... List of shells from which excitations are to be considered.
    + excitationToShells    ::Array{Shell,1}    
        ... List of shells to which excitations are to be considered.
    + printTransitions      ::Bool      
        ... Print the bound-bound line list that the opacity is built from -- wavelength, oscillator strength in
            both gauges and the lower-level energy -- so that one can see which lines dominate a given bin.

        NB an energy shift is deliberately NOT a field of this scheme.  It adjusts how the computed lines are
        INTERPRETED, not which lines are computed, and therefore belongs to the simulation side, where
        Cascade.ExpansionOpacities.transitionEnergyShift provides it.  The former meanEnergyShift here was read
        nowhere and merely duplicated that field, inviting the two to disagree.
"""
struct   ExpansionOpacityScheme  <:  Cascade.AbstractCascadeScheme
    multipoles              ::Array{EmMultipole}  
    minPhotonEnergy         ::Float64 
    maxPhotonEnergy         ::Float64   
    NoExcitations           ::Int64
    excitationFromShells    ::Array{Shell,1}
    excitationToShells      ::Array{Shell,1}
    printTransitions        ::Bool
end


"""
`Cascade.ExpansionOpacityScheme()`  ... constructor for an 'default' instance of a Cascade.ExpansionOpacityScheme.
"""
function ExpansionOpacityScheme()
    ExpansionOpacityScheme([E1], 0., 1.0, 1, Shell[], Shell[], false)
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

        This scheme always computes the DIRECT electron-impact excitation ImpactExc() and carries no list of
        processes: there is nothing to choose from.  The resonant channel -- capture into a doubly-excited state
        with subsequent re-autoionization into an excited level of the same ion, which yields the same final
        state -- belongs to Cascade.ElectronExcitationScheme.

    + fromShells             ::Array{Shell,1}    
        ... List of shells from which impact-excitations are to be considered.
    + toShells               ::Array{Shell,1}    
        ... List of shells into which impact-excitations are to be considered, including possibly already occupied shells.
    + electronEnergies       ::Array{Float64,1}                
        ... List of electron energies for which this electron-impact excitation scheme is to be calculated.  The
            collision strengths are computed AT these energies; a Cascade.EieRateCoefficients simulation later
            interpolates over them, so three or more well-spread values are needed for any rate coefficient.
    + lValues                ::Array{Int64,1}
        ... Orbital angular momentum values of the free electron, i.e. the partial waves that are summed over.
            maxKappa is taken as maximum(lValues)+1.  NOTE that a truncated sum does not merely lower the collision
            strength, it can invert its energy dependence, and it does so silently; read the `convergence` column.
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
    ImpactExcitationScheme(Shell[], Shell[], Float64[], Int64[], 0, 0., 0. )
end


# `Base.string(scheme::ImpactExcitationScheme)`  ... provides a String notation for the variable scheme::ImpactExcitationScheme.
function Base.string(scheme::ImpactExcitationScheme)
    sa = "Direct electron-impact excitation (scheme):"
    return( sa )
end


# `Base.show(io::IO, scheme::ImpactExcitationScheme)`  ... prepares a proper printout of the scheme::ImpactExcitationScheme.
function Base.show(io::IO, scheme::ImpactExcitationScheme)
    sa = Base.string(scheme);                 print(io, sa, "\n")
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
        ... RESERVED FOR THE DIRECT electron-impact ionization channel; NOT IMPLEMENTED, and there is no
            Cascade.perform for it.  Recorded here so that the gap is visible rather than surprising.

            The obstacle is not effort but shape.  JAC's ImpactIonization module is semi-empirical, not
            amplitude-based: its Settings are an AbstractEmpiricalSettings, it is driven from an
            Empirical.Computation rather than an Atomic.Computation, it takes a Basis rather than two multiplets,
            and it has no Line type at all.  Its elementary datum is an ImpactIonization.CrossSection per
            (subshell, impact energy), obtained from a BEB, BED or relativistic model.  A cascade built on it
            would therefore store cross sections rather than lines, and the cascade approaches (AverageSCA, SCA)
            would influence only the SCF basis behind the binding energies -- not what they mean elsewhere.

            That is a workable design, but it is a different one, and it should be settled deliberately rather
            than by analogy.  The INDIRECT channel, excitation-autoionization, is implemented and lives in
            Cascade.ElectronIonizationScheme.

    + processes             ::Array{Basics.AbstractProcess,1} 
        ... List of the atomic processes that are supported and should be included into the cascade.
    + electronEnergies      ::Array{Float64,1}                
        ... List of electron energies for which this scheme is to be calculated.
"""
struct   ImpactIonizationScheme  <:  Cascade.AbstractCascadeScheme
    processes               ::Array{Basics.AbstractProcess,1}
    electronEnergies        ::Array{Float64,1}
end


"""
`Cascade.ImpactIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.ImpactIonizationScheme.
"""
function ImpactIonizationScheme()
    ## ImpactExc is what this constructor has always named; only the missing () is repaired here, which is what
    ## made it raise.  Whether an IONIZATION scheme should default to the impact-EXCITATION process is a
    ## separate question: Basics has no impact-ionization tag at all (the nearest is Coulion, i.e. Coulomb
    ## ionization by ion impact), so choosing one is a decision for the maintainer and not a repair.
    ImpactIonizationScheme([ImpactExc()], Float64[] )
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
    println(io, "fromShells:                  $(scheme.fromShells)  ")
    println(io, "toShells:                    $(scheme.toShells)  ")
    println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
    println(io, "lValues:                     $(scheme.lValues)  ")
    println(io, "NoFreeElectronEnergies:      $(scheme.NoFreeElectronEnergies)  ")
    println(io, "maxFreeElectronEnergy:       $(scheme.maxFreeElectronEnergy)  ")
    println(io, "electronEnergyShift:         $(scheme.electronEnergyShift)  ")
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
    ## These two tests referred to bare `photonEnergies` / `electronEnergies` instead of the fields of the
    ## given scheme, so that ANY printout of a PhotoIonizationScheme -- and hence of any Cascade.Computation
    ## carrying one -- raised an UndefVarError.  Every example does println(wa), so the scheme could not be
    ## used in the documented way at all.
    if  length(scheme.photonEnergies) > 0  &&  length(scheme.electronEnergies) > 0
        error("Only photon or electron energies can be specified.")
    end
end




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


# Cascade.AbstractCascadeApproach, Cascade.AverageSCA, Cascade.SCA, Cascade.UserMCA moved to
# Basics.AbstractCascadeApproach / Basics.AverageSCA / Basics.SCA / Basics.UserMCA (module-Basics-inc-abstract.jl),
# to break a circular dependency with DecayYield.Settings; re-exported here via `using ..Basics` (already imported
# above) so both bare (AverageSCA()) and self-qualified (Cascade.AverageSCA()) usage throughout this module's
# included files keep working unchanged.





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
    ## ImpactExcitation.Settings was missing from this Union, so a Cascade.Step could not even be constructed
    ## for an electron-impact excitation cascade -- the scheme was unreachable however its steps were built.
    settings           ::Union{PhotoEmission.Settings, AutoIonization.Settings, PhotoIonization.Settings, PhotoExcitation.Settings,
                                PhotoRecombination.Settings, ImpactExcitation.Settings}
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
    
    if  !isnothing(initialConfigs)   &&   !isnothing(initialMultiplets)   
        error("Only initialConfigs=..  'xor'  initialMultiplets=..  can be specified for a Cascade.Computation().")
    end

    if  isnothing(name)                   namex                 = comp.name                    else  namex = name                              end 
    if  isnothing(nuclearModel)           nuclearModelx         = comp.nuclearModel            else  nuclearModelx = nuclearModel              end 
    if  isnothing(grid)                   gridx                 = comp.grid                    else  gridx = grid                              end 
    if  isnothing(asfSettings)            asfSettingsx          = comp.asfSettings             else  asfSettingsx = asfSettings                end 
    if  isnothing(scheme)                 schemex               = comp.scheme                  else  schemex = scheme                          end 
    if  isnothing(approach)               approachx             = comp.approach                else  approachx = approach                      end 
    if  isnothing(initialConfigs)         initialConfigsx       = comp.initialConfigs          else  initialConfigsx = initialConfigs          end 
    if  isnothing(initialMultiplets)      initialMultipletsx    = comp.initialMultiplets       else  initialMultipletsx = initialMultiplets    end
    
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
## STORED CASCADE DATA (.jld):  a format stamp, and a reader that fails with an explanation.
##
## WHAT THESE FILES ARE, stated plainly because it is easy to assume otherwise.  A cascade `.jld` is written by
## JLD2 from the `results` dictionary, which holds the WHOLE OBJECT GRAPH -- Multiplet -> Level -> Basis -> CsfR
## -> Orbital -> Radial.Grid, plus the process Lines.  JAC itself never reads one back: every "JLD2.load(...)"
## in this module is help text printed for the user, who loads the file in their own session.
##
## THE CONSEQUENCE, measured on 14-Aug-2026 rather than assumed: a file written on 5-Aug-2026 CANNOT be loaded
## by today's JAC.  Not because of anything in the cascade, but because Radial.Grid lost a field and
## Radial.MeshNone ceased to exist in the Radial rework (cfc1848, d2722fe, 968d269).  JLD2 reconstructs what it
## can and then dies with a MethodError about a ReconstructedMutable, six frames deep, saying nothing about the
## cause.  ANY struct change anywhere in that graph does this, and there have already been several.
##
## So a stored cascade is tied to the JAC that wrote it.  The stamp below does not change that -- serialising an
## object graph and then changing the objects cannot be made safe by a version number.  What it does is make the
## failure legible: Cascade.checkDataFile says which JAC wrote the file and what is missing, instead of leaving
## the user with a JLD2 internal error.  The durable fix would be to store DATA (indices, energies, rates)
## rather than objects; that is a design decision about what these files are for, recorded and not taken here.
#######################################################################################################################################


"""
`Cascade.GBL_CASCADE_DATA_FORMAT`
    ... version of the layout of the `results` dictionary that Cascade writes to a .jld file. Increase it
        whenever a KEY of that dictionary is added, removed or renamed. It says nothing about the structs
        stored underneath, which are the part that actually breaks; see the note above.
"""
const GBL_CASCADE_DATA_FORMAT = 1


"""
`Cascade.dataFormatStamp()`
    ... returns a Dict{String,Any} identifying the writer of a cascade .jld file; it is merged into `results`
        by Cascade.writeDataFile below, under the key "data format:".
"""
function dataFormatStamp()
    return( Dict{String,Any}("format version" => Cascade.GBL_CASCADE_DATA_FORMAT,
                             "JAC version"    => Cascade.jacVersion(),
                             "written"        => string(Dates.now())) )
end


"""
`Cascade.jacVersion()`
    ... returns the version of the JenaAtomicCalculator package as a String, or "unknown" if it cannot be
        determined; used only to stamp stored cascade data.
"""
function jacVersion()
    try
        wa = pkgversion(@__MODULE__)
        return( isnothing(wa) ? "unknown" : string(wa) )
    catch
        return( "unknown" )
    end
end


"""
`Cascade.writeDataFile(filename::String, results::Dict{String,Any})`
    ... stamps `results` with Cascade.dataFormatStamp() and writes it to `filename` with JLD2; nothing is
        returned. Every cascade that stores its results goes through here, so that the stamp cannot be
        forgotten at one of the fourteen call sites.
"""
function writeDataFile(filename::String, results::Dict{String,Any})
    results = Base.merge( results, Dict("data format:" => Cascade.dataFormatStamp()) )
    ## jldopen rather than JLD2.save, which routes through FileIO and rejects the .jld extension; this is what
    ## the JLD2.@save macro these call sites used before expands to, so the file layout is unchanged -- one
    ## entry named "results".
    JLD2.jldopen(filename, "w") do f
        f["results"] = results
    end
    return( nothing )
end


"""
`Cascade.checkDataFile(filename::String; stream::IO=stdout)`
    ... reports what a stored cascade .jld file contains and whether THIS JAC can read it; true is returned if
        the file loaded, false otherwise. Nothing is thrown: a file that cannot be read is the normal case this
        function exists for, and it prints which JAC wrote it and which type failed, instead of leaving a JLD2
        internal error.
"""
function checkDataFile(filename::String; stream::IO=stdout)
    println(stream, "\n  Stored cascade data:  $filename")
    if  !isfile(filename)   println(stream, "  >>> no such file.");    return( false )    end
    println(stream, "  size: " * string(round(filesize(filename)/1024/1024, digits=2)) * " MB")
    local results
    try
        ## jldopen rather than JLD2.load: the latter goes through FileIO, which prints its own "Fatal error"
        ## banner to stdout BEFORE throwing -- which is exactly the noise this function exists to replace.
        JLD2.jldopen(filename, "r") do f
            results = haskey(f, "results") ? f["results"] : Dict{String,Any}(k => f[k] for k in keys(f))
        end
    catch  e
        println(stream, "  >>> THIS FILE CANNOT BE READ BY THIS VERSION OF JAC.")
        println(stream, "  >>> " * first(sprint(showerror, e), 400))
        println(stream, "  >>> A cascade .jld stores whole objects (Multiplet, Level, Basis, Orbital, " *
                        "Radial.Grid, ...), so it can only be\n" *
                        "  >>> opened by a JAC whose struct definitions still match the ones that wrote it. " *
                        "Re-run the cascade,\n  >>> or check out the JAC version named in the file's stamp " *
                        "(if it has one; files written before 14-Aug-2026 do not).")
        return( false )
    end
    wa = results
    if  wa isa Dict
        if  haskey(wa, "data format:")
            st = wa["data format:"]
            println(stream, "  written by JAC " * string(get(st, "JAC version", "?")) *
                            " on " * string(get(st, "written", "?")) *
                            ",  data format " * string(get(st, "format version", "?")))
            if  get(st, "format version", 0) != Cascade.GBL_CASCADE_DATA_FORMAT
                println(stream, "  >>> this JAC writes data format $(Cascade.GBL_CASCADE_DATA_FORMAT); the keys may differ.")
            end
        else
            println(stream, "  >>> no format stamp: written before 14-Aug-2026.")
        end
        println(stream, "  keys:")
        for  k in sort(collect(keys(wa)))
            v = wa[k];   sa = v isa AbstractArray ? "$(typeof(v)), $(length(v)) entries" : string(typeof(v))
            println(stream, "     " * rpad(k, 30) * sa)
        end
    else
        println(stream, "  contents: $(typeof(wa))")
    end
    println(stream, "  >>> loaded successfully.")
    return( true )
end

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

include("module-Cascade-inc-approaches.jl")
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
include("module-Cascade-inc-resonant-ionization.jl")
include("module-Cascade-inc-stepwise-decay.jl")
include("module-Cascade-inc-simulations.jl")


end # module


