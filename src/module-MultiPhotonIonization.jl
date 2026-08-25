
"""
`module  JenaAtomicCalculator.MultiPhotonIonization`  
    ... a submodel of JAC for the ionization of an atom or ion by the absorption of SEVERAL photons in a single
        process, i.e. for those multi-photon transitions whose final state carries a FREE ELECTRON.

        THE BOUNDARY, and it is drawn from the other side as well.  MultiPhotonTransition is exclusively
        BOUND-BOUND and says so in its own docstring: "As soon as the final state carries a free electron the
        process belongs to MultiPhotonIonization."  PhotonScattering is for processes with a photon on BOTH sides
        of the reaction.  Here the photons go in and an electron comes out, which is neither.

        WHAT THE ESCAPING ELECTRON CHANGES, relative to two-photon absorption between bound levels:
        + the final state is a continuum partial wave, so there is a PARTIAL-WAVE SUM to truncate, and the
          different final waves are distinct states and add INCOHERENTLY;
        + the ejected-electron energy is NOT free: at one photon energy, energy conservation fixes it at
          epsilon = 2*omega - I_P.  There is therefore no energy-sharing scan of the kind a two-photon EMISSION
          line needs, and omega is the variable to scan instead.  (Sharing returns for two-electron processes,
          and differently for their sequential and simultaneous routes.)

        UNITS.  A two-photon process needs two photons at once, so its rate goes as the SQUARE of the photon flux,
        W [1/s] = sigma^(2) F^2 with F in photons cm^-2 s^-1.  A two-photon cross section is therefore NOT an
        area: it carries cm^4 s, and is quoted in GM (1 GM = 1e-50 cm^4 s), after Maria Goeppert-Mayer, who
        predicted two-photon absorption in her 1931 thesis.  One atomic unit, a0^4 t_au, is 1.8966 GM.  The
        single-beam factor is a convention and is pinned as in MultiPhotonTransition, by requiring the
        monochromatic result to agree with the bichromatic one as omega_1 -> omega_2; every table says so.

        STRUCTURE.  Every process is selected by a scheme <: AbstractMultiPhotonIonizationScheme, following
        MultiPhotonTransition and Cascade: the scheme types live in this file and the algorithms one-per-file in
        the `-inc-` includes.  Two-photon one-electron ionization is implemented; three-photon and two-electron
        processes have a place in the hierarchy and no code.
"""
module MultiPhotonIonization

using  Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron,
       ..Nuclear, ..Radial, ..TableStrings


"""
`abstract type  MultiPhotonIonization.AbstractMultiPhotonIonizationScheme`  
    ... defines an abstract type to distinguish the multi-photon ionization processes.

    + struct TwoPhotonOneElectronScheme    ... two photons absorbed, ONE electron ejected.
"""
abstract type  AbstractMultiPhotonIonizationScheme       end


"""
`struct  MultiPhotonIonization.TwoPhotonOneElectronScheme  <:  AbstractMultiPhotonIonizationScheme`  
    ... two photons of the same energy are absorbed and one electron is ejected; the electron energy follows from
        energy conservation, epsilon = 2*omega - I_P, so omega is the variable to scan.

    + omegas             ::Array{Float64,1}   ... photon energies to be computed [in the selected energy unit].
    + multipoles         ::Array{EmMultipole,1}  ... multipoles of the radiation field to be included.
    + resonanceTolerance ::Float64            ... a computation is REFUSED if an intermediate state comes this
                                                   close to the energy shell; see below.
"""
struct   TwoPhotonOneElectronScheme  <:  MultiPhotonIonization.AbstractMultiPhotonIonizationScheme
    omegas               ::Array{Float64,1}
    multipoles           ::Array{EmMultipole,1}
    resonanceTolerance   ::Float64
end


"""
`MultiPhotonIonization.TwoPhotonOneElectronScheme()`  ... constructor for a default TwoPhotonOneElectronScheme.
"""
function TwoPhotonOneElectronScheme()
    TwoPhotonOneElectronScheme(Float64[], [E1], 1.0e-3)
end


# `Base.show(io::IO, scheme::MultiPhotonIonization.TwoPhotonOneElectronScheme)`  ... prepares a proper printout.
function Base.show(io::IO, scheme::MultiPhotonIonization.TwoPhotonOneElectronScheme)
    println(io, "omegas:                   $(scheme.omegas)  ")
    println(io, "multipoles:               $(scheme.multipoles)  ")
    println(io, "resonanceTolerance:       $(scheme.resonanceTolerance)  ")
end


"""
`struct  MultiPhotonIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of a multi-photon ionization computation.  The PROCESS itself is
        chosen by the scheme, so that a further process is added by adding a scheme rather than by adding fields here.

    + scheme          ::MultiPhotonIonization.AbstractMultiPhotonIonizationScheme   ... which process is computed.
    + gauges          ::Array{UseGauge,1}   ... gauges to be included; the agreement between them is the only internal
                          check a second-order amplitude has, so both are computed by default.
    + printBefore     ::Bool                ... print all selected lines before the computation starts.
    + lineSelection   ::LineSelection       ... specifies the selected levels, if any.
"""
struct Settings  <:  AbstractProcessSettings
    scheme            ::MultiPhotonIonization.AbstractMultiPhotonIonizationScheme
    gauges            ::Array{UseGauge,1}
    printBefore       ::Bool
    lineSelection     ::LineSelection
end


"""
`MultiPhotonIonization.Settings()`  
    ... constructor for a default set of multi-photon ionization settings; a `settings::MultiPhotonIonization.Settings`
        is returned, which selects two-photon one-electron ionization in both gauges.
"""
function Settings()
    Settings( TwoPhotonOneElectronScheme(), [UseCoulomb, UseBabushkin], false, LineSelection() )
end


"""
`MultiPhotonIonization.Settings(set::MultiPhotonIonization.Settings;`
    
        scheme=..,              gauges=..,              printBefore=..,         lineSelection=..)
                    
    ... constructor for modifying the given MultiPhotonIonization.Settings by 'overwriting' the previously selected
        parameters; a `settings::MultiPhotonIonization.Settings` is returned.
"""
function Settings(set::MultiPhotonIonization.Settings;
    scheme::Union{Nothing,MultiPhotonIonization.AbstractMultiPhotonIonizationScheme}=nothing,
    gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing)

    if  isnothing(scheme)           schemex         = set.scheme         else   schemex         = scheme         end
    if  isnothing(gauges)           gaugesx         = set.gauges         else   gaugesx         = gauges         end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection  else   lineSelectionx  = lineSelection  end

    Settings( schemex, gaugesx, printBeforex, lineSelectionx)
end


# `Base.show(io::IO, settings::MultiPhotonIonization.Settings)`  ... prepares a proper printout of the settings.
function Base.show(io::IO, settings::MultiPhotonIonization.Settings)
    println(io, "scheme:                   $(settings.scheme)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
end


include("module-MultiPhotonIonization-inc-2p-one-electron-hydrogenic.jl")

end # module
