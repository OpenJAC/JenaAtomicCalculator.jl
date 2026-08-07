
"""
`module  JAC.MultiPhotonTransition`
... a submodel of JAC that contains all methods for computing multi-photon EXCITATION and DECAY between BOUND
    atomic levels, i.e. the absorption or emission of two (and, later, three) photons in a single transition.

    RENAMED from `MultiPhotonDeExcitation` on 06-Aug-2026. The old name described only half of what the module
    does -- it has always carried absorption processes as well -- and a name that is wrong for half of its
    content is exactly what misleads a later reader.

    THE BOUNDARY, which must not be blurred. This module is exclusively BOUND-BOUND: initial and final states
    are levels of a `Multiplet`, the total transition energy is fixed by `E_i - E_f`, and an emission line is
    resolved by energy SHARINGS between the photons rather than by scanning a photon energy. As soon as the
    final state carries a free electron the process belongs to `MultiPhotonIonization`, `MultiPhotonDoubleIon`
    or the two-colour modules; and one photon with TWO electrons jumping is `TwoElectronOnePhoton`.

    STRUCTURE. Every process is selected by a scheme <: AbstractMultiPhotonScheme, following the same pattern as
    `Cascade` and `Plasma`: all scheme types and their constructors live in this file, while the algorithms live
    one-scheme-per-file in the `-inc-` includes below, and `computeLines` dispatches on the scheme.
"""
module MultiPhotonTransition

using Printf, QuadGK, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear,
                        ..PhotoEmission, ..Radial, ..TableStrings


"""
`abstract type MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... defines an abstract type to distinguish the observables that may be requested for a multi-photon process;
        by far not all of them are meaningful for every scheme. See also the JAC User Guide:

    + struct EnergyDiffCs
        ... energy-differential cross section / rate for initially unpolarized atoms (emission schemes).
    + struct TotalAlpha0
        ... the scalar (K = 0) part of the two-photon transition strength (absorption schemes).
    + struct TotalCsLinear
        ... total cross section for linearly-polarized incident radiation.
    + struct TotalCsRightCircular
        ... total cross section for right-circularly polarized incident radiation.
    + struct TotalCsUnpolarized
        ... total cross section for unpolarized incident radiation.
    + struct TotalCsDensityMatrix
        ... total cross section for incident radiation of arbitrary polarization, given by the 2x2 photon
            density matrix, i.e. by the Stokes parameters in `Settings.stokes`.

    ONE CALCULATION, EVERY OBSERVABLE. These are not independent computations. The two-photon amplitude is
    decomposed into its irreducible parts of rank K, where K couples the two photon multipoles
    (K in oplus(J_f, J_i)); every polarization observable is then a fixed algebraic combination of the SAME
    |M_K|^2. The amplitudes are therefore evaluated once, and the polarization cases cost nothing extra:

        TotalAlpha0            the K = 0 term alone -- the only one surviving for J_i = J_f = 0, i.e. exactly
                               the classic S -> S two-photon transition
        TotalCsLinear          Stokes (P1,P2,P3) = (1,0,0)
        TotalCsRightCircular   Stokes (0,0,1)
        TotalCsUnpolarized     Stokes (0,0,0)
        TotalCsDensityMatrix   the general case, of which the three above are special values

    Their RATIOS are fixed by angular algebra alone and are therefore a check on the implementation that is
    independent of any overall normalization.
"""
abstract type  AbstractMultiPhotonProperty       end


"""
`struct  MultiPhotonTransition.EnergyDiffCs  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the energy-differential cross sections/rates of a multi-photon emission process.
"""
struct   EnergyDiffCs  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

## Note the `io` argument: until 06-Aug-2026 all of these were declared `Base.show(property::X)` -- without io --
## while their bodies called `print(io, ...)`, so every one of them raised an UndefVarError if ever shown.
function Base.show(io::IO, property::EnergyDiffCs)
    print(io, "energy-differential cross sections for multi-photon emission", "\n")
end


"""
`struct  MultiPhotonTransition.TotalAlpha0  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the scalar (K = 0) part of the two-photon transition strength; see the units note at
        MultiPhotonTransition.Settings.
"""
struct   TotalAlpha0  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

function Base.show(io::IO, property::TotalAlpha0)
    print(io, "scalar (K = 0) two-photon transition strength alpha_0", "\n")
end


"""
`struct  MultiPhotonTransition.TotalCsLinear  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the total (absorption) cross section for linearly-polarized incident radiation,
        i.e. Stokes (1,0,0).
"""
struct   TotalCsLinear  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

function Base.show(io::IO, property::TotalCsLinear)
    print(io, "total cross section for linearly-polarized radiation", "\n")
end


"""
`struct  MultiPhotonTransition.TotalCsRightCircular  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the total (absorption) cross section for right-circularly polarized incident radiation,
        i.e. Stokes (0,0,1).
"""
struct   TotalCsRightCircular  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

function Base.show(io::IO, property::TotalCsRightCircular)
    print(io, "total cross section for right-circularly polarized radiation", "\n")
end


"""
`struct  MultiPhotonTransition.TotalCsUnpolarized  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the total (absorption) cross section for unpolarized incident radiation, i.e. Stokes (0,0,0).
"""
struct   TotalCsUnpolarized  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

function Base.show(io::IO, property::TotalCsUnpolarized)
    print(io, "total cross section for unpolarized radiation", "\n")
end


"""
`struct  MultiPhotonTransition.TotalCsDensityMatrix  <:  MultiPhotonTransition.AbstractMultiPhotonProperty`
    ... to request the total (absorption) cross section for incident radiation of ARBITRARY polarization, as
        described by the 2x2 photon density matrix.

        DEFINED 06-Aug-2026; it had been listed in the docstring of AbstractMultiPhotonProperty since the module
        was written but never declared, so requesting it was an UndefVarError. The density matrix of a photon is
        fixed by its three Stokes parameters, which are taken from `Settings.stokes` and carried by the existing
        `Basics.ExpStokes` type -- the same type `PhotoExcitation.Settings` already uses. Linear, right-circular
        and unpolarized light are the special values (1,0,0), (0,0,1) and (0,0,0) of that vector, so this
        property is the general case and the other three are conveniences.
"""
struct   TotalCsDensityMatrix  <:  MultiPhotonTransition.AbstractMultiPhotonProperty      end

function Base.show(io::IO, property::TotalCsDensityMatrix)
    print(io, "total cross section for radiation with a user-defined density matrix (Stokes parameters)", "\n")
end


#################################################################################################################
#################################################################################################################


"""
`abstract type MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... defines an abstract type to distinguish the different multi-photon processes; see also:

    + struct TwoPhotonEmissionScheme
        ... two-photon emission (decay) between two bound levels, often for initially-unpolarized atoms.
    + struct TwoPhotonAbsorptionScheme
        ... two-photon absorption (excitation) with monochromatic and equally-polarized photons, i.e. both
            photons from the same beam.
    + struct TwoPhotonAbsorptionBichromaticScheme
        ... two-photon absorption by photons of two different, well-defined frequencies, usually from two
            different (oriented) beams.
    + struct ThreePhotonEmissionScheme
        ... three-photon emission; NOT yet implemented, see module-MultiPhotonTransition-inc-3p.jl.
    + struct ThreePhotonAbsorptionScheme
        ... three-photon absorption; NOT yet implemented, see module-MultiPhotonTransition-inc-3p.jl.

    Each scheme carries ITS OWN parameters. That is deliberate: the previous design kept `NoEnergySharings` in
    Settings while `TwoPhotonEmission` carried its own `noSharings`, so the same quantity had two sources of
    truth and nothing said which one won.
"""
abstract type  AbstractMultiPhotonScheme       end


"""
`struct  MultiPhotonTransition.TwoPhotonEmissionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... a scheme for two-photon emission between two bound levels.

    + properties   ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
        ... observables to be computed for this scheme.
    + noSharings   ::Int64
        ... number of energy sharings between the two emitted photons, placed at the zeros of a Gauss-Legendre
            quadrature so that the differential rates integrate to the total rate.
"""
struct   TwoPhotonEmissionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme
    properties          ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
    noSharings          ::Int64
end


"""
`MultiPhotonTransition.TwoPhotonEmissionScheme()`  ... constructor for the default values.
"""
function TwoPhotonEmissionScheme()
    TwoPhotonEmissionScheme( AbstractMultiPhotonProperty[EnergyDiffCs()], 4 )
end


"""
`MultiPhotonTransition.TwoPhotonEmissionScheme(scheme::MultiPhotonTransition.TwoPhotonEmissionScheme;`

        properties=..,      noSharings=..)

    ... the standard JAC keyword copy-constructor.
"""
function TwoPhotonEmissionScheme(scheme::MultiPhotonTransition.TwoPhotonEmissionScheme;
    properties::Union{Nothing,Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}}=nothing,
    noSharings::Union{Nothing,Int64}=nothing)

    if  isnothing(properties)   propertiesx = scheme.properties   else   propertiesx = properties   end
    if  isnothing(noSharings)   noSharingsx = scheme.noSharings   else   noSharingsx = noSharings   end

    TwoPhotonEmissionScheme(propertiesx, noSharingsx)
end


function Base.string(scheme::MultiPhotonTransition.TwoPhotonEmissionScheme)
    return( "Two-photon emission between two bound levels:" )
end


function Base.show(io::IO, scheme::MultiPhotonTransition.TwoPhotonEmissionScheme)
    println(io, Base.string(scheme))
    println(io, "properties:                 $(scheme.properties)  ")
    println(io, "noSharings:                 $(scheme.noSharings)  ")
end


"""
`struct  MultiPhotonTransition.TwoPhotonAbsorptionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... a scheme for two-photon absorption with monochromatic, equally-polarized photons from the SAME beam.
        Both photons then carry omega = (E_f - E_i)/2.

    + properties   ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
        ... observables to be computed for this scheme.
"""
struct   TwoPhotonAbsorptionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme
    properties          ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
end


"""
`MultiPhotonTransition.TwoPhotonAbsorptionScheme()`  ... constructor for the default values.
"""
function TwoPhotonAbsorptionScheme()
    TwoPhotonAbsorptionScheme( AbstractMultiPhotonProperty[TotalAlpha0(), TotalCsLinear()] )
end


"""
`MultiPhotonTransition.TwoPhotonAbsorptionScheme(scheme::MultiPhotonTransition.TwoPhotonAbsorptionScheme;`

        properties=..)

    ... the standard JAC keyword copy-constructor.
"""
function TwoPhotonAbsorptionScheme(scheme::MultiPhotonTransition.TwoPhotonAbsorptionScheme;
    properties::Union{Nothing,Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}}=nothing)

    if  isnothing(properties)   propertiesx = scheme.properties   else   propertiesx = properties   end

    TwoPhotonAbsorptionScheme(propertiesx)
end


function Base.string(scheme::MultiPhotonTransition.TwoPhotonAbsorptionScheme)
    return( "Two-photon absorption with monochromatic photons (from the same beam):" )
end


function Base.show(io::IO, scheme::MultiPhotonTransition.TwoPhotonAbsorptionScheme)
    println(io, Base.string(scheme))
    println(io, "properties:                 $(scheme.properties)  ")
end


"""
`struct  MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... a scheme for two-photon absorption by photons of two different frequencies, usually from two different
        (oriented) beams.

    + properties   ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
        ... observables to be computed for this scheme.
    + omegaLess    ::Float64
        ... energy of the photon with the SMALLER frequency; the larger one follows from energy conservation,
            omegaMore = (E_f - E_i) - omegaLess.

    THE BICHROMATIC CASE IS THE UNAMBIGUOUS ONE, and that matters beyond its own physics: with two
    distinguishable beams the rate is W = sigma^(2) * F_1 * F_2 with no combinatorial factor to argue about, so
    it is the reference against which the single-beam convention of TwoPhotonAbsorptionScheme is fixed.
"""
struct   TwoPhotonAbsorptionBichromaticScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme
    properties          ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
    omegaLess           ::Float64
end


"""
`MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme()`  ... constructor for the default values.
"""
function TwoPhotonAbsorptionBichromaticScheme()
    TwoPhotonAbsorptionBichromaticScheme( AbstractMultiPhotonProperty[TotalCsUnpolarized()], 0.05 )
end


"""
`MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme(scheme::MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme;`

        properties=..,      omegaLess=..)

    ... the standard JAC keyword copy-constructor.
"""
function TwoPhotonAbsorptionBichromaticScheme(scheme::MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme;
    properties::Union{Nothing,Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}}=nothing,
    omegaLess::Union{Nothing,Float64}=nothing)

    if  isnothing(properties)   propertiesx = scheme.properties   else   propertiesx = properties   end
    if  isnothing(omegaLess)    omegaLessx  = scheme.omegaLess    else   omegaLessx  = omegaLess    end

    TwoPhotonAbsorptionBichromaticScheme(propertiesx, omegaLessx)
end


function Base.string(scheme::MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme)
    return( "Two-photon absorption with bi-chromatic photons (from different beams):" )
end


function Base.show(io::IO, scheme::MultiPhotonTransition.TwoPhotonAbsorptionBichromaticScheme)
    println(io, Base.string(scheme))
    println(io, "properties:                 $(scheme.properties)  ")
    println(io, "omegaLess:                  $(scheme.omegaLess)  ")
end


"""
`struct  MultiPhotonTransition.ThreePhotonEmissionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... a scheme for three-photon emission between two bound levels. NOT YET IMPLEMENTED; the scheme exists so
        that the hierarchy is complete and so that selecting it fails with a clear message rather than silently.

    + properties   ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
    + noSharings   ::Int64   ... sharings of the total energy among THREE photons (a two-dimensional simplex).
"""
struct   ThreePhotonEmissionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme
    properties          ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
    noSharings          ::Int64
end


"""
`MultiPhotonTransition.ThreePhotonEmissionScheme()`  ... constructor for the default values.
"""
function ThreePhotonEmissionScheme()
    ThreePhotonEmissionScheme( AbstractMultiPhotonProperty[EnergyDiffCs()], 4 )
end


function Base.string(scheme::MultiPhotonTransition.ThreePhotonEmissionScheme)
    return( "Three-photon emission between two bound levels (NOT yet implemented):" )
end


function Base.show(io::IO, scheme::MultiPhotonTransition.ThreePhotonEmissionScheme)
    println(io, Base.string(scheme))
    println(io, "properties:                 $(scheme.properties)  ")
    println(io, "noSharings:                 $(scheme.noSharings)  ")
end


"""
`struct  MultiPhotonTransition.ThreePhotonAbsorptionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme`
    ... a scheme for three-photon absorption. NOT YET IMPLEMENTED; see ThreePhotonEmissionScheme.

    + properties   ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
"""
struct   ThreePhotonAbsorptionScheme  <:  MultiPhotonTransition.AbstractMultiPhotonScheme
    properties          ::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1}
end


"""
`MultiPhotonTransition.ThreePhotonAbsorptionScheme()`  ... constructor for the default values.
"""
function ThreePhotonAbsorptionScheme()
    ThreePhotonAbsorptionScheme( AbstractMultiPhotonProperty[TotalCsUnpolarized()] )
end


function Base.string(scheme::MultiPhotonTransition.ThreePhotonAbsorptionScheme)
    return( "Three-photon absorption (NOT yet implemented):" )
end


function Base.show(io::IO, scheme::MultiPhotonTransition.ThreePhotonAbsorptionScheme)
    println(io, Base.string(scheme))
    println(io, "properties:                 $(scheme.properties)  ")
end


#################################################################################################################
#################################################################################################################


"""
`struct  MultiPhotonTransition.Settings  <:  AbstractProcessSettings`
    ... defines the settings for computing multi-photon excitation and decay between bound levels.

    + scheme             ::MultiPhotonTransition.AbstractMultiPhotonScheme
        ... the multi-photon process to be computed; carries its own parameters.
    + multipoles         ::Array{EmMultipole,1}   ... multipoles of the radiation field to be included.
    + gauges             ::Array{UseGauge,1}      ... gauges to be included into the computations.
    + intermediateStates ::Union{Multiplet, Array{AtomicState.GreenChannel,1}}
        ... the intermediate states |nu> over which the second-order sum runs; see the note below.
    + selfTolerance      ::Float64
        ... intermediate levels whose energy denominator |E_i + omega - E_nu| falls below this value are
            excluded; see the note below.
    + stokes             ::ExpStokes
        ... Stokes parameters of the incident radiation, used by TotalCsDensityMatrix. (0,0,0) is unpolarized.
    + printBefore        ::Bool                   ... print all selected lines before their evaluation.
    + calcOverview       ::Bool                   ... run in overview mode; see the note below.
    + lineSelection      ::LineSelection          ... selected lines, if any.
    + photonEnergyShift  ::Float64                ... overall shift of all initial-to-final transition energies.

    THE INTERMEDIATE-STATE SUM accepts BOTH representations JAC offers, dispatched internally:
      * a mean-field `Multiplet` -- cheap, and what `TwoElectronOnePhoton` uses for the same kind of
        second-order sum; the right choice for surveys across many atoms;
      * an `Array{AtomicState.GreenChannel,1}` -- the accurate route, cf. Fritzsche, Approximate Atomic Green
        Functions, Molecules 26 (2021).
    Until 06-Aug-2026 the absorption path used the first while the emission path referred to `settings.green`, a
    field that had been commented out of Settings -- so the emission path could not run at all.

    SELF-TOLERANCE is not cosmetic. The denominator (E_i + omega - E_nu) vanishes when an intermediate level is
    resonant, which is the physical boundary between NON-RESONANT and RESONANT multi-photon absorption. The
    guard removes the exactly-singular term; `calcOverview` shows the merely dangerous ones. The same convention
    and default (1.0e-8) is used by `TwoElectronOnePhoton`.

    CALCOVERVIEW = true is a cheap first pass that stops BEFORE the expensive part: it builds the multiplets,
    prints them with their indices and symmetries, and ranks the intermediate levels by their contribution
        | <f|O(mp2)|nu> <nu|O(mp1)|i> / (E_i + omega - E_nu) |
    together with each denominator -- then returns without forming any amplitude or property. It answers three
    questions that cannot be guessed from a configuration list: which intermediate levels actually matter, which
    denominators are near-resonant (so that perturbation theory is in doubt), and whether the intermediate basis
    is large enough -- if the largest contributions sit on the highest levels included, the sum is truncated.

    UNITS (settled 06-Aug-2026; nothing had been defined or documented before). The generalized two-photon cross
    section is given in cm^4 s, defined by W [1/s] = sigma^(2) * F^2 with F the photon flux density in
    photons cm^-2 s^-1, and is additionally reported in GM (1 GM = 1e-50 cm^4 s), the unit the two-photon
    literature quotes. `alpha_0` is the K = 0 component of the same quantity and carries the same units -- NOT
    the `cm^4/Ws` of the previous docstring, which is dimensionally inconsistent with a rate quadratic in the
    photon flux. For a SINGLE beam the two photons are indistinguishable, and whether the rate carries an extra
    factor 2 relative to the two-beam case is pure convention; it is fixed here by requiring the monochromatic
    result to agree with the bichromatic one, W = sigma^(2) * F_1 * F_2, in the limit omega_1 -> omega_2.
"""
struct Settings  <:  AbstractProcessSettings
    scheme              ::MultiPhotonTransition.AbstractMultiPhotonScheme
    multipoles          ::Array{EmMultipole,1}
    gauges              ::Array{UseGauge,1}
    intermediateStates  ::Union{Multiplet, Array{AtomicState.GreenChannel,1}}
    selfTolerance       ::Float64
    stokes              ::ExpStokes
    printBefore         ::Bool
    calcOverview        ::Bool
    lineSelection       ::LineSelection
    photonEnergyShift   ::Float64
end


"""
`MultiPhotonTransition.Settings()`  ... constructor for the default values.
"""
function Settings()
    Settings( TwoPhotonEmissionScheme(), EmMultipole[E1], UseGauge[Basics.UseCoulomb], Multiplet(), 1.0e-8,
              Basics.ExpStokes(), false, false, LineSelection(), 0. )
end


"""
`MultiPhotonTransition.Settings(set::MultiPhotonTransition.Settings;`

        scheme=..,              multipoles=..,          gauges=..,          intermediateStates=..,
        selfTolerance=..,       stokes=..,              printBefore=..,     calcOverview=..,
        lineSelection=..,       photonEnergyShift=..)

    ... the standard JAC keyword copy-constructor, which this Settings previously lacked.
"""
function Settings(set::MultiPhotonTransition.Settings;
    scheme::Union{Nothing,MultiPhotonTransition.AbstractMultiPhotonScheme}=nothing,
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,     gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    intermediateStates::Union{Nothing,Multiplet,Array{AtomicState.GreenChannel,1}}=nothing,
    selfTolerance::Union{Nothing,Float64}=nothing,               stokes::Union{Nothing,ExpStokes}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                    calcOverview::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,         photonEnergyShift::Union{Nothing,Float64}=nothing)

    if  isnothing(scheme)              schemex             = set.scheme              else  schemex             = scheme              end
    if  isnothing(multipoles)          multipolesx         = set.multipoles          else  multipolesx         = multipoles          end
    if  isnothing(gauges)              gaugesx             = set.gauges              else  gaugesx             = gauges             end
    if  isnothing(intermediateStates)  intermediateStatesx = set.intermediateStates  else  intermediateStatesx = intermediateStates  end
    if  isnothing(selfTolerance)       selfTolerancex      = set.selfTolerance       else  selfTolerancex      = selfTolerance       end
    if  isnothing(stokes)              stokesx             = set.stokes              else  stokesx             = stokes              end
    if  isnothing(printBefore)         printBeforex        = set.printBefore         else  printBeforex        = printBefore         end
    if  isnothing(calcOverview)        calcOverviewx       = set.calcOverview        else  calcOverviewx       = calcOverview        end
    if  isnothing(lineSelection)       lineSelectionx      = set.lineSelection       else  lineSelectionx      = lineSelection       end
    if  isnothing(photonEnergyShift)   photonEnergyShiftx  = set.photonEnergyShift   else  photonEnergyShiftx  = photonEnergyShift   end

    Settings( schemex, multipolesx, gaugesx, intermediateStatesx, selfTolerancex, stokesx, printBeforex,
              calcOverviewx, lineSelectionx, photonEnergyShiftx )
end


# `Base.show(io::IO, settings::MultiPhotonTransition.Settings)`
#		... prepares a proper printout of the variable settings::MultiPhotonTransition.Settings.
function Base.show(io::IO, settings::MultiPhotonTransition.Settings)
    println(io, "scheme:                   $(settings.scheme)  ")
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "intermediateStates:       $(MultiPhotonTransition.stringIntermediateStates(settings.intermediateStates))  ")
    println(io, "selfTolerance:            $(settings.selfTolerance)  ")
    println(io, "stokes:                   $(settings.stokes)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "calcOverview:             $(settings.calcOverview)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:        $(settings.photonEnergyShift)  ")
end


"""
`MultiPhotonTransition.stringIntermediateStates(mp::Multiplet)`  or
`MultiPhotonTransition.stringIntermediateStates(gChannels::Array{AtomicState.GreenChannel,1})`
    ... to describe the intermediate-state representation in ONE line, without printing every level of it; a
        String is returned. The previous `Base.show` printed the whole gMultiplet, which for any realistic
        intermediate basis buried everything else in the settings.
"""
function stringIntermediateStates(mp::Multiplet)
    return( "mean-field multiplet with $(length(mp.levels)) level(s)" )
end

function stringIntermediateStates(gChannels::Array{AtomicState.GreenChannel,1})
    nl = 0;   for ch in gChannels    nl = nl + length(ch.gMultiplet.levels)    end
    return( "Green expansion with $(length(gChannels)) channel(s) and $nl level(s) in total" )
end


"""
`MultiPhotonTransition.intermediateLevels(mp::Multiplet, Jsym::LevelSymmetry)`  or
`MultiPhotonTransition.intermediateLevels(gChannels::Array{AtomicState.GreenChannel,1}, Jsym::LevelSymmetry)`
    ... to return those intermediate levels |nu> that carry the total symmetry Jsym and hence contribute to the
        second-order sum for that symmetry; an Array{Level,1} is returned, empty if the basis does not span it.

        THIS IS THE ONE PLACE where the two intermediate-state representations differ, and isolating it here is
        what lets every amplitude routine in this module be written once rather than twice. A mean-field
        Multiplet must be FILTERED by symmetry, whereas a Green expansion is already organised into channels of
        definite symmetry and only needs the matching channel selected.
"""
function intermediateLevels(mp::Multiplet, Jsym::LevelSymmetry)
    return( filter(lev -> LevelSymmetry(lev.J, lev.parity) == Jsym, mp.levels) )
end

function intermediateLevels(gChannels::Array{AtomicState.GreenChannel,1}, Jsym::LevelSymmetry)
    levels = Level[]
    for  ch in gChannels
        if  Jsym == ch.symmetry    append!(levels, ch.gMultiplet.levels)    end
    end
    return( levels )
end


#################################################################################################################
#################################################################################################################


include("module-MultiPhotonTransition-inc-2p-emission.jl")
include("module-MultiPhotonTransition-inc-2p-absorption.jl")
include("module-MultiPhotonTransition-inc-2p-bichromatic.jl")
include("module-MultiPhotonTransition-inc-3p.jl")


end # module
