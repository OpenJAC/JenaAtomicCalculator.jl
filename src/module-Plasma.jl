
"""
`module  JAC.Plasma`  
... a submodel of JAC that contains all methods to set-up and process (simple) plasma computations 
    and simulations.
"""
module Plasma

using  Dates, JLD2, LinearAlgebra, Printf, SpecialFunctions, ..ResultKeys
using  ..AngularMomentum, ..Atomic, ..AtomicState, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..DielectronicRecombination,
       ..Distribution, ..Empirical, ..Hamiltonian, ..ImpactExcitation, ..ImpactIonization,
       ..ManyElectron, ..Nuclear, ..PhotoRecombination, ..Radial, ..RadialIntegrals,
       ..Semiempirical, ..TableStrings, ..FormFactor, ..PhotoEmission, ..PhotoIonization, ..AutoIonization, ..SelfConsistent,
       ..ThreeBodyRecombination


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
    + struct KineticChargeBalanceScheme
        ... the cross-charge-state sibling of CollisionalRadiativeScheme: a kinetic (rate-based, non-LTE)
            population balance across several adjacent charge states, with PhotoEmission + ImpactExcitation
            within each stage and ImpactIonization + PhotoRecombination + ThreeBodyRecombination between
            adjacent stages. Each of the five legs independently selects an Empirical (tabulated/fitted) or
            a structure-based (ab-initio/semi-empirical) treatment -- see Plasma.AbstractRateTreatment.
    + struct SuperConfigurationOpacityScheme
        ... to build one SuperConfigurations.Group per reference configuration, report its Boltzmann
            population at a given plasma temperature, build a SuperConfigurations.TransitionArray for every
            E1-connected pair with its oscillator-strength sum, and -- if opacityOmegas is given -- sum
            every array's Gaussian onto that frequency grid into a photo-absorption cross section
            kappa(omega). All three stages are built; see that struct's own docstring for what even that
            still is not (a macroscopic opacity in cm^2/g).
    + struct SuperConfigurationIonizationScheme
        ... to compute the EQUILIBRIUM (LTE) population across several ionization stages of one species,
            each stage described by its own SuperConfigurations.Group's rather than by diagonalized levels
            -- the natural extension of SuperConfigurationOpacityScheme ACROSS charge states instead of
            within one. See that struct's own docstring for the grand-canonical Saha-Boltzmann equations
            and their self-consistent electron-density solve.
    + struct SuperConfigurationMixtureOpacityScheme
        ... to combine several SPECIES (each its own Nuclear.Model, own number fraction, and now possibly
            several of its own ionization STAGES) into one coupled multi-species Saha-Boltzmann balance and
            one PLANCK MEAN opacity at a shared LTE temperature -- the natural extension of
            SuperConfigurationOpacityScheme ACROSS a composition AND of SuperConfigurationIonizationScheme
            ACROSS species, both axes at once. OPTIONALLY (`boundFreeApprox`, 18-Oct-2026) also sums a
            bound-free (photoionization) contribution between every pair of adjacent stages, reusing
            `Empirical.photoionizationCrossSection`, and OPTIONALLY (`freeFreeIncluded`, 19-Oct-2026) also
            sums a free-free (inverse-bremsstrahlung) contribution per stage, via
            `Empirical.freeFreeAbsorptionCoefficient` -- once EITHER is on, the ROSSELAND mean, previously
            withheld as ill-posed, becomes well-defined and is reported too. See that struct's own
            docstring for the coupled balance, the two-stage mixing rule, both continuum additions, and the
            Planck/Rosseland weight functions.
"""
abstract type  AbstractPlasmaScheme       end


"""
`struct  Plasma.AverageAtomScheme  <:  Plasma.AbstractPlasmaScheme`  
    ... a struct to perform an average-atom computation.

    + nMax                  ::Int64                ... maximum principal quantum number for the subshells of the AA model.
    + lMax                  ::Int64                ... maximum orbital quantum number for the subshells of the AA model.
    + scField               ::AbstractScField      ... self-consistent field of the AA model, e.g. Basics.AaDFSField().
    + calcPhotoionizationCs ::Bool                 ... True, if photoionization cross sections are to be calculated.
    + calcFormFactor        ::Bool                 ... True, if the form factor need to be calculated.
    + calcScatteringFactor  ::Bool                 ... True, if scattering factor need to be calculated.
    + piSubshells           ::Array{Subshell,1}    ... Bound subshells to be included into the photoionization cross sections.
    + omegas                ::Array{Float64,1}     ... energies for the photoionization & scattering factors [in a.u.].
    + qValues               ::Array{Float64,1}     ... q-values for calculating the form factor [in a.u.].
    + calcConductivity      ::Bool                 ... True, if the Kubo-Greenwood conductivity and the optical
                                                       constants are to be calculated.
    + condOmegas            ::Array{Float64,1}     ... photon energies for the conductivity [in a.u.];  carry them
                                                       WELL PAST the deepest ionization threshold and log-spaced,
                                                       since the Kramers-Kronig integral runs to infinity.
    + condBroadening        ::Float64              ... width [a.u.] with which the conductivity lines of the box
                                                       are folded into a spectrum;  a model parameter, see
                                                       Plasma.computeConductivitySpectrum.
    + condCollisionRate     ::Float64              ... electron-ion collision rate gamma [a.u.] of the Drude
                                                       free-electron term.  Zero leaves that term out;  a
                                                       POSITIVE value is used as given;  a NEGATIVE value means
                                                       COMPUTE it from the average atom's own scattering phase
                                                       shifts through `Plasma.dcConductivity`, which is the only
                                                       setting that leaves nothing in the optical chain supplied
                                                       from outside (see that function for what it is worth).
    + calcMeanOpacities     ::Bool                 ... True, if the Planck and Rosseland mean opacities are to be
                                                       formed from the conductivity spectrum.  NEEDS
                                                       calcConductivity, and needs condOmegas to span roughly
                                                       0.1 to 30 kT, since that is where the two thermal weights
                                                       live;  both are checked.
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
    calcConductivity        ::Bool
    condOmegas              ::Array{Float64,1}
    condBroadening          ::Float64
    condCollisionRate       ::Float64
    calcMeanOpacities       ::Bool
end  


"""
`Plasma.AverageAtomScheme(nMax::Int64, lMax::Int64, scField::AbstractScField, calcPhotoionizationCs::Bool,
                           calcFormFactor::Bool, calcScatteringFactor::Bool, piSubshells::Array{Subshell,1},
                           omegas::AbstractArray{<:Real,1}, qValues::AbstractArray{<:Real,1})`
    ... the nine-argument constructor, kept so that every call written before the conductivity existed still
        works;  it leaves the conductivity out.  The two meshes are accepted as ANY real arrays and converted
        here, since a mesh written with whole numbers -- `qValues = [q  for q = 1:10]` -- is an `Array{Int64,1}`:
        the struct's own constructor converts such a list silently, but a typed method of the same arity
        shadows it, so a narrower signature here would refuse the very calls this constructor exists to keep
        working.  A `scheme::AverageAtomScheme` is returned.
"""
function AverageAtomScheme(nMax::Int64, lMax::Int64, scField::AbstractScField, calcPhotoionizationCs::Bool,
                           calcFormFactor::Bool, calcScatteringFactor::Bool, piSubshells::Array{Subshell,1},
                           omegas::AbstractArray{<:Real,1}, qValues::AbstractArray{<:Real,1})
    AverageAtomScheme( nMax, lMax, scField, calcPhotoionizationCs, calcFormFactor, calcScatteringFactor,
                       piSubshells, Float64.(omegas), Float64.(qValues), false, Float64[], 0., 0., false )
end


"""
`struct  Plasma.AverageAtomOutcome`
    ... the complete result of one average-atom computation, as named fields rather than as string keys into a
        Dict.

    **THE DIVISION OF THE FIELDS IS THE SCHEME'S OWN DIVISION.**  The cell block is ALWAYS filled, because a
    converged cell always has it;  every optional property is `Union{Nothing,...}` and is `nothing` exactly when
    its switch was not set.  That distinction is the reason for the type:  a string-keyed Dict cannot tell
    "this was not computed" from "this came out zero", and a conductivity of zero and an absent conductivity
    mean entirely different things to whoever reads the number afterwards.

    + temperature          ::Float64   ... the temperature [a.u.] the cell was converged at.
    + massDensity          ::Float64   ... the mass density [g/cm^3].
    + radiusWS             ::Float64   ... the Wigner-Seitz cell radius [a.u.];  NOT the B-spline box, which is
                                           several times larger and must be -- see SelfConsistent's own note.
    + iCell                ::Int64     ... the grid index of the cell boundary, i.e. the last point with
                                           r <= radiusWS.  Every cell integral runs 1:iCell.
    + chemicalMu           ::Float64   ... the chemical potential [a.u.], fixed by neutrality of the cell.
    + meanCharge           ::Float64   ... Z*, the mean ionization, i.e. the free-electron count of the cell.
    + noBoundElectrons     ::Float64   ... the bound-electron count;  noBound + meanCharge = Z exactly.
    + radialDensity        ::Array{Float64,1}        ... the total electron density, carrying Z over the cell.
    + orbitals             ::Dict{Subshell, Orbital} ... the converged bound orbitals.
    + couplingParameter    ::Float64   ... the COUPLING PARAMETER Gamma = Z*^2/(R_WS kT):  the ratio of the
                                           typical Coulomb energy between neighbouring ions to their thermal
                                           energy.  Gamma << 1 is a weakly coupled (ideal) plasma, Gamma >> 1 a
                                           strongly coupled one, Gamma ~ 1 the warm-dense regime this model is
                                           for.
    + degeneracyParameter  ::Float64   ... the DEGENERACY PARAMETER Theta = kT/E_F with E_F the free-electron
                                           Fermi energy:  Theta << 1 is a degenerate electron gas, Theta >> 1 a
                                           classical one.  **Gamma and Theta say WHERE in the (rho,T) plane a
                                           number came from, so they are always reported and never optional** --
                                           a result quoted without them can be misused in a regime the model
                                           does not describe.
    + qValues              ::Union{Nothing, Array{Float64,1}}  ... momentum transfers of the form factor.
    + formFactors          ::Union{Nothing, Array{Float64,1}}  ... F(q), integrated over the cell.
    + piOmegas             ::Union{Nothing, Array{Float64,1}}  ... photon energies of the cross sections.
    + piCrossSections      ::Union{Nothing, Any}               ... bound-free cross sections per subshell.
    + scatteringFactor     ::Union{Nothing, Array{Float64,1}}  ... the imaginary scattering factor f_2.
    + massAttenuation      ::Union{Nothing, Array{Float64,1}}  ... mu/rho [cm^2/g].
    + condOmegas           ::Union{Nothing, Array{Float64,1}}  ... photon energies of the conductivity.
    + conductivity         ::Union{Nothing, Array{Float64,1}}  ... Re sigma(omega) [a.u.].
    + condLines            ::Union{Nothing, Any}               ... the individual box transitions and their
                                                                   strengths, classified by kind.
    + condSumRules         ::Union{Nothing, Any}               ... the per-channel sum rules.
    + optical              ::Union{Nothing, Any}               ... eps(omega), n + ik, the absorption coefficient.
    + dcConductivity       ::Union{Nothing, Any}               ... sigma(0) and the collision rate, from the
                                                                   scattering phase shifts.
    + meanOpacities        ::Union{Nothing, Any}               ... the Planck and Rosseland mean opacities.
"""
struct  AverageAtomOutcome
    temperature             ::Float64
    massDensity             ::Float64
    radiusWS                ::Float64
    iCell                   ::Int64
    chemicalMu              ::Float64
    meanCharge              ::Float64
    noBoundElectrons        ::Float64
    radialDensity           ::Array{Float64,1}
    orbitals                ::Dict{Subshell, Orbital}
    couplingParameter       ::Float64
    degeneracyParameter     ::Float64
    qValues                 ::Union{Nothing, Array{Float64,1}}
    formFactors             ::Union{Nothing, Array{Float64,1}}
    piOmegas                ::Union{Nothing, Array{Float64,1}}
    piCrossSections         ::Union{Nothing, Any}
    scatteringFactor        ::Union{Nothing, Array{Float64,1}}
    massAttenuation         ::Union{Nothing, Array{Float64,1}}
    condOmegas              ::Union{Nothing, Array{Float64,1}}
    conductivity            ::Union{Nothing, Array{Float64,1}}
    condLines               ::Union{Nothing, Any}
    condSumRules            ::Union{Nothing, Any}
    optical                 ::Union{Nothing, Any}
    dcConductivity          ::Union{Nothing, Any}
    meanOpacities           ::Union{Nothing, Any}
end


"""
`Plasma.AverageAtomOutcome(temperature::Float64, massDensity::Float64, radiusWS::Float64, iCell::Int64,
                            chemicalMu::Float64, meanCharge::Float64, noBoundElectrons::Float64,
                            radialDensity::Array{Float64,1}, orbitals::Dict{Subshell, Orbital})`
    ... constructor for the CELL alone, with every optional property left as `nothing`;  the coupling and
        degeneracy parameters are derived here rather than passed, since they follow from the cell:

            Gamma = Z*^2 / (R_WS kT)            ... coupling, the ion-sphere form
            Theta = kT / E_F,   E_F = (3 pi^2 n_e)^(2/3) / 2,   n_e = Z* / (4pi/3 R_WS^3)

        An `outcome::Plasma.AverageAtomOutcome` is returned.
"""
function AverageAtomOutcome(temperature::Float64, massDensity::Float64, radiusWS::Float64, iCell::Int64,
                            chemicalMu::Float64, meanCharge::Float64, noBoundElectrons::Float64,
                            radialDensity::Array{Float64,1}, orbitals::Dict{Subshell, Orbital})
    volume = 4pi/3 * radiusWS^3
    gamma  = (radiusWS > 0. && temperature > 0.) ? meanCharge^2 / (radiusWS * temperature) : 0.
    nE     = volume > 0. ? meanCharge / volume : 0.
    eFermi = nE > 0. ? (3pi^2 * nE)^(2/3) / 2 : 0.
    theta  = eFermi > 0. ? temperature / eFermi : 0.

    AverageAtomOutcome( temperature, massDensity, radiusWS, iCell, chemicalMu, meanCharge, noBoundElectrons,
                        radialDensity, orbitals, gamma, theta,
                        nothing, nothing, nothing, nothing, nothing, nothing,
                        nothing, nothing, nothing, nothing, nothing, nothing, nothing )
end


# `Base.show(io::IO, outcome::Plasma.AverageAtomOutcome)`  ... prepares a proper printout of the outcome.
function Base.show(io::IO, outcome::Plasma.AverageAtomOutcome)
    println(io, "Average-atom outcome:")
    println(io, "  temperature:              $(outcome.temperature) a.u. = " *
                "$(Defaults.convertUnits("energy: from atomic to eV", outcome.temperature)) eV  ")
    println(io, "  mass density:             $(outcome.massDensity) g/cm^3  ")
    println(io, "  Wigner-Seitz radius:      $(outcome.radiusWS) a.u.   (cell boundary at grid index " *
                "iCell = $(outcome.iCell))  ")
    println(io, "  chemical potential:       $(outcome.chemicalMu) a.u.  ")
    println(io, "  mean charge Z*:           $(outcome.meanCharge)  ")
    println(io, "  bound electrons:          $(outcome.noBoundElectrons)  ")
    println(io, "  coupling parameter Gamma    = $(outcome.couplingParameter)  " *
                "($(outcome.couplingParameter < 1. ? "weakly" : "strongly") coupled)  ")
    println(io, "  degeneracy parameter Theta  = $(outcome.degeneracyParameter)  " *
                "($(outcome.degeneracyParameter < 1. ? "degenerate" : "classical") electrons)  ")
    wa = String[]
    isnothing(outcome.formFactors)     || push!(wa, "form factors")
    isnothing(outcome.piCrossSections) || push!(wa, "photoionization cross sections")
    isnothing(outcome.scatteringFactor)|| push!(wa, "scattering factor f_2")
    isnothing(outcome.massAttenuation) || push!(wa, "mass attenuation mu/rho")
    isnothing(outcome.conductivity)    || push!(wa, "conductivity + optical constants")
    isnothing(outcome.dcConductivity)  || push!(wa, "DC conductivity")
    isnothing(outcome.meanOpacities)   || push!(wa, "Planck + Rosseland mean opacities")
    println(io, "  optional properties:      " * (isempty(wa) ? "none (the cell alone)" : join(wa, ", ")) * "  ")
end


"""
`Plasma.AverageAtomScheme(scheme::Plasma.AverageAtomScheme;`

        `nMax=..,                lMax=..,                 scField=..,`
        
        `calcPhotoionizationCs=..,  calcFormFactor=..,    calcScatteringFactor=..,`
        
        `piSubshells=..,         omegas=..,               qValues=..,`
        
        `calcConductivity=..,    condOmegas=..,           condBroadening=..,       condCollisionRate=..)`
                        
    ... constructor for modifying the given Plasma.AverageAtomScheme by 'overwriting' the previously selected
        parameters.  Everything not named is carried over unchanged, so a scheme is built by naming only what
        differs from the default -- which is what makes a thirteen-field scheme readable at the call site:

            scheme = Plasma.AverageAtomScheme(Plasma.AverageAtomScheme();
                         nMax = 8, lMax = 5, scField = Basics.AaDFSField(),
                         calcConductivity = true, condOmegas = [ 10^x  for x = -2.0:0.08:2.6 ],
                         condBroadening = 0.05, condCollisionRate = -1.0 )

        **THE STRUCTURE IS NEVER OPTIONAL AND HAS NO FLAG.**  A converged cell always yields the chemical
        potential, the bound/free split, every subshell's energy and fractional occupation, the mean charge and
        the radial density;  the `calc...` switches add the OPTIONAL properties on top of that, each at its own
        cost.  So a scheme with every switch false is a legitimate and useful computation, not an empty one.

        **EACH SWITCH NEEDS ITS OWN MESH, and that dependency is enforced rather than assumed:** 
        `calcFormFactor` needs `qValues`, `calcPhotoionizationCs` needs `omegas` (and reads `piSubshells`), and
        `calcConductivity` needs `condOmegas`.  A switch set without its mesh raises at construction time, where
        the mistake is, rather than returning an empty array from deep inside `perform`.

        **EVERY MESH AND MODEL PARAMETER TAKES WHOLE NUMBERS AS WELL AS DECIMALS**, so `qValues = [1, 2, 3]`
        and `condCollisionRate = -1` are accepted exactly as `[1., 2., 3.]` and `-1.0` are;  the conversion to
        Float64 happens on the way into the struct and the stored values are identical.
"""
function AverageAtomScheme(scheme::Plasma.AverageAtomScheme;
    nMax::Union{Nothing,Int64}=nothing,                          lMax::Union{Nothing,Int64}=nothing, 
    scField::Union{Nothing,AbstractScField}=nothing,             calcPhotoionizationCs::Union{Nothing,Bool}=nothing, 
    calcFormFactor::Union{Nothing,Bool}=nothing,                 calcScatteringFactor::Union{Nothing,Bool}=nothing, 
    piSubshells::Union{Nothing,Array{Subshell,1}}=nothing,       omegas::Union{Nothing,AbstractArray{<:Real,1}}=nothing, 
    qValues::Union{Nothing,AbstractArray{<:Real,1}}=nothing,     calcConductivity::Union{Nothing,Bool}=nothing, 
    condOmegas::Union{Nothing,AbstractArray{<:Real,1}}=nothing,  condBroadening::Union{Nothing,Real}=nothing, 
    condCollisionRate::Union{Nothing,Real}=nothing,      calcMeanOpacities::Union{Nothing,Bool}=nothing)
    
    if  isnothing(nMax)                   nMaxx                  = scheme.nMax                   else   nMaxx                  = nMax                   end 
    if  isnothing(lMax)                   lMaxx                  = scheme.lMax                   else   lMaxx                  = lMax                   end 
    if  isnothing(scField)                scFieldx               = scheme.scField                else   scFieldx               = scField                end 
    if  isnothing(calcPhotoionizationCs)  calcPhotoionizationCsx = scheme.calcPhotoionizationCs  else   calcPhotoionizationCsx = calcPhotoionizationCs  end 
    if  isnothing(calcFormFactor)         calcFormFactorx        = scheme.calcFormFactor         else   calcFormFactorx        = calcFormFactor         end 
    if  isnothing(calcScatteringFactor)   calcScatteringFactorx  = scheme.calcScatteringFactor   else   calcScatteringFactorx  = calcScatteringFactor   end 
    if  isnothing(piSubshells)            piSubshellsx           = scheme.piSubshells            else   piSubshellsx           = piSubshells            end 
    if  isnothing(omegas)                 omegasx                = scheme.omegas                 else   omegasx                = omegas                 end 
    if  isnothing(qValues)                qValuesx               = scheme.qValues                else   qValuesx               = qValues                end 
    if  isnothing(calcConductivity)       calcConductivityx      = scheme.calcConductivity       else   calcConductivityx      = calcConductivity       end 
    if  isnothing(condOmegas)             condOmegasx            = scheme.condOmegas             else   condOmegasx            = condOmegas             end 
    if  isnothing(condBroadening)         condBroadeningx        = scheme.condBroadening         else   condBroadeningx        = condBroadening         end 
    if  isnothing(condCollisionRate)      condCollisionRatex     = scheme.condCollisionRate      else   condCollisionRatex     = condCollisionRate      end 
    if  isnothing(calcMeanOpacities)      calcMeanOpacitiesx     = scheme.calcMeanOpacities      else   calcMeanOpacitiesx     = calcMeanOpacities      end 

    # THE DEPENDENCIES ARE CHECKED HERE, where the mistake is, and not inside perform() where the symptom is.
    # A switch without its mesh used to return an empty array from deep in the computation, which reads as
    # "this property is zero" rather than as "you forgot the mesh".
    if  calcFormFactorx         &&  isempty(qValuesx)
        error("Plasma.AverageAtomScheme(): calcFormFactor = true needs qValues, the momentum transfers [a.u.] "  *
              "at which the form factor is wanted;  none were given.")                                              end
    if  calcPhotoionizationCsx  &&  isempty(omegasx)
        error("Plasma.AverageAtomScheme(): calcPhotoionizationCs = true needs omegas, the photon energies "       *
              "[a.u.] at which the cross section is wanted;  none were given.")                                     end
    if  calcScatteringFactorx   &&  isempty(omegasx)
        error("Plasma.AverageAtomScheme(): calcScatteringFactor = true needs omegas, since f_2 follows from "     *
              "sigma_bf(omega);  none were given.")                                                                 end
    if  calcConductivityx       &&  isempty(condOmegasx)
        error("Plasma.AverageAtomScheme(): calcConductivity = true needs condOmegas, the photon energies [a.u.] " *
              "of the conductivity spectrum.  Carry them WELL PAST the deepest ionization threshold and "        *
              "log-spaced, since the Kramers-Kronig integral runs to infinity;  none were given.")                  end
    if  calcMeanOpacitiesx      &&  !calcConductivityx
        error("Plasma.AverageAtomScheme(): calcMeanOpacities = true needs calcConductivity = true -- the two "    *
              "thermal means are functionals OF the absorption spectrum, so there is nothing to average "       *
              "without it.")                                                                                       end

    AverageAtomScheme( nMaxx, lMaxx, scFieldx, calcPhotoionizationCsx, calcFormFactorx, calcScatteringFactorx,
                       piSubshellsx, omegasx, qValuesx, calcConductivityx, condOmegasx, condBroadeningx,
                       condCollisionRatex, calcMeanOpacitiesx )
end


"""
`Plasma.AverageAtomScheme()`  ... constructor for an 'default' instance of a Plasma.AverageAtomScheme.
"""
function AverageAtomScheme()
    AverageAtomScheme( 1, 0, Basics.AaHSField(), false, false, false, Subshell[], Float64[], Float64[],
                       false, Float64[], 0., 0., false )
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
    println(io, "calcConductivity:           $(scheme.calcConductivity)  ")
    if  scheme.calcConductivity
        println(io, "condOmegas:                 $(length(scheme.condOmegas)) photon energies  ")
        println(io, "condBroadening:             $(scheme.condBroadening)  ")
        println(io, "condCollisionRate:          $(scheme.condCollisionRate)  ")
    end
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
`struct  Plasma.ChargeStage`
    ... one charge state ("ionization stage") entering a Plasma.KineticChargeBalanceScheme: identified by two reference
        configurations of a shared electron count, a ground configuration and one excited configuration, from which a
        small 2-level reference Multiplet is built (mean-field orbitals + CI) whenever the scheme needs real computed
        levels for that stage (cf. Plasma.AbstractRateTreatment). Keeping exactly two configurations per stage is a
        deliberate simplification -- enough to give the within-stage processes (PhotoEmission, ImpactExcitation) a real
        transition to compute, without the larger shell-generation machinery Plasma.CollisionalRadiativeScheme uses for
        a single ion.

    + groundConfig    ::Configuration   ... the stage's ground configuration; also the one used by the cross-stage
                                               processes (ImpactIonization, PhotoRecombination, ThreeBodyRecombination),
                                               which connect ground levels only.
    + excitedConfig   ::Configuration   ... one excited configuration of the same charge state, giving
                                               PhotoEmission/ImpactExcitation a within-stage transition.
"""
struct  ChargeStage
    groundConfig      ::Configuration
    excitedConfig     ::Configuration
end


"""
`abstract type Plasma.AbstractRateTreatment`
    ... distinguishes, for ONE leg (one process) of a Plasma.KineticChargeBalanceScheme, whether that rate is obtained
        from one of Empirical's tabulated/fitted approximations or from JAC's own ab-initio/semi-empirical structure
        computation (SCF + CI, the full process module). Every leg (PhotoEmission, ImpactExcitation, ImpactIonization,
        PhotoRecombination, ThreeBodyRecombination) carries its own AbstractRateTreatment, so a single scheme can mix
        empirical and structure-based legs freely ("mixed use") -- this is the one thing the whole scheme is built to
        make visible. This type is deliberately kept separate from Empirical.AbstractEmpiricalApproximation (it is NOT
        a subtype of it): the two routes are different in kind, not two flavours of the same thing.

    + Plasma.EmpiricalTreatment        ... wraps one Empirical.AbstractEmpiricalApproximation singleton.
    + Plasma.StructureBasedTreatment   ... selects the real computation; carries the one knob that differs between
                                             models, an ImpactIonization.AbstractModel, read only by the
                                             ImpactIonization/ThreeBodyRecombination legs and ignored by the other three.
"""
abstract type  AbstractRateTreatment                                             end


"""
`struct  Plasma.EmpiricalTreatment  <:  Plasma.AbstractRateTreatment`
    ... selects one of Empirical's own tabulated/fitted approximations for one leg of a Plasma.KineticChargeBalanceScheme;
        no SCF or CI is performed for a stage whose every leg uses this treatment. The quantity computed is otherwise
        identical in kind (a plasma rate coefficient, or a spontaneous rate for PhotoEmission) to what
        Plasma.StructureBasedTreatment returns for the same leg -- only the INPUT to that same rate equation differs,
        which is the comparison this scheme exists to make legible.

    + approx   ::Empirical.AbstractEmpiricalApproximation   ... e.g. Empirical.ScaledHydrogenic(), Empirical.VanRegemorter1962(),
                                                                   Empirical.Lotz1967().
"""
struct  EmpiricalTreatment  <:  AbstractRateTreatment
    approx     ::Empirical.AbstractEmpiricalApproximation
end


"""
`struct  Plasma.StructureBasedTreatment  <:  Plasma.AbstractRateTreatment`
    ... selects JAC's own ab-initio/semi-empirical structure computation for one leg of a Plasma.KineticChargeBalanceScheme.
        Despite the name, this is NOT one uniform level of rigor across the five processes: PhotoEmission and
        ImpactExcitation are intrinsically ab-initio (there is no other JAC route for them); ImpactIonization and
        ThreeBodyRecombination are semi-empirical, since eiiModel itself selects a BEB/BED/... cross-section model over
        REAL computed levels; PhotoRecombination is ab-initio, a genuine structure-based cross section. A scheme prints
        which is which (see Plasma.displayKineticChargeBalanceTreatments); this field exists only to carry the one knob
        that the EII/TBR legs still need.

    + eiiModel   ::ImpactIonization.AbstractModel   ... read only by the ImpactIonization and ThreeBodyRecombination legs;
                                                           ignored otherwise.
"""
struct  StructureBasedTreatment  <:  AbstractRateTreatment
    eiiModel     ::ImpactIonization.AbstractModel
end


"""
`Plasma.StructureBasedTreatment()`  ... constructor with the default eiiModel = ImpactIonization.BEBmodel().
"""
function StructureBasedTreatment()
    StructureBasedTreatment( ImpactIonization.BEBmodel() )
end


# `Base.show(io::IO, treatment::EmpiricalTreatment)`  ... provides a String notation for treatment::EmpiricalTreatment.
function Base.show(io::IO, treatment::EmpiricalTreatment)
    print(io, "empirical ($(typeof(treatment.approx)))")
end


# `Base.show(io::IO, treatment::StructureBasedTreatment)`  ... provides a String notation for treatment::StructureBasedTreatment.
function Base.show(io::IO, treatment::StructureBasedTreatment)
    print(io, "structure-based (eiiModel = $(typeof(treatment.eiiModel)))")
end


"""
`struct  Plasma.KineticChargeBalanceScheme  <:  Plasma.AbstractPlasmaScheme`
    ... a kinetic (rate-based, non-LTE) population balance across SEVERAL adjacent charge states of one element -- the
        cross-charge-state sibling to Plasma.CollisionalRadiativeScheme, which is confined by design to one fixed charge
        state. Every chargeStages[k] is connected to itself by PhotoEmission (radiative decay, excited -> ground) and
        ImpactExcitation (collisional excitation/de-excitation, both directions, by detailed balance from one
        computed/estimated rate), and to chargeStages[k+1] by ImpactIonization (up, ground -> ground) and by
        PhotoRecombination + ThreeBodyRecombination (down, ground -> ground, added into the same matrix entry).
        DielectronicRecombination is deliberately NOT included here (an explicit scope decision, not an oversight --
        its resonant capture-then-stabilize mechanism does not reduce to a single rate coefficient the way the other
        four processes do).

        EACH of the five legs carries its OWN Plasma.AbstractRateTreatment (peTreatment, ieTreatment, eiiTreatment,
        rrTreatment, tbrTreatment), so a run can mix a cheap Empirical formula for some processes with JAC's own
        structure computation for others, in any combination -- "mixed use". A run in which every leg is
        Plasma.EmpiricalTreatment never calls a mean-field/CI computation at all and completes in well under a second;
        the moment ANY leg is Plasma.StructureBasedTreatment, a small reference Multiplet (ground + excited
        configuration, mean-field orbitals + frozen-orbital CI) is generated for EVERY stage -- not only the stage(s)
        that leg touches -- since this is cheap here (one mean-field pass per stage) and a per-leg, per-stage decision
        of which levels to generate would add real bookkeeping for little saved time. An Empirical leg still uses its
        own Configuration-level formula regardless, so this costs extra computation, never wrong physics.

    + chargeStages          ::Array{Plasma.ChargeStage,1}   ... the charge states, ordered from the LOWEST to the
                                                                   HIGHEST charge (fewest to most electrons removed);
                                                                   chargeStages[k] and chargeStages[k+1] must differ by
                                                                   exactly one electron in their groundConfig.
    + electronTemperatures ::Array{Float64,1}               ... electron temperatures [K] at which the balance is solved.
    + peTreatment           ::Plasma.AbstractRateTreatment  ... within-stage PhotoEmission.
    + ieTreatment           ::Plasma.AbstractRateTreatment  ... within-stage ImpactExcitation.
    + eiiTreatment          ::Plasma.AbstractRateTreatment  ... cross-stage ImpactIonization (up).
    + rrTreatment           ::Plasma.AbstractRateTreatment  ... cross-stage PhotoRecombination (down).
    + tbrTreatment          ::Plasma.AbstractRateTreatment  ... cross-stage ThreeBodyRecombination (down).
"""
struct  KineticChargeBalanceScheme  <:  Plasma.AbstractPlasmaScheme
    chargeStages           ::Array{Plasma.ChargeStage,1}
    electronTemperatures   ::Array{Float64,1}
    peTreatment            ::Plasma.AbstractRateTreatment
    ieTreatment            ::Plasma.AbstractRateTreatment
    eiiTreatment           ::Plasma.AbstractRateTreatment
    rrTreatment            ::Plasma.AbstractRateTreatment
    tbrTreatment           ::Plasma.AbstractRateTreatment
end


"""
`Plasma.KineticChargeBalanceScheme()`  ... constructor for a default instance of a Plasma.KineticChargeBalanceScheme.
"""
function KineticChargeBalanceScheme()
    KineticChargeBalanceScheme( Plasma.ChargeStage[], [1.0e4],
                                EmpiricalTreatment(Empirical.ScaledHydrogenic()),
                                EmpiricalTreatment(Empirical.VanRegemorter1962()),
                                EmpiricalTreatment(Empirical.Lotz1967()),
                                EmpiricalTreatment(Empirical.ScaledHydrogenic()),
                                EmpiricalTreatment(Empirical.Lotz1967()) )
end


"""
`Plasma.KineticChargeBalanceScheme(scheme::Plasma.KineticChargeBalanceScheme;`

        chargeStages=..,          electronTemperatures=..,
        peTreatment=..,           ieTreatment=..,           eiiTreatment=..,          rrTreatment=..,        tbrTreatment=..)

    ... constructor for modifying the given Plasma.KineticChargeBalanceScheme by 'overwriting' the previously selected
        parameters.
"""
function KineticChargeBalanceScheme(scheme::Plasma.KineticChargeBalanceScheme;
    chargeStages::Union{Nothing,Array{Plasma.ChargeStage,1}}=nothing,
    electronTemperatures::Union{Nothing,Array{Float64,1}}=nothing,
    peTreatment::Union{Nothing,Plasma.AbstractRateTreatment}=nothing,
    ieTreatment::Union{Nothing,Plasma.AbstractRateTreatment}=nothing,
    eiiTreatment::Union{Nothing,Plasma.AbstractRateTreatment}=nothing,
    rrTreatment::Union{Nothing,Plasma.AbstractRateTreatment}=nothing,
    tbrTreatment::Union{Nothing,Plasma.AbstractRateTreatment}=nothing)

    if  isnothing(chargeStages)           chargeStagesx          = scheme.chargeStages          else   chargeStagesx          = chargeStages          end
    if  isnothing(electronTemperatures)   electronTemperaturesx  = scheme.electronTemperatures  else   electronTemperaturesx  = electronTemperatures  end
    if  isnothing(peTreatment)            peTreatmentx           = scheme.peTreatment           else   peTreatmentx           = peTreatment            end
    if  isnothing(ieTreatment)            ieTreatmentx           = scheme.ieTreatment           else   ieTreatmentx           = ieTreatment            end
    if  isnothing(eiiTreatment)           eiiTreatmentx          = scheme.eiiTreatment          else   eiiTreatmentx          = eiiTreatment           end
    if  isnothing(rrTreatment)            rrTreatmentx           = scheme.rrTreatment           else   rrTreatmentx           = rrTreatment            end
    if  isnothing(tbrTreatment)           tbrTreatmentx          = scheme.tbrTreatment          else   tbrTreatmentx          = tbrTreatment           end

    KineticChargeBalanceScheme( chargeStagesx, electronTemperaturesx, peTreatmentx, ieTreatmentx, eiiTreatmentx, rrTreatmentx, tbrTreatmentx )
end


# `Base.show(io::IO, scheme::KineticChargeBalanceScheme)`  ... prepares a proper printout of the scheme::KineticChargeBalanceScheme.
function Base.show(io::IO, scheme::KineticChargeBalanceScheme)
    println(io, "chargeStages:          $(length(scheme.chargeStages)) stage(s)  ")
    println(io, "electronTemperatures:  $(scheme.electronTemperatures)  ")
    println(io, "peTreatment:           $(scheme.peTreatment)  ")
    println(io, "ieTreatment:           $(scheme.ieTreatment)  ")
    println(io, "eiiTreatment:          $(scheme.eiiTreatment)  ")
    println(io, "rrTreatment:           $(scheme.rrTreatment)  ")
    println(io, "tbrTreatment:          $(scheme.tbrTreatment)  ")
end


"""
`struct  Plasma.SuperConfigurationOpacityScheme  <:  Plasma.AbstractPlasmaScheme`
    ... a type for a superconfiguration-based plasma computation: ONE Gaussian `SuperConfigurations.Group`
        per reference configuration (a shared orbital set is built once, from ALL of `refConfigs` together),
        the Boltzmann population each group carries at the plasma temperature
        `computation.settings.temperature`, ONE `SuperConfigurations.TransitionArray` for every PAIR of
        `refConfigs` connected by a single-electron E1 transition (opposite parity), and -- if
        `opacityOmegas` is non-empty -- the SUM of every array's own Gaussian line profile on that
        frequency grid: a per-absorber photo-absorption cross section kappa(omega), in atomic units (a0^2).

        STAGE 3 OF 3. Stage 1 was the population weighting alone; Stage 2 added the per-array transition
        moment (mean energy, width, oscillator-strength sum `gf`) that Stage 1 lacked. This stage folds
        every array's own Gaussian -- whose shape was CALIBRATED, not assumed, against a direct
        diagonalization (skewness +0.28, excess kurtosis +1.11 on a real case; see
        `SuperConfigurations.TransitionArray`'s own docstring) -- onto one shared frequency grid, weighted
        by its lower configuration's population, into the single summed spectrum the scheme is named for.
        `opacityOmegas` left empty (the default) skips this stage exactly as before; Stages 1 and 2 do not
        need it. WHAT THIS STILL IS NOT: a macroscopic opacity in cm^2/g -- that needs an ion number density
        and mass this scheme is not given, a multiplication the caller supplies.

    + refConfigs      ::Array{Configuration,1}   ... the configurations whose `SuperConfigurations.Group`
                                                        this scheme builds, one group per configuration.
    + opacityOmegas   ::Array{Float64,1}          ... photon energies [a.u.] at which to sum the Stage-3
                                                        spectrum; empty (the default) skips Stage 3.
"""
struct  SuperConfigurationOpacityScheme  <:  Plasma.AbstractPlasmaScheme
    refConfigs        ::Array{Configuration,1}
    opacityOmegas     ::Array{Float64,1}
end


"""
`Plasma.SuperConfigurationOpacityScheme()`  ... constructor for a standard (empty) instance.
"""
function SuperConfigurationOpacityScheme()
    SuperConfigurationOpacityScheme( Configuration[], Float64[] )
end


"""
`Plasma.SuperConfigurationOpacityScheme(refConfigs::Array{Configuration,1})`
    ... constructor with `opacityOmegas` defaulted to empty, i.e. Stages 1 and 2 only -- every existing
        caller from before Stage 3 was added keeps working unchanged.
"""
function SuperConfigurationOpacityScheme(refConfigs::Array{Configuration,1})
    SuperConfigurationOpacityScheme( refConfigs, Float64[] )
end


# `Base.show(io::IO, scheme::Plasma.SuperConfigurationOpacityScheme)`  ... prepares a proper printout.
function Base.show(io::IO, scheme::Plasma.SuperConfigurationOpacityScheme)
    println(io, "refConfigs:      $(scheme.refConfigs)  ")
    println(io, "opacityOmegas:   $(scheme.opacityOmegas)  ")
end


"""
`struct  Plasma.SuperConfigurationIonizationScheme  <:  Plasma.AbstractPlasmaScheme`
    ... a type for the EQUILIBRIUM (LTE) population across several ionization stages of ONE species, each
        stage described by `SuperConfigurations.Group`'s rather than by diagonalized levels -- the natural
        extension of `SuperConfigurationOpacityScheme` ACROSS charge states instead of within one.

        THE PHYSICS IS THE GRAND-CANONICAL REFORMULATION OF THE SAHA EQUATION, the SAME one
        `Plasma.SahaBoltzmannScheme` already uses and the SAME equations example-Jc.jl independently
        verified by hand: every (stage q, group g) pair carries a Boltzmann weight
        `g.degeneracy * exp(-(mu_e*q + g.meanEnergy - Eref)/kT)` in ONE shared sum over every stage at
        once, where `mu_e` is the electron chemical potential and `q` the stage's charge (nuclear charge
        minus its own electron count). `mu_e` is tied to the electron density by the standard non-degenerate
        ideal-electron-gas relation, and the two are found together by a SELF-CONSISTENT fixed-point
        iteration on charge neutrality -- reimplemented here rather than calling `SahaBoltzmannScheme`'s own
        code, which is NOT in the test suite and carries isotopic-mixture and ionization-potential-
        depression machinery this scheme deliberately does not need.

        REUSES, NOT REIMPLEMENTS, EVERY STAGE'S OWN STRUCTURE: one `SuperConfigurations.Group` per
        configuration, exactly as `SuperConfigurationOpacityScheme`'s Stage 1 already builds them -- no
        diagonalization, no new structure physics.

        SCOPE, DELIBERATELY NARROW, the same choice `SuperConfigurations` itself made for the same reason:
        ONE species (no isotopic mixture), vacuum ionization energies (no ionization-potential-depression
        model), and every stage considered must carry AT LEAST ONE electron -- a fully stripped bare
        nucleus needs no `Group` at all and is not handled here.

    + refConfigsPerStage   ::Array{Array{Configuration,1},1}   ... one list of reference configurations per
                                                                     ionization stage, ORDERED from the
                                                                     LEAST to the MOST ionized stage
                                                                     considered.
"""
struct  SuperConfigurationIonizationScheme  <:  Plasma.AbstractPlasmaScheme
    refConfigsPerStage     ::Array{Array{Configuration,1},1}
end


"""
`Plasma.SuperConfigurationIonizationScheme()`  ... constructor for a standard (empty) instance.
"""
function SuperConfigurationIonizationScheme()
    SuperConfigurationIonizationScheme( Array{Configuration,1}[] )
end


# `Base.show(io::IO, scheme::Plasma.SuperConfigurationIonizationScheme)`  ... prepares a proper printout.
function Base.show(io::IO, scheme::Plasma.SuperConfigurationIonizationScheme)
    println(io, "refConfigsPerStage:   $(scheme.refConfigsPerStage)  ")
end


"""
`struct  Plasma.SuperConfigurationMixtureOpacityScheme  <:  Plasma.AbstractPlasmaScheme`
    ... a type for a MIXTURE of species, each possibly spread over SEVERAL of its own ionization stages,
        combined at a SHARED LTE temperature `computation.settings.temperature` into one PLANCK MEAN
        opacity over `opacityOmegas`. This is C9's own named remainder, item C20 in the competition list:
        "opacity for a realistic mixture, not one species" -- EXTENDED 18-Oct-2026 to also cover the stage
        axis `SuperConfigurationIonizationScheme` covers for one species alone, so a mixture no longer has
        to be given a single hand-picked charge state per species.

        EACH SPECIES CARRIES ITS OWN `Nuclear.Model` AND OWN GRID PER STAGE, built automatically
        (`Basics.recommendedGrid`) rather than shared from `computation` -- a mixture is, by definition, more
        than one element (and, now, each element possibly more than one charge state), so
        `computation.nuclearModel`/`computation.grid` (sized for ONE species at ONE stage) cannot be reused
        without breaking Rule 12's box-matching requirement for whichever (species,stage) is not the one
        `computation` was built for. `computation.nuclearModel`/`computation.grid` are therefore UNUSED by
        this scheme, the same way `computation.refConfigs` is already unused by
        `SuperConfigurationOpacityScheme`'s own `refConfigs` field.

        THE IONIZATION BALANCE IS SOLVED FIRST, COUPLED ACROSS EVERY SPECIES AT ONCE -- the genuinely new
        physics here, not a wiring exercise: `SuperConfigurationIonizationScheme`'s own charge-neutrality
        equation assumes ONE species supplies every electron, which is wrong the moment a species is a
        trace component of the mixture -- its own ionization balance then depends on an electron density
        set mostly by everyone ELSE. The SAME grand-canonical fixed-point iteration is reused, with the
        charge-neutrality sum lifted to range over every (species, stage) pair:

            n_e  =  n_ion,total * SUM_species  numberFractions[s] * SUM_stages_i  charge[s,i] * p[s,i]

        where `p[s,i]` is species `s`'s OWN Boltzmann population of stage `i`, normalized to 1 WITHIN that
        species (`SUM_i p[s,i] = 1`) at the SHARED chemical potential `mu_e` every species sees -- each
        species keeps its own, purely numerical energy reference for the exponentials (Eref_s), which never
        needs to agree ACROSS species since it cancels out of that species' own normalization. A SINGLE-
        STAGE species has `p[s,1] = 1` identically, whatever `mu_e` is -- which is exactly why the
        backward-compatible constructor below reproduces the ORIGINAL (17-Oct-2026) scheme's numbers
        unchanged: the balance above is trivial when there is only one state to belong to.

        THE MIXING RULE builds on that balance in two stages: WITHIN a species, its own stages' Stage-3
        `kappa_{s,i}(omega)` (each built with EXACTLY `SuperConfigurationOpacityScheme`'s own machinery) are
        weighted by that stage's OWN population `p[s,i]`; ACROSS species, by `numberFractions`:

            kappa_s(omega)    =  SUM_i  p[s,i] * kappa_{s,i}(omega)
            kappa_mix(omega)  =  SUM_s  numberFractions[s] * kappa_s(omega)

        THE PLANCK MEAN is the standard dimensionless-weight form, written in terms of `x = omega/kT`:

            kappa_P  =  INTEGRAL  kappa_mix(x) * w_P(x) dx,     w_P(x) = (15/pi^4) * x^3/(e^x-1)

        normalized so `INTEGRAL w_P(x) dx = 1` over `x in (0,infinity)` -- the standard check this scheme's
        own test uses: a FREQUENCY-INDEPENDENT `kappa_mix` must return exactly that same constant as its
        Planck mean, since a weighted average of a constant is that constant. The integral is evaluated by
        the trapezoidal rule over `opacityOmegas` itself, so the grid must be wide and fine enough to
        resolve `w_P`'s own peak near `x ~ 3-4` AND every species' own line widths -- not merely wide
        enough for the lines alone, as a single-species `opacityOmegas` choice would be.

        **THE ROSSELAND MEAN IS WITHHELD UNLESS `boundFreeApprox` IS GIVEN**, and this was a finding, not
        an oversight: a Rosseland mean is a HARMONIC average, dominated by whichever frequencies are LEAST
        opaque, and a `kappa(omega)` built from E1 lines ALONE -- no bound-free or free-free continuum --
        leaves every frequency between lines a true gap with nothing sensible to divide into. Measured, not
        assumed, 17-Oct-2026: sweeping the numerical floor `Plasma.rosselandMeanOpacity` needs just to avoid
        dividing by machine-zero moved the returned value by MORE THAN ELEVEN ORDERS OF MAGNITUDE (0.049
        down to 7.7e-13) on the SAME line-only spectrum -- proof the number had no content, set by an
        implementation detail rather than physics. The PLANCK mean, an ARITHMETIC (kappa-weighted) average,
        has no such gap-sensitivity and remains reported unconditionally.

        **BOUND-FREE, 18-Oct-2026 (`boundFreeApprox`, optional, default `nothing` = off):** when given an
        `Empirical.AbstractEmpiricalApproximation` (e.g. `Empirical.ScaledHydrogenic()`, the default choice
        once turned on), `Plasma.perform` ALSO sums a photoionization contribution between every pair of
        ADJACENT stages of a species, reusing the ALREADY-VALIDATED `Empirical.photoionizationCrossSection`
        -- Kramers' (1923) formula scaled by the ionized shell's own binding energy and a bound-free Gaunt
        factor, the SAME machinery `Cascade.BoundFreeOpacity` already relies on for exactly this purpose
        elsewhere in JAC. No new ab initio machinery is needed: that function already takes a bare
        `(initialConfig, finalConfig)` pair, which is exactly what `refConfigsPerSpeciesStages`' own
        adjacent-stage configurations already are, so this is substantially a WIRING addition on top of
        machinery this project has already audited (see [[project_empirical_module]]-equivalent history),
        NOT a new physical capability built from scratch. For every GROUP `g` in stage `i` and GROUP `h` in
        stage `i+1` whose configurations differ by exactly one electron in exactly one shell (checked via
        `Basics.extractFromConfigurations(Basics.OccupationDifference(), ...)`, the SAME compatibility test
        `Empirical.photoionizationCrossSection` itself applies), the cross section is weighted by `g`'s own
        ABSOLUTE population within the species -- `stagePopulations[s][i] * (g's population WITHIN stage i)`
        -- and summed; a pair whose configurations are not a genuine single-electron ionization channel is
        simply SKIPPED, not an error, exactly as a same-parity pair is skipped (not an error) in the
        bound-bound Stage 2 above. Reported SEPARATELY from the bound-bound spectrum (own ResultKeys), so a
        user always sees bound-bound, bound-free, and their sum as three distinct, clearly labelled things
        -- never silently merged. Once this exists, `Plasma.rosselandMeanOpacity` is called on the SUMMED
        (bound-bound + bound-free) spectrum and reported.

        **FREE-FREE (inverse bremsstrahlung), 19-Oct-2026 (`freeFreeIncluded`, optional `Bool`, default
        `false`):** when `true`, `Plasma.perform` ALSO sums, for every (species,stage), the free-free
        (inverse-bremsstrahlung) absorption coefficient of a free electron scattering off that stage's own
        ion charge, `Empirical.freeFreeAbsorptionCoefficient` -- the classical Kramers cross section times
        the quantum Gaunt-factor correction, read from Sutherland's (1998) thermally-averaged table and
        derived, step for step, against the typeset equations of Karzas & Latter (1961); see that function's
        own docstring for the full derivation and provenance. Unlike bound-free, this needs no PAIR of
        adjacent stages -- a stage's own charge alone sets the Coulomb field a passing electron feels -- so
        it is available even for a single-stage species, PROVIDED a real electron density exists, i.e.
        `useNumberDensity`/`density` must be set exactly as bound-free's multi-stage case already requires;
        `freeFreeIncluded = true` with no density raises, naming what is missing. A neutral stage (charge 0)
        contributes exactly zero, by construction of the formula's own Z^2 factor, not as a special case.
        Reported under its OWN new key, SuperConfigurationMixtureFreeFreeOpacity, never merged into the
        bound-free number; the TOTAL spectrum/Planck/Rosseland keys now mean "bound-bound plus whichever of
        bound-free/free-free are switched on", and are present whenever EITHER is active.

        **RESONANCES (autoionizing lines), 20-Oct-2026 (`resonanceContinuumGrid`, optional
        `Union{Nothing,Radial.Grid}`, default `nothing` = off):** every stage's own excited configurations
        already include some that lie ABOVE the adjacent, more-ionized stage's ground energy --
        `Plasma.classifyResonantGroups` flags these RESONANT -- and are reached from that stage's own ground
        configuration by an ORDINARY bound-bound E1 transition, exactly like any other line this scheme
        already builds. What makes such a line different is its DECAY: besides radiating, it can
        AUTOIONIZE into the adjacent stage's continuum, which gives it a finite lifetime and hence a
        LORENTZIAN (not delta-function) contribution to its width, on top of the GAUSSIAN statistical
        spread (`TransitionArray.variance`) every line already carries. The two combine into a Voigt
        profile, `Plasma.voigtLineProfile` -- the Faddeeva-function closed form, verified (not assumed)
        against a direct numerical convolution of the two component profiles to `1e-12` relative, against
        its own normalization integral to `1e-10`, and against both single-component limits
        (`gammaFWHM -> 0` reduces to `Plasma.gaussianLineProfile`; `sigma -> 0` reduces to a pure
        Lorentzian), before being written into this scheme.

        The width itself is `Plasma.resonanceAutoionizationWidth` -- a REPRESENTATIVE-LEVEL width, reusing
        JAC's own, already-validated `AutoIonization` process UNCHANGED (no new decay-rate formula), applied
        to the SAME ground-to-resonance `TransitionArray` this scheme's own Stage 2 already builds for that
        pair. Each such line is therefore corrected IN PLACE, from the Gaussian its own Stage 2/3 would
        otherwise give it to the Voigt profile its finite lifetime actually implies -- the correction is
        `gfEff * (voigtLineProfile(...) - gaussianLineProfile(...))`, added directly into that stage's own
        `kappa(omega)` before it is weighted into the species and mixture totals, so every other line (every
        line with no open autoionization channel) is completely unaffected and the whole treatment is OFF,
        byte-identically, whenever `resonanceContinuumGrid` is `nothing`. No new `ResultKeys` category is
        needed for the bound-bound spectrum itself -- it simply contains better-shaped lines where they
        apply -- but the resonances actually corrected are reported separately for inspection, under
        `SuperConfigurationMixtureResonanceCorrections`.

    + nuclearModels              ::Array{Nuclear.Model,1}                      ... one nuclear model per species.
    + refConfigsPerSpeciesStages ::Array{Array{Array{Configuration,1},1},1}
        ... one species' own `refConfigsPerStage` -- EXACTLY `SuperConfigurationIonizationScheme`'s own field
            shape, one `Array{Configuration,1}` per ionization stage of that species -- matching
            `nuclearModels` and `numberFractions` by index. A species with only ONE stage (the common case
            before the coupled balance below existed) is simply a length-1 outer array; the backward-
            compatible constructor below builds exactly that from the OLD one-stage-per-species call shape.
    + numberFractions            ::Array{Float64,1}                            ... each species' number
                                                                        fraction of the mixture; expected to
                                                                        sum to 1, checked, not enforced silently.
    + opacityOmegas              ::Array{Float64,1}                           ... the SHARED photon-energy
                                                                        grid [a.u.] every species' Stage-3
                                                                        spectrum and the Planck mean are built on.
    + boundFreeApprox            ::Union{Nothing,Empirical.AbstractEmpiricalApproximation}
        ... `nothing` (default) turns the bound-free contribution OFF, reproducing every pre-18-Oct-2026
            result unchanged; `Empirical.ScaledHydrogenic()` turns it ON. See above for what it adds.
    + freeFreeIncluded           ::Bool
        ... `false` (default) turns the free-free contribution OFF; `true` turns it ON. See above for what
            it adds. Unlike `boundFreeApprox` there is only ONE trusted free-free formula available, so this
            is a plain switch rather than a choice among `Empirical.AbstractEmpiricalApproximation`s.
    + resonanceContinuumGrid    ::Union{Nothing,Radial.Grid}
        ... `nothing` (default) turns the RESONANCE (autoionizing-line) treatment OFF, reproducing every
            pre-20-Oct-2026 result unchanged. A `Radial.Grid` turns it ON and is the grid
            `AutoIonization.computeLines` uses for the ejected electron's continuum orbital -- NOT auto-
            derived (Rule 12: `Basics.recommendedGrid` is sized for BOUND orbitals only), so the CALLER must
            supply one already checked adequate for the ejected electron's kinetic energy in every channel
            the mixture opens. **ONE GRID SERVES EVERY SPECIES, STAGE AND RESONANT GROUP** -- a known
            simplification, stated rather than hidden: a mixture whose channels span a wide range of
            ejected-electron energies may need a wider or finer grid than any single channel alone would,
            and this field offers no per-channel adjustment. See `Plasma.resonanceAutoionizationWidth`'s own
            docstring for the full Rule-12 discussion this field's caller must satisfy.
    + resonanceAugerSettings    ::AutoIonization.Settings
        ... the settings passed to every `AutoIonization` call the resonance treatment makes; unused when
            `resonanceContinuumGrid` is `nothing`. Default `AutoIonization.Settings()`.
"""
struct  SuperConfigurationMixtureOpacityScheme  <:  Plasma.AbstractPlasmaScheme
    nuclearModels              ::Array{Nuclear.Model,1}
    refConfigsPerSpeciesStages ::Array{Array{Array{Configuration,1},1},1}
    numberFractions            ::Array{Float64,1}
    opacityOmegas              ::Array{Float64,1}
    boundFreeApprox            ::Union{Nothing,Empirical.AbstractEmpiricalApproximation}
    freeFreeIncluded           ::Bool
    resonanceContinuumGrid     ::Union{Nothing,Radial.Grid}
    resonanceAugerSettings     ::AutoIonization.Settings
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme()`  ... constructor for a standard (empty) instance.
"""
function SuperConfigurationMixtureOpacityScheme()
    SuperConfigurationMixtureOpacityScheme( Nuclear.Model[], Array{Array{Configuration,1},1}[], Float64[], Float64[],
                                              nothing, false, nothing, AutoIonization.Settings() )
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                  refConfigsPerSpecies::Array{Array{Configuration,1},1},
                                                  numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1})`
    ... BACKWARD-COMPATIBLE constructor for the ORIGINAL, single-stage-per-species call shape (17-Oct-2026):
        each species' flat `Array{Configuration,1}` is wrapped into a length-1 outer array and
        `boundFreeApprox`/`freeFreeIncluded` are set to `nothing`/`false`, so every pre-existing call site
        keeps compiling and -- since a single-stage species has a TRIVIAL ionization balance
        (`stagePopulations = [1.0]` whatever the temperature or density, there being only one state to
        belong to) AND both continuum contributions are off by default -- keeps returning EXACTLY the
        result it always did, unaffected by the coupled multi-species balance or either continuum addition
        the general constructor now also allows. A scheme::SuperConfigurationMixtureOpacityScheme is returned.
"""
function SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                   refConfigsPerSpecies::Array{Array{Configuration,1},1},
                                                   numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1})
    SuperConfigurationMixtureOpacityScheme( nuclearModels, [ [refs]  for refs in refConfigsPerSpecies ],
                                              numberFractions, opacityOmegas, nothing, false, nothing, AutoIonization.Settings() )
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                  refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                  numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1})`
    ... BACKWARD-COMPATIBLE constructor for the 18-Oct-2026 species+stage call shape WITHOUT a 5th argument
        (`example-Jh.jl` branch b's own call shape): `boundFreeApprox`/`freeFreeIncluded` are set to
        `nothing`/`false`, so every such call site keeps compiling and both continuum contributions stay
        off unless explicitly requested through the full 6-argument form. A
        scheme::SuperConfigurationMixtureOpacityScheme is returned.
"""
function SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                   refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                   numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1})
    SuperConfigurationMixtureOpacityScheme( nuclearModels, refConfigsPerSpeciesStages, numberFractions, opacityOmegas,
                                              nothing, false, nothing, AutoIonization.Settings() )
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                  refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                  numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                  boundFreeApprox::Empirical.AbstractEmpiricalApproximation)`
    ... BACKWARD-COMPATIBLE constructor for the 18-Oct-2026 bound-free-only 5-argument call shape (the one
        `TestFrames.testModule_SuperConfigurationMixtureOpacity`'s own check (8) uses): `freeFreeIncluded` is
        set to `false`, so this call shape keeps compiling and free-free stays off unless explicitly
        requested through the full 6-argument form. A scheme::SuperConfigurationMixtureOpacityScheme is returned.
"""
function SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                   refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                   numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                   boundFreeApprox::Empirical.AbstractEmpiricalApproximation)
    SuperConfigurationMixtureOpacityScheme( nuclearModels, refConfigsPerSpeciesStages, numberFractions, opacityOmegas,
                                              boundFreeApprox, false, nothing, AutoIonization.Settings() )
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                  refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                  numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                  boundFreeApprox::Union{Nothing,Empirical.AbstractEmpiricalApproximation},
                                                  freeFreeIncluded::Bool)`
    ... BACKWARD-COMPATIBLE constructor for the 19-Oct-2026 full 6-argument call shape (the one
        `TestFrames.testModule_SuperConfigurationMixtureOpacity`'s own check (9), free-free, uses): this
        shape USED TO BE the struct's own plain positional constructor, before `resonanceContinuumGrid`/
        `resonanceAugerSettings` were added 20-Oct-2026 and pushed it to 8 fields; this overload keeps it
        compiling, with the resonance treatment OFF (`nothing`, `AutoIonization.Settings()`). A
        scheme::SuperConfigurationMixtureOpacityScheme is returned.
"""
function SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                   refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                   numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                   boundFreeApprox::Union{Nothing,Empirical.AbstractEmpiricalApproximation},
                                                   freeFreeIncluded::Bool)
    SuperConfigurationMixtureOpacityScheme( nuclearModels, refConfigsPerSpeciesStages, numberFractions, opacityOmegas,
                                              boundFreeApprox, freeFreeIncluded, nothing, AutoIonization.Settings() )
end


"""
`Plasma.SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                  refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                  numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                  boundFreeApprox::Union{Nothing,Empirical.AbstractEmpiricalApproximation},
                                                  freeFreeIncluded::Bool, resonanceContinuumGrid::Radial.Grid)`
    ... CONVENIENCE constructor for turning the 20-Oct-2026 resonance treatment ON with the default
        `AutoIonization.Settings()`; equivalent to the full 8-argument struct constructor with
        `resonanceAugerSettings = AutoIonization.Settings()`. A scheme::SuperConfigurationMixtureOpacityScheme
        is returned.
"""
function SuperConfigurationMixtureOpacityScheme(nuclearModels::Array{Nuclear.Model,1},
                                                   refConfigsPerSpeciesStages::Array{Array{Array{Configuration,1},1},1},
                                                   numberFractions::Array{Float64,1}, opacityOmegas::Array{Float64,1},
                                                   boundFreeApprox::Union{Nothing,Empirical.AbstractEmpiricalApproximation},
                                                   freeFreeIncluded::Bool, resonanceContinuumGrid::Radial.Grid)
    SuperConfigurationMixtureOpacityScheme( nuclearModels, refConfigsPerSpeciesStages, numberFractions, opacityOmegas,
                                              boundFreeApprox, freeFreeIncluded, resonanceContinuumGrid, AutoIonization.Settings() )
end


# `Base.show(io::IO, scheme::Plasma.SuperConfigurationMixtureOpacityScheme)`  ... prepares a proper printout.
function Base.show(io::IO, scheme::Plasma.SuperConfigurationMixtureOpacityScheme)
    println(io, "nuclearModels:               $(scheme.nuclearModels)  ")
    println(io, "refConfigsPerSpeciesStages:  $(scheme.refConfigsPerSpeciesStages)  ")
    println(io, "numberFractions:             $(scheme.numberFractions)  ")
    println(io, "opacityOmegas:               $(scheme.opacityOmegas)  ")
    println(io, "boundFreeApprox:             $(scheme.boundFreeApprox)  ")
    println(io, "freeFreeIncluded:            $(scheme.freeFreeIncluded)  ")
    println(io, "resonanceContinuumGrid:      $(scheme.resonanceContinuumGrid)  ")
    println(io, "resonanceAugerSettings:      $(scheme.resonanceAugerSettings)  ")
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
        nuclearM = Nuclear.Model(18., FermiNucleus())
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
include("module-Plasma-inc-average-atom-optical.jl")
include("module-Plasma-inc-line-shifts.jl")
include("module-Plasma-inc-saha-boltzmann-mixture.jl")
include("module-Plasma-inc-satellite-diagnostic.jl")
include("module-Plasma-inc-collisional-radiative.jl")
include("module-Plasma-inc-kinetic-charge-balance.jl")
include("module-Plasma-inc-superconfiguration-opacity.jl")
include("module-Plasma-inc-superconfiguration-ionization.jl")
include("module-Plasma-inc-superconfiguration-mixture.jl")

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
