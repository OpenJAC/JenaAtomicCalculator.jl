
"""
`module  JAC.Empirical`  
... a submodel of JAC that contains all methods to set-up and process (simple) empirical computations such as 
    electron-impact ionization, charge exchange, etc.
"""
module Empirical


using  ..AtomicState, ..Basics, ..Continuum, ..Defaults, ..Distribution, ..Radial, ..ManyElectron, ..Nuclear,
       ..InteractionStrength, ..ImpactIonization, ..PeriodicTable, ..Semiempirical, ..SelfConsistent, SpecialFunctions

"""
`abstract type Empirical.AbstractEmpiricalApproximation`
    ... defines an abstract and a number of singleton types to distinguish the empirical approximations that are
        applied for estimating cross sections, rates and (plasma) rate coefficients within this module. Which
        approximations are supported by a given function follows from its methods (multiple dispatch); an
        unsupported combination raises an informative error.

    + ScaledHydrogenic     ... hydrogenic (Kramers-type) formulas, scaled by tabulated binding energies;
                               cf. Empirical.scaledBindingEnergy.
    + UsingJAC             ... simple mean-field (DFS) computations of energies and one-electron amplitudes.
    + GivenEinsteinA       ... an externally given Einstein-A value (multipole, energy, rate); cf. below.
    + Lotz1967             ... Lotz's (1967) formula for electron-impact ionization cross sections
                               [W. Lotz, Z. Physik 206, 205 (1967)].
    + VanRegemorter1962    ... Van Regemorter's (1962) formula for electron-impact excitation of optically
                               allowed transitions [H. Van Regemorter, ApJ 136, 906 (1962)].
    + Bohr1913             ... Bohr's classical stopping power of a plasma, with the Coulomb logarithm
                               ln(1.123 m v^3 / (e^2 omega_p))  [N. Bohr, Philos. Mag. 25, 10 (1913)].
    + Bethe1931            ... Bethe's quantum stopping power of a plasma, with the Coulomb logarithm
                               ln(2 E / (hbar omega_p)) for an electron projectile
                               [H. Bethe, Ann. Phys. (Leipzig) 397, 325 (1930)].
    + KozmaFranson1992     ... the piecewise electron loss function of Kozma & Fransson, quantal above and
                               classical below E = 14 eV [C. Kozma & C. Fransson, ApJ 390, 602 (1992), Eqs. (1-3);
                               after R. Schunk & P. Hays, Planet. Space Sci. 19, 113 (1971)].
    + Axelrod1980          ... Axelrod's relativistic plasma energy loss, i.e. the Bethe-type loss with the plasma
                               energy hbar omega_p as effective ionization potential [T. Axelrod, PhD thesis, UC
                               Santa Cruz (1980); as transcribed by P. Milne et al., ApJS 124, 503 (1999), Eqs. (1,3)].
    + ADK1986              ... the (quasiclassical) tunneling ionization rate of Ammosov, Delone & Krainov for an
                               atom or ion in a strong, quasi-static electric field [M. Ammosov, N. Delone &
                               V. Krainov, Sov. Phys. JETP 64, 1191 (1986); as summarized, e.g., by J. Bauer &
                               P. Mulser, Phys. Rev. A 59, 569 (1999) [arXiv:physics/9802042], Eq. (10)].
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

struct     ADK1986                       <:  AbstractEmpiricalApproximation  end


export  AbstractEmpiricalApproximation, GivenEinsteinA, ScaledHydrogenic, UsingJAC,
        Bohr1913, Bethe1931, Axelrod1980, KozmaFranson1992, Lotz1967, VanRegemorter1962, ADK1986


"""
`abstract type Empirical.AbstractStoppingProjectile`
    ... defines an abstract and a number of types for the projectile whose stopping power is requested. The projectile
        enters the stopping formulas only through its charge z (as z^2 in the prefactor), its mass M (via the velocity
        v = sqrt(2 E/M) for a given kinetic energy E) and its quantum statistics: an electron projectile is
        indistinguishable from the target electrons, which limits the maximum energy transfer to E/2 and reduces the
        Coulomb logarithm by ln 2 relative to a heavy projectile of the same velocity.

    + ElectronProjectile   ... an electron (z = -1, M = 1).
    + PositronProjectile   ... a positron (z = +1, M = 1); it is distinguishable from the target electrons but is
                               treated with the same Coulomb logarithms as an electron -- the differences appear
                               in the relativistic corrections and are beyond the accuracy of these estimates
                               (cf. Milne et al. 1999, who sum electron and positron losses in just this way).
    + IonProjectile        ... a bare (structureless) ion of charge z and mass M [in units of m_e]; for a proton,
                               IonProjectile(1.0, Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U). Electron capture
                               and loss of a dressed ion are not accounted for.
"""
abstract type  AbstractStoppingProjectile                                    end
struct     ElectronProjectile            <:  AbstractStoppingProjectile      end
struct     PositronProjectile            <:  AbstractStoppingProjectile      end


"""
`struct  Empirical.IonProjectile  <:  AbstractStoppingProjectile`
    ... a bare (structureless) ion as projectile of a stopping-power calculation.

    + z           ::Float64   ... charge of the ion [in units of e].
    + M           ::Float64   ... mass of the ion [in units of m_e].
"""
struct  IonProjectile  <:  AbstractStoppingProjectile
    z             ::Float64
    M             ::Float64
end


# `Base.show(io::IO, p::ElectronProjectile)`  ... provides a String notation for the variable p::ElectronProjectile.
function Base.show(io::IO, p::ElectronProjectile)
    print(io, "an electron (z = -1, M = m_e)")
end


# `Base.show(io::IO, p::PositronProjectile)`  ... provides a String notation for the variable p::PositronProjectile.
function Base.show(io::IO, p::PositronProjectile)
    print(io, "a positron (z = +1, M = m_e)")
end


# `Base.show(io::IO, p::IonProjectile)`  ... provides a String notation for the variable p::IonProjectile.
function Base.show(io::IO, p::IonProjectile)
    print(io, "a bare ion (z = $(p.z), M = $(p.M) m_e)")
end


"""
`abstract type Empirical.AbstractStoppingMaterial`
    ... defines an abstract and a number of types for the material in which the projectile is stopped. In the
        (Bohr/Bethe-type) stopping formulas of this module, a material enters only through the density of the
        electrons that take up the energy and through one characteristic energy of these electrons: the plasma
        energy hbar omega_p = hbar sqrt(4 pi n_e e^2/m_e) for free electrons, and the mean excitation energy
        Ibar for bound electrons. Nuclear stopping and radiative (bremsstrahlung) losses belong to different
        mechanisms and are *not* included in any of these materials.

    + FreeElectronGas      ... the free electrons of a (fully ionized) plasma with electron density ne.
    + NeutralAtomGas       ... the bound electrons of a neutral, monoatomic gas of atoms (Z, A) with number
                               density natom; the mean excitation energy is estimated as
                               Ibar = 9.1 Z (1 + 1.9 Z^(-2/3)) eV (Segre 1977; Roy & Reed 1968).
    + PartiallyIonizedGas  ... a monoatomic gas of atoms (Z, A) with number density natom, of which chie electrons
                               per atom are ionized: (Z - chie) natom bound electrons with Ibar plus chie natom
                               free electrons with hbar omega_p; cf. Milne et al. (1999, Eq. 3).
"""
abstract type  AbstractStoppingMaterial                                      end


"""
`struct  Empirical.FreeElectronGas  <:  AbstractStoppingMaterial`
    ... the free electrons of a (fully ionized) plasma.

    + ne          ::Float64   ... electron density [a.u.].
"""
struct  FreeElectronGas  <:  AbstractStoppingMaterial
    ne            ::Float64
end


"""
`struct  Empirical.NeutralAtomGas  <:  AbstractStoppingMaterial`
    ... the bound electrons of a neutral, monoatomic gas.

    + Z           ::Int64     ... nuclear charge of the atoms.
    + A           ::Float64   ... atomic mass number of the atoms; only used to express the mass stopping power.
    + natom       ::Float64   ... atom density [a.u.].
"""
struct  NeutralAtomGas  <:  AbstractStoppingMaterial
    Z             ::Int64
    A             ::Float64
    natom         ::Float64
end


"""
`struct  Empirical.PartiallyIonizedGas  <:  AbstractStoppingMaterial`
    ... a partially ionized, monoatomic gas: (Z - chie) bound electrons per atom plus chie free electrons per atom.

    + Z           ::Int64     ... nuclear charge of the atoms.
    + A           ::Float64   ... atomic mass number of the atoms; only used to express the mass stopping power.
    + chie        ::Float64   ... number of free electrons per atom (ionization fraction), 0 <= chie <= Z.
    + natom       ::Float64   ... atom density [a.u.].
"""
struct  PartiallyIonizedGas  <:  AbstractStoppingMaterial
    Z             ::Int64
    A             ::Float64
    chie          ::Float64
    natom         ::Float64
end


# `Base.show(io::IO, m::FreeElectronGas)`  ... provides a String notation for the variable m::FreeElectronGas.
function Base.show(io::IO, m::FreeElectronGas)
    print(io, "the free electrons of a plasma with n_e = $(m.ne) [a.u.]")
end


# `Base.show(io::IO, m::NeutralAtomGas)`  ... provides a String notation for the variable m::NeutralAtomGas.
function Base.show(io::IO, m::NeutralAtomGas)
    print(io, "the bound electrons of a neutral atom gas with Z = $(m.Z), A = $(m.A) and n_atom = $(m.natom) [a.u.]")
end


# `Base.show(io::IO, m::PartiallyIonizedGas)`  ... provides a String notation for the variable m::PartiallyIonizedGas.
function Base.show(io::IO, m::PartiallyIonizedGas)
    print(io, "a partially ionized atom gas with Z = $(m.Z), A = $(m.A), chi_e = $(m.chie) free electrons per atom " *
              "and n_atom = $(m.natom) [a.u.]")
end


export  AbstractStoppingProjectile, ElectronProjectile, PositronProjectile, IonProjectile,
        AbstractStoppingMaterial, FreeElectronGas, NeutralAtomGas, PartiallyIonizedGas
   

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
include("module-Empirical-inc-tunneling-ionization.jl")

end # module
