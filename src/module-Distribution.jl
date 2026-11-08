
"""
`module  JAC.Distribution`  
	... a submodel of JAC that contains all methods for computing different distributions, such as the Maxwell-Boltzmann
        (normalized) velocity and energy distributions, Planck's spectral black-body photon-energy distributions
        and various others.
"""
module Distribution


using  Printf, ..Basics, ..Defaults


"""
`abstract type Distribution.AbstractLineProfile` 
    ... defines an abstract and a number of singleton types to deal with different line profiles, for instance,
        in a plasma.

    + GaussianProfile       ... assumes a Gaussian line profile L^Gaussian (omega)
    + LorentzianProfile     ... assumes a Lorentzian line profile L^Lorentzian (omega)
"""
abstract type  AbstractLineProfile                                   end
struct     GaussianProfile               <:  AbstractLineProfile     end
struct     LorentzianProfile             <:  AbstractLineProfile     end

export  AbstractLineProfile, GaussianProfile, LorentzianProfile

#################################################################################################################################
#################################################################################################################################


"""
`abstract type Distribution.AbstractPhotonDistribution` 
    ... defines an abstract and a number of singleton types to deal with different electron and photon-flux
        distributions in a plasma.

    + PhotonBremsstrahlungThin      ... to apply a free-free bremsstrahlung photon distribution in hot, optically-thin plasma.
    + PhotonBremsstrahlungThick     ... to apply a free-free bremsstrahlung photon distribution in hot, optically-thick plasma.
    + PhotonDilute                  ... to apply a dilute photon distribution.
    + PhotonPlanck                  ... to apply Planck's black-body photon distribution.
    + PhotonPowerLaw                ... to apply a power-law photon distribution.
    + PhotonVacuumField             ... to apply a vacuum (zero-intensity) photon distribution.
"""
abstract type  AbstractPhotonDistribution                            end

    
export  AbstractPhotonDistribution, PhotonBremsstrahlungThin, PhotonBremsstrahlungThick, PhotonDilute, PhotonPlanck,
        PhotonPowerLaw, PhotonVacuumField


"""
`struct  PhotonBremsstrahlungThin  <:  Distribution.AbstractPhotonDistribution`  
    ... to apply a free-free bremsstrahlungs photon spectrum of a hot, optically-thin plasma 
        to photoexcitation, photonionization and recombination processes.

    + T           ::Float64   ... temperature
"""
struct   PhotonBremsstrahlungThin  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
end


function Base.show(io::IO, dist::PhotonBremsstrahlungThin)
    sa = "Free-free bremsstrahlungs distribution of a hot, optically-thin plasma at temperature T = $(dist.T)"
    print(io, sa)
end
        

"""
`struct  PhotonBremsstrahlungThick  <:  Distribution.AbstractPhotonDistribution`  
    ... to apply a free-free bremsstrahlungs photon spectrum of a hot, optically-thick plasma 
        to photoexcitation, photonionization and recombination processes.

    + T           ::Float64   ... temperature
"""
struct   PhotonBremsstrahlungThick  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
end


function Base.show(io::IO, dist::PhotonBremsstrahlungThick)
    sa = "Free-free bremsstrahlungs distribution of a hot, optically-thick plasma at temperature T = $(dist.T)"
    print(io, sa)
end
   

"""
`struct  PhotonDilute  <:  Distribution.AbstractPhotonDistribution`  
    ... to apply a dilute photon spectrum to photoexcitation, photonionization and recombination processes. 
    
    + T           ::Float64   ... temperature
    + w           ::Float64   ... geometric dilution factor
"""
struct   PhotonDilute  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
    w             ::Float64
end


function Base.show(io::IO, dist::PhotonDilute)
    sa = "Dilute photon distribution of a plasma at temperature T = $(dist.T)"
    print(io, sa)
end
   

"""
`struct  PhotonPlanck  <:  Distribution.AbstractPlasmaDistribution`  
    ... to apply Planck's black-body photon spectrum to photoexcitation, photonionization and recombination processes.

    + T           ::Float64   ... temperature
"""
struct   PhotonPlanck  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
end


function Base.show(io::IO, dist::PhotonPlanck)
    sa = "Planck's black-body photon distribution of a hot, optically-thick plasma at temperature T = $(dist.T)"
    print(io, sa)
end
   

"""
`struct  PhotonPowerLaw  <:  Distribution.AbstractPhotonDistribution`  
    ... to apply a dilute photon spectrum n(omega; T) = A * (omega/omega0)^(-p) / ( exp(omega/dist.T) - 1.0 )
        to photoexcitation, photonionization and recombination processes. The constants A, omega0 can be combined but 
        for the price that the units are less obvious.

    + T           ::Float64   ... temperature
    + A           ::Float64   ... Constant amplitude    [a0^-3 Hz^-1]
    + p           ::Float64   ... decay parameter; 1 < p < 2 ... high-energy dominated; p > 2 ... low-energy dominated.
    + omegaMin    ::Float64   ... minimum omega for which the power law applies.
    + omegaMax    ::Float64   ... maximum omega for which the power law applies.
"""
struct   PhotonPowerLaw  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
    A             ::Float64 
    p             ::Float64 
    omegaMin      ::Float64 
    omegaMax      ::Float64 
end


function Base.show(io::IO, dist::PhotonPowerLaw)
    sa = "Power-law photon distribution with parameters A = $(dist.A), p = $(dist.p), omegaMin = $(dist.omegaMin) " * 
         "omegaMax = $(dist.omegaMax) of a plasma at temperature T = $(dist.T)"
    print(io, sa)
    
end
   

"""
`struct  PhotonVacuumField  <:  Distribution.AbstractPhotonDistribution`  
    ... to apply a vacuum (zero-intensity) photon spectrum to photoexcitation, photonionization and 
        recombination processes. This vacuum field is suitable to extract spontaneous rates and 
        rate coefficients.

    + T           ::Float64   ... temperature
"""
struct   PhotonVacuumField  <:  Distribution.AbstractPhotonDistribution
    T             ::Float64
end


function Base.show(io::IO, dist::PhotonVacuumField)
    sa = "Vacuum (zero-intensity) photon distribution, independent of the temperature T = $(dist.T)"
    print(io, sa)
end


#################################################################################################################################
#################################################################################################################################

"""
`abstract type Distribution.AbstractElectronDistribution` 
    ... defines an abstract and a number of singleton types to deal with different electron distributions in a plasma.

    + ElectronBiMaxwell             ... to apply a non-relativistic Bi-Maxwellian electron distribution.
    + ElectronFermiDirac            ... to apply a non-relativistic electron distribution in a degenerate plasma.
    + ElectronMaxwell               ... to apply a non-relativistic Maxwell-Boltzmann electron distribution.
"""
abstract type  AbstractElectronDistribution                            end

    
export  AbstractElectronDistribution, ElectronBiMaxwell, ElectronFermiDirac, ElectronMaxwell

   
   
"""
`struct  ElectronBiMaxwell  <:  Distribution.AbstractElectronDistribution`  
    ... to apply a bi-Maxwellian electron distribution to photoexcitation, photonionization and recombination processes.

    + Tpar        ::Float64   ... parallel temperature
    + Tperp       ::Float64   ... perpendicular temperature
"""
struct   ElectronBiMaxwell   <:  Distribution.AbstractElectronDistribution
    Tpar          ::Float64
    Tperp         ::Float64 
end


function Base.show(io::IO, dist::ElectronBiMaxwell)
    sa = "Bi-Maxwellian electron distribution at the temperatures T_par = $(dist.Tpar) and T_perp = $(dist.Tperp)"
    print(io, sa)
end
   

"""
`struct  ElectronFermiDirac  <:  Distribution.AbstractElectronDistribution`  
    ... to apply a Fermi-Dirac electron distribution in plasma simulations to excitation and recombination processes.

    + T           ::Float64   ... temperature
    + muChem      ::Float64   ... chemical potential mu.
"""
struct   ElectronFermiDirac  <:  Distribution.AbstractElectronDistribution  
    T             ::Float64
    muChem        ::Float64 
end


function Base.show(io::IO, dist::ElectronFermiDirac)
    sa = "Fermi-Dirac electron distribution of a generated plasma at temperature T = $(dist.T)"
    print(io, sa)
end


"""
`struct  ElectronMaxwell  <:  Distribution.AbstractElectronDistribution`  
    ... to apply a Maxwell-Boltzmann electron distribution in plasma simulations to excitation and recombination processes.

    + T           ::Float64   ... temperature
"""
struct   ElectronMaxwell  <:  Distribution.AbstractElectronDistribution 
    T             ::Float64
end


function Base.show(io::IO, dist::ElectronMaxwell)
    sa = "Maxwell-Boltzmann electron distribution of a generated plasma at temperature T = $(dist.T)"
    print(io, sa)
end

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

# Distribution.ElectronMaxwell

"""
`Distribution.electronEnergyDistribution(dist::Distribution.ElectronMaxwell, energy::Float64)`  
    ... to provide the Maxwell-Boltzmann electron energy distribution of given temperature at electron energy.
        A value::Float64 is returned. 
"""
function  electronEnergyDistribution(dist::Distribution.ElectronMaxwell, energy::Float64)
    wa = pi * dist.T^3
    fe = 2 * sqrt(energy / wa) * exp( -energy/dist.T )
    return( fe )
end


"""
`Distribution.electronVelocityDistribution(dist::Distribution.ElectronMaxwell, velocity::Float64)`  
    ... to provide the Maxwell-Boltzmann electron velocity distribution of given temperature at electron velocity.
        A value::Float64 is returned. 
"""
function  electronVelocityDistribution(dist::Distribution.ElectronMaxwell, velocity::Float64)
    wa = 2*pi * dist.T
    fe = (1/wa)^1.5  * exp( - velocity^2 / (2*dist.T) )
    return( fe )
end


#################################################################################################################################
#################################################################################################################################

# Distribution.ElectronFermiDirac

"""
`Distribution.electronEnergyDistribution(dist::Distribution.ElectronFermiDirac, energy::Float64)`  
    ... to provide the Fermi-Dirac electron energy distribution at the given energy, i.e. the spin-degenerate
        free-electron density of states sqrt(2 E)/pi^2 times the Fermi occupation of that energy.  A
        value::Float64 is returned, with the dimension of a number density per unit energy.

        THIS IS AN ABSOLUTE DENSITY AND NOT A PROBABILITY, which is the one thing to know before using it:  it
        integrates to the electron number density n, whereas the two ElectronMaxwell methods of the same names
        integrate to ONE.  The two families are therefore not interchangeable, and the relation between them is
        f_FD -> f_Maxwell * n in the non-degenerate limit, which `TestFrames.testModule_Distribution` asserts.
        A caller wanting a rate coefficient <sigma v> from a degenerate distribution must divide by n itself.

        The prefactor was 2^(3/2)/pi^2 until 08-Nov-2026, i.e. twice this, and so counted the two spin states
        once too often.  Four independent checks agreed on the factor and none of them on the old value:  the
        Jacobian against this module's own velocity form, the textbook density of states
        (2m)^(3/2) sqrt(E) / (2 pi^2 hbar^3), the T -> 0 density (2 mu)^(3/2)/(3 pi^2), and the non-degenerate
        limit above.  Nothing in JAC constructs an ElectronFermiDirac, so no published number moves.
"""
function  electronEnergyDistribution(dist::Distribution.ElectronFermiDirac, energy::Float64)
    wa = (energy - dist.muChem) / dist.T
    fe = sqrt(2.0) / (1*pi^2) * sqrt(energy) / (exp(wa) + 1.0)
    return( fe )
end


"""
`Distribution.electronVelocityDistribution(dist::Distribution.ElectronFermiDirac, velocity::Float64)`  
    ... to provide the Fermi-Dirac velocity distribution of given temperature at electron velocity, with regard to d^3 v,
        i.e. 2/(2 pi)^3 times the Fermi occupation of the state of that speed.  A value::Float64 is returned.

        LIKE ITS ENERGY PARTNER THIS IS AN ABSOLUTE DENSITY, not a probability:  integrated over d^3 v it gives the
        electron number density n, while the two ElectronMaxwell methods integrate to ONE.  Taken through the
        Jacobian f_E(E) dE = f_v(v) 4 pi v^2 dv at E = v^2/2 it yields sqrt(2 E)/pi^2 per unit energy, which is the
        spin-degenerate free-electron density of states, and that agreement with electronEnergyDistribution above
        was established on 08-Nov-2026 by correcting THAT method, which had carried twice this weight.
"""
function  electronVelocityDistribution(dist::Distribution.ElectronFermiDirac, velocity::Float64)
    wa = (velocity^2 / 2.0 - dist.muChem) / dist.T
    fe = 1.0 / (4 * pi^3) / ( exp(wa) + 1.0)
    return( fe )
end


#################################################################################################################################
#################################################################################################################################

# Distribution.PhotonDilute

"""
`Distribution.photonNumberDensity(dist::Distribution.PhotonDilute, omega::Float64)`  
    ... to provide the photon number density of a diluite Planckian black-body field of given temperature at frequency omega.
        A value::Float64 is returned. 
"""
function  photonNumberDensity(dist::Distribution.PhotonDilute, omega::Float64)
    nPhoton = dist.w  * Distribution.photonNumberDensity(Distribution.PhotonPlanck(dist.T), omega)
    return( nPhoton )
end


#################################################################################################################################
#################################################################################################################################

# Distribution.PhotonPlanck

"""
`Distribution.photonNumberDensity(dist::Distribution.PhotonPlanck, omega::Float64)`  
    ... to provide the photon number density of a Planckian black-body field of given temperature at frequency omega.
        A value::Float64 is returned. 
"""
function  photonNumberDensity(dist::Distribution.PhotonPlanck, omega::Float64)
    c = Defaults.getDefaults("speed of light: c")
    nPhoton = omega^2 / (pi^2 * c^3) / ( exp(omega/dist.T) - 1.0 )
    return( nPhoton )
end

#################################################################################################################################
#################################################################################################################################

# Distribution.PhotonPowerLaw

"""
`Distribution.photonNumberDensity(dist::Distribution.PhotonPowerLaw, omega::Float64)`  
    ... to provide the photon number density of a dilute photon field of given temperature at frequency omega.
        A value::Float64 is returned. 
"""
function  photonNumberDensity(dist::Distribution.PhotonPowerLaw, omega::Float64)
    nPhoton = 0.
    if  dist.omegaMin  <= omega  <=  dist.omegaMax
        nPhoton = dist.A * omega^(-dist.p) * exp(- omega/dist.T)
    end
    return( nPhoton )
end

#################################################################################################################################
#################################################################################################################################

# Distribution.PhotonVacuumField

"""
`Distribution.photonNumberDensity(dist::Distribution.PhotonVacuumField, omega::Float64)`  
    ... to provide the photon number density of a vacuum photon field , i.e. a zero field.
        A value::Float64 is returned. 
"""
function  photonNumberDensity(dist::Distribution.PhotonVacuumField, omega::Float64)
    return( 0. )
end


end # module
