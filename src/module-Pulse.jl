
"""
`module  JAC.Pulse`  
... a submodel of JAC that contains all structs and methods to deal with time-dependent pulses of the em field.
"""
module Pulse


using   ..Basics, ..Defaults, ..Radial, GSL, SpecialFunctions

export  AbstractEnvelope, AbstractBeam, AbstractPulse


"""
`abstract type Pulse.AbstractEnvelope` 
    ... defines an abstract type to comprise various envelopes of (possible) laser pulses in terms of their shape,  
        pulse duration or number of cycles, etc.

    + InfiniteEnvelope         ... to represent an infinte (plane-wave) pulse.
    + RectangularEnvelope      ... to represent a finite rectangular pulse.
    + SinSquaredEnvelope       ... to represent a finite sin^2 pulse.
    + GaussianEnvelope         ... to represent a Gaussian light pulse.
"""
abstract type  AbstractEnvelope  end


"""
`struct  Pulse.InfiniteEnvelope  <: Pulse.AbstractEnvelope`   ... to represent an infinte (plane-wave) pulse.
"""
struct   InfiniteEnvelope        <: Pulse.AbstractEnvelope    end

function Base.string(env::InfiniteEnvelope)
    sa = "infinite pulse"
    return( sa )
end

function Base.show(io::IO, env::InfiniteEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.RectangularEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a finte rectangular pulse.

    + cycles      ::Int64     ... Number of cycles of the pulse.
"""
struct   RectangularEnvelope       <: Pulse.AbstractEnvelope
    cycles        ::Int64
end

function Base.string(env::RectangularEnvelope)
    sa = "rectangular pulse of $(env.cycles) cycles"
    return( sa )
end

function Base.show(io::IO, env::RectangularEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.SinSquaredEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a finite sin^2 pulse.

    + cycles      ::Int64     ... Number of cycles of the pulse.
"""
struct   SinSquaredEnvelope        <: Pulse.AbstractEnvelope
    cycles      ::Int64
end


function Base.string(env::SinSquaredEnvelope)
    sa = "sin^2 pulse of $(env.cycles) cycles"
    return( sa )
end

function Base.show(io::IO, env::SinSquaredEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`struct  Pulse.GaussianEnvelope  <: Pulse.AbstractEnvelope`   ... to represent a (infinite) Gaussian pulse.

    + fwhm        ::Float64     ... FWHM which is often taken as pulse duration
"""
struct   GaussianEnvelope        <: Pulse.AbstractEnvelope
    fwhm          ::Float64
end


function Base.string(env::GaussianEnvelope)
    sa = "Gaussian pulse with pulse duration or FWHM = $(env.fwhm)"
    return( sa )
end

function Base.show(io::IO, env::GaussianEnvelope)
    sa = string(env);       print(io, sa, "\n")
end


"""
`abstract type Pulse.AbstractBeam` 
    ... defines an abstract type to comprise various basic laser pulses that are characterized in terms of their amplitude, 
        frequency, carrier-envelope phase, etc. In general, the basic beam properties are independent of the (pulse) envelope 
        and the polarization properties which are handled and communicated separately (to and within the program).

    + PlaneWaveBeam            ... to represent a plane-wave beam.
    + BesselBeam               ... to represent a Bessel beam (not yet).
"""
abstract type  AbstractBeam  end


"""
`struct  Pulse.PlaneWaveBeam  <: AbstractBeam`   
    ... to represent a plane-wave beam with given (real) amplitude, frequency and carrier-envelope phase.

    + A0            ::Float64            ... (Constant) Amplitude of the light pulse. 
    + omega         ::Float64            ... Central frequency. 
    + cep           ::Float64            ... Carrier-envelope phase. 
"""
struct   PlaneWaveBeam  <: AbstractBeam
    A0              ::Float64
    omega           ::Float64
    cep             ::Float64
end


"""
`Pulse.PlaneWaveBeam()`  ... constructor for an `empty` instance of Pulse.PlaneWaveBeam().
"""
function PlaneWaveBeam()
    PlaneWaveBeam( 0., 0., 0. )
end


function Base.string(beam::PlaneWaveBeam)
    sa = "Plane-wave beam/pulse with amplitude A0=$(beam.A0), frequency omega=$(beam.omega) a.u. and carrier-envelope phase cep=$(beam.cep)"
    return( sa )
end

function Base.show(io::IO, beam::PlaneWaveBeam)
    sa = string(beam);       print(io, sa, "\n")
end


"""
`abstract type Pulse.AbstractPulseParameter` 
    ... defines an abstract type to deal with the (experimental) specification of light pulses of type Pulse.GeneratedPulse.
        These (experimental) GeneratedPulse's are always converted into (computational) Gaussian, SinSquared, ... pulses
        before they are used in Liouville time-evolutions and simulations. These parameters are usually defined as singletons
        but can provide also a few additional parameters (like the ellipticity of elliptical-polarized pulse).
        The list of these parameters can be readily extended, whereas the conversion from the GeneratedPulse's to the 
        computational pulses is usually nontrivial and need to be done step-by-step.

    + GaussianShape            ... to represent a Gaussian-shaped pulse with a given width.
    + SinSquaredShape          ... to represent a sin^2-shaped pulse of given period (a.u.).
    + LinearPolarized          ... to represent a linearly-polarized pulse.
    + EllipticalPolarized      ... to represent a elliptically-polarized pulse with given ellipticity.
    + LeftCircular             ... to represent a left-circular polarized pulse.
    + RightCircular            ... to represent a right-circular polarized pulse.
"""
abstract type  AbstractPulseParameter  end
struct   GaussianShape          <:  Pulse.AbstractPulseParameter   end
struct   SinSquaredShape        <:  Pulse.AbstractPulseParameter   end
struct   LinearPolarized        <:  Pulse.AbstractPulseParameter   end
struct   EllipticalPolarized    <:  Pulse.AbstractPulseParameter   end
struct   LeftCircular           <:  Pulse.AbstractPulseParameter   end
struct   RightCircular          <:  Pulse.AbstractPulseParameter   end


"""
`abstract type Pulse.AbstractPulse` 
    ... defines an abstract type to comprise different full laser pulses that are characterized in terms of their amplitude, 
        frequency, carrier-envelope phase, polarization as well their occurrence time (with regard to a reference time refTime).
        Here, we shall distinguish between an GeneratedPulse that is specified in terms of experimentally accessible parameters
        and specifications as well as various pulses (Gaussian, ...) that are used in computations. Each GeneratedPulse need
        first to be converted into the internal (computational) format before it can be applied in simulations.

    + FelPulse                 ... to represent an FEL pulse in experimentally simple terms and parameters.
    + GeneratedPulse           ... to represent a pulse in experimentally accessible terms and parameters.
    + Gaussian                 ... to represent a Gaussian pulse
    + SinSquared               ... to represent a SinSquared pulse
"""
abstract type  AbstractPulse  end


"""
`struct  Pulse.FelPulse   <:  Pulse.AbstractPulse`  
    ... defines an FEL pulse with experimentally accessible properties; these pulses need first to the converted into a 
        proper computational pulse before they can be utilized in Liouville time evolutions.

    + shape            ::String        ... Shape (GaussianSimplified) of the pulse.
    + omega            ::Float64       ... Central frequency [eV] of the em pulse.
    + intensity        ::Float64       ... Peak intensity [W/cm^2] of the pulse.
    + fwhm             ::Float64       ... full-width-half-maximum in time [fs] of pulse. 
    + timeDelay        ::Float64       ... Propagation time [fs] of the pulse before its peak arrives at the atomic target.
    + energySpread     ::Float64       ... Energy spread [eV] of the em pulse.
"""
struct FelPulse   <:  Pulse.AbstractPulse
    shape              ::String
    omega              ::Float64 
    intensity          ::Float64 
    fwhm               ::Float64 
    timeDelay          ::Float64
    energySpread       ::Float64
end 


"""
`Pulse.FelPulse()`  ... constructor for an `empty` instance of Pulse.FelPulse().
"""
function FelPulse()
    FelPulse( "NoShape", 0., 0., 0., 0., 0. )
end


# `Base.show(io::IO, pulse::Pulse.FelPulse)`  ... prepares a proper printout of the variable pulse::Pulse.FelPulse.
function Base.show(io::IO, pulse::Pulse.FelPulse) 
    println(io, "shape:         $(pulse.shape)  ")
    println(io, "omega:         $(pulse.omega)  ")
    println(io, "intensity:     $(pulse.intensity)  ")
    println(io, "fwhm:          $(pulse.fwhm)       ")
    println(io, "timeDelay:     $(pulse.timeDelay)  ")
    println(io, "energySpread:  $(pulse.energySpread)  ")
end


"""
`struct  Pulse.GeneratedPulse   <:  Pulse.AbstractPulse`  
    ... defines a (full) GeneratedPulse light pulse in terms of experimentally accessible parameters, keywords, etc.

    + propagation      ::SolidAngle                     ... Propagation direction of the pulse as described by the unit vector u = u(Omega).
    + omega            ::Float64                        ... Central frequency of the em pulse (in atomic units)
    + intensity        ::Float64                        ... Peak intensity of the em pulse (in W/cm^2).
    + timeDelay        ::Float64                        
        ... Initial time-delay with regard to the reference time (a.u.); the detailed meaning of this time-delay depends
            on the shape of the pulse.
    + multipoles       ::Array{EmMultipole,1}           ... Multipoles of the em field to be included in the description of the light field.
    + shape            ::Pulse.AbstractPulseParameter   ... to characterize the shape/envelope of the light pulse.
    + polarization     ::Pulse.AbstractPulseParameter   ... to characterize the polarization of the light pulse.
"""
struct GeneratedPulse   <:  Pulse.AbstractPulse
    propagation        ::SolidAngle 
    omega              ::Float64 
    intensity          ::Float64 
    timeDelay          ::Float64
    multipoles         ::Array{EmMultipole,1}
    shape              ::Pulse.AbstractPulseParameter 
    polarization       ::Pulse.AbstractPulseParameter
end 


"""
`Pulse.GeneratedPulse()`  ... constructor for an `empty` instance of Pulse.GeneratedPulse().
"""
function GeneratedPulse()
    GeneratedPulse( SolidAngle(0., 0.), 0., 0., 0., EmMultipole[], Pulse.GaussianShape, Pulse.LinearPolarized  )
end


# `Base.show(io::IO, pulse::Pulse.GeneratedPulse)`  ... prepares a proper printout of the variable pulse::Pulse.GeneratedPulse.
function Base.show(io::IO, pulse::Pulse.GeneratedPulse) 
    println(io, "propagation:        $(pulse.propagation)  ")
    println(io, "omega:              $(pulse.omega)  ")
    println(io, "intensity:          $(pulse.intensity)  ")
    println(io, "timeDelay:          $(pulse.timeDelay)  ")
    println(io, "multipoles:         $(pulse.multipoles)  ")
    println(io, "shape:              $(pulse.shape)  ")
    println(io, "polarization:       $(pulse.polarization)  ")
end


"""
`struct  Pulse.Gaussian   <:  Pulse.AbstractPulse`  
    ... defines a (full) Gaussian light pulse with well-defined frequency, carrier-envelope phase, polarization as well their 
        occurrence time (with regard to a reference time refTime).

    + propagation      ::SolidAngle                 ... Propagation direction of the pulse as described by the unit vector u = u(Omega).
    + omega            ::Float64                    ... Central frequency of the em pulse.
    + timeDelay        ::Float64                    ... Initial time-delay of the peak center with regard to the reference time (a.u.).
    + multipoles       ::Array{EmMultipole,1}       ... Multipoles of the em field to be included in the description of the light field.
    + envelope         ::Pulse.GaussianEnvelope     ... Gaussian envelope (function) of the light pulse.
    + polarization     ::Pulse.AbstractPolarization ... Polarization of the light pulse with typically well-defined 
                                                        g_+1 and g_-1 coefficients; not yet.
"""
struct Gaussian   <:  Pulse.AbstractPulse
    propagation        ::SolidAngle 
    omega              ::Float64 
    timeDelay          ::Float64 
    multipoles         ::Array{EmMultipole,1}
    envelope           ::Pulse.GaussianEnvelope
    # polarization     ::Pulse.AbstractPolarization
end 


"""
`Pulse.Gaussian()`  ... constructor for an `empty` instance of Pulse.Gaussian().
"""
function Gaussian()
    Gaussian( SolidAngle(0., 0.), 0., 0., EmMultipole[], Pulse.GaussianEnvelope(0.)  )
end


# `Base.show(io::IO, pulse::Pulse.Gaussian)`  ... prepares a proper printout of the variable pulse::Pulse.Gaussian.
function Base.show(io::IO, pulse::Pulse.Gaussian) 
    println(io, "propagation:        $(pulse.propagation)  ")
    println(io, "omega:              $(pulse.omega)  ")
    println(io, "timeDelay:          $(pulse.timeDelay)  ")
    println(io, "multipoles:         $(pulse.multipoles)  ")
    println(io, "envelope:           $(pulse.envelope)  ")
end


"""
`struct  Pulse.GaussianSimplified   <:  Pulse.AbstractPulse`  
    ... defines a full but simplified and linearly-polarized Gaussian light pulse with well-defined frequency, 
        peak amplitude, fwhm as well as energy spread, such as they occur at FEL sources. These simplified pulses 
        always propagate along the z-axis and are polarized along the x-axis. The full pulse is internal 
        restricted to 5*fwhm.

    + omega            ::Float64       ... Central frequency (a.u.) of the em pulse.
    + A0               ::Float64       ... peak amplitude of the pulse.
    + fwhm             ::Float64       ... full-width-half-maximum in time (a.u.) of the Gaussian pulse. 
    + timeDelay        ::Float64       ... Propagation time (a.u.) of the pulse before its peak arrives at the atomic target.
    + energySpread     ::Float64       ... Energy spread (a.u.) of the em pulse.
"""
struct GaussianSimplified   <:  Pulse.AbstractPulse
    omega              ::Float64 
    A0                 ::Float64 
    fwhm               ::Float64 
    timeDelay          ::Float64
    energySpread       ::Float64
end 


"""
`Pulse.GaussianSimplified()`  ... constructor for an `empty` instance of Pulse.GaussianSimplified().
"""
function GaussianSimplified()
    Gaussian( 0., 0., 0., 0., 0.)
end


# `Base.show(io::IO, pulse::Pulse.GaussianSimplified)`  ... prepares a proper printout of the variable pulse::Pulse.GaussianSimplified.
function Base.show(io::IO, pulse::Pulse.GaussianSimplified) 
    println(io, "omega:         $(pulse.omega)  ")
    println(io, "A0:            $(pulse.A0)  ")
    println(io, "fwhm:          $(pulse.fwhm)  ")
    println(io, "timeDelay:     $(pulse.timeDelay)  ")
    println(io, "energySpread:  $(pulse.energySpread)  ")
end


"""
`struct  Pulse.SinSquared   <:  Pulse.AbstractPulse`  
    ... defines a (full) sin^2 light pulse with well-defined frequency, carrier-envelope phase, polarization as well their 
        occurrence time (with regard to a reference time refTime).

    + propagation      ::SolidAngle                 ... Propagation direction of the pulse as described by the unit vector u = u(Omega).
    + omega            ::Float64                    ... Central frequency of the em pulse.
    + timeDelay        ::Float64                    ... Initial time-delay of the peak center with regard to the reference time (a.u.).
    + multipoles       ::Array{EmMultipole,1}       ... Multipoles of the em field to be included in the description of the light field.
    + envelope         ::Pulse.SinSquaredEnvelope   ... Gaussian envelope (function) of the light pulse.
    + polarization     ::Pulse.AbstractPolarization ... Polarization of the light pulse with typically well-defined 
                                                        g_+1 and g_-1 coefficients; not yet.
"""
struct SinSquared   <:  Pulse.AbstractPulse
    propagation        ::SolidAngle 
    omega              ::Float64 
    timeDelay          ::Float64 
    multipoles         ::Array{EmMultipole,1}
    envelope           ::Pulse.SinSquaredEnvelope
    # polarization     ::Pulse.AbstractPolarization
end 


"""
`Pulse.SinSquared()`  ... constructor for an `empty` instance of Pulse.SinSquared().
"""
function SinSquared()
    SinSquared( SolidAngle(0., 0.), 0., 0., EmMultipole[], Pulse.SinSquaredEnvelope(0.)  )
end


# `Base.show(io::IO, pulse::Pulse.SinSquared)`  ... prepares a proper printout of the variable pulse::Pulse.SinSquared.
function Base.show(io::IO, pulse::Pulse.SinSquared) 
    println(io, "propagation:        $(pulse.propagation)  ")
    println(io, "omega:              $(pulse.omega)  ")
    println(io, "timeDelay:          $(pulse.timeDelay)  ")
    println(io, "multipoles:         $(pulse.multipoles)  ")
    println(io, "envelope:           $(pulse.envelope)  ")
end


###################################################################################################################
###################################################################################################################
###################################################################################################################



"""
`Pulse.computeFieldAmplitude(intensity::Float64, omega::Float64)`  
    ... compute the field amplitude from the given (maximum) intensity [in a.u.] and the central frequency [in a.u.]
        of the light field: A0 = sqrt( 8pi * alpha * intensity) / omega.
"""
function computeFieldAmplitude(intensity::Float64, omega::Float64)
    wa = sqrt(8 * pi * Defaults.getDefaults("alpha") * intensity) / omega
    return( wa )
end


"""
`Pulse.convertPulse(expPulse::Pulse.FelPulse)`  
    ... converts the experimentally specified FEL pulse into a proper computational pulse that can be utilized for 
        the Liouville time evolution.
"""
function convertPulse(expPulse::Pulse.FelPulse)
    
    omega        = Defaults.convertUnits("energy: from eV to atomic",        expPulse.omega)
    energySpread = Defaults.convertUnits("energy: from eV to atomic",        expPulse.energySpread)
    A0           = Defaults.convertUnits("intensity: from W/cm^2 to atomic", expPulse.intensity)
    fwhm         = expPulse.fwhm      * 1.0e-15 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    timeDelay    = expPulse.timeDelay * 1.0e-15 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    
    if    expPulse.shape == "GaussianSimplified"     
        compPulse = Pulse.GaussianSimplified(omega, A0, fwhm, timeDelay, energySpread);
    else  error("Unknown pulse shape = $(expPulse.shape) ")
    end 
    
    return( compPulse )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(infinite); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;                            
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -10:10
        if  plus   wb = Basics.diracDelta((s-1)*beam.omega + energyp - initialEn + Up, 1.0e-3)
        else       wb = Basics.diracDelta((s+1)*beam.omega + energyp - initialEn + Up, 1.0e-3)       
        end
        #
        if  wb != 0.
            wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * (phiCep - lambda*phip)) * wb
        end
    end
    if  plus   wa = wa * 2pi * beam.A0 * exp(-im*phiCep)
    else       wa = wa * 2pi * beam.A0 * exp(im*phiCep)
    end
end
    
    return( wa )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   Tp = 2pi * envelope.cycles / beam.omega
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / ( sqrt(2) * beam.omega );   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -20:20
        if  plus   wb = (s-1)*beam.omega + energyp - initialEn + Up
        else       wb = (s+1)*beam.omega + energyp - initialEn + Up       
        end
        #
        if  wb != 0.
            wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * phaseConstant) / wb * ( exp(im*wb*Tp) - 1 )
        end
    end
    if  plus   wa = wa * (-im) * exp(-im * a * sin(phaseConstant)) * beam.A0 * exp(-im*phiCep)
    else       wa = wa * (-im) * exp(-im * a * sin(phaseConstant)) * beam.A0 * exp(im*phiCep)
    end
    end
    
    return( wa )
end


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a sine-squared pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   np = envelope.cycles;   Tp = 2pi * np / beam.omega
    omega = beam.omega
    phiCep = beam.cep;
    sinSqrArg = 0.5 * omega / np;
    lambda = Basics.determinePolarizationLambda(polarization)
    
    p = sqrt(2.0*energyp)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = 1.0
    if polarization != Basics.RightCircular() && polarization != Basics.LeftCircular()
        epsilon = polarization.ellipticity
    end
    
    A0eps = beam.A0/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 100*np
    else             orderGL = 1000
    end
    gaussLegendre = Radial.GridGL(Radial.GridGaussLegendreFinite(),0.0,Tp,orderGL)
    tgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = energyp*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = beam.A0 * exp(-im * phiCep) * wa
    else        wa = beam.A0 * exp(im * phiCep) * wa
    end
    end
    
    return( wa )
end    


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int )
if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    #Define Gauss-Hermite grid, convergence is typically good for orderGH = 10000
    orderGH = 10000
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = transformFac * beam.A0 * exp(-im * phiCep) * wa
    else        wa = transformFac * beam.A0 * exp(im * phiCep) * wa
    end
end
    
    return( wa )
end  


"""
`Pulse.pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderGH::Int64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA) [+/-; omega; f^(rectangular); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeIntegral(plus::Bool, envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderGH::Int64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        if  plus    integrand = exp( -im * ( ( initialEn + beam.omega ) * t - SVolkov ) )
        else        integrand = exp( -im * ( ( initialEn - beam.omega ) * t - SVolkov ) )
        end
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    if  plus    wa = transformFac * beam.A0 * exp(-im * phiCep) * wa
    else        wa = transformFac * beam.A0 * exp(im * phiCep) * wa
    end
end
    
    return( wa )
end  


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for an infinite pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.InfiniteEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
if orderSFA == 0
    wa = 0. * im;                            
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    #
    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -10:10
        wb = Basics.diracDelta(s*beam.omega + energyp - initialEn + Up, 1.0e-3)
        #
        if  wb != 0.    wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * (phiCep - lambda*phip)) * wb     
        end
    end
    wa = wa * 4pi * Up
end
return( wa )
end


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies]  for a rectangular pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.RectangularEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   Tp = 2pi * envelope.cycles / beam.omega
    phiCep = beam.cep;       a = beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) / beam.omega;   Up = beam.A0^2 / 4
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip

    # Compute the summation over the Bessel functions first; start with value for s = 0
    for  s = -20:20
        wb = s*beam.omega + energyp - initialEn + Up
        #
        if  wb != 0.    wa = wa + GSL.sf_bessel_Jn(s, a) * exp(im*s * phaseConstant) / wb * ( exp(im*wb*Tp) - 1 )   
        end
    end
    wa = wa * (-im) * 2 * Up * exp(-im * a * sin(phaseConstant))
    end
    
    return( wa )
end


"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for a sine-squared pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.SinSquaredEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   np = envelope.cycles;   Tp = 2pi * np / beam.omega
    omega = beam.omega
    phiCep = beam.cep
    sinSqrArg = 0.5 * omega / np
    lambda = Basics.determinePolarizationLambda(polarization)
    
    p = sqrt(2.0*energyp)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = 1.0
    if polarization != Basics.RightCircular() && polarization != Basics.LeftCircular()
        epsilon = polarization.ellipticity
    end
    
    A0eps = beam.A0/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 100*np
    else             orderGL = 1000
    end
    gaussLegendre = Radial.GridGL(Radial.GridGaussLegendreFinite(),0.0,Tp,orderGL)
    tgrid = gaussLegendre.t
    weights = gaussLegendre.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = energyp*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        integrand = sin( sinSqrArg * t )^4 * ( cos(omega*t+phiCep)^2 + epsilon^2 * sin(omega*t+phiCep)^2 ) * exp( -im * ( initialEn * t - SVolkov ) )
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    wa = A0eps^2 * wa
    end
    
    return( wa )
end



"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)` 
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies] for a gaussian pulse with given 
        parameters; an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.GaussianEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;  Tp = envelope.fwhm
    phiCep = beam.cep;   Up = beam.A0^2 / 4
    a1 = 0.25 * Up * sqrt(pi/log(4)) * Tp
    a2 = 0.125 * Tp * sqrt(pi/log(2)) * beam.A0 * sqrt(2*energyp) * sin(thetap) / sqrt(2) * exp( -beam.omega^2 * Tp^2 / log( 65536 ) )
    lambda = Basics.determinePolarizationLambda(polarization)
    phaseConstant = phiCep - lambda*phip
    
    transformFac = Tp/(2*sqrt(2*log(2))) #factor from (linear) variable transformation tau -> t in order to allow Gauss-Hermite integration ( weight function e^(-t^2) )
    
    #Define Gauss-Hermite grid, convergence is typically good for ......
    orderGH = 10000
    gaussHermite = Radial.GridGH(orderGH)
    tgrid = gaussHermite.t
    weights = gaussHermite.wt
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGH
        t = transformFac * tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        SVolkov = energyp * t + a1 * erf(2 * sqrt(log(4)) * t/Tp) 
                    + a2 * ( exp( -im*phaseConstant ) * erf( (im*Tp^2*beam.omega + 8*log(2)*t)/(4*Tp*sqrt(log(2))) )
                                - exp( im*phaseConstant ) * erf( (im*Tp^2*beam.omega - 8*log(2)*t)/(4*Tp*sqrt(log(2))) ) )
        
        #Compute integrand at gridpoint t
        integrand = exp( -im * ( initialEn * t - SVolkov ) )
        
        #Gauss-Hermite sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    wa = 0.5 * transformFac * beam.A0^2 * wa
    end
    
    return( wa )
end  

"""
`Pulse.pulseShapeQuadIntegral(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)`  
    ... evaluates the pulse-shape integral F^(orderSFA)_2[f^(infinite); A; angles & energies]  for all pulses for which no analytical
        expression is so easily available.an ntg::Complex{Float64} is retured.
"""
function pulseShapeQuadIntegral(envelope::Pulse.AbstractEnvelope, beam::AbstractBeam, polarization::Basics.AbstractPolarization,
                                    thetap::Float64, phip::Float64, energyp::Float64, initialEn::Float64, orderSFA::Int)
    if orderSFA == 0
    wa = 0. * im;   
    # Collect parameters that are specific to a given pulse envelope
    if       typeof(envelope) == SinSquaredEnvelope
    elseif   typeof(envelope) == GaussianEnvelope
    end
    #
    # Determine first the integrant; the timeGrid must still be adapted to thos
    timeGrid = [0.1i for i = 1:10]
    for   t  in  timeGrid 
        phase = Pulse.volkovPhase(t, envelope) ## , ...)
    end
    #
    error("Not yet implemented.")
end
return( wa )
end



###################################################################################################################
###################################################################################################################
###################################################################################################################

end # module
