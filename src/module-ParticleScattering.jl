
"""
`module  JAC.ParticleScattering`  
... a submodel of JAC that contains all methods for computing scattering amplitudes and cross sections
    with regard to "transitions" between some initial and final-state multiplets. It covers different types
    of particles (electrons, protons, ...), different types of scattering processes (elastic, ...), 
    different types of incoming beams (plane-wave, Bessel, Laguerre-Gauss, ...) and different treatments
    (nonrelativistic, relativistic, ...).
"""
module ParticleScattering


using  Printf, GSL, SpecialFunctions,
        ..AngularMomentum, ..Basics, ..Beam, ..Bsplines, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron,
        ..Nuclear, ..Radial, ..RadialIntegrals, ..SpinAngular, ..TableStrings


"""
`abstract type ParticleScattering.AbstractProcessType` 
    ... defines an abstract type to distinguish different scattering processes of particles with atoms and ions; see also:
    
    + struct ParticleScattering.ElasticElectronNR    ... to model the elastic electron scattering.
    + struct ParticleScattering.InelasticElectronNR  ... to model the inelastic electron scattering (not yet).
"""
abstract type  AbstractProcessType                                            end
struct   InelasticElectronNR    <:  ParticleScattering.AbstractProcessType    end
    

"""
`struct  ParticleScattering.ElasticElectronNR      <:  ParticleScattering.AbstractProcessType`  
    ... to model the elastic electron scattering and to compute different cross sections.

    + calcd2SigmaHeadon       ::Bool       ... to compute the double-differential (head-on) scattering cross sections.
    + calcd2SigmaBorn         ::Bool       ... to compute the double-differential (Born) scattering cross sections.
    + calcd2SigmaMacroscopic  ::Bool       ... to compute the double-differential (macroscopic) scattering cross sections
                                               by integrating over all impact parameters.
    + calcBdependentAmps      ::Bool       ... to compute the impact-parameter b-dependent scattering amplitudes for a 
                                               given set of b-Vectors.
"""
struct ElasticElectronNR      <:  ParticleScattering.AbstractProcessType
    calcd2SigmaHeadon         ::Bool
    calcd2SigmaBorn           ::Bool
    calcd2SigmaMacroscopic    ::Bool
    calcBdependentAmps        ::Bool
end 


"""
`ParticleScattering.ElasticElectronNR()`  ... constructor for the default ParticleScattering.ElasticElectronNR.
"""
function ElasticElectronNR()
    ElasticElectronNR(true, false, false, false)
end


# `Base.show(io::IO, proc::ParticleScattering.ElasticElectronNR)`  ... printout of the variable pw::ParticleScattering.ElasticElectronNR.
function Base.show(io::IO, proc::ParticleScattering.ElasticElectronNR) 
    print(io, "ElasticElectronNR[calcd2SigmaHeadon=$(proc.calcd2SigmaHeadon),  calcd2SigmaBorn=$(proc.calcd2SigmaBorn),  " * 
              "calcd2SigmaMacroscopic=$(proc.calcd2SigmaMacroscopic),  calcBdependentAmps=$(proc.calcBdependentAmps)]")
end
    

"""
`struct  ParticleScattering.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters in calculating the scattering amplitudes and cross sections
        for selected scattering event.

    + processType         ::ParticleScattering.AbstractProcessType
    + beamType            ::Beam.AbstractBeamType
    + polarization        ::Basics.AbstractPolarization
    + impactEnergies      ::Array{Float64,1}
    + polarThetas         ::Array{Float64,1}
    + polarPhis           ::Array{Float64,1}
    + bVectors            ::Vector{Vector{Float64}}
        ... Provides a list of 2d cartesian vectors for the impact parameter (vector) b.
    + printBefore         ::Bool               ... True, if all energies and events are printed before their evaluation.
    + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
    + epsPartialWave      ::Float64            
        ... accuracy criterion; all partial waves (l, kappa) with 
            abs(contribution_l) + abs(contribution_l+1) < epsPartialWave are neglected;
"""
struct Settings  <:  AbstractProcessSettings
    processType           ::ParticleScattering.AbstractProcessType
    beamType              ::Beam.AbstractBeamType
    polarization          ::Basics.AbstractPolarization
    impactEnergies        ::Array{Float64,1}
    polarThetas           ::Array{Float64,1}
    polarPhis             ::Array{Float64,1}
    bVectors              ::Vector{Vector{Float64}}
    printBefore           ::Bool 
    lineSelection         ::LineSelection 
    epsPartialWave        ::Float64
end 


"""
`ParticleScattering.Settings()`  ... constructor for the default ParticleScattering.Settings.
"""
function Settings()
    Settings(ElasticElectronNR(), Beam.PlaneWave(), LinearPolarization(), Float64[], Float64[], Float64[], Vector{Float64}[], false, LineSelection(), 2.)
end


"""
`ParticleScattering.Settings(set::ParticleScattering.Settings;`

        processType=..,         beamType=..,                polarization=..,          
        impactEnergies=..,      polarThetas=..,             polarPhis=..,           bVectors=..,
        printBefore=..,         lineSelection=..,           epsPartialWave=.. )
                    
    ... constructor for modifying the given ParticleScattering.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ParticleScattering.Settings; 
    processType::Union{Nothing,ParticleScattering.AbstractProcessType}=nothing,
    beamType::Union{Nothing,Beam.AbstractBeamType}=nothing,
    polarization::Union{Nothing,Basics.AbstractPolarization}=nothing,
    impactEnergies::Union{Nothing,Array{Float64,1}}=nothing,    polarThetas::Union{Nothing,Array{Float64,1}}=nothing,
    polarPhis::Union{Nothing,Array{Float64,1}}=nothing,         bVectors::Union{Nothing,Vector{Vector{Float64}}}=nothing, 
    printBefore::Union{Nothing,Bool}=nothing, 
    lineSelection::Union{Nothing,LineSelection}=nothing,        epsPartialWave::Union{Nothing,Float64}=nothing)  
    
    if  isnothing(processType)      processTypex    = set.processType       else  processTypex    = processType       end
    if  isnothing(beamType)         beamTypex       = set.beamType          else  beamTypex       = beamType          end
    if  isnothing(polarization)     polarizationx   = set.polarization      else  polarizationx   = polarization      end
    if  isnothing(impactEnergies)   impactEnergiesx = set.impactEnergies    else  impactEnergiesx = impactEnergies    end
    if  isnothing(polarThetas)      polarThetasx    = set.polarThetas       else  polarThetasx    = polarThetas       end
    if  isnothing(polarPhis)        polarPhisx      = set.polarPhis         else  polarPhisx      = polarPhis         end
    if  isnothing(bVectors)         bVectorsx       = set.bVectors          else  bVectorsx       = bVectors          end
    if  isnothing(printBefore)      printBeforex    = set.printBefore       else  printBeforex    = printBefore       end 
    if  isnothing(lineSelection)    lineSelectionx  = set.lineSelection     else  lineSelectionx  = lineSelection     end 
    if  isnothing(epsPartialWave)   epsPartialWavex = set.epsPartialWave    else  epsPartialWavex = epsPartialWave    end

    Settings( processTypex, beamTypex, polarizationx, impactEnergiesx, polarThetasx, polarPhisx, bVectorsx, 
              printBeforex, lineSelectionx, epsPartialWavex )
end


# `Base.show(io::IO, settings::ParticleScattering.Settings)`  
#       ... prepares a proper printout of the variable settings::ParticleScattering.Settings.
function Base.show(io::IO, settings::ParticleScattering.Settings) 
    println(io, "processType:           $(settings.processType)  ")
    println(io, "beamType:              $(settings.beamType)  ")
    println(io, "polarization:          $(settings.polarization)  ")
    println(io, "impactEnergies :       $(settings.impactEnergies )  ")
    println(io, "polarThetas:           $(settings.polarThetas)  ")
    println(io, "polarPhis:             $(settings.polarPhis)  ")
    println(io, "bVectors:              $(settings.bVectors)  ")
    println(io, "printBefore:           $(settings.printBefore)  ")
    println(io, "lineSelection:         $(settings.lineSelection)  ")
    println(io, "epsPartialWave:        $(settings.epsPartialWave)  ")
end
    

"""
`struct  ParticleScattering.PartialWaveNR`  
    ... defines a type to represent a single partial-wave amplitude.

    + l                   ::Int64      ... OAM of the partial wave.
    + phase               ::Float64    ... phase shift of partial-wave.
    + amplitude           ::ComplexF64 ... partial-wave amplitude
"""
struct PartialWaveNR
    l                     ::Int64   
    phase                 ::Float64  
    amplitude             ::ComplexF64 
end 


"""
`ParticleScattering.PartialWaveNR()`  ... constructor for the default ParticleScattering.PartialWaveNR.
"""
function PartialWaveNR()
    PartialWaveNR(0, 0., 0.)
end


# `Base.show(io::IO, pw::ParticleScattering.PartialWaveNR)`  ... printout of the variable pw::ParticleScattering.PartialWaveNR.
function Base.show(io::IO, pw::ParticleScattering.PartialWaveNR) 
    println(io, "l:               $(pw.l)  ")
    println(io, "phase:           $(pw.phase)  ")
    println(io, "amplitude:       $(pw.amplitude)  ")
end
    

"""
`struct  ParticleScattering.EventNR`  
    ... defines a type to collect data & results for a (no-relativistic) scattering "event".

    + processType    ::ParticleScattering.AbstractProcessType
    + beamType       ::Beam.AbstractBeamType
    + initialLevel   ::Level           ... initial-(state) level
    + finalLevel     ::Level           ... final-(state) level
    + impactEnergy   ::Float64         ... Energy of the (incoming) particle.
    + theta          ::Float64         ... polar theta of partial-wave amplitude.
    + phi            ::Float64         ... polar phi of partial-wave amplitude.
    + d2SigmaHeadon  ::Float64         ... Double-differential CS: d^2 sigma / d Omega for head-on coliisions
    + d2SigmaBorn    ::Float64         ... Double-differential CS: d^2 sigma / d Omega for head-on Born coliisions
    + d2SigmaMacros  ::Float64         ... Double-differential CS: d^2 sigma / d Omega for a macroscopic target,
                                           i.e. by integrating over all impact parameters b.
    + dSigmadEHeadon ::Float64         ... Energy-differential CS: d sigma / d epsilon
    + partialWaves   ::Array{ParticleScattering.PartialWaveNR,1}
"""
struct EventNR
    processType      ::ParticleScattering.AbstractProcessType
    beamType         ::Beam.AbstractBeamType
    initialLevel     ::Level  
    finalLevel       ::Level 
    impactEnergy     ::Float64  
    theta            ::Float64
    phi              ::Float64
    d2SigmaHeadon    ::Float64 
    d2SigmaBorn      ::Float64 
    d2SigmaMacros    ::Float64
    dSigmadEHeadon   ::Float64
    partialWaves     ::Array{ParticleScattering.PartialWaveNR,1}
end 


"""
`ParticleScattering.EventNR()`  ... constructor for the default ParticleScattering.EventNR.
"""
function EventNR()
    EventNR(ElasticElectronNR(), Beam.PlaneWave(), Level(), Level(), 0., 0., 0., 0., 0., 0., 0., PartialWaveNR[])
end


# `Base.show(io::IO, event::ParticleScattering.EventNR)`  ... printout of the variable event::ParticleScattering.EventNR.
function Base.show(io::IO, event::ParticleScattering.EventNR) 
    println(io, "processType:        $(event.processType)  ")
    println(io, "beamType:           $(event.beamType)  ")
    println(io, "initialLevel:       $(event.initialLevel)  ")
    println(io, "finalLevel:         $(event.finalLevel)  ")
    println(io, "impactEnergy:       $(event.impactEnergy)  ")
    println(io, "theta:              $(event.theta)  ")
    println(io, "phi:                $(event.phi)  ")
    println(io, "d2SigmaHeadon:      $(event.d2SigmaHeadon)  ")
    println(io, "d2SigmaBorn:        $(event.d2SigmaBorn)  ")
    println(io, "d2SigmaMacros:      $(event.d2SigmaMacros)  ")
    println(io, "dSigmadEHeadon:     $(event.dSigmadEHeadon)  ")
    println(io, "partialWaves:       $(event.partialWaves)  ")
end


"""
`ParticleScattering.amplitude(processType::ElasticElectronNR, beamType::Beam.PlaneWave, 
                                l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)`
    ... to compute the partial-wave (l) amplitude for the nonrelativistic elastic electron scattering of plane-wave beam 
        at given theta and phi. An  amplitude::ComplexF64 is returned.
"""
function amplitude(processType::ElasticElectronNR, beamType::Beam.PlaneWave, 
                   l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)
    amplitude = ComplexF64(0.)
    if  beamType.kx == beamType.ky == 0.
        amplitude = amplitude + (2*l + 1) / (4*pi) * exp( im*lPhase ) * sin(lPhase) * GSL.sf_legendre_Pl_e(l, cos(theta)).val
        amplitude = 4pi * amplitude / beamType.kz
    else  error("stop a")
    end
    
    return( amplitude )
end


"""
`ParticleScattering.amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                                l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)`
    ... to compute the partial-wave (l) amplitude for the nonrelativistic elastic electron scattering of Bessel beam 
        at given theta and phi. An  amplitude::ComplexF64 is returned.
"""
function amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                   l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)
    amplitude = ComplexF64(0.)
    # 
    if  l >= abs(beamType.mOAM)
        # Only partial amplitudes with l >= abs(mOAM) are nonzero
        k         = beamType.kz / cos(beamType.openingAngle)
        wa        = GSL.sf_legendre_Plm(l, abs(beamType.mOAM), cos(beamType.openingAngle));   wc = 1.0
        if   beamType.mOAM < 0
            m  = abs(beamType.mOAM)
            wc = (-1)^m * factorial(big(l-m)) / factorial(big(l+m))     
            wc = Float64(wc)
        end
        wb        = 4pi * (-1.0im)^beamType.mOAM / k * exp( im*lPhase ) * (-1.)^beamType.mOAM                               *
                    sqrt( (2*l + 1) * factorial( big(l - beamType.mOAM) ) / (4*pi * factorial( big(l + beamType.mOAM) ) ) ) * 
                    wa * wc * AngularMomentum.sphericalYlm(l, beamType.mOAM, theta, phi) * sin(lPhase)
        amplitude = amplitude + ComplexF64(wb)
    end
    
    return( amplitude )
end


"""
`ParticleScattering.amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, nu::Int64,  
                              l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)`
    ... to compute the partial-wave (l) amplitude for the nonrelativistic elastic electron scattering of Bessel beam 
        at given theta and phi as well as at given impact parameter (b, phib). An  amplitude::ComplexF64 is returned.
"""
function amplitude(processType::ElasticElectronNR, beamType::Beam.BesselBeam, nu::Int64,  
                   l::Int64, lPhase::Float64, theta::Float64, phi::Float64, grid::Radial.Grid)
    
    amplitude = ComplexF64(0.)
    #
    if  l >= abs(beamType.mOAM + nu)
        # Only partial amplitudes with l >= abs(mOAM) are nonzero
        wa        = GSL.sf_legendre_Plm(l, abs(beamType.mOAM + nu), cos(beamType.openingAngle));   wc = 1.0
        if   beamType.mOAM + nu < 0
            m  = abs(beamType.mOAM + nu)
            wc = (-1)^m * factorial(big(l-m)) / factorial(big(l+m))     
            wc = Float64(wc)
        end
        wb        = exp( im*lPhase ) * (-1.)^(beamType.mOAM+nu)           *
                    sqrt( (2*l + 1) * factorial( big(l - beamType.mOAM - nu) ) / (4*pi * factorial( big(l + beamType.mOAM - nu) ) ) ) * 
                    wa * wc * AngularMomentum.sphericalYlm(l, beamType.mOAM+nu, theta, phi) * sin(lPhase)
        amplitude = amplitude + ComplexF64(wb)
    end
    
    return( amplitude )
end


"""
`ParticleScattering.computeAmplitudesProperties(processType::ElasticElectronNR, event::ParticleScattering.EventNR, 
                                                nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                settings::ParticleScattering.Settings; printout::Bool=true)` 
    ... to compute all amplitudes and properties of the given event; an event::ParticleScattering.EventNR is returned 
        for which the amplitudes and properties are now evaluated. For elastic scattering, the initial and final levels are 
        expected to be the same.
"""
function computeAmplitudesProperties(processType::ElasticElectronNR, event::ParticleScattering.EventNR, nm::Nuclear.Model, 
                                     grid::Radial.Grid, nrContinuum::Int64, settings::ParticleScattering.Settings; printout::Bool=true) 
    newPws   = ParticleScattering.PartialWaveNR[];   contSettings = Continuum.Settings(false, nrContinuum)
    d2SigmaHeadon = 0.;   d2SigmaBorn = 0.;   d2SigmaMacros = 0.
    # The nuclear potential and the B-spline primitives depend only on (nm, grid), which are constant here, while the
    # partial-wave loops below call for one continuum orbital each -- 31 for the head-on cross sections and, where the
    # b-dependent amplitudes are switched on, 11 x 1001 of them. Building both once and passing them keeps that work
    # out of the loops; omitting the keywords would reproduce the same numbers, only more slowly.
    nucPot     = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    
    if  false  ## processType.calcd2SigmaHeadon
        # Compute the head-on scattering cross sections
        amp1 = amp2 = amp3 = amp4 = amp5 = amp6 = amp7 = amp8 = 1.0e6;   maxamp = 0.
        totalAmp = ComplexF64(0.)
        for  l = 0:1000
            kappa      = -l - 1;    
            cOrbital, lPhase  = Continuum.generateOrbitalForLevel(event.impactEnergy, Subshell(101, kappa), event.finalLevel,
                                                                  nm, grid, contSettings; nuclearPot=nucPot,
                                                                  primitives=primitives)
            amplitude  = ParticleScattering.amplitude(event.processType, event.beamType, l, lPhase, event.theta, event.phi, grid)
            totalAmp   = totalAmp + amplitude
            push!( newPws, ParticleScattering.PartialWaveNR(l, lPhase, amplitude) )
            amp1 = amp2;   amp2 = amp3;   amp3 = amp4;   amp4 = amp5;   amp5 = amp6;   amp6 = amp7;   amp7 = amp8;   
            amp8 = abs(amplitude)^2;   maxamp = max(maxamp, amp8)
            if  l > 10  &&  (amp1 + amp2 + amp3 + amp4 + amp5 + amp6 + amp7 + amp8) / maxamp < settings.epsPartialWave    break     end
        end
        d2SigmaHeadon  = conj(totalAmp) * totalAmp
    end
    
    if  processType.calcd2SigmaHeadon
        # Compute the head-on scattering cross sections
        totalAmp = ComplexF64(0.)
        for  l = 0:30
            kappa      = -l - 1;    
            cOrbital, lPhase  = Continuum.generateOrbitalForLevel(event.impactEnergy, Subshell(101, kappa), event.finalLevel,
                                                                  nm, grid, contSettings; nuclearPot=nucPot,
                                                                  primitives=primitives)
            amplitude  = ParticleScattering.amplitude(event.processType, event.beamType, l, lPhase, event.theta, event.phi, grid)
            totalAmp   = totalAmp + amplitude
            push!( newPws, ParticleScattering.PartialWaveNR(l, lPhase, amplitude) )
        end
        d2SigmaHeadon  = conj(totalAmp) * totalAmp
    end
    
    if  processType.calcd2SigmaBorn
        # Compute the macroscopic scattering cross sections, integrated over all impact parameters
        println("\n\n>> Compute (Born) double-differential cross sections. ... \n\n")
        amp1 = amp2 = amp3 = amp4 = amp5 = amp6 = amp7 = amp8 = 1.0e6;   maxamp = 0.
        totalAmp = ComplexF64(0.);   k = event.beamType.kz / cos(event.beamType.openingAngle)
        for  l = 0:1000
            #== Compute the lPhaseBorn by means of JAC orbitals
            kappa      = -l - 1;
            potZero = Radial.Potential("zero potential", zeros(grid.NoPoints), grid)
            potDFS  = Basics.computePotential(Basics.DFSField(0.7), grid, event.finalLevel) 
            cOrbital, lPhase  = Continuum.generateOrbitalLocalPotential(event.impactEnergy, Subshell(101, kappa), potZero, contSettings)
            wx = Float64[];   mtp = length(cOrbital.P)   
            for  m = 1:mtp 
                wy = (cOrbital.P[m]^2 + cOrbital.Q[m]^2) * potDFS.Zr[m] / grid.r[m];    push!(wx, -wy )
            end  ==#
            # Compute the lPhaseBorn by means of spherical bessel functions
            wfl2 = Float64[];  mtp = grid.NoPoints-20;      wx = 0.;  kr = 0.
            nuclearPotential  = Nuclear.nuclearPotential(nm, grid)
            potDFS            = Basics.computePotential(Basics.DFSField(0.7), grid, event.finalLevel) 
            pot               = Basics.add(nuclearPotential, potDFS)
            for  m = 1:mtp   kr = k * grid.r[m];    wx = kr * SpecialFunctions.sphericalbesselj(l, kr); 
                push!(wfl2, - wx^2 * pot.Zr[m] / grid.r[m])   
            end
            wint = - RadialIntegrals.V0(wfl2, mtp, grid::Radial.Grid) / k
            #
            #
            lPhaseBorn = atan(wint)
            amplitude  = ParticleScattering.amplitude(event.processType, event.beamType, l, lPhaseBorn, event.theta, event.phi, grid)
            totalAmp   = totalAmp + amplitude
            push!( newPws, ParticleScattering.PartialWaveNR(l, lPhaseBorn, amplitude) )
            amp1 = amp2;   amp2 = amp3;   amp3 = amp4;   amp4 = amp5;   amp5 = amp6;   amp6 = amp7;   amp7 = amp8;   
            amp8 = abs(amplitude)^2;   maxamp = max(maxamp, amp8)
            if  l > 10  &&  (amp1 + amp2 + amp3 + amp4 + amp5 + amp6 + amp7 + amp8) / maxamp < settings.epsPartialWave    break     end
        end
        d2SigmaBorn = conj(totalAmp) * totalAmp    
    end
    
    if  processType.calcd2SigmaMacroscopic
        # Compute the macroscopic scattering cross sections, integrated over all impact parameters
        d2SigmaMacros = ParticleScattering.crossSectionMacroscopic(event.processType, event.beamType, newPws, event.theta, event.phi, grid)        
    end
    
    if  processType.calcBdependentAmps
        # Compute b-vector dependent scattering amplitudes; they are generated an printed in a neat table but not (yet)
        # brought to the final outcome of the event. ... Compute the head-on scattering amplitudes
        nb = length(settings.bVectors)
        wb = NamedTuple{(:bx, :by, :b, :phib, :amp), Tuple{Float64, Float64, Float64, Float64, ComplexF64}}[]
        #
        for  bVector in settings.bVectors
            bx, by   = bVector;   b = sqrt(bx^2 + by^2);   phib = angle(bx + by*im)
            k        = event.beamType.kz / cos(event.beamType.openingAngle)
            krho     = k * sin(event.beamType.openingAngle)
            nuAmp    = ComplexF64(0.)
            for  nu = -5:5
                amp1 = amp2 = amp3 = amp4 = 1.0e6;    lAmp = ComplexF64(0.)
                for  l = 0:1000
                    kappa      = -l - 1;    
                    cOrbital, lPhase  = Continuum.generateOrbitalForLevel(event.impactEnergy, Subshell(101, kappa), event.finalLevel,
                                                                          nm, grid, contSettings; nuclearPot=nucPot,
                                                                          primitives=primitives)
                    amp  = ParticleScattering.amplitude(event.processType, event.beamType, nu, l, lPhase, event.theta, event.phi, grid)
                    lAmp = lAmp + amp
                    amp1 = amp2;   amp2 = amp3;   amp3 = amp4;   amp4 = abs(amp)^2
                    if  l > 5  &&  (amp1 + amp2 + amp3 + amp4) < settings.epsPartialWave    break     end
                end
                nuAmp = nuAmp + (-1.0im)^nu * GSL.sf_bessel_Jnu(nu, krho*b) * exp( - im*nu+phib )
            end
            nuAmp = nuAmp * 2 * (-1.0im)^event.beamType.mOAM / k 
            println("******* nuAmp = $nuAmp ")

            push!(wb, (bx=bx, by=by, b=b, phib=phib, amp=nuAmp) )
        end
        #
        sym = LevelSymmetry(event.initialLevel.J, event.initialLevel.parity)
        println("\n > Compute b-vector dependent scattering amplitudes for $nb impact parameters with ..." *
                "\n      J^P = " * string(sym)                    * "             ... initial level "      *
                "\n      " * @sprintf("%.6e", event.impactEnergy) * "    ... impact energy [Hartree]  \n" )
        #
        ParticleScattering.displayBdependentAmplitudes(stdout, wb)
    end
    
    newEvent = ParticleScattering.EventNR(event.processType, event.beamType, event.initialLevel, event.finalLevel, event.impactEnergy, 
                                          event.theta, event.phi, d2SigmaHeadon, d2SigmaBorn, d2SigmaMacros, 0., newPws)

    return( newEvent )
end


"""
`ParticleScattering.computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::ParticleScattering.Settings; output=true, printout::Bool=true)`  
    ... to compute the particle scattering amplitudes and all properties as requested by the given settings. A list of 
        events::Array{ParticleScattering.Events} is returned.
"""
function  computeEvents(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::ParticleScattering.Settings; output=true, printout::Bool=true)
    println("") 
    printstyled("ParticleScattering.computeEvents(): The computation of Auger rates and properties starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------ \n", color=:light_green)
    println("")
    #
    if !( typeof(settings.processType) in [ElasticElectronNR] )   error("stop a")    end
    #
    # Distinguish different sub-procedures ... if relativistic scattering will be considered in the future
    events = ParticleScattering.determineEventsNR(finalMultiplet, initialMultiplet, settings)
    # Display all selected events before the computations start
    if  settings.printBefore    ParticleScattering.displayEvents(settings.processType, events)    end  
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  event in events   maxEnergy = max(maxEnergy, event.impactEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newEvents = ParticleScattering.EventNR[]
    for  event in events
        newEvent = ParticleScattering.computeAmplitudesProperties(settings.processType, event, nm, grid, nrContinuum, settings) 
        push!( newEvents, newEvent)
    end
    # Print all results to screen
    ParticleScattering.displayAmplitudes(stdout, settings.processType, newEvents, settings)
    ParticleScattering.displayCrossSections(stdout, settings.processType, newEvents, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ParticleScattering.displayAmplitudes(iostream, settings.processType, newEvents, settings)
                       ParticleScattering.displayCrossSections(iostream, settings.processType, newEvents, settings)    end

    if    output    return( newEvents )
    else            return( nothing )
    end
end


"""
`ParticleScattering.crossSectionMacroscopic(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                                pws::Array{ParticleScattering.PartialWaveNR,1}, theta::Float64, phi::Float64, grid::Radial.Grid)`
    ... to compute the (macroscopic) triple-differential cross section d2 Sigma / d Omega, intergrated over all the 
        impact parameters. A d2SigmaMacros::Float64 is returned.
"""
function crossSectionMacroscopic(processType::ElasticElectronNR, beamType::Beam.BesselBeam, 
                                 pws::Array{ParticleScattering.PartialWaveNR,1}, theta::Float64, phi::Float64, grid::Radial.Grid)
    #
    function Theta(l::Int64, m::Int64, theta::Float64)
        println(">>>>> l = $l,  m = $m, theta = $theta ")
        if  abs(m) > l   wx =  ComplexF64(0.) 
        elseif  m >= 0   wx = (-1.0)^m * sqrt( (2*l+1) * factorial(big(l-m)) / (2 * factorial(big(l+m)) ) ) * 
                              GSL.sf_legendre_Plm(l, m, cos(theta))
        else             wx = (-1.0)^abs(m) * Theta(l, abs(m), theta)
        end
        return( wx )
    end
    #
    d2SigmaMacros = 0.
    #
    println(">> Compute (macroscopic) double-differential cross sections based on $(length(pws)) partial waves.")
    #
    wb = 0.;   k = beamType.kz / cos(beamType.openingAngle);   kperp = k * sin(beamType.openingAngle)
    #
    wa = ComplexF64(0.)
    for  pwa in pws,    pwb in pws
        for  m = -pwa.l:pwa.l 
            wa = wa + exp( im*(pwa.phase - pwb.phase) ) * sin(pwa.phase) * sin(pwb.phase) * 
                      Theta(pwa.l, m, theta) * Theta(pwb.l, m, theta) * 
                      Theta(pwa.l, m, beamType.openingAngle) * Theta(pwb.l, m, beamType.openingAngle) 
        end
    end
    #
    wa = Float64( wa.re )
    d2SigmaMacros = 16 * pi^2 * wa / (k * beamType.kz) 
    #
    #==  First attempt to calculate macroscopic cross sections ... which should not be m-dependent !!
    for  nu = -10:50
        wa = ComplexF64(0.)
        for  pw in pws 
            l = pw.l;   lPhase = pw.phase
            if  l >= abs(beamType.mOAM + nu)
                # Only partial amplitudes with l >= abs(mOAM+nu) are nonzero
                wx        = GSL.sf_legendre_Plm(l, abs(beamType.mOAM+nu), cos(beamType.openingAngle));   wc = 1.0
                if   beamType.mOAM + nu < 0
                    m  = abs(beamType.mOAM + nu)
                    wc = (-1)^m * factorial(big(l-m)) / factorial(big(l+m))     
                    wc = Float64(wc)
                end
                wx = exp( im*lPhase ) * (-1.)^(beamType.mOAM+nu)                               *
                     sqrt( (2*l + 1) * factorial( big(l - beamType.mOAM - nu) ) / (4*pi * factorial( big(l + beamType.mOAM - nu) ) ) ) * 
                     wx * wc * AngularMomentum.sphericalYlm(l, beamType.mOAM + nu, theta, phi) * sin(lPhase)
                wa = wa + ComplexF64(wx)
            end
        end
        wb = wb + 4 / k^2 * abs(wa)^2
        @show nu, wb, 4 / k^2 * abs(wa)^2
    end
    #
    d2SigmaMacros = wb / kperp ==#
    
    println(">>>> d2SigmaMacros = $d2SigmaMacros ")
    
    return( d2SigmaMacros )
end


"""
`ParticleScattering.determineEventsNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)`  
    ... to determine a list of ParticleScattering.Event's for transitions between levels from the initial- and final-state multiplets, and  
        by taking into account the particular selections and settings for this computation; an Array{ParticleScattering.Event,1} is returned. 
        Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineEventsNR(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ParticleScattering.Settings)
    events = ParticleScattering.EventNR[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                for en in settings.impactEnergies
                    enau    = Defaults.convertUnits("energy: to atomic", en)
                    newBeam = Beam.redefineEnergy(enau, settings.beamType)
                    for  theta in settings.polarThetas,   phi in settings.polarPhis
                        push!(events, ParticleScattering.EventNR(settings.processType, newBeam, iLevel, fLevel, enau, 
                                                                 theta, phi, 0., 0., 0., 0., ParticleScattering.PartialWaveNR[] ) )
                    end
                end
            end
        end
    end
    return( events )
end
            

"""
`ParticleScattering.displayBdependentAmplitudes(stream::IO, 
                    wb::Array{NamedTuple{(:bx, :by, :b, :phib, :amp), Tuple{Float64, Float64, Float64, Float64, ComplexF64}},1})`  
    ... to display a a neat table of impact-parameter dependent amplitudes; nothing is returned otherwise.
"""
function  displayBdependentAmplitudes(stream::IO, 
                 wb::Array{NamedTuple{(:bx, :by, :b, :phib, :amp), Tuple{Float64, Float64, Float64, Float64, ComplexF64}},1})
    nx = 90
    println(stream, " ")
    println(stream, "  Selected b-dependent scattering amplitudes:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(12, "b_x [a_o]"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "b_y [a_o]"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, " b [a_o] "; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "  phi_b  "; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(30, "(complex) amplitude"; na=2);                      sb = sb * TableStrings.hBlank(20)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  w in wb
        sa  = "    " * @sprintf("%.3e", w.bx) * "     " * @sprintf("%.3e", w.by)   * "     "
        sa  = sa     * @sprintf("%.3e", w.b)  * "     " * @sprintf("%.3e", w.phib) * "     "
        sa  = sa     * @sprintf("%.8e", w.amp.re)  * "  " * @sprintf("%.8e", w.amp.im) * "  "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`ParticleScattering.displayEvents(processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1})`  
    ... to display a list of events and partial waves that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayEvents(processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1})
    nx = 106
    println(" ")
    println("  Selected elastic electron scattering events:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(18, "No partial waves"; na=2);                         sb = sb * TableStrings.hBlank(20)           
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    for  event in events
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                        fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(20, string(typeof(event.beamType)); na=2 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))         * "    "
        sa = sa * TableStrings.center(18, string( length(event.partialWaves)); na=2 )
        println( sa )
    end
    println("  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`ParticleScattering.displayAmplitudes(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                        settings::ParticleScattering.Settings)`  
    ... to display a list of events and partial waves that have been selected due to the prior settings. A neat table of all selected 
        transitions, energies, angles and partial-wave amplitudes is printed but nothing is returned otherwise.
"""
function  displayAmplitudes(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                            settings::ParticleScattering.Settings)
    nx = 158
    println(stream, " ")
    println(stream, "  Selected elastic electron scattering events:  Scattering amplitudes")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=6)
    sa = sa * TableStrings.center(10, "theta"; na=2);                                    sb = sb * TableStrings.hBlank(12)               
    sa = sa * TableStrings.center(10, "phi";   na=3);                                    sb = sb * TableStrings.hBlank(13)               
    sa = sa * TableStrings.center( 3, "l";     na=2);                                    sb = sb * TableStrings.hBlank( 3)               
    sa = sa * TableStrings.center(30, "Partial-wave amplitudes";  na=2);                  
    sb = sb * TableStrings.center(30, "Re-Amplitude-Im";          na=2)           
    sa = sa * TableStrings.center(10, "phase"; na=2);                                    sb = sb * TableStrings.hBlank(12)               
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    nsa = 0
    for  event in events
        if  length(event.partialWaves) == 0     continue    end
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                       fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(18, string(typeof(event.beamType)); na=1 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=1)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))  * "     ";     nsa = length(sa)
        sa = sa * @sprintf("%.3e", event.theta)                           * "  "
        sa = sa * @sprintf("%.3e", event.phi)                             * "     "
        sa = sa * string(event.partialWaves[1].l)                         * "   " 
        sa = sa * @sprintf("% 1.6e", event.partialWaves[1].amplitude.re)    * "   " 
        sa = sa * @sprintf("% 1.6e", event.partialWaves[1].amplitude.im)    * "   " 
        println(stream,  sa)
        for  j in 2:length(event.partialWaves)
            sa = " "^nsa
            sa = sa * @sprintf("%.3e", event.theta)                             * "  "
            sa = sa * @sprintf("%.3e", event.phi)                               * "     "
            sa = sa * string(event.partialWaves[j].l)                           * "   " 
            sa = sa * @sprintf("% 1.6e", event.partialWaves[j].amplitude.re)    * "   " 
            sa = sa * @sprintf("% 1.6e", event.partialWaves[j].amplitude.im)    * "   " 
            sa = sa * @sprintf("% 1.6e", event.partialWaves[j].phase)           * "   " 
            println(stream,  sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end


"""
`ParticleScattering.displayCrossSections(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                            settings::ParticleScattering.Settings)`  
    ... to display a list of events and associated triple-differential and energy-differential cross sections. 
        A neat table of all selected transitions, energies and cross sections is printed but nothing is returned otherwise.
"""
function  displayCrossSections(stream::IO, processType::ElasticElectronNR, events::Array{ParticleScattering.EventNR,1},
                                settings::ParticleScattering.Settings)
    nx = 208
    println(stream, " ")
    println(stream, "  Selected elastic electron scattering events:  Differential cross sections for $(settings.beamType)")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sc = TableStrings.inUnits("cross section") 
    sa = sa * TableStrings.center(20, "Beam type"; na=2);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(16, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(16, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Impact energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=6)
    sa = sa * TableStrings.center(10, "theta"; na=2);                                    sb = sb * TableStrings.hBlank(12)               
    sa = sa * TableStrings.center(10, "phi";   na=0);                                    sb = sb * TableStrings.hBlank(10)               
    sa = sa * TableStrings.center(24, "(head-on) d^2 sigma/d Om";  na=1); 
    sb = sb * TableStrings.center(24, sc * "/sr" ;  na=1)           
    sa = sa * TableStrings.center(24, "(Born) d^2 sigma/d Om";  na=1); 
    sb = sb * TableStrings.center(24, sc * "/sr" ;  na=1)           
    sa = sa * TableStrings.center(24, "(macros.) d^2 sigma/d Om";  na=1); 
    sb = sb * TableStrings.center(24, sc * "/sr" ;  na=1)           
    sa = sa * TableStrings.center(24, "(head-on) d sigma / dE";  na=1);                  
    sb = sb * TableStrings.center(24, sc * "/dE" ;  na=1)           
    sa = sa * TableStrings.center( 8, "kz";   na=2);                                     sb = sb * TableStrings.hBlank(10)               
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    nsa = 0;    d2Convert = Defaults.convertUnits("cross section: from atomic to predefined unit", 1.0) 
                deConvert = d2Convert / Defaults.convertUnits("energy: from atomic to predefined unit", 1.0)
    for  event in events
        sa  = "  ";    isym = LevelSymmetry( event.initialLevel.J, event.initialLevel.parity)
                        fsym = LevelSymmetry( event.finalLevel.J,   event.finalLevel.parity)
        sa = sa * TableStrings.center(18, string(typeof(event.beamType)); na=1 ) 
        sa = sa * TableStrings.center(18, TableStrings.levels_if(event.initialLevel.index, event.finalLevel.index); na=1)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", event.impactEnergy))  * "     ";     nsa = length(sa)
        sa = sa * @sprintf("%.3e", event.theta)                   * "  "
        sa = sa * @sprintf("%.3e", event.phi)                     * "        "
        sa = sa * @sprintf("%.6e", event.d2SigmaHeadon * d2Convert)     * "             " 
        sa = sa * @sprintf("%.6e", event.d2SigmaBorn   * d2Convert)     * "             " 
        sa = sa * @sprintf("%.6e", event.d2SigmaMacros * d2Convert)     * "             " 
        sa = sa * @sprintf("%.6e", event.dSigmadEHeadon * deConvert)    * "    " 
        sa = sa * @sprintf("%.3e", event.beamType.kz)             * "   " 
        println(stream,  sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #
    return( nothing )
end

"""
`ParticleScattering.extractCrossSections(events::Array{ParticleScattering.EventNR,1}, energy, phi)`  
    ... to extract all triple and energy-differential scattering cross sections for (constant) energy and phi which are 
        provided by the events. A name tuple (thetas::Array{Float64,1}, d2SigmaHeadons::Array{Float64,1}, 
                                              d2SigmaMacross::Array{Float64,1}, dSigmadEHeadons::Array{Float64,1})
        is returned which provide the corresponding values.
"""
function  extractCrossSections(events::Array{ParticleScattering.EventNR,1}, energy, phi)
    thetas = Float64[];    d2SigmaHeadons = Float64[];    d2SigmaBorns = Float64[];    d2SigmaMacross = Float64[];    dSigmadEHeadons = Float64[];  
    enau   = Defaults.convertUnits("energy: to atomic", energy)
    for event in events
        if   enau == event.impactEnergy   &&   phi == event.phi 
            push!(thetas, event.theta);   push!(d2SigmaHeadons, event.d2SigmaHeadon);   push!(d2SigmaBorns, event.d2SigmaBorn)   
            push!(d2SigmaMacross, event.d2SigmaMacros);   push!(dSigmadEHeadons, event.dSigmadEHeadon)
        end
    end
    
    wa = (thetas=thetas, d2SigmaHeadons=d2SigmaHeadons, d2SigmaBorns=d2SigmaBorns, d2SigmaMacross=d2SigmaMacross, 
          dSigmadEHeadons=dSigmadEHeadons)
    return(wa)
end    

end # module
