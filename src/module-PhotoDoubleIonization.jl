
"""
`module  JAC.PhotoDoubleIonization`  
... a submodel of JAC that contains all methods for computing photo-double ionization properties between some initial 
    and final-state multiplets. This modules implements a simplified model for the photo-double ionization by treating
    one of the ionized electrons as semi-bound in some high-n shells and by using a (resonant) Green function multiplet
    for second-order computations. From a physics viewpoint, this model becomes better as higher n is chosen and as 
    larger the number of sharings is. The use of a semi-bound electron ensures that most parts of the photo-double
    ionization is similar to the (standard) photoionization of atoms.
"""
module PhotoDoubleIonization


using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron, ..PhotoEmission, 
                ..TableStrings

"""
`struct  PhotoDoubleIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing photo-double ionization lines.

    + multipoles              ::Array{EmMultipole}      ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge}         ... Specifies the gauges to be included into the computations.
    + quasiShells             ::Array{Shell,1}          
        ... Non-relativistic shells that emulate the quasi-free part of the two-electron continuum.
    + photonEnergies          ::Array{Float64,1}        ... List of photon energies [in user-selected units].  
    + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
    + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
    + calcDifferentialCs      ::Bool                    ... True, if the energy-differential cs are to be calculated and false otherwise.  
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + eeInteraction           ::AbstractEeInteraction   ... Type of the electron-electron interaction in the second-order treatment.
    + gMultiplet              ::Multiplet               ... Mean-field multiplet of intermediate levels in the computations.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole}
    gauges                    ::Array{UseGauge}
    quasiShells               ::Array{Shell,1} 
    photonEnergies            ::Array{Float64,1} 
    NoEnergySharings          ::Int64         
    maxKappa                  ::Int64 
    calcDifferentialCs        ::Bool 
    printBefore               ::Bool
    lineSelection             ::LineSelection
    eeInteraction             ::AbstractEeInteraction
    gMultiplet                ::Multiplet
end 


"""
`PhotoDoubleIonization.Settings()`  ... constructor for the default values of PhotoDoubleIonization line computations
"""
function Settings()
    Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Shell[], Float64[], 0, 0, false, false, 
                LineSelection(), CoulombInteraction(), Multiplet())
end


"""
`PhotoDoubleIonization.Settings(set::PhotoDoubleIonization.Settings;`

        multipoles=..,          gauges=..,                  quasiShells=..,          
        photonEnergies=..,      NoEnergySharings=..,     
        maxKappa=..,            calcDifferentialCs..,       printBefore=..,             lineSelection=..,       
        eeInteraction=..,       gMultiplet=..)
                    
    ... constructor for modifying the given PhotoDoubleIonization.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoDoubleIonization.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    quasiShells::Union{Nothing,Array{Shell,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                NoEnergySharings::Union{Nothing,Int64}=nothing,       
    maxKappa::Union{Nothing,Int64}=nothing,                                 calcDifferentialCs::Union{Nothing,Bool}=nothing,      
    printBefore::Union{Nothing,Bool}=nothing,                               lineSelection::Union{Nothing,LineSelection}=nothing, 
    eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,            gMultiplet::Union{Nothing,Multiplet}=nothing)  
    
    if  multipoles         == nothing   multipolesx         = set.multipoles         else  multipolesx         = multipoles          end 
    if  gauges             == nothing   gaugesx             = set.gauges             else  gaugesx             = gauges              end 
    if  quasiShells        == nothing   quasiShellsx        = set.quasiShells        else  quasiShellsx        = quasiShells         end 
    if  photonEnergies     == nothing   photonEnergiesx     = set.photonEnergies     else  photonEnergiesx     = photonEnergies      end 
    if  NoEnergySharings   == nothing   NoEnergySharingsx   = set.electronEnergies   else  NoEnergySharingsx   = NoEnergySharings    end 
    if  maxKappa           == nothing   maxKappax           = set.maxKappa           else  maxKappasx          = maxKappa            end 
    if  calcDifferentialCs == nothing   calcDifferentialCsx = set.calcDifferentialCs else  calcDifferentialCsx = calcDifferentialCs  end 
    if  printBefore        == nothing   printBeforex        = set.printBefore        else  printBeforex        = printBefore         end 
    if  lineSelection      == nothing   lineSelectionx      = set.lineSelection      else  lineSelectionx      = lineSelection       end 
    if  eeInteraction      == nothing   eeInteractionx      = set.eeInteraction      else  eeInteractionx      = eeInteraction       end 
    if  gMultiplet         == nothing   gMultipletx         = set.gMultiplet         else  gMultipletx         = gMultiplet          end 

    Settings( multipolesx, gaugesx, quasiShellsx, photonEnergiesx, NoEnergySharingsx, maxKappax, calcDifferentialCsx, 
                printBeforex, lineSelectionx, eeInteractionx, gMultipletx)
end


# `Base.show(io::IO, settings::PhotoDoubleIonization.Settings)`  
#   ... prepares a proper printout of the variable settings::PhotoDoubleIonization.Settings.
function Base.show(io::IO, settings::PhotoDoubleIonization.Settings) 
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "quasiShells:              $(settings.quasiShells)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "NoEnergySharings:         $(settings.NoEnergySharings)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "calcDifferentialCs:       $(settings.calcDifferentialCs)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "eeInteraction:            $(settings.eeInteraction)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
end


"""
`struct  PhotoDoubleIonization.Channel`  
    ... defines a type for a PhotoDoubleIonization channel to help characterize a single multipole and scattering (continuum) state 
        with one semi-bound and one free electron.

    + multipole      ::EmMultipole          ... Multipole of the photon absorption.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + quasiSubshell  ::Subshell             ... (Bound) subshell that represents a quasi-free electron.
    + xSymmetry      ::LevelSymmetry        ... angular momentum and parity due to the coupling of the quasi-bound subshell.
    + kappa          ::Int64                ... partial-wave of the free electron
    + tSymmetry      ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + phase          ::Float64              ... phase of the partial wave
    + amplitude      ::Complex{Float64}     ... PhotoDoubleIonization amplitude associated with the given channel.
"""
struct  Channel
    multipole        ::EmMultipole
    gauge            ::EmGauge
    quasiSubshell    ::Subshell  
    xSymmetry        ::LevelSymmetry
    kappa            ::Int64
    tSymmetry        ::LevelSymmetry
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


# `Base.show(io::IO, channel::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of sharing::PhotoDoubleIonization.Channel.
function Base.show(io::IO, channel::PhotoDoubleIonization.Channel) 
    println(io, "multipole:              $(channel.multipole)  ")
    println(io, "gauge:                  $(channel.gauge)  ")
    println(io, "quasiSubshell:          $(channel.quasiSubshell)  ")
    println(io, "xSymmetry:              $(channel.xSymmetry)  ")
    println(io, "kappa:                  $(channel.kappa)  ")
    println(io, "tSymmetry:              $(channel.tSymmetry)  ")
    println(io, "phase:                  $(channel.phase)  ")
    println(io, "amplitude:              $(channel.amplitude)  ")
end


"""
`struct  PhotoDoubleIonization.Sharing`  
    ... defines a type for a PhotoDoubleIonization sharing to help characterize energy sharing between the two emitted electrons.

    + omega          ::Float64         ... Energy of the incident photon
    + epsilon1       ::Float64         ... Energy of (free) electron 1.
    + epsilon2       ::Float64         ... Energy of (free) electron 2.
    + weight         ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
    + differentialCs ::EmProperty      ... differential cross section of this energy sharing.
    + channels       ::Array{PhotoDoubleIonization.Channel,1}  ... List of PhotoDoubleIonization channels of this line.
"""
struct  Sharing
    omega            ::Float64
    epsilon1         ::Float64
    epsilon2         ::Float64
    weight           ::Float64
    differentialCs   ::EmProperty
    channels         ::Array{PhotoDoubleIonization.Channel,1}
end


# `Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing)`  ... prepares a proper printout of sharing::PhotoDoubleIonization.Sharing.
function Base.show(io::IO, sharing::PhotoDoubleIonization.Sharing) 
    println(io, "omega:                  $(sharing.omega)  ")
    println(io, "epsilon1:               $(sharing.epsilon1)  ")
    println(io, "epsilon2:               $(sharing.epsilon2)  ")
    println(io, "weight:                 $(sharing.weight)  ")
    println(io, "differentialCs:         $(sharing.differentialCs)  ")
    println(io, "channels:               $(sharing.channels)  ")
end



"""
`struct  PhotoDoubleIonization.Line`  ... defines a type for a photo-double ionization line that may include the definition of channels.

    + initialLevel   ::Level                  ... initial-(state) level
    + finalLevel     ::Level                  ... final-(state) level
    + electronEnergy ::Float64                ... Energy of the (outgoing free) electron.
    + photonEnergy   ::Float64                ... Energy of the absorbed photon.
    + crossSection   ::EmProperty             ... Cross section for this photo-double ionization.
    + sharings       ::Array{PhotoDoubleIonization.Sharing,1}  ... List of PhotoDoubleIonization.Sharings of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    photonEnergy     ::Float64
    crossSection     ::EmProperty
    sharings         ::Array{PhotoDoubleIonization.Sharing,1}
end


"""
`PhotoDoubleIonization.Line()`  
    ... constructor an empty PhotoDoubleIonization line.
"""
function Line()
    Line(Level(), Level(), 0., 0., EmProperty(0., 0.), PhotoDoubleIonization.Sharing[])
end


# `Base.show(io::IO, line::PhotoDoubleIonization.Line)`  ... prepares a proper printout of the variable line::PhotoDoubleIonization.Line.
function Base.show(io::IO, line::PhotoDoubleIonization.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "electronEnergy:    $(line.electronEnergy)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "sharings:          $(line.sharings)  ")
end



"""
`PhotoDoubleIonization.amplitude(kind::String, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                                    gMultiplet::Multiplet, grid::Radial.Grid; display::Bool=false, printout::Bool=false)`  
    ... to compute the photo-double ionization amplitude  
    
                <alpha_f J_f || O^(Mp, absorption) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>  
            +   <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, absorption) || alpha_i J_i> 
            
        absorption amplitude for the interaction with the photon field of frequency omega, multipolarity Mp and gauge. 
        A value::ComplexF64 is returned. The amplitude value is printed to screen if display=true.
"""
function amplitude(kind::String, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                    gMultiplet::Multiplet, grid::Radial.Grid; display::Bool=false, printout::Bool=false)
    #
    # Always ensure the same subshell list for all initial, intermediate and final levels
    subshells  = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
    @show  initialLevel.basis.subshells, finalLevel.basis.subshells, subshells
    iLevel     = Level(initialLevel, subshells)
    fLevel     = Level(finalLevel, subshells)
    nMultiplet = Multiplet(gMultiplet, subshells)
    
    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity);    eni = iLevel.energy
    nn = length(nMultiplet.levels[1].basis.csfs)
    
    if  printout   printstyled("Compute photo-double $(Mp) ionizations amplitude for the transition [$(iLevel.index)-$(fLevel.index)] ... ", 
                                color=:light_green)    end
    amplitude = ComplexF64(0.)
    #
    for  r = 1:nf
        symr = LevelSymmetry(fLevel.basis.csfs[r].J, fLevel.basis.csfs[r].parity);      if  symr != symf    continue    end
        for  s = 1:ni
            syms = LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity);  if  syms != symi    continue    end
            for  nLevel in nMultiplet.levels
                symn = LevelSymmetry(nLevel.J, nLevel.parity);    enn = nLevel.energy;   @show symn
                #
                #   Compute <alpha_f J_f || O^(Mp, kind) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
                if  symn != symi     continue    end
                for  t = 1:nn
                    if  nLevel.mc[t] == 0.  continue    end
                    Vee       = ManyElectron.matrixElement_Vee(CoulombInteraction(), nLevel.basis, t, iLevel.basis, s, grid)
                    OMp       = ManyElectron.matrixElement_Mab(Mp, gauge, omega, fLevel.basis, r, nLevel.basis, t, grid)
                    amplitude = amplitude + fLevel.mc[r] * OMp * nLevel.mc[t]^2 * Vee * iLevel.mc[s] / (eni - enn)
                    @show  t, eni, (eni - enn), amplitude
                end
                #
                #   Compute <alpha_f J_f || V^(e-e) || alpha_n J_f> <alpha_n J_f || O^(Mp, kind) || alpha_i J_i>
                if  symn != symf     continue    end
                for  t = 1:nn
                    if  nLevel.mc[t] == 0.  continue    end
                    OMp       = ManyElectron.matrixElement_Mab(Mp, gauge, omega, nLevel.basis, t, iLevel.basis, s, grid)
                    Vee       = ManyElectron.matrixElement_Vee(CoulombInteraction(), fLevel.basis, r, nLevel.basis, t, grid)
                    amplitude = amplitude + fLevel.mc[r] * Vee * nLevel.mc[t]^2 * OMp * iLevel.mc[s] / (eni - omega - enn)
                    @show  t, eni, (eni + omega - enn), amplitude
                end
            end
        end
    end
    if  printout   printstyled("done. \n", color=:light_green)    end
    
    if  display  
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " PhotoDouble^($Mp, $kind) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end
    
    return( amplitude )
end



"""
`PhotoDoubleIonization.computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, quasiOrbitals::Dict{Subshell, Orbital},
                        grid::Radial.Grid, nrContinuum::Int64, settings::PhotoDoubleIonization.Settings; printout::Bool=false)`  
    ... to compute all amplitudes and properties of the given line; a line::PhotoDoubleIonization.Line is returned for which the amplitudes and 
        properties are now evaluated.
"""
function  computeAmplitudesProperties(line::PhotoDoubleIonization.Line, nm::Nuclear.Model, quasiOrbitals::Dict{Subshell, Orbital}, 
                                        grid::Radial.Grid, nrContinuum::Int64, settings::PhotoDoubleIonization.Settings; printout::Bool=false)
    newSharings = PhotoDoubleIonization.Sharing[];   contSettings = Continuum.Settings(false, nrContinuum);    csC = 0.;    csB = 0.
    for  sharing in line.sharings
        newChannels = PhotoDoubleIonization.Channel[];    dcsC = 0.;    dcsB = 0.   
        for channel in sharing.channels
            subshells = line.initialLevel.basis.subshells;  
            push!(subshells, channel.quasiSubshell);    push!(subshells, Subshell(101, channel.kappa))
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshells)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshells[1:end-2])
            newfLevel = Basics.generateLevelWithExtraElectron(quasiOrbitals[channel.quasiSubshell], channel.xSymmetry, newfLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(sharing.epsilon2, Subshell(101, channel.kappa), 
                                                                    line.finalLevel, nm, grid, contSettings)
            #
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.tSymmetry, newfLevel)
            newChannel = PhotoDoubleIonization.Channel(channel.multipole, channel.gauge, channel.quasiSubshell, channel.xSymmetry,
                                                        channel.kappa, channel.tSymmetry, phase, 0.)
            amplitude  = PhotoDoubleIonization.amplitude("absorption", channel.multipole, channel.gauge, sharing.omega, newcLevel,
                                                            newiLevel, settings.gMultiplet, grid, display=true, printout=true)           
            push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.quasiSubshell, newChannel.xSymmetry, 
                                                        newChannel.kappa, newChannel.tSymmetry, newChannel.phase, amplitude) )
            if       channel.gauge == Basics.Coulomb     dcsC = dcsC + abs(amplitude)^2
            elseif   channel.gauge == Basics.Babushkin   dcsB = dcsB + abs(amplitude)^2
            elseif   channel.gauge == Basics.Magnetic    dcsB = dcsB + abs(amplitude)^2;   dcsC = dcsC + abs(amplitude)^2
            end
        end
        dcs = Basics.EmProperty(dcsC, dcsB)
        push!(newSharings, PhotoDoubleIonization.Sharing(sharing.omega, sharing.epsilon1, sharing.epsilon2, sharing.weight, dcs, newChannels) )
    end
    Ji2 = Basics.twice(line.initialLevel.J)
    csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
    csFactor     = csFactor / 2.   # Not fully clear, arises likely from the Rydberg normalization
    ##  Correct for energy normalization 
    ##  if  line.electronEnergy < 2.0   csFactor = csFactor * (line.electronEnergy/2.0)^1.5     end
    crossSection = EmProperty(csFactor * csC, csFactor * csB)
    newLine = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, 
                                    crossSection, newSharings)
    return( newLine )
end



"""
`PhotoDoubleIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::PhotoDoubleIonization.Settings; output::Bool=true)`  
    ... to compute the photo-double ionization transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{PhotoDoubleIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::PhotoDoubleIonization.Settings; output::Bool=true)
    println("")
    printstyled("PhotoDoubleIonization.computeLines(): The computation of photo-double ionization properties starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    # Generate orbitals for all quasi-subshells
    quasiSubshells = Basics.generateSubshellList(settings.quasiShells)
    meanPot        = Basics.computePotentialDFS(grid, finalMultiplet.levels[1].basis)
    quasiOrbitals  = Basics.generateOrbitalsForPotential(grid, meanPot, quasiSubshells)
    #
    lines = PhotoDoubleIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoDoubleIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoDoubleIonization.Line[]
    for  line in lines
        println("\n>> Calculate photo-double ionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
        newLine = PhotoDoubleIonization.computeAmplitudesProperties(line, nm, quasiOrbitals, grid, nrContinuum, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoDoubleIonization.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoDoubleIonization.displayResults(iostream, newLines, settings)     end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoDoubleIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)`  
    ... to determine a list of PhotoDoubleIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{PhotoDoubleIonization.Line,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoDoubleIonization.Settings)
    lines = PhotoDoubleIonization.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                # Add lines for all photon energies
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    energy   = omega_au - (fLevel.energy - iLevel.energy)
                    @show omega, energy
                    if  energy < 0.    continue   end  
                    sharings = PhotoDoubleIonization.determineSharingsAndChannels(fLevel, iLevel, omega, energy, settings) 
                    @show sharings
                    push!( lines, PhotoDoubleIonization.Line(iLevel, fLevel, energy, omega_au, EmProperty(0., 0.), sharings) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoDoubleIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                                                    settings::PhotoDoubleIonization.Settings)`  
    ... to determine a list of PhotoDoubleIonization Sharing's and Channel's for a transitions from the initial to final level 
        and by taking into account the particular settings of for this computation; an Array{PhotoDoubleIonization.Sharing,1} 
        is returned.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, omega::Float64, energy::Float64,
                                        settings::PhotoDoubleIonization.Settings)
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    subshellList = Basics.generateSubshellList(settings.quasiShells)
    if  Basics.UseCoulomb  in  settings.gauges   gaugeM = Basics.UseCoulomb    else   gaugeM = Basics.UseBabushkin    end
        
    sharings  = PhotoDoubleIonization.Sharing[]
    eSharings = Basics.determineEnergySharings(omega - energy, settings.NoEnergySharings) 
    for  es in eSharings
        epsilon1  = es[1];    epsilon2 = es[2];    weight = es[3];    channels = PhotoDoubleIonization.Channel[]
        for  mp in settings.multipoles
            for  subsh in subshellList
                xSymmetries = AngularMomentum.allowedTotalSymmetries(symf, subsh.kappa)
                for  symx in xSymmetries
                    tSymmetries = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
                    for  symt in tSymmetries
                        kappaList = AngularMomentum.allowedKappaSymmetries(symt, symx)
                        for  kappa in kappaList
                            for  gauge in settings.gauges
                                # Include further restrictions if appropriate
                                if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                                    push!(channels, PhotoDoubleIonization.Channel(mp, Basics.Coulomb,   subsh, symx, kappa, symt, 0., Complex(0.)) )
                                elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                                    push!(channels, PhotoDoubleIonization.Channel(mp, Basics.Babushkin, subsh, symx, kappa, symt, 0., Complex(0.)) )  
                                elseif string(mp)[1] == 'M'  &&   gauge == gaugeM                               
                                    push!(channels, PhotoDoubleIonization.Channel(mp, Basics.Magnetic,  subsh, symx, kappa, symt, 0., Complex(0.)) ) 
                                end
                            end
                        end 
                    end
                end
            end
        end
        push!(sharings, PhotoDoubleIonization.Sharing(omega, epsilon1, epsilon2, weight, EmProperty(0., 0.), channels) )
    end
    return( sharings )  
end


"""
`PhotoDoubleIonization.displayLines(iostream::IO, lines::Array{PhotoDoubleIonization.Line,1})`  
    ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(iostream::IO, lines::Array{PhotoDoubleIonization.Line,1})
    #
    # First, print lines and sharings
    nx = 94
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & sharings:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
    sb = sb * TableStrings.flushleft(54, "  i -- f         omega        epsilon_1    epsilon_2"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        #
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega))    * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "   "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Second, print lines and channles
    nx = 130
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & channels:   ... channels are shown only for the first sharing")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(100, "Channels (all energies in " * TableStrings.inUnits("energy") * ")" ; na=5);              
    sb = sb * TableStrings.flushleft(100, "Multipole  Gauge      quasi-Subshell   J^P_x   kappa   -->    J^P_t"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        for (ic, ch) in enumerate(line.sharings[1].channels)
            if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * string(ch.multipole) * "         " * string(ch.gauge)[1:3] * "            " 
            sb = sb * string(ch.quasiSubshell) * "       " * string(ch.xSymmetry) * "   "
            sb = sb * string(Subshell(1,ch.kappa))[end-4:end] * "   -->    "
            sb = sb * string(ch.tSymmetry)
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
