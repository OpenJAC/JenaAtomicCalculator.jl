
"""
`module  PhotoExcitationAutoion`  
... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross sections and rates.
"""
module PhotoExcitationAutoion 

using Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..PhotoEmission, 
                ..PhotoIonization, ..Radial, ..TableStrings

"""
`struct  PhotoExcitationAutoion.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing photon-impact excitation-autoionization pathways 
        |i(N)>  --> |m(N)>  --> |f(N-1)>.

    + multipoles              ::Array{Basics.EmMultipole,1} ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{Basics.UseGauge,1}    ... Specifies the gauges to be included into the computations.
    + calcPartialCs           ::Bool                        ... Calculate the partial excitation-autoionization cross sections
                                                                sigma(i-e-f), and where the direct photoionization is neglected.
    + calcAngular             ::Bool                        ... Calculate the angular distribution of the resonantly emitted Auger electrons
                                                                at selected solid angles.
    + calcFano                ::Bool                        ... Calculate the Fano parameters of individual cs resonances sigma(i-e-f).
    + printBefore             ::Bool                        ... True, if all energies and lines are printed before their evaluation.
    + incidentStokes          ::ExpStokes                   ... Stokes parameters of the incident radiation.
    + solidAngles             ::Array{SolidAngle,1}         ... List of solid angles [(theta_1, pho_1), ...].  
    + electronEnergyShift     ::Float64                     ... An overall energy shift for all electron energies.
    + maxKappa                ::Int64                       ... Maximum kappa value of partial waves to be included for free electrons.
    + pathwaySelection        ::PathwaySelection            ... Specifies the selected levels/pathways, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                ::Array{Basics.EmMultipole,1}
    gauges                    ::Array{Basics.UseGauge,1} 
    calcPartialCs             ::Bool
    calcAngular               ::Bool
    calcFano                  ::Bool
    printBefore               ::Bool  
    incidentStokes            ::ExpStokes
    solidAngles               ::Array{SolidAngle,1}
    electronEnergyShift       ::Float64 
    maxKappa                  ::Int64 
    pathwaySelection          ::PathwaySelection  
end 


"""
`PhotoExcitationAutoion.Settings()`
    ... constructor for the default values of photon-impact excitation-autoionizaton settings.
"""
function Settings()
    Settings( Basics.EmMultipole[], UseGauge[], false, false, false, false, Basics.ExpStokes(), SolidAngle[], 0., 0, PathwaySelection())
end


"""
`PhotoExcitationAutoion.Settings(set::PhotoExcitationAutoion.Settings; keywords...)`
    ... keyword constructor to modify selected fields of a PhotoExcitationAutoion.Settings object;
        for fields not given, the values of set are used.
"""
function Settings(set::PhotoExcitationAutoion.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,           gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
        calcPartialCs::Union{Nothing,Bool}=nothing,                        calcAngular::Union{Nothing,Bool}=nothing,
        calcFano::Union{Nothing,Bool}=nothing,                             printBefore::Union{Nothing,Bool}=nothing,
        incidentStokes::Union{Nothing,ExpStokes}=nothing,                  solidAngles::Union{Nothing,Array{SolidAngle,1}}=nothing,
        electronEnergyShift::Union{Nothing,Float64}=nothing,               maxKappa::Union{Nothing,Int64}=nothing,
        pathwaySelection::Union{Nothing,PathwaySelection}=nothing)
    if  isnothing(multipoles)            multipolsx           = set.multipoles          else   multipolsx           = multipoles          end
    if  isnothing(gauges)                gaugesx              = set.gauges              else   gaugesx              = gauges              end
    if  isnothing(calcPartialCs)         calcPartialCsx       = set.calcPartialCs       else   calcPartialCsx       = calcPartialCs       end
    if  isnothing(calcAngular)           calcAngularx         = set.calcAngular         else   calcAngularx         = calcAngular         end
    if  isnothing(calcFano)              calcFanox            = set.calcFano            else   calcFanox            = calcFano            end
    if  isnothing(printBefore)           printBeforex         = set.printBefore         else   printBeforex         = printBefore         end
    if  isnothing(incidentStokes)        incidentStokesx      = set.incidentStokes      else   incidentStokesx      = incidentStokes      end
    if  isnothing(solidAngles)           solidAnglesx         = set.solidAngles         else   solidAnglesx         = solidAngles         end
    if  isnothing(electronEnergyShift)   electronEnergyShiftx = set.electronEnergyShift else   electronEnergyShiftx = electronEnergyShift end
    if  isnothing(maxKappa)              maxKappax            = set.maxKappa            else   maxKappax            = maxKappa            end
    if  isnothing(pathwaySelection)      pathwaySelectionx    = set.pathwaySelection    else   pathwaySelectionx    = pathwaySelection    end

    Settings( multipolsx, gaugesx, calcPartialCsx, calcAngularx, calcFanox, printBeforex,
              incidentStokesx, solidAnglesx, electronEnergyShiftx, maxKappax, pathwaySelectionx )
end


# `Base.show(io::IO, settings::PhotoExcitationAutoion.Settings)`  
# 	... prepares a proper printout of the variable settings::PhotoExcitationAutoion.Settings.  
function Base.show(io::IO, settings::PhotoExcitationAutoion.Settings) 
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "calcPartialCs:           $(settings.calcPartialCs)  ")
    println(io, "calcAngular:             $(settings.calcAngular)  ")
    println(io, "calcFano:                $(settings.calcFano)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "incidentStokes:          $(settings.incidentStokes)  ")
    println(io, "solidAngles:             $(settings.solidAngles)  ")
    println(io, "electronEnergyShift:     $(settings.electronEnergyShift)  ")
    println(io, "maxKappa:                $(settings.maxKappa)  ")
    println(io, "pathwaySelection:        $(settings.pathwaySelection)  ")
end


"""
`struct  PhotoExcitationAutoion.Pathway`  
    ... defines a type for a photon-impact excitation pathway that may include the definition of different 
        excitation and autoionization channels and their corresponding amplitudes.

    + initialLevel        ::Level           ... initial-(state) level
    + intermediateLevel   ::Level           ... intermediate-(state) level
    + finalLevel          ::Level           ... final-(state) level
    + excitEnergy         ::Float64         ... photon excitation energy of this pathway
    + electronEnergy      ::Float64         ... energy of the (finally outgoing, scattered) electron
    + partialCs           ::EmProperty      ... partial cross section sigma(i-e-f) of this pathway
    + qFano               ::EmProperty      ... Fano-q parameter of a resonance (i-e-f)
    + excitChannels       ::Array{PhotoEmission.Channel,1}      ... List of excitation channels of this pathway.
    + augerChannels       ::Array{AutoIonization.Channel,1}     ... List of Auger channels of this pathway.
    + photoChannels       ::Array{PhotoIonization.PartialWave,1} ... Partial waves of the photoionization step.
"""
struct  Pathway
    initialLevel          ::Level
    intermediateLevel     ::Level
    finalLevel            ::Level
    excitEnergy           ::Float64
    electronEnergy        ::Float64
    partialCs             ::EmProperty
    qFano                 ::EmProperty 
    excitChannels         ::Array{PhotoEmission.Channel,1}  
    augerChannels         ::Array{AutoIonization.Channel,1}
    photoChannels         ::Array{PhotoIonization.PartialWave,1}
end 


"""
`PhotoExcitationAutoion.Pathway()`  
    ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified initial, intermediate 
        and final level.
"""
function Pathway()
    Pathway(Level(), Level(), Level(), 0., 0., EmProperty(0., 0.), EmProperty(0., 0.), PhotoEmission.Channel[], 
            AutoIonization.Channel[], PhotoIonization.PartialWave[] )
end


# `Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway)`  
#		... prepares a proper printout of the variable pathway::PhotoExcitationAutoion.Pathway.
function Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "excitEnergy                 $(pathway.excitEnergy)  ") 
    println(io, "electronEnergy              $(pathway.electronEnergy)  ")
    println(io, "partialCs:                  $(pathway.partialCs)  ")
    println(io, "qFano:                      $(pathway.qFano)  ")
    println(io, "excitChannels:              $(pathway.excitChannels)  ")
    println(io, "augerChannels:              $(pathway.augerChannels)  ")
    println(io, "photoChannels:              $(pathway.photoChannels)  ")
end



"""
`PhotoExcitationAutoion.computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, nm::Nuclear.Model, grid::Radial.Grid, 
                                                    nrContinuum::Int64, settings::PhotoExcitationAutoion.Settings)` 
    ... to compute all amplitudes and properties of the given pathway; a line::PhotoExcitationAutoion.Pathway is returned for which 
        the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                        settings::PhotoExcitationAutoion.Settings)
    # Compute all excitation channels
    neweChannels = PhotoEmission.Channel[]
    for eChannel in pathway.excitChannels
        amplitude   = PhotoEmission.amplitude(Absorption(), eChannel.multipole, eChannel.gauge, pathway.excitEnergy, 
                                                pathway.intermediateLevel, pathway.initialLevel, grid)
        push!( neweChannels, PhotoEmission.Channel( eChannel.multipole, eChannel.gauge, amplitude))
    end
    # Compute all AutoIonization decay channels
    newaChannels = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)
    for aChannel in pathway.augerChannels
        newnLevel   = Basics.generateLevelWithSymmetryReducedBasis(pathway.intermediateLevel, pathway.intermediateLevel.basis.subshells)
        newnLevel   = Basics.generateLevelWithExtraSubshell(Subshell(101, aChannel.kappa), newnLevel)
        newfLevel   = Basics.generateLevelWithSymmetryReducedBasis(pathway.finalLevel, pathway.finalLevel.basis.subshells)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(pathway.electronEnergy, Subshell(101, aChannel.kappa), newfLevel, nm, grid, contSettings)
        newcLevel   = Basics.generateLevelWithExtraElectron(cOrbital, aChannel.symmetry, newfLevel)
        newcChannel = AutoIonization.Channel( aChannel.kappa, aChannel.symmetry, phase, Complex(0.))
        amplitude = 1.0
        ## amplitude   = AutoIonization.amplitude("Coulomb", aChannel, newnLevel, newcLevel, grid)
        push!( newaChannels, AutoIonization.Channel( aChannel.kappa, aChannel.symmetry, phase, amplitude))
    end
    # Compute all photoionization channels
    ## The photoionization step now carries PARTIAL WAVES, each holding the channels of one kappa and, in
    ## each channel, one amplitude per multipole holding both gauges.  The amplitude here is still the stub
    ## 1.0im that this module has always used -- the real call is commented out two lines below, as it was.
    newpChannels = PhotoIonization.PartialWave[]
    for  pw in pathway.photoChannels
        ## amplitude = PhotoIonization.amplitude("photoionization", ma.multipole, gauge, pw.kappa, pw.phase,
        ##                                       pathway.excitEnergy, pathway.intermediateLevel, pathway.initialLevel, grid)
        newChannels = PhotoIonization.Channel[]
        for  ch in pw.channels
            newAmps = MultipoleAmplitude[]
            for  ma in ch.amplitudes    push!( newAmps, MultipoleAmplitude(ma.multipole, EmPropertyC(1.0im)) )    end
            push!( newChannels, PhotoIonization.Channel(ch.symmetry, newAmps) )
        end
        push!( newpChannels, PhotoIonization.PartialWave(pw.kappa, pw.energy, pw.phase, newChannels) )
    end
    #
    partialCs = EmProperty(-1., -1.)
    qFano     = EmProperty(-2., -2.)
    pathway = PhotoExcitationAutoion.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, pathway.excitEnergy, 
                                                pathway.electronEnergy, partialCs, qFano, neweChannels, newaChannels, newpChannels)
    return( pathway )
end



"""
`PhotoExcitationAutoion.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                        nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoExcitation.Settings; output=true)`  
    ... to compute the photo-excitation-autoionization amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{PhotoExcitationAutoion.Lines} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, 
                            grid::Radial.Grid, settings::PhotoExcitationAutoion.Settings; output=true)
    println("")
    printstyled("PhotoExcitationAutoion.computePathways(): The computation of photo-excitation-autoionization amplitudes starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    pathways = PhotoExcitationAutoion.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoExcitationAutoion.displayPathways(stdout, pathways, settings)    end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  pathway in pathways   maxEnergy = max(maxEnergy, pathway.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newPathways = PhotoExcitationAutoion.Pathway[]
    for  pathway in pathways
        push!( newPathways, PhotoExcitationAutoion.computeAmplitudesProperties(pathway, nm, grid, nrContinuum, settings) )
    end
    # Print all results to screen
    if  settings.calcPartialCs      PhotoExcitationAutoion.displayPartialCs(stdout, newPathways)         end
    if  settings.calcAngular        PhotoExcitationAutoion.displayAngular(stdout, newPathways)           end
    if  settings.calcFano           PhotoExcitationAutoion.displayFanoParameters(stdout, newPathways)    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    
        if  settings.calcPartialCs      PhotoExcitationAutoion.displayPartialCs(iostream, newPathways)         end
        if  settings.calcAngular        PhotoExcitationAutoion.displayAngular(iostream, newPathways)           end
        if  settings.calcFano           PhotoExcitationAutoion.displayFanoParameters(iostream, newPathways)    end
    end
    #
    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`PhotoExcitationAutoion.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                            settings::PhotoExcitationAutoion.Settings)`  
    ... to determine a list of photoexcitation-autoionization pathways between the levels from the given initial-, intermediate- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; an 
        Array{PhotoExcitationAutoion.Line,1} is returned. Apart from the level specification, all physical properties are set to zero 
        during the initialization process.  
"""
function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                            settings::PhotoExcitationAutoion.Settings)
    pathways = PhotoExcitationAutoion.Pathway[]
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                    eEnergy = nLevel.energy - iLevel.energy
                    aEnergy = nLevel.energy - fLevel.energy
                    pEnergy = fLevel.energy - iLevel.energy
                    if  eEnergy < 0.   ||   aEnergy < 0    continue    end
                    rSettings = PhotoEmission.Settings( settings.multipoles, settings.gauges, false, false, CorePolarization(), LineSelection(), 0., 0., 0., false)
                    eChannels = PhotoEmission.determineChannels(nLevel, iLevel, rSettings) 
                    aSettings = AutoIonization.Settings( false, false, false, false, LineSelection(), 0., 0., 0., settings.maxKappa, CoulombInteraction(), Multiplet())
                    aChannels = AutoIonization.determineChannels(fLevel, nLevel, aSettings) 
                    pSettings = PhotoIonization.Settings( settings.multipoles, settings.gauges, [pEnergy], false, false, false, false, 
                                                            LineSelection(), ExpStokes())
                    pChannels = PhotoIonization.determineChannels(fLevel, iLevel, pSettings) 
                    push!( pathways, PhotoExcitationAutoion.Pathway(iLevel, nLevel, fLevel, eEnergy, aEnergy, EmProperty(0., 0.), EmProperty(0., 0.), 
                                                                    eChannels, aChannels, pChannels) )
                end
            end
        end
    end
    return( pathways )
end


"""
`PhotoExcitationAutoion.displayPathways(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1}, settings::PhotoExcitationAutoion.Settings)`  
    ... to display a list of pathways and channels that have been selected due to the prior settings. A neat table of 
        all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPathways(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1}, settings::PhotoExcitationAutoion.Settings)
    nx = 158
    println(stream, " ")
    println(stream, "  Selected photo-excitation-autoionization pathways to be calculated:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=0);              
    sb = sb * TableStrings.center(30, "  e--i         e--f    "; na=0)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=2)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4) 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))    * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "    "
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  ech in pathway.excitChannels
            for  ach in pathway.augerChannels
                push!( kappaMultipoleSymmetryList, (ach.kappa, ech.multipole, ech.gauge, ach.symmetry) )
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        if      length(wa) == 0                       println(stream,  sa ) 
        elseif  length(wa) > 0    sb = sa * wa[1];    println(stream,  sb )                                 end  
        for  i = 2:length(wa)   sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )    end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    println(stream, "\n  + Given solid angles:  $(settings.solidAngles) ")
    if settings.calcPartialCs  println(stream, "  + Calculate the partial cross sections sigma(i-e-f) at the given solid angles.")           end
    if settings.calcAngular    println(stream, "  + Calculate the angular distribution of the Auger electrons at the given solid angles.")   end
    if settings.calcFano       println(stream, "  + Calculate the Fano parameters of the cross section resonance sigma(i-e-f).")             end
    #
    return( nothing )
end


"""
`PhotoExcitationAutoion.displayAngular(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but 
        nothing is returned otherwise.
"""
function  displayAngular(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})
    nx = 135
    println(stream, " ")
    println(stream, "  Angular parameters of electrons following the photo-excitation & autoionizationof resonances: ... to be adapted !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=0);              
    sb = sb * TableStrings.center(30, "  e--i         e--f    "; na=0)
    sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(30, "Cou -- Cross sections -- Bab"; na=2);       
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section")*"          "*
                                            TableStrings.inUnits("cross section"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=2)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))    * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "    "
        multipoles = EmMultipole[]
        for  ech in pathway.excitChannels
            multipoles = push!( multipoles, ech.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Coulomb))     * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Babushkin))   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end



"""
`PhotoExcitationAutoion.displayFanoParameters(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but 
        nothing is returned otherwise.
"""
function  displayFanoParameters(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})
    nx = 135
    println(stream, " ")
    println(stream, "  Fano parameters of photo-excited Auger resonances: ... to be adapted !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=0);              
    sb = sb * TableStrings.center(30, "  e--i         e--f    "; na=0)
    sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(30, "Cou -- Cross sections -- Bab"; na=2);       
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section")*"          "*
                                        TableStrings.inUnits("cross section"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=2)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))    * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "    "
        multipoles = EmMultipole[]
        for  ech in pathway.excitChannels
            multipoles = push!( multipoles, ech.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Coulomb))     * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Babushkin))   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end



"""
`PhotoExcitationAutoion.displayPartialCs(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but 
        nothing is returned otherwise.
"""
function  displayPartialCs(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})
    nx = 135
    println(stream, " ")
    println(stream, "  Partial photo-excitation & autoionization cross sections: ... to be adapted !!")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i   --   e   --   f"; na=1);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=0);              
    sb = sb * TableStrings.center(30, "  e--i         e--f    "; na=0)
    sa = sa * TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(30, "Cou -- Cross sections -- Bab"; na=2);       
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section")*"          "*
                                        TableStrings.inUnits("cross section"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=2)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))    * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.electronEnergy)) * "    "
        multipoles = EmMultipole[]
        for  ech in pathway.excitChannels
            multipoles = push!( multipoles, ech.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=3)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Coulomb))     * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.partialCs.Babushkin))   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
