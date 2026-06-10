
"""
`module  JAC.ResonantInelastic`  
... a submodel of JAC that contains all methods for computing Resonant-Inelastic X-ray Scattering (RIXS) observables,
    based on initial, intermediate and final-state multiplets.
"""
module ResonantInelastic


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..PhotoEmission, ..PhotoExcitation, 
      ..Radial, ..TableStrings

"""
`struct  ResonantInelastic.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing RIXS lines.

    + multipoles              ::Array{EmMultipoles}     ... Specifies the (radiat. field) multipoles to be included.
    + gauges                  ::Array{UseGauge}         ... Gauges to be included into the computations.
    + calcRixsCI              ::Bool                    
        ... True, if a modified configuration-interaction (CI) scheme is applied in order to compute RIXS  amplitude,
            and false otherwise.
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before comput.
    + ciEnhancement           ::Float64                 ... (external) enhancement factor for the RixsCI method.      
    + width                   ::Float64                 ... width of intermediate levels [user-defined energies].      
    + pathwaySelection        ::PathwaySelection        ... Specifies the selected pathways, if any.
    + omegaInInterval         ::NamedTuple{(:lower, :upper), Tuple{Float64, Float64}}
        ... A NamedTuple (lower=.., upper=..) to indicate the lower and upper incoming photon energies.
    + omegaOutInterval        ::NamedTuple{(:lower, :upper), Tuple{Float64, Float64}}
        ... A NamedTuple (lower=.., upper=..) to indicate the lower and upper outgoing photon energies.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge}
    calcRixsCI                ::Bool                    
    printBefore               ::Bool 
    ciEnhancement             ::Float64
    width                     ::Float64       
    pathwaySelection          ::PathwaySelection 
    omegaInInterval           ::NamedTuple{(:lower, :upper), Tuple{Float64, Float64}}
    omegaOutInterval          ::NamedTuple{(:lower, :upper), Tuple{Float64, Float64}}
end 


"""
`ResonantInelastic.Settings()`  ... constructor for the default values of radiative line computations
"""
function Settings()
    Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], false, false, 1.0, 0., PathwaySelection(), 
             (lower=0., upper=0.), (lower=0., upper=0.) )
end


"""
`ResonantInelastic.Settings(set::ResonantInelastic.Settings;`

        multipoles::=..,        gauges=..,                calcRixsCI=..,            printBefore=..,
        ciEnhancement=..,       width=..,                 pathwaySelection=..,      
        omegaInInterval=..,       omegaOutInterval=..) 
                    
    ... constructor for modifying the given ResonantInelastic.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ResonantInelastic.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    gauges::Union{Nothing,Array{UseGauge}}=nothing,
    calcRixsCI::Union{Nothing,Bool}=nothing,                    printBefore::Union{Nothing,Bool}=nothing,
    ciEnhancement::Union{Nothing,Float64}=nothing,              width::Union{Nothing,Float64}=nothing,      
    pathwaySelection::Union{Nothing,PathwaySelection}=nothing, 
    omegaInInterval::Union{Nothing,NamedTuple{(:lower, :upper), Tuple{Float64, Float64}} }=nothing,            
    omegaOutInterval::Union{Nothing,NamedTuple{(:lower, :upper), Tuple{Float64, Float64}} }=nothing)
    
    if  multipoles          == nothing   multipolesx          = set.multipoles              else  multipolesx          = multipoles            end 
    if  gauges              == nothing   gaugesx              = set.gauges                  else  gaugesx              = gauges                end 
    if  calcRixsCI          == nothing   calcRixsCIx          = set.calcRixsCI              else  calcRixsCIx          = calcRixsCI            end 
    if  printBefore         == nothing   printBeforex         = set.printBefore             else  printBeforex         = printBefore           end 
    if  ciEnhancement       == nothing   ciEnhancementx       = set.ciEnhancement           else  ciEnhancementx       = ciEnhancement         end 
    if  width               == nothing   widthx               = set.width                   else  widthx               = width                 end 
    if  pathwaySelection    == nothing   pathwaySelectionx    = set.pathwaySelection        else  pathwaySelectionx    = pathwaySelection      end 
    if  omegaInInterval     == nothing   omegaInIntervalx     = set.omegaInInterval         else  omegaInIntervalx     = omegaInInterval       end 
    if  omegaOutInterval    == nothing   omegaOutIntervalx    = set.omegaOutInterval        else  omegaOutIntervalx    = omegaOutInterval      end 
    
    Settings( multipolesx, gaugesx, calcRixsCIx, printBeforex, ciEnhancementx, widthx, pathwaySelectionx, omegaInIntervalx, omegaOutIntervalx)
end


# `Base.show(io::IO, settings::ResonantInelastic.Settings)`  ... prepares a proper printout of the variable settings::ResonantInelastic.Settings.
function Base.show(io::IO, settings::ResonantInelastic.Settings) 
    println(io, "multipoles:             $(settings.multipoles)  ")
    println(io, "gauges:                 $(settings.gauges)  ")
    println(io, "calcRixsCI:             $(settings.calcRixsCI)  ")
    println(io, "printBefore:            $(settings.printBefore)  ")
    println(io, "ciEnhancement:          $(settings.ciEnhancement)  ")
    println(io, "width:                  $(settings.width)  ")
    println(io, "pathwaySelection:       $(settings.pathwaySelection)  ")
    println(io, "omegaInInterval:        $(settings.omegaInInterval)  ")
    println(io, "omegaOutInterval:       $(settings.omegaOutInterval)  ")
end


"""
`struct  ResonantInelastic.Pathway`  
    ... defines a type for a photoexcitation-photoemission pathway, based on a number of 
        radiative line that may include the definition of sublines and their corresponding amplitudes.

    + initialLevel      ::Level                             ... initial-(state) level
    + intermediateLevel ::Level                             ... intermediate-(state) level
    + finalLevel        ::Level                             ... final-(state) level
    + omegaIn           ::Float64                           ... Transition frequency of excitation.
    + omegaOut          ::Float64                           ... Transition frequency of emission.
    + intermediateGamma ::Float64                           ... Widths of the intermediate level
    + productF          ::EmProperty                        ... F (product of amplitudes is assumed to be real)
    + relativeCS        ::EmProperty                        ... relative CS = |F|^2
    + excitationChannels::Array{PhotoEmission.Channel,1}    ... List of incoming channels.
    + emissionChannels  ::Array{PhotoEmission.Channel,1}    ... List of outgoing RIXS channels.
"""
struct  Pathway
    initialLevel        ::Level
    intermediateLevel   ::Level
    finalLevel          ::Level
    omegaIn             ::Float64 
    omegaOut            ::Float64 
    intermediateGamma   ::Float64 
    productF            ::EmProperty 
    relativeCS          ::EmProperty
    excitationChannels  ::Array{PhotoEmission.Channel,1}
    emissionChannels    ::Array{PhotoEmission.Channel,1}
end 


# `Base.show(io::IO, pathway::ResonantInelastic.Pathway)`  ... prepares a proper printout of the variable pathway::ResonantInelastic.Pathway
function Base.show(io::IO, pathway::ResonantInelastic.Pathway) 
    println(io, "initialLevel:         $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:    $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:           $(pathway.finalLevel)  ")
    println(io, "omegaIn:              $(pathway.omegaIn)  ")
    println(io, "omegaOut:             $(pathway.omegaOut)  ")
    println(io, "intermediateGamma:    $(pathway.intermediateGamma)  ")
    println(io, "productF:             $(pathway.productF)  ")
    println(io, "relativeCS:           $(pathway.relativeCS)  ")
    println(io, "channels:             $(pathway.channels)  ")
end


"""
`ResonantInelastic.computeAmplitudesProperties(pathway::ResonantInelastic.Pathway, grid::Radial.Grid, 
                                               settings::ResonantInelastic.Settings)` 
    ... to compute all amplitudes and properties of the given RIXS pathway; a newPathway::ResonantInelastic.Pathway 
        is returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(pathway::ResonantInelastic.Pathway, grid::Radial.Grid, 
                                      settings::ResonantInelastic.Settings)
    productF = EmProperty(0.);   relativeCS = EmProperty(0.)
    ## println(">> pathway: $(pathway.initialLevel.index)--$(pathway.intermediateLevel.index)--$(pathway.finalLevel.index) ...")
    #
    Defaults.setDefaults("relativistic subshell list", pathway.intermediateLevel.basis.subshells; printout=false)
    intermediateLevel = deepcopy(pathway.intermediateLevel)
    initialLevel      = deepcopy(pathway.initialLevel)
    
    #  The particular application requires to decide, which levels need to be "mixed with enhanced factor"
    
    ## Farid + Ruben (June 2025): Interference in the RIXS out-channel
    ## if settings.calcRixsCI    intermediateLevel =  Basics.modifyLevelMixing(intermediateLevel, settings.ciEnhancement)   end 
    
    ## Ruben + Johan (Nov 2025): RIXS signals with double excitations
    if settings.calcRixsCI       initialLevel      =  Basics.modifyLevelMixing(initialLevel, settings.ciEnhancement)        end 
    
    # Compute the amplitudes of the excitation and emission channels separately
    newExChannels = PhotoEmission.Channel[];    newEmChannels = PhotoEmission.Channel[]
    for exChannel in pathway.excitationChannels
        amplitude    = PhotoEmission.amplitude("emission", exChannel.multipole, exChannel.gauge, pathway.omegaIn, 
                                               intermediateLevel, initialLevel, grid, display=false, printout=false)
        newExChannel = PhotoEmission.Channel( exChannel.multipole, exChannel.gauge, amplitude)
        push!( newExChannels, newExChannel)
    end
    
    for emChannel in pathway.emissionChannels
        amplitude    = PhotoEmission.amplitude("emission", emChannel.multipole, emChannel.gauge, pathway.omegaOut, 
                                               intermediateLevel, pathway.finalLevel, grid, display=false, printout=false)
        newEmChannel = PhotoEmission.Channel( emChannel.multipole, emChannel.gauge, amplitude)
        push!( newEmChannels, newEmChannel)
    end

    # Choose a proper width of the intermediate level from the settings
    width = Defaults.convertUnits("energy: to atomic", settings.width)
    
    # Combine the amplitudes of the proper gauges to form the relativeCS
    ampCou = ampBab = ComplexF64(0.);    ampEx = ampEm = 0.
    for exChannel in newExChannels
        for emChannel in newEmChannels
            ampEx = ampEx + abs(exChannel.amplitude)  # Just to understand which part is always zero
            ampEm = ampEm + abs(emChannel.amplitude)  # Just to understand which part is always zero
            if  exChannel.gauge == Basics.Coulomb  &&  emChannel.gauge == Basics.Babushkin     continue    end
            if  emChannel.gauge == Basics.Coulomb  &&  exChannel.gauge == Basics.Babushkin     continue    end
            if      exChannel.gauge == Basics.Coulomb    &&  emChannel.gauge == Basics.Coulomb
                ampCou = ampCou + emChannel.amplitude * conj(exChannel.amplitude)
            elseif  exChannel.gauge == Basics.Babushkin  &&  emChannel.gauge == Basics.Babushkin
                ampBab = ampBab + emChannel.amplitude * conj(exChannel.amplitude)
            else
                ampCou = ampCou + emChannel.amplitude * conj(exChannel.amplitude)
                ampBab = ampBab + emChannel.amplitude * conj(exChannel.amplitude)
            end
        end
    end

    @show ampCou, ampBab, ampEx, ampEm
    productF   = EmProperty(    ampCou.re,     ampBab.re )
    relativeCS = EmProperty(abs(ampCou)^2, abs(ampBab)^2)
    pathway    = ResonantInelastic.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, 
                                            pathway.omegaIn, pathway.omegaOut, width, productF, relativeCS, 
                                            newExChannels, newEmChannels)
    return( pathway )
end


"""
`ResonantInelastic.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                   grid::Radial.Grid, settings::ResonantInelastic.Settings; output=true)`  
    ... to compute the resonant-inelastic x-ray scattering observables as characterized by the associated settings.
        A list of pathways::Array{ResonantInelastic.Pathway} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                          grid::Radial.Grid, settings::ResonantInelastic.Settings; output=true) 
    # Define a common subshell list for all three multiplets
    subshellList = Basics.generate("subshells: ordered list for two bases", intermediateMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    subshellList = Basics.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, intermediateMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    println("")
    printstyled("ResonantInelastic.computeLines(): The computation of the pathway properties starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
    println("")
    pathways = ResonantInelastic.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    # Display all selected pathways before the computations start
    if  settings.printBefore    ResonantInelastic.displayPathways(stdout, pathways)    end
    # Calculate all amplitudes and requested properties
    newPathways = ResonantInelastic.Pathway[]
    for  pathway in pathways
        newPathway = ResonantInelastic.computeAmplitudesProperties(pathway, grid, settings) 
        push!( newPathways, newPathway)
    end
    # Print all results to screen
    ResonantInelastic.displayResults(stdout, newPathways, settings)
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ResonantInelastic.displayResults(iostream, newPathways, settings)   end
    #
    if    output    return( newPathways )
    else            return( nothing )
    end
end


"""
`ResonantInelastic.determineExcitationChannels(intermediateLevel::Level, initialLevel::Level, 
                                               settings::ResonantInelastic.Settings)` 
    ... to determine a list of PhotoExcitation.Channel's for a photoexcitation transitions from the initial to an 
        intermediate level, and by taking into account the particular settings of for this computation;  
        an Array{PhotoExcitation.Channel,1} is returned.
"""
function determineExcitationChannels(intermediateLevel::Level, initialLevel::Level, settings::ResonantInelastic.Settings)
    channels = PhotoEmission.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symn, mp, symi)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)           push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end
    
    return( channels )  
end


"""
`ResonantInelastic.determineEmissionChannels(finalLevel::Level, intermediateLevel::Level, settings::ResonantInelastic.Settings)` 
    ... to determine a list of PhotoEmission.Channel for the photon transitions from the intermediate and to a final level, and by 
        taking into account the particular settings of for this computation; an Array{PhotoEmission.Channel,1} is returned.
"""
function determineEmissionChannels(finalLevel::Level, intermediateLevel::Level, settings::ResonantInelastic.Settings)
    channels = PhotoEmission.Channel[];   
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symn, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)           push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end

    return( channels )  
end


"""
`ResonantInelastic.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                     settings::ResonantInelastic.Settings)`  
    ... to determine a list of RIXS pathways between the levels from the given initial-, intermediate- and final-state 
        multiplets, and by taking into account the particular selections and settings for this computation; 
        an Array{ResonantInelastic.Pathway,1} is returned. Apart from the level specification, all physical properties 
        are set to zero during the initialization process.  
"""
function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                            settings::ResonantInelastic.Settings)
    pathways   = ResonantInelastic.Pathway[]
    omegaInau  = ( lower = Defaults.convertUnits("energy: to atomic", settings.omegaInInterval.lower), 
                   upper = Defaults.convertUnits("energy: to atomic", settings.omegaInInterval.upper) )
    omegaOutau = ( lower = Defaults.convertUnits("energy: to atomic", settings.omegaOutInterval.lower), 
                   upper = Defaults.convertUnits("energy: to atomic", settings.omegaOutInterval.upper) )
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                    omegaIn  = nLevel.energy - iLevel.energy
                    omegaOut = nLevel.energy - fLevel.energy
                    if  omegaIn  < omegaInau.lower   ||  omegaInau.upper  < omegaIn           ||
                        omegaOut < omegaOutau.lower  ||  omegaOutau.upper < omegaOut    continue    end
                    exChannels = ResonantInelastic.determineExcitationChannels(nLevel, iLevel, settings) 
                    emChannels = ResonantInelastic.determineEmissionChannels(  fLevel, nLevel, settings)
                    if  length(exChannels) * length(emChannels) > 0
                        push!( pathways, ResonantInelastic.Pathway(iLevel, nLevel, fLevel, omegaIn, omegaOut, 0., 
                                              EmProperty(0., 0.),  EmProperty(0., 0.), exChannels, emChannels) )
                    end
                end
            end
        end
    end
    return( pathways )
end


"""
`ResonantInelastic.displayPathways(stream::IO, pathways::Array{ResonantInelastic.Pathway,1})`  
    ... to display a list of pathways and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPathways(stream::IO, pathways::Array{ResonantInelastic.Pathway,1})
    nx = 180
    println(stream, " ")
    println(stream, "  Selected RIXS pathways:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(23, "Levels"; na=4);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=4);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=3);
    sa = sa * TableStrings.center(26, "Energies  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(26, "omega-in        omega-out"; na=5)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                       msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                       fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                  pathway.finalLevel.index); na=5)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.omegaIn))    * "   "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pathway.omegaOut))   * "    "
        multipolesSymmetryList = Tuple{EmMultipole,LevelSymmetry,EmMultipole,EmGauge}[]
        for  exChannel in pathway.excitationChannels
            for  emChannel in pathway.emissionChannels
                if  exChannel.gauge == Basics.Coulomb  &&  emChannel.gauge == Basics.Babushkin    continue    end
                if  emChannel.gauge == Basics.Coulomb  &&  exChannel.gauge == Basics.Babushkin    continue    end
                if      exChannel.gauge == Basics.Coulomb   ||  emChannel.gauge == Basics.Coulomb      gauge = Basics.Coulomb
                elseif  exChannel.gauge == Basics.Babushkin ||  emChannel.gauge == Basics.Babushkin    gauge = Basics.Babushkin
                else                                                                                   gauge = Basics.Magnetic
                end
                push!( multipolesSymmetryList, (exChannel.multipole, msym, emChannel.multipole, gauge) )
            end
        end
        wa = TableStrings.multipolesSymmetryTupels(105, multipolesSymmetryList)
        if  length(wa) > 0    sb = sa * wa[1];    println(stream,  sb )    end  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "\n>> A total of $(length(pathways)) RIXS pathways will be calculated. \n")
    #
    return( nothing )
end


"""
`ResonantInelastic.displayResults(stream::IO, pathways::Array{ResonantInelastic.Pathway,1},
                                  settings::ResonantInelastic.Settings)`  
    ... to list all results, energies, relative amplitudes, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(stream::IO, pathways::Array{ResonantInelastic.Pathway,1},  settings::ResonantInelastic.Settings)
    nx = 196
    #
    println(stream, " ")
    println(stream, "  RIXS energies, widths, relative F = Re (<f|Mout|m> <m|Min|i>) values" * 
                    "  and relative cross sections (without omega/Gamma-dependent denominator):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2);
    sa = sa * TableStrings.center(50, "Energies  " * TableStrings.inUnits("energy");        na=4);              
    sb = sb * TableStrings.center(50, " omega_in     omega_out    in - out     width(m)  "; na=1)
    sa = sa * TableStrings.center(10, "Multipoles"; na=6)        
    sb = sb * TableStrings.center(10, "in;  out  "; na=5)
    sa = sa * TableStrings.center(36, " relative |F|^2 "; na=3);   
    sb = sb * TableStrings.center(36, "  Cou-Re (<f|Mout|m> <m|Min|i>)-Bab"; na=3)
    sa = sa * TableStrings.center(36, " relative CS "; na=2);   
    sb = sb * TableStrings.center(36, "  Cou-|<f|Mout|m>  <m|Min|i>|^2-Bab"; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";      isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                           pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        wm = pathway.omegaIn - pathway.omegaOut
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.omegaIn))           * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.omegaOut))          * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wm))                        * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.intermediateGamma)) * "   "
        exMultipoles = EmMultipole[];   emMultipoles = EmMultipole[]
        for  ch in pathway.excitationChannels   push!( exMultipoles, ch.multipole)   end;   exMultipoles = unique(exMultipoles)
        for  ch in pathway.emissionChannels     push!( emMultipoles, ch.multipole)   end;   emMultipoles = unique(emMultipoles)
        mpString = TableStrings.multipoleList(exMultipoles) * ";  " * TableStrings.multipoleList(emMultipoles) * "               "
        sa = sa * TableStrings.flushleft(18, mpString[1:18];  na=0)
        sa = sa * @sprintf("% 2.6e", pathway.productF.Coulomb)     * "       "
        sa = sa * @sprintf("% 2.6e", pathway.productF.Babushkin)   * "      "
        sa = sa * @sprintf("% 2.6e", pathway.relativeCS.Coulomb)   * "       "
        sa = sa * @sprintf("% 2.6e", pathway.relativeCS.Babushkin) * " "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


end # module

