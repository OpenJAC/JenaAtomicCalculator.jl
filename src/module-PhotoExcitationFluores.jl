
"""
`module  JAC.PhotoExcitationFluores`  
... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross sections and rates.
"""
module PhotoExcitationFluores 


using Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..PhotoEmission, ..Statistical, ..TableStrings

"""
`struct  PhotoExcitationFluores.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing photon-impact excitation-autoionization 
        pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

    + multipoles              ::Array{EmMultipole,1}   ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                  ::Array{UseGauge,1}      ... Specifies the gauges to be included into the computations.
    + calcPhotonDm            ::Bool                   ... Calculate the reduced density matrix of the fluorescence photon.
    + calcAngular             ::Bool                   ... Calculate the angular distribution of the fluorescence photon.
    + calcStokes              ::Bool                   ... Calculate the Stokes parameters of the fluorescence photon.
                                                            In all these (three) cases, the atoms and incident plane-wave photons are
                                                            assumed to be unpolarized initially.
    + printBefore             ::Bool                   ... True, if all energies and lines are printed before their evaluation.
    + incidentStokes          ::ExpStokes              ... Stokes parameters of the incident radiation.
    + solidAngles             ::Array{SolidAngle,1}    ... List of solid angles [(theta_1, pho_1), ...].  
    + photonEnergyShift       ::Float64                ... An overall energy shift for all photon energies.
    + pathwaySelection        ::PathwaySelection       ... Specifies the selected levels/pathways, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge,1} 
    calcPhotonDm              ::Bool
    calcAngular               ::Bool
    calcStokes                ::Bool 
    printBefore               ::Bool  
    incidentStokes            ::ExpStokes
    solidAngles               ::Array{SolidAngle,1}
    photonEnergyShift         ::Float64  
    pathwaySelection          ::PathwaySelection
    end 


"""
`PhotoExcitationFluores.Settings()`  
    ... constructor for the default values of photon-impact excitation-autoionizaton settings.
"""
function Settings()
    Settings( Basics.EmMultipole[], UseGauge[], false, false, false, false, ExpStokes(), SolidAngle[], 0., PathwaySelection() )
end


"""
`PhotoExcitationFluores.Settings(set::PhotoExcitationFluores.Settings; keywords...)`
    ... keyword constructor to modify selected fields of a PhotoExcitationFluores.Settings object;
        for fields not given, the values of set are used.
"""
function Settings(set::PhotoExcitationFluores.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,           gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
        calcPhotonDm::Union{Nothing,Bool}=nothing,                         calcAngular::Union{Nothing,Bool}=nothing,
        calcStokes::Union{Nothing,Bool}=nothing,                           printBefore::Union{Nothing,Bool}=nothing,
        incidentStokes::Union{Nothing,ExpStokes}=nothing,                  solidAngles::Union{Nothing,Array{SolidAngle,1}}=nothing,
        photonEnergyShift::Union{Nothing,Float64}=nothing,                 pathwaySelection::Union{Nothing,PathwaySelection}=nothing)
    if  isnothing(multipoles)          multipolsx          = set.multipoles         else   multipolsx          = multipoles         end
    if  isnothing(gauges)              gaugesx             = set.gauges             else   gaugesx             = gauges             end
    if  isnothing(calcPhotonDm)        calcPhotonDmx       = set.calcPhotonDm       else   calcPhotonDmx       = calcPhotonDm       end
    if  isnothing(calcAngular)         calcAngularx        = set.calcAngular        else   calcAngularx        = calcAngular        end
    if  isnothing(calcStokes)          calcStokesx         = set.calcStokes         else   calcStokesx         = calcStokes         end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end
    if  isnothing(incidentStokes)      incidentStokesx     = set.incidentStokes     else   incidentStokesx     = incidentStokes     end
    if  isnothing(solidAngles)         solidAnglesx        = set.solidAngles        else   solidAnglesx        = solidAngles        end
    if  isnothing(photonEnergyShift)   photonEnergyShiftx  = set.photonEnergyShift  else   photonEnergyShiftx  = photonEnergyShift  end
    if  isnothing(pathwaySelection)    pathwaySelectionx   = set.pathwaySelection   else   pathwaySelectionx   = pathwaySelection   end

    Settings( multipolsx, gaugesx, calcPhotonDmx, calcAngularx, calcStokesx, printBeforex,
              incidentStokesx, solidAnglesx, photonEnergyShiftx, pathwaySelectionx )
end


# `Base.show(io::IO, settings::PhotoExcitationFluores.Settings)`
# 	 ... prepares a proper printout of the variable settings::PhotoExcitationFluores.Settings.  
function Base.show(io::IO, settings::PhotoExcitationFluores.Settings) 
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "calcPhotonDm:            $(settings.calcPhotonDm)  ")
    println(io, "calcAngular:             $(settings.calcAngular)  ")
    println(io, "calcStokes:              $(settings.calcStokes)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "incidentStokes:          $(settings.incidentStokes)  ")
    println(io, "solidAngles:             $(settings.solidAngles)  ")
    println(io, "photonEnergyShift:       $(settings.photonEnergyShift)  ")
    println(io, "pathwaySelection:        $(settings.pathwaySelection)  ")
end


"""
`struct  PhotoExcitationFluores.Pathway`  
    ... defines a type for a photon-impact excitation pathway that may include the definition of different excitation and 
        autoionization channels and their corresponding amplitudes.

    + initialLevel        ::Level                  ... initial-(state) level
    + intermediateLevel   ::Level                  ... intermediate-(state) level
    + finalLevel          ::Level                  ... final-(state) level
    + excitEnergy         ::Float64                ... photon excitation energy of this pathway
    + fluorEnergy         ::Float64                ... photon fluorescence energy of this pathway
    + crossSection        ::EmProperty             ... total cross section of this pathway
    + excitChannels       ::Array{MultipoleAmplitude,1}     ... excitation amplitudes, one per multipole, both gauges.
    + fluorChannels       ::Array{MultipoleAmplitude,1}     ... fluorescence amplitudes, one per multipole, both gauges.
"""
struct  Pathway
    initialLevel          ::Level
    intermediateLevel     ::Level
    finalLevel            ::Level
    excitEnergy           ::Float64
    fluorEnergy           ::Float64
    crossSection          ::EmProperty
    excitChannels         ::Array{MultipoleAmplitude,1}
    fluorChannels         ::Array{MultipoleAmplitude,1}
end 


"""
`PhotoExcitationFluores.Pathway()`  
    ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified initial, 
        intermediate and final level.
"""
function Pathway()
    Pathway(Level(), Level(), Level(), 0., 0., 0., false, PhotoExcitationFluores.Channel[] )
end


# `Base.show(io::IO, pathway::PhotoExcitationFluores.Pathway)` 
#		 ... prepares a proper printout of the variable pathway::PhotoExcitationFluores.Pathway.
function Base.show(io::IO, pathway::PhotoExcitationFluores.Pathway) 
    println(io, "initialLevel:               $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:                 $(pathway.finalLevel)  ")
    println(io, "excitEnergy                 $(pathway.excitEnergy)  ") 
    println(io, "fluorEnergy                 $(pathway.fluorEnergy)  ")
    println(io, "crossSection:               $(pathway.crossSection)  ")
    println(io, "excitChannels:              $(pathway.excitChannels)  ")
    println(io, "fluorChannels:              $(pathway.fluorChannels)  ")
end


"""
`PhotoExcitationFluores.computeAmplitudesProperties(pathway::PhotoExcitationFluores.Pathway, grid::Radial.Grid, 
                                                        settings::PhotoExcitationFluores.Settings)` 
    ... to compute all amplitudes and properties of the given line; a pathway::PhotoExcitationFluores.Pathway is 
        returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(pathway::PhotoExcitationFluores.Pathway, grid::Radial.Grid, settings::PhotoExcitationFluores.Settings)
    # Compute all excitation channels
    neweChannels = MultipoleAmplitude[]
    for  ma in pathway.excitChannels
        mp = ma.multipole
        if  string(mp)[1] == 'E'
            ampC = PhotoEmission.amplitude(Absorption(), mp, Basics.Coulomb,   pathway.excitEnergy,
                                            pathway.intermediateLevel, pathway.initialLevel, grid)
            ampB = PhotoEmission.amplitude(Absorption(), mp, Basics.Babushkin, pathway.excitEnergy,
                                            pathway.intermediateLevel, pathway.initialLevel, grid)
            amp  = EmPropertyC(ampC, ampB)
        else
            ampM = PhotoEmission.amplitude(Absorption(), mp, Basics.Magnetic,  pathway.excitEnergy,
                                            pathway.intermediateLevel, pathway.initialLevel, grid)
            amp  = EmPropertyC(ampM)
        end
        push!( neweChannels, MultipoleAmplitude(mp, amp))
    end
    # Compute all fluorescence channels
    newfChannels = MultipoleAmplitude[]
    for  ma in pathway.fluorChannels
        mp = ma.multipole
        if  string(mp)[1] == 'E'
            ampC = PhotoEmission.amplitude(Emission(), mp, Basics.Coulomb,   pathway.fluorEnergy,
                                            pathway.finalLevel, pathway.intermediateLevel, grid)
            ampB = PhotoEmission.amplitude(Emission(), mp, Basics.Babushkin, pathway.fluorEnergy,
                                            pathway.finalLevel, pathway.intermediateLevel, grid)
            amp  = EmPropertyC(ampC, ampB)
        else
            ampM = PhotoEmission.amplitude(Emission(), mp, Basics.Magnetic,  pathway.fluorEnergy,
                                            pathway.finalLevel, pathway.intermediateLevel, grid)
            amp  = EmPropertyC(ampM)
        end
        push!( newfChannels, MultipoleAmplitude(mp, amp))
    end
    crossSection = EmProperty(-1., -1.)
    pathway = PhotoExcitationFluores.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, pathway.excitEnergy, 
                                                pathway.fluorEnergy, crossSection, true, neweChannels, newfChannels)
    return( pathway )
end


"""
`PhotoExcitationFluores.computePhotonDm(pathway::PhotoExcitationFluores.Pathway, settings::PhotoExcitationFluores.Settings)`  
    ... computes the density matrix of the fluorescence photon, in Coulomb gauge, for the given pathway and for all
        selected solid angles.  An `dmArray::Array{ComplexF64,2}` of shape (n_solidAngles, 6) is returned, each row
        holding (theta, phi, rho(1,1), rho(1,-1), rho(-1,1), rho(-1,-1)).

        THE EXCITED ENSEMBLE IS NOW COMPUTED rather than assumed.  What decides the polarization of the fluorescence is
        the alignment left behind by the excitation step, and that is taken from
        `Statistical.computeTensorsPhotoExcitation` on the pathway's excitation channels, so this module and
        `PhotoExcitation` describe the same excited atom.  Its sublevel density matrix is recovered with
        `Statistical.densityMatrix`.

        THE PHOTON AMPLITUDE IS ASSEMBLED EXPLICITLY, with no Racah reduction, for the same reason as in the excitation
        step: at these angular momenta the sums cost nothing and their phases can be checked.  A photon emitted into
        (theta,phi) with helicity lambda is described by rotating the multipole operator out of the beam frame,

            M(M_e -> M_f, lambda)  =  SUM_(L)  i^(L+p) lambda^p A_L  conj(D^L_(M_e-M_f, lambda)(phi,theta,0))
                                                <J_f M_f, L (M_e-M_f) | J_e M_e>

        the projection being fixed at mu = M_e - M_f by the Clebsch-Gordan coefficient, and the photon density matrix
        follows by summing over the unobserved final sublevels,

            rho^gamma_(lambda lambda')  =  SUM_(M_f) SUM_(M_e,M_e') rho_e(M_e,M_e') M(M_e,M_f,lambda) conj(M(M_e',M_f,lambda')) .
"""
function  computePhotonDm(pathway::PhotoExcitationFluores.Pathway, settings::PhotoExcitationFluores.Settings)
    Je = pathway.intermediateLevel.J;             Jf = pathway.finalLevel.J
    Jex = AngularMomentum.oneJ(Je);               Jfx = AngularMomentum.oneJ(Jf)
    MeList = AngularMomentum.m_values(Je);        MfList = AngularMomentum.m_values(Jf)

    # The alignment left by the excitation step, and the sublevel density matrix it stands for
    eKey    = Basics.LevelKey( LevelSymmetry(Je, pathway.intermediateLevel.parity), pathway.intermediateLevel.index,
                               pathway.intermediateLevel.energy, 1.)
    tensors = Statistical.computeTensorsPhotoExcitation( Int64(round(2*Jex)), pathway.initialLevel.J, eKey,
                                                         pathway.excitChannels, settings.incidentStokes, Basics.Coulomb)
    rhoE    = length(tensors) == 0 ? Dict{Tuple{AngularM64,AngularM64},ComplexF64}() : Statistical.densityMatrix(tensors)

    # The emission amplitude into (theta,phi) with helicity lambda
    function emissionAmplitude(Me::AngularM64, Mf::AngularM64, lambda::Int64, solidAngle::SolidAngle)
        mu = AngularMomentum.oneM(Me) - AngularMomentum.oneM(Mf)
        wa = ComplexF64(0.)
        for  ch in pathway.fluorChannels
            L = ch.multipole.L;    p = ch.multipole.electric ? 1 : 0
            if  abs(mu) > L    continue    end
            cg = AngularMomentum.ClebschGordan(Jfx, AngularMomentum.oneM(Mf), 1.0*L, mu, Jex, AngularMomentum.oneM(Me))
            if  cg == 0.    continue    end
            wa = wa + (1.0im)^(L + p) * lambda^p * ch.amplitude.Coulomb * cg *
                      conj( AngularMomentum.Wigner_DFunction(L, mu, 1.0*lambda, solidAngle.phi, solidAngle.theta, 0.0) )
        end

        return( wa )
    end

    function photonDmElement(solidAngle::SolidAngle, lambda::Int64, lambdap::Int64)
        wa = ComplexF64(0.)
        for  Mf in MfList
            for  Me in MeList
                for  Mep in MeList
                    if  !haskey(rhoE, (Me,Mep))    continue    end
                    wa = wa + rhoE[(Me,Mep)] * emissionAmplitude(Me, Mf, lambda, solidAngle) *
                              conj( emissionAmplitude(Mep, Mf, lambdap, solidAngle) )
                end
            end
        end

        return( 2pi * wa )
    end

    dmArray = zeros(ComplexF64, length(settings.solidAngles), 6)
    for  (s, solidAngle) in enumerate(settings.solidAngles)
        dmArray[s,1] = solidAngle.theta;                    dmArray[s,2] = solidAngle.phi
        dmArray[s,3] = photonDmElement(solidAngle,  1,  1); dmArray[s,4] = photonDmElement(solidAngle,  1, -1)
        dmArray[s,5] = photonDmElement(solidAngle, -1,  1); dmArray[s,6] = photonDmElement(solidAngle, -1, -1)
    end

    return( dmArray )
end


"""
`PhotoExcitationFluores.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                        grid::Radial.Grid, settings::PhotoExcitationFluores.Settings; output=true)`  
    ... to compute the photo-excitation-fluorescence amplitudes and all properties as requested by the given settings. A list of
        pathways::Array{PhotoExcitationFluores.Pathway,1} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                            settings::PhotoExcitationFluores.Settings; output=true)
    println("")
    printstyled("PhotoExcitationFluores.computePathways(): The computation of photo-excitation-fluorescence amplitudes starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    pathways = PhotoExcitationFluores.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    # Display all selected pathways before the computations start
    if  settings.printBefore    PhotoExcitationFluores.displayPathways(stdout, pathways, settings)    end
    # Calculate all amplitudes and requested properties
    newPathways = PhotoExcitationFluores.Pathway[]
    for  pathway in pathways
        newPathway = PhotoExcitationFluores.computeAmplitudesProperties(pathway, grid, settings) 
        push!( newPathways, newPathway)
    end
    # Print all results to screen
    PhotoExcitationFluores.displayResults(stdout, newPathways, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoExcitationFluores.displayResults(iostream, newPathways, settings)     end
    #
    # Calculate the photon density matrix of the fluorescence photon if needed and display all requested properties
    # These calculations are done in turn for all selected pathways
    if  settings.calcPhotonDm  ||  settings.calcAngular  ||  settings.calcStokes
        for  pathway in newPathways
            dmArray = PhotoExcitationFluores.computePhotonDm(pathway, settings) 
            PhotoExcitationFluores.displayPhotonDm(stdout, pathway, dmArray, settings)
            if  printSummary   PhotoExcitationFluores.displayPhotonDm(stdout, pathway, dmArray, settings)     end
        end
    end
    
    if    output    return( pathways )
    else            return( nothing )
    end
end


"""
`PhotoExcitationFluores.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                            settings::PhotoExcitationFluores.Settings)`  
    ... to determine a list of dielectronic-recombination pathways between the levels from the given initial-, intermediate- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{PhotoExcitationFluores.Pathway,1} is returned. Apart from the level specification, all physical properties are 
        set to zero during the initialization process.  
"""
function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                            settings::PhotoExcitationFluores.Settings)
    pathways = PhotoExcitationFluores.Pathway[]
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)
                    eEnergy = nLevel.energy - iLevel.energy + settings.photonEnergyShift
                    fEnergy = nLevel.energy - fLevel.energy + settings.photonEnergyShift
                    if  eEnergy < 0.   ||   fEnergy < 0.    continue    end
                    rSettings = PhotoEmission.Settings( settings.multipoles, settings.gauges, false, false, CorePolarization(), LineSelection(), 0., 0., 0., false)
                    eChannels = PhotoEmission.determineChannels(nLevel, iLevel, rSettings) 
                    fChannels = PhotoEmission.determineChannels(fLevel, nLevel, rSettings) 
                    push!( pathways, PhotoExcitationFluores.Pathway(iLevel, nLevel, fLevel, eEnergy, fEnergy, EmProperty(0., 0.), 
                                                                    true, eChannels, fChannels) )
                end
            end
        end
    end
    return( pathways )
end


"""
`PhotoExcitationFluores.displayPhotonDm(stream::IO, pathway::PhotoExcitationFluores.Pathway, dmArray::Array{ComplexF64,2}, 
                                        settings::PhotoExcitationFluores.Settings)`  
    ... to display all requested properties for the given pathway that are related to the photon density matrix of the fluorescence
        photon. One or a few neat tables are printed but nothing is returned otherwise.
"""
function  displayPhotonDm(stream::IO, pathway::PhotoExcitationFluores.Pathway, dmArray::Array{ComplexF64,2}, settings::PhotoExcitationFluores.Settings)
    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
    msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
    fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
    #
    println(stream, "\n\n  Properties of the fluorescence photon related to its density matrix for pathway:" * 
            TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, pathway.finalLevel.index) *
                        "\n  --------------------------------------------------------------------------------")
    println(stream, " ")
    println(stream, "  + Symmetries of pathway: " * TableStrings.symmetries_imf(isym, msym, fsym))
    println(stream, "  + Excitation and fluorescence energies of pathway " * TableStrings.inUnits("energy") * ":  " *
            @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))             * ",  " *
            @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.fluorEnergy))  )
    #
    # Display the density matrix of the fluorescence photon if requested
    if  settings.calcPhotonDm
        nx = 119
        println(stream, "\n  + Photon density matrix:")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "    theta       phi               rho(1,1)               rho(1,-1)              rho(-1,1)              rho(-1,-1) ";   
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        for  i = 1:size(dmArray, 1)
            sa = "    " * @sprintf("%.4e", dmArray[i,1].re) * "  " * @sprintf("%.4e", dmArray[i,2].re) * "  " * 
                            @sprintf("%.4e", dmArray[i,3].re) * " "  * @sprintf("%.4e", dmArray[i,3].im) * "  " * 
                            @sprintf("%.4e", dmArray[i,4].re) * " "  * @sprintf("%.4e", dmArray[i,4].im) * "  " * 
                            @sprintf("%.4e", dmArray[i,5].re) * " "  * @sprintf("%.4e", dmArray[i,5].im) * "  " * 
                            @sprintf("%.4e", dmArray[i,6].re) * " "  * @sprintf("%.4e", dmArray[i,6].im)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    # Display the angular distribution of the fluorescence photon if requested
    if  settings.calcAngular
        nx = 44
        println(stream, "\n  + Angular distribution:")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "    theta       phi           W(theta, phi)   ";   println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        for  i = 1:size(dmArray, 1)
            sa = "    " * @sprintf("%.4e", dmArray[i,1].re) * "  " * @sprintf("%.4e", dmArray[i,2].re) * "    " * 
                            @sprintf("%.4e", dmArray[i,3].re + dmArray[i,6].re)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    # Display the Stokes parameters of the fluorescence photon if requested
    if  settings.calcStokes
        nx = 67
        println(stream, "\n  + Stokes parameters:")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "    theta       phi                P_1         P_2         P_3";   println(stream, sa)
        println(stream, "  ", TableStrings.hLine(nx))
        for  i = 1:size(dmArray, 1)
            p1 = -(dmArray[i,4] + dmArray[i,5]) / (dmArray[i,3] + dmArray[i,6])
            p2 =  (dmArray[i,5] - dmArray[i,4]) / (dmArray[i,3] + dmArray[i,6]) * (0.0+1.0im)
            p3 =  (dmArray[i,4] - dmArray[i,5]) / (dmArray[i,3] + dmArray[i,6])
            sa = "    " * @sprintf("%.4e", dmArray[i,1].re) * "  " * @sprintf("%.4e", dmArray[i,2].re) * "    " * 
                            @sprintf("%.4e", p1.re)           * "  " * @sprintf("%.4e", p2.re)           * "  "   *  @sprintf("%.4e", p3.re)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    
    return( nothing )
end


"""
`PhotoExcitationFluores.displayPathways(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, settings::PhotoExcitationFluores.Settings)`  
    ... to list the energies and multipoles of the selected pathways. A neat table is printed but nothing is returned otherwise.
"""
function  displayPathways(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, settings::PhotoExcitationFluores.Settings)
    nx = 110
    println(stream, " ")
    println(stream, "  Partial excitation & fluorescence (emission) amplitudes to be calculated:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  e  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  e  --  f"; na=2);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(30, "  e--i         e--f  "; na=2)
    sa = sa * TableStrings.flushleft(30, "Multipoles"; na=0);        
    sb = sb * TableStrings.flushleft(30, "i--e           e--f  "; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))      * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.fluorEnergy))      * "     "
        #
        multipoles = EmMultipole[]
        for  ech in pathway.excitChannels
            multipoles = push!( multipoles, ech.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)            * "                   "
        sa = sa * TableStrings.flushleft(14, mpString[1:14];  na=1)
        #
        multipoles = EmMultipole[]
        for  fch in pathway.fluorChannels
            multipoles = push!( multipoles, fch.multipole)
        end  
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)             * "                   "
        sa = sa * TableStrings.flushleft(14, mpString[1:14];  na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    println(stream, "\n  + Given solid angles:  $(settings.solidAngles) ")
    if settings.calcPhotonDm  println(stream, "  + Calculate the photon density matrix of the fluorescence photon at the given solid angles.")    end
    if settings.calcAngular   println(stream, "  + Calculate the angular distribution of the fluorescence photon at the given solid angles.")     end
    if settings.calcStokes    println(stream, "  + Calculate the Stokes parameters of the fluorescence photon at the given solid angles.")        end
    #
    return( nothing )
end


"""
`PhotoExcitationFluores.displayResults(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, 
                                        settings::PhotoExcitationFluores.Settings)`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, settings::PhotoExcitationFluores.Settings)
    nx = 146
    println(stream, " ")
    println(stream, "  Partial excitation & fluorescence cross sections:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    ";   sb = "    "
    sa = sa * TableStrings.center(23, "Levels"; na=2);            sb = sb * TableStrings.center(23, "i  --  e  --  f"; na=2);          
    sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * TableStrings.center(23, "i  --  e  --  f"; na=2);
    sa = sa * TableStrings.center(30, " Energies  " * TableStrings.inUnits("energy"); na=4);              
    sb = sb * TableStrings.center(30, "  e--i         e--f  "; na=2)
    sa = sa * TableStrings.flushleft(30, "Multipoles"; na=0);        
    sb = sb * TableStrings.flushleft(30, "i--e           e--f  "; na=2)
    sa = sa * TableStrings.center(30, "Cou -- cross sections -- Bab";        na=2)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "      " * 
                                        TableStrings.inUnits("cross section"); na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  pathway in pathways
        sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                        msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                        fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
        sa = sa * TableStrings.center(23, TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                            pathway.finalLevel.index); na=3)
        sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.excitEnergy))      * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pathway.fluorEnergy))      * "     "
        #
        multipoles = EmMultipole[]
        for  ech in pathway.excitChannels
            multipoles = push!( multipoles, ech.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)            * "                   "
        sa = sa * TableStrings.flushleft(14, mpString[1:14];  na=1)
        #
        multipoles = EmMultipole[]
        for  fch in pathway.fluorChannels
            multipoles = push!( multipoles, fch.multipole)
        end  
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles)             * "                   "
        sa = sa * TableStrings.flushleft(14, mpString[1:14];  na=5)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.crossSection.Coulomb))     * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", pathway.crossSection.Babushkin))   * "  "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
