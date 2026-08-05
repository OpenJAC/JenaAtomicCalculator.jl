
# Functions and methods for cascade simulations

"""
`Cascade.addLevels(levelsA::Array{Cascade.Level,1}, levelsB::Array{Cascade.Level,1})` 
    ... adds two sets of levels so that each levels occurs only 'once' in the list; in practice, however, this 'addition' requires also 
        that the parent and daughter processes are added properly so that all information is later available for the simulations.
        It is assumed here that all daughter and parent (processes) appear only once if levels from different data sets
        (Cascade.DecayData, Cascade.PhotoIonData) are added to each other. A message is issued about the number of levels before and 
        after this 'addition', and how many of the levels have been modified by this method. Note that all relative occucations are 
        set to zero in this addition; a newlevels::Array{Cascade.Level,1} is returned.
"""
function  addLevels(levelsA::Array{Cascade.Level,1}, levelsB::Array{Cascade.Level,1})
    nA = length(levelsA);   nB = length(levelsB);    nmod = 0;    nnew = 0;    newlevels = Cascade.Level[];  appendedB = falses(nB)
    
    # First all levels from levels A but taking additional parents and daughters into accout
    for  levA in levelsA
        parents   = deepcopy(levA.parents);     daughters = deepcopy(levA.daughters);   
        for  (i,levB) in enumerate(levelsB)
            if levA == levB     appendedB[i] = true
                for p in levB.parents     push!(parents,   p)   end
                for d in levB.daughters   push!(daughters, d)   end
                if  length(parents) > length(levA.parents)  ||  length(daughters) > length(levA.daughters)  nmod = nmod + 1   end
                break
            end
        end
        push!(newlevels, Cascade.Level(levA.energy, levA.J, levA.parity, levA.NoElectrons, levA.majorConfig, 0., parents, daughters) )
    end
    
    # Append those levels from levelsB that are not yet appended
    for  (i,levB) in enumerate(levelsB)
        if appendedB[i]
        else    push!(newlevels, levB);     nnew = nnew + 1
        end 
    end
    nN = length(newlevels)
    println("> Append $nnew (new) levels to $nA levels results in a total of $nN levels (with $nmod modified levels) in the list.")
    return( newlevels )
end


"""
`Cascade.assignOccupation!(levels::Array{Cascade.Level,1}, property::AbstractSimulationProperty)` 
    ... assigns the occupation due to the given property
"""
function assignOccupation!(levels::Array{Cascade.Level,1}, property::AbstractSimulationProperty)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    
    if  typeof(property) in [Cascade.IonDistribution, Cascade.PhotonIntensities]
        # Assign the relative occupation of the levels in this list due to the given initial occupation
        for  pair  in  property.initialOccupations
            levels[pair[1]].relativeOcc = pair[2]
        end
        println("> Assign an initial occupation for given level numbers.")
        if printSummary  println(iostream, "> Assign an initial occupation for given level numbers.")   end
    end

    return( nothing )
end


"""
`Cascade.combineEnergiesIntensities(w1::Float64, w1enInts::Array{Tuple{Float64,Float64},1}, 
                                    w2::Float64, w2enInts::Array{Tuple{Float64,Float64},1})` 
    ... combines w1 * w1enInts + w2 * w2enInts; a newEnergiesInts::Array{Tuple{Float64,Float64},1} is returned.
"""
function combineEnergiesIntensities(w1::Float64, w1enInts::Array{Tuple{Float64,Float64},1}, 
                                    w2::Float64, w2enInts::Array{Tuple{Float64,Float64},1})
    newEnergiesInts = Tuple{Float64,Float64}[]
    for  enInt in w1enInts   push!(newEnergiesInts, (enInt[1], w1*enInt[2]))    end
    for  enInt in w2enInts   push!(newEnergiesInts, (enInt[1], w2*enInt[2]))    end
    
    return( newEnergiesInts )
end


"""
`Cascade.displayExpansionOpacities(stream::IO, sc::String, property::Cascade.ExpansionOpacities, 
                                    energyInterval::Tuple{Float64, Float64}, kappas::Array{Basics.EmProperty,1})` 
    ... displays the expansion opacities in a neat table. Nothing is returned.
"""
function displayExpansionOpacities(stream::IO, sc::String, property::Cascade.ExpansionOpacities, 
                                    energyInterval::Tuple{Float64, Float64}, kappas::Array{Basics.EmProperty,1})
    nx = 63
    println(stream, " ")
    sa = "  Expansion opacities:  $sc       ... are evaluated for the following parameters: \n" *
        "\n    + level population                   = $(property.levelPopulation)    " *
        "\n    + opacityDependence                  = $(property.opacityDependence)    " *
        "\n    + ion density [ions/cm^3]            = $(property.ionDensity)    " *
        "\n    + plasma temperature [K]             = $(property.temperature)   " *
        "\n    + expansion/observation time [sec]   = $(property.expansionTime) " *
        "\n    + binning [Hartree]                  = $(property.opacityDependence.binning) " *
        "\n    + energy interval [Hartree]          = $(energyInterval) " *
        "\n    + energy shift  [Hartree]            = $(property.transitionEnergyShift) \n"
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    sb = TableStrings.inUnits("energy")
    if  typeof(property.opacityDependence) == Cascade.TemperatureOpacityDependence    sb = "[dim-less]"     end
    sa = "  "
    sa = sa * TableStrings.center(20, "Values " * sb; na=1)        
    sa = sa * TableStrings.center(36, "Cou -- kappa^(expansion) [cm^2/g] -- Bab";      na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for (i,value) in enumerate(property.dependencyValues)
        sa = "       " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", value)) * 
                "         " * @sprintf("%.6e", kappas[i].Coulomb) * "        " * @sprintf("%.6e", kappas[i].Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})` 
    ... displays the (current or final) ion distribution in a neat table. Nothing is returned.
"""
function displayIonDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1})
    minElectrons = 1000;   maxElectrons = 0;   totalProb = 0.;    nx = 31
    for  level in levels   minElectrons = min(minElectrons, level.NoElectrons);   maxElectrons = max(maxElectrons, level.NoElectrons)   end
    println(stream, " ")
    println(stream, "  (Final) Ion distribution for the cascade:  $sc ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "No. electrons"; na=4)        
    sa = sa * TableStrings.center(10,"Rel. occ.";      na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for n = maxElectrons:-1:minElectrons
        sa = "             " * string(n);   sa = sa[end-10:end];   prob = 0.
        for  level in levels    if  n == level.NoElectrons   prob = prob + level.relativeOcc    end    end
        sa = sa * "         " * @sprintf("%.5e", prob)
        println(stream, sa)
        totalProb = totalProb + prob
    end
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  Total distributed probability:  " * @sprintf("%.5e", totalProb)
    println(stream, sa)

    return( nothing )
end


"""
`Cascade.displayFinalLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1}, finalConfigs::Array{Configuration,1})` 
    ... displays the (current or final) level distribution in a neat table. Only those levels with a non-zero 
        occupation are displayed here. Nothing is returned.
"""
function displayFinalLevelDistribution(stream::IO, sc::String, levels::Array{Cascade.Level,1}, finalConfigs::Array{Configuration,1})
    minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels));    nx = 69
    for  i = 1:length(levels)
        minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
        energies[i]  = levels[i].energy   
    end
    enIndices = sortperm(energies, rev=true)
    # Now printout the results
    println(stream, " ")
    println(stream, "  (Final) Level distribution for the cascade:  $sc")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
    sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
    sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
    sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Rel. occ.";                                na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for n = maxElectrons:-1:minElectrons
        sa = "            " * string(n);        sa  = sa[end-10:end]
        for  en in enIndices
            saa = "            " * string(en);  saa = saa[end-12:end]
            if  n == levels[en].NoElectrons  ##    &&  levels[en].relativeOcc > 0
                sx = "    " * string( LevelSymmetry(levels[en].J, levels[en].parity) )                       * "           "
                sb = sa * saa * sx[1:15]
                sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "      "
                sb = sb * @sprintf("%.5e", levels[en].relativeOcc) 
                sa = "           "
                if       length(finalConfigs) == 0                       println(stream, sb)
                elseif   levels[en].majorConfig  in  finalConfigs        println(stream, sb)
                end
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}; extended::Bool=false)` 
    ... displays all defined levels  in a neat table, together with their No. of electrons, symmetry, level energy, 
        current (relative) population as well as analogue information about their parents and daughter levels. This 
        enables one to recognize (and perhaps later add) missing parent and daughter levels. Nothing is returned.
"""
function displayLevelTree(stream::IO, levels::Array{Cascade.Level,1}; extended::Bool=false)
    minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels));    nx = 179;    ny = 65
    for  i = 1:length(levels)
        minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
        energies[i]  = levels[i].energy   
    end
    enIndices = sortperm(energies, rev=true)
    # Now printout the results
    println(stream, " ")
    println(stream, "* Level tree of this cascade:")
    println(stream, " ")
    if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end
    sa = " "
    sa = sa * TableStrings.center( 6, "No. e-"; na=2)        
    sa = sa * TableStrings.center( 6, "Lev-No"; na=2)        
    sa = sa * TableStrings.center( 6, "J^P"          ; na=2);               
    sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(10, "Rel. occ.";                                na=3)
    if  extended
        sb = "Parents P(A/R: No_e, sym, energy) and Daughters D(A/R: No_e, sym, energy);  all energies in " * TableStrings.inUnits("energy")
        sa = sa * TableStrings.flushleft(100, sb; na=2)
    end
    # 
    println(stream, sa)
    if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end
    for n = maxElectrons:-1:minElectrons
        sa = "            " * string(n);     sa  = sa[end-5:end]
        for  en in enIndices
            saa = "         " * string(en);  saa = saa[end-8:end]
            if  n == levels[en].NoElectrons
                sx = "    " * string( LevelSymmetry(levels[en].J, levels[en].parity) )                       * "           "
                sb = sa * saa * sx[1:15]
                sb = sb * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "    "
                sb = sb * @sprintf("%.4e", levels[en].relativeOcc)                                           * "  "
                if extended
                pProcessSymmetryEnergyList = Tuple{AtomicProcess,Int64,LevelSymmetry,Float64}[]
                dProcessSymmetryEnergyList = Tuple{Basics.AtomicProcess,Int64,LevelSymmetry,Float64}[]
                for  p in levels[en].parents
                    idx = p.index
                    if      p.process == Basics.Auger()         lev = p.lineSet.linesA[idx].initialLevel
                    elseif  p.process == Basics.Radiative()     lev = p.lineSet.linesR[idx].initialLevel
                    elseif  p.process == Basics.Photo()         lev = p.lineSet.linesP[idx].initialLevel
                    else    error("stop a")    end
                    push!( pProcessSymmetryEnergyList, (p.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                end
                for  d in levels[en].daughters
                    idx = d.index
                    if      d.process == Basics.Auger()         lev = d.lineSet.linesA[idx].finalLevel
                    elseif  d.process == Basics.Radiative()     lev = d.lineSet.linesR[idx].finalLevel
                    elseif  d.process == Basics.Photo()         lev = d.lineSet.linesP[idx].finalLevel
                    else    error("stop b")    end
                    push!( dProcessSymmetryEnergyList, (d.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                end
                wa = TableStrings.processSymmetryEnergyTupels(120, pProcessSymmetryEnergyList, "P")
                if  length(wa) > 0    sc = sb * wa[1];    println(stream,  sc )    else    println(stream,  sb )   end  
                for  i = 2:length(wa)
                    sc = TableStrings.hBlank( length(sb) ) * wa[i];    println(stream,  sc )
                end
                wa = TableStrings.processSymmetryEnergyTupels(120, dProcessSymmetryEnergyList, "D")
                for  i = 1:length(wa)
                    sc = TableStrings.hBlank( length(sb) ) * wa[i];    println(stream,  sc )
                end
                else    println(stream,  sb )
                end  ## extended
                sa = "      "
            end
        end
    end
    if  extended    println(stream, "  ", TableStrings.hLine(nx))  else    println(stream, "  ", TableStrings.hLine(ny))  end

    return( nothing )
end


"""
`Cascade.displayIntensities(stream::IO, property::PhotonIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})` 
    ... displays the (tuples of) energiesIntensities in a neat table. Nothing is returned.
"""
function displayIntensities(stream::IO, property::PhotonIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})
    nx = 40
    sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.minPhotonEnergy))
    sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.maxPhotonEnergy))
    println(stream, " ")
    println(stream, "* Energies & (relative) photon intensities between " * sMinEn * " and "  * sMaxEn * 
                        TableStrings.inUnits("energy") * ":  ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Rel. Intensity"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    totalIntensity = 0.
    for  enInt in  energiesIntensities
        sa = "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", enInt[1])) * "         " * @sprintf("%.3e", enInt[2])
        totalIntensity = totalIntensity + enInt[2]
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "  Total (relative) intensity:  " * @sprintf("%.3e", totalIntensity))

    return( nothing )
end


"""
`Cascade.displayIntensities(stream::IO, property::ElectronIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})`
    ... displays the (tuples of) energiesIntensities of the emitted electrons in a neat table. Nothing is
        returned.
"""
function displayIntensities(stream::IO, property::ElectronIntensities, energiesIntensities::Array{Tuple{Float64,Float64},1})
    nx = 40
    sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.minElectronEnergy))
    sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.maxElectronEnergy))
    println(stream, " ")
    println(stream, "* Energies & (relative) electron intensities between " * sMinEn * " and "  * sMaxEn *
                        TableStrings.inUnits("energy") * ":  ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Rel. Intensity"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    totalIntensity = 0.
    for  enInt in  energiesIntensities
        sa = "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", enInt[1])) * "         " * @sprintf("%.3e", enInt[2])
        totalIntensity = totalIntensity + enInt[2]
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "  Total (relative) intensity:  " * @sprintf("%.3e", totalIntensity))

    return( nothing )
end


"""
`Cascade.displayPhotoAbsorptionSpectrum(stream::IO, pEnergies::Array{Float64,1}, crossSections::Array{EmProperty,1},
                                        property::Cascade.PhotoAbsorptionSpectrum)` 
    ... displays the photoabsorption cross sections a neat table. Nothing is returned.
"""
function displayPhotoAbsorptionSpectrum(stream::IO, pEnergies::Array{Float64,1}, crossSections::Array{EmProperty,1},
                                        property::Cascade.PhotoAbsorptionSpectrum)
    # Photon energies enter via the property
    if  length(pEnergies) != length(crossSections)  error("stop a")    end
    #
    nx = 46
    println(stream, " ")
    println(stream, "* Absorption cross sections:  ")
    println(stream, " ")
    println(stream, "  Absorption cross sections are determined for the given photon energies and for levels \n  with the" *
                    " initial population $(property.initialOccupations) \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "   * TableStrings.inUnits("energy"); na=7)
    sa = sa * TableStrings.center(10, "Total CS " * TableStrings.inUnits("cross section"); na=11)
    println(stream, sa)
    sa = "                      Coulomb       Babushkin"
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  (i, cs)  in  enumerate(crossSections)
        sa = "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", pEnergies[i]))        * "     "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.Coulomb))   * "   " 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.Babushkin)) * "         "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


#==
"""
`Cascade.displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Cascade.AbsorptionCrossSection,1}, settings::Cascade.SimulationSettings)` 
    ... displays the photoabsorption cross sections a neat table. Nothing is returned.
"""
function displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Cascade.AbsorptionCrossSection,1}, settings::Cascade.SimulationSettings)
    nx = 83
    println(stream, " ")
    println(stream, "* Absorption cross sections:  ")
    println(stream, " ")
    sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", settings.minPhotonEnergy))
    sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", settings.maxPhotonEnergy))
    println(stream, "  Absorption cross sections are determined for photon energies between " * sMinEn * " and "  *
                    sMaxEn * TableStrings.inUnits("energy") * " as well as for levels \n  with the initial population " *
                    "$(settings.initialOccupations) \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Ionization CS " * TableStrings.inUnits("cross section"); na=11)
    sa = sa * TableStrings.center(10, "Excitation CS " * TableStrings.inUnits("cross section"); na=11)
    println(stream, sa)
    sa = "                      Coulomb       Babushkin         Coulomb       Babushkin"
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  cs in  crossSections
        sa = "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", cs.photonEnergy)) * "     "
        if     cs.ionizationCS == Basics.EmProperty(0.)       sa = sa * "                                         "
        else   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.ionizationCS.Coulomb))   * "   " 
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.ionizationCS.Babushkin)) * "         "
        end
        #
        if     cs.ionizationCS == Basics.EmProperty(0.)       sa = sa * "                                         "
        else   sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.excitationCS.Coulomb))   * "   " 
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.excitationCS.Babushkin)) * "    "
        end
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Basics.ScalarProperty{EmProperty},1}, 
                                        property::Cascade.PhotoAbsorptionSpectrum)` 
    ... displays the photoabsorption cross sections a neat table. Nothing is returned.
"""
function displayPhotoAbsorptionSpectrum(stream::IO, crossSections::Array{Basics.ScalarProperty{EmProperty},1}, 
                                        property::Cascade.PhotoAbsorptionSpectrum)
    nx = 46
    println(stream, " ")
    println(stream, "* Absorption cross sections:  ")
    println(stream, " ")
    sMinEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.photonEnergies[1]))
    sMaxEn = @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", property.photonEnergies[end]))
    println(stream, "  Absorption cross sections are determined for photon energies between " * sMinEn * " and "  *
                    sMaxEn * TableStrings.inUnits("energy") * " as well as for levels \n  with the initial population " *
                    "$(property.initialOccupations) \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "   * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Total CS " * TableStrings.inUnits("cross section"); na=11)
    println(stream, sa)
    sa = "                      Coulomb       Babushkin"
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  cs in  crossSections
        sa = "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", cs.arg)) * "     "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.value.Coulomb))   * "   " 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs.value.Babushkin)) * "         "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end
==#


"""
`Cascade.displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1}, settings::Cascade.SimulationSettings)` 
    ... displays the (initial) relative occupation of the levels in a neat table; an error message is issued if the population is
        given for those levels in the settings, which do not exist in the present simulation. Nothing is returned.
"""
function displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1}, settings::Cascade.SimulationSettings)
    nx = 69
    println(stream, " ")
    println(stream, "* Initial level occupation:  ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
    sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
    sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
    sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
    println(stream, sa)
    #
    for  initialOcc in  settings.initialOccupations
        idx = initialOcc[1];   occ = initialOcc[2]
        if   idx < 1   ||   idx > length(levels)       error("In appropriate choice of initial occupation; idx = $idx")    end
        level = levels[idx]
        sa = "            " * string(level.NoElectrons);                                  sa  = sa[end-10:end]
        sb = "            " * string(idx);                                                sb  = sb[end-12:end]
        sc = "    " * string( LevelSymmetry(level.J, level.parity) )  * "           ";    sc  = sc[1:15]
        sd = sa * sb * sc  * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy))
        sd = sd * "      " * @sprintf("%.5e", level.relativeOcc) 
        println(stream, sd)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1})` 
    ... displays the (initial) relative occupation of the levels in a neat table. Nothing is returned.
"""
function displayRelativeOccupation(stream::IO, levels::Array{Cascade.Level,1})
    nx = 69
    println(stream, " ")
    println(stream, "* Initial level occupation:  ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
    sa = sa * TableStrings.center( 8, "Lev-No"; na=2)        
    sa = sa * TableStrings.center( 6, "J^P"          ; na=3);               
    sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  (levi, level) in enumerate(levels)
        sa = "            " * string(level.NoElectrons);                                  sa  = sa[end-10:end]
        sb = "            " * string(levi);                                               sb  = sb[end-12:end]
        sc = "    " * string( LevelSymmetry(level.J, level.parity) )  * "           ";    sc  = sc[1:15]
        sd = sa * sb * sc  * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", level.energy))
        sd = sd * "      " * @sprintf("%.5e", level.relativeOcc) 
        println(stream, sd)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.extractOccupation(levels::Array{Cascade.Level,1}, groundConfigs::Array{Configuration,1})` 
    ... determines the total occupation of the levels in (one of) the groundConfigs. A occ::Float64 is returned.
"""
function extractOccupation(levels::Array{Cascade.Level,1}, groundConfigs::Array{Configuration,1})
    #
    wocc = 0.
    for level in levels
        if  Basics.extractConfiguration(Basics.LeadingConfiguration(), level) in groundConfigs   wocc = wocc + level.relativeOcc     end
    end

    return( wocc )
end



"""
`Cascade.extractPhotoExcitationData(dataDicts::Array{Dict{String,Any},1})` 
    ... returns the available photoexcitation data.
"""
function extractPhotoExcitationData(dataDicts::Array{Dict{String,Any},1})
    photoexcitationData = Cascade.Data[]
    for data  in  dataDicts       results = data["results"]
        if  haskey(results, "photoexcitation lines:")
            linesE = results["photoexcitation lines:"]
            push!(photoexcitationData, Cascade.Data{PhotoExcitation.Line}(linesE))   
        end
    end
    
    return( photoexcitationData )
end



"""
`Cascade.extractPhotoIonizationData(dataDicts::Array{Dict{String,Any},1})` 
    ... returns the available photoionization data.
"""
function extractPhotoIonizationData(dataDicts::Array{Dict{String,Any},1})
    photoionizationData = Cascade.Data[]
    for data  in  dataDicts       results = data["results"]
        if  haskey(results, "photoionization lines:")  
            linesP = results["photoionization lines:"]
            push!(photoionizationData, Cascade.Data{PhotoIonization.Line}(linesP))  
        end
    end
    
    return( photoionizationData )
end


"""
`Cascade.extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)` 
    ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
        propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
        with its parent lines (which may populate the level) and the daughter lines (to which the pobability may decay). 
        A levelList::Array{Cascade.Level,1} is returned.
"""
function extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    levels = Cascade.Level[]
    print("> Extract and sort the list of levels for the given decay data ... ")
    if printSummary     print(iostream, "> Extract and sort the list of levels for the given decay data ... ")     end
        
    for cData in data
        #
        if  typeof(cData) == Cascade.Data{PhotoEmission.Line}
            linesR = cData.lines
            for  (i,line)  in  enumerate(linesR)
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[],
                                        [ Cascade.LineReference(linesR, Basics.Radiative(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesR, Basics.Radiative(), i)], 
                                        Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)  
            end
            #
        elseif  typeof(cData) == Cascade.Data{AutoIonization.Line}
            linesA = cData.lines
            for  (i,line)  in  enumerate(linesA)
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[], 
                                        [ Cascade.LineReference(linesA, Basics.Auger(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesA, Basics.Auger(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)
            end
            #
        elseif  typeof(cData) == Cascade.Data{PhotoIonization.Line}
            linesP = cData.lines
            for  (i,line)  in  enumerate(linesP)
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(linesP, Basics.Photo(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesP, Basics.Photo(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)
            end
            #
            #
        elseif  typeof(cData) == Cascade.Data{PhotoExcitation.Line}
            linesE = cData.lines
            for  (i,line)  in  enumerate(linesE)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(linesE, Basics.Radiative(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        line.finalLevel.relativeOcc, [ Cascade.LineReference(linesE, Basics.Radiative(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)  
            end
        else  error("stop a")
        end
    end
    
    # Sort the levels by energy in reversed order
    energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
    enIndices = sortperm(energies, rev=true)
    newlevels = Cascade.Level[]
    for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
    
    println("a total of $(length(newlevels)) levels were found.")
    if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
    
    return( newlevels )
end


#==
"""
`Cascade.extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)` 
    ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
        propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
        with its parent lines (which may populate the level) and the daughter lines (to which the pobability may decay). 
        A levelList::Array{Cascade.Level,1} is returned.
"""
function extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream");   found = false
    photonEnergy = settings.initialPhotonEnergy;   newlevels = Cascade.Level[]
    for dataset in data
        println(">>> photonEnergy = $(dataset.photonEnergy) ")
        if  dataset.photonEnergy == photonEnergy       found     = true
            newlevels = Cascade.extractLevels(dataset, settings);    break
        end
    end
    
    if  !found  error("No proper photo-ionizing data set (Cascade.PhotoIonData) found for the photon energy $photonEnergy ")  end
    return( newlevels )
end
==#


#==
"""
`Cascade.extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)` 
    ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
        propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
        with its parent lines (which may populate the level) and the daughter lines (to which the pobability may decay). 
        A levelList::Array{Cascade.Level,1} is returned.
"""
function extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    levels = Cascade.Level[]
    print("> Extract and sort the list of levels for the given excitation data ... ")
    if printSummary     print(iostream, "> Extract and sort the list of levels for the given excitation data ... ")     end
    
    for  i = 1:length(data.linesE)
        line = data.linesE[i]
        iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(data, Basics.Radiative(), i)] ) 
        Cascade.pushLevels!(levels, iLevel)  
        fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                line.finalLevel.relativeOcc, [ Cascade.LineReference(data, Basics.Radiative(), i)], Cascade.LineReference[] ) 
        Cascade.pushLevels!(levels, fLevel)  
    end
    
    # Sort the levels by energy in reversed order
    energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
    enIndices = sortperm(energies, rev=true)
    newlevels = Cascade.Level[]
    for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
    
    println("a total of $(length(newlevels)) levels were found.")
    if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
    
    return( newlevels )
end
==#

#==
"""
`Cascade.extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)` 
    ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
        propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
        with its parent lines (which may populate the level) and the daughter lines (to which the pobability may decay). 
        A levelList::Array{Cascade.Level,1} is returned.
"""
function extractLevels(data::Array{Cascade.Data,1}, settings::Cascade.SimulationSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    levels = Cascade.Level[]
    print("> Extract and sort the list of levels for the given photo-ionization data ... ")
    if printSummary     print(iostream, "> Extract and sort the list of levels for the given photo-ionization data ... ")     end

    for  i = 1:length(data.linesP)
        line = data.linesP[i]
        iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(data, Basics.Photo(), i)] ) 
        Cascade.pushLevels!(levels, iLevel)  
        fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                line.finalLevel.relativeOcc, [ Cascade.LineReference(data, Basics.Photo(), i)], Cascade.LineReference[] ) 
        Cascade.pushLevels!(levels, fLevel)  
    end
    
    # Sort the levels by energy in reversed order
    energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
    enIndices = sortperm(energies, rev=true)
    newlevels = Cascade.Level[]
    for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
    
    println("a total of $(length(newlevels)) levels were found.")
    if printSummary     println(iostream, "a total of $(length(newlevels)) levels were found.")     end
    
    return( newlevels )
end
==#


"""
`Cascade.findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})` 
    ... find the index of the given level within the given list of levels; an idx::Int64 is returned and an error message is 
        issued if the level is not found in the list.
"""
function findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})
    for  k = 1:length(levels)
        if  level.energy == levels[k].energy  &&   level.J == levels[k].J   &&   level.parity == levels[k].parity   &&
            level.NoElectrons == levels[k].NoElectrons
            kk = k;   return( kk )
        end
    end
    error("findLevelIndex():  No index was found;\n   level = $(level) ")
end


"""
`Cascade.interpolateIonizationCS(photonEnergy::Float64, ionizationCS::Array{Basics.ScalarProperty{EmProperty},1})` 
    ... interpolates (or extrapolates) the ionization cross sections as defined by ionizationCS for the given photonEnergy.
        If photonEnergy is outside the photon energies from ionizationCS, simply the cross section from the nearest energy
        is returned; if photonEnergy lays between two photon energies from ionizationCS, a simple linear interpolation
        rules is applied here. A cs::Basics.EmProperty is returned.
"""
function interpolateIonizationCS(photonEnergy::Float64, ionizationCS::Array{Basics.ScalarProperty{EmProperty},1})
    imin = imax = 0
    for  (i, ionCS)  in  enumerate(ionizationCS)
        if  ionCS.arg <= photonEnergy   imin = i    end
    end
    for  (i, ionCS)  in  enumerate(ionizationCS)
        if  ionCS.arg >  photonEnergy   imax = i;   break    end
    end
    #
    if       imin == 0  &&  imax == 1        return(ionizationCS[1].value)
    elseif   imax == 0                       return(ionizationCS[end].value)
    elseif   imax - imin == 1
        deltaEnergy = photonEnergy - ionizationCS[imin].arg
        totalEnergy = ionizationCS[imax].arg - ionizationCS[imin].arg 
        cs          = ionizationCS[imin].value + deltaEnergy/totalEnergy * (ionizationCS[imax].value - ionizationCS[imin].value)
        return( cs )
    else  error("stop b")    
    end
end


"""
`Cascade.perform(simulation::Cascade.Simulation`  
    ... to simulate a cascade decay (and excitation) from the given data. Different computational methods and different properties of 
        the ionic system, such as the ion distribution or final-level distribution can be derived and displayed from these simulations. 
        Of course, the details of these simulations strongly depend on the atomic processes and data that have been generated before by 
        performing a computation::Cascade.Computation. The results of all individual steps are printed to screen but nothing is 
        returned otherwise.

`Cascade.perform(simulation::Cascade.Simulation; output=true)`   
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the method and 
        specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function perform(simulation::Cascade.Simulation; output::Bool=false)
    ## The property (and, where it matters, the method) select what is simulated BY DISPATCH; this used to be a
    ## chain of ten `typeof(simulation.property) == ...` tests ending in error("stop b"). That chain is why
    ## Cascade.ElectronIntensities stayed dead for years: it had no branch, so a request for it fell straight
    ## through and returned nothing, silently. Under dispatch an unsupported combination hits the explicit
    ## fallback below and says so. Cascade.PiRateCoefficients is still in that position and now reports it.
    wa = Cascade.simulate(simulation.property, simulation.method, simulation)

    ## NOTE (05-Aug-2026): the individual branches used to merge their own keys into `results`
    ## ("energies/intensities:", "alpha^RR:", "relaxPercentage:", ...), but the closing block then RESET
    ## `results` to a fresh Dict before filling in name/property/data -- so none of those keys ever reached a
    ## caller. They were aliases for `wa` in any case, and the behaviour that callers actually saw is kept
    ## unchanged here.
    if  output
        results = Dict{String, Any}()
        results = Base.merge( results, Dict("name:"         => simulation.name) )
        results = Base.merge( results, Dict("property:"     => simulation.property) )
        results = Base.merge( results, Dict("data:"         => wa) )
        return( results )
    end

    return( nothing )
end


"""
`Cascade.simulate(property::Cascade.AbstractSimulationProperty, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`
    ... carries out the simulation that the given property asks for, by the given method. A method is defined
        for every supported (property, method) combination; this fallback catches the rest and names what is
        missing, rather than letting the request pass silently.
"""
function simulate(property::Cascade.AbstractSimulationProperty, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    error("No simulation is implemented for property $(typeof(property)) with method $(typeof(method)). " *
          "Supported: PhotoAbsorptionSpectrum, DrRateCoefficients, RrRateCoefficients, ExpansionOpacities and " *
          "RosselandOpacities with any method; IonDistribution, FinalLevelDistribution, PhotonIntensities, " *
          "ElectronIntensities and MeanRelaxationTime with Cascade.ProbPropagation().")
end


"""
`Cascade.simulate(property::Cascade.PhotoAbsorptionSpectrum, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates the photo-absorption spectrum.
"""
function simulate(property::Cascade.PhotoAbsorptionSpectrum, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    if    haskey(simulation.computationData[1]["results"], "photoionization lines:")
            linesP = simulation.computationData[1]["results"]["photoionization lines:"]
    else  linesP = PhotoIonization.Line[]
    end
    if    haskey(simulation.computationData[1]["results"], "photoexcitation lines:")
            linesE = simulation.computationData[1]["results"]["photoexcitation lines:"]
    else  linesE = PhotoExcitation.Line[]
    end
    # Display the line data if appropriate
    if  simulation.settings.printTree
        PhotoIonization.displayLineData(stdout, linesP)
        PhotoExcitation.displayLineData(stdout, linesE)
    else
        println(">>>> Set settings.printTree to list all line data explicitly.")
    end

    return( Cascade.simulatePhotoAbsorptionSpectrum(simulation, linesP, linesE) )
end


"""
`Cascade.simulate(property::Cascade.IonDistribution, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the final ion (charge-state) distribution.
"""
function simulate(property::Cascade.IonDistribution, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateIonDistribution(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.FinalLevelDistribution, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the final level distribution.
"""
function simulate(property::Cascade.FinalLevelDistribution, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateFinalLevelDistribution(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.PhotonIntensities, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the emitted photon (fluorescence) spectrum.
"""
function simulate(property::Cascade.PhotonIntensities, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulatePhotonIntensities(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.ElectronIntensities, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the emitted electron (Auger) spectrum.
"""
function simulate(property::Cascade.ElectronIntensities, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateElectronIntensities(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.MeanRelaxationTime, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the mean relaxation time of the cascade.
"""
function simulate(property::Cascade.MeanRelaxationTime, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateMeanRelaxationTime(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.DrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates dielectronic-recombination rate coefficients.
"""
function simulate(property::Cascade.DrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateDrRateCoefficients(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.RrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates radiative-recombination rate coefficients.
"""
function simulate(property::Cascade.RrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    linesR = simulation.computationData[1]["results"]["photo-recombination line data:"].linesR
    return( Cascade.simulateRrRateCoefficients(linesR, simulation) )
end


"""
`Cascade.simulate(property::Cascade.ExpansionOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates expansion opacities.
"""
function simulate(property::Cascade.ExpansionOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    photoExcData = Cascade.extractPhotoExcitationData(simulation.computationData)
    return( Cascade.simulateExpansionOpacities(photoExcData, simulation.name, simulation.property, printout=true) )
end


"""
`Cascade.simulate(property::Cascade.RosselandOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates Rosseland opacities.
"""
function simulate(property::Cascade.RosselandOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    photoExcData = Cascade.extractPhotoExcitationData(simulation.computationData)
    return( Cascade.simulateRosselandOpacities(photoExcData, simulation) )
end


"""
`Cascade.propagateOccupationInTime!(levels::Array{Cascade.Level,1}, dt::Float64)` 
    ... propagates the occupation of the levels by dt in time. 
"""
function propagateOccupationInTime!(levels::Array{Cascade.Level,1}, dt::Float64)
    #
    relativeOcc = zeros(length(levels));    relativeLoss = zeros(length(levels));    down = 0.0;    loss = 0.0
    for (i, level) in  enumerate(levels)
        # Cycle through all daughters of level and 'shift' the part of occupation down the daughter levels
        lossFactor = 0.;    occ = level.relativeOcc
        for  daughter in level.daughters
            idx = daughter.index
            if      daughter.process == Basics.Radiative()     line = daughter.lineSet.linesR[idx];  rate = line.photonRate.Coulomb
            elseif  daughter.process == Basics.Auger()         line = daughter.lineSet.linesA[idx];  rate = line.totalRate
            else    error("stop b; process = $(daughter.process) ")
            end
            downFactor = (1.0 - exp(-rate*dt));         
            newLevel   = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                        line.finalLevel.basis.NoElectrons, 0., Cascade.LineReference[], Cascade.LineReference[] )
            kk         = Cascade.findLevelIndex(newLevel, levels)
            if lossFactor + downFactor < 1.0    lossFactor = lossFactor + downFactor
            else                                downFactor = 1.0 - lossFactor;      lossFactor = 1.0
            end
            relativeOcc[kk] = relativeOcc[kk] + downFactor * occ;   down = down + downFactor * occ
        end
        levels[i].relativeOcc = levels[i].relativeOcc - lossFactor * occ;   loss = loss + lossFactor * occ
    end
    for i = 1:length(levels)  levels[i].relativeOcc = levels[i].relativeOcc + relativeOcc[i]    end
    
    return( nothing )
end


"""
`Cascade.propagateProbability!(levels::Array{Cascade.Level,1}; 
                                collectPhotonIntensities::Bool=false, collectElectronIntensities::Bool=false)` 
    ... propagates the relative level occupation through the levels of the cascade until no further change occur in the 
        relative level occupation. The argument levels is modified during the propagation, but nothing is returned otherwise.
"""
function propagateProbability!(levels::Array{Cascade.Level,1}; 
                                collectPhotonIntensities::Bool=false, collectElectronIntensities::Bool=false)
    # Inititalize arrays to be returned
    if      collectPhotonIntensities  && collectElectronIntensities      error("stop a")
    elseif  collectPhotonIntensities  || collectElectronIntensities
            energiesIntensities = Tuple{Float64,Float64}[]
    end
    
    n = 0
    println("\n*  Probability propagation through $(length(levels)) levels of the cascade:")
    while true ## n < 2  ## 
        n = n + 1;    totalProbability = 0.
        print("    $n-th round ... ")
        relativeOcc = zeros(length(levels));    relativeLoss = zeros(length(levels))
        for  level in levels
            if   level.relativeOcc > 0.   && length(level.daughters) > 0
                # A level with relative occupation > 0 has still 'daughter' levels; collect all excitation/decay rates for this level
                # Here, an excitation cross section is formally treated as a rate as it is assumed that the initial levels of
                # the photoionization process cannot decay by photon emission or autoionization processes.
                prob  = level.relativeOcc;   rates = zeros(length(level.daughters))
                for  (i,daughter) in  enumerate(level.daughters)
                    idx = daughter.index
                    if      daughter.process == Basics.Radiative()     rates[i] = daughter.lines[idx].photonRate.Coulomb
                    elseif  daughter.process == Basics.Auger()         rates[i] = daughter.lines[idx].totalRate
                    elseif  daughter.process == Basics.Photo()         rates[i] = daughter.lines[idx].crossSection.Coulomb
                    else    error("stop a; process = $(daughter.process) ")
                    end
                end
                totalRate = sum(rates)
                if      totalRate <  0.    error("stop b")
                elseif  totalRate == 0.    # do nothing
                else
                    # Shift the relative occupation to the 'daughter' levels due to the different ionization and decay pathes
                    for  (i,daughter) in  enumerate(level.daughters)
                        idx = daughter.index
                        if      daughter.process == Basics.Radiative()     line = daughter.lines[idx]
                        elseif  daughter.process == Basics.Auger()         line = daughter.lines[idx]
                        elseif  daughter.process == Basics.Photo()         line = daughter.lines[idx]
                        else    error("stop b; process = $(daughter.process) ")
                        end
                        major    = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                        newLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                                    line.finalLevel.basis.NoElectrons, major, 0., Cascade.LineReference[], Cascade.LineReference[] )
                        kk    = Cascade.findLevelIndex(newLevel, levels)
                        relativeOcc[kk] = relativeOcc[kk] + prob * rates[i] / totalRate
                        #
                        if  collectPhotonIntensities   && daughter.process == Basics.Radiative()
                            push!(energiesIntensities, (level.energy - levels[kk].energy, prob * rates[i] / totalRate))     end
                        if  collectElectronIntensities && daughter.process == Basics.Auger()
                            push!(energiesIntensities, (level.energy - levels[kk].energy, prob * rates[i] / totalRate))     end
                    end
                    level.relativeOcc = 0.;   totalProbability = totalProbability + prob
                end
            end
        end
        for i = 1:length(levels)  levels[i].relativeOcc = levels[i].relativeOcc + relativeOcc[i]    end
        println("has propagated a total of $totalProbability level occupation.")
        # Cycle once more if the relative occupation has still changed
        if  totalProbability == 0.    break    end
    end

    if  collectPhotonIntensities   ||   collectElectronIntensities   return(energiesIntensities)
    else                                                             return( nothing )
    end
end


"""
`Cascade.pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)` 
    ... push's the information of newLevel of levels. This is the standard 'push!(levels, newLevel)' if newLevel is not yet 
        including in levels, and the proper modification of the parent and daughter lines of this level otherwise. The argument 
        levels::Array{Cascade.Level,1} is modified and nothing is returned otherwise.
"""
function pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)
    for  i = 1:length(levels)
        if  newLevel.energy == levels[i].energy  &&  newLevel.J == levels[i].J  &&  newLevel.parity == levels[i].parity
            append!(levels[i].parents,   newLevel.parents)
            append!(levels[i].daughters, newLevel.daughters)
            return( nothing )
        end
    end
    push!( levels, newLevel)
    return( nothing )
end


"""
`Cascade.reviewData(simulation::Cascade.Simulation; ascendingOrder::Bool=false)` 
    ... reviews and displays the (computation) data for the given simulation; these data contains the name of the data set, 
        its initial and generated multiplets for the various blocks of (some part of the ionization and/or decay) cascade as 
        well as all the line data [lineR, linesA, lineP, ...]. From these data, this function also generates and returns
        the level tree that is to be used in the subsequent simulations, and where levels are odered in `ascending' order
        if selected.
"""
function reviewData(simulation::Cascade.Simulation; ascendingOrder::Bool=false)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    dataDicts = simulation.computationData;     settings = simulation.settings;     allLevels = Cascade.Level[]
    
    # Loop through all (computation) data set and display the major results
    for  (i,data) in  enumerate(dataDicts)
        results     = data["results"]
        multiplets  = results["initial multiplets:"]
        gMultiplets = results["generated multiplets:"]
        nlev = 0;    for multiplet in multiplets     nlev  = nlev  + length(multiplet.levels)     end
        nglev = 0;   for multiplet in gMultiplets    nglev = nglev + length(multiplet.levels)     end
        println("\n* $i) Data dictionary for cascade computation:   $(results["name"])  with  $nlev initial and  $nglev generated levels") 
        println(  "  ===========================================")
        
        Cascade.displayLevels(stdout, multiplets, sa="initial ")
        if  printSummary 
            println(iostream, "\n* $i) Data dictionary for cascade computation:   $(results["name"])  with  $nlev initial and  $nglev generated levels") 
            println(iostream,   "  ===========================================")
            Cascade.displayLevels(iostream, multiplets,  sa="initial ")
            Cascade.displayLevels(iostream, gMultiplets, sa="generated ")        
        end
        #
        if      haskey(results, "cascade data:")             lineData = results["cascade data:"]
        else    error("stop a")
        end
        
        @show typeof(lineData)
        
        #==
        if      haskey(results,"decay line data:")                          lineData = results["decay line data:"]
        elseif  haskey(results,"photo-ionizing line data:")                 lineData = results["photo-ionizing line data:"]
        elseif  haskey(results,"photo-excited line data:")                  lineData = results["photo-excited line data:"]
        elseif  haskey(results,"hollow-ion line data:")                     lineData = results["hollow-ion line data:"]
        else    error("stop a")
        end ==#
        levels = Cascade.extractLevels(lineData, settings)
        allLevels = Cascade.addLevels(allLevels, levels)
    end
    
    allLevels = Cascade.sortByEnergy(allLevels, ascendingOrder=ascendingOrder)
    Cascade.assignOccupation!(allLevels, simulation.property)
    if  simulation.settings.printTree      Cascade.displayLevelTree(stdout, allLevels, extended=false)     end
    if  simulation.settings.printLongTree  Cascade.displayLevelTree(stdout, allLevels, extended=true)      end
    
    return( allLevels )
end


"""
`Cascade.simulateDrRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... Determines and prints the DR resonance strength and (plasma) rate coefficients for all resonance levels.
        Nothing is returned.
"""
function simulateDrRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    resonances = DielectronicRecombination.CaptureLine[]
    rSelection = simulation.property.resonanceSelection
    @show length(levels)
    #
    # Collect the information about all resonances
    for  level in levels
        for daughter in level.daughters
            if  daughter.process != Basics.Auger();                            continue   end
            aLine   = daughter.lines[daughter.index]
            if  aLine.finalLevel.index != simulation.property.initialLevelNo   continue   end
            #
            dJ      = aLine.initialLevel.J;         iJ      = aLine.finalLevel.J
            dM      = aLine.initialLevel.M;         iM      = aLine.finalLevel.M
            dParity = aLine.initialLevel.parity;    iParity = aLine.finalLevel.parity
            dIndex  = aLine.initialLevel.index;     iIndex  = aLine.finalLevel.index
            dEnergy = aLine.initialLevel.energy;    iEnergy = aLine.finalLevel.energy
            captureRate = aLine.totalRate
            augerRate   = Cascade.computeTotalAugerRate(level)
            photonRate  = Cascade.computeTotalPhotonRate(level)
            strength    = Basics.EmProperty(0.)
            #
            iLevel      = ManyElectron.Level(iJ, iM, iParity, iIndex, iEnergy, 0., false, ManyElectron.Basis(), Float64[])
            dLevel      = ManyElectron.Level(dJ, dM, dParity, dIndex, dEnergy, 0., false, ManyElectron.Basis(), Float64[])
            es          = Defaults.convertUnits("energy: to atomic", simulation.property.electronEnergyShift)
            en          = dLevel.energy-iLevel.energy + es;    if  en < 0.  continue   end
            #
            wa          = Defaults.convertUnits("kinetic energy to wave number: atomic units", en)
            if   augerRate + photonRate.Babushkin == 0.
                strength = EmProperty(0.)
            elseif  DielectronicRecombination.isResonanceToBeExcluded(aLine.initialLevel, aLine.finalLevel, rSelection)
                @show DielectronicRecombination.isResonanceToBeExcluded(aLine.initialLevel, aLine.finalLevel, rSelection)
                # Set the strength to zero, if the initial (resonance) level of aLine is not selected explicitly
                strength = EmProperty(0.)
            else
                ## CORRECTED 05-Aug-2026, on two counts.
                ## (i) The factor 2 that stood here is spurious. Tu et al., Phys. Plasmas 23, 053301 (2016), Eq. (1)
                ##     give S = (g_d/2g_i) * pi^2 hbar^3/(m_e E_res) * A_r A_a/(sum A_r + sum A_a); with k^2 = 2E in
                ##     atomic units, pi^2/k^2 * g_d/g_i is identically pi^2/E * g_d/(2g_i), so the 2 is ALREADY
                ##     contained in k^2 and including it again made every Cascade DR strength twice too large. The
                ##     DielectronicRecombination module has always had it (correctly) commented out; this copy did not.
                ## (ii) Both gauges were divided by the BABUSHKIN total width, so the Coulomb strength used a
                ##     mismatched denominator. Each gauge now uses its own.
                wa       = pi*pi / (wa*wa) * captureRate *
                            ((Basics.twice(dJ) + 1) / (Basics.twice(iJ) + 1))
                sC       = augerRate + photonRate.Coulomb   == 0. ? 0. : wa * photonRate.Coulomb   / (augerRate + photonRate.Coulomb)
                sB       = augerRate + photonRate.Babushkin == 0. ? 0. : wa * photonRate.Babushkin / (augerRate + photonRate.Babushkin)
                strength = EmProperty(sC, sB)
            end
            ## The captureRate of the individual (i,m) channel is not tracked by the cascade, so it is left at 0.;
            ## the two TOTAL widths are what the strength and the printout need, and both are known here.
            newResonance = DielectronicRecombination.CaptureLine(iLevel, dLevel, en, 0., augerRate, photonRate, strength,
                                                                 AutoIonization.Channel[])
            push!(resonances, newResonance)
        end
    end
    
    # Add contributions for the high-n shell if requested; in this case, each level and its daughters are 
    # tested for having an orbital with principal quantum number nDetailed. If this is the case, scaled resonances
    # are added for all n = nDetailed+1 : nMax, ie. resonances with scaled energies and rates.
    # At present, a simple (nDetailed/n)^beta with beta = 1.1 is applied
    beta = 1.1
    if  simulation.property.nDetailed < simulation.property.nMax
        @warn("This feature nDetailed < nMax has been implemented but never tested; first check the individual n-contributions." * 
                "i.e. the nEnergy and nStrength below ... and how they contribute.")
        for  level in levels
            for daughter in level.daughters
                if  daughter.process != Basics.Auger();                            continue   end
                aLine   = daughter.lineSet.linesA[daughter.index]
                if  aLine.finalLevel.index != simulation.property.initialLevelNo   continue   end
                #
                dJ      = aLine.initialLevel.J;         iJ      = aLine.finalLevel.J
                dM      = aLine.initialLevel.M;         iM      = aLine.finalLevel.M
                dParity = aLine.initialLevel.parity;    iParity = aLine.finalLevel.parity
                dIndex  = aLine.initialLevel.index;     iIndex  = aLine.finalLevel.index
                dEnergy = aLine.initialLevel.energy;    iEnergy = aLine.finalLevel.energy
                # Determine of whether the initial level has an electron with principal quantum number nDetailed
                # scale the contribution if the basis has such a subshell for all n > nDetailed
                if  !hasSubshell(simulation.property.nDetailed, aLine.initialLevel.basis.subshells)  continue   end
                captureRate = aLine.totalRate
                augerRate   = Cascade.computeTotalAugerRate(level)
                photonRate  = Cascade.computeTotalPhotonRate(level)
                strength    = Basics.EmProperty(0.)
                #
                wa          = Defaults.convertUnits("kinetic energy to wave number: atomic units", en)
                wa          = pi*pi / (wa*wa) * captureRate  * 2 * # factor 2 is not really clear.
                                ((Basics.twice(dJ) + 1) / (Basics.twice(iJ) + 1)) /
                                (augerRate + photonRate.Babushkin)
                strength     = EmProperty(wa * photonRate.Coulomb, wa * photonRate.Babushkin)
                #
                Zeff         = sqrt( 2*simulation.property.nDetailed^2 * 
                                        Basics.extractMeanEnergy(simulation.property.nDetailed, aLine.initialLevel.basis) )
                for  n = simulation.property.nDetailed+1:simulation.property.nMax
                    nEnergy     = dEnergy - Zeff^2 / 2 * (1/n^2 - 1/simulation.property.nDetailed^2)
                    nStrength   = (simulation.property.nDetailed / n)^beta * strength
                    iLevel      = ManyElectron.Level(iJ, iM, iParity, iIndex, iEnergy, 0., false, ManyElectron.Basis(), Float64[])
                    dLevel      = ManyElectron.Level(dJ, dM, dParity, nIndex, dEnergy, 0., false, ManyElectron.Basis(), Float64[])
                    es          = Defaults.convertUnits("energy: to atomic", simulation.property.electronEnergyShift)
                    en          = dLevel.energy-iLevel.energy + es;    if  en < 0.  continue   end
                    newResonance = DielectronicRecombination.CaptureLine(iLevel, dLevel, en, 0., augerRate, photonRate,
                                                                         nFactor * nStrength, AutoIonization.Channel[])
                    push!(resonances, newResonance)
                    @show n, augerRate, photonRate, nStrength
                end
            end
        end
    end
    
    # Printout the resonance strength
    settings = DielectronicRecombination.Settings(DielectronicRecombination.Settings(), calcRateAlpha=true, temperatures=simulation.property.temperatures)
    DielectronicRecombination.displayResults(stdout, resonances, DielectronicRecombination.PhotonLine[], settings)
    DielectronicRecombination.displayRateCoefficients(stdout, resonances, settings)
    wb = DielectronicRecombination.extractRateCoefficients(resonances, settings)

    return( wb )
end


"""
`Cascade.simulateExpansionOpacities(photoexcitationData::Array{Cascade.Data,1}, name::String, 
                                    property::Cascade.ExpansionOpacities; printout::Bool=true)` 
    ... runs through all excitation lines, sums up their contributions and form a (list of) expansion opacities for the given 
        parameters. Nothing is returned.
"""
function simulateExpansionOpacities(photoexcitationData::Array{Cascade.Data,1}, name::String, 
                                    property::Cascade.ExpansionOpacities; printout::Bool=true)
    #
    function lambda_over_dlambda(opacityDependence::AbstractOpacityDependence, lineOmega::Float64, kT::Float64, depValue::Float64)
        # Calculates the value lambda / Delta lambda for the given binning and depValue, for which the opacity need to be
        # determined; the binning is assumed in Hartree (for frequency- and temperature-normalized dependence) and in
        # nm for wavelength-dependent opacities; an value = 0. is returned if the lineOmega does not fall into the (binning)
        # interval
        wa = 0.;  halfBinning = opacityDependence.binning / 2.
        #
        if     typeof(opacityDependence) == FrequencyOpacityDependence
            # Binning, depValue and lineOmega are all in Hartree and readily to compare; ratio need to be inverted, when compared
            # with wavelength.
            if  depValue - halfBinning < lineOmega < depValue + halfBinning             wa = 2* halfBinning / lineOmega        end
        elseif typeof(opacityDependence) == WavelengthOpacityDependence
            # Binning in nm, lineOmega & depValue in Hartree are first converted into nm ... and ratio is determined in nm as well
            lineOmega_nm = convertUnits("energy: from atomic to Angstrom", lineOmega) / 10.
            depValue_nm  = convertUnits("energy: from atomic to Angstrom", depValue)  / 10.
            if  depValue_nm - halfBinning < lineOmega_nm < depValue_nm + halfBinning    wa = lineOmega_nm / (2* halfBinning)   end
        elseif typeof(opacityDependence) == TemperatureOpacityDependence
            # Binning and lineOmega are in Hartree, depValue in [u] and need to be converted; 
            # ratio still need to be inverted, when compared with wavelength.
            if  depValue * kT - halfBinning < lineOmega < depValue + halfBinning * kT   wa = 2* halfBinning / lineOmega        end
        else   error("stop a")
        end
        
        return( wa )
    end
    #
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    values     = property.dependencyValues;                
    exptime    = property.expansionTime;                   
    dependence = property.opacityDependence;                   
    exptime_au = exptime / convertUnits("time: from atomic to sec", 1.0)
    T          = property.temperature;                                
    rho        = property.ionDensity
    eshift     = property.transitionEnergyShift
    ne         = 1.0 ## number density [1/a_o^3] ?? 
    NoValues   = length(values);                                      
    kappas     = Basics.EmProperty[];    for  i = 1:NoValues     push!(kappas, Basics.EmProperty(0.))   end
    #
    # Determine c in cm/s
    alpha   = Defaults.getDefaults("alpha")
    c_in_SI = Defaults.getDefaults("speed of light: c") * convertUnits("length: from atomic to fm", 1.0) * 1.0e-13 /
                convertUnits("time: from atomic to sec", 1.0)
    factor  = 1.0 / (rho * c_in_SI * exptime)
    kT      = convertUnits("temperature: from Kelvin to (Hartree) units", T)
    A_au    = convertUnits("length: from fm to atomic", 1.0e5)
    #
    if length(photoexcitationData) == 0     error("No photoexcitationData provided.")     end
    #
    minEnergy = 1000.;   maxEnergy = 0.
    for  excData  in photoexcitationData
        for  line  in excData.linesE
            if  minEnergy > line.omega     minEnergy = line.omega   end
            if  maxEnergy < line.omega     maxEnergy = line.omega   end
                omega       = line.omega + eshift
                fosc        = line.oscStrength
                g0          = Basics.twice(line.initialLevel.J) + 1;   ge = Basics.twice(line.finalLevel.J) + 1
                lambda_au   = convertUnits("energy: from atomic to Angstrom", omega) * A_au
            for  ivalue = 1:NoValues
                lmd_over_dl    = lambda_over_dlambda(dependence, omega, kT, values[ivalue])
                if  lmd_over_dl == 0.  continue  end
                tau_Cou        = pi * alpha * ne * lambda_au * exptime_au * ge / g0 * fosc.Coulomb   * exp(-omega/kT)
                tau_Bab        = pi * alpha * ne * lambda_au * exptime_au * ge / g0 * fosc.Babushkin * exp(-omega/kT)
                term_Cou       = factor * lmd_over_dl * (1.0 - exp(-tau_Cou))
                term_Bab       = factor * lmd_over_dl * (1.0 - exp(-tau_Bab))
                kappas[ivalue] = kappas[ivalue] + Basics.EmProperty(term_Cou, term_Bab)
            end
        end
    end
    #
    if  printout
        Cascade.displayExpansionOpacities(stdout, name, property, (minEnergy,maxEnergy), kappas)     
        if  printSummary   
            Cascade.displayExpansionOpacities(iostream, name, property, (minEnergy,maxEnergy), kappas)      end
            
            return( nothing )
    else  return( kappas )
    end

end


"""
`Cascade.simulateFinalLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
        propagation, it runs through all levels and propagates the probability until no level probability changes anymore. 
        Nothing is returned.
"""
function simulateFinalLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Specify and display the initial (relative) occupation
    if length(simulation.property.initialOccupations) > 0
            Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
    else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
    end
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    Cascade.propagateProbability!(levels, collectPhotonIntensities=false)  
    #
    finalDist = Cascade.displayFinalLevelDistribution(stdout, simulation.name, levels, simulation.property.finalConfigs)     
    if  printSummary   Cascade.displayFinalLevelDistribution(iostream, simulation.name, levels, simulation.property.finalConfigs)      end

    return( finalDist )
end


"""
`Cascade.simulateIonDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
        propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
        level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
"""
function simulateIonDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Specify and display the initial (relative) occupation
    if length(simulation.property.initialOccupations) > 0
            Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
    else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
    end
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    Cascade.propagateProbability!(levels, collectPhotonIntensities=false)  
    #
    ionDist = Cascade.displayIonDistribution(stdout, simulation.name, levels)     
    if  printSummary   Cascade.displayIonDistribution(iostream, simulation.name, levels)      end

    return( ionDist )
end


"""
`Cascade.simulateLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
        propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
        level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
"""
function simulateLevelDistribution(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    Cascade.propagateProbability!(levels)   
    #
    if  typeof(simulation.property) == Cascade.IonDistribution   
        Cascade.displayIonDistribution(stdout, simulation.name, levels)     
        if  printSummary   Cascade.displayIonDistribution(iostream, simulation.name, levels)      end
    end
    if  typeof(simulation.property) == Cascade.FinalLevelDistribution    
        Cascade.displayLevelDistribution(stdout, simulation.name, levels)   
        if  printSummary   Cascade.displayLevelDistribution(iostream, simulation.name, levels)    end
    end

    return( nothing )
end


"""
`Cascade.simulateMeanRelaxationTime(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... determine the mean relaxation time until 70%, 80%, 90%  of the initially occupied levels decay down to
        the ground configurations. An relaxTimes::Array{Float64,1} is returned that contains the
        mean relaxation times for 70%, 80%, 90%, ...
"""
function simulateMeanRelaxationTime(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    relaxTimes = zeros(3);     relaxPercentage = [0.7, 0.8, 0.9];
    time = 0.0;     dt = simulation.property.timeStep;   nx = 0
    # Specify and display the initial (relative) occupation
    if length(simulation.property.initialOccupations) > 0
            Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
    else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
    end
    Cascade.displayRelativeOccupation(stdout, levels)
    # Determine the smallest and largest rate for the given cascade tree
    minRate = 1.0e100;    maxRate = 0.0
    for  level in levels
        for daughter in level.daughters
            if      daughter.process == Basics.Auger()
                aLine  = daughter.lineSet.linesA[daughter.index]
                @show aLine.totalRate
                if  minRate > aLine.totalRate > 0.  minRate = aLine.totalRate   end
                if  maxRate < aLine.totalRate       maxRate = aLine.totalRate   end
            elseif  daughter.process == Basics.Radiative()
                rLine  = daughter.lineSet.linesR[daughter.index]
                @show rLine.photonRate.Coulomb
                if  minRate > rLine.photonRate.Coulomb > 0.  minRate = rLine.photonRate.Coulomb   end
                if  maxRate < rLine.photonRate.Coulomb       maxRate = rLine.photonRate.Coulomb   end
            else    error("stop a")
            end
        end
    end
    println(">> Simulate cascade with time step $(dt) for minRate = $minRate and  maxRate = $maxRate; " *
            "minRate*timeStep = $(minRate*dt) ")
    wocc = Cascade.extractOccupation(levels, simulation.property.groundConfigs)
    @show wocc
    #
    # Now propagate the occupation in time
    goon = true;    pr70 = true;    pr80 = true;    pr90 = true
    while  goon
        time = time + dt;   nx = nx + 1
        Cascade.propagateOccupationInTime!(levels, dt)
        wocc = Cascade.extractOccupation(levels, simulation.property.groundConfigs)
        if  rem(nx,10000) == 0     println(">>> $nx) time = $time [a.u.]   ground conf. occupation = $wocc")   end
        for  i = 1:length(relaxPercentage) 
            if wocc < relaxPercentage[i]    relaxTimes[i] = time    end
        end
        if  pr70  &&  wocc > 0.70   pr70 = false;   println("70%  at  time = $(relaxTimes[1])")   end
        if  pr80  &&  wocc > 0.80   pr80 = false;   println("80%  at  time = $(relaxTimes[2])")   end
        if  pr90  &&  wocc > 0.90   pr90 = false;   println("90%  at  time = $(relaxTimes[3])")   end
        if  wocc > 0.91     goon = false    end
    end

    return( relaxPercentage, relaxTimes )
end
"""
`Cascade.simulatePhotoAbsorptionSpectrum(simulation::Cascade.Simulation, 
                                            linesP::Array{PhotoIonization.Line,1}, linesE::Array{PhotoExcitation.Line,1})` 
    ... cycle through all lines and (incident photon) energies to derive the overall photo-absorption spectrum.
        The procedure interpolates the photoionization and 'adds' the photoexcitation cross sections to obtain the 
        total photoabsorption CS. A linear interpolation is used inside of the energy interval, for which photoionization 
        lines and cross sections have been calculated before. No extrapolation of cross sections is done here.
        It is also assumed that the same initial levels (indices) appear in the photoionization and photoexcitation lines.
        Moreover, energy units must be one of "eV", "Kayser", "Hartree"].
        All absorption cross sections are displayed in a neat table and are returned as lists.
"""
function simulatePhotoAbsorptionSpectrum(simulation::Cascade.Simulation, 
                                         linesP::Array{PhotoIonization.Line,1}, linesE::Array{PhotoExcitation.Line,1})
    if  !(getDefaults("unit: energy") in ["eV", "Kayser", "Hartree"])   
        error("\nFor a photo-absorption spectrum, the energy units must be one of eV, Kayser, Hartree;  units = " *
                getDefaults("unit: energy") )             
    end
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    paProperty = simulation.property;   pEnergies = Float64[];      crossSections = Basics.EmProperty[]
    for  en in paProperty.photonEnergies    push!(pEnergies, Defaults.convertUnits("energy: to atomic", en))    end
    
    # First determine and display all initial levels, which contribute to the partial or total photoabsorption cs
    initialLevels  = ManyElectron.Level[];   initialWeights = Float64[]
    initialIndices = Int64[];   for tp in paProperty.initialOccupations     push!(initialIndices, tp[1])  end
    for  occ in  paProperty.initialOccupations
        notYet = true
        for  line  in  linesP       
            if  occ[1] == line.initialLevel.index   &&   notYet   
                push!(initialLevels, line.initialLevel);   push!(initialWeights, occ[2])  
                notYet = false 
            end
        end  
        #
        for  line  in  linesE       
            if  occ[1] == line.initialLevel.index   &&   notYet   
                push!(initialLevels, line.initialLevel);   push!(initialWeights, occ[2])  
                notYet = false 
            end
        end    
    end
    #
    println(stdout, "\n  Initial levels, for which cross section data contribute to the photoabsorption cross section")
    println(stdout, "\n  ", TableStrings.hLine(55))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Weight"   ; na=3);               
    println(stdout, sa);    println(stdout, sb);    println(stdout, "  ", TableStrings.hLine(55))
    for  (i, initialLevel)  in  enumerate(initialLevels)
        sa  = "  ";    sym = LevelSymmetry( initialLevel.J, initialLevel.parity )
        sa = sa * TableStrings.center(10, TableStrings.level(initialLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", initialLevel.energy)) * "    "
        sa = sa * @sprintf("%.6e", initialWeights[i]) * "    "
        println(stdout, sa)
    end
    
    println(stdout, "\n  Number of (original) photoionization lines = $(length(linesP)) " *
                    "\n  Number of (original) photoexcitation lines = $(length(linesE)) \n ")
                    
    if  length(paProperty.shells) != 0
        # Reduce the number of linesP and linesE if partial cross sections are requested;
        # analyze the population of the (leading configuation of the) initial, final levels and exclude 
        # all levels with iocc != focc + 1 for just one of the shells
        newLinesP = PhotoIonization.Line[];    newLinesE = PhotoExcitation.Line[]
        for  line in linesP
            iConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
            fConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
            diffInsideShells = 0;   diffOutsideShells = 0;   addLine = true;   diff = 0
            for  (sh, occ) in  iConf.shells
                if  sh  in  paProperty.shells    
                    diff = occ - fConf.shells[sh]
                    if      diff == 1    diffInsideShells = diffInsideShells + 1
                    elseif  diff >  1    ||  diff < 0       addLine = false
                    end     
                else    diffOutsideShells = diffOutsideShells + abs(occ - fConf.shells[sh]) 
                end
            end
            if  addLine  &&  diffInsideShells == 1  &&   diffOutsideShells == 0
                push!(newLinesP, line)  
            end
        end
        linesP = newLinesP
    
        for  line in linesE
            iConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
            fConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
            diffInsideShells = 0;   diffOutsideShells = 0;   addLine = true;   diff = 0
            for  (sh, occ) in  iConf.shells
                if  sh  in  paProperty.shells    
                    diff = occ - fConf.shells[sh]
                    if      diff == 1    diffInsideShells = diffInsideShells + 1
                    elseif  diff >  1    ||  diff < 0       addLine = false
                    end     
                else    diffOutsideShells = diffOutsideShells + abs(occ - fConf.shells[sh]) 
                end
            end
            if  addLine  &&  diffInsideShells == 1  &&   diffOutsideShells == 0
                push!(newLinesE, line)  
            end
        end
        linesE = newLinesE
    
        println(stdout, "\n  Number of (reduced) photoionization lines = $(length(linesP)) " *
                        "\n  Number of (reduced) photoexcitation lines = $(length(linesE)) \n ")
    end
    
    # Define an empty array of proper size
    for pEnergy  in pEnergies     push!(crossSections, Basics.EmProperty(0.))    end
    
    # Collect photoionization cross section data for all photon energies and initial levels involved
    if  paProperty.includeIonization
        for  (p, pEnergy)  in  enumerate(pEnergies)
            cs = Basics.EmProperty(0.)
            for  (i, initialLevel)  in  enumerate(initialLevels)
                # The selection of individual subshells has been considered above
                cs = cs + initialWeights[i] * PhotoIonization.interpolateCrossSection(linesP, pEnergy, initialLevel)
                #== 
                if  length(paProperty.shells) != 0
                    # Add cross section data only if they refer to shells
                    error("aa: not yet implemented")
                else
                    cs = cs + initialWeights[i] * PhotoIonization.interpolateCrossSection(linesP, pEnergy, initialLevel)
                end ==#
            end 
            crossSections[p] = crossSections[p] + cs
        end
    end
    
    # Add, if requested, all photoexitation cross section data for all photon energies and initial levels involved
    gam    = Defaults.convertUnits("energy: to atomic", paProperty.resonanceWidth)
    
    if  paProperty.includeExcitation
        for  (p, pEnergy)  in  enumerate(pEnergies)
            cs = Basics.EmProperty(0.)
            for  (i, initialLevel)  in  enumerate(initialLevels)
                cs = cs + initialWeights[i] * paProperty.csScaling * 
                          PhotoExcitation.estimateCrossSection(linesE, pEnergy, gam, initialLevel)
                #==
                if  length(paProperty.shells) != 0
                    # Add cross section data only if they refer to shells
                    error("bb: not yet implemented")
                else
                    cs = cs + initialWeights[i] * paProperty.csScaling * 
                                PhotoExcitation.estimateCrossSection(linesE, pEnergy, gam, initialLevel)
                end   ==#
            end 
            crossSections[p] = crossSections[p] + cs
        end
    end        
    
    
    # Display the total or partial cross sections in tabular form
    Cascade.displayPhotoAbsorptionSpectrum(stdout, pEnergies, crossSections, paProperty)
    if  printSummary   Cascade.displayPhotoAbsorptionSpectrum(iostream, pEnergies, crossSections, paProperty)    end

    return( pEnergies, crossSections )
end    


"""
`Cascade.simulatePhotonIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
        propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
        level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
"""
function simulatePhotonIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Specify and display the initial (relative) occupation
    if length(simulation.property.initialOccupations) > 0
            Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations) 
    else    Cascade.specifyInitialOccupation!(levels, simulation.property.leadingConfigs) 
    end
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    prop = simulation.property
    energiesInts = Cascade.propagateProbability!(levels, collectPhotonIntensities=true)
    energiesInts = Cascade.truncateEnergiesIntensities(energiesInts, prop.minPhotonEnergy, prop.maxPhotonEnergy)
    #
    Cascade.displayIntensities(stdout, simulation.property, energiesInts)
    if  printSummary   Cascade.displayIntensities(iostream, simulation.property, energiesInts)      end

    return( energiesInts )
end


"""
`Cascade.simulateElectronIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)`
    ... simulates the (relative) intensities of the ELECTRONS emitted during the cascade, i.e. the Auger
        spectrum, as a function of the electron energy. This is the exact counterpart of
        Cascade.simulatePhotonIntensities: the same probability propagation is run through the cascade tree,
        but the emitted energy is collected at the Auger rather than at the radiative steps. For a K-shell
        hole in a light element the electron spectrum carries most of the decay, since the fluorescence yield
        is small, so this is usually the more informative of the two.

        The propagation machinery for this already existed -- Cascade.propagateProbability! has taken a
        collectElectronIntensities keyword all along and collects (energy, intensity) at every
        Basics.Auger() daughter -- but nothing ever called it: Cascade.perform(simulation) had no branch for
        ElectronIntensities and there was no simulate function to reach. Only the wiring was missing.

        An  energiesInts::Array{Tuple{Float64,Float64},1}  of (energy, relative intensity) is returned.
"""
function simulateElectronIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Specify and display the initial (relative) occupation
    Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations)
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    prop = simulation.property
    energiesInts = Cascade.propagateProbability!(levels, collectElectronIntensities=true)
    energiesInts = Cascade.truncateEnergiesIntensities(energiesInts, prop.minElectronEnergy, prop.maxElectronEnergy)
    #
    Cascade.displayIntensities(stdout, simulation.property, energiesInts)
    if  printSummary   Cascade.displayIntensities(iostream, simulation.property, energiesInts)      end

    return( energiesInts )
end


"""
`Cascade.simulateRosselandOpacities(photoexcitationData::Array{Cascade.Data,1}, simulation::Cascade.Simulation)` 
    ... runs through all excitation lines, sums up their contributions and form a (list of) Rosseland opacities, based on
        the expansion opacities, for the given parameters. Nothing is returned.
"""
function simulateRosselandOpacities(photoexcitationData::Array{Cascade.Data,1}, simulation::Cascade.Simulation)
    ulist, wlist = FastGaussQuadrature.gausslaguerre(8);            kappaList    = Float64[]
    
    opacityDependence = simulation.property.opacityDependence
    for  rho in simulation.property.ionDensities
        for  T in simulation.property.temperatures
            property=Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(), opacityDependence, rho, T, 
                                simulation.property.expansionTime, simulation.property.transitionEnergyShift, ulist)
                                
            kappas = Cascade.simulateExpansionOpacities(photoexcitationData, "expansion opacity for rho = $rho & T=$T", 
                                                        property, printout=false)    
            #
            # Form the u-integral of the Rosseland opacity
            rosseland = Basics.EmProperty(0.)
            for (i,u)  in  enumerate(ulist)
                rosseland = rosseland + 15.0/ (4*pi^4) * u^4 * wlist[i] / ( (1.0 - exp(-u))^2 ) * kappas[i]
            end
            #
            sa = "> Rosseland opacity for rho = " * @sprintf("%.3e",rho) * "[g/cm^3]  &  T="    * @sprintf("%.3e",T) *
                    " [K]  is  kappa^Rosseland [cm^2/g] = " * @sprintf("%.5e",rosseland.Coulomb)   * " [Coulomb]  "  *
                                                            @sprintf("%.5e",rosseland.Babushkin) * " [Babushkin]"
            println(sa)
        end
    end
    
    return( nothing )
end


"""
`Cascade.simulateRrRateCoefficients(lines::Array{PhotoRecombination.Line,1}, simulation)` 
    ... Integrates over all selected cross sections in order to determine the RR plasma rate coefficients for a Maxwellian
        distribution.
"""
function simulateRrRateCoefficients(lines::Array{PhotoRecombination.Line,1}, simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Check consistency of data
    isym = LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity);   @show isym
    ## for  line  in  lines
    ##     if  LevelSymmetry(line.initialLevel.J, line.initialLevel.parity) != isym   error("error a")    end
    ## end
    
    # Convert units into cm^3 / s
    factor  = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 * 
                Defaults.convertUnits("rate: from atomic to 1/s", 1.0) 
    
    temperatures = simulation.property.temperatures
    alphaRR      = EmProperty[]
    #
    for temp  in  temperatures
        temp_au  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
        wa       = EmProperty(0., 0.)
        #
        ## @warn "Cross sections not yet properly set."
        for  line  in  lines
            if  LevelSymmetry(line.initialLevel.J, line.initialLevel.parity) != 
                LevelSymmetry(AngularJ64(1), Basics.plus)                                        continue    end
            # Determine cross section of this line
            cs  = EmProperty(0., 0.)
            if   length(simulation.property.finalConfigurations) > 0
                # If finalConfigurations are given, only their contributions are counted and all final levels just
                # refer to these configurations
                finalConfigurations = Basics.extractConfigurations(Basics.FromBasis(), line.finalLevel.basis)
                if  !(finalConfigurations[1]  in  simulation.property.finalConfigurations)           continue    end
                if   simulation.property.finalLevelSelection.active  &&  
                    !(line.finalLevel.index  in  simulation.property.finalLevelSelection.indices)    continue    end
                @show finalConfigurations
            else  println("No configuration/level selection.")
            end 
            pcs = PhotoRecombination.computeCrossSectionForMultipoles(simulation.property.multipoles, line)
            cs  = cs + pcs
            #==
            # Determine cross section if only one final level contributes; for test purposes
            if   line.finalLevel.index == 1
                    wb = PhotoRecombination.crossSectionKramers(line.electronEnergy, 26.0::Float64, (1,50))
                    ## wb = PhotoRecombination.crossSectionStobbe(line.electronEnergy, 26.0::Float64)
                    cs = EmProperty(wb, wb)
            else  cs = EmProperty(0., 0.)
            end  ==#
            wa = wa + 2*2 / sqrt(2pi) / temp_au^(3/2) * factor * line.electronEnergy * 
                        exp(-line.electronEnergy / temp_au ) * line.weight * cs
        end
        push!(alphaRR, wa)
    end
    
    # Display the RR plasma rate coefficients
    PhotoRecombination.displayRateCoefficients(stdout, isym, simulation.property.temperatures, alphaRR)

    return( alphaRR )
end


"""
`Cascade.sortByEnergy(levels::Array{Cascade.Level,1}; ascendingOrder::Bool=false)` 
    ... sorts all levels by energy and assigns the occupation as given by the simulation
"""
function sortByEnergy(levels::Array{Cascade.Level,1}; ascendingOrder::Bool=false)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    newlevels = Cascade.Level[]
    #
    
    # Sort the levels by energy in reversed order
    energies  = zeros(length(levels));       for  i = 1:length(levels)   energies[i]  = levels[i].energy   end
    if     ascendingOrder   enIndices = sortperm(energies, rev=true)
    else                    enIndices = sortperm(energies, rev=false)
    end
    #
    newlevels = Cascade.Level[]
    for i = 1:length(enIndices)   ix = enIndices[i];    push!(newlevels, levels[ix])    end
    
    println("> Sort a total of $(length(newlevels)) levels, and to which all level numbers refer below. " *
            "Here all charged states are considered together in the overall cascade.")
    if printSummary     
        println(iostream, "> Sort a total of $(length(newlevels)) levels, and to which all level numbers refer below. " *
                            "Here all charged states are considered together in the overall cascade.")
    end

    return( newlevels )
end


"""
`Cascade.specifyInitialOccupation!(levels::Array{Cascade.Level,1}, initialOccupations::Array{Tuple{Int64,Float64},1})` 
    ... specifies the initial occupation of levels for the given relOccupation; it modifies the occupation 
        but returns nothing otherwise.
"""
function specifyInitialOccupation!(levels::Array{Cascade.Level,1}, initialOccupations::Array{Tuple{Int64,Float64},1})
    #
    for  initialOcc in  initialOccupations
        idx = initialOcc[1];   occ = initialOcc[2]
        if   idx < 1   ||   idx > length(levels)       error("In appropriate choice of initial occupation; idx = $idx")    end
        levels[idx].relativeOcc = occ
    end

    return( nothing )
end


"""
`Cascade.specifyInitialOccupation!(levels::Array{Cascade.Level,1}, leadingConfigs::Array{Configuration,1})` 
    ... specifies the initial occupation of levels for the given leadingConfigs; it modifies the occupation 
        but returns nothing otherwise.
"""
function specifyInitialOccupation!(levels::Array{Cascade.Level,1}, leadingConfigs::Array{Configuration,1})
    #
    nx = 0
    for  lev = 1:length(levels)
        if  Basics.extractConfiguration(Basics.LeadingConfiguration(), levels[lev])  in  leadingConfigs   nx = nx + 1    end
    end
    if  nx == 0     error("Inappropriate selection of leading configurations for the given set of cascade levels.")     end
    #
    # Now distribute the occupation
    for  lev = 1:length(levels)
        if  Basics.extractConfiguration(Basics.LeadingConfiguration(), levels[lev]) in leadingConfigs   levels[lev].relativeOcc = 1/nx    end
    end

    return( nothing )
end


"""
`Cascade.truncateEnergiesIntensities(energiesInts::Array{Tuple{Float64,Float64},1}, minPhotonEnergy::Float64, maxPhotonEnergy::Float64)` 
    ... reduces and truncates the energies & intensities  energiesInts; 'reduce' hereby refer to omit all intensity < 1.0e-8,
        while 'truncate' omits all energies outside of the interval [minPhotonEnergy, miaxPhotonEnergy]. An 
        newEnergiesInts::Array{Tuple{Float64,Float64},1} is returned.
"""
function truncateEnergiesIntensities(energiesInts::Array{Tuple{Float64,Float64},1}, minPhotonEnergy::Float64, maxPhotonEnergy::Float64)
    # Firt, truncate contributions to given energy range
    w1EnergiesInts = Tuple{Float64,Float64}[];   we = Float64[]
    for  enInt in energiesInts
        if  minPhotonEnergy <= enInt[1] <= maxPhotonEnergy    push!(w1EnergiesInts, enInt);     push!(we, enInt[1])   end
    end
    # Add contributions with equal energies contributions
    w2EnergiesInts = Tuple{Float64,Float64}[];  wasConsidered = falses(length(we))
    for  (en, enInt) in  enumerate(w1EnergiesInts)
        if      wasConsidered[en]   continue
        else    
            wa = findall(isequal(enInt[1]), we);    tInt = 0. 
            for  a in wa
                tInt = tInt + w1EnergiesInts[a][2];     wasConsidered[a] = true
            end
            push!(w2EnergiesInts, (enInt[1], tInt))
        end
    end
    # Finally, omit all small contributions
    newEnergiesInts = Tuple{Float64,Float64}[]
    for  enInt in w2EnergiesInts
        if enInt[2] > 1.0e-8    push!(newEnergiesInts, enInt)   end
    end
    
    return( newEnergiesInts )
end
