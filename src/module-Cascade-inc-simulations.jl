
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
                pProcessSymmetryEnergyList = Tuple{Basics.AtomicProcess,Int64,LevelSymmetry,Float64}[]
                dProcessSymmetryEnergyList = Tuple{Basics.AtomicProcess,Int64,LevelSymmetry,Float64}[]
                for  p in levels[en].parents
                    ## A Cascade.LineReference already carries the line list of its OWN process, so no
                    ## per-process selection is needed; the former p.lineSet.linesA/linesR/linesP referred to
                    ## an earlier layout in which one reference pointed at a SET of lists.
                    idx = p.index
                    lev = p.lines[idx].initialLevel
                    push!( pProcessSymmetryEnergyList, (p.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                end
                for  d in levels[en].daughters
                    idx = d.index
                    lev = d.lines[idx].finalLevel
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
    println(stream, "* Energies & photon yields between " * sMinEn * " and "  * sMaxEn *
                        TableStrings.inUnits("energy") * ":  ")
    println(stream, "  The yield is the number of photons emitted PER INITIAL ION at that energy: an absolute,")
    println(stream, "  dimensionless quantity, not a relative one and not a field intensity in W/cm^2.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(13, "Number / ion"; na=2)
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
    println(stream, "  Total number emitted per ion:  " * @sprintf("%.3e", totalIntensity))

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
    println(stream, "* Energies & electron yields between " * sMinEn * " and "  * sMaxEn *
                        TableStrings.inUnits("energy") * ":  ")
    println(stream, "  The yield is the number of electrons emitted PER INITIAL ION at that energy: an absolute,")
    println(stream, "  dimensionless quantity, not a relative one.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(13, "Number / ion"; na=2)
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
    println(stream, "  Total number emitted per ion:  " * @sprintf("%.3e", totalIntensity))

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
          "ElectronIntensities and RelaxationCurve with Cascade.ProbPropagation().")
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
`Cascade.simulate(property::Cascade.RelaxationCurve, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates the mean relaxation time of the cascade.
"""
function simulate(property::Cascade.RelaxationCurve, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateRelaxationCurve(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.ParticleCoincidences, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)`   ... simulates a coincidence spectrum.
"""
function simulate(property::Cascade.ParticleCoincidences, method::Cascade.ProbPropagation,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateParticleCoincidences(levels, simulation) )
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
`Cascade.simulate(property::Cascade.EaCrossSections, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`
    ... simulates the excitation-autoionization contribution to the electron-impact ionization cross section from
        the ImpactExcitation.Line's and AutoIonization.Line's of a previous ElectronIonizationScheme computation.
        Only those excited levels that carry at least one Auger line contribute, i.e. exactly the autoionizing
        ones.  An Array{Basics.ScalarProperty{Float64},1} of (impact energy, sigma^EA) is returned.
"""
function simulate(property::Cascade.EaCrossSections, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    linesE = simulation.computationData[1]["results"]["impact-excitation lines:"]
    linesA = simulation.computationData[1]["results"]["autoionization lines:"]
    if  length(linesE) == 0     error("Cascade.EaCrossSections: the cascade data carry no impact-excitation lines.")   end
    if  length(linesA) == 0
        error("Cascade.EaCrossSections: the cascade data carry no autoionization lines, so NO excited level of " *
              "this cascade autoionizes and the EA cross section is zero by construction.  Widen " *
              "ElectronIonizationScheme.excitationToShells until levels above the ionization threshold are reached.")
    end
    ## The autoionizing levels are those that appear as the INITIAL level of an Auger line.  Their branching ratio
    ## is taken as 1: this scheme computes no radiative rates, so the result is an upper bound (see the docstring).
    autoIonizing = unique([ (l.initialLevel.index, l.initialLevel.energy)  for l in linesA ])
    #
    energies = length(property.electronEnergies) > 0 ? property.electronEnergies :
                      sort(unique([Defaults.convertUnits("energy: from atomic", l.initialElectronEnergy) for l in linesE]))
    results  = Basics.ScalarProperty{Float64}[]
    for  en  in  energies
        en_au = Defaults.convertUnits("energy: to atomic", en)
        cs    = 0.
        for  l  in  linesE
            if  abs(l.initialElectronEnergy - en_au) / max(en_au, 1.0e-10) > 1.0e-6      continue    end
            if  !( (l.finalLevel.index, l.finalLevel.energy)  in  autoIonizing )          continue    end
            cs = cs + l.crossSection
        end
        push!( results, Basics.ScalarProperty(en_au, cs) )
    end
    #
    Cascade.displayEaCrossSections(stdout, results, property)
    if  printSummary   Cascade.displayEaCrossSections(iostream, results, property)    end
    #
    return( results )
end


"""
`Cascade.displayEaCrossSections(stream::IO, results::Array{Basics.ScalarProperty{Float64},1},
                                property::Cascade.EaCrossSections)`
    ... displays the excitation-autoionization cross sections in a neat table; nothing is returned.
"""
function displayEaCrossSections(stream::IO, results::Array{Basics.ScalarProperty{Float64},1},
                                property::Cascade.EaCrossSections)
    nx = 64
    println(stream, " ")
    println(stream, "  Excitation-autoionization (EA) contribution to the electron-impact ionization:")
    println(stream, " ")
    println(stream, "    Summed over the excited levels that autoionize, with a branching ratio of 1;")
    println(stream, "    this is an UPPER BOUND, since no radiative decay of those levels is computed here.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "      Impact energy         sigma^EA        ")
    println(stream, "      " * TableStrings.inUnits("energy") * "               " *
                    TableStrings.inUnits("cross section"))
    println(stream, "  ", TableStrings.hLine(nx))
    for  r  in  results
        sa = "      " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", r.arg)) * "        " *
                        @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", r.value))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Cascade.simulate(property::Cascade.EieRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`
    ... simulates electron-impact excitation rate coefficients and effective collision strengths from the
        ImpactExcitation.Line's of a previous cascade computation.  The lines carry a cross section and a
        collision strength for each incident energy; this method interpolates the collision strengths over
        energy and integrates them against a Maxwellian at the requested temperatures.  An
        Array{ImpactExcitation.RateCoefficients,1} is returned.
"""
function simulate(property::Cascade.EieRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    linesE = simulation.computationData[1]["results"]["impact-excitation lines:"]
    if  length(linesE) < 3
        error("Cascade.EieRateCoefficients needs the collision strengths at three or more electron energies; " *
              "the given cascade data carry only $(length(linesE)) line(s).  Widen ImpactExcitationScheme.electronEnergies.")
    end
    ## The aggregation itself lives in ImpactExcitation; only the temperatures are taken from the property.
    ## numElectronEnergies must match the number of DISTINCT incident energies actually present in the lines:
    ## ImpactExcitation.groupLines reshapes the line list into (energies x transitions) and refuses anything
    ## that does not divide evenly.  The Settings default of 6 has nothing to do with what the cascade computed,
    ## so it is taken from the data here.
    nEnergies = length(unique([l.initialElectronEnergy  for l in linesE]))
    settings  = ImpactExcitation.Settings(ImpactExcitation.Settings(), calcRateCoefficient=true,
                                          temperatures=property.temperatures, numElectronEnergies=nEnergies)
    rates    = ImpactExcitation.computeRateCoefficients(linesE, settings)
    ImpactExcitation.displayResults(stdout, rates)
    if  printSummary    ImpactExcitation.displayResults(iostream, rates)    end
    #
    return( rates )
end


"""
`Cascade.simulate(property::Cascade.RrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates radiative-recombination rate coefficients.
"""
function simulate(property::Cascade.RrRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    ## Cascade.Data carries its lines in the field `lines`; `.linesR` has not existed since that struct was
    ## generalised, so this path raised a FieldError on every call and the RR rate coefficients could never
    ## be simulated at all.
    linesR = simulation.computationData[1]["results"]["photo-recombination line data:"].lines
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
    ## Exact for one time step, at ANY dt: the probability that a level has decayed at all after dt is
    ## 1 - exp(-totalRate*dt), and that amount is shared among its daughters in proportion to their individual
    ## rates. The previous formulation applied 1 - exp(-rate*dt) to each daughter separately and capped the
    ## running sum at 1, which is only valid while every rate*dt << 1: for a larger step the first daughter
    ## in the list absorbed the whole population and the remaining branches received nothing, so the branching
    ## ratios were silently destroyed. That did not matter while the caller used a tiny fixed step, but it is
    ## fatal on a logarithmic time grid, where dt deliberately becomes large compared with the fast rates.
    relativeOcc = zeros(length(levels))
    for (i, level) in  enumerate(levels)
        occ = level.relativeOcc
        if  occ <= 0.  ||  length(level.daughters) == 0    continue    end
        # Total decay rate out of this level, and the rate of each individual branch
        rates = zeros(length(level.daughters))
        for  (k, daughter) in  enumerate(level.daughters)
            line = daughter.lines[daughter.index]
            if      daughter.process == Basics.Radiative()     rates[k] = line.photonRate.Coulomb
            elseif  daughter.process == Basics.Auger()         rates[k] = line.totalRate
            else    error("stop b; process = $(daughter.process) ")
            end
        end
        totalRate = sum(rates)
        if  totalRate <= 0.    continue    end
        pLeave = (1.0 - exp(-totalRate*dt)) * occ
        #
        for  (k, daughter) in  enumerate(level.daughters)
            line       = daughter.lines[daughter.index]
            major      = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
            newLevel   = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity,
                                        line.finalLevel.basis.NoElectrons, major, 0., Cascade.LineReference[], Cascade.LineReference[] )
            kk         = Cascade.findLevelIndex(newLevel, levels)
            relativeOcc[kk] = relativeOcc[kk] + pLeave * rates[k] / totalRate
        end
        levels[i].relativeOcc = levels[i].relativeOcc - pLeave
    end
    for i = 1:length(levels)  levels[i].relativeOcc = levels[i].relativeOcc + relativeOcc[i]    end

    return( nothing )
end


"""
`Cascade.propagateProbability!(levels::Array{Cascade.Level,1})`
    ... propagates the relative level occupation through the levels of the cascade until no further change
        occurs. The argument levels is modified during the propagation.

        A  fluxes::Array{Cascade.LineFlux,1}  is returned: one record per line that probability has flowed
        through, carrying the process, the energy of the emitted particle, the flux itself, and the indices of
        the initial and final level. Observables are NOT collected here -- each simulation property derives
        what it needs afterwards, see Cascade.extractIntensities. Until 05-Aug-2026 this function took
        collectPhotonIntensities/collectElectronIntensities keywords, could honour only one at a time
        (error("stop a")), and so had to be re-run in full for a second observable; a photon-electron
        coincidence could not be expressed at all.
"""
function propagateProbability!(levels::Array{Cascade.Level,1})
    fluxes = Cascade.LineFlux[]

    n = 0
    println("\n*  Probability propagation through $(length(levels)) levels of the cascade:")
    while true
        n = n + 1;    totalProbability = 0.
        print("    $n-th round ... ")
        relativeOcc = zeros(length(levels));    relativeLoss = zeros(length(levels))
        for  (li, level) in  enumerate(levels)
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
                        flux  = prob * rates[i] / totalRate
                        relativeOcc[kk] = relativeOcc[kk] + flux
                        ## Record the flux through this line, whatever its process; the properties sort it out later.
                        push!(fluxes, Cascade.LineFlux(daughter.process, level.energy - levels[kk].energy, flux, li, kk))
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

    return( fluxes )
end


"""
`Cascade.extractIntensities(fluxes::Array{Cascade.LineFlux,1}, process::Basics.AbstractProcess)`
    ... extracts, from the flux records of a cascade propagation, the (energy, intensity) tuples of all
        emissions belonging to the given process -- Basics.Radiative() for the photon spectrum,
        Basics.Auger() for the electron spectrum. An  Array{Tuple{Float64,Float64},1}  is returned.
"""
function extractIntensities(fluxes::Array{Cascade.LineFlux,1}, process::Basics.AbstractProcess)
    energiesIntensities = Tuple{Float64,Float64}[]
    for  lf in fluxes
        if  lf.process == process    push!(energiesIntensities, (lf.energy, lf.flux))    end
    end

    return( energiesIntensities )
end


"""
`Cascade.matchesGate(gate::Cascade.ElectronGate, emission::Tuple{Basics.AbstractProcess,Float64})`
    ... true if the emission is an electron within the gate's energy window. A Bool is returned.
"""
function matchesGate(gate::Cascade.ElectronGate, emission::Tuple{Basics.AbstractProcess,Float64})
    return( emission[1] == Basics.Auger()  &&  gate.minEnergy <= emission[2] <= gate.maxEnergy )
end


"""
`Cascade.matchesGate(gate::Cascade.PhotonGate, emission::Tuple{Basics.AbstractProcess,Float64})`
    ... true if the emission is a photon within the gate's energy window. A Bool is returned.
"""
function matchesGate(gate::Cascade.PhotonGate, emission::Tuple{Basics.AbstractProcess,Float64})
    return( emission[1] == Basics.Radiative()  &&  gate.minEnergy <= emission[2] <= gate.maxEnergy )
end


"""
`Cascade.walkPathways!(levels::Array{Cascade.Level,1}, index::Int64, prob::Float64,
                       emissions::Array{Tuple{Basics.AbstractProcess,Float64},1}, cutoff::Float64, pathways::Array)`
    ... walks the cascade tree from the level with the given index, follows every branch while accumulating the
        probability along it, and records one entry per completed PATHWAY: the ordered list of emissions
        (process and energy of each emitted particle), the index of the level the pathway ends in, and the
        probability of the pathway as a whole. Nothing is returned; `pathways` is appended to.

        This is deliberately not what Cascade.propagateProbability! does. The propagation yields the total flux
        through each line, which is all a singles spectrum needs, but it forgets by which route the probability
        arrived. A coincidence is a statement about routes -- the JOINT probability that one particle was
        emitted AND then another -- so pathways have to be kept intact. Cascade.Level's `daughters` links are
        what make this possible at all.

        Branches whose probability falls below `cutoff` are terminated: pathway counts multiply with cascade
        depth, and a deep cascade is not otherwise enumerable.
"""
function walkPathways!(levels::Array{Cascade.Level,1}, index::Int64, prob::Float64,
                       emissions::Array{Tuple{Basics.AbstractProcess,Float64},1}, cutoff::Float64, pathways::Array)
    level = levels[index]
    if  prob < cutoff  ||  length(level.daughters) == 0
        push!(pathways, (emissions, index, prob));    return( nothing )
    end
    #
    rates = zeros(length(level.daughters))
    for  (i,daughter) in  enumerate(level.daughters)
        idx = daughter.index
        if      daughter.process == Basics.Radiative()     rates[i] = daughter.lines[idx].photonRate.Coulomb
        elseif  daughter.process == Basics.Auger()         rates[i] = daughter.lines[idx].totalRate
        elseif  daughter.process == Basics.Photo()         rates[i] = daughter.lines[idx].crossSection.Coulomb
        else    error("stop a; process = $(daughter.process) ")
        end
    end
    totalRate = sum(rates)
    if  totalRate <= 0.
        push!(pathways, (emissions, index, prob));    return( nothing )
    end
    #
    for  (i,daughter) in  enumerate(level.daughters)
        idx      = daughter.index;    line = daughter.lines[idx]
        major    = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
        newLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity,
                                  line.finalLevel.basis.NoElectrons, major, 0., Cascade.LineReference[], Cascade.LineReference[] )
        kk       = Cascade.findLevelIndex(newLevel, levels)
        branch   = prob * rates[i] / totalRate
        if  branch < cutoff    continue    end
        newEmissions = copy(emissions);    push!(newEmissions, (daughter.process, level.energy - levels[kk].energy))
        Cascade.walkPathways!(levels, kk, branch, newEmissions, cutoff, pathways)
    end

    return( nothing )
end


"""
`Cascade.simulateParticleCoincidences(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)`
    ... simulates a coincidence measurement: the spectrum of the observed particle, restricted to those decay
        pathways in which every gate is satisfied. Gates are matched in emission order, and the observed
        particle is the emission FOLLOWING the last gate; a Cascade.ChargeStateGate among the gates is instead
        applied to the charge state that the pathway ends in. All energies are in ATOMIC UNITS, as everywhere
        else in these properties. An  Array{Tuple{Float64,Float64},1}  of (energy, intensity) is returned.
"""
function simulateParticleCoincidences(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    prop = simulation.property
    if  typeof(prop.observed) == Cascade.ChargeStateGate
        error("Cascade.ParticleCoincidences: a ChargeStateGate cannot be the OBSERVED quantity -- that is an " *
              "ion distribution, not a spectrum. Use Cascade.IonDistribution instead.")
    end
    Cascade.specifyInitialOccupation!(levels, prop.initialOccupations)
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    cutoff   = 1.0e-10
    pathways = Any[]
    for  (li, level) in  enumerate(levels)
        if  level.relativeOcc > 0.
            Cascade.walkPathways!(levels, li, level.relativeOcc,
                                  Tuple{Basics.AbstractProcess,Float64}[], cutoff, pathways)
        end
    end
    println("\n*  Coincidence analysis: $(length(pathways)) decay pathways enumerated (probability cutoff $cutoff).")
    #
    seqGates    = filter(g -> typeof(g) != Cascade.ChargeStateGate, prop.gates)
    chargeGates = filter(g -> typeof(g) == Cascade.ChargeStateGate, prop.gates)
    #
    energiesInts = Tuple{Float64,Float64}[];   nAccepted = 0
    for  (emissions, finalIndex, pathProb) in pathways
        ok = true
        for  g in chargeGates
            if  levels[finalIndex].NoElectrons != g.NoElectrons    ok = false    end
        end
        if  !ok    continue    end
        ## The gates are matched as an ordered SUBSEQUENCE of the emissions, not against adjacent positions.
        ## That is what a coincidence experiment actually measures: between two detected Auger electrons the
        ## ion may well have emitted a photon, and a detector counting electrons is indifferent to it. Strict
        ## adjacency would discard exactly those pathways. With no gates at all, every emission is therefore a
        ## candidate, and the result reduces to the ungated singles spectrum -- which is the check that caught
        ## this: adjacency gave only the FIRST emission and so recovered just 0.95 of the 1.91 electrons that
        ## this cascade emits per ion.
        pos = 0
        for  g in seqGates
            found = 0
            for  j = pos+1:length(emissions)
                if  Cascade.matchesGate(g, emissions[j])    found = j;   break    end
            end
            if  found == 0    ok = false;   break    end
            pos = found
        end
        if  !ok    continue    end
        ## every later emission that matches the observation window contributes
        for  j = pos+1:length(emissions)
            if  Cascade.matchesGate(prop.observed, emissions[j])
                push!(energiesInts, (emissions[j][2], pathProb));    nAccepted = nAccepted + 1
            end
        end
    end
    println("*  $nAccepted pathways satisfy all gates and contribute to the coincidence spectrum.")
    #
    energiesInts = Cascade.truncateEnergiesIntensities(energiesInts, -1.0e10, 1.0e10)
    Cascade.displayCoincidences(stdout, prop, energiesInts)
    if  printSummary   Cascade.displayCoincidences(iostream, prop, energiesInts)      end

    return( energiesInts )
end


"""
`Cascade.displayCoincidences(stream::IO, property::Cascade.ParticleCoincidences,
                             energiesIntensities::Array{Tuple{Float64,Float64},1})`
    ... displays the coincidence spectrum in a neat table. Nothing is returned.
"""
function displayCoincidences(stream::IO, property::Cascade.ParticleCoincidences,
                             energiesIntensities::Array{Tuple{Float64,Float64},1})
    nx = 40
    println(stream, " ")
    println(stream, "* Coincidence spectrum of the observed $(typeof(property.observed)), " *
                    "with $(length(property.gates)) gate(s):  ")
    for  (k,g) in enumerate(property.gates)   println(stream, "    gate $k:    $g")   end
    println(stream, "    observed:  $(property.observed)")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(16, "Energy "        * TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(13, "Number / ion"; na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    totalIntensity = 0.
    for  enInt in  energiesIntensities
        sa = "     " * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", enInt[1])) *
                       "         " * @sprintf("%.3e", enInt[2])
        totalIntensity = totalIntensity + enInt[2]
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "  Total coincidence counts per ion:  " * @sprintf("%.3e", totalIntensity))

    return( nothing )
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
                aLine   = daughter.lines[daughter.index]
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
    ## alpha^DR is the observable this whole scheme exists for, yet it used to go to the screen only, while the
    ## intermediate Auger and radiative tables were written to the summary file.  Send it there as well.
    if  printSummary
        DielectronicRecombination.displayResults(iostream, resonances, DielectronicRecombination.PhotonLine[], settings)
        DielectronicRecombination.displayRateCoefficients(iostream, resonances, settings)
    end
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
    Cascade.propagateProbability!(levels)  
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
    Cascade.propagateProbability!(levels)  
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
`Cascade.simulateRelaxationCurve(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)` 
    ... determine the mean relaxation time until 70%, 80%, 90%  of the initially occupied levels decay down to
        the ground configurations. An relaxTimes::Array{Float64,1} is returned that contains the
        mean relaxation times for 70%, 80%, 90%, ...
"""
function simulateRelaxationCurve(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    relaxPercentage = [0.7, 0.8, 0.9];    relaxTimes = zeros(length(relaxPercentage))
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
                aLine  = daughter.lines[daughter.index]
                if  minRate > aLine.totalRate > 0.  minRate = aLine.totalRate   end
                if  maxRate < aLine.totalRate       maxRate = aLine.totalRate   end
            elseif  daughter.process == Basics.Radiative()
                rLine  = daughter.lines[daughter.index]
                if  minRate > rLine.photonRate.Coulomb > 0.  minRate = rLine.photonRate.Coulomb   end
                if  maxRate < rLine.photonRate.Coulomb       maxRate = rLine.photonRate.Coulomb   end
            else    error("stop a")
            end
        end
    end

    ## The time grid is LOGARITHMIC, not uniform. A cascade spans an enormous range of rates -- for the Mg
    ## K-hole tree of example-Fb.jl the fastest channel is 4.9e-2 a.u. and the slowest 4.4e-12 a.u., a
    ## stiffness ratio of 1e10 -- because a single metastable level decays only through a strongly forbidden
    ## channel. Reaching the last percent of the population therefore takes ~1e11 a.u., which a fixed step of
    ## 1e-3 a.u. would need 1e14 iterations to cover. Stepping geometrically from just below the fastest
    ## lifetime to well beyond the slowest needs a few thousand points instead, and resolves every decade
    ## equally well -- which is what a relaxation curve spanning ten decades actually calls for.
    tFirst        = 0.01 / maxRate
    tLast         = 20.0 / minRate
    pointsPerDec  = 40
    nSteps        = max(2, round(Int, pointsPerDec * log10(tLast/tFirst)))
    println(">> Relaxation curve on a logarithmic time grid: minRate = $minRate, maxRate = $maxRate a.u., " *
            "i.e. lifetimes from $(round(1/maxRate, digits=2)) to $(round(1/minRate, sigdigits=3)) a.u.")
    println(">> $nSteps steps from t = $(round(tFirst, sigdigits=3)) to $(round(tLast, sigdigits=3)) a.u. " *
            "($pointsPerDec per decade); property.timeStep is not used by this grid.")

    times = Float64[];   occupations = Float64[]
    tPrev = 0.0
    for  n = 0:nSteps
        t  = tFirst * 10.0^(n * log10(tLast/tFirst) / nSteps)
        Cascade.propagateOccupationInTime!(levels, t - tPrev);    tPrev = t
        wocc = Cascade.extractOccupation(levels, simulation.property.groundConfigs)
        push!(times, t);    push!(occupations, wocc)
        for  i = 1:length(relaxPercentage)
            if  relaxTimes[i] == 0.  &&  wocc >= relaxPercentage[i]    relaxTimes[i] = t    end
        end
    end
    #
    Cascade.displayRelaxationCurve(stdout, relaxPercentage, relaxTimes, times, occupations)
    if  printSummary   Cascade.displayRelaxationCurve(iostream, relaxPercentage, relaxTimes, times, occupations)   end

    return( relaxPercentage, relaxTimes )
end


"""
`Cascade.displayRelaxationCurve(stream::IO, relaxPercentage::Array{Float64,1}, relaxTimes::Array{Float64,1},
                                times::Array{Float64,1}, occupations::Array{Float64,1})`
    ... displays the relaxation curve, i.e. the fraction of the initial population that has reached the ground
        configurations as a function of time, together with the times at which the requested percentiles are
        crossed. Nothing is returned.
"""
function displayRelaxationCurve(stream::IO, relaxPercentage::Array{Float64,1}, relaxTimes::Array{Float64,1},
                                times::Array{Float64,1}, occupations::Array{Float64,1})
    auToFs = 2.4188843265e-2      # 1 a.u. of time in femtoseconds
    println(stream, " ")
    println(stream, "* Relaxation curve: population that has reached the ground configuration(s) ")
    println(stream, " ")
    println(stream, "  ------------------------------------------------------------")
    println(stream, "      time [a.u.]        time [fs]         rel. occupation     ")
    println(stream, "  ------------------------------------------------------------")
    for  (i,t) in  enumerate(times)
        ## print a readable subset: every fourth grid point, plus the ends
        if  i == 1  ||  i == length(times)  ||  rem(i,4) == 0
            println(stream, "     " * @sprintf("%.4e", t) * "      " * @sprintf("%.4e", t*auToFs) *
                            "        " * @sprintf("%.6f", occupations[i]))
        end
    end
    println(stream, "  ------------------------------------------------------------")
    println(stream, " ")
    for  (i,p) in  enumerate(relaxPercentage)
        if  relaxTimes[i] > 0.
            println(stream, "  $(round(Int, 100p))% relaxed at t = " * @sprintf("%.4e", relaxTimes[i]) *
                            " a.u. = " * @sprintf("%.4e", relaxTimes[i]*auToFs) * " fs")
        else
            println(stream, "  $(round(Int, 100p))% relaxed:  NOT reached within the simulated time range.")
        end
    end

    return( nothing )
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
    fluxes       = Cascade.propagateProbability!(levels)
    energiesInts = Cascade.extractIntensities(fluxes, Basics.Radiative())
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

        The propagation machinery for this already existed -- Cascade.propagateProbability! had taken a
        collectElectronIntensities keyword all along, and collected (energy, intensity) at every
        Basics.Auger() daughter -- but nothing ever called it: Cascade.perform(simulation) had no branch for
        ElectronIntensities and there was no simulate function to reach. Only the wiring was missing. Since
        05-Aug-2026 the propagation records the flux through every line instead, and the Auger emissions are
        selected from those records by Cascade.extractIntensities.

        An  energiesInts::Array{Tuple{Float64,Float64},1}  of (energy, relative intensity) is returned.
"""
function simulateElectronIntensities(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    # Specify and display the initial (relative) occupation
    Cascade.specifyInitialOccupation!(levels, simulation.property.initialOccupations)
    Cascade.displayRelativeOccupation(stdout, levels)
    #
    prop = simulation.property
    fluxes       = Cascade.propagateProbability!(levels)
    energiesInts = Cascade.extractIntensities(fluxes, Basics.Auger())
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
    isym = LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity)
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
            ## Select the initial symmetry from the DATA, not from a hard-wired J^P = 1+.  The 1+ was a
            ## leftover from one particular test case and silently discarded every line of any ion whose
            ## initial level has a different symmetry -- e.g. every closed-shell (0+) ground state, for which
            ## the rate coefficient then came out as exactly zero.  isym is what the table header announces.
            if  LevelSymmetry(line.initialLevel.J, line.initialLevel.parity) != isym               continue    end
            # Determine cross section of this line
            cs  = EmProperty(0., 0.)
            if   length(simulation.property.finalConfigurations) > 0
                # If finalConfigurations are given, only their contributions are counted and all final levels just
                # refer to these configurations
                finalConfigurations = Basics.extractConfigurations(Basics.FromBasis(), line.finalLevel.basis)
                if  !(finalConfigurations[1]  in  simulation.property.finalConfigurations)           continue    end
                if   simulation.property.finalLevelSelection.active  &&  
                    !(line.finalLevel.index  in  simulation.property.finalLevelSelection.indices)    continue    end
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
    ## ... and to the summary file, for the same reason as for alpha^DR above.
    if  printSummary
        PhotoRecombination.displayRateCoefficients(iostream, isym, simulation.property.temperatures, alphaRR)
    end

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
