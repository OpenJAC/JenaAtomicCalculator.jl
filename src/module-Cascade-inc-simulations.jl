
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
`Cascade.binningUnit(dependence::Cascade.AbstractOpacityDependence)`  
    ... returns the unit in which the binning of the given opacity dependence is specified; a String is returned.
        The binning is in nm for a wavelength dependence and in Hartree otherwise, a mixture that is easily
        mistaken and that was previously printed as [Hartree] in all three cases.
"""
function binningUnit(dependence::Cascade.AbstractOpacityDependence)
    if      typeof(dependence) == Cascade.WavelengthOpacityDependence     return( "[nm]" )
    else                                                                  return( "[Hartree]" )
    end
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
        "\n    + ion number density [1/cm^3]        = $(property.ionNumberDensity)    " *
        "\n    + mass density [g/cm^3]              = $(property.massDensity)    " *
        "\n    + plasma temperature [K]             = $(property.temperature)   " *
        "\n    + expansion/observation time [sec]   = $(property.expansionTime) " *
        "\n    + binning                            = $(property.opacityDependence.binning) $(Cascade.binningUnit(property.opacityDependence)) " *
        "\n    + energy interval [Hartree]          = $(energyInterval) " *
        "\n    + energy shift  [Hartree]            = $(property.transitionEnergyShift) \n"
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    ## The bin centres are handed over as photon energies in all three cases, but a wavelength-dependent
    ## opacity is read against wavelength -- printing it in eV made the standard plot of the literature
    ## impossible to compare with. The binning likewise was labelled [Hartree] whatever the dependence.
    isWavelength = typeof(property.opacityDependence) == Cascade.WavelengthOpacityDependence
    sb = TableStrings.inUnits("energy")
    if  typeof(property.opacityDependence) == Cascade.TemperatureOpacityDependence    sb = "[dim-less]"     end
    if  isWavelength                                                                  sb = "[Angstrom]"    end
    sa = "  "
    sa = sa * TableStrings.center(20, "Values " * sb; na=1)        
    sa = sa * TableStrings.center(36, "Cou -- kappa^(expansion) [cm^2/g] -- Bab";      na=2)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for (i,value) in enumerate(property.dependencyValues)
        wx = isWavelength ?  convertUnits("energy: from atomic to Angstrom", value)  :
                             Defaults.convertUnits("energy: from atomic", value)
        sa = "       " * @sprintf("%.6e", wx) * 
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
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(linesE, Basics.PhotoExc(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesE, Basics.PhotoExc(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)  
            end
        elseif  typeof(cData) == Cascade.Data{PhotoRecombination.Line}
            linesR = cData.lines
            for  (i,line)  in  enumerate(linesR)
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(linesR, Basics.Rec(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesR, Basics.Rec(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)  
            end
        elseif  typeof(cData) == Cascade.Data{ImpactExcitation.Line}
            linesI = cData.lines
            for  (i,line)  in  enumerate(linesI)
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
                iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                        major, line.initialLevel.relativeOcc, Cascade.LineReference[], [ Cascade.LineReference(linesI, Basics.ImpactExc(), i)] ) 
                Cascade.pushLevels!(levels, iLevel)  
                major  = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
                fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                        major, line.finalLevel.relativeOcc, [ Cascade.LineReference(linesI, Basics.ImpactExc(), i)], Cascade.LineReference[] ) 
                Cascade.pushLevels!(levels, fLevel)  
            end
        else
            error("Cascade.extractLevels(): no branch for $(typeof(cData)).  A Cascade.Data{T} is built for " *
                  "PhotoEmission, AutoIonization, PhotoIonization, PhotoExcitation, PhotoRecombination and " *
                  "ImpactExcitation lines; any other T has to be given its own branch here before a simulation " *
                  "can use it.  Silently dropping the lines would leave the level list incomplete.")
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

`Cascade.perform(simulation::Cascade.Simulation; output::Bool=false)`   
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
          "Supported: PhotoAbsorptionSpectrum, DrRateCoefficients, RrRateCoefficients, EaCrossSections, " *
          "ResonantIonizationStrengths, ExpansionOpacities and " *
          "MeanOpacities with any method; IonDistribution, FinalLevelDistribution, PhotonIntensities, " *
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
`Cascade.simulate(property::Cascade.ResonantIonizationStrengths, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates the resonant contribution to electron-impact ionization.
"""
function simulate(property::Cascade.ResonantIonizationStrengths, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateResonantIonizationStrengths(levels, simulation) )
end


"""
`Cascade.simulate(property::Cascade.PiRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`
    ... simulates the photoionization rate per ion by folding the cascade's cross sections with the given photon
        fields.  An Array{Basics.EmProperty,1} in 1/s is returned.
"""
function simulate(property::Cascade.PiRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    return( Cascade.simulatePiRateCoefficients(simulation) )
end


"""
`Cascade.simulate(property::Cascade.EiiRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`
    ... simulates the electron-impact ionization plasma rate coefficients, summed over whichever of the resonant and
        excitation-autoionization channels the cascade data contain.  An Array{Basics.EmProperty,1} in cm^3/s is
        returned.
"""
function simulate(property::Cascade.EiiRateCoefficients, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    levels = Cascade.reviewData(simulation, ascendingOrder=true)
    return( Cascade.simulateEiiRateCoefficients(levels, simulation) )
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
`Cascade.simulate(property::Cascade.MeanOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)`   ... simulates mean (Rosseland or Planck) opacities.
"""
function simulate(property::Cascade.MeanOpacities, method::Cascade.AbstractSimulationMethod,
                  simulation::Cascade.Simulation)
    photoExcData = Cascade.extractPhotoExcitationData(simulation.computationData)
    return( Cascade.simulateMeanOpacities(photoExcData, simulation) )
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
`Cascade.findCascadeLevel(levels::Array{Cascade.Level,1}, mLevel::ManyElectron.Level, NoElectrons::Int64)`
    ... finds, among `levels`, the cascade level that corresponds to the many-electron level `mLevel` of an ion with
        `NoElectrons` electrons; a level::Union{Cascade.Level,Nothing} is returned, and nothing if there is no match.

        The match is on the electron number together with the energy, J and parity, since a Cascade.Level carries no
        reference back to the ManyElectron.Level it was built from.  The energy tolerance is relative and loose
        (1e-10), because both sides come from the same computation and differ only by round-trip through the data file.
"""
function findCascadeLevel(levels::Array{Cascade.Level,1}, mLevel::ManyElectron.Level, NoElectrons::Int64)
    for  level in levels
        if  level.NoElectrons == NoElectrons   &&   level.J == mLevel.J   &&   level.parity == mLevel.parity   &&
            abs(level.energy - mLevel.energy) <= 1.0e-10 * max(abs(mLevel.energy), 1.0)
            return( level )
        end
    end

    return( nothing )
end


"""
`Cascade.extractIonizingResonances(levels::Array{Cascade.Level,1}, property::Cascade.AbstractSimulationProperty)`
    ... determines and prints the resonance strengths of the two resonant electron-capture channels of electron-impact
        ionization, together with the recombination strength of the SAME resonances so that the competition between the
        two is visible.  Nothing is returned.

        A resonance is any level with one electron MORE than the initial one that has an Auger line back to the initial
        level; that Auger rate is the capture rate, by detailed balance.  An intermediate of the sequential route is
        admitted when it has the initial electron number AND carries at least one Auger line of its own, i.e. when the
        computation itself found it to be autoionizing.
"""
function extractIonizingResonances(levels::Array{Cascade.Level,1}, property::Cascade.AbstractSimulationProperty;
                                   strict::Bool=true)
    resonances = Cascade.IonizingResonance[]
    if  length(levels) == 0    return( resonances )    end
    es    = Defaults.convertUnits("energy: to atomic", property.electronEnergyShift)
    #
    # IDENTIFY THE INITIAL LEVEL BY ELECTRON NUMBER AND ENERGY, NOT BY ITS INDEX ALONE.  A ManyElectron.Level's index
    # counts within ITS OWN multiplet, so level 1 of the excited block carries the same index as level 1 of the initial
    # block; testing the index alone made every ordinary Auger to an excited level look like a capture, and produced a
    # set of spurious "resonances" at the wrong energies entirely.  The resonances are the levels with the LARGEST
    # electron number, and the initial ion is one electron below them; its ground level is the lowest of those.
    nMax  = maximum( level.NoElectrons  for level in levels )
    nIni  = nMax - 1
    iniLs = filter(level -> level.NoElectrons == nIni, levels)
    if  length(iniLs) == 0
        if  !strict    return( resonances )    end
        error("Cascade.extractIonizingResonances(): no level with $nIni electrons was found, so the ion the " *
              "capture starts from is not in the data.  This property needs a computation that requested one of the " *
              "resonant channels of Cascade.ElectronIonizationScheme.")
    end
    iniLevel = iniLs[ argmin([level.energy for level in iniLs]) ]
    #
    for  level in levels
        if  level.NoElectrons != nMax                                          continue   end
        # find the Auger line back to the initial GROUND level; its rate is the capture rate, by detailed balance
        captureRate = 0.;   twoJi = 0
        for  daughter in level.daughters
            if  daughter.process != Basics.Auger()                             continue   end
            aLine = daughter.lines[daughter.index]
            if  aLine.finalLevel.J      != iniLevel.J                          continue   end
            if  aLine.finalLevel.parity != iniLevel.parity                     continue   end
            if  abs(aLine.finalLevel.energy - iniLevel.energy) > 1.0e-8 * max(abs(iniLevel.energy), 1.0)  continue  end
            captureRate = captureRate + aLine.totalRate;    twoJi = Basics.twice(aLine.finalLevel.J)
        end
        if  captureRate == 0.                                                  continue   end
        en = level.energy - iniLevel.energy + es
        if  en <= 0.                                                           continue   end
        #
        augerD  = Cascade.computeTotalAugerRate(level)
        photonD = Cascade.computeTotalPhotonRate(level)
        gammaDb = augerD + photonD.Babushkin;      gammaDc = augerD + photonD.Coulomb
        if  gammaDb <= 0.                                                      continue   end
        twoJd   = Basics.twice(level.J)
        #
        # the SEQUENTIAL route: every Auger daughter that lands on a level which itself autoionizes.  Each gauge is
        # carried through with its OWN total width, since the widths sit in the denominators of the branching ratios.
        sSeqB = 0.;   sSeqC = 0.
        for  daughter in level.daughters
            if  daughter.process != Basics.Auger()                             continue   end
            aLine = daughter.lines[daughter.index]
            nLevel = Cascade.findCascadeLevel(levels, aLine.finalLevel, level.NoElectrons - 1)
            if  nLevel === iniLevel                                            continue   end
            if  nLevel === nothing                                             continue   end
            augerN  = Cascade.computeTotalAugerRate(nLevel)
            photonN = Cascade.computeTotalPhotonRate(nLevel)
            gammaNb = augerN + photonN.Babushkin;  gammaNc = augerN + photonN.Coulomb
            if  augerN <= 0.  ||  gammaNb <= 0.                                continue   end
            sSeqB = sSeqB + ResonantImpactIonization.sequentialStrength(en, captureRate, twoJd, twoJi,
                                                                        aLine.totalRate, gammaDb, augerN, gammaNb)
            sSeqC = sSeqC + ResonantImpactIonization.sequentialStrength(en, captureRate, twoJd, twoJi,
                                                                        aLine.totalRate, gammaDc, augerN, gammaNc)
        end
        # the SIMULTANEOUS route, only if a double-Auger probability was supplied
        if  property.dblAugerProbability <= 0.
            sSim = Basics.EmProperty(0.)
        else
            sSim = Basics.EmProperty(ResonantImpactIonization.simultaneousStrength(en, captureRate, twoJd, twoJi,
                                            property.dblAugerProbability * augerD, gammaDc),
                                     ResonantImpactIonization.simultaneousStrength(en, captureRate, twoJd, twoJi,
                                            property.dblAugerProbability * augerD, gammaDb))
        end
        # the RECOMBINATION strength of the same resonance, for comparison
        s0   = ResonantImpactIonization.resonanceStrength(en, captureRate, twoJd, twoJi)
        sDr  = Basics.EmProperty(s0 * photonD.Coulomb / gammaDc, s0 * photonD.Babushkin / gammaDb)
        #
        push!( resonances, Cascade.IonizingResonance(level, en, captureRate, augerD, photonD,
                                                     Basics.EmProperty(sSeqC, sSeqB), sSim, sDr) )
    end

    return( resonances )
end


"""
`Cascade.simulateResonantIonizationStrengths(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)`
    ... displays the resonance strengths of the resonant electron-capture channels of electron-impact ionization,
        beside the dielectronic-recombination strength of the same resonances.  Nothing is returned.
"""
function simulateResonantIonizationStrengths(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    property   = simulation.property
    resonances = Cascade.extractIonizingResonances(levels, property)
    #
    nRes   = 0
    sumSeq = 0.;   sumSim = 0.;   sumDr = Basics.EmProperty(0.)
    ## Every line goes to BOTH streams: the summary file is what the test suite compares against, so a table that
    ## reached only stdout could not be regression-tested at all.
    sayBoth = function(line::String)
        println(line);   if  printSummary   println(iostream, line)   end
    end
    sayBoth("\n  Resonance strengths of the resonant electron-capture channels of electron-impact ionization:")
    sayBoth("\n    resonance         E(res) [eV]     S(sequential)     S(simultaneous)    S(recombination)   sum of branchings")
    sayBoth("  " * "-"^116)
    for  r  in  resonances
        sumSeq = sumSeq + r.sequential.Babushkin;   sumSim = sumSim + r.simultaneous.Babushkin
        sumDr  = sumDr + r.recombination;           nRes   = nRes + 1
        gammaD = r.augerRate + r.photonRate.Babushkin
        sa = "    " * TableStrings.level(nRes) * "  " * string(LevelSymmetry(r.level.J, r.level.parity)) * "  "
        sayBoth(sa * "    " * @sprintf("%10.4f", Defaults.convertUnits("energy: from atomic to eV", r.electronEnergy)) *
                "     " * @sprintf("%.6e", r.sequential.Babushkin) * "      " * @sprintf("%.6e", r.simultaneous.Babushkin) *
                "      " * @sprintf("%.6e", r.recombination.Babushkin) * "        " *
                @sprintf("%.8f", (r.augerRate + r.photonRate.Babushkin)/gammaD))
    end
    sayBoth("  " * "-"^116)
    sayBoth("    TOTAL over $nRes resonances        " * @sprintf("%.6e", sumSeq) * "      " * @sprintf("%.6e", sumSim) *
            "      " * @sprintf("%.6e", sumDr.Babushkin))
    sayBoth("\n    Strengths are energy-integrated, in atomic units.  The last column is the sum of ALL branchings of the")
    sayBoth("    resonance and must be 1 to machine precision: it checks the arithmetic, NOT that every decay route was")
    sayBoth("    generated.  S(recombination) is the dielectronic-recombination strength of the SAME resonances, on the")
    sayBoth("    same footing, so that the competition between recombination and ionization can be read off directly.")
    if  property.dblAugerProbability <= 0.
        sayBoth("    S(simultaneous) is identically zero because dblAugerProbability was left at 0.; see the property's")
        sayBoth("    docstring on why the cascade does not choose that number for you.")
    end

    return( nothing )
end


"""
`Cascade.simulateEiiRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)`
    ... forms the electron-impact ionization plasma rate coefficient alpha^EII (T) by adding whichever of its channels
        the given cascade data contain, and reports the breakdown rather than only the sum.  An
        Array{Basics.EmProperty,1} of length length(property.temperatures) is returned, in cm^3/s.
"""
function simulateEiiRateCoefficients(levels::Array{Cascade.Level,1}, simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    property = simulation.property
    temps    = property.temperatures
    if  length(temps) == 0
        error("Cascade.EiiRateCoefficients: no temperature was given, so there is nothing to form.  alpha^EII is a " *
              "function of temperature; set the `temperatures` field of the property.")
    end
    #
    # ---- the RESONANT half ---------------------------------------------------------------------------------------
    # The Boltzmann factor and the conversion to cm^3/s are taken from DielectronicRecombination.computeRateCoefficient
    # itself, rather than written out again here, so that the IONIZATION and the RECOMBINATION rate coefficients of the
    # same resonances cannot drift apart.  Only the STRENGTH differs between the two: it is the double-autoionization
    # strength here and the radiative one there.  The dummy levels below are inert -- that function reads only the
    # resonance energy and the strength.
    ## strict=false: a cascade that computed ONLY the impact-excitation channel carries no resonance at all, and for
    ## this property that is a legitimate input rather than an error -- it sums whichever channels are present.  The
    ## alternative, wrapping the call in a bare `catch` that substitutes an empty list, would also have swallowed a
    ## genuine failure and reported "no resonances" for it.
    resonances = Cascade.extractIonizingResonances(levels, property, strict=false)
    dummy    = ManyElectron.Level(AngularJ64(0), AngularM64(0), Basics.plus, 0, 0., 0., false, ManyElectron.Basis(), Float64[])
    alphaRes = Basics.EmProperty[];    alphaDr = Basics.EmProperty[]
    for  temp  in  temps
        wa = Basics.EmProperty(0.);    wb = Basics.EmProperty(0.)
        for  r  in  resonances
            sIon  = r.sequential + r.simultaneous
            cIon  = DielectronicRecombination.CaptureLine(dummy, dummy, r.electronEnergy, 0., r.augerRate, r.photonRate,
                                                          sIon, AutoIonization.PartialWave[])
            cRec  = DielectronicRecombination.CaptureLine(dummy, dummy, r.electronEnergy, 0., r.augerRate, r.photonRate,
                                                          r.recombination, AutoIonization.PartialWave[])
            wa = wa + DielectronicRecombination.computeRateCoefficient(cIon, temp)
            wb = wb + DielectronicRecombination.computeRateCoefficient(cRec, temp)
        end
        push!(alphaRes, wa);    push!(alphaDr, wb)
    end
    #
    # ---- the EXCITATION-AUTOIONIZATION half ----------------------------------------------------------------------
    csEnergies, csValues = Cascade.extractEaCrossSections(simulation)
    alphaEa, weightAbove = Cascade.foldWithMaxwellian(csEnergies, csValues, temps)
    #
    # ---- the DIRECT half, semi-empirically ---------------------------------------------------------------------
    if  property.directCharge > 0.
        alphaDir, skipped = Cascade.directIonizationAlpha(property.directCharge, property.directConfig, temps)
    else
        alphaDir = zeros(length(temps));    skipped = String[]
    end
    #
    total = Basics.EmProperty[ alphaRes[i] + Basics.EmProperty(alphaEa[i] + alphaDir[i])  for i = 1:length(temps) ]
    #
    Cascade.displayEiiRateCoefficients(stdout, temps, alphaRes, alphaEa, alphaDir, skipped, alphaDr, total,
                                       weightAbove, resonances, csEnergies, property)
    if  printSummary
        Cascade.displayEiiRateCoefficients(iostream, temps, alphaRes, alphaEa, alphaDir, skipped, alphaDr, total,
                                           weightAbove, resonances, csEnergies, property)
    end

    return( total )
end


"""
`Cascade.displayEiiRateCoefficients(stream::IO, temperatures::Array{Float64,1}, alphaRes::Array{Basics.EmProperty,1},
                                    alphaEa::Array{Float64,1}, alphaDr::Array{Basics.EmProperty,1},
                                    total::Array{Basics.EmProperty,1}, weightAbove::Array{Float64,1},
                                    resonances::Array{Cascade.IonizingResonance,1}, csEnergies::Array{Float64,1},
                                    property::Cascade.EiiRateCoefficients)`
    ... displays the electron-impact ionization rate coefficients channel by channel, together with the two diagnostics
        that say whether the numbers may be used: the position of the maximum of alpha^res(T), which is fixed
        analytically by the resonance energies, and the share of the Maxwellian weight lying above the last computed
        impact energy.  Nothing is returned.
"""
function displayEiiRateCoefficients(stream::IO, temperatures::Array{Float64,1}, alphaRes::Array{Basics.EmProperty,1},
                                    alphaEa::Array{Float64,1}, alphaDir::Array{Float64,1}, skipped::Array{String,1},
                                    alphaDr::Array{Basics.EmProperty,1},
                                    total::Array{Basics.EmProperty,1}, weightAbove::Array{Float64,1},
                                    resonances::Array{Cascade.IonizingResonance,1}, csEnergies::Array{Float64,1},
                                    property::Cascade.EiiRateCoefficients)
    nx = 118
    println(stream, " ")
    println(stream, "  Electron-impact ionization plasma rate coefficients  alpha^EII (T),  Babushkin gauge:")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "        T [K]        kT [eV]        alpha(resonant)     alpha(exc-autoion)     alpha(direct)" *
                    "        alpha(TOTAL)   ")
    println(stream, "                                       [cm^3/s]             [cm^3/s]             [cm^3/s]  " *
                    "          [cm^3/s]     ")
    println(stream, "  ", "-"^nx)
    for  i = 1:length(temperatures)
        kT = Defaults.convertUnits("energy: from atomic to eV",
                                   Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temperatures[i]))
        sd = property.directCharge > 0. ? @sprintf("%.6e", alphaDir[i]) : "  not avail. "
        sa = "     " * @sprintf("%.4e", temperatures[i]) * "    " * @sprintf("%10.3f", kT) *
             "       " * @sprintf("%.6e", alphaRes[i].Babushkin) * "        " * @sprintf("%.6e", alphaEa[i]) *
             "        " * sd * "      " * @sprintf("%.6e", total[i].Babushkin)
        println(stream, sa)
    end
    println(stream, "  ", "-"^nx)
    println(stream, " ")
    if  property.directCharge > 0.
        println(stream, "    alpha(direct) IS A SEMI-EMPIRICAL ESTIMATE AND NOT A CASCADE RESULT.  It is the Lotz (1967)")
        println(stream, "    cross section, summed over the occupied subshells of $(property.directConfig) at Z = " *
                        @sprintf("%.1f", property.directCharge) * " and folded with the same Maxwellian, using")
        println(stream, "    tabulated binding energies.  Expect tens of per cent, worse near threshold and for")
        println(stream, "    near-neutral ions -- it is NOT of the same quality as the two computed channels beside it.")
        println(stream, "    It is included because omitting it is worse: for a light ion the direct channel carries 98%")
        println(stream, "    or more of the rate, so a total without it is wrong by a factor of order 100, and a number")
        println(stream, "    wrong by 100 misleads where one wrong by 30% does not.")
        if  length(skipped) > 0
            println(stream, "    SUBSHELLS OMITTED from the direct sum; the reason each gave, verbatim:")
            for  sa  in  skipped     println(stream, "       " * sa)     end
            println(stream, "    Their contribution is missing from alpha(direct), which is therefore a lower bound.")
        end
    else
        println(stream, "    alpha(direct) is NOT AVAILABLE and is shown as absent rather than as zero.  There is no")
        println(stream, "    Cascade.perform for ImpactIonizationScheme, so no cascade produces the direct lines; for a")
        println(stream, "    neutral or near-neutral target the direct channel is normally the LARGEST of the three, and")
        println(stream, "    the TOTAL above is therefore a lower bound on the ionization rate, not the ionization rate.")
        println(stream, "    Set directCharge and directConfig to add a semi-empirical Lotz estimate of it.")
    end
    println(stream, " ")
    #
    # ---- the DR comparison, free of charge: it is the other fate of the very same resonances -------------------
    if  length(resonances) > 0
        println(stream, "    For comparison, the DIELECTRONIC RECOMBINATION rate coefficient of the SAME resonances,")
        println(stream, "    i.e. the competing fate of each capture, on the same footing:")
        println(stream, " ")
        println(stream, "          T [K]          alpha^DR [cm^3/s]      alpha^res(ion) / alpha^DR")
        for  i = 1:length(temperatures)
            ratio = alphaDr[i].Babushkin == 0. ? 0. : alphaRes[i].Babushkin / alphaDr[i].Babushkin
            println(stream, "       " * @sprintf("%.4e", temperatures[i]) * "        " *
                            @sprintf("%.6e", alphaDr[i].Babushkin) * "            " * @sprintf("%12.4f", ratio))
        end
        println(stream, " ")
    end
    #
    # ---- diagnostic 1: WHERE THE MAXIMUM MUST LIE.  alpha^res(T) ~ T^(-3/2) exp(-E/T) has d(ln alpha)/dT = 0 at
    # kT = 2E/3 exactly.  BE CLEAR ABOUT WHAT THIS TESTS.  It tests the FOLD -- the exponent, the sign of the
    # exponential and the Kelvin-to-Hartree conversion -- and it is falsifiable there: a T^(-1/2) in place of
    # T^(-3/2) moves the true maximum to kT = 2E and the bracket test fails at once.  It does NOT test that the
    # resonances were correctly identified, because both sides of the comparison are built from the SAME energies:
    # had the level-index collision of 22-Aug-2026 put the resonances at 24 eV instead of 320, E would have been
    # 24 eV, the prediction 16 eV, and the curve would have peaked obediently at 16 eV.  What guards against THAT
    # is the line below printing E in eV beside the temperature grid, where a reader who knows the ionization
    # threshold of the ion will see 24 eV and stop.  That is a reporting virtue, not an automatic check.
    if  length(resonances) > 0
        sTot = sum(r.sequential.Babushkin + r.simultaneous.Babushkin  for r in resonances)
        eBar = sTot > 0. ? sum((r.sequential.Babushkin + r.simultaneous.Babushkin) * r.electronEnergy
                                for r in resonances) / sTot :
                           sum(r.electronEnergy for r in resonances) / length(resonances)
        kTpk = Defaults.convertUnits("energy: from atomic to eV", 2*eBar/3)
        eBev = Defaults.convertUnits("energy: from atomic to eV", eBar)
        kTs  = [Defaults.convertUnits("energy: from atomic to eV",
                    Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", t))  for t in temperatures]
        imax = argmax([a.Babushkin for a in alphaRes])
        println(stream, "    SHAPE CHECK on alpha^res(T).  With alpha ~ T^(-3/2) E exp(-E/T) the maximum sits at kT = 2E/3")
        println(stream, "    exactly.  This checks the FOLD, not the resonances: it fails if the exponent, the sign of the")
        println(stream, "    exponential or the temperature conversion is wrong, but NOT if the resonances sit at the wrong")
        println(stream, "    energies, since prediction and curve are built from the same ones.  Read E below against the")
        println(stream, "    ionization threshold you expect -- that comparison is yours to make, not the code's.")
        println(stream, "      strength-weighted mean resonance energy   E     = " * @sprintf("%10.3f", eBev) * " eV")
        println(stream, "      predicted maximum at                      kT    = " * @sprintf("%10.3f", kTpk) * " eV")
        println(stream, "      largest tabulated value falls at          kT    = " * @sprintf("%10.3f", kTs[imax]) * " eV")
        ## The test is that the tabulated maximum is one of the two grid points BRACKETING the predicted one, and
        ## not that it lies numerically close to it.  alpha(T) is unimodal, so it is monotone on either side of the
        ## true maximum; no grid point can then exceed the bracketing point on its own side, and the discrete
        ## argmax must be one of the two.  That statement is exact and holds on ANY grid, whereas comparing
        ## |kT_max - kT_predicted| against a tolerance merely measures how coarse the grid is: the first version of
        ## this check used a 50% tolerance and reported INCONSISTENT for a perfectly sound four-point grid whose
        ## neighbouring points sat at 86 and 431 eV around a predicted 204 eV.
        if      kTpk < kTs[1]
            println(stream, "      -> the predicted maximum lies BELOW the whole temperature grid; widen it downwards.")
            if  imax != 1                     println(stream, "         INCONSISTENT: the tabulated maximum is not at the lowest temperature.")   end
        elseif  kTpk > kTs[end]
            println(stream, "      -> the predicted maximum lies ABOVE the whole temperature grid; widen it upwards.")
            if  imax != length(kTs)           println(stream, "         INCONSISTENT: the tabulated maximum is not at the highest temperature.")   end
        else
            lo = findlast(k -> k <= kTpk, kTs);     hi = findfirst(k -> k >= kTpk, kTs)
            println(stream, "      bracketing grid points                          " * @sprintf("%10.3f", kTs[lo]) *
                            " and " * @sprintf("%10.3f", kTs[hi]) * " eV")
            if  imax == lo  ||  imax == hi
                println(stream, "      -> consistent: the tabulated maximum is one of the two points bracketing kT = 2E/3.")
            else
                println(stream, "      -> INCONSISTENT: the tabulated maximum lies outside the bracket, which alpha(T) being")
                println(stream, "         unimodal forbids.  Check the resonance energies before using these numbers; the")
                println(stream, "         strengths may belong to the wrong levels.")
            end
        end
        println(stream, " ")
    end
    #
    # ---- diagnostic 2: how much of the Maxwellian the computed impact energies actually cover -------------------
    if  length(csEnergies) > 0
        eMax = Defaults.convertUnits("energy: from atomic to eV", csEnergies[end])
        println(stream, "    TRUNCATION of the excitation-autoionization integral.  sigma^EA was computed at " *
                        "$(length(csEnergies)) impact")
        println(stream, "    energies, the largest being " * @sprintf("%.3f", eMax) * " eV, and is taken as zero above it.  " *
                        "The share of the Maxwellian")
        println(stream, "    weight lying beyond that energy, (1+x)exp(-x) with x = E_max/kT, is:")
        for  i = 1:length(temperatures)
            flag = weightAbove[i] > 0.1 ? "   <-- alpha(exc-autoion) is a LOWER BOUND here; compute more energies" : ""
            println(stream, "       T = " * @sprintf("%.4e", temperatures[i]) * " K :  " *
                            @sprintf("%8.4f", weightAbove[i]) * flag)
        end
    else
        println(stream, "    No impact-excitation lines are present in these cascade data, so alpha(exc-autoion) is")
        println(stream, "    structurally zero rather than small.  Add Basics.ImpactExcAuto() to the scheme's processes")
        println(stream, "    to compute that channel.")
    end
    println(stream, "  ", "-"^nx)

    return( nothing )
end


"""
`Cascade.simulatePiRateCoefficients(simulation::Cascade.Simulation)`
    ... folds the photoionization cross sections of the cascade data with each of the given photon fields and reports
        the photoionization rate per ion.  An Array{Basics.EmProperty,1} in 1/s, one entry per field, is returned.
"""
function simulatePiRateCoefficients(simulation::Cascade.Simulation)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    property = simulation.property
    if  length(property.photonDistributions) == 0
        error("Cascade.PiRateCoefficients: no photon distribution was given, so there is nothing to fold with.  " *
              "Set photonDistributions, e.g. [Distribution.PhotonPlanck(kT)] with kT in atomic units.")
    end
    results = simulation.computationData[1]["results"]
    if  !haskey(results, "photoionization lines:")
        error("Cascade.PiRateCoefficients: these cascade data carry no photoionization lines.  This property needs a " *
              "computation of Cascade.PhotoIonizationScheme; example-Fd.jl branch a is the smallest one.")
    end
    lines = results["photoionization lines:"]
    if  length(lines) < 2
        error("Cascade.PiRateCoefficients: a fold over photon energy needs at least two computed energies; the given " *
              "data carry $(length(lines)) line(s).  Widen PhotoIonizationScheme.photonEnergies.")
    end
    #
    # sigma^PI(omega), summed over the final levels and over the selected initial level(s)
    energies = sort(unique([l.photonEnergy  for l in lines]))
    csC = Float64[];    csB = Float64[]
    for  om  in  energies
        c1 = 0.;   b1 = 0.
        for  l  in  lines
            if  abs(l.photonEnergy - om) > 1.0e-6 * max(om, 1.0e-10)                                     continue   end
            if  property.initialLevelNo != 0  &&  l.initialLevel.index != property.initialLevelNo        continue   end
            c1 = c1 + l.crossSection.Coulomb;    b1 = b1 + l.crossSection.Babushkin
        end
        push!(csC, c1);    push!(csB, b1)
    end
    if  sum(csC) + sum(csB) == 0.
        error("Cascade.PiRateCoefficients: every cross section is zero for initialLevelNo = " *
              "$(property.initialLevelNo).  Either that level carries no photoionization line in these data, or its " *
              "index does not exist; initialLevelNo = 0 sums over all initial levels.")
    end
    #
    rates = Basics.EmProperty[];    edges = Float64[]
    for  dist  in  property.photonDistributions
        rC, eC = Cascade.foldWithPhotonField(energies, csC, dist)
        rB, _  = Cascade.foldWithPhotonField(energies, csB, dist)
        push!(rates, Basics.EmProperty(rC, rB));    push!(edges, eC)
    end
    #
    Cascade.displayPiRateCoefficients(stdout, energies, rates, edges, property)
    if  printSummary   Cascade.displayPiRateCoefficients(iostream, energies, rates, edges, property)   end

    return( rates )
end


"""
`Cascade.foldWithPhotonField(energies::Array{Float64,1}, values::Array{Float64,1},
                             dist::Distribution.AbstractPhotonDistribution)`
    ... folds a tabulated cross section with a photon field, R = INT d(omega) n(omega) c sigma(omega), by the
        trapezoidal rule over the tabulated energies and with sigma taken as ZERO outside them.  A tuple
        (rate::Float64 in 1/s, edgeShare::Float64) is returned, where edgeShare is the fraction of the integral
        contributed by the two OUTERMOST intervals together.  That fraction is the diagnostic: if the integrand is
        still large at the ends of the computed range, the range is too narrow for this field and the rate is a lower
        bound.  It makes no assumption about the shape of the field, which a closed-form tail estimate would.
"""
function foldWithPhotonField(energies::Array{Float64,1}, values::Array{Float64,1},
                             dist::Distribution.AbstractPhotonDistribution)
    cLight = Defaults.getDefaults("speed of light: c")
    n      = length(energies)
    if  n < 2    return( (0., 1.) )    end
    contrib = Float64[]
    for  i = 1:n-1
        f1 = values[i]   * cLight * Distribution.photonNumberDensity(dist, energies[i])
        f2 = values[i+1] * cLight * Distribution.photonNumberDensity(dist, energies[i+1])
        push!(contrib, 0.5 * (f1 + f2) * (energies[i+1] - energies[i]))
    end
    total = sum(contrib)
    edge  = total == 0. ? 1. : (contrib[1] + contrib[end]) / total
    rate  = Defaults.convertUnits("rate: from atomic to 1/s", total)

    return( (rate, edge) )
end


"""
`Cascade.displayPiRateCoefficients(stream::IO, energies::Array{Float64,1}, rates::Array{Basics.EmProperty,1},
                                   edges::Array{Float64,1}, property::Cascade.PiRateCoefficients)`
    ... displays the photoionization rate per ion for each photon field, with the range folded over and the share of
        the integral carried by the outermost intervals.  Nothing is returned.
"""
function displayPiRateCoefficients(stream::IO, energies::Array{Float64,1}, rates::Array{Basics.EmProperty,1},
                                   edges::Array{Float64,1}, property::Cascade.PiRateCoefficients)
    nx = 118
    eMin = Defaults.convertUnits("energy: from atomic to eV", energies[1])
    eMax = Defaults.convertUnits("energy: from atomic to eV", energies[end])
    println(stream, " ")
    println(stream, "  Photoionization rate per ion  R^PI = INT d(omega) n(omega) c sigma^PI(omega):")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "     photon field                                              R^PI (Coulomb)   R^PI (Babushkin)" *
                    "   edge share")
    println(stream, "                                                                    [1/s]             [1/s]     ")
    println(stream, "  ", "-"^nx)
    for  i = 1:length(rates)
        ## Each field prints a whole sentence describing itself, with its temperature in ATOMIC units.  Truncating
        ## that to fit the column removed exactly what distinguishes two fields of the same kind -- two Planck
        ## entries became identical labels against different numbers.  The type and the temperature in eV are what
        ## the reader needs, so they are built here rather than taken from the sentence.
        dist = property.photonDistributions[i]
        lab  = replace(string(typeof(dist)), "JenaAtomicCalculator." => "", "Distribution." => "")
        if  hasproperty(dist, :T)
            lab = lab * @sprintf("  kT = %.1f eV", Defaults.convertUnits("energy: from atomic to eV", dist.T))
        end
        if  hasproperty(dist, :w)    lab = lab * @sprintf(",  w = %.3e", dist.w)    end
        sa = "     " * rpad(lab, 54)
        println(stream, sa * @sprintf("%.6e", rates[i].Coulomb) * "    " * @sprintf("%.6e", rates[i].Babushkin) *
                        "     " * @sprintf("%8.4f", edges[i]))
    end
    println(stream, "  ", "-"^nx)
    println(stream, " ")
    println(stream, "    THIS IS A RATE [1/s] AND NOT A RATE COEFFICIENT.  The convolution already carries the photon")
    println(stream, "    number density of the field, so it needs no further multiplication by a density; multiply by")
    println(stream, "    the ION number density for a volumetric rate.  The electron density does not enter at all.")
    println(stream, " ")
    println(stream, "    THE FOLD IS OVER THE COMPUTED ENERGIES ONLY, " * @sprintf("%.3f", eMin) * " to " *
                    @sprintf("%.3f", eMax) * " eV, with sigma^PI taken as ZERO")
    println(stream, "    outside them: no extrapolation to threshold and none to high energy.  The last column is the")
    println(stream, "    share of the integral carried by the two OUTERMOST intervals together, which is the honest")
    println(stream, "    test of whether that range suits the field -- a large share means the integrand is still big")
    println(stream, "    where the data stop, and the rate is then a LOWER BOUND.  Widening it means recomputing the")
    println(stream, "    cascade with more photonEnergies; it cannot be repaired at this stage.")
    println(stream, "    Resonant photoabsorption is NOT included here; use Cascade.PhotoAbsorptionSpectrum for that.")
    println(stream, "  ", "-"^nx)

    return( nothing )
end


"""
`Cascade.directIonizationAlpha(Z::Float64, conf::Configuration, temperatures::Array{Float64,1})`
    ... estimates the DIRECT electron-impact ionization rate coefficient semi-empirically, by summing the Lotz rate of
        Empirical.impactIonizationPlasmaAlpha over every occupied subshell of `conf`.  A tuple
        (alphas::Array{Float64,1} in cm^3/s, skipped::Array{String,1}) is returned, where `skipped` names any subshell
        the empirical binding-energy tables could not supply, together with the reason -- those subshells contribute
        nothing and the caller must be able to say which, since a silently dropped inner shell would lower the total
        without any sign of it.

        THIS IS NOT A CASCADE COMPUTATION and is not of the same kind as the other channels.  It is a fit: Lotz (1967),
        folded with a Maxwellian, using tabulated binding energies.  Its accuracy is tens of per cent at best and worse
        near threshold and for near-neutral ions.  It is here because the alternative is worse -- without it a total
        that omits the direct channel is wrong by ORDERS OF MAGNITUDE for a light ion, where direct ionization carries
        98% or more of the rate, and a number wrong by 100 is more misleading than one wrong by 30%.
"""
function directIonizationAlpha(Z::Float64, conf::Configuration, temperatures::Array{Float64,1})
    alphas = Float64[];   skipped = String[]
    conv   = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    shells = sort(collect(keys(conf.shells)), by = sh -> (sh.n, sh.l))
    ## Empirical.impactIonizationPlasmaAlpha reads the nuclear charge from the global defaults rather than taking it as
    ## an argument, so it is set here and restored afterwards; leaving it changed would silently alter whatever the
    ## caller does next.
    oldZ = Defaults.getDefaults("nuclear: charge")
    Defaults.setDefaults("nuclear: charge", Z)
    try
        for  temp  in  temperatures
            tAu = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
            wa  = 0.
            for  sh  in  shells
                occ = conf.shells[sh];      if  occ < 1    continue    end
                d   = copy(conf.shells);    d[sh] = occ - 1
                fConf = Configuration(d, conf.NoElectrons - 1)
                try
                    wa = wa + redirect_stdout(devnull) do
                             Empirical.impactIonizationPlasmaAlpha(Distribution.ElectronMaxwell(tAu), conf, fConf)
                         end
                catch  ex
                    ## An UndefVarError or a MethodError is a fault in THIS code, not a gap in the tables, and must
                    ## not be recorded as a physics limitation -- the first version of this catch reported exactly
                    ## that, turning a missing `using ..Distribution` into "the empirical tables do not cover them".
                    ## Those two are re-thrown; only a genuine failure of the empirical routine is skipped, and the
                    ## reason it gave is reported verbatim rather than interpreted.
                    if  ex isa UndefVarError  ||  ex isa MethodError    rethrow(ex)    end
                    sa = "$(sh): " * first(split(sprint(showerror, ex), "\n"))
                    if  !(sa in skipped)    push!(skipped, sa)    end
                end
            end
            push!(alphas, conv * wa)
        end
    finally
        Defaults.setDefaults("nuclear: charge", oldZ)
    end

    return( (alphas, skipped) )
end


"""
`Cascade.extractEaCrossSections(simulation::Cascade.Simulation)`
    ... collects the excitation-autoionization cross section sigma^EA(E) from the impact-excitation and autoionization
        lines of the cascade data, summing over those excited levels that actually autoionize.  A tuple
        (energies::Array{Float64,1}, values::Array{Float64,1}) in ATOMIC UNITS is returned, sorted by energy and empty
        if the data carry no impact-excitation lines at all -- which is the normal case for a purely resonant
        computation and is not an error.
"""
function extractEaCrossSections(simulation::Cascade.Simulation)
    results = simulation.computationData[1]["results"]
    if  !haskey(results, "impact-excitation lines:")     return( (Float64[], Float64[]) )   end
    linesE = results["impact-excitation lines:"]
    linesA = haskey(results, "autoionization lines:") ? results["autoionization lines:"] : AutoIonization.Line[]
    if  length(linesE) == 0  ||  length(linesA) == 0     return( (Float64[], Float64[]) )   end
    ## An excited level counts as autoionizing exactly if it appears as the INITIAL level of an Auger line; its
    ## branching ratio is taken as 1, as in Cascade.EaCrossSections, so this half is an UPPER BOUND.
    autoIonizing = unique([ (l.initialLevel.index, l.initialLevel.energy)  for l in linesA ])
    energies     = sort(unique([l.initialElectronEnergy  for l in linesE]))
    values       = Float64[]
    for  en  in  energies
        cs = 0.
        for  l  in  linesE
            if  abs(l.initialElectronEnergy - en) / max(en, 1.0e-10) > 1.0e-6           continue    end
            if  !( (l.finalLevel.index, l.finalLevel.energy)  in  autoIonizing )         continue    end
            cs = cs + l.crossSection
        end
        push!(values, cs)
    end

    return( (energies, values) )
end


"""
`Cascade.foldWithMaxwellian(energies::Array{Float64,1}, values::Array{Float64,1}, temperatures::Array{Float64,1})`
    ... folds a tabulated cross section sigma(E) with a Maxwellian electron distribution at each temperature,

            alpha(T)  =  4/sqrt(2pi) * T^(-3/2) * INT sigma(E) * E * exp(-E/T) dE ,

        the SAME prefactor that DielectronicRecombination.computeRateCoefficient carries, since sqrt(8/pi) and
        4/sqrt(2pi) are identically equal; a delta-like sigma therefore reproduces the isolated-resonance formula
        exactly.  The integral is a trapezoidal sum over the tabulated energies and is TRUNCATED above the largest of
        them.  A tuple (alphas::Array{Float64,1} in cm^3/s, weightAbove::Array{Float64,1}) is returned, where
        weightAbove[i] = (1 + x) exp(-x) with x = E_max/T_i is the exact share of the Maxwellian weight lying beyond
        the last tabulated energy -- the number that says whether the truncation matters at that temperature.
"""
function foldWithMaxwellian(energies::Array{Float64,1}, values::Array{Float64,1}, temperatures::Array{Float64,1})
    alphas = Float64[];    weightAbove = Float64[]
    conv   = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 *
             Defaults.convertUnits("rate: from atomic to 1/s", 1.0)
    for  temp  in  temperatures
        tAu = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
        if  length(energies) < 2   push!(alphas, 0.);   push!(weightAbove, length(energies) == 0 ? 1. : 0.);  continue  end
        wa = 0.
        for  i = 1:length(energies)-1
            f1 = values[i]   * energies[i]   * exp(-energies[i]/tAu)
            f2 = values[i+1] * energies[i+1] * exp(-energies[i+1]/tAu)
            wa = wa + 0.5 * (f1 + f2) * (energies[i+1] - energies[i])
        end
        push!(alphas, conv * 4/sqrt(2pi) * tAu^(-3/2) * wa)
        x = energies[end]/tAu;     push!(weightAbove, (1. + x) * exp(-x))
    end

    return( (alphas, weightAbove) )
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
                                                                 AutoIonization.PartialWave[])
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
                                                                         nFactor * nStrength, AutoIonization.PartialWave[])
                    push!(resonances, newResonance)
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
        parameters; `printout` decides whether the result is also displayed, but never whether it is returned. An 
        `Array{Basics.EmProperty,1}` in [cm^2/g], one entry per `property.dependencyValues`, is returned.
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
        ## The Eastman & Pinto factor is lambda/Delta-lambda, i.e. the INVERSE relative bin width. Since
        ## |Delta-lambda|/lambda = |Delta-omega|/omega exactly, the same factor omega/Delta-omega must be used
        ## in every dependence -- the frequency and temperature branches returned Delta-omega/omega instead,
        ## which is wrong by (omega/Delta-omega)^2, typically four orders of magnitude.
        if     typeof(opacityDependence) == FrequencyOpacityDependence
            # Binning, depValue and lineOmega are all in Hartree and readily to compare.
            if  depValue - halfBinning < lineOmega < depValue + halfBinning             wa = lineOmega / (2* halfBinning)      end
        elseif typeof(opacityDependence) == WavelengthOpacityDependence
            # Binning in nm, lineOmega & depValue in Hartree are first converted into nm ... and ratio is determined in nm as well
            lineOmega_nm = convertUnits("energy: from atomic to Angstrom", lineOmega) / 10.
            depValue_nm  = convertUnits("energy: from atomic to Angstrom", depValue)  / 10.
            if  depValue_nm - halfBinning < lineOmega_nm < depValue_nm + halfBinning    wa = lineOmega_nm / (2* halfBinning)   end
        elseif typeof(opacityDependence) == TemperatureOpacityDependence
            # Binning and lineOmega are in Hartree, depValue in [u = omega/kT] and is converted here. The upper bound
            # read  depValue + halfBinning * kT,  i.e. the kT was attached to the wrong factor and the bin was not even
            # centred on depValue * kT.
            if  depValue * kT - halfBinning < lineOmega < depValue * kT + halfBinning   wa = lineOmega / (2* halfBinning)      end
        else   error("Unsupported opacityDependence :: $(typeof(opacityDependence)).")
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
    rho        = property.massDensity          ## [g/cm^3], enters the 1/(rho c t) prefactor
    eshift     = property.transitionEnergyShift
    ## The ion number density enters the Sobolev optical depth and must be in ATOMIC units there (1/a_0^3).
    ## It used to be hard-wired as `ne = 1.0` with a "??" comment, so the input density never reached tau and
    ## the absolute scale of every opacity was arbitrary.
    a0_in_cm   = convertUnits("length: from atomic to fm", 1.0) * 1.0e-13
    nion_au    = property.ionNumberDensity * a0_in_cm^3          ## [1/cm^3] -> [1/a_0^3]
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
    ## The LTE population of each LOWER level, normalised over the levels that actually occur in this line
    ## list.  The former expression used  ge/g0 * exp(-omega/kT), which is not a population of the lower level
    ## at all: there is no partition function, omega is the TRANSITION energy rather than the excitation energy
    ## of the absorbing level, and the statistical-weight ratio is inverted.  property.levelPopulation is now
    ## dispatched on, so that the field means something.
    if  typeof(property.levelPopulation) != Basics.BoltzmannLevelPopulation
        error("Cascade.ExpansionOpacities presently supports levelPopulation = BoltzmannLevelPopulation() only; " *
              "got $(property.levelPopulation).")
    end
    levelEnergies = Float64[];    levelWeights = Float64[]
    for  excData in photoexcitationData,  line in excData.lines
        for  lev  in  [line.initialLevel, line.finalLevel]
            if  !any(abs.(levelEnergies .- lev.energy) .< 1.0e-12)
                push!(levelEnergies, lev.energy);   push!(levelWeights, Basics.twice(lev.J) + 1.0)
            end
        end
    end
    eGround      = minimum(levelEnergies)
    partitionFct = sum( levelWeights .* exp.(-(levelEnergies .- eGround) ./ kT) )
    nLower(lev)  = nion_au * (Basics.twice(lev.J) + 1.0) * exp(-(lev.energy - eGround)/kT) / partitionFct
    #
    minEnergy = 1000.;   maxEnergy = 0.
    for  excData  in photoexcitationData
        ## Cascade.Data{T} carries its lines in `lines`; `.linesE` belonged to the retired
        ## Cascade.ExcitationData and raised a FieldError here (the same class as the RR `.linesR` bug).
        for  line  in excData.lines
            if  minEnergy > line.omega     minEnergy = line.omega   end
            if  maxEnergy < line.omega     maxEnergy = line.omega   end
                omega       = line.omega + eshift
                fosc        = line.oscStrength
                lambda_au   = convertUnits("energy: from atomic to Angstrom", omega) * A_au
                nl          = nLower(line.initialLevel)
            for  ivalue = 1:NoValues
                lmd_over_dl    = lambda_over_dlambda(dependence, omega, kT, values[ivalue])
                if  lmd_over_dl == 0.  continue  end
                ## Sobolev optical depth  tau_l = pi alpha n_l lambda_l t_exp f_l   (atomic units)
                tau_Cou        = pi * alpha * nl * lambda_au * exptime_au * fosc.Coulomb
                tau_Bab        = pi * alpha * nl * lambda_au * exptime_au * fosc.Babushkin
                term_Cou       = factor * lmd_over_dl * (1.0 - exp(-tau_Cou))
                term_Bab       = factor * lmd_over_dl * (1.0 - exp(-tau_Bab))
                kappas[ivalue] = kappas[ivalue] + Basics.EmProperty(term_Cou, term_Bab)
            end
        end
    end
    #
    ## `printout` controls PRINTING only. It used to control the return value as well, so the printing path
    ## handed back nothing and perform(::Simulation) filled results["data:"] with nothing -- the opacities
    ## could be read on screen but not retrieved, unlike every other simulation property.
    if  printout
        Cascade.displayExpansionOpacities(stdout, name, property, (minEnergy,maxEnergy), kappas)
        if  printSummary
            Cascade.displayExpansionOpacities(iostream, name, property, (minEnergy,maxEnergy), kappas)      end
    end

    return( kappas )
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
`Cascade.simulateMeanOpacities(photoexcitationData::Array{Cascade.Data,1}, simulation::Cascade.Simulation)`
    ... runs through all excitation lines, assembles the spectral opacity kappa(u_i) on the eight Gauss-Laguerre
        nodes of the temperature-normalised variable u = omega/kT, and reduces it to the MEAN named by
        simulation.property.opacityMean -- RosselandMean() or PlanckMean(). One value in [cm^2/g] is printed
        for each (rho, T) pair, in both gauges. An `Array{Basics.EmProperty,1}` carrying those means, in the order
        printed -- mass densities outermost, temperatures innermost -- is returned.
"""
function simulateMeanOpacities(photoexcitationData::Array{Cascade.Data,1}, simulation::Cascade.Simulation)
    ulist, weights = Cascade.rosselandWeights(8)
    property       = simulation.property

    if  length(property.ionNumberDensities) != length(property.massDensities)
        error("Cascade.MeanOpacities needs one ion NUMBER density per MASS density; got " *
              "$(length(property.ionNumberDensities)) and $(length(property.massDensities)).")
    end

    means = Basics.EmProperty[]
    for  (id, rho)  in  enumerate(property.massDensities)
        nion = property.ionNumberDensities[id]
        for  T in property.temperatures
            ## Assemble kappa_nu as the SUM over the requested contributions; each returns one value per node.
            kappas = Basics.EmProperty[ Basics.EmProperty(0.)  for i = 1:length(ulist) ]
            for  contribution  in  property.contributions
                wa = Cascade.spectralOpacityContribution(contribution, photoexcitationData, property,
                                                         nion, rho, T, ulist)
                for  i = 1:length(ulist)     kappas[i] = kappas[i] + wa[i]     end
            end
            #
            mean, note = Cascade.applyOpacityMean(property.opacityMean, kappas, weights)
            sa = "> " * Cascade.nameOfOpacityMean(property.opacityMean) * " opacity for rho = " *
                    @sprintf("%.3e",rho) * " [g/cm^3],  n_ion = " * @sprintf("%.3e",nion) * " [1/cm^3]  &  T = " *
                    @sprintf("%.3e",T) * " [K]  is  kappa [cm^2/g] = " * @sprintf("%.5e",mean.Coulomb) *
                    " [Coulomb]  " * @sprintf("%.5e",mean.Babushkin) * " [Babushkin]"
            println(sa)
            if  note != ""      println("    ++ " * note)      end
            push!(means, mean)
        end
    end

    return( means )
end


"""
`Cascade.rosselandWeights(noNodes::Int64)`
    ... returns the Gauss-Laguerre nodes u_i of the temperature-normalised variable u = omega/kT together with
        the NORMALISED Rosseland weights

            w(u) = 15/(4 pi^4) u^4 e^-u / (1 - e^-u)^2 ,        int_0^inf w(u) du = 1 .

        Gauss-Laguerre supplies the factor e^-u, so only 15/(4 pi^4) u^4 / (1 - e^-u)^2 is evaluated against
        its weights. The prefactor is EXACT rather than fitted, because
        int_0^inf u^4 e^u/(e^u-1)^2 du = 4 int_0^inf u^3/(e^u-1) du = 4 Gamma(4) zeta(4) = 4 pi^4/15.
        Measured 14-Aug-2026: eight nodes give sum_i w_i = 1 to 2.8e-6, sixteen to 1.8e-11. Eight is ample
        here, and that residue is the accuracy floor of every mean built on these weights -- the grey and
        Thomson known-answer checks in TestFrames.testMethod_Opacities reproduce their exact values to
        precisely this 2.8e-6 and to nothing better, which is how one knows the deviation is the quadrature
        and not the physics.

        A tuple (nodes::Array{Float64,1}, weights::Array{Float64,1}) is returned.
"""
function rosselandWeights(noNodes::Int64)
    ulist, wlist = FastGaussQuadrature.gausslaguerre(noNodes)
    weights      = [ 15.0/(4*pi^4) * ulist[i]^4 * wlist[i] / (1.0 - exp(-ulist[i]))^2   for i = 1:noNodes ]
    return( (ulist, weights) )
end


"""
`Cascade.spectralOpacityContribution(contribution::Cascade.BoundBoundOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})`
    ... returns the bound-bound (line) contribution to the spectral opacity at the given nodes, by handing the
        work to Cascade.simulateExpansionOpacities. An Array{Basics.EmProperty,1} in [cm^2/g] is returned.
"""
function spectralOpacityContribution(contribution::Cascade.BoundBoundOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})
    expansion = Cascade.ExpansionOpacities(contribution.levelPopulation, property.opacityDependence,
                                           nion, rho, T, property.expansionTime,
                                           property.transitionEnergyShift, ulist)
    return( Cascade.simulateExpansionOpacities(photoexcitationData, "spectral opacity for rho = $rho & T = $T",
                                               expansion, printout=false) )
end


"""
`Cascade.spectralOpacityContribution(contribution::Cascade.BoundFreeOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})`
    ... returns the bound-free contribution kappa_nu^bf = n_ion sum_shells sigma^PI(omega) / rho at the given
        nodes, with the cross sections from Empirical.photoionizationCrossSection in the ScaledHydrogenic
        approximation. The nodes u_i are converted to photon energies omega_i = u_i kT, and the cross
        sections from atomic units (a_0^2) to cm^2. Being a continuum, this contribution carries NO
        lambda/Delta-lambda bin factor. An Array{Basics.EmProperty,1} in [cm^2/g] is returned; the estimate
        is gauge-independent, so both components carry the same value.

        Empirical.photoionizationCrossSection reads the nuclear charge from the GLOBAL Defaults, which is
        1.0 unless set. It is therefore set here from contribution.nuclearCharge and restored afterwards,
        so that a forgotten setDefaults cannot silently produce hydrogen cross sections for another element.
"""
function spectralOpacityContribution(contribution::Cascade.BoundFreeOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})
    kT       = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", T)
    omegas   = [ u * kT   for u in ulist ]                                          ## [Hartree]
    a0_in_cm = Defaults.convertUnits("length: from atomic to fm", 1.0) * 1.0e-13     ## [cm]
    kappas   = Basics.EmProperty[ Basics.EmProperty(0.)  for i = 1:length(ulist) ]
    #
    zOld = Defaults.getDefaults("nuclear: charge")
    Defaults.setDefaults("nuclear: charge", contribution.nuclearCharge)
    try
        for  (iConf, fConf)  in  contribution.ionizationPairs
            css = Empirical.photoionizationCrossSection(omegas, iConf, fConf, Empirical.ScaledHydrogenic())
            for  i = 1:length(ulist)
                wa = nion * css[i] * a0_in_cm^2 / rho                                ## [cm^2/g]
                kappas[i] = kappas[i] + Basics.EmProperty(wa, wa)
            end
        end
    finally
        Defaults.setDefaults("nuclear: charge", zOld)
    end
    return( kappas )
end


"""
`Cascade.spectralOpacityContribution(contribution::Cascade.ScatteringOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})`
    ... returns the Thomson-scattering contribution kappa^sc = n_e sigma_T / rho, which is INDEPENDENT of
        frequency and therefore identical in every bin, and gauge-independent so that both components of the
        EmProperty carry the same value. An Array{Basics.EmProperty,1} in [cm^2/g] is returned.
"""
function spectralOpacityContribution(contribution::Cascade.ScatteringOpacity, photoexcitationData::Array{Cascade.Data,1},
                                     property::Cascade.MeanOpacities, nion::Float64, rho::Float64, T::Float64,
                                     ulist::Array{Float64,1})
    sigmaThomson = 6.6524587321e-25                              ## [cm^2]
    kappaSc      = contribution.electronNumberDensity * sigmaThomson / rho       ## [cm^2/g]
    return( Basics.EmProperty[ Basics.EmProperty(kappaSc, kappaSc)  for i = 1:length(ulist) ] )
end


"""
`Cascade.nameOfOpacityMean(mean::Cascade.AbstractOpacityMean)`
    ... returns the name of the given mean for printout. A String is returned.
"""
nameOfOpacityMean(mean::Cascade.RosselandMean)  = "Rosseland"
nameOfOpacityMean(mean::Cascade.PlanckMean)     = "Planck"


"""
`Cascade.applyOpacityMean(mean::Cascade.PlanckMean, kappas::Array{Basics.EmProperty,1}, weights::Array{Float64,1})`
    ... forms the ARITHMETIC, weight-normalised mean sum_i w_i kappa_i of the given spectral opacities.
        A tuple (value::Basics.EmProperty, note::String) is returned; note is empty for this mean, which is
        finite whatever the spectral opacity looks like.
"""
function applyOpacityMean(mean::Cascade.PlanckMean, kappas::Array{Basics.EmProperty,1}, weights::Array{Float64,1})
    wa = Basics.EmProperty(0.)
    for  i = 1:length(kappas)     wa = wa + weights[i] * kappas[i]     end
    return( (wa, "") )
end


"""
`Cascade.applyOpacityMean(mean::Cascade.RosselandMean, kappas::Array{Basics.EmProperty,1}, weights::Array{Float64,1})`
    ... forms the HARMONIC, weight-normalised mean  1 / sum_i (w_i / kappa_i)  of the given spectral
        opacities. A tuple (value::Basics.EmProperty, note::String) is returned.

        An EMPTY BIN (kappa_i = 0) makes the harmonic mean vanish identically, and that is physics rather
        than a numerical accident: the Rosseland mean is dominated by the most transparent frequencies, so a
        single perfectly transparent window short-circuits it. Rather than silently returning zero -- or,
        worse, substituting an arithmetic mean, as JAC did before 14-Aug-2026 -- the vanishing case is
        reported through `note`, naming how many bins are empty and what is missing: a continuum
        contribution (electron scattering, bound-free) that a bound-bound line list alone cannot supply.
"""
function applyOpacityMean(mean::Cascade.RosselandMean, kappas::Array{Basics.EmProperty,1}, weights::Array{Float64,1})
    nEmptyCou = count(k -> k.Coulomb    <= 0., kappas)
    nEmptyBab = count(k -> k.Babushkin  <= 0., kappas)
    note      = ""
    if  nEmptyCou > 0  ||  nEmptyBab > 0
        note = "the Rosseland mean VANISHES: $(max(nEmptyCou,nEmptyBab)) of $(length(kappas)) bins carry no " *
               "opacity at all. A harmonic mean is dominated by its most transparent bin, so this is the " *
               "correct value for the given spectral opacity -- and it says that a bound-bound line list " *
               "alone cannot support a Rosseland mean. Add a continuum contribution (electron scattering, " *
               "bound-free), or ask for a PlanckMean(), which stays finite on an incomplete line list."
        return( (Basics.EmProperty(0.), note) )
    end
    sumCou = sum( weights[i] / kappas[i].Coulomb    for i = 1:length(kappas) )
    sumBab = sum( weights[i] / kappas[i].Babushkin  for i = 1:length(kappas) )
    return( (Basics.EmProperty(1.0/sumCou, 1.0/sumBab), note) )
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
