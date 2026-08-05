
# Functions and methods for cascade computation


function Base.isless(x::Cascade.AbsorptionCrossSection, y::Cascade.AbsorptionCrossSection)
    return x.photonEnergy < y.photonEnergy
end

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

    
"""
`Basics.extractConfiguration(theme::Basics.LeadingConfiguration, cLevel::Cascade.Level)`  
    ... extract the leading configuration of the given level; a conf::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.LeadingConfiguration, cLevel::Cascade.Level)
    # First extract the right ManyElectron.Level that is to be analyzed
    level = Level()
    if length(cLevel.parents) > 0   
        parent = cLevel.parents[1]
        if      parent.process == Basics.Auger()       level = parent.lines[parent.index].finalLevel
        elseif  parent.process == Basics.Radiative()   level = parent.lines[parent.index].finalLevel
        else    error("stop a")
        end
    elseif length(cLevel.daugthers) > 0   
        daugth = cLevel.daugthers[1]
        if      daugth.process == Basics.Auger()       level = daugth.lines[daugth.index].initialLevel
        elseif  daugth.process == Basics.Radiative()   level = daugth.lines[daugth.index].initialLevel
        else    error("stop b")
        end
    else        error("stop c")
    end

    allConfs = Basics.extractConfigurations(Basics.FromBasis(), level.basis)
    weights  = zeros(length(allConfs))
    for  (ia, allConf) in  enumerate(allConfs)
        for (ic, csf) in enumerate(level.basis.csfs)
            if  allConf == Basics.extractConfiguration(Basics.FromBasis(), level.basis, csf)     
                weights[ia] = weights[ia] + level.mc[ic]
            end
        end
    end
    # Determine index of maximum and return the corresponding configuration
    wx   = findmax(weights)
    conf = allConfs[ wx[2] ]
    
    return( conf )
end

    
"""
`Basics.isSimilar()`  ... returns true if two instances are similar to each other, and false otherwise.

+ `(keya::LevelKey, keyb::LevelKey, relAcc::Float64)`  
    ... returns true if two level keys refer to the same level, i.e. level with the same symmetry and
        if the relative energy abs( (E_a - E_b)/E_a ) < relAcc. It returns false otherwise.
"""
function Basics.isSimilar(keya::LevelKey, keyb::LevelKey, relAcc::Float64)
    if  keya.sym == keyb.sym   &&   abs( (keya.energy - keyb.energy)/keya.energy ) < relAcc    return(true)
    else                                                                                       return(false)
    end
end


"""
`Basics.perform(comp::Cascade.Computation)`  
    ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
        various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
        to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
        otherwise.

`Basics.perform(comp::Cascade.Computation; output::Bool=true, outputToFile::Bool=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Cascade.Computation; output::Bool=false, outputToFile::Bool=true,
                        outputDirectory::String="")
    ## Reject an unrealizable cascade approach HERE, before any orbital or amplitude work is done; a cascade
    ## that cannot succeed should say so in the first second rather than in the middle of some step.
    Cascade.validateConditions(comp.approach)
    Cascade.perform(comp.scheme, comp::Cascade.Computation, output=output, outputToFile=outputToFile,
                    outputDirectory=outputDirectory)
end



"""
`Basics.perform(comp::Cascade.Simulation)`  
    ... to set-up and perform a cascade computation that starts from a given set of initial configurations and proceeds via 
        various steps until a given number of electrons has been removed or the decay stops at some stable levels with regard 
        to the given atomic processes. The results of all individual steps are printed to screen but nothing is returned 
        otherwise.

`Basics.perform(comp::Cascade.Simulation; output=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the type 
        and specifications of the cascade but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Cascade.Simulation; output::Bool=false)
    Cascade.perform(comp::Cascade.Simulation, output=output)
end


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Cascade.computeDecayProbabilities(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                    linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)` 
    ... to compute the decay probabilities for all pairs and triples of subshells; these probabilities only depend
        on the holes in different subshells and are sumed over all final levels that share the same subshell occupation.
        The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeDecayProbabilities(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                    linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)
    #
    # Determine the leading configuration and shells of the initial level
    subshList     = Basics.extractRelativisticSubshellList(outcome.level)
    relConfigs    = Basics.extractConfigurations(Basics.RelativisticConfigurations(), outcome.level.basis, outcome.level.J)
    holeSubshells = Basics.extractFromConfiguration(Basics.OpenSubshells(), relConfigs[1]);   holeSubshell = holeSubshells[1]
    #
    # Initialize and fill dictionaries for collecting the radiative and Auger rates
    level = outcome.level;    sym = LevelSymmetry(level.J, level.parity)
    if  subshList[1] != Subshell("1s_1/2")    error("stop a")     end
    println("\n Probabilities are determined for the initial hole $holeSubshell of a level with symmetry $sym " *
            "and configurations \n  $(relConfigs)")
    subshEnergies  = Dict{Subshell,Float64}()
    rProbabilities = Dict{Subshell,Float64}()
    aProbabilities = Dict{Tuple{Subshell,Subshell},Float64}()
    for  (i,subshi)  in  enumerate(subshList)
        subshEnergies[subshi]  = outcome.level.basis.orbitals[subshi].energy
        rProbabilities[subshi] = 0.
        for  (j,subshj)  in  enumerate(subshList)
            if  j >= i  aProbabilities[(subshi, subshj)] = 0.   end 
        end
    end
    #
    @show subshList, relConfigs, holeSubshell
    @show subshEnergies
    #
    # Calculate the total rate and convert the dictionaries into probabilities
    # First, identify the level key of the given level also in the lists of radiative and Auger lines
    level = outcome.level;    levelKey = LevelKey( LevelSymmetry(level.J, level.parity), level.index, level.energy, 0.)
    similarKey = LevelKey();  rateR = Basics.EmProperty(0.);    rateA = 0.;   NoPhotonLines = 0;   NoAugerLines = 0
    for  line in linesR
        compareKey = LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
        if   Basics.isSimilar(levelKey, compareKey, 1.0e-3)    println("** compareKey = $compareKey");   similarKey = deepcopy(compareKey)    end
    end
    if   similarKey == LevelKey()    @warn("No similar level found !")   end
    
    for  line in linesR
        if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            rateR = rateR + line.photonRate;   NoPhotonLines = NoPhotonLines + 1  
            # Extract the subshells that make the difference between the initial and final levels; terminate if NO 1s_1/2 occurs
            confi = Basics.extractConfiguration(Basics.LeadingConfigurationR(), line.initialLevel)
            conff = Basics.extractConfiguration(Basics.LeadingConfigurationR(), line.finalLevel)
            occDiffs = Basics.extractFromConfigurations(Basics.OccupationDifference(), confi, conff);   subshList = Subshell[]
            for  diff in occDiffs
                # holeSubshell must differ by -1
                if diff[1] == holeSubshell
                    if diff[2] == -1     else   error("stop b")    end
                end
                if     diff[2] ==  1     push!(subshList, diff[1])                                end       
                if     diff[2] ==  2     push!(subshList, diff[1]);    push!(subshList, diff[1])  end       
            end
            if   length(subshList)  != 1   error("stop c")    end
            rProbabilities[ subshList[1] ] = rProbabilities[ subshList[1] ] + 
                                                (line.photonRate.Babushkin + line.photonRate.Coulomb) / 2
        end
    end
    @show rProbabilities
    #
    for  line in linesA
        if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            rateA = rateA + line.totalRate;   NoAugerLines = NoAugerLines + 1    
            # Extract the subshells that make the difference between the initial and final levels; terminate if NO 1s_1/2 occurs
            confi = Basics.extractConfiguration(Basics.LeadingConfigurationR(), line.initialLevel)
            conff = Basics.extractConfiguration(Basics.LeadingConfigurationR(), line.finalLevel)
            occDiffs = Basics.extractFromConfigurations(Basics.OccupationDifference(), confi, conff);   subshList = Subshell[]
            for  diff in occDiffs
                # holeSubshell must differ by -1
                if diff[1] == holeSubshell
                    if diff[2] == -1     else   error("stop d")    end
                end
                if     diff[2] ==  1     push!(subshList, diff[1])                                end       
                if     diff[2] ==  2     push!(subshList, diff[1]);    push!(subshList, diff[1])  end       
            end
            if   length(subshList)  != 2   error("stop e")    end
            aProbabilities[ (subshList[1], subshList[2]) ] = aProbabilities[ (subshList[1], subshList[2]) ] + line.totalRate
        end
    end
    @show aProbabilities
    #
    # Add all rates, convert to probabilities and print these dictionaries into a neat table
    rate = 0.
    for (k,v) in rProbabilities   rate = rate + v   end
    for (k,v) in aProbabilities   rate = rate + v   end
    @show outcome.rateR + outcome.rateA, rate
    #
    for (k,v) in rProbabilities   rProbabilities[k] = v/(rate + 1.0e-20)   end
    for (k,v) in aProbabilities   aProbabilities[k] = v/(rate + 1.0e-20)   end
    #
    Cascade.displayDecayProbabilities(stdout, outcome, rProbabilities, aProbabilities, settings)
    # Generate an additional output file if printDataX & printDataY = true
    printDataX, iodataX = Defaults.getDefaults("data-X flag/stream")
    printDataY, iodataY = Defaults.getDefaults("data-Y flag/stream")
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary  Cascade.displayDecayProbabilities(iostream, outcome, rProbabilities, aProbabilities, settings)  end
    if  printDataX    Cascade.dumpDecayProbabilities(iodataX, outcome, subshEnergies, rProbabilities, settings)       end
    if  printDataY    Cascade.dumpDecayProbabilities(iodataY, outcome, subshEnergies, aProbabilities, settings)       end
        
    return( nothing )
end


"""
`Cascade.computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                    linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)` 
    ... to compute the flourescence and Auger yields for a single decay yield outcome as specified by the corresponding
        level; an outcome::DecayYield.Outcome is returned in which all physical parameters are now specified for the given
        decay-yield.
"""
function computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                    linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)
    # Identify the level key of the given level also in the lists of radiative and Auger lines
    level = outcome.level;    levelKey = LevelKey( LevelSymmetry(level.J, level.parity), level.index, level.energy, 0.)
    @show levelKey
    similarKey = LevelKey();  rateR = Basics.EmProperty(0.);    rateA = 0.;   NoPhotonLines = 0;   NoAugerLines = 0
    for  line in linesR
        compareKey = LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
        @show compareKey
        if   Basics.isSimilar(levelKey, compareKey, 1.0e-3)    println("** compareKey = $compareKey");   similarKey = deepcopy(compareKey)    end
    end
    if   similarKey == LevelKey()    @warn("No similar level found !")   end
    
    for  line in linesR
        if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            rateR = rateR + line.photonRate;   NoPhotonLines = NoPhotonLines + 1    
        end
    end
    for  line in linesA
        if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            rateA = rateA + line.totalRate;   NoAugerLines = NoAugerLines + 1    
        end
    end
    
    omegaR     = rateR / (rateR + rateA + 1.0e-20);   omegaA = rateA / (rateR + rateA + 1.0e-20)  # add 1.0e-20 to avoid division by 0.
    newOutcome = DecayYield.Outcome(level, NoPhotonLines, NoAugerLines, rateR, rateA, omegaR, omegaA)
    return( newOutcome )
end


"""
`Cascade.computeTotalAugerRate(level::Cascade.Level)` 
    ... computes the total Auger rate of level as given by its daugther levels; a rate::Float64 is returned.
"""
function computeTotalAugerRate(level::Cascade.Level)
    rate = 0.
    for daugther in level.daugthers
        if  daugther.process != Basics.Auger();      continue   end
        aLine   = daugther.lines[daugther.index]
        rate    = rate + aLine.totalRate
    end
    return( rate )
end


"""
`Cascade.computeTotalPhotonRate(level::Cascade.Level)` 
    ... computes the total photon (radiative) rate of level as given by its daugther levels; a rate::EmProperty is returned.
"""
function computeTotalPhotonRate(level::Cascade.Level)
    rate = Basics.EmProperty(0.)
    for daugther in level.daugthers
        if  daugther.process != Basics.Radiative();      continue   end
        rLine   = daugther.lines[daugther.index]
        rate    = rate + rLine.photonRate
    end
    
    return( rate )
end



"""
`Cascade.displayBlocks(stream::IO, blockList::Array{Cascade.Block,1}; sa::String="")` 
    ... group & display the blocks of the cascade with same No. of electrons; this blocks are displayed with the
        minimum and maximum energy of each multiplet. The optional sa::String can be used to display some details
        about the given blocks. nothing is returned.
"""
function displayBlocks(stream::IO, blockList::Array{Cascade.Block,1}; sa::String="")
    #
    nx = 150
    println(stream, "\n* Configuration 'blocks' (multiplets): " * sa * "\n")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "      No.   Configurations                                                                       " *
                    "      No. CSF  ",
                    "      Range of total energies " * TableStrings.inUnits("energy") ) 
    println(stream, "  ", TableStrings.hLine(nx))
    i = 0
    for  block  in  blockList
        i = i + 1;    
        sa = "   " * TableStrings.flushright( 6, string(i); na=2)
        sb = " ";         for conf  in blockList[i].confs   sb = sb * string(conf, true) * ", "    end
        en = Float64[];   for level in  block.multiplet.levels    push!(en, level.energy)    end
        minEn = minimum(en);   minEn = Defaults.convertUnits("energy: from atomic", minEn)
        maxEn = maximum(en);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
        sa = sa * TableStrings.flushleft(87, sb[1:end-2]; na=2) 
        sb = "          " * string( length(block.multiplet.levels) )
        sa = sa * sb[end-9:end] * "        "
        sa = sa * TableStrings.flushleft(30, string( round(minEn)) * " ... " * string( round(maxEn)); na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`Cascade.displayDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, rProbabilities::Dict{Subshell,Float64},
                                    aProbabilities::Dict{Tuple{Subshell,Subshell},Float64}, settings::DecayYield.Settings)`  
    ... displays the -- radiative and Auger -- decay probabilities in a neat table and a format close to the geant4 input
        files. However, here the subshells are still displayed in the standard form. If suitable output files are selected,
        these -- radiative and Auger -- decay probabilities are also dumped independently into two ASCII files by using a
        format very similar to GEANT4. Overall this procedure is rather specific in that the level of each outcome 
        is given by a single relativistic configuration and with just a single core hole (subshell). The many-electron rates 
        are then brought back to a single-particle subshell notation. 
        A neat table is printed to the given stream but nothing is returned otherwise.
"""
function displayDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, rProbabilities::Dict{Subshell,Float64},
                                    aProbabilities::Dict{Tuple{Subshell,Subshell},Float64}, settings::DecayYield.Settings)
    
    #
    # Determine the leading configuration and shells of the initial level
    subshList     = Basics.extractRelativisticSubshellList(outcome.level)
    relConfigs    = Basics.extractConfigurations(Basics.RelativisticConfigurations(), outcome.level.basis, outcome.level.J)
    holeSubshells = Basics.extractFromConfiguration(Basics.OpenSubshells(), relConfigs[1]);   holeSubshell = holeSubshells[1]
    @show subshList, relConfigs, holeSubshell
    
    nx = 28
    println(stream, " ")
    println(stream, "  Fluorescence and Auger decay probabilities:")
    println(stream, " ")
    println(stream, "    + Configuration: $(string(relConfigs[1]))   ")
    println(stream, "    + Hole subshell: $holeSubshell ")
    println(stream, "    + Level:         $(outcome.level.index) with symmetry $(LevelSymmetry( outcome.level.J, outcome.level.parity)) ")
    println(stream, "    + $(outcome.NoPhotonLines)           ... number of photon lines. ")
    println(stream, "    + $(outcome.NoAugerLines )           ... number of Auger lines. ")
    if    settings.approach isa Union{Basics.AverageSCA, Basics.SCA}
    println(stream, "    + Approach:  $(settings.approach)  ... all decay probabilities only in Babushkin gauge. ")
    end
    #
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = " "
    sa = sa * TableStrings.center(10, "Shell"; na=2);     
    sa = sa * TableStrings.center(16, "p (radiative)";   na=4)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    for subsh in subshList
        if  haskey(rProbabilities, subsh)
            if  rProbabilities[subsh] == 0.  continue    end
            sa  = "    $subsh     ";    sa = sa * @sprintf("%.6e", rProbabilities[subsh]);     println(stream, sa )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    #  
    sa = " ";       nx = 38
    sa = sa * TableStrings.center(10, "Shell"; na=0) * TableStrings.center(10, "Shell"; na=2);     
    sa = sa * TableStrings.center(16, "p (Auger)";   na=4)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    for subsha in subshList
        for subshb in subshList
            if  haskey(aProbabilities, (subsha, subshb) )
                if  aProbabilities[(subsha, subshb)] == 0.  continue    end
                sa  = "    $subsha    $subshb       ";    sa = sa * @sprintf("%.6e", aProbabilities[(subsha, subshb)]);   println(stream, sa )
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end


"""
`Cascade.dumpGeant4Index(subsh::Subshell)`  
    ... returns the (integer) index of subshell in the Geant4 convention; an i::Int64 is returned
"""
function dumpGeant4Index(subsh::Subshell)
    idx = 0
    if      subsh == Subshell("1s_1/2")       idx = 1
    elseif  subsh == Subshell("2s_1/2")       idx = 3
    elseif  subsh == Subshell("2p_1/2")       idx = 5
    elseif  subsh == Subshell("2p_3/2")       idx = 6
        
    elseif  subsh == Subshell("3s_1/2")       idx = 8
    elseif  subsh == Subshell("3p_1/2")       idx = 10
    elseif  subsh == Subshell("3p_3/2")       idx = 11
    elseif  subsh == Subshell("3d_3/2")       idx = 13
    elseif  subsh == Subshell("3d_5/2")       idx = 14
        
    elseif  subsh == Subshell("4s_1/2")       idx = 16
    elseif  subsh == Subshell("4p_1/2")       idx = 18
    elseif  subsh == Subshell("4p_3/2")       idx = 19
    elseif  subsh == Subshell("4d_3/2")       idx = 21
    elseif  subsh == Subshell("4d_5/2")       idx = 22
    elseif  subsh == Subshell("4f_5/2")       idx = 24
    elseif  subsh == Subshell("4f_7/2")       idx = 25
        
    elseif  subsh == Subshell("5s_1/2")       idx = 27
    elseif  subsh == Subshell("5p_1/2")       idx = 29
    elseif  subsh == Subshell("5p_3/2")       idx = 30
    elseif  subsh == Subshell("5d_3/2")       idx = 32
    elseif  subsh == Subshell("5d_5/2")       idx = 33
    elseif  subsh == Subshell("5f_5/2")       idx = 35
    elseif  subsh == Subshell("5f_7/2")       idx = 36
    elseif  subsh == Subshell("5g_7/2")       idx = 38
    elseif  subsh == Subshell("5g_9/2")       idx = 39
        
    elseif  subsh == Subshell("6s_1/2")       idx = 41
    elseif  subsh == Subshell("6p_1/2")       idx = 43
    elseif  subsh == Subshell("6p_3/2")       idx = 44
    elseif  subsh == Subshell("6d_3/2")       idx = 46
    elseif  subsh == Subshell("6d_5/2")       idx = 47
    elseif  subsh == Subshell("6f_5/2")       idx = 49
    elseif  subsh == Subshell("6f_7/2")       idx = 50
    elseif  subsh == Subshell("6g_7/2")       idx = 52
    elseif  subsh == Subshell("6g_9/2")       idx = 53
    elseif  subsh == Subshell("6h_9/2")       idx = 55
    elseif  subsh == Subshell("6h_11/2")      idx = 56
        
    elseif  subsh == Subshell("7s_1/2")       idx = 58
    elseif  subsh == Subshell("7p_1/2")       idx = 60
    elseif  subsh == Subshell("7p_3/2")       idx = 61
    else    error("stop a; subsh = $subsh")
    end
    
    return( idx )
end


"""
`Cascade.dumpDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, subshEnergies::Dict{Subshell,Float64},
                                rProbabilities::Dict{Subshell,Float64},settings::DecayYield.Settings)`  
    ... dumps the radiative decay probabilities to a selected data file in geant4 form; this is caused by the 
        geant4 boolean of the DecayYield.Settings; this procedure is rather specific in that a single relativistic 
        configuration with a single core hole (subshell) is assumed and that the many-electron rates have been brought back 
        to a single-particle subshell notation. A neat table is printed to the given stream but nothing is returned otherwise.
"""
function dumpDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, subshEnergies::Dict{Subshell,Float64},
                                rProbabilities::Dict{Subshell,Float64},settings::DecayYield.Settings)
    #
    # Determine the leading configuration and shells of the initial level
    subshList     = Basics.extractRelativisticSubshellList(outcome.level)
    relConfigs    = Basics.extractConfigurations(Basics.RelativisticConfigurations(), outcome.level.basis, outcome.level.J)
    holeSubshells = Basics.extractFromConfiguration(Basics.OpenSubshells(), relConfigs[1]);   holeSubshell = holeSubshells[1]
    
    # Dump probabilities for the given holeSubshell to the data file
    idx = Cascade.dumpGeant4Index(holeSubshell)
    sa  = string(idx) * "              "
    sa  = sa[1:13] * sa[1:13] * sa[1:13] * "        " * string(holeSubshell);      println(stream, sa)
    for subsh in subshList
        if  haskey(rProbabilities, subsh)
            if  subsh == holeSubshell        continue    end
            if  rProbabilities[subsh] == 0.  continue    end
            idy = Cascade.dumpGeant4Index(subsh)
            sa  = "$idy            ";    sa = sa[1:13] * @sprintf("%.6e", rProbabilities[subsh])
            en  = subshEnergies[holeSubshell] - subshEnergies[subsh]
            en  = Defaults.convertUnits("energy: from atomic to eV", en) * 1.0e-6 ## energies in MeV
            sa  = sa * " " * @sprintf("%.6e", abs(en));     println(stream, sa)
        end
    end
    sa  = string(-1) * "              "
    sa  = sa[1:13] * sa[1:13] * sa[1:13];      println(stream, sa)

    return( nothing )
end


"""
`Cascade.dumpDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, subshEnergies::Dict{Subshell,Float64},
                                aProbabilities::Dict{Tuple{Subshell,Subshell},Float64},settings::DecayYield.Settings)`  
    ... dumps the Auger decay probabilities to a selected data file in geant4 form; this is caused by the 
        geant4 boolean of the DecayYield.Settings; this procedure is rather specific in that a single relativistic 
        configuration with a single core hole (subshell) is assumed and that the many-electron rates have been brought back 
        to a single-particle subshell notation. A neat table is printed to the given stream but nothing is returned otherwise.
"""
function dumpDecayProbabilities(stream::IO, outcome::DecayYield.Outcome, subshEnergies::Dict{Subshell,Float64},
                                aProbabilities::Dict{Tuple{Subshell,Subshell},Float64},settings::DecayYield.Settings)
    #
    # Determine the leading configuration and shells of the initial level
    subshList     = Basics.extractRelativisticSubshellList(outcome.level)
    relConfigs    = Basics.extractConfigurations(Basics.RelativisticConfigurations(), outcome.level.basis, outcome.level.J)
    holeSubshells = Basics.extractFromConfiguration(Basics.OpenSubshells(), relConfigs[1]);   holeSubshell = holeSubshells[1]
    
    # Dump probabilities for the given holeSubshell to the data file
    idx = Cascade.dumpGeant4Index(holeSubshell)
    sa  = string(idx) * "              "
    sa  = sa[1:13] * sa[1:13] * sa[1:13] * sa[1:13] * "        " * string(holeSubshell);      println(stream, sa)
    for subsha in subshList
        for subshb in subshList
            if  haskey(aProbabilities, (subsha, subshb) )
                if  subsha == holeSubshell  ||  subshb == holeSubshell   continue    end
                if  aProbabilities[(subsha, subshb)] == 0.               continue    end
                ida = Cascade.dumpGeant4Index(subsha);    idb = Cascade.dumpGeant4Index(subshb)
                sa  = "$ida            $idb            ";    sa = sa[1:26] * @sprintf("%.6e", aProbabilities[(subsha, subshb)])
                en  = subshEnergies[holeSubshell] - subshEnergies[subsha] - subshEnergies[subshb]
                en  = Defaults.convertUnits("energy: from atomic to eV", en) * 1.0e-6 ## energies in MeV
                sa  = sa * " " * @sprintf("%.6e", abs(en));     println(stream, sa)
            end
        end
    end
    sa  = string(-1) * "              "
    sa  = sa[1:13] * sa[1:13] * sa[1:13] * sa[1:13];      println(stream, sa)

    return( nothing )
end


"""
`Cascade.displayLevels(stream::IO, multiplets::Array{Multiplet,1}; sa::String="")`  
    ... display on stream the initial configurations as well as the calculated levels for all initial multiplets.
"""
function displayLevels(stream::IO, multiplets::Array{Multiplet,1}; sa::String="")
    nx = 44
    println(stream, " ")
    println(stream, "* Configurations and levels for all given " * sa * "multiplets of the cascade, relative to the lowest:")
    for  multiplet  in multiplets
        println(stream, "  ")
        confList = Basics.extractConfigurations(Basics.FromBasis(), multiplet.levels[1].basis)
         for  conf in confList
            println(stream, "  $conf")
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "    Level  J Parity          Energy " * TableStrings.inUnits("energy") ) 
        println(stream, "  ", TableStrings.hLine(nx))
        for  i = 1:length(multiplet.levels)
            lev = multiplet.levels[i]
            en  = lev.energy - multiplet.levels[1].energy;    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = "   "  * TableStrings.level(i) * "     " * string(LevelSymmetry(lev.J, lev.parity)) * "     "
            @printf(stream, "%s %.15e %s", sc, en_requested, "\n")
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    return( nothing )
end


"""
`Cascade.displaySteps(stream::IO, steps::Array{Cascade.Step,1}; sa::String="")` 
    ... displays all predefined steps in a neat table and supports to delete individual steps from the list.
        sa::String can be used to display details about the given steps
"""
function displaySteps(stream::IO, steps::Array{Cascade.Step,1}; sa::String="")
    nx = 170
    println(stream, " ")
    println(stream, "* Steps that are defined for the current " * sa * "cascade due to the given approach:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center( 9, "Step-No"; na=2)
    sa = sa * TableStrings.flushleft(11, "Process"; na=1)
    sa = sa * TableStrings.flushleft(55, "Initial:  No CSF, configuration(s)"; na=4)
    sa = sa * TableStrings.flushleft(55, "Final:  No CSF, configuration(s)"; na=4)
    sa = sa * TableStrings.flushleft(40, "Energies from ... to in " * TableStrings.inUnits("energy"); na=4)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  i = 1:length(steps)
        sa = " " * TableStrings.flushright( 7, string(i); na=5)
        sa = sa  * TableStrings.flushleft( 11, string(steps[i].process); na=1)
        sb = "";   for conf in steps[i].initialConfigs   sb = sb * string(conf,true) * ", "    end
        sa = sa  * TableStrings.flushright( 5, string( length(steps[i].initialMultiplet.levels[1].basis.csfs) )*", "; na=0) 
        sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
        sb = "";   for conf in steps[i].finalConfigs     sb = sb * string(conf,true) * ", "    end
        sa = sa  * TableStrings.flushright( 5, string( length(steps[i].finalMultiplet.levels[1].basis.csfs) )*", "; na=0) 
        sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
        minEn = 1000.;   maxEn = -1000.;
        for  p = 1:length(steps[i].initialMultiplet.levels),  q = 1:length(steps[i].finalMultiplet.levels)
            minEn = min(minEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
            maxEn = max(maxEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
        end
        minEn = Defaults.convertUnits("energy: from atomic", minEn);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
        sa = sa * string( round(minEn)) * " ... " * string( round(maxEn))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
end


"""
`Cascade.generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}; 
                        sa::String="", printout::Bool=true)`  
    ... generate all block::Cascade.Block's, that need to be computed for this cascade, and compute also the corresponding multiplets.
        The different cascade approches enables one to realized follow different strategies how these block are selected and computed. 
        A blockList::Array{Cascade.Block,1} is returned.
"""
function generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}; 
                        sa::String="", printout::Bool=true)
    blockList = Cascade.Block[]
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    ## The approach is no longer interrogated by name here. It is resolved into its conditions, and those
    ## conditions decide what happens -- which is also what the six lines printed by Cascade.displayApproach
    ## report. Announcing and doing therefore come from ONE source and cannot drift apart, as they had: this
    ## block used to print "all blocks (multiplets) are generated from single-CSF levels and without any
    ## configuration mixing" while the code below was running full configuration interaction.
    cd = Cascade.conditions(comp.approach)
    if  printout
        println("\n* Generate blocks " * sa)
        Cascade.displayApproach(stdout, comp.approach, comp.asfSettings)
        if  printSummary
            println(iostream, "\n* Generate blocks " * sa)
            Cascade.displayApproach(iostream, comp.approach, comp.asfSettings)
        end
    end
    #
    ## ChargeStateOrbitals needs one self-consistent field per charge state, prepared once for all blocks.
    ## GlobalOrbitals uses the initial multiplet's own orbitals, which the caller has already passed in and
    ## which are precisely "the mean field of the initial ion" that the approach asks for. MultipletOrbitals
    ## prepares nothing here: every block runs its own self-consistent field in the loop below.
    orbitalSets = Dict{Int64, Dict{Subshell, Orbital}}()
    if  cd.bound == Cascade.ChargeStateOrbitals()
        orbitalSets = Cascade.generateBoundOrbitals(comp.approach, comp, confs; printout=printout)
    end
    #
    for  confa  in confs
        print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
        if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
        if      cd.bound == Cascade.MultipletOrbitals()
            multiplet = SelfConsistent.performSCF([confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
        elseif  cd.bound == Cascade.GlobalOrbitals()
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa], initalOrbitals, comp.nuclearModel, comp.grid,
                                                                Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings); printout=false)
        else
            multiplet = Hamiltonian.performCIwithFrozenOrbitals([confa], orbitalSets[confa.NoElectrons], comp.nuclearModel, comp.grid,
                                                                Cascade.asfSettingsForApproach(comp.approach, comp.asfSettings); printout=false)
        end

        # Shift the total energies of all levels if requested for the StepwiseDecayScheme
        if  typeof(comp.scheme) == StepwiseDecayScheme   &&   haskey(comp.scheme.chargeStateShifts, confa.NoElectrons)
            energyShift = comp.scheme.chargeStateShifts[confa.NoElectrons]
            multiplet   = Basics.shiftTotalEnergies(multiplet, energyShift)
            print("shift all levels by $energyShift [a.u.] ... ")
        end
        push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
        println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
    end

    return( blockList )
end


"""
`Cascade.generateConfigurationList(multiplets::Array{Multiplet,1}, further::Int64, NoShake::Int64)`  
    ... generates all possible (decay) configurations with up to further holes and with NoShake displacements with regard
        to the given multiplets. First, all configuratons are generated for which the hole is either moved 'outwards' or 
        is moved and a second 'outer' hole is created; this step is repated further + 2 times to make sure that all relevant
        configurations are met. From the generated list, however, only those configurations are kept eventually with 
        up to further holes, when compared to the configurations of the given multiplets. A confList::Array{Configuration,1} 
        is returned.
"""
function generateConfigurationList(multiplets::Array{Multiplet,1}, further::Int64, NoShake::Int64)
    # Determine all (different) configurations from multiplets
    confList = Configuration[]
    for mp  in  multiplets   
        cfList = Basics.extractConfigurations(Basics.FromBasis(), mp.levels[1].basis)
        for  cf in cfList   if  cf in confList   nothing   else   push!(confList, cf)      end      end
    end
    cList = copy(confList);   initialNoElectrons = multiplets[1].levels[1].basis.NoElectrons
    # First, move and generate new 'outer' hole without displacements
    for  fur = 1:further+1
        newConfList = Configuration[]
        for conf  in cList
            wa = Basics.generateConfigurations(Basics.ForPhotoEmission(), [conf]);   append!(newConfList, wa)
            wa = Basics.generateConfigurations(Basics.ForAutoIonization(), [conf]);  append!(newConfList, wa)
        end
        newConfList = unique(newConfList)
        cList = newConfList
        append!(confList, newConfList)
    end
    # Make sure that only configurations with up to further holes are returned
    newConfList = Configuration[]
    for   conf in confList   
        if  conf.NoElectrons + further >= initialNoElectrons   push!(newConfList, conf)    end
    end
    # Add further shake-displacements if appropriate
    newConfList = unique(newConfList)
    return( newConfList )
end


"""
`Cascade.modifySteps(stepList::Array{Cascade.Step,1})` 
    ... allows the user to modify the steps, for instance, by deleting selected steps of the cascade or by modifying the settings of
        one or several steps. A newStepList::Array{Cascade.Step,1} for which the transition data are eventually computed.
"""
function modifySteps(stepList::Array{Cascade.Step,1})
    #
    newStepList = Cascade.Step[]
    #
    println("\n* Here, modify the individual steps explicitly in the code, if needed, ...... and just do it !!")
    # 
    #  Delete individual steps from stepList
    #  if  i in [1,2,5, ...] modify the particular settings, etc.
    for  i = 1:length(stepList)
        step = stepList[i]
        #
        if  i in []
            println("  Modify step $i :")
            newStep = Cascade.Step(step.process, step.settings, step.initialConfs, step.finalConfs, step.initialMultiplet, step.initialMultiplet)
            push!(newStepList, newStep)
        else
            push!(newStepList, step)
        end
    end
    #
    # wa = [1,2,3]
    # delete from list
    #
    println("\n  A total of $(length(newStepList)) steps are still defined in the cascade.")
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   println(iostream, "\n* A total of $(length(newStepList)) steps are still defined in the cascade.")    end      
    
    return( newStepList )
end
