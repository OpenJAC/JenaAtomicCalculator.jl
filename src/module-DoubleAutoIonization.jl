
"""
`module  JAC.DoubleAutoIonization`  
... a submodel of JAC that contains all methods for computing double Auger and autoionization amplitudes and rates.
"""
module DoubleAutoIonization


using Printf, ..AngularMomentum, ..AtomicState, ..AutoIonization, ..Basics, ..Bsplines, ..ManyElectron, ..Nuclear, ..Radial,
      ..TableStrings


"""
`struct  DoubleAutoIonization.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the settings in estimating double-Auger and autoionization rates.

    + green                   ::Array{GreenChannel,1}   ... Precalculated and user-specified Green function of the ion.
    + NoEnergySharings        ::Int64                   ... Number of energy sharings that are used in the computations for each line.
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before their evaluation.
    + minAugerEnergy          ::Float64                 ... Minimum energy of free (Auger) electrons to be included.
    + maxAugerEnergy          ::Float64                 ... Maximum energy of free (Auger) electrons to be included.
    + maxKappa                ::Int64                   ... Maximum kappa value of partial waves to be included.
    + operator                ::AbstractEeInteraction   ... Auger operator that is to be used for evaluating the Auger amplitudes: 
                                                            allowed values are: CoulombInteraction(), BreitInteraction(), ...
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractProcessSettings
    green                     ::Array{GreenChannel,1}
    NoEnergySharings          ::Int64 
    printBefore               ::Bool 
    minAugerEnergy            ::Float64
    maxAugerEnergy            ::Float64
    maxKappa                  ::Int64  
    operator                  ::AbstractEeInteraction
    lineSelection             ::LineSelection 
end 


"""
`DoubleAutoIonization.Settings()`  ... constructor for the default values of DoubleAutoIonization line computations
"""
function Settings()
    Settings(GreenChannel[], 0, false, 0., 10e5, 100, CoulombInteraction(), LineSelection())
end


"""
`DoubleAutoIonization.Settings(set::DoubleAutoIonization.Settings;`

    green..., NoEnergySharings..., printBefore..., minAugerEnergy..., maxAugerEnergy..., maxKappa..., operator...,
    lineSelection...)

    ... constructor for modifying the given DoubleAutoIonization.Settings by 'overwriting' the previously selected parameters.
        A settings::DoubleAutoIonization.Settings is returned.
"""
function Settings(set::DoubleAutoIonization.Settings;
                        green::Union{Nothing,Array{GreenChannel,1}} = nothing,
                        NoEnergySharings::Union{Nothing,Int64} = nothing,
                        printBefore::Union{Nothing,Bool} = nothing,
                        minAugerEnergy::Union{Nothing,Float64} = nothing,
                        maxAugerEnergy::Union{Nothing,Float64} = nothing,
                        maxKappa::Union{Nothing,Int64} = nothing,
                        operator::Union{Nothing,AbstractEeInteraction} = nothing,
                        lineSelection::Union{Nothing,LineSelection} = nothing )

    if  isnothing(green)              greenx            = set.green              else   greenx            = green              end
    if  isnothing(NoEnergySharings)   NoEnergySharingsx = set.NoEnergySharings   else   NoEnergySharingsx = NoEnergySharings   end
    if  isnothing(printBefore)        printBeforex      = set.printBefore        else   printBeforex      = printBefore        end
    if  isnothing(minAugerEnergy)     minAugerEnergyx   = set.minAugerEnergy     else   minAugerEnergyx   = minAugerEnergy     end
    if  isnothing(maxAugerEnergy)     maxAugerEnergyx   = set.maxAugerEnergy     else   maxAugerEnergyx   = maxAugerEnergy     end
    if  isnothing(maxKappa)           maxKappax         = set.maxKappa           else   maxKappax         = maxKappa           end
    if  isnothing(operator)           operatorx         = set.operator           else   operatorx         = operator           end
    if  isnothing(lineSelection)      lineSelectionx    = set.lineSelection      else   lineSelectionx    = lineSelection      end

    Settings( greenx, NoEnergySharingsx, printBeforex, minAugerEnergyx, maxAugerEnergyx, maxKappax, operatorx,
              lineSelectionx )
end


# `Base.show(io::IO, settings::DoubleAutoIonization.Settings)`  ... prepares a proper printout of the variable settings::DoubleAutoIonization.Settings.
function Base.show(io::IO, settings::DoubleAutoIonization.Settings) 
    println(io, "green:                         (settings.green)  ")
    println(io, "NoEnergySharings:             $(settings.NoEnergySharings)  ")
    println(io, "printBefore:                  $(settings.printBefore)  ")
    println(io, "minAugerEnergy:               $(settings.minAugerEnergy)  ")
    println(io, "maxAugerEnergy:               $(settings.maxAugerEnergy)  ")
    println(io, "maxKappa:                     $(settings.maxKappa)  ")
    println(io, "operator:                     $(settings.operator)  ")
    println(io, "lineSelection:                $(settings.lineSelection)  ")
end


"""
`struct  DoubleAutoIonization.ReducedChannel`  
    ... defines a type for a DoubleAutoIonization (reduced) channel to help characterize a scattering (continuum) state of many 
        electron-states with two free electrons.

    + energy1        ::Float64              ... energy of the partial wave_1
    + kappa1         ::Int64                ... partial-wave of the free electron_1
    + phase1         ::Float64              ... phase of the partial wave_1
    + xSymmetry      ::LevelSymmetry        ... angular momentum and parity after coupling electron1.
    + energy2        ::Float64              ... energy of the partial wave_2
    + kappa2         ::Int64                ... partial-wave of the free electron_2
    + phase2         ::Float64              ... phase of the partial wave_2
    + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
"""
struct  ReducedChannel
    energy1          ::Float64 
    kappa1           ::Int64 
    phase1           ::Float64
    xSymmetry        ::LevelSymmetry 
    energy2          ::Float64 
    kappa2           ::Int64   
    phase2           ::Float64 
    amplitude        ::Complex{Float64}
end


# `Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel)`  
#       ... prepares a proper printout of the variable cannel::DoubleAutoIonization.ReducedChannel.
function Base.show(io::IO, channel::DoubleAutoIonization.ReducedChannel) 
    ## sa = "reduced DA channel for J^P=(channel.tSymmetry)  with " * 
    ##      "kappa1=$(channel.kappa1), energy1=$(channel.energy1), phase1=$(channel.phase1), " *
    ##      "kappa2=$(channel.kappa2), energy2=$(channel.energy2), phase2=$(channel.phase2), amp=$(channel.amplitude) "
    ## println(io, sa)
    println(io, "multipole:                    $(channel.multipole)  ")
    println(io, "gauge:                        $(channel.gauge)  ")
    println(io, "omega:                        $(channel.omega)  ")
    println(io, "energy1:                      $(channel.energy1)  ")
    println(io, "kappa1:                       $(channel.kappa1)  ")
    println(io, "phase1:                       $(channel.phase1)  ")
    println(io, "xSymmetry:                    $(channel.xSymmetry)  ")
    println(io, "energy2:                      $(channel.energy2)  ")
    println(io, "kappa2:                       $(channel.kappa2)  ")
    println(io, "phase2:                       $(channel.phase2)  ")
    println(io, "amplitude:                    $(channel.amplitude)  ")
end


"""
`struct  DoubleAutoIonization.Sharing`  
    ... defines a type for a DoubleAutoIonization sharing to help characterize energy sharing between the two emitted electrons.

    + epsilon1         ::Float64         ... Energy of (free) electron 1.
    + epsilon2         ::Float64         ... Energy of (free) electron 2.
    + weight           ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
    + differentialRate ::Float64         ... differential cross section of this energy sharing.
    + channels         ::Array{DoubleAutoIonization.ReducedChannel,1}  ... List of DoubleAutoIonization (reduced) channels of this line.
"""
struct  Sharing
    epsilon1           ::Float64
    epsilon2           ::Float64
    weight             ::Float64 
    differentialRate   ::Float64
    channels           ::Array{DoubleAutoIonization.ReducedChannel,1}
end


# `Base.show(io::IO, sharing::DoubleAutoIonization.Sharing)`  ... prepares a proper printout of the variable sharing::DoubleAutoIonization.Sharing.
function Base.show(io::IO, sharing::DoubleAutoIonization.Sharing) 
    println(io, "epsilon1:               $(sharing.epsilon1)  ")
    println(io, "epsilon2:               $(sharing.epsilon2)  ")
    println(io, "differentialCs:         $(sharing.differentialCs)  ")
    println(io, "channels:               $(sharing.channels)  ")
end


"""
`struct  DoubleAutoIonization.Line`  ... defines a type for a double Auger line that includes sharings and their reduced amplitudes.

    + initialLevel   ::Level            ... initial-(state) level
    + finalLevel     ::Level            ... final-(state) level
    + totalRate      ::Float64          ... Total rate of this line.
    + sharings       ::Array{DoubleAutoIonization.Sharing,1}  ... List of DoubleAutoIonization sharings of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    totalRate        ::Float64
    sharings         ::Array{DoubleAutoIonization.Sharing,1}
end 


# `Base.show(io::IO, line::DoubleAutoIonization.Line)`  ... prepares a proper printout of the variable line::DoubleAutoIonization.Line.
    function Base.show(io::IO, line::DoubleAutoIonization.Line) 
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "totalRate:              $(line.totalRate)  ")
    println(io, "sharings:               $(line.sharings)  ")
end


"""
`DoubleAutoIonization.amplitude(kind::String, channel::DoubleAutoIonization.ReducedChannel, 
                                continuumLevel::Level, green::Array{GreenChannel,1}, initialLevel::Level, grid::Radial.Grid)`  
    ... to compute the kind = (autoionization) amplitude  
        <((alpha_f J_f, epsilon1 kappa1) J_x, epsilon2 kappa2) J_i || A || alpha_i J_i>  with second-order operator
    
            A  =  V^(e-e)  V^(e-e)  =  Sum_n  (V^(e-e)  |n >< n|  V^(e-e)) / (E_i - E_n)
                
        due to the electron-electron interaction V^(e-e) for the given final and initial level. Of course, this amplitude 
        depends on the partial wave of the outgoing electrons. A value::ComplexF64 is returned.
        A present, the CoulombInteraction() is hard-coded for the electron-electron interaction.
"""
function amplitude(kind::String, channel::DoubleAutoIonization.ReducedChannel, 
                    continuumLevel::Level, green::Array{GreenChannel,1}, initialLevel::Level, grid::Radial.Grid)
    if      kind in [ "double-Auger"]
    #--------------------------------
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);      symc = LevelSymmetry(continuumLevel.J, continuumLevel.parity)
        eni  = initialLevel.energy
        amplitude = ComplexF64(0.);
        # Compute the terms for <((alpha_f J_f, epsilon1 kappa1) J_x, epsilon2 kappa2) J_i || A || alpha_i J_i>
        for greenChannel in green
            if  greenChannel.symmetry == symi == symc
                println(">>> Apply Green function channel with symmetry $(greenChannel.symmetry)")
                for  (nlev, level)  in  enumerate(greenChannel. gMultiplet.levels)
                    ## if  nlev > 2  println(">>> Skip nlev > 2 ... need to be corrected !!!!");    break     end
                    am_ni     = AutoIonization.amplitude(CoulombInteraction(), channel.kappa1, channel.phase1, level, initialLevel, grid; printout=false)
                    am_fn     = AutoIonization.amplitude(CoulombInteraction(), channel.kappa2, channel.phase2, continuumLevel, level, grid; printout=false)
                    amplitude = amplitude  +  am_fn * am_ni / (eni - level.energy)
                    println(">> DAI for level $nlev  with  am_fn(A) = $am_fn  and  am_ni(A) = $am_ni")
                end 
            else    ## println(">>> Skip Green function channel with symmetry $(greenChannel.symmetry)")
            end
        end
        #
    else    error("stop b")
    end
    
    return( amplitude )
end


"""
`DoubleAutoIonization.computeAmplitudesProperties(line::DoubleAutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                    settings::DoubleAutoIonization.Settings)`  
    ... to compute all amplitudes and properties of the given line; a line::DoubleAutoIonization.Line is returned 
        for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(line::DoubleAutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::DoubleAutoIonization.Settings)
    newSharings  = DoubleAutoIonization.Sharing[]
    contSettings = Continuum.Settings(false, grid.NoPoints-50);    trate = 0.
    # The nuclear potential and the B-spline primitives depend only on (nm, grid) and are constant here, while the loop
    # below generates TWO continuum orbitals per channel and per sharing; building them once keeps that work out of it.
    nucPot       = Nuclear.nuclearPotential(nm, grid)
    primitives   = Bsplines.generatePrimitives(grid)
    symi         = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    #
    for sharing in line.sharings
        newChannels = DoubleAutoIonization.ReducedChannel[];    drate = 0.
        for (ich, ch)  in  enumerate(sharing.channels)
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
            # Generate two continuum orbitals as associated with this sharing; these orbitals move each in the potential of the final level
            cOrbital1, phase1 = Continuum.generateOrbitalForLevel(ch.energy1, Subshell(101, ch.kappa1), newfLevel, nm, grid, contSettings;
                                                                   nuclearPot=nucPot, primitives=primitives)
            cOrbital2, phase2 = Continuum.generateOrbitalForLevel(ch.energy2, Subshell(102, ch.kappa2), newfLevel, nm, grid, contSettings;
                                                                   nuclearPot=nucPot, primitives=primitives)
            # Expand these continuum orbitals in the orbital set as given by the Green function
            greenOrbs  = settings.green[1].gMultiplet.levels[1].basis.orbitals
            bOrbitals1 = Basics.expandOrbital(cOrbital1, greenOrbs, 3.0e-5, grid, printout=true)
            bOrbitals2 = Basics.expandOrbital(cOrbital2, greenOrbs, 3.0e-5, grid, printout=true)
            # Loop through all sub-orbital contribution in the expansion of the continuum orbitals
            newChannel = DoubleAutoIonization.ReducedChannel(ch.energy1, ch.kappa1, phase1, ch.xSymmetry, 
                                                                ch.energy2, ch.kappa2, phase2, ComplexF64(0.) )
            amplitude  = ComplexF64(0.)
            for  bOrb1  in  bOrbitals1
                if   bOrb1.subshell  in  newfLevel.basis.subshells        continue    end
                for  bOrb2  in  bOrbitals2
                    if   bOrb2.subshell  in  newfLevel.basis.subshells    continue    end
                    if   bOrb1.subshell < bOrb2.subshell
                        # Only  subshell_1 < subshell_2  need to be considered since the 'reversed' case is part of another channel 
                        # Perhaps, one needs to multiply with 2 here (??)
                        symf = LevelSymmetry(newfLevel.J, newfLevel.parity) 
                        newcLevel  = Basics.generateLevelWithExtraTwoElectrons(bOrb1, ch.xSymmetry, bOrb2, symi, newfLevel)
                        amplitude  = amplitude + 
                                        DoubleAutoIonization.amplitude("double-Auger", newChannel, newcLevel, settings.green, newiLevel, grid)
                    elseif   bOrb1.subshell  ==  bOrb2.subshell     error("stop a")
                    end
                end
            end
            push!( newChannels, ReducedChannel(ch.energy1, ch.kappa1, phase1, ch.xSymmetry,
                                                ch.energy2, ch.kappa2, phase2, amplitude) )
            drate = drate  +  2pi * abs(amplitude)^2
        end
        # Calculate the differential decay rates
        push!( newSharings, DoubleAutoIonization.Sharing( sharing.epsilon1, sharing.epsilon2, sharing.weight, drate, newChannels) )
        trate = trate + sharing.weight * drate
    end
    # Add the total decay rate
    newLine = DoubleAutoIonization.Line( line.initialLevel, line.finalLevel, trate, newSharings)
    return( newLine )

end


"""
`DoubleAutoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::DoubleAutoIonization.Settings; output=true)`  
    ... to compute the double-Auger transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{DoubleAutoIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                        settings::DoubleAutoIonization.Settings; output=true)
    println("")
    printstyled("DoubleAutoIonization.computeLines(): The computation of double-Auger rates starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    lines = DoubleAutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    DoubleAutoIonization.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = DoubleAutoIonization.Line[]
    for  line in lines
        newLine = DoubleAutoIonization.computeAmplitudesProperties(line, nm, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    DoubleAutoIonization.displayTotalRates(stdout, newLines, settings)
    DoubleAutoIonization.displayDifferentialRates(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   DoubleAutoIonization.displayTotalRates(iostream, newLines, settings)
                        DoubleAutoIonization.displayDifferentialRates(iostream, newLines, settings)     end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`DoubleAutoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::DoubleAutoIonization.Settings)` 
    ... to determine a list of DoubleAutoIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{DoubleAutoIonization.Line,1} is 
        returned. Apart from the level specification and sharing, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::DoubleAutoIonization.Settings)
    lines = DoubleAutoIonization.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy     = iLevel.energy - fLevel.energy
                if  energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end  
                sharings = DoubleAutoIonization.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                push!( lines, DoubleAutoIonization.Line(iLevel, fLevel, 0., sharings) )
            end
        end
    end
    return( lines )
end


"""
`DoubleAutoIonization.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::DoubleAutoIonization.Settings)`  
    ... to determine a list of DoubleAutoIonization Sharing's and Channel's for a transitions from the initial to final level and by taking into 
        account the particular settings of for this computation; an Array{DoubleAutoIonization.Sharing,1} is returned.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::DoubleAutoIonization.Settings)
    sharings = DoubleAutoIonization.Sharing[];    eSharings = Basics.determineEnergySharings(energy, settings.NoEnergySharings) 
    for  es in eSharings
        epsilon1  = es[1];    epsilon2 = es[2];    weight = es[3] 
        channels  = DoubleAutoIonization.ReducedChannel[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        couplings = AngularMomentum.allowedDoubleKappaCouplingSequence(symf, symi, settings.maxKappa)
        for  coupling  in  couplings
            push!(channels, ReducedChannel(epsilon1, coupling[1], 0., coupling[2], epsilon2, coupling[3], 0., Complex(0.)) ) 
        end
        push!(sharings, DoubleAutoIonization.Sharing(epsilon1, epsilon2, weight, 0., channels) )
    end
    return( sharings )  
end


"""
`DoubleAutoIonization.displayDifferentialRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)`  
    ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayDifferentialRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)
    #
    # First, print lines and sharings
    nx = 110
    println(stream, " ")
    println(stream, "  Energy-differential rates of selected double-Auger lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=3);              
    sb = sb * TableStrings.flushleft(38, "  i -- f        epsilon_1     epsilon_2"; na=3)
    sa = sa * TableStrings.center(14, "Weight"; na=0);                            sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(16, "Diff. rate"; na=3)      
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        #
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1))           * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2))           * "    "
            sb = sb * @sprintf("%.4e",                                              sharing.weight)              * "     "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic",   sharing.differentialRate))   * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`DoubleAutoIonization.displayLines(stream::IO, lines::Array{DoubleAutoIonization.Line,1})`  
    ... to display a list of lines, sharings and channels that have been selected due to the prior settings. A neat table 
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{DoubleAutoIonization.Line,1})
    #
    # First, print lines and sharings
    nx = 82
    println(stream, " ")
    println(stream, "  Selected double-Auger lines & sharings:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(54, "Energies (all in )" * TableStrings.inUnits("energy") * ")"; na=5);              
    sb = sb * TableStrings.flushleft(54, "   i -- f       epsilon_1     epsilon_2"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        #
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon1)) * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.epsilon2)) * "    "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Second, print lines and channles
    nx = 123
    println(stream, " ")
    println(stream, "  Selected photo-double ionization lines & channels:   ... channels are shown only for the first sharing")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(77, "Channels (all energies in " * TableStrings.inUnits("energy") * ")" ; na=5);              
    sb = sb * TableStrings.flushleft(77, "(epsilon, kappa)_1   -->    J^P_x  -->   (epsilon, kappa)_2   -->  J^P_t = J^P_i"; na=5)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        for (ic, ch) in enumerate(line.sharings[1].channels)
            if  ic == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            skappa1 = "  " * string(ch.kappa1);    skappa2 = "  " * string(ch.kappa2);    sxsym = "      " * string(ch.xSymmetry)
            sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy1)) * ", " * skappa1[end-2:end] * ")    -->"
            sb = sb * sxsym[end-8:end] * "  -->   "
            sb = sb * " (" * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", ch.energy2)) * ", " * skappa2[end-2:end] * ")    -->     "
            sb = sb * string(isym) 
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`DoubleAutoIonization.displayTotalRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)`  
    ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayTotalRates(stream::IO, lines::Array{DoubleAutoIonization.Line,1}, settings::DoubleAutoIonization.Settings)
    #
    # First, print lines and sharings
    nx = 70
    println(stream, " ")
    println(stream, "  Total rates of selected double-Auger lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(21)
    sa = sa * TableStrings.center(12, "i--Energy--f"; na=3)               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(16, "Total rate"; na=4)      
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))        * "       "
        sb = sa * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate))  * "     "
        println(stream,  sb )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
