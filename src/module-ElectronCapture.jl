
"""
`module  JAC.ElectronCapture`  
... a submodel of JAC that contains all methods for computing electron-capture properties between some initial and final-state 
    multiplets.
"""
module ElectronCapture


## ..AutoIonization is imported because the capture amplitude IS the complex conjugate of the Auger amplitude; the
## module is included immediately before this one in JenaAtomicCalculator.jl, so the dependency is well ordered.
using  Printf, ..AngularMomentum, ..AutoIonization, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength,
                ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

"""
`struct  ElectronCapture.Settings`  ... defines a type for the details and parameters of computing electron-capture lines.

    + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
    + minCaptureEnergy        ::Float64                      ... Minimum energy of free (Auger) electrons to be captured.
    + maxCaptureEnergy        ::Float64                      ... Maximum energy of free (Auger) electrons to be captured.
    + maxKappa                ::Int64                        ... Maximum kappa value of partial waves to be included.
    + operator                ::AbstractEeInteraction        ... Auger/capture operator that is to be used for evaluating the capture amplitudes: 
                                                                    allowed values are: CoulombInteraction(), BreitInteraction(), ...
"""
struct Settings  <:  AbstractProcessSettings
    printBefore               ::Bool 
    lineSelection             ::LineSelection 
    minCaptureEnergy          ::Float64
    maxCaptureEnergy          ::Float64
    maxKappa                  ::Int64
    operator                  ::AbstractEeInteraction 
    calcAlignment             ::Bool
end 


"""
`ElectronCapture.Settings()`  ... constructor for the default values of ElectronCapture line computations
"""
function Settings()
    Settings(false, LineSelection(), 0., 1.0e5, 2, CoulombInteraction(), true)
end


"""
`ElectronCapture.Settings(set::ElectronCapture.Settings;`

        printBefore=..,       minCaptureEnergy=..,       maxCaptureEnergy=..,       maxKappa=..,       operator=..)
                    
    ... constructor for modifying the given ElectronCapture.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ElectronCapture.Settings;    
    printBefore::Union{Nothing,Bool}=nothing,               lineSelection::Union{Nothing,LineSelection}=nothing, 
    minCaptureEnergy::Union{Nothing,Float64}=nothing,       maxCaptureEnergy::Union{Nothing,Float64}=nothing,
    maxKappa::Union{Nothing,Int64}=nothing,                 operator::Union{Nothing,AbstractEeInteraction}=nothing,
    calcAlignment::Union{Nothing,Bool}=nothing)  
    
    if  isnothing(printBefore)        printBeforex      = set.printBefore       else  printBeforex      = printBefore        end 
    if  isnothing(lineSelection)      lineSelectionx    = set.lineSelection     else  lineSelectionx    = lineSelection      end 
    if  isnothing(minCaptureEnergy)   minCaptureEnergyx = set.minCaptureEnergy  else  minCaptureEnergyx = minCaptureEnergy   end 
    if  isnothing(maxCaptureEnergy)   maxCaptureEnergyx = set.maxCaptureEnergy  else  maxCaptureEnergyx = maxCaptureEnergy   end 
    if  isnothing(maxKappa)           maxKappax         = set.maxKappa          else  maxKappax         = maxKappa           end 
    if  isnothing(operator)           operatorx         = set.operator          else  operatorx         = operator           end 
    if  isnothing(calcAlignment)      calcAlignmentx    = set.calcAlignment     else  calcAlignmentx    = calcAlignment      end 

    Settings( printBeforex, lineSelectionx, minCaptureEnergyx, maxCaptureEnergyx, maxKappax, operatorx, calcAlignmentx)
end


# `Base.show(io::IO, settings::ElectronCapture.Settings)`  ... prepares a proper printout of the variable settings::ElectronCapture.Settings.
function Base.show(io::IO, settings::ElectronCapture.Settings) 
    println(io, "printBefore:                   $(settings.printBefore)  ")
    println(io, "lineSelection:                 $(settings.lineSelection)  ")
    println(io, "minCaptureEnergy:              $(settings.minCaptureEnergy)  ")
    println(io, "maxCaptureEnergy:              $(settings.maxCaptureEnergy)  ")
    println(io, "maxKappa:                      $(settings.maxKappa)  ")
    println(io, "operator:                      $(settings.operator)  ")
    println(io, "calcAlignment:                 $(settings.calcAlignment)  ")
end


"""
`struct  Channel`   
    ... defines a type for a ElectronCapture channel to help characterize a scattering (continuum) state of many 
        electron-states with a single free electron.

    + kappa          ::Int64                ... partial-wave of the free electron
    + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + phase          ::Float64              ... phase of the partial wave
    + amplitude      ::Complex{Float64}     ... electron-capture amplitude associated with the given channel.
"""
struct  Channel
    kappa            ::Int64
    symmetry         ::LevelSymmetry
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  Line`  
    ... defines a type for a ElectronCapture line that may include the definition of sublines and their 
        corresponding amplitudes.

    + initialLevel   ::Level           ... initial-(state) level
    + finalLevel     ::Level           ... final-(state) level
    + electronEnergy ::Float64         ... Energy of the (incoming free) electron.
    + totalRate      ::Float64         ... Capture rate P_cap(0 -> d) of this line.
    + alignment      ::Array{Float64,1}  ... Alignment parameters A_20, A_40, ... of the captured (doubly excited)
                                             state; even ranks only, k <= 2J_d, and empty when J_d < 1.
    + channels       ::Array{ElectronCapture.Channel,1}  ... List of ElectronCapture channels of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    totalRate        ::Float64
    alignment        ::Array{Float64,1}
    channels         ::Array{ElectronCapture.Channel,1}
end 


"""
`ElectronCapture.Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)`  
    ... constructor for an ElectronCapture line between a specified initial and final level; a
        line::ElectronCapture.Line with no channels and no alignment is returned.
"""
function Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)
    Line(initialLevel, finalLevel, 0., totalRate, Float64[], ElectronCapture.Channel[])
end


# `Base.show(io::IO, line::ElectronCapture.Line)`  ... prepares a proper printout of the variable line::ElectronCapture.Line.
function Base.show(io::IO, line::ElectronCapture.Line) 
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "electronEnergy:         $(line.electronEnergy)  ")
    println(io, "totalRate:              $(line.totalRate)  ")
    println(io, "channels:               $(line.channels)  ")
end


"""
`ElectronCapture.amplitude(kind::String, channel::ElectronCapture.Channel, finalLevel::Level, continuumLevel::Level,
                            grid::Radial.Grid; printout::Bool=true)`  
    ... to compute, for the interelectronic interaction `kind` -- CoulombInteraction(), BreitInteraction(), CoulombBreit()
        or CoulombGaunt() -- the ElectronCapture amplitude 
        <alpha_f J_f || O^(capture, kind) || (alpha_i J_i, kappa) J_f> = <(alpha_i J_i, kappa) J_f || O^(Auger, kind) || alpha_f J_f>^*
        due to the interelectronic interaction for the given final and initial level. A value::ComplexF64 is returned.
"""
function amplitude(kind::AbstractEeInteraction, channel::ElectronCapture.Channel, finalLevel::Level, continuumLevel::Level,
                   grid::Radial.Grid; printout::Bool=true)
    amplitude = conj(AutoIonization.amplitude(kind, channel.kappa, channel.phase, continuumLevel::Level, finalLevel::Level, grid; printout=printout))
    
    return( amplitude )
end


"""
`ElectronCapture.computeAmplitudesProperties(line::ElectronCapture.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                settings::ElectronCapture.Settings; printout::Bool=true)` 
    ... to compute all amplitudes and properties of the given line; a line::ElectronCapture.Line is returned for which the amplitudes 
        and properties are now evaluated.
"""
function computeAmplitudesProperties(line::ElectronCapture.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                        settings::ElectronCapture.Settings; printout::Bool=true) 
    newChannels = ElectronCapture.Channel[];   contSettings = Continuum.Settings(false, nrContinuum);   rate = 0.
    ## BOTH levels must be expressed in ONE ordered subshell list before anything is contracted; taking the final
    ## level's own list, as this function did until 02-Sep-2026, leaves the initial ion's 1s appended at the end and
    ## the contraction raises "Improper subshell order".  This mirrors AutoIonization.computeAmplitudesProperties.
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    redILevel    = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    newfLevel    = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
    for channel in line.channels
        ## The free electron approaches the INITIAL ion, so its orbital is generated in that field, and the
        ## scattering state (ion + electron) carries the symmetry of the resonance it forms.
        ## The empty 101 subshell goes onto the RESONANCE and the captured electron onto the ION -- two DIFFERENT
        ## levels.  Putting both on one level appends 101s twice and the contraction refuses.
        newdLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newfLevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa),
                                                            redILevel, nm, grid, contSettings)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, redILevel)
        newChannel = ElectronCapture.Channel(channel.kappa, channel.symmetry, phase, 0.)
        amplitude  = ElectronCapture.amplitude(settings.operator, newChannel, newdLevel, newcLevel, grid, printout=printout)
        rate       = rate + conj(amplitude) * amplitude
        push!( newChannels, ElectronCapture.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude) )
    end
    ## Fritzsche, Kabachnik & Surzhykov, PRA 78, 032703 (2008), Eq. (3):
    ##
    ##     P_cap(0 -> d) = 2pi / (2(2J_0+1)) * SUM_lj |<a_d J_d || V || a_0 J_0, lj : J_d>|^2
    ##                   = (2J_d+1) / (2(2J_0+1)) * P_A(d -> 0)
    ##
    ## WHERE THE (2J_d+1) COMES FROM, since this is a Rule 18 site and was got wrong once.  AutoIonization forms its
    ## Auger rate as 2pi * SUM |A|^2 with NO 1/(2J_d+1), and that is CORRECT: for autoionization the scattering state
    ## (f + eps kappa) is coupled to exactly J_i, so the full and the reduced matrix element coincide and no such
    ## factor arises.  Comparing with Eq. (3) therefore fixes the relation between the two normalizations,
    ## <...>_paper = sqrt(2J_d+1) * A_JAC, and the capture rate must carry the (2J_d+1) explicitly.
    ##
    ## The ALIGNMENT below is untouched by this: Eq. (4) carries the same amplitudes in its numerator and in N, so
    ## the factor cancels exactly and A_k0 is the same either way.  Only the RATE changes.
    totalRate = 2pi * rate * (Basics.twice(line.finalLevel.J) + 1) / (2 * (Basics.twice(line.initialLevel.J) + 1))
    newLine   = ElectronCapture.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, Float64[], newChannels)
    if  settings.calcAlignment    newLine = ElectronCapture.Line(newLine.initialLevel, newLine.finalLevel,
                                                newLine.electronEnergy, newLine.totalRate,
                                                ElectronCapture.alignmentParameters(newLine), newLine.channels)   end
    
    return( newLine )
end


"""
`ElectronCapture.alignmentParameters(line::ElectronCapture.Line)`  
    ... computes the reduced alignment parameters A_k0(alpha_d J_d) of the doubly excited state that the resonant
        capture produces, from the capture amplitudes already stored on the line's channels; an
        alignment::Array{Float64,1} with the even ranks k = 2, 4, ... <= 2J_d is returned, and it is empty for
        J_d < 1, where no alignment exists.

        THIS IS THE OBSERVABLE OF THIS MODULE, and the reason the module has one at all. A capture rate cannot be
        measured on its own -- the resonance always decays -- but the ALIGNMENT of the state it produces can be, and
        it belongs to the capture step alone: in a pure two-step process every angular and polarization property of
        what follows is fixed by the alignment of the intermediate state times the decay amplitudes (Balashov,
        Grum-Grzhimailo & Kabachnik, *Polarization and Correlation Phenomena in Atomic Collisions*, Sec. 4.3.2).
        The DR resonance STRENGTH is deliberately NOT computed here: it contains the branching ratio
        A_r/(A_a + A_r), so it is not a capture-only quantity, and it is what DielectronicRecombination produces.

        Following Fritzsche, Kabachnik & Surzhykov, PRA 78, 032703 (2008), Eq. (4),

            A_k0 = N^-1 [J_d]^(1/2) SUM_(l l' j j')  (-1)^(J_d + J_0 - 1/2) [l,l',j,j']^(1/2) <l 0 l' 0 | k 0>
                        x { j   l   1/2 }   { j    J_d  J_0 }
                          { l'  j'  k   }   { J_d  j'   k   }
                        x <a_d J_d || V || a_0 J_0, l  j  : J_d>
                        x <a_d J_d || V || a_0 J_0, l' j' : J_d>^*

        with N = SUM_(lj) |<a_d J_d || V || a_0 J_0, l j : J_d>|^2 and [a,b,...] = (2a+1)(2b+1)... The parameters
        are nonzero only for EVEN k, because an unpolarized beam on an unpolarized ion can align the resonance but
        not orient it.
"""
function alignmentParameters(line::ElectronCapture.Line)
    Jd = Basics.twice(line.finalLevel.J) / 2.0;    J0 = Basics.twice(line.initialLevel.J) / 2.0
    if  Jd < 1.0    return( Float64[] )    end
    ## the partial wave (l, j) of each channel follows from its Dirac quantum number kappa
    lOf(kappa) = kappa < 0  ?  -kappa - 1  :  kappa
    jOf(kappa) = abs(kappa) - 0.5
    norm = 0.
    for  ch in line.channels    norm = norm + abs(ch.amplitude)^2    end
    if  norm == 0.    return( Float64[] )    end

    alignment = Float64[]
    for  k = 2:2:Int(floor(2*Jd))
        wa = 0.0 + 0.0im
        for  cha in line.channels,  chb in line.channels
            l  = lOf(cha.kappa);   j  = jOf(cha.kappa)
            lp = lOf(chb.kappa);   jp = jOf(chb.kappa)
            cg = AngularMomentum.ClebschGordan(AngularJ64(l), AngularM64(0), AngularJ64(lp), AngularM64(0),
                                               AngularJ64(k), AngularM64(0))
            if  cg == 0.    continue    end
            w6a = AngularMomentum.Wigner_6j(j, l, 1//2, lp, jp, k)
            w6b = AngularMomentum.Wigner_6j(j, Jd, J0, Jd, jp, k)
            if  w6a == 0.  ||  w6b == 0.    continue    end
            wa  = wa + sqrt( (2l+1)*(2lp+1)*(2j+1)*(2jp+1) ) * cg * w6a * w6b * cha.amplitude * conj(chb.amplitude)
        end
        wa = (-1)^Int(round(Jd + J0 - 0.5)) * sqrt(2*Jd + 1) * wa / norm
        push!(alignment, real(wa))
    end

    return( alignment )
end


"""
`ElectronCapture.determineChannels(finalLevel::Level, initialLevel::Level, settings::ElectronCapture.Settings)`  
    ... determines the partial waves of the captured electron that couple the initial ion to the doubly excited
        final level; an Array{ElectronCapture.Channel,1} with all amplitudes still zero is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::ElectronCapture.Settings)
    channels = ElectronCapture.Channel[]
    symi     = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf     = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    for  kappa in AngularMomentum.allowedKappaSymmetries(symi, symf)
        if  abs(kappa) > settings.maxKappa      continue    end
        push!( channels, ElectronCapture.Channel(kappa, symf, 0., Complex(0.)) )
    end

    return( channels )
end


"""
`ElectronCapture.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ElectronCapture.Settings)`  
    ... selects the capture lines by forming all pairs of an initial (ion) and a final (doubly excited) level that
        pass the line selection and whose capture energy lies inside [minCaptureEnergy, maxCaptureEnergy]; an
        Array{ElectronCapture.Line,1} with all amplitudes still zero is returned.

        The captured electron's energy is E(final) - E(initial) with the SIGN of an electron that is bound into a
        resonance: it must be POSITIVE, since the doubly excited level lies above the initial ion by exactly the
        kinetic energy the free electron brings in.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ElectronCapture.Settings)
    lines = ElectronCapture.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)      continue    end
            energy = fLevel.energy - iLevel.energy
            if  energy < settings.minCaptureEnergy  ||  energy > settings.maxCaptureEnergy    continue    end
            channels = ElectronCapture.determineChannels(fLevel, iLevel, settings)
            if  length(channels) == 0                                               continue    end
            push!( lines, ElectronCapture.Line(iLevel, fLevel, energy, 0., Float64[], channels) )
        end
    end

    return( lines )
end


"""
`ElectronCapture.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                              grid::Radial.Grid, settings::ElectronCapture.Settings; output=true)`  
    ... computes the capture rates and alignment parameters of all selected lines; an
        Array{ElectronCapture.Line,1} is returned if output=true, and nothing otherwise.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                      settings::ElectronCapture.Settings; output=true)
    println("")
    printstyled("ElectronCapture.computeLines(): The computation of capture rates and alignments starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = ElectronCapture.determineLines(finalMultiplet, initialMultiplet, settings)
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    if  settings.printBefore    ElectronCapture.displayLines(stdout, lines)    end
    maxEnergy   = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    newLines    = ElectronCapture.Line[]
    for  line in lines
        push!( newLines, ElectronCapture.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) )
    end
    ElectronCapture.displayRates(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   ElectronCapture.displayRates(iostream, newLines, settings)    end

    if  output    return( newLines )   else    return( nothing )    end
end


"""
`ElectronCapture.displayLines(stream::IO, lines::Array{ElectronCapture.Line,1})`  
    ... lists the selected capture lines and their partial waves before the amplitudes are computed; a neat table
        is printed but nothing is returned otherwise.
"""
function displayLines(stream::IO, lines::Array{ElectronCapture.Line,1})
    nx = 105
    println(stream, " ");   println(stream, "  Selected electron-capture lines:");    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);            sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);            sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy e_capt"; na=3)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.flushleft(37, "List of partial waves"; na=4)
    sb = sb * TableStrings.flushleft(37, "kappa                    "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        sa = sa * TableStrings.flushleft(37, join([string(ch.kappa) * "  " for ch in line.channels]); na=2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`ElectronCapture.displayRates(stream::IO, lines::Array{ElectronCapture.Line,1}, settings::ElectronCapture.Settings)`  
    ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned 
        otherwise.
"""
function  displayRates(stream::IO, lines::Array{ElectronCapture.Line,1}, settings::ElectronCapture.Settings)
    nx = 130
    println(stream, " ")
    println(stream, "  Electron-capture rates and alignment of the captured state: \n")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Electron energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(14, "Capture rate"; na=2);       
    sb = sb * TableStrings.center(14, TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.flushleft(32, "Alignment A_20, A_40, ..."; na=2)
    sb = sb * TableStrings.flushleft(32, "of the captured state"; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.totalRate))              * "    "
        ## No alignment is printed where none exists: a J_d < 1 resonance cannot be aligned, whatever the capture.
        if  length(line.alignment) == 0    sa = sa * "-- (J_d < 1)"
        else                               sa = sa * join([@sprintf("%9.5f", a) for a in line.alignment], "  ")
        end
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
