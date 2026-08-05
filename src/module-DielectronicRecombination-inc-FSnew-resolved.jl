
##
## RE-STRUCTURED, FINE-STRUCTURE RESOLVED DIELECTRONIC RECOMBINATION  (04/05-Aug-2026)
##
## Dielectronic recombination is physically a TWO-step process, and this file represents it as such:
##
##      capture          i + e^-  -->  m        described by a CaptureLine, one per (i,m)
##      stabilization    m        -->  f + hv   described by a PhotonLine,  one per (m,f)
##
## The older route in module-DielectronicRecombination-inc-FS-resolved.jl instead carries three overlapping
## structs -- Pathway (i,m,f), Passage (i,m) and Resonance (i,m) -- together with TWO parallel implementations
## of the same aggregation over m. Those two copies had drifted apart and were broken in complementary ways:
## the pathway route ASSIGNED the total Auger width where it had to SUM it and produced NaN for intermediate
## levels without any open channel, while the passage route had both of these right but crashed outright on a
## threadid()/nthreads() mismatch. Neither failure was visible because the two were never run side by side.
##
## The essential change here is therefore not the naming but WHERE the coupling between the two steps lives.
## The two total widths
##
##      Gamma_a(m) = sum_i A_a(m --> i)         Gamma_r(m) = sum_f A_r(m --> f)
##
## are properties of the intermediate level m alone. They are accumulated ONCE, in setTotalRates(), and then
## STORED on every CaptureLine. Nothing downstream re-derives them by scanning a list, which is what removes
## the entire class of error found above.
##
## The resonance strength then factorizes cleanly. With
##
##      C(i,m) = pi^2 / k^2 * A_a(m --> i) * (2J_m+1) / (2J_i+1)        [ k = wave number of the captured e^- ]
##
## one has, per gauge,      S(i,m) = C(i,m) * Gamma_r(m) / ( Gamma_a(m) + Gamma_r(m) )
##
## and the intensity of the individual DR satellite line (i,m,f) follows by combining the two objects:
##
##      I(i,m,f) = S(i,m) * A_r(m --> f) / Gamma_r(m)
##
## which is what settings.calcPhotonSpectrum makes available.
##
## Names carrying a "New" suffix (displayResultsNew, ...) only do so because the un-suffixed name is still
## taken by the old route; all three files share ONE module namespace. They lose the suffix when
## module-DielectronicRecombination-inc-FS-resolved.jl retires.
##


"""
`struct  DielectronicRecombination.CaptureLine`
    ... defines a type for a dielectronic capture line, i.e. the first step  i + e^-  -->  m  of dielectronic
        recombination, together with the two total widths of the intermediate level m and the resonance strength
        that follows from them.

    + initialLevel      ::Level             ... initial-(state) level i of the recombining ion.
    + intermediateLevel ::Level             ... intermediate-(state), autoionizing level m.
    + electronEnergy    ::Float64           ... energy of the captured electron, E_m - E_i + electronEnergyShift.
    + captureRate       ::Float64           ... Auger rate A_a(m --> i) of THIS channel alone.
    + totalAugerRate    ::Float64           ... Gamma_a(m) = sum over ALL initial levels i.
    + totalPhotonRate   ::EmProperty        ... Gamma_r(m) = sum over ALL final levels f.
    + resonanceStrength ::EmProperty        ... DR resonance strength of this (i,m) resonance.
    + captureChannels   ::Array{AutoIonization.Channel,1}  ... List of capture channels.
"""
struct  CaptureLine
    initialLevel        ::Level
    intermediateLevel   ::Level
    electronEnergy      ::Float64
    captureRate         ::Float64
    totalAugerRate      ::Float64
    totalPhotonRate     ::EmProperty
    resonanceStrength   ::EmProperty
    captureChannels     ::Array{AutoIonization.Channel,1}
end


"""
`DielectronicRecombination.CaptureLine()`  ... constructor for an 'empty' instance of a CaptureLine.
"""
function CaptureLine()
    em = EmProperty(0., 0.)
    CaptureLine(Level(), Level(), 0., 0., 0., em, em, AutoIonization.Channel[])
end


# `Base.show(io::IO, line::DielectronicRecombination.CaptureLine)`  ... prepares a proper printout of a CaptureLine.
function Base.show(io::IO, line::DielectronicRecombination.CaptureLine)
    println(io, "initialLevel:               $(line.initialLevel)  ")
    println(io, "intermediateLevel:          $(line.intermediateLevel)  ")
    println(io, "electronEnergy:             $(line.electronEnergy)  ")
    println(io, "captureRate:                $(line.captureRate)  ")
    println(io, "totalAugerRate:             $(line.totalAugerRate)  ")
    println(io, "totalPhotonRate:            $(line.totalPhotonRate)  ")
    println(io, "resonanceStrength:          $(line.resonanceStrength)  ")
    println(io, "captureChannels:            $(line.captureChannels)  ")
end


"""
`struct  DielectronicRecombination.PhotonLine`
    ... defines a type for the radiative stabilization  m  -->  f + hv,  the second step of dielectronic
        recombination. A PhotonLine depends on (m,f) only; it knows nothing about the initial level, since the
        radiative rate does not. The intensity with which a given photon line appears in a DR spectrum follows
        by combination with a CaptureLine, see the header of this file.

    + intermediateLevel ::Level             ... intermediate-(state), autoionizing level m.
    + finalLevel        ::Level             ... final-(state) level f of the recombined ion.
    + photonEnergy      ::Float64           ... energy of the emitted photon, E_m - E_f + photonEnergyShift.
    + photonRate        ::EmProperty        ... radiative rate A_r(m --> f) of this line.
    + photonChannels    ::Array{PhotoEmission.Channel,1}  ... List of photon channels.
"""
struct  PhotonLine
    intermediateLevel   ::Level
    finalLevel          ::Level
    photonEnergy        ::Float64
    photonRate          ::EmProperty
    photonChannels      ::Array{PhotoEmission.Channel,1}
end


"""
`DielectronicRecombination.PhotonLine()`  ... constructor for an 'empty' instance of a PhotonLine.
"""
function PhotonLine()
    PhotonLine(Level(), Level(), 0., EmProperty(0., 0.), PhotoEmission.Channel[])
end


# `Base.show(io::IO, line::DielectronicRecombination.PhotonLine)`  ... prepares a proper printout of a PhotonLine.
function Base.show(io::IO, line::DielectronicRecombination.PhotonLine)
    println(io, "intermediateLevel:          $(line.intermediateLevel)  ")
    println(io, "finalLevel:                 $(line.finalLevel)  ")
    println(io, "photonEnergy:               $(line.photonEnergy)  ")
    println(io, "photonRate:                 $(line.photonRate)  ")
    println(io, "photonChannels:             $(line.photonChannels)  ")
end


"""
`DielectronicRecombination.determineCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                                                 initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)`
    ... to determine a list of capture lines (i,m) for which the amplitudes are subsequently computed; apart from the level
        specification and the electron energy, all physical properties are set to zero here.
        An Array{DielectronicRecombination.CaptureLine,1} is returned.

        The finalMultiplet is required although no final level enters a CaptureLine: settings.pathwaySelection selects
        level TRIPLES (i,m,f), and the (i,m) pair must be kept whenever SOME final level makes the triple selected.
        Testing for existence over f reproduces the selection semantics of the old route exactly.
"""
function  determineCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                settings::DielectronicRecombination.Settings)
    captureLines        = DielectronicRecombination.CaptureLine[]
    electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            eEnergy = nLevel.energy - iLevel.energy + electronEnergyShift
            if  eEnergy < 0.    continue    end
            ## Keep the (i,m) pair if any final level makes the triple selected
            isSelected = false
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)   isSelected = true;   break   end
            end
            if  !isSelected     continue    end
            cChannels = DielectronicRecombination.determineCaptureChannels(nLevel, iLevel, settings)
            push!( captureLines, DielectronicRecombination.CaptureLine(iLevel, nLevel, eEnergy, 0., 0., EmProperty(0., 0.),
                                                                       EmProperty(0., 0.), cChannels) )
        end
    end
    return( captureLines )
end


"""
`DielectronicRecombination.determinePhotonLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                                                initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)`
    ... to determine a list of photon lines (m,f) for which the amplitudes are subsequently computed; apart from the level
        specification and the photon energy, all physical properties are set to zero here.
        An Array{DielectronicRecombination.PhotonLine,1} is returned.

        As in determineCaptureLines, the initialMultiplet enters only through settings.pathwaySelection, which selects
        level triples; the (m,f) pair is kept whenever SOME initial level makes the triple selected.
"""
function  determinePhotonLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                               settings::DielectronicRecombination.Settings)
    photonLines       = DielectronicRecombination.PhotonLine[]
    photonEnergyShift = Defaults.convertUnits("energy: to atomic", settings.photonEnergyShift)
    #
    for  nLevel  in  intermediateMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            pEnergy = nLevel.energy - fLevel.energy + photonEnergyShift
            if  pEnergy < 0.    continue    end
            isSelected = false
            for  iLevel  in  initialMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)   isSelected = true;   break   end
            end
            if  !isSelected     continue    end
            pChannels = DielectronicRecombination.determinePhotonChannels(fLevel, nLevel, settings)
            if  length(pChannels) == 0      continue    end
            push!( photonLines, DielectronicRecombination.PhotonLine(nLevel, fLevel, pEnergy, EmProperty(0., 0.), pChannels) )
        end
    end
    return( photonLines )
end


"""
`DielectronicRecombination.computeCaptureAmplitudes(captureLine::DielectronicRecombination.CaptureLine, nm::Nuclear.Model,
                                    grid::Radial.Grid, nrContinuum::Int64, settings::DielectronicRecombination.Settings)`
    ... to compute the Auger amplitudes and the capture rate A_a(m --> i) of the given capture line; a new
        captureLine::DielectronicRecombination.CaptureLine is returned in which captureRate and captureChannels are filled.
        The two total widths and the resonance strength remain zero here; they are set later by setTotalRates().
"""
function  computeCaptureAmplitudes(captureLine::DielectronicRecombination.CaptureLine, nm::Nuclear.Model, grid::Radial.Grid,
                                   nrContinuum::Int64, settings::DielectronicRecombination.Settings)
    rateA             = 0.
    newcChannels      = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum)
    initialLevel      = deepcopy(captureLine.initialLevel)
    intermediateLevel = deepcopy(captureLine.intermediateLevel)
    ## NOTE: this sets a GLOBAL default and is therefore not thread-safe; the same is true of the corresponding
    ## loops in the old route. It is harmless as long as JAC runs single-threaded, but it is the reason why the
    ## amplitude loops here must not be assumed safe under -t N without further work.
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    #
    for  cChannel in captureLine.captureChannels
        newnLevel = Basics.generateLevelWithSymmetryReducedBasis(intermediateLevel, intermediateLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, newnLevel.basis.subshells)
        newnLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, cChannel.kappa), newnLevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(captureLine.electronEnergy, Subshell(101, cChannel.kappa),
                                                            newiLevel, nm, grid, contSettings)
        newcLevel   = Basics.generateLevelWithExtraElectron(cOrbital, cChannel.symmetry, newiLevel)
        amplitude   = AutoIonization.amplitude(settings.augerOperator, cChannel, newcLevel, newnLevel, grid)
        rateA       = rateA + conj(amplitude) * amplitude
        push!( newcChannels, AutoIonization.Channel( cChannel.kappa, cChannel.symmetry, phase, amplitude) )
    end
    captureRate = 2pi * rateA
    #
    return( DielectronicRecombination.CaptureLine(captureLine.initialLevel, captureLine.intermediateLevel,
                                                  captureLine.electronEnergy, captureRate, 0., EmProperty(0., 0.),
                                                  EmProperty(0., 0.), newcChannels) )
end


"""
`DielectronicRecombination.computePhotonAmplitudes(photonLine::DielectronicRecombination.PhotonLine, grid::Radial.Grid)`
    ... to compute the photon amplitudes and the radiative rate A_r(m --> f) of the given photon line; a new
        photonLine::DielectronicRecombination.PhotonLine is returned in which photonRate and photonChannels are filled.
"""
function  computePhotonAmplitudes(photonLine::DielectronicRecombination.PhotonLine, grid::Radial.Grid)
    finalLevel        = deepcopy(photonLine.finalLevel)
    intermediateLevel = deepcopy(photonLine.intermediateLevel)
    Defaults.setDefaults("relativistic subshell list", intermediateLevel.basis.subshells; printout=false)
    #
    newpChannels = PhotoEmission.Channel[];    rateC = 0.;    rateB = 0.
    for  pChannel in photonLine.photonChannels
        amplitude = PhotoEmission.amplitude(Emission(), pChannel.multipole, pChannel.gauge, photonLine.photonEnergy,
                                            finalLevel, intermediateLevel, grid, display=false, printout=false)
        push!( newpChannels, PhotoEmission.Channel( pChannel.multipole, pChannel.gauge, amplitude) )
        if       pChannel.gauge == Basics.Coulomb     rateC = rateC + abs(amplitude)^2
        elseif   pChannel.gauge == Basics.Babushkin   rateB = rateB + abs(amplitude)^2
        elseif   pChannel.gauge == Basics.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
        end
    end
    ##
    ## THE RADIATIVE RATE IS THE EINSTEIN A COEFFICIENT, with exactly the prefactor that
    ## PhotoEmission.computeAmplitudesProperties uses for the same transition and the same sum over
    ## |PhotoEmission.amplitude|^2:
    ##       wa = 8pi * alpha * omega / (2J_upper + 1)      [ J_upper = J_m, the emitting level ]
    ##
    ## CORRECTED 05-Aug-2026. Until then both DR routes used
    ##       wa = 8pi * alpha * omega / (2J_m + 1) * (2J_f + 1)  ;  wa = wa / pi
    ## i.e. they carried a spurious (2J_f+1) together with a 1/pi that had been introduced while fitting one test
    ## case ("modified for test with Xe^53+ (March/2024)"). The two together gave
    ##       A_DR / A_Einstein = (2J_f + 1) / pi
    ## which is never unity. Measured for Li-like carbon over six independent transitions, the ratio was 0.63661
    ## for every J_f = 1/2 final level and 1.27322 for J_f = 3/2, against 2/pi = 0.63662 and 4/pi = 1.27324.
    ##
    ## Since Gamma_r(m) = sum_f A_r(m --> f) IS the total radiative width, the Einstein value is the required one.
    ## Why this survived: the effect is suppressed exactly where the one available benchmark sits. Radiatively
    ## dominated resonances (Gamma_r >> Gamma_a, K-LL in Xe^53+) have a branching ratio near unity and barely
    ## respond to a rescaling of Gamma_r, and their dominant J_f = 1 final levels give 3/pi = 0.955 anyway.
    ## Auger-dominated resonances (Gamma_a >> Gamma_r, K-LL in carbon) have S proportional to Gamma_r and took the
    ## full error. Absolute DR strengths computed before this date are therefore unreliable in the Auger-dominated
    ## regime, and mildly affected in the radiative one.
    ##
    wa = 8.0pi * Defaults.getDefaults("alpha") * photonLine.photonEnergy / (Basics.twice(photonLine.intermediateLevel.J) + 1)
    photonRate = EmProperty(wa * rateC, wa * rateB)
    #
    return( DielectronicRecombination.PhotonLine(photonLine.intermediateLevel, photonLine.finalLevel,
                                                 photonLine.photonEnergy, photonRate, newpChannels) )
end


"""
`DielectronicRecombination.setTotalRates(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                         photonLines::Array{DielectronicRecombination.PhotonLine,1})`
    ... to accumulate, for every intermediate level m, the two total widths

            Gamma_a(m) = sum_i A_a(m --> i)          Gamma_r(m) = sum_f A_r(m --> f)

        and to store them, together with the resulting resonance strength, on each capture line. A new
        Array{DielectronicRecombination.CaptureLine,1} is returned.

        This is deliberately a separate pass rather than something folded into the amplitude loops: it is the one
        place where the capture and the photon side couple, and it is exactly the step that both implementations of
        the old route got wrong -- once by assigning where it had to sum, once by producing NaN. Keeping it alone
        and short makes it readable and lets the two identities it establishes be tested directly.

        No `!` in the name: CaptureLine is immutable, so a new array is built rather than the old one modified.
"""
function  setTotalRates(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                        photonLines::Array{DielectronicRecombination.PhotonLine,1})
    ## Accumulate both widths, keyed by the index of the intermediate level
    totalAuger  = Dict{Int64,Float64}()
    totalPhoton = Dict{Int64,EmProperty}()
    for  cLine in captureLines
        m = cLine.intermediateLevel.index
        totalAuger[m]  = get(totalAuger,  m, 0.) + cLine.captureRate
    end
    for  pLine in photonLines
        m = pLine.intermediateLevel.index
        totalPhoton[m] = get(totalPhoton, m, EmProperty(0., 0.)) + pLine.photonRate
    end
    #
    newCaptureLines = DielectronicRecombination.CaptureLine[]
    for  cLine in captureLines
        m         = cLine.intermediateLevel.index
        gammaA    = get(totalAuger,  m, 0.)
        gammaR    = get(totalPhoton, m, EmProperty(0., 0.))
        ## C(i,m) = pi^2 / k^2 * A_a(m --> i) * (2J_m+1)/(2J_i+1)
        wavenb    = Defaults.convertUnits("kinetic energy to wave number: atomic units", cLine.electronEnergy)
        factor    = pi*pi / (wavenb*wavenb) * cLine.captureRate *
                    ((Basics.twice(cLine.intermediateLevel.J) + 1) / (Basics.twice(cLine.initialLevel.J) + 1))
        ## S = C(i,m) * Gamma_r / (Gamma_a + Gamma_r), per gauge. An intermediate level with NO open channel at all
        ## has zero total width; it contributes nothing and must be reported as ZERO. A NaN here would propagate
        ## through every later sum, in particular into the rate coefficient alpha(T).
        totC      = gammaA + gammaR.Coulomb
        totB      = gammaA + gammaR.Babushkin
        if  totC == 0.   sC = 0.   else   sC = factor * gammaR.Coulomb   / totC   end
        if  totB == 0.   sB = 0.   else   sB = factor * gammaR.Babushkin / totB   end
        push!( newCaptureLines, DielectronicRecombination.CaptureLine(cLine.initialLevel, cLine.intermediateLevel,
                                    cLine.electronEnergy, cLine.captureRate, gammaA, gammaR, EmProperty(sC, sB),
                                    cLine.captureChannels) )
    end
    return( newCaptureLines )
end


"""
`DielectronicRecombination.computeCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                            initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::DielectronicRecombination.Settings; output::Bool=true)`
    ... to compute dielectronic recombination in the re-structured, fine-structure resolved representation. This is the
        driver of this file and the analogue of DielectronicRecombination.computePathways of the old route.
        A tuple (captureLines, photonLines) is returned if output = true, and nothing otherwise; photonLines is empty
        unless settings.calcPhotonSpectrum.
"""
function  computeCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                              nm::Nuclear.Model, grid::Radial.Grid, settings::DielectronicRecombination.Settings; output::Bool=true)
    println("")
    printstyled("DielectronicRecombination.computeCaptureLines(): The computation of DR capture and photon lines starts now ... \n",
                color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    ## The same two consistency checks as the old route; both are reused unchanged
    DielectronicRecombination.checkConsistentMultiplets(finalMultiplet, intermediateMultiplet, initialMultiplet)
    DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet, intermediateMultiplet, initialMultiplet)
    #
    captureLines = DielectronicRecombination.determineCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    photonLines  = DielectronicRecombination.determinePhotonLines( finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
    if  settings.printBefore
        DielectronicRecombination.displayCaptureLines(stdout, captureLines)
        DielectronicRecombination.displayPhotonLines(stdout, photonLines)
    end
    ## Determine the maximum continuum energy and check the grid against it
    maxEnergy = 0.;   for  cLine in captureLines   maxEnergy = max(maxEnergy, cLine.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    ## The photon side is computed FIRST: Gamma_r(m) is needed for every resonance strength, whether or not the
    ## individual photon lines are afterwards retained. settings.calcPhotonSpectrum controls RETENTION, not work.
    ## Both loops write into a preallocated vector BY INDEX; indexing a per-thread buffer with threadid() is what
    ## crashed the old passage route, because nthreads() counts only the :default pool while threadid() also
    ## counts the :interactive one.
    newPhotonLines = Vector{DielectronicRecombination.PhotonLine}(undef, length(photonLines))
    @threads for  p in eachindex(photonLines)
        newPhotonLines[p] = DielectronicRecombination.computePhotonAmplitudes(photonLines[p], grid)
    end
    #
    newCaptureLines = Vector{DielectronicRecombination.CaptureLine}(undef, length(captureLines))
    @threads for  c in eachindex(captureLines)
        newCaptureLines[c] = DielectronicRecombination.computeCaptureAmplitudes(captureLines[c], nm, grid, nrContinuum, settings)
    end
    #
    newCaptureLines = DielectronicRecombination.setTotalRates(newCaptureLines, newPhotonLines)
    #
    ## Print all results to screen and, if requested, to the summary file
    DielectronicRecombination.displayResultsNew(stdout, newCaptureLines, newPhotonLines, settings)
    DielectronicRecombination.displayRateCoefficientsNew(stdout, newCaptureLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        DielectronicRecombination.displayResultsNew(iostream, newCaptureLines, newPhotonLines, settings)
        DielectronicRecombination.displayRateCoefficientsNew(iostream, newCaptureLines, settings)
    end
    #
    if  !settings.calcPhotonSpectrum    newPhotonLines = DielectronicRecombination.PhotonLine[]    end
    if  output    return( (newCaptureLines, newPhotonLines) )
    else          return( nothing )
    end
end


"""
`DielectronicRecombination.displayCaptureLines(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1})`
    ... to display a list of capture lines and channels that have been selected due to the prior settings. A neat table
        is printed but nothing is returned otherwise.
"""
function  displayCaptureLines(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1})
    nx = 120
    println(stream, " ")
    println(stream, "  Selected dielectronic capture lines  i + e^-  -->  m:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(16, "Levels"; na=4);            sb = sb * TableStrings.center(16, "i  --  m"; na=4)
    sa = sa * TableStrings.center(16, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(16, "i  --  m"; na=3)
    sa = sa * TableStrings.center(18, "Energies  " * TableStrings.inUnits("energy"); na=5)
    sb = sb * TableStrings.center(18, "electron     "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of kappas and total symmetries"; na=4)
    sb = sb * TableStrings.flushleft(57, "partial (total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  cLine in captureLines
        sa  = "  ";    isym = LevelSymmetry( cLine.initialLevel.J,      cLine.initialLevel.parity)
                       msym = LevelSymmetry( cLine.intermediateLevel.J, cLine.intermediateLevel.parity)
        sa = sa * TableStrings.center(16, TableStrings.levels_if(cLine.initialLevel.index, cLine.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(16, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", cLine.electronEnergy))  * "     "
        kappaSymmetryList = Tuple{Int64,LevelSymmetry}[]
        for  cChannel in cLine.captureChannels    push!( kappaSymmetryList, (cChannel.kappa, cChannel.symmetry) )    end
        sa = sa * TableStrings.kappaSymmetryTupels(85, kappaSymmetryList)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">> A total of $(length(captureLines)) dielectronic capture lines will be calculated.")
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayPhotonLines(stream::IO, photonLines::Array{DielectronicRecombination.PhotonLine,1})`
    ... to display a list of photon lines and channels that have been selected due to the prior settings. A neat table
        is printed but nothing is returned otherwise.
"""
function  displayPhotonLines(stream::IO, photonLines::Array{DielectronicRecombination.PhotonLine,1})
    nx = 110
    println(stream, " ")
    println(stream, "  Selected radiative stabilization lines  m  -->  f + hv:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(16, "Levels"; na=4);            sb = sb * TableStrings.center(16, "m  --  f"; na=4)
    sa = sa * TableStrings.center(16, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(16, "m  --  f"; na=3)
    sa = sa * TableStrings.center(18, "Energies  " * TableStrings.inUnits("energy"); na=5)
    sb = sb * TableStrings.center(18, "photon       "; na=5)
    sa = sa * TableStrings.flushleft(47, "List of multipoles and gauges"; na=4)
    sb = sb * TableStrings.flushleft(47, "partial (multipole, gauge)   "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  pLine in photonLines
        sa  = "  ";    msym = LevelSymmetry( pLine.intermediateLevel.J, pLine.intermediateLevel.parity)
                       fsym = LevelSymmetry( pLine.finalLevel.J,        pLine.finalLevel.parity)
        sa = sa * TableStrings.center(16, TableStrings.levels_if(pLine.intermediateLevel.index, pLine.finalLevel.index); na=4)
        sa = sa * TableStrings.center(16, TableStrings.symmetries_if(msym, fsym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pLine.photonEnergy))  * "     "
        mpGaugeList = Tuple{EmMultipole,EmGauge}[]
        for  pChannel in pLine.photonChannels    push!( mpGaugeList, (pChannel.multipole, pChannel.gauge) )    end
        sa = sa * TableStrings.multipoleGaugeTupels(65, mpGaugeList)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">> A total of $(length(photonLines)) radiative stabilization lines will be calculated.")
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayResultsNew(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            photonLines::Array{DielectronicRecombination.PhotonLine,1},
                            settings::DielectronicRecombination.Settings)`
    ... to list the total Auger and radiative rates and the resonance strengths of all capture lines and, if
        settings.calcPhotonSpectrum, the individual DR satellite lines that follow by combining the two line types.
        Neat tables are printed but nothing is returned otherwise.

        The first table is deliberately identical in layout and column content to displayResults(stream, resonances,
        settings) of the old route, so that the two implementations can be compared line by line.
"""
function  displayResultsNew(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            photonLines::Array{DielectronicRecombination.PhotonLine,1},
                            settings::DielectronicRecombination.Settings)
    nx = 160
    println(stream, " ")
    println(stream, "  Total Auger rates, radiative rates and resonance strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(42, "Auger rate     Cou -- rad. rates -- Bab"; na=1)
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=1)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=6)
    sa = sa * TableStrings.center(30, "Cou -- res. strength -- Bab"; na=3)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=2)
    sa = sa * TableStrings.center(18, "Widths Gamma_m"; na=2)
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=6)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  cLine in captureLines
        sa  = "";      isym = LevelSymmetry( cLine.initialLevel.J,      cLine.initialLevel.parity)
                       msym = LevelSymmetry( cLine.intermediateLevel.J, cLine.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(cLine.initialLevel.index, cLine.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", cLine.electronEnergy))                   * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", cLine.totalAugerRate))                     * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", cLine.totalPhotonRate.Coulomb))            * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", cLine.totalPhotonRate.Babushkin))          * "        "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", cLine.resonanceStrength.Coulomb))      * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", cLine.resonanceStrength.Babushkin))    * "     "
        wa = cLine.totalAugerRate + cLine.totalPhotonRate.Coulomb
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa))                                     * "   "
        wa = cLine.totalAugerRate + cLine.totalPhotonRate.Babushkin
        sa = sa * "(" * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa)) * ")"                         * "   "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    ## The DR satellite spectrum: intensity of the individual line (i,m,f) is  S(i,m) * A_r(m,f) / Gamma_r(m)
    if  settings.calcPhotonSpectrum   &&   length(photonLines) > 0
        nx = 152
        println(stream, " ")
        println(stream, "  DR satellite spectrum,  intensity = S(i,m) * A_r(m,f) / Gamma_r(m):")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(23, "Levels"; na=2);           sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2)
        sa = sa * TableStrings.center(23, "J^P symmetries"; na=0);   sb = sb * TableStrings.center(23, "i  --  m  --  f"; na=2)
        sa = sa * TableStrings.center(26, "Energies  " * TableStrings.inUnits("energy"); na=4)
        sb = sb * TableStrings.center(26, "electron        photon "; na=1)
        sa = sa * TableStrings.center(26, "A_r(m,f)  " * TableStrings.inUnits("rate"); na=2)
        sb = sb * TableStrings.center(26, " Cou -- rate -- Bab";  na=2)
        sa = sa * TableStrings.center(30, "Cou -- intensity -- Bab"; na=3)
        sb = sb * TableStrings.center(30, TableStrings.inUnits("strength"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  cLine in captureLines
            for  pLine in photonLines
                if  pLine.intermediateLevel.index != cLine.intermediateLevel.index    continue    end
                gammaR = cLine.totalPhotonRate
                if  gammaR.Coulomb   == 0.   iC = 0.   else   iC = cLine.resonanceStrength.Coulomb   * pLine.photonRate.Coulomb   / gammaR.Coulomb    end
                if  gammaR.Babushkin == 0.   iB = 0.   else   iB = cLine.resonanceStrength.Babushkin * pLine.photonRate.Babushkin / gammaR.Babushkin  end
                if  iC == 0.  &&  iB == 0.   continue    end
                sa  = "";      isym = LevelSymmetry( cLine.initialLevel.J,      cLine.initialLevel.parity)
                               msym = LevelSymmetry( cLine.intermediateLevel.J, cLine.intermediateLevel.parity)
                               fsym = LevelSymmetry( pLine.finalLevel.J,        pLine.finalLevel.parity)
                sa = sa * TableStrings.center(23, TableStrings.levels_imf(cLine.initialLevel.index, cLine.intermediateLevel.index,
                                                                          pLine.finalLevel.index); na=3)
                sa = sa * TableStrings.center(23, TableStrings.symmetries_imf(isym, msym, fsym);  na=3)
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", cLine.electronEnergy))       * "   "
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", pLine.photonEnergy))         * "     "
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pLine.photonRate.Coulomb))     * "  "
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", pLine.photonRate.Babushkin))   * "      "
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", iC))                       * "  "
                sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", iB))                       * "   "
                println(stream, sa)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end


"""
`DielectronicRecombination.computeRateCoefficientNew(captureLine::DielectronicRecombination.CaptureLine, temp::Float64)`
    ... to compute the DR plasma rate coefficient alpha^DR (T) of the given capture line for the temperature temp [K],
        in the isolated-resonance (delta-like) approximation. An alphaDR::EmProperty is returned in cm^3/s.
"""
function  computeRateCoefficientNew(captureLine::DielectronicRecombination.CaptureLine, temp::Float64)
    temp_au = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
    factor  = 4 / sqrt(2pi) * temp_au^(-3/2) * captureLine.electronEnergy * exp(- captureLine.electronEnergy/temp_au)
    alphaDR = factor * captureLine.resonanceStrength
    ## Convert into cm^3 / s
    factor  = Defaults.convertUnits("length: from atomic to fm", 1.0)^3 * 1.0e-39 *
                Defaults.convertUnits("rate: from atomic to 1/s", 1.0)
    alphaDR = factor * alphaDR
    #
    return( alphaDR )
end


"""
`DielectronicRecombination.displayRateCoefficientsNew(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                                      settings::DielectronicRecombination.Settings)`
    ... to list, if settings.calcRateAlpha, the DR plasma rate coefficients for the selected temperatures. Both the
        individual and the total rate coefficients are printed; nothing is returned otherwise.
"""
function  displayRateCoefficientsNew(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                     settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures)
    if  !settings.calcRateAlpha  ||  ntemps == 0     return( nothing )     end
    ## Unlike the old displayRateCoefficients, which silently showed only min(ntemps,7) columns while
    ## extractRateCoefficients returned all of them, every requested temperature is printed here; the table is
    ## simply split into blocks of at most seven columns.
    nblocks = cld(ntemps, 7)
    for  nb = 1:nblocks
        nt1 = 7*(nb-1) + 1;    nt2 = min(7*nb, ntemps)
        nx  = 54 + 17 * (nt2 - nt1 + 1)
        println(stream, " ")
        println(stream, "  DR rate coefficients for delta-like resonances [cm^3/s]:        ... all results in Babushkin gauge")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(14, "Energy"   ; na=2)
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=2)
        for  nt = nt1:nt2
            sa = sa * TableStrings.center(14, "T = " * @sprintf("%.2e", settings.temperatures[nt]); na=3)
            sb = sb * TableStrings.center(14, "[K]"; na=3)
        end
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  cLine in captureLines
            sa  = "";      isym = LevelSymmetry( cLine.initialLevel.J,      cLine.initialLevel.parity)
                           msym = LevelSymmetry( cLine.intermediateLevel.J, cLine.intermediateLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(cLine.initialLevel.index, cLine.intermediateLevel.index); na=4)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
            sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", cLine.electronEnergy))  * "      "
            for  nt = nt1:nt2
                alphaDR = DielectronicRecombination.computeRateCoefficientNew(cLine, settings.temperatures[nt])
                sa = sa * @sprintf("%.4e", alphaDR.Babushkin)  * "       "
            end
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        sa = "       alpha^DR (T, i; Coulomb gauge):                      "
        sb = "       alpha^DR (T, i; Babushkin gauge):                    "
        for  nt = nt1:nt2
            wa = EmProperty(0., 0.)
            for  cLine in captureLines   wa = wa + DielectronicRecombination.computeRateCoefficientNew(cLine, settings.temperatures[nt])   end
            sa = sa * @sprintf("%.4e", wa.Coulomb)    * "       "
            sb = sb * @sprintf("%.4e", wa.Babushkin)  * "       "
        end
        println(stream, sa);    println(stream, sb)
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end


"""
`DielectronicRecombination.extractRateCoefficientsNew(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                                      settings::DielectronicRecombination.Settings)`
    ... to extract the total DR plasma rate coefficients, summed over all capture lines, for each of the selected
        temperatures. An Array{EmProperty,1} of length length(settings.temperatures) is returned, and an empty array
        if settings.calcRateAlpha is false or no temperature is given.
"""
function  extractRateCoefficientsNew(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                     settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures);     alphaDRs = EmProperty[]
    if  !settings.calcRateAlpha  ||  ntemps == 0     return( EmProperty[] )     end
    for  nt = 1:ntemps
        wa = EmProperty(0., 0.)
        for  cLine in captureLines   wa = wa + DielectronicRecombination.computeRateCoefficientNew(cLine, settings.temperatures[nt])   end
        push!( alphaDRs, wa)
    end
    #
    return( alphaDRs )
end
