
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
## Names carrying a "New" suffix (displayResults, ...) only do so because the un-suffixed name is still
## taken by the old route; all three files share ONE module namespace. They lose the suffix when
## module-DielectronicRecombination-inc-FS-resolved.jl retires.
##














#####################################################################################################################
## THE PHYSICAL FORM, built BESIDE the flat one above.  This module is a COMPOSITE: it has no channel type of its
## own, and its two line types are built entirely out of the other modules' --
##     CaptureLine.capturePartialWaves ::Array{AutoIonization.PartialWave,1}
##     PhotonLine.photonAmplitudes     ::Array{MultipoleAmplitude,1}
## so its line types are built from AutoIonization.PartialWave and
## MultipoleAmplitude (fd7e6de).
##
## WHERE EmPropertyC BELONGS HERE, AND WHERE IT DOES NOT.  The photon AMPLITUDES are complex and gauge-paired, so
## they are EmPropertyC -- and they arrive already so, inside the embedded PhotoEmission type, rather than as a
## field declared here.  The four observables of this module are RATES and STRENGTHS: captureRate and
## totalAugerRate are real and gauge-free (the Auger operator has no gauge), totalPhotonRate and
## resonanceStrength are real and gauge-dependent.  Float64 and EmProperty are therefore already right for all
## four, and an EmPropertyC would be wrong -- they are not complex quantities.
##
## What this buys is not a smaller struct but the DELETION of a third copy of PhotoEmission's gauge machinery:
## determinePhotonChannels reproduces its gauge loop and hasMagnetic flag verbatim, and computePhotonAmplitudes
## reproduces its three-way rateC/rateB/Magnetic accumulation.  Below, the first becomes a CALL into PhotoEmission
## and the second becomes one abs2 on an EmPropertyC.
#####################################################################################################################


"""
`struct  DielectronicRecombination.CaptureLine`
    ... as DielectronicRecombination.CaptureLine, but carrying partial waves instead of flat Auger channels.

    + initialLevel        ::Level             ... initial-(state) level i of the recombining ion.
    + intermediateLevel   ::Level             ... intermediate-(state), autoionizing level m.
    + electronEnergy      ::Float64           ... energy of the captured electron.
    + captureRate         ::Float64           ... Auger rate A_a(m --> i) of THIS channel alone.
    + totalAugerRate      ::Float64           ... Gamma_a(m) = sum over ALL initial levels i.
    + totalPhotonRate     ::EmProperty        ... Gamma_r(m) = sum over ALL final levels f.
    + resonanceStrength   ::EmProperty        ... DR resonance strength of this (i,m) resonance.
    + capturePartialWaves ::Array{AutoIonization.PartialWave,1}  ... partial waves of the captured electron.
"""
struct  CaptureLine
    initialLevel        ::Level
    intermediateLevel   ::Level
    electronEnergy      ::Float64
    captureRate         ::Float64
    totalAugerRate      ::Float64
    totalPhotonRate     ::EmProperty
    resonanceStrength   ::EmProperty
    capturePartialWaves ::Array{AutoIonization.PartialWave,1}
end


"""
`struct  DielectronicRecombination.PhotonLine`
    ... as DielectronicRecombination.PhotonLine, but carrying one entry per MULTIPOLE, each holding both gauges,
        instead of one flat channel per (multipole, gauge).

    + intermediateLevel ::Level             ... intermediate-(state), autoionizing level m.
    + finalLevel        ::Level             ... final-(state) level f of the recombined ion.
    + photonEnergy      ::Float64           ... energy of the emitted photon.
    + photonRate        ::EmProperty        ... radiative rate A_r(m --> f) of this line.
    + photonAmplitudes  ::Array{MultipoleAmplitude,1}  ... one entry per multipole.
"""
struct  PhotonLine
    intermediateLevel   ::Level
    finalLevel          ::Level
    photonEnergy        ::Float64
    photonRate          ::EmProperty
    photonAmplitudes    ::Array{MultipoleAmplitude,1}
end










"""
`DielectronicRecombination.setTotalRates(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                         photonLines::Array{DielectronicRecombination.PhotonLine,1},
                                         empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to accumulate, for every intermediate level m, the two total widths

            Gamma_a(m) = sum_i A_a(m --> i)          Gamma_r(m) = sum_f A_r(m --> f)  [ + hydrogenic ]

        and to store them, together with the resulting resonance strength, on each capture line. A new
        Array{DielectronicRecombination.CaptureLine,1} is returned.

        This is deliberately a separate pass rather than something folded into the amplitude loops: it is the one
        place where the capture and the photon side couple, and it is exactly the step that both implementations of
        the old route got wrong -- once by assigning where it had to sum, once by producing NaN. Keeping it alone
        and short makes it readable and lets the two identities it establishes be tested directly.

        It is also the natural home of the HydrogenicCorrections: the decay of the captured Rydberg spectator is a
        contribution to Gamma_r(m) and to nothing else, so it belongs where Gamma_r(m) is formed. The old route
        instead evaluated it inside the per-(i,m) amplitude loop, which recomputed the same m-only quantity once for
        every initial level. Note the consequence for the internal identity: with hydrogenic corrections active,

            sum_f A_r(m --> f)  over the PhotonLines   +   hydrogenic(m)   =   totalPhotonRate(m),

        i.e. the plain sum over PhotonLines no longer reproduces the stored width. The added amount is printed per
        intermediate level so that it can always be subtracted back out.

        No `!` in the name: CaptureLine is immutable, so a new array is built rather than the old one modified.
"""
function  setTotalRates(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                        photonLines::Array{DielectronicRecombination.PhotonLine,1},
                        empTreatment::DielectronicRecombination.EmpiricalTreatment)
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
    ## The hydrogenic correction depends on m ALONE -- evaluate it once per intermediate level, not once per (i,m).
    ## Its re-autoionization cutoff needs the height of m above the LOWEST initial level, since the landing state
    ## can Auger back to any initial level below it; that is the LARGEST electronEnergy among the capture lines
    ## sharing this m.
    hydrogenic = Dict{Int64,Float64}()
    if  empTreatment.doHydrogenicCorrections
        maxEe = Dict{Int64,Float64}()
        for  cLine in captureLines
            m = cLine.intermediateLevel.index
            maxEe[m] = max(get(maxEe, m, -Inf), cLine.electronEnergy)
        end
        seen = Tuple{Int64,Int64,Int64,Float64,Float64,Int64,Float64}[]
        for  cLine in captureLines
            m = cLine.intermediateLevel.index
            if  haskey(hydrogenic, m)    continue    end
            rate, ni, li, nBound, dropped =
                DielectronicRecombination.computeHydrogenicPhotonRate(cLine.intermediateLevel, maxEe[m], empTreatment)
            explicitR      = get(totalPhoton, m, EmProperty(0., 0.)).Coulomb
            hydrogenic[m]  = rate
            totalPhoton[m] = get(totalPhoton, m, EmProperty(0., 0.)) + rate
            push!(seen, (m, ni, li, rate, explicitR, nBound, dropped))
        end
        DielectronicRecombination.displayHydrogenicCorrections(stdout, seen, empTreatment)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    DielectronicRecombination.displayHydrogenicCorrections(iostream, seen, empTreatment)    end
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
                                    cLine.capturePartialWaves) )
    end
    return( newCaptureLines )
end


"""
`DielectronicRecombination.computeCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                            initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::DielectronicRecombination.Settings; output::Bool=true)`
    ... to compute dielectronic recombination in the re-structured, fine-structure resolved representation. This is the
        driver of this file and the analogue of `computePathways` of the old route.
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
    ## All four corrections are implemented on this route: HydrogenicCorrections (with MaximumlCorrection capping the
    ## l of the Rydberg shell it acts on) inside setTotalRates where Gamma_r(m) is formed, RydbergTailCorrection and
    ## ResonanceWindowCorrection afterwards. An unknown correction type is refused rather than ignored -- silently
    ## dropping one returns a plausible number that is simply too small, which is the failure mode this whole
    ## re-structuring exists to eliminate.
    for  correction in settings.corrections
        if  !(typeof(correction) in [DielectronicRecombination.HydrogenicCorrections,
                                     DielectronicRecombination.MaximumlCorrection,
                                     DielectronicRecombination.RydbergTailCorrection,
                                     DielectronicRecombination.ResonanceWindowCorrection])
            error("\n\nDielectronicRecombination.computeCaptureLines():  STOP -- unknown correction type "     *
                  "$(typeof(correction)),\nrequested as\n\n    $correction\n"                                  *
                  ">>> Implemented are HydrogenicCorrections, MaximumlCorrection, RydbergTailCorrection and \n" *
                  "    ResonanceWindowCorrection.\n")
        end
    end
    ## Two consistency checks on the given multiplets; determineEmpiricalTreatment reads the shell structure of all
    ## three multiplets and is therefore run only once they are known to be mutually consistent.
    DielectronicRecombination.checkConsistentMultiplets(finalMultiplet, intermediateMultiplet, initialMultiplet)
    DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet, intermediateMultiplet, initialMultiplet)
    empTreatment = DielectronicRecombination.determineEmpiricalTreatment(finalMultiplet, intermediateMultiplet, nm,
                                                                         initialMultiplet, settings)
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
    ## Display-only, and the same for every line; set ONCE here rather than from inside the threaded loops below.
    Defaults.setStandardSubshellList(intermediateMultiplet.levels[1].basis.subshells; printout=false)
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
    ## Both are built ONCE for the whole computation and handed down to every capture line and every partial wave:
    ## the nuclear potential depends only on the nuclear model and the grid, the B-spline basis only on the grid --
    ## and the grid is fixed here, since Continuum.gridConsistency above is called once for the maximum energy.
    ## Both are READ-ONLY below (nothing writes into primitives.bsplinesL/bsplinesS), so sharing them across the
    ## threads of the following loop introduces no race; it also spares each thread its own copy of the basis.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    @threads for  c in eachindex(captureLines)
        newCaptureLines[c] = DielectronicRecombination.computeCaptureAmplitudes(captureLines[c], nm, grid, nrContinuum,
                                                                                settings; nuclearPot=nuclearPot,
                                                                                primitives=primitives)
    end
    #
    newCaptureLines = DielectronicRecombination.setTotalRates(newCaptureLines, newPhotonLines, empTreatment)
    #
    ## The Rydberg tail is built from the FINISHED explicit lines and is itself finished on construction, so it never
    ## re-enters setTotalRates. That separation is deliberate: the extrapolated lines carry the level identity of the
    ## shell they were scaled from, and an index-keyed aggregation would fold the whole series back onto that level.
    tailLines, tailShells, tailReport =
        DielectronicRecombination.computeRydbergTailLines(newCaptureLines, empTreatment)
    if  empTreatment.doRydbergTailCorrection
        DielectronicRecombination.displayRydbergTail(stdout, tailLines, tailShells, tailReport, empTreatment)
    end
    ## The energy window is applied to explicit and extrapolated lines alike
    newCaptureLines, windowReport = DielectronicRecombination.applyResonanceWindow(newCaptureLines, empTreatment)
    tailLines,       windowReport2 = DielectronicRecombination.applyResonanceWindow(tailLines, empTreatment)
    if  windowReport != ""    println(stdout, windowReport * " (explicit)");   println(stdout, windowReport2 * " (tail)")   end
    #
    ## Print all results to screen and, if requested, to the summary file
    DielectronicRecombination.displayResults(stdout, newCaptureLines, newPhotonLines, settings)
    DielectronicRecombination.displayRateCoefficients(stdout, newCaptureLines, settings, tailLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        if  empTreatment.doRydbergTailCorrection
            DielectronicRecombination.displayRydbergTail(iostream, tailLines, tailShells, tailReport, empTreatment)
        end
        DielectronicRecombination.displayResults(iostream, newCaptureLines, newPhotonLines, settings)
        DielectronicRecombination.displayRateCoefficients(iostream, newCaptureLines, settings, tailLines)
    end
    ## The tail lines join the returned list, so that every downstream consumer -- the rate coefficients, the Cascade
    ## simulations, the satellite diagnostic -- sees the full series rather than only the part that fitted in memory.
    append!(newCaptureLines, tailLines)
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
        ## One entry per partial wave; the total symmetry is that of the level which autoionizes -- the
        ## INTERMEDIATE level here -- and is the same for all of them.
        symn = LevelSymmetry(cLine.intermediateLevel.J, cLine.intermediateLevel.parity)
        for  pw in cLine.capturePartialWaves    push!( kappaSymmetryList, (pw.kappa, symn) )    end
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
        ## One row per (multipole, GAUGE): an electric multipole is one amplitude holding two gauges, a
        ## magnetic one has no gauge freedom.  The gauge is a property of this table, not of the amplitude.
        for  ma in pLine.photonAmplitudes
            if  string(ma.multipole)[1] == 'E'
                push!( mpGaugeList, (ma.multipole, Basics.Coulomb) );   push!( mpGaugeList, (ma.multipole, Basics.Babushkin) )
            else
                push!( mpGaugeList, (ma.multipole, Basics.Magnetic) )
            end
        end
        sa = sa * TableStrings.multipoleGaugeTupels(65, mpGaugeList)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">> A total of $(length(photonLines)) radiative stabilization lines will be calculated.")
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayHydrogenicCorrections(stream::IO,
                            seen::Array{Tuple{Int64,Int64,Int64,Float64,Float64,Int64,Float64},1},
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to list, for every intermediate level m, the hydrogenic rate that has been ADDED to the explicitly computed
        Gamma_r(m) in order to account for the radiative stabilization of the captured Rydberg spectator itself.
        A neat table is printed but nothing is returned otherwise.

        Two columns carry the message. `fraction` is how much of the corrected width the correction itself supplies:
        a few per cent means the explicit final-state set does most of the work and the estimate merely tidies up
        the tail, whereas anything approaching unity means the answer is being carried by a hydrogenic model with an
        adjustable Z_eff, and n^(final) should be raised instead of trusting it. `dropped` is the rate removed by
        the re-autoionization cutoff at n^(bound), i.e. spectator decays that emit a photon but leave the ion still
        above threshold; adding it back gives the upper bound of the bracket described on computeHydrogenicPhotonRate.
"""
function  displayHydrogenicCorrections(stream::IO,
                                       seen::Array{Tuple{Int64,Int64,Int64,Float64,Float64,Int64,Float64},1},
                                       empTreatment::DielectronicRecombination.EmpiricalTreatment)
    nx = 122
    println(stream, " ")
    println(stream, "  Hydrogenic corrections to the total photon rate  Gamma_r(m):")
    println(stream, " ")
    ## Print the EFFECTIVE range, i.e. after the internal cap at n_R - 1; nHydrogenic alone would mislead, since
    ## a value above the Rydberg shell (the Be-like gold app passes 22 against n_R = 19) contributes nothing.
    nRmax = 0;    for  entry in seen    nRmax = max(nRmax, entry[2])    end
    println(stream, "  + Rydberg decay  n_R, l_R  -->  n_f, l_R +- 1   summed over  " *
                    "$(empTreatment.nFinal+1) <= n_f <= $(min(empTreatment.nHydrogenic, nRmax-1))" *
                    (empTreatment.nHydrogenic > nRmax-1 ?
                        "   (nHydrogenic = $(empTreatment.nHydrogenic) capped at n_R - 1)" : ""))
    println(stream, "  + effective charge Z_eff = $(empTreatment.hydrogenicEffectiveZ),   " *
                    "rate scaling = $(empTreatment.hydrogenicRateScaling),   l_max = $(empTreatment.maximum_l)")
    println(stream, "  + these decays reach shells that are NOT represented in the final multiplet " *
                    "(explicit up to n = $(empTreatment.nFinal));")
    println(stream, "    without them Gamma_r(m) is too small and every resonance strength with it.")
    println(stream, "  + cut at n^(bound), below which the ion can Auger back: only decays that land BELOW the")
    println(stream, "    autoionization threshold complete the recombination.  'dropped' is what the cut removed;")
    println(stream, "    adding it back gives the upper bound, this table the lower one.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level m"; na=3);      sb = sb * TableStrings.center(10, "       "; na=3)
    sa = sa * TableStrings.center(12, "Rydberg"; na=3);      sb = sb * TableStrings.center(12, "n_R   l_R"; na=3)
    sa = sa * TableStrings.center( 9, "n^(bound)"; na=2);    sb = sb * TableStrings.center( 9, "         "; na=2)
    sa = sa * TableStrings.center(18, "Gamma_r explicit"; na=3)
    sb = sb * TableStrings.center(18, TableStrings.inUnits("rate"); na=3)
    sa = sa * TableStrings.center(18, "hydrogenic added"; na=3)
    sb = sb * TableStrings.center(18, TableStrings.inUnits("rate"); na=3)
    sa = sa * TableStrings.center(18, "dropped (unbound)"; na=3)
    sb = sb * TableStrings.center(18, TableStrings.inUnits("rate"); na=3)
    sa = sa * TableStrings.center(14, "fraction"; na=2);     sb = sb * TableStrings.center(14, "added/total"; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    maxFraction = 0.;   sumKept = 0.;   sumDropped = 0.
    for  (m, ni, li, rate, explicitR, nBound, dropped)  in  seen
        total    = explicitR + rate
        fraction = total == 0.  ?  0.  :  rate / total
        maxFraction = max(maxFraction, fraction)
        sumKept     = sumKept + rate;    sumDropped = sumDropped + dropped
        sa  = "  " * TableStrings.center(10, string(m); na=3)
        sa  = sa * TableStrings.center(12, "$ni    $li"; na=3)
        sa  = sa * TableStrings.center( 9, string(nBound); na=2)
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", explicitR)) * "        "
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", rate))      * "        "
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", dropped))   * "        "
        sa  = sa * @sprintf("%.3f", fraction)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    if  sumDropped > 0.
        println(stream, ">>> The re-autoionization cut removed " *
                        "$(round(100*sumDropped/(sumKept+sumDropped), digits=1)) % of the summed hydrogenic rate.")
        println(stream, "    Those decays do emit a photon, but leave the ion above threshold, where it can Auger " *
                        "back; a full")
        println(stream, "    radiative cascade would recover the part of them that stabilizes on a later step. " *
                        "Treat as a lower bound.")
    end
    if      maxFraction > 0.5
        println(stream, ">>> WARNING: the hydrogenic estimate supplies up to $(round(maxFraction*100, digits=1)) % of " *
                        "Gamma_r. The resonance strengths are then")
        println(stream, "    governed by a model rate with an adjustable Z_eff rather than by the computed structure; " *
                        "raise n^(final).")
    elseif  maxFraction > 0.
        println(stream, ">>> The hydrogenic estimate supplies at most $(round(maxFraction*100, digits=1)) % of Gamma_r.")
    end
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayResults(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            photonLines::Array{DielectronicRecombination.PhotonLine,1},
                            settings::DielectronicRecombination.Settings)`
    ... to list the total Auger and radiative rates and the resonance strengths of all capture lines and, if
        settings.calcPhotonSpectrum, the individual DR satellite lines that follow by combining the two line types.
        Neat tables are printed but nothing is returned otherwise.

        The first table is deliberately identical in layout and column content to displayResults(stream, resonances,
        settings) of the old route, so that the two implementations can be compared line by line.
"""
function  displayResults(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
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
`DielectronicRecombination.computeRateCoefficient(captureLine::DielectronicRecombination.CaptureLine, temp::Float64)`
    ... to compute the DR plasma rate coefficient alpha^DR (T) of the given capture line for the temperature temp [K],
        in the isolated-resonance (delta-like) approximation. An alphaDR::EmProperty is returned in cm^3/s.
"""
function  computeRateCoefficient(captureLine::DielectronicRecombination.CaptureLine, temp::Float64)
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
`DielectronicRecombination.displayRateCoefficients(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            settings::DielectronicRecombination.Settings,
                            tailLines::Array{DielectronicRecombination.CaptureLine,1}=DielectronicRecombination.CaptureLine[])`
    ... to list, if settings.calcRateAlpha, the DR plasma rate coefficients for the selected temperatures. Both the
        individual and the total rate coefficients are printed; nothing is returned otherwise.

        When a Rydberg tail was extrapolated, alpha(T) is reported SPLIT into the explicitly computed part and the
        extrapolated one. That split is the number which says whether any of the extrapolation mattered: a tail
        contributing a per cent or two can be quoted with the total, whereas one contributing half of alpha(T) means
        the result is being carried by a scaling law and the explicit calculation needs more shells. The split is also
        strongly temperature dependent -- the tail sits at higher electron energies than the explicit resonances only
        marginally, but it multiplies the level density -- so a single number would hide the effect.
"""
function  displayRateCoefficients(stream::IO, captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                  settings::DielectronicRecombination.Settings,
                                  tailLines::Array{DielectronicRecombination.CaptureLine,1}=DielectronicRecombination.CaptureLine[])
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
                alphaDR = DielectronicRecombination.computeRateCoefficient(cLine, settings.temperatures[nt])
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
            for  cLine in captureLines   wa = wa + DielectronicRecombination.computeRateCoefficient(cLine, settings.temperatures[nt])   end
            sa = sa * @sprintf("%.4e", wa.Coulomb)    * "       "
            sb = sb * @sprintf("%.4e", wa.Babushkin)  * "       "
        end
        println(stream, sa);    println(stream, sb)
        ## The explicit / extrapolated split, in Babushkin gauge, plus the share the tail carries
        if  length(tailLines) > 0
            sc = "       ... of which EXTRAPOLATED tail (Babushkin):        "
            sd = "       ... tail share of alpha^DR:                        "
            for  nt = nt1:nt2
                wexp = EmProperty(0., 0.);   wtail = EmProperty(0., 0.)
                for  cLine in captureLines   wexp  = wexp  + DielectronicRecombination.computeRateCoefficient(cLine, settings.temperatures[nt])   end
                for  cLine in tailLines      wtail = wtail + DielectronicRecombination.computeRateCoefficient(cLine, settings.temperatures[nt])   end
                share = wexp.Babushkin + wtail.Babushkin == 0.  ?  0.  :  wtail.Babushkin / (wexp.Babushkin + wtail.Babushkin)
                sc = sc * @sprintf("%.4e", wtail.Babushkin)  * "       "
                sd = sd * @sprintf("%13.3f", share)          * "    "
            end
            println(stream, sc);    println(stream, sd)
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream, "  NOTE: the totals above cover the EXPLICIT lines only; add the tail row for the full alpha^DR.")
            println(stream, "        A large tail share means the answer rests on the n^(-p) scaling rather than on the")
            println(stream, "        computed structure -- compute more Rydberg shells rather than raise nMax further.")
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end


"""
`DielectronicRecombination.extractRateCoefficients(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                                      settings::DielectronicRecombination.Settings)`
    ... to extract the total DR plasma rate coefficients, summed over all capture lines, for each of the selected
        temperatures. An Array{EmProperty,1} of length length(settings.temperatures) is returned, and an empty array
        if settings.calcRateAlpha is false or no temperature is given.
"""
function  extractRateCoefficients(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                     settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures);     alphaDRs = EmProperty[]
    if  !settings.calcRateAlpha  ||  ntemps == 0     return( EmProperty[] )     end
    for  nt = 1:ntemps
        wa = EmProperty(0., 0.)
        for  cLine in captureLines   wa = wa + DielectronicRecombination.computeRateCoefficient(cLine, settings.temperatures[nt])   end
        push!( alphaDRs, wa)
    end
    #
    return( alphaDRs )
end


##
## ---- Shared helpers, MOVED verbatim from module-DielectronicRecombination-inc-FS-resolved.jl (05-Aug-2026) ----
## Two consistency checks on the given multiplets, the two channel-determination routines, and the resonance
## selection predicate used by module-Cascade-inc-simulations.jl. They are route-independent, and are MOVED
## rather than copied so that this file no longer depends on the old one and that file can be retired simply
## by removing its include.
##


"""
`DielectronicRecombination.checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                                     initialMultiplet::Multiplet)`  
    ... to check that the given initial-, intermediate- and final-state levels and multiplets are consistent to each other and
        to avoid later problems with the computations. An error message is issued if an inconsistency occurs,
        and nothing is returned otherwise.
"""
function  checkConsistentMultiplets(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    initialSubshells      = initialMultiplet.levels[1].basis.subshells;             ni = length(initialSubshells)
    intermediateSubshells = intermediateMultiplet.levels[1].basis.subshells
    finalSubshells        = finalMultiplet.levels[1].basis.subshells;               nf = length(finalSubshells)
    
    if initialSubshells[1:end] == intermediateSubshells[1:ni]   &&
        intermediateSubshells[1:nf]  == finalSubshells
    else
        error("\nThe order of subshells must be equal for the initial-, intermediate and final states, and the same \n"  *
                "subshells must occur in the definition of the intermediate and final states. Only the initial states \n" *
                "may have FEWER subshells; this limitation arises from the angular coefficients.\n\n"                     *
                "    initial      = $initialSubshells \n"                                                                 *
                "    intermediate = $intermediateSubshells \n"                                                            *
                "    final        = $finalSubshells \n\n"                                                                 *
                ">>> The remedy is to write EVERY configuration over the SAME set of subshells, giving explicit ZERO \n"  *
                "    occupations wherever a shell is empty. For K-LL dielectronic recombination of an H-like ion: \n\n"   *
                "        initialConfigs      = [Configuration(\"1s^1 2s^0 2p^0\")] \n"                                    *
                "        intermediateConfigs = [Configuration(\"1s^0 2s^2 2p^0\"), Configuration(\"1s^0 2s^1 2p^1\")] \n" *
                "        finalConfigs        = [Configuration(\"1s^2 2s^0 2p^0\"), Configuration(\"1s^1 2s^1 2p^0\")] \n\n" *
                "    Writing Configuration(\"2s^2\") in place of Configuration(\"1s^0 2s^2 2p^0\") drops 1s from the \n"  *
                "    subshell list of the intermediate states and produces exactly this error.\n")
    end
        
    return( nothing )
end


"""
`DielectronicRecombination.checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                                      initialMultiplet::Multiplet)`  
    ... to check (and analyze) that all high nl orbitals in these multiplets are properly represented on the given grid.
        The function prints for each symmetry block kappa the high-nl orbitals and checks that they are all bound.
        An error message is issued if this is not the case, and nothing is returned otherwise.
"""
function  checkOrbitalRepresentation(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet)
    
    println("\n  Check the energies and orbital representation:")
    println("\n  ----------------------------------------------")
    subshells = initialMultiplet.levels[1].basis.subshells;     basis = initialMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Initial occupied subshells:      $(subshells)")
    
    subshells = intermediateMultiplet.levels[1].basis.subshells;     basis = intermediateMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.
            error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Intermediate occupied subshells: $(subshells)")
    
    subshells = finalMultiplet.levels[1].basis.subshells;     basis = finalMultiplet.levels[1].basis
    for sub in subshells   
        if  basis.orbitals[sub].energy >= 0.    error("$sub orbital not bound; enlarge the box size !")     end
    end
    println("  > Final occupied subshells:        $(subshells)")
        
    return( nothing )
end






"""
`DielectronicRecombination.isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)`  
    returns true, if level is to be excluded from the valid resonances, and false otherwise.
    It returns false if the ResonanceSelection() is inactive or if level belongs to the selected resoances.
    It is true only of ResonanceSelection() is active but the level does not belong to the selected resonances.
"""
function  isResonanceToBeExcluded(level::Level, refLevel, rSelection::ResonanceSelection)
    wa = false
    if !rSelection.active    return( wa )
        # This is the standard case if no additional limitations are specified by the user
    else
        # Analyze of whether level belongs to the selected resonances; it determines of whether level has 
        # electrons in either toShells or intoShells, and if there is one electron less in the fromShells
        fromwb    = false;    towb    = false;    intowb    = true
        confList  = Basics.extractConfigurations(Basics.FromBasis(), level.basis)
        occShells = Basics.extractShellList(confList) 
        for  shell in rSelection.toShells   
            if  shell in occShells   towb   = true;   break    end
        end
        for  shell in occShells   
            if  !(shell in rSelection.intoShells)   intowb = false;   break    end
        end
        # Compare the occupation of the given resonance level with those of the refLevel; there should be (at least)
        # one electron more in the shells of leadingConfig than in the same shells of level
        refConfig = Basics.extractConfiguration(Basics.LeadingConfiguration(), refLevel)
        levConfig = confList[1];   NoElectrons = 0
        for  (k,v) in refConfig.shells
            NoElectrons = NoElectrons + levConfig.shells[k]
        end
        if  NoElectrons < refConfig.NoElectrons     fromwb = true   end
        
        wa = !(fromwb &&  towb  && intowb)
    end
    
    return( wa )
end


"""
`DielectronicRecombination.computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64, Zeff::Float64)`
    ... to compute the non-relativistic electric-dipole rate for the transition from shell ni,li --> nf,lf of a
        hydrogenic ion with effective charge Zeff. The recursion formulas by Infeld and Hull (1951) are used
        together with the absorption oscillator strength. This makes the overall formulation/computation rather
        obscure, unfortunately. Uses SpecialFunctions.logfactorial. A rate::Float64 [a.u.] is returned.
        This procedure has been worked out by Stefan Schippers (2023).

        Moved here unchanged from module-DielectronicRecombination-inc-FS-resolved.jl when that route was retired
        (05-Aug-2026); the numerics are deliberately untouched.
"""
function  computeHydrogenicRate(ni::Int64, li::Int64, nf::Int64, lf::Int64,  Zeff::Float64)
    # Compute A(n,l) coeffient in the recursion formulas
    function computeA(n::Int64, l::Int64)
        A = 0.0
        if n>l  &&  n*l > 0     A = sqrt( (n+l) * (n-l)) / (n*l)     end
        return ( A )
    end
    # Compute I(n,l; n',l') integral in the recursion formulas
    function computeI(n::Int64, l::Int64, np::Int64, lp::Int64)
        wi = 0.0
        if       l>=n  ||  l<0  || lp>n  ||  lp>=np  ||  lp<0  || abs(l-lp)!=1    wi = 0.0
        elseif   l == lp-1
            if  lp == n
                wi = (n+2)*log(4*n*np)+(np-n-2)*log(np-n)-(np+n+2)*log(np+n) + 
                        0.5 * ( SpecialFunctions.logfactorial(np+n) - SpecialFunctions.logfactorial(np-n-1) -
                                SpecialFunctions.logfactorial(2*n-1) )
                wi = 0.25* exp(wi)
            else
                wi = (2*lp+1)   * computeA(np, lp+1) * computeI(n,lp, np,lp+1) + computeA(n, lp+1) * computeI(n,lp+1, np,lp)
                wi = wi / (2*lp * computeA(n, lp))
            end
        elseif   l == lp+1
            wi = computeA(np, l+1) * computeI(n,l, np,l+1) + (2*l+1) * computeA(n, l+1) * computeI(n,l+1, np,l)
            wi = wi / (2*l * computeA(np, l))
        else
            error("Unexpected set of quantum number n=$n l=$l  np=$np  lp=$lp ")
        end
        return ( wi )
    end
    # Compute absorption oscillator strength
    function computeOsc(n::Int64, l::Int64, np::Int64, lp::Int64)
        # This oscillator strength is used in absorption
        wx = 0.0
        if       l == lp+1   wx = (1/n^2 - 1/np^2) * (lp+1) / (2*(lp+1) +1) * computeI(n,lp+1, np,lp)^2
        elseif   l == lp-1   wx = (1/n^2 - 1/np^2) *  lp    / (2*(lp-1) +1) * computeI(n,lp-1, np,lp)^2
        else     error("stop a")
        end
        return( wx / 3.0 )
    end
    #
    rate = 0.0;
    if  abs(li-lf)!=1  ||  ni <= nf  ||  li >= ni  ||  li<0  ||  lf >= nf  ||  lf<0   return( rate )   
    elseif  nf > 40    error("Don't use a recursive scheme ... but make a new implementation for ni = $ni ")
    end
    #
    if      lf == (li + 1)     rate = (2*(li+1) + 1) / (2*li+1) * computeOsc(nf, li+1, ni, li)
    elseif  lf == (li - 1)     rate = (2*(li-1) + 1) / (2*li+1) * computeOsc(nf, li-1, ni, li)
    else    error("stop b")
    end
    rate = rate * Zeff^4 / 2.0 * Defaults.getDefaults("alpha")^3 * (1/nf^2 - 1/ni^2)^2

    return( rate )
end


"""
`DielectronicRecombination.determineEmpiricalTreatment(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                                    nm::Nuclear.Model, initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)` 
    ... to determine an instance of empiricalTreatment::EmpiricalTreatment that is (internally) applied to simplify
        the use of "corrections" to the DR strenghts. This data structure summarizes all parameters that help introduce
        several empirical corrections. The procedure is simple but slightly sophisticated as we wish to support "missing"
        parameters in the individual corrections as well as the knowledge that can be derived internally. 
        The definition of EmpiricalTreatment() can readily be extended as the need arises from the user side.

        Moved here from module-DielectronicRecombination-inc-FS-resolved.jl when that route was retired (05-Aug-2026),
        with three cosmetic repairs and no change to any derived quantity: a leftover @show of the shell lists was
        removed, a dead accumulator (occ, permanently zero) was dropped from the core-shell loop, and a missing
        maximum_l now resolves to 1000 = no limit rather than to 1. The old value of 1 contradicted both the default
        set a few lines above it and the meaning of the field, and would have silently suppressed every Rydberg shell
        with l > 1 for a user who wrote MaximumlCorrection(missing) to mean "no limit".
"""
function determineEmpiricalTreatment(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, 
                  nm::Nuclear.Model, initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)
    # First specify all parameters in turn
    doRydbergTailCorrection = doHydrogenicCorrections = doMaximumlCorrection = doResonanceWindowCorrection = false
    nHydrogenic            = nMax             = 0
    maximum_l              = 1000;    lMaxTail = -1;    nExponent = 3.0
    hydrogenicEffectiveZ   = hydrogenicRateScaling = 1.0
    resonanceEnergyMin     = 0.;    resonanceEnergyMax = 1.0e8;    tailEffectiveZ = 0.

    for correction in settings.corrections
        if  typeof(correction) == RydbergTailCorrection
            doRydbergTailCorrection = true;    nMax                   = correction.nMax
                                               tailEffectiveZ         = correction.effectiveZ
                                               lMaxTail               = correction.lMax
                                               nExponent              = correction.nExponent
        end
        if  typeof(correction) == HydrogenicCorrections
            doHydrogenicCorrections = true;    nHydrogenic            = correction.nHydrogenic
                                               hydrogenicEffectiveZ   = correction.effectiveZ
                                               hydrogenicRateScaling  = correction.rateScaling
        end
        if  typeof(correction) == MaximumlCorrection
            doMaximumlCorrection = true;       maximum_l = correction.maximum_l
        end
        if  typeof(correction) == ResonanceWindowCorrection
            doResonanceWindowCorrection = true;     resonanceEnergyMin = correction.energyMin
                                                    resonanceEnergyMax = correction.energyMax
        end
    end

    # Extract the n-shell quantum numbers from the various multiplets
    NoCoreElectrons          = initialMultiplet.levels[1].basis.NoElectrons
    basis                    = intermediateMultiplet.levels[1].basis
    intermediateConfs        = Basics.extractConfigurations(Basics.FromBasis(), basis)
    intermediateShellList    = Basics.extractNonrelativisticShellList(intermediateMultiplet.levels[1].basis.subshells)
    finalShellList           = Basics.extractNonrelativisticShellList(finalMultiplet.levels[1].basis.subshells)
    
    intermediateNs = Int64[];   for  shell  in  intermediateShellList    push!(intermediateNs, shell.n)   end
    finalNs        = Int64[];   for  shell  in  finalShellList           push!(finalNs,        shell.n)   end
    coreNs         = Int64[]
    for  conf  in  intermediateConfs     ne = 0
        for  shell in intermediateShellList
            if haskey(conf.shells, shell)   
                ne = ne + conf.shells[shell];    push!(coreNs, shell.n)
            end
            if  ne >= NoCoreElectrons    break   end
        end
    end
    #
    nCore          = maximum(coreNs)
    if  typeof(nHydrogenic)     == Missing    nHydrogenic = 0   end
    ## n^(lowest-captured) is the first n-shell that re-appears in the intermediate levels AFTER a gap; for a
    ## Be-like 1s^2 2s nl calculation with n = 1,2,3 and 19 explicit, this is 19.
    afterGap = false;   nLowestCaptured = 0
    for  i = 1:100
        if     i in intermediateNs && afterGap   nLowestCaptured = i;   break
        elseif !(i in intermediateNs)            afterGap = true
        end
    end
    ## n^(final) = the highest n treated explicitly in the final levels, EXCLUDING the captured Rydberg shell.
    ##
    ## CORRECTED 05-Aug-2026, see the note on computeHydrogenicPhotonRate. The old route set nFinal =
    ## maximum(finalNs) over the whole final basis. But a DR final level is |core> nl with the SAME captured
    ## electron, so the Rydberg shell is always present in that basis and nFinal came out equal to
    ## n^(lowest-captured) -- for the Be-like gold app, 19 rather than the intended 3. That contradicts the
    ## hierarchy this struct documents,
    ##
    ##      n^(core)  <  n^(final)  <  n^(hydrogenic)  <  n^(lowest-captured),
    ##
    ## and it made the correction range nFinal+1 : nHydrogenic start ABOVE the Rydberg shell, where the
    ## hydrogenic rate is identically zero.
    if  nLowestCaptured > 0    finalNs = filter(n -> n < nLowestCaptured, finalNs)    end
    if  isempty(finalNs)
        error("\n\nDielectronicRecombination.determineEmpiricalTreatment():  STOP -- after removing the captured \n" *
              "Rydberg shell (n >= $nLowestCaptured), no n-shell is left in the final levels, so n^(final) is \n"     *
              "undefined.\n>>> The final configurations must contain the low-n shells into which the resonance \n"    *
              "    stabilizes, not only the shell the electron was captured into.\n")
    end
    nFinal         = maximum(finalNs)

    if  typeof(nMax)                    == Missing    nMax                  = 0                          end
    if  typeof(tailEffectiveZ)          == Missing    tailEffectiveZ        = nm.Z - nCore               end
    if  typeof(lMaxTail)                == Missing    lMaxTail              = -1                         end
    if  typeof(nExponent)               == Missing    nExponent             = 3.0                        end
    if  typeof(hydrogenicEffectiveZ)    == Missing    hydrogenicEffectiveZ  = nm.Z - nCore               end
    if  typeof(hydrogenicRateScaling)   == Missing    hydrogenicRateScaling = 1.0                        end
    if  typeof(maximum_l)               == Missing    maximum_l             = 1000                       end
    if  typeof(resonanceEnergyMin)      == Missing    resonanceEnergyMin    = 0.                         end
    if  typeof(resonanceEnergyMax)      == Missing    resonanceEnergyMax    = 1.0e7                      end
    ## A RydbergTailCorrection that does not reach beyond the explicit shells is a user error, not a no-op:
    ## it reads as "extrapolate" while doing nothing at all.
    if  doRydbergTailCorrection  &&  nMax <= nLowestCaptured
        error("\n\nDielectronicRecombination.determineEmpiricalTreatment():  STOP -- RydbergTailCorrection was \n" *
              "requested with nMax = $nMax, but the highest explicitly computed Rydberg shell is already \n"       *
              "n^(lowest-captured) = $nLowestCaptured, so there is nothing above it to extrapolate to.\n"          *
              ">>> Raise nMax above $nLowestCaptured, or drop the correction.\n")
    end

    empTreatment = DielectronicRecombination.EmpiricalTreatment(doRydbergTailCorrection, doHydrogenicCorrections,
                                             doMaximumlCorrection, doResonanceWindowCorrection, nCore, nFinal, nHydrogenic,
                                             nLowestCaptured, nMax, maximum_l, lMaxTail, nExponent, hydrogenicEffectiveZ,
                                             hydrogenicRateScaling, tailEffectiveZ, resonanceEnergyMin, resonanceEnergyMax)
    return( empTreatment )
end


"""
`DielectronicRecombination.computeHydrogenicPhotonRate(intermediateLevel::Level, electronEnergy::Float64,
                            empTreatment::EmpiricalTreatment)`
    ... to estimate the radiative rate by which the captured Rydberg electron itself stabilizes the resonance, i.e.
        the sum of the hydrogenic rates A(n_R, l_R  -->  n_f, l_R +- 1) over all n^(final) < n_f <= n^(hydrogenic).
        These are decay channels into shells that are NOT represented explicitly in the final multiplet, so the
        explicitly computed Gamma_r(m) = sum_f A_r(m --> f) misses them entirely and is too small.

        This is a genuine physical decay channel, not a fudge: for high-n capture the spectator electron carries
        almost no binding, its decay is hydrogenic to good accuracy, and it stabilizes the ion below the
        autoionization threshold exactly as a core transition does. What IS adjustable is how well it is modeled --
        effectiveZ and rateScaling are sensitivity knobs, and a computation that leans hard on them is telling you
        that n^(final) was chosen too small, not that the physics is uncertain.

        RE-AUTOIONIZATION CUTOFF (added 05-Aug-2026). Not every such decay actually completes the recombination.
        The strength formula S = C * Gamma_r/(Gamma_a+Gamma_r) treats each radiative decay as TERMINATING the
        process, which is true only if the state one lands in lies below the autoionization threshold. The core
        stays excited during a spectator decay, so after n_R --> n_f the ion sits at

            E(landing) - E(threshold)  =  E_electron  +  B(n_R)  -  B(n_f),        B(n) = Zeff^2 / (2n^2),

        with E_electron the capture energy of this resonance. Where that is still positive the ion can Auger
        straight back to the initial level: the photon was emitted and the electron then left again, so counting
        the rate in Gamma_r overestimates the DR yield. The sum is therefore cut at the largest n_f satisfying
        B(n_f) > E_electron + B(n_R), i.e.

            n_f  <=  n^(bound)  =  floor( 1 / sqrt( 2*E_electron/Zeff^2 + 1/n_R^2 ) ).

        Two properties worth knowing. First, this is automatically consistent with the capture step: for the
        resonance group whose core excitation is what makes n_R the LOWEST capturable shell, E_electron is small
        and n^(bound) comes out as n_R - 1, so nothing is removed. Only groups built on a larger core excitation
        -- a 2p_3/2 rather than a 2p_1/2 hole, say -- are cut back, which is exactly right, since those sit far
        above threshold and the spectator must fall much further to get underneath it. Second, the cut makes the
        result a LOWER bound rather than the truth: a landing state above threshold may still radiate a second
        time -- the core decaying, or the spectator falling further -- and stabilize then. Resolving that needs a
        radiative cascade with branching at every step; the uncut sum is the upper bound and this one the lower,
        and both are reported so the bracket is visible.

        A tuple (rate::Float64 [a.u.], nRydberg::Int64, lRydberg::Int64, nBound::Int64, dropped::Float64 [a.u.])
        is returned, where `dropped` is the rate removed by the cutoff; all entries are zero if no hydrogenic
        corrections were requested.

        THE SUMMATION RANGE WAS WRONG IN THE OLD ROUTE and is corrected here (05-Aug-2026). Two independent
        defects combined to make the whole correction identically zero:

        (1) the live loop in computeAmplitudesProperties(::Passage) ran  nf = nFinal+1 : nHydrogenic, i.e. UPWARD
            from the explicit final shells, whereas computeHydrogenicRate returns exactly 0.0 unless nf < n_R;
        (2) nFinal was taken as the largest n in the final basis, which always contains the captured Rydberg shell
            itself, so nFinal = n^(lowest-captured) and the range started above n_R even in principle.

        The evidence that (1) is a slip and not a convention: the comment sitting directly above that loop reads
        "Compute and add hydrogenic rates A(ni,li --> nDetailed < n <= ni-1, li +- 1)", i.e. it describes the
        bound ni-1 that the code does not use; and the second, disabled copy of the same block in computeResonances
        -- the one carrying the note "This code need to be adapted" -- does loop  nf = nDetailed+1 : ni-1.
        So the correct range was written down twice and executed never.

        With both repaired, the range is  n^(final) < n_f <= min(n^(hydrogenic), n_R - 1)  and reads as the
        documented hierarchy intends: the shells below the Rydberg electron that are NOT represented explicitly.
"""
function  computeHydrogenicPhotonRate(intermediateLevel::Level, electronEnergy::Float64,
                                      empTreatment::EmpiricalTreatment)
    if  !empTreatment.doHydrogenicCorrections    return( (0., 0, 0, 0, 0.) )    end
    ## Identify the single Rydberg shell of this intermediate level
    rydbergSubshs = Basics.extractRydbergSubshellList(intermediateLevel, empTreatment.nLowestCaptured-1, 1.0e-1)
    rydbergShells = Basics.extractNonrelativisticShellList(rydbergSubshs)
    if  length(rydbergShells) != 1
        error("\n\nDielectronicRecombination.computeHydrogenicPhotonRate():  STOP -- the intermediate level "        *
              "$(intermediateLevel.index) \nhas $(length(rydbergShells)) Rydberg shells above n = "                  *
              "$(empTreatment.nLowestCaptured-1), namely $rydbergShells, \nbut the hydrogenic correction models "   *
              "the decay of exactly ONE spectator electron.\n>>> Either the intermediate configurations place "      *
              "more than one electron into high-n shells -- in which case this \n    correction does not apply -- "  *
              "or n^(lowest-captured) = $(empTreatment.nLowestCaptured) was mis-derived because the explicit \n"     *
              "    n-shells of the intermediate levels leave no gap below the captured shell.\n")
    end
    ni = rydbergShells[1].n;    li = rydbergShells[1].l
    ## MaximumlCorrection(lmax) suppresses the correction for the high-l Rydberg shells altogether
    if  li > empTreatment.maximum_l    return( (0., ni, li, 0, 0.) )    end
    #
    ## n^(bound): the deepest shell the spectator must reach for the landing state to be below the autoionization
    ## threshold. B(n) = Zeff^2/(2n^2) in atomic units; see the docstring for the derivation.
    Zeff   = empTreatment.hydrogenicEffectiveZ
    wa     = 2 * electronEnergy / (Zeff*Zeff)  +  1 / (ni*ni)
    nBound = wa <= 0.  ?  ni - 1  :  Int64(floor(1 / sqrt(wa)))
    #
    ## The Rydberg electron can only decay DOWNWARD, so the sum runs to n_R - 1; nHydrogenic caps it from above
    ## and n^(bound) cuts away the steps that leave the ion able to autoionize again.
    ## CORRECTED 05-Aug-2026 -- see the note above.
    rate = 0.;   dropped = 0.
    for  nf = empTreatment.nFinal+1 : min(empTreatment.nHydrogenic, ni-1)
        wb = DielectronicRecombination.computeHydrogenicRate(ni, li, nf, li-1, Zeff) +
             DielectronicRecombination.computeHydrogenicRate(ni, li, nf, li+1, Zeff)
        if  nf <= nBound    rate    = rate    + wb
        else                dropped = dropped + wb
        end
    end
    return( (rate * empTreatment.hydrogenicRateScaling, ni, li, nBound, dropped * empTreatment.hydrogenicRateScaling) )
end


"""
`DielectronicRecombination.measureRydbergExponent(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to MEASURE the exponent of the n-scaling of the capture rates from the user's own data, whenever two or more
        Rydberg shells were computed explicitly. Returns a tuple (exponent::Float64, nShells::Int64, details::String);
        exponent is NaN when fewer than two shells are available, in which case nothing can be measured.

        WHY THIS EXISTS. The whole Rydberg-tail extrapolation rests on a single assumption, A_a ~ n^(-p) with p = 3,
        which follows from the normalization of a Rydberg orbital near the core -- the capture amplitude samples the
        wavefunction where the core is, and |psi_nl(r->0)|^2 scales as n^(-3). That is textbook, but it is an
        assumption about THIS calculation only until it has been checked against it, and a computation that includes
        two Rydberg shells contains the check for free:

            p  =  ln( W(n1,l) / W(n2,l) )  /  ln( n2 / n1 ),        W(n,l) = sum_m A_a(m) (2J_m+1)  at shell n, momentum l

        The weighted sum W rather than a single A_a, because the individual levels of a shell need not correspond
        one-to-one between two shells, whereas the (2J+1)-weighted total over a given l does.
"""
function  measureRydbergExponent(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                 empTreatment::DielectronicRecombination.EmpiricalTreatment)
    ## Collect W(n,l) = sum_m A_a(m) * (2J_m+1) over the Rydberg shells actually present, and COUNT the levels that
    ## enter each. The count is the guard: the exponent is only meaningful if the two shells carry the same set of
    ## levels, and unequal counts mean they do not -- either because the shells interleave with something else, or
    ## because CI mixing left some levels without a single identifiable Rydberg shell.
    wnl = Dict{Tuple{Int64,Int64},Float64}();    cnl = Dict{Tuple{Int64,Int64},Int64}()
    for  cLine in captureLines
        nl = DielectronicRecombination.rydbergShellOf(cLine.intermediateLevel, empTreatment)
        if  nl === nothing    continue    end
        n, l = nl
        wnl[(n,l)] = get(wnl, (n,l), 0.) + cLine.captureRate * (Basics.twice(cLine.intermediateLevel.J) + 1)
        cnl[(n,l)] = get(cnl, (n,l), 0) + 1
    end
    nShells = sort(unique([k[1] for k in keys(wnl)]))
    if  length(nShells) < 2    return( (NaN, length(nShells), "") )    end
    #
    ## Every adjacent pair of shells, every l common to both, gives one estimate of p
    exponents = Float64[];    sa = ""
    for  i = 1:length(nShells)-1
        n1 = nShells[i];    n2 = nShells[i+1]
        for  l  in  sort(unique([k[2] for k in keys(wnl)]))
            w1 = get(wnl, (n1,l), 0.);    w2 = get(wnl, (n2,l), 0.)
            c1 = get(cnl, (n1,l), 0);     c2 = get(cnl, (n2,l), 0)
            if  w1 <= 0.  ||  w2 <= 0.    continue    end
            p  = log(w1/w2) / log(n2/n1)
            sa = sa * "\n    n = $n1 -> $n2,  l = $l :   W = " * @sprintf("%.4e", w1) * " -> " * @sprintf("%.4e", w2) *
                      "    p = " * @sprintf("%.3f", p) * "   ($c1 vs $c2 levels)"
            ## An unequal number of levels means the two shells are NOT the same series sampled twice, so their
            ## ratio measures the difference in level content and not the n-scaling. Report it and exclude it.
            if  c1 != c2
                sa = sa * "  <== EXCLUDED, unequal level content"
            else
                push!(exponents, p)
            end
        end
    end
    if  length(exponents) == 0    return( (NaN, length(nShells), "") )    end
    return( (sum(exponents)/length(exponents), length(nShells), sa) )
end


"""
`DielectronicRecombination.rydbergShellOf(level::Level, empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to return the (n,l) of the single Rydberg shell occupied in the given level, or `nothing` if the level has no
        Rydberg electron above n^(lowest-captured) - 1 or has more than one, in which case the Rydberg-series picture
        does not apply to it and it must simply be left alone rather than guessed at.
"""
function  rydbergShellOf(level::Level, empTreatment::DielectronicRecombination.EmpiricalTreatment)
    if  empTreatment.nLowestCaptured <= 1    return( nothing )    end
    rydbergSubshs = Basics.extractRydbergSubshellList(level, empTreatment.nLowestCaptured-1, 1.0e-1)
    rydbergShells = Basics.extractNonrelativisticShellList(rydbergSubshs)
    if  length(rydbergShells) != 1    return( nothing )    end
    return( (rydbergShells[1].n, rydbergShells[1].l) )
end


"""
`DielectronicRecombination.fitLProfile(wByL::Dict{Int64,Float64})`
    ... to fit  ln W(l) = a + b*l  by least squares over the computed orbital angular momenta, so that the l values
        which EXIST at a Rydberg shell but were never computed can be estimated. Returns a tuple
        (a::Float64, b::Float64, maxResidual::Float64, ok::Bool).

        `ok` is false -- and no l extrapolation must then be attempted -- when fewer than three l values are available
        or when the fit does not DECAY (b >= 0). Both cases mean the data do not constrain a profile, and inventing one
        would put an unbounded amount of strength into shells nobody has looked at. Refusing is the honest answer.

        The exponential ansatz is the crude part of this whole correction and is labelled as such wherever it is used:
        the capture rate falls with l because the centrifugal barrier keeps the electron away from the core, which is
        not exactly an exponential in l. What makes it tolerable is that it is fitted to the user's own numbers, that
        the residuals are reported, and that the SHARE of the total carried by the extrapolated l is reported too --
        if that share is large, the answer is to compute more l, not to trust the fit.
"""
function  fitLProfile(wByL::Dict{Int64,Float64})
    ls = sort([l for (l,w) in wByL  if w > 0.])
    if  length(ls) < 3    return( (0., 0., 0., false) )    end
    xs = Float64[l for l in ls];    ys = Float64[log(wByL[l]) for l in ls]
    nn = length(xs);    sx = sum(xs);   sy = sum(ys);   sxx = sum(xs.*xs);   sxy = sum(xs.*ys)
    den = nn*sxx - sx*sx
    if  den == 0.    return( (0., 0., 0., false) )    end
    b   = (nn*sxy - sx*sy) / den
    a   = (sy - b*sx) / nn
    maxRes = 0.;   for  i = 1:nn    maxRes = max(maxRes, abs(ys[i] - (a + b*xs[i])))   end
    return( (a, b, maxRes, b < 0.) )
end


"""
`DielectronicRecombination.computeRydbergTailLines(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to estimate the DR capture lines of the Rydberg shells  n^(lowest-captured) < n <= nMax  that are too numerous
        to be computed explicitly, by extrapolating the series that WAS computed. A tuple
        (tailLines::Array{CaptureLine,1}, tailShells::Array{Int64,1}, report::String) is returned, where tailShells[k]
        is the principal quantum number of tailLines[k] -- a CaptureLine does not carry it, and recovering it
        afterwards from the energy would be guesswork; the lines are fully formed -- their two total
        widths and their resonance strength are set here -- so that they never pass through setTotalRates and cannot
        contaminate the aggregation of the explicit levels, which is exactly how the previous implementation of this
        idea went wrong.

        For each shell n and each computed l, with n0 = n^(lowest-captured) and B(n) = Zeff^2/(2n^2):

            E_e(n)      = E_e(n0) + B(n0) - B(n)                         exact Rydberg shift, no ansatz
            A_a(n)      = A_a(n0) * (n0/n)^p                             p measured where possible, else 3
            Gamma_a(n)  = Gamma_a(n0) * (n0/n)^p
            Gamma_r(n)  = Gamma_r,core(n0) + hydrogenic spectator(n,l)    core part n-INDEPENDENT, spectator computed
            S(n)        = pi^2/k^2 * A_a * (2J_m+1)/(2J_i+1) * Gamma_r/(Gamma_a+Gamma_r)

        Holding the core radiative rate fixed while scaling the Auger rate is the physical content: the stabilizing
        core transition does not know where the spectator sits, whereas the capture amplitude samples the Rydberg
        wavefunction at the core and therefore carries the n^(-3). This is also what makes the tail matter -- the
        branching ratio Gamma_r/(Gamma_a+Gamma_r) creeps towards 1 as n grows, so the strength per shell falls much
        more slowly than the capture rate until the two widths cross.
"""
function  computeRydbergTailLines(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                                  empTreatment::DielectronicRecombination.EmpiricalTreatment)
    tailLines = DielectronicRecombination.CaptureLine[];    tailShells = Int64[];    report = ""
    if  !empTreatment.doRydbergTailCorrection    return( (tailLines, tailShells, report) )    end
    n0   = empTreatment.nLowestCaptured
    Zeff = empTreatment.tailEffectiveZ
    bind(n) = Zeff*Zeff / (2.0*n*n)
    #
    ## (1) The n-exponent: measure it if two shells are present, otherwise use the requested/default value
    pMeasured, nShells, details = DielectronicRecombination.measureRydbergExponent(captureLines, empTreatment)
    pUsed = empTreatment.nExponent
    if  isnan(pMeasured)
        report = report * "\n  + n-scaling exponent p = " * @sprintf("%.2f", pUsed) * " ASSUMED (only $nShells Rydberg " *
                          "shell computed; two are needed to measure it)."
    else
        report = report * "\n  + n-scaling exponent:  MEASURED p = " * @sprintf("%.3f", pMeasured) *
                          "   vs   assumed p = " * @sprintf("%.2f", pUsed) * details
        if  abs(pMeasured - pUsed) > 0.3
            report = report * "\n    >>> WARNING: measured and assumed exponent differ by " *
                              @sprintf("%.2f", abs(pMeasured-pUsed)) * ". The extrapolation uses the ASSUMED value; " *
                              "\n        set nExponent to the measured one if the deviation is believed."
        end
    end
    #
    ## (2) Group the explicit lines of the reference shell n0 by their Rydberg l
    byL = Dict{Int64,Array{DielectronicRecombination.CaptureLine,1}}()
    for  cLine in captureLines
        nl = DielectronicRecombination.rydbergShellOf(cLine.intermediateLevel, empTreatment)
        if  nl === nothing  ||  nl[1] != n0    continue    end
        byL[nl[2]] = push!(get(byL, nl[2], DielectronicRecombination.CaptureLine[]), cLine)
    end
    if  length(byL) == 0
        report = report * "\n  + NO capture line of the reference shell n0 = $n0 could be identified; no tail built."
        return( (tailLines, tailShells, report) )
    end
    lComputed = sort(collect(keys(byL)))
    report = report * "\n  + reference shell n0 = $n0 with l = $lComputed, extrapolated to $(n0+1) <= n <= $(empTreatment.nMax)"
    #
    ## (3) Extrapolate each computed l in n. Statistical weights, level identities and gauges are inherited from the
    ##     reference line, which is what makes these lines representable as ordinary CaptureLines at all.
    ## The CORE radiative rate of a reference line is its explicit Gamma_r minus whatever hydrogenic spectator rate
    ## was added to it there. It depends on the reference line ALONE, so it is formed once per line here rather than
    ## once per (line, shell) inside the loop below -- the subtraction is cheap but the extractRydbergSubshellList
    ## call hidden inside computeHydrogenicPhotonRate is not.
    refs = Tuple{Int64,DielectronicRecombination.CaptureLine,EmProperty}[]
    for  l  in  lComputed
        for  cLine  in  byL[l]
            hydRef, _, _, _, _ = DielectronicRecombination.computeHydrogenicPhotonRate(cLine.intermediateLevel,
                                                                                cLine.electronEnergy, empTreatment)
            push!( refs, (l, cLine, cLine.totalPhotonRate - EmProperty(hydRef, hydRef)) )
        end
    end
    for  n = n0+1 : empTreatment.nMax
        scale = (n0/n)^pUsed
        for  (l, cLine, coreR)  in  refs
            eEnergy = cLine.electronEnergy + bind(n0) - bind(n)
            if  eEnergy <= 0.    continue    end       ## no longer a resonance: the shell has dropped below threshold
            ## the spectator's own radiative decay IS recomputed at the new shell -- it is n-dependent
            gammaR  = coreR + DielectronicRecombination.hydrogenicRateAt(n, l, eEnergy, empTreatment)
            push!( tailLines, DielectronicRecombination.buildTailLine(cLine, eEnergy, cLine.captureRate * scale,
                                                                      cLine.totalAugerRate * scale, gammaR) )
            push!( tailShells, n )
        end
    end
    #
    ## (4) The l extrapolation: estimate the l that exist at these shells but were never computed
    tailL, shellL, reportL = DielectronicRecombination.computeLTailLines(byL, lComputed, empTreatment, pUsed)
    append!(tailLines, tailL);    append!(tailShells, shellL);    report = report * reportL
    #
    return( (tailLines, tailShells, report) )
end


"""
`DielectronicRecombination.buildTailLine(refLine::DielectronicRecombination.CaptureLine, electronEnergy::Float64,
                            captureRate::Float64, gammaA::Float64, gammaR::EmProperty)`
    ... to assemble one extrapolated CaptureLine with its resonance strength already formed. The strength uses exactly
        the expression of setTotalRates -- written out here rather than shared, because these lines deliberately do NOT
        pass through that aggregation: they carry the level identity of their reference line, so an index-keyed sum
        would fold the whole extrapolated series back onto the reference level. That is precisely the defect that made
        the previous implementation of this idea inflate Gamma_a for every member of the series.
"""
function  buildTailLine(refLine::DielectronicRecombination.CaptureLine, electronEnergy::Float64,
                        captureRate::Float64, gammaA::Float64, gammaR::EmProperty)
    wavenb = Defaults.convertUnits("kinetic energy to wave number: atomic units", electronEnergy)
    factor = pi*pi / (wavenb*wavenb) * captureRate *
             ((Basics.twice(refLine.intermediateLevel.J) + 1) / (Basics.twice(refLine.initialLevel.J) + 1))
    totC   = gammaA + gammaR.Coulomb;      totB = gammaA + gammaR.Babushkin
    sC     = totC == 0.  ?  0.  :  factor * gammaR.Coulomb   / totC
    sB     = totB == 0.  ?  0.  :  factor * gammaR.Babushkin / totB
    return( DielectronicRecombination.CaptureLine(refLine.initialLevel, refLine.intermediateLevel, electronEnergy,
                                                  captureRate, gammaA, gammaR, EmProperty(sC, sB),
                                                  AutoIonization.PartialWave[]) )
end


"""
`DielectronicRecombination.hydrogenicRateAt(n::Int64, l::Int64, electronEnergy::Float64,
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to evaluate the hydrogenic spectator rate for a Rydberg electron at shell (n,l) directly from the quantum
        numbers, i.e. without needing a Level to read them off. Same summation range and same re-autoionization cutoff
        as computeHydrogenicPhotonRate; an EmProperty is returned, gauge-independent because the estimate is.
"""
function  hydrogenicRateAt(n::Int64, l::Int64, electronEnergy::Float64,
                           empTreatment::DielectronicRecombination.EmpiricalTreatment)
    if  !empTreatment.doHydrogenicCorrections  ||  l > empTreatment.maximum_l    return( EmProperty(0., 0.) )    end
    Zeff   = empTreatment.hydrogenicEffectiveZ
    wa     = 2 * electronEnergy / (Zeff*Zeff)  +  1 / (n*n)
    nBound = wa <= 0.  ?  n - 1  :  Int64(floor(1 / sqrt(wa)))
    rate   = 0.
    for  nf = empTreatment.nFinal+1 : min(empTreatment.nHydrogenic, n-1, nBound)
        rate = rate + DielectronicRecombination.computeHydrogenicRate(n, l, nf, l-1, Zeff) +
                      DielectronicRecombination.computeHydrogenicRate(n, l, nf, l+1, Zeff)
    end
    rate = rate * empTreatment.hydrogenicRateScaling
    return( EmProperty(rate, rate) )
end


"""
`DielectronicRecombination.computeLTailLines(byL::Dict{Int64,Array{DielectronicRecombination.CaptureLine,1}},
                            lComputed::Array{Int64,1}, empTreatment::DielectronicRecombination.EmpiricalTreatment,
                            pUsed::Float64)`
    ... to estimate the contribution of the orbital angular momenta that EXIST at the Rydberg shells but were never
        computed at the reference shell. A tuple (lines::Array{CaptureLine,1}, shells::Array{Int64,1},
        report::String) is returned.

        At shell n the orbital angular momentum runs to n-1, whereas the explicit calculation carries only the few l
        the user could afford. Those missing l are not negligible in general: the capture rate falls with l, but the
        number of states rises as (2l+1), and the two only cancel once the centrifugal barrier bites.

        This is the most assumption-laden part of the correction and it is fenced accordingly. The weighted capture
        strength W(l) = sum_m A_a(m) (2J_m+1) is fitted by ln W = a + b l over the computed l, and the fit is used only
        if it DECAYS and rests on at least three points; otherwise no l extrapolation happens and the report says so.
        The extrapolation stops at lMax, at n-1, or once a term has fallen below 1e-4 of the l = 0 term.
"""
function  computeLTailLines(byL::Dict{Int64,Array{DielectronicRecombination.CaptureLine,1}},
                            lComputed::Array{Int64,1}, empTreatment::DielectronicRecombination.EmpiricalTreatment,
                            pUsed::Float64)
    lines = DielectronicRecombination.CaptureLine[];    shells = Int64[]
    lMax  = empTreatment.lMaxTail
    if  lMax <= maximum(lComputed)
        return( (lines, shells,
                 "\n  + l extrapolation NOT requested (lMax = $lMax <= largest computed l = $(maximum(lComputed)))") )
    end
    ## Weighted capture strength per computed l, at the reference shell
    wByL = Dict{Int64,Float64}()
    for  (l, cls)  in  byL
        wByL[l] = sum( cl.captureRate * (Basics.twice(cl.intermediateLevel.J) + 1)  for cl in cls )
    end
    a, b, maxRes, ok = DielectronicRecombination.fitLProfile(wByL)
    if  !ok
        return( (lines, shells, "\n  + l extrapolation REFUSED: the fit of ln W(l) over l = $lComputed is either based on " *
                        "fewer than\n    three points or does not decay. Extrapolating it would put an unbounded " *
                        "amount of strength\n    into shells that were never examined; compute more l instead.") )
    end
    report = "\n  + l extrapolation:  ln W(l) = " * @sprintf("%.3f", a) * " + " * @sprintf("%.3f", b) *
             "*l   (max residual " * @sprintf("%.3f", maxRes) * " in ln W)"
    #
    ## Represent each extrapolated l by the reference line of the LARGEST computed l, rescaled so that its weighted
    ## capture strength matches the fit. Its level identity is wrong in detail -- there is no computed level at that l
    ## -- but the quantities that enter the strength (A_a, the two widths, the weights) are all carried explicitly.
    lRef    = maximum(lComputed);    refLine = byL[lRef][1]
    wRef    = wByL[lRef]
    n0      = empTreatment.nLowestCaptured;    Zeff = empTreatment.tailEffectiveZ
    bind(n) = Zeff*Zeff / (2.0*n*n)
    hydRef, _, _, _, _ = DielectronicRecombination.computeHydrogenicPhotonRate(refLine.intermediateLevel,
                                                                               refLine.electronEnergy, empTreatment)
    coreR   = refLine.totalPhotonRate - EmProperty(hydRef, hydRef)
    w0      = exp(a)
    nAdded  = 0;    wExtrap = 0.
    for  n = n0 : empTreatment.nMax
        scale = (n0/n)^pUsed
        for  l = lRef+1 : min(lMax, n-1)
            wl = exp(a + b*l)
            if  wl < 1.0e-4 * w0    break    end
            if  n == n0    wExtrap = wExtrap + wl    end
            eEnergy = refLine.electronEnergy + bind(n0) - bind(n)
            if  eEnergy <= 0.    continue    end
            ## Scale the reference line so its weighted capture strength equals the fitted W(l), then apply n^(-p)
            aRate  = refLine.captureRate    * (wl/wRef) * scale
            gammaA = refLine.totalAugerRate * (wl/wRef) * scale
            gammaR = coreR + DielectronicRecombination.hydrogenicRateAt(n, l, eEnergy, empTreatment)
            push!( lines, DielectronicRecombination.buildTailLine(refLine, eEnergy, aRate, gammaA, gammaR) )
            push!( shells, n );    nAdded = nAdded + 1
        end
    end
    wComputed = sum(values(wByL))
    share     = wComputed + wExtrap == 0.  ?  0.  :  wExtrap / (wComputed + wExtrap)
    report = report * "\n    $nAdded extrapolated (n,l) lines;  at n0 the extrapolated l carry " *
                      @sprintf("%.1f", 100*share) * " % of the weighted capture strength"
    if  share > 0.3
        report = report * "\n    >>> WARNING: the extrapolated l carry more than 30 % of the total. The fit is then " *
                          "\n        doing more work than the calculation; compute more l at n0 rather than trust it."
    end
    return( (lines, shells, report) )
end


"""
`DielectronicRecombination.applyResonanceWindow(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to keep only those capture lines whose electron energy falls inside [energyMin, energyMax], for instance the
        range actually scanned in an EBIT or storage-ring measurement. A tuple (lines, report::String) is returned.
        The number of lines dropped is always reported: a silently shortened list looks exactly like a calculation
        that found fewer resonances.
"""
function  applyResonanceWindow(captureLines::Array{DielectronicRecombination.CaptureLine,1},
                               empTreatment::DielectronicRecombination.EmpiricalTreatment)
    if  !empTreatment.doResonanceWindowCorrection    return( (captureLines, "") )    end
    kept = filter(cl -> empTreatment.resonanceEnergyMin <= cl.electronEnergy <= empTreatment.resonanceEnergyMax,
                  captureLines)
    nDrop = length(captureLines) - length(kept)
    sa = "\n>>> ResonanceWindowCorrection: kept $(length(kept)) of $(length(captureLines)) capture lines with " *
         @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", empTreatment.resonanceEnergyMin)) * " <= E_e <= " *
         @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", empTreatment.resonanceEnergyMax)) * " " *
         TableStrings.inUnits("energy") * ";  $nDrop dropped."
    return( (kept, sa) )
end




"""
`DielectronicRecombination.displayRydbergTail(stream::IO, tailLines::Array{DielectronicRecombination.CaptureLine,1},
                            tailShells::Array{Int64,1}, report::String,
                            empTreatment::DielectronicRecombination.EmpiricalTreatment)`
    ... to summarize the extrapolated Rydberg tail: how it was built, what each shell contributes, where the two total
        widths cross, and whether the sum has converged by nMax. Nothing is returned.

        The shell-by-shell table is the one to read, and the last two columns say whether to believe the total.
        If the LAST shell still carries a noticeable share, nMax was set too low and the sum is still climbing.
        The Gamma_a/Gamma_r column says which regime each shell is in: while it is >> 1 the resonance is Auger
        dominated and its strength is nearly independent of n, so truncating there discards a lot; once it falls
        below 1 the strength has begun its n^(-3) descent and the remaining tail is genuinely small.
"""
function  displayRydbergTail(stream::IO, tailLines::Array{DielectronicRecombination.CaptureLine,1},
                             tailShells::Array{Int64,1}, report::String,
                             empTreatment::DielectronicRecombination.EmpiricalTreatment)
    println(stream, " ")
    println(stream, "  Rydberg-tail extrapolation of the DR resonances:")
    println(stream, report)
    if  length(tailLines) == 0    println(stream, " ");   return( nothing )    end
    #
    ## Aggregate per shell
    strengthOf = Dict{Int64,Float64}();   gaOf = Dict{Int64,Float64}();   grOf = Dict{Int64,Float64}()
    countOf    = Dict{Int64,Int64}();     energyOf = Dict{Int64,Float64}()
    for  k in eachindex(tailLines)
        n = tailShells[k];   cl = tailLines[k]
        strengthOf[n] = get(strengthOf, n, 0.) + cl.resonanceStrength.Coulomb
        gaOf[n]       = max(get(gaOf, n, 0.), cl.totalAugerRate)
        grOf[n]       = max(get(grOf, n, 0.), cl.totalPhotonRate.Coulomb)
        countOf[n]    = get(countOf, n, 0) + 1
        energyOf[n]   = cl.electronEnergy
    end
    shells = sort(collect(keys(strengthOf)))
    total  = sum(values(strengthOf))
    #
    nx = 108
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center( 8, "shell n"; na=3);    sb = sb * TableStrings.center( 8, "       "; na=3)
    sa = sa * TableStrings.center( 8, "lines"; na=3);      sb = sb * TableStrings.center( 8, "     "; na=3)
    sa = sa * TableStrings.center(16, "E_e  " * TableStrings.inUnits("energy"); na=3)
    sb = sb * TableStrings.center(16, "                "; na=3)
    sa = sa * TableStrings.center(18, "strength (Coulomb)"; na=3)
    sb = sb * TableStrings.center(18, TableStrings.inUnits("strength"); na=3)
    sa = sa * TableStrings.center(14, "share of tail"; na=3);   sb = sb * TableStrings.center(14, "             "; na=3)
    sa = sa * TableStrings.center(14, "Gamma_a/Gamma_r"; na=2); sb = sb * TableStrings.center(14, "               "; na=2)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    ## The crossover is located by scanning DOWNWARD for the last shell that is still Auger dominated, not by
    ## taking the first shell below 1. The lowest row of this table can hold only the l-extrapolated lines of the
    ## reference shell -- a handful of high-l lines with tiny capture rates and hence a tiny ratio -- and a
    ## forward scan would then report the crossover at the very first shell, which is exactly backwards.
    crossover = 0
    for  k = length(shells):-1:1
        n = shells[k];    ratio = grOf[n] == 0.  ?  Inf  :  gaOf[n] / grOf[n]
        if  ratio >= 1.0    crossover = (k == length(shells)) ? 0 : shells[k+1];    break    end
    end
    for  n in shells
        ratio = grOf[n] == 0.  ?  Inf  :  gaOf[n] / grOf[n]
        sa  = "  " * TableStrings.center( 8, string(n); na=3)
        sa  = sa * TableStrings.center( 8, string(countOf[n]); na=3)
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energyOf[n]))   * "       "
        sa  = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", strengthOf[n])) * "         "
        sa  = sa * @sprintf("%.4f", total == 0. ? 0. : strengthOf[n]/total) * "          "
        sa  = sa * (isinf(ratio) ? "     inf" : @sprintf("%8.2f", ratio))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    lastShare = total == 0. ? 0. : strengthOf[shells[end]] / total
    println(stream, ">>> Tail total (Coulomb) = " * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", total)) *
                    " " * TableStrings.inUnits("strength") * " from $(length(tailLines)) extrapolated lines.")
    if  crossover == 0
        println(stream, ">>> Gamma_a > Gamma_r at EVERY shell up to nMax = $(empTreatment.nMax): the series is still")
        println(stream, "    Auger dominated there, so the strength per shell has not yet started to fall as n^(-3)")
        println(stream, "    and the shells above nMax still matter. Raise nMax.")
    else
        println(stream, ">>> Gamma_a = Gamma_r near n = $crossover; above it the strength falls as n^(-3) and the sum converges.")
    end
    if  lastShare > 0.05
        println(stream, ">>> WARNING: the LAST shell (n = $(shells[end])) still carries " *
                        @sprintf("%.1f", 100*lastShare) * " % of the tail -- the sum has not converged; raise nMax.")
    else
        println(stream, ">>> The last shell carries " * @sprintf("%.2f", 100*lastShare) * " % of the tail; the sum has converged.")
    end
    #
    return( nothing )
end


#####################################################################################################################
## The physical form: the core, the bridge back to the flat one, and the drivers.  Everything below is ADDITIVE.
#####################################################################################################################


"""
`DielectronicRecombination.determineCaptureChannels(intermediateLevel::Level, initialLevel::Level,
                                                          settings::DielectronicRecombination.Settings)`
    ... as determineCaptureChannels, but returning partial waves; an
        Array{AutoIonization.PartialWave,1} is returned with all amplitudes still zero.

        Note that this is NOT AutoIonization.determineChannels: that one applies settings.maxKappa from an
        AutoIonization.Settings, which a DielectronicRecombination.Settings does not carry, and the flat
        determineCaptureChannels applies no kappa restriction at all. The selection rule itself --
        allowedKappaSymmetries(symi, symn) -- is the same, and is reproduced here unrestricted, exactly as the
        flat version has it.
"""
function determineCaptureChannels(intermediateLevel::Level, initialLevel::Level,
                                        settings::DielectronicRecombination.Settings)
    partialWaves = AutoIonization.PartialWave[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
    for  kappa in AngularMomentum.allowedKappaSymmetries(symi, symn)
        push!( partialWaves, AutoIonization.PartialWave(kappa, 0., 0., Complex(0.)) )
    end
    return( partialWaves )
end


"""
`DielectronicRecombination.determinePhotonChannels(finalLevel::Level, intermediateLevel::Level,
                                                         settings::DielectronicRecombination.Settings)`
    ... as determinePhotonChannels, but returning one entry per MULTIPOLE; an
        Array{MultipoleAmplitude,1} is returned with all amplitudes still zero.

        THIS IS A CALL, NOT A COPY.  The flat determinePhotonChannels reproduces PhotoEmission's gauge loop and
        its hasMagnetic flag verbatim -- a third copy of the same machinery, after PhotoEmission's own and
        PhotoRecombination's.  Here the selection rule is asked of the module that owns it. The only reason a
        wrapper exists at all is that PhotoEmission.determineChannelsClaude expects a PhotoEmission.Settings.
        NOTE the Claude suffix: it survives ONLY because PhotoEmission is not retired yet, and goes when it is,
        so the multipoles have to be handed over.
"""
function determinePhotonChannels(finalLevel::Level, intermediateLevel::Level,
                                       settings::DielectronicRecombination.Settings)
    peSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=settings.multipoles,
                                        gauges=settings.gauges)
    return( PhotoEmission.determineChannelsClaude(finalLevel, intermediateLevel, peSettings) )
end


"""
`DielectronicRecombination.computeCaptureAmplitudes(captureLine::DielectronicRecombination.CaptureLine,
        nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, settings::DielectronicRecombination.Settings;
        nuclearPot::Union{Nothing,Radial.Potential}=nothing, primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... as computeCaptureAmplitudes, but on partial waves; a CaptureLine with the Auger amplitudes and the
        capture rate filled is returned. The two total widths and the resonance strength remain zero here.
"""
function  computeCaptureAmplitudes(captureLine::DielectronicRecombination.CaptureLine, nm::Nuclear.Model,
                                         grid::Radial.Grid, nrContinuum::Int64, settings::DielectronicRecombination.Settings;
                                         nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                         primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    rateA             = 0.
    newPartialWaves   = AutoIonization.PartialWave[];   contSettings = Continuum.Settings(false, nrContinuum)
    initialLevel      = deepcopy(captureLine.initialLevel)
    intermediateLevel = deepcopy(captureLine.intermediateLevel)
    redNLevel = Basics.generateLevelWithSymmetryReducedBasis(intermediateLevel, intermediateLevel.basis.subshells)
    newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, redNLevel.basis.subshells)
    ## The one total symmetry of the scattering state: that of the level which autoionizes, i.e. the
    ## INTERMEDIATE level here.  See the note at AutoIonization.PartialWave.
    symn      = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)

    for  pw in captureLine.capturePartialWaves
        newnLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, pw.kappa), redNLevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(captureLine.electronEnergy, Subshell(101, pw.kappa),
                                                            newiLevel, nm, grid, contSettings;
                                                            nuclearPot=nuclearPot, primitives=primitives)
        newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, symn, newiLevel)
        ## The computed phase, restored on 14-Aug-2026 and as AutoIonization's own driver has always done:
        ## AutoIonization.amplitude multiplies by exp(-im*phase), which the capture path used to drop.
        amplitude = AutoIonization.amplitude(settings.augerOperator, pw.kappa, phase, newcLevel, newnLevel, grid)
        rateA     = rateA + conj(amplitude) * amplitude
        push!( newPartialWaves, AutoIonization.PartialWave(pw.kappa, captureLine.electronEnergy, phase, amplitude) )
    end
    captureRate = 2pi * rateA
    #
    return( DielectronicRecombination.CaptureLine(captureLine.initialLevel, captureLine.intermediateLevel,
                                                        captureLine.electronEnergy, captureRate, 0.,
                                                        EmProperty(0., 0.), EmProperty(0., 0.), newPartialWaves) )
end


"""
`DielectronicRecombination.computePhotonAmplitudes(photonLine::DielectronicRecombination.PhotonLine,
                                                         grid::Radial.Grid)`
    ... as computePhotonAmplitudes, but on the physical form; a PhotonLine with the amplitudes and the
        radiative rate filled is returned.

        THIS IS WHERE EmPropertyC EARNS ITS PLACE IN THIS MODULE.  The flat version carries a third copy of
        PhotoEmission's three-way accumulation --
              if  gauge == Coulomb  rateC += ...;  elseif Babushkin  rateB += ...;  elseif Magnetic  BOTH
        -- which here becomes one line, because abs2 of an EmPropertyC is an EmProperty and a magnetic amplitude
        has equal components and so enters both by itself.

        The Einstein-A prefactor is reproduced verbatim, including the correction of 05-Aug-2026: it is
        8pi * alpha * omega / (2J_m + 1), with NO (2J_f+1) and NO 1/pi.
"""
function  computePhotonAmplitudes(photonLine::DielectronicRecombination.PhotonLine, grid::Radial.Grid)
    finalLevel        = deepcopy(photonLine.finalLevel)
    intermediateLevel = deepcopy(photonLine.intermediateLevel)
    newAmplitudes     = MultipoleAmplitude[];    rate = EmProperty(0., 0.)
    for  ma in photonLine.photonAmplitudes
        mp = ma.multipole
        if  string(mp)[1] == 'E'
            ampC = PhotoEmission.amplitude(Emission(), mp, Basics.Coulomb,   photonLine.photonEnergy,
                                            finalLevel, intermediateLevel, grid, display=false, printout=false)
            ampB = PhotoEmission.amplitude(Emission(), mp, Basics.Babushkin, photonLine.photonEnergy,
                                            finalLevel, intermediateLevel, grid, display=false, printout=false)
            amp  = EmPropertyC(ampC, ampB)
        else
            ampM = PhotoEmission.amplitude(Emission(), mp, Basics.Magnetic,  photonLine.photonEnergy,
                                            finalLevel, intermediateLevel, grid, display=false, printout=false)
            amp  = EmPropertyC(ampM)
        end
        rate = rate + abs2(amp)
        push!( newAmplitudes, MultipoleAmplitude(mp, amp) )
    end
    wa = 8.0pi * Defaults.getDefaults("alpha") * photonLine.photonEnergy /
                 (Basics.twice(photonLine.intermediateLevel.J) + 1)
    photonRate = wa * rate
    #
    return( DielectronicRecombination.PhotonLine(photonLine.intermediateLevel, photonLine.finalLevel,
                                                       photonLine.photonEnergy, photonRate, newAmplitudes) )
end




"""
`DielectronicRecombination.determineCaptureLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                                    initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)`
`DielectronicRecombination.determinePhotonLines(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet,
                                    initialMultiplet::Multiplet, settings::DielectronicRecombination.Settings)`
    ... to determine the capture and photon lines for which amplitudes are subsequently computed; the selection,
        the energy shifts and the "keep the pair if SOME third level makes the triple selected" rule are
        reproduced exactly.
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
            isSelected = false
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelTriple(iLevel, nLevel, fLevel, settings.pathwaySelection)   isSelected = true;   break   end
            end
            if  !isSelected     continue    end
            pws = DielectronicRecombination.determineCaptureChannels(nLevel, iLevel, settings)
            push!( captureLines, DielectronicRecombination.CaptureLine(iLevel, nLevel, eEnergy, 0., 0.,
                                                                             EmProperty(0., 0.), EmProperty(0., 0.), pws) )
        end
    end
    return( captureLines )
end

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
            amps = DielectronicRecombination.determinePhotonChannels(fLevel, nLevel, settings)
            if  length(amps) == 0      continue    end
            push!( photonLines, DielectronicRecombination.PhotonLine(nLevel, fLevel, pEnergy,
                                                                           EmProperty(0., 0.), amps) )
        end
    end
    return( photonLines )
end
