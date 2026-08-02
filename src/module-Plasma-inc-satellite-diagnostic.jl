
# Functions and methods for scheme::Plasma.SatelliteDiagnosticScheme computations


"""
`Plasma.reconcileTemperatures(drTemps::Array{Float64,1}, ieTemps::Array{Float64,1})`
    ... reconciles the two independent temperature grids of scheme.drSettings and scheme.ieSettings into
        a single, shared grid [K] for the R(Te) ratio. If one is empty, the other is used. If both are
        non-empty and differ, the longer one is used and a warning is printed to the summary/log stream.
        If both are empty, an error is raised. An Array{Float64,1} is returned.
"""
function reconcileTemperatures(drTemps::Array{Float64,1}, ieTemps::Array{Float64,1})
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if      isempty(drTemps)  &&  isempty(ieTemps)
        error("Plasma.SatelliteDiagnosticScheme: both drSettings.temperatures and ieSettings.temperatures " *
              "are empty; a shared temperature grid [K] must be given in at least one of them.")
    elseif  isempty(drTemps)                        return( ieTemps )
    elseif  isempty(ieTemps)                        return( drTemps )
    elseif  drTemps == ieTemps                      return( drTemps )
    else
        sa = "Plasma.SatelliteDiagnosticScheme: drSettings.temperatures and ieSettings.temperatures differ; " *
             "using the longer of the two ($(length(drTemps) >= length(ieTemps) ? "drSettings" : "ieSettings"))."
        @warn(sa);    println(sa)
        if  printSummary   println(iostream, sa)   end
        return( length(drTemps) >= length(ieTemps) ? drTemps : ieTemps )
    end
end


"""
`Plasma.perform(scheme::Plasma.SatelliteDiagnosticScheme, computation::Plasma.Computation; output::Bool=true)`
    ... to compute the dielectronic-satellite-to-parent-line intensity-ratio Te diagnostic R(Te) for one
        recombining ion. For output=true, a dictionary is returned from which the relevant results can be
        accessed by proper keys.
"""
function  perform(scheme::Plasma.SatelliteDiagnosticScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end

    nm  = computation.nuclearModel
    grid = computation.grid

    temperatures = Plasma.reconcileTemperatures(scheme.drSettings.temperatures, scheme.ieSettings.temperatures)

    # -------------------------------------------------------------------------------------------------------
    # (1) Satellites:  DR capture fromShells --> toShells, + intoShells, stabilizing via decayShells.
    #     scheme.intoShells is treated TOGETHER, exactly as Basics.ForDielectronicRecombination already does.
    # -------------------------------------------------------------------------------------------------------
    theme = Basics.ForDielectronicRecombination(scheme.fromShells, scheme.toShells, scheme.intoShells, scheme.decayShells)
    (intermediateConfs, finalConfs) = Basics.generateConfigurations(theme, computation.refConfigs)

    # DielectronicRecombination.Settings()'s own default has gauges=UseGauge[] (EMPTY) -- with no gauge
    # selected, computePathways silently gives photonRate=EmProperty(0.,0.) for every resonance (verified
    # today: reproducibly zero, both gauges, regardless of grid size/asfSettings/printBefore, none of
    # which were the actual cause). Only override gauges here if the caller left it empty, so an explicit
    # choice in scheme.drSettings is still respected.
    gauges     = isempty(scheme.drSettings.gauges) ? [UseCoulomb, UseBabushkin] : scheme.drSettings.gauges
    drSettings = DielectronicRecombination.Settings(scheme.drSettings; calcRateAlpha=true, temperatures=temperatures,
                                                    gauges=gauges)

    drComp  = Atomic.Computation(Atomic.Computation(), name="SatelliteDiagnosticScheme: DR satellites", grid=grid,
                                 nuclearModel=nm, initialConfigs=computation.refConfigs,
                                 intermediateConfigs=intermediateConfs, finalConfigs=finalConfs,
                                 asfSettings=computation.asfSettings, processSettings=drSettings )
    drResults  = Basics.perform(drComp; output=true)
    pathways   = drResults["dielectronic recombination pathways:"]
    resonances = DielectronicRecombination.computeResonances(pathways, drSettings)

    # A resonance with BOTH augerRate==0 and photonRate==0 (no decay channel at all within the computed
    # scope -- a real, if unusual, edge case, distinct from a resonance with only one channel zero) gives
    # a genuine 0/0 = NaN resonanceStrength; since NaN propagates through addition, a single such
    # resonance would otherwise poison the entire satelliteAlpha sum. Skipped here, not folded in as 0,
    # so it is not silently misrepresented as "computed and negligible".
    satelliteAlpha = Dict{Float64, EmProperty}( t => EmProperty(0., 0.) for t in temperatures )
    for  resonance in resonances
        if  !isfinite(resonance.resonanceStrength.Coulomb)  ||  !isfinite(resonance.resonanceStrength.Babushkin)
            continue
        end
        for  t in temperatures
            satelliteAlpha[t] = satelliteAlpha[t] + DielectronicRecombination.computeRateCoefficient(resonance, t)
        end
    end

    if output    results = Base.merge( results, Dict("satellite pathways:"    => pathways,
                                                       "satellite resonances:" => resonances,
                                                       "satellite alpha_DR:"   => satelliteAlpha) )   end

    # -------------------------------------------------------------------------------------------------------
    # (2) Parent line:  the SAME core promotion fromShells-->toShells, with NO spectator electron.
    # -------------------------------------------------------------------------------------------------------
    bareExcitedConfs = Basics.generateConfigurations(ExciteElectrons(1, scheme.fromShells, scheme.toShells), computation.refConfigs)
    bareExcitedConfs = Basics.extractConfigurations(ContractShells(), bareExcitedConfs)

    # (2a) Radiative decay of the bare excited level, back to refConfigs. Routed through
    #      Atomic.Computation/Basics.perform, exactly like the DR and IE steps, for the same proven
    #      orbital-/basis-reconciliation between independently-generated multiplets.
    #
    #      NO competing autoionization channel is computed here, deliberately: autoionization requires
    #      a final state with ONE FEWER bound electron (the ion after ejecting one), but the bare
    #      excited level (fromShells-->toShells) has the SAME electron count as computation.refConfigs
    #      -- an ordinary valence excitation, not an ionizing process. A genuine autoionization channel
    #      only opens up when EI excites an INNER-SHELL electron, placing the excited level in the
    #      continuum of the NEXT-lower charge state (a different, much less commonly studied DR regime
    #      than the valence-shell case this scheme targets) -- an earlier attempt to compute this via
    #      AutoIonization.computeLines(bareExcitedConfs, refConfigs, ...) was dimensionally wrong for
    #      exactly this reason (same electron count on both sides) and gave a meaningless empty-channel
    #      result, not a genuine parity/angular-momentum forbiddenness.
    #      Branching fraction is therefore fixed at 1 (pure radiative) -- valid for the common,
    #      non-inner-shell case this scheme is designed for.
    @warn("Plasma.SatelliteDiagnosticScheme: the parent line's upper level is assumed to decay purely " *
          "radiatively (branching fraction = 1); autoionizing electron-impact excitations (which occur " *
          "only when an INNER-SHELL electron is excited) are not covered by this scheme.")

    ieLineSelection = scheme.ieSettings.lineSelection
    lineSelection    = LineSelection(ieLineSelection.active;
                                     indexPairs=[ (f, i) for (i, f) in ieLineSelection.indexPairs ])

    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(); lineSelection=lineSelection)
    photoComp     = Atomic.Computation(Atomic.Computation(), name="SatelliteDiagnosticScheme: parent line (radiative)",
                                       grid=grid, nuclearModel=nm, initialConfigs=bareExcitedConfs,
                                       finalConfigs=computation.refConfigs, asfSettings=computation.asfSettings,
                                       processSettings=photoSettings )
    photoResults  = Basics.perform(photoComp; output=true)
    photoLines    = photoResults["radiative lines:"]

    branchingFraction = EmProperty(1., 1.)

    # (2b) Population of the bare excited level via electron-impact excitation.
    ieSettings = ImpactExcitation.Settings(scheme.ieSettings; calcRateCoefficient=true, temperatures=temperatures,
                                           printBefore=false)   # false always -- module-ImpactExcitation.jl displayLines() bug

    ieComp    = Atomic.Computation(Atomic.Computation(), name="SatelliteDiagnosticScheme: parent-line excitation", grid=grid,
                                   nuclearModel=nm, initialConfigs=computation.refConfigs, finalConfigs=bareExcitedConfs,
                                   asfSettings=computation.asfSettings, processSettings=ieSettings )
    ieResults = Basics.perform(ieComp; output=true)
    ieLines, ieRates = ieResults["impact-excitation lines:"]

    parentAlpha = Dict{Float64, Float64}( rc.temperatures[i] => rc.alphas[i] for rc in ieRates for i in eachindex(rc.temperatures) )

    if output    results = Base.merge( results, Dict("parent radiative lines:"     => photoLines,
                                                       "parent branching fraction:" => branchingFraction,
                                                       "parent excitation lines:"   => ieLines,
                                                       "parent excitation alpha:"   => parentAlpha) )   end

    # -------------------------------------------------------------------------------------------------------
    # (3) Combine:  R(Te) = satellite alpha_DR / (parent excitation alpha * parent radiative branching fraction)
    # -------------------------------------------------------------------------------------------------------
    ratios = Dict{Float64, EmProperty}()
    for  t in temperatures
        if  haskey(parentAlpha, t)  &&  isfinite(parentAlpha[t])  &&  parentAlpha[t] > 0.  &&
            isfinite(satelliteAlpha[t].Coulomb)  &&  isfinite(satelliteAlpha[t].Babushkin)
            denom = parentAlpha[t] * branchingFraction
            ratios[t] = satelliteAlpha[t] / denom
        end
    end

    Plasma.displaySatelliteDiagnosticResults(stdout, temperatures, satelliteAlpha, parentAlpha, branchingFraction, ratios)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displaySatelliteDiagnosticResults(iostream, temperatures, satelliteAlpha, parentAlpha,
                                                                 branchingFraction, ratios)   end

    if output    results = Base.merge( results, Dict("R(Te):" => ratios) )   end

    println("SatelliteDiagnosticScheme computation complete ...")
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end


"""
`Plasma.displaySatelliteDiagnosticResults(stream::IO, temperatures, satelliteAlpha, parentAlpha, branchingFraction, ratios)`
    ... prints a table of the satellite/parent rate coefficients and the combined R(Te) ratio.
"""
function  displaySatelliteDiagnosticResults(stream::IO, temperatures::Array{Float64,1},
                                            satelliteAlpha::Dict{Float64,EmProperty}, parentAlpha::Dict{Float64,Float64},
                                            branchingFraction::EmProperty, ratios::Dict{Float64,EmProperty})
    println(stream, "\n  Satellite-to-parent-line intensity-ratio Te diagnostic:")
    println(stream, "  parent radiative branching fraction (Cou/Bab) = $(branchingFraction)")
    println(stream, "  ", TableStrings.hLine(90))
    println(stream, "    Te [K]        alpha_DR (satellites)         alpha_exc (parent)            R(Te) = ratio")
    println(stream, "  ", TableStrings.hLine(90))
    for  t in temperatures
        aDR  = get(satelliteAlpha, t, EmProperty(0., 0.))
        aExc = get(parentAlpha, t, 0.)
        r    = get(ratios, t, EmProperty(0., 0.))
        println(stream, "    $t        $aDR        $aExc        $r")
    end
    println(stream, "  ", TableStrings.hLine(90))
    return( nothing )
end
