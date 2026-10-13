
"""
`module  JAC.DecayYield`  
... a submodel of JAC that contains all methods for computing fluorescence and Auger yields for some level(s).
"""
module DecayYield


using Printf, JenaAtomicCalculator, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

"""
`struct  DecayYield.Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing fluorescence and Auger yields.

    + approach                 ::Basics.AbstractCascadeApproach  ... Cascade approach applied for the yields:
                                                    Basics.AverageSCA(), Basics.SCA() or Basics.UserMCA().
    + printBefore              ::Bool           ... True if the selected levels is to be printed before the actual computations start.
    + geant4                   ::Bool           ... True, if decay probabilities for the Geant code are to be generated and printed.
    + levelSelection           ::LevelSelection ... Specifies the selected levels, if any.
    + decayShells              ::Array{Shell,1} ... Additional shells to be included among Cascade.StepwiseDecayScheme's decayShells,
                                                    on top of those already occurring in the initial configuration(s); needed, for
                                                    instance, to communicate a spectator shell that is completely absent from a
                                                    double-core-hole initial configuration (e.g. "1s^2 2p^5" needs 2s added
                                                    explicitly so the cascade knows this shell can become populated during decay).
"""
struct Settings  <:  AbstractPropertySettings
    approach                   ::Basics.AbstractCascadeApproach
    printBefore                ::Bool
    geant4                     ::Bool
    levelSelection             ::LevelSelection
    decayShells                ::Array{Shell,1}
end


"""
`DecayYield.Settings()`
    ... constructor for an `empty` instance of DecayYield.Settings for the computation of fluorescence and Auger yields.
"""
function Settings()
    Settings(Basics.AverageSCA(), false, false, LevelSelection(), Shell[] )
end


"""
`DecayYield.Settings(set::DecayYield.Settings;`

        approach=.., printBefore=.., geant4=.., levelSelection=.., decayShells=..)

    ... keyword copy-constructor for re-defining selected values of a settings::DecayYield.Settings.
"""
function Settings(set::DecayYield.Settings;
        approach::Union{Nothing,Basics.AbstractCascadeApproach}=nothing,   printBefore::Union{Nothing,Bool}=nothing,
        geant4::Union{Nothing,Bool}=nothing,               levelSelection::Union{Nothing,LevelSelection}=nothing,
        decayShells::Union{Nothing,Array{Shell,1}}=nothing)
    if  isnothing(approach)         approachx       = set.approach       else   approachx       = approach       end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(geant4)           geant4x         = set.geant4         else   geant4x         = geant4         end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end
    if  isnothing(decayShells)      decayShellsx    = set.decayShells    else   decayShellsx    = decayShells    end

    Settings( approachx, printBeforex, geant4x, levelSelectionx, decayShellsx )
end


# `Base.show(io::IO, settings::DecayYield.Settings)`  ... prepares a proper printout of the variable settings::DecayYield.Settings.
function Base.show(io::IO, settings::DecayYield.Settings)
    println(io, "approach:                 $(settings.approach)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "geant4:                   $(settings.geant4)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
    println(io, "decayShells:              $(settings.decayShells)  ")
end


"""
`struct  DecayYield.Outcome`  
    ... defines a type to keep the outcome of a fluorescence and Auger yield computation as well other results.

    + level                     ::Level              ... Atomic level to which the outcome refers to.
    + NoPhotonLines             ::Int64              ... Number of photon decay lines
    + NoAugerLines              ::Int64              ... Number of Auger decay lines
    + rateR                     ::EmProperty         ... Total fluorescence rate
    + rateA                     ::Float64            ... Total Auger rate
    + omegaR                    ::EmProperty         ... fluorescence yield
    + omegaA                    ::EmProperty         ... Auger yield
"""
struct Outcome 
    level                       ::Level 
    NoPhotonLines               ::Int64
    NoAugerLines                ::Int64 
    rateR                       ::EmProperty
    rateA                       ::Float64
    omegaR                      ::EmProperty
    omegaA                      ::EmProperty 
end 


"""
`DecayYield.Outcome()`  
    ... constructor for an `empty` instance of DecayYield.Outcome for the computation of fluorescence and Auger yields.
"""
function Outcome()
    Outcome(Level(), 0, 0, EmProperty(0.), 0., EmProperty(0.), EmProperty(0.))
end


# `Base.show(io::IO, outcome::DecayYield.Outcome)`  ... prepares a proper printout of the variable outcome::DecayYield.Outcome.
function Base.show(io::IO, outcome::DecayYield.Outcome) 
    println(io, "level:                $(outcome.level)  ")
    println(io, "NoPhotonLines:        $(outcome.NoPhotonLines)  ")
    println(io, "NoAugerLines:         $(outcome.NoAugerLines)  ")
    println(io, "rateR:                $(outcome.rateR)  ")
    println(io, "rateA:                $(outcome.rateA)  ")
    println(io, "omegaR:               $(outcome.omegaR)  ")
    println(io, "omegaA:               $(outcome.omegaA)  ")
end


"""
`DecayYield.computeOutcomes(initialConfigs::Array{Configuration,1},asfSettings::AsfSettings,
                            multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings)`
    ... to compute (as selected) the fluorescence and Auger yields for the levels of the given multiplet and as specified
        by the given settings. The results are printed in neat tables to screen but nothing is returned otherwise.
"""
function computeOutcomes(initialConfigs::Array{Configuration,1}, asfSettings::AsfSettings,
                            multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::DecayYield.Settings; output=true)
    println("")
    printstyled("DecayYield.computeOutcomes(): The computation of the fluorescence & Auger yields starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    # `multiplet` (the caller's own, independent standalone SCF) is deliberately NOT used to build the
    # outcomes below -- see the cascadeLevels comment further down (FIXED 13-Oct-2026, item 65).
    #
    # Overwrite the default settings for the generation of the continuum orbitals, if needed
    ## setDefaults("method: continuum, asymptotic Coulomb")
    ## setDefaults("method: normalization, pure sine")
    #
    cApproach   = settings.approach
    decayShells = union( Basics.extractNonrelativisticShellList(initialConfigs), settings.decayShells )
    # NOTE: Cascade.StepwiseDecayScheme.decayShells (passed below) is NOT consulted by the currently-active
    # Basics.ForStepwiseDecay-based decay-configuration generator, which decides which shells may participate
    # purely via haskey(conf.shells, shell) on the CONFIGURATION itself (Basics.ValenceShells()). A shell that
    # is genuinely absent from the initial configuration (e.g. a double-core-hole's missing 2s in "1s 2p^6")
    # therefore must be injected here directly, as an explicit zero-occupation entry, or the cascade will
    # never consider it -- exactly mirroring what Cascade's own now-superseded (Oct-2025) decay-configuration
    # generator used to do internally.
    decayInitialConfigs = initialConfigs
    if  !isempty(settings.decayShells)
        decayInitialConfigs = Configuration[]
        for  conf  in  initialConfigs
            nshells = copy(conf.shells)
            for  shell  in  settings.decayShells    if  !haskey(nshells, shell)    nshells[shell] = 0    end    end
            push!(decayInitialConfigs, Configuration(nshells, conf.NoElectrons))
        end
    end
    #
    println("\nPerform a cascade computation for all radiative decay channels of the levels from the initial configurations:")
    wa = JenaAtomicCalculator.Cascade.Computation(JenaAtomicCalculator.Cascade.Computation(), name="photon lines", nuclearModel=nm, grid=grid, asfSettings=asfSettings,
                                scheme=JenaAtomicCalculator.Cascade.StepwiseDecayScheme([Radiative()], 0, Dict{Int64,Float64}(), 0, decayShells, Shell[], Shell[]),
                                approach=cApproach, initialConfigs=decayInitialConfigs)
    wb = perform(wa, output=true, outputToFile=false);   linesR = wb[PhotoEmission.Settings]
    # The levels DecayYield reports on must be identified from THIS SAME cascade computation, not from an
    # independent standalone SCF: linesR/linesA's line.initialLevel energies are whatever `cApproach` actually
    # produced for "generation 0" (e.g. Basics.SCA()'s ChargeStateOrbitals() shares ONE mean field, built by
    # generateMeanFieldOrbitals, across every configuration of the same electron count anywhere in the
    # cascade -- a deliberate speed/accuracy trade-off that does not, in general, reproduce a configuration's
    # own dedicated SCF energy). Matching an EXTERNALLY-computed `outcome.level` against those lines by a
    # fixed relative-energy tolerance (Cascade.computeDecayYieldOutcome's Basics.isSimilar, tol 1e-3) silently
    # failed for exactly this reason (FIXED 13-Oct-2026, item 65): measured 2.1e-3 apart on a Mg K-hole, over
    # twice the tolerance, for nominally "the same" level.
    genMultiplets = wb[ResultKeys.GeneratedMultiplets]
    cascadeLevels = DecayYield.extractGenerationZeroLevels(genMultiplets, decayInitialConfigs)
    #
    outcomes = DecayYield.determineOutcomes(cascadeLevels, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    DecayYield.displayOutcomes(outcomes)    end
    #
    println("\nPerform a cascade computation for all (single-electron) Auger decay channels of the levels from the initial configurations:")
    wa = JenaAtomicCalculator.Cascade.Computation(JenaAtomicCalculator.Cascade.Computation(), name="Auger lines", nuclearModel=nm, grid=grid, asfSettings=asfSettings,
                                scheme=JenaAtomicCalculator.Cascade.StepwiseDecayScheme([Auger()], 1, Dict{Int64,Float64}(), 0, decayShells, Shell[], Shell[]),
                                approach=cApproach, initialConfigs=decayInitialConfigs)
    wb = perform(wa, output=true, outputToFile=false);   linesA = wb[AutoIonization.Settings]
    #
    # Calculate all amplitudes and requested properties
    newOutcomes = DecayYield.Outcome[]
    for  outcome in outcomes
        newOutcome = JenaAtomicCalculator.Cascade.computeDecayYieldOutcome(outcome, linesR, linesA, settings)
        push!( newOutcomes, newOutcome)
        #
        if settings.geant4
            JenaAtomicCalculator.Cascade.computeDecayProbabilities(newOutcome, linesR, linesA, settings)
        end
    end
    #
    # Print all results to screen
    DecayYield.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    DecayYield.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`DecayYield.extractGenerationZeroLevels(genMultiplets::Array{Multiplet,1}, decayInitialConfigs::Array{Configuration,1})`
    ... to pick out, from a Cascade.StepwiseDecayScheme's own `generated multiplets:` result, the multiplet(s) that
        represent "generation 0" -- the initial configuration(s) themselves, computed under the SAME mean-field
        approach (`cApproach`) that also produced the radiative and Auger lines these levels will later be matched
        against. `Cascade.generateBlocks` appends the initial configurations to the end of its own configuration
        list and builds one block per configuration in that same order, so the last `length(decayInitialConfigs)`
        entries of `genMultiplets` are exactly the generation-0 blocks, in order; this is verified rather than
        assumed, by comparing each block's own extracted configuration against the one it is expected to match.
        A combined `multiplet::Multiplet` of all their levels is returned.
"""
function  extractGenerationZeroLevels(genMultiplets::Array{Multiplet,1}, decayInitialConfigs::Array{Configuration,1})
    n = length(decayInitialConfigs)
    if  length(genMultiplets) < n
        error("DecayYield.extractGenerationZeroLevels(): $(length(genMultiplets)) generated multiplets, fewer than "  *
              "the $n initial configuration(s); the cascade did not generate a generation-0 block for each of them.")
    end
    genZero = genMultiplets[end-n+1:end]
    levels  = Level[]
    for  (mp, conf)  in  zip(genZero, decayInitialConfigs)
        mpConfs = Basics.extractConfigurations(Basics.FromMultiplet(), mp)
        if  length(mpConfs) != 1  ||  mpConfs[1] != conf
            error("DecayYield.extractGenerationZeroLevels(): expected the cascade's generation-0 block to be the "  *
                  "single configuration $conf, but its multiplet's own configuration(s) are $mpConfs instead; the "  *
                  "generation-0 block could no longer be identified by position and must not be guessed at.")
        end
        append!(levels, mp.levels)
    end
    return( Multiplet("generation 0 (cascade-internal)", levels) )
end


"""
`DecayYield.determineOutcomes(multiplet::Multiplet, settings::DecayYield.Settings)`
    ... to determine a list of Outcomes's for the fluorescence and Auger yields of the levels of the given
        multiplet. It takes into account the particular selections and settings. An Array{DecayYield.Outcome,1} is
        returned. Apart from the level specification, all physical properties are set to zero during the
        initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::DecayYield.Settings)
    outcomes = DecayYield.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, DecayYield.Outcome(level, 0, 0, EmProperty(0.), 0., EmProperty(0.), EmProperty(0.)) )
        end
    end
    return( outcomes )
end


"""
`DecayYield.displayOutcomes(outcomes::Array{DecayYield.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{DecayYield.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected DecayYield levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`DecayYield.displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1}, settings::DecayYield.Settings)`  
    ... to display the energies, M_ms and F-parameters, etc. for the selected levels. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{DecayYield.Outcome,1}, settings::DecayYield.Settings)
    nx = 149
    println(stream, " ")
    println(stream, "  Fluorescence and Auger decay yields:")
    println(stream, " ")
    println(stream, "    + No_R, rate_R, omega_R  ... number of fluorescence lines, total fluorescence rate and yield. ")
    println(stream, "    + No_A, rate_A, omega_A  ... number of Auger lines, total Auger rate and Auger yield. ")
    if    settings.approach isa Union{Basics.AverageSCA, Basics.SCA}
    println(stream, "    + Approach:  $(settings.approach)  ... all fluorescence rates and yields only in Babushkin gauge. ")
    end
    #
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=3)              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center( 9, "No_R"    ; na=0);                          sb = sb * TableStrings.hBlank(9)
    sa = sa * TableStrings.center(24, "Cou -- rate_R  -- Bab"; na=1);             sb = sb * TableStrings.hBlank(25)
    sa = sa * TableStrings.center(24, "Cou -- omega_R -- Bab"; na=3);             sb = sb * TableStrings.hBlank(27)
    sa = sa * TableStrings.center( 6, "No_A"    ; na=3);                          sb = sb * TableStrings.hBlank(9)
    sa = sa * TableStrings.center(10, "rate_A"; na=2);                            sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(24, "Cou -- omega_A -- Bab"; na=3);             sb = sb * TableStrings.hBlank(27)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.level.energy
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))    * "     "
        sa = sa * string(outcome.NoPhotonLines, pad=4)                                      * "    "
        sa = sa * @sprintf("%.4e", outcome.rateR.Coulomb)                                   * " "
        sa = sa * @sprintf("%.4e", outcome.rateR.Babushkin)                                 * "    "
        sa = sa * @sprintf("%.4e", outcome.omegaR.Coulomb)                                  * " "
        sa = sa * @sprintf("%.4e", outcome.omegaR.Babushkin)                                * "     "
        sa = sa * string(outcome.NoAugerLines, pad=4)                                       * "    "
        sa = sa * @sprintf("%.4e", outcome.rateA)                                           * "    "
        sa = sa * @sprintf("%.4e", outcome.omegaA.Coulomb)                                  * " "
        sa = sa * @sprintf("%.4e", outcome.omegaA.Babushkin)                                * "  "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n")
    #
    return( nothing )
end


end # module
