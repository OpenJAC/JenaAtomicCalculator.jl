
# Functions and methods for scheme::Plasma.CollisionalRadiativeScheme computations


"""
`Plasma.generateCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)`
    ... generates the (small) level set for computation.refConfigs, following the same generation policy
        as Plasma.generateIonLevelData (single excitations from all shells of refConfigs; double & higher
        excitations, if scheme.NoExcitations > 1, from the same shells -- refConfigs is given explicitly by
        the caller here, so no Plasma.determineValenceShells electron-count heuristic is needed or used).
        All configurations share ONE mean-field orbital basis, and every level is then extracted via
        Hamiltonian.performCIwithFrozenOrbitals against that shared basis -- this is what keeps the levels
        properly reconciled for the later PhotoEmission/ImpactExcitation calls (no separate SCF per config).
        A Multiplet, sorted by energy with levels re-indexed 1..N (Basics.sortByEnergy), is returned.
"""
function generateCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)
    nm = computation.nuclearModel;   grid = computation.grid;   refConfigs = computation.refConfigs

    fromShells = Basics.extractShellList(refConfigs)
    toShells   = Basics.generateShellList(refConfigs, scheme.upperShellNo, [0,1,2,3])

    sConfigs  = Configuration[];    dtConfigs = Configuration[]
    if  scheme.NoExcitations >= 1    sConfigs  = Basics.generateConfigurations(refConfigs, fromShells, toShells, 1)                  end
    if  scheme.NoExcitations >  1    dtConfigs = Basics.generateConfigurations(refConfigs, fromShells, toShells, scheme.NoExcitations)  end
    allConfigs = unique( append!(sConfigs, dtConfigs) )

    println(">> Plasma.CollisionalRadiativeScheme: generating the level set from configurations \n   $allConfigs")

    name       = "CollisionalRadiativeScheme mean multiplet"
    mfSettings = AtomicState.MeanFieldSettings()
    repBasis   = Representation(name, nm, grid, allConfigs, MeanFieldBasis(mfSettings))
    repOutput  = generate(repBasis, output=true)
    orbitals   = repOutput["mean-field basis"].orbitals

    levels = Level[]
    for  conf in allConfigs
        mp = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, nm, grid, computation.asfSettings; printout=false)
        append!(levels, mp.levels)
    end

    repMultiplet = Basics.sortByEnergy( Multiplet(name, levels) )
    Basics.displayLevels(stdout, [repMultiplet], N=100)
    println(">> Plasma.CollisionalRadiativeScheme: generated $(length(repMultiplet.levels)) levels.")

    return( repMultiplet )
end


"""
`Plasma.collisionalRadiativeLevelFingerprint(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)`
    ... builds a small Dict{String,Any} of PLAIN, comparable values (Float64/Int64/String only, never whole
        structs) that together identify the physics inputs of Plasma.generateCollisionalRadiativeLevels for
        the given scheme/computation: element, charge state, reference configuration(s), NoExcitations,
        upperShellNo, the grid's own recipe scalars (rnt, h, hp, NoPoints -- NOT its large r/rp/wr arrays,
        which would risk truncated/misleading string comparisons), and asfSettings (safe to stringify whole,
        since none of its own fields are large arrays, unlike Radial.Grid). Used both when writing a
        level-cache file (stored alongside the data) and when validating a candidate cache file for reuse --
        an EXACT match on every key is required, not a sufficiency/`>=`-style check.
"""
function collisionalRadiativeLevelFingerprint(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)
    nm  = computation.nuclearModel;   grid = computation.grid;   asf = computation.asfSettings
    # Base.string(::AsfSettings) is deliberately unimplemented in JAC ("Not yet implemented."), so the
    # struct cannot be fingerprinted as a whole; string(typeof(...)) is used for its enum-like interaction/
    # QED fields (safe and generic -- works for any JAC type without needing a custom Base.string method).
    # jjLS/levelSelectionCI/shellSequenceScf/frozenSubshells are deliberately NOT included -- a documented
    # simplification, not a guarantee that every AsfSettings field is covered.
    return  Dict{String,Any}(
        "Z"             => trunc(nm.Z, digits=3),
        "A"             => trunc(nm.mass, digits=3),
        "electronCount" => sum( c.NoElectrons for c in computation.refConfigs ),
        "refConfigs"    => join( sort(string.(computation.refConfigs)), ";" ),
        "NoExcitations" => scheme.NoExcitations,
        "upperShellNo"  => scheme.upperShellNo,
        "grid_rnt"      => grid.rnt,   "grid_h" => grid.h,   "grid_hp" => grid.hp,   "grid_NoPoints" => grid.NoPoints,
        "asf_eeInteractionCI" => string(typeof(asf.eeInteractionCI)),
        "asf_qedModel"        => string(typeof(asf.qedModel)),
        "asf_accuracyScf"     => asf.accuracyScf,
        "asf_maxIterationsScf" => asf.maxIterationsScf,
        # Version tag for the level-generation ALGORITHM itself, which no settings-derived key can capture:
        # the settings can be identical while the code that turns them into levels has changed, in which case
        # every key above matches and a stale cache is silently reused -- the fix then looks like a no-op.
        # This bit us once: until 04-Aug-2026 Hamiltonian.performCIwithFrozenOrbitals ignored eeInteractionCI
        # and returned single-CSF levels, so "asf_eeInteractionCI" read CoulombInteraction both before and
        # after the fix. BUMP THIS STRING whenever generateCollisionalRadiativeLevels changes how it builds
        # levels; that invalidates every existing cache file, which is exactly what should happen.
        "generationMethod"    => "performCI-with-frozen-orbitals-2026-08-04" )
end


"""
`Plasma.readEvaluateCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)`
    ... tries each file in scheme.levelsFilenames, in order, and returns the cached Multiplet from the first
        one whose stored fingerprint EXACTLY matches collisionalRadiativeLevelFingerprint(scheme,
        computation) -- not a sufficiency/`>=`-style check, unlike Plasma.readEvaluateIonLevelData, since
        most of the relevant parameters here (grid, asfSettings, ...) have no natural "richer implies
        superset" ordering. missing is returned if no candidate matches (or scheme.levelsFilenames is
        empty), so the caller can fall back to a fresh Plasma.generateCollisionalRadiativeLevels call.
"""
function readEvaluateCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation)
    fingerprint = Plasma.collisionalRadiativeLevelFingerprint(scheme, computation)
    for  filename in scheme.levelsFilenames
        if  !isfile(filename)
            println(">>> Plasma.CollisionalRadiativeScheme: level-cache file $filename not found; skipped.")
            continue
        end
        data = JLD2.load(filename)
        if  data["fingerprint"] == fingerprint
            println(">>> Plasma.CollisionalRadiativeScheme: reusing level set from $filename.")
            return( data["repMultiplet"] )
        else
            println(">>> Plasma.CollisionalRadiativeScheme: $filename does not match the current level " *
                    "settings; skipped.")
        end
    end
    return( missing )
end


"""
`Plasma.writeCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                        repMultiplet::Multiplet)`
    ... writes repMultiplet, together with its collisionalRadiativeLevelFingerprint(scheme, computation), to a
        new, auto-named, never-overwritten file inside scheme.cacheDirectory (created if missing; ""
        writes into the current directory). The filename is printed so it can be adopted into a later run's
        scheme.levelsFilenames.
"""
function writeCollisionalRadiativeLevels(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                         repMultiplet::Multiplet)
    if  scheme.cacheDirectory != ""    mkpath(scheme.cacheDirectory)    end
    nm = computation.nuclearModel
    sa = string(round(Int, nm.Z));   sb = string(round(Int, nm.mass))
    filename    = joinpath(scheme.cacheDirectory, "newCRLevelsZ" * sa * "A" * sb * "-" * string(Dates.now())[1:13] * ".jld")
    fingerprint = Plasma.collisionalRadiativeLevelFingerprint(scheme, computation)
    JLD2.@save filename fingerprint repMultiplet
    println(">>> Plasma.CollisionalRadiativeScheme: level set written to $filename -- add this filename to " *
            "levelsFilenames in a later run to reuse it.")
    return( filename )
end


"""
`Plasma.collisionalRadiativeRatesFingerprint(scheme::Plasma.CollisionalRadiativeScheme, levelFingerprint::Dict{String,Any})`
    ... extends the given levelFingerprint with the ImpactExcitation-specific settings that also determine
        the raw collision data (maxKappa, maxEnergyMultiplier, numElectronEnergies, operator) -- copying in
        the level fingerprint keeps the rates-cache file self-contained, so validating it does not require
        separately cross-checking a level-cache file. Deliberately does NOT include ieSettings.temperatures:
        the cached data (ImpactExcitation's raw, pre-thermal-average lines) is temperature-independent by
        construction, so any temperature list can be evaluated from it afterward at negligible cost.
"""
function collisionalRadiativeRatesFingerprint(scheme::Plasma.CollisionalRadiativeScheme, levelFingerprint::Dict{String,Any})
    fp = deepcopy(levelFingerprint)
    fp["maxKappa"]            = scheme.ieSettings.maxKappa
    fp["maxEnergyMultiplier"] = scheme.ieSettings.maxEnergyMultiplier
    fp["numElectronEnergies"] = scheme.ieSettings.numElectronEnergies
    fp["operator"]            = string(typeof(scheme.ieSettings.operator))
    return( fp )
end


"""
`Plasma.readEvaluateCollisionalRadiativeRates(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                              levelFingerprint::Dict{String,Any})`
    ... tries each file in scheme.ratesFilenames, in order, and returns the cached raw ImpactExcitation
        lines from the first one whose stored fingerprint EXACTLY matches
        collisionalRadiativeRatesFingerprint(scheme, levelFingerprint). missing is returned if no candidate
        matches (or scheme.ratesFilenames is empty), so the caller can fall back to a fresh
        ImpactExcitation.computeLines call.
"""
function readEvaluateCollisionalRadiativeRates(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                               levelFingerprint::Dict{String,Any})
    fingerprint = Plasma.collisionalRadiativeRatesFingerprint(scheme, levelFingerprint)
    for  filename in scheme.ratesFilenames
        if  !isfile(filename)
            println(">>> Plasma.CollisionalRadiativeScheme: rates-cache file $filename not found; skipped.")
            continue
        end
        data = JLD2.load(filename)
        if  data["fingerprint"] == fingerprint
            println(">>> Plasma.CollisionalRadiativeScheme: reusing ImpactExcitation collision data from $filename.")
            return( data["ieLines"] )
        else
            println(">>> Plasma.CollisionalRadiativeScheme: $filename does not match the current level/rate " *
                    "settings; skipped.")
        end
    end
    return( missing )
end


"""
`Plasma.writeCollisionalRadiativeRates(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                       levelFingerprint::Dict{String,Any}, ieLines::Array{ImpactExcitation.Line,1})`
    ... writes ieLines (ImpactExcitation's raw, pre-thermal-average lines), together with its
        collisionalRadiativeRatesFingerprint(scheme, levelFingerprint), to a new, auto-named,
        never-overwritten file inside scheme.cacheDirectory. The filename is printed so it can be adopted
        into a later run's scheme.ratesFilenames.
"""
function writeCollisionalRadiativeRates(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation,
                                        levelFingerprint::Dict{String,Any}, ieLines::Array{ImpactExcitation.Line,1})
    if  scheme.cacheDirectory != ""    mkpath(scheme.cacheDirectory)    end
    nm = computation.nuclearModel
    sa = string(round(Int, nm.Z));   sb = string(round(Int, nm.mass))
    filename    = joinpath(scheme.cacheDirectory, "newCRRatesZ" * sa * "A" * sb * "-" * string(Dates.now())[1:13] * ".jld")
    fingerprint = Plasma.collisionalRadiativeRatesFingerprint(scheme, levelFingerprint)
    JLD2.@save filename fingerprint ieLines
    println(">>> Plasma.CollisionalRadiativeScheme: ImpactExcitation collision data written to $filename -- add " *
            "this filename to ratesFilenames in a later run to reuse it (any scheme.ieSettings.temperatures " *
            "list will work, no recomputation needed).")
    return( filename )
end


"""
`Plasma.estimateCollisionalRadiativeCost(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation;`

        referencePair::Tuple{Int64,Int64}=(1,2))

    ... estimates, WITHOUT running ImpactExcitation for every level pair, whether a Plasma.
        CollisionalRadiativeScheme computation is tractable. Generates the real level set (mean-field basis
        + per-configuration CI, always cheap) to get the real level count N, times ONE representative
        ImpactExcitation pair (referencePair, by energy-sorted index) at the scheme's own scheme.ieSettings,
        and extrapolates the total collisional-rate cost via N(N-1)/2 pair-count scaling -- the dominant
        cost for any level set beyond a handful of levels (found 02-Aug-2026: a 45-level/990-pair He-like
        carbon case extrapolated to ~7h; a 1879-level/1.76M-pair Ar3+ case -- single excitations from ALL
        shells of a 3s^2 3p^3 reference configuration, including the closed 1s/2s/2p core -- extrapolated
        to ~17000h, genuinely intractable rather than merely slow). This is a ROUGH order-of-magnitude
        estimate, not a precise prediction: per-pair cost varies with channel count and CSF-basis size, and
        only one representative pair is timed here. Meant to be checked BEFORE attempting a full run, not
        as a replacement for one. A named tuple (levels, pairs, levelGenerationTime, referencePairTime,
        estimatedTotalIETime) is returned; a human-readable summary is always printed.

    + scheme        ::Plasma.CollisionalRadiativeScheme
        ... the scheme (NoExcitations, upperShellNo, ieSettings) whose cost is to be estimated.
    + computation   ::Plasma.Computation   ... provides nuclearModel, grid, refConfigs, asfSettings.
    + referencePair ::Tuple{Int64,Int64}   ... which (by-energy-sorted) level indices to time; default (1,2).
"""
function estimateCollisionalRadiativeCost(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation;
                                          referencePair::Tuple{Int64,Int64}=(1,2))
    nm = computation.nuclearModel;   grid = computation.grid

    levelGenerationTime = @elapsed  repMultiplet = Plasma.generateCollisionalRadiativeLevels(scheme, computation)
    N     = length(repMultiplet.levels)
    pairs = div(N*(N-1), 2)

    if  N < maximum(referencePair)
        error("Plasma.estimateCollisionalRadiativeCost: referencePair=$referencePair is out of range for the " *
              "generated $N-level set; choose indices <= $N (or leave the default).")
    end

    ieSettingsRef = ImpactExcitation.Settings(scheme.ieSettings; calcRateCoefficient=true, printBefore=false,
                                              lineSelection=LineSelection(true; indexPairs=[referencePair]))
    ImpactExcitation.computeLines(repMultiplet, repMultiplet, nm, grid, ieSettingsRef; output=true)    # warm-up / compile
    referencePairTime    = @elapsed  ImpactExcitation.computeLines(repMultiplet, repMultiplet, nm, grid, ieSettingsRef; output=true)
    estimatedTotalIETime = pairs * referencePairTime

    println("\n  Plasma.CollisionalRadiativeScheme cost estimate:")
    println("  ", TableStrings.hLine(80))
    println("    generated levels                          N        = $N")
    println("    level pairs                               N(N-1)/2 = $pairs")
    println("    level-generation time (real)                       = $(round(levelGenerationTime, digits=2)) s")
    println("    reference pair $referencePair ImpactExcitation time (real) = $(round(referencePairTime, digits=3)) s")
    println("    extrapolated total ImpactExcitation time           = $(round(estimatedTotalIETime, digits=1)) s " *
            "= $(round(estimatedTotalIETime/60, digits=2)) min = $(round(estimatedTotalIETime/3600, digits=2)) h")
    if  estimatedTotalIETime > 3600.
        println("    >>> WARNING: estimated cost exceeds 1 hour -- consider a smaller level set (fewer shells, " *
                "lower NoExcitations/upperShellNo, or restricting which shells may be singly excited) before " *
                "attempting a full run.")
    end
    println("  ", TableStrings.hLine(80))

    return( (levels=N, pairs=pairs, levelGenerationTime=levelGenerationTime, referencePairTime=referencePairTime,
             estimatedTotalIETime=estimatedTotalIETime) )
end


"""
`Plasma.solveCollisionalRadiativeBalance(matrixA::Array{Float64,2})`
    ... solves the steady-state collisional-radiative balance dn_i/dt = 0 for the given N x N rate matrix
        (matrixA[i,j] = total rate FROM level i TO level j for i != j; matrixA[i,i] = -(total loss rate
        out of level i)). Since the rows of matrixA sum to zero by construction (population conservation),
        this system is singular; exactly one of the N balance equations is redundant and is replaced by the
        normalization constraint sum(n) = 1 to obtain a well-posed, uniquely solvable system. The equation
        that gets sacrificed is the GROUND STATE's (index 1, the lowest level, since repMultiplet is sorted
        ascending in energy) -- see the discussion below. A relative-population vector (summing to 1) is
        returned.
"""
function solveCollisionalRadiativeBalance(matrixA::Array{Float64,2})
    N = size(matrixA, 1)
    fullMatrix = collect(transpose(matrixA))
    ##
    ## Which equation to sacrifice is NOT arbitrary in floating-point arithmetic, even though all N choices
    ## are equivalent in exact arithmetic. Whichever level's equation is replaced loses its own balance
    ## relation and is then determined only indirectly, through its appearance as a source term in the other
    ## equations -- i.e. through a near-singular combination that amplifies roundoff. Replacing the LAST
    ## equation is therefore the worst possible choice: level N is the highest-lying, least-populated level,
    ## and at low temperatures it can be an inner-shell-excited level whose every feeding rate has underflowed
    ## to exact 0.0, so that its true population is exactly zero. The indirect determination then returns
    ## roundoff-scale noise of arbitrary sign -- the source of the (physically impossible) tiny NEGATIVE
    ## populations, e.g. -2.89e-27, seen for the decoupled levels 11-18 of boron-like Ne5+ (example-Je.jl,
    ## branch b) at low Te. Sacrificing the GROUND STATE's equation instead is strictly better on both counts:
    ## it always carries an O(1) share of the total population (not necessarily the single largest share --
    ## for Ne5+ the 2p_3/2 level 2 holds ~0.62 against the ground state's ~0.38 -- but never a roundoff-scale
    ## one), so it is well determined by the normalization constraint alone, and every weakly-populated level
    ## keeps its own, well-conditioned balance equation (for a decoupled level j that equation reduces to
    ## -loss_j * n_j = 0, which returns a clean zero instead of signed noise).
    fullMatrix[1, :] .= 1.0
    rhs               = zeros(N);   rhs[1] = 1.0
    populations       = fullMatrix \ rhs
    ##
    ## Residual cleanup, applied to NEGATIVE components ONLY. Positive populations are never touched, however
    ## small: in a coronal-limit balance the genuinely-populated levels span many orders of magnitude (for
    ## He-like carbon, example-Je.jl branch a, the n=2 levels sit near ~1e-16 at Te=1e5 K purely because of the
    ## exp(-dE/kT) Boltzmann factor), and flooring on |population| would silently erase that real physics
    ## instead of merely removing noise. A negative population, by contrast, is unphysical whatever its size.
    ## Note the reformulation above is what actually fixes the problem -- a decoupled level now comes out of
    ## the solve as an exact zero, so this loop is a safety net, not the mechanism. It does still do one
    ## necessary job: that exact zero is often a NEGATIVE zero (-0.0), which compares equal to 0. but which
    ## @sprintf("%.4e", ...) renders as "-0.0000e+00" in the results table, i.e. as an apparent negative
    ## population. Assigning 0. normalizes the sign of zero away.
    popFloor = 1.0e-13 * maximum(abs, populations)
    for  i = 1:N
        if  populations[i] < 0.  &&  abs(populations[i]) < popFloor    populations[i] = 0.    end
        if  populations[i] == 0.                                       populations[i] = 0.    end
    end
    if  any(populations .< 0.)
        @warn("Plasma.solveCollisionalRadiativeBalance: the balance returned NEGATIVE level populations well " *
              "above the roundoff floor (most negative: $(minimum(populations))). This is unphysical and " *
              "indicates an inconsistent rate matrix, not a numerical artefact -- treat these populations as " *
              "invalid.")
    end
    ##
    ## Renormalize: zeroing the roundoff components perturbs the sum only at the popFloor level, but the
    ## returned vector is documented to sum to 1, so restore that exactly.
    return( populations / sum(populations) )
end


"""
`Plasma.perform(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation; output::Bool=true)`
    ... to compute the (relative) collisional-radiative population balance among the levels generated for
        one ion (computation.refConfigs). For output=true, a dictionary is returned from which the relevant
        results can be accessed by proper keys.
"""
function  perform(scheme::Plasma.CollisionalRadiativeScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end

    nm  = computation.nuclearModel;   grid = computation.grid

    if  !computation.settings.useNumberDensity
        error("Plasma.CollisionalRadiativeScheme: computation.settings.useNumberDensity must be true, with " *
              "settings.density interpreted as the ELECTRON density Ne [cm^-3] -- the mass-density convention " *
              "used by other Plasma schemes is not supported here, since converting it to Ne requires " *
              "ionization-state information this scheme does not have.")
    end
    Ne = computation.settings.density

    # -------------------------------------------------------------------------------------------------------
    # (1) Generate the level set -- or reuse a cached one from scheme.levelsFilenames if it matches exactly.
    # -------------------------------------------------------------------------------------------------------
    levelFingerprint = Plasma.collisionalRadiativeLevelFingerprint(scheme, computation)
    repMultiplet      = Plasma.readEvaluateCollisionalRadiativeLevels(scheme, computation)
    if  ismissing(repMultiplet)
        repMultiplet = Plasma.generateCollisionalRadiativeLevels(scheme, computation)
        Plasma.writeCollisionalRadiativeLevels(scheme, computation, repMultiplet)
    end
    N = length(repMultiplet.levels)

    @warn("Plasma.CollisionalRadiativeScheme: autoionizing decay of any generated level that might lie above " *
          "the next charge state's ionization threshold is NOT computed by this scheme; all generated levels " *
          "are treated as purely bound (radiative + collisional only). This is a real limitation for a level " *
          "set that reaches that far in energy -- keep scheme.NoExcitations/upperShellNo modest enough that it " *
          "does not, or treat any near-threshold result with caution.")

    # -------------------------------------------------------------------------------------------------------
    # (2) Radiative rate matrix (PhotoEmission.determineLines already restricts this to downward, omega>0,
    #     transitions on its own -- no explicit lineSelection needed). gauges is forced to both Coulomb and
    #     Babushkin regardless of what scheme.peSettings carries: PhotoEmission.Settings()'s own default is
    #     gauges=[UseCoulomb] only (single-gauge), but this driver always reports both gauges side by side
    #     (see step 4), so Babushkin must actually be computed, not silently left at its zero default.
    # -------------------------------------------------------------------------------------------------------
    peSettings = PhotoEmission.Settings(scheme.peSettings; printBefore=false, gauges=[UseCoulomb, UseBabushkin])
    peLines    = PhotoEmission.computeLines(repMultiplet, repMultiplet, grid, peSettings; output=true)
    rateFactor = Defaults.convertUnits("rate: from atomic to 1/s", 1.0)

    # -------------------------------------------------------------------------------------------------------
    # (3) Collisional rate matrix: ImpactExcitation computed for upward pairs only (index i < f, since
    #     repMultiplet is sorted ascending in energy); de-excitation obtained by detailed balance from the
    #     SAME (symmetric) effective collision strength Omega(Te), not by a second, separately-computed
    #     ImpactExcitation call -- this also avoids the exp(+dE/kT) overflow/NaN risk a naive alpha-based
    #     detailed-balance formula would have far off threshold.
    # -------------------------------------------------------------------------------------------------------
    upwardPairs = [ (i, f)  for i in 1:N  for f in i+1:N ]
    ieLineSelection = LineSelection(true; indexPairs=upwardPairs)
    ieSettings = ImpactExcitation.Settings(scheme.ieSettings; calcRateCoefficient=true, lineSelection=ieLineSelection,
                                           printBefore=false)   # false always -- module-ImpactExcitation.jl displayLines() bug
    temperatures = ieSettings.temperatures

    # Reuse cached raw (pre-thermal-average) lines from scheme.ratesFilenames if they match exactly; the
    # thermal average into rate coefficients is always redone below, cheaply, for the CURRENT temperatures
    # list -- this is what makes any new scheme.ieSettings.temperatures free once the rates cache is valid.
    ieLines = Plasma.readEvaluateCollisionalRadiativeRates(scheme, computation, levelFingerprint)
    if  ismissing(ieLines)
        ieLines, _ = ImpactExcitation.computeLines(repMultiplet, repMultiplet, nm, grid, ieSettings; output=true)
        Plasma.writeCollisionalRadiativeRates(scheme, computation, levelFingerprint, ieLines)
    end
    ieRates = ImpactExcitation.computeRateCoefficients(ieLines, ieSettings)

    if output    results = Base.merge( results, Dict("levels:" => repMultiplet, "radiative lines:" => peLines,
                                                       "excitation lines:" => ieLines, "excitation rates:" => ieRates) )   end

    # -------------------------------------------------------------------------------------------------------
    # (4) Assemble and solve the balance, once per temperature and once per gauge. The population balance is
    #     level-resolved (a vector of N fractions), not a single scalar per temperature -- both gauges' full
    #     vectors are kept side by side (as a NamedTuple) rather than forced through EmProperty, which is
    #     built for scalar rates/cross sections, not vectors.
    # -------------------------------------------------------------------------------------------------------
    popVectors = Dict{Float64, NamedTuple{(:Coulomb,:Babushkin), Tuple{Array{Float64,1},Array{Float64,1}}}}()
    for  temp in temperatures
        tIndex = findfirst(==(temp), temperatures)
        cMatrix = zeros(N, N);   bMatrix = zeros(N, N)

        for  line in peLines
            i = line.initialLevel.index;   f = line.finalLevel.index
            cRate = rateFactor * line.photonRate.Coulomb;    bRate = rateFactor * line.photonRate.Babushkin
            if  isfinite(cRate)   cMatrix[i, f] += cRate;   cMatrix[i, i] -= cRate   end
            if  isfinite(bRate)   bMatrix[i, f] += bRate;   bMatrix[i, i] -= bRate   end
        end
        for  rc in ieRates
            i = rc.initialLevel.index;   f = rc.finalLevel.index
            excRate   = rc.alphas[tIndex] * Ne
            gUpper    = Float64(Basics.twice(rc.finalLevel.J) + 1)
            deExcRate = 8.6291269e-6 * rc.effOmegas[tIndex] / (gUpper * sqrt(temp)) * Ne
            if  isfinite(excRate)
                cMatrix[i, f] += excRate;   cMatrix[i, i] -= excRate
                bMatrix[i, f] += excRate;   bMatrix[i, i] -= excRate
            end
            if  isfinite(deExcRate)
                cMatrix[f, i] += deExcRate;   cMatrix[f, f] -= deExcRate
                bMatrix[f, i] += deExcRate;   bMatrix[f, f] -= deExcRate
            end
        end

        popVectors[temp] = (Coulomb = Plasma.solveCollisionalRadiativeBalance(cMatrix),
                            Babushkin = Plasma.solveCollisionalRadiativeBalance(bMatrix))
    end

    Plasma.displayCollisionalRadiativeResults(stdout, repMultiplet, temperatures, Ne, popVectors)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   Plasma.displayCollisionalRadiativeResults(iostream, repMultiplet, temperatures, Ne, popVectors)   end

    if output    results = Base.merge( results, Dict("populations:" => popVectors) )   end

    println("CollisionalRadiativeScheme computation complete ...")
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end


"""
`Plasma.displayCollisionalRadiativeResults(stream::IO, repMultiplet, temperatures, Ne, popVectors)`
    ... prints a table of the relative level populations (Coulomb/Babushkin), one row per (level,
        temperature) pair, following JAC's usual table style (2-line header, truncated scientific
        notation, energies/temperatures in the usual units).
"""
function  displayCollisionalRadiativeResults(stream::IO, repMultiplet::Multiplet, temperatures::Array{Float64,1}, Ne::Float64,
                                             popVectors::Dict)
    nx = 100
    println(stream, "\n  Collisional-radiative relative level populations (Ne = $Ne cm^-3):")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "    " * TableStrings.flushleft(9, "Level") * TableStrings.flushleft(9, "J^P") *
                  TableStrings.flushright(16, "Energy") * TableStrings.flushright(14, "Temperature") *
                  TableStrings.flushright(22, "Population") * TableStrings.flushright(22, "Population")
    sb = "    " * TableStrings.flushleft(9, "")      * TableStrings.flushleft(9, "") *
                  TableStrings.flushright(16, TableStrings.inUnits("energy")) * TableStrings.flushright(14, "[K]") *
                  TableStrings.flushright(22, "(Coulomb)") * TableStrings.flushright(22, "(Babushkin)")
    println(stream, sa);   println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    for  level in repMultiplet.levels
        energy = Defaults.convertUnits("energy: from atomic", level.energy)
        jp     = "$(level.J) $(string(level.parity))"
        for  temp in temperatures
            pv = popVectors[temp]
            sc = "    " * TableStrings.flushleft(9, string(level.index)) * TableStrings.flushleft(9, jp) *
                          TableStrings.flushright(16, @sprintf("%.4e", energy)) *
                          TableStrings.flushright(14, @sprintf("%.4e", temp)) *
                          TableStrings.flushright(22, @sprintf("%.4e", pv.Coulomb[level.index])) *
                          TableStrings.flushright(22, @sprintf("%.4e", pv.Babushkin[level.index]))
            println(stream, sc)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end
