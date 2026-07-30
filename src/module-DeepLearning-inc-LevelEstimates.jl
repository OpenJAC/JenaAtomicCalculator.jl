
# Date structures, functions and methods for (level) estimates of levels that are missing in the NIST database


"""
`struct  DeepLearning.LsjLevel`
    ... defines a type for collecting information about a single LSJ with two different energies, a calculated or
        predicted energy and the associated NIST, if this assignment is possible.

    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity J of the LSJ-level
    + thisEnergy        ::Float64             ... excitation energy, either compute or predicted.
    + nistEnergy        ::Float64             ... excitation energy, either compute or predicted.
"""
struct  LsjLevel
    conf                ::Configuration
    L                   ::AngularJ64
    S                   ::AngularJ64
    J                   ::AngularJ64
    parity              ::Basics.Parity
    thisEnergy          ::Float64
    nistEnergy          ::Float64
end


# `Base.show(io::IO, level::DeepLearning.LsjLevel)`  ... prepares a proper printout of the variable level::DeepLearning.LsjLevel
function Base.show(io::IO, level::DeepLearning.LsjLevel)
    println(io, "conf:        $(level.conf)  ")
    println(io, "L:           $(level.L)  ")
    println(io, "S:           $(level.S)  ")
    println(io, "J:           $(level.J)  ")
    println(io, "parity:      $(level.parity)  ")
    println(io, "thisEnergy:  $(level.thisEnergy)  ")
    println(io, "nistEnergy:  $(level.nistEnergy)  ")
end


"""
`struct  DeepLearning.NistLevel`
    ... defines a type for collecting information about a single level from the NIST Level database.

    + conf              ::Configuration       ... Configuration of the LSJ-level
    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity of the LSJ-level
    + energy            ::Float64             ... total energy, compute with the given atomic model.
"""
struct  NistLevel
    conf                ::Configuration
    L                   ::AngularJ64
    S                   ::AngularJ64
    J                   ::AngularJ64
    parity              ::Basics.Parity
    energy              ::Float64
end


# `Base.show(io::IO, level::DeepLearning.NistLevel)`  ... prepares a proper printout of the variable level::DeepLearning.NistLevel
function Base.show(io::IO, level::DeepLearning.NistLevel)
    println(io, "conf:        $(level.conf)  ")
    println(io, "L:           $(level.L)  ")
    println(io, "S:           $(level.S)  ")
    println(io, "J:           $(level.J)  ")
    println(io, "parity:      $(level.parity)  ")
    println(io, "energy:      $(level.energy)  ")
end


"""
`struct  DeepLearning.LevelEstimationRequest  <:  DeepLearning.AbstractNeuralNetworkRequest`
    ... to define a (deep-learning) request for estimating the level energies for a given set of configurations.

    + configs       ::Array{Configuration,1}      ... List of configuration for which levels are to be estimated.
    + nistLevels    ::Array{NistLevel,1}          ... Available list of NIST levels used for comparison.
"""
struct   LevelEstimationRequest  <:  DeepLearning.AbstractNeuralNetworkRequest
    configs         ::Array{Configuration,1}
    nistLevels      ::Array{NistLevel,1}
end


"""
`DeepLearning.LevelEstimationRequest()`  ... constructor for an 'default' instance of a DeepLearning.LevelEstimationRequest.
"""
function LevelEstimationRequest()
    LevelEstimationRequest( Configuration[], NistLevel[])
end


# `Base.string(request::LevelEstimationRequest)`  ... provides a String notation for the variable request::LevelEstimationRequest.
function Base.string(request::LevelEstimationRequest)
    sa = "NN request for level estimation:"
    return( sa )
end


# `Base.show(io::IO, request::LevelEstimationRequest)`  ... prepares a proper printout of the request::LevelEstimationRequest.
function Base.show(io::IO, request::LevelEstimationRequest)
    sa = Base.string(request);        print(io, sa, "\n")
    println(io, "configs:             $(request.configs)  ")
    nShow = min(3, length(request.nistLevels))
    println(io, "nistLevels[1:$nShow]:     $(request.nistLevels[1:nShow])  ")
end



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
`DeepLearning.checkConfiguration(conf::Configuration, am::AtomicFeatures.AtomicModel)`
    ... checks that the given configuration can be encoded into/used by the given atomic model;
        a boolian value of true/false is returned. If false is returned, a warning is printed about
        the reasons, usually concerning the maximum n quantum numbers involved.
"""
function checkConfiguration(conf::Configuration, am::AtomicFeatures.AtomicModel)
    shells = Basics.extractShellList([conf]);   val = true

    # No shell should be larger than nMax and lMax
    for  shell in shells
        if  shell.n > am.nMax  ||  shell.l > am.lMax   val = false
            @warn "$conf has shells > $(Shell(am.nMax, am.lMax)) "
            break
        end
    end

    return( val )
end


"""
`DeepLearning.computeMultiplet(conf::Configuration, am::AtomicFeatures.AtomicModel)`
    ... computes the level multiplet as associated with the given configuration, (frozen) set of orbitals and asfSettings;
        these settings just control the CI computations for the levels of a single configuration.
"""
function computeMultiplet(conf::Configuration, am::AtomicFeatures.AtomicModel)
    orbitals  = DeepLearning.extractOrbitals(conf, am)
    multiplet = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, am.nuclearModel, am.grid, am.asfSettings;
                                                        printout=true)
    return( multiplet )
end


"""
`DeepLearning.computeLSjjMultiplet(multiplet::Multiplet)`
    ... computes the jj-to-LS expansion of all levels in `multiplet`, composed from LSjj's own (exported)
        building blocks -- LSjj.generateNonrelativisticCsfList, LSjj.BasisNR, LSjj.expandCsfRintoNonrelativisticBasis,
        LSjj.LevelNR, LSjj.MultipletNR -- the same pieces that LSjj.expandLevelsIntoLS uses internally; that
        function itself only prints its result and always returns nothing, so it cannot be reused directly here.
        A multipletNR::LSjj.MultipletNR is returned, or nothing if multiplet's basis does not have a standard
        subshell list.
"""
function computeLSjjMultiplet(multiplet::Multiplet)
    if  !Basics.isStandardSubshellList(multiplet.levels[1].basis)
        @warn "DeepLearning.computeLSjjMultiplet(): non-standard subshell list, no jj-LS expansion possible."
        return( nothing )
    end

    shellList = Basics.extractNonrelativisticShellList(multiplet.levels[1].basis.subshells)
    confList  = Basics.extractConfigurations(Basics.FromBasis(), multiplet.levels[1].basis)
    csfsNR    = LSjj.generateNonrelativisticCsfList(confList, shellList)
    basisNR   = LSjj.BasisNR(multiplet.levels[1].basis.NoElectrons, shellList, csfsNR)
    ncsfs     = length(basisNR.csfs)

    mcVectors = Dict{Int64, Array{Float64,1}}()
    for  levelR in multiplet.levels    mcVectors[levelR.index] = zeros(ncsfs)    end

    for  r = 1:length(multiplet.levels[1].basis.csfs)
        csfR       = multiplet.levels[1].basis.csfs[r]
        conf       = Basics.extractConfiguration(Basics.FromBasis(), multiplet.levels[1].basis, csfR)
        openShells = Basics.extractFromConfiguration(Basics.OpenShellNumber(), conf)
        mcCsfR     = LSjj.expandCsfRintoNonrelativisticBasis(openShells, csfR, multiplet.levels[1].basis, basisNR)
        for  levelR in multiplet.levels
            mcVectors[levelR.index] = mcVectors[levelR.index] + levelR.mc[r] * mcCsfR
        end
    end

    levelsNR = LSjj.LevelNR[]
    for  levelR in multiplet.levels
        push!( levelsNR, LSjj.LevelNR(levelR.J, levelR.parity, levelR.index, levelR.energy, basisNR, mcVectors[levelR.index]) )
    end

    return( LSjj.MultipletNR("LS-expanded " * multiplet.name, levelsNR) )
end


"""
`DeepLearning.extractLeadingLSTerm(levelNR::LSjj.LevelNR)`
    ... determines the (L,S) term of the leading CsfNR (largest |mixing coefficient|^2) of levelNR, for the
        closed-shell / one-open-shell / several-open-shells cases relevant to the low-lying Ar-ion levels
        considered here: zero open shells -> L=S=0; one open shell -> that shell's own (shellL,shellS); two
        or more open shells -> the LAST open shell's cumulative (shellLX,shellSX) -- the standard JAC
        convention for progressively-coupled antisymmetric shell states (mirroring how subshellX in the
        relativistic CsfR already represents a cumulative intermediate coupling). levelNR.J/.parity already
        give the total J/parity directly and are unaffected by this LS decomposition. A tuple
        (L::AngularJ64, S::AngularJ64, weight::Float64) is returned.
"""
function extractLeadingLSTerm(levelNR::LSjj.LevelNR)
    idx    = sortperm(abs.(levelNR.mc).^2, rev=true)
    r      = idx[1];   weight = abs(levelNR.mc[r])^2
    csfNR  = levelNR.basis.csfs[r]
    openShellQN = LSjj.extractOpenShellQNfromCsfNR(csfNR, levelNR.basis)

    if      length(openShellQN) == 0
        L = AngularJ64(0);   S = AngularJ64(0)
    elseif  length(openShellQN) == 1
        qn = first(values(openShellQN))
        L  = AngularJ64(qn[5], 2);   S = AngularJ64(qn[6], 2)
    else
        lastQn = nothing
        for  shell in levelNR.basis.shells
            if  haskey(openShellQN, shell)   lastQn = openShellQN[shell]   end
        end
        L = AngularJ64(lastQn[7], 2);   S = AngularJ64(lastQn[8], 2)
    end

    return( (L, S, weight) )
end


"""
`DeepLearning.extractLsjLevels(conf::Configuration, multiplet::Multiplet, nistLevels::Array{NistLevel,1})`
    ... extracts all LsjLevels that are associated with the given configuration. It also assigns the correct
        excitation energy from the NIST tables, if available. The procedure performs three steps:
        (i)    the jj-LS expansion of all levels in the multiplet (DeepLearning.computeLSjjMultiplet);
        (ii)   the extraction of each level's leading (L,S) term (DeepLearning.extractLeadingLSTerm);
        (iii)  the assignment of the matching NIST energy, if a NIST level with the same configuration
               and (L,S,J,parity) exists.
        A list levels::Array{LsjLevel,1} is returned, one entry per level in multiplet (in the same
        order), with nistEnergy set to 0. where no matching NIST level was found.
"""
function extractLsjLevels(conf::Configuration, multiplet::Multiplet, nistLevels::Array{NistLevel,1})
    lsjLevels = LsjLevel[]

    # Extract the NIST levels that belong to the given configuration
    relevantLevels = NistLevel[]
    for level in nistLevels
        if  conf == level.conf   push!(relevantLevels, level)   end
    end

    # Perform the jj-LS expansion of all levels in the multiplet
    multipletNR = DeepLearning.computeLSjjMultiplet(multiplet)

    for  (i, level)  in  enumerate(multiplet.levels)
        if  multipletNR === nothing
            L = AngularJ64(0);   S = AngularJ64(0)
        else
            L, S, _ = DeepLearning.extractLeadingLSTerm(multipletNR.levels[i])
        end

        nistEnergy = 0.
        for  rLevel in relevantLevels
            # Compare via Basics.twice(...) rather than raw AngularJ64 == : AngularJ64 is not auto-reduced
            # (e.g. J=1 could be stored as (num=1,den=1) or (num=2,den=2) depending on where it came from),
            # so a direct struct == could spuriously fail for integer L/S/J values.
            if  Basics.twice(rLevel.L) == Basics.twice(L)  &&  Basics.twice(rLevel.S) == Basics.twice(S)  &&
                Basics.twice(rLevel.J) == Basics.twice(level.J)  &&  rLevel.parity == level.parity
                nistEnergy = rLevel.energy;   break
            end
        end

        push!(lsjLevels, LsjLevel(conf, L, S, level.J, level.parity, level.energy, nistEnergy))
    end

    return( lsjLevels )
end


"""
`DeepLearning.displayLsjLevels(conf::Configuration, lsjLevels::Array{LsjLevel,1})`
    ... displays (and lists) all given LsjLevels in a neat format; nothing is returned.
"""
function displayLsjLevels(conf::Configuration, lsjLevels::Array{LsjLevel,1})
    nx = 105
    println(" ")
    println("  Comparison of calculated/predicted (excitation) energies with NIST energies for $conf:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(14, "level  "       ; na=0);          sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(18, "^(2S+1) L_J^P" ; na=2);          sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "This energy"   ; na=3);
    sb = sb * TableStrings.center(10, "eV"            ; na=3)
    sa = sa * TableStrings.center(10, "NIST energy"   ; na=3)
    sb = sb * TableStrings.center(10, "eV"            ; na=2)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    #
    for  (i, level) in enumerate(lsjLevels)
        sym = LevelSymmetry( level.J, level.parity)
        sm  = round(Int64, Basics.twice(level.S)+1)
        sL  = round(Int64, Basics.twice(level.L)/2.)
        si  = "    " * string(i)
        sy  = string(sym) * "         "
        sa  = "   " * si[end-4:end] * "   " * "^" * string(sm) * " " * string(sL) * "_" * sy[1:8]
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", level.thisEnergy))    * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", level.nistEnergy))    * "   "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end


"""
`DeepLearning.extractNistConfigurations(nistLevels::Array{NistLevel,1})`
    ... extracts all configurations that are defined by the given list of NIST levels.
"""
function extractNistConfigurations(nistLevels::Array{NistLevel,1})
    configs = Configuration[]

    for level  in  nistLevels   push!(configs, level.conf)  end
    configs = unique(configs)

    return( configs )
end


"""
`DeepLearning.parseNistConfiguration(sa::AbstractString, coreConfig::String="")`
    ... converts a dot-separated NIST-style configuration string, such as "3s2.3p5", into a proper
        JAC Configuration. Each dot-separated token is a shell followed by its occupation (e.g. "3p5"
        means shell 3p with occupation 5); a token without a trailing digit (e.g. "3s") means
        occupation 1. Some NIST strings additionally include a parenthesized PARENT-TERM annotation
        between shells (e.g. "3s2.3p3.(4S).3d", the "(4S)" indicating the ⁴S parent term of the 3p^3
        core coupled to the added 3d electron) -- such tokens are skipped, since a plain Configuration
        carries occupations only, not parent-coupling terms; this is a deliberate Stage-1
        simplification (occasionally two parent-coupling variants of the same configuration could
        share the same overall (L,S,J), which this parser cannot distinguish). NIST configuration
        strings only ever list the VALENCE shells, never the closed inner core -- if `coreConfig` is
        given (e.g. "[Ne]"), it is prepended so the returned Configuration matches the FULL
        configuration used in an actual atomic-structure computation (needed so that NIST levels can
        later be matched against computed levels by Configuration equality). A conf::Configuration
        is returned.
"""
function parseNistConfiguration(sa::AbstractString, coreConfig::String="")
    tokens = split(strip(sa), '.');   parts = String[]

    for  tok in tokens
        if  startswith(tok, '(')  &&  endswith(tok, ')')   continue    end
        m = match(r"^(\d+)([a-zA-Z])(\d*)$", tok)
        if  m === nothing   error("DeepLearning.parseNistConfiguration(): cannot parse shell token '$tok' in '$sa'.")   end
        n = m.captures[1];   l = lowercase(m.captures[2]);   occStr = m.captures[3]
        occ = occStr == "" ? 1 : parse(Int64, occStr)
        push!(parts, "$n$l^$occ")
    end

    valencePart = join(parts, " ")
    fullPart    = coreConfig == "" ? valencePart : coreConfig * " " * valencePart

    return( Configuration(fullPart) )
end


"""
`DeepLearning.extractNistLevels(filenames::Array{String,1}, coreConfig::String="")`
    ... extracts all available NIST levels from a (number of) ASCII-files, in which these data are encoded.
        Each file is expected to have one header line, followed by whitespace-separated columns
        `Configuration  L  S  J  Level(eV)` (the exact header text is not checked, the first line of
        every file is simply skipped). `coreConfig` (e.g. "[Ne]" for the Ar ions considered here) is
        prepended to every parsed configuration (see DeepLearning.parseNistConfiguration) so that the
        returned NistLevel.conf matches the FULL configuration of an actual computed level, not just
        its valence part -- required for DeepLearning.extractLsjLevels's Configuration-equality match
        to ever succeed. Parity is not given in the files and is derived from the parsed (full)
        configuration via Basics.extractFromConfiguration(Basics.GetParity(), conf). Energies are
        converted from eV to atomic units. A list nistLevels::Array{NistLevel,1} is returned.
"""
function extractNistLevels(filenames::Array{String,1}, coreConfig::String="")
    nistLevels = NistLevel[]

    for filename  in  filenames
        rawLines = readlines(filename)
        for  (i, rawLine) in enumerate(rawLines)
            if  i == 1   continue    end
            line = strip(rawLine)
            if  line == ""   continue    end
            tokens = split(line)
            if  length(tokens) < 5   continue    end

            conf   = DeepLearning.parseNistConfiguration(tokens[1], coreConfig)
            L      = AngularJ64( round(Int64, 2*parse(Float64, tokens[2])), 2 )
            S      = AngularJ64( round(Int64, 2*parse(Float64, tokens[3])), 2 )
            J      = AngularJ64( round(Int64, 2*parse(Float64, tokens[4])), 2 )
            parity = Basics.extractFromConfiguration(Basics.GetParity(), conf)
            energy = Defaults.convertUnits("energy: to atomic", parse(Float64, tokens[5]))

            push!(nistLevels, NistLevel(conf, L, S, J, parity, energy))
        end
    end

    return( nistLevels )
end


"""
`DeepLearning.extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)`
    ... extracts for the given configuration the orbitals from the atomic model.
"""
function extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)
    q = round(Int64, am.nuclearModel.Z - conf.NoElectrons)

    if      q == 1   orbitals = am.orbitals01
    elseif  q == 2   orbitals = am.orbitals02
    elseif  q == 3   orbitals = am.orbitals03
    else    error("stop a")
    end

    return( orbitals )
end


"""
`DeepLearning.generateAtomicModelForLE_Arn4()`
    ... generates an atomic model for the training, test and feature extraction for the (level) estimation of atomic levels.
        Here, in particular, we consider the levels of Ar^+, Ar^2+ and Ar^3+ and configurations with maximum nMax = 4 shells.
        All details are hard-coded and follow some 'script-like' style. We expect that a particular generateAtomicModel...()
        is designed for each neutral network, which we shall train and consider.

        An atomicModel::AtomicFeatures.AtomicModel is returned.
"""
function generateAtomicModelForLE_Arn4()
    nMax = 4;   lMax = 3
    nm          = Nuclear.Model(18.)
    grid        = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    asfSettings = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))

    # Generate all shells and subshells of the model
    shells    = Basics.generateShellList(1, nMax, lMax)
    subshells = Basics.generateSubshellList(shells)

    # Generate mean-field orbitals separately for each charge state q=1,2,3 (Ar+, Ar2+, Ar3+), each from
    # its own reference configuration with 1-electron excitations across the full shell list, so that
    # every subshell in `subshells` has an available mean-field orbital for that charge state.
    mfSettings = AtomicState.MeanFieldSettings(Basics.DFSField(1.0))

    refConfigs01 = Basics.generateConfigurations(Basics.ExciteElectrons(1, shells, shells), [Configuration("[Ne] 3s^2 3p^5")])
    meanField01  = AtomicState.Representation("Ar+ mean field",  nm, grid, refConfigs01, AtomicState.MeanFieldBasis(mfSettings))
    mfrep01      = Basics.generate(meanField01; output=true)
    orbitals01   = mfrep01["mean-field basis"].orbitals

    refConfigs02 = Basics.generateConfigurations(Basics.ExciteElectrons(1, shells, shells), [Configuration("[Ne] 3s^2 3p^4")])
    meanField02  = AtomicState.Representation("Ar2+ mean field", nm, grid, refConfigs02, AtomicState.MeanFieldBasis(mfSettings))
    mfrep02      = Basics.generate(meanField02; output=true)
    orbitals02   = mfrep02["mean-field basis"].orbitals

    refConfigs03 = Basics.generateConfigurations(Basics.ExciteElectrons(1, shells, shells), [Configuration("[Ne] 3s^2 3p^3")])
    meanField03  = AtomicState.Representation("Ar3+ mean field", nm, grid, refConfigs03, AtomicState.MeanFieldBasis(mfSettings))
    mfrep03      = Basics.generate(meanField03; output=true)
    orbitals03   = mfrep03["mean-field basis"].orbitals

    atomicModel = AtomicFeatures.AtomicModel(nMax, lMax, grid, asfSettings, nm, shells, subshells,
                                             orbitals01, orbitals02, orbitals03)

    return( atomicModel )
end


"""
`DeepLearning.generateFeatureVectors(confs::Array{Configuration,1}, am::AtomicFeatures.AtomicModel, probability::Float64,
                                     nistLevels::Array{NistLevel,1})`
    ... generates a list of -- training and test -- feature vectors which are associated with the given configurations.
        These feature vectors can be written out and used for the training and testing of NN as needed by the
        applied DeepLearning toolboxes (outside of JAC). The procedure assumes that the atomic model am has
        been fully determined; it simply loops through the (NIST) configurations and extracts all useful
        feature x-vectors (which also contain the NIST energy as y-vector). Two level lists
        (trainXyVectors::Array{AtomicFeatures.XyVector,1}, testXyVectors::Array{AtomicFeatures.XyVector,1})
        are returned.
"""
function generateFeatureVectors(confs::Array{Configuration,1}, am::AtomicFeatures.AtomicModel, probability::Float64,
                                nistLevels::Array{NistLevel,1})

    trainXyVectors = AtomicFeatures.XyVector[];   testXyVectors = AtomicFeatures.XyVector[]

    for conf in confs
        if  !DeepLearning.checkConfiguration(conf, am)   continue    end

        # Compute the level multiplet in the given atomic model and assign the LSJ-levels
        multiplet = DeepLearning.computeMultiplet(conf, am)
        lsjLevels = DeepLearning.extractLsjLevels(conf, multiplet, nistLevels)
        DeepLearning.displayLsjLevels(conf, lsjLevels)

        # Generate and select the feature x- (and y-) vectors
        xyVectors = DeepLearning.generateFeatureVectors(conf, am, multiplet, lsjLevels)
        wxy       = DeepLearning.selectFeatureVectors(xyVectors, probability)
        append!(trainXyVectors, wxy[1]);    append!(testXyVectors, wxy[2])
    end

    return( (trainXyVectors, testXyVectors) )
end


"""
`DeepLearning.generateFeatureVectors(conf::Configuration, am::AtomicFeatures.AtomicModel, multiplet::Multiplet,
                                     lsjLevels::Array{LsjLevel,1})`
    ... generates a list of feature vectors which are associated with the given configuration.
        This generation is based on 'script-like' code to enable the user to play with different assumptions
        and settings about the relevance of different features. Usually, the length of the feature vectors
        increases very rapidly with the number of "physical features" that are included into the training
        and extraction of data.

        A list of xyVectors::AtomicFeatures.XyVector is returned.
"""
function generateFeatureVectors(conf::Configuration, am::AtomicFeatures.AtomicModel, multiplet::Multiplet, lsjLevels::Array{LsjLevel,1})
    xyVectors = AtomicFeatures.XyVector[]

    orbitals  = DeepLearning.extractOrbitals(conf, am)

    # Append for each level in multiplet all desired features to the xVector
    for  (i, level)  in  enumerate(multiplet.levels)
        xVector = Float64[]

        # Add shell occupations of the configuration
        append!(xVector,  AtomicFeatures.extractShellOccupations(am.shells, conf) )

        # Add mean occupation of the given LSJ level
        append!(xVector,  AtomicFeatures.extractMeanOccupationNumbers(am.subshells, level) )

        # Compute and add LSJ quantum numbers of level
        L = lsjLevels[i].L;    S = lsjLevels[i].S;    J = lsjLevels[i].J
        if  J != level.J  error("stop a")   end
        append!(xVector,  [Basics.twice(L)/2.0, Basics.twice(S)/2.0, Basics.twice(J)/2.0] )

        # Compute and add intermediate coupling quantum numbers of the n-leading CSF of the levels
        # (with zeros, if < n CSF are defined for the given level
        append!(xVector,  AtomicFeatures.extractIntermediateQN(2, level) )

        # Add <r^k> of the given subshells
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells,  1, orbitals, am.grid) )
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells,  2, orbitals) )

        # Add F^k (a,b) = R^k(a,b,a,b) of the given subshells
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  0, orbitals, am.grid) )
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  1, orbitals) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  2, orbitals, am.grid) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  4, orbitals, am.grid) )

        # Add G^k (a,b) = R^k(a,b,b,a) of the given subshells
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  0, orbitals, am.grid) )
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  1, orbitals) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  2, orbitals, am.grid) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  4, orbitals, am.grid) )


        # Append the newxVector to the desired list
        push!(xyVectors, AtomicFeatures.XyVector(L, S, J, level.parity, level.energy, lsjLevels[i].nistEnergy, xVector))
    end

    println(">> $(length(xyVectors)) feature vectors have been generated for configuration $conf ")

    return( xyVectors )
end


"""
`DeepLearning.run(request::DeepLearning.LevelEstimationRequest, applic::DeepLearning.Application; output::Bool=true)`
    ... to run a deep-learning level estimation (request) to estimate level energies for given configurations,
        and based on a given atomic model. The predicted energies are compared explicitly with the data availalbe
        from the NIST Level database. A dict::Dict{String, Any} is returned if output=true, and nothing otherwise.

        NOTE (Stage 1 status): this function needs a *trained* neural network (module-DeepLearningNetworks.jl,
        not yet implemented) to actually predict energies; it is not part of Stage 1 (feature-vector generation)
        and is left as a placeholder for Stage 2.
"""
function  run(request::DeepLearning.LevelEstimationRequest, applic::DeepLearning.Application; output::Bool=true)
    error("DeepLearning.run(): not yet implemented -- needs a trained network (Stage 2, module-DeepLearningNetworks.jl).")
end



"""
`DeepLearning.selectFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, probability::Float64)`
    ... divides the set of feature vectors into two sets: for training and for testing due to the given
        probability for selecting test data (probability <= 0.2). Two list of feature vectors
        (trainXyVectors::Array{AtomicFeatures.XyVector,1}, testXyVectors::Array{AtomicFeatures.XyVector,1})
        are returned.
"""
function selectFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, probability::Float64)
    trainXyVectors = AtomicFeatures.XyVector[];    testXyVectors = AtomicFeatures.XyVector[]
    wrn = rand( length(xyVectors) );   nw = 0

    for xyVector in xyVectors
        if   xyVector.nistEnergy == 0.   continue    end
        nw = nw + 1
        if   wrn[nw] < probability    push!(testXyVectors, xyVector)
        else                          push!(trainXyVectors, xyVector)
        end
    end

    return( (trainXyVectors, testXyVectors) )
end


"""
`DeepLearning.writeFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, filename::String)`
    ... writes out to a file(name) the xyVectors in a format suitable for the training and test of neural networks:
        one row per xyVector, its xVector entries space-separated, followed by its nistEnergy. Nothing is returned.
"""
function writeFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, filename::String)
    ioFeatures = open(filename, "w")

    for  xyVector in  xyVectors
        line = join( [@sprintf("%.8e", v) for v in xyVector.xVector], "  " )
        line = line * "        " * @sprintf("%.8e", xyVector.nistEnergy) * "\n"
        write(ioFeatures, line )
    end

    close(ioFeatures)

    return(nothing)
end
