
#####################################################################################################################
#####  Second-order treatment of the Q space  #######################################################################
#####################################################################################################################

"""
`struct  Hamiltonian.QPartition`
    ... records how the Q space of one step was divided, and whether the division was justified. Returned beside the
        multiplet by `Hamiltonian.performCIwithSecondOrderQ` so that a caller can act on the numbers rather than on
        the printed report.

    + promoted       ::Array{Configuration,1}          ... configurations put INTO the CI and diagonalized exactly.
    + folded         ::Array{Configuration,1}          ... configurations folded in through second order only.
    + discarded      ::Array{Configuration,1}          ... configurations dropped entirely.
    + weights        ::Array{Tuple{Configuration,Float64},1}  ... |c|^2 of every configuration of the Q space, the
                                                           MAXIMUM over the target levels (a configuration that matters
                                                           to any one level must enter the basis they all share),
                                                           sorted descending. A LIST and not a Dict because the order
                                                           is the point: a caller reads this to see WHICH configurations
                                                           carry the weight, and a Dict would throw that away.
    + deltaE         ::Array{Float64,1}                ... second-order contribution of the folded part, per level.
    + sumC2Folded    ::Array{Float64,1}                ... Sum |c|^2 over the folded part, per level. THE GATE.
    + isJustified    ::Bool                            ... whether that sum stayed below `promoteAbove` for every
                                                           level, i.e. whether the folded remainder really is a weak
                                                           perturbation. False means the energies are not to be used.
    + unsound        ::Array{Tuple{Configuration,Float64},1}  ... configurations carrying an INTRUDER level: one that
                                                           came out of the enlarged CI below the reference level while
                                                           holding almost none of the reference, with that eigenvalue.
                                                           Empty is the normal case. A non-empty list means the CI has
                                                           been handed a state that displaces the reference rather than
                                                           corrects it -- see the guard in `Hamiltonian.secondOrderBlock`.
"""
struct  QPartition
    promoted            ::Array{Configuration,1}
    folded              ::Array{Configuration,1}
    discarded           ::Array{Configuration,1}
    weights             ::Array{Tuple{Configuration,Float64},1}
    deltaE              ::Array{Float64,1}
    sumC2Folded         ::Array{Float64,1}
    isJustified         ::Bool
    unsound             ::Array{Tuple{Configuration,Float64},1}
end


# `Base.show(io::IO, part::Hamiltonian.QPartition)`  ... prepares a proper printout of part::Hamiltonian.QPartition.
function Base.show(io::IO, part::Hamiltonian.QPartition)
    println(io, "QPartition:  $(length(part.promoted)) promoted, $(length(part.folded)) folded in, " *
                "$(length(part.discarded)) discarded;  justified = $(part.isJustified)" *
                (isempty(part.unsound) ? "" : ";  $(length(part.unsound)) UNSOUND promotions"))
end


"""
`Hamiltonian.matrixElement(basis::Basis, r::Int64, s::Int64, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::AsfSettings, potential::Radial.Potential, cache::InteractionStrength.XLCache)`
    ... computes the single Hamiltonian matrix element between the CSF with indices r and s of the given basis, under
        the electron-electron interaction and QED model named by the settings. This is the body of the inner loop of
        `Hamiltonian.setupMatrix`, taken out so that a caller which needs only a BLOCK of the matrix -- the P rows
        against one configuration of Q, say -- can have it without allocating the whole symmetry block. A
        value::Float64 is returned.
"""
function matrixElement(basis::Basis, r::Int64, s::Int64, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::AsfSettings, potential::Radial.Potential, cache::InteractionStrength.XLCache)
    if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    return( 0. )    end
    subshellList = basis.subshells
    opa  = SpinAngular.OneParticleOperator(0, plus)
    waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)
    opa  = SpinAngular.TwoParticleOperator(0, plus)
    waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)
    #
    me = 0.
    for  coeff in waG1
        me = me + coeff.T * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
        if  settings.qedModel != NoneQed()
            me = me + InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nm,
                                                      settings.qedModel, potential, grid)
        end
    end
    for  coeff in waG2
        if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                        basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, cache)
        end
        if  typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                      basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                      settings.eeInteractionCI, cache)
        end
    end

    return( me )
end


"""
`Hamiltonian.groupCsfsByConfiguration(basis::Basis, indices::Array{Int64,1})`
    ... groups the given CSF indices of the basis by their NON-RELATIVISTIC configuration. Grouping is by the
        nl-configuration and not by the relativistic one on purpose: the jj partners of one nl shell then stay
        together, so a threshold can never keep `3p_^2 3p` while dropping `3p_ 3p^2` and leave the expansion
        unbalanced in the spin-orbit splitting.

        THE GROUPING IS BY A CANONICAL STRING AND NOT BY THE Configuration ITSELF. The original reason was a defect:
        `Configuration` satisfied `==` while `hash` fell back to object identity, so a `Dict` keyed by one put 426
        CSFs into 425 groups -- one per CSF -- with no error anywhere. THAT DEFECT IS FIXED (12-Sep-2026; see
        `Base.hash(conf::Basics.Configuration, h::UInt)` in `module-Basics-inc-first.jl`) and a Dict keyed by a
        Configuration is now sound. The String key is kept because it is also SORTABLE and READABLE, which is what
        the reports below need, not because it must be. A tuple
        `(groups, confOf)::Tuple{Dict{String,Array{Int64,1}}, Dict{String,Configuration}}` is returned, the second
        map holding one representative Configuration per key for reporting.
"""
function groupCsfsByConfiguration(basis::Basis, indices::Array{Int64,1})
    groups = Dict{String,Array{Int64,1}}();    confOf = Dict{String,Configuration}()
    for  idx in indices
        conf = Basics.extractConfiguration(Basics.FromBasis(), basis, basis.csfs[idx])
        key  = Hamiltonian.configurationKey(conf)
        if  haskey(groups, key)   push!(groups[key], idx)
        else                      groups[key] = [idx];    confOf[key] = conf
        end
    end

    return( (groups, confOf) )
end


"""
`Hamiltonian.configurationKey(conf::Configuration)`
    ... builds a canonical String for the given (non-relativistic) configuration -- the occupied shells in a fixed
        order -- so that configurations which compare equal give the same key, and so that the key can be sorted and
        printed. A key::String is returned.
"""
function configurationKey(conf::Configuration)
    return( join(sort([ string(sh) * "^" * string(w)  for (sh,w) in conf.shells  if w > 0 ]), " ") )
end


"""
`Hamiltonian.rankQspace(basis::Basis, pIndices::Array{Int64,1}, qGroups::Dict{String,Array{Int64,1}},
                         energies::Array{Float64,1}, vectors::Array{Float64,2}, nm::Nuclear.Model, grid::Radial.Grid,
                         settings::AsfSettings, cache::InteractionStrength.XLCache,
                         targetLevels::Array{Int64,1})`
    ... computes, for every configuration of the Q space, the fraction of the atomic state function it carries and
        the second-order energy it would contribute. The P space must already have been diagonalized, its energies
        and eigenvectors being passed in; the Q space is visited ONE CONFIGURATION AT A TIME and its CSFs are never
        held, so the memory is set by the P space and by the largest single configuration, not by |Q|.

        For each Q-CSF `q` and each P level `i` the weight is `|c|^2 = (<Psi_i|V|q> / D)^2` with the Epstein-Nesbet
        denominator `D = <q|H|q> - E_i`: the Q-CSF's OWN diagonal element against the CI eigenvalue, so that no
        configuration-average energy is needed and no H0/V partitioning has to be declared. The weight of a
        configuration is the MAXIMUM over the TARGET levels named by `targetLevels`, since one basis must serve
        them all -- but only them: a correlation level of P that nobody asked for must not enlarge the basis.

        A tuple `(weights, deltaE, c2)::Tuple{Dict{String,Float64}, Dict{String,Array{Float64,1}},
        Dict{String,Array{Float64,1}}}` is returned, keyed as `Hamiltonian.groupCsfsByConfiguration` keys its groups
        and holding per configuration the level-maximum weight, the per-level second-order energy and the per-level
        weight.
"""
function rankQspace(basis::Basis, pIndices::Array{Int64,1}, qGroups::Dict{String,Array{Int64,1}},
                     energies::Array{Float64,1}, vectors::Array{Float64,2}, nm::Nuclear.Model, grid::Radial.Grid,
                     settings::AsfSettings, cache::InteractionStrength.XLCache, targetLevels::Array{Int64,1})
    potential = Nuclear.nuclearPotential(nm, grid)
    nLev      = length(energies)
    weights   = Dict{String,Float64}()
    deltaE    = Dict{String,Array{Float64,1}}()
    c2        = Dict{String,Array{Float64,1}}()

    for  (key, qIdx)  in  qGroups
        dEconf = zeros(nLev);    c2conf = zeros(nLev)
        for  q in qIdx
            # the coupling of this one Q-CSF to every CSF of P, and its own diagonal element
            vRow = [ Hamiltonian.matrixElement(basis, p, q, nm, grid, settings, potential, cache)  for p in pIndices ]
            hqq  =   Hamiltonian.matrixElement(basis, q, q, nm, grid, settings, potential, cache)
            for  i = 1:nLev
                viq = sum( vectors[k,i] * vRow[k]  for k = 1:length(pIndices) )
                dd  = hqq - energies[i]
                # A vanishing denominator is an intruder rather than a divergence: the configuration is
                # degenerate with the level and belongs in P.  Give it an infinite weight so that it is
                # promoted by any threshold, and never divide by it.
                if  abs(dd) < 1.0e-10   c2conf[i] = Inf;    dEconf[i] = 0.;    continue    end
                c2conf[i] = c2conf[i] + (viq/dd)^2
                dEconf[i] = dEconf[i] - viq^2/dd
            end
        end
        weights[key] = maximum( c2conf[i]  for i in targetLevels );    deltaE[key] = dEconf;    c2[key] = c2conf
    end

    return( (weights, deltaE, c2) )
end


"""
`Hamiltonian.secondOrderBlock(sym::LevelSymmetry, basis::Basis, pBasis::Basis,
                              refConfigs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid,
                              settings::AsfSettings, treatment::Basics.SecondOrder; printout::Bool=true)`
    ... performs the CI of ONE symmetry block of the given basis with its Q space treated by the thresholds of the
        treatment, and is the body of `Hamiltonian.performCIwithSecondOrderQ`. Every
        configuration whose weight exceeds `promoteAbove` enters the CI and is diagonalized exactly, those below
        `discardBelow` are dropped, and the remainder is folded in through the second-order effective Hamiltonian.
        The P space is taken to be the CSFs belonging to the reference configurations.

        THE STEP REFUSES TO STAY SILENT WHEN ITS OWN PRECONDITION FAILS. `Sum |c|^2` over the folded part says
        whether that part really is a weak perturbation; where it is not, the computation still returns -- the
        multiplet is a legitimate CI in the space actually used -- but `isJustified` is false, a warning is printed,
        and the configurations that would have to be promoted are NAMED, in descending order with their running
        sum, so that a second attempt succeeds by reading the output rather than by guessing a threshold.

        A tuple `(levels::Array{Level,1}, partition::Hamiltonian.QPartition)` is returned.
"""
function secondOrderBlock(sym::LevelSymmetry, basis::Basis, pBasis::Basis, refConfigs::Array{Configuration,1},
                           nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                           treatment::Basics.SecondOrder; printout::Bool=true)
    cache     = InteractionStrength.XLCache()
    potential = Nuclear.nuclearPotential(nm, grid)

    # -- only the CSFs of THIS symmetry take part; the Hamiltonian is block diagonal in J^P, so a block may be
    #    partitioned and perturbed on its own, and the levels it returns carry its own J and parity.
    allIdx            = [ k for k = 1:length(basis.csfs)
                              if basis.csfs[k].J == sym.J  &&  basis.csfs[k].parity == sym.parity ]
    (groups, confOf)  = Hamiltonian.groupCsfsByConfiguration(basis, allIdx)
    # P IS THE PREVIOUS STEP'S SPACE, NOT THE REFERENCE.  A perturbative step sits at the end of a chain of
    # variational ones, and everything those built belongs in P -- the correlation they describe was decided by
    # physics and must not be re-decided by a threshold.  Taking P from the reference alone put the whole valence
    # correlation into Q, where the ranking promoted most of it back: the same answer for the wrong reason, and
    # dependent on the threshold rather than on the model.
    (pGroups, _)      = Hamiltonian.groupCsfsByConfiguration(pBasis, collect(1:length(pBasis.csfs)))
    refKeys           = collect(keys(pGroups))
    pIndices = Int64[];    qGroups = Dict{String,Array{Int64,1}}()
    for  (key, idx)  in  groups
        if  key in refKeys   append!(pIndices, idx)   else   qGroups[key] = idx   end
    end
    sort!(pIndices)
    if  length(pIndices) == 0
        error("Hamiltonian.performCIwithSecondOrderQ(): no CSF of this basis belongs to the reference " *
              "configurations, so there is no P space to perturb about.")
    end
    if  printout
        println("\n>> [second order, $(string(sym))] P = $(length(pIndices)) CSF of $(length(refKeys)) configuration(s) inherited from the previous step;  " *
                "Q = $(sum(length(v) for v in values(qGroups); init=0)) CSF in $(length(qGroups)) configurations.")
    end

    # -- the CI in P alone, which defines the zeroth order
    hPP = [ Hamiltonian.matrixElement(basis, r, s, nm, grid, settings, potential, cache)
            for r in pIndices, s in pIndices ]
    hPP = (hPP + transpose(hPP)) / 2
    fP  = Hamiltonian.diagonalizeCiMatrix(hPP, LevelSelection())
    ePl = fP.values;    cPl = hcat(fP.vectors...)

    # -- WHICH LEVELS DOES THE WEIGHT ANSWER FOR?  P now carries the whole previous variational space, so most of
    #    its levels are correlation states nobody asked for -- and taking the maximum over ALL of them lets an
    #    obscure level force a promotion the wanted levels do not need. Measured on Cl III: the maximum came from
    #    P level 5 at 3.6e-06 against level 1's 9.0e-07, and 9 of 11 configurations were promoted, leaving 875 CSF
    #    of 908 and no saving at all. The targets are therefore chosen by WEIGHT ON THE REFERENCE, the same rule
    #    and the same 0.5 as SelfConsistent.selectTargetLevelsEOL -- selection by index is unstable for exactly
    #    the reason given there, that a correlation level sinking below the reference silently redirects it.
    refKeysTarget = [ Hamiltonian.configurationKey(conf)  for conf in refConfigs ]
    refRows       = [ k  for (k, p) in enumerate(pIndices)
                          if Hamiltonian.configurationKey(
                                 Basics.extractConfiguration(Basics.FromBasis(), basis, basis.csfs[p])) in refKeysTarget ]
    targetLevels  = Int64[]
    for  i = 1:length(ePl)
        w = isempty(refRows) ? 1.0 : sum( cPl[k,i]^2  for k in refRows )
        if  w >= 0.5    push!(targetLevels, i)    end
    end
    if  isempty(targetLevels)
        @warn("Hamiltonian.secondOrderBlock(): no level of P carries a weight of 0.5 on the reference " *
              "configurations, so every level is treated as a target and the selection is CONSERVATIVE " *
              "(it will promote more than it needs).")
        targetLevels = collect(1:length(ePl))
    end
    if  printout
        println(">>   $(length(targetLevels)) of $(length(ePl)) levels in P carry the reference and set the weights.")
    end

    # -- rank every Q configuration, one at a time
    (weights, deltaE, c2) = Hamiltonian.rankQspace(basis, pIndices, qGroups, ePl, cPl, nm, grid, settings, cache,
                                                   targetLevels)

    # -- route each configuration by its own weight
    promoted = Configuration[];    folded = Configuration[];    discarded = Configuration[]
    promotedKeys = String[];        foldedKeys = String[]
    for  (key, w)  in  weights
        if      w >  treatment.promoteAbove    push!(promoted,  confOf[key]);    push!(promotedKeys, key)
        elseif  w >= treatment.discardBelow    push!(folded,    confOf[key]);    push!(foldedKeys,   key)
        else                                   push!(discarded, confOf[key])
        end
    end

    # -- the CI in P plus everything promoted, exactly
    newIdx = deepcopy(pIndices)
    for  key in promotedKeys    append!(newIdx, qGroups[key])    end
    sort!(newIdx)
    hNN = [ Hamiltonian.matrixElement(basis, r, s, nm, grid, settings, potential, cache)
            for r in newIdx, s in newIdx ]
    hNN = (hNN + transpose(hNN)) / 2

    # -- and the folded remainder, through the effective Hamiltonian of Gaigalas Eq. (37) with the symmetrised
    #    denominator of his Eq. (30), in Epstein-Nesbet form: each CSF's own diagonal element in place of a
    #    configuration average.  The correction is symmetric in (a,b) by construction, so the matrix stays
    #    Hermitian and a single diagonalization delivers energies and mixing coefficients together.
    hDiag = [ hNN[a,a]  for a = 1:length(newIdx) ]
    for  key in foldedKeys
        for  q in qGroups[key]
            vRow = [ Hamiltonian.matrixElement(basis, p, q, nm, grid, settings, potential, cache)  for p in newIdx ]
            hqq  =   Hamiltonian.matrixElement(basis, q, q, nm, grid, settings, potential, cache)
            for  a = 1:length(newIdx),  b = 1:length(newIdx)
                da = hqq - hDiag[a];    db = hqq - hDiag[b]
                if  abs(da) < 1.0e-10  ||  abs(db) < 1.0e-10    continue    end
                hNN[a,b] = hNN[a,b] - vRow[a]*vRow[b] * 0.5 * (1/da + 1/db)
            end
        end
    end
    fN  = Hamiltonian.diagonalizeCiMatrix(hNN, LevelSelection())

    # -- THE PROMOTION GUARD, AND IT MUST LOOK AT THE RESULT AND NOT AT THE INPUT.  `Sum |c|^2` polices the
    #    FOLDED remainder and is blind to the promoted part, which is where a bad configuration does the real
    #    damage: folding one costs a small error, PROMOTING one puts it into the variational CI.
    #      A FIRST VERSION TESTED THE DIAGONAL ELEMENTS of the promoted configurations and did NOT fire, because
    #    a variational minimum can lie far below every diagonal -- the spurious level is made by OFF-DIAGONAL
    #    mixing.  The honest test is therefore on the eigenvalues: an excitation OUT OF THE CORE must COST
    #    energy, so a level that carries almost none of the reference and yet lies BELOW the reference level is
    #    not a correction to it but an intruder.
    #      Found on Fe XV 3s^2 (12-Sep-2026): a 2p-hole configuration landed 10 Ha below the 3s^2 reference,
    #    which cannot be -- the hole must cost some 28 Ha.  The cause is that a variational layer's correlation
    #    orbitals are CONTRACTED pseudo-orbitals (its "4s" came out at <r> = 0.698 a.u., on top of the 3s
    #    valence) and this step reuses them to build core-hole CSFs, where they are unphysical.
    refRowsN = [ k  for (k, idx) in enumerate(newIdx)
                     if Hamiltonian.configurationKey(
                            Basics.extractConfiguration(Basics.FromBasis(), basis, basis.csfs[idx])) in refKeysTarget ]
    unsound  = Tuple{Configuration,Float64}[]
    if  !isempty(refRowsN)
        wRef  = [ sum(fN.vectors[i][k]^2  for k in refRowsN; init=0.0)  for i = 1:length(fN.values) ]
        iRef  = argmax(wRef)
        for  i = 1:length(fN.values)
            if  fN.values[i] < fN.values[iRef]  &&  wRef[i] < 0.1
                # name the configuration this intruder is built on: the one carrying most of its weight
                best = nothing;  bw = -1.0
                for  (key, idx)  in  qGroups
                    w = sum( fN.vectors[i][k]^2  for (k, ix) in enumerate(newIdx) if ix in idx; init=0.0 )
                    if  w > bw    bw = w;   best = confOf[key]   end
                end
                isnothing(best)  ||  push!(unsound, (best, fN.values[i]))
            end
        end
    end

    # -- the gate: is the folded part really a weak perturbation?
    nLev        = length(targetLevels)
    sumC2Folded = zeros(nLev);    dETotal = zeros(nLev)
    for  key in foldedKeys
        for (j, i) in enumerate(targetLevels)
            sumC2Folded[j] += c2[key][i];   dETotal[j] += deltaE[key][i]
        end
    end
    isJustified = all(sumC2Folded .<= treatment.promoteAbove)

    # The level vectors are expanded back onto the WHOLE basis -- zero on every CSF that was folded in or
    # discarded -- so that a level from here is indistinguishable from one performCI would return, and every
    # downstream routine that reads level.basis.csfs keeps working.
    levels = Level[]
    for  i = 1:length(fN.values)
        vector = zeros( length(basis.csfs) )
        for  (k, idx)  in  enumerate(newIdx)    vector[idx] = fN.vectors[i][k]    end
        push!(levels, Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, fN.values[i], 0., true,
                             basis, vector) )
    end
    weightList = sort( [ (confOf[key], w)  for (key, w) in weights ], by = x -> -x[2] )
    partition  = QPartition(promoted, folded, discarded, weightList, dETotal, sumC2Folded, isJustified, unsound)

    if  printout    Hamiltonian.displayQPartition(stdout, partition, treatment, length(newIdx))    end

    return( (levels, partition) )
end


"""
`Hamiltonian.performCIwithSecondOrderQ(basis::Basis, pBasis::Basis, refConfigs::Array{Configuration,1},
                                       nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                                       treatment::Basics.SecondOrder; printout::Bool=true)`
    ... performs the CI of the given basis, symmetry block by symmetry block, with the Q space of each treated by
        the thresholds of the treatment; see `Hamiltonian.secondOrderBlock` for what happens inside one block. The
        levels of every block are merged and sorted by energy, exactly as `Hamiltonian.performCI` does, so that the
        multiplet returned here is indistinguishable from one that routine would give.

        P IS `pBasis`, THE SPACE THE PREVIOUS STEP BUILT, and Q is what this basis adds beyond it. A perturbative
        step is therefore always a LEAF: it reports what the excluded excitations are still worth and is never
        inherited by a later variational step, so a promotion is internal to it -- some Q configurations simply
        treated better -- and needs no carrying forward.

        AND Q IS FINITE. It is the set of CSFs the CURRENT ORBITAL SET can build but the previous step excluded,
        not the complete second-order sum over an infinite virtual space. A correction from here says what this
        orbital set is still missing through the excitations left out; it does not say what a larger orbital set
        would add. That is why a new layer ENLARGES Q rather than shrinking it, and why the sequence of these
        corrections across layers is a convergence measure for the variational space.

        A tuple `(multiplet::Multiplet, partitions::Array{Hamiltonian.QPartition,1})` is returned, one partition
        per symmetry block, in the order the blocks were computed.
"""
function performCIwithSecondOrderQ(basis::Basis, pBasis::Basis, refConfigs::Array{Configuration,1},
                                    nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                                    treatment::Basics.SecondOrder; printout::Bool=true)
    symList = LevelSymmetry[]
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if  !(sym in symList)   push!(symList, sym)   end
    end

    blockMp = Multiplet[];    partitions = Hamiltonian.QPartition[]
    for  sym in symList
        (levels, part) = Hamiltonian.secondOrderBlock(sym, basis, pBasis, refConfigs, nm, grid, settings,
                                                      treatment; printout=printout)
        push!(blockMp, Multiplet(string(sym) * "+", levels));    push!(partitions, part)
    end
    mp = Basics.sortByEnergy( Basics.merge(blockMp) )

    return( (mp, partitions) )
end


"""
`Hamiltonian.displayQPartition(stream::IO, partition::Hamiltonian.QPartition, treatment::Basics.SecondOrder,
                                nCsfUsed::Int64)`
    ... reports how the Q space was divided and, where the fold-in was not justified, NAMES the configurations that
        would have to be promoted to make it so. Nothing is returned.
"""
function displayQPartition(stream::IO, partition::Hamiltonian.QPartition, treatment::Basics.SecondOrder,
                            nCsfUsed::Int64)
    println(stream, ">>   promoted into P  (|c|^2 > $(treatment.promoteAbove)) : " *
                    "$(length(partition.promoted)) configurations;  the CI uses $nCsfUsed CSF")
    println(stream, ">>   folded in                              : $(length(partition.folded)) configurations")
    println(stream, ">>   discarded        (|c|^2 < $(treatment.discardBelow)) : $(length(partition.discarded)) configurations")
    for  i = 1:length(partition.sumC2Folded)
        println(stream, @sprintf(">>   level %d:  Sum |c|^2 over the folded part = %.4e   second-order contribution = %.6e Ha",
                                 i, partition.sumC2Folded[i], partition.deltaE[i]))
    end

    if  !isempty(partition.unsound)
        println(stream, ">>   *** WARNING: $(length(partition.unsound)) level(s) of the CI lie BELOW the reference " *
                        "level while carrying almost none of the reference.  An excitation out of the core must")
        println(stream, ">>   COST energy, so these are INTRUDERS rather than corrections, and any level ordering " *
                        "taken from this step is unsafe.  This is not a threshold question: a variational")
        println(stream, ">>   layer's correlation orbitals are CONTRACTED pseudo-orbitals, and reusing them to " *
                        "build core-hole configurations is what produces it.  See priority item 21.")
        for  (conf, e)  in  sort(partition.unsound, by = x -> x[2])[1:min(5,end)]
            println(stream, @sprintf(">>     intruder at %.6f Ha, mostly   ", e) * string(conf))
        end
    end
    if  partition.isJustified
        println(stream, ">>   the folded remainder is a WEAK perturbation; second order is justified here.")
    else
        println(stream, ">>   *** WARNING: the folded remainder is NOT a weak perturbation.  Sum |c|^2 reaches " *
                        @sprintf("%.4e", maximum(partition.sumC2Folded)) * " against a promotion threshold of " *
                        "$(treatment.promoteAbove), so second order is NOT justified for this step and the")
        println(stream, ">>   energies above should not be used.  PROMOTE THE FOLLOWING and run again -- the running " *
                        "sum shows how far down the list is needed:")
        foldedKeys = [ Hamiltonian.configurationKey(c)  for c in partition.folded ]
        acc        = 0.
        for  (conf, w)  in  partition.weights
            Hamiltonian.configurationKey(conf) in foldedKeys   ||   continue
            acc = acc + w
            println(stream, @sprintf(">>     |c|^2 = %.4e   running sum %.4e   ", w, acc) * string(conf))
            if  acc > maximum(partition.sumC2Folded) - treatment.promoteAbove    break    end
        end
    end

    return( nothing )
end
