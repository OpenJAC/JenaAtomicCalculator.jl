
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
`Hamiltonian.foldQGroup!(hNN::Array{Float64,2}, newIdx::Array{Int64,1}, hDiag::Array{Float64,1}, basis::Basis,
                         qIndices::Array{Int64,1}, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                         potential::Radial.Potential, cache::InteractionStrength.XLCache,
                         perturbation::Basics.DiagonalCSFPerturbation)`
    ... folds ONE Q-space configuration's CSFs (`qIndices`) into `hNN`, the growing CI matrix of the retained
        (P + promoted) basis -- the body of `Hamiltonian.secondOrderBlock`'s own fold loop, taken out so that
        HOW a folded configuration's own internal coupling is treated can be dispatched on `perturbation`; see
        `Basics.AbstractQSpacePerturbation`. This method is `Basics.DiagonalCSFPerturbation`'s: every CSF of
        `qIndices` is folded as its own, isolated intermediate state, using its own diagonal Hamiltonian
        element as the Epstein-Nesbet denominator -- exactly what this loop has always done. Mutates `hNN` in
        place. Nothing is returned.
"""
function foldQGroup!(hNN::Array{Float64,2}, newIdx::Array{Int64,1}, hDiag::Array{Float64,1}, basis::Basis,
                      qIndices::Array{Int64,1}, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                      potential::Radial.Potential, cache::InteractionStrength.XLCache,
                      perturbation::Basics.DiagonalCSFPerturbation)
    for  q in qIndices
        vRow = [ Hamiltonian.matrixElement(basis, p, q, nm, grid, settings, potential, cache)  for p in newIdx ]
        hqq  =   Hamiltonian.matrixElement(basis, q, q, nm, grid, settings, potential, cache)
        for  a = 1:length(newIdx),  b = 1:length(newIdx)
            da = hqq - hDiag[a];    db = hqq - hDiag[b]
            if  abs(da) < 1.0e-10  ||  abs(db) < 1.0e-10    continue    end
            hNN[a,b] = hNN[a,b] - vRow[a]*vRow[b] * 0.5 * (1/da + 1/db)
        end
    end

    return( nothing )
end


"""
`Hamiltonian.foldQGroup!(hNN::Array{Float64,2}, newIdx::Array{Int64,1}, hDiag::Array{Float64,1}, basis::Basis,
                         qIndices::Array{Int64,1}, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                         potential::Radial.Potential, cache::InteractionStrength.XLCache,
                         perturbation::Basics.BlockDiagonalizedPerturbation)`
    ... the `Basics.BlockDiagonalizedPerturbation` method: `qIndices`' own CSF block is built and diagonalized
        FIRST (`Hamiltonian.matrixElement` to build it, `Hamiltonian.diagonalizeCiMatrix` to diagonalize it --
        the same two primitives `Hamiltonian.secondOrderBlock` already uses elsewhere), then the coupling
        vector to every retained CSF of `newIdx` is transformed into this configuration's OWN eigenbasis, and
        the SAME Gaigalas-Eq.-(37)-style fold is applied summed over its eigenstates -- their eigenvalues in
        place of each bare CSF's diagonal element -- rather than over its bare CSFs one at a time. This
        captures correlation WITHIN `qIndices` that the diagonal method cannot see. Mutates `hNN` in place.
        Nothing is returned.
"""
function foldQGroup!(hNN::Array{Float64,2}, newIdx::Array{Int64,1}, hDiag::Array{Float64,1}, basis::Basis,
                      qIndices::Array{Int64,1}, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                      potential::Radial.Potential, cache::InteractionStrength.XLCache,
                      perturbation::Basics.BlockDiagonalizedPerturbation)
    nq   = length(qIndices)
    hQQ  = [ Hamiltonian.matrixElement(basis, qIndices[i], qIndices[j], nm, grid, settings, potential, cache)
             for i = 1:nq, j = 1:nq ]
    hQQ  = (hQQ + transpose(hQQ)) / 2
    fQ   = Hamiltonian.diagonalizeCiMatrix(hQQ, LevelSelection())
    eQ   = fQ.values;    cQ = hcat(fQ.vectors...)
    vP   = [ Hamiltonian.matrixElement(basis, p, qIndices[i], nm, grid, settings, potential, cache)
             for p in newIdx, i = 1:nq ]
    vEig = vP * cQ

    for  k = 1:nq
        hqq = eQ[k]
        for  a = 1:length(newIdx),  b = 1:length(newIdx)
            da = hqq - hDiag[a];    db = hqq - hDiag[b]
            if  abs(da) < 1.0e-10  ||  abs(db) < 1.0e-10    continue    end
            hNN[a,b] = hNN[a,b] - vEig[a,k]*vEig[b,k] * 0.5 * (1/da + 1/db)
        end
    end

    return( nothing )
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
    #    denominator of his Eq. (30).  The correction is symmetric in (a,b) by construction, so the matrix stays
    #    Hermitian and a single diagonalization delivers energies and mixing coefficients together.  HOW each
    #    folded configuration's own CSFs are treated -- isolated (Epstein-Nesbet) or diagonalized together first
    #    -- is dispatched on treatment.qSpacePerturbation; see Hamiltonian.foldQGroup! and
    #    Basics.AbstractQSpacePerturbation.
    hDiag = [ hNN[a,a]  for a = 1:length(newIdx) ]
    for  key in foldedKeys
        Hamiltonian.foldQGroup!(hNN, newIdx, hDiag, basis, qGroups[key], nm, grid, settings, potential, cache,
                                 treatment.qSpacePerturbation)
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
    # THE TEST IS AGAINST `justifyBelow` AND NOT AGAINST `promoteAbove`, fixed 13-Sep-2026.  Judging the fold-in by
    # the same number that decides what gets folded made the gate LOOSEN as more was folded, so the one setting
    # that most needs policing -- promoteAbove = 1.0, "promote nothing, fold everything" -- was exactly the setting
    # at which it could never fire.  Measured on Fe VII [Ar] 3d^2: the step announced "a WEAK perturbation; second
    # order is justified here" while returning a 3F term inverted and displaced thirtyfold, with Sum |c|^2 = 0.048
    # printed beside it and not judged.  See priority item 28.
    isJustified = all(sumC2Folded .<= treatment.justifyBelow)

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
    # THE SPREAD OF THE CORRECTION ACROSS THE TARGET LEVELS, and for a SPECTRUM it decides more than the size does.
    # Core correlation is largely a common shift, which cancels out of every excitation energy; what survives is
    # the part that DIFFERS between levels.  On Fe VII the correction was -0.28 to -0.30 Ha -- a spread of some
    # 3300 cm^-1 across a term spanning 2329 cm^-1, which is how a 5 % fold-in inverted a multiplet.
    if  length(partition.deltaE) > 1
        spread = maximum(partition.deltaE) - minimum(partition.deltaE)
        println(stream, @sprintf(">>   second-order contribution: %+.6e to %+.6e Ha over the target levels, a SPREAD of ",
                                 minimum(partition.deltaE), maximum(partition.deltaE)) *
                        @sprintf("%.4e Ha = %.1f cm^-1.", spread,
                                 Defaults.convertUnits("energy: from atomic to Kayser", spread)))
        println(stream, ">>   For a SPECTRUM it is that spread and not the size that matters: a common shift cancels " *
                        "out of every excitation energy, and what is left over does not.")
    end
    if  partition.isJustified
        println(stream, ">>   the folded remainder is a WEAK perturbation; second order is justified here.")
    else
        println(stream, ">>   *** WARNING: the folded remainder is NOT a weak perturbation.  Sum |c|^2 reaches " *
                        @sprintf("%.4e", maximum(partition.sumC2Folded)) * " against the justification bound " *
                        "$(treatment.justifyBelow), so second order is NOT justified for this step and the")
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


#####################################################################################################################
#####  Second-order treatment of a SHAPE-CROSSING Q space (AtomicState.RasStep.shapedTo)  ##########################
#####################################################################################################################

"""
`struct  Hamiltonian.ShapeCrossingPartition`
    ... what one shape-crossing second-order junction (`AtomicState.RasStep.shapedTo`) is worth, returned by
        `Hamiltonian.secondOrderShapeCrossingCorrection`.

    + targetLevels   ::Array{Int64,1}    ... which levels of P (by index into the ladder's own priorMultiplet)
        this correction was computed for -- those carrying weight >= 0.5 on the reference configurations, the
        SAME rule `Hamiltonian.secondOrderBlock` uses, for the same reason: P carries the whole accumulated
        ladder, most of whose levels are correlation states nobody asked for.
    + deltaE         ::Array{Float64,1}  ... the second-order energy shift for each of `targetLevels`, in the
        SAME order.
    + newQIndices    ::Array{Int64,1}    ... which CSFs of Q's OWN basis are genuinely NEW -- not already
        present (by non-relativistic configuration key) in P -- and therefore the ones actually summed over.
        Q's basis always re-includes `refConfigs` and typically re-includes much of what P already built
        (`Basics.generateBasis` was written for a P-like space, not a bare "new excitations only" one), and
        coupling P against one of THOSE would be a near-diagonal, not a genuine off-diagonal, term.
    + coupling       ::Array{Float64,2}  ... the `<P level i|H-E|Q CSF s>` coupling, P's levels (rows) against
        `newQIndices` (columns, in the SAME order), for a caller who wants more than the target levels alone.
"""
struct  ShapeCrossingPartition
    targetLevels    ::Array{Int64,1}
    deltaE          ::Array{Float64,1}
    newQIndices     ::Array{Int64,1}
    coupling        ::Array{Float64,2}
end


function Base.show(io::IO, part::Hamiltonian.ShapeCrossingPartition)
    println(io, "ShapeCrossingPartition: second-order correction for $(length(part.targetLevels)) target " *
                "level(s) of P against $(length(part.newQIndices)) genuinely new CSF(s) of Q.")
    for  (j, i)  in  enumerate(part.targetLevels)
        println(io, @sprintf("  P level %d:  deltaE = %+.6e Ha  (%+.4f mHa)", i, part.deltaE[j], 1000*part.deltaE[j]))
    end
end


"""
`Hamiltonian.foldNearDegenerateQCluster(ePLevel::Float64, hqqCluster::Array{Float64,1}, couplingCluster::Array{Float64,1})`
    ... folds a P target level together with a cluster of Q-CSFs whose Epstein-Nesbet denominator with it is too
        small for ordinary non-degenerate perturbation theory to be trusted (see
        `secondOrderShapeCrossingCorrection`'s own classification, priority item 5's dense-Q half). Builds the
        small `(1+n) x (1+n)` matrix with `ePLevel` and `hqqCluster` on the diagonal and `couplingCluster` in the
        first row/column -- Q-Q off-diagonals are zero, since no such matrix element is ever computed by the
        caller -- diagonalizes it, and returns the shift of the eigenvalue closest to `ePLevel`, relative to
        `ePLevel` itself. Validated against a hand-built toy case before being wired in: reduces to the ordinary
        `-coupling^2/dd` sum to machine precision in the well-separated limit, matches the closed-form 2x2
        degenerate-PT result for a single flagged CSF, and matches a direct external diagonalization for a
        multi-CSF cluster. A `shift::Float64` is returned.
"""
function foldNearDegenerateQCluster(ePLevel::Float64, hqqCluster::Array{Float64,1}, couplingCluster::Array{Float64,1})
    n = length(hqqCluster)
    M = zeros(1+n, 1+n)
    M[1,1] = ePLevel
    for  s = 1:n
        M[1,1+s] = couplingCluster[s];   M[1+s,1] = couplingCluster[s]
        M[1+s,1+s] = hqqCluster[s]
    end
    vals = eigvals(Symmetric(M))
    iMin = argmin(abs.(vals .- ePLevel))

    return( vals[iMin] - ePLevel )
end


"""
`Hamiltonian.secondOrderShapeCrossingCorrection(pMp::Multiplet, qMp::Multiplet, refConfigs::Array{Configuration,1},
                                                 grid::Radial.Grid, nm::Nuclear.Model, settings::AsfSettings)`
    ... the shape-crossing counterpart of `Hamiltonian.secondOrderBlock`, for a Q space built by a
        `AtomicState.RasStep.shapedTo` layer (`BasicsAZ.retagShapedSubshells`/`generateOrbitalsForShapedBasis`
        -- its own, separately SCF-refined `Basis`, never merged into P). P and Q are biorthogonalized via
        `BiOrthogonal.computeTransformation(BiOrthogonal.ShapeCrossingJunction(), pMp, qMp, grid)`, and every
        P-Q coupling then goes through `ManyElectron.matrixElement_HminusE` on the RESULTING biorthonormal
        bases, summed over P's own CI mixing coefficients -- never a bare, unweighted CSF-pair value, which
        Phase 0's toy validation (`work/clIII-corevalence/phase0-toy/validate2.jl`; see the module docstring
        of `BiOrthogonal.jl`) showed is not a well-defined quantity on its own.

        ONLY Q'S GENUINELY NEW CSFs ENTER THE SUM -- those whose non-relativistic configuration is NOT already
        present in P. `Basics.generateBasis` (used to build Q in `BasicsAZ.runRasStep!`) was written for a
        P-LIKE space and always re-includes `refConfigs`, and typically much of P's own valence-correlation
        excitation set besides (it has no "new excitations only" mode); treating Q's WHOLE CSF space, or its
        own full diagonalization, as "the Q space" therefore couples P against configurations P ALREADY HAS.
        For such a pair the biorthogonal overlap after transformation is exactly 1 (same spectators, same
        shape, no reshaping at that label), so the coupling is really a DIAGONAL-like term, equal to roughly
        the CSF's own total energy (`ManyElectron.matrixElement_HminusE` carries no `-E<r|s>` term -- see its
        own docstring -- which is harmless only where the overlap is exactly 0, as it is for a genuinely new
        CSF). Measured while building this function: including P's own CSFs on the Q side gave a "correction"
        of +5.7e7 Ha; restricting to the genuinely new CSFs, as this function now does, is what removes it.

        For each of Q's new CSFs `s` the Epstein-Nesbet denominator is ITS OWN diagonal Hamiltonian element,
        computed on Q's OWN, untransformed, internally-orthonormal basis (a single-basis quantity, exactly
        `Hamiltonian.rankQspace`'s own per-CSF Epstein-Nesbet formula) -- NOT on the biorthogonal-transformed
        one, since the transform is not norm-preserving in general (`BiOrthogonal.jl`'s own module docstring)
        and the denominator needs Q's true, physical self-energy. The second-order shift for target P level
        `i` is `- sum_s |<Psi_Pi|H-E|s>|^2 / (h_ss - E_Pi)`, the same sign convention as `Hamiltonian.rankQspace`
        and `Hamiltonian.foldQGroup!` -- a state lying ABOVE the reference (`h_ss > E_Pi`, the normal case)
        LOWERS it, as a correlation correction must.

        ONLY THE TARGET LEVELS OF P ARE CORRECTED -- those carrying weight >= 0.5 on `refConfigs`, the SAME
        rule `Hamiltonian.secondOrderBlock` uses and for the same reason.

        THE DENOMINATOR GUARD IS A HARD EXCLUSION AT `abs(dd) < 1.0e-10` (literal degeneracy, UNCHANGED) PLUS
        A QUASI-DEGENERATE FOLD, added 05-Nov-2026 -- filed as priority item 5, 04-Nov-2026, after a dense Q
        space (Ti III 3d^2 -> 4d, 5371 genuinely-new CSFs against 9 target levels) gave a correction swinging
        -386 to -537 mHa, LARGER than the whole NIST term span, with no warning printed anywhere. For each
        target level, every surviving Q-CSF `s` is classified by `abs(coupling[i,s]/dd)`: non-degenerate PT is
        trusted (`ratioThreshold`, default 0.4, below which a Taylor expansion in the coupling is a fair
        description) and stays in the ordinary `-coupling^2/dd` sum; above it, the CSF is FLAGGED and its
        whole flagged cluster for that level is folded instead -- `foldNearDegenerateQCluster` diagonalizes
        the small matrix of the target level plus its flagged Q-CSFs (Q-Q off-diagonals zero, since this
        function never computes them) and takes the eigenvalue closest to the unperturbed level as the
        correction, which stays bounded as `dd -> 0` where the ordinary sum diverges.

        THE FOLD IS VALIDATED AND HARMLESS, BUT DOES NOT REPAIR THE TI III CASE IT WAS BUILT FOR -- state this
        plainly rather than letting the mechanism's correctness stand in for the fix having worked. Measured
        05-Nov-2026: on the real Ti III 3d^2 run, `abs(coupling[i,s]/dd)` never exceeds `ratioThreshold` for a
        SINGLE one of the 9x5371 pairs, so nothing is ever flagged, the `@warn` below never fires, and the
        nine `deltaE` values are UNCHANGED from before this fold existed (-386.2 to -537.4 mHa, still larger
        than the whole NIST term span). The diagnosis behind this fold -- a few outlier near-degenerate
        denominators -- was WRONG for this case: the blowup is a DENSE ACCUMULATION of many only-moderately-
        small denominators across those 5371 CSFs, none individually near-degenerate enough to trip a ratio
        test, which is a different failure mode that a per-term classifier cannot catch by construction. The
        dense-Q half of item 5 is UNFIXED and is held on the challenges list (see `src/` `CLAUDE.md`'s working
        lists) rather than carried here as an open priority item; this fold stays in `src/` because it is
        correct, toy-validated (see the helper's own docstring) and a byte-identical no-op everywhere tested
        (Cl III's sparse p^3 Q space, branch g: nothing is ever flagged there either, and the five dated
        `deltaE` values reproduce exactly), so it is ready for whatever genuinely isolated near-degeneracy
        turns up in a future case -- it is simply not that case.

        P MAY CARRY MORE SUBSHELLS THAN Q -- priority item 5's second half, also closed 04-Nov-2026. A prior
        RAS layer can deepen P's own ladder beyond the fixed subshell set Q's own `shapedTo` step touches
        (branch h's n3+n4 ladder is exactly this). Q's subshells are now remapped by (n,kappa) SEARCH and
        then RE-EXPANDED onto P's (possibly longer) subshell list via `ManyElectron.Basis(basis, subshells)`,
        the already-existing zero-padding re-expansion constructor, rather than requiring the two lists to
        match length-for-length and position-for-position. This is exact, not approximate: every CSF of Q
        genuinely has zero electrons on the extra slots, and `ManyElectron.matrixElement_HminusE`'s angular
        coefficients are driven by CSF OCCUPATION, never by subshell-list length, so a slot no CSF occupies
        is never looked up. A no-op (bit-identical) when P and Q already agree, as on branch g.

        A `Hamiltonian.ShapeCrossingPartition` is returned.
"""
function secondOrderShapeCrossingCorrection(pMp::Multiplet, qMp::Multiplet, refConfigs::Array{Configuration,1},
                                             grid::Radial.Grid, nm::Nuclear.Model, settings::AsfSettings)
    # WHICH OF Q'S OWN CSFs ARE GENUINELY NEW, decided on the ORIGINAL (untransformed) bases -- configuration
    # keys do not depend on orbitals at all, so the biorthogonal transform below has nothing to do with this.
    pBasis0 = pMp.levels[1].basis;    qBasis0 = qMp.levels[1].basis
    (pGroups, _) = Hamiltonian.groupCsfsByConfiguration(pBasis0, collect(1:length(pBasis0.csfs)))
    pKeys        = Set(keys(pGroups))
    newQIndices  = Int64[ s  for s = 1:length(qBasis0.csfs)
                             if !(Hamiltonian.configurationKey(
                                      Basics.extractConfiguration(Basics.FromBasis(), qBasis0, qBasis0.csfs[s])) in pKeys) ]
    if  isempty(newQIndices)
        error("Hamiltonian.secondOrderShapeCrossingCorrection(): every CSF of Q's basis already belongs to a " *
              "configuration P has -- there is nothing genuinely new to couple against, which likely means " *
              "this step's *To lists never reached beyond what P already built.")
    end

    # THE DENOMINATOR, on Q's OWN untransformed basis -- a single-basis diagonal element, so no cross-basis issue.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    op         = Basics.HminusE(settings.eeInteractionCI)
    hqq        = [ ManyElectron.matrixElement_HminusE(op, qBasis0, s, qBasis0, s, grid, nuclearPot)  for s in newQIndices ]

    # THE NUMERATOR needs the biorthogonal transform, since P's and Q's reshaped orbital genuinely differ.
    newPMp, newQMp = JenaAtomicCalculator.BiOrthogonal.computeTransformation(
                         JenaAtomicCalculator.BiOrthogonal.ShapeCrossingJunction(), pMp, qMp, grid)
    pBasis = newPMp.levels[1].basis;    qBasis = newQMp.levels[1].basis

    # `ManyElectron.matrixElement_HminusE` generates its angular coefficients from `leftBasis.subshells` ALONE
    # (see its own code) and then looks each `coeff.a`/`coeff.b` up in BOTH orbital dicts -- correct wherever the
    # two bases share the SAME Subshell objects label for label, which every OTHER BiOrthogonal consumer does,
    # since nothing there ever differs in `.shape`. A shape-crossing P/Q pair is the first case where it does:
    # pBasis's reshaped label is `SCFShape`, qBasis's own is `CoreValenceCorrelationShape`, so a direct call
    # raises a KeyError -- qBasis's dict genuinely holds the orbital, just under the OTHER tag. Fixed by
    # PRESENTING qBasis's orbitals under pBasis's own subshell objects for this call only (same physical
    # orbitals, same (n,kappa) label), matched by (n,kappa) SEARCH rather than assumed position: this is a
    # pure relabelling for the angular machinery's benefit, not a change to any radial content.

    # P MAY CARRY MORE SUBSHELLS THAN Q (priority item 5, 04-Nov-2026: a prior RAS layer deepened P's own
    # ladder beyond the fixed core-valence excitation space Q's own step built) -- so Q's subshell list is
    # remapped label-for-label, NOT forced to `pBasis.subshells` wholesale, and then RE-EXPANDED onto the
    # (possibly longer) `pBasis.subshells` via `ManyElectron.Basis(basis, subshells)`, the already-existing
    # zero-padding re-expansion constructor. This is a no-op when the two already agree (its own first line
    # returns `basis` unchanged), so the equal-length case (branch g) is untouched; it is exact rather than
    # approximate on the extra slots, since every CSF of Q genuinely has ZERO electrons there, and the
    # occupation-driven coefficient generation never needs an orbital at a subshell no CSF occupies.
    remapToP(sh::Subshell) = pBasis.subshells[ findfirst(s -> (s.n,s.kappa) == (sh.n,sh.kappa), pBasis.subshells) ]
    if  any( findfirst(s -> (s.n,s.kappa) == (sh.n,sh.kappa), pBasis.subshells) === nothing  for sh in qBasis.subshells )
        error("Hamiltonian.secondOrderShapeCrossingCorrection(): some (n,kappa) label of Q's subshells " *
              "$(qBasis.subshells) is absent from P's $(pBasis.subshells) altogether -- P must carry at " *
              "least every label Q's own excitation space touches, even if P carries others besides.")
    end
    qSubshellsPresented     = [ remapToP(sh)  for sh in qBasis.subshells ]
    qOrbitalsPresented      = Dict{Subshell,Orbital}( remapToP(sh) => orb  for (sh, orb) in qBasis.orbitals )
    qCoreSubshellsPresented = [ remapToP(sh)  for sh in qBasis.coreSubshells ]
    qBasisPresented = Basis(qBasis.isDefined, qBasis.NoElectrons, qSubshellsPresented, qBasis.csfs,
                             qCoreSubshellsPresented, qOrbitalsPresented)
    qBasis = ManyElectron.Basis(qBasisPresented, pBasis.subshells)

    nP = length(pBasis.csfs);    nNewQ = length(newQIndices)
    M  = [ ManyElectron.matrixElement_HminusE(op, pBasis, r, qBasis, newQIndices[j], grid, nuclearPot)
           for r = 1:nP, j = 1:nNewQ ]

    # THE FULL LEVEL-TO-CSF COUPLING is one matrix product away once M is in hand -- no new physics here, just
    # P's own mixing coefficients (unchanged by the transform's own counter-rotation, which acts on newPMp's
    # mc, not on pMp's own).
    mcP      = hcat( [ lv.mc  for lv in newPMp.levels ]... )
    coupling = mcP' * M

    refKeysTarget = [ Hamiltonian.configurationKey(conf)  for conf in refConfigs ]
    targetLevels  = Int64[]
    for  i = 1:length(pMp.levels)
        w = sum( pMp.levels[i].mc[r]^2  for r = 1:length(pBasis0.csfs)
                 if  Hamiltonian.configurationKey(
                         Basics.extractConfiguration(Basics.FromBasis(), pBasis0, pBasis0.csfs[r])) in refKeysTarget )
        if  w >= 0.5    push!(targetLevels, i)    end
    end
    if  isempty(targetLevels)
        @warn("Hamiltonian.secondOrderShapeCrossingCorrection(): no level of P carries a weight of 0.5 on the " *
              "reference configurations, so every level is treated as a target and the selection is " *
              "CONSERVATIVE (it will promote more than it needs).")
        targetLevels = collect(1:length(pMp.levels))
    end

    deltaE = zeros(length(targetLevels))
    anatomy = Dict{Int64,Tuple{Int64,Float64,Int64,Float64,Float64,Float64}}()   # see the block that fills it
    ratioThreshold = 0.4   # |coupling/dd| above this is folded rather than summed, item 5 -- see the module docstring
    nFolded = 0;    worstRatio = 0.0;    worstDenominator = Inf
    for  (j, i)  in  enumerate(targetLevels)
        ePLevel      = pMp.levels[i].energy
        unflaggedIdx = Int64[];    flaggedIdx = Int64[]
        for  s = 1:nNewQ
            dd = hqq[s] - ePLevel
            abs(dd) < 1.0e-10   &&   continue                       # literal degeneracy, UNCHANGED hard exclusion
            if  abs(coupling[i,s] / dd) > ratioThreshold
                push!(flaggedIdx, s)
                if  abs(coupling[i,s]/dd) > worstRatio
                    worstRatio = abs(coupling[i,s]/dd);    worstDenominator = dd
                end
            else
                push!(unflaggedIdx, s)
            end
        end
        for  s  in  unflaggedIdx
            deltaE[j] -= coupling[i,s]^2 / (hqq[s] - ePLevel)
        end
        if  !isempty(flaggedIdx)
            nFolded += length(flaggedIdx)
            deltaE[j] += foldNearDegenerateQCluster(ePLevel, hqq[flaggedIdx], coupling[i,flaggedIdx])
        end
        # THE ANATOMY OF THE SUM, so that a blowup cannot arrive silently.  A per-term ratio test says whether
        # any ONE denominator is dangerous;  it cannot say whether the TOTAL is made of a few large terms or of
        # thousands of moderate ones, and those two call for opposite remedies.  So the share carried by the
        # largest terms and the number of terms needed to reach nine tenths of the total are recorded here:  a
        # handful of terms carrying most of it is an outlier problem, which foldNearDegenerateQCluster exists
        # for;  a tail in which no term matters much and the count does the work is a VOLUME effect, which no
        # reclassification of individual terms can reach.
        if  !isempty(unflaggedIdx)
            terms = sort!( [ -coupling[i,s]^2 / (hqq[s] - ePLevel)  for s in unflaggedIdx ],
                           by = abs, rev = true )
            tot   = sum(abs, terms)
            n90   = 0;    acc = 0.0
            for  t  in  terms
                acc += abs(t);    n90 += 1
                acc >= 0.9 * tot   &&   break
            end
            # AND THE DOMINANT TERM'S OWN COUPLING AND DENOMINATOR, because the two shapes of outlier call
            # for opposite remedies and only these two numbers separate them:  a SMALL denominator is a
            # near-degeneracy, which folding repairs;  a LARGE coupling at an ordinary denominator is a state
            # that does not belong in a perturbative Q space at all, and no denominator treatment reaches it.
            top  = argmax( [ abs(-coupling[i,s]^2 / (hqq[s] - ePLevel))  for s in unflaggedIdx ] )
            sTop = unflaggedIdx[top]
            dds  = sort( [ abs(hqq[s] - ePLevel)  for s in unflaggedIdx ] )   # Statistics is not imported here,
            ddMed = dds[ (length(dds)+1) >> 1 ]                                   # so the median is taken directly
            anatomy[j] = (length(terms), abs(terms[1])/tot, n90, coupling[i,sTop], hqq[sTop] - ePLevel, ddMed)
        end
    end

    # AND THE GUARD THAT THE TI III CASE ASKED FOR (challenge S34):  a second-order correction LARGER THAN THE
    # SPAN OF THE LEVELS IT CORRECTS is a breakdown of the expansion whatever its individual terms look like,
    # because the quantity being corrected has become smaller than its own correction.  Measured on Ti III
    # 3d^2 -> 4d: deltaE between -386 and -537 mHa against a computed term span of a few tens of mHa, so some
    # twenty times over -- and nothing printed a word about it.  The span is taken over the TARGET levels, the
    # same ones deltaE is reported for, so the comparison needs nothing from outside this function.
    if  length(targetLevels) > 1
        eTargets = [ pMp.levels[i].energy  for i in targetLevels ]
        span     = maximum(eTargets) - minimum(eTargets)
        worst    = maximum(abs, deltaE)
        if  span > 0.0  &&  worst > span
            @warn("Hamiltonian.secondOrderShapeCrossingCorrection(): the second-order correction reaches " *
                  "$(round(worst*1000, digits=2)) mHa while the target levels span only " *
                  "$(round(span*1000, digits=2)) mHa, i.e. a factor of $(round(worst/span, digits=1)).  A " *
                  "correction larger than the structure it corrects is NOT a usable number:  the expansion " *
                  "has broken down, and no choice of denominator floor or fold repairs that.  The anatomy " *
                  "below says which kind of breakdown it is.  See challenge S34.")
            for  (j, i)  in  enumerate(targetLevels)
                haskey(anatomy, j) || continue
                (nT, topShare, n90, cTop, ddTop, ddMed) = anatomy[j]
                @info("   level $i (J = $(pMp.levels[i].J)):  deltaE = " *
                      "$(round(deltaE[j]*1000, digits=2)) mHa from $nT terms;  the largest carries " *
                      "$(round(topShare*100, digits=2)) % and $n90 terms are needed for 90 %;  that " *
                      "largest term has coupling $(round(cTop, digits=5)) Ha over a denominator of " *
                      "$(round(ddTop, digits=5)) Ha (median |dd| = $(round(ddMed, digits=4)) Ha) -- " *
                      (n90 > 0.05*nT ? "a VOLUME effect, so truncating or shifting the whole sum is the only " *
                                       "lever;  reclassifying individual terms cannot reach it"
                                     : "an OUTLIER effect, which foldNearDegenerateQCluster is the mechanism for"))
            end
        end
    end
    if  nFolded > 0
        @warn("Hamiltonian.secondOrderShapeCrossingCorrection(): $nFolded of $(length(targetLevels)*nNewQ) " *
              "P-level/Q-CSF denominator(s) were too small for ordinary perturbation theory (worst " *
              "|coupling/dd| = $worstRatio at dd = $worstDenominator Ha) and were FOLDED (quasi-degenerate " *
              "diagonalization, see foldNearDegenerateQCluster) rather than summed term by term -- see " *
              "priority item 5.  Q-Q coupling WITHIN a flagged cluster is not included, since this function " *
              "never computes it;  if a flagged cluster is a large fraction of Q, that omission may itself " *
              "become the dominant remaining error, which would point at real multi-reference work rather " *
              "than this fold.  Inspect `hqq` and `coupling` directly before trusting `deltaE` here.")
    end

    return( ShapeCrossingPartition(targetLevels, deltaE, newQIndices, coupling) )
end
