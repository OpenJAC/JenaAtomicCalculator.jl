
"""
`struct  AdjustedLevels.RadialParameter`
    ... one fitted scaling parameter for ONE automatically-discovered integral class -- the sibling of
        `AdjustedLevels.ScalingParameter` for the orthogonal route: `fitScaleFactor` always returns exactly
        one, uniform across every rank; `fitOrthogonalParameters` returns one PER multipole rank it finds
        among this block's two-body coefficients (see `AdjustedLevels.buildParameterMatrices`).

    + parameterClass   ::String    ... which integral class this factor multiplies, e.g. "rank k=1".
    + value             ::Float64  ... the fitted factor itself; NOT 1.0 by construction unless that class
                                        turned out to be exactly undetermined by the given targets.
"""
struct  RadialParameter
    parameterClass      ::String
    value               ::Float64
end


"""
`struct  AdjustedLevels.OrthogonalParameter`
    ... one ORTHOGONAL combination of `RadialParameter`s -- a direction in parameter space, found from the
        measured targets themselves (see `AdjustedLevels.fitOrthogonalParameters`'s own docstring for how),
        not stored or universal. `value` is the fitted SHIFT along this direction, away from the uniform,
        unadjusted (every `RadialParameter` = 1) starting point -- NOT an absolute scale like
        `ScalingParameter.value`, since a direction is only meaningful as a combination of DEVIATIONS.

    + label         ::String                         ... an arbitrary but stable name, e.g. "mode 1".
    + value         ::Float64                         ... the fitted shift along this direction.
    + composition   ::Array{Tuple{String,Float64},1}  ... which `RadialParameter`s make up this direction,
                                                            and with what weight; together the weights form
                                                            one column of an orthogonal (unit-norm) matrix.
    + determinacy   ::Float64                          ... how strongly the given targets pin this direction
                                                             down -- the matching singular value of the
                                                             sensitivity matrix. Near zero means the targets
                                                             cannot tell this combination from no change at
                                                             all, and `value` was accordingly left at 0.0.
"""
struct  OrthogonalParameter
    label           ::String
    value           ::Float64
    composition     ::Array{Tuple{String,Float64},1}
    determinacy     ::Float64
end


"""
`struct  AdjustedLevels.FittedMultipletOrthogonal`
    ... the result of `AdjustedLevels.fitOrthogonalParameters`. The sibling of `AdjustedLevels.FittedMultiplet`
        for the orthogonal route -- a DIFFERENT type, deliberately, so that a reader (or a method dispatch)
        can tell at a glance whether one or several parameters were fit, without inspecting array lengths.

    + multiplet             ::ManyElectron.Multiplet         ... the single-symmetry-block multiplet from the
                                                                   FITTED, re-diagonalized Hamiltonian.
    + radialParameters      ::Array{RadialParameter,1}        ... the fitted factors, one per integral class,
                                                                    translated back to the familiar language.
    + orthogonalParameters  ::Array{OrthogonalParameter,1}     ... what was actually fit -- the orthogonal
                                                                     directions, each with its own determinacy.
    + fitTargets            ::Array{FitTarget,1}                ... the measured levels this fit was matched
                                                                      against.
    + rmsDeviation          ::Float64                           ... root-mean-square deviation, fitted vs
                                                                      measured, over the fit targets, in cm^-1,
                                                                      from the EXACT re-diagonalized spectrum
                                                                      (not the linearized estimate).
    + source                ::String                           ... a dated citation for the measured data used.
"""
struct  FittedMultipletOrthogonal
    multiplet               ::ManyElectron.Multiplet
    radialParameters         ::Array{RadialParameter,1}
    orthogonalParameters     ::Array{OrthogonalParameter,1}
    fitTargets               ::Array{FitTarget,1}
    rmsDeviation              ::Float64
    source                    ::String
end


"""
`Base.show(io::IO, fm::AdjustedLevels.FittedMultipletOrthogonal)`
    ... ALWAYS prints the "fitted, not ab initio" disclosure first, exactly as the single-parameter route's
        `Base.show(io::IO, ::AdjustedLevels.FittedMultiplet)` does, plus the orthogonal breakdown and each
        direction's determinacy, so a near-zero determinacy is visible without inspecting the struct by hand.
        Nothing is returned.
"""
function  Base.show(io::IO, fm::FittedMultipletOrthogonal)
    println(io, "FittedMultipletOrthogonal -- FITTED, NOT AB INITIO.")
    println(io, "  source: $(fm.source)")
    for  rp in fm.radialParameters
        println(io, "  radial parameter [$(rp.parameterClass)] = $(rp.value)")
    end
    for  op in fm.orthogonalParameters
        flag = op.determinacy < 1.0e-6 ? "  (UNDETERMINED by these targets)" : ""
        println(io, "  orthogonal $(op.label): shift = $(op.value), determinacy = $(op.determinacy)$flag")
    end
    println(io, "  fit against $(length(fm.fitTargets)) measured level(s), rms deviation = $(fm.rmsDeviation) cm^-1")
    println(io, "  multiplet: $(fm.multiplet.name), $(length(fm.multiplet.levels)) level(s)")
end


"""
`AdjustedLevels.buildParameterMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)`
    ... the orthogonal route's generalization of `AdjustedLevels.buildOneAndTwoBodyMatrices`: instead of
        summing every two-body (Coulomb) contribution into ONE lumped matrix, each contribution is kept in
        the matrix belonging to its own multipole RANK (`coeff.nu`) -- the same grouping Cowan's RCE uses to
        separate F^k (direct) from G^k (exchange), since for two orbitals of different l these occupy
        disjoint rank sets automatically; no direct/exchange bookkeeping is needed beyond the rank itself.
        A (idxCsf, oneBody, parameterClasses, matrices) tuple is returned: idxCsf::Array{Int64,1} as in
        `buildOneAndTwoBodyMatrices`; oneBody::Array{Float64,2} the shared one-body matrix; parameterClasses
        ::Array{String,1} a label per discovered rank; matrices::Array{Array{Float64,2},1} the matching
        upper-triangular two-body matrices, one per class, same convention as `buildOneAndTwoBodyMatrices`.
"""
function  buildParameterMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)
    idxCsf = [idx for idx = 1:length(basis.csfs) if basis.csfs[idx].J == sym.J && basis.csfs[idx].parity == sym.parity]
    n = length(idxCsf)
    potential = Nuclear.nuclearPotential(nm, grid)
    oneBody = zeros(Float64, n, n)
    twoBodyByRank = Dict{Int64,Array{Float64,2}}()

    for  r = 1:n
        for  s = r:n
            subshellList = basis.subshells
            opa  = SpinAngular.OneParticleOperator(0, plus)
            waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[idxCsf[r]], basis.csfs[idxCsf[s]], subshellList)
            opa  = SpinAngular.TwoParticleOperator(0, plus)
            waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[idxCsf[r]], basis.csfs[idxCsf[s]], subshellList)

            me1 = 0.0
            for  coeff in waG1
                me1 = me1 + coeff.T * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
            end
            oneBody[r,s] = me1

            for  coeff in waG2
                contribution = coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                          basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                if  !haskey(twoBodyByRank, coeff.nu)   twoBodyByRank[coeff.nu] = zeros(Float64, n, n)   end
                twoBodyByRank[coeff.nu][r,s] += contribution
            end
        end
    end

    ranks = sort(collect(keys(twoBodyByRank)))
    parameterClasses = [ "rank k=$(k)" for k in ranks ]
    matrices         = [ twoBodyByRank[k] for k in ranks ]

    return( (idxCsf, oneBody, parameterClasses, matrices) )
end


"""
`AdjustedLevels.fitOrthogonalParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                          targets::Array{Tuple{Int64,Float64},1}, source::String;
                                          determinacyThreshold::Float64=1.0e-6)`
    ... the orthogonal (Uylings/Raassen-style) route: discovers the integral classes present in `sym`'s block
        via `buildParameterMatrices`, fits one scaling parameter PER class against `targets` (same convention
        as `fitScaleFactor`'s own `targets`), and reports the fit in terms of the ORTHOGONAL combinations the
        targets can actually determine. THE FIT ITSELF IS A FIRST-ORDER (LINEARIZED) LEAST-SQUARES STEP: a
        sensitivity matrix J[target,class] = d(energy)/d(parameter), evaluated once at every class = 1 via the
        ab initio eigenvectors (no extra diagonalization), is decomposed by `LinearAlgebra.svd`; the fitted
        correction is the MINIMUM-NORM solution of `J * delta = measured - ab initio`, which leaves any
        direction the targets cannot determine untouched (determinacy ~ 0) rather than guessing. The spectrum
        is then re-diagonalized EXACTLY at the resulting parameters, so `rmsDeviation` reflects the real fit,
        not the linear estimate. This is honest for the SMALL-correction regime a spectroscopic fit usually
        is in; it is not a general nonlinear optimizer, and a determinacy near zero is reported rather than
        hidden. A fm::AdjustedLevels.FittedMultipletOrthogonal is returned.
"""
function  fitOrthogonalParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                    targets::Array{Tuple{Int64,Float64},1}, source::String;
                                    determinacyThreshold::Float64=1.0e-6)
    isempty(targets)   &&   error("AdjustedLevels.fitOrthogonalParameters(): need at least one target level to fit against.")
    idxCsf, oneBody, parameterClasses, matrices = buildParameterMatrices(sym, basis, nm, grid)
    nParam = length(matrices)
    nParam == 0   &&   error("AdjustedLevels.fitOrthogonalParameters(): no two-body (Coulomb) contributions were found for this block.")
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)

    # ab initio (every parameter = 1): one diagonalization gives both the reference spectrum and the eigenvectors
    # the sensitivity matrix needs.
    totalAbInitio = copy(oneBody)
    for  mtx in matrices   totalAbInitio .+= mtx   end
    eig0       = Hamiltonian.diagonalizeCiMatrix(totalAbInitio, LevelSelection())
    order0     = sortperm(eig0.values)
    energies0  = eig0.values[order0]
    relKayser0 = (energies0 .- energies0[1]) ./ auPerKayser

    nTargets  = length(targets)
    J         = zeros(Float64, nTargets, nParam)
    residual0 = zeros(Float64, nTargets)
    for  (row, (idx, measured)) in enumerate(targets)
        vec = eig0.vectors[order0[idx]]
        for  p = 1:nParam   J[row, p] = (vec' * Symmetric(matrices[p]) * vec) / auPerKayser   end
        residual0[row] = measured - relKayser0[idx]
    end

    # Minimum-norm least-squares correction, by hand through the SVD: a direction with a near-zero singular
    # value is left at delta = 0 rather than amplifying measurement-sized noise into an enormous "fit".
    F = svd(J)
    sMax = isempty(F.S) ? 0.0 : maximum(F.S)
    deltaC = zeros(Float64, nParam)
    for  k = 1:length(F.S)
        F.S[k] <= determinacyThreshold * max(sMax, 1.0)   &&   continue
        deltaC .+= ( (F.U[:,k]' * residual0) / F.S[k] ) .* F.V[:,k]
    end
    cVec = ones(Float64, nParam) .+ deltaC

    # The real fit: re-diagonalize EXACTLY at cVec, so rmsDeviation is not the linearized estimate.
    totalFitted = copy(oneBody)
    for  p = 1:nParam   totalFitted .+= cVec[p] .* matrices[p]   end
    eig1       = Hamiltonian.diagonalizeCiMatrix(totalFitted, LevelSelection())
    order1     = sortperm(eig1.values)
    energies1  = eig1.values[order1]
    relKayser1 = (energies1 .- energies1[1]) ./ auPerKayser

    s = 0.0
    for  (idx, measured) in targets   s = s + (relKayser1[idx] - measured)^2   end
    rmsBest = sqrt(s/length(targets))

    levels = Level[]
    for  (pos, ev) in enumerate(order1)
        vector = zeros( length(basis.csfs) );   ns = 0
        for  r = 1:length(basis.csfs)
            if  LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = eig1.vectors[ev][ns]   end
        end
        push!( levels, Level(sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, pos, energies1[pos], 0., true, basis, vector) )
    end

    radialParameters = [ RadialParameter(parameterClasses[p], cVec[p]) for p = 1:nParam ]

    orthogonalParameters = OrthogonalParameter[]
    for  k = 1:length(F.S)
        composition = [ (parameterClasses[p], F.V[p,k]) for p = 1:nParam ]
        determined  = F.S[k] > determinacyThreshold * max(sMax, 1.0)
        shiftAlong  = determined ? (F.U[:,k]' * residual0) / F.S[k] : 0.0
        push!( orthogonalParameters, OrthogonalParameter("mode $k", shiftAlong, composition, F.S[k]) )
    end

    fitTargets = [ FitTarget(idx, measured, relKayser1[idx]) for (idx, measured) in targets ]

    return( FittedMultipletOrthogonal(Multiplet(string(sym)*"+ (fitted, orthogonal)", levels), radialParameters,
                                        orthogonalParameters, fitTargets, rmsBest, source) )
end


"""
`AdjustedLevels.mergeIntoMultiplet(fm::AdjustedLevels.FittedMultipletOrthogonal, abInitioMultiplet::ManyElectron.Multiplet)`
    ... the same merge step as the single-parameter route's `mergeIntoMultiplet(::FittedMultiplet, ...)`, for a
        `FittedMultipletOrthogonal` result: `fm`'s own symmetry block replaces the matching block of
        `abInitioMultiplet`, every other block is kept unchanged, and the result is re-sorted by energy. A
        mp::ManyElectron.Multiplet is returned.
"""
function  mergeIntoMultiplet(fm::FittedMultipletOrthogonal, abInitioMultiplet::ManyElectron.Multiplet)
    isempty(fm.multiplet.levels)   &&   error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletOrthogonal has no levels.")
    sym = LevelSymmetry(fm.multiplet.levels[1].J, fm.multiplet.levels[1].parity)
    otherLevels = [lev for lev in abInitioMultiplet.levels if LevelSymmetry(lev.J, lev.parity) != sym]

    length(otherLevels) + length(fm.multiplet.levels) != length(abInitioMultiplet.levels)   &&
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletOrthogonal's block has a different number of levels than " *
              "the matching block in abInitioMultiplet -- they must come from the same Basis.")

    allLevels = vcat(otherLevels, fm.multiplet.levels)
    mp = Multiplet(abInitioMultiplet.name * "+fitted", allLevels)

    return( Basics.sortByEnergy(mp) )
end
