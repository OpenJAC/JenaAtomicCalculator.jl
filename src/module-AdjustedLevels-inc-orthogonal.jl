
"""
`struct  AdjustedLevels.RadialParameter`
    ... one fitted scaling parameter for ONE automatically-discovered integral class -- the sibling of
        `AdjustedLevels.ScalingParameter` for the orthogonal route: `fitScaleFactor` always returns exactly
        one, uniform across every rank and every shell pair; `fitOrthogonalParameters` returns one PER
        (operator type, shell pair, rank) class it finds among this block's two-body coefficients -- Cowan's
        own F^k (direct)/G^k (exchange)/R^k (configuration-interaction) separation, not multipole rank alone
        (see `AdjustedLevels.twoBodyParameterClass` and `AdjustedLevels.buildParameterMatrices`).

    + parameterClass   ::String    ... which integral class this factor multiplies, e.g. "F^2(4f,4f)" or
                                        "G^2(2p,3d)" -- see `AdjustedLevels.twoBodyParameterClass`.
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
    + degreesOfFreedom      ::Int64                             ... length(fitTargets) - length(radialParameters).
                                                                      THE OVERFITTING GUARD: zero or negative
                                                                      means every class was fit with no spare
                                                                      data to check it against (zero) or more
                                                                      classes than data at all (negative) --
                                                                      a legitimate fit (the SVD cutoff still
                                                                      protects an individual undetermined
                                                                      direction) but not one whose per-class
                                                                      values should be read as independently
                                                                      confirmed. `fitOrthogonalParameters`
                                                                      prints a warning when this is <= 0.
    + source                ::String                           ... a dated citation for the measured data used.
"""
struct  FittedMultipletOrthogonal
    multiplet               ::ManyElectron.Multiplet
    radialParameters         ::Array{RadialParameter,1}
    orthogonalParameters     ::Array{OrthogonalParameter,1}
    fitTargets               ::Array{FitTarget,1}
    rmsDeviation              ::Float64
    degreesOfFreedom          ::Int64
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
    println(io, "  degrees of freedom: $(fm.degreesOfFreedom)")
    fm.degreesOfFreedom <= 0   &&
        println(io, "  WARNING: at or under the identifiability boundary -- read each class's own determinacy " *
                     "before trusting its value in isolation.")
    println(io, "  multiplet: $(fm.multiplet.name), $(length(fm.multiplet.levels)) level(s)")
end


"""
`AdjustedLevels.twoBodyParameterClass(a::Subshell, b::Subshell, c::Subshell, d::Subshell, nu::Int64)`
    ... to classify ONE two-body coefficient `coeff.V * XL_Coulomb(nu, orbital[a], orbital[b], orbital[c],
        orbital[d], ...)` into Cowan's own RCE integral classes, by the SHELL (n,l) of each involved orbital
        -- not the (possibly split, j-dependent) subshell itself, so that a p_1/2/p_3/2 or f_5/2/f_7/2 pair
        shares ONE parameter exactly as Cowan's nonrelativistic scheme does, rather than fitting a spin-orbit
        splitting of the Slater parameter itself that nobody asks for. Three cases, read directly off which
        subshells coincide (the same convention `SpinAngular.computeCoefficients` itself uses to build a
        direct term as `(sha,shb,sha,shb)` and an exchange term as `(sha,shb,shb,sha)`):

        - DIRECT, `a==c && b==d`: the density-density (F^k) term -- labelled `"F^k(shellA,shellB)"`, shells
          sorted so the label does not depend on which orbital happened to be listed first.
        - EXCHANGE, `a==d && b==c && a!=b`: the antisymmetrization (G^k) term -- `"G^k(shellA,shellB)"`,
          shells sorted likewise; a==b falls through to the DIRECT case above instead (equivalent electrons
          have no separate exchange term, exactly as `SpinAngular`'s own same-shell coefficient has none).
        - OTHERWISE: a genuine configuration-interaction (R^k) term connecting two different orbital
          occupations -- `"R^k(shellA,shellB;shellC,shellD)"`, canonicalized over the integral's own 4-index
          symmetry `R^k(ab,cd) = R^k(ba,dc) = R^k(cd,ab) = R^k(dc,ba)` by taking the lexicographically
          smallest of the four equivalent labellings, so the SAME physical integral is never split into
          several distinct classes merely because of which CSF pair it was read off from.

        A value::String is returned.
"""
function  twoBodyParameterClass(a::Subshell, b::Subshell, c::Subshell, d::Subshell, nu::Int64)
    shellOf(sh::Subshell) = string(sh.n) * Basics.shellNotation( Basics.subshell_l(sh) )

    if  a == c  &&  b == d
        sa, sb = shellOf(a), shellOf(b)
        lo, hi = sa <= sb ? (sa, sb) : (sb, sa)
        return( "F^$(nu)($(lo),$(hi))" )
    elseif  a == d  &&  b == c  &&  a != b
        sa, sb = shellOf(a), shellOf(b)
        lo, hi = sa <= sb ? (sa, sb) : (sb, sa)
        return( "G^$(nu)($(lo),$(hi))" )
    else
        tup1 = (shellOf(a), shellOf(b), shellOf(c), shellOf(d))
        tup2 = (shellOf(b), shellOf(a), shellOf(d), shellOf(c))
        tup3 = (shellOf(c), shellOf(d), shellOf(a), shellOf(b))
        tup4 = (shellOf(d), shellOf(c), shellOf(b), shellOf(a))
        best = minimum((tup1, tup2, tup3, tup4))
        return( "R^$(nu)($(best[1]),$(best[2]);$(best[3]),$(best[4]))" )
    end
end


"""
`AdjustedLevels.buildParameterMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)`
    ... the orthogonal route's generalization of `AdjustedLevels.buildOneAndTwoBodyMatrices`: instead of
        summing every two-body (Coulomb) contribution into ONE lumped matrix, each contribution is kept in
        the matrix belonging to its own (operator type, shell pair, rank) class, via
        `AdjustedLevels.twoBodyParameterClass` -- Cowan's own F^k/G^k/R^k separation, which rank alone cannot
        give: two different shell pairs can share the same rank (a p-d pair's F^2 and G^2 both exist at
        k=2), and a correlation model with more than one active shell pair would otherwise have them fitted
        as if they were one physical quantity. A CLOSED subshell is never separately parametrized (Cowan's
        RCE fits the open shell(s) only, see this function's own comment on why), and a class whose matrix
        turns out to be a multiple of the identity on this block -- the rank-0 monopole within one fixed-
        occupation open shell, most often -- is folded into `oneBody` too, since no splitting could ever
        determine it. Even a single OPEN shell with l>=1 still separates into several classes, not one: its
        own relativistic j-splitting (e.g. 4f_5/2 against 4f_7/2) already gives distinct F^k, G^k and R^k
        integrals at a given rank, which a single rank-only bucket could not have told apart either.
        A (idxCsf, oneBody, parameterClasses, matrices) tuple is returned: idxCsf::Array{Int64,1} as in
        `buildOneAndTwoBodyMatrices`; oneBody::Array{Float64,2} the shared one-body matrix; parameterClasses
        ::Array{String,1} a label per discovered class; matrices::Array{Array{Float64,2},1} the matching
        upper-triangular two-body matrices, one per class, same convention as `buildOneAndTwoBodyMatrices`.
"""
function  buildParameterMatrices(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid)
    idxCsf = [idx for idx = 1:length(basis.csfs) if basis.csfs[idx].J == sym.J && basis.csfs[idx].parity == sym.parity]
    n = length(idxCsf)
    potential = Nuclear.nuclearPotential(nm, grid)
    oneBody = zeros(Float64, n, n)
    twoBodyByClass = Dict{String,Array{Float64,2}}()

    # CLOSED subshells (full occupation in EVERY CSF of this block) are NEVER separately parametrized here,
    # matching Cowan's own RCE exactly: it fits Slater parameters among the OPEN shell(s) only and treats a
    # closed core as part of the fixed, frozen background -- never as a source of its own adjustable integral,
    # whether that closed shell is interacting with itself, with another closed shell, or with an open one.
    # (MEASURED, not merely argued: a core-open EXCHANGE term can still vary CSF to CSF within one fixed open
    # shell's own internal coupling -- the closed shell does not simply average it away -- so this is a scope
    # choice tracking standard RCE practice, not a claim that every excluded term is itself numerically inert.)
    # Left unfiltered, an explicit closed core (as any `[Xe] ...`-style reference carries) turns into hundreds
    # of such classes -- academic clutter, not capability.
    closedSubshell = Dict{Subshell,Bool}()
    for  (i, sh) in enumerate(basis.subshells)
        capacity = Basics.subshell_2j(sh) + 1
        closedSubshell[sh] = all( basis.csfs[idxCsf[r]].occupation[i] == capacity  for  r = 1:n )
    end
    allOpen(a,b,c,d) = !closedSubshell[a] && !closedSubshell[b] && !closedSubshell[c] && !closedSubshell[d]

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

            for  coeff in waG2
                contribution = coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                          basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                if  allOpen(coeff.a, coeff.b, coeff.c, coeff.d)
                    class = twoBodyParameterClass(coeff.a, coeff.b, coeff.c, coeff.d, coeff.nu)
                    if  !haskey(twoBodyByClass, class)   twoBodyByClass[class] = zeros(Float64, n, n)   end
                    twoBodyByClass[class][r,s] += contribution
                else
                    me1 = me1 + contribution
                end
            end
            oneBody[r,s] = me1
        end
    end

    # AMONG THE SURVIVING, OPEN-SHELL-ONLY classes, one further, numerical (not scope-based) simplification:
    # a class whose matrix is a multiple of the IDENTITY on this block -- same value on every diagonal, zero
    # off it -- shifts every level by the SAME amount and so cancels EXACTLY from any splitting regardless of
    # its scale factor. The rank-0 (monopole) direct term within a single fixed-occupation open shell is
    # exactly this case (F^0 is the coupling-independent average Coulomb repulsion). Folded into `oneBody` at
    # its full, UNSCALED value, same reasoning as the closed-shell exclusion above, one level down.
    for  (class, M) in collect(twoBodyByClass)
        diagVals  = [ M[i,i] for i = 1:n ]
        scale     = max(1.0, maximum(abs.(diagVals)))
        offDiag   = n > 1 ? maximum( abs(M[i,j])  for i = 1:n  for j = i+1:n ) : 0.0
        isInert   = (maximum(diagVals) - minimum(diagVals)) < 1.0e-10*scale   &&   offDiag < 1.0e-10*scale
        if  isInert
            oneBody .+= M
            delete!(twoBodyByClass, class)
        end
    end

    classes          = sort(collect(keys(twoBodyByClass)))
    parameterClasses = classes
    matrices         = [ twoBodyByClass[cl] for cl in classes ]

    return( (idxCsf, oneBody, parameterClasses, matrices) )
end


"""
`AdjustedLevels.fitOrthogonalParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                          targets::Array{Tuple{Int64,Float64},1}, source::String;
                                          determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                          tolerance::Float64=1.0e-6)`
    ... the orthogonal (Uylings/Raassen-style) route: discovers the integral classes present in `sym`'s block
        via `buildParameterMatrices`, fits one scaling parameter PER class against `targets` (same convention
        as `fitScaleFactor`'s own `targets`), and reports the fit in terms of the ORTHOGONAL combinations the
        targets can actually determine. THE FIT IS A GAUSS-NEWTON ITERATION, not a single linearized step
        (that single-step version was tried first and MEASURED to leave a real, growing residual even in an
        idealised recovery test with exactly as many targets as free parameters -- one linear step is only the
        leading-order correction to a genuinely nonlinear eigenvalue problem). At each iteration: diagonalize
        EXACTLY at the current parameters, build the sensitivity matrix J[target,class] = d(energy)/d(parameter)
        from the CURRENT eigenvectors (no extra diagonalization beyond the one already needed), take the
        MINIMUM-NORM correction via `LinearAlgebra.svd` (a direction the targets cannot determine is left
        untouched rather than amplifying noise), and re-center there for the next round. Iteration stops when
        the rms stops improving by more than `tolerance`, or after `maxIterations` -- AND THE BEST ITERATE SEEN
        is what is returned, never just the last one: an undamped Gauss-Newton step is not guaranteed to
        improve the fit (measured directly: without this guard, one real case got WORSE after iterating, 175
        to 323 cm^-1), so every candidate point is compared against the best found so far. MEASURED RANGE OF
        VALIDITY: on an idealised recovery test (synthetic targets generated from this route's own model, so a
        perfect fit exists in principle), this converges to the true parameters for corrections up to roughly
        30% per class; beyond that the iteration can stall in a poorly-conditioned region well short of the
        solution (a slow, near-linear crawl, not an oscillation) rather than diverging outright. A correction
        this large is exactly the sign that the underlying ab initio calculation needs more correlation
        (more configurations), not a better fit -- see this module's own docstring. `value` on each returned
        `OrthogonalParameter` is the TOTAL accumulated shift along that direction, and `determinacy` is read at
        the CONVERGED point, not the starting one -- so a mode whose shift is reported but whose determinacy is
        still near zero signals the fit has not actually converged along it.
        THE ITERATION IS WARM-STARTED from the single-factor (uniform) optimum over these same matrices, not
        from all-ones: "every rank equal" is itself a point inside this route's own parameter space, so this
        route's rms can never be worse than `fitScaleFactor`'s on the same targets -- it WAS, measured, before
        this warm start was added, since a plain Gauss-Newton iteration can fail to reach even that nested
        point within `maxIterations` when starting far from it (confirmed on a real case: Pr IV 4f^2, J=2,
        where the orthogonal route returned a worse rms than the single-factor route on identical targets --
        mathematically impossible for the true optimum of a strictly richer, nested model).
        A fm::AdjustedLevels.FittedMultipletOrthogonal is returned.
"""
function  fitOrthogonalParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                    targets::Array{Tuple{Int64,Float64},1}, source::String;
                                    determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                    tolerance::Float64=1.0e-6)
    isempty(targets)   &&   error("AdjustedLevels.fitOrthogonalParameters(): need at least one target level to fit against.")
    idxCsf, oneBody, parameterClasses, matrices = buildParameterMatrices(sym, basis, nm, grid)
    nParam = length(matrices)
    nParam == 0   &&   error("AdjustedLevels.fitOrthogonalParameters(): no two-body (Coulomb) contributions were found for this block.")
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    nTargets    = length(targets)

    # WARM START from the UNIFORM (single-factor-equivalent) optimum over these SAME matrices -- cheap, since
    # it costs only re-diagonalizations, no rebuild. "every rank equal" is a point INSIDE this route's own
    # parameter space, so without this the richer model could (and was measured to) return a WORSE rms than
    # the simpler, nested one it contains as a special case -- a mathematical impossibility for the true
    # optimum, and a sign the iteration was starting too far away to reach it within `maxIterations`.
    function  rmsAtUniform(c::Float64)
        total = copy(oneBody);   for  mtx in matrices   total .+= c .* mtx   end
        eig = Hamiltonian.diagonalizeCiMatrix(total, LevelSelection())
        rel = (sort(eig.values) .- sort(eig.values)[1]) ./ auPerKayser
        s = 0.0;   for  (idx, measured) in targets   s = s + (rel[idx] - measured)^2   end
        return( sqrt(s/length(targets)) )
    end
    cUniform, _ = goldenSectionMinimize(rmsAtUniform, 0.5, 1.5)

    cVec    = fill(cUniform, nParam)
    rmsPrev = Inf
    bestCVec = copy(cVec);   bestRms = Inf
    bestEigC = nothing;   bestOrderC = nothing;   bestRelKayserC = Float64[]

    for  iteration = 1:maxIterations
        totalC = copy(oneBody);   for  p = 1:nParam   totalC .+= cVec[p] .* matrices[p]   end
        eigC   = Hamiltonian.diagonalizeCiMatrix(totalC, LevelSelection())
        orderC = sortperm(eigC.values)
        energiesC  = eigC.values[orderC]
        relKayserC = (energiesC .- energiesC[1]) ./ auPerKayser

        J = zeros(Float64, nTargets, nParam);   residual = zeros(Float64, nTargets)
        for  (row, (idx, measured)) in enumerate(targets)
            vec = eigC.vectors[orderC[idx]]
            for  p = 1:nParam   J[row, p] = (vec' * Symmetric(matrices[p]) * vec) / auPerKayser   end
            residual[row] = measured - relKayserC[idx]
        end
        rmsNow = sqrt(sum(residual.^2)/nTargets)

        # TRACK THE BEST ITERATE SEEN, and return that one -- Gauss-Newton without damping is NOT guaranteed
        # to improve every step (confirmed by measurement: an undamped run made a real case worse, 175 -> 323
        # cm^-1, before this guard was added), so trusting the LAST point blindly is unsafe.
        if  rmsNow < bestRms
            bestRms = rmsNow;   bestCVec = copy(cVec)
            bestEigC = eigC;    bestOrderC = orderC;   bestRelKayserC = relKayserC
        end

        ( rmsNow < tolerance   ||   ( iteration > 1   &&   rmsNow >= rmsPrev - tolerance * max(rmsPrev, 1.0) ) )   &&   break
        rmsPrev = rmsNow

        F    = svd(J)
        sMax = isempty(F.S) ? 0.0 : maximum(F.S)
        deltaC = zeros(Float64, nParam)
        for  k = 1:length(F.S)
            F.S[k] <= determinacyThreshold * max(sMax, 1.0)   &&   continue
            deltaC .+= ( (F.U[:,k]' * residual) / F.S[k] ) .* F.V[:,k]
        end
        cVec = cVec .+ deltaC
    end

    cVec   = bestCVec;   eigC = bestEigC;   orderC = bestOrderC;   relKayserC = bestRelKayserC
    rmsBest = bestRms

    levels = Level[]
    for  (pos, ev) in enumerate(orderC)
        vector = zeros( length(basis.csfs) );   ns = 0
        for  r = 1:length(basis.csfs)
            if  LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = eigC.vectors[ev][ns]   end
        end
        push!( levels, Level(sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, pos, eigC.values[orderC[pos]], 0., true, basis, vector) )
    end

    radialParameters = [ RadialParameter(parameterClasses[p], cVec[p]) for p = 1:nParam ]

    # Determinacy and the orthogonal breakdown are read AT THE CONVERGED point, from a fresh sensitivity matrix --
    # not carried over from the first iteration's, which described a point the fit has since moved away from.
    Jfinal = zeros(Float64, nTargets, nParam);   residualFinal = zeros(Float64, nTargets)
    for  (row, (idx, measured)) in enumerate(targets)
        vec = eigC.vectors[orderC[idx]]
        for  p = 1:nParam   Jfinal[row, p] = (vec' * Symmetric(matrices[p]) * vec) / auPerKayser   end
        residualFinal[row] = measured - relKayserC[idx]
    end
    Ffinal    = svd(Jfinal)
    sMaxFinal = isempty(Ffinal.S) ? 0.0 : maximum(Ffinal.S)
    totalShift = cVec .- ones(nParam)

    orthogonalParameters = OrthogonalParameter[]
    for  k = 1:length(Ffinal.S)
        composition = [ (parameterClasses[p], Ffinal.V[p,k]) for p = 1:nParam ]
        push!( orthogonalParameters, OrthogonalParameter("mode $k", dot(Ffinal.V[:,k], totalShift), composition, Ffinal.S[k]) )
    end

    fitTargets = [ FitTarget(idx, measured, relKayserC[idx]) for (idx, measured) in targets ]

    return( FittedMultipletOrthogonal(Multiplet(string(sym)*"+ (fitted, orthogonal)", levels), radialParameters,
                                        orthogonalParameters, fitTargets, rmsBest, nTargets - nParam, source) )
end


"""
`struct  AdjustedLevels.FittedBlock`
    ... one symmetry block's share of a `AdjustedLevels.FittedMultipletOrthogonalMultiBlock` fit: the radial
        parameters are shared GLOBALLY across every block of such a fit, but each block keeps its own CSFs,
        its own re-diagonalized levels and its own measured targets.

    + sym          ::LevelSymmetry              ... which J^P block this is.
    + multiplet    ::ManyElectron.Multiplet      ... this block's re-diagonalized, fitted levels.
    + fitTargets   ::Array{FitTarget,1}          ... the measured levels THIS block was matched against.
"""
struct  FittedBlock
    sym          ::LevelSymmetry
    multiplet    ::ManyElectron.Multiplet
    fitTargets   ::Array{FitTarget,1}
end


"""
`struct  AdjustedLevels.FittedMultipletOrthogonalMultiBlock`
    ... the result of `AdjustedLevels.fitOrthogonalParametersMultiBlock`: ONE shared set of radial parameters,
        fit SIMULTANEOUSLY against measured levels drawn from SEVERAL symmetry blocks of the SAME basis --
        the capability this module's own docstring names as missing from both single-block routes. An
        `F^k(4f,4f)`-type class means the SAME physical radial integral in every block it appears in, so
        pooling blocks turns a borderline, exactly-determined (or worse) per-block fit into a comfortably
        over-determined one, at no extra configuration-space cost -- see `degreesOfFreedom`.

    + blocks                ::Array{FittedBlock,1}            ... one entry per symmetry block, in the order
                                                                     given to `fitOrthogonalParametersMultiBlock`.
    + radialParameters       ::Array{RadialParameter,1}         ... the fitted factors, GLOBAL across every
                                                                       block -- the same class in two blocks is
                                                                       the SAME parameter, not two independent
                                                                       ones.
    + orthogonalParameters   ::Array{OrthogonalParameter,1}      ... the orthogonal breakdown of the GLOBAL
                                                                        parameter space, read at the converged
                                                                        point, over EVERY target in EVERY
                                                                        block together -- the sibling of
                                                                        `FittedMultipletOrthogonal`'s own field.
    + rmsDeviation            ::Float64                           ... root-mean-square deviation over every
                                                                         target in every block together, in
                                                                         cm^-1.
    + degreesOfFreedom        ::Int64                             ... total targets (all blocks) minus the
                                                                         number of GLOBAL parameters -- THE
                                                                         OVERFITTING GUARD, see
                                                                         `AdjustedLevels.FittedMultipletOrthogonal`'s
                                                                         own field of the same name.
    + source                  ::String                            ... a dated citation for the measured data.
"""
struct  FittedMultipletOrthogonalMultiBlock
    blocks                ::Array{FittedBlock,1}
    radialParameters       ::Array{RadialParameter,1}
    orthogonalParameters    ::Array{OrthogonalParameter,1}
    rmsDeviation             ::Float64
    degreesOfFreedom         ::Int64
    source                   ::String
end


"""
`Base.show(io::IO, fm::AdjustedLevels.FittedMultipletOrthogonalMultiBlock)`
    ... ALWAYS prints the "fitted, not ab initio" disclosure first, the GLOBAL parameters, the overfitting-
        guard warning when `degreesOfFreedom <= 0`, and then each block's own share of the fit. Nothing is
        returned.
"""
function  Base.show(io::IO, fm::FittedMultipletOrthogonalMultiBlock)
    println(io, "FittedMultipletOrthogonalMultiBlock -- FITTED, NOT AB INITIO.")
    println(io, "  source: $(fm.source)")
    for  rp in fm.radialParameters
        println(io, "  radial parameter [$(rp.parameterClass)] = $(rp.value)  (shared across every block below)")
    end
    for  op in fm.orthogonalParameters
        flag = op.determinacy < 1.0e-6 ? "  (UNDETERMINED by these targets)" : ""
        println(io, "  orthogonal $(op.label): shift = $(op.value), determinacy = $(op.determinacy)$flag")
    end
    nTargetsTotal = sum( length(b.fitTargets) for b in fm.blocks )
    println(io, "  fit against $(nTargetsTotal) measured level(s) across $(length(fm.blocks)) block(s), " *
                 "rms deviation = $(fm.rmsDeviation) cm^-1")
    println(io, "  degrees of freedom: $(fm.degreesOfFreedom)")
    fm.degreesOfFreedom <= 0   &&
        println(io, "  WARNING: at or under the identifiability boundary -- read each class's own determinacy " *
                     "before trusting its value in isolation.")
    for  b in fm.blocks
        println(io, "  block $(b.sym): $(length(b.fitTargets)) target(s), $(length(b.multiplet.levels)) level(s)")
    end
end


"""
`AdjustedLevels.fitOrthogonalParametersMultiBlock(blocks::Array{Tuple{LevelSymmetry,Array{Tuple{Int64,Float64},1}},1},
                                                    basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, source::String;
                                                    determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                                    tolerance::Float64=1.0e-6)`
    ... the multi-block generalization of `AdjustedLevels.fitOrthogonalParameters`: ONE shared parameter per
        (type, shell pair, rank) class -- see `AdjustedLevels.twoBodyParameterClass` -- fit SIMULTANEOUSLY
        against targets drawn from SEVERAL symmetry blocks of the SAME `basis`, e.g.
        `[(sym2, targets2), (sym4, targets4)]` for two different J^P blocks of one open-shell reference.
        THIS IS WHAT MAKES CLASS-SEPARATED FITTING ACTUALLY WORTH HAVING, not merely correct: a class shared
        by several blocks is the SAME physical radial integral in every one of them, so pooling the blocks'
        targets turns an exactly-determined (or worse) per-block fit into a comfortably over-determined one,
        without adding a single CSF -- see `AdjustedLevels.FittedMultipletOrthogonalMultiBlock`'s own
        `degreesOfFreedom` field. The iteration is otherwise the SAME Gauss-Newton + SVD scheme
        `fitOrthogonalParameters` uses (warm-started the same way, from the single uniform factor minimizing
        the COMBINED rms over every block), just with each block's own Jacobian rows and residuals stacked
        together before the one, shared SVD solve each iteration -- a class absent from a given block simply
        contributes a zero column there. A fm::AdjustedLevels.FittedMultipletOrthogonalMultiBlock is returned.
"""
function  fitOrthogonalParametersMultiBlock(blocks::Array{Tuple{LevelSymmetry,Array{Tuple{Int64,Float64},1}},1},
                                              basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, source::String;
                                              determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                              tolerance::Float64=1.0e-6)
    isempty(blocks)   &&   error("AdjustedLevels.fitOrthogonalParametersMultiBlock(): need at least one block.")
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    nBlocks = length(blocks)

    idxCsfB   = Vector{Array{Int64,1}}(undef, nBlocks)
    oneBodyB  = Vector{Array{Float64,2}}(undef, nBlocks)
    classesB  = Vector{Array{String,1}}(undef, nBlocks)
    matricesB = Vector{Array{Array{Float64,2},1}}(undef, nBlocks)
    for  ib = 1:nBlocks
        sym, _ = blocks[ib]
        idxCsfB[ib], oneBodyB[ib], classesB[ib], matricesB[ib] = buildParameterMatrices(sym, basis, nm, grid)
    end

    globalClasses = sort(unique(vcat(classesB...)))
    nParam = length(globalClasses)
    nParam == 0   &&
        error("AdjustedLevels.fitOrthogonalParametersMultiBlock(): no two-body (Coulomb) contributions were found in any block.")
    classIndex = Dict{String,Int64}( cl => k for (k,cl) in enumerate(globalClasses) )

    # Each block's own two-body matrices, reindexed onto the GLOBAL class list; a class absent from a block
    # simply has no entry here (equivalent to a zero column for that block in the Jacobian below).
    globalMatricesB = Vector{Dict{Int64,Array{Float64,2}}}(undef, nBlocks)
    for  ib = 1:nBlocks
        d = Dict{Int64,Array{Float64,2}}()
        for  (loc, cl) in enumerate(classesB[ib])   d[classIndex[cl]] = matricesB[ib][loc]   end
        globalMatricesB[ib] = d
    end

    nTargetsTotal = sum( length(blocks[ib][2]) for ib = 1:nBlocks )

    # Warm start: ONE uniform factor minimizing the COMBINED rms over every block's own targets together --
    # the multi-block counterpart of fitOrthogonalParameters' own single-block warm start, for the same reason.
    function  rmsAtUniformAll(c::Float64)
        s = 0.0
        for  ib = 1:nBlocks
            total = copy(oneBodyB[ib]);   for  mtx in matricesB[ib]   total .+= c .* mtx   end
            eig = Hamiltonian.diagonalizeCiMatrix(total, LevelSelection())
            rel = (sort(eig.values) .- sort(eig.values)[1]) ./ auPerKayser
            for  (idx, measured) in blocks[ib][2]   s = s + (rel[idx] - measured)^2   end
        end
        return( sqrt(s/nTargetsTotal) )
    end
    cUniform, _ = goldenSectionMinimize(rmsAtUniformAll, 0.5, 1.5)

    cVec = fill(cUniform, nParam)
    rmsPrev = Inf
    bestCVec = copy(cVec);   bestRms = Inf
    bestEigB = Vector{Any}(undef, nBlocks);   bestOrderB = Vector{Any}(undef, nBlocks)
    bestRelB = Vector{Array{Float64,1}}(undef, nBlocks)

    for  iteration = 1:maxIterations
        eigB   = Vector{Any}(undef, nBlocks);   orderB = Vector{Any}(undef, nBlocks)
        relB   = Vector{Array{Float64,1}}(undef, nBlocks)
        J        = zeros(Float64, nTargetsTotal, nParam)
        residual = zeros(Float64, nTargetsTotal)
        row = 0

        for  ib = 1:nBlocks
            total = copy(oneBodyB[ib])
            for  (gIdx, mtx) in globalMatricesB[ib]   total .+= cVec[gIdx] .* mtx   end
            eig   = Hamiltonian.diagonalizeCiMatrix(total, LevelSelection())
            order = sortperm(eig.values)
            relKayser = (eig.values[order] .- eig.values[order][1]) ./ auPerKayser
            eigB[ib] = eig;   orderB[ib] = order;   relB[ib] = relKayser

            for  (idx, measured) in blocks[ib][2]
                row += 1
                vec = eig.vectors[order[idx]]
                for  (gIdx, mtx) in globalMatricesB[ib]   J[row, gIdx] = (vec' * Symmetric(mtx) * vec) / auPerKayser   end
                residual[row] = measured - relKayser[idx]
            end
        end
        rmsNow = sqrt(sum(residual.^2)/nTargetsTotal)

        if  rmsNow < bestRms
            bestRms = rmsNow;   bestCVec = copy(cVec)
            bestEigB = eigB;    bestOrderB = orderB;   bestRelB = relB
        end

        ( rmsNow < tolerance   ||   ( iteration > 1   &&   rmsNow >= rmsPrev - tolerance * max(rmsPrev, 1.0) ) )   &&   break
        rmsPrev = rmsNow

        F    = svd(J)
        sMax = isempty(F.S) ? 0.0 : maximum(F.S)
        deltaC = zeros(Float64, nParam)
        for  k = 1:length(F.S)
            F.S[k] <= determinacyThreshold * max(sMax, 1.0)   &&   continue
            deltaC .+= ( (F.U[:,k]' * residual) / F.S[k] ) .* F.V[:,k]
        end
        cVec = cVec .+ deltaC
    end

    cVec = bestCVec;   eigB = bestEigB;   orderB = bestOrderB;   relB = bestRelB
    rmsBest = bestRms
    radialParameters = [ RadialParameter(globalClasses[p], cVec[p]) for p = 1:nParam ]

    # Determinacy and the orthogonal breakdown, read AT THE CONVERGED (best) point -- same discipline as the
    # single-block route: a mode's determinacy describes where the fit ENDED, not an earlier iterate.
    Jfinal = zeros(Float64, nTargetsTotal, nParam);   residualFinal = zeros(Float64, nTargetsTotal);   row = 0
    for  ib = 1:nBlocks
        for  (idx, measured) in blocks[ib][2]
            row += 1
            vec = eigB[ib].vectors[orderB[ib][idx]]
            for  (gIdx, mtx) in globalMatricesB[ib]   Jfinal[row, gIdx] = (vec' * Symmetric(mtx) * vec) / auPerKayser   end
            residualFinal[row] = measured - relB[ib][idx]
        end
    end
    Ffinal = svd(Jfinal)
    totalShift = cVec .- ones(nParam)
    orthogonalParameters = OrthogonalParameter[]
    for  k = 1:length(Ffinal.S)
        composition = [ (globalClasses[p], Ffinal.V[p,k]) for p = 1:nParam ]
        push!( orthogonalParameters, OrthogonalParameter("mode $k", dot(Ffinal.V[:,k], totalShift), composition, Ffinal.S[k]) )
    end

    fittedBlocks = FittedBlock[]
    for  ib = 1:nBlocks
        sym, targets = blocks[ib]
        levels = Level[]
        for  (pos, ev) in enumerate(orderB[ib])
            vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs)
                if  LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = eigB[ib].vectors[ev][ns]   end
            end
            push!( levels, Level(sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, pos, eigB[ib].values[orderB[ib][pos]], 0., true, basis, vector) )
        end
        fitTargets = [ FitTarget(idx, measured, relB[ib][idx]) for (idx, measured) in targets ]
        push!( fittedBlocks, FittedBlock(sym, Multiplet(string(sym)*"+ (fitted, multi-block)", levels), fitTargets) )
    end

    return( FittedMultipletOrthogonalMultiBlock(fittedBlocks, radialParameters, orthogonalParameters,
                                                  rmsBest, nTargetsTotal - nParam, source) )
end


"""
`AdjustedLevels.mergeIntoMultiplet(fm::AdjustedLevels.FittedMultipletOrthogonalMultiBlock, abInitioMultiplet::ManyElectron.Multiplet)`
    ... the multi-block generalization of the other `mergeIntoMultiplet` methods: EVERY block carried by `fm`
        replaces its own matching block of `abInitioMultiplet` in one call, every OTHER symmetry block is
        kept unchanged, and the result is re-sorted by energy. A mp::ManyElectron.Multiplet is returned.
"""
function  mergeIntoMultiplet(fm::FittedMultipletOrthogonalMultiBlock, abInitioMultiplet::ManyElectron.Multiplet)
    isempty(fm.blocks)   &&   error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletOrthogonalMultiBlock has no blocks.")
    syms = [ b.sym for b in fm.blocks ]
    length(unique(syms)) == length(syms)   ||
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletOrthogonalMultiBlock carries the same symmetry block twice.")

    otherLevels  = [ lev for lev in abInitioMultiplet.levels if !(LevelSymmetry(lev.J, lev.parity) in syms) ]
    fittedLevels = vcat( [ b.multiplet.levels for b in fm.blocks ]... )

    length(otherLevels) + length(fittedLevels) != length(abInitioMultiplet.levels)   &&
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletOrthogonalMultiBlock's blocks do not account for " *
              "exactly the levels of the matching blocks in abInitioMultiplet -- they must come from the same Basis.")

    mp = Multiplet(abInitioMultiplet.name * "+fitted", vcat(otherLevels, fittedLevels))

    return( Basics.sortByEnergy(mp) )
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


### GROUPED (Slater-vs-CI) route: the physically-motivated middle ground between ONE uniform factor and
### one factor PER Cowan class, added 23-Oct-2026 on the maintainer's own request while scoping priority
### item 91. Once a correlating configuration is added, `fitOrthogonalParameters`'s full per-class split
### routinely needs more classes than any realistic set of measured levels can determine (measured: 20
### classes against 2-13 targets for a correlated open d/f-shell block); `fitScaleFactor`'s single factor,
### at the other extreme, forces ONE number to carry both the Slater-integral size (which standard RCE
### practice corrects by roughly 0.80-0.85) and the configuration-interaction coupling strength (which
### published Cowan-code analyses instead scale separately, typically by a SMALLER factor, around 0.70) --
### two different physical quantities with different, well-documented empirical corrections, forced to
### share one number. `twoBodyParameterGroup` sorts every class `AdjustedLevels.twoBodyParameterClass`
### finds into exactly those two Cowan families: F^k/G^k (Slater, within or between the ORIGINAL orbitals)
### against R^k (configuration-interaction, connecting different configurations) -- never a class count
### picked without a physical reason, exactly the discipline the maintainer asked this route to keep.

"""
`AdjustedLevels.twoBodyParameterGroup(className::String)`
    ... to sort one `AdjustedLevels.twoBodyParameterClass` label into Cowan's own TWO scaling families --
        `"Slater (F^k/G^k)"` for a direct or exchange integral (the class name starts with `"F^"` or
        `"G^"`), `"Configuration-interaction (R^k)"` for a class connecting different orbital occupations
        (the class name starts with `"R^"`). This is the grouping `AdjustedLevels.fitGroupedParameters` and
        `AdjustedLevels.fitGroupedParametersMultiBlock` use, matching published Cowan-code practice of
        scaling the two families by different empirical factors rather than one shared number. A
        value::String is returned.
"""
function  twoBodyParameterGroup(className::String)

    return(  startswith(className, "R^")  ?  "Configuration-interaction (R^k)"  :  "Slater (F^k/G^k)"  )
end


"""
`AdjustedLevels.groupParameterMatrices(classes::Array{String,1}, matrices::Array{Array{Float64,2},1})`
    ... to fold `AdjustedLevels.buildParameterMatrices`'s per-class two-body matrices into the TWO Cowan
        families `AdjustedLevels.twoBodyParameterGroup` sorts every class into, by plain matrix ADDITION --
        since every class in one group shares the SAME intended scale factor, the matrix a single shared
        factor would multiply is exactly the sum of the matrices the individual classes would have used. A
        (groupLabels, groupMatrices, memberLists)::Tuple{Array{String,1},Array{Array{Float64,2},1},
        Array{Array{String,1},1}} is returned: `groupLabels` the (at most two) group names found, sorted;
        `groupMatrices` the matching summed two-body matrices, same upper-triangular convention as
        `buildParameterMatrices`; `memberLists` which raw classes were folded into each group, so a reader
        of the fitted result can see exactly what a group label stands for rather than trusting the name
        alone.
"""
function  groupParameterMatrices(classes::Array{String,1}, matrices::Array{Array{Float64,2},1})
    summed  = Dict{String,Array{Float64,2}}()
    members = Dict{String,Array{String,1}}()
    for  (cl, mtx) in zip(classes, matrices)
        g = twoBodyParameterGroup(cl)
        if  !haskey(summed, g)   summed[g] = zeros(size(mtx));   members[g] = String[]   end
        summed[g]  .+= mtx
        push!(members[g], cl)
    end
    groupLabels   = sort(collect(keys(summed)))
    groupMatrices = [ summed[g]         for g in groupLabels ]
    memberLists   = [ sort(members[g])  for g in groupLabels ]

    return( (groupLabels, groupMatrices, memberLists) )
end


"""
`struct  AdjustedLevels.GroupedParameter`
    ... one fitted scaling parameter for one of the TWO Cowan families `AdjustedLevels.twoBodyParameterGroup`
        sorts every integral class into -- the sibling of `AdjustedLevels.RadialParameter` for the grouped
        route, carrying the member classes too so a group label is never read as a bare, unexplained name.

    + groupLabel      ::String          ... which family this factor multiplies, `"Slater (F^k/G^k)"` or
                                              `"Configuration-interaction (R^k)"` -- see
                                              `AdjustedLevels.twoBodyParameterGroup`.
    + value            ::Float64        ... the fitted factor itself; NOT 1.0 by construction unless this
                                              group turned out to be exactly undetermined by the given targets.
    + memberClasses    ::Array{String,1} ... the raw `AdjustedLevels.twoBodyParameterClass` labels folded
                                               into this one group, e.g. `["F^0(3d,3d)", "F^2(3d,3d)", ...]`.
"""
struct  GroupedParameter
    groupLabel         ::String
    value              ::Float64
    memberClasses       ::Array{String,1}
end


"""
`struct  AdjustedLevels.FittedMultipletGrouped`
    ... the result of `AdjustedLevels.fitGroupedParameters`. The sibling of `AdjustedLevels.FittedMultipletOrthogonal`
        for the grouped (Slater-vs-CI) route -- a DIFFERENT type, deliberately, matching that module
        convention, so a reader can tell at a glance which of the three class-aware routes produced a given
        result.

    + multiplet             ::ManyElectron.Multiplet         ... the single-symmetry-block multiplet from the
                                                                    FITTED, re-diagonalized Hamiltonian.
    + groupedParameters      ::Array{GroupedParameter,1}       ... the fitted factors, one PER FAMILY (at
                                                                      most two), translated back to the
                                                                      familiar language.
    + orthogonalParameters   ::Array{OrthogonalParameter,1}     ... what was actually fit -- the orthogonal
                                                                       directions over the (at most two)
                                                                       grouped parameters, each with its own
                                                                       determinacy; the sibling of
                                                                       `FittedMultipletOrthogonal`'s own field.
    + fitTargets             ::Array{FitTarget,1}                ... the measured levels this fit was matched
                                                                        against.
    + rmsDeviation           ::Float64                           ... root-mean-square deviation, fitted vs
                                                                        measured, over the fit targets, in
                                                                        cm^-1, from the EXACT re-diagonalized
                                                                        spectrum (not the linearized estimate).
    + degreesOfFreedom       ::Int64                             ... length(fitTargets) - length(groupedParameters);
                                                                        see `FittedMultipletOrthogonal`'s own
                                                                        field of the same name for what this
                                                                        guards against.
    + source                 ::String                           ... a dated citation for the measured data used.
"""
struct  FittedMultipletGrouped
    multiplet               ::ManyElectron.Multiplet
    groupedParameters         ::Array{GroupedParameter,1}
    orthogonalParameters       ::Array{OrthogonalParameter,1}
    fitTargets                  ::Array{FitTarget,1}
    rmsDeviation                 ::Float64
    degreesOfFreedom              ::Int64
    source                        ::String
end


"""
`Base.show(io::IO, fm::AdjustedLevels.FittedMultipletGrouped)`
    ... ALWAYS prints the "fitted, not ab initio" disclosure first, exactly as the other two routes'
        `Base.show` methods do, plus which raw classes each group actually contains. Nothing is returned.
"""
function  Base.show(io::IO, fm::FittedMultipletGrouped)
    println(io, "FittedMultipletGrouped -- FITTED, NOT AB INITIO.")
    println(io, "  source: $(fm.source)")
    for  gp in fm.groupedParameters
        println(io, "  grouped parameter [$(gp.groupLabel)] = $(gp.value)   " *
                     "($(length(gp.memberClasses)) class(es): $(join(gp.memberClasses, ", ")))")
    end
    for  op in fm.orthogonalParameters
        flag = op.determinacy < 1.0e-6 ? "  (UNDETERMINED by these targets)" : ""
        println(io, "  orthogonal $(op.label): shift = $(op.value), determinacy = $(op.determinacy)$flag")
    end
    println(io, "  fit against $(length(fm.fitTargets)) measured level(s), rms deviation = $(fm.rmsDeviation) cm^-1")
    println(io, "  degrees of freedom: $(fm.degreesOfFreedom)")
    fm.degreesOfFreedom <= 0   &&
        println(io, "  WARNING: at or under the identifiability boundary -- read each group's own determinacy " *
                     "before trusting its value in isolation.")
    println(io, "  multiplet: $(fm.multiplet.name), $(length(fm.multiplet.levels)) level(s)")
end


"""
`struct  AdjustedLevels.FittedMultipletGroupedMultiBlock`
    ... the result of `AdjustedLevels.fitGroupedParametersMultiBlock`: ONE shared pair of GROUP factors, fit
        SIMULTANEOUSLY against measured levels drawn from SEVERAL symmetry blocks of the SAME basis -- the
        grouped-route sibling of `AdjustedLevels.FittedMultipletOrthogonalMultiBlock`, and the form this
        route is actually meant to be used in: at most two global parameters make pooling several blocks'
        targets the easiest way in this module to reach a genuinely positive `degreesOfFreedom`.

    + blocks                ::Array{FittedBlock,1}            ... one entry per symmetry block, in the order
                                                                     given to `fitGroupedParametersMultiBlock`.
    + groupedParameters      ::Array{GroupedParameter,1}        ... the fitted factors, GLOBAL across every
                                                                       block -- the same family in two blocks
                                                                       is the SAME parameter, not two
                                                                       independent ones.
    + orthogonalParameters   ::Array{OrthogonalParameter,1}      ... the orthogonal breakdown of the GLOBAL
                                                                        (at most two-dimensional) parameter
                                                                        space, read at the converged point,
                                                                        over EVERY target in EVERY block
                                                                        together.
    + rmsDeviation            ::Float64                           ... root-mean-square deviation over every
                                                                         target in every block together, in
                                                                         cm^-1.
    + degreesOfFreedom        ::Int64                             ... total targets (all blocks) minus the
                                                                         number of GLOBAL groups (at most two)
                                                                         -- see `FittedMultipletOrthogonal`'s
                                                                         own field of the same name.
    + source                  ::String                            ... a dated citation for the measured data.
"""
struct  FittedMultipletGroupedMultiBlock
    blocks                 ::Array{FittedBlock,1}
    groupedParameters        ::Array{GroupedParameter,1}
    orthogonalParameters      ::Array{OrthogonalParameter,1}
    rmsDeviation                ::Float64
    degreesOfFreedom             ::Int64
    source                       ::String
end


"""
`Base.show(io::IO, fm::AdjustedLevels.FittedMultipletGroupedMultiBlock)`
    ... ALWAYS prints the "fitted, not ab initio" disclosure first, the GLOBAL (at most two) group factors,
        the overfitting-guard warning when `degreesOfFreedom <= 0`, and then each block's own share of the
        fit. Nothing is returned.
"""
function  Base.show(io::IO, fm::FittedMultipletGroupedMultiBlock)
    println(io, "FittedMultipletGroupedMultiBlock -- FITTED, NOT AB INITIO.")
    println(io, "  source: $(fm.source)")
    for  gp in fm.groupedParameters
        println(io, "  grouped parameter [$(gp.groupLabel)] = $(gp.value)  (shared across every block below)   " *
                     "($(length(gp.memberClasses)) class(es): $(join(gp.memberClasses, ", ")))")
    end
    for  op in fm.orthogonalParameters
        flag = op.determinacy < 1.0e-6 ? "  (UNDETERMINED by these targets)" : ""
        println(io, "  orthogonal $(op.label): shift = $(op.value), determinacy = $(op.determinacy)$flag")
    end
    nTargetsTotal = sum( length(b.fitTargets) for b in fm.blocks )
    println(io, "  fit against $(nTargetsTotal) measured level(s) across $(length(fm.blocks)) block(s), " *
                 "rms deviation = $(fm.rmsDeviation) cm^-1")
    println(io, "  degrees of freedom: $(fm.degreesOfFreedom)")
    fm.degreesOfFreedom <= 0   &&
        println(io, "  WARNING: at or under the identifiability boundary -- read each group's own determinacy " *
                     "before trusting its value in isolation.")
    for  b in fm.blocks
        println(io, "  block $(b.sym): $(length(b.fitTargets)) target(s), $(length(b.multiplet.levels)) level(s)")
    end
end


"""
`AdjustedLevels.fitGroupedParametersMultiBlock(blocks::Array{Tuple{LevelSymmetry,Array{Tuple{Int64,Float64},1}},1},
                                                 basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, source::String;
                                                 determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                                 tolerance::Float64=1.0e-6)`
    ... the grouped-route sibling of `AdjustedLevels.fitOrthogonalParametersMultiBlock`: instead of one
        parameter PER Cowan class (routinely more classes than any realistic set of measured levels can
        determine once a correlating configuration is present), this fits exactly ONE parameter per Cowan
        FAMILY -- `AdjustedLevels.twoBodyParameterGroup`'s "Slater (F^k/G^k)" and "Configuration-interaction
        (R^k)" -- SIMULTANEOUSLY against targets drawn from SEVERAL symmetry blocks of the SAME `basis`, same
        calling convention as `fitOrthogonalParametersMultiBlock`. Each block's classes are folded into its
        (at most two) group matrices via `AdjustedLevels.groupParameterMatrices` BEFORE the global class list
        is built, so the Gauss-Newton solve below never sees more than two unknowns regardless of how rich
        the correlating configuration's own class list is. The iteration is otherwise the SAME Gauss-Newton +
        SVD scheme `fitOrthogonalParametersMultiBlock` uses (warm-started from the single uniform factor
        minimizing the combined rms over every block, the best iterate tracked rather than trusting the last
        one), just over the grouped matrices. A fm::AdjustedLevels.FittedMultipletGroupedMultiBlock is
        returned.
"""
function  fitGroupedParametersMultiBlock(blocks::Array{Tuple{LevelSymmetry,Array{Tuple{Int64,Float64},1}},1},
                                           basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, source::String;
                                           determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                           tolerance::Float64=1.0e-6)
    isempty(blocks)   &&   error("AdjustedLevels.fitGroupedParametersMultiBlock(): need at least one block.")
    auPerKayser = 1.0 / Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
    nBlocks = length(blocks)

    idxCsfB   = Vector{Array{Int64,1}}(undef, nBlocks)
    oneBodyB  = Vector{Array{Float64,2}}(undef, nBlocks)
    groupsB   = Vector{Array{String,1}}(undef, nBlocks)
    matricesB = Vector{Array{Array{Float64,2},1}}(undef, nBlocks)
    membersB  = Vector{Array{Array{String,1},1}}(undef, nBlocks)
    for  ib = 1:nBlocks
        sym, _ = blocks[ib]
        idxCsfB[ib], oneBodyB[ib], rawClasses, rawMatrices = buildParameterMatrices(sym, basis, nm, grid)
        groupsB[ib], matricesB[ib], membersB[ib] = groupParameterMatrices(rawClasses, rawMatrices)
    end

    globalGroups = sort(unique(vcat(groupsB...)))
    nParam = length(globalGroups)
    nParam == 0   &&
        error("AdjustedLevels.fitGroupedParametersMultiBlock(): no two-body (Coulomb) contributions were found in any block.")
    groupIndex = Dict{String,Int64}( g => k for (k,g) in enumerate(globalGroups) )

    globalMemberLists = [ String[] for _ = 1:nParam ]
    for  ib = 1:nBlocks
        for  (loc, g) in enumerate(groupsB[ib])   append!(globalMemberLists[groupIndex[g]], membersB[ib][loc])   end
    end
    globalMemberLists = [ sort(unique(ml)) for ml in globalMemberLists ]

    # Each block's own (at most two) group matrices, reindexed onto the GLOBAL group list; a group absent
    # from a block simply has no entry here (equivalent to a zero column for that block in the Jacobian).
    globalMatricesB = Vector{Dict{Int64,Array{Float64,2}}}(undef, nBlocks)
    for  ib = 1:nBlocks
        d = Dict{Int64,Array{Float64,2}}()
        for  (loc, g) in enumerate(groupsB[ib])   d[groupIndex[g]] = matricesB[ib][loc]   end
        globalMatricesB[ib] = d
    end

    nTargetsTotal = sum( length(blocks[ib][2]) for ib = 1:nBlocks )

    function  rmsAtUniformAll(c::Float64)
        s = 0.0
        for  ib = 1:nBlocks
            total = copy(oneBodyB[ib]);   for  mtx in matricesB[ib]   total .+= c .* mtx   end
            eig = Hamiltonian.diagonalizeCiMatrix(total, LevelSelection())
            rel = (sort(eig.values) .- sort(eig.values)[1]) ./ auPerKayser
            for  (idx, measured) in blocks[ib][2]   s = s + (rel[idx] - measured)^2   end
        end
        return( sqrt(s/nTargetsTotal) )
    end
    cUniform, _ = goldenSectionMinimize(rmsAtUniformAll, 0.5, 1.5)

    # Diagonalize EVERY block at a GIVEN parameter vector, once, returning everything a caller could need:
    # the per-block eigensystems (for the Jacobian) and the combined rms over every block's own targets --
    # used for the warm start, the per-iteration Jacobian build, AND the backtracking line search below, so
    # "evaluate the fit at this cVec" has exactly one implementation.
    function  evaluateAtVec(cv::Array{Float64,1})
        eigB = Vector{Any}(undef, nBlocks);   orderB = Vector{Any}(undef, nBlocks)
        relB = Vector{Array{Float64,1}}(undef, nBlocks);   s = 0.0
        for  ib = 1:nBlocks
            total = copy(oneBodyB[ib])
            for  (gIdx, mtx) in globalMatricesB[ib]   total .+= cv[gIdx] .* mtx   end
            eig   = Hamiltonian.diagonalizeCiMatrix(total, LevelSelection())
            order = sortperm(eig.values)
            relKayser = (eig.values[order] .- eig.values[order][1]) ./ auPerKayser
            eigB[ib] = eig;   orderB[ib] = order;   relB[ib] = relKayser
            for  (idx, measured) in blocks[ib][2]   s = s + (relKayser[idx] - measured)^2   end
        end
        return( (eigB, orderB, relB, sqrt(s/nTargetsTotal)) )
    end

    cVec = fill(cUniform, nParam)
    eigB, orderB, relB, rmsNow = evaluateAtVec(cVec)
    bestCVec = copy(cVec);   bestRms = rmsNow
    bestEigB = eigB;   bestOrderB = orderB;   bestRelB = relB

    for  iteration = 1:maxIterations
        J        = zeros(Float64, nTargetsTotal, nParam)
        residual = zeros(Float64, nTargetsTotal)
        row = 0
        for  ib = 1:nBlocks
            for  (idx, measured) in blocks[ib][2]
                row += 1
                vec = eigB[ib].vectors[orderB[ib][idx]]
                for  (gIdx, mtx) in globalMatricesB[ib]   J[row, gIdx] = (vec' * Symmetric(mtx) * vec) / auPerKayser   end
                residual[row] = measured - relB[ib][idx]
            end
        end

        F    = svd(J)
        sMax = isempty(F.S) ? 0.0 : maximum(F.S)
        deltaC = zeros(Float64, nParam)
        for  k = 1:length(F.S)
            F.S[k] <= determinacyThreshold * max(sMax, 1.0)   &&   continue
            deltaC .+= ( (F.U[:,k]' * residual) / F.S[k] ) .* F.V[:,k]
        end
        all(deltaC .== 0.0)   &&   break   # the Jacobian has nothing left to say; stop rather than loop idle

        # BACKTRACKING LINE SEARCH, added for this (at-most-two-parameter) route because an UNDAMPED full
        # Gauss-Newton step was MEASURED to overshoot and then OSCILLATE on a real correlated basis (Ti III
        # 3d^2+3p^4.3d^4, pooled J=0+2+4: rms 377.9 -> 964.4 -> 395.4 -> 790.2 cm^-1 ... never settling) --
        # worse than the per-class route's own documented "stall" near its 30%-per-class validity edge, this
        # genuinely diverged. Halving the step is the standard remedy and costs only re-diagonalizations
        # (already the dominant cost here, with at most two parameters to search over).
        step = 1.0;   trialVec = cVec .+ step .* deltaC
        trialEigB, trialOrderB, trialRelB, trialRms = evaluateAtVec(trialVec)
        nHalvings = 0
        while  trialRms >= rmsNow  &&  nHalvings < 10
            step = step / 2;   trialVec = cVec .+ step .* deltaC
            trialEigB, trialOrderB, trialRelB, trialRms = evaluateAtVec(trialVec);   nHalvings += 1
        end
        ( trialRms >= rmsNow   &&   nHalvings >= 10 )   &&   break   # no improving step found at all; stop here

        if  trialRms < bestRms
            bestRms = trialRms;   bestCVec = copy(trialVec)
            bestEigB = trialEigB;   bestOrderB = trialOrderB;   bestRelB = trialRelB
        end
        converged = abs(rmsNow - trialRms) < tolerance * max(rmsNow, 1.0)
        cVec = trialVec;   eigB = trialEigB;   orderB = trialOrderB;   relB = trialRelB;   rmsNow = trialRms
        converged   &&   break
    end

    cVec = bestCVec;   eigB = bestEigB;   orderB = bestOrderB;   relB = bestRelB
    rmsBest = bestRms
    groupedParameters = [ GroupedParameter(globalGroups[p], cVec[p], globalMemberLists[p]) for p = 1:nParam ]

    Jfinal = zeros(Float64, nTargetsTotal, nParam);   residualFinal = zeros(Float64, nTargetsTotal);   row = 0
    for  ib = 1:nBlocks
        for  (idx, measured) in blocks[ib][2]
            row += 1
            vec = eigB[ib].vectors[orderB[ib][idx]]
            for  (gIdx, mtx) in globalMatricesB[ib]   Jfinal[row, gIdx] = (vec' * Symmetric(mtx) * vec) / auPerKayser   end
            residualFinal[row] = measured - relB[ib][idx]
        end
    end
    Ffinal = svd(Jfinal)
    totalShift = cVec .- ones(nParam)
    orthogonalParameters = OrthogonalParameter[]
    for  k = 1:length(Ffinal.S)
        composition = [ (globalGroups[p], Ffinal.V[p,k]) for p = 1:nParam ]
        push!( orthogonalParameters, OrthogonalParameter("mode $k", dot(Ffinal.V[:,k], totalShift), composition, Ffinal.S[k]) )
    end

    fittedBlocks = FittedBlock[]
    for  ib = 1:nBlocks
        sym, targets = blocks[ib]
        levels = Level[]
        for  (pos, ev) in enumerate(orderB[ib])
            vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs)
                if  LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = eigB[ib].vectors[ev][ns]   end
            end
            push!( levels, Level(sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, pos, eigB[ib].values[orderB[ib][pos]], 0., true, basis, vector) )
        end
        fitTargets = [ FitTarget(idx, measured, relB[ib][idx]) for (idx, measured) in targets ]
        push!( fittedBlocks, FittedBlock(sym, Multiplet(string(sym)*"+ (fitted, grouped)", levels), fitTargets) )
    end

    return( FittedMultipletGroupedMultiBlock(fittedBlocks, groupedParameters, orthogonalParameters,
                                               rmsBest, nTargetsTotal - nParam, source) )
end


"""
`AdjustedLevels.fitGroupedParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                       targets::Array{Tuple{Int64,Float64},1}, source::String;
                                       determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                       tolerance::Float64=1.0e-6)`
    ... the single-block form of `AdjustedLevels.fitGroupedParametersMultiBlock`, for a caller who only has
        targets from ONE symmetry block -- a thin wrapper calling the multi-block route with exactly one
        block and unwrapping its result, so the Gauss-Newton mechanics exist in exactly one place. Same
        `targets` convention as `AdjustedLevels.fitScaleFactor`. A fm::AdjustedLevels.FittedMultipletGrouped
        is returned.
"""
function  fitGroupedParameters(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                 targets::Array{Tuple{Int64,Float64},1}, source::String;
                                 determinacyThreshold::Float64=1.0e-6, maxIterations::Int64=20,
                                 tolerance::Float64=1.0e-6)
    isempty(targets)   &&   error("AdjustedLevels.fitGroupedParameters(): need at least one target level to fit against.")
    fmMulti = fitGroupedParametersMultiBlock([(sym, targets)], basis, nm, grid, source;
                 determinacyThreshold=determinacyThreshold, maxIterations=maxIterations, tolerance=tolerance)
    b = fmMulti.blocks[1]

    return( FittedMultipletGrouped(b.multiplet, fmMulti.groupedParameters, fmMulti.orthogonalParameters,
                                     b.fitTargets, fmMulti.rmsDeviation, fmMulti.degreesOfFreedom, fmMulti.source) )
end


"""
`AdjustedLevels.mergeIntoMultiplet(fm::AdjustedLevels.FittedMultipletGrouped, abInitioMultiplet::ManyElectron.Multiplet)`
    ... the grouped-route form of `mergeIntoMultiplet(::FittedMultipletOrthogonal, ...)`: `fm`'s own symmetry
        block replaces the matching block of `abInitioMultiplet`, every other block is kept unchanged, and
        the result is re-sorted by energy. A mp::ManyElectron.Multiplet is returned.
"""
function  mergeIntoMultiplet(fm::FittedMultipletGrouped, abInitioMultiplet::ManyElectron.Multiplet)
    isempty(fm.multiplet.levels)   &&   error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletGrouped has no levels.")
    sym = LevelSymmetry(fm.multiplet.levels[1].J, fm.multiplet.levels[1].parity)
    otherLevels = [lev for lev in abInitioMultiplet.levels if LevelSymmetry(lev.J, lev.parity) != sym]

    length(otherLevels) + length(fm.multiplet.levels) != length(abInitioMultiplet.levels)   &&
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletGrouped's block has a different number of levels than " *
              "the matching block in abInitioMultiplet -- they must come from the same Basis.")

    allLevels = vcat(otherLevels, fm.multiplet.levels)
    mp = Multiplet(abInitioMultiplet.name * "+fitted", allLevels)

    return( Basics.sortByEnergy(mp) )
end


"""
`AdjustedLevels.mergeIntoMultiplet(fm::AdjustedLevels.FittedMultipletGroupedMultiBlock, abInitioMultiplet::ManyElectron.Multiplet)`
    ... the grouped-route form of `mergeIntoMultiplet(::FittedMultipletOrthogonalMultiBlock, ...)`: EVERY
        block carried by `fm` replaces its own matching block of `abInitioMultiplet` in one call, every
        OTHER symmetry block is kept unchanged, and the result is re-sorted by energy. A
        mp::ManyElectron.Multiplet is returned.
"""
function  mergeIntoMultiplet(fm::FittedMultipletGroupedMultiBlock, abInitioMultiplet::ManyElectron.Multiplet)
    isempty(fm.blocks)   &&   error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletGroupedMultiBlock has no blocks.")
    syms = [ b.sym for b in fm.blocks ]
    length(unique(syms)) == length(syms)   ||
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletGroupedMultiBlock carries the same symmetry block twice.")

    otherLevels  = [ lev for lev in abInitioMultiplet.levels if !(LevelSymmetry(lev.J, lev.parity) in syms) ]
    fittedLevels = vcat( [ b.multiplet.levels for b in fm.blocks ]... )

    length(otherLevels) + length(fittedLevels) != length(abInitioMultiplet.levels)   &&
        error("AdjustedLevels.mergeIntoMultiplet(): the FittedMultipletGroupedMultiBlock's blocks do not account for " *
              "exactly the levels of the matching blocks in abInitioMultiplet -- they must come from the same Basis.")

    mp = Multiplet(abInitioMultiplet.name * "+fitted", vcat(otherLevels, fittedLevels))

    return( Basics.sortByEnergy(mp) )
end
